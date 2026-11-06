"""Discovery, freezing, and sealed evaluation for V3 fragment rules."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any

from matching.feature_fragment_benchmark.mirror import UnionFind
from matching.feature_fragment_benchmark.statistics import extract_features, feature_center

from .parsing import ce_program_key
from .statistics import (
    benjamini_hochberg,
    estimate_conditional_lr,
    evaluate_primary_gate,
    holm_adjust,
)


def _canonical(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def exact_stratum_key(record: Mapping[str, Any]) -> tuple[Any, ...]:
    acquisition = record["acquisition"]
    return (
        acquisition.get("instrument_domain_id"),
        acquisition.get("instrument_vendor"),
        acquisition.get("instrument_model"),
        acquisition.get("analyzer"),
        acquisition.get("activation"),
        acquisition.get("collision_cell_region"),
        ce_program_key(acquisition.get("collision_energy") or {}),
        acquisition.get("polarity"),
        acquisition.get("adduct"),
        acquisition.get("charge"),
    )


def exact_stratum_id(record: Mapping[str, Any]) -> str:
    return "stratum-v3:" + hashlib.sha256(_canonical(exact_stratum_key(record))).hexdigest()[:24]


def _majority(values: Sequence[bool]) -> bool | None:
    positives = sum(values)
    negatives = len(values) - positives
    if positives == negatives:
        return None
    return positives > negatives


def _source_equal_feature_states(
    records: Sequence[Mapping[str, Any]],
) -> dict[str, bool | None]:
    by_source: dict[str, dict[str, list[Mapping[str, Any]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    all_features: set[str] = set()
    feature_sets: dict[str, set[str]] = {}
    for record in records:
        spectrum_id = str(record["spectrum_id"])
        features = extract_features(record)
        feature_sets[spectrum_id] = features
        all_features.update(features)
        by_source[str(record["source"]["source_id"])][
            str(record["independence"]["acquisition_replication_group_id"])
        ].append(record)
    result: dict[str, bool | None] = {}
    for feature in all_features:
        source_votes: list[bool] = []
        source_indeterminate = False
        for acquisition_groups in by_source.values():
            acquisition_votes: list[bool] = []
            for group in acquisition_groups.values():
                scan_vote = _majority(
                    [feature in feature_sets[str(row["spectrum_id"])] for row in group]
                )
                if scan_vote is not None:
                    acquisition_votes.append(scan_vote)
            source_vote = _majority(acquisition_votes) if acquisition_votes else None
            if source_vote is None:
                source_indeterminate = True
            else:
                source_votes.append(source_vote)
        identity_vote = _majority(source_votes) if source_votes else None
        result[feature] = None if source_indeterminate and identity_vote is None else identity_vote
    return result


def build_identity_observations(
    records: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Collapse spectra to one source-equal identity vote per exact stratum."""
    grouped: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[
            (
                str(record["compound"]["connectivity"]),
                exact_stratum_id(record),
            )
        ].append(record)
    observations = []
    for (identity, stratum_id), members in sorted(grouped.items()):
        formulas = {str(row["compound"].get("formula")) for row in members}
        classes = Counter(
            str(row["compound"]["chemical_class"]["leaf_class"])
            for row in members
        )
        if len(formulas) != 1 or not classes:
            continue
        maximum = max(classes.values())
        class_winners = sorted(name for name, count in classes.items() if count == maximum)
        if len(class_winners) != 1:
            continue
        representative = min(members, key=lambda row: str(row["spectrum_id"]))
        source_ids = sorted({str(row["source"]["source_id"]) for row in members})
        laboratory_ids = sorted(
            {str(row["source"]["laboratory_id"]) for row in members}
        )
        source_cluster = (
            source_ids[0]
            if len(source_ids) == 1
            else "source_bundle:"
            + hashlib.sha256(_canonical(source_ids)).hexdigest()[:16]
        )
        laboratory_cluster = (
            laboratory_ids[0]
            if len(laboratory_ids) == 1
            else "laboratory_bundle:"
            + hashlib.sha256(_canonical(laboratory_ids)).hexdigest()[:16]
        )
        observations.append(
            {
                "observation_id": f"{identity}|{stratum_id}",
                "identity_id": identity,
                "formula_id": next(iter(formulas)),
                "chemical_class": class_winners[0],
                "chemical_class_provenance": representative["compound"][
                    "chemical_class"
                ],
                "stratum_id": stratum_id,
                "stratum_key": list(exact_stratum_key(representative)),
                "feature_states": _source_equal_feature_states(members),
                "source_ids": source_ids,
                "laboratory_ids": laboratory_ids,
                "source_id": source_cluster,
                "laboratory_id": laboratory_cluster,
                "record_count": len(members),
                "acquisition_replication_group_count": len(
                    {
                        row["independence"]["acquisition_replication_group_id"]
                        for row in members
                    }
                ),
                "positive_eligible": all(
                    row["quantitative_eligibility"]["positive_lr_eligible"]
                    for row in members
                ),
                "negative_eligible": all(
                    row["quantitative_eligibility"]["negative_lr_eligible"]
                    for row in members
                ),
                "representative_record": representative,
            }
        )
    return observations


def _rule_rows(
    observations: Sequence[Mapping[str, Any]], feature_id: str
) -> list[dict[str, Any]]:
    rows = []
    for observation in observations:
        state = observation["feature_states"].get(feature_id, False)
        determinate = state is not None
        rows.append(
            {
                "identity_id": observation["identity_id"],
                "formula_id": observation["formula_id"],
                "source_id": observation["source_id"],
                "laboratory_id": observation["laboratory_id"],
                "source_ids": observation["source_ids"],
                "laboratory_ids": observation["laboratory_ids"],
                "present": bool(state),
                "positive_eligible": bool(
                    determinate and observation["positive_eligible"]
                ),
                "negative_eligible": bool(
                    determinate and observation["negative_eligible"]
                ),
            }
        )
    return rows


def _rule_statistics(
    observations: Sequence[Mapping[str, Any]],
    *,
    chemical_class: str,
    feature_id: str,
    bootstrap_replicates: int,
    seed: int,
) -> dict[str, Any]:
    target_observations = [
        row for row in observations if row["chemical_class"] == chemical_class
    ]
    background_observations = [
        row for row in observations if row["chemical_class"] != chemical_class
    ]
    target = _rule_rows(target_observations, feature_id)
    background = _rule_rows(background_observations, feature_id)
    lr = estimate_conditional_lr(
        target,
        background,
        bootstrap_replicates=bootstrap_replicates,
        seed=seed,
    )
    eligible_target = [row for row in target if row["positive_eligible"]]
    eligible_background = [row for row in background if row["positive_eligible"]]
    return {
        **lr,
        "target_formula_families": len(
            {row["formula_id"] for row in eligible_target}
        ),
        "background_formula_families": len(
            {row["formula_id"] for row in eligible_background}
        ),
        "target_sources": sorted(
            {source for row in eligible_target for source in row["source_ids"]}
        ),
        "background_sources": sorted(
            {source for row in eligible_background for source in row["source_ids"]}
        ),
        "target_laboratories": sorted(
            {lab for row in eligible_target for lab in row["laboratory_ids"]}
        ),
        "background_laboratories": sorted(
            {lab for row in eligible_background for lab in row["laboratory_ids"]}
        ),
    }


def _candidate_id(chemical_class: str, stratum_id: str, feature_id: str) -> str:
    return "FFV3-" + hashlib.sha256(
        _canonical([chemical_class, stratum_id, feature_id])
    ).hexdigest()[:20].upper()


def _freeze_fragment_pathway_groups(
    rules: list[dict[str, Any]], observations: Sequence[Mapping[str, Any]], threshold: float
) -> None:
    by_context: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for rule in rules:
        by_context[(rule["chemical_class"]["name"], rule["stratum_id"])].append(rule)
    for (chemical_class, stratum_id), context_rules in by_context.items():
        ids = [rule["rule_id"] for rule in context_rules]
        union = UnionFind(ids)
        target = [
            row
            for row in observations
            if row["chemical_class"] == chemical_class
            and row["stratum_id"] == stratum_id
        ]
        presence = {
            rule["rule_id"]: {
                row["identity_id"]
                for row in target
                if row["feature_states"].get(rule["feature_id"]) is True
            }
            for rule in context_rules
        }
        for index, left in enumerate(context_rules):
            for right in context_rules[index + 1 :]:
                a, b = presence[left["rule_id"]], presence[right["rule_id"]]
                similarity = len(a & b) / len(a | b) if a | b else 0.0
                if similarity >= threshold:
                    union.union(left["rule_id"], right["rule_id"])
        groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        by_rule = {rule["rule_id"]: rule for rule in context_rules}
        for rule_id in ids:
            groups[union.find(rule_id)].append(by_rule[rule_id])
        for members in groups.values():
            group_id = "fragpath:" + hashlib.sha256(
                _canonical(sorted(rule["rule_id"] for rule in members))
            ).hexdigest()[:20]
            representative = min(
                members,
                key=lambda rule: (
                    rule["benjamini_hochberg_q"],
                    -rule["discovery_statistics"]["lr_positive"]["point"],
                    rule["rule_id"],
                ),
            )
            for rule in members:
                rule["fragment_pathway_group_id"] = group_id
                rule["frozen_group_representative"] = (
                    rule["rule_id"] == representative["rule_id"]
                )
                rule["fragment_group_policy"] = {
                    "discovery_target_identity_jaccard_minimum": threshold,
                    "representative_frozen_before_test": True,
                    "query_time_maximum_absolute_lr_selection": False,
                    "maximum_runtime_contributions_from_group": 1,
                }


def discover_rules(
    observations: Sequence[Mapping[str, Any]],
    *,
    minimum_target_identities: int = 20,
    minimum_background_identities: int = 20,
    minimum_target_present: int = 5,
    maximum_q_value: float = 0.05,
    minimum_lr_positive: float = 2.0,
    maximum_rules_per_class_stratum: int = 12,
    fragment_group_jaccard: float = 0.85,
) -> list[dict[str, Any]]:
    """Discover class-level rules from the discovery partition only."""
    by_stratum: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for observation in observations:
        by_stratum[observation["stratum_id"]].append(observation)
    candidates: list[dict[str, Any]] = []
    for stratum_id, rows in sorted(by_stratum.items()):
        classes = sorted({row["chemical_class"] for row in rows})
        for chemical_class in classes:
            target = [row for row in rows if row["chemical_class"] == chemical_class]
            background = [row for row in rows if row["chemical_class"] != chemical_class]
            if len(target) < minimum_target_identities or len(background) < minimum_background_identities:
                continue
            counts = Counter(
                feature
                for row in target
                for feature, state in row["feature_states"].items()
                if state is True
            )
            for feature_id, present in counts.items():
                if present < minimum_target_present:
                    continue
                stats = _rule_statistics(
                    rows,
                    chemical_class=chemical_class,
                    feature_id=feature_id,
                    bootstrap_replicates=0,
                    seed=151503,
                )
                if stats.get("status") != "estimated":
                    continue
                representative = target[0]["representative_record"]
                observation_type = feature_id.split(":", 1)[0]
                candidate = {
                    "schema_version": "3.0.0",
                    "rule_id": _candidate_id(chemical_class, stratum_id, feature_id),
                    "status": "discovery_candidate_research_only",
                    "chemical_class": {
                        "name": chemical_class,
                        "ontology": "ChemOnt",
                        "provenance": target[0]["chemical_class_provenance"],
                    },
                    "feature_id": feature_id,
                    "observation_type": observation_type,
                    "product_mz": feature_center(feature_id)
                    if observation_type == "product_ion"
                    else None,
                    "neutral_loss_da": feature_center(feature_id)
                    if observation_type == "neutral_loss"
                    else None,
                    "minimum_relative_intensity": 0.01,
                    "mass_tolerance": {
                        "absolute_da": 0.01,
                        "ppm": 20.0,
                        "policy": "larger_of_absolute_or_ppm",
                    },
                    "stratum_id": stratum_id,
                    "stratum_key": target[0]["stratum_key"],
                    "applicability": {
                        "instrument_domain_id": representative["acquisition"][
                            "instrument_domain_id"
                        ],
                        "instrument_vendor": representative["acquisition"][
                            "instrument_vendor"
                        ],
                        "instrument_model": representative["acquisition"][
                            "instrument_model"
                        ],
                        "analyzer": representative["acquisition"]["analyzer"],
                        "activation": representative["acquisition"]["activation"],
                        "collision_cell_region": representative["acquisition"][
                            "collision_cell_region"
                        ],
                        "collision_energy": representative["acquisition"][
                            "collision_energy"
                        ],
                        "polarity": representative["acquisition"]["polarity"],
                        "adduct": representative["acquisition"]["adduct"],
                        "charge": representative["acquisition"]["charge"],
                    },
                    "background_definition": "all eligible non-target identities in the identical exact acquisition stratum",
                    "discovery_statistics": stats,
                    "selection_p": stats["fisher_exact_two_sided_p"],
                    "benjamini_hochberg_q": None,
                    "evidence_stage_output": "class_level_msms_support_only",
                    "identity_upgrade_allowed": False,
                    "formula_upgrade_allowed": False,
                    "isomer_resolution_allowed": False,
                    "standard_confirmation_allowed": False,
                    "evidence_grade": "G2",
                    "runtime_influence": False,
                }
                candidates.append(candidate)
    adjusted = benjamini_hochberg([row["selection_p"] for row in candidates])
    for row, q_value in zip(candidates, adjusted):
        row["benjamini_hochberg_q"] = q_value
    retained = [
        row
        for row in candidates
        if row["benjamini_hochberg_q"] is not None
        and row["benjamini_hochberg_q"] <= maximum_q_value
        and row["discovery_statistics"]["lr_positive"]["point"]
        >= minimum_lr_positive
    ]
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in retained:
        grouped[(row["chemical_class"]["name"], row["stratum_id"])].append(row)
    limited = []
    for rows in grouped.values():
        rows.sort(
            key=lambda row: (
                row["benjamini_hochberg_q"],
                -row["discovery_statistics"]["lr_positive"]["point"],
                row["rule_id"],
            )
        )
        limited.extend(rows[:maximum_rules_per_class_stratum])
    _freeze_fragment_pathway_groups(limited, observations, fragment_group_jaccard)
    return sorted(limited, key=lambda row: row["rule_id"])


def calibrate_rules(
    rules: Sequence[Mapping[str, Any]], observations: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    """Freeze one discovery representative per fragment pathway group."""
    frozen = []
    for source in rules:
        rule = deepcopy(source)
        if not rule.get("frozen_group_representative"):
            rule["status"] = "fragment_pathway_nonrepresentative_not_tested"
            rule["selected_for_sealed_test"] = False
            frozen.append(rule)
            continue
        rows = [row for row in observations if row["stratum_id"] == rule["stratum_id"]]
        stats = _rule_statistics(
            rows,
            chemical_class=rule["chemical_class"]["name"],
            feature_id=rule["feature_id"],
            bootstrap_replicates=0,
            seed=151504,
        )
        failures = []
        if stats.get("status") != "estimated":
            failures.append("calibration_lr_unavailable")
        else:
            if stats["target_total"] < 10:
                failures.append("calibration_target_identities_lt_10")
            if stats["background_total"] < 10:
                failures.append("calibration_background_identities_lt_10")
            if stats["lr_positive"]["point"] <= 1.0:
                failures.append("calibration_lr_direction_failed")
        rule["calibration_statistics"] = stats
        rule["calibration_gate_failures"] = failures
        rule["selected_for_sealed_test"] = not failures
        rule["status"] = (
            "calibrated_frozen_research_only"
            if not failures
            else "calibration_failed_research_only"
        )
        frozen_material = {
            key: rule[key]
            for key in (
                "rule_id",
                "chemical_class",
                "feature_id",
                "observation_type",
                "product_mz",
                "neutral_loss_da",
                "minimum_relative_intensity",
                "mass_tolerance",
                "stratum_id",
                "stratum_key",
                "applicability",
                "background_definition",
                "fragment_pathway_group_id",
                "fragment_group_policy",
            )
        }
        rule["frozen_parameter_sha256"] = hashlib.sha256(
            _canonical(frozen_material)
        ).hexdigest()
        frozen.append(rule)
    return sorted(frozen, key=lambda row: row["rule_id"])


def evaluate_sealed_rules(
    calibrated: Sequence[Mapping[str, Any]],
    observations: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Evaluate the frozen rule family once and apply Holm correction."""
    evaluated = []
    for index, source in enumerate(calibrated):
        if not source.get("selected_for_sealed_test"):
            continue
        rule = deepcopy(source)
        rows = [row for row in observations if row["stratum_id"] == rule["stratum_id"]]
        stats = _rule_statistics(
            rows,
            chemical_class=rule["chemical_class"]["name"],
            feature_id=rule["feature_id"],
            bootstrap_replicates=2000,
            seed=151503 + index * 2,
        )
        rule["sealed_test_statistics"] = stats
        evaluated.append(rule)
    holm = holm_adjust(
        [
            row["sealed_test_statistics"].get("fisher_exact_two_sided_p")
            for row in evaluated
        ]
    )
    for rule, adjusted in zip(evaluated, holm):
        stats = rule["sealed_test_statistics"]
        gate = evaluate_primary_gate(
            stats,
            target_formula_count=int(stats.get("target_formula_families", 0)),
            background_formula_count=int(stats.get("background_formula_families", 0)),
            target_sources=stats.get("target_sources", []),
            background_sources=stats.get("background_sources", []),
            target_laboratories=stats.get("target_laboratories", []),
            background_laboratories=stats.get("background_laboratories", []),
            holm_adjusted_p=adjusted,
        )
        rule["holm_adjusted_p"] = adjusted
        rule["sealed_test_gate"] = gate
        rule["evidence_grade"] = gate["evidence_grade"]
        rule["status"] = "sealed_test_evaluated_research_only"
        rule["runtime_influence"] = False
    return evaluated
