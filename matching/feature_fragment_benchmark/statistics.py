"""Conditional feature-fragment discovery, calibration, and sealed testing.

All occurrence calculations collapse spectra to a two-dimensional identity
within one exact acquisition stratum.  Product-ion and neutral-loss bins are
never pooled across instruments, activation methods, CE units, CE modes,
adducts, polarities, or charges.  Scores are log likelihood ratios, not fixed
fragment points.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from scipy.stats import fisher_exact, spearmanr

from .mirror import UnionFind
from .split import collision_energy_bin, stratum_id, stratum_key


DEFAULT_FEATURE_BIN_DA = 0.01
DEFAULT_MIN_RELATIVE_INTENSITY = 0.01


def _feature_id(kind: str, bin_index: int, bin_width: float) -> str:
    return f"{kind}:{bin_index}:{bin_width:.8g}"


def _parse_feature_id(feature_id: str) -> tuple[str, int, float]:
    kind, raw_bin, raw_width = feature_id.split(":")
    return kind, int(raw_bin), float(raw_width)


def feature_center(feature_id: str) -> float:
    _, bin_index, width = _parse_feature_id(feature_id)
    return round(bin_index * width, 8)


def extract_features(
    record: Mapping[str, Any],
    *,
    bin_width_da: float = DEFAULT_FEATURE_BIN_DA,
    minimum_relative_intensity: float = DEFAULT_MIN_RELATIVE_INTENSITY,
) -> set[str]:
    precursor = float(record["acquisition"]["precursor_mz"])
    features: set[str] = set()
    for peak in record["peaks"]:
        mz, relative = float(peak[0]), float(peak[1])
        if relative < minimum_relative_intensity or mz < 40.0 or mz >= precursor - 1.0:
            continue
        product_bin = int(round(mz / bin_width_da))
        features.add(_feature_id("product_ion", product_bin, bin_width_da))
        loss = precursor - mz
        if 10.0 <= loss <= 500.0:
            loss_bin = int(round(loss / bin_width_da))
            features.add(_feature_id("neutral_loss", loss_bin, bin_width_da))
    return features


def _strict_absence_eligible(record: Mapping[str, Any], feature_id: str) -> bool:
    quality = record.get("quality") or {}
    required_true = (
        "centroid_peak_list",
        "positive_finite_peaks",
        "at_least_three_peaks",
        "collision_energy_verified",
        "instrument_domain_verified",
        "reported_mass_range_present",
        "reported_mass_range_contains_precursor",
    )
    if any(quality.get(field) is not True for field in required_true):
        return False
    mass_range = record["acquisition"].get("reported_mass_range")
    if not mass_range:
        return False
    kind, _, _ = _parse_feature_id(feature_id)
    expected = feature_center(feature_id)
    if kind == "neutral_loss":
        expected = float(record["acquisition"]["precursor_mz"]) - expected
    return float(mass_range["minimum"]) <= expected <= float(mass_range["maximum"])


def _exploratory_absence_eligible(record: Mapping[str, Any], feature_id: str) -> bool:
    quality = record.get("quality") or {}
    if any(
        quality.get(field) is not True
        for field in (
            "centroid_peak_list",
            "positive_finite_peaks",
            "at_least_three_peaks",
            "collision_energy_verified",
            "instrument_domain_verified",
        )
    ):
        return False
    kind, _, _ = _parse_feature_id(feature_id)
    expected = feature_center(feature_id)
    precursor = float(record["acquisition"]["precursor_mz"])
    if kind == "neutral_loss":
        expected = precursor - expected
    return 40.0 <= expected < precursor - 1.0


def build_identity_observations(
    records: Sequence[Mapping[str, Any]],
    *,
    bin_width_da: float = DEFAULT_FEATURE_BIN_DA,
    minimum_relative_intensity: float = DEFAULT_MIN_RELATIVE_INTENSITY,
) -> list[dict[str, Any]]:
    """Collapse mirror spectra and technical repeats to identity×stratum rows."""
    grouped: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[(record["compound"]["connectivity"], stratum_id(record))].append(record)
    observations: list[dict[str, Any]] = []
    for (identity, sid), members in grouped.items():
        by_mirror: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
        for member in members:
            by_mirror[member["independence"]["mirror_cluster_id"]].append(member)
        mirror_deduplicated_records = [
            sorted(bucket, key=lambda row: row["spectrum_id"])[0]
            for bucket in by_mirror.values()
        ]
        by_fragmentation_group: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
        for member in mirror_deduplicated_records:
            group_id = str(
                member["independence"].get("independent_fragmentation_group_id")
                or member["independence"]["mirror_cluster_id"]
            )
            by_fragmentation_group[group_id].append(member)
        independent_groups = [
            sorted(bucket, key=lambda row: row["spectrum_id"])
            for _, bucket in sorted(by_fragmentation_group.items())
        ]
        independent_records = [group[0] for group in independent_groups]
        feature_votes: Counter[str] = Counter()
        for group in independent_groups:
            within_group_votes: Counter[str] = Counter()
            for member in group:
                within_group_votes.update(
                    extract_features(
                        member,
                        bin_width_da=bin_width_da,
                        minimum_relative_intensity=minimum_relative_intensity,
                    )
                )
            within_group_required = len(group) // 2 + 1
            group_features = {
                feature
                for feature, votes in within_group_votes.items()
                if votes >= within_group_required
            }
            feature_votes.update(group_features)
        required_votes = len(independent_groups) // 2 + 1
        consensus = {
            feature for feature, votes in feature_votes.items() if votes >= required_votes
        }
        class_votes = Counter(
            member["compound"]["chemical_class"]["leaf_class"] for member in members
        )
        leaf_class = min(
            (name for name, count in class_votes.items() if count == max(class_votes.values())),
            default="",
        )
        representative = sorted(members, key=lambda row: row["spectrum_id"])[0]
        observations.append(
            {
                "observation_id": f"{identity}|{sid}",
                "identity": identity,
                "formula": representative["compound"]["formula"],
                "chemical_class": leaf_class,
                "chemical_class_provenance": representative["compound"]["chemical_class"],
                "stratum_id": sid,
                "stratum_key": list(stratum_key(representative)),
                "features": consensus,
                "record_count": len(members),
                "mirror_cluster_count": len(by_mirror),
                "mirror_deduplicated_record_count": len(mirror_deduplicated_records),
                "independent_fragmentation_group_count": len(independent_groups),
                "sources": sorted({member["source"]["source_id"] for member in members}),
                "laboratories": sorted(
                    {member["source"]["laboratory_id"] for member in members}
                ),
                "research_ids": sorted({member["source"]["research_id"] for member in members}),
                "precursor_mz": float(representative["acquisition"]["precursor_mz"]),
                "representative_record": representative,
                "independent_records": independent_records,
                "independent_groups": independent_groups,
            }
        )
    return observations


def _eligible(observation: Mapping[str, Any], feature_id: str, *, strict: bool) -> bool:
    predicate = _strict_absence_eligible if strict else _exploratory_absence_eligible
    groups = observation.get("independent_groups") or [
        [record] for record in observation["independent_records"]
    ]
    group_votes = 0
    for group in groups:
        within_group = sum(1 for record in group if predicate(record, feature_id))
        if within_group >= len(group) // 2 + 1:
            group_votes += 1
    return group_votes >= len(groups) // 2 + 1


def _katz_interval(
    events_target: int,
    total_target: int,
    events_background: int,
    total_background: int,
    *,
    z: float = 1.959963984540054,
) -> dict[str, float]:
    a = events_target + 0.5
    nt = total_target + 1.0
    c = events_background + 0.5
    nb = total_background + 1.0
    risk_target = a / nt
    risk_background = c / nb
    ratio = risk_target / risk_background
    variance = max(0.0, 1.0 / a - 1.0 / nt + 1.0 / c - 1.0 / nb)
    half = z * math.sqrt(variance)
    log_ratio = math.log(ratio)
    return {
        "point": ratio,
        "log_point": log_ratio,
        "lower_95": math.exp(log_ratio - half),
        "upper_95": math.exp(log_ratio + half),
        "log_half_width_95": half,
        "method": "Katz log-risk-ratio interval with Haldane-Anscombe correction",
    }


def _frequency_stats(
    target: Sequence[Mapping[str, Any]],
    background: Sequence[Mapping[str, Any]],
    feature_id: str,
    *,
    strict: bool,
) -> dict[str, Any]:
    eligible_target = [row for row in target if _eligible(row, feature_id, strict=strict)]
    eligible_background = [
        row for row in background if _eligible(row, feature_id, strict=strict)
    ]
    target_present = sum(feature_id in row["features"] for row in eligible_target)
    background_present = sum(feature_id in row["features"] for row in eligible_background)
    nt, nb = len(eligible_target), len(eligible_background)
    if nt and nb:
        lr_positive = _katz_interval(target_present, nt, background_present, nb)
        lr_negative = _katz_interval(nt - target_present, nt, nb - background_present, nb)
        p_value = float(
            fisher_exact(
                [
                    [target_present, nt - target_present],
                    [background_present, nb - background_present],
                ],
                alternative="two-sided",
            ).pvalue
        )
    else:
        lr_positive = None
        lr_negative = None
        p_value = None
    return {
        "absence_denominator_basis": (
            "reported_scan_mass_range"
            if strict
            else "exploratory_precursor_bounded_centroid_spectrum"
        ),
        "target_present_identities": target_present,
        "target_total_identities": nt,
        "background_present_identities": background_present,
        "background_total_identities": nb,
        "target_formula_families": len({row["formula"] for row in eligible_target}),
        "background_formula_families": len(
            {row["formula"] for row in eligible_background}
        ),
        "target_sources": sorted({x for row in eligible_target for x in row["sources"]}),
        "background_sources": sorted(
            {x for row in eligible_background for x in row["sources"]}
        ),
        "target_laboratories": sorted(
            {x for row in eligible_target for x in row["laboratories"]}
        ),
        "background_laboratories": sorted(
            {x for row in eligible_background for x in row["laboratories"]}
        ),
        "lr_positive": lr_positive,
        "lr_negative": lr_negative,
        "fisher_exact_two_sided_p": p_value,
    }


def _bh_adjust(rows: list[dict[str, Any]], p_key: str) -> None:
    valid = [(index, float(row[p_key])) for index, row in enumerate(rows) if row[p_key] is not None]
    valid.sort(key=lambda item: item[1])
    running = 1.0
    m = len(valid)
    adjusted: dict[int, float] = {}
    for rank_from_end, (index, p_value) in enumerate(reversed(valid), start=1):
        rank = m - rank_from_end + 1
        running = min(running, p_value * m / rank)
        adjusted[index] = running
    for index, row in enumerate(rows):
        row["benjamini_hochberg_q"] = adjusted.get(index)


def _candidate_id(
    chemical_class: str, sid: str, feature_id: str, minimum_relative_intensity: float
) -> str:
    material = json.dumps(
        [chemical_class, sid, feature_id, minimum_relative_intensity],
        ensure_ascii=False,
        separators=(",", ":"),
    )
    return "FFV2-" + hashlib.sha256(material.encode("utf-8")).hexdigest()[:20].upper()


def _assign_fragment_groups(
    rules: list[dict[str, Any]], observations: Sequence[Mapping[str, Any]], threshold: float
) -> None:
    by_context: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for rule in rules:
        by_context[(rule["chemical_class"]["name"], rule["stratum_id"])].append(rule)
    for (chemical_class, sid), context_rules in by_context.items():
        rule_ids = [rule["rule_id"] for rule in context_rules]
        union = UnionFind(rule_ids)
        target_rows = [
            row
            for row in observations
            if row["chemical_class"] == chemical_class and row["stratum_id"] == sid
        ]
        presence_sets = {
            rule["rule_id"]: {
                row["identity"] for row in target_rows if rule["feature_id"] in row["features"]
            }
            for rule in context_rules
        }
        for i, left in enumerate(context_rules):
            for right in context_rules[i + 1 :]:
                a = presence_sets[left["rule_id"]]
                b = presence_sets[right["rule_id"]]
                union_size = len(a | b)
                jaccard = len(a & b) / union_size if union_size else 0.0
                if jaccard >= threshold:
                    union.union(left["rule_id"], right["rule_id"])
        members_by_root: dict[str, list[str]] = defaultdict(list)
        for rule_id in rule_ids:
            members_by_root[union.find(rule_id)].append(rule_id)
        group_by_rule: dict[str, str] = {}
        for members in members_by_root.values():
            material = "\n".join(sorted(members))
            group_id = "fraggrp:" + hashlib.sha256(material.encode("utf-8")).hexdigest()[:16]
            for rule_id in members:
                group_by_rule[rule_id] = group_id
        for rule in context_rules:
            rule["independent_fragmentation_group"] = group_by_rule[rule["rule_id"]]
            rule["fragment_group_policy"] = {
                "construction": f"discovery target-identity Jaccard >= {threshold}",
                "runtime_combination": "at_most_one_frozen_group_contribution",
                "representative_selection": "frozen_before_test_not_dynamic_max_lr",
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
    minimum_relative_intensity: float = DEFAULT_MIN_RELATIVE_INTENSITY,
    fragment_group_jaccard: float = 0.85,
) -> list[dict[str, Any]]:
    """Discover candidate rules from discovery observations only."""
    by_stratum: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in observations:
        by_stratum[row["stratum_id"]].append(row)
    discovered: list[dict[str, Any]] = []
    for sid, stratum_rows in sorted(by_stratum.items()):
        by_class: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
        for row in stratum_rows:
            by_class[row["chemical_class"]].append(row)
        for chemical_class, target in sorted(by_class.items()):
            background = [row for row in stratum_rows if row["chemical_class"] != chemical_class]
            if len(target) < minimum_target_identities or len(background) < minimum_background_identities:
                continue
            feature_counts = Counter(feature for row in target for feature in row["features"])
            candidates: list[dict[str, Any]] = []
            for feature_id, present_count in feature_counts.items():
                if present_count < minimum_target_present:
                    continue
                exploratory = _frequency_stats(target, background, feature_id, strict=False)
                strict = _frequency_stats(target, background, feature_id, strict=True)
                if (
                    exploratory["target_total_identities"] < minimum_target_identities
                    or exploratory["background_total_identities"] < minimum_background_identities
                    or exploratory["lr_positive"] is None
                ):
                    continue
                center = feature_center(feature_id)
                kind, _, bin_width = _parse_feature_id(feature_id)
                representative_ce = target[0]["representative_record"]["acquisition"][
                    "collision_energy"
                ]
                rule_id = _candidate_id(
                    chemical_class, sid, feature_id, minimum_relative_intensity
                )
                candidates.append(
                    {
                        "rule_id": rule_id,
                        "status": "discovery_candidate_shadow",
                        "evidence_stage_output": "class_level_msms_support_only",
                        "identity_upgrade_allowed": False,
                        "formula_upgrade_allowed": False,
                        "isomer_resolution_allowed": False,
                        "standard_confirmation_allowed": False,
                        "chemical_class": {
                            "name": chemical_class,
                            "ontology": "ChemOnt",
                            "provenance": target[0]["chemical_class_provenance"],
                        },
                        "feature_id": feature_id,
                        "observation_type": kind,
                        "product_mz": center if kind == "product_ion" else None,
                        "neutral_loss_da": center if kind == "neutral_loss" else None,
                        "mass_tolerance": {
                            "absolute_da": max(bin_width, 0.01),
                            "ppm": 20.0,
                            "match_policy": "larger_of_absolute_or_ppm",
                        },
                        "minimum_relative_intensity": minimum_relative_intensity,
                        "stratum_id": sid,
                        "stratum_key": target[0]["stratum_key"],
                        "applicability": {
                            "instrument_domain_id": target[0]["stratum_key"][0],
                            "activation": target[0]["stratum_key"][1],
                            "activation_evidence": target[0]["representative_record"][
                                "acquisition"
                            ]["activation_evidence"],
                            "psi_beam_type_cid_confirmed": target[0][
                                "representative_record"
                            ]["acquisition"].get("psi_beam_type_cid_confirmed"),
                            "psi_hcd_cv_confirmed": target[0]["representative_record"][
                                "acquisition"
                            ].get("psi_hcd_cv_confirmed"),
                            "activation_architecture_interpretation": target[0][
                                "representative_record"
                            ]["acquisition"].get(
                                "activation_architecture_interpretation"
                            ),
                            "collision_cell_region": target[0]["stratum_key"][2],
                            "collision_energy_unit": target[0]["stratum_key"][3],
                            "collision_energy_mode": target[0]["stratum_key"][4],
                            "collision_energy_bin": target[0]["stratum_key"][5],
                            "collision_energy_minimum": representative_ce["minimum"],
                            "collision_energy_maximum": representative_ce["maximum"],
                            "collision_energy_unit_provenance": representative_ce[
                                "unit_provenance"
                            ],
                            "cross_vendor_or_cross_model_ce_equivalence_assumed": False,
                            "adduct": target[0]["stratum_key"][6],
                            "polarity": target[0]["stratum_key"][7],
                            "charge": target[0]["stratum_key"][8],
                        },
                        "background_definition": {
                            "included": "other ChemOnt leaf classes in the identical acquisition stratum",
                            "candidate_specific_mass_condition": "precursor permits the product-ion or neutral-loss observation",
                            "excluded": "target ChemOnt leaf class",
                        },
                        "discovery_statistics": {
                            "exploratory": exploratory,
                            "strict_reported_range": strict,
                        },
                        "selection_p": exploratory["fisher_exact_two_sided_p"],
                        "source_provenance": {
                            "data_origin": "MassBank public experimental spectra",
                            "source_ids": sorted(
                                {x for row in target + background for x in row["sources"]}
                            ),
                            "statistical_origin": "identity-collapsed discovery split",
                        },
                        "positive_evidence_limit": "class support within the exact acquisition stratum only",
                        "negative_evidence_limit": "absence is interpretable only under the saved denominator basis; isolation window is unreported",
                        "evidence_grade": "provisional_shadow",
                    }
                )
            _bh_adjust(candidates, "selection_p")
            retained = [
                rule
                for rule in candidates
                if rule["benjamini_hochberg_q"] is not None
                and rule["benjamini_hochberg_q"] <= maximum_q_value
                and rule["discovery_statistics"]["exploratory"]["lr_positive"]["point"]
                >= minimum_lr_positive
            ]
            retained.sort(
                key=lambda rule: (
                    rule["benjamini_hochberg_q"],
                    -rule["discovery_statistics"]["exploratory"]["lr_positive"]["point"],
                    rule["rule_id"],
                )
            )
            discovered.extend(retained[:maximum_rules_per_class_stratum])
    _assign_fragment_groups(discovered, observations, fragment_group_jaccard)
    return discovered


def _evaluate_rule(
    rule: Mapping[str, Any], observations: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    rows = [row for row in observations if row["stratum_id"] == rule["stratum_id"]]
    target = [row for row in rows if row["chemical_class"] == rule["chemical_class"]["name"]]
    background = [row for row in rows if row["chemical_class"] != rule["chemical_class"]["name"]]
    return {
        "exploratory": _frequency_stats(target, background, rule["feature_id"], strict=False),
        "strict_reported_range": _frequency_stats(
            target, background, rule["feature_id"], strict=True
        ),
    }


def calibrate_rules(
    discovered_rules: Sequence[Mapping[str, Any]],
    calibration_observations: Sequence[Mapping[str, Any]],
    *,
    minimum_target_identities: int = 10,
    minimum_background_identities: int = 10,
) -> list[dict[str, Any]]:
    """Evaluate fixed discovery candidates and freeze the test rule set."""
    calibrated: list[dict[str, Any]] = []
    for source_rule in discovered_rules:
        rule = json.loads(json.dumps(source_rule))
        statistics = _evaluate_rule(rule, calibration_observations)
        exploratory = statistics["exploratory"]
        unmet = []
        if exploratory["target_total_identities"] < minimum_target_identities:
            unmet.append("calibration_target_identities")
        if exploratory["background_total_identities"] < minimum_background_identities:
            unmet.append("calibration_background_identities")
        if exploratory["lr_positive"] is None or exploratory["lr_positive"]["point"] <= 1.0:
            unmet.append("calibration_lr_direction")
        rule["calibration_statistics"] = statistics
        rule["calibration_gate_unmet"] = unmet
        rule["selected_for_sealed_test"] = not unmet
        rule["status"] = (
            "calibrated_frozen_shadow" if not unmet else "calibration_failed_shadow"
        )
        rule["frozen_parameter_sha256"] = hashlib.sha256(
            json.dumps(
                {
                    key: rule[key]
                    for key in (
                        "rule_id",
                        "feature_id",
                        "observation_type",
                        "product_mz",
                        "neutral_loss_da",
                        "mass_tolerance",
                        "minimum_relative_intensity",
                        "stratum_id",
                        "applicability",
                        "chemical_class",
                        "background_definition",
                        "independent_fragmentation_group",
                        "fragment_group_policy",
                        "positive_evidence_limit",
                        "negative_evidence_limit",
                    )
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()
        calibrated.append(rule)
    return calibrated


def _empirical_gate(statistics: Mapping[str, Any]) -> tuple[str, list[str]]:
    strict = statistics["strict_reported_range"]
    unmet = []
    if strict["target_total_identities"] < 30:
        unmet.append("test_target_identities_lt_30")
    if strict["background_total_identities"] < 30:
        unmet.append("test_background_identities_lt_30")
    if strict["target_formula_families"] < 10:
        unmet.append("test_target_formulas_lt_10")
    if strict["background_formula_families"] < 10:
        unmet.append("test_background_formulas_lt_10")
    if len(strict["target_sources"]) < 3 or len(strict["background_sources"]) < 3:
        unmet.append("test_sources_lt_3")
    if len(strict["target_laboratories"]) < 2 or len(strict["background_laboratories"]) < 2:
        unmet.append("test_laboratories_lt_2")
    lr = strict["lr_positive"]
    if lr is None or lr["lower_95"] <= 1.0:
        unmet.append("lr_positive_interval_not_above_1")
    if lr is None or lr["log_half_width_95"] > math.log(2.0):
        unmet.append("lr_positive_interval_too_wide")
    return ("empirically_validated_exact_domain" if not unmet else "provisional_shadow", unmet)


def test_frozen_rules(
    calibrated_rules: Sequence[Mapping[str, Any]],
    test_observations: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Evaluate calibration-selected rules once on a sealed test split."""
    tested: list[dict[str, Any]] = []
    for source_rule in calibrated_rules:
        if not source_rule.get("selected_for_sealed_test"):
            continue
        rule = json.loads(json.dumps(source_rule))
        statistics = _evaluate_rule(rule, test_observations)
        grade, unmet = _empirical_gate(statistics)
        rule["sealed_test_statistics"] = statistics
        rule["evidence_grade"] = grade
        rule["empirical_validation_gate_unmet"] = unmet
        rule["status"] = "sealed_test_evaluated_shadow"
        rule["conditional_evidence_score"] = {
            "present": (
                statistics["exploratory"]["lr_positive"] or {}
            ).get("log_point"),
            "absent": (
                statistics["exploratory"]["lr_negative"] or {}
            ).get("log_point"),
            "scale": "natural_log_likelihood_ratio",
            "fixed_points_used": False,
        }
        tested.append(rule)
    return tested


def evaluate_ce_dependence(
    tested_rules: Sequence[Mapping[str, Any]],
    observations: Sequence[Mapping[str, Any]],
    *,
    minimum_shared_identities: int = 20,
    minimum_ce_settings: int = 3,
) -> list[dict[str, Any]]:
    """Evaluate same-unit, same-mode CE series for frozen rule observations."""
    families: dict[tuple[Any, ...], list[Mapping[str, Any]]] = defaultdict(list)
    for rule in tested_rules:
        app = rule["applicability"]
        if app["collision_energy_mode"] != "fixed":
            continue
        key = (
            rule["chemical_class"]["name"],
            rule["feature_id"],
            app["instrument_domain_id"],
            app["activation"],
            app["collision_cell_region"],
            app["collision_energy_unit"],
            app["collision_energy_mode"],
            app["adduct"],
            app["polarity"],
            app["charge"],
        )
        families[key].append(rule)
    results = []
    rows_by_stratum: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in observations:
        rows_by_stratum[row["stratum_id"]].append(row)
    for key, rules in families.items():
        unique_by_ce = {float(rule["applicability"]["collision_energy_bin"]): rule for rule in rules}
        if len(unique_by_ce) < minimum_ce_settings:
            continue
        identities_by_ce = {}
        for ce_value, rule in unique_by_ce.items():
            identities_by_ce[ce_value] = {
                row["identity"]
                for row in rows_by_stratum[rule["stratum_id"]]
                if row["chemical_class"] == rule["chemical_class"]["name"]
            }
        shared = set.intersection(*identities_by_ce.values()) if identities_by_ce else set()
        rates = []
        for ce_value, rule in sorted(unique_by_ce.items()):
            target_rows = {
                row["identity"]: row
                for row in rows_by_stratum[rule["stratum_id"]]
                if row["chemical_class"] == rule["chemical_class"]["name"]
                and row["identity"] in shared
            }
            present = sum(rule["feature_id"] in row["features"] for row in target_rows.values())
            rates.append(
                {
                    "ce_value": ce_value,
                    "present_identities": present,
                    "shared_identity_total": len(shared),
                    "frequency": present / len(shared) if shared else None,
                }
            )
        if len(shared) >= minimum_shared_identities:
            correlation = spearmanr(
                [row["ce_value"] for row in rates],
                [row["frequency"] for row in rates],
            )
            rho, p_value = float(correlation.statistic), float(correlation.pvalue)
            status = "ce_series_evaluated"
        else:
            rho, p_value = None, None
            status = "insufficient_shared_identities"
        results.append(
            {
                "series_id": "CE-" + hashlib.sha256(
                    json.dumps(key, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
                ).hexdigest()[:20],
                "context": list(key),
                "ce_setting_count": len(unique_by_ce),
                "shared_identity_count": len(shared),
                "minimum_shared_identities": minimum_shared_identities,
                "rates": rates,
                "spearman_rho": rho,
                "spearman_p": p_value,
                "status": status,
                "cross_unit_conversion_applied": False,
            }
        )
    return results
