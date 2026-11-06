"""Per-class, per-acquisition-stratum quantitative-readiness accounting."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from .metadata import recompute_v31_eligibility


LANES = (
    "current_positive",
    "current_negative",
    "v31_positive",
    "v31_negative",
    "positive_metadata_rescue_ceiling",
    "negative_metadata_rescue_ceiling",
)


def _present(value: Any) -> bool:
    text = str(value or "").strip()
    return bool(text) and not any(
        token in text.lower()
        for token in ("unreported", "unknown", "unresolved", "not_reported")
    )


def _leaf_class(record: Mapping[str, Any]) -> str | None:
    chemical_class = (record.get("compound") or {}).get("chemical_class") or {}
    leaf = chemical_class.get("leaf_class")
    return str(leaf).strip() if _present(leaf) else None


def _identity(record: Mapping[str, Any]) -> str | None:
    value = (record.get("compound") or {}).get("connectivity")
    return str(value).strip() if _present(value) else None


def _ce_stratum(program: Mapping[str, Any], bin_width: float) -> str:
    unit = str(program.get("unit") or "unreported")
    mode = str(program.get("mode") or "unreported")
    if mode == "fixed" and isinstance(program.get("center"), (int, float)):
        value = float(program["center"])
        lower = math.floor(value / bin_width) * bin_width
        upper = lower + bin_width
        return f"fixed:[{lower:g},{upper:g}):{unit}"
    if mode in {"stepped", "composite", "CES", "MSE", "bbCID", "all_ions"}:
        material = {
            "mode": mode,
            "unit": unit,
            "values": program.get("values"),
            "weights": program.get("weights"),
            "equation": program.get("equation"),
        }
        digest = hashlib.sha256(
            json.dumps(material, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()[:16]
        return f"{mode}:program:{digest}:{unit}"
    material = {
        "mode": mode,
        "unit": unit,
        "minimum": program.get("minimum"),
        "maximum": program.get("maximum"),
        "equation": program.get("equation"),
    }
    digest = hashlib.sha256(
        json.dumps(material, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()[:16]
    return f"{mode}:unpooled:{digest}:{unit}"


def acquisition_stratum(
    record: Mapping[str, Any], *, fixed_ce_bin_width: float
) -> dict[str, Any]:
    acquisition = record.get("acquisition") or {}
    ce = acquisition.get("collision_energy") or {}
    source = record.get("source") or {}
    license_data = source.get("license") or {}
    return {
        "instrument_domain_id": acquisition.get("instrument_domain_id")
        or "unresolved_instrument_domain",
        "instrument_vendor": acquisition.get("instrument_vendor"),
        "instrument_model": acquisition.get("instrument_model"),
        "activation": acquisition.get("activation"),
        "collision_cell_region": acquisition.get("collision_cell_region"),
        "polarity": acquisition.get("polarity"),
        "adduct": acquisition.get("adduct"),
        "charge": acquisition.get("charge"),
        "collision_energy_unit": ce.get("unit"),
        "collision_energy_mode": ce.get("mode"),
        "collision_energy_stratum": _ce_stratum(ce, fixed_ce_bin_width),
        "license_tier": license_data.get("tier") or "excluded",
    }


def _stratum_key(stratum: Mapping[str, Any]) -> str:
    return json.dumps(stratum, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _new_bucket() -> dict[str, Any]:
    return {
        "record_count": 0,
        "identities": set(),
        "formulas": set(),
        "sources": set(),
        "laboratories": set(),
        "acquisition_groups": set(),
        "mirror_clusters": set(),
    }


def _add(bucket: dict[str, Any], record: Mapping[str, Any]) -> None:
    bucket["record_count"] += 1
    compound = record.get("compound") or {}
    source = record.get("source") or {}
    independence = record.get("independence") or {}
    values = {
        "identities": compound.get("connectivity"),
        "formulas": compound.get("formula"),
        "sources": independence.get("source_origin_id") or source.get("source_id"),
        "laboratories": independence.get("laboratory_id")
        or source.get("laboratory_id"),
        "acquisition_groups": independence.get("acquisition_replication_group_id"),
        "mirror_clusters": independence.get("mirror_cluster_id"),
    }
    for name, value in values.items():
        if _present(value):
            bucket[name].add(str(value))


def _merge(buckets: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    result = _new_bucket()
    for bucket in buckets:
        result["record_count"] += int(bucket["record_count"])
        for name in (
            "identities",
            "formulas",
            "sources",
            "laboratories",
            "acquisition_groups",
            "mirror_clusters",
        ):
            result[name].update(bucket[name])
    return result


def _counts(bucket: Mapping[str, Any]) -> dict[str, int]:
    return {
        "record_count": int(bucket["record_count"]),
        "identity_count": len(bucket["identities"]),
        "formula_family_count": len(bucket["formulas"]),
        "source_origin_count": len(bucket["sources"]),
        "laboratory_count": len(bucket["laboratories"]),
        "acquisition_group_count": len(bucket["acquisition_groups"]),
        "mirror_cluster_count": len(bucket["mirror_clusters"]),
    }


def _gap(
    target: Mapping[str, int],
    background: Mapping[str, int],
    thresholds: Mapping[str, int],
    planning_required_per_arm: int,
) -> dict[str, Any]:
    required_target_identities = max(
        int(thresholds["minimum_target_identities"]), planning_required_per_arm
    )
    required_background_identities = max(
        int(thresholds["minimum_background_identities"]), planning_required_per_arm
    )
    requirements = {
        "target_identities": (
            target["identity_count"],
            required_target_identities,
        ),
        "background_identities": (
            background["identity_count"],
            required_background_identities,
        ),
        "target_formula_families": (
            target["formula_family_count"],
            int(thresholds["minimum_target_formula_families"]),
        ),
        "background_formula_families": (
            background["formula_family_count"],
            int(thresholds["minimum_background_formula_families"]),
        ),
        "target_sources": (
            target["source_origin_count"],
            int(thresholds["minimum_target_sources"]),
        ),
        "background_sources": (
            background["source_origin_count"],
            int(thresholds["minimum_background_sources"]),
        ),
        "target_laboratories": (
            target["laboratory_count"],
            int(thresholds["minimum_target_laboratories"]),
        ),
        "background_laboratories": (
            background["laboratory_count"],
            int(thresholds["minimum_background_laboratories"]),
        ),
    }
    shortfalls = {
        name: max(0, required - observed)
        for name, (observed, required) in requirements.items()
    }
    return {
        "required": {name: required for name, (_, required) in requirements.items()},
        "shortfalls": shortfalls,
        "passes_count_and_precision_planning_gate": not any(shortfalls.values()),
    }


def _identity_class_contract(
    records: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, str], dict[str, list[str]], int]:
    by_identity: dict[str, set[str]] = defaultdict(set)
    missing = 0
    for record in records:
        identity = _identity(record)
        leaf = _leaf_class(record)
        if identity is None or leaf is None:
            missing += 1
            continue
        by_identity[identity].add(leaf)
    valid = {
        identity: next(iter(classes))
        for identity, classes in by_identity.items()
        if len(classes) == 1
    }
    conflicts = {
        identity: sorted(classes)
        for identity, classes in by_identity.items()
        if len(classes) > 1
    }
    return valid, conflicts, missing


def build_eligibility_gap_matrix(
    records: Sequence[Mapping[str, Any]],
    *,
    thresholds: Mapping[str, int],
    planning_required_per_arm: int,
    fixed_ce_bin_width: float = 5.0,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    valid_classes, conflicts, missing_class_records = _identity_class_contract(records)
    cells: dict[tuple[str, str], dict[str, Any]] = {}
    strata: dict[str, dict[str, Any]] = {}
    overall_lane_records = Counter()
    overall_lane_identities: dict[str, set[str]] = defaultdict(set)
    derivations = Counter()
    excluded_conflict_records = 0
    for record in records:
        identity = _identity(record)
        leaf = valid_classes.get(identity or "")
        if leaf is None:
            if identity in conflicts:
                excluded_conflict_records += 1
            continue
        stratum = acquisition_stratum(record, fixed_ce_bin_width=fixed_ce_bin_width)
        key = _stratum_key(stratum)
        strata[key] = stratum
        cell = cells.setdefault(
            (key, leaf),
            {
                "all": _new_bucket(),
                "lanes": {lane: _new_bucket() for lane in LANES},
                "missing_fields": Counter(),
                "derivation_events": Counter(),
            },
        )
        _add(cell["all"], record)
        for field in (record.get("metadata_completeness") or {}).get(
            "missing_fields", []
        ):
            cell["missing_fields"][str(field)] += 1
        for event in (record.get("v31_metadata_derivation") or {}).get("events", []):
            if event.get("status") == "derived":
                event_id = str(event.get("derivation_id"))
                cell["derivation_events"][event_id] += 1
                derivations[event_id] += 1
        original = record.get("quantitative_eligibility") or {}
        recalculated = recompute_v31_eligibility(record)
        lane_flags = {
            "current_positive": original.get("positive_lr_eligible") is True,
            "current_negative": original.get("negative_lr_eligible") is True,
            "v31_positive": recalculated["positive_lr_eligible"],
            "v31_negative": recalculated["negative_lr_eligible"],
            "positive_metadata_rescue_ceiling": recalculated[
                "positive_ceiling_if_precursor_selection_is_verified"
            ],
            "negative_metadata_rescue_ceiling": recalculated[
                "negative_ceiling_if_isolation_is_verified"
            ],
        }
        for lane, enabled in lane_flags.items():
            if enabled:
                _add(cell["lanes"][lane], record)
                overall_lane_records[lane] += 1
                overall_lane_identities[lane].add(str(identity))

    by_stratum: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for (stratum_key, leaf), cell in cells.items():
        by_stratum[stratum_key][leaf] = cell
    rows: list[dict[str, Any]] = []
    pass_counts = Counter()
    for stratum_key, class_cells in sorted(by_stratum.items()):
        for leaf, cell in sorted(class_cells.items()):
            lane_results: dict[str, Any] = {}
            other_cells = [
                other for other_leaf, other in class_cells.items() if other_leaf != leaf
            ]
            for lane in LANES:
                target = _counts(cell["lanes"][lane])
                background = _counts(
                    _merge([other["lanes"][lane] for other in other_cells])
                )
                gap = _gap(target, background, thresholds, planning_required_per_arm)
                if gap["passes_count_and_precision_planning_gate"]:
                    pass_counts[lane] += 1
                lane_results[lane] = {
                    "target": target,
                    "background": background,
                    "gap": gap,
                }
            rows.append(
                {
                    "cell_id": hashlib.sha256(
                        f"{stratum_key}|{leaf}".encode("utf-8")
                    ).hexdigest()[:24],
                    "chemical_class": {
                        "ontology": "ChemOnt",
                        "level": "leaf_class",
                        "name": leaf,
                    },
                    "acquisition_stratum": strata[stratum_key],
                    "audit_record_counts": _counts(cell["all"]),
                    "metadata_missing_field_record_counts": dict(
                        sorted(cell["missing_fields"].items())
                    ),
                    "metadata_derivation_event_counts": dict(
                        sorted(cell["derivation_events"].items())
                    ),
                    "evidence_lanes": lane_results,
                    "runtime_effect": False,
                }
            )
    summary = {
        "record_count": len(records),
        "two_dimensional_identity_count": len(
            {identity for record in records if (identity := _identity(record))}
        ),
        "formula_family_count": len(
            {
                str((record.get("compound") or {}).get("formula"))
                for record in records
                if _present((record.get("compound") or {}).get("formula"))
            }
        ),
        "frozen_leaf_class_count": len(set(valid_classes.values())),
        "identity_class_conflict_count": len(conflicts),
        "identity_class_conflicts": conflicts,
        "records_missing_identity_or_leaf_class": missing_class_records,
        "records_excluded_for_identity_class_conflict": excluded_conflict_records,
        "cell_count": len(rows),
        "lane_record_counts": {
            lane: int(overall_lane_records[lane]) for lane in LANES
        },
        "lane_identity_counts": {
            lane: len(overall_lane_identities[lane]) for lane in LANES
        },
        "cells_passing_count_and_precision_planning_gate": {
            lane: int(pass_counts[lane]) for lane in LANES
        },
        "metadata_derivation_event_counts": dict(sorted(derivations.items())),
        "planning_required_identities_per_target_and_background_arm": planning_required_per_arm,
        "runtime_effect": False,
    }
    return rows, summary


def collision_energy_summary(
    records: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    groups: dict[tuple[Any, ...], dict[str, Any]] = {}
    for record in records:
        acquisition = record.get("acquisition") or {}
        ce = acquisition.get("collision_energy") or {}
        key = (
            acquisition.get("instrument_domain_id") or "unresolved_instrument_domain",
            acquisition.get("instrument_model"),
            acquisition.get("activation"),
            ce.get("unit"),
            ce.get("mode"),
            acquisition.get("polarity"),
        )
        group = groups.setdefault(
            key,
            {
                "record_count": 0,
                "identities": set(),
                "sources": set(),
                "centers": [],
                "programs": set(),
            },
        )
        group["record_count"] += 1
        identity = _identity(record)
        if identity:
            group["identities"].add(identity)
        source_id = (record.get("source") or {}).get("source_id")
        if _present(source_id):
            group["sources"].add(str(source_id))
        if isinstance(ce.get("center"), (int, float)):
            group["centers"].append(float(ce["center"]))
        program = {
            "values": ce.get("values"),
            "weights": ce.get("weights"),
            "equation": ce.get("equation"),
        }
        group["programs"].add(
            json.dumps(program, sort_keys=True, separators=(",", ":"))
        )
    rows = []
    for key, group in sorted(groups.items(), key=lambda item: str(item[0])):
        domain, model, activation, unit, mode, polarity = key
        centers = group["centers"]
        rows.append(
            {
                "instrument_domain_id": domain,
                "instrument_model": model,
                "activation": activation,
                "collision_energy_unit": unit,
                "collision_energy_mode": mode,
                "polarity": polarity,
                "record_count": group["record_count"],
                "identity_count": len(group["identities"]),
                "source_count": len(group["sources"]),
                "distinct_native_program_count": len(group["programs"]),
                "center_minimum": min(centers) if centers else None,
                "center_maximum": max(centers) if centers else None,
                "cross_unit_conversion_applied": False,
                "ce_response_model_fitted": False,
                "reason": "no frozen rule family has an eligible identity denominator",
            }
        )
    return {
        "group_count": len(rows),
        "cross_unit_conversion_applied": False,
        "ce_response_model_count": 0,
        "groups": rows,
    }


class _UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        root_left, root_right = self.find(left), self.find(right)
        if root_left != root_right:
            self.parent[root_right] = root_left


def _view_components(
    records: Sequence[Mapping[str, Any]], paths: Sequence[tuple[str, ...]]
) -> dict[str, Any]:
    if not records:
        return {
            "status": "infeasible_no_positive_lr_eligible_records",
            "record_count": 0,
            "connected_component_count": 0,
            "component_size_counts": {},
            "three_nonempty_zero_overlap_partitions_structurally_possible": False,
        }
    union = _UnionFind(len(records))
    first_for_value: dict[tuple[int, str], int] = {}
    for index, record in enumerate(records):
        for path_index, path in enumerate(paths):
            value: Any = record
            for part in path:
                value = value.get(part) if isinstance(value, Mapping) else None
            if not _present(value):
                continue
            key = (path_index, str(value))
            if key in first_for_value:
                union.union(index, first_for_value[key])
            else:
                first_for_value[key] = index
    sizes = Counter(union.find(index) for index in range(len(records)))
    size_counts = Counter(sizes.values())
    return {
        "status": "structural_check_complete",
        "record_count": len(records),
        "connected_component_count": len(sizes),
        "component_size_counts": {
            str(size): count for size, count in sorted(size_counts.items())
        },
        "three_nonempty_zero_overlap_partitions_structurally_possible": len(sizes)
        >= 3,
        "claim_limit": "this component count is necessary but not sufficient for a class-stratified discovery/calibration/test split",
    }


def split_isolation_readiness(
    records: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    eligible = [
        record
        for record in records
        if recompute_v31_eligibility(record)["positive_lr_eligible"]
    ]
    views = {
        "identity_and_formula": _view_components(
            eligible,
            (("compound", "connectivity"), ("compound", "formula")),
        ),
        "source_and_laboratory": _view_components(
            eligible,
            (
                ("independence", "source_origin_id"),
                ("independence", "laboratory_id"),
            ),
        ),
        "mirror_and_raw_acquisition": _view_components(
            eligible,
            (
                ("independence", "mirror_cluster_id"),
                ("independence", "acquisition_replication_group_id"),
            ),
        ),
    }
    return {
        "positive_lr_eligible_record_count": len(eligible),
        "views": views,
        "origin_isolated_split_attempted": False,
        "reason": (
            "no_positive_lr_eligible_records"
            if not eligible
            else "phase_3_1_reports_readiness_without_consuming_a_sealed_split"
        ),
    }
