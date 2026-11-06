"""Deterministic identity/formula/mirror-isolated benchmark splits."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from .mirror import UnionFind


SPLITS = ("discovery", "calibration", "test")


def collision_energy_bin(record: Mapping[str, Any]) -> str:
    ce = record["acquisition"]["collision_energy"]
    if ce["mode"] == "fixed":
        return f"{float(ce['value']):.8g}"
    return f"{float(ce['minimum']):.8g}-{float(ce['maximum']):.8g}"


def stratum_key(record: Mapping[str, Any]) -> tuple[Any, ...]:
    acquisition = record["acquisition"]
    ce = acquisition["collision_energy"]
    return (
        acquisition["instrument_domain_id"],
        acquisition["activation"],
        acquisition["collision_cell_region"],
        ce["unit"],
        ce["mode"],
        collision_energy_bin(record),
        acquisition["adduct"],
        acquisition["polarity"],
        int(acquisition["charge"]),
    )


def stratum_id(record: Mapping[str, Any]) -> str:
    material = json.dumps(stratum_key(record), separators=(",", ":"))
    return "stratum:" + hashlib.sha256(material.encode("utf-8")).hexdigest()[:20]


def _component_map(records: Sequence[Mapping[str, Any]]) -> dict[str, list[str]]:
    ids = [str(record["spectrum_id"]) for record in records]
    union = UnionFind(ids)
    first_seen: dict[tuple[str, str], str] = {}
    for record in records:
        spectrum_id = str(record["spectrum_id"])
        independence = record["independence"]
        keys = (
            ("identity", str(independence["two_dimensional_identity_id"])),
            ("formula", str(independence["formula_family_id"])),
            ("mirror", str(independence["mirror_cluster_id"])),
        )
        for key in keys:
            previous = first_seen.setdefault(key, spectrum_id)
            union.union(previous, spectrum_id)
    components: dict[str, list[str]] = defaultdict(list)
    for spectrum_id in ids:
        components[union.find(spectrum_id)].append(spectrum_id)
    return dict(components)


def _stable_tie(seed: str, value: str) -> str:
    return hashlib.sha256(f"{seed}|{value}".encode("utf-8")).hexdigest()


def build_isolated_split(
    records: Sequence[Mapping[str, Any]],
    *,
    ratios: Mapping[str, float] | None = None,
    seed: str = "v15-feature-fragment-benchmark-v2",
) -> dict[str, Any]:
    """Assign whole identity/formula/mirror components to three fixed splits."""
    if not records:
        raise ValueError("cannot split an empty corpus")
    ratios = dict(ratios or {"discovery": 0.60, "calibration": 0.20, "test": 0.20})
    if set(ratios) != set(SPLITS) or any(float(ratios[name]) <= 0 for name in SPLITS):
        raise ValueError("ratios must contain positive discovery/calibration/test values")
    total_ratio = sum(float(ratios[name]) for name in SPLITS)
    ratios = {name: float(ratios[name]) / total_ratio for name in SPLITS}

    record_by_id = {str(record["spectrum_id"]): record for record in records}
    if len(record_by_id) != len(records):
        raise ValueError("duplicate spectrum_id in split input")
    components = _component_map(records)
    total_records = len(records)
    total_by_stratum = Counter(stratum_id(record) for record in records)
    targets = {name: total_records * ratios[name] for name in SPLITS}
    target_by_stratum = {
        name: {key: count * ratios[name] for key, count in total_by_stratum.items()}
        for name in SPLITS
    }
    assigned_counts = Counter()
    assigned_by_stratum: dict[str, Counter[str]] = {
        name: Counter() for name in SPLITS
    }
    assignment: dict[str, dict[str, str]] = {}
    component_rows: list[dict[str, Any]] = []

    ordered_components = []
    for root, members in components.items():
        stratum_counts = Counter(stratum_id(record_by_id[member]) for member in members)
        stable_member_key = min(members)
        ordered_components.append(
            (len(members), _stable_tie(seed, stable_member_key), root, members, stratum_counts)
        )
    ordered_components.sort(key=lambda item: (-item[0], item[1]))

    for size, _, root, members, component_strata in ordered_components:
        choices: list[tuple[float, str]] = []
        for split_name in SPLITS:
            projected = assigned_counts[split_name] + size
            overall_fill = projected / max(targets[split_name], 1.0)
            stratum_fill = 0.0
            stratum_weight = 0
            for key, count in component_strata.items():
                target = target_by_stratum[split_name][key]
                after = assigned_by_stratum[split_name][key] + count
                stratum_fill += (after / max(target, 1.0)) * count
                stratum_weight += count
            mean_stratum_fill = stratum_fill / max(stratum_weight, 1)
            score = overall_fill + 0.10 * mean_stratum_fill
            choices.append((score, split_name))
        _, selected = min(choices, key=lambda item: (item[0], SPLITS.index(item[1])))
        component_material = "\n".join(sorted(members))
        component_id = "component:" + hashlib.sha256(
            component_material.encode("utf-8")
        ).hexdigest()[:20]
        for member in members:
            assignment[member] = {
                "split": selected,
                "component_id": component_id,
                "reason": "identity_formula_mirror_union_component",
            }
        assigned_counts[selected] += size
        assigned_by_stratum[selected].update(component_strata)
        component_rows.append(
            {
                "component_id": component_id,
                "split": selected,
                "record_count": size,
                "formula_count": len(
                    {
                        record_by_id[member]["compound"]["formula"] for member in members
                    }
                ),
                "identity_count": len(
                    {
                        record_by_id[member]["compound"]["connectivity"] for member in members
                    }
                ),
                "member_hash": hashlib.sha256(component_material.encode("utf-8")).hexdigest(),
            }
        )

    isolation = validate_split_isolation(records, assignment)
    seal_material = [
        [spectrum_id, row["split"], row["component_id"]]
        for spectrum_id, row in sorted(assignment.items())
    ]
    test_material = sorted(
        spectrum_id for spectrum_id, row in assignment.items() if row["split"] == "test"
    )
    return {
        "schema_version": "2.0.0",
        "seed": seed,
        "ratios": ratios,
        "isolation_dimensions": [
            "two_dimensional_identity",
            "normalized_formula_family",
            "mirror_cluster",
        ],
        "source_isolation": {
            "enforced": False,
            "reason": "single-source exact instrument domains prevent three-way source holdout",
        },
        "assignment": assignment,
        "components": component_rows,
        "split_record_counts": dict(sorted(assigned_counts.items())),
        "split_manifest_sha256": hashlib.sha256(
            json.dumps(seal_material, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "sealed_test_spectrum_ids_sha256": hashlib.sha256(
            json.dumps(test_material, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "isolation_audit": isolation,
        "split_feasible": isolation["zero_forbidden_overlap"]
        and all(assigned_counts[name] > 0 for name in SPLITS),
    }


def validate_split_isolation(
    records: Sequence[Mapping[str, Any]],
    assignment: Mapping[str, Mapping[str, str]],
) -> dict[str, Any]:
    dimensions = {
        "identity": lambda record: record["compound"]["connectivity"],
        "formula": lambda record: record["compound"]["formula"],
        "mirror": lambda record: record["independence"]["mirror_cluster_id"],
    }
    sets: dict[str, dict[str, set[str]]] = {
        dimension: {name: set() for name in SPLITS} for dimension in dimensions
    }
    sources = {name: set() for name in SPLITS}
    missing_assignments = []
    for record in records:
        spectrum_id = str(record["spectrum_id"])
        row = assignment.get(spectrum_id)
        if not row:
            missing_assignments.append(spectrum_id)
            continue
        split_name = str(row["split"])
        if split_name not in SPLITS:
            raise ValueError(f"invalid split name for {spectrum_id}: {split_name}")
        for dimension, getter in dimensions.items():
            sets[dimension][split_name].add(str(getter(record)))
        sources[split_name].add(str(record["source"]["source_id"]))
    overlaps: dict[str, dict[str, int]] = {}
    for dimension in dimensions:
        overlaps[dimension] = {}
        for i, left in enumerate(SPLITS):
            for right in SPLITS[i + 1 :]:
                overlaps[dimension][f"{left}__{right}"] = len(
                    sets[dimension][left] & sets[dimension][right]
                )
    forbidden_total = sum(
        value for dimension in overlaps.values() for value in dimension.values()
    )
    source_overlaps = {}
    for i, left in enumerate(SPLITS):
        for right in SPLITS[i + 1 :]:
            source_overlaps[f"{left}__{right}"] = len(sources[left] & sources[right])
    return {
        "missing_assignment_count": len(missing_assignments),
        "forbidden_overlap_counts": overlaps,
        "forbidden_overlap_total": forbidden_total,
        "zero_forbidden_overlap": forbidden_total == 0 and not missing_assignments,
        "source_overlap_counts_reported_not_forbidden": source_overlaps,
    }
