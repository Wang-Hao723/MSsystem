"""Source- and laboratory-isolated deterministic benchmark splitting."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from matching.feature_fragment_benchmark.mirror import UnionFind


def _canonical(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _token(dimension: str, value: Any) -> str:
    return f"{dimension}:{value}"


def _components(
    records: Sequence[Mapping[str, Any]], dimensions: Sequence[str]
) -> list[list[str]]:
    ids = [str(row["spectrum_id"]) for row in records]
    union = UnionFind(ids)
    seen: dict[str, str] = {}
    for row in records:
        spectrum_id = str(row["spectrum_id"])
        source = row["source"]
        independence = row["independence"]
        values = {
            "identity": independence.get("two_dimensional_identity_id"),
            "formula": independence.get("formula_family_id"),
            "mirror": independence.get("mirror_cluster_id"),
            "acquisition": independence.get("acquisition_replication_group_id"),
            "source": source.get("source_id"),
            "laboratory": source.get("laboratory_id"),
        }
        for dimension in dimensions:
            value = values.get(dimension)
            if value is None:
                continue
            key = _token(dimension, value)
            if key in seen:
                union.union(spectrum_id, seen[key])
            else:
                seen[key] = spectrum_id
    groups: dict[str, list[str]] = defaultdict(list)
    for spectrum_id in ids:
        groups[union.find(spectrum_id)].append(spectrum_id)
    return sorted((sorted(group) for group in groups.values()), key=lambda row: row[0])


def _select_holdout(
    components: Sequence[Sequence[str]], *, target_count: int, seed: str
) -> set[str]:
    if len(components) < 2 or target_count <= 0:
        return set()
    ranked = sorted(
        components,
        key=lambda group: (
            hashlib.sha256(_canonical([seed, list(group)])).hexdigest(),
            group[0],
        ),
    )
    selected: list[Sequence[str]] = []
    count = 0
    for index, group in enumerate(ranked):
        if index == len(ranked) - 1 and not selected:
            break
        before = abs(target_count - count)
        after = abs(target_count - (count + len(group)))
        if after <= before or not selected:
            selected.append(group)
            count += len(group)
        if count >= target_count:
            break
    if len(selected) == len(ranked):
        selected.pop()
    return {spectrum_id for group in selected for spectrum_id in group}


def _sets_by_split(
    records: Sequence[Mapping[str, Any]], assignment: Mapping[str, str]
) -> dict[str, dict[str, set[str]]]:
    result: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in records:
        split = assignment[str(row["spectrum_id"])]
        source = row["source"]
        independence = row["independence"]
        for dimension, value in {
            "identity": independence.get("two_dimensional_identity_id"),
            "formula": independence.get("formula_family_id"),
            "mirror": independence.get("mirror_cluster_id"),
            "acquisition": independence.get("acquisition_replication_group_id"),
            "source": source.get("source_id"),
            "laboratory": source.get("laboratory_id"),
        }.items():
            if value is not None:
                result[split][dimension].add(str(value))
    return result


def _overlap_audit(
    records: Sequence[Mapping[str, Any]], assignment: Mapping[str, str]
) -> dict[str, Any]:
    values = _sets_by_split(records, assignment)
    pairs = (("discovery", "calibration"), ("discovery", "test"), ("calibration", "test"))
    dimensions = ("identity", "formula", "mirror", "acquisition", "source", "laboratory")
    counts: dict[str, dict[str, int]] = {}
    for dimension in dimensions:
        counts[dimension] = {
            f"{left}__{right}": len(values[left][dimension] & values[right][dimension])
            for left, right in pairs
        }
    forbidden = []
    for dimension in ("identity", "formula", "mirror", "acquisition"):
        for pair, count in counts[dimension].items():
            if count:
                forbidden.append(f"{dimension}:{pair}:{count}")
    for dimension in ("source", "laboratory"):
        for pair in ("discovery__test", "calibration__test"):
            if counts[dimension][pair]:
                forbidden.append(f"{dimension}:{pair}:{counts[dimension][pair]}")
    return {
        "overlap_counts": counts,
        "forbidden_overlap_failures": forbidden,
        "zero_forbidden_overlap": not forbidden,
    }


def build_origin_isolated_split(
    records: Sequence[Mapping[str, Any]],
    *,
    seed: str = "v15-feature-fragment-benchmark-v3-origin-isolated",
    test_fraction: float = 0.2,
    calibration_fraction_within_development: float = 0.25,
) -> dict[str, Any]:
    """Create a sealed source/lab holdout and isolated dev partitions."""
    if not records:
        return {
            "schema_version": "3.0.0",
            "split_feasible": False,
            "failure_reasons": ["no_quantitative_records"],
            "assignment": {},
            "split_record_counts": {},
            "isolation_audit": {
                "overlap_counts": {},
                "forbidden_overlap_failures": [],
                "zero_forbidden_overlap": False,
            },
        }
    ids = {str(row["spectrum_id"]) for row in records}
    test_components = _components(
        records,
        ("identity", "formula", "mirror", "acquisition", "source", "laboratory"),
    )
    test_ids = _select_holdout(
        test_components,
        target_count=max(1, round(len(records) * test_fraction)),
        seed=seed + ":test",
    )
    development = [row for row in records if str(row["spectrum_id"]) not in test_ids]
    dev_components = _components(
        development, ("identity", "formula", "mirror", "acquisition")
    )
    calibration_ids = _select_holdout(
        dev_components,
        target_count=max(1, round(len(development) * calibration_fraction_within_development))
        if development
        else 0,
        seed=seed + ":calibration",
    )
    assignment = {
        spectrum_id: (
            "test"
            if spectrum_id in test_ids
            else "calibration"
            if spectrum_id in calibration_ids
            else "discovery"
        )
        for spectrum_id in sorted(ids)
    }
    counts = Counter(assignment.values())
    audit = _overlap_audit(records, assignment)
    failures = []
    if len(test_components) < 2 or not test_ids or len(test_ids) == len(records):
        failures.append("source_laboratory_isolated_test_partition_infeasible")
    if len(dev_components) < 2 or not calibration_ids:
        failures.append("isolated_discovery_calibration_partition_infeasible")
    if any(counts.get(name, 0) == 0 for name in ("discovery", "calibration", "test")):
        failures.append("one_or_more_required_splits_empty")
    failures.extend(audit["forbidden_overlap_failures"])
    seal = hashlib.sha256(
        _canonical([[spectrum_id, assignment[spectrum_id]] for spectrum_id in sorted(assignment)])
    ).hexdigest()
    test_seal = hashlib.sha256(_canonical(sorted(test_ids))).hexdigest()
    return {
        "schema_version": "3.0.0",
        "seed": seed,
        "split_feasible": not failures,
        "failure_reasons": sorted(set(failures)),
        "assignment": assignment,
        "split_record_counts": dict(sorted(counts.items())),
        "development_test_component_count": len(test_components),
        "development_component_count": len(dev_components),
        "split_manifest_sha256": seal,
        "sealed_test_spectrum_ids_sha256": test_seal,
        "isolation_audit": audit,
    }
