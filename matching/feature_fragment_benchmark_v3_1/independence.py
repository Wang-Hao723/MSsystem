"""Multi-key mirror isolation for the V3.1 research amendment."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from matching.feature_fragment_benchmark.mirror import (
    UnionFind,
    greedy_peak_cosine,
    sparse_peak_signature,
)


def canonical_usi_core(value: Any) -> str | None:
    """Return the five-part spectrum locator, excluding interpretation text."""
    raw = str(value or "").strip()
    parts = raw.split(":", 5)
    if len(parts) < 5 or parts[0].lower() != "mzspec":
        return None
    collection, run, index_type, index_value = parts[1:5]
    if not collection or not run or index_type not in {"scan", "index", "nativeId"} or not index_value:
        return None
    return ":".join(("mzspec", collection, run, index_type, index_value))


def origin_raw_locator(record: Mapping[str, Any]) -> str | None:
    source = record.get("source") or {}
    raw = record.get("raw_acquisition") or {}
    values = (
        source.get("origin_dataset_id"),
        raw.get("raw_file_id"),
        raw.get("scan_id"),
    )
    if not all(str(value or "").strip() for value in values):
        return None
    payload = json.dumps(values, ensure_ascii=False, separators=(",", ":"))
    return "origin:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _link_group(
    rows: Sequence[Mapping[str, Any]],
    union: UnionFind,
    edge_type: str,
    edge_sets: dict[str, set[tuple[str, str]]],
) -> None:
    if len(rows) < 2:
        return
    anchor = str(rows[0]["spectrum_id"])
    for row in rows[1:]:
        other = str(row["spectrum_id"])
        edge = tuple(sorted((anchor, other)))
        edge_sets[edge_type].add(edge)
        union.union(anchor, other)


def cluster_mirrors_v31(
    records: Sequence[Mapping[str, Any]],
    *,
    near_calibration: Mapping[str, Any] | None = None,
    mz_tolerance: float = 0.01,
) -> dict[str, Any]:
    """Union separate exact lineage/content edges, then calibrated near edges."""
    ids = [str(row["spectrum_id"]) for row in records]
    union = UnionFind(ids)
    groups: dict[str, dict[str, list[Mapping[str, Any]]]] = {
        "canonical_usi_core": defaultdict(list),
        "origin_raw_locator": defaultdict(list),
        "exact_peak_hash": defaultdict(list),
    }
    for row in records:
        source = row.get("source") or {}
        usi = canonical_usi_core(source.get("usi"))
        origin = origin_raw_locator(row)
        peak_hash = str((row.get("independence") or {}).get("exact_peak_hash") or "")
        if usi:
            groups["canonical_usi_core"][usi].append(row)
        if origin:
            groups["origin_raw_locator"][origin].append(row)
        if peak_hash:
            groups["exact_peak_hash"][peak_hash].append(row)
    edge_sets: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for edge_type, values in groups.items():
        for rows in values.values():
            _link_group(rows, union, edge_type, edge_sets)

    calibration = dict(near_calibration or {})
    near_allowed = bool(
        calibration.get("status") == "calibrated"
        and calibration.get("near_mirror_links_allowed") is True
        and calibration.get("threshold") is not None
        and calibration.get("minimum_bilateral_coverage") is not None
    )
    near_candidates: list[dict[str, Any]] = []
    if near_allowed:
        signature_groups: dict[tuple[Any, ...], list[Mapping[str, Any]]] = defaultdict(list)
        for row in records:
            if "peaks" not in row:
                raise ValueError("peaks_required_for_calibrated_near_mirror_evaluation")
            acquisition = row.get("acquisition") or {}
            compound = row.get("compound") or {}
            signature_groups[
                (
                    compound.get("connectivity"),
                    acquisition.get("polarity"),
                    acquisition.get("adduct"),
                    sparse_peak_signature(row["peaks"]),
                )
            ].append(row)
        for rows in signature_groups.values():
            for index, left in enumerate(rows):
                for right in rows[index + 1 :]:
                    left_id = str(left["spectrum_id"])
                    right_id = str(right["spectrum_id"])
                    if union.find(left_id) == union.find(right_id):
                        continue
                    cosine, coverage_left, coverage_right = greedy_peak_cosine(
                        left["peaks"], right["peaks"], mz_tolerance=mz_tolerance
                    )
                    linked = bool(
                        cosine >= float(calibration["threshold"])
                        and min(coverage_left, coverage_right)
                        >= float(calibration["minimum_bilateral_coverage"])
                    )
                    near_candidates.append(
                        {
                            "left": left_id,
                            "right": right_id,
                            "cosine": cosine,
                            "coverage_left": coverage_left,
                            "coverage_right": coverage_right,
                            "linked": linked,
                        }
                    )
                    if linked:
                        edge = tuple(sorted((left_id, right_id)))
                        edge_sets["calibrated_near_mirror"].add(edge)
                        union.union(left_id, right_id)

    components: dict[str, list[str]] = defaultdict(list)
    for spectrum_id in ids:
        components[union.find(spectrum_id)].append(spectrum_id)
    ordered = sorted((sorted(members) for members in components.values()), key=lambda row: row[0])
    assignment: dict[str, str] = {}
    clusters: list[dict[str, Any]] = []
    all_exact_edges = set().union(
        edge_sets["canonical_usi_core"],
        edge_sets["origin_raw_locator"],
        edge_sets["exact_peak_hash"],
    )
    for members in ordered:
        cluster_id = "mirror31:" + members[0]
        member_set = set(members)
        cluster_edge_types = sorted(
            edge_type
            for edge_type, edges in edge_sets.items()
            if any(left in member_set and right in member_set for left, right in edges)
        )
        clusters.append(
            {
                "mirror_cluster_id": cluster_id,
                "members": members,
                "member_count": len(members),
                "edge_types": cluster_edge_types,
            }
        )
        for spectrum_id in members:
            assignment[spectrum_id] = cluster_id
    return {
        "schema_version": "3.1.0",
        "exact_edge_count_by_type": {
            edge_type: len(edge_sets[edge_type])
            for edge_type in (
                "canonical_usi_core",
                "origin_raw_locator",
                "exact_peak_hash",
            )
        },
        "unique_exact_edge_count": len(all_exact_edges),
        "near_duplicate_edge_count": len(edge_sets["calibrated_near_mirror"]),
        "component_count": len(clusters),
        "component_size_counts": dict(
            sorted(Counter(row["member_count"] for row in clusters).items())
        ),
        "near_calibration": calibration
        or {
            "status": "not_provided",
            "threshold": None,
            "near_mirror_links_allowed": False,
        },
        "near_candidate_status": "evaluated_with_labelled_calibration"
        if near_allowed
        else "not_evaluated_without_labelled_calibration",
        "near_candidates": near_candidates,
        "assignment": assignment,
        "clusters": clusters,
    }
