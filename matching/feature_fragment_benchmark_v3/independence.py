"""Mirror detection with labelled-threshold calibration and safe fallbacks."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from matching.feature_fragment_benchmark.mirror import (
    UnionFind,
    greedy_peak_cosine,
    sparse_peak_signature,
)


def calibrate_near_mirror_threshold(
    labelled_pairs: Sequence[Mapping[str, Any]],
    *,
    minimum_positive_pairs: int = 20,
    minimum_negative_pairs: int = 20,
    maximum_false_positive_rate: float = 0.01,
    minimum_bilateral_coverage: float = 0.8,
) -> dict[str, Any]:
    """Choose a cosine threshold only from labelled cross-repository pairs."""
    positives = [row for row in labelled_pairs if bool(row.get("is_duplicate"))]
    negatives = [row for row in labelled_pairs if not bool(row.get("is_duplicate"))]
    if len(positives) < minimum_positive_pairs or len(negatives) < minimum_negative_pairs:
        return {
            "status": "insufficient_labelled_pairs",
            "threshold": None,
            "minimum_bilateral_coverage": minimum_bilateral_coverage,
            "positive_pair_count": len(positives),
            "negative_pair_count": len(negatives),
            "near_mirror_links_allowed": False,
        }
    candidates = sorted(
        {
            float(row["cosine"])
            for row in labelled_pairs
            if float(row.get("coverage_left", 0.0)) >= minimum_bilateral_coverage
            and float(row.get("coverage_right", 0.0)) >= minimum_bilateral_coverage
        },
        reverse=True,
    )
    selected = None
    for threshold in candidates:
        false_positive = sum(
            float(row.get("cosine", 0.0)) >= threshold
            and float(row.get("coverage_left", 0.0)) >= minimum_bilateral_coverage
            and float(row.get("coverage_right", 0.0)) >= minimum_bilateral_coverage
            for row in negatives
        )
        fpr = false_positive / len(negatives)
        if fpr > maximum_false_positive_rate:
            continue
        true_positive = sum(
            float(row.get("cosine", 0.0)) >= threshold
            and float(row.get("coverage_left", 0.0)) >= minimum_bilateral_coverage
            and float(row.get("coverage_right", 0.0)) >= minimum_bilateral_coverage
            for row in positives
        )
        sensitivity = true_positive / len(positives)
        candidate = {
            "threshold": threshold,
            "false_positive_rate": fpr,
            "sensitivity": sensitivity,
        }
        if selected is None or (sensitivity, threshold) > (
            selected["sensitivity"],
            selected["threshold"],
        ):
            selected = candidate
    if selected is None:
        return {
            "status": "no_threshold_meets_false_positive_gate",
            "threshold": None,
            "minimum_bilateral_coverage": minimum_bilateral_coverage,
            "positive_pair_count": len(positives),
            "negative_pair_count": len(negatives),
            "near_mirror_links_allowed": False,
        }
    return {
        "status": "calibrated",
        **selected,
        "maximum_false_positive_rate": maximum_false_positive_rate,
        "minimum_bilateral_coverage": minimum_bilateral_coverage,
        "positive_pair_count": len(positives),
        "negative_pair_count": len(negatives),
        "near_mirror_links_allowed": True,
    }


def cluster_mirrors(
    records: Sequence[Mapping[str, Any]],
    *,
    near_calibration: Mapping[str, Any] | None = None,
    mz_tolerance: float = 0.01,
) -> dict[str, Any]:
    """Cluster exact mirrors and calibrated near mirrors.

    Without a labelled calibration, suspected near duplicates are reported but
    cannot create links.  Exact peak-content hashes always share a component.
    """
    ids = [str(row["spectrum_id"]) for row in records]
    union = UnionFind(ids)
    by_hash: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in records:
        by_hash[str(row["independence"]["exact_peak_hash"])].append(row)
    exact_edges = 0
    for rows in by_hash.values():
        anchor = str(rows[0]["spectrum_id"])
        for row in rows[1:]:
            union.union(anchor, str(row["spectrum_id"]))
            exact_edges += 1

    calibration = dict(near_calibration or {})
    near_allowed = bool(
        calibration.get("status") == "calibrated"
        and calibration.get("near_mirror_links_allowed")
        and calibration.get("threshold") is not None
    )
    if not near_allowed:
        components: dict[str, list[str]] = defaultdict(list)
        for spectrum_id in ids:
            components[union.find(spectrum_id)].append(spectrum_id)
        ordered = sorted(
            (sorted(members) for members in components.values()), key=lambda x: x[0]
        )
        assignment: dict[str, str] = {}
        clusters = []
        for members in ordered:
            cluster_id = "mirror:" + members[0]
            clusters.append(
                {
                    "mirror_cluster_id": cluster_id,
                    "members": members,
                    "member_count": len(members),
                }
            )
            for spectrum_id in members:
                assignment[spectrum_id] = cluster_id
        return {
            "schema_version": "3.0.0",
            "exact_duplicate_edge_count": exact_edges,
            "near_duplicate_edge_count": 0,
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
            "near_candidate_status": "not_evaluated_without_labelled_calibration",
            "suspected_near_mirrors": [],
            "assignment": assignment,
            "clusters": clusters,
        }
    signatures: dict[tuple[Any, ...], list[Mapping[str, Any]]] = defaultdict(list)
    for row in records:
        acquisition = row["acquisition"]
        key = (
            row["compound"]["connectivity"],
            acquisition.get("polarity"),
            acquisition.get("adduct"),
            sparse_peak_signature(row["peaks"]),
        )
        signatures[key].append(row)
    suspected: list[dict[str, Any]] = []
    accepted_near = 0
    for rows in signatures.values():
        if len(rows) < 2:
            continue
        for index, left in enumerate(rows):
            for right in rows[index + 1 :]:
                if (
                    left["independence"]["exact_peak_hash"]
                    == right["independence"]["exact_peak_hash"]
                ):
                    continue
                cosine, coverage_left, coverage_right = greedy_peak_cosine(
                    left["peaks"], right["peaks"], mz_tolerance=mz_tolerance
                )
                if cosine < 0.9 or min(coverage_left, coverage_right) < 0.7:
                    continue
                row = {
                    "left": left["spectrum_id"],
                    "right": right["spectrum_id"],
                    "cosine": cosine,
                    "coverage_left": coverage_left,
                    "coverage_right": coverage_right,
                    "cross_repository": left["source"].get("repository")
                    != right["source"].get("repository"),
                    "linked": False,
                }
                if (
                    near_allowed
                    and cosine >= float(calibration["threshold"])
                    and min(coverage_left, coverage_right)
                    >= float(calibration["minimum_bilateral_coverage"])
                ):
                    union.union(str(left["spectrum_id"]), str(right["spectrum_id"]))
                    row["linked"] = True
                    accepted_near += 1
                suspected.append(row)

    components: dict[str, list[str]] = defaultdict(list)
    for spectrum_id in ids:
        components[union.find(spectrum_id)].append(spectrum_id)
    ordered = sorted((sorted(members) for members in components.values()), key=lambda x: x[0])
    assignment: dict[str, str] = {}
    clusters = []
    for members in ordered:
        cluster_id = "mirror:" + members[0]
        clusters.append({"mirror_cluster_id": cluster_id, "members": members, "member_count": len(members)})
        for spectrum_id in members:
            assignment[spectrum_id] = cluster_id
    return {
        "schema_version": "3.0.0",
        "exact_duplicate_edge_count": exact_edges,
        "near_duplicate_edge_count": accepted_near,
        "component_count": len(clusters),
        "component_size_counts": dict(sorted(Counter(row["member_count"] for row in clusters).items())),
        "near_calibration": calibration
        or {
            "status": "not_provided",
            "threshold": None,
            "near_mirror_links_allowed": False,
        },
        "suspected_near_mirrors": suspected,
        "assignment": assignment,
        "clusters": clusters,
    }
