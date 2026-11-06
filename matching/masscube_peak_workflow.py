"""Deterministic feature-level orchestration for MassCube-derived evidence.

This module combines chromatographic peak quality with conservative candidate
in-source-fragment (ISF) relations.  It does not delete features, consume
compound labels, or alter library candidates.
"""

from __future__ import annotations

import bisect
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from matching.masscube_peak_evidence import assess_isf_pair, evaluate_peak_quality


MZ_TOLERANCE_DA = 0.01
RT_TOLERANCE_MIN = 0.05
CORRELATION_MIN = 0.7
MIN_COMMON_POINTS = 5
CHILD_TO_PARENT_HEIGHT_MAX_EXCLUSIVE = 3.0


def _finite_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return parsed if math.isfinite(parsed) else None


def _json_scalar(value: Any) -> str | int | float | bool | None:
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return str(value)


def _missing_quality() -> dict[str, Any]:
    return {
        "status": "missing_trace",
        "n_points": 0,
        "gaussian_similarity": None,
        "noise_score": None,
        "asymmetry_factor": None,
        "flags": ["missing_trace"],
        "apex_rt_minutes": None,
    }


def _normalized_features(features: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    seen_ids: set[Any] = set()
    for position, feature in enumerate(features):
        if not isinstance(feature, Mapping):
            raise ValueError(f"feature at position {position} is not a mapping")
        feature_id = _json_scalar(feature.get("id"))
        try:
            hash(feature_id)
        except TypeError as exc:
            raise ValueError(f"feature at position {position} has an invalid id") from exc
        if feature_id is None or feature_id in seen_ids:
            raise ValueError(f"feature id at position {position} is missing or duplicated")
        seen_ids.add(feature_id)
        mz = _finite_float(feature.get("mz"))
        rt = _finite_float(feature.get("rt_minutes"))
        height = _finite_float(feature.get("height"))
        if mz is None or mz <= 0 or rt is None or height is None or height < 0:
            raise ValueError(f"feature {feature_id!r} has invalid mz, rt_minutes, or height")
        normalized.append(
            {
                "id": feature_id,
                "mz": mz,
                "rt_minutes": rt,
                "height": height,
                "is_isotope": bool(feature.get("is_isotope", False)),
            }
        )
    return normalized


def _record_tic(record: Mapping[str, Any]) -> float:
    tic = _finite_float(record.get("tic"))
    return tic if tic is not None and tic >= 0 else -math.inf


def _spectrum_order(record: Mapping[str, Any]) -> str:
    spectrum_id = _json_scalar(record.get("spectrum_id"))
    return "" if spectrum_id is None else str(spectrum_id)


def _select_parent_spectra(
    records: Sequence[Mapping[str, Any]], valid_feature_ids: set[Any]
) -> dict[Any, dict[str, Any]]:
    grouped: dict[Any, list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        if not isinstance(record, Mapping):
            continue
        feature_id = _json_scalar(record.get("feature_id"))
        if feature_id in valid_feature_ids:
            grouped[feature_id].append(record)

    selected: dict[Any, dict[str, Any]] = {}
    for feature_id, candidates in grouped.items():
        # Highest TIC wins.  Equal TICs use the lexicographically smallest
        # spectrum_id so selection is stable across input order.
        best = min(candidates, key=lambda item: (-_record_tic(item), _spectrum_order(item)))
        selected[feature_id] = {
            "feature_id": feature_id,
            "peaks": best.get("peaks", []),
            "tic": None if _record_tic(best) == -math.inf else _record_tic(best),
            "spectrum_id": _json_scalar(best.get("spectrum_id")),
        }
    return selected


def _valid_fragment_masses(peaks: Any) -> list[float]:
    try:
        peak_records = list(peaks)
    except TypeError:
        return []
    masses: list[float] = []
    for peak in peak_records:
        try:
            mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError, IndexError, OverflowError):
            continue
        if math.isfinite(mz) and mz > 0 and math.isfinite(intensity) and intensity > 0:
            masses.append(mz)
    return masses


def annotate_feature_evidence(
    features: Sequence[Mapping[str, Any]],
    ms2_records: Sequence[Mapping[str, Any]],
    traces: Mapping[Any, Mapping[str, Sequence[Any]]],
    *,
    polarity: str | None = None,
) -> dict[str, Any]:
    """Annotate peak quality and candidate ISF relations without deletion.

    Parent features are processed by descending m/z, matching MassCube's
    outer traversal.  One measured MS2 is
    fixed per parent by highest TIC, with ``spectrum_id`` as the deterministic
    tie-break.  Fragment masses query a sorted child-m/z index through
    ``bisect``; the workflow never performs a full feature-pair by spectrum-
    peak scan.

    Every relation also needs five simultaneous positive MS1 measurements
    in both features; zero-filled XIC padding contributes no signal support.
    """
    if not isinstance(traces, Mapping):
        raise ValueError("traces must be a mapping keyed by feature id")
    normalized = _normalized_features(features)
    feature_by_id = {feature["id"]: feature for feature in normalized}
    selected_spectra = _select_parent_spectra(ms2_records, set(feature_by_id))

    output_features: list[dict[str, Any]] = []
    trace_by_id: dict[Any, Any] = {}
    for feature in normalized:
        trace = traces.get(feature["id"])
        trace_by_id[feature["id"]] = trace
        if isinstance(trace, Mapping):
            quality = evaluate_peak_quality(trace.get("rt", []), trace.get("intensity", []))
        else:
            quality = _missing_quality()
        output_features.append(
            {
                **feature,
                "quality": quality,
                "candidate_isf_parent_ids": [],
                "automatic_remove": False,
            }
        )
    output_by_id = {feature["id"]: feature for feature in output_features}

    child_candidates = [feature for feature in normalized if not feature["is_isotope"]]
    child_candidates.sort(key=lambda feature: (feature["mz"], str(feature["id"])))
    child_masses = [feature["mz"] for feature in child_candidates]

    parents = [
        feature
        for feature in normalized
        if not feature["is_isotope"] and feature["id"] in selected_spectra
    ]
    parents.sort(key=lambda feature: (-feature["mz"], str(feature["id"])))

    relations: list[dict[str, Any]] = []
    considered_pairs: set[tuple[Any, Any]] = set()
    candidate_prefilter_count = 0
    assessed_pair_count = 0
    for parent in parents:
        spectrum = selected_spectra[parent["id"]]
        peaks = spectrum["peaks"]
        for fragment_mz in _valid_fragment_masses(peaks):
            left = bisect.bisect_left(child_masses, fragment_mz - MZ_TOLERANCE_DA)
            right = bisect.bisect_right(child_masses, fragment_mz + MZ_TOLERANCE_DA)
            for child in child_candidates[left:right]:
                pair = (parent["id"], child["id"])
                if pair in considered_pairs:
                    continue
                considered_pairs.add(pair)
                candidate_prefilter_count += 1
                if child["mz"] >= parent["mz"]:
                    continue
                if abs(child["rt_minutes"] - parent["rt_minutes"]) > RT_TOLERANCE_MIN:
                    continue
                if child["height"] >= CHILD_TO_PARENT_HEIGHT_MAX_EXCLUSIVE * parent["height"]:
                    continue
                assessed_pair_count += 1
                evidence = assess_isf_pair(
                    parent["mz"],
                    child["mz"],
                    peaks,
                    trace_by_id[parent["id"]],
                    trace_by_id[child["id"]],
                    parent_polarity=polarity,
                    child_polarity=polarity,
                    mz_tolerance_da=MZ_TOLERANCE_DA,
                    rt_tolerance_min=RT_TOLERANCE_MIN,
                    correlation_min=CORRELATION_MIN,
                    min_common_points=MIN_COMMON_POINTS,
                )
                if evidence["supported"]:
                    relations.append(
                        {
                            "parent_id": parent["id"],
                            "child_id": child["id"],
                            "parent_spectrum_id": spectrum["spectrum_id"],
                            **evidence,
                            "ambiguous_parent": False,
                        }
                    )

    parent_ids_by_child: dict[Any, list[Any]] = defaultdict(list)
    for relation in relations:
        parent_ids_by_child[relation["child_id"]].append(relation["parent_id"])
    ambiguous_children = {
        child_id for child_id, parent_ids in parent_ids_by_child.items() if len(parent_ids) > 1
    }
    for relation in relations:
        relation["ambiguous_parent"] = relation["child_id"] in ambiguous_children
    for child_id, parent_ids in parent_ids_by_child.items():
        output_by_id[child_id]["candidate_isf_parent_ids"] = parent_ids

    quality_ok = sum(feature["quality"]["status"] == "ok" for feature in output_features)
    missing_trace = sum(
        feature["quality"]["status"] == "missing_trace" for feature in output_features
    )
    return {
        "features": output_features,
        "isf_relations": relations,
        "summary": {
            "feature_count": len(output_features),
            "quality_ok_count": quality_ok,
            "quality_missing_trace_count": missing_trace,
            "selected_parent_spectrum_count": len(selected_spectra),
            "candidate_prefilter_count": candidate_prefilter_count,
            "assessed_pair_count": assessed_pair_count,
            "isf_relation_count": len(relations),
            "ambiguous_child_count": len(ambiguous_children),
            "automatic_remove_count": 0,
        },
    }


__all__ = [
    "annotate_feature_evidence",
    "MZ_TOLERANCE_DA",
    "RT_TOLERANCE_MIN",
    "CORRELATION_MIN",
    "MIN_COMMON_POINTS",
    "CHILD_TO_PARENT_HEIGHT_MAX_EXCLUSIVE",
]
