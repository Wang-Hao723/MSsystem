"""Chromatographic peak-shape and candidate in-source-fragment evidence.

Peak-shape formulas in this module are adapted from MassCube 1.0.22
(`feature_evaluation.py` and `feature_grouping.py`), authored by Huaxu Yu and
distributed under the MIT License.  MassCube 1.0.22 is the version used by
the benchmark paper (doi:10.1038/s41467-025-60640-5).

Copyright (c) 2023 Huaxu Yu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

The ISF evidence gate intentionally differs from MassCube's fallback: fewer
than five common scans and constant traces fail closed instead of receiving a
correlation of 1.0. At least five common scans must also have measured positive
signal in both features; padded zeros cannot turn a lone spike into a peak.
The result is candidate evidence only and never requests automatic removal.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np


def _finite_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return parsed if math.isfinite(parsed) else None


def _peak_flags(rt: np.ndarray, intensity: np.ndarray) -> list[str]:
    flags: list[str] = []
    if rt.ndim != 1 or intensity.ndim != 1:
        return ["non_1d_input"]
    if rt.size != intensity.size:
        flags.append("length_mismatch")
        return flags
    if rt.size == 0:
        flags.append("missing_trace")
        return flags
    if rt.size < 5:
        flags.append("too_few_points")
    if not np.all(np.isfinite(rt)) or not np.all(np.isfinite(intensity)):
        flags.append("nonfinite_values")
        return flags
    if np.any(intensity < 0):
        flags.append("negative_intensity")
    if rt.size > 1 and np.any(np.diff(rt) <= 0):
        flags.append("non_increasing_rt")
    if np.all(intensity == 0):
        flags.append("zero_intensity")
    elif np.all(intensity == intensity[0]):
        flags.append("constant_intensity")
    if intensity.size and int(np.argmax(intensity)) in (0, intensity.size - 1):
        flags.append("endpoint_apex")
    return flags


def _masscube_gaussian_similarity(rt: np.ndarray, intensity: np.ndarray) -> float:
    """MassCube 1.0.22 half-height width estimate and cosine formula."""
    amplitude = float(np.max(intensity))
    apex_idx = int(np.argmax(intensity))
    mean = float(rt[apex_idx])

    right_width = float(rt[-1] - mean)
    left_width = float(mean - rt[0])
    for idx in range(apex_idx, intensity.size):
        if intensity[idx] < amplitude / 2:
            right_width = float(rt[idx] - mean)
            break
    for idx in range(apex_idx, 0, -1):
        if intensity[idx] < amplitude / 2:
            left_width = float(mean - rt[idx])
            break

    sigma = (right_width + left_width) / 2.355
    # Eligible traces have an interior apex and strictly increasing RT, so
    # sigma and both cosine norms are positive without warning suppression.
    fitted = amplitude * np.exp(-0.5 * ((rt - mean) / sigma) ** 2)
    value = np.dot(intensity, fitted) / (
        np.linalg.norm(intensity) * np.linalg.norm(fitted)
    )
    return float(value)


def _masscube_noise_score(intensity: np.ndarray) -> float:
    retained = intensity[intensity > np.max(intensity) * 0.05]
    if retained.size < 5:
        return 0.0
    signs = np.sign(np.diff(retained))
    counter = -1
    for idx in range(1, signs.size):
        if signs[idx] != signs[idx - 1]:
            counter += 1
    return 0.0 if counter == -1 else float(counter / (retained.size - 2))


def _masscube_asymmetry_factor(intensity: np.ndarray) -> float:
    apex_idx = int(np.argmax(intensity))
    below = intensity < 0.1 * intensity[apex_idx]
    left = np.where(below[:apex_idx])[0]
    right = np.where(below[apex_idx:])[0]
    left_idx = 0 if left.size == 0 else int(left[-1])
    right_idx = intensity.size - 1 if right.size == 0 else int(right[0] + apex_idx)
    return float((right_idx - apex_idx) / (apex_idx - left_idx))


def evaluate_peak_quality(rt_minutes: Sequence[float], intensities: Sequence[float]) -> dict[str, Any]:
    """Evaluate one chromatographic peak without mutating the input arrays.

    Valid traces reproduce MassCube 1.0.22's three raw metrics.  Invalid or
    degenerate traces return JSON-safe ``None`` metrics and explicit flags.
    """
    try:
        rt = np.asarray(rt_minutes, dtype=float)
        intensity = np.asarray(intensities, dtype=float)
    except (TypeError, ValueError, OverflowError):
        rt = np.asarray([], dtype=float)
        intensity = np.asarray([], dtype=float)
        flags = ["non_numeric_input"]
    else:
        flags = _peak_flags(rt, intensity)

    n_points = int(rt.size) if rt.ndim == 1 else 0
    apex_rt: float | None = None
    if (
        rt.ndim == intensity.ndim == 1
        and rt.size == intensity.size
        and rt.size > 0
        and np.all(np.isfinite(rt))
        and np.all(np.isfinite(intensity))
    ):
        apex_rt = float(rt[int(np.argmax(intensity))])

    result: dict[str, Any] = {
        "status": "ok" if not flags else "invalid_peak",
        "n_points": n_points,
        "gaussian_similarity": None,
        "noise_score": None,
        "asymmetry_factor": None,
        "flags": flags,
        "apex_rt_minutes": apex_rt,
    }
    if flags:
        return result

    gaussian = _masscube_gaussian_similarity(rt, intensity)
    noise = _masscube_noise_score(intensity)
    asymmetry = _masscube_asymmetry_factor(intensity)
    metrics = (gaussian, noise, asymmetry)
    if not all(math.isfinite(value) for value in metrics):
        result["status"] = "invalid_peak"
        result["flags"] = ["nonfinite_metric"]
        return result
    result.update(
        gaussian_similarity=gaussian,
        noise_score=noise,
        asymmetry_factor=asymmetry,
    )
    return result


def _base_isf_result() -> dict[str, Any]:
    return {
        "supported": False,
        "status": "unsupported",
        "reasons": [],
        "matched_fragment_mz": None,
        "mass_error_da": None,
        "common_points": 0,
        "common_positive_points": 0,
        "correlation": None,
        "apex_rt_delta_min": None,
        "evidence_level": "candidate_isf_relation",
        "automatic_remove": False,
    }


def _parse_trace(trace: Any, label: str, min_points: int) -> tuple[dict[str, Any] | None, list[str]]:
    reasons: list[str] = []
    if not isinstance(trace, Mapping):
        return None, [f"{label}_trace_missing"]
    try:
        rt = np.asarray(trace.get("rt", []), dtype=float)
        intensity = np.asarray(trace.get("intensity", []), dtype=float)
    except (TypeError, ValueError, OverflowError):
        return None, [f"{label}_trace_non_numeric"]
    if rt.ndim != 1 or intensity.ndim != 1 or rt.size != intensity.size:
        return None, [f"{label}_trace_shape_invalid"]
    if rt.size < min_points:
        reasons.append(f"{label}_trace_too_short")
    if rt.size == 0:
        return None, reasons
    if not np.all(np.isfinite(rt)) or not np.all(np.isfinite(intensity)):
        reasons.append(f"{label}_trace_nonfinite")
    if np.any(intensity < 0):
        reasons.append(f"{label}_trace_negative_intensity")
    if rt.size > 1 and np.any(np.diff(rt) <= 0):
        reasons.append(f"{label}_trace_non_increasing_rt")
    if np.all(np.isfinite(intensity)) and (
        np.all(intensity == 0) or np.all(intensity == intensity[0])
    ):
        reasons.append(f"{label}_trace_constant")

    scan_ids = trace.get("scan_ids")
    if scan_ids is not None:
        try:
            scan_ids = list(scan_ids)
        except TypeError:
            reasons.append(f"{label}_scan_ids_invalid")
            scan_ids = None
        else:
            try:
                unique_scan_count = len(set(scan_ids))
            except TypeError:
                unique_scan_count = -1
            if len(scan_ids) != rt.size or unique_scan_count != len(scan_ids):
                reasons.append(f"{label}_scan_ids_invalid")
    if reasons:
        return None, reasons
    return {"rt": rt, "intensity": intensity, "scan_ids": scan_ids}, []


def _common_intensities(parent: dict[str, Any], child: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, int, str | None]:
    if parent["scan_ids"] is not None and child["scan_ids"] is not None:
        parent_keys = parent["scan_ids"]
        child_keys = child["scan_ids"]
    else:
        parent_keys = [round(float(value), 8) for value in parent["rt"]]
        child_keys = [round(float(value), 8) for value in child["rt"]]
        if len(set(parent_keys)) != len(parent_keys) or len(set(child_keys)) != len(child_keys):
            return np.array([]), np.array([]), 0, "ambiguous_rounded_rt"

    child_by_key = {key: idx for idx, key in enumerate(child_keys)}
    parent_values: list[float] = []
    child_values: list[float] = []
    for idx, key in enumerate(parent_keys):
        child_idx = child_by_key.get(key)
        if child_idx is not None:
            parent_values.append(float(parent["intensity"][idx]))
            child_values.append(float(child["intensity"][child_idx]))
    return (
        np.asarray(parent_values, dtype=float),
        np.asarray(child_values, dtype=float),
        len(parent_values),
        None,
    )


def assess_isf_pair(
    parent_mz: float,
    child_mz: float,
    parent_ms2_peaks: Sequence[Sequence[float]],
    parent_trace: Mapping[str, Sequence[Any]],
    child_trace: Mapping[str, Sequence[Any]],
    *,
    parent_polarity: str | None = None,
    child_polarity: str | None = None,
    mz_tolerance_da: float = 0.01,
    rt_tolerance_min: float = 0.05,
    correlation_min: float = 0.7,
    min_common_points: int = 5,
) -> dict[str, Any]:
    """Assess evidence that a lower-mass feature is an in-source fragment.

    Only measured parent-MS2 mass evidence and simultaneous MS1 scans are
    used. At least ``min_common_points`` simultaneous scans must have positive
    measured intensity in both traces. Intensities are never interpolated,
    and no identity, library score, or truth label is accepted by this API.
    """
    result = _base_isf_result()
    reasons: list[str] = []
    parent_mass = _finite_float(parent_mz)
    child_mass = _finite_float(child_mz)
    mz_tol = _finite_float(mz_tolerance_da)
    rt_tol = _finite_float(rt_tolerance_min)
    corr_min = _finite_float(correlation_min)
    if parent_mass is None or parent_mass <= 0:
        reasons.append("parent_mz_invalid")
    if child_mass is None or child_mass <= 0:
        reasons.append("child_mz_invalid")
    if parent_mass is not None and child_mass is not None and child_mass >= parent_mass:
        reasons.append("child_not_lower_mass")
    if mz_tol is None or mz_tol < 0:
        reasons.append("mz_tolerance_invalid")
    if rt_tol is None or rt_tol < 0:
        reasons.append("rt_tolerance_invalid")
    if corr_min is None or not -1 <= corr_min <= 1:
        reasons.append("correlation_threshold_invalid")
    if not isinstance(min_common_points, int) or isinstance(min_common_points, bool) or min_common_points < 5:
        reasons.append("min_common_points_invalid")

    if parent_polarity is not None and child_polarity is not None:
        polarity_aliases = {"+": "positive", "pos": "positive", "-": "negative", "neg": "negative"}
        normalized_parent = str(parent_polarity).strip().lower()
        normalized_child = str(child_polarity).strip().lower()
        normalized_parent = polarity_aliases.get(normalized_parent, normalized_parent)
        normalized_child = polarity_aliases.get(normalized_child, normalized_child)
        if not normalized_parent or not normalized_child or normalized_parent != normalized_child:
            reasons.append("polarity_conflict")

    if reasons:
        result.update(status="invalid_input", reasons=reasons)
        return result

    assert child_mass is not None and mz_tol is not None and rt_tol is not None and corr_min is not None
    valid_fragments: list[tuple[float, float]] = []
    try:
        peaks = list(parent_ms2_peaks)
    except TypeError:
        peaks = []
    for peak in peaks:
        try:
            fragment_mz = float(peak[0])
            fragment_intensity = float(peak[1])
        except (TypeError, ValueError, IndexError, OverflowError):
            continue
        if (
            math.isfinite(fragment_mz)
            and math.isfinite(fragment_intensity)
            and fragment_mz > 0
            and fragment_intensity > 0
        ):
            valid_fragments.append((fragment_mz, fragment_intensity))
    matches = [peak for peak in valid_fragments if abs(peak[0] - child_mass) <= mz_tol]
    if not matches:
        reasons.append("child_mass_absent_from_parent_ms2")
    else:
        matched = min(matches, key=lambda peak: (abs(peak[0] - child_mass), -peak[1], peak[0]))
        result["matched_fragment_mz"] = float(matched[0])
        result["mass_error_da"] = float(matched[0] - child_mass)

    parent, parent_reasons = _parse_trace(parent_trace, "parent", min_common_points)
    child, child_reasons = _parse_trace(child_trace, "child", min_common_points)
    reasons.extend(parent_reasons)
    reasons.extend(child_reasons)
    if parent is None or child is None:
        result.update(status="insufficient_trace_evidence", reasons=reasons)
        return result

    parent_apex_rt = float(parent["rt"][int(np.argmax(parent["intensity"]))])
    child_apex_rt = float(child["rt"][int(np.argmax(child["intensity"]))])
    apex_delta = abs(parent_apex_rt - child_apex_rt)
    result["apex_rt_delta_min"] = float(apex_delta)
    if apex_delta > rt_tol:
        reasons.append("apex_rt_outside_tolerance")

    parent_common, child_common, common_points, common_error = _common_intensities(parent, child)
    result["common_points"] = int(common_points)
    common_positive_points = int(np.count_nonzero((parent_common > 0) & (child_common > 0)))
    result["common_positive_points"] = common_positive_points
    if common_error is None and common_positive_points < min_common_points:
        reasons.append("too_few_common_positive_points")
    if common_error is not None:
        reasons.append(common_error)
    elif common_points < min_common_points:
        reasons.append("too_few_common_points")
    elif np.all(parent_common == parent_common[0]) or np.all(child_common == child_common[0]):
        reasons.append("constant_common_intensity")
    else:
        correlation = float(np.corrcoef(parent_common, child_common)[0, 1])
        if math.isfinite(correlation):
            result["correlation"] = correlation
            if correlation <= corr_min:
                reasons.append("correlation_not_above_threshold")
        else:
            reasons.append("correlation_nonfinite")

    result["reasons"] = reasons
    if not reasons:
        result["supported"] = True
        result["status"] = "supported"
    return result


__all__ = ["evaluate_peak_quality", "assess_isf_pair"]
