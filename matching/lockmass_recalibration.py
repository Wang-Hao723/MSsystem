"""Candidate-independent lock-mass audit for precursor recalibration.

The module measures an explicitly declared vendor reference ion in MS1 scans.
It produces a shadow correction artifact only; chemical library matches are
never used to estimate the correction.
"""

from __future__ import annotations

import bisect
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

from matching.library_quality import compound_identity_key


SCHEMA_VERSION = "1.0.0"
MODULE_ID = "candidate_independent_lockmass_recalibration"
AGILENT_NEGATIVE_REFERENCE_MZ = 112.985587
AGILENT_POSITIVE_REFERENCE_MZ = 121.050873


def _number(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _retention_time_minutes(spectrum: Mapping[str, Any]) -> float | None:
    scans = ((spectrum.get("scanList") or {}).get("scan") or [])
    if not scans:
        return None
    raw = scans[0].get("scan start time")
    value = _number(raw)
    if value is None:
        return None
    unit = " ".join((
        str(getattr(raw, "unit_info", "") or ""),
        str(getattr(raw, "unit_name", "") or ""),
    )).lower()
    if "millisecond" in unit:
        return value / 60000.0
    if "second" in unit and "minute" not in unit:
        return value / 60.0
    return value


def summarize_lockmass_observations(
    observations: Iterable[Mapping[str, Any]],
    *,
    task_id: str,
    source_mzml: Path,
    reference_mz: float,
    ms1_scan_count: int,
    search_window_da: float,
    minimum_intensity: float,
    minimum_observations: int = 50,
    minimum_presence_fraction: float = 0.25,
    maximum_mad_ppm: float = 5.0,
) -> dict[str, Any]:
    rows = [dict(item) for item in observations]
    ppm_values = [float(item["ppm_error"]) for item in rows]
    presence = len(rows) / ms1_scan_count if ms1_scan_count else 0.0
    base = {
        "schema_version": SCHEMA_VERSION,
        "task_id": task_id,
        "module_id": MODULE_ID,
        "status": "not_qualified",
        "deployment_state": "shadow",
        "runtime_influence": False,
        "ranking_changed": False,
        "confidence_changed": False,
        "claim_scope_changed": False,
        "candidate_independent": True,
        "chemical_identity_used_for_fit": False,
        "source_mzml": str(source_mzml),
        "reference": {
            "reference_mz": reference_mz,
            "polarity": "NEG" if reference_mz == AGILENT_NEGATIVE_REFERENCE_MZ else "POS",
            "vendor": "Agilent",
            "role": "continuous_internal_reference_mass",
            "source": "declared_vendor_reference_ion_not_library_candidate",
        },
        "selection": {
            "ms_level": 1,
            "search_window_da": search_window_da,
            "minimum_intensity": minimum_intensity,
            "minimum_observations": minimum_observations,
            "minimum_presence_fraction": minimum_presence_fraction,
            "maximum_mad_ppm": maximum_mad_ppm,
        },
        "summary": {
            "ms1_scan_count": ms1_scan_count,
            "qualified_reference_observations": len(rows),
            "presence_fraction": round(presence, 6),
        },
        "rt_bins": [],
        "correction": None,
        "production_authorized": False,
        "interpretation": (
            "A qualified correction is suitable for shadow reprocessing only. "
            "Production activation requires a frozen standard/control benchmark."
        ),
    }
    if not ppm_values:
        base["summary"]["reason"] = "no_reference_observations"
        return base
    median_ppm = statistics.median(ppm_values)
    mad_ppm = statistics.median(abs(value - median_ppm) for value in ppm_values)
    ordered = sorted(ppm_values)
    q10 = ordered[round((len(ordered) - 1) * 0.10)]
    q90 = ordered[round((len(ordered) - 1) * 0.90)]
    base["summary"].update({
        "median_error_ppm": round(median_ppm, 6),
        "mad_error_ppm": round(mad_ppm, 6),
        "q10_error_ppm": round(q10, 6),
        "q90_error_ppm": round(q90, 6),
        "median_observed_mz": round(
            statistics.median(float(item["observed_mz"]) for item in rows), 9
        ),
        "median_intensity": round(
            statistics.median(float(item["intensity"]) for item in rows), 6
        ),
    })
    bins: dict[int, list[float]] = defaultdict(list)
    for item in rows:
        rt = _number(item.get("retention_time_min"))
        if rt is not None:
            bins[int(rt // 2) * 2].append(float(item["ppm_error"]))
    for start in sorted(bins):
        values = bins[start]
        local_median = statistics.median(values)
        local_mad = statistics.median(abs(value - local_median) for value in values)
        base["rt_bins"].append({
            "rt_start_min": float(start),
            "rt_end_min": float(start + 2),
            "observation_count": len(values),
            "median_error_ppm": round(local_median, 6),
            "mad_error_ppm": round(local_mad, 6),
        })
    reasons = []
    if len(rows) < minimum_observations:
        reasons.append("insufficient_reference_observations")
    if presence < minimum_presence_fraction:
        reasons.append("insufficient_reference_presence")
    if mad_ppm > maximum_mad_ppm:
        reasons.append("reference_peak_not_stable")
    if reasons:
        base["summary"]["reason_codes"] = reasons
        return base
    base["status"] = "qualified_shadow"
    base["correction"] = {
        "model": "single_reference_multiplicative_mass_axis",
        "estimated_error_ppm": round(median_ppm, 6),
        "formula": "corrected_mz = observed_mz / (1 + estimated_error_ppm * 1e-6)",
        "applies_to": "precursor_mz_shadow_reprocessing",
        "fragment_axis_effect_at_1000_mz_da": round(abs(median_ppm) * 0.001, 6),
        "runtime_influence": False,
    }
    return base


def audit_mzml_lockmass(
    mzml_path: Path,
    *,
    task_id: str,
    reference_mz: float = AGILENT_NEGATIVE_REFERENCE_MZ,
    search_window_da: float = 0.006,
    minimum_intensity: float = 1000.0,
    minimum_observations: int = 50,
    minimum_presence_fraction: float = 0.25,
    maximum_mad_ppm: float = 5.0,
) -> dict[str, Any]:
    try:
        from pyteomics import mzml
    except ImportError as exc:
        return {
            "schema_version": SCHEMA_VERSION,
            "task_id": task_id,
            "module_id": MODULE_ID,
            "status": "unavailable",
            "deployment_state": "shadow",
            "runtime_influence": False,
            "candidate_independent": True,
            "summary": {"reason": f"pyteomics_unavailable: {exc}"},
        }
    observations: list[dict[str, Any]] = []
    ms1_scan_count = 0
    with mzml.read(str(mzml_path)) as reader:
        for spectrum in reader:
            if spectrum.get("ms level") != 1:
                continue
            ms1_scan_count += 1
            mz_array = spectrum.get("m/z array")
            intensity_array = spectrum.get("intensity array")
            if mz_array is None or intensity_array is None:
                continue
            mzs = mz_array.tolist()
            intensities = intensity_array.tolist()
            start = bisect.bisect_left(mzs, reference_mz - search_window_da)
            stop = bisect.bisect_right(mzs, reference_mz + search_window_da)
            if stop <= start:
                continue
            index = max(range(start, stop), key=lambda item: intensities[item])
            intensity = _number(intensities[index])
            observed_mz = _number(mzs[index])
            if (
                intensity is None
                or observed_mz is None
                or intensity < minimum_intensity
            ):
                continue
            observations.append({
                "observed_mz": observed_mz,
                "intensity": intensity,
                "ppm_error": (observed_mz - reference_mz) / reference_mz * 1e6,
                "retention_time_min": _retention_time_minutes(spectrum),
            })
    return summarize_lockmass_observations(
        observations,
        task_id=task_id,
        source_mzml=mzml_path,
        reference_mz=reference_mz,
        ms1_scan_count=ms1_scan_count,
        search_window_da=search_window_da,
        minimum_intensity=minimum_intensity,
        minimum_observations=minimum_observations,
        minimum_presence_fraction=minimum_presence_fraction,
        maximum_mad_ppm=maximum_mad_ppm,
    )


def corrected_mz(observed_mz: Any, artifact: Mapping[str, Any]) -> float | None:
    value = _number(observed_mz)
    if value is None or artifact.get("status") != "qualified_shadow":
        return None
    correction = artifact.get("correction") or {}
    ppm = _number(correction.get("estimated_error_ppm"))
    if ppm is None:
        return None
    return value / (1.0 + ppm * 1e-6)


def build_selected_match_counterfactual(
    results: Iterable[Mapping[str, Any]], artifact: Mapping[str, Any]
) -> dict[str, Any]:
    """Audit fixed selected matches before/after correction without reranking."""
    rows = []
    for row_index, row in enumerate(results):
        match = row.get("match") if isinstance(row, Mapping) else None
        if not isinstance(match, Mapping):
            continue
        query = _number(match.get("mass_query_value"))
        library = _number(match.get("mass_library_value"))
        corrected = corrected_mz(query, artifact)
        if query is None or library is None or corrected is None or corrected <= 0:
            continue
        before = abs(query - library) / query * 1e6
        after = abs(corrected - library) / corrected * 1e6
        rows.append({
            "row_index": row_index,
            "spectrum_id": row.get("spectrum_id"),
            "identity_key": compound_identity_key(match),
            "compound_name": match.get("compound_name"),
            "existing_confidence": match.get("confidence"),
            "mass_match_basis": match.get("mass_match_basis"),
            "observed_query_mass": query,
            "corrected_query_mass": corrected,
            "library_mass": library,
            "before_abs_ppm": round(before, 6),
            "after_abs_ppm": round(after, 6),
            "selected_candidate_fixed": True,
            "runtime_influence": False,
        })
    high_identities: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for item in rows:
        if item["existing_confidence"] == "High" and item["identity_key"]:
            high_identities[item["identity_key"]].append(item)
    before_values = [item["before_abs_ppm"] for item in rows]
    after_values = [item["after_abs_ppm"] for item in rows]
    return {
        "status": "completed_shadow" if rows else "not_applicable",
        "candidate_set_fixed": True,
        "retrieval_rerun": False,
        "ranking_changed": False,
        "confidence_changed": False,
        "rows": rows,
        "summary": {
            "selected_matches_evaluated": len(rows),
            "before_within_5ppm": sum(value <= 5 for value in before_values),
            "after_within_5ppm": sum(value <= 5 for value in after_values),
            "before_within_10ppm": sum(value <= 10 for value in before_values),
            "after_within_10ppm": sum(value <= 10 for value in after_values),
            "before_median_abs_ppm": (
                round(statistics.median(before_values), 6) if before_values else None
            ),
            "after_median_abs_ppm": (
                round(statistics.median(after_values), 6) if after_values else None
            ),
            "existing_high_unique_identities": len(high_identities),
            "existing_high_unique_identities_after_within_10ppm": sum(
                min(item["after_abs_ppm"] for item in items) <= 10
                for items in high_identities.values()
            ),
            "interpretation": (
                "Fixed-candidate mass-axis counterfactual; it does not estimate "
                "new retrieval recall, Top-1 accuracy, or High promotions."
            ),
        },
    }


__all__ = [
    "AGILENT_NEGATIVE_REFERENCE_MZ",
    "AGILENT_POSITIVE_REFERENCE_MZ",
    "audit_mzml_lockmass",
    "build_selected_match_counterfactual",
    "corrected_mz",
    "summarize_lockmass_observations",
]
