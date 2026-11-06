"""Link OpenMS chromatographic features to DDA spectra and build consensus MS2."""
from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from typing import Any

import numpy as np

from utils.spectrum import PeakProcessingConfig, process_peaks


def _number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _charge_state(value: Any) -> int | None:
    """Return an unsigned precursor charge state from vendor/OpenMS metadata.

    mzML selected-ion charge is commonly stored as a positive magnitude even
    in negative mode, whereas OpenMS may serialize the same state as ``-1``.
    Polarity is carried separately, so comparing signs would reject valid NEG
    links.  Only a genuine state mismatch such as |z| 1 versus |z| 2 is a
    charge conflict here.
    """

    number = _number(value)
    if number in (None, 0):
        return None
    rounded = int(round(abs(number)))
    return rounded if rounded > 0 else None


def resolve_search_precursor_mz(
    spectrum: dict[str, Any],
    *,
    maximum_link_delta_ppm: float = 50.0,
    absolute_floor_da: float = 0.02,
) -> dict[str, Any]:
    """Choose the measured mass axis used for high-resolution library search.

    Vendor mzML commonly stores the quadrupole isolation target as the MS2
    ``selected ion m/z``.  That acquisition setpoint is traceable, but it is
    not necessarily the measured precursor centroid.  Once an MS2 event has a
    unique OpenMS continuous-MS1 feature link, the feature centroid is the
    appropriate HRMS retrieval/scoring mass.  The selected-ion value remains
    unchanged in the result for raw-data navigation and audit.

    The function fails closed to the selected-ion value when the feature link
    is absent, ambiguous, non-finite, charge-inconsistent, or implausibly far
    away.  It never infers a chemical identity or adduct.
    """

    selected = _number(spectrum.get("precursor_mz"))
    base = {
        "selected_ion_mz": selected,
        "search_precursor_mz": selected,
        "search_precursor_mz_source": "ms2_selected_ion_mz_fallback",
        "search_precursor_mz_reliable": selected is not None and selected > 0,
        "lcms_feature_used_for_search": False,
        "search_precursor_delta_ppm": None,
        "fallback_reason": "continuous_ms1_feature_not_linked",
    }
    if selected is None or selected <= 0:
        base["fallback_reason"] = "selected_ion_mz_missing"
        return base
    if spectrum.get("lcms_feature_link_status") != "linked_openms_feature":
        base["fallback_reason"] = str(
            spectrum.get("lcms_feature_link_status")
            or "continuous_ms1_feature_not_linked"
        )
        return base
    feature_mz = _number(spectrum.get("lcms_feature_mz"))
    if feature_mz is None or feature_mz <= 0:
        base["fallback_reason"] = "linked_feature_mz_missing"
        return base

    selected_charge = _charge_state(spectrum.get("precursor_charge"))
    feature_charge = _charge_state(spectrum.get("lcms_feature_charge"))
    if (
        selected_charge is not None
        and feature_charge is not None
        and selected_charge != feature_charge
    ):
        base["fallback_reason"] = "linked_feature_charge_conflict"
        return base

    delta_da = feature_mz - selected
    delta_ppm = delta_da / feature_mz * 1e6
    allowed_da = max(
        float(absolute_floor_da),
        feature_mz * float(maximum_link_delta_ppm) / 1e6,
    )
    if abs(delta_da) > allowed_da:
        base["fallback_reason"] = "linked_feature_mass_delta_implausible"
        base["search_precursor_delta_ppm"] = round(delta_ppm, 6)
        return base

    return {
        **base,
        "search_precursor_mz": feature_mz,
        "search_precursor_mz_source": (
            "measured_continuous_ms1_feature_centroid"
        ),
        "search_precursor_mz_reliable": True,
        "lcms_feature_used_for_search": True,
        "search_precursor_delta_ppm": round(delta_ppm, 6),
        "fallback_reason": "",
        "lcms_feature_id": spectrum.get("lcms_feature_id"),
    }


def attach_ms2_to_lcms_features(
    spectra: list[dict],
    features: list[dict],
    *,
    mz_tolerance_ppm: float = 10.0,
    mz_floor_da: float = 0.003,
    rt_padding_minutes: float = 0.0,
) -> dict[str, Any]:
    """Attach each MS2 precursor to at most one chromatographic MS1 feature."""

    usable = [
        feature
        for feature in features or []
        if (_number(feature.get("mz")) or 0) > 0
        and _number(feature.get("rt_minutes")) is not None
    ]
    usable.sort(key=lambda feature: float(feature["mz"]))
    mz_values = [float(feature["mz"]) for feature in usable]
    linked = 0
    ambiguous = 0
    no_mass = 0
    no_rt = 0
    charge_conflict_candidates = 0

    for spectrum in spectra or []:
        precursor = _number(spectrum.get("precursor_mz"))
        rt = _number(spectrum.get("retention_time"))
        if precursor is None or precursor <= 0:
            spectrum["lcms_feature_link_status"] = "missing_precursor_mz"
            no_mass += 1
            continue
        if rt is None:
            spectrum["lcms_feature_link_status"] = "missing_retention_time"
            no_rt += 1
            continue
        tolerance = max(float(mz_floor_da), precursor * mz_tolerance_ppm * 1e-6)
        start = bisect_left(mz_values, precursor - tolerance)
        end = bisect_right(mz_values, precursor + tolerance)
        candidates = []
        for feature in usable[start:end]:
            precursor_charge = _charge_state(spectrum.get("precursor_charge"))
            feature_charge = _charge_state(feature.get("charge"))
            if (
                precursor_charge is not None
                and feature_charge is not None
                and precursor_charge != feature_charge
            ):
                charge_conflict_candidates += 1
                continue
            apex = float(feature["rt_minutes"])
            rt_start = _number(feature.get("rt_start_minutes"))
            rt_end = _number(feature.get("rt_end_minutes"))
            if rt_start is None:
                rt_start = apex
            if rt_end is None:
                rt_end = apex
            if not (
                rt_start - rt_padding_minutes
                <= rt
                <= rt_end + rt_padding_minutes
            ):
                continue
            mz_error = abs(float(feature["mz"]) - precursor) / tolerance
            rt_scale = max(
                rt_padding_minutes,
                (rt_end - rt_start) / 2.0,
                1e-6,
            )
            rt_error = abs(apex - rt) / rt_scale
            candidates.append(
                (
                    mz_error + 0.25 * rt_error,
                    -float(feature.get("intensity") or 0.0),
                    str(feature.get("feature_id") or ""),
                    feature,
                )
            )
        if not candidates:
            spectrum["lcms_feature_link_status"] = "no_feature_in_mz_rt_window"
            continue
        candidates.sort(key=lambda item: item[:3])
        if (
            len(candidates) > 1
            and abs(candidates[1][0] - candidates[0][0]) <= 0.05
        ):
            ambiguous += 1
            spectrum["lcms_feature_link_status"] = "ambiguous_best_feature"
            spectrum["lcms_feature_candidate_count"] = len(candidates)
            continue
        feature = candidates[0][3]
        spectrum.update(
            {
                "lcms_feature_id": feature.get("feature_id"),
                "lcms_feature_mz": feature.get("mz"),
                "lcms_feature_rt_minutes": feature.get("rt_minutes"),
                "lcms_feature_rt_start_minutes": feature.get("rt_start_minutes"),
                "lcms_feature_rt_end_minutes": feature.get("rt_end_minutes"),
                "lcms_feature_intensity": feature.get("intensity"),
                "lcms_feature_charge": feature.get("charge"),
                "lcms_isotope_trace_count": feature.get("isotope_trace_count"),
                "lcms_adduct": feature.get("adduct"),
                "lcms_neutral_mass": feature.get("neutral_mass"),
                "lcms_adduct_group_id": feature.get("adduct_group_id"),
                "lcms_feature_link_status": "linked_openms_feature",
                "lcms_feature_candidate_count": len(candidates),
                "lcms_feature_mz_error_ppm": round(
                    (precursor - float(feature["mz"]))
                    / float(feature["mz"])
                    * 1e6,
                    4,
                ),
                "lcms_feature_rt_error_seconds": round(
                    (rt - float(feature["rt_minutes"])) * 60.0,
                    4,
                ),
            }
        )
        linked += 1

    return {
        "method": "OpenMS MS1 feature-to-MS2 assignment (m/z + feature RT bounds)",
        "feature_count": len(usable),
        "ms2_spectra": len(spectra or []),
        "linked_ms2_spectra": linked,
        "linked_fraction_pct": (
            round(100.0 * linked / len(spectra), 3) if spectra else None
        ),
        "ambiguous_ms2_spectra": ambiguous,
        "missing_precursor_mz": no_mass,
        "missing_retention_time": no_rt,
        "charge_conflict_candidates_rejected": charge_conflict_candidates,
        "mz_tolerance_ppm": float(mz_tolerance_ppm),
        "mz_floor_da": float(mz_floor_da),
        "rt_padding_minutes": float(rt_padding_minutes),
    }


def _consensus_peaks(
    members: list[dict],
    precursor_mz: float,
    *,
    bin_width_da: float = 0.02,
    occurrence_fraction: float = 0.30,
    top_n: int = 70,
) -> list[list[float]]:
    config = PeakProcessingConfig(
        preset="network_rich",
        relative_intensity_threshold=0.01,
        precursor_exclusion_da=0.5,
        top_n=top_n,
        upper_mass_bound_enabled=True,
    )
    cleaned = []
    for member in members:
        peaks, _ = process_peaks(
            member.get("raw_peaks") or member.get("query_peaks") or [],
            precursor_mz,
            config=config,
            precursor_charge=member.get("precursor_charge"),
        )
        if not peaks:
            continue
        maximum = max(float(peak[1]) for peak in peaks)
        if maximum <= 0:
            continue
        cleaned.append(
            [(float(mz), float(intensity) / maximum) for mz, intensity in peaks]
        )
    if not cleaned:
        return []
    if len(cleaned) == 1:
        return [[mz, intensity] for mz, intensity in cleaned[0]]

    observations: dict[int, list[tuple[float, float, int]]] = defaultdict(list)
    for member_index, peaks in enumerate(cleaned):
        best_by_bin: dict[int, tuple[float, float]] = {}
        for mz, intensity in peaks:
            bin_index = int(round(mz / bin_width_da))
            previous = best_by_bin.get(bin_index)
            if previous is None or intensity > previous[1]:
                best_by_bin[bin_index] = (mz, intensity)
        for bin_index, (mz, intensity) in best_by_bin.items():
            observations[bin_index].append((mz, intensity, member_index))

    required = max(1, int(math.ceil(len(cleaned) * occurrence_fraction)))
    output = []
    for values in observations.values():
        if len({value[2] for value in values}) < required:
            continue
        intensity_sum = sum(value[1] for value in values)
        mz = (
            sum(value[0] * value[1] for value in values) / intensity_sum
            if intensity_sum > 0
            else float(np.median([value[0] for value in values]))
        )
        output.append([mz, intensity_sum / len(cleaned)])
    output.sort(key=lambda peak: peak[1], reverse=True)
    return sorted(output[:top_n], key=lambda peak: peak[0])


def build_lcms_consensus_spectra(
    results: list[dict],
    *,
    top_n: int = 70,
    bin_width_da: float = 0.02,
    occurrence_fraction: float = 0.30,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build condition-safe spectrum artifacts for linked LC-MS features.

    Until precursor-purity and apex-EIC weights are both measured, the
    architecture requires a representative spectrum rather than a synthetic
    consensus. Unknown acquisition conditions never form a shared group.
    These artifacts are read-only evidence and do not change library ranking.
    """

    grouped: dict[str, list[tuple[int, dict]]] = defaultdict(list)
    for index, row in enumerate(results or []):
        feature_id = str(row.get("lcms_feature_id") or "")
        if feature_id:
            grouped[feature_id].append((index, row))

    artifacts = []
    repeated = 0
    mixed_condition_features = 0
    for feature_id, indexed in sorted(grouped.items()):
        condition_groups: dict[tuple[Any, ...], list[tuple[int, dict]]] = (
            defaultdict(list)
        )
        for index, row in indexed:
            ce = _number(row.get("collision_energy"))
            unit = str(row.get("collision_energy_unit") or "").strip()
            mode = str(row.get("collision_energy_mode") or "").strip().lower()
            activation = str(row.get("activation_method") or "").strip()
            configuration = str(
                row.get("instrument_configuration_id") or ""
            ).strip()
            if (
                ce is None
                or not unit
                or mode not in {"fixed", "single"}
                or not activation
                or not configuration
            ):
                condition_key = (
                    "unqualified",
                    str(row.get("spectrum_id") or index),
                )
            else:
                condition_key = (
                    "recorded_condition",
                    activation,
                    round(ce, 8),
                    unit,
                    mode,
                    configuration,
                    str(row.get("ion_mode") or "").upper(),
                    abs(int(_charge_state(row.get("precursor_charge")) or 0)),
                )
            condition_groups[condition_key].append((index, row))

        if len(condition_groups) > 1:
            mixed_condition_features += 1
        ranked_groups = sorted(
            condition_groups.items(),
            key=lambda item: (
                    -max(len(
                        row.get("query_peaks") or row.get("raw_peaks") or []
                    ) for _, row in item[1]),
                -len(item[1]),
                min(index for index, _ in item[1]),
                repr(item[0]),
            ),
        )
        for condition_index, (condition_key, members) in enumerate(
            ranked_groups, start=1
        ):
            rows = [row for _, row in members]
            if len(rows) > 1:
                repeated += 1
            feature_apex = _number(rows[0].get("lcms_feature_rt_minutes"))
            representative_index, representative = min(
                members,
                key=lambda item: (
                    abs(
                        (_number(item[1].get("retention_time")) or 0.0)
                        - (feature_apex or 0.0)
                    ) if feature_apex is not None else 0.0,
                    -len(
                        item[1].get("query_peaks")
                        or item[1].get("raw_peaks")
                        or []
                    ),
                    item[0],
                ),
            )
            peaks = [
                [float(peak[0]), float(peak[1])]
                for peak in (
                    representative.get("query_peaks")
                    or representative.get("raw_peaks")
                    or []
                )[:top_n]
                if len(peak) >= 2
            ]
            if not peaks:
                continue
            artifact_id = f"RS:{feature_id}:{condition_index}"
            artifact = {
                "consensus_spectrum_id": artifact_id,
                "spectrum_artifact_type": "representative_spectrum",
                "lcms_feature_id": feature_id,
                "primary_for_feature": condition_index == 1,
                "precursor_mz": _number(representative.get("precursor_mz")),
                "retention_time_minutes": _number(
                    representative.get("retention_time")
                ),
                "member_count": len(rows),
                "member_spectrum_ids": [
                    row.get("spectrum_id") for row in rows
                ],
                "representative_spectrum_id": representative.get(
                    "spectrum_id"
                ),
                "condition_key": list(condition_key),
                "condition_qualified": condition_key[0] == (
                    "recorded_condition"
                ),
                "peaks": peaks,
                "peak_count": len(peaks),
                "method": (
                    "candidate-independent apex-nearest representative; "
                    "consensus disabled because precursor-purity and EIC "
                    "quality weights are unavailable"
                ),
                "ranking_influence": False,
            }
            artifacts.append(artifact)
            for index, row in members:
                row["lcms_consensus_spectrum_id"] = artifact_id
                row["lcms_consensus_member_count"] = len(rows)
                row["lcms_consensus_peak_count"] = len(peaks)
                row["lcms_consensus_is_representative"] = (
                    index == representative_index
                )
                if index == representative_index:
                    row["lcms_consensus_peaks"] = peaks

    return artifacts, {
        "method": "condition-safe feature-level representative MS2",
        "linked_feature_groups": len(grouped),
        "consensus_spectra": 0,
        "actual_consensus_spectra": 0,
        "representative_spectra": len(artifacts),
        "multi_spectrum_consensus": 0,
        "single_spectrum_consensus": 0,
        "multi_spectrum_representative_groups": repeated,
        "single_spectrum_representatives": len(artifacts) - repeated,
        "mixed_condition_features": mixed_condition_features,
        "unknown_conditions_never_merged": True,
        "ranking_influence": False,
        "top_n": int(top_n),
        "bin_width_da": None,
        "minimum_occurrence_fraction": None,
    }
