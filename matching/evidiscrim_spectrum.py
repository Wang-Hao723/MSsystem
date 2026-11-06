"""Versioned matchms spectrum semantics for the EviDiscrim benchmark."""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

import matchms
import numpy as np
from matchms import Spectrum
from matchms.filtering import (
    normalize_intensities,
    require_minimum_number_of_peaks,
    select_by_relative_intensity,
)
from matchms.similarity import CosineHungarian


MATCHMS_PIPELINE_ID = "evidiscrim-matchms-v1"


def _number(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _coalesce_peaks(peaks: Any) -> tuple[np.ndarray, np.ndarray]:
    grouped: dict[float, float] = defaultdict(float)
    if not isinstance(peaks, Sequence) or isinstance(peaks, (str, bytes)):
        return np.array([], dtype=float), np.array([], dtype=float)
    for raw in peaks:
        if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
            continue
        if len(raw) < 2:
            continue
        mz = _number(raw[0])
        intensity = _number(raw[1])
        if mz is None or intensity is None or mz <= 0 or intensity <= 0:
            continue
        grouped[round(mz, 8)] += intensity
    ordered = sorted(grouped.items())
    return (
        np.asarray([item[0] for item in ordered], dtype=float),
        np.asarray([item[1] for item in ordered], dtype=float),
    )


def prepare_matchms_spectrum(
    peaks: Any,
    *,
    precursor_mz: float | None = None,
    minimum_relative_intensity: float = 0.01,
    minimum_peak_count: int = 2,
) -> Spectrum | None:
    """Apply the frozen benchmark peak pipeline in a fixed order."""

    mz, intensities = _coalesce_peaks(peaks)
    if not len(mz):
        return None
    metadata: dict[str, Any] = {
        "processing_pipeline": MATCHMS_PIPELINE_ID,
    }
    if precursor_mz is not None and math.isfinite(float(precursor_mz)):
        metadata["precursor_mz"] = float(precursor_mz)
    spectrum = Spectrum(mz=mz, intensities=intensities, metadata=metadata)
    spectrum = normalize_intensities(spectrum)
    spectrum = select_by_relative_intensity(
        spectrum,
        intensity_from=max(0.0, float(minimum_relative_intensity)),
        intensity_to=1.0,
    )
    spectrum = require_minimum_number_of_peaks(
        spectrum, n_required=max(1, int(minimum_peak_count))
    )
    return spectrum


def processed_peak_pairs(spectrum: Spectrum | None) -> list[list[float]]:
    if spectrum is None:
        return []
    return [
        [round(float(mz), 8), round(float(intensity), 10)]
        for mz, intensity in zip(
            spectrum.peaks.mz,
            spectrum.peaks.intensities,
            strict=True,
        )
    ]


def spectrum_fingerprint(spectrum: Spectrum | None) -> str:
    payload = {
        "pipeline_id": MATCHMS_PIPELINE_ID,
        "matchms_version": matchms.__version__,
        "peaks": processed_peak_pairs(spectrum),
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def whole_spectrum_hungarian_cosine(
    observed_peaks: Any,
    predicted_peaks: Any,
    *,
    precursor_mz: float | None,
    mz_tolerance_da: float = 0.02,
    minimum_relative_intensity: float = 0.01,
    minimum_peak_count: int = 2,
) -> dict[str, Any]:
    observed = prepare_matchms_spectrum(
        observed_peaks,
        precursor_mz=precursor_mz,
        minimum_relative_intensity=minimum_relative_intensity,
        minimum_peak_count=minimum_peak_count,
    )
    predicted = prepare_matchms_spectrum(
        predicted_peaks,
        precursor_mz=precursor_mz,
        minimum_relative_intensity=minimum_relative_intensity,
        minimum_peak_count=minimum_peak_count,
    )
    if observed is None or predicted is None:
        return {
            "status": "unavailable_insufficient_processed_peaks",
            "score": None,
            "matched_peaks": 0,
            "pipeline_id": MATCHMS_PIPELINE_ID,
            "matchms_version": matchms.__version__,
        }
    score = CosineHungarian(tolerance=float(mz_tolerance_da)).pair(
        predicted, observed
    )
    return {
        "status": "scored",
        "score": round(float(score["score"]), 10),
        "matched_peaks": int(score["matches"]),
        "pipeline_id": MATCHMS_PIPELINE_ID,
        "matchms_version": matchms.__version__,
        "mz_tolerance_da": float(mz_tolerance_da),
        "minimum_relative_intensity": float(minimum_relative_intensity),
        "observed_peak_count": len(observed.peaks.mz),
        "predicted_peak_count": len(predicted.peaks.mz),
    }


def aligned_predicted_observed_intensities(
    predicted: Spectrum,
    observed: Spectrum,
    *,
    mz_tolerance_da: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Map every predicted peak to at most one observed peak for calibration."""

    candidates: list[tuple[float, float, int, int]] = []
    for pred_index, (pred_mz, pred_intensity) in enumerate(zip(
        predicted.peaks.mz, predicted.peaks.intensities, strict=True
    )):
        for obs_index, (obs_mz, obs_intensity) in enumerate(zip(
            observed.peaks.mz, observed.peaks.intensities, strict=True
        )):
            distance = abs(float(pred_mz) - float(obs_mz))
            if distance <= mz_tolerance_da:
                candidates.append((
                    distance,
                    -float(pred_intensity) * float(obs_intensity),
                    pred_index,
                    obs_index,
                ))
    assigned_pred: set[int] = set()
    assigned_obs: set[int] = set()
    targets = np.zeros(len(predicted.peaks.mz), dtype=float)
    for _, _, pred_index, obs_index in sorted(candidates):
        if pred_index in assigned_pred or obs_index in assigned_obs:
            continue
        assigned_pred.add(pred_index)
        assigned_obs.add(obs_index)
        targets[pred_index] = float(observed.peaks.intensities[obs_index])
    return np.asarray(predicted.peaks.intensities, dtype=float), targets


def spectrum_contract(policy: Mapping[str, Any]) -> dict[str, Any]:
    processing = policy.get("spectrum_processing") or {}
    return {
        "pipeline_id": MATCHMS_PIPELINE_ID,
        "matchms_version": matchms.__version__,
        "normalization": "base_peak_to_one",
        "peak_matching": "CosineHungarian_optimal_one_to_one",
        "mz_tolerance_da": float(processing.get("mz_tolerance_da", 0.02)),
        "minimum_relative_intensity": float(
            processing.get("minimum_relative_intensity", 0.01)
        ),
        "minimum_peak_count": int(processing.get("minimum_peak_count", 2)),
        "precursor_peak_removal": False,
        "modified_cosine_mass_shift": False,
    }


__all__ = [
    "MATCHMS_PIPELINE_ID",
    "aligned_predicted_observed_intensities",
    "prepare_matchms_spectrum",
    "processed_peak_pairs",
    "spectrum_contract",
    "spectrum_fingerprint",
    "whole_spectrum_hungarian_cosine",
]
