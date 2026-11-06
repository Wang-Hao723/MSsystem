"""Condition-aware pair features for same-formula MS/MS discrimination.

These features intentionally preserve two forms of isomer evidence that a
single broad Modified Cosine score can blur:

* strict versus broad product-ion agreement;
* strict versus broad neutral-loss agreement;
* relative-intensity reproducibility among matched ions, including cases with
  no candidate-exclusive m/z value.

The module is dependency-light and deterministic so the same implementation is
used for offline training, locked evaluation, and production inference.
"""

from __future__ import annotations

import math
from typing import Any, Iterable

import numpy as np


TOLERANCES = (0.01, 0.02, 0.05, 0.10)


def _metric_names(prefix: str, tolerance: float) -> tuple[str, ...]:
    suffix = f"{tolerance:.2f}".replace(".", "p")
    return (
        f"{prefix}_cosine_{suffix}",
        f"{prefix}_matches_{suffix}",
        f"{prefix}_query_coverage_{suffix}",
        f"{prefix}_reference_coverage_{suffix}",
        f"{prefix}_matched_intensity_corr_{suffix}",
        f"{prefix}_log_intensity_mae_{suffix}",
    )


FEATURE_NAMES = (
    "query_peak_count",
    "reference_peak_count",
    "peak_count_ratio",
    "query_entropy",
    "reference_entropy",
    "entropy_abs_delta",
    "collision_energy_comparable",
    "collision_energy_delta_scaled",
    "instrument_family_match",
    *(
        name
        for prefix in ("product", "neutral_loss")
        for tolerance in TOLERANCES
        for name in _metric_names(prefix, tolerance)
    ),
    "modified_cosine_0p10",
    "strict_product_fraction",
    "strict_loss_fraction",
    "product_loss_score_agreement",
    "strict_broad_intensity_stability",
)


def _finite_number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def prepare_spectrum(
    peaks: Iterable[Iterable[Any]] | None,
    precursor_mz: Any,
    *,
    relative_intensity_threshold: float = 0.01,
    top_n: int = 100,
    precursor_exclusion_da: float = 0.5,
) -> np.ndarray:
    """Return a rich, normalized product-ion view for pair discrimination."""
    precursor = _finite_number(precursor_mz) or 0.0
    by_mz: dict[float, float] = {}
    for peak in peaks or []:
        try:
            mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if (
            math.isfinite(mz)
            and math.isfinite(intensity)
            and mz > 0.0
            and intensity > 0.0
            and (precursor <= 0.0 or mz <= precursor + 1.5)
            and (
                precursor <= 0.0
                or abs(mz - precursor) > precursor_exclusion_da
            )
        ):
            # Coalescing exact duplicates prevents a malformed source spectrum
            # from contributing the same ion more than once.
            by_mz[mz] = by_mz.get(mz, 0.0) + intensity
    if not by_mz:
        return np.empty((0, 2), dtype=np.float64)
    spectrum = np.asarray(sorted(by_mz.items()), dtype=np.float64)
    maximum = float(spectrum[:, 1].max())
    spectrum = spectrum[
        spectrum[:, 1] >= maximum * max(0.0, relative_intensity_threshold)
    ]
    if len(spectrum) > top_n:
        strongest = np.argsort(spectrum[:, 1], kind="stable")[::-1][:top_n]
        spectrum = spectrum[strongest]
        spectrum = spectrum[np.argsort(spectrum[:, 0], kind="stable")]
    norm = float(np.linalg.norm(spectrum[:, 1]))
    if norm > 0.0:
        spectrum[:, 1] /= norm
    return spectrum


def _entropy(intensities: np.ndarray) -> float:
    if intensities.size == 0:
        return 0.0
    total = float(intensities.sum())
    if total <= 0.0:
        return 0.0
    probabilities = intensities / total
    value = -float(np.sum(probabilities * np.log(probabilities + 1e-12)))
    maximum = math.log(len(probabilities)) if len(probabilities) > 1 else 1.0
    return value / maximum if maximum > 0.0 else 0.0


def _candidate_pairs(
    left_mz: np.ndarray,
    left_intensity: np.ndarray,
    right_mz: np.ndarray,
    right_intensity: np.ndarray,
    tolerance: float,
    shifts: tuple[float, ...] = (0.0,),
) -> list[tuple[float, int, int]]:
    candidates: list[tuple[float, int, int]] = []
    for shift in shifts:
        differences = np.abs(left_mz[:, None] - (right_mz[None, :] + shift))
        left_indices, right_indices = np.where(differences <= tolerance)
        candidates.extend(
            (
                float(left_intensity[left] * right_intensity[right]),
                int(left),
                int(right),
            )
            for left, right in zip(left_indices, right_indices)
        )
    # Stable index suffixes make equal-intensity ties reproducible.
    candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
    return candidates


def _greedy_matches(
    left: np.ndarray,
    right: np.ndarray,
    tolerance: float,
    *,
    shifts: tuple[float, ...] = (0.0,),
) -> list[tuple[int, int]]:
    if len(left) == 0 or len(right) == 0:
        return []
    used_left: set[int] = set()
    used_right: set[int] = set()
    matches: list[tuple[int, int]] = []
    for _, left_index, right_index in _candidate_pairs(
        left[:, 0], left[:, 1], right[:, 0], right[:, 1], tolerance, shifts
    ):
        if left_index in used_left or right_index in used_right:
            continue
        used_left.add(left_index)
        used_right.add(right_index)
        matches.append((left_index, right_index))
    return matches


def _pearson(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 2 or len(right) < 2:
        return 0.0
    left_std = float(np.std(left))
    right_std = float(np.std(right))
    if left_std <= 1e-12 or right_std <= 1e-12:
        return 0.0
    value = float(np.corrcoef(left, right)[0, 1])
    return value if math.isfinite(value) else 0.0


def _match_metrics(
    left: np.ndarray,
    right: np.ndarray,
    tolerance: float,
    *,
    shifts: tuple[float, ...] = (0.0,),
) -> tuple[float, float, float, float, float, float]:
    matches = _greedy_matches(left, right, tolerance, shifts=shifts)
    if not matches:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 4.0
    left_indices = np.asarray([item[0] for item in matches], dtype=int)
    right_indices = np.asarray([item[1] for item in matches], dtype=int)
    left_intensity = left[left_indices, 1]
    right_intensity = right[right_indices, 1]
    score = float(np.sum(left_intensity * right_intensity))
    query_coverage = float(np.sum(left_intensity * left_intensity))
    reference_coverage = float(np.sum(right_intensity * right_intensity))
    correlation = _pearson(
        np.log1p(left_intensity * 100.0),
        np.log1p(right_intensity * 100.0),
    )
    left_relative = left_intensity / max(float(left_intensity.sum()), 1e-12)
    right_relative = right_intensity / max(float(right_intensity.sum()), 1e-12)
    log_mae = float(np.mean(np.abs(
        np.log(left_relative + 1e-6) - np.log(right_relative + 1e-6)
    )))
    return (
        max(0.0, min(1.0, score)),
        float(len(matches)),
        max(0.0, min(1.0, query_coverage)),
        max(0.0, min(1.0, reference_coverage)),
        max(-1.0, min(1.0, correlation)),
        min(log_mae, 4.0),
    )


def _neutral_losses(spectrum: np.ndarray, precursor_mz: float) -> np.ndarray:
    if len(spectrum) == 0 or precursor_mz <= 0.0:
        return np.empty((0, 2), dtype=np.float64)
    losses = precursor_mz - spectrum[:, 0]
    keep = (losses > 0.0) & (losses <= 500.0)
    output = np.column_stack((losses[keep], spectrum[keep, 1]))
    return output[np.argsort(output[:, 0], kind="stable")]


def spectral_pair_feature_dict(
    query_peaks: Iterable[Iterable[Any]] | None,
    reference_peaks: Iterable[Iterable[Any]] | None,
    *,
    query_precursor_mz: Any,
    reference_precursor_mz: Any,
    query_collision_energy: Any = None,
    reference_collision_energy: Any = None,
    query_collision_energy_unit: str | None = None,
    reference_collision_energy_unit: str | None = None,
    query_instrument_family: str | None = None,
    reference_instrument_family: str | None = None,
) -> dict[str, float]:
    query_pmz = _finite_number(query_precursor_mz) or 0.0
    reference_pmz = _finite_number(reference_precursor_mz) or 0.0
    query = prepare_spectrum(query_peaks, query_pmz)
    reference = prepare_spectrum(reference_peaks, reference_pmz)
    query_losses = _neutral_losses(query, query_pmz)
    reference_losses = _neutral_losses(reference, reference_pmz)
    features: dict[str, float] = {
        "query_peak_count": min(len(query) / 100.0, 1.0),
        "reference_peak_count": min(len(reference) / 100.0, 1.0),
        "peak_count_ratio": (
            min(len(query), len(reference)) / max(len(query), len(reference))
            if len(query) and len(reference) else 0.0
        ),
        "query_entropy": _entropy(query[:, 1]),
        "reference_entropy": _entropy(reference[:, 1]),
        "entropy_abs_delta": abs(
            _entropy(query[:, 1]) - _entropy(reference[:, 1])
        ),
    }
    query_ce = _finite_number(query_collision_energy)
    reference_ce = _finite_number(reference_collision_energy)
    query_ce_unit = str(query_collision_energy_unit or "").strip().lower()
    reference_ce_unit = str(reference_collision_energy_unit or "").strip().lower()
    ce_comparable = bool(
        query_ce is not None
        and reference_ce is not None
        and query_ce_unit
        and query_ce_unit == reference_ce_unit
    )
    features["collision_energy_comparable"] = float(ce_comparable)
    features["collision_energy_delta_scaled"] = (
        min(abs(query_ce - reference_ce) / 50.0, 4.0)
        if ce_comparable else 0.0
    )
    query_instrument = str(query_instrument_family or "").strip().lower()
    reference_instrument = str(reference_instrument_family or "").strip().lower()
    features["instrument_family_match"] = float(
        bool(query_instrument)
        and query_instrument == reference_instrument
    )

    for prefix, left, right in (
        ("product", query, reference),
        ("neutral_loss", query_losses, reference_losses),
    ):
        for tolerance in TOLERANCES:
            values = _match_metrics(left, right, tolerance)
            for name, value in zip(_metric_names(prefix, tolerance), values):
                features[name] = value

    shift = query_pmz - reference_pmz
    shifts = (0.0,) if abs(shift) <= 0.10 else (0.0, shift)
    features["modified_cosine_0p10"] = _match_metrics(
        query, reference, 0.10, shifts=shifts
    )[0]
    product_strict = features["product_cosine_0p02"]
    product_broad = features["product_cosine_0p10"]
    loss_strict = features["neutral_loss_cosine_0p02"]
    loss_broad = features["neutral_loss_cosine_0p10"]
    features["strict_product_fraction"] = product_strict / max(product_broad, 1e-6)
    features["strict_loss_fraction"] = loss_strict / max(loss_broad, 1e-6)
    features["product_loss_score_agreement"] = 1.0 - min(
        abs(product_broad - loss_broad), 1.0
    )
    features["strict_broad_intensity_stability"] = 1.0 - min(
        abs(features["product_log_intensity_mae_0p02"]
            - features["product_log_intensity_mae_0p10"]) / 4.0,
        1.0,
    )
    return {name: float(features.get(name, 0.0)) for name in FEATURE_NAMES}


def spectral_pair_features(*args: Any, **kwargs: Any) -> np.ndarray:
    values = spectral_pair_feature_dict(*args, **kwargs)
    return np.asarray([values[name] for name in FEATURE_NAMES], dtype=np.float32)
