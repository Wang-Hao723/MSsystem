"""Read-only multi-metric spectral-similarity evidence for public Top-K rows.

The production search score is intentionally left untouched.  This module
computes several independently interpretable matchms similarities only after
candidate selection so users can see whether a match is supported by direct
product ions, precursor-shifted ions, or neutral-loss patterns.
"""

from __future__ import annotations

import math
from typing import Any, Iterable

import numpy as np
from matchms import Spectrum
from matchms.similarity import (
    CosineHungarian,
    ModifiedCosineHungarian,
    NeutralLossesCosine,
)
from utils.spectrum import process_peaks, resolve_peak_processing_config


PEAK_VIEW_PRESETS = ("library_match", "sensitive", "network_rich")


def _spectrum(
    peaks: Iterable[Iterable[Any]],
    precursor_mz: Any,
) -> tuple[Spectrum | None, int]:
    """Return a valid, sorted matchms spectrum and its retained peak count."""
    cleaned: list[tuple[float, float]] = []
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
        ):
            cleaned.append((mz, intensity))
    if not cleaned:
        return None, 0

    # matchms requires sorted m/z values.  Coalesce exact duplicates so a
    # malformed legacy spectrum cannot make the evidence layer fail.
    by_mz: dict[float, float] = {}
    for mz, intensity in cleaned:
        by_mz[mz] = by_mz.get(mz, 0.0) + intensity
    ordered = sorted(by_mz.items())
    try:
        pmz = float(precursor_mz)
    except (TypeError, ValueError):
        pmz = 0.0
    if not math.isfinite(pmz) or pmz <= 0.0:
        return None, len(ordered)
    return (
        Spectrum(
            mz=np.asarray([item[0] for item in ordered], dtype=float),
            intensities=np.asarray([item[1] for item in ordered], dtype=float),
            metadata={"precursor_mz": pmz},
        ),
        len(ordered),
    )


def _metric(similarity: Any, query: Spectrum, library: Spectrum) -> dict[str, Any]:
    """Convert matchms' structured scalar into stable JSON evidence."""
    try:
        value = similarity.pair(query, library)
        score = float(value["score"])
        matches = int(value["matches"])
        if not math.isfinite(score):
            raise ValueError("non-finite similarity")
        return {
            "status": "available",
            "score": round(max(0.0, min(1.0, score)), 6),
            "matched_peaks": max(0, matches),
        }
    except Exception as exc:
        return {
            "status": "evaluation_failed",
            "score": None,
            "matched_peaks": None,
            "error": str(exc)[:160],
        }


def build_multimetric_shadow(
    query_peaks: Iterable[Iterable[Any]],
    library_peaks: Iterable[Iterable[Any]],
    *,
    query_precursor_mz: Any,
    library_precursor_mz: Any,
    production_modified_cosine: Any = None,
    spectral_entropy: Any = None,
    reverse_coverage: Any = None,
    fragment_tolerance_da: float = 0.10,
    strict_fragment_tolerance_da: float = 0.02,
) -> dict[str, Any]:
    """Build shadow-only similarity evidence for one selected candidate."""
    evidence: dict[str, Any] = {
        "schema_version": 1,
        "status": "unavailable",
        "method": "matchms Top-K multi-metric shadow audit",
        "runtime_influence": False,
        "score_contribution": 0.0,
        "ranking_changed": False,
        "fragment_tolerance_da": float(fragment_tolerance_da),
        "strict_fragment_tolerance_da": float(strict_fragment_tolerance_da),
        "production_modified_cosine": production_modified_cosine,
        "spectral_entropy": spectral_entropy,
        "reverse_coverage": reverse_coverage,
    }
    query, query_n = _spectrum(query_peaks, query_precursor_mz)
    library, library_n = _spectrum(library_peaks, library_precursor_mz)
    evidence["query_peak_count"] = query_n
    evidence["library_peak_count"] = library_n
    if query is None or library is None:
        evidence["status"] = "insufficient_spectrum_or_precursor"
        return evidence

    tol = max(0.001, float(fragment_tolerance_da))
    strict_tol = max(0.001, min(tol, float(strict_fragment_tolerance_da)))
    evidence.update({
        "status": "available",
        "direct_cosine_hungarian": _metric(
            CosineHungarian(tolerance=tol), query, library
        ),
        "modified_cosine_hungarian": _metric(
            ModifiedCosineHungarian(tolerance=tol), query, library
        ),
        "strict_modified_cosine_hungarian": _metric(
            ModifiedCosineHungarian(tolerance=strict_tol), query, library
        ),
        "neutral_loss_cosine": _metric(
            NeutralLossesCosine(tolerance=tol), query, library
        ),
    })

    score_values: list[float] = []
    for key in (
        "production_modified_cosine",
        "spectral_entropy",
        "reverse_coverage",
    ):
        try:
            value = float(evidence.get(key))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            score_values.append(value)
    for key in (
        "direct_cosine_hungarian",
        "modified_cosine_hungarian",
        "neutral_loss_cosine",
    ):
        try:
            value = float(evidence[key]["score"])
        except (TypeError, ValueError, KeyError):
            continue
        if math.isfinite(value):
            score_values.append(value)
    evidence["metric_spread"] = (
        round(max(score_values) - min(score_values), 6)
        if len(score_values) >= 2 else None
    )
    evidence["interpretation"] = (
        "Diagnostic comparison only; no metric in this object changes "
        "production score, ranking, or confidence."
    )
    return evidence


def build_peak_view_shadow(
    query_raw_peaks: Iterable[Iterable[Any]],
    library_raw_peaks: Iterable[Iterable[Any]],
    *,
    query_precursor_mz: Any,
    library_precursor_mz: Any,
    query_precursor_charge: Any = None,
    library_precursor_charge: Any = None,
    fragment_tolerance_da: float = 0.10,
) -> dict[str, Any]:
    """Compare the same candidate across three symmetric peak views.

    The production Top-15/5% view is intentionally preserved.  The Top-40
    and Top-70 views expose whether discarded low-abundance ions would support
    or weaken the candidate, but they cannot change candidate selection or a
    confidence tier until an identity-disjoint calibration approves a policy.
    """
    result: dict[str, Any] = {
        "schema_version": 1,
        "status": "unavailable",
        "method": "symmetric query/library peak-view shadow",
        "runtime_influence": False,
        "score_contribution": 0.0,
        "ranking_changed": False,
        "production_view": "library_match",
        "views": {},
    }
    tolerance = max(0.001, float(fragment_tolerance_da))
    available_scores: dict[str, float] = {}
    for preset in PEAK_VIEW_PRESETS:
        config = resolve_peak_processing_config(preset)
        query_peaks, query_audit = process_peaks(
            query_raw_peaks,
            query_precursor_mz,
            config=config,
            precursor_charge=query_precursor_charge,
        )
        library_peaks, library_audit = process_peaks(
            library_raw_peaks,
            library_precursor_mz,
            config=config,
            precursor_charge=library_precursor_charge,
        )
        query, query_n = _spectrum(query_peaks, query_precursor_mz)
        library, library_n = _spectrum(library_peaks, library_precursor_mz)
        if query is None or library is None:
            metric = {
                "status": "insufficient_spectrum_or_precursor",
                "score": None,
                "matched_peaks": None,
            }
        else:
            metric = _metric(
                ModifiedCosineHungarian(tolerance=tolerance),
                query,
                library,
            )
            if metric.get("score") is not None:
                available_scores[preset] = float(metric["score"])
        result["views"][preset] = {
            "config": config.as_dict(),
            "query_peak_count": query_n,
            "library_peak_count": library_n,
            "modified_cosine_hungarian": metric,
            "query_retained_pct": query_audit.get("retained_pct"),
            "library_retained_pct": library_audit.get("retained_pct"),
        }

    if available_scores:
        result["status"] = "available"
        best_view = max(
            PEAK_VIEW_PRESETS,
            key=lambda key: available_scores.get(key, float("-inf")),
        )
        production_score = available_scores.get("library_match")
        result["highest_similarity_view"] = best_view
        result["highest_similarity"] = round(available_scores[best_view], 6)
        result["sensitive_gain_over_production_view"] = (
            round(
                available_scores["sensitive"] - production_score,
                6,
            )
            if production_score is not None and "sensitive" in available_scores
            else None
        )
        result["network_rich_gain_over_production_view"] = (
            round(
                available_scores["network_rich"] - production_score,
                6,
            )
            if production_score is not None and "network_rich" in available_scores
            else None
        )
        result["view_score_range"] = (
            round(max(available_scores.values()) - min(available_scores.values()), 6)
            if len(available_scores) > 1 else 0.0
        )
    result["interpretation"] = (
        "Counterfactual peak-retention audit only. Broader views are not "
        "calibrated to rerank candidates or promote High evidence."
    )
    return result


__all__ = [
    "PEAK_VIEW_PRESETS",
    "build_multimetric_shadow",
    "build_peak_view_shadow",
]
