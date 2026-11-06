"""Contrastive empirical expert evidence (shadow-only).

The candidate rules were mined with target-class prevalence, background
prevalence and compound-level counting.  They are more informative than raw
within-class frequency, but source-family overlap and unitless collision energy
mean they are not yet calibrated for production ranking.  This module exposes
the evidence without changing any score, rank or confidence label.
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

from .fragmentation_rules import normalize_classification

_RULE_PATH = (
    Path(__file__).resolve().parent.parent
    / "kb_v3"
    / "empirical_rules_candidates_v2.json"
)
_NON_CLASSES = {"", "unknown", "other", "unclassified", "none", "nan"}


@lru_cache(maxsize=1)
def _load_rules() -> dict:
    try:
        payload = json.loads(_RULE_PATH.read_text(encoding="utf-8"))
    except Exception:
        return {}
    if not isinstance(payload, dict):
        return {}
    return payload


def _polarity(value: str | None, adduct: str | None) -> str | None:
    text = str(value or "").strip().lower()
    if text in {"pos", "positive", "+"}:
        return "pos"
    if text in {"neg", "negative", "-"}:
        return "neg"
    try:
        from utils.adduct import adduct_polarity
        return {"POS": "pos", "NEG": "neg"}.get(adduct_polarity(adduct))
    except Exception:
        return None


def _closest_match(values: list[float], target: float, tolerance: float) -> float | None:
    if not values:
        return None
    closest = min(values, key=lambda value: abs(value - target))
    return closest if abs(closest - target) <= tolerance else None


def _evidence_weight(rule: dict) -> float:
    """Bounded descriptive weight; never consumed by production ranking."""
    likelihood_ratio = max(float(rule.get("likelihood_ratio") or 1.0), 1.0)
    prevalence = max(float(
        rule.get("prevalence")
        or rule.get("prevalence_within_adduct")
        or 0.0
    ), 0.0)
    source_families = max(int(rule.get("n_source_families") or 0), 0)
    dominance = min(max(float(rule.get("max_source_dominance") or 1.0), 0.0), 1.0)
    return (
        math.log1p(likelihood_ratio)
        * prevalence
        * math.log1p(source_families)
        * max(1.0 - dominance, 0.05)
    )


def evaluate_contrastive_shadow(
    query_peaks: list[list[float]],
    query_pmz: float,
    candidate_classification: str,
    *,
    query_polarity: str | None = None,
    candidate_adduct: str | None = None,
    query_charge: int | None = None,
    mz_tolerance: float = 0.03,
    relative_intensity_floor: float = 0.01,
) -> dict[str, Any]:
    """Return class-contrastive evidence with an explicit zero-influence guard."""
    payload = _load_rules()
    classes = payload.get("classes") or {}
    normalized_class = normalize_classification(candidate_classification)
    normalized_class = str(normalized_class or "").strip().lower()
    polarity = _polarity(query_polarity, candidate_adduct)
    if (
        not payload
        or normalized_class in _NON_CLASSES
        or normalized_class not in classes
        or polarity not in {"pos", "neg"}
    ):
        return {
            "status": "not_applicable",
            "candidate_class": normalized_class or "unknown",
            "polarity": polarity,
            "runtime_influence": False,
            "score_contribution": 0.0,
            "matched_ions": [],
            "matched_losses": [],
        }

    valid_peaks = [
        (float(peak[0]), float(peak[1]))
        for peak in (query_peaks or [])
        if len(peak) >= 2 and float(peak[0]) > 0 and float(peak[1]) > 0
    ]
    base_peak = max((intensity for _, intensity in valid_peaks), default=0.0)
    mz_values = [
        mz for mz, intensity in valid_peaks
        if base_peak <= 0 or intensity >= base_peak * relative_intensity_floor
    ]
    rules = (classes.get(normalized_class) or {}).get(polarity) or {}
    selected_rules = rules
    adduct_specific = False
    adduct_text = str(candidate_adduct or "").strip()
    if adduct_text:
        stratum = (rules.get("adduct_strata") or {}).get(adduct_text)
        if stratum and int(stratum.get("n_independent_compounds") or 0) >= 20:
            selected_rules = {
                **rules,
                "diagnostic_ions": stratum.get("diagnostic_ions") or [],
                "neutral_losses": stratum.get("neutral_losses") or [],
            }
            adduct_specific = True

    matched_ions = []
    for rule in selected_rules.get("diagnostic_ions") or []:
        target = float(rule.get("mz") or 0.0)
        observed = _closest_match(mz_values, target, mz_tolerance)
        if observed is not None:
            matched_ions.append({
                "expected_mz": round(target, 6),
                "observed_mz": round(observed, 6),
                "prevalence": rule.get("prevalence") or rule.get(
                    "prevalence_within_adduct"
                ),
                "background_prevalence": rule.get("background_prevalence"),
                "likelihood_ratio": rule.get("likelihood_ratio"),
                "n_independent_compounds": rule.get("n_independent_compounds"),
                "n_source_families": rule.get("n_source_families"),
                "max_source_dominance": rule.get("max_source_dominance"),
                "descriptive_weight": round(_evidence_weight(rule), 6),
            })

    try:
        multiply_charged = abs(int(query_charge or 1)) > 1
    except (TypeError, ValueError):
        multiply_charged = False
    observed_losses = (
        [float(query_pmz) - mz for mz in mz_values if 0 < float(query_pmz) - mz <= 400]
        if query_pmz and not multiply_charged else []
    )
    matched_losses = []
    for rule in selected_rules.get("neutral_losses") or []:
        target = float(rule.get("loss_da") or 0.0)
        observed = _closest_match(observed_losses, target, mz_tolerance)
        if observed is not None:
            matched_losses.append({
                "expected_loss_da": round(target, 6),
                "observed_loss_da": round(observed, 6),
                "prevalence": rule.get("prevalence") or rule.get(
                    "prevalence_within_adduct"
                ),
                "background_prevalence": rule.get("background_prevalence"),
                "likelihood_ratio": rule.get("likelihood_ratio"),
                "n_independent_compounds": rule.get("n_independent_compounds"),
                "n_source_families": rule.get("n_source_families"),
                "max_source_dominance": rule.get("max_source_dominance"),
                "descriptive_weight": round(_evidence_weight(rule), 6),
            })

    evidence_weight = sum(
        item["descriptive_weight"] for item in matched_ions + matched_losses
    )
    return {
        "status": "shadow_only_uncalibrated",
        "candidate_class": normalized_class,
        "polarity": polarity,
        "adduct_specific": adduct_specific,
        "runtime_influence": False,
        "score_contribution": 0.0,
        "descriptive_evidence_weight": round(evidence_weight, 6),
        "n_class_compounds": int(rules.get("n_independent_compounds") or 0),
        "matched_ions": sorted(
            matched_ions,
            key=lambda item: item["descriptive_weight"],
            reverse=True,
        )[:10],
        "matched_losses": sorted(
            matched_losses,
            key=lambda item: item["descriptive_weight"],
            reverse=True,
        )[:10],
        "limitations": [
            "Source labels are not guaranteed independent because public libraries republish spectra.",
            "Collision-energy values in the mined database are unitless and were not used.",
            "Evidence is descriptive until compound-split target-decoy calibration is complete.",
        ],
        "rule_set": "empirical_rules_candidates_v2",
        "rule_rows_used": int((payload.get("meta") or {}).get("rows_used") or 0),
    }


EXPERT_CONTRASTIVE_AVAILABLE = True
