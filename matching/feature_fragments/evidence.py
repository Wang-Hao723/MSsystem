"""Conditional likelihood-ratio evidence for feature fragments.

No fixed point awards are used.  A present or eligible absent observation is
converted to a conditional natural-log likelihood ratio using Jeffreys-smoothed
target and background frequencies.  Correlated rules share an independent
fragmentation group and contribute at most one representative observation.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


JEFFREYS_ALPHA = 0.5
JEFFREYS_BETA = 0.5


def feature_fragment_score_standard() -> dict[str, Any]:
    """Return the stable, non-probabilistic feature-fragment score contract."""

    return {
        "schema_version": "mssystem-v15-feature-fragment-score/1.0",
        "primary_unit": "natural_log_likelihood_ratio_nats",
        "one_independent_fragment": (
            "ln[P(observation|target_class,acquisition_condition)/"
            "P(observation|background,acquisition_condition)]"
        ),
        "two_or_more_independent_fragments": (
            "sum of the per-group natural-log likelihood ratios"
        ),
        "correlated_rule_handling": (
            "within one independent_fragmentation_group retain only the rule "
            "with the largest absolute natural-log likelihood ratio"
        ),
        "count_alone_has_fixed_points": False,
        "collision_energy_requirement": (
            "same reported unit, same fixed mode, and query value inside the "
            "rule minimum/maximum; no cross-unit or cross-instrument conversion"
        ),
        "absence_requirement": (
            "negative evidence is scored only when the rule authorizes absence "
            "and every declared spectrum-quality field is true"
        ),
        "score_is_probability": False,
        "identity_rank_influence": False,
        "production_score_contribution": 0.0,
    }


def _score_units(log_likelihood_ratio: float) -> dict[str, float]:
    value = float(log_likelihood_ratio)
    return {
        "score_nats": value,
        "score_bits": value / math.log(2.0),
        "score_decibans": 10.0 * value / math.log(10.0),
    }


def _rule_standard(rule: Mapping[str, Any]) -> dict[str, Any]:
    expected = (
        {"product_mz": rule.get("product_mz")}
        if rule.get("observation_type") == "product_ion"
        else {"neutral_loss_da": rule.get("neutral_loss_da")}
    )
    return {
        "observation_type": rule.get("observation_type"),
        **expected,
        "chemical_classes": list(rule.get("chemical_classes") or []),
        "polarity": rule.get("polarity"),
        "adducts": list(rule.get("adducts") or []),
        "charges": list(rule.get("charges") or []),
        "instrument_domain_ids": list(rule.get("instrument_domain_ids") or []),
        "activation": rule.get("activation"),
        "collision_energy_domain": dict(rule.get("collision_energy") or {}),
        "minimum_relative_intensity": rule.get("minimum_relative_intensity"),
        "mass_tolerance": dict(rule.get("mass_tolerance") or {}),
    }


def mass_window(mz: float, tolerance: Mapping[str, Any]) -> float:
    """Return the rule-specific maximum of the absolute and ppm windows."""
    return max(
        float(tolerance["absolute_da"]),
        abs(float(mz)) * float(tolerance["ppm"]) / 1_000_000.0,
    )


def conditional_log_likelihood_ratio(
    *,
    target_present: int,
    target_total: int,
    background_present: int,
    background_total: int,
    observed_present: bool,
    alpha: float = JEFFREYS_ALPHA,
    beta: float = JEFFREYS_BETA,
) -> dict[str, float]:
    """Calculate a smoothed conditional LR from independent identity counts."""
    if not (0 <= target_present <= target_total):
        raise ValueError("invalid target counts")
    if not (0 <= background_present <= background_total):
        raise ValueError("invalid background counts")
    if target_total <= 0 or background_total <= 0:
        raise ValueError("target and background totals must be positive")
    target_probability = (target_present + alpha) / (target_total + alpha + beta)
    background_probability = (background_present + alpha) / (
        background_total + alpha + beta
    )
    if observed_present:
        numerator = target_probability
        denominator = background_probability
    else:
        numerator = 1.0 - target_probability
        denominator = 1.0 - background_probability
    likelihood_ratio = numerator / denominator
    return {
        "target_probability": target_probability,
        "background_probability": background_probability,
        "likelihood_ratio": likelihood_ratio,
        "log_likelihood_ratio": math.log(likelihood_ratio),
    }


def _ce_applicability(
    rule_ce: Mapping[str, Any], query_ce: Mapping[str, Any] | None
) -> tuple[bool, str]:
    if not query_ce:
        return False, "collision_energy_missing"
    if query_ce.get("unit") != rule_ce.get("unit"):
        return False, "collision_energy_unit_mismatch"
    if query_ce.get("mode") != rule_ce.get("mode"):
        return False, "collision_energy_mode_mismatch"
    if rule_ce.get("mode") != "fixed":
        return False, "composite_energy_not_supported_by_fixed_rule"
    value = query_ce.get("value")
    if value is None:
        return False, "collision_energy_value_missing"
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return False, "collision_energy_value_invalid"
    if not float(rule_ce["minimum"]) <= numeric <= float(rule_ce["maximum"]):
        return False, "collision_energy_outside_rule_domain"
    return True, "applicable"


def _normalized_peaks(peaks: Sequence[Sequence[float]]) -> list[tuple[float, float]]:
    cleaned: list[tuple[float, float]] = []
    for peak in peaks:
        if not peak or len(peak) < 2:
            continue
        try:
            mz, intensity = float(peak[0]), float(peak[1])
        except (TypeError, ValueError):
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity >= 0:
            cleaned.append((mz, intensity))
    base = max((intensity for _, intensity in cleaned), default=0.0)
    if base <= 0:
        return []
    return [(mz, intensity / base) for mz, intensity in cleaned]


def _observation(
    rule: Mapping[str, Any],
    peaks: Sequence[tuple[float, float]],
    precursor_mz: float | None,
) -> dict[str, Any]:
    if rule["observation_type"] == "product_ion":
        expected = float(rule["product_mz"])
    else:
        if precursor_mz is None:
            return {"observable": False, "reason": "precursor_mz_missing"}
        expected = float(precursor_mz) - float(rule["neutral_loss_da"])
        if expected <= 0:
            return {"observable": False, "reason": "neutral_loss_exceeds_precursor"}
    window = mass_window(expected, rule["mass_tolerance"])
    intensity_floor = float(rule["minimum_relative_intensity"])
    matches = [
        (mz, relative)
        for mz, relative in peaks
        if relative >= intensity_floor and abs(mz - expected) <= window
    ]
    best = min(matches, key=lambda item: abs(item[0] - expected)) if matches else None
    return {
        "observable": True,
        "present": best is not None,
        "expected_mz": expected,
        "mass_window_da": window,
        "matched_mz": best[0] if best else None,
        "matched_relative_intensity": best[1] if best else None,
    }


def _absence_eligible(
    rule: Mapping[str, Any], quality: Mapping[str, Any] | None
) -> tuple[bool, list[str]]:
    negative = rule["negative_evidence"]
    if negative.get("absence_interpretable") is not True:
        return False, ["rule_does_not_authorize_absence_evidence"]
    quality = quality or {}
    missing = [
        field
        for field in negative.get("required_quality_fields", [])
        if quality.get(field) is not True
    ]
    return not missing, missing


def _not_applicable(rule: Mapping[str, Any], reason: str) -> dict[str, Any]:
    return {
        "rule_id": rule["rule_id"],
        "status": "not_applicable",
        "reason": reason,
        "independent_fragmentation_group": rule["independent_fragmentation_group"],
        "rule_standard": _rule_standard(rule),
        "score_nats": 0.0,
        "score_bits": 0.0,
        "score_decibans": 0.0,
        "log_likelihood_ratio": 0.0,
    }


def evaluate_feature_fragment_evidence(
    *,
    peaks: Sequence[Sequence[float]],
    precursor_mz: float | None,
    context: Mapping[str, Any],
    catalog: Mapping[str, Any],
) -> dict[str, Any]:
    """Evaluate auditable class-level MS/MS evidence in shadow mode.

    ``context`` must provide an independently assigned chemical class.  A
    class inferred by this module cannot be fed back as its own premise.
    """
    normalized = _normalized_peaks(peaks)
    base_output = {
        "module": "v15_feature_fragments",
        "status": "available",
        "evidence_stage": "msms_rule_support",
        "runtime_influence": False,
        "identity_rank_influence": False,
        "score_contribution": 0.0,
        "exact_identity_confirmed": False,
        "standard_confirmed": False,
        "isomer_discriminating": False,
        "confirmation_status": "not_confirmed_by_feature_fragments",
        "feature_fragment_score_standard": feature_fragment_score_standard(),
        "informative_independent_group_count": 0,
        "present_independent_fragment_count": 0,
        "eligible_absence_independent_group_count": 0,
        "cumulative_evidence_by_independent_count": [],
        "combined_score_nats": 0.0,
        "combined_score_bits": 0.0,
        "combined_score_decibans": 0.0,
        "query_rule_context": {
            "chemical_class": context.get("chemical_class"),
            "polarity": context.get("polarity"),
            "adduct": context.get("adduct"),
            "charge": context.get("charge"),
            "instrument_domain_id": context.get("instrument_domain_id"),
            "activation": context.get("activation"),
            "collision_energy": dict(context.get("collision_energy") or {}),
        },
    }
    if not normalized:
        return {
            **base_output,
            "status": "unavailable",
            "reason": "no_valid_intensity_bearing_product_ions",
            "rule_evidence": [],
            "combined_log_likelihood_ratio": 0.0,
            "likelihood_ratio": 1.0,
        }
    if context.get("classification_reliable") is not True:
        return {
            **base_output,
            "status": "not_applicable",
            "reason": "independent_chemical_class_not_established",
            "rule_evidence": [],
            "combined_log_likelihood_ratio": 0.0,
            "likelihood_ratio": 1.0,
        }
    if context.get("classification_source") == "feature_fragments":
        return {
            **base_output,
            "status": "not_applicable",
            "reason": "circular_classification_source",
            "rule_evidence": [],
            "combined_log_likelihood_ratio": 0.0,
            "likelihood_ratio": 1.0,
        }

    chemical_class = str(context.get("chemical_class") or "")
    polarity = str(context.get("polarity") or "").lower()
    adduct = str(context.get("adduct") or "")
    domain_id = str(context.get("instrument_domain_id") or "")
    activation = str(context.get("activation") or "")
    try:
        charge = int(context.get("charge"))
    except (TypeError, ValueError):
        charge = 0

    evaluated: list[dict[str, Any]] = []
    for rule in catalog["rules"]:
        if chemical_class not in rule["chemical_classes"]:
            evaluated.append(_not_applicable(rule, "chemical_class_mismatch"))
            continue
        if polarity != rule["polarity"]:
            evaluated.append(_not_applicable(rule, "polarity_mismatch"))
            continue
        if adduct not in rule["adducts"]:
            evaluated.append(_not_applicable(rule, "adduct_mismatch"))
            continue
        if charge not in rule["charges"]:
            evaluated.append(_not_applicable(rule, "charge_mismatch"))
            continue
        if domain_id not in rule["instrument_domain_ids"]:
            evaluated.append(_not_applicable(rule, "instrument_domain_mismatch"))
            continue
        if activation != rule["activation"]:
            evaluated.append(_not_applicable(rule, "activation_mismatch"))
            continue
        ce_ok, ce_reason = _ce_applicability(
            rule["collision_energy"], context.get("collision_energy")
        )
        if not ce_ok:
            evaluated.append(_not_applicable(rule, ce_reason))
            continue
        observed = _observation(rule, normalized, precursor_mz)
        if not observed.get("observable"):
            evaluated.append(_not_applicable(rule, observed["reason"]))
            continue

        present = bool(observed["present"])
        absence_ok, absence_missing = _absence_eligible(rule, context.get("quality"))
        if not present and not absence_ok:
            evaluated.append(
                {
                    "rule_id": rule["rule_id"],
                    "status": "uninformative_absence",
                    "reason": "absence_quality_gate_failed",
                    "missing_quality_fields": absence_missing,
                    "independent_fragmentation_group": rule[
                        "independent_fragmentation_group"
                    ],
                    "observation": observed,
                    "rule_standard": _rule_standard(rule),
                    "score_nats": 0.0,
                    "score_bits": 0.0,
                    "score_decibans": 0.0,
                    "log_likelihood_ratio": 0.0,
                }
            )
            continue
        occurrence = rule["occurrence_statistics"]
        background = rule["background_statistics"]
        lr = conditional_log_likelihood_ratio(
            target_present=int(occurrence["present"]),
            target_total=int(occurrence["total"]),
            background_present=int(background["present"]),
            background_total=int(background["total"]),
            observed_present=present,
        )
        evaluated.append(
            {
                "rule_id": rule["rule_id"],
                "status": "present" if present else "eligible_absence",
                "independent_fragmentation_group": rule[
                    "independent_fragmentation_group"
                ],
                "observation": observed,
                "frequency_basis": {
                    "target": occurrence,
                    "background": background,
                },
                **lr,
                **_score_units(float(lr["log_likelihood_ratio"])),
                "rule_standard": _rule_standard(rule),
                "collision_energy_applicability": {
                    "status": "applicable",
                    "query": dict(context.get("collision_energy") or {}),
                    "rule": dict(rule.get("collision_energy") or {}),
                    "cross_unit_conversion_applied": False,
                },
                "positive_evidence_limitations": rule["positive_evidence_limitations"],
                "negative_evidence_limitations": rule["negative_evidence"][
                    "limitations"
                ],
                "evidence_level": rule["evidence_level"],
            }
        )

    informative = [
        item for item in evaluated if item["status"] in {"present", "eligible_absence"}
    ]
    representatives: dict[str, dict[str, Any]] = {}
    for item in informative:
        group = item["independent_fragmentation_group"]
        previous = representatives.get(group)
        if previous is None or abs(item["log_likelihood_ratio"]) > abs(
            previous["log_likelihood_ratio"]
        ):
            representatives[group] = item
    representative_items = sorted(
        representatives.values(),
        key=lambda item: (
            str(item["independent_fragmentation_group"]),
            str(item["rule_id"]),
        ),
    )
    combined = sum(item["log_likelihood_ratio"] for item in representative_items)
    cumulative_steps = []
    cumulative = 0.0
    for count, item in enumerate(representative_items, 1):
        contribution = float(item["log_likelihood_ratio"])
        cumulative += contribution
        cumulative_steps.append({
            "independent_evidence_count": count,
            "independent_fragmentation_group": item[
                "independent_fragmentation_group"
            ],
            "selected_rule_id": item["rule_id"],
            "observation_status": item["status"],
            "increment_score_nats": contribution,
            "increment_score_bits": contribution / math.log(2.0),
            "increment_score_decibans": 10.0 * contribution / math.log(10.0),
            "cumulative_score_nats": cumulative,
            "cumulative_score_bits": cumulative / math.log(2.0),
            "cumulative_score_decibans": 10.0 * cumulative / math.log(10.0),
            "cumulative_likelihood_ratio": math.exp(
                max(min(cumulative, 700.0), -700.0)
            ),
        })
    present_representatives = [
        item for item in representative_items if item["status"] == "present"
    ]
    absent_representatives = [
        item for item in representative_items if item["status"] == "eligible_absence"
    ]
    return {
        **base_output,
        "reason": "conditional_lr_shadow_evidence",
        "rule_evidence": evaluated,
        "independent_group_representatives": sorted(representatives),
        "independent_group_representative_evidence": representative_items,
        "informative_independent_group_count": len(representative_items),
        "present_independent_fragment_count": len(present_representatives),
        "eligible_absence_independent_group_count": len(absent_representatives),
        "cumulative_evidence_by_independent_count": cumulative_steps,
        "combined_log_likelihood_ratio": combined,
        **{
            f"combined_{key}": value
            for key, value in _score_units(combined).items()
        },
        "likelihood_ratio": math.exp(max(min(combined, 700.0), -700.0)),
        "interpretation": (
            "Class-conditional MS/MS support only; exact identity, formula, "
            "regioisomer, stereoisomer, and standard confirmation remain unchanged."
        ),
    }
