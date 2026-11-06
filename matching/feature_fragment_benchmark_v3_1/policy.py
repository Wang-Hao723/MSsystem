"""V3.1 evidence eligibility rules layered over immutable V3 records."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


POSITIVE_ONLY_MODES = {
    "stepped",
    "CES",
    "ramp",
    "rolling",
    "dynamic",
    "MSE",
    "bbCID",
    "all_ions",
}

CONTINUOUS_OR_VARIABLE_MODES = {
    "ramp",
    "rolling",
    "dynamic",
    "MSE",
    "bbCID",
    "all_ions",
}


def _finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(
        float(value)
    )


def composite_program_weight_status(program: Mapping[str, Any]) -> dict[str, Any]:
    """Return whether a CE program can enter quantitative V3.1 analysis.

    A continuous or variable program additionally needs a pre-frozen,
    same-model program and detection-floor validation.  Numeric endpoint values
    alone do not describe dwell weights for a ramp.
    """
    mode = str(program.get("mode") or "unresolved")
    values: Sequence[Any] = program.get("values") or []
    weights: Sequence[Any] = program.get("weights") or []
    if mode == "fixed":
        complete = bool(
            len(values) == 1
            and len(weights) == 1
            and _finite_number(values[0])
            and _finite_number(weights[0])
            and float(weights[0]) > 0
        )
        return {
            "status": "complete" if complete else "incomplete",
            "quantitative_program_eligible": complete,
            "mode": mode,
            "reason": None if complete else "fixed_value_or_weight_incomplete",
        }
    numeric_weights = bool(
        len(values) >= 2
        and len(weights) == len(values)
        and all(_finite_number(value) and float(value) >= 0 for value in values)
        and all(_finite_number(weight) and float(weight) >= 0 for weight in weights)
        and sum(float(weight) for weight in weights) > 0
    )
    if mode in CONTINUOUS_OR_VARIABLE_MODES:
        calibration = program.get("quantitative_program_calibration") or {}
        calibrated = bool(
            calibration.get("status")
            == "validated_same_model_program_detection_floor"
            and calibration.get("frozen_before_discovery") is True
        )
        complete = numeric_weights and calibrated
        return {
            "status": "complete" if complete else "incomplete",
            "quantitative_program_eligible": complete,
            "mode": mode,
            "numeric_weights_complete": numeric_weights,
            "same_model_program_detection_floor_validated": calibrated,
            "reason": None
            if complete
            else "continuous_or_variable_program_not_prevalidated",
        }
    if mode in {"stepped", "CES"}:
        return {
            "status": "complete" if numeric_weights else "incomplete",
            "quantitative_program_eligible": numeric_weights,
            "mode": mode,
            "numeric_weights_complete": numeric_weights,
            "reason": None if numeric_weights else "composite_numeric_weights_incomplete",
        }
    return {
        "status": "incomplete",
        "quantitative_program_eligible": False,
        "mode": mode,
        "reason": "collision_energy_program_unresolved_or_unsupported",
    }


def _peak_qc(record: Mapping[str, Any], qc_policy: Mapping[str, Any]) -> dict[str, Any]:
    peak_policy = qc_policy["peak_presence"]
    relative_threshold = float(peak_policy["minimum_relative_intensity"])
    peaks = record.get("peaks") or []
    retained = [
        peak
        for peak in peaks
        if len(peak) >= 2
        and _finite_number(peak[0])
        and _finite_number(peak[1])
        and float(peak[1]) >= relative_threshold
    ]
    minimum_gate = peak_policy["minimum_deduplicated_peak_count_gate"]
    gate_active = minimum_gate.get("active") is True
    minimum = minimum_gate.get("minimum")
    gate_passed = bool(
        not gate_active
        or (_finite_number(minimum) and len(retained) >= int(float(minimum)))
    )
    return {
        "relative_intensity_threshold": relative_threshold,
        "mass_tolerance": deepcopy(peak_policy["mass_tolerance"]),
        "observed_peak_count": len(peaks),
        "peak_count_at_or_above_relative_threshold": len(retained),
        "minimum_peak_gate_active": gate_active,
        "minimum_peak_gate_passed": gate_passed,
        "minimum_peak_gate_value": minimum,
        "input_peaks_modified": False,
    }


def _precursor_purity(
    record: Mapping[str, Any], qc_policy: Mapping[str, Any]
) -> dict[str, Any]:
    acquisition = record.get("acquisition") or {}
    reported = acquisition.get("precursor_purity") or {}
    policy = qc_policy["precursor_purity"]
    threshold = policy.get("minimum_fraction")
    fraction = reported.get("fraction")
    verified = bool(
        policy.get("threshold_status") == "calibrated_and_frozen"
        and _finite_number(threshold)
        and reported.get("status") == "verified"
        and _finite_number(fraction)
        and float(fraction) >= float(threshold)
    )
    return {
        "status": "verified_at_or_above_frozen_threshold" if verified else "unknown_or_unqualified",
        "reported_status": reported.get("status") or "not_reported",
        "reported_fraction": fraction,
        "frozen_minimum_fraction": threshold,
        "threshold_status": policy.get("threshold_status"),
        "eligible_for_lr_negative": verified,
    }


def apply_v31_policy(
    record: Mapping[str, Any], *, qc_policy: Mapping[str, Any]
) -> dict[str, Any]:
    """Apply only stricter V3.1 gates; never promote a V3-ineligible record."""
    normalized = deepcopy(record)
    acquisition = normalized.get("acquisition") or {}
    program = acquisition.get("collision_energy") or {}
    if program.get("conversion_applied") is not False:
        raise ValueError("collision_energy_conversion_forbidden_v3_1")
    mode = str(program.get("mode") or "unresolved")
    program_status = composite_program_weight_status(program)
    peak_qc = _peak_qc(normalized, qc_policy)
    purity = _precursor_purity(normalized, qc_policy)
    base = normalized.get("quantitative_eligibility") or {}
    base_positive = base.get("positive_lr_eligible") is True
    base_negative = base.get("negative_lr_eligible") is True
    activation_location = acquisition.get("activation_location")
    positive_only = bool(mode in POSITIVE_ONLY_MODES or activation_location == "in_source")
    shared_or_deconvoluted = bool(
        acquisition.get("shared_fragment_channel") is True
        or acquisition.get("deconvoluted_fragment_channel") is True
    )
    positive = bool(
        base_positive
        and program_status["quantitative_program_eligible"]
        and peak_qc["minimum_peak_gate_passed"]
    )
    negative = bool(
        base_negative
        and positive
        and not positive_only
        and not shared_or_deconvoluted
        and purity["eligible_for_lr_negative"]
    )
    reasons = set(base.get("ineligible_reason_codes") or [])
    if not program_status["quantitative_program_eligible"]:
        reasons.add("v31_collision_energy_program_not_quantitatively_eligible")
    if mode in {"stepped", "CES"} and program_status["status"] != "complete":
        reasons.add("v31_composite_numeric_weights_incomplete")
    if mode in CONTINUOUS_OR_VARIABLE_MODES and program_status["status"] != "complete":
        reasons.add("v31_variable_program_detection_floor_not_validated")
    if not peak_qc["minimum_peak_gate_passed"]:
        reasons.add("v31_minimum_peak_gate_failed")
    if positive_only:
        reasons.add("v31_multi_energy_or_in_source_lr_negative_forbidden")
    if shared_or_deconvoluted:
        reasons.add("v31_shared_or_deconvoluted_peak_lr_negative_forbidden")
    if not purity["eligible_for_lr_negative"]:
        reasons.add("v31_precursor_purity_unavailable_for_lr_negative")
    normalized["schema_version"] = "3.1.0"
    normalized["quantitative_eligibility"] = {
        **base,
        "positive_lr_eligible": positive,
        "negative_lr_eligible": negative,
        "positive_only_composite_or_in_source": bool(positive and positive_only),
        "ineligible_reason_codes": sorted(reasons),
    }
    normalized["v3_1_amendment"] = {
        "amendment_id": "v15-instrument-conditioned-feature-fragment-v3.1-20260827",
        "base_schema_version": record.get("schema_version"),
        "base_positive_lr_eligible": base_positive,
        "base_negative_lr_eligible": base_negative,
        "collision_energy_program": program_status,
        "multi_energy_or_in_source_positive_only": positive_only,
        "spectral_qc": peak_qc,
        "precursor_purity": purity,
        "shared_or_deconvoluted_fragment_channel": shared_or_deconvoluted,
        "nce_to_ev_conversion_applied": False,
        "nce_to_ev_derived_value": None,
        "observed_fragment_evidence_allowed": bool(peak_qc["observed_peak_count"]),
        "quantitative_evidence_grade_ceiling": "G4"
        if positive and negative
        else "G3"
        if positive
        else "G0",
        "observed_evidence_grade_ceiling": "G1"
        if peak_qc["observed_peak_count"]
        else "G0",
    }
    return normalized
