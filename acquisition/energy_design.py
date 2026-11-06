"""Instrument-domain-aware collision-energy scouting for EviAcquire-MS.

The planner deliberately does not infer a universally sufficient collision
energy from the energy used in the original run.  A single observed energy
cannot reveal an energy-response curve, and values in eV, NCE, volts and
percent are not interchangeable across platforms.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from utils.collision_energy import canonical_collision_energy_unit


def _number(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) and parsed >= 0 else None


def _numbers(values: Any) -> tuple[float, ...]:
    if not isinstance(values, (list, tuple, set)):
        return ()
    parsed = {_number(value) for value in values}
    return tuple(sorted(value for value in parsed if value is not None and value > 0))


def classify_instrument_domain(
    row: Mapping[str, Any], context: Mapping[str, Any] | None
) -> tuple[str, str]:
    context = context or {}
    vendor = context.get("vendor_conversion") or {}
    source = vendor.get("source") or {}
    acquisition = source.get("acquisition") or {}
    fragments = [
        row.get("instrument_type"),
        row.get("instrument"),
        context.get("source_raw_filename"),
        vendor.get("format"),
        source.get("format"),
        acquisition.get("instrument"),
        acquisition.get("acquisition_mode"),
    ]
    text = " ".join(str(value) for value in fragments if value).strip()
    normalized = text.lower().replace("–", "-")
    is_agilent = "agilent" in normalized or "masshunter" in normalized
    is_qtof = any(token in normalized for token in ("q-tof", "qtof", "quadrupole time"))
    # Some mzML records expose only ``tof`` plus an Agilent model.  In an MS2
    # acquisition this still identifies the Agilent Q-TOF domain conservatively.
    if not is_qtof and is_agilent and "tof" in normalized:
        is_qtof = True
    if is_agilent and is_qtof:
        return "agilent_qtof", text or "Agilent Q-TOF"
    if text:
        return "unqualified_instrument_domain", text
    return "unknown_instrument_domain", "not reported"


def _observed_energies(row: Mapping[str, Any]) -> tuple[float, ...]:
    observed = _numbers(row.get("precursor_feature_collision_energies"))
    current = _number(row.get("collision_energy"))
    if current is not None and current > 0:
        observed = tuple(sorted(set(observed) | {current}))
    return observed


def _stage(
    stage_id: str,
    energies: tuple[float, ...],
    acquire_when: str,
    objective: str,
) -> dict[str, Any]:
    return {
        "stage_id": stage_id,
        "collision_energies": energies,
        "acquire_when": acquire_when,
        "objective": objective,
    }


def build_collision_energy_recommendation(
    row: Mapping[str, Any],
    policy: Mapping[str, Any],
    *,
    context: Mapping[str, Any] | None = None,
    ignore_observed_energy: bool = False,
) -> dict[str, Any]:
    """Return a staged CE recommendation without cross-instrument conversion."""

    configured = policy.get("collision_energy_design") or {}
    raw_unit = row.get("collision_energy_unit")
    unit = canonical_collision_energy_unit(raw_unit)
    observed = () if ignore_observed_energy else _observed_energies(row)
    if ignore_observed_energy:
        current = None
    else:
        current = (
            observed[-1]
            if len(observed) == 1
            else _number(row.get("collision_energy"))
        )
    domain, domain_detail = classify_instrument_domain(row, context)

    if len(observed) == 1:
        evidence_status = "single_energy_only_no_response_curve"
    elif len(observed) > 1:
        evidence_status = "multiple_energies_observed_response_curve_not_yet_evaluated"
    else:
        evidence_status = "no_observed_energy_response_curve"

    profile = configured.get("agilent_qtof_ev") or {}
    if domain == "agilent_qtof" and unit == "eV" and profile:
        initial = _numbers(profile.get("initial_scout_energies_eV"))
        low = _numbers(profile.get("conditional_low_energy_eV"))
        midpoint = _numbers(profile.get("conditional_midpoint_energy_eV"))
        high = _numbers(profile.get("conditional_high_energy_eV"))
        ceiling = _number(profile.get("automatic_ceiling_eV"))
        stages = (
            _stage(
                "initial_scout",
                initial,
                "first energy-resolved targeted reacquisition",
                "Bracket low-, intermediate- and high-energy fragmentation without assuming any one energy is sufficient.",
            ),
            _stage(
                "labile_low_energy_refinement",
                low,
                "only if precursor survival is already below 0.50 at the lowest scout energy or low-energy diagnostic onset is unresolved",
                "Avoid overfragmentation and retain precursor-to-primary-fragment information.",
            ),
            _stage(
                "transition_midpoint_refinement",
                midpoint,
                "only if precursor survival crosses 0.50 between adjacent scout energies or fragment novelty changes inside that bracket",
                "Refine the energy-response transition instead of adding every intermediate energy by default.",
            ),
            _stage(
                "resistant_high_energy_extension",
                high,
                "only if precursor survival remains above 0.50 at 40 eV or informative product-ion count remains below the preregistered minimum",
                "Test resistant precursors above 30 eV while keeping high-energy acquisition conditional.",
            ),
        )
        design = {
            "profile_id": str(profile.get("profile_id") or "agilent_qtof_ev_sequential"),
            "strategy": "sequential_energy_resolved_scout_and_refine",
            "instrument_domain": domain_detail,
            "collision_energy_unit": "eV",
            "existing_collision_energy": current,
            "observed_collision_energies": observed,
            "evidence_status": evidence_status,
            "stages": stages,
            "response_metrics": tuple(profile.get("response_metrics") or ()),
            "automatic_ceiling": ceiling,
            "numeric_schedule_available": bool(initial),
            "universal_complete_fragmentation_assumed": False,
            "cross_instrument_conversion_applied": False,
            "method_owner_review_required": True,
        }
        return {
            "initial_energies": initial,
            "collision_energy_unit": "eV",
            "design": design,
            "rationale": (
                "10/20/40 eV is an initial Agilent Q-TOF scouting bracket, not a universal 'complete fragmentation' claim.",
                "5, 30 or 60 eV is acquired only after the measured energy-response curve indicates a labile, transition or resistant precursor.",
                "The highest automatic recommendation is 60 eV; any further extension requires explicit method-owner review.",
            ),
        }

    unsupported = configured.get("unsupported_or_unknown_domain") or {}
    stages = (
        _stage(
            "method_owner_defined_scout",
            (),
            "after instrument family, activation mode and collision-energy unit are verified",
            "Define a platform-native energy bracket without converting eV, NCE, volts or percent.",
        ),
    )
    design = {
        "profile_id": str(
            unsupported.get("profile_id") or "manual_instrument_native_energy_design"
        ),
        "strategy": "method_owner_defined_energy_resolved_scout",
        "instrument_domain": domain_detail,
        "collision_energy_unit": unit or "instrument_native_unspecified",
        "existing_collision_energy": current,
        "observed_collision_energies": observed,
        "evidence_status": evidence_status,
        "stages": stages,
        "response_metrics": tuple(unsupported.get("response_metrics") or ()),
        "automatic_ceiling": None,
        "numeric_schedule_available": False,
        "universal_complete_fragmentation_assumed": False,
        "cross_instrument_conversion_applied": False,
        "method_owner_review_required": True,
    }
    return {
        "initial_energies": (),
        "collision_energy_unit": unit or "instrument_native_unspecified",
        "design": design,
        "rationale": (
            "No numeric collision-energy schedule is emitted because the instrument/activation domain or unit is not qualified.",
            "eV, NCE, collision voltage and percent are preserved as distinct domains and are never numerically converted here.",
        ),
    }
