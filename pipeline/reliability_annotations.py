"""Identity-free reliability annotations for feature-level result display.

These helpers add provenance and review status only. They do not alter query
peaks, library scores, candidate order, confidence, or EviGate decisions.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from matching.rule_quality import collision_energy_value_unit


def _finite(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def ce_comparability(
    query: Mapping[str, Any],
    reference: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Describe recorded CE compatibility without changing a score."""

    reference = reference or {}
    query_value, query_unit = collision_energy_value_unit({
        "value": query.get("collision_energy"),
        "unit": query.get("collision_energy_unit"),
    })
    reference_value, reference_unit = collision_energy_value_unit({
        "value": reference.get("lib_collision_energy"),
        "unit": reference.get("lib_collision_energy_unit"),
    })
    query_mode = str(query.get("collision_energy_mode") or "").strip().lower()
    reference_mode = str(
        reference.get("lib_collision_energy_mode") or ""
    ).strip().lower()
    query_quantities = list(query.get("collision_energy_quantities") or [])

    if query_mode == "conflict":
        status = "conflict"
    elif query_value is None:
        status = "query_unknown"
    elif reference_value is None:
        status = "reference_unknown"
    elif not query_unit or not reference_unit or query_unit != reference_unit:
        status = "incompatible_encoding"
    elif query_mode in {"stepped", "ramp", "ramped"}:
        # A single scalar cannot stand in for a stepped or ramped program.
        status = "incompatible_encoding"
    elif reference_mode in {"stepped", "ramp", "ramped"}:
        status = "incompatible_encoding"
    elif abs(float(query_value) - float(reference_value)) <= 1e-9:
        status = "same_recorded_condition"
    else:
        status = "different_recorded_condition"

    return {
        "status": status,
        "query": {
            "value": query_value,
            "unit": query_unit or None,
            "mode": query_mode or "unknown",
            "quantities": query_quantities,
            "source": query.get("acquisition_metadata_source") or None,
        },
        "reference": {
            "value": reference_value,
            "unit": reference_unit or None,
            "mode": reference_mode or "unknown",
            "raw": reference.get("lib_collision_energy_raw"),
            "source_accession": reference.get("lib_source_accession"),
        },
        "score_multiplier": 1.0,
        "ranking_influence": False,
        "governance_influence": False,
        "same_recorded_condition_is_not_physical_equivalence": True,
    }


def fragmentation_annotation(row: Mapping[str, Any]) -> dict[str, Any]:
    """Classify query information content for display and audit only."""

    effective_count = len(row.get("query_peaks") or [])
    if effective_count <= 1:
        status = "insufficient_fragments"
    elif bool(row.get("limited_fragment_list")) or effective_count < 5:
        status = "low_information"
    else:
        status = "information_present"
    return {
        "status": status,
        "effective_fragment_count": effective_count,
        "score_multiplier": 1.0,
        "ranking_influence": False,
        "governance_influence": False,
    }


def background_evidence(row: Mapping[str, Any]) -> dict[str, Any]:
    """Expose B1/B2/B3 evidence without inventing evidence or deleting peaks."""

    blank_used = bool(row.get("blank_used"))
    b3 = "unknown"
    b3_reason = "no_comparable_blank_provided"
    if blank_used:
        role = str(row.get("blank_evidence_role") or "unspecified")
        if role == "unspecified":
            b3_reason = "blank_role_or_method_comparability_unverified"
        elif bool(row.get("blank_feature_matched")):
            # Existing blank matching is retained as evidence, but the new
            # strict three-part exclusion rule is not claimed as validated.
            b3_reason = "blank_signal_observed_strict_rule_unvalidated"
        else:
            b3_reason = "blank_signal_not_observed_detection_power_unverified"

    return {
        "B1_cross_precursor": "unknown",
        "B2_incoherent_chromatography": "unknown",
        "B3_blank_confirmed": b3,
        "B3_reason": b3_reason,
        "background_hint": False,
        "action": "retain",
        "new_background_excluded_peak_count": 0,
        "new_background_soft_downweighted_peak_count": 0,
        "raw_peaks_preserved": True,
        "score_multiplier": 1.0,
        "ranking_influence": False,
        "governance_influence": False,
        "rule_version": "mssystem-v16-background-three-state-1.0.0",
    }
