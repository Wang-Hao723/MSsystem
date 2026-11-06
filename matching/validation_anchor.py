"""Truth-anchor integrity checks for validation and reference-standard QC.

The functions in this module are deliberately outside candidate retrieval and
ranking.  They compare a declared validation target with molecular-formula,
polarity, adduct and already observed precursor events.  A suggested alternate
anchor is an audit lead only; truth records are never rewritten automatically.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from matching.library_quality import formula_exact_mass


PROTON_MASS = 1.007276466621
SODIUM_ION_MASS = 22.989218
POTASSIUM_ION_MASS = 38.963158
AMMONIUM_ION_MASS = 18.033823
CHLORIDE_SHIFT = 34.968853
FORMATE_SHIFT = 44.998201
ACETATE_SHIFT = 59.013851
WATER_MASS = 18.010564684

DEFAULT_POSITIVE_ADDUCTS = (
    "[M]+",
    "[M+H]+",
    "[M+Na]+",
    "[M+NH4]+",
    "[M+K]+",
    "[M-H2O+H]+",
    "[M+2H]2+",
)
DEFAULT_NEGATIVE_ADDUCTS = (
    "[M]-",
    "[M-H]-",
    "[M+Cl]-",
    "[M+FA-H]-",
    "[M+Hac-H]-",
    "[M-2H]2-",
)


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def normalize_polarity(value: Any) -> str | None:
    text = str(value or "").strip().upper()
    if text in {"POS", "POSITIVE", "+", "1", "+1"}:
        return "POS"
    if text in {"NEG", "NEGATIVE", "-", "-1"}:
        return "NEG"
    return None


def precursor_tolerance(mz: float, *, da: float = 0.02, ppm: float = 10.0) -> float:
    return max(float(da), abs(float(mz)) * float(ppm) / 1e6)


def _adduct_mz(neutral_mass: float, label: str) -> tuple[float, int] | None:
    formulas = {
        "[M]+": (neutral_mass, 1),
        "[M+H]+": (neutral_mass + PROTON_MASS, 1),
        "[M+NA]+": (neutral_mass + SODIUM_ION_MASS, 1),
        "[M+NH4]+": (neutral_mass + AMMONIUM_ION_MASS, 1),
        "[M+K]+": (neutral_mass + POTASSIUM_ION_MASS, 1),
        "[M-H2O+H]+": (neutral_mass - WATER_MASS + PROTON_MASS, 1),
        "[M+H-H2O]+": (neutral_mass - WATER_MASS + PROTON_MASS, 1),
        "[M+2H]2+": ((neutral_mass + 2 * PROTON_MASS) / 2, 2),
        "[M]-": (neutral_mass, 1),
        "[M-H]-": (neutral_mass - PROTON_MASS, 1),
        "[M+CL]-": (neutral_mass + CHLORIDE_SHIFT, 1),
        "[M+FA-H]-": (neutral_mass + FORMATE_SHIFT, 1),
        "[M+HCOO]-": (neutral_mass + FORMATE_SHIFT, 1),
        "[M+HAC-H]-": (neutral_mass + ACETATE_SHIFT, 1),
        "[M+CH3COO]-": (neutral_mass + ACETATE_SHIFT, 1),
        "[M-2H]2-": ((neutral_mass - 2 * PROTON_MASS) / 2, 2),
    }
    result = formulas.get(str(label or "").strip().upper().replace(" ", ""))
    if result is None or result[0] <= 0:
        return None
    return result


def common_adduct_hypotheses(
    neutral_mass: float,
    *,
    positive_adducts: Sequence[str] = DEFAULT_POSITIVE_ADDUCTS,
    negative_adducts: Sequence[str] = DEFAULT_NEGATIVE_ADDUCTS,
) -> list[dict[str, Any]]:
    """Return deterministic, mode-labelled theoretical precursor hypotheses."""

    output: list[dict[str, Any]] = []
    for polarity, labels in (("POS", positive_adducts), ("NEG", negative_adducts)):
        for label in labels:
            calculated = _adduct_mz(float(neutral_mass), label)
            if calculated is None:
                continue
            mz, charge = calculated
            output.append(
                {
                    "adduct": str(label),
                    "polarity": polarity,
                    "charge_state": charge,
                    "theoretical_mz": round(mz, 8),
                }
            )
    return output


def _first_number(record: Mapping[str, Any], fields: Sequence[str]) -> float | None:
    for field in fields:
        value = _finite(record.get(field))
        if value is not None:
            return value
    return None


def _event_projection(event: Mapping[str, Any], fallback_polarity: str | None) -> dict[str, Any] | None:
    mz = _first_number(
        event,
        (
            "observed_mz",
            "representative_search_precursor_mz",
            "representative_precursor_mz",
            "search_precursor_mz",
            "precursor_mz",
            "mz",
        ),
    )
    if mz is None or mz <= 0:
        return None
    return {
        "event_id": str(
            event.get("event_id")
            or event.get("entity_id")
            or event.get("spectrum_id")
            or ""
        ),
        "observed_mz": mz,
        "polarity": normalize_polarity(
            event.get("polarity") or event.get("ion_mode") or fallback_polarity
        ),
        "charge_state": (
            int(round(abs(charge)))
            if (charge := _finite(event.get("charge_state") or event.get("precursor_charge")))
            not in (None, 0)
            else None
        ),
        "source": str(event.get("source") or event.get("event_source") or "observed_event"),
    }


def _matches(
    observed_mz: float,
    theoretical_mz: float,
    *,
    tolerance_da: float,
    tolerance_ppm: float,
) -> tuple[bool, float, float]:
    delta = float(observed_mz) - float(theoretical_mz)
    ppm = delta / float(theoretical_mz) * 1e6
    return (
        abs(delta)
        <= precursor_tolerance(
            theoretical_mz,
            da=tolerance_da,
            ppm=tolerance_ppm,
        ),
        delta,
        ppm,
    )


def audit_validation_anchor(
    truth: Mapping[str, Any],
    observed_events: Iterable[Mapping[str, Any]],
    *,
    precursor_tolerance_da: float = 0.02,
    precursor_tolerance_ppm: float = 10.0,
    formula_mass_tolerance_da: float = 0.005,
    formula_mass_tolerance_ppm: float = 10.0,
    positive_adducts: Sequence[str] = DEFAULT_POSITIVE_ADDUCTS,
    negative_adducts: Sequence[str] = DEFAULT_NEGATIVE_ADDUCTS,
) -> dict[str, Any]:
    """Audit one validation anchor without changing truth or predictions.

    Flexible input keys are accepted so the same routine can cover CSV truth,
    reference-standard libraries and saved entity artifacts.  The return value
    always includes a fail-closed ``metric_eligible`` flag.
    """

    declared_formula = str(
        truth.get("formula") or truth.get("molecular_formula") or ""
    ).strip()
    intrinsic_formula_charge = None
    formula = declared_formula
    if formula.endswith("+") or formula.endswith("-"):
        intrinsic_formula_charge = 1 if formula.endswith("+") else -1
        formula = formula[:-1].strip()
    formula_mass = formula_exact_mass(formula) if formula else None
    supplied_neutral = _first_number(
        truth,
        ("neutral_mass", "monoisotopic_mass", "exact_mass"),
    )
    neutral_mass = formula_mass if formula_mass is not None else supplied_neutral
    polarity = normalize_polarity(
        truth.get("polarity") or truth.get("mode") or truth.get("ion_mode")
    )
    declared_adduct = str(
        truth.get("declared_adduct")
        or truth.get("selected_adduct")
        or truth.get("adduct")
        or ""
    ).strip()
    declared_mz = _first_number(
        truth,
        ("declared_mz", "truth_mz", "theoretical_mz", "mz", "precursor_mz"),
    )
    declared_mz_source = "declared"
    if declared_mz is None and neutral_mass is not None and declared_adduct:
        calculated = _adduct_mz(neutral_mass, declared_adduct)
        if calculated is not None:
            declared_mz = calculated[0]
            declared_mz_source = "derived_from_formula_and_declared_adduct"

    projected_events = [
        projected
        for event in observed_events
        if isinstance(event, Mapping)
        and (projected := _event_projection(event, polarity)) is not None
    ]
    hypotheses = (
        common_adduct_hypotheses(
            neutral_mass,
            positive_adducts=positive_adducts,
            negative_adducts=negative_adducts,
        )
        if neutral_mass is not None
        else []
    )

    formula_mass_check: dict[str, Any] | None = None
    formula_mass_conflict = False
    if formula_mass is not None and supplied_neutral is not None:
        delta = supplied_neutral - formula_mass
        tolerance = max(
            formula_mass_tolerance_da,
            formula_mass * formula_mass_tolerance_ppm / 1e6,
        )
        formula_mass_conflict = abs(delta) > tolerance
        formula_mass_check = {
            "formula_exact_mass": round(formula_mass, 8),
            "supplied_neutral_mass": supplied_neutral,
            "delta_da": round(delta, 8),
            "tolerance_da": tolerance,
            "consistent": not formula_mass_conflict,
        }

    declared_matches: list[dict[str, Any]] = []
    if declared_mz is not None:
        for hypothesis in hypotheses:
            matched, delta, ppm = _matches(
                declared_mz,
                float(hypothesis["theoretical_mz"]),
                tolerance_da=precursor_tolerance_da,
                tolerance_ppm=precursor_tolerance_ppm,
            )
            if matched:
                declared_matches.append(
                    {
                        **hypothesis,
                        "declared_delta_da": round(delta, 8),
                        "declared_error_ppm": round(ppm, 4),
                    }
                )
    declared_matches.sort(
        key=lambda item: (abs(item["declared_delta_da"]), item["adduct"])
    )
    same_mode_declared = [
        item for item in declared_matches if polarity and item["polarity"] == polarity
    ]
    opposite_mode_declared = [
        item for item in declared_matches if polarity and item["polarity"] != polarity
    ]
    nominal_resolution = None
    nominal_same_mode_matches: list[dict[str, Any]] = []
    if declared_mz is not None and not same_mode_declared:
        if abs(declared_mz - round(declared_mz)) < 1e-10:
            nominal_resolution = 1.0
        elif abs(declared_mz * 10 - round(declared_mz * 10)) < 1e-9:
            nominal_resolution = 0.1
        if nominal_resolution is not None:
            for hypothesis in hypotheses:
                if hypothesis["polarity"] != polarity:
                    continue
                delta = declared_mz - float(hypothesis["theoretical_mz"])
                if abs(delta) <= nominal_resolution / 2 + precursor_tolerance_da:
                    nominal_same_mode_matches.append(
                        {
                            **hypothesis,
                            "declared_delta_da": round(delta, 8),
                            "inferred_declared_mz_resolution_da": nominal_resolution,
                        }
                    )
            nominal_same_mode_matches.sort(
                key=lambda item: (abs(item["declared_delta_da"]), item["adduct"])
            )

    observed_at_declared: list[dict[str, Any]] = []
    if declared_mz is not None:
        for event in projected_events:
            if event["polarity"] not in {None, polarity}:
                continue
            matched, delta, ppm = _matches(
                event["observed_mz"],
                declared_mz,
                tolerance_da=precursor_tolerance_da,
                tolerance_ppm=precursor_tolerance_ppm,
            )
            if matched:
                observed_at_declared.append(
                    {
                        **event,
                        "delta_to_declared_da": round(delta, 8),
                        "error_to_declared_ppm": round(ppm, 4),
                    }
                )

    formula_compatible_events: list[dict[str, Any]] = []
    for event in projected_events:
        for hypothesis in hypotheses:
            if event["polarity"] not in {None, hypothesis["polarity"]}:
                continue
            charge = event.get("charge_state")
            if charge is not None and charge != hypothesis["charge_state"]:
                continue
            matched, delta, ppm = _matches(
                event["observed_mz"],
                float(hypothesis["theoretical_mz"]),
                tolerance_da=precursor_tolerance_da,
                tolerance_ppm=precursor_tolerance_ppm,
            )
            if matched:
                formula_compatible_events.append(
                    {
                        **event,
                        "adduct": hypothesis["adduct"],
                        "hypothesis_polarity": hypothesis["polarity"],
                        "theoretical_mz": hypothesis["theoretical_mz"],
                        "delta_da": round(delta, 8),
                        "mass_error_ppm": round(ppm, 4),
                    }
                )
    formula_compatible_events.sort(
        key=lambda item: (
            item["hypothesis_polarity"] != polarity,
            abs(item["delta_da"]),
            item["adduct"],
            item["event_id"],
        )
    )
    alternate_events = [
        item
        for item in formula_compatible_events
        if item["hypothesis_polarity"] == polarity
        and not any(
            item["event_id"] == declared["event_id"]
            and abs(item["observed_mz"] - declared["observed_mz"]) < 1e-12
            for declared in observed_at_declared
        )
    ]

    adduct_polarity = None
    if declared_adduct:
        calculated = _adduct_mz(neutral_mass, declared_adduct) if neutral_mass else None
        if calculated is not None:
            adduct_polarity = "POS" if declared_adduct.rstrip().endswith("+") else "NEG"
    adduct_polarity_conflict = bool(
        polarity and adduct_polarity and polarity != adduct_polarity
    )
    declared_adduct_mass_conflict = False
    if declared_adduct and declared_mz is not None and neutral_mass is not None:
        expected = _adduct_mz(neutral_mass, declared_adduct)
        if expected is not None:
            declared_adduct_mass_conflict = not _matches(
                declared_mz,
                expected[0],
                tolerance_da=precursor_tolerance_da,
                tolerance_ppm=precursor_tolerance_ppm,
            )[0]
    intrinsic_charge_adduct_conflict = bool(
        intrinsic_formula_charge is not None
        and declared_adduct
        and declared_adduct not in ({"[M]+"} if intrinsic_formula_charge > 0 else {"[M]-"})
    )

    if not formula or formula_mass is None:
        status = "invalid_formula_or_formula_mass"
        reason = "formula_cannot_generate_neutral_exact_mass"
    elif polarity is None:
        status = "unresolved_polarity"
        reason = "declared_polarity_is_missing_or_invalid"
    elif formula_mass_conflict:
        status = "formula_neutral_mass_conflict"
        reason = "supplied_neutral_mass_conflicts_with_formula_exact_mass"
    elif intrinsic_charge_adduct_conflict:
        status = "intrinsic_formula_charge_adduct_review_required"
        reason = "charged_formula_is_paired_with_an_additional_adduct_declaration"
    elif declared_mz is None:
        status = "missing_declared_anchor"
        reason = "neither_declared_mz_nor_supported_declared_adduct_is_available"
    elif adduct_polarity_conflict:
        status = "declared_adduct_polarity_conflict"
        reason = "declared_adduct_charge_conflicts_with_declared_polarity"
    elif declared_adduct_mass_conflict:
        status = "declared_adduct_mass_conflict"
        reason = "declared_mz_conflicts_with_declared_adduct_and_formula"
    elif not same_mode_declared and opposite_mode_declared:
        status = "opposite_polarity_formula_anchor"
        reason = "declared_mz_matches_only_opposite_polarity_common_adduct"
    elif not same_mode_declared and nominal_same_mode_matches:
        status = "nominal_mz_not_exact_hrms_anchor"
        reason = "declared_mz_precision_is_insufficient_for_exact_precursor_scoring"
    elif not same_mode_declared:
        status = "formula_adduct_anchor_unresolved"
        reason = "declared_mz_matches_no_common_same_polarity_formula_adduct"
    elif observed_at_declared:
        status = "valid_observed_anchor"
        reason = "formula_polarity_adduct_and_observed_precursor_are_consistent"
    elif alternate_events:
        status = "declared_anchor_missing_alternate_observed"
        reason = "no_observed_event_at_declared_mz_but_same_mode_formula_adduct_exists"
    else:
        status = "formula_anchor_not_observed"
        reason = "formula_consistent_declared_anchor_has_no_observed_precursor_event"

    correction_candidates = [
        {
            "observed_mz": item["observed_mz"],
            "adduct": item["adduct"],
            "polarity": item["hypothesis_polarity"],
            "event_id": item["event_id"],
            "mass_error_ppm": item["mass_error_ppm"],
            "action": "manual_truth_and_source_review_required",
        }
        for item in alternate_events[:5]
    ]
    return {
        "schema_version": "mssystem-validation-anchor/1.0.0",
        "status": status,
        "reason": reason,
        "metric_eligible": status == "valid_observed_anchor",
        "formula": declared_formula,
        "neutral_formula_used_for_mass": formula,
        "intrinsic_formula_charge": intrinsic_formula_charge,
        "formula_exact_mass": round(formula_mass, 8) if formula_mass is not None else None,
        "formula_mass_check": formula_mass_check,
        "declared_polarity": polarity,
        "declared_adduct": declared_adduct or None,
        "declared_mz": declared_mz,
        "declared_mz_source": declared_mz_source if declared_mz is not None else None,
        "declared_formula_adduct_matches": declared_matches,
        "nominal_same_mode_adduct_matches": nominal_same_mode_matches,
        "inferred_declared_mz_resolution_da": nominal_resolution,
        "observed_event_count": len(projected_events),
        "observed_at_declared_count": len(observed_at_declared),
        "observed_at_declared": observed_at_declared,
        "formula_compatible_observed_events": formula_compatible_events,
        "correction_candidates": correction_candidates,
        "automatic_truth_rewrite_allowed": False,
        "blind_window_widening_used": False,
        "truth_used_for_retrieval_or_ranking": False,
        "tolerances": {
            "precursor_da": float(precursor_tolerance_da),
            "precursor_ppm": float(precursor_tolerance_ppm),
            "formula_mass_da": float(formula_mass_tolerance_da),
            "formula_mass_ppm": float(formula_mass_tolerance_ppm),
        },
    }


def audit_validation_anchor_batch(
    truth_records: Iterable[Mapping[str, Any]],
    observed_events_by_record: Mapping[str, Iterable[Mapping[str, Any]]],
    **audit_kwargs: Any,
) -> dict[str, Any]:
    """Apply the fail-closed anchor gate to an explicitly keyed validation set."""

    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, truth_record in enumerate(truth_records):
        if not isinstance(truth_record, Mapping):
            raise ValueError(f"truth record {index} is not an object")
        record_id = str(truth_record.get("record_id") or "").strip()
        if not record_id:
            raise ValueError(f"truth record {index} has no record_id")
        if record_id in seen:
            raise ValueError(f"duplicate validation record_id: {record_id}")
        seen.add(record_id)
        truth = truth_record.get("truth")
        if truth is None:
            truth = {
                key: value
                for key, value in truth_record.items()
                if key != "record_id"
            }
        if not isinstance(truth, Mapping):
            raise ValueError(f"truth payload for {record_id} is not an object")
        events = observed_events_by_record.get(record_id, ())
        audit = audit_validation_anchor(
            truth,
            events,
            **audit_kwargs,
        )
        records.append(
            {
                "record_id": record_id,
                "metric_eligible": bool(audit["metric_eligible"]),
                "audit": audit,
            }
        )
    unknown_event_keys = sorted(
        str(key) for key in observed_events_by_record if str(key) not in seen
    )
    statuses = Counter(record["audit"]["status"] for record in records)
    eligible_ids = [
        record["record_id"] for record in records if record["metric_eligible"]
    ]
    correction_queue = [
        {
            "record_id": record["record_id"],
            "status": record["audit"]["status"],
            "reason": record["audit"]["reason"],
            "declared_mz": record["audit"]["declared_mz"],
            "declared_polarity": record["audit"]["declared_polarity"],
            "correction_candidates": record["audit"]["correction_candidates"],
            "action": "manual_truth_and_source_review_required",
        }
        for record in records
        if not record["metric_eligible"]
    ]
    return {
        "schema_version": "mssystem-validation-anchor-batch/1.0.0",
        "summary": {
            "records": len(records),
            "metric_eligible": len(eligible_ids),
            "excluded_from_accuracy_metrics": len(records) - len(eligible_ids),
            "status_counts": dict(sorted(statuses.items())),
            "unknown_observed_event_record_keys": len(unknown_event_keys),
        },
        "metric_eligible_record_ids": eligible_ids,
        "manual_correction_queue": correction_queue,
        "unknown_observed_event_record_ids": unknown_event_keys,
        "records": records,
        "automatic_truth_rewrite_allowed": False,
        "truth_used_for_retrieval_or_ranking": False,
    }


def require_metric_eligible_anchors(batch: Mapping[str, Any]) -> None:
    """Stop an accuracy calculation when unresolved anchors remain included."""

    summary = batch.get("summary") or {}
    excluded = int(summary.get("excluded_from_accuracy_metrics") or 0)
    unknown = int(summary.get("unknown_observed_event_record_keys") or 0)
    if excluded or unknown:
        raise ValueError(
            "validation anchor gate failed: "
            f"excluded={excluded}, unknown_event_keys={unknown}"
        )


__all__ = [
    "DEFAULT_NEGATIVE_ADDUCTS",
    "DEFAULT_POSITIVE_ADDUCTS",
    "audit_validation_anchor",
    "audit_validation_anchor_batch",
    "common_adduct_hypotheses",
    "normalize_polarity",
    "precursor_tolerance",
    "require_metric_eligible_anchors",
]
