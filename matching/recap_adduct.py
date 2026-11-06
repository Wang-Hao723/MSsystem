"""Auditable adduct belief states for formula-conditional isomer ranking.

RECAP-Iso must not turn a library candidate's adduct label into query truth.
This module therefore separates two evidence states:

* an explicit, polarity-consistent query annotation, which is independent of
  the retrieved connectivity candidates; and
* a formula-conditional mass hypothesis, which is shared by every
  connectivity in the retrieved molecular-formula family.

The second state is sufficient to decide whether a same-formula ranker is being
applied to a compatible adduct stratum.  It is not independent confirmation of
the molecular formula or adduct and must never promote or preserve High by
itself.  All return values make that limitation machine-readable.
"""

from __future__ import annotations

import math
import re
from typing import Any, Iterable, Mapping

from matching.library_quality import formula_exact_mass
from utils.adduct import NEUTRAL_SHIFTS, adduct_polarity, ion_mz


_SPECIFIC_ADDUCTS_BY_POLARITY = {
    "POS": (
        "[M]+",
        "[M+H]+",
        "[M+Na]+",
        "[M+K]+",
        "[M+NH4]+",
        "[M+H-H2O]+",
        "[M+H-NH3]+",
    ),
    "NEG": (
        "[M]-",
        "[M-H]-",
        "[M+Cl]-",
        "[M+FA-H]-",
        "[M+Hac-H]-",
        "[M-H-H2O]-",
    ),
}


def _compact(value: Any) -> str:
    return re.sub(r"\s+", "", str(value or "")).upper()


_ADDUCT_KEY_BY_ALIAS: dict[str, str] = {}
for _polarity, _labels in _SPECIFIC_ADDUCTS_BY_POLARITY.items():
    for _label in _labels:
        _parsed = NEUTRAL_SHIFTS.get(_label)
        if _parsed is None:
            continue
        _identity = (
            round(float(_parsed[0]), 6),
            abs(int(_parsed[1])),
            _polarity,
        )
        for _alias, _alias_parsed in NEUTRAL_SHIFTS.items():
            _alias_identity = (
                round(float(_alias_parsed[0]), 6),
                abs(int(_alias_parsed[1])),
                adduct_polarity(_alias),
            )
            if _alias_identity == _identity:
                _ADDUCT_KEY_BY_ALIAS[_compact(_alias)] = _compact(_label)

# Common public-library spellings whose chemistry is explicit but whose text
# is absent from the legacy neutral-shift table.  These are normalization
# aliases only; they do not add a new mass hypothesis.
_ADDUCT_KEY_BY_ALIAS.update({
    "[M+HCOOH-H]-": _compact("[M+FA-H]-"),
    "[M+CH3COOH-H]-": _compact("[M+Hac-H]-"),
    "[M]+*": _compact("[M]+"),
})


def canonical_adduct_key(value: Any) -> str:
    """Return a specific canonical adduct key without polarity-only guessing."""
    compact = _compact(value)
    if not compact or compact in {
        "POSITIVE", "POS", "+", "NEGATIVE", "NEG", "-",
    }:
        return ""
    return _ADDUCT_KEY_BY_ALIAS.get(compact, compact)


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _query_polarity(value: Any, charge: Any = None) -> str:
    text = str(value or "").strip().upper()
    if text in {"POS", "POSITIVE", "+"}:
        return "POS"
    if text in {"NEG", "NEGATIVE", "-"}:
        return "NEG"
    charge_value = _finite(charge)
    if charge_value is not None and charge_value > 0:
        return "POS"
    if charge_value is not None and charge_value < 0:
        return "NEG"
    return ""


def _base_context(status: str, reason: str) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "status": status,
        "reason": reason,
        "resolved_adduct": "",
        "adduct_key": "",
        "evidence_source": "none",
        "independent_query_evidence": False,
        "formula_conditioned": False,
        "connectivity_neutral": True,
        "ranking_eligible": False,
        "high_precision_gate_eligible": False,
        "high_promotion_allowed": False,
        "not_global_adduct_probability": True,
        "hypotheses": [],
    }


def resolve_isopair_adduct_context(
    *,
    formula: Any,
    query_precursor_mz: Any,
    query_adduct: Any = None,
    query_adduct_reliable: bool | None = None,
    query_polarity: Any = None,
    query_charge: Any = None,
    candidates: Iterable[Mapping[str, Any]] = (),
    tolerance_da: float = 0.02,
    tolerance_ppm: float = 10.0,
) -> dict[str, Any]:
    """Resolve the adduct stratum that IsoPair-MS may audit.

    A formula-conditioned resolution is deliberately local to the already
    retrieved formula family.  It can activate within-family ranking, but it
    cannot validate that family, expand candidate recall, or support High.
    """
    polarity = _query_polarity(query_polarity, query_charge)
    explicit_key = canonical_adduct_key(query_adduct)
    if explicit_key and query_adduct_reliable is not False:
        explicit_polarity = adduct_polarity(str(query_adduct or ""))
        if polarity and explicit_polarity and explicit_polarity != polarity:
            context = _base_context("unresolved", "explicit_adduct_polarity_conflict")
            context.update({
                "query_polarity": polarity,
                "reported_query_adduct": str(query_adduct or ""),
            })
            return context
        return {
            **_base_context("resolved", "explicit_query_adduct"),
            "resolved_adduct": str(query_adduct or "").strip(),
            "adduct_key": explicit_key,
            "evidence_source": "explicit_query_metadata",
            "query_polarity": polarity or explicit_polarity or "",
            "reported_query_adduct": str(query_adduct or ""),
            "independent_query_evidence": True,
            "ranking_eligible": True,
            "high_precision_gate_eligible": True,
            "not_global_adduct_probability": False,
            "hypotheses": [{
                "adduct": str(query_adduct or "").strip(),
                "adduct_key": explicit_key,
                "status": "explicit",
            }],
        }

    if query_adduct_reliable is False and str(query_adduct or "").strip():
        rejected_reason = "explicit_query_adduct_marked_unreliable"
    else:
        rejected_reason = "query_adduct_not_explicit"

    precursor = _finite(query_precursor_mz)
    neutral_mass = formula_exact_mass(formula)
    if precursor is None or precursor <= 0.0:
        return _base_context("unresolved", "invalid_query_precursor")
    if neutral_mass is None:
        return _base_context("unresolved", "formula_exact_mass_unavailable")
    if polarity not in _SPECIFIC_ADDUCTS_BY_POLARITY:
        return _base_context("unresolved", "query_polarity_unavailable")

    tolerance = max(float(tolerance_da), precursor * float(tolerance_ppm) * 1e-6)
    candidate_support: dict[str, set[str]] = {}
    for candidate in candidates:
        candidate_key = canonical_adduct_key(candidate.get("adduct"))
        connectivity = str(candidate.get("inchikey") or "").strip().upper()[:14]
        if candidate_key and len(connectivity) == 14:
            candidate_support.setdefault(candidate_key, set()).add(connectivity)

    hypotheses: list[dict[str, Any]] = []
    for label in _SPECIFIC_ADDUCTS_BY_POLARITY[polarity]:
        predicted = ion_mz(neutral_mass, label)
        if predicted is None:
            continue
        error_da = abs(float(predicted) - precursor)
        error_ppm = error_da / precursor * 1e6
        key = canonical_adduct_key(label)
        hypotheses.append({
            "adduct": label,
            "adduct_key": key,
            "predicted_mz": round(float(predicted), 8),
            "error_da": round(error_da, 8),
            "error_ppm": round(error_ppm, 4),
            "within_tolerance": error_da <= tolerance,
            "candidate_connectivity_support": len(candidate_support.get(key, set())),
        })
    hypotheses.sort(key=lambda item: (item["error_da"], item["adduct_key"]))
    compatible = [item for item in hypotheses if item["within_tolerance"]]
    if not compatible:
        context = _base_context("unresolved", "no_formula_mass_consistent_adduct")
        context.update({
            "query_polarity": polarity,
            "formula_exact_mass": round(neutral_mass, 8),
            "mass_tolerance_da": round(tolerance, 8),
            "explicit_adduct_status": rejected_reason,
            "hypotheses": hypotheses,
        })
        return context
    if len(compatible) != 1:
        context = _base_context("unresolved", "multiple_formula_mass_consistent_adducts")
        context.update({
            "query_polarity": polarity,
            "formula_exact_mass": round(neutral_mass, 8),
            "mass_tolerance_da": round(tolerance, 8),
            "explicit_adduct_status": rejected_reason,
            "hypotheses": hypotheses,
        })
        return context

    selected = compatible[0]
    return {
        **_base_context("resolved", "unique_formula_mass_consistent_adduct"),
        "resolved_adduct": selected["adduct"],
        "adduct_key": selected["adduct_key"],
        "evidence_source": "retrieved_formula_exact_mass_and_query_precursor",
        "query_polarity": polarity,
        "formula_exact_mass": round(neutral_mass, 8),
        "mass_tolerance_da": round(tolerance, 8),
        "mass_error_da": selected["error_da"],
        "mass_error_ppm": selected["error_ppm"],
        "explicit_adduct_status": rejected_reason,
        "formula_conditioned": True,
        "ranking_eligible": True,
        "candidate_connectivity_support": selected[
            "candidate_connectivity_support"
        ],
        "hypotheses": hypotheses,
    }


RECAP_ADDUCT_AVAILABLE = True
