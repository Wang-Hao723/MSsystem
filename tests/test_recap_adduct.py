"""Safety tests for RECAP-Iso adduct belief states."""

from matching.library_quality import formula_exact_mass
from matching.recap_adduct import (
    canonical_adduct_key,
    resolve_isopair_adduct_context,
)
from utils.adduct import ion_mz


FORMULA = "C10H12O2"


def _candidate(inchikey, adduct):
    return {"inchikey": inchikey, "adduct": adduct}


def test_explicit_query_adduct_is_independent_and_high_gate_eligible():
    context = resolve_isopair_adduct_context(
        formula=FORMULA,
        query_precursor_mz=165.0,
        query_adduct="M+H",
        query_adduct_reliable=True,
        query_polarity="POS",
    )

    assert context["status"] == "resolved"
    assert context["reason"] == "explicit_query_adduct"
    assert context["adduct_key"] == canonical_adduct_key("[M+H]+")
    assert context["independent_query_evidence"] is True
    assert context["formula_conditioned"] is False
    assert context["high_precision_gate_eligible"] is True
    assert context["not_global_adduct_probability"] is False


def test_public_library_adduct_spellings_are_canonicalized_by_chemistry():
    assert canonical_adduct_key("[M+HCOOH-H]-") == canonical_adduct_key(
        "[M+FA-H]-"
    )
    assert canonical_adduct_key("[M+CH3COOH-H]-") == canonical_adduct_key(
        "[M+Hac-H]-"
    )
    assert canonical_adduct_key("[M]+*") == canonical_adduct_key("[M]+")


def test_unique_formula_mass_hypothesis_is_ranking_only():
    precursor = ion_mz(formula_exact_mass(FORMULA), "[M+H]+")
    candidates = [
        _candidate("AAAAAAAAAAAAAA-BBBBBBBBBB-C", "[M+H]+"),
        _candidate("DDDDDDDDDDDDDD-EEEEEEEEEE-F", "M+H"),
    ]

    context = resolve_isopair_adduct_context(
        formula=FORMULA,
        query_precursor_mz=precursor,
        query_adduct="",
        query_polarity="POS",
        candidates=candidates,
    )

    assert context["status"] == "resolved"
    assert context["reason"] == "unique_formula_mass_consistent_adduct"
    assert context["resolved_adduct"] == "[M+H]+"
    assert context["candidate_connectivity_support"] == 2
    assert context["formula_conditioned"] is True
    assert context["connectivity_neutral"] is True
    assert context["ranking_eligible"] is True
    assert context["high_precision_gate_eligible"] is False
    assert context["high_promotion_allowed"] is False
    assert context["not_global_adduct_probability"] is True


def test_candidate_labels_cannot_manufacture_a_mass_hypothesis():
    candidates = [
        _candidate("AAAAAAAAAAAAAA-BBBBBBBBBB-C", "[M+H]+"),
        _candidate("DDDDDDDDDDDDDD-EEEEEEEEEE-F", "[M+H]+"),
    ]

    context = resolve_isopair_adduct_context(
        formula=FORMULA,
        query_precursor_mz=101.0,
        query_adduct="",
        query_polarity="POS",
        candidates=candidates,
    )

    assert context["status"] == "unresolved"
    assert context["reason"] == "no_formula_mass_consistent_adduct"
    assert context["ranking_eligible"] is False


def test_explicit_adduct_polarity_conflict_fails_closed():
    context = resolve_isopair_adduct_context(
        formula=FORMULA,
        query_precursor_mz=165.0,
        query_adduct="[M+H]+",
        query_adduct_reliable=True,
        query_polarity="NEG",
    )

    assert context["status"] == "unresolved"
    assert context["reason"] == "explicit_adduct_polarity_conflict"
    assert context["ranking_eligible"] is False
    assert context["high_precision_gate_eligible"] is False


def test_unreliable_explicit_label_is_not_used_as_query_truth():
    precursor = ion_mz(formula_exact_mass(FORMULA), "[M+Na]+")

    context = resolve_isopair_adduct_context(
        formula=FORMULA,
        query_precursor_mz=precursor,
        query_adduct="[M+H]+",
        query_adduct_reliable=False,
        query_polarity="POS",
    )

    assert context["status"] == "resolved"
    assert context["resolved_adduct"] == "[M+Na]+"
    assert context["explicit_adduct_status"] == (
        "explicit_query_adduct_marked_unreliable"
    )
    assert context["formula_conditioned"] is True
    assert context["high_precision_gate_eligible"] is False
