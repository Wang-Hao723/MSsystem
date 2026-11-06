"""High-precision and identity-scope tests for IsoPair-MS calibration."""

from matching.honest_calibration import annotate


IK_A1 = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
IK_A2 = "AAAAAAAAAAAAAA-CCCCCCCCCC-D"
IK_B = "DDDDDDDDDDDDDD-EEEEEEEEEE-F"


def _candidate(inchikey, cosine, *, rank=None):
    return {
        "compound_name": inchikey[:4],
        "inchikey": inchikey,
        "formula": "C10H12O2",
        "smiles": "CCOc1ccccc1",
        "has_structure": True,
        "adduct": "[M+H]+",
        "cosine": cosine,
        "final_score": 0.72,
        "rank_score": 0.72 if rank is None else rank,
        "reverse_match": 0.90,
        "mass_ppm": 3.0,
        "classification": "flavonoid",
        "source_tier": "experimental",
    }


def _isopair_evidence(*, support, changed=False, high_gate=True):
    return {
        "schema_version": 1,
        "status": "applied" if changed else "clear_baseline_winner_no_rerank",
        "evaluated": True,
        "high_precision_gate_applicable": high_gate,
        "high_precision_support": support,
        "high_eligible": support,
        "high_promotion_allowed": False,
        "runtime_influence": changed,
        "ranking_changed": changed,
        "selected_connectivity": IK_A1[:14],
        "identity_scope": "2d_connectivity_only",
        "stereochemistry_claimed": False,
        "adduct_context": {
            "formula_conditioned": not high_gate,
            "high_precision_gate_eligible": high_gate,
        },
    }


def test_supported_existing_high_is_limited_to_2d_connectivity_scope():
    best = _candidate(IK_A1, 0.95)
    alternative = _candidate(IK_B, 0.75, rank=0.50)
    evidence = _isopair_evidence(support=True)
    best.update({
        "isopair_evidence": evidence,
        "isopair_high_precision_support": True,
        "isopair_high_eligible": True,
        "top_candidates": [best.copy(), alternative],
    })

    result = annotate(best)

    assert result["compound_confidence"] == "High"
    assert result["two_dimensional_connectivity_confidence"] == "High"
    assert result["exact_structure_confidence"] == "Medium"
    assert result["identification_scope"] == "exact_2d_connectivity"
    assert result["honesty_tier"] == "isopair_connectivity_high"
    assert result["isopair_high_precision_gate_result"] == "supported"


def test_unsupported_existing_high_is_vetoed_to_medium():
    best = _candidate(IK_A1, 0.95)
    alternative = _candidate(IK_B, 0.75, rank=0.50)
    evidence = _isopair_evidence(support=False)
    best.update({
        "isopair_evidence": evidence,
        "isopair_high_precision_support": False,
        "top_candidates": [best.copy(), alternative],
    })

    result = annotate(best)

    assert result["compound_confidence"] == "Medium"
    assert result["isopair_high_precision_gate_required"] is True
    assert result["isopair_high_precision_gate_result"] == "vetoed_to_medium"


def test_reranked_near_tie_remains_ambiguous_and_cannot_be_high():
    best = _candidate(IK_A1, 0.91, rank=0.70)
    alternative = _candidate(IK_B, 0.90, rank=0.69)
    evidence = _isopair_evidence(support=False, changed=True)
    best.update({
        "isopair_evidence": evidence,
        "isopair_high_precision_support": False,
        "top_candidates": [best.copy(), alternative],
    })

    result = annotate(best)

    assert result["isopair_reranked_2d_connectivity"] is True
    assert result["isomer_ambiguous"] is True
    assert result["compound_confidence"] == "Medium"
    assert result["identification_scope"] == "isomer_group"


def test_same_connectivity_variant_preserves_high_at_2d_scope_only():
    best = _candidate(IK_A1, 0.95, rank=0.72)
    stereo_variant = _candidate(IK_A2, 0.94, rank=0.71)
    other_connectivity = _candidate(IK_B, 0.70, rank=0.50)
    evidence = _isopair_evidence(support=True)
    best.update({
        "isopair_evidence": evidence,
        "isopair_high_precision_support": True,
        "top_candidates": [best.copy(), stereo_variant, other_connectivity],
    })

    result = annotate(best)

    assert result["isomer_ambiguous"] is True
    assert result["isomer_ambiguity_type"] == "same_connectivity_variant_group"
    assert result["compound_confidence"] == "Medium"
    assert result["exact_structure_confidence"] == "Medium"
    assert result["two_dimensional_connectivity_confidence"] == "High"
    assert result["identification_scope"] == "exact_2d_connectivity"
    assert result["honesty_tier"] == "connectivity_high_stereo_unresolved"
    assert result["isomer_group_evidence"][
        "two_dimensional_connectivity_resolved"
    ] is True


def test_formula_conditioned_rerank_is_evaluated_without_high_gate():
    best = _candidate(IK_A1, 0.91, rank=0.70)
    alternative = _candidate(IK_B, 0.90, rank=0.69)
    evidence = _isopair_evidence(
        support=False, changed=True, high_gate=False
    )
    best.update({
        "isopair_evidence": evidence,
        "isopair_high_precision_support": False,
        "top_candidates": [best.copy(), alternative],
    })

    result = annotate(best)

    assert result["isopair_reranked_2d_connectivity"] is True
    assert result["isomer_group_evidence"]["isopair_ms"]["evaluated"] is True
    assert result["isopair_high_precision_gate_required"] is False
    assert result["isopair_high_precision_gate_result"] == (
        "ranking_only_not_applicable"
    )
    assert result["compound_confidence"] == "Medium"


def test_formula_conditioned_clear_winner_does_not_preserve_or_veto_high():
    best = _candidate(IK_A1, 0.95)
    alternative = _candidate(IK_B, 0.75, rank=0.50)
    evidence = _isopair_evidence(
        support=False, changed=False, high_gate=False
    )
    best.update({
        "isopair_evidence": evidence,
        "isopair_high_precision_support": False,
        "top_candidates": [best.copy(), alternative],
    })

    result = annotate(best)

    assert result["compound_confidence"] == "High"
    assert result["honesty_tier"] == "confident"
    assert result["isopair_high_precision_gate_required"] is False
    assert result["isopair_high_precision_gate_result"] == (
        "ranking_only_not_applicable"
    )


def test_formula_conditioned_rerank_is_medium_even_without_serialized_rival():
    best = _candidate(IK_A1, 0.95)
    evidence = _isopair_evidence(
        support=False, changed=True, high_gate=False
    )
    best.update({
        "isopair_evidence": evidence,
        "isopair_high_precision_support": False,
        "top_candidates": [best.copy()],
    })

    result = annotate(best)

    assert result["isomer_ambiguous"] is False
    assert result["isopair_reranked_2d_connectivity"] is True
    assert result["compound_confidence"] == "Medium"
    assert result["isopair_high_precision_gate_required"] is False
