"""Regression tests for expert-score provenance and applicability semantics."""

from matching.expert_scorer import (
    ChemicalExpertScorer,
    expert_evaluate,
    finalize_expert_evidence,
)
from matching.fragment_feature_library import (
    classify_compound_from_fragments,
    get_profile,
)
from matching.fragmentation_rules import normalize_classification
from matching.expert_v3 import score_and_reason as expert_v3_score


FLAVONOID_SMILES = "c1cc2c(cc1O)C(=O)C=C(O2)c3ccc(O)cc3"
GLYCOSIDE_SMILES = "c1ccccc1OC2OC(CO)C(O)C(O)C2O"


def test_fusion_explains_positive_expert_with_zero_base_dimensions():
    details = {
        "expert_total_score": 0.0,
        "neutral_loss_score": 0.0,
        "rule_coverage_score": 0.0,
        "fg_signature_score": 0.0,
        "feature_library_score": 0.0,
        "structure_fragment_score": 0.0,
        "diagnostic_ion_score": 0.0,
        "bond_energy_score": 0.0,
        "isomer_differentiation_score": 0.0,
        "classification_reliable": True,
        "knowledge_applicability": {"neutral_loss_applicable": True},
        "neutral_loss_analysis": {},
        "diagnostic_ion_analysis": {},
        "rule_coverage_analysis": {
            "note": "unsupported_classification",
        },
        "fg_signature_analysis": {},
        "structure_fragment_analysis": {
            "consistency_note": "0/5 fragments structurally plausible",
            "total_fragments_observed": 5,
        },
        "bond_energy_analysis": {
            "note": "labile_bonds_are_prior_only_without_observed_neutral_loss",
        },
        "feature_library_analysis": {
            "note": "no_profile_for_Unknown",
        },
        "isomer_differentiation_analysis": {
            "note": "no_isomer_rules_for_class",
        },
        "pathway_align": {"score": 0.625},
    }

    finalized = finalize_expert_evidence(details, fused_score=0.625)

    assert finalized["fused_expert_score"] == 0.0
    assert "compound_pathway" in finalized["score_fusion"]["suppressed_sources"]
    assert finalized["score_fusion"]["base_dimensions_all_zero"] is True
    assert finalized["score_fusion"]["components"]["base_8d"]["score"] == 0.0
    assert (
        finalized["score_fusion"]["components"]["compound_pathway"]["score"]
        == 0.625
    )
    assert (
        finalized["score_fusion"]["components"]["compound_pathway"]
        ["eligible_for_fusion"] is False
    )


def test_ml_class_does_not_unlock_class_rules_without_structure_confirmation():
    result = expert_evaluate(
        query_peaks=[[171.1386, 100.0], [153.1270, 30.0], [127.1110, 20.0]],
        query_pmz=255.2330,
        candidate_smiles="CCCCCCCCCCCCCCCC(=O)O",
        candidate_classification="lipid",
        query_polarity="NEG",
        candidate_adduct="[M-H]-",
        candidate_classification_source="ml",
        candidate_classification_confidence=0.99,
    )

    provenance = result["classification_provenance"]
    assert provenance["effective_classification"] == "lipid"
    assert provenance["reliable"] is False
    assert provenance["reliability_reason"] == "untrusted_source:ml"
    assert result["analysis_class"] == "Unknown"
    assert (
        result["dimension_assessments"]["rule_coverage"]["status"]
        == "not_applicable"
    )


def test_v3_requires_explicit_independent_class_confirmation():
    result = expert_v3_score(
        [[153.0188, 100.0], [137.0239, 60.0], [121.0290, 30.0]],
        271.0600,
        "flavonoid",
        "[M+H]+",
        query_polarity="POS",
    )

    assert result["v3_score"] == 0.0
    assert result["classification_reliable"] is False
    assert result["reason"] == "classification_unconfirmed"


def test_ml_class_is_confirmed_when_independent_smiles_rule_agrees():
    result = expert_evaluate(
        query_peaks=[[153.0188, 100.0], [137.0239, 60.0], [121.0290, 30.0]],
        query_pmz=271.0600,
        candidate_smiles=FLAVONOID_SMILES,
        candidate_classification="flavonoid",
        query_polarity="POS",
        candidate_adduct="[M+H]+",
        candidate_classification_source="ml",
        candidate_classification_confidence=0.99,
    )

    provenance = result["classification_provenance"]
    assert provenance["structure_class"] == "flavonoid"
    assert provenance["reliable"] is True
    assert (
        provenance["reliability_reason"]
        == "stored_class_agrees_with_smiles_structure_rule"
    )
    assert result["analysis_class"] == "flavonoid"


def test_candidate_adduct_fallback_is_explicitly_audited():
    result = expert_evaluate(
        query_peaks=[[163.0606, 100.0], [145.0500, 50.0], [127.0395, 20.0]],
        query_pmz=301.1000,
        candidate_smiles=GLYCOSIDE_SMILES,
        candidate_classification="glycoside",
        query_polarity="NEG",
        query_adduct="",
        candidate_adduct="[M-H]-",
        candidate_classification_source="smarts",
        candidate_classification_confidence=1.0,
    )

    applicability = result["knowledge_applicability"]
    assert applicability["adduct"] == "[M-H]-"
    assert applicability["adduct_source"] == "candidate_hypothesis"
    assert "candidate_adduct_used_as_hypothesis" in applicability["warnings"]


def test_generic_glycoside_profile_does_not_vote_in_fragment_only_classifier():
    assert get_profile("glycoside") is not None
    classified = classify_compound_from_fragments(
        [145.9, 129.0, 111.0, 87.0, 101.1, 73.0],
        tolerance=0.5,
    )
    assert classified["class_name"] == "sugar"


def test_xanthone_structure_output_has_a_supported_broad_parent():
    assert normalize_classification("xanthone") == "benzenoid"


def test_insufficient_input_is_unavailable_instead_of_false_zero_evidence():
    result = expert_evaluate(
        query_peaks=[],
        query_pmz=300.0,
        candidate_smiles=FLAVONOID_SMILES,
        candidate_classification="flavonoid",
    )

    assert result["evaluation_status"] == "unavailable"
    assert {
        assessment["status"]
        for assessment in result["dimension_assessments"].values()
    } == {"unavailable"}


def test_unknown_class_dynamic_weights_stay_normalized():
    result = ChemicalExpertScorer()._compute_expert_scores(
        nl_analysis={"explanation_rate": 1.0, "priority_weighted_score": 1.0},
        di_analysis={"detection_rate": 1.0, "matched_ions": []},
        rc_analysis={"classification": "Unknown", "coverage_rate": 1.0},
        fg_analysis={"overall_signature_score": 1.0},
        sf_analysis={"score": 1.0},
        be_analysis={"score": 1.0},
        ff_analysis={"score": 1.0},
        id_analysis={"score": 1.0},
    )

    assert result["total"] == 1.0
