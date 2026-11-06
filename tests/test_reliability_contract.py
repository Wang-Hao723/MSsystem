from __future__ import annotations

import json

import pytest

from matching.v15_flagship_ranker import load_v15_flagship_policy
from pipeline.reliability_contract import build_reliability_contract


def _contract() -> dict:
    return build_reliability_contract(
        target_defaults={"target_screening_enabled": False},
        expected_defaults={"enable_expected_compounds": False},
        flagship_policy=load_v15_flagship_policy(),
    )


def test_contract_records_current_formula_and_rejects_obsolete_weights():
    contract = _contract()
    effective = contract["weighted_score"]["effective_ranking_components"]

    assert effective == {
        "cosine": 0.40,
        "entropy": 0.30,
        "mass": 0.10,
        "fragment": 0.10,
        "intensity": 0.05,
        "expert": 0.0,
    }
    assert [effective[key] for key in ("cosine", "entropy", "mass")] != [
        0.35,
        0.25,
        0.15,
    ]
    assert contract["weighted_score"]["effective_ranking_weight_sum"] == pytest.approx(
        0.95
    )
    assert contract["weighted_score"]["configured_weight_sum"] == pytest.approx(1.0)
    assert contract["weighted_score"]["expert"] == {
        "evidence_computed": True,
        "ranking_influence_enabled": False,
        "effective_input_when_disabled": 0.0,
    }


def test_contract_exposes_gate_post_score_and_shadow_boundaries():
    contract = _contract()

    assert contract["reverse_match"]["role"] == "gate_only"
    assert contract["reverse_match"]["ranking_weight"] == 0.0
    assert contract["fragment_tolerance_da"] == {
        "production_similarity_components": 0.10,
        "multimetric_shadow": 0.10,
        "multimetric_shadow_strict": 0.02,
        "isomer_diagnostic_shadow": 0.02,
    }

    post_score = {stage["name"]: stage for stage in contract["post_score_stages"]}
    assert post_score["classification_consistency"]["enabled"] is False
    assert post_score["classification_consistency"]["mode"] == "soft"
    assert post_score["strict_evidence_gate"]["enabled"] is False

    assert contract["shadow_stages"] == [
        {
            "name": "multimetric_top_k_audit",
            "enabled": True,
            "execution_position": "after_candidate_selection",
            "ranking_influence": False,
            "acceptance_influence": False,
            "confidence_influence": False,
        }
    ]
    assert contract["all_assertions_satisfied"] is True


def test_contract_distinguishes_scan_search_from_production_identity_selector():
    contract = _contract()
    selector = contract["production_identity_selector"]

    assert selector["scan_score_components"] == {
        "cosine": 0.35,
        "entropy": 0.25,
        "reverse_match": 0.15,
        "mass_accuracy": 0.10,
        "fragment_coverage": 0.10,
        "intensity_agreement": 0.05,
    }
    assert selector["scan_score_weight_sum"] == pytest.approx(1.0)
    assert selector["entity_score"] == {
        "spectral_fraction": 0.70,
        "reciprocal_rank_fraction": 0.30,
    }
    conditioned = selector["conditioned_multilibrary"]
    assert conditioned["conditioned_consensus_rule"] == (
        "maximum_per_compatible_provenance_stratum_then_median"
    )
    assert conditioned["best_reference_fraction"] == 0.5
    assert conditioned["consensus_fraction"] == 0.5
    assert conditioned["unconditioned_consensus_affects_ranking"] is False


def test_contract_uses_detached_caller_defaults_and_is_json_serializable():
    target = {"target_screening_enabled": False, "datasets": ("curated",)}
    expected = {"enable_expected_compounds": False, "datasets": []}
    contract = build_reliability_contract(
        target_defaults=target,
        expected_defaults=expected,
    )

    target["target_screening_enabled"] = True
    expected["datasets"].append("later-mutation")

    assert contract["caller_defaults"] == {
        "target": {
            "target_screening_enabled": False,
            "datasets": ["curated"],
        },
        "expected": {
            "enable_expected_compounds": False,
            "datasets": [],
        },
        "source": "caller",
    }
    json.dumps(contract, allow_nan=False)


def test_contract_rejects_non_json_defaults():
    with pytest.raises(TypeError, match="target_defaults must be JSON-serializable"):
        build_reliability_contract(
            target_defaults={"bad": object()},
            expected_defaults={},
        )
