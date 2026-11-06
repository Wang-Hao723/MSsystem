from __future__ import annotations

import copy

from matching.v15_decision_assurance import attach_v15_decision_assurance


def _payload() -> dict:
    return {
        "metadata": {},
        "invariants": {},
        "entities": [
            {
                "entity_id": "entity-1",
                "winner": {
                    "identity_key": "ik14:AAAAAAAAAAAAAA",
                    "entity_rank": 1,
                    "entity_score": 0.91,
                },
                "candidates": [
                    {
                        "identity_key": "ik14:AAAAAAAAAAAAAA",
                        "entity_rank": 1,
                        "entity_score": 0.91,
                    }
                ],
                "decision": {
                    "status": "selected",
                    "scope": "2d_connectivity_candidate",
                    "query_spectrum_release_gate": {
                        "status": "failed_audit_only",
                        "deployment_state": "audit_only_failed_locked_transfer_guardrail",
                        "affects_release_scope": False,
                        "failed_gates": ["insufficient_quality_qualified_ce_bins"],
                        "observed_quality_qualified_ce_bins": 1,
                    },
                    "cross_platform_abstention_gate": {
                        "status": "passed",
                        "failed_gates": [],
                        "external_independent_validation_complete": False,
                    },
                },
                "candidate_recall_support": {
                    "status": "evaluated",
                    "target_total_connectivities": 40,
                    "expanded_candidates": [{"identity_key": "std2d:BBBBBBBBBBBBBB"}],
                },
            }
        ],
    }


def test_decision_assurance_fails_closed_without_project_qvalue():
    original = _payload()
    output = attach_v15_decision_assurance(
        original,
        calibration_state={
            "status": "unavailable",
            "production_calibration_id": None,
            "method_id": "mssystem-project-target-decoy-v1",
            "q_value": None,
            "reason": "no_production_project_target_decoy_calibration_registered",
        },
    )
    assurance = output["entities"][0]["decision_assurance"]

    assert assurance["status"] == "released_without_project_qvalue"
    assert assurance["candidate_recall"]["expanded_candidate_count"] == 1
    assert assurance["query_spectrum_quality"]["status"] == (
        "failed_audit_only"
    )
    assert assurance["query_spectrum_quality"]["affects_release_scope"] is False
    assert assurance["project_target_decoy_calibration"]["q_value"] is None
    assert output["metadata"]["entity_qvalues_exposed"] == 0
    assert output["invariants"][
        "uncalibrated_flagship_entities_expose_qvalues"
    ] is False


def test_decision_assurance_does_not_change_rank_score_or_decision():
    original = _payload()
    before = copy.deepcopy(original["entities"][0])
    output = attach_v15_decision_assurance(
        original,
        calibration_state={
            "status": "unavailable",
            "production_calibration_id": None,
            "q_value": None,
        },
    )
    after = output["entities"][0]

    assert after["winner"] == before["winner"]
    assert after["candidates"] == before["candidates"]
    assert after["decision"] == before["decision"]
    assert "decision_assurance" not in original["entities"][0]
