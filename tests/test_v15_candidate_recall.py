from __future__ import annotations

import copy

import pytest

from matching.v15_candidate_recall import (
    ALGORITHM_ID,
    attach_v15_candidate_recall,
    load_v15_candidate_recall_policy,
    semantic_sha256,
)


def _base_payload() -> dict:
    candidate = {
        "identity_key": "std2d:AAAAAAAAAAAAAA",
        "entity_score": 0.81,
        "entity_rank": 1,
        "candidate": {"formula": "C10H12O2"},
    }
    return {
        "metadata": {
            "task_id": "recall-test",
            "source_semantic_sha256": "a" * 64,
        },
        "invariants": {},
        "entities": [
            {
                "entity_id": "entity-1",
                "decision": {
                    "status": "selected",
                    "scope": "2d_connectivity_candidate",
                },
                "winner": copy.deepcopy(candidate),
                "candidates": [copy.deepcopy(candidate)],
            }
        ],
    }


def _recall_bundle(base: dict) -> dict:
    policy = load_v15_candidate_recall_policy()
    return {
        "metadata": {
            "base_source_semantic_sha256": "a" * 64,
            "base_artifact_sha256": semantic_sha256(base),
            "policy_sha256": semantic_sha256(policy),
            "library_snapshot_id": "library-test",
        },
        "entities": [
            {
                "entity_id": "entity-1",
                "target_total_connectivities": 40,
                "expanded_candidate_count": 1,
                "expanded_candidates": [
                    {
                        "identity_key": "std2d:BBBBBBBBBBBBBB",
                        "connectivity_identity": "BBBBBBBBBBBBBB",
                        "retrieval_rank": 1,
                        "retrieval_rrf_utility": 0.031,
                        "distinct_condition_support": 2,
                        "spectral_view_support": ["library_match", "sensitive"],
                    }
                ],
            }
        ],
    }


def test_v15_candidate_recall_policy_is_locked_to_selected_rrf40_variant():
    policy = load_v15_candidate_recall_policy()

    assert policy["algorithm_id"] == ALGORITHM_ID
    assert policy["candidate_pool"]["target_total_connectivities"] == 40
    assert policy["candidate_pool"]["maximum_additions_per_entity"] == 40
    assert [view["name"] for view in policy["spectral_views"]] == [
        "library_match",
        "sensitive",
    ]
    assert policy["v15_integration"]["neutral_loss_view_selected"] is False


def test_candidate_recall_attachment_cannot_change_primary_rank_score_or_decision():
    base = _base_payload()
    final = attach_v15_candidate_recall(
        base,
        _recall_bundle(base),
        base_payload=base,
    )
    entity = final["entities"][0]
    support = entity["candidate_recall_support"]

    assert entity["winner"] == base["entities"][0]["winner"]
    assert entity["candidates"] == base["entities"][0]["candidates"]
    assert entity["decision"] == base["entities"][0]["decision"]
    assert support["expanded_candidates"][0]["can_replace_winner"] is False
    assert support["expanded_candidates"][0][
        "can_receive_2d_connectivity_decision"
    ] is False
    assert final["metadata"]["candidate_recall_target_total_connectivities"] == 40
    assert final["invariants"]["candidate_recall_support_affects_primary_ranking"] is False


def test_candidate_recall_cannot_be_rebound_to_changed_base_ranking():
    base = _base_payload()
    bundle = _recall_bundle(base)
    changed = copy.deepcopy(base)
    changed["entities"][0]["winner"]["entity_score"] = 0.79

    with pytest.raises(ValueError, match="another base ranking"):
        attach_v15_candidate_recall(
            changed,
            bundle,
            base_payload=changed,
        )
