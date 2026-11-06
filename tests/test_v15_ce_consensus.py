from __future__ import annotations

import copy
import hashlib
import json

import pytest

from matching.v15_ce_consensus import (
    apply_v15_ce_evidence,
    build_v15_ce_evidence_bundle,
)
from matching.v15_flagship_ranker import load_v15_flagship_policy


def _candidate(identity: str, formula: str, score: float, rank: int) -> dict:
    return {
        "identity_key": f"std2d:{identity}",
        "candidate": {
            "formula": formula,
            "compound_name": f"display-{identity}",
            "potential_contaminant": False,
        },
        "entity_score": score,
        "entity_rank": rank,
        "vetoed": False,
        "score_is_probability": False,
    }


def _payload(*, margin: float = 0.04) -> dict:
    leader = _candidate("AAAAAAAAAAAAAA", "C10H12O2", 0.80, 1)
    challenger = _candidate(
        "BBBBBBBBBBBBBB", "C10H12O2", 0.80 - margin, 2
    )
    return {
        "metadata": {
            "task_id": "task-v15",
            "source_row_count": 2,
            "source_semantic_sha256": "a" * 64,
        },
        "policy": load_v15_flagship_policy(),
        "invariants": {},
        "entities": [
            {
                "entity_id": "entity-1",
                "member_result_indices": [0, 1],
                "member_spectrum_ids": ["scan-1", "scan-2"],
                "distinct_collision_energy_bins": 2,
                "collision_energy_bins": [
                    {"value": 20.0, "unit": "eV"},
                    {"value": 30.0, "unit": "eV"},
                ],
                "decision": {
                    "status": "selected_isomer_group",
                    "scope": "molecular_formula_isomer_group",
                },
                "winner": leader,
                "leading_candidate_when_abstained": None,
                "candidates": [leader, challenger],
            }
        ],
    }


class _Retriever:
    library_snapshot_id = "sqlite-test-snapshot"
    policy_sha256 = "b" * 64

    def __init__(self, *, platform_delta: float = 0.02) -> None:
        self.calls = 0
        self.platform_delta = platform_delta

    def retrieve_entity(self, entity_id, scans, *, existing_candidates, detection_method_id):
        self.calls += 1
        comparison = {
            "challenger_identity_key": "std2d:BBBBBBBBBBBBBB",
            "qualified": True,
            "median_score_delta": 0.03,
            "reference_platform_audit": {
                "paired_observation_count": 6,
                "common_provenance_stratum_count": 1,
                "distinct_query_collision_energy_condition_count": 2,
                "spectral_view_count": 2,
                "unique_provenance_evidence_profile_count": 1,
                "median_score_delta": self.platform_delta,
                "observation_win_fraction": 5 / 6,
                "unique_provenance_profile_win_fraction": 1.0,
                "reference_collision_energy_missing_treated_as_match": False,
            },
            "reference_condition_audit": {
                "paired_observation_count": 0,
                "cross_unit_conversion_applied": False,
            },
        }
        return {
            "expanded_candidate_count": 3,
            "scan_audit": [],
            "source_paired_consensus": {
                "status": "qualified_switch",
                "baseline_identity_key": "std2d:AAAAAAAAAAAAAA",
                "selected_identity_key": "std2d:BBBBBBBBBBBBBB",
                "comparisons": [comparison],
                "ground_truth_used_at_inference": False,
            },
        }


def _rows() -> list[dict]:
    return [
        {"spectrum_id": "scan-1", "collision_energy": 20.0},
        {"spectrum_id": "scan-2", "collision_energy": 30.0},
    ]


def test_qualified_platform_anchored_ce_switch_preserves_base_scores():
    base = _payload()
    retriever = _Retriever()
    evidence = build_v15_ce_evidence_bundle(
        base, _rows(), retriever=retriever, detection_method_id="qtof-dda"
    )
    final = apply_v15_ce_evidence(base, evidence)

    entity = final["entities"][0]
    assert retriever.calls == 1
    assert entity["winner"]["identity_key"] == "std2d:BBBBBBBBBBBBBB"
    assert entity["winner"]["entity_score"] == 0.76
    assert entity["winner"]["base_conditioned_multilibrary_rank"] == 2
    assert entity["winner"]["entity_rank"] == 1
    assert entity["decision"]["status"] == "selected_candidate_only"
    assert entity["decision"]["blockers"] == [
        "same_formula_connectivity_requires_orthogonal_confirmation"
    ]
    assert final["metadata"]["entities_reranked_by_source_paired_ce_gate"] == 1
    assert final["invariants"]["missing_reference_ce_treated_as_exact_match"] is False


def test_negative_platform_anchor_retains_v15_baseline():
    base = _payload()
    evidence = build_v15_ce_evidence_bundle(
        base, _rows(), retriever=_Retriever(platform_delta=-0.02)
    )
    final = apply_v15_ce_evidence(base, evidence)

    assert final["entities"][0]["winner"]["identity_key"] == "std2d:AAAAAAAAAAAAAA"
    gate = final["entities"][0]["source_paired_ce_evidence"]["production_gate"]
    assert gate["status"] == "baseline_retained"
    assert "platform_median_delta_below_gate" in gate["failed_gates"]


def test_far_margin_entity_skips_expensive_retrieval():
    base = _payload(margin=0.20)
    retriever = _Retriever()
    evidence = build_v15_ce_evidence_bundle(base, _rows(), retriever=retriever)

    assert retriever.calls == 0
    assert evidence["metadata"]["entities_retrieved"] == 0
    assert evidence["entities"][0]["production_gate"]["reason"] == "v15_margin_above_gate"


def test_evidence_cannot_be_rebound_to_another_base_ranking():
    base = _payload()
    evidence = build_v15_ce_evidence_bundle(
        base, _rows(), retriever=_Retriever()
    )
    changed = copy.deepcopy(base)
    changed["entities"][0]["candidates"][0]["entity_score"] = 0.79

    with pytest.raises(ValueError, match="another base ranking"):
        apply_v15_ce_evidence(changed, evidence)


def test_ce_gate_reuses_bound_top40_candidate_recall_without_second_retrieval():
    base = _payload()
    retriever = _Retriever()
    recall = retriever.retrieve_entity(
        "entity-1",
        _rows(),
        existing_candidates=base["entities"][0]["candidates"],
        detection_method_id="qtof-dda",
    )
    recall["entity_id"] = "entity-1"
    base_hash = hashlib.sha256(
        json.dumps(
            base,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode("utf-8")
    ).hexdigest()
    bundle = {
        "metadata": {
            "base_source_semantic_sha256": "a" * 64,
            "base_artifact_sha256": base_hash,
            "detection_method_id": "qtof-dda",
            "library_snapshot_id": "sqlite-test-snapshot",
            "policy_sha256": "c" * 64,
        },
        "entities": [recall],
    }

    evidence = build_v15_ce_evidence_bundle(
        base,
        _rows(),
        candidate_recall_bundle=bundle,
        detection_method_id="qtof-dda",
    )

    assert retriever.calls == 1
    assert evidence["metadata"]["candidate_recall_bundle_reused"] is True
    assert evidence["metadata"]["entities_retrieved"] == 1
