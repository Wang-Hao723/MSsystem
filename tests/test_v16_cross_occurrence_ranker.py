from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline.v16_cross_occurrence_ranker import (
    ALGORITHM_ID,
    apply_consensus_to_payload,
    build_consensus_plan,
    load_policy,
    project_entity,
)
from web import v16_cross_occurrence_api as api


X = "AAAAAAAAAAAAAA"
Y = "BBBBBBBBBBBBBB"


def _event(task, *, baseline, direct, formula="C10H10O2", mz=101.0, rt=2.0):
    order = [baseline] + [value for value in (X, Y) if value != baseline]
    return {
        "event_id": f"{task}::entity-{task}",
        "task_id": task,
        "entity_id": f"entity-{task}",
        "polarity": "NEG",
        "measured_precursor_mz": mz,
        "retention_time_min": rt,
        "baseline_formula": formula,
        "baseline_top1_connectivity": baseline,
        "native_same_formula_candidates": [
            {"connectivity": value, "formula": formula, "rank": rank}
            for rank, value in enumerate(order, 1)
        ],
        "member_result_indices": [0],
        "direct_top1_connectivity": direct,
        "candidate_scores": {X: 0.4, Y: 0.6},
        "baseline_score_available": True,
    }


def _candidate(identity, rank, score, formula="C10H10O2"):
    return {
        "identity_key": f"std2d:{identity}",
        "entity_rank": rank,
        "entity_score": score,
        "vetoed": False,
        "candidate": {
            "formula": formula,
            "inchikey": f"{identity}-UHFFFAOYSA-N",
        },
    }


def _payload(task, baseline):
    other = Y if baseline == X else X
    candidates = [_candidate(baseline, 1, 0.7), _candidate(other, 2, 0.6)]
    return {
        "metadata": {
            "task_id": task,
            "algorithm_id": "mssystem-v16-replicate-audit-multilibrary-1.0.0",
            "policy_sha256": "test",
            "identity_resolution_status": "available",
            "entity_semantic_sha256": "before",
        },
        "policy": {},
        "entities": [
            {
                "entity_id": f"entity-{task}",
                "member_result_indices": [0],
                "representative_search_precursor_mz": 101.0,
                "representative_retention_time": 2.0,
                "ion_modes": ["NEG"],
                "winner": dict(candidates[0]),
                "leading_candidate_when_abstained": None,
                "candidates": candidates,
                "decision": {
                    "status": "selected_isomer_group",
                    "scope": "molecular_formula_isomer_group",
                    "winner_margin": 0.1,
                    "ranking_lead_identity_key": f"std2d:{baseline}",
                },
                "identity_resolution": {
                    "status": "unresolved_same_formula_isomer_group",
                    "reportable_scope": "molecular_formula_isomer_group",
                    "ranking_lead_identity_key": f"std2d:{baseline}",
                    "competitive_candidate_identity_keys": [
                        f"std2d:{baseline}", f"std2d:{other}"
                    ],
                },
                "scan_top1_audit_baseline": {"identity_key": f"std2d:{baseline}"},
                "v14_audit_baseline": {"identity_key": f"std2d:{baseline}"},
            }
        ],
        "invariants": {},
    }


def test_locked_policy_is_controlled_top5_same_formula_opt_in():
    policy = load_policy()
    assert policy["algorithm_id"] == ALGORITHM_ID
    assert policy["deployment_state"] == "controlled_opt_in"
    assert policy["candidate_scoring"]["maximum_native_candidates"] == 5
    assert policy["candidate_scoring"]["expanded_candidates_can_replace_winner"] is False
    assert policy["peer_rule"]["target_direct_top1_must_equal_peer_winner"] is True


def test_consensus_switch_requires_two_external_votes_and_target_agreement():
    plan = build_consensus_plan(
        [
            _event("t1", baseline=X, direct=Y),
            _event("t2", baseline=Y, direct=Y),
            _event("t3", baseline=Y, direct=Y),
        ],
        load_policy(),
    )
    proposals = {row["task_id"]: row for row in plan["proposals"]}
    assert plan["applied_count"] == 1
    assert proposals["t1"]["applied"] is True
    assert proposals["t1"]["proposed_top1_connectivity"] == Y
    assert proposals["t1"]["same_formula_peer_count"] == 2
    assert proposals["t2"]["applied"] is False


def test_target_is_excluded_and_formula_conflict_abstains():
    target_excluded = build_consensus_plan(
        [
            _event("t1", baseline=X, direct=Y),
            _event("t2", baseline=X, direct=X),
            _event("t3", baseline=X, direct=X),
        ],
        load_policy(),
    )
    assert not next(row for row in target_excluded["proposals"] if row["task_id"] == "t1")["applied"]

    formula_locked = build_consensus_plan(
        [
            _event("t1", baseline=X, direct=Y, formula="C10H10O2"),
            _event("t2", baseline=Y, direct=Y, formula="C11H12O2"),
            _event("t3", baseline=Y, direct=Y, formula="C11H12O2"),
        ],
        load_policy(),
    )
    row = next(item for item in formula_locked["proposals"] if item["task_id"] == "t1")
    assert row["applied"] is False
    assert "insufficient_same_formula_peer_tasks" in row["failure_reasons"]


def test_unscored_target_or_peer_cannot_drive_consensus():
    target_unscored = _event("t1", baseline=X, direct=Y)
    target_unscored["baseline_score_available"] = False
    plan = build_consensus_plan(
        [
            target_unscored,
            _event("t2", baseline=Y, direct=Y),
            _event("t3", baseline=Y, direct=Y),
        ],
        load_policy(),
    )
    row = next(item for item in plan["proposals"] if item["task_id"] == "t1")
    assert row["applied"] is False
    assert "target_baseline_direct_neutral_score_unavailable" in row["failure_reasons"]

    peer_unscored = _event("t2", baseline=Y, direct=Y)
    peer_unscored["baseline_score_available"] = False
    plan = build_consensus_plan(
        [
            _event("t1", baseline=X, direct=Y),
            peer_unscored,
            _event("t3", baseline=Y, direct=Y),
        ],
        load_policy(),
    )
    row = next(item for item in plan["proposals"] if item["task_id"] == "t1")
    assert row["applied"] is False
    assert row["vote_count"] == 1


def test_projection_rejects_more_than_native_top5():
    entity = _payload("t1", X)["entities"][0]
    entity["candidates"] = [
        _candidate(chr(65 + index) * 14, index + 1, 1.0 - index / 10)
        for index in range(6)
    ]
    entity["winner"] = dict(entity["candidates"][0])
    with pytest.raises(ValueError, match="Top-5"):
        project_entity("t1", entity)


def test_apply_changes_order_not_scores_or_decision_scope():
    payload = _payload("t1", X)
    lead_semantics = "ranking_lead_within_unresolved_isomer_group"
    payload["entities"][0]["winner"]["reporting_semantics"] = lead_semantics
    payload["entities"][0]["candidates"][0]["reporting_semantics"] = lead_semantics
    proposal = build_consensus_plan(
        [
            _event("t1", baseline=X, direct=Y),
            _event("t2", baseline=Y, direct=Y),
            _event("t3", baseline=Y, direct=Y),
        ],
        load_policy(),
    )["proposals"]
    updated, changed = apply_consensus_to_payload(
        payload, proposal, cohort_id="XCR-test", policy=load_policy()
    )
    entity = updated["entities"][0]
    assert changed == 1
    assert entity["winner"]["identity_key"] == f"std2d:{Y}"
    assert [row["identity_key"] for row in entity["candidates"]] == [f"std2d:{Y}", f"std2d:{X}"]
    assert sorted(row["entity_score"] for row in entity["candidates"]) == [0.6, 0.7]
    assert entity["candidates"][0]["reporting_semantics"] == lead_semantics
    assert "reporting_semantics" not in entity["candidates"][1]
    assert entity["decision"]["scope"] == "molecular_formula_isomer_group"
    assert entity["cross_occurrence_consensus"]["candidate_score_changed"] is False


class _Retriever:
    def __init__(self, path: Path):
        self.db_path = path
        self.library_snapshot_id = "sqlite-test"
        self.policy = {
            "spectral_views": [
                {"name": "library_match"},
                {"name": "sensitive"},
            ]
        }


def test_service_activates_with_backups_then_rolls_back(tmp_path, monkeypatch):
    policy = load_policy()
    task_ids = ["t1", "t2", "t3"]
    originals = {}
    for task_id, baseline in zip(task_ids, [X, Y, Y]):
        artifact = tmp_path / f"{task_id}_v16_flagship.json"
        artifact.write_text(json.dumps(_payload(task_id, baseline)), encoding="utf-8")
        (tmp_path / f"{task_id}_results.json").write_text(
            json.dumps([{"instrument_type": "qtof"}]), encoding="utf-8"
        )
        originals[task_id] = artifact.read_bytes()
    database = tmp_path / "library.sqlite"
    database.write_bytes(b"read-only-test-library")

    def fake_score(event, rows, retriever, selected_policy):
        assert selected_policy == policy
        direct = Y
        return {
            "direct_top1_connectivity": direct,
            "candidate_scores": {X: 0.4, Y: 0.6},
            "scored_pair_count": 2,
            "baseline_score_available": True,
            "single_candidate_shortcut": False,
        }

    monkeypatch.setattr(api, "score_event_direct_neutral", fake_score)
    service = api.CrossOccurrenceService(
        task_directory=lambda: tmp_path,
        task_lookup=lambda task_id: {
            "id": task_id,
            "status": "completed",
            "params": {"rt_method_id": "LC-A"},
        },
        artifact_path=lambda directory, task_id: directory / f"{task_id}_v16_flagship.json",
        artifact_current=lambda path: True,
        retriever_factory=lambda: _Retriever(database),
        runtime_snapshot=lambda: "v16-test",
        runtime_idle=lambda: True,
    )
    receipt = service.compute(
        api.CrossOccurrenceRequest(
            task_ids=task_ids,
            same_lc_method_confirmed=True,
            activate=True,
        )
    )
    assert receipt["status"] == "activated"
    assert receipt["applied_count"] == 1
    active = json.loads((tmp_path / "t1_v16_flagship.json").read_text())
    assert active["entities"][0]["winner"]["identity_key"] == f"std2d:{Y}"
    assert Path(receipt["artifacts"]["t1"]["backup"]).is_file()

    result = service.rollback(
        api.CrossOccurrenceRollbackRequest(
            cohort_id=receipt["cohort_id"], confirm_rollback=True
        )
    )
    assert result["status"] == "rolled_back"
    assert (tmp_path / "t1_v16_flagship.json").read_bytes() == originals["t1"]


def test_service_requires_same_method_confirmation(tmp_path):
    database = tmp_path / "library.sqlite"
    database.write_bytes(b"test")
    service = api.CrossOccurrenceService(
        task_directory=lambda: tmp_path,
        task_lookup=lambda task_id: {"status": "completed", "params": {}},
        artifact_path=lambda directory, task_id: directory / task_id,
        artifact_current=lambda path: True,
        retriever_factory=lambda: _Retriever(database),
        runtime_snapshot=lambda: "v16-test",
        runtime_idle=lambda: True,
    )
    with pytest.raises(api.HTTPException, match="same-LC-method"):
        service.compute(
            api.CrossOccurrenceRequest(
                task_ids=["t1", "t2", "t3"],
                same_lc_method_confirmed=False,
            )
        )


def test_runtime_top5_replay_reproduces_frozen_gain(tmp_path):
    from eval.v16_cross_occurrence_consensus_exploration.validate_runtime_top5 import run

    result = run(tmp_path / "runtime_top5_replay.json")
    assert result["event_count"] == 184
    assert (result["baseline_successes"], result["variant_successes"]) == (131, 140)
    assert (result["gained"], result["lost"]) == (9, 0)
    assert result["identity_cluster"]["bootstrap_95pct_interval"][0] > 0
