import asyncio
from copy import deepcopy
import json
import os
from pathlib import Path

import pytest

from runtime.task_reconciliation import (
    reconciled_completed_task,
    verify_completed_task_artifacts,
)


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")


def _completed_run_fixture(tmp_path: Path) -> dict:
    task_id = "task-reconcile-test"
    run_id = "run-7"
    semantic_hash = "sha256:retained-source-semantic-hash"
    source = tmp_path / f"{task_id}_results.json"
    flagship = tmp_path / f"{task_id}_flagship.json"
    compound = tmp_path / f"{task_id}_compound.json"

    results = [
        {"spectrum_id": "scan-1", "match": {"score": 0.91, "name": "A"}},
        {"spectrum_id": "scan-2", "match": {"score": 0.73, "name": "B"}},
    ]
    _write_json(source, results)
    _write_json(
        flagship,
        {
            "metadata": {
                "task_id": task_id,
                "source_row_count": len(results),
                "source_semantic_sha256": semantic_hash,
            },
            "rows": [{"spectrum_id": "scan-1"}],
        },
    )
    source_stat = source.stat()
    flagship_stat = flagship.stat()
    _write_json(
        compound,
        {
            "metadata": {
                "task_id": task_id,
                "source_artifact": str(source),
                "source_artifact_size_bytes": source_stat.st_size,
                "source_artifact_mtime_ns": source_stat.st_mtime_ns,
                "primary_identity_artifact": str(flagship),
                "primary_identity_artifact_size_bytes": flagship_stat.st_size,
                "primary_identity_artifact_mtime_ns": flagship_stat.st_mtime_ns,
                "source_row_count": len(results),
                "source_semantic_sha256": semantic_hash,
            },
            "rows": [{"feature": "feature-1"}],
        },
    )

    return {
        "id": task_id,
        "status": "failed",
        "state_version": 7,
        "active_run_id": run_id,
        "result_run_id": run_id,
        "total_spectra": len(results),
        "progress": 86,
        "error": "late status overwrite",
        "message": "worker exit was observed after artifact commit",
        "results": results,
        "compound_conclusion_summary": {
            "artifact": str(compound),
            "source_row_count": len(results),
            "source_semantic_sha256": semantic_hash,
        },
        "flagship_summary": {
            "artifact": str(flagship),
            "source_row_count": len(results),
            "source_semantic_sha256": semantic_hash,
        },
        "state_history": [
            {
                "version": 6,
                "from": "running",
                "to": "finalizing",
                "at": "2026-09-22T07:59:00+00:00",
                "reason": "analysis_finished",
                "run_id": run_id,
            },
            {
                "version": 7,
                "from": "finalizing",
                "to": "completed",
                "at": "2026-09-22T08:00:00+00:00",
                "reason": "v11_flagship_and_supporting_artifacts_committed",
                "run_id": run_id,
            },
        ],
    }


def test_verified_completed_run_is_reconciled_without_changing_results_or_scores(tmp_path):
    task = _completed_run_fixture(tmp_path)
    original_history = deepcopy(task["state_history"])
    original_results = task["results"]
    original_scores = [row["match"]["score"] for row in original_results]

    evidence = verify_completed_task_artifacts(task, tmp_path)
    repaired = reconciled_completed_task(task, evidence)

    assert evidence["source_row_count"] == 2
    assert evidence["score_rank_or_evidence_changes"] is False
    assert repaired["status"] == "completed"
    assert repaired["state_version"] == 8
    assert repaired["progress"] == 100
    assert repaired["end_time"] == original_history[-1]["at"]
    assert repaired["error"] is None
    assert repaired["results"] is original_results
    assert [row["match"]["score"] for row in repaired["results"]] == original_scores
    assert repaired["compound_conclusion_summary"] == task["compound_conclusion_summary"]
    assert repaired["flagship_summary"] == task["flagship_summary"]
    assert repaired["state_history"][:-1] == original_history
    assert repaired["state_history"][-1]["prior_error"] == task["error"]
    assert repaired["state_history"][-1]["prior_message"] == task["message"]
    assert repaired["state_history"][-1]["verification"] is evidence
    assert repaired["state_reconciliations"] == [repaired["state_history"][-1]]
    assert task["status"] == "failed"
    assert task["error"] == "late status overwrite"
    assert task["state_history"] == original_history


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("reason", "some_other_completion"),
        ("run_id", "another-run"),
        ("version", 6),
    ],
)
def test_completed_ledger_mismatch_is_rejected(tmp_path, field, value):
    task = _completed_run_fixture(tmp_path)
    task["state_history"][-1][field] = value

    with pytest.raises(ValueError, match="Completed ledger"):
        verify_completed_task_artifacts(task, tmp_path)


def test_source_file_change_invalidates_registered_stat_binding(tmp_path):
    task = _completed_run_fixture(tmp_path)
    source = tmp_path / f"{task['id']}_results.json"
    source.write_text(source.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="source_artifact_(size_bytes|mtime_ns)"):
        verify_completed_task_artifacts(task, tmp_path)


def test_task_change_after_verification_is_rejected(tmp_path):
    task = _completed_run_fixture(tmp_path)
    evidence = verify_completed_task_artifacts(task, tmp_path)
    task["message"] = "a newer worker observation"

    with pytest.raises(ValueError, match="Task changed during verification"):
        reconciled_completed_task(task, evidence)


def test_same_size_same_mtime_artifact_change_is_rejected(tmp_path):
    task = _completed_run_fixture(tmp_path)
    evidence = verify_completed_task_artifacts(task, tmp_path)
    source = tmp_path / f"{task['id']}_results.json"
    stat = source.stat()
    changed = source.read_text(encoding="utf-8").replace("0.91", "0.92", 1)
    source.write_text(changed, encoding="utf-8")
    os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert source.stat().st_size == stat.st_size

    with pytest.raises(ValueError, match="Artifact changed"):
        reconciled_completed_task(task, evidence)


def test_api_persistence_failure_restores_failed_task(monkeypatch):
    from runtime import task_reconciliation
    from web import server

    task_id = "task-reconcile-persist-failure"
    original = {
        "id": task_id,
        "status": "failed",
        "state_version": 3,
        "error": "retained failure evidence",
        "message": "retained worker message",
        "results": [{"spectrum_id": "scan-1", "match": {"score": 0.88}}],
    }
    repaired = {**original, "status": "completed", "state_version": 4, "error": None}
    broadcasts = []

    monkeypatch.setitem(server._tasks, task_id, original)
    monkeypatch.setattr(
        task_reconciliation,
        "verify_completed_task_artifacts",
        lambda task, directory: {"task_id": task_id, "verified": True},
    )
    monkeypatch.setattr(
        task_reconciliation,
        "reconciled_completed_task",
        lambda task, evidence: repaired,
    )
    monkeypatch.setattr(server, "_save_tasks", lambda: False)
    monkeypatch.setattr(server, "_sse_broadcast", lambda *args: broadcasts.append(args))

    request = server.TaskStateReconciliationRequest(
        expected_state_version=original["state_version"], dry_run=False
    )
    with pytest.raises(server.HTTPException) as caught:
        asyncio.run(server.reconcile_completed_task_state(task_id, request))

    assert caught.value.status_code == 500
    assert server._tasks[task_id] is original
    assert server._tasks[task_id]["status"] == "failed"
    assert server._tasks[task_id]["error"] == "retained failure evidence"
    assert broadcasts == []
