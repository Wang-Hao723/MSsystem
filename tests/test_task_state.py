from __future__ import annotations

import pytest

from runtime.task_state import (
    InvalidTaskTransition,
    TaskState,
    initialize_task_state,
    is_active_task_state,
    transition_task,
)


def test_task_state_happy_path_is_versioned_and_auditable():
    task = {"id": "task-1", "status": "pending", "upload_time": "t0"}
    initialize_task_state(task)
    task["active_run_id"] = "run-1"
    transition_task(task, TaskState.QUEUED, reason="user_started", at="t1")
    transition_task(task, TaskState.RUNNING, reason="worker_started", at="t2")
    transition_task(task, TaskState.FINALIZING, reason="matching_done", at="t3")
    transition_task(task, TaskState.COMPLETED, reason="artifacts_committed", at="t4")

    assert task["status"] == "completed"
    assert task["state_version"] == 5
    assert task["run_attempt"] == 1
    assert [event["to"] for event in task["state_history"]] == [
        "pending",
        "queued",
        "running",
        "finalizing",
        "completed",
    ]
    assert task["state_history"][-1]["run_id"] == "run-1"


def test_illegal_transition_is_rejected():
    task = {"id": "task-1", "status": "pending"}
    with pytest.raises(InvalidTaskTransition, match="pending -> completed"):
        transition_task(task, TaskState.COMPLETED, reason="invalid")


def test_cancel_and_retry_path_is_explicit():
    task = {"id": "task-1", "status": "pending"}
    transition_task(task, TaskState.QUEUED, reason="start")
    transition_task(task, TaskState.CANCELLING, reason="cancel_requested")
    transition_task(task, TaskState.CANCELLED, reason="worker_acknowledged")
    transition_task(task, TaskState.QUEUED, reason="retry")
    assert task["run_attempt"] == 2
    assert is_active_task_state(task["status"])
