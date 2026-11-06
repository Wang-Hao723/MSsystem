"""Explicit, auditable task lifecycle transitions for V11."""

from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum
from typing import Any


class TaskState(str, Enum):
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    FINALIZING = "finalizing"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLING = "cancelling"
    CANCELLED = "cancelled"


ACTIVE_TASK_STATES = frozenset(
    {
        TaskState.QUEUED,
        TaskState.RUNNING,
        TaskState.FINALIZING,
        TaskState.CANCELLING,
    }
)
TERMINAL_TASK_STATES = frozenset(
    {TaskState.COMPLETED, TaskState.FAILED, TaskState.CANCELLED}
)
ALLOWED_TRANSITIONS = {
    TaskState.PENDING: {TaskState.QUEUED, TaskState.CANCELLED, TaskState.FAILED},
    # PENDING is also the durable restart-recovery state.  A process restart
    # destroys executor workers but not managed uploads, so active work can be
    # returned to PENDING and re-queued after the requested library is ready.
    TaskState.QUEUED: {
        TaskState.PENDING, TaskState.RUNNING, TaskState.CANCELLING,
        TaskState.FAILED,
    },
    TaskState.RUNNING: {
        TaskState.PENDING, TaskState.FINALIZING, TaskState.CANCELLING,
        TaskState.FAILED,
    },
    TaskState.FINALIZING: {
        TaskState.PENDING, TaskState.COMPLETED, TaskState.FAILED,
    },
    TaskState.COMPLETED: {TaskState.QUEUED},
    TaskState.FAILED: {TaskState.QUEUED},
    TaskState.CANCELLING: {TaskState.CANCELLED, TaskState.FAILED},
    TaskState.CANCELLED: {TaskState.QUEUED},
}


class InvalidTaskTransition(ValueError):
    pass


class TaskCancellationRequested(RuntimeError):
    pass


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def task_state(value: Any) -> TaskState:
    if isinstance(value, TaskState):
        return value
    try:
        return TaskState(str(value or "").strip().lower())
    except ValueError as exc:
        raise InvalidTaskTransition(f"Unknown task state: {value!r}") from exc


def is_active_task_state(value: Any) -> bool:
    try:
        return task_state(value) in ACTIVE_TASK_STATES
    except InvalidTaskTransition:
        return False


def initialize_task_state(task: dict[str, Any], *, at: str | None = None) -> dict[str, Any]:
    state = task_state(task.get("status") or TaskState.PENDING.value)
    timestamp = at or task.get("upload_time") or _now_iso()
    task.setdefault("state_version", 1)
    task.setdefault("run_attempt", 0)
    task.setdefault(
        "state_history",
        [
            {
                "version": int(task["state_version"]),
                "from": None,
                "to": state.value,
                "at": timestamp,
                "reason": "task_registered_or_legacy_state_initialized",
                "run_id": task.get("active_run_id"),
            }
        ],
    )
    return task


def transition_task(
    task: dict[str, Any],
    target: TaskState | str,
    *,
    reason: str,
    at: str | None = None,
    run_id: str | None = None,
) -> dict[str, Any]:
    initialize_task_state(task, at=at)
    current = task_state(task.get("status"))
    target_state = target if isinstance(target, TaskState) else task_state(target)
    if current == target_state:
        return task
    if target_state not in ALLOWED_TRANSITIONS[current]:
        raise InvalidTaskTransition(
            f"Illegal task transition {current.value} -> {target_state.value}"
        )
    timestamp = at or _now_iso()
    version = int(task.get("state_version") or 0) + 1
    effective_run_id = run_id or task.get("active_run_id")
    task["status"] = target_state.value
    task["state_version"] = version
    history = list(task.get("state_history") or [])
    history.append(
        {
            "version": version,
            "from": current.value,
            "to": target_state.value,
            "at": timestamp,
            "reason": str(reason),
            "run_id": effective_run_id,
        }
    )
    task["state_history"] = history[-100:]
    if target_state == TaskState.QUEUED:
        task["queued_time"] = timestamp
        task["run_attempt"] = int(task.get("run_attempt") or 0) + 1
    elif target_state == TaskState.RUNNING:
        task["start_time"] = timestamp
    elif target_state in TERMINAL_TASK_STATES:
        task["end_time"] = timestamp
    return task
