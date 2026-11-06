"""Native enum values and saved JSON strings describe the same task state."""
import importlib.util
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[1] / "runtime/task_state.py"
SPEC = importlib.util.spec_from_file_location("verified_task_state_input", SOURCE)
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


@pytest.mark.parametrize("state", list(M.TaskState))
def test_enum_and_saved_string_have_identical_state_and_activity(state):
    assert M.task_state(state) is M.task_state(state.value) is state
    assert M.is_active_task_state(state) == M.is_active_task_state(state.value)
    assert M.is_active_task_state(state) == (state in M.ACTIVE_TASK_STATES)


@pytest.mark.parametrize("state", list(M.TaskState))
def test_enum_initialization_preserves_serialized_history_value(state):
    task = {"status": state}
    M.initialize_task_state(task, at="2026-09-15T00:00:00+00:00")
    assert task["state_history"][0]["to"] == state.value


@pytest.mark.parametrize("source,target", [(s,t) for s, targets in M.ALLOWED_TRANSITIONS.items() for t in targets])
def test_allowed_transition_from_enum_keeps_json_status_string(source, target):
    enum_task, string_task = {"status": source}, {"status": source.value}
    for task in (enum_task, string_task):
        M.transition_task(task, target, reason="test", at="2026-09-15T00:00:00+00:00")
    assert enum_task == string_task and enum_task["status"] == target.value
    assert type(enum_task["status"]) is str


@pytest.mark.parametrize("value", [None, True, 1, "TaskState.RUNNING", "unknown"])
def test_unrecognized_states_still_fail(value):
    with pytest.raises(M.InvalidTaskTransition):
        M.task_state(value)
    assert M.is_active_task_state(value) is False
