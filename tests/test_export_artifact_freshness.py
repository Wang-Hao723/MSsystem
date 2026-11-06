"""A saved task must not export conclusions through stale SQL metadata."""
import json

import pytest
from fastapi import HTTPException


@pytest.fixture
def context(monkeypatch, tmp_path):
    import web.server as server
    tid = "persisted-export-freshness"
    task = {"id": tid, "filename": "test.mzML", "status": "completed",
            "compound_conclusion_summary": {"stale": True}}
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    (tmp_path / f"{tid}_results.json").write_text("[]")
    primary = tmp_path / "primary.json"
    compound = tmp_path / "compound.json"
    primary.write_text("{}")
    compound.write_text("{}")
    monkeypatch.setattr(server, "_flagship_artifact_path", lambda *a: primary)
    monkeypatch.setattr(server, "_compound_conclusion_artifact_path", lambda *a: compound)
    monkeypatch.setattr(server, "_get_task_results", lambda *a: [])
    monkeypatch.setattr(server, "_apply_final_compound_count_summary", lambda *a: None)
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    state = {"primary": True, "compound": True, "events": []}

    def primary_current(*args):
        state["events"].append("check_primary")
        return state["primary"]

    def compound_current(*args):
        state["events"].append("check_compound")
        return state["compound"]

    def reader(*args, **kwargs):
        state["events"].append("read_compound")
        return {"metadata": {"current": True}, "compounds": [{"compound_id": "current"}], "has_more": False}

    monkeypatch.setattr(server, "_flagship_artifact_is_current", primary_current)
    monkeypatch.setattr(server, "_compound_conclusion_artifact_is_current", compound_current)
    monkeypatch.setattr(server, "_read_compound_conclusion_artifact", reader)
    return server, tid, task, state, primary, compound


def test_stale_primary_cannot_hide_behind_current_compound(context, monkeypatch):
    server, tid, task, state, *_ = context
    state["primary"] = False

    def invalid_sql(*args):
        raise ValueError("SQL metadata source hash mismatch")

    monkeypatch.setattr(server, "_materialize_current_v15_flagship", invalid_sql)
    with pytest.raises(HTTPException) as error:
        server._compound_conclusions_for_export(tid, task)
    assert error.value.status_code == 500
    assert "SQL metadata source hash mismatch" in error.value.detail
    assert "read_compound" not in state["events"]
    assert state["events"] == ["check_primary"]
    assert task["compound_conclusion_summary"] == {"stale": True}


def test_rebuild_verifies_primary_then_compound_before_reading(context, monkeypatch):
    server, tid, task, state, *_ = context
    state.update(primary=False, compound=False)

    def rebuild_primary(*args):
        state["events"].append("rebuild_primary")
        state["primary"] = True
        return {"primary_current": True}

    def rebuild_compound(*args, **kwargs):
        state["events"].append("rebuild_compound")
        assert state["primary"] is True
        state["compound"] = True
        return {"current": True}

    monkeypatch.setattr(server, "_materialize_current_v15_flagship", rebuild_primary)
    monkeypatch.setattr(server, "_write_compound_conclusion_artifact", rebuild_compound)
    assert server._compound_conclusions_for_export(tid, task) == [{"compound_id": "current"}]
    assert state["events"] == ["check_primary", "rebuild_primary", "check_primary",
                               "check_compound", "rebuild_compound", "check_compound", "read_compound"]
    assert task["flagship_summary"] == {"primary_current": True}


@pytest.mark.parametrize("failed_stage", ["primary", "compound"])
def test_rebuild_success_return_does_not_replace_freshness_check(context, monkeypatch, failed_stage):
    server, tid, task, state, *_ = context
    state[failed_stage] = False
    monkeypatch.setattr(server, "_materialize_current_v15_flagship", lambda *a: {"pretends_success": True})
    monkeypatch.setattr(server, "_write_compound_conclusion_artifact", lambda *a, **k: {"pretends_success": True})
    with pytest.raises(HTTPException) as error:
        server._compound_conclusions_for_export(tid, task)
    assert error.value.status_code == 500
    assert "remains stale" in error.value.detail
    assert "read_compound" not in state["events"]


def test_compound_rebuild_failure_does_not_return_old_summary(context, monkeypatch):
    server, tid, task, state, *_ = context
    state["compound"] = False

    def broken(*args, **kwargs):
        raise RuntimeError("primary SQL dependency could not be bound")

    monkeypatch.setattr(server, "_write_compound_conclusion_artifact", broken)
    with pytest.raises(HTTPException) as error:
        server._compound_scope_summary_for_export(tid, task)
    assert error.value.status_code == 500
    assert task["compound_conclusion_summary"] == {"stale": True}
    assert "read_compound" not in state["events"]


def test_missing_compound_file_never_returns_legacy_none(context):
    server, tid, task, state, _, compound = context
    compound.unlink()
    with pytest.raises(HTTPException) as error:
        server._compound_conclusions_for_export(tid, task)
    assert error.value.status_code == 500


def test_compound_read_error_is_not_silently_replaced_by_scan_rows(context, monkeypatch):
    server, tid, task, state, *_ = context

    def unreadable(*args, **kwargs):
        raise json.JSONDecodeError("invalid artifact", "{", 0)

    monkeypatch.setattr(server, "_read_compound_conclusion_artifact", unreadable)
    with pytest.raises(HTTPException) as error:
        server._compound_conclusions_for_export(tid, task)
    assert error.value.status_code == 500


def test_missing_current_metadata_does_not_reuse_old_summary(context, monkeypatch):
    server, tid, task, state, *_ = context
    monkeypatch.setattr(server, "_read_compound_conclusion_artifact", lambda *a, **k: {"compounds": []})
    with pytest.raises(HTTPException) as error:
        server._compound_scope_summary_for_export(tid, task)
    assert error.value.status_code == 500
    assert task["compound_conclusion_summary"] == {"stale": True}


def test_in_memory_legacy_task_is_the_only_none_fallback(context):
    server, tid, task, state, *_ = context
    (server.TASK_DIR / f"{tid}_results.json").unlink()
    assert server._compound_conclusions_for_export(tid, task) is None
    assert server._compound_scope_summary_for_export(tid, task) == {"stale": True}
    assert state["events"] == []
