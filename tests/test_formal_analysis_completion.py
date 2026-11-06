from copy import deepcopy

import pytest

from runtime.analysis_completion import formal_output_check, assert_formal_outputs


def summaries():
    return {"primary_identity": {"status": "available", "production_identity_selector": True},
            "identification_records": {"status": "available"},
            "claim_contract": {"status": "available"},
            "compound_conclusions": {"status": "available", "main_result_count": 0}}


def test_zero_candidate_formal_run_can_complete():
    check = formal_output_check(summaries(), primary_artifact_current=True)
    assert check["status"] == "passed"
    assert_formal_outputs(check)


@pytest.mark.parametrize("module", ["primary_identity", "identification_records", "claim_contract", "compound_conclusions"])
def test_failed_required_module_blocks_completion_even_with_available_audit(module):
    value = summaries()
    value[module] = {"status": "failed", "error": "missing registered dependency"}
    value["v11_audit"] = {"status": "available", "production_identity_selector": True}
    check = formal_output_check(value, primary_artifact_current=True)
    with pytest.raises(RuntimeError, match=module):
        assert_formal_outputs(check)
    assert check["audit_artifact_substitution_allowed"] is False


def test_stale_or_audit_only_artifact_blocks_completion():
    check = formal_output_check(summaries(), primary_artifact_current=False)
    with pytest.raises(RuntimeError, match="absent or stale"):
        assert_formal_outputs(check)
    value = summaries()
    value["primary_identity"]["production_identity_selector"] = False
    with pytest.raises(RuntimeError, match="production identity selector"):
        assert_formal_outputs(formal_output_check(value, primary_artifact_current=True))


def test_malformed_primary_has_structured_failure():
    value = summaries()
    value["primary_identity"] = []
    check = formal_output_check(value, primary_artifact_current=False)
    assert check["status"] == "failed"
    with pytest.raises(RuntimeError, match="primary_identity"):
        assert_formal_outputs(check)


def test_missing_dependency_marks_native_task_failed_before_parsing(monkeypatch, tmp_path):
    from web import server
    from runtime import verified_sql_retrieval
    task = {"id": "formal-missing", "filename": "sample.mgf", "file_path": str(tmp_path / "sample.mgf"),
            "status": "queued", "params": {}, "results": []}
    monkeypatch.setattr(server, "_tasks", {task["id"]: deepcopy(task)})
    monkeypatch.setattr(server, "_task_cancel_events", {})
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    monkeypatch.setattr(server, "_sse_broadcast", lambda *a, **kw: None)
    monkeypatch.setattr(server, "_parse_spectrum_file", lambda *a, **kw: pytest.fail("parse must not run"))
    monkeypatch.setattr(verified_sql_retrieval, "sql_metadata_dependency_paths", lambda *a: [tmp_path / "missing_bridge.jsonl"])
    server._run_analysis(task["id"])
    actual = server._tasks[task["id"]]
    assert actual["status"] == "failed"
    assert "missing_bridge.jsonl" in actual["error"]
    assert actual["primary_dependency_preflight"]["status"] == "failed"


def test_formal_check_precedes_native_completed_transition():
    import ast
    from pathlib import Path
    source = Path(__file__).parents[1] / "web/server.py"
    module = ast.parse(source.read_text())
    worker = next(n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == "_run_analysis")
    checks = [n.lineno for n in ast.walk(worker) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "assert_formal_outputs"]
    completed = [n.lineno for n in ast.walk(worker) if isinstance(n, ast.Attribute) and n.attr == "COMPLETED"]
    assert len(checks) == 1 and completed and checks[0] < min(completed)
