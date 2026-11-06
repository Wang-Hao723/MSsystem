import ast
from pathlib import Path

import pytest

from runtime.analysis_completion import assert_formal_outputs, formal_output_check


def _available_outputs():
    return {
        "primary_identity": {"status": "available", "production_identity_selector": True},
        "identification_records": {"status": "available"},
        "claim_contract": {"status": "available"},
        "compound_conclusions": {"status": "available"},
    }


@pytest.mark.parametrize("failed_output", ["primary_identity", "claim_contract", "compound_conclusions"])
def test_required_formal_failure_stops_before_completed_side_effects(failed_output):
    outputs = _available_outputs()
    outputs[failed_output] = {"status": "failed", "error": "synthetic failure"}
    completed_side_effects = []
    with pytest.raises(RuntimeError, match=failed_output):
        check = formal_output_check(outputs, primary_artifact_current=True)
        assert_formal_outputs(check)
        completed_side_effects.extend(["TaskState.COMPLETED", "SSE completed"])
    assert completed_side_effects == []


def test_server_places_formal_gate_before_completed_transition_and_sse():
    source = Path(__file__).parents[1] / "web/server.py"
    module = ast.parse(source.read_text(encoding="utf-8"))
    worker = next(node for node in module.body if isinstance(node, ast.FunctionDef) and node.name == "_run_analysis")
    gate_lines = [
        node.lineno for node in ast.walk(worker)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        and node.func.id == "assert_formal_outputs"
    ]
    completed_transition_lines = []
    completed_sse_lines = []
    for node in ast.walk(worker):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Name) and node.func.id == "_transition_task":
            if any(isinstance(arg, ast.Attribute) and arg.attr == "COMPLETED" for arg in node.args):
                completed_transition_lines.append(node.lineno)
        if isinstance(node.func, ast.Name) and node.func.id == "_sse_broadcast":
            for arg in node.args:
                if isinstance(arg, ast.Dict):
                    pairs = zip(arg.keys, arg.values)
                    if any(isinstance(key, ast.Constant) and key.value == "status"
                           and isinstance(value, ast.Constant) and value.value == "completed"
                           for key, value in pairs):
                        completed_sse_lines.append(node.lineno)
    assert len(gate_lines) == 1
    assert completed_transition_lines and completed_sse_lines
    assert gate_lines[0] < min(completed_transition_lines + completed_sse_lines)
