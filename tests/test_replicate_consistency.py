from __future__ import annotations

import json

from fastapi.testclient import TestClient

from pipeline.replicate_consistency import (
    build_replicate_consistency,
    parse_replicate_filename,
)
from web import server


def _compound(identity: str, name: str, *, high: bool, main: bool = True) -> dict:
    return {
        "identity_key": identity,
        "identity_scope": "isomer_group" if identity.startswith("iso") else "2d_connectivity",
        "formula": "C6H12O6" if identity.startswith("iso") else "C10H12O2",
        "display_name": name,
        "lead_candidate": {
            "compound_name": name,
            "rank_score": 0.9 if high else 0.7,
            "cosine": 0.95 if high else 0.8,
        },
        "conclusion": {
            "key": (
                "isomer" if high and identity.startswith("iso")
                else "high" if high
                else "candidate" if main
                else "insufficient"
            ),
            "counted_as_sample_compound": main,
            "high_confidence_group_detected": high,
        },
    }


def _write_artifact(task_dir, task_id: str, compounds: list[dict]) -> None:
    from pipeline.compound_conclusions import compound_conclusion_artifact_path
    compound_conclusion_artifact_path(task_dir, task_id).write_text(
        json.dumps({"metadata": {}, "compounds": compounds}),
        encoding="utf-8",
    )


def test_explicit_terminal_repeat_number_is_required():
    parsed = parse_replicate_filename("F001_NEG_3.mzML")
    assert parsed == {
        "filename": "F001_NEG_3.mzML",
        "base": "F001_NEG",
        "base_key": "f001_neg",
        "replicate_number": 3,
        "extension": ".mzml",
    }
    assert parse_replicate_filename("F001_NEG.mzML") is None
    assert parse_replicate_filename("F001_NEG_three.mzML") is None


def test_three_injection_consistency_separates_per_run_high_from_robust_high(
    tmp_path,
):
    params = {"db_mode": "full", "ion_mode": "auto", "blank_filter_policy": "soft_evidence"}
    tasks = {
        "t1": {"filename": "F001_NEG_1.mzML", "status": "completed", "params": params},
        "t2": {"filename": "F001_NEG_2.mzML", "status": "completed", "params": params},
        "t3": {"filename": "F001_NEG_3.mzML", "status": "completed", "params": params},
        "different": {
            "filename": "F001_NEG_4.mzML",
            "status": "completed",
            "params": {**params, "db_mode": "small"},
        },
    }
    _write_artifact(tmp_path, "t1", [
        _compound("A", "A", high=True),
        _compound("iso:B", "B lead", high=True),
        _compound("D", "D", high=False),
    ])
    _write_artifact(tmp_path, "t2", [
        _compound("A", "A", high=True),
        _compound("C", "C", high=True),
        _compound("D", "D", high=False),
    ])
    _write_artifact(tmp_path, "t3", [
        _compound("A", "A", high=True),
        _compound("iso:B", "B alternate lead", high=True),
    ])
    _write_artifact(tmp_path, "different", [_compound("X", "X", high=True)])

    payload = build_replicate_consistency("t1", tasks, tmp_path)

    assert payload["status"] == "available"
    assert payload["injection_count"] == 3
    assert payload["robustness_threshold"] == 2
    assert [item["high_group_count"] for item in payload["injections"]] == [2, 2, 2]
    assert payload["summary"] == {
        "high_union_count": 3,
        "robust_high_count": 2,
        "all_injection_high_count": 1,
        "single_injection_high_count": 1,
        "main_union_count": 4,
        "robust_main_count": 3,
    }
    by_identity = {item["identity_key"]: item for item in payload["compounds"]}
    assert by_identity["A"]["high_injection_count"] == 3
    assert by_identity["A"]["all_injection_high"] is True
    assert by_identity["iso:B"]["high_injection_count"] == 2
    assert by_identity["iso:B"]["robust_high"] is True
    assert by_identity["C"]["single_injection_high_only"] is True
    assert payload["ranking_changed"] is False
    assert payload["confidence_changed"] is False


def test_task_without_compatible_repeat_returns_honest_unavailable(tmp_path):
    tasks = {
        "single": {
            "filename": "ordinary_sample.mzML",
            "status": "completed",
            "params": {},
        }
    }
    payload = build_replicate_consistency("single", tasks, tmp_path)
    assert payload["status"] == "unavailable"
    assert payload["reason"] == "filename_has_no_explicit_replicate_suffix"
    assert payload["ranking_changed"] is False


def test_replicate_consistency_api_exposes_projection(monkeypatch, tmp_path):
    params = {"db_mode": "full"}
    tasks = {
        "api1": {
            "id": "api1", "filename": "sample_1.mzML",
            "status": "completed", "params": params,
        },
        "api2": {
            "id": "api2", "filename": "sample_2.mzML",
            "status": "completed", "params": params,
        },
    }
    _write_artifact(tmp_path, "api1", [_compound("A", "A", high=True)])
    _write_artifact(tmp_path, "api2", [_compound("A", "A", high=True)])
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(server, "_save_tasks", lambda: None)

    with TestClient(server.app) as client:
        monkeypatch.setattr(server, "_tasks", tasks)
        response = client.get("/api/tasks/api1/replicate-consistency")

    assert response.status_code == 200
    payload = response.json()
    assert payload["status"] == "available"
    assert payload["summary"]["robust_high_count"] == 1
    assert payload["ranking_changed"] is False
