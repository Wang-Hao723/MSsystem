from __future__ import annotations

import json

from fastapi.testclient import TestClient

from web import server


def _row(scan: str, energy: float) -> dict:
    return {
        "spectrum_id": scan,
        "precursor_mz": 301.0712,
        "retention_time": 1.25,
        "collision_energy": energy,
        "ion_mode": "NEG",
        "precursor_feature_id": "event-1",
        "lcms_feature_id": "feature-1",
        "match": {
            "compound_name": "Candidate",
            "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            "formula": "C15H10O7",
            "final_score": 0.73,
            "rank_score": 0.73,
            "confidence": "High",
        },
    }


def test_entity_api_materializes_candidate_independent_graph(tmp_path, monkeypatch):
    task_id = "entity-task"
    legacy_path = tmp_path / f"{task_id}_results.json"
    legacy_path.write_text(
        json.dumps([_row("scan=1", 20.0), _row("scan=2", 40.0)]),
        encoding="utf-8",
    )
    before = legacy_path.read_bytes()
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {
            task_id: {
                "id": task_id,
                "status": "completed",
                "filename": "sample.mzML",
            }
        },
    )
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    client = TestClient(server.app)

    response = client.get(f"/api/tasks/{task_id}/entities?node_limit=2")

    assert response.status_code == 200
    payload = response.json()
    assert payload["summary"]["scan_count"] == 2
    assert payload["summary"]["precursor_event_count"] == 1
    assert payload["summary"]["same_injection_extra_scans_collapsed"] == 1
    assert payload["summary"]["candidate_independent"] is True
    assert payload["feature_consensus_runtime_influence"] is False
    assert payload["nodes_total"] == 4
    assert len(payload["nodes"]) == 2
    assert payload["nodes_has_more"] is True
    assert payload["legacy_results_unchanged"] is True
    assert legacy_path.read_bytes() == before
    assert (tmp_path / f"{task_id}_entity_bundle.json").is_file()


def test_entity_api_can_refuse_lazy_materialization(tmp_path, monkeypatch):
    task_id = "entity-task-no-cache"
    (tmp_path / f"{task_id}_results.json").write_text(
        json.dumps([_row("scan=1", 20.0)]), encoding="utf-8"
    )
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {
            task_id: {
                "id": task_id,
                "status": "completed",
                "filename": "sample.mzML",
            }
        },
    )
    client = TestClient(server.app)

    response = client.get(
        f"/api/tasks/{task_id}/entities?materialize=false"
    )

    assert response.status_code == 404
    assert "not been materialized" in response.json()["detail"]
