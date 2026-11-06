from __future__ import annotations

from fastapi.testclient import TestClient

from web import server


def test_minimal_task_projection_has_a_fixed_small_field_set(monkeypatch):
    task = {
        "id": "task-pos",
        "filename": "sample.mzML",
        "status": "completed",
        "ion_mode": "POS",
        "upload_time": "2026-08-31T00:00:00Z",
        "progress": 100,
        "results": [{"raw_peaks": [[100.0, 1.0]]}] * 1000,
        "flagship_summary": {"large": "x" * 100_000},
        "untargeted_workflow": {"large": "y" * 100_000},
    }
    monkeypatch.setattr(server, "_tasks", {task["id"]: task})
    monkeypatch.setattr(
        server,
        "_full_queue_snapshot",
        lambda: {"queued": [], "reorder_supported": True},
    )

    response = TestClient(server.app).get(
        "/api/tasks?status=completed&projection=minimal"
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["projection"] == "minimal"
    assert payload["count"] == 1
    assert set(payload["tasks"][0]) == {
        "id",
        "filename",
        "status",
        "ion_mode",
        "upload_time",
    }
    assert len(response.content) < 500
