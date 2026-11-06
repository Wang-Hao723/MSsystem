from __future__ import annotations

import threading

from fastapi.testclient import TestClient

from runtime.task_state import initialize_task_state
from web import server


def test_pending_task_can_be_cancelled_and_state_history_is_persistable(monkeypatch):
    task = {
        "id": "pending-task",
        "filename": "sample.mzML",
        "status": "pending",
        "progress": 0,
        "upload_time": "2026-08-03T00:00:00+00:00",
    }
    initialize_task_state(task)
    monkeypatch.setattr(server, "_tasks", {"pending-task": task})
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    monkeypatch.setattr(server, "_sse_broadcast", lambda *_args, **_kwargs: None)
    client = TestClient(server.app)

    cancelled = client.post("/api/tasks/pending-task/cancel")
    history = client.get("/api/tasks/pending-task/state-history")

    assert cancelled.status_code == 200
    assert cancelled.json()["status"] == "cancelled"
    assert cancelled.json()["cancellation_requested"] is False
    assert history.status_code == 200
    assert [event["to"] for event in history.json()["history"]] == [
        "pending",
        "cancelled",
    ]


def test_running_task_cancellation_sets_cooperative_event(monkeypatch):
    task = {
        "id": "running-task",
        "filename": "sample.mzML",
        "status": "running",
        "progress": 50,
        "upload_time": "2026-08-03T00:00:00+00:00",
    }
    initialize_task_state(task)
    events: dict[str, threading.Event] = {}
    monkeypatch.setattr(server, "_tasks", {"running-task": task})
    monkeypatch.setattr(server, "_task_cancel_events", events)
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    monkeypatch.setattr(server, "_sse_broadcast", lambda *_args, **_kwargs: None)

    response = TestClient(server.app).post("/api/tasks/running-task/cancel")

    assert response.status_code == 200
    assert response.json()["status"] == "cancelling"
    assert events["running-task"].is_set()


def test_finalizing_task_refuses_unsafe_cancellation(monkeypatch):
    task = {
        "id": "finalizing-task",
        "filename": "sample.mzML",
        "status": "finalizing",
        "progress": 98,
        "upload_time": "2026-08-03T00:00:00+00:00",
    }
    initialize_task_state(task)
    monkeypatch.setattr(server, "_tasks", {"finalizing-task": task})

    response = TestClient(server.app).post("/api/tasks/finalizing-task/cancel")

    assert response.status_code == 409
    assert "committing result artifacts" in response.json()["detail"]
