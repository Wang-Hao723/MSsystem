from __future__ import annotations

import threading
from pathlib import Path

from fastapi.testclient import TestClient

from runtime.reorderable_executor import ReorderableSerialExecutor
from runtime.task_state import initialize_task_state
from web import server


def test_reorderable_serial_executor_moves_only_waiting_work():
    executor = ReorderableSerialExecutor(thread_name_prefix="queue-test")
    started = threading.Event()
    release = threading.Event()
    execution_order: list[str] = []

    def first(task_id: str):
        execution_order.append(task_id)
        started.set()
        assert release.wait(timeout=3)

    def record(task_id: str):
        execution_order.append(task_id)
        return task_id

    first_future = executor.submit(first, "first")
    assert started.wait(timeout=2)
    second_future = executor.submit(record, "second")
    third_future = executor.submit(record, "third")

    assert executor.snapshot() == {
        "running": "first",
        "queued": ["second", "third"],
    }
    assert executor.reorder("third", "front")["queued"] == [
        "third", "second",
    ]
    # The running task is never displaced by a queue edit.
    assert executor.snapshot()["running"] == "first"

    release.set()
    first_future.result(timeout=3)
    assert third_future.result(timeout=3) == "third"
    assert second_future.result(timeout=3) == "second"
    executor.shutdown(wait=True)
    assert execution_order == ["first", "third", "second"]


def _queued_task(task_id: str, upload_time: str) -> dict:
    task = {
        "id": task_id,
        "filename": f"{task_id}.mzML",
        "file_path": f"/tmp/{task_id}.mzML",
        "status": "queued",
        "progress": 0,
        "message": "Queued for Full-library serial analysis...",
        "upload_time": upload_time,
        "results": [],
        "params": {"db_mode": "full", "ion_mode": "POS"},
    }
    initialize_task_state(task)
    return task


def test_queue_api_persists_order_and_exposes_positions(monkeypatch):
    executor = ReorderableSerialExecutor(thread_name_prefix="queue-api-test")
    started = threading.Event()
    release = threading.Event()

    def blocker(_task_id: str):
        started.set()
        assert release.wait(timeout=3)

    running_future = executor.submit(blocker, "running")
    assert started.wait(timeout=2)
    queued_one = executor.submit(lambda task_id: task_id, "q1")
    queued_two = executor.submit(lambda task_id: task_id, "q2")

    tasks = {
        "q1": _queued_task("q1", "2026-08-24T01:00:00+00:00"),
        "q2": _queued_task("q2", "2026-08-24T02:00:00+00:00"),
    }
    events: list[tuple[str, dict]] = []
    monkeypatch.setattr(server, "_executor_full", executor)
    monkeypatch.setattr(server, "_tasks", tasks)
    monkeypatch.setattr(server, "_save_tasks", lambda **_kwargs: True)
    monkeypatch.setattr(
        server,
        "_sse_broadcast",
        lambda event_type, payload: events.append((event_type, payload)),
    )

    client = TestClient(server.app)
    try:
        before = client.get("/api/tasks/queue")
        assert before.status_code == 200, before.text
        assert before.json()["order"] == ["q1", "q2"]

        moved = client.post(
            "/api/tasks/queue/reorder",
            json={"task_id": "q2", "action": "front"},
        )
        assert moved.status_code == 200, moved.text
        assert moved.json()["order"] == ["q2", "q1"]
        assert tasks["q2"]["queue_order"] < tasks["q1"]["queue_order"]
        assert events[-1][0] == "queue_reordered"

        listed = client.get("/api/tasks")
        assert listed.status_code == 200, listed.text
        rows = {row["id"]: row for row in listed.json()["tasks"]}
        assert rows["q2"]["queue_position"] == 1
        assert rows["q1"]["queue_position"] == 2
        assert rows["q2"]["queue_reorderable"] is True
    finally:
        client.close()
        release.set()
        running_future.result(timeout=3)
        queued_two.result(timeout=3)
        queued_one.result(timeout=3)
        executor.shutdown(wait=True)


def test_task_page_offers_queue_reorder_controls():
    root = Path(__file__).resolve().parents[1]
    script = (root / "web/static/js/tasks.js").read_text(encoding="utf-8")
    html = (root / "web/static/index.html").read_text(encoding="utf-8")

    assert "/tasks/queue/reorder" in script
    assert "⏫ 置顶" in script
    assert "queue_position" in script
    assert "Tasks.reorderTask" in script
    assert "等待任务可随时调整先后顺序" in html
