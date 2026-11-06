from __future__ import annotations

import threading
from pathlib import Path

from fastapi.testclient import TestClient

from runtime.task_state import initialize_task_state
from web import server


def _isolate_task_storage(monkeypatch, tmp_path):
    upload_dir = tmp_path / "uploads"
    task_dir = tmp_path / "tasks"
    upload_dir.mkdir()
    task_dir.mkdir()
    monkeypatch.setattr(server, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASK_DIR", task_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks_index.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_task_cancel_events", {})
    monkeypatch.setattr(server, "_restart_recovery_task_ids", set())
    monkeypatch.setattr(server, "_task_delete_previews", {})
    monkeypatch.setattr(server, "_task_recycle_purge_previews", {})
    monkeypatch.setattr(server, "_sse_broadcast", lambda *_args, **_kwargs: None)
    return upload_dir, task_dir


def _task(task_id: str, status: str, file_path: Path) -> dict:
    task = {
        "id": task_id,
        "filename": file_path.name,
        "file_path": str(file_path),
        "status": status,
        "progress": 0,
        "upload_time": "2026-08-12T00:00:00+00:00",
    }
    initialize_task_state(task)
    return task


def test_queued_task_is_atomically_removed_with_owned_files(monkeypatch, tmp_path):
    upload_dir, task_dir = _isolate_task_storage(monkeypatch, tmp_path)
    sample = upload_dir / "queued.mzML"
    sample.write_text("<mzML/>", encoding="utf-8")
    result = task_dir / "queued-task_results.json"
    result.write_text("[]", encoding="utf-8")
    artifact_dir = task_dir / "queued-task_artifacts"
    artifact_dir.mkdir()
    (artifact_dir / "partial.tmp").write_text("partial", encoding="utf-8")
    server._tasks["queued-task"] = _task("queued-task", "queued", sample)
    server._task_cancel_events["queued-task"] = threading.Event()

    response = TestClient(server.app).delete("/api/tasks/queued-task")

    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["deleted"] == "queued-task"
    assert payload["previous_status"] == "queued"
    assert payload["removed_from_queue"] is True
    assert payload["recoverable"] is True
    assert payload["deletion_id"]
    assert "queued-task" not in server._tasks
    assert "queued-task" not in server._task_cancel_events
    assert not sample.exists()
    assert not result.exists()
    assert not artifact_dir.exists()
    recycled = tmp_path / "task_recycle_bin" / payload["deletion_id"] / "queued-task"
    assert (recycled / "task.json").is_file()
    assert (recycled / "manifest.json").is_file()

    # A Future already submitted to the serial executor becomes a harmless
    # no-op after the atomic registry removal.
    assert server._run_analysis("queued-task") is None
    assert "queued-task" not in server._tasks


def test_running_task_requires_safe_stop_before_delete(monkeypatch, tmp_path):
    upload_dir, _ = _isolate_task_storage(monkeypatch, tmp_path)
    sample = upload_dir / "running.mzML"
    sample.write_text("<mzML/>", encoding="utf-8")
    server._tasks["running-task"] = _task("running-task", "running", sample)

    response = TestClient(server.app).delete("/api/tasks/running-task")

    assert response.status_code == 409, response.text
    assert "请先点击“停止”" in response.json()["detail"]
    assert server._tasks["running-task"]["status"] == "running"
    assert sample.exists()


def test_recycled_queued_task_can_be_restored_without_auto_running(
    monkeypatch, tmp_path,
):
    upload_dir, task_dir = _isolate_task_storage(monkeypatch, tmp_path)
    sample = upload_dir / "recover-me.mzML"
    sample.write_text("<mzML/>", encoding="utf-8")
    result = task_dir / "recoverable_results.json"
    result.write_text("[]", encoding="utf-8")
    server._tasks["recoverable"] = _task("recoverable", "queued", sample)

    client = TestClient(server.app)
    deleted = client.delete("/api/tasks/recoverable")
    assert deleted.status_code == 200, deleted.text
    deletion_id = deleted.json()["deletion_id"]

    restored = client.post(
        f"/api/task-recycle-bin/{deletion_id}/recoverable/restore"
    )

    assert restored.status_code == 200, restored.text
    assert restored.json()["status"] == "pending"
    assert server._tasks["recoverable"]["status"] == "pending"
    assert sample.is_file()
    assert result.is_file()


def test_recycle_bin_list_exposes_project_metadata_size_and_pagination(
    monkeypatch, tmp_path,
):
    upload_dir, _ = _isolate_task_storage(monkeypatch, tmp_path)
    client = TestClient(server.app)
    for index in range(2):
        sample = upload_dir / f"project-{index}.mzML"
        sample.write_bytes(b"<mzML/>" + bytes([index]) * 32)
        task_id = f"recycled-{index}"
        server._tasks[task_id] = _task(task_id, "completed", sample)
        deleted = client.delete(f"/api/tasks/{task_id}")
        assert deleted.status_code == 200, deleted.text

    first = client.get("/api/task-recycle-bin?limit=1&offset=0")
    assert first.status_code == 200, first.text
    payload = first.json()
    assert payload["total"] == 2
    assert payload["has_more"] is True
    assert payload["summary"]["active_task_count"] == 2
    assert payload["summary"]["size_bytes"] > 0
    assert payload["summary"]["file_count"] >= 6
    assert len(payload["entries"]) == 1
    entry = payload["entries"][0]
    assert entry["filename"].startswith("project-")
    assert entry["previous_status"] == "completed"
    assert entry["recoverable"] is True
    assert entry["safe_to_purge"] is True
    assert entry["size_bytes"] > 0

    second = client.get("/api/task-recycle-bin?limit=1&offset=1")
    assert second.status_code == 200, second.text
    assert second.json()["has_more"] is False
    assert second.json()["entries"][0]["task_id"] != entry["task_id"]


def test_empty_recycle_bin_uses_state_bound_preview_and_writes_receipt(
    monkeypatch, tmp_path,
):
    upload_dir, _ = _isolate_task_storage(monkeypatch, tmp_path)
    client = TestClient(server.app)
    first_sample = upload_dir / "first.mzML"
    first_sample.write_text("<mzML/>", encoding="utf-8")
    server._tasks["first"] = _task("first", "completed", first_sample)
    assert client.delete("/api/tasks/first").status_code == 200

    preview = client.post("/api/task-recycle-bin/purge/preview")
    assert preview.status_code == 200, preview.text
    stale_token = preview.json()["confirmation_token"]
    assert preview.json()["task_count"] == 1
    assert preview.json()["typed_confirmation_required"] is False

    # A new deletion changes the recycle state and invalidates the preview.
    second_sample = upload_dir / "second.mzML"
    second_sample.write_text("<mzML/>", encoding="utf-8")
    server._tasks["second"] = _task("second", "completed", second_sample)
    assert client.delete("/api/tasks/second").status_code == 200
    rejected = client.request(
        "DELETE",
        "/api/task-recycle-bin",
        json={"confirmation_token": stale_token},
    )
    assert rejected.status_code == 409
    assert list((tmp_path / "task_recycle_bin").glob("*/*/manifest.json"))

    fresh = client.post("/api/task-recycle-bin/purge/preview")
    assert fresh.status_code == 200, fresh.text
    removed = client.request(
        "DELETE",
        "/api/task-recycle-bin",
        json={"confirmation_token": fresh.json()["confirmation_token"]},
    )
    assert removed.status_code == 200, removed.text
    result = removed.json()
    assert result["removed_task_count"] == 2
    assert result["removed_size_bytes"] > 0
    assert result["failure_count"] == 0
    assert result["remaining_task_count"] == 0
    assert Path(result["receipt_path"]).is_file()
    assert not list((tmp_path / "task_recycle_bin").glob("*/*/manifest.json"))
    assert client.get("/api/task-recycle-bin").json()["total"] == 0


def test_literature_worker_cannot_recreate_artifact_after_task_deletion(
    monkeypatch, tmp_path,
):
    upload_dir, task_dir = _isolate_task_storage(monkeypatch, tmp_path)
    sample = upload_dir / "late-worker.mzML"
    sample.write_text("<mzML/>", encoding="utf-8")
    task_id = "late-worker"
    server._tasks[task_id] = _task(task_id, "completed", sample)

    def finish_after_deletion(*_args, **_kwargs):
        with server._tasks_lock:
            server._tasks.pop(task_id, None)
        return {
            "groups_processed": 1,
            "total_unique_ambiguous_groups": 1,
            "status_counts": {"ok": 1},
        }

    monkeypatch.setattr(
        server,
        "enrich_isomer_groups_with_pubmed",
        finish_after_deletion,
    )
    server._run_task_isomer_literature(task_id, results=[])

    assert not (task_dir / f"{task_id}_isomer_literature.json").exists()


def test_batch_delete_requires_fresh_state_bound_preview(monkeypatch, tmp_path):
    upload_dir, _ = _isolate_task_storage(monkeypatch, tmp_path)
    sample = upload_dir / "guarded.mzML"
    sample.write_text("<mzML/>", encoding="utf-8")
    server._tasks["guarded"] = _task("guarded", "completed", sample)
    client = TestClient(server.app)

    legacy_direct = client.post("/api/tasks/batch/delete", json=["guarded"])
    assert legacy_direct.status_code == 422

    preview = client.post(
        "/api/tasks/batch/delete/preview", json={"task_ids": ["guarded"]}
    )
    assert preview.status_code == 200, preview.text
    confirmation = preview.json()
    assert confirmation["selects_entire_registry"] is True
    assert confirmation["typed_confirmation_required"] is True

    server._tasks["guarded"]["state_version"] += 1
    rejected = client.post("/api/tasks/batch/delete", json={
        "task_ids": ["guarded"],
        "confirmation_token": confirmation["confirmation_token"],
        "confirmation_phrase": confirmation["confirmation_phrase"],
    })
    assert rejected.status_code == 409
    assert sample.is_file()
    assert "guarded" in server._tasks


def test_batch_delete_removes_queue_and_requests_running_stop(monkeypatch, tmp_path):
    upload_dir, _ = _isolate_task_storage(monkeypatch, tmp_path)
    for task_id, status in (
        ("queued-task", "queued"),
        ("completed-task", "completed"),
        ("running-task", "running"),
        ("finalizing-task", "finalizing"),
    ):
        sample = upload_dir / f"{task_id}.mzML"
        sample.write_text("<mzML/>", encoding="utf-8")
        server._tasks[task_id] = _task(task_id, status, sample)

    client = TestClient(server.app)
    task_ids = [
            "queued-task",
            "completed-task",
            "running-task",
            "finalizing-task",
    ]
    preview = client.post(
        "/api/tasks/batch/delete/preview", json={"task_ids": task_ids}
    )
    assert preview.status_code == 200, preview.text
    confirmation = preview.json()
    response = client.post("/api/tasks/batch/delete", json={
        "task_ids": task_ids,
        "confirmation_token": confirmation["confirmation_token"],
        "confirmation_phrase": confirmation["confirmation_phrase"],
    })

    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["deleted"] == ["queued-task", "completed-task"]
    assert payload["deleted_count"] == 2
    assert payload["stop_requested"] == ["running-task"]
    assert payload["stop_requested_count"] == 1
    assert payload["failed"][0]["id"] == "finalizing-task"

    assert "queued-task" not in server._tasks
    assert "completed-task" not in server._tasks
    assert not (upload_dir / "queued-task.mzML").exists()
    assert not (upload_dir / "completed-task.mzML").exists()
    recycle_root = tmp_path / "task_recycle_bin" / payload["deletion_id"]
    assert (recycle_root / "queued-task" / "manifest.json").is_file()
    assert (recycle_root / "completed-task" / "manifest.json").is_file()

    running = server._tasks["running-task"]
    assert running["status"] == "cancelling"
    assert running["state_history"][-1]["reason"] == "batch_delete_requested_safe_stop"
    assert server._task_cancel_events["running-task"].is_set()
    assert (upload_dir / "running-task.mzML").exists()

    assert server._tasks["finalizing-task"]["status"] == "finalizing"
    assert (upload_dir / "finalizing-task.mzML").exists()


def test_task_table_exposes_state_appropriate_actions():
    tasks_js = (
        Path(__file__).resolve().parents[1] / "web/static/js/tasks.js"
    ).read_text(encoding="utf-8")

    assert "Tasks.cancelTask('${t.id}')" in tasks_js
    assert "移出分析队列并放入回收站" in tasks_js
    assert "停止后可再次移除" in tasks_js
    assert "stop_requested_count" in tasks_js
    assert "['completed', 'failed', 'cancelled']" in tasks_js
    assert "/tasks/batch/delete/preview" in tasks_js
    assert "confirmation_token" in tasks_js
    assert "任务回收站" in tasks_js
    assert "/task-recycle-bin/purge/preview" in tasks_js
    assert "API.del('/task-recycle-bin'" in tasks_js
    assert "beginEmptyRecycleBin" in tasks_js
    assert "永久清空" in tasks_js
    single_delete = tasks_js.split("async function deleteTask", 1)[1].split(
        "async function cancelTask", 1
    )[0]
    batch_delete = tasks_js.split("async function batchDelete", 1)[1].split(
        "async function toggleRecycleBin", 1
    )[0]
    assert "confirm(" not in single_delete
    assert "prompt(" not in batch_delete
    assert "confirm(" not in batch_delete
