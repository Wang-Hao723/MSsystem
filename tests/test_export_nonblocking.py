from __future__ import annotations

import inspect
import threading

from fastapi.testclient import TestClient

from web import server


def test_large_excel_export_does_not_block_health_event_loop(monkeypatch):
    # This test isolates report scheduling. Real cold library verification can
    # take minutes regardless of an export; readiness integrity has separate
    # integration checks and must not determine this five-second event gate.
    monkeypatch.setattr(server, "_current_runtime_manifest_status",
                        lambda: {"release_ready": True})
    monkeypatch.setattr(server, "library_metadata_repair_status",
                        lambda *_: {"status": "verified"})
    monkeypatch.setattr(server, "sql_metadata_repair_status",
                        lambda: {"status": "verified"})
    monkeypatch.setattr(server, "source_compound_mass_status",
                        lambda *_: {"status": "verified"})
    task_id = "slow-export-task"
    started = threading.Event()
    release = threading.Event()
    task = {
        "id": task_id,
        "filename": "sample.mzML",
        "status": "completed",
        "progress": 100,
        "upload_time": "2026-08-03T00:00:00+00:00",
    }

    def slow_workbook(*_args, **_kwargs):
        started.set()
        assert release.wait(timeout=5)
        return b"xlsx"

    monkeypatch.setattr(server, "_tasks", {task_id: task})
    monkeypatch.setattr(server, "_build_excel_export_subprocess", slow_workbook)
    export_response = {}

    def request_export():
        export_response["response"] = TestClient(server.app).get(
            f"/api/tasks/{task_id}/download/excel"
        )

    worker = threading.Thread(target=request_export)
    worker.start()
    assert started.wait(timeout=5)
    health = TestClient(server.app).get("/api/health")
    release.set()
    worker.join(timeout=5)

    assert health.status_code == 200
    assert health.json()["status"] in {"ok", "degraded"}
    assert export_response["response"].status_code == 200
    assert export_response["response"].content == b"xlsx"


def test_task_metadata_and_status_never_lazy_load_disk_results(
    monkeypatch, tmp_path
):
    task_id = "disk-heavy-task"
    task = {
        "id": task_id,
        "filename": "large_sample.mzML",
        "status": "completed",
        "progress": 100,
        "message": "complete",
        "upload_time": "2026-08-23T00:00:00+00:00",
        "results": [],
        "presentation_calibration_version": "stale-version",
    }
    (tmp_path / f"{task_id}_results.json").write_text(
        '[{"large":"artifact"}]', encoding="utf-8"
    )

    def fail_if_loaded(_task_id):
        raise AssertionError("metadata polling must not load result artifacts")

    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(server, "_tasks", {task_id: task})
    monkeypatch.setattr(server, "_get_task_results", fail_if_loaded)

    client = TestClient(server.app)
    metadata = client.get(f"/api/tasks/{task_id}")
    status = client.get(f"/api/tasks/{task_id}/status")

    assert metadata.status_code == 200
    assert status.status_code == 200
    assert "results" not in metadata.json()
    assert task["results"] == []


def test_presentation_refresh_discards_persisted_rows_from_registry(
    monkeypatch, tmp_path
):
    task_id = "refresh-heavy-task"
    task = {
        "id": task_id,
        "filename": "sample.mzML",
        "status": "completed",
        "progress": 100,
        "message": "old",
        "results": [{"spectrum_id": "scan-1", "match": None}],
        "presentation_calibration_version": "stale-version",
    }
    (tmp_path / f"{task_id}_results.json").write_text(
        "[]", encoding="utf-8"
    )

    def refresh(view, _rows):
        view["presentation_calibration_version"] = (
            server.PRESENTATION_CALIBRATION_VERSION
        )
        view["message"] = "refreshed"

    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(server, "_tasks", {task_id: task})
    monkeypatch.setattr(server, "_project_claims_for_task", lambda _id, rows: rows)
    monkeypatch.setattr(server, "_refresh_completed_task_presentation", refresh)
    monkeypatch.setattr(server, "_save_tasks", lambda: True)

    refreshed = server._ensure_task_presentation_current(task_id)

    assert refreshed["message"] == "refreshed"
    assert task["results"] == []
    assert (
        task["presentation_calibration_version"]
        == server.PRESENTATION_CALIBRATION_VERSION
    )


def test_batch_zip_is_disk_spooled_and_never_copied_with_getvalue():
    source = inspect.getsource(server.batch_export)
    assert "tempfile.SpooledTemporaryFile" in source
    assert "max_size=8 * 1024 * 1024" in source
    assert "iter([buf.getvalue()])" not in source
    assert "stream_archive()" in source
