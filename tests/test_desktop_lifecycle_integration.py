"""Lightweight desktop/recycle integration, using in-memory synthetic tasks only."""
import asyncio
import json
from datetime import datetime, timezone
from types import SimpleNamespace

import pytest
from web import server
from runtime.storage_maintenance import StorageMaintenanceManager


def test_shutdown_snapshot_is_live_and_does_not_materialize_results(monkeypatch):
    monkeypatch.setattr(server, "_tasks", {
        "one":{"status":"running"}, "two":{"status":"queued"}, "three":{"status":"completed"}
    })
    monkeypatch.setattr(server, "_current_runtime_manifest_status", lambda: pytest.fail("heavy runtime probe"))
    value = asyncio.run(server.desktop_shutdown_status())
    assert value["task_counts"] == {"running":1,"queued":1,"completed":1}
    assert value["task_state_unknown"] is False
    assert value["project_root"] == str(server.PROJECT_ROOT)
    assert value["runtime_root"] == str(server.RUNTIME_ROOT)


def test_unrecognized_shutdown_state_is_explicit(monkeypatch):
    monkeypatch.setattr(server, "_tasks", {"one":{"status":"future-task-state"}})
    assert asyncio.run(server.desktop_shutdown_status())["task_state_unknown"] is True


@pytest.mark.parametrize("enabled", [False, True])
def test_existing_recycle_ui_gets_saved_retention_policy(monkeypatch, enabled):
    inventory = {"active_entries":[], "active_total":0,"restored_total":0,
                 "active_size_bytes":0,"active_file_count":0,"unsafe_active_entries":[]}
    monkeypatch.setattr(server, "_task_recycle_inventory", lambda: inventory)
    manager = SimpleNamespace(load_policy=lambda: {"recycle":{
        "automatic_cleanup_enabled":enabled,"retention_days":30,"last_enabled_at":None
    }})
    monkeypatch.setattr(server.app.state, "storage_maintenance_manager", manager)
    result = asyncio.run(server.list_task_recycle_bin(limit=100, offset=0))
    assert result["automatic_permanent_deletion"] is enabled
    assert result["retention_days"] == 30


def test_expired_subset_uses_real_server_inventory_and_preserves_other_files(monkeypatch, tmp_path):
    """Exercise real manifest/audit validation, not a mocked inventory schema."""
    for name, relative in [("OUTPUT_DIR","."),("UPLOAD_DIR","uploads"),("TASK_DIR","tasks")]:
        path = tmp_path / relative
        path.mkdir(exist_ok=True)
        monkeypatch.setattr(server, name, path)
    monkeypatch.setattr(server, "_tasks", {})
    protected = tmp_path / "uploads/keep.mzML"
    protected.write_bytes(b"synthetic retained original")
    roots = []
    for task_id, deleted in [("expired-task","2026-08-01T00:00:00+00:00"),("recent-task","2026-10-23T00:00:00+00:00")]:
        root = tmp_path / "task_recycle_bin/generation" / task_id
        (root / "files").mkdir(parents=True)
        recycled = root / "files/fixture.mzML"
        recycled.write_bytes(b"synthetic recycled payload")
        (root / "task.json").write_text(json.dumps({"id":task_id,"filename":task_id+".mzML","status":"completed"}))
        (root / "manifest.json").write_text(json.dumps({
            "schema_version":1,"deletion_id":"generation","task_id":task_id,
            "deleted_at":deleted,"recoverable":True,"failures":[],
            "moved":[{"original_path":str(tmp_path/"uploads"/(task_id+".mzML")),"recycle_path":str(recycled)}]
        }))
        server._append_task_deletion_audit({"event":"single_task_recycled","deletion_id":"generation","task_ids":[task_id]})
        roots.append(root)
    manager = StorageMaintenanceManager(
        output_dir=tmp_path, tasks_snapshot=server._storage_tasks_snapshot,
        recycle_inventory_unlocked=server._task_recycle_inventory_unlocked,
        recycle_lock=server._task_recycle_lock, cache_lock=server._report_export_lock,
        append_audit=server._append_task_deletion_audit,
        clock=lambda:datetime(2026,10,24,tzinfo=timezone.utc),
    )
    preview = manager.preview("recycle")
    assert not preview["blocked_reasons"] and preview["item_count"] == 1
    assert preview["items"][0]["task_id"] == "expired-task"
    result = manager.clean(preview["token"], confirmed=True)
    assert result["removed_count"] == 1 and result["failure_count"] == 0
    assert not roots[0].exists() and roots[1].is_dir()
    assert protected.read_bytes() == b"synthetic retained original"
    remaining = server._task_recycle_inventory_unlocked()
    assert remaining["active_total"] == 1 and not remaining["structural_errors"]
    assert remaining["active_entries"][0]["safe_to_purge"] is True
