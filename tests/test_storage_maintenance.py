from __future__ import annotations

import asyncio
import hashlib
import json
import os
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from runtime.storage_maintenance import (
    StorageMaintenanceError,
    StorageMaintenanceManager,
)
from web.storage_management import install_storage_management
from web.storage_management import _automatic_worker


class Clock:
    def __init__(self, value: datetime):
        self.value = value

    def __call__(self) -> datetime:
        return self.value

    def advance(self, **kwargs) -> None:
        self.value += timedelta(**kwargs)


class RecycleFixture:
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.root = output_dir / "task_recycle_bin"
        self.rows: dict[tuple[str, str], dict] = {}
        self.calls = 0

    def add(
        self,
        deletion_id: str,
        task_id: str,
        deleted_at: datetime,
        *,
        restored: bool = False,
        safe: bool = True,
    ) -> Path:
        task_root = self.root / deletion_id / task_id
        task_root.mkdir(parents=True)
        (task_root / "manifest.json").write_text(
            json.dumps({"deleted_at": deleted_at.isoformat()}), encoding="utf-8"
        )
        (task_root / "task.json").write_text(
            json.dumps({"id": task_id, "filename": f"{task_id}.mzML"}),
            encoding="utf-8",
        )
        data = task_root / "files" / "payload.bin"
        data.parent.mkdir()
        data.write_bytes(task_id.encode("utf-8") * 8)
        stamp = deleted_at.timestamp()
        for path in sorted(task_root.rglob("*"), reverse=True):
            os.utime(path, (stamp, stamp), follow_symlinks=False)
        os.utime(task_root, (stamp, stamp), follow_symlinks=False)
        self.rows[(deletion_id, task_id)] = {
            "deleted_at": deleted_at.isoformat(),
            "restored": restored,
            "safe": safe,
        }
        return task_root

    @staticmethod
    def _stats(root: Path) -> tuple[int, int, str, bool]:
        size = 0
        files = 0
        linked = False
        rows = []
        stack = [root]
        while stack:
            path = stack.pop()
            stat = path.lstat()
            relative = str(path.relative_to(root)) or "."
            if path.is_symlink():
                linked = True
                kind = "link"
            elif path.is_file():
                kind = "file"
                size += stat.st_size
                files += 1
            else:
                kind = "dir"
                stack.extend(path.iterdir())
            rows.append((relative, kind, stat.st_size, stat.st_mtime_ns))
        signature = hashlib.sha256(
            json.dumps(sorted(rows), separators=(",", ":")).encode()
        ).hexdigest()
        return size, files, signature, linked

    def inventory(self) -> dict:
        self.calls += 1
        active = []
        restored_entries = []
        unsafe = []
        for (deletion_id, task_id), metadata in sorted(self.rows.items()):
            task_root = self.root / deletion_id / task_id
            if not task_root.exists() and not task_root.is_symlink():
                continue
            size, files, signature, linked = self._stats(task_root)
            safe = bool(metadata["safe"] and not linked)
            entry = {
                "deletion_id": deletion_id,
                "task_id": task_id,
                "filename": f"{task_id}.mzML",
                "deleted_at": metadata["deleted_at"],
                "restored": metadata["restored"],
                "size_bytes": size,
                "file_count": files,
                "tree_signature": signature,
                "manifest_sha256": hashlib.sha256(
                    (task_root / "manifest.json").read_bytes()
                ).hexdigest(),
                "safe_to_purge": safe,
                "purge_block_reasons": [] if safe else ["unsafe_fixture"],
            }
            if metadata["restored"]:
                restored_entries.append(entry)
            else:
                active.append(entry)
                if not safe:
                    unsafe.append(entry)
        signature = hashlib.sha256(
            json.dumps(
                [
                    (item["deletion_id"], item["task_id"], item["tree_signature"])
                    for item in active
                ],
                separators=(",", ":"),
            ).encode()
        ).hexdigest()
        return {
            "active_entries": active,
            "restored_entries": restored_entries,
            "unsafe_active_entries": unsafe,
            "structural_errors": [],
            "signature": signature,
        }


def _manager(tmp_path: Path, *, now: datetime | None = None):
    clock = Clock(now or datetime(2026, 10, 24, tzinfo=timezone.utc))
    recycle = RecycleFixture(tmp_path)
    tasks: dict[str, dict] = {}
    lock = threading.RLock()
    cache_lock = threading.RLock()
    audits = []
    manager = StorageMaintenanceManager(
        output_dir=tmp_path,
        tasks_snapshot=lambda: dict(tasks),
        recycle_inventory_unlocked=recycle.inventory,
        recycle_lock=lock,
        cache_lock=cache_lock,
        append_audit=audits.append,
        clock=clock,
        open_file_detector=lambda: set(),
    )
    return manager, recycle, tasks, clock, cache_lock, audits


def _old_cache(root: Path, name: str, when: datetime) -> Path:
    path = root / "export_tmp" / name
    path.mkdir(parents=True, exist_ok=True)
    child = path / "export.part"
    child.write_bytes(b"temporary-export")
    stamp = when.timestamp()
    for item in path.rglob("*"):
        os.utime(item, (stamp, stamp), follow_symlinks=False)
    os.utime(path, (stamp, stamp))
    return path


def test_defaults_are_readable_disabled_and_do_not_scan_automatically(tmp_path):
    manager, recycle, *_ = _manager(tmp_path)

    policy = manager.load_policy()
    automatic = manager.run_auto_cleanup()

    assert policy["recycle"] == {
        "automatic_cleanup_enabled": False,
        "retention_days": 30,
        "last_enabled_at": None,
    }
    assert policy["cache"]["automatic_cleanup_enabled"] is False
    assert policy["cache"]["retention_days"] == 7
    assert automatic == {"status": "disabled", "results": []}
    assert recycle.calls == 0


def test_first_enable_grants_full_retention_period_to_existing_recycle_items(tmp_path):
    manager, recycle, _, clock, *_ = _manager(tmp_path)
    recycle.add("delete-old", "old-task", clock() - timedelta(days=100))

    before = manager.preview("recycle")
    assert before["item_count"] == 1
    manager.save_policy(
        recycle_enabled=True,
        recycle_retention_days=30,
        cache_enabled=False,
        cache_retention_days=7,
        confirm_auto_delete=True,
    )

    grace = manager.preview("recycle")
    assert grace["item_count"] == 0
    clock.advance(days=30)
    expired = manager.preview("recycle")
    assert expired["item_count"] == 1


def test_enabling_auto_delete_requires_explicit_confirmation(tmp_path):
    manager, *_ = _manager(tmp_path)

    with pytest.raises(StorageMaintenanceError, match="显式确认"):
        manager.save_policy(
            recycle_enabled=True,
            recycle_retention_days=30,
            cache_enabled=False,
            cache_retention_days=7,
            confirm_auto_delete=False,
        )

    assert not manager.policy_path.exists()


def test_shortening_enabled_retention_requires_explicit_confirmation(tmp_path):
    manager, *_ = _manager(tmp_path)
    manager.save_policy(
        recycle_enabled=True,
        recycle_retention_days=90,
        cache_enabled=True,
        cache_retention_days=30,
        confirm_auto_delete=True,
    )

    with pytest.raises(StorageMaintenanceError, match="显式确认"):
        manager.save_policy(
            recycle_enabled=True,
            recycle_retention_days=30,
            cache_enabled=True,
            cache_retention_days=30,
            confirm_auto_delete=False,
        )
    assert manager.load_policy()["recycle"]["retention_days"] == 90

    updated = manager.save_policy(
        recycle_enabled=True,
        recycle_retention_days=30,
        cache_enabled=True,
        cache_retention_days=7,
        confirm_auto_delete=True,
    )
    assert updated["recycle"]["retention_days"] == 30
    assert updated["cache"]["retention_days"] == 7


def test_persisted_policy_rejects_unknown_version_types_and_retention(tmp_path):
    manager, *_ = _manager(tmp_path)
    manager.policy_path.parent.mkdir(parents=True)
    invalid_payloads = [
        {
            "schema_version": 2,
            "recycle": {},
            "cache": {},
        },
        {
            "schema_version": 1,
            "recycle": {
                "automatic_cleanup_enabled": "false",
                "retention_days": 30,
                "last_enabled_at": None,
            },
            "cache": {
                "automatic_cleanup_enabled": False,
                "retention_days": 7,
                "last_enabled_at": None,
            },
            "last_cleanup": None,
            "last_auto_cleanup_date": None,
        },
        {
            "schema_version": 1,
            "recycle": {
                "automatic_cleanup_enabled": False,
                "retention_days": 31,
                "last_enabled_at": None,
            },
            "cache": {
                "automatic_cleanup_enabled": False,
                "retention_days": 7,
                "last_enabled_at": None,
            },
            "last_cleanup": None,
            "last_auto_cleanup_date": None,
        },
    ]
    for payload in invalid_payloads:
        manager.policy_path.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(StorageMaintenanceError, match="策略"):
            manager.load_policy()


def test_first_cache_enable_also_grants_a_full_retention_period(tmp_path):
    manager, _, _, clock, *_ = _manager(tmp_path)
    old = _old_cache(tmp_path, "old-export", clock() - timedelta(days=100))
    assert manager.preview("cache")["item_count"] == 1

    manager.save_policy(
        recycle_enabled=False,
        recycle_retention_days=30,
        cache_enabled=True,
        cache_retention_days=7,
        confirm_auto_delete=True,
    )

    assert manager.preview("cache")["item_count"] == 0
    assert old.exists()
    clock.advance(days=7)
    assert manager.preview("cache")["item_count"] == 1


def test_recycle_cleanup_removes_only_expired_safe_subset_and_preserves_data(tmp_path):
    manager, recycle, _, clock, _, audits = _manager(tmp_path)
    old = recycle.add("delete-old", "old-task", clock() - timedelta(days=40))
    recent = recycle.add("delete-new", "new-task", clock() - timedelta(days=5))
    restored = recycle.add(
        "delete-restored", "restored-task", clock() - timedelta(days=90),
        restored=True,
    )
    protected = [
        tmp_path / "tmp_uploads" / "raw.mzML",
        tmp_path / "tasks" / "task_results.json",
        tmp_path / "cache" / "active.bin",
        tmp_path / "models" / "model.bin",
        tmp_path / "logs" / "service.log",
    ]
    for path in protected:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"protected")

    preview = manager.preview("recycle")
    result = manager.clean(preview["token"], confirmed=True)

    assert preview["item_count"] == 1
    assert preview["items"][0]["task_id"] == "old-task"
    assert result["removed_count"] == 1
    assert not old.exists()
    assert recent.exists()
    assert restored.exists()
    assert all(path.read_bytes() == b"protected" for path in protected)
    assert Path(result["receipt_path"]).is_file()
    assert audits[-1]["event"] == "storage_maintenance_cleanup"


def test_recycle_preview_is_rejected_after_tree_tampering_and_token_is_consumed(tmp_path):
    manager, recycle, _, clock, *_ = _manager(tmp_path)
    task_root = recycle.add("delete-old", "old-task", clock() - timedelta(days=40))
    preview = manager.preview("recycle")
    (task_root / "files" / "payload.bin").write_bytes(b"tampered")

    with pytest.raises(StorageMaintenanceError, match="变化"):
        manager.clean(preview["token"], confirmed=True)
    with pytest.raises(StorageMaintenanceError, match="失效"):
        manager.clean(preview["token"], confirmed=True)
    assert task_root.exists()


def test_active_or_restored_task_id_cannot_be_purged(tmp_path):
    manager, recycle, tasks, clock, *_ = _manager(tmp_path)
    task_root = recycle.add("delete-old", "same-task", clock() - timedelta(days=40))
    tasks["same-task"] = {"id": "same-task", "status": "completed"}

    preview = manager.preview("recycle")

    assert preview["token"] is None
    assert preview["item_count"] == 0
    assert any("task_id_exists" in reason for reason in preview["blocked_reasons"])
    assert task_root.exists()


def test_recycle_only_purges_safe_entries_absent_from_registry(tmp_path):
    manager, recycle, tasks, clock, *_ = _manager(tmp_path)
    removable = recycle.add(
        "delete-removable", "removable-task", clock() - timedelta(days=40)
    )
    registered = recycle.add(
        "delete-registered", "registered-task", clock() - timedelta(days=5)
    )
    tasks["registered-task"] = {
        "id": "registered-task", "status": "future-launcher-state"
    }

    preview = manager.preview("recycle")
    result = manager.clean(preview["token"], confirmed=True)

    assert [item["task_id"] for item in preview["items"]] == ["removable-task"]
    assert result["removed_count"] == 1
    assert not removable.exists()
    assert registered.exists()


def test_recycle_generation_symlink_cannot_escape_allowlisted_root(tmp_path):
    now = datetime(2026, 10, 24, tzinfo=timezone.utc)
    outside_task = tmp_path / "outside" / "same-task"
    outside_task.mkdir(parents=True)
    (outside_task / "manifest.json").write_text("{}", encoding="utf-8")
    (outside_task / "task.json").write_text("{}", encoding="utf-8")
    recycle_root = tmp_path / "task_recycle_bin"
    recycle_root.mkdir()
    (recycle_root / "delete-link").symlink_to(outside_task.parent, target_is_directory=True)
    inventory = {
        "active_entries": [{
            "deletion_id": "delete-link",
            "task_id": "same-task",
            "filename": "same-task.mzML",
            "deleted_at": (now - timedelta(days=40)).isoformat(),
            "size_bytes": 4,
            "file_count": 2,
            "tree_signature": "trusted-callback-signature",
            "manifest_sha256": "trusted-callback-manifest",
            "safe_to_purge": True,
            "purge_block_reasons": [],
        }],
        "restored_entries": [],
        "unsafe_active_entries": [],
        "structural_errors": [],
        "signature": "trusted-callback-inventory",
    }
    manager = StorageMaintenanceManager(
        output_dir=tmp_path,
        tasks_snapshot=lambda: {},
        recycle_inventory_unlocked=lambda: inventory,
        recycle_lock=threading.RLock(),
        cache_lock=threading.RLock(),
        clock=lambda: now,
        open_file_detector=lambda: set(),
    )

    preview = manager.preview("recycle")

    assert preview["token"] is None
    assert preview["item_count"] == 0
    assert any("recycle_path_unavailable" in item for item in preview["blocked_reasons"])
    assert (outside_task / "manifest.json").exists()


def test_cache_cleanup_is_allowlisted_and_preserves_active_cache_and_results(tmp_path):
    manager, _, _, clock, *_ = _manager(tmp_path)
    old = _old_cache(tmp_path, "old-export", clock() - timedelta(days=10))
    recent = _old_cache(tmp_path, "new-export", clock() - timedelta(days=1))
    protected = [
        tmp_path / "cache" / "active.cache",
        tmp_path / "tasks" / "result.json",
        tmp_path / "tmp_uploads" / "raw.mzML",
    ]
    for path in protected:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"protected")

    preview = manager.preview("cache")
    result = manager.clean(preview["token"], confirmed=True)

    assert preview["item_count"] == 1
    assert preview["items"][0]["name"] == "old-export"
    assert result["removed_count"] == 1
    assert not old.exists()
    assert recent.exists()
    assert all(path.exists() for path in protected)


def test_cache_cleanup_refuses_active_tasks_busy_export_and_unknown_open_files(tmp_path):
    manager, _, tasks, clock, cache_lock, _ = _manager(tmp_path)
    old = _old_cache(tmp_path, "old-export", clock() - timedelta(days=10))
    tasks["running"] = {"id": "running", "status": "running"}
    active = manager.preview("cache")
    assert active["token"] is None
    assert "active_tasks_present" in active["blocked_reasons"]

    tasks.clear()
    cache_lock.acquire()
    try:
        busy = manager.preview("cache")
    finally:
        cache_lock.release()
    # RLock is reentrant in the same thread; use an open-file guard as the
    # deterministic in-process export safety assertion.
    manager.open_file_detector = lambda: {old / "export.part"}
    opened = manager.preview("cache")
    assert opened["token"] is None
    assert "export_tmp_contains_open_file" in opened["blocked_reasons"]
    manager.open_file_detector = lambda: None
    unknown = manager.preview("cache")
    assert unknown["token"] is None
    assert "open_file_detection_unavailable" in unknown["blocked_reasons"]
    assert busy["item_count"] == 1


@pytest.mark.parametrize(
    ("status", "bucket"),
    [
        ("converting", "converting"),
        ("stopping", "stopping"),
        ("cancel_requested", "cancel_requested"),
        ("future-launcher-state", "unknown:future-launcher-state"),
        (None, "unknown:missing"),
    ],
)
def test_nonterminal_and_unknown_task_states_block_cache_and_auto_cleanup(
    tmp_path, status, bucket
):
    manager, _, tasks, clock, *_ = _manager(tmp_path)
    _old_cache(tmp_path, "old-export", clock() - timedelta(days=10))
    tasks["task"] = {"id": "task", "status": status}
    manager.save_policy(
        recycle_enabled=False,
        recycle_retention_days=30,
        cache_enabled=True,
        cache_retention_days=7,
        confirm_auto_delete=True,
    )

    active = manager.active_tasks()
    preview = manager.preview("cache")
    automatic = manager.run_auto_cleanup()

    assert active == {"count": 1, "by_status": {bucket: 1}}
    assert preview["token"] is None
    assert "active_tasks_present" in preview["blocked_reasons"]
    assert automatic["status"] == "skipped_active_tasks"


def test_terminal_task_states_allow_cache_preview(tmp_path):
    manager, _, tasks, clock, *_ = _manager(tmp_path)
    _old_cache(tmp_path, "old-export", clock() - timedelta(days=10))
    tasks.update({
        status: {"id": status, "status": status}
        for status in ("completed", "failed", "cancelled")
    })

    assert manager.active_tasks() == {"count": 0, "by_status": {}}
    assert manager.preview("cache")["token"] is not None


def test_cache_lock_held_by_another_thread_blocks_preview(tmp_path):
    manager, _, _, clock, cache_lock, _ = _manager(tmp_path)
    _old_cache(tmp_path, "old-export", clock() - timedelta(days=10))
    entered = threading.Event()
    release = threading.Event()

    def holder():
        with cache_lock:
            entered.set()
            release.wait(5)

    thread = threading.Thread(target=holder)
    thread.start()
    assert entered.wait(2)
    try:
        preview = manager.preview("cache")
    finally:
        release.set()
        thread.join(2)

    assert preview["token"] is None
    assert "cache_lock_unavailable_or_busy" in preview["blocked_reasons"]


def test_cache_clean_rechecks_lock_and_open_files_after_preview(tmp_path):
    manager, _, _, clock, cache_lock, _ = _manager(tmp_path)
    old = _old_cache(tmp_path, "old-export", clock() - timedelta(days=10))
    preview = manager.preview("cache")
    entered = threading.Event()
    release = threading.Event()

    def holder():
        with cache_lock:
            entered.set()
            release.wait(5)

    thread = threading.Thread(target=holder)
    thread.start()
    assert entered.wait(2)
    try:
        with pytest.raises(StorageMaintenanceError, match="导出锁"):
            manager.clean(preview["token"], confirmed=True)
    finally:
        release.set()
        thread.join(2)
    assert old.exists()

    preview = manager.preview("cache")
    manager.open_file_detector = lambda: {old / "export.part"}
    with pytest.raises(StorageMaintenanceError, match="占用状态"):
        manager.clean(preview["token"], confirmed=True)
    assert old.exists()


def test_symlink_and_path_escape_are_never_cache_cleanup_targets(tmp_path):
    manager, _, _, clock, *_ = _manager(tmp_path)
    outside = tmp_path / "outside"
    outside.mkdir()
    secret = outside / "raw.mzML"
    secret.write_bytes(b"raw")
    export_tmp = tmp_path / "export_tmp"
    export_tmp.mkdir()
    (export_tmp / "escape").symlink_to(outside, target_is_directory=True)
    old_stamp = (clock() - timedelta(days=10)).timestamp()
    os.utime(export_tmp / "escape", (old_stamp, old_stamp), follow_symlinks=False)

    preview = manager.preview("cache")

    assert preview["token"] is None
    assert any("symbolic_link" in reason for reason in preview["blocked_reasons"])
    assert secret.read_bytes() == b"raw"


def test_cache_preview_detects_inode_or_content_tampering(tmp_path):
    manager, _, _, clock, *_ = _manager(tmp_path)
    old = _old_cache(tmp_path, "old-export", clock() - timedelta(days=10))
    preview = manager.preview("cache")
    (old / "export.part").write_bytes(b"changed-after-preview")

    with pytest.raises(StorageMaintenanceError, match="变化"):
        manager.clean(preview["token"], confirmed=True)
    assert old.exists()


def test_api_contract_policy_validation_single_use_token_and_app_state(tmp_path):
    recycle = RecycleFixture(tmp_path)
    now = datetime.now(timezone.utc)
    recycle.add("delete-old", "old-task", now - timedelta(days=40))
    tasks = {}
    app = FastAPI()
    manager = install_storage_management(
        app,
        output_dir=tmp_path,
        tasks_snapshot=lambda: dict(tasks),
        recycle_inventory_unlocked=recycle.inventory,
        recycle_lock=threading.RLock(),
        cache_lock=threading.RLock(),
    )
    manager.open_file_detector = lambda: set()

    with TestClient(app) as client:
        overview = client.get("/api/storage-management")
        assert overview.status_code == 200
        assert overview.json()["policy"]["recycle"]["retention_days"] == 30
        invalid_path = client.put("/api/storage-management/policy", json={
            "recycle_automatic_cleanup_enabled": False,
            "recycle_retention_days": 30,
            "cache_automatic_cleanup_enabled": False,
            "cache_retention_days": 7,
            "confirm_auto_delete": False,
            "output_dir": "/tmp/escape",
        })
        assert invalid_path.status_code == 422
        preview = client.post(
            "/api/storage-management/preview", json={"target": "recycle"}
        )
        assert preview.status_code == 200
        token = preview.json()["token"]
        cleaned = client.post("/api/storage-management/clean", json={
            "token": token, "confirmed": True,
        })
        assert cleaned.status_code == 200
        replay = client.post("/api/storage-management/clean", json={
            "token": token, "confirmed": True,
        })
        assert replay.status_code == 409

    assert app.state.storage_maintenance_manager is manager


def test_history_api_returns_records_sorted_limited_and_without_paths(tmp_path):
    runtime = tmp_path / "runtime"
    runtime.mkdir()
    receipts = [
        (
            "storage_maintenance_recycle_20261023T010203Z-abcdef12.json",
            {
                "schema_version": 1,
                "target": "recycle",
                "finished_at": "2026-10-23T01:02:03+00:00",
                "requested_items": [{
                    "deletion_id": "delete-1", "task_id": "task-1",
                    "filename": "sample-one.mzML",
                }],
                "removed": [{
                    "deletion_id": "delete-1", "task_id": "task-1", "size_bytes": 12,
                }],
                "failures": [],
            },
        ),
        (
            "storage_maintenance_cache_20261024T010203Z-abcdef34.json",
            {
                "schema_version": 1,
                "target": "cache",
                "finished_at": "2026-10-24T01:02:03Z",
                "requested_items": [{"name": "export-1"}],
                "removed": [],
                "failures": [{"name": "export-1", "error": "file busy"}],
            },
        ),
    ]
    for name, payload in receipts:
        (runtime / name).write_text(json.dumps(payload), encoding="utf-8")
    app = FastAPI()
    install_storage_management(
        app,
        output_dir=tmp_path,
        tasks_snapshot=lambda: {},
        recycle_inventory_unlocked=lambda: {},
        recycle_lock=threading.RLock(),
        cache_lock=threading.RLock(),
    )

    with TestClient(app) as client:
        response = client.get("/api/storage-management/history", params={"limit": 1})
        too_small = client.get("/api/storage-management/history", params={"limit": 0})
        too_large = client.get("/api/storage-management/history", params={"limit": 21})

    assert response.status_code == 200
    assert too_small.status_code == 422
    assert too_large.status_code == 422
    body = response.json()
    assert list(body) == ["records", "rejected_count"]
    assert body["records"] == [{
        "target": "cache",
        "finished_at": "2026-10-24T01:02:03+00:00",
        "removed_count": 0,
        "removed_bytes": 0,
        "failure_count": 1,
        "removed_items": [],
        "failed_items": [{"name": "export-1", "error": "file busy"}],
    }]
    assert str(tmp_path) not in json.dumps(body)


def test_history_without_receipts_returns_empty_records(tmp_path):
    manager, *_ = _manager(tmp_path)

    assert manager.cleanup_history() == {"records": [], "rejected_count": 0}


def test_history_rejects_symlink_oversize_malformed_and_target_mismatch(tmp_path):
    manager, *_ = _manager(tmp_path)
    manager.runtime_dir.mkdir(parents=True)
    outside = tmp_path / "outside.json"
    outside.write_text("{}", encoding="utf-8")
    (manager.runtime_dir / "storage_maintenance_cache_20261020T000000Z-aaaa0001.json").symlink_to(outside)
    (manager.runtime_dir / "storage_maintenance_cache_20261020T000001Z-aaaa0002.json").write_bytes(
        b"x" * (1024 * 1024 + 1)
    )
    (manager.runtime_dir / "storage_maintenance_cache_20261020T000002Z-aaaa0003.json").write_text(
        "not-json", encoding="utf-8"
    )
    (manager.runtime_dir / "storage_maintenance_cache_20261020T000003Z-aaaa0004.json").write_text(
        json.dumps({
            "schema_version": 1,
            "target": "recycle",
            "finished_at": "2026-10-20T00:00:03Z",
        }),
        encoding="utf-8",
    )
    (manager.runtime_dir / "unrelated.json").write_text(
        json.dumps({"absolute_path": "/private/data"}), encoding="utf-8"
    )

    history = manager.cleanup_history()

    assert history == {"records": [], "rejected_count": 4}
    assert outside.read_text(encoding="utf-8") == "{}"


def test_startup_worker_does_not_immediately_delete_historical_items(tmp_path):
    recycle = RecycleFixture(tmp_path)
    old = recycle.add(
        "delete-old", "old-task",
        datetime.now(timezone.utc) - timedelta(days=100),
    )
    app = FastAPI()
    manager = install_storage_management(
        app,
        output_dir=tmp_path,
        tasks_snapshot=lambda: {},
        recycle_inventory_unlocked=recycle.inventory,
        recycle_lock=threading.RLock(),
        cache_lock=threading.RLock(),
    )
    manager.save_policy(
        recycle_enabled=True,
        recycle_retention_days=7,
        cache_enabled=False,
        cache_retention_days=7,
        confirm_auto_delete=True,
    )

    with TestClient(app):
        assert old.exists()
        assert manager.load_policy()["last_cleanup"] is None


def test_automatic_worker_stop_waits_for_inflight_receipt_finalization(tmp_path):
    receipt = tmp_path / "automatic-receipt.json"
    started = threading.Event()
    release = threading.Event()

    class FakeManager:
        def run_auto_cleanup(self):
            started.set()
            assert release.wait(2)
            receipt.write_text('{"finished": true}', encoding="utf-8")

    async def scenario():
        stop_event = asyncio.Event()
        worker = asyncio.create_task(
            _automatic_worker(FakeManager(), stop_event, interval_seconds=0.001)
        )
        assert await asyncio.to_thread(started.wait, 1)
        stop_event.set()
        await asyncio.sleep(0)
        assert not worker.done()
        release.set()
        await asyncio.wait_for(worker, timeout=1)

    asyncio.run(scenario())
    assert json.loads(receipt.read_text(encoding="utf-8")) == {"finished": True}


def test_automatic_worker_logs_cleanup_exception(caplog):
    async def scenario():
        loop = asyncio.get_running_loop()
        stop_event = asyncio.Event()

        class FailingManager:
            def run_auto_cleanup(self):
                loop.call_soon_threadsafe(stop_event.set)
                raise RuntimeError("synthetic cleanup failure")

        await _automatic_worker(
            FailingManager(), stop_event, interval_seconds=0.001
        )

    caplog.set_level("ERROR", logger="web.storage_management")
    asyncio.run(scenario())
    assert "automatic storage maintenance failed" in caplog.text
    assert "synthetic cleanup failure" in caplog.text
