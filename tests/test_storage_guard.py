from pathlib import Path

import pytest

from utils import storage_guard


def test_dynamic_reserve_uses_larger_of_floor_and_ratio(monkeypatch, tmp_path):
    usage = shutil_usage(total=1000 * storage_guard.GIB, used=900 * storage_guard.GIB)
    monkeypatch.setattr(storage_guard.shutil, "disk_usage", lambda _path: usage)
    monkeypatch.setenv("MSSYSTEM_MIN_FREE_GIB", "25")
    monkeypatch.setenv("MSSYSTEM_MIN_FREE_RATIO", "0.05")

    status = storage_guard.storage_status(tmp_path, 40 * storage_guard.GIB)

    assert status["reserve_bytes"] == 50 * storage_guard.GIB
    assert status["ok"] is True
    with pytest.raises(storage_guard.InsufficientStorageError):
        storage_guard.require_storage_capacity(
            tmp_path, 51 * storage_guard.GIB, "test write"
        )


def test_generated_backup_retention_is_narrow(tmp_path):
    database = tmp_path / "local_library.db"
    backups = []
    for index in range(5):
        path = tmp_path / f"local_library.db.prebuild_20260719T12000{index}.bak"
        path.write_bytes(str(index).encode())
        path.touch()
        backups.append(path)
    protected = tmp_path / "local_library.db.pre_literature_sync.bak"
    protected.write_bytes(b"milestone")

    removed = storage_guard.prune_generated_backups(database, keep=3)

    assert len(removed) == 2
    assert len(list(tmp_path.glob("local_library.db.prebuild_*.bak"))) == 3
    assert protected.exists()


def shutil_usage(*, total: int, used: int):
    return storage_guard.shutil._ntuple_diskusage(total, used, total - used)
