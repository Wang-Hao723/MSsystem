from __future__ import annotations

import threading
import time

from web import server


def _isolate_runtime_status(monkeypatch, tmp_path) -> None:
    manifest_path = tmp_path / "config" / "runtime_manifest.json"
    manifest_path.parent.mkdir()
    manifest_path.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(server, "RUNTIME_SNAPSHOT_ID", "startup-snapshot")
    monkeypatch.setattr(server, "RUNTIME_SOURCE_TREE_SHA256", "source-hash")
    monkeypatch.setattr(server, "_runtime_snapshot_lock", threading.RLock())
    monkeypatch.setattr(
        server, "_runtime_manifest_status_check_lock", threading.Lock()
    )
    monkeypatch.setattr(
        server,
        "_runtime_snapshot_cache",
        {"signature": None, "value": "startup-snapshot"},
    )
    monkeypatch.setattr(
        server,
        "_runtime_manifest_status_cache",
        {"signature": None, "value": None},
    )
    monkeypatch.setattr(
        server, "_runtime_manifest_dependency_signature", lambda _root: ("stable",)
    )
    monkeypatch.setattr(
        server, "_runtime_snapshot_id_from_manifest", lambda _root: "new-snapshot"
    )
    monkeypatch.setattr(
        server, "_runtime_manifest_source_tree_sha256", lambda _root: "source-hash"
    )


def test_slow_health_check_does_not_block_snapshot_and_is_single_flight(
    monkeypatch, tmp_path
):
    _isolate_runtime_status(monkeypatch, tmp_path)
    freshness_started = threading.Event()
    release_freshness = threading.Event()
    calls_lock = threading.Lock()
    calls = 0

    def slow_freshness(_root, _manifest_path):
        nonlocal calls
        with calls_lock:
            calls += 1
        freshness_started.set()
        assert release_freshness.wait(timeout=5)
        return {
            "status": "ready",
            "fresh": True,
            "release_ready": True,
            "manifest_source_tree_sha256": "source-hash",
            "current_source_tree_sha256": "source-hash",
        }

    monkeypatch.setattr(server, "_runtime_manifest_freshness", slow_freshness)
    statuses: list[dict] = []
    first_health = threading.Thread(
        target=lambda: statuses.append(server._current_runtime_manifest_status())
    )
    first_health.start()
    assert freshness_started.wait(timeout=2)

    second_health = threading.Thread(
        target=lambda: statuses.append(server._current_runtime_manifest_status())
    )
    second_health.start()
    snapshot: list[str | None] = []
    snapshot_done = threading.Event()

    def read_snapshot() -> None:
        snapshot.append(server._current_runtime_snapshot_id())
        snapshot_done.set()

    snapshot_reader = threading.Thread(target=read_snapshot)
    snapshot_reader.start()
    assert snapshot_done.wait(timeout=0.5), "health preflight blocked snapshot identity"
    time.sleep(0.05)
    with calls_lock:
        assert calls == 1

    release_freshness.set()
    first_health.join(timeout=2)
    second_health.join(timeout=2)
    snapshot_reader.join(timeout=2)
    assert not first_health.is_alive()
    assert not second_health.is_alive()
    assert snapshot == ["new-snapshot"]
    assert len(statuses) == 2
    assert all(status["release_ready"] is True for status in statuses)
    with calls_lock:
        assert calls == 1


def test_failed_health_check_never_reports_release_ready(monkeypatch, tmp_path):
    _isolate_runtime_status(monkeypatch, tmp_path)

    monkeypatch.setattr(
        server,
        "_runtime_manifest_freshness",
        lambda _root, _manifest_path: {
            "status": "invalid",
            "fresh": False,
            "release_ready": False,
            "manifest_source_tree_sha256": "source-hash",
            "current_source_tree_sha256": "source-hash",
        },
    )

    status = server._current_runtime_manifest_status()

    assert status["status"] == "invalid"
    assert status["fresh"] is False
    assert status["release_ready"] is False
    assert status["restart_required"] is False
