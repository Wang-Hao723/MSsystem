from __future__ import annotations

from web import server


def test_progress_persistence_is_throttled_and_skips_backup_rotation(monkeypatch):
    calls = []

    def fake_save_tasks(**kwargs):
        calls.append(kwargs)
        return True

    monkeypatch.setattr(server, "_save_tasks", fake_save_tasks)
    server._task_progress_last_persisted.clear()

    for tick in range(1000):
        assert server._save_task_progress_throttled(
            "task-a", now_monotonic=tick / 100.0
        )

    assert calls == [{"backup_previous": False}]

    assert server._save_task_progress_throttled(
        "task-a", now_monotonic=10.0
    )
    assert calls == [
        {"backup_previous": False},
        {"backup_previous": False},
    ]


def test_forced_progress_save_and_failed_save_retry_immediately(monkeypatch):
    outcomes = iter([False, True, True])
    calls = []

    def fake_save_tasks(**kwargs):
        calls.append(kwargs)
        return next(outcomes)

    monkeypatch.setattr(server, "_save_tasks", fake_save_tasks)
    server._task_progress_last_persisted.clear()

    assert not server._save_task_progress_throttled(
        "task-b", now_monotonic=1.0
    )
    assert server._save_task_progress_throttled(
        "task-b", now_monotonic=1.1
    )
    assert server._save_task_progress_throttled(
        "task-b", force=True, now_monotonic=1.2
    )
    assert calls == [
        {"backup_previous": False},
        {"backup_previous": False},
        {"backup_previous": False},
    ]
