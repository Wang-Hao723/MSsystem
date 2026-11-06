"""Keep in-process API tests away from the user's live persistence files."""
import pytest


@pytest.fixture(autouse=True)
def isolated_server_persistence(monkeypatch, tmp_path):
    # TestClient executes the real startup hook. Its ledger lock must belong to
    # the test, not the running desktop instance. Individual tests may override
    # these paths with their own fixtures; production locking is not bypassed.
    from web import server

    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks_index.json")
    monkeypatch.setattr(server, "VENDOR_CONVERSION_LEDGER_FILE",
                        tmp_path / "runtime" / "vendor_conversion_ledger.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_restart_recovery_task_ids", set())
