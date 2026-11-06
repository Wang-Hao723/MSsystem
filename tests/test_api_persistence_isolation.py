from fastapi.testclient import TestClient

from web import server


def test_default_api_test_persistence_is_isolated_from_the_desktop_instance(tmp_path):
    assert server.TASKS_FILE.resolve().is_relative_to(tmp_path.resolve())
    assert server.VENDOR_CONVERSION_LEDGER_FILE.resolve().is_relative_to(tmp_path.resolve())
    assert server.TASKS_FILE.resolve() != (server.PROJECT_ROOT / "output/tasks_index.json").resolve()
    with TestClient(server.app):
        assert server._vendor_conversion_service_lock_handle is not None
        assert server._vendor_conversion_service_lock_path.is_relative_to(tmp_path.resolve())
        assert server._tasks == {}
    assert server._vendor_conversion_service_lock_handle is None
    assert server.TASKS_FILE.is_file()
