"""Execute the actual API loader's source without touching the live library."""
import ast
import asyncio
from pathlib import Path
import threading

import pytest


def _api_loader(tmp_path, callback):
    from web import server
    source = Path(__file__).resolve().parents[1] / "web/server.py"
    node, = [n for n in ast.parse(source.read_text()).body
             if isinstance(n, ast.AsyncFunctionDef) and n.name == "library_load"]
    node.decorator_list = []
    database = tmp_path / "fixture.db"
    database.touch()
    async def recovery(mode):
        return {"resumed": [], "deferred": [], "failed": []}
    scope = {**server.__dict__, "DB_PATH": database, "DB_PATH_SMALL": database,
             "_lib_loaded": True, "_lib_records": ["old"], "_mz_index": {1: [0]},
             "_lib_name": "old", "_lib_spectra_count": 1, "_lib_load_time": 1.0,
             "_library_state_lock": threading.RLock(), "_tasks_lock": threading.RLock(),
             "_tasks": {}, "load_library": callback, "_resume_recovered_tasks": recovery}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(source), "exec"), scope)
    return scope, scope["library_load"]


def test_failed_verification_clears_previous_public_library(tmp_path):
    from fastapi import HTTPException
    def fail(path):
        raise ValueError("changed source binding")
    scope, load = _api_loader(tmp_path, fail)
    with pytest.raises(HTTPException) as error:
        asyncio.run(load("full"))
    assert error.value.status_code == 500 and "changed source binding" in error.value.detail
    assert scope["_lib_loaded"] is False
    assert scope["_lib_records"] is None and scope["_mz_index"] is None
    assert scope["_lib_name"] == "" and scope["_lib_spectra_count"] == 0


def test_successful_verification_publishes_one_generation(tmp_path):
    records, index = ["one", "two"], {10: [0, 1]}
    def accepted(path):
        assert scope["_lib_loaded"] is False and scope["_lib_records"] is None
        return records, index
    scope, load = _api_loader(tmp_path, accepted)
    response = asyncio.run(load("full"))
    assert response.loaded and response.spectra_count == 2
    assert scope["_lib_loaded"] is True and scope["_lib_records"] is records
    assert scope["_mz_index"] is index and scope["_lib_spectra_count"] == 2


def test_active_analysis_rejection_preserves_existing_library(tmp_path):
    from fastapi import HTTPException
    def forbidden(path):
        raise AssertionError("Must reject before loading")
    scope, load = _api_loader(tmp_path, forbidden)
    scope["_tasks"] = {"running": {"status": "running"}}
    with pytest.raises(HTTPException) as error:
        asyncio.run(load("full"))
    assert error.value.status_code == 409
    assert scope["_lib_loaded"] is True and scope["_lib_records"] == ["old"]
