"""Both web entrypoints bind the verified loader; failed reload clears the UI."""
import ast
from pathlib import Path
from types import SimpleNamespace

import pytest


APP = Path(__file__).resolve().parents[1] / "web/app.py"


def _session_loader(loader):
    tree = ast.parse(APP.read_text())
    node, = [n for n in tree.body if isinstance(n, ast.FunctionDef)
             and n.name == "_load_library_into_session"]
    state = SimpleNamespace(records=["old"], idx={"old": [0]}, lib_loaded=True,
                            lib_name="old", lib_spectra_count=1, lib_load_time=1.0)
    times = iter([100.0, 101.25])
    from runtime.legacy_library_session import LibraryCoordinator
    coordinator = LibraryCoordinator(runtime_binding=lambda: ())
    scope = {"st": SimpleNamespace(session_state=state), "load_library": loader,
             "time": SimpleNamespace(time=lambda: next(times)),
             "_library_coordinator": coordinator, "_tasks_lock": coordinator.lock}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(APP), "exec"), scope)
    return state, scope["_load_library_into_session"]


def test_web_surfaces_import_same_verified_callable():
    from runtime.library_metadata_repairs import load_library_with_verified_repairs
    from web import server
    tree = ast.parse(APP.read_text())
    imports = [n for n in tree.body if isinstance(n, ast.ImportFrom)
               and n.module == "runtime.library_metadata_repairs"]
    scope = {}
    exec(compile(ast.Module(body=imports, type_ignores=[]), str(APP), "exec"), scope)
    assert scope["load_library"] is server.load_library is load_library_with_verified_repairs
    assert not any(isinstance(n, ast.ImportFrom) and n.module == "matching"
                   and any(a.name == "load_library" for a in n.names) for n in tree.body)
    assert sum(isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
               and n.func.id == "_load_library_into_session" for n in ast.walk(tree)) == 2


def test_successful_load_publishes_actual_library_and_count():
    records = ["a", "b"]
    index = {10: [0, 1]}
    calls = []
    def loader(path):
        calls.append(path)
        assert not state.lib_loaded and state.records is None and state.idx is None
        return records, index
    state, load = _session_loader(loader)
    assert load("selected.db", "Full") == (2, 1.25)
    assert calls == ["selected.db"]
    assert state.lib_loaded and state.records is records and state.idx is index
    assert state.lib_name == "Full" and state.lib_spectra_count == 2
    assert state.lib_load_time == round(1.25, 1)


def test_failed_reload_cannot_leave_previous_library_marked_loaded():
    def fail(path):
        raise ValueError("source binding differs")
    state, load = _session_loader(fail)
    with pytest.raises(ValueError, match="source binding differs"):
        load("changed.db", "Full")
    assert state.lib_loaded is False and state.records is None and state.idx is None
    assert state.lib_name is None and state.lib_spectra_count == 0 and state.lib_load_time == 0.0
