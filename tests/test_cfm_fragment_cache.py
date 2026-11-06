import hashlib
import json
import subprocess

from utils import cfm_fragment


def _cache_file(tmp_path, smiles, ion="-"):
    digest = hashlib.md5(smiles.encode("utf-8")).hexdigest()[:16]
    return tmp_path / f"v2_d1_{ion}_{digest}.json"


def test_empty_negative_cache_is_retried_for_interactive_request(tmp_path, monkeypatch):
    smiles = "CCO"
    monkeypatch.setattr(cfm_fragment, "CACHE_DIR", tmp_path)
    cache_file = _cache_file(tmp_path, smiles)
    cache_file.write_text("[]", encoding="utf-8")
    predicted = [{"mass": 31.018, "smiles": "CO"}]
    calls = []

    def fake_run(*args, **kwargs):
        calls.append((args, kwargs))
        return predicted

    monkeypatch.setattr(cfm_fragment, "run_fraggraph", fake_run)
    assert cfm_fragment.get_fragments(smiles, ionization="-", cache_only=False) == predicted
    assert len(calls) == 1
    assert json.loads(cache_file.read_text(encoding="utf-8")) == predicted


def test_empty_cache_remains_nonblocking_for_cache_only_request(tmp_path, monkeypatch):
    smiles = "CCO"
    monkeypatch.setattr(cfm_fragment, "CACHE_DIR", tmp_path)
    _cache_file(tmp_path, smiles).write_text("[]", encoding="utf-8")

    def must_not_run(*args, **kwargs):
        raise AssertionError("cache_only must not invoke Docker")

    monkeypatch.setattr(cfm_fragment, "run_fraggraph", must_not_run)
    assert cfm_fragment.get_fragments(smiles, ionization="-", cache_only=True) == []


def test_fraggraph_timeout_removes_exact_named_container(monkeypatch):
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        if command[:2] == ["docker", "run"]:
            raise subprocess.TimeoutExpired(command, timeout=kwargs["timeout"])
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(cfm_fragment.subprocess, "run", fake_run)
    assert cfm_fragment.run_fraggraph("CCO", timeout=1) == []
    run_command = calls[0][0]
    container_name = run_command[run_command.index("--name") + 1]
    assert container_name.startswith("mssystem-cfmid-")
    assert calls[1][0] == ["docker", "rm", "-f", container_name]
