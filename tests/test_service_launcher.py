"""Process-control regressions use only temporary files and synthetic HTTP data."""
from __future__ import annotations

import errno
import io
import json
import os
import signal
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest

from scripts import service_launcher as launcher


ROOT = Path(launcher.__file__).resolve().parents[1]
FAKE_SERVER = """
import json, sys
from http.server import BaseHTTPRequestHandler, HTTPServer
port, root = int(sys.argv[1]), sys.argv[2]
class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        payload = dict(status='ok', project_root=root, service_version='16.0.0', active_port=port)
        self.send_response(200)
        self.end_headers()
        self.wfile.write(json.dumps(payload).encode())
    def log_message(self, *args):
        pass
print('synthetic server started', flush=True)
HTTPServer(('127.0.0.1', port), Handler).serve_forever()
"""


def free_port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def test_service_survives_launcher_exit_and_refuses_duplicate(tmp_path):
    """The parent exits before we probe the child; this is not a Popen mock."""
    port = free_port()
    command = [sys.executable, "-u", "-c", FAKE_SERVER, str(port), str(tmp_path)]
    parent = (
        "import sys; from pathlib import Path; from scripts.service_launcher import start_service; "
        f"sys.exit(start_service(Path({str(tmp_path)!r}), {command!r}, '127.0.0.1', {port}, wait_seconds=10))"
    )
    pid_file = tmp_path / "output/runtime/mssystem_v16.pid"
    result = None
    try:
        result = subprocess.run([sys.executable, "-c", parent], cwd=ROOT,
                                capture_output=True, text=True, timeout=20)
        assert result.returncode == 0, result.stdout + result.stderr
        pid = int(pid_file.read_text())
        assert launcher.pid_is_live(pid_file)
        assert os.getsid(pid) == pid
        assert os.getpgid(pid) == pid
        assert launcher.check_health(tmp_path, "127.0.0.1", port, "fastapi")["status"] == "ok"
        assert (tmp_path / "output/runtime/mssystem_v16.endpoint").read_text() == f"127.0.0.1 {port}\n"
        log = (tmp_path / "output/runtime/mssystem_v16_service.log").read_text()
        assert "synthetic server started" in log
        assert launcher.start_service(tmp_path, command, "127.0.0.1", port) == 2
        assert int(pid_file.read_text()) == pid
        assert (tmp_path / "output/runtime/mssystem_v16_service.log").read_text() == log
    finally:
        if pid_file.exists():
            pid = int(pid_file.read_text())
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass


def test_startup_failure_preserves_error_log(tmp_path, capsys):
    result = launcher.start_service(tmp_path,
        [sys.executable, "-u", "-c", "import sys; print('synthetic startup error', file=sys.stderr); sys.exit(7)"],
        "127.0.0.1", free_port(), wait_seconds=5)
    assert result == 3
    assert "exit 7" in capsys.readouterr().err
    state = tmp_path / "output/runtime"
    assert "synthetic startup error" in (state / "mssystem_v16_service.log").read_text()
    assert launcher.pid_is_live(state / "mssystem_v16.pid") is False


def test_not_ready_keeps_owned_process_but_does_not_report_success(tmp_path, capsys):
    pid_file = tmp_path / "output/runtime/mssystem_v16.pid"
    try:
        result = launcher.start_service(tmp_path,
            [sys.executable, "-u", "-c", "import time; time.sleep(20)"],
            "127.0.0.1", free_port(), wait_seconds=.05)
        assert result == 4
        assert launcher.pid_is_live(pid_file)
        captured = capsys.readouterr()
        assert "not yet healthy" in captured.err
        assert "is ready" not in captured.out
    finally:
        if pid_file.exists():
            pid = int(pid_file.read_text())
            try:
                os.kill(pid, signal.SIGTERM)
                os.waitpid(pid, 0)
            except (ProcessLookupError, ChildProcessError):
                pass


def test_live_pid_is_not_replaced(tmp_path, monkeypatch):
    state = tmp_path / "output/runtime"
    state.mkdir(parents=True)
    pid_file = state / "mssystem_v16.pid"
    pid_file.write_text(str(os.getpid()))
    monkeypatch.setattr(launcher.subprocess, "Popen", lambda *a, **kw: pytest.fail("must not spawn"))
    assert launcher.start_service(tmp_path, ["unused"], "127.0.0.1", 8516) == 2
    assert pid_file.read_text() == str(os.getpid())


def test_concurrent_launcher_lock_prevents_spawn(tmp_path, monkeypatch):
    def busy(*args):
        raise BlockingIOError(errno.EAGAIN, "already locked")
    monkeypatch.setattr(launcher.fcntl, "flock", busy)
    monkeypatch.setattr(launcher.subprocess, "Popen", lambda *a, **kw: pytest.fail("must not spawn"))
    assert launcher.start_service(tmp_path, ["unused"], "127.0.0.1", 8516) == 2
    assert not (tmp_path / "output/runtime/mssystem_v16.pid").exists()


def health_opener(monkeypatch, body, seen=None):
    class Opener:
        def open(self, url, timeout):
            if seen is not None:
                seen.append((url, timeout))
            return io.BytesIO(body)
    def build(handler):
        assert handler.proxies == {}
        return Opener()
    monkeypatch.setattr(launcher.urllib.request, "build_opener", build)


@pytest.mark.parametrize("change", [
    {"status": "error"}, {"project_root": "/other/checkout"},
    {"active_port": 8515}, {"service_version": "15.0.0"},
])
def test_health_rejects_other_instance_or_failure(tmp_path, monkeypatch, change):
    payload = dict(status="ok", project_root=str(tmp_path), active_port=8516, service_version="16.0.0")
    payload.update(change)
    health_opener(monkeypatch, json.dumps(payload).encode())
    with pytest.raises(ValueError):
        launcher.check_health(tmp_path, "127.0.0.1", 8516, "fastapi")


@pytest.mark.parametrize("host,address", [("0.0.0.0", "127.0.0.1"), ("::", "[::1]"), ("localhost", "localhost")])
def test_probe_uses_bounded_local_direct_request(tmp_path, monkeypatch, host, address):
    payload = dict(status="ok", project_root=str(tmp_path), active_port=8516, service_version="16.0.0")
    seen = []
    health_opener(monkeypatch, json.dumps(payload).encode(), seen)
    assert launcher.check_health(tmp_path, host, 8516, "fastapi") == payload
    assert seen == [(f"http://{address}:8516/api/health", 30.0)]


def test_streamlit_has_its_own_health_endpoint(tmp_path, monkeypatch):
    seen = []
    health_opener(monkeypatch, b"ok", seen)
    assert launcher.check_health(tmp_path, "127.0.0.1", 8526, "streamlit")["status"] == "ok"
    assert seen[0][0].endswith(":8526/_stcore/health")
    health_opener(monkeypatch, b"not ready")
    with pytest.raises(ValueError):
        launcher.check_health(tmp_path, "127.0.0.1", 8526, "streamlit")


def test_failed_atomic_publication_keeps_original(tmp_path, monkeypatch):
    target = tmp_path / "process.pid"
    target.write_text("original\n")
    def fail(*args):
        raise OSError("synthetic publication failure")
    monkeypatch.setattr(launcher.os, "replace", fail)
    with pytest.raises(OSError):
        launcher.atomic_text(target, "replacement\n")
    assert target.read_text() == "original\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["process.pid"]


@pytest.mark.parametrize("value", ["", "not-a-pid", "0", "-1", "1"])
def test_invalid_pid_never_targets_a_process_group(tmp_path, monkeypatch, value):
    path = tmp_path / "process.pid"
    path.write_text(value)
    monkeypatch.setattr(launcher.os, "kill", lambda *args: pytest.fail("must not signal"))
    assert launcher.pid_is_live(path) is False


def test_status_cli_fails_when_health_check_fails(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["service_launcher", "--project-root", str(tmp_path), "--check"])
    def fail(*args):
        raise OSError("synthetic unavailable endpoint")
    monkeypatch.setattr(launcher, "check_health", fail)
    assert launcher.main() == 4
    assert "not healthy" in capsys.readouterr().err


def test_shell_defaults_to_detached_and_keeps_debug_option():
    source = (ROOT / "run.sh").read_text()
    assert "FOREGROUND=0" in source
    assert "--foreground)" in source
    assert "scripts/service_launcher.py" in source
    assert "show_status\n        exit $?" in source
    assert 'curl -fsS "http://${status_host}:${status_port}/api/health" 2>/dev/null || true' not in source


@pytest.mark.parametrize("delay,wait_seconds,expected", [(1.2,5.0,0),(3.0,.8,4)])
def test_one_slow_health_request_uses_the_remaining_startup_deadline(tmp_path, delay, wait_seconds, expected):
    """Real delayed HTTP proves readiness is neither rejected at 1s nor unbounded."""
    port = free_port()
    server = FAKE_SERVER.replace("import json, sys", "import json, sys, time")
    server = server.replace("    def do_GET(self):", "    def do_GET(self):\n        print('health-request', flush=True)\n        time.sleep(float(sys.argv[3]))")
    command = [sys.executable, "-u", "-c", server, str(port), str(tmp_path), str(delay)]
    pid_file = tmp_path / "output/runtime/mssystem_v16.pid"
    started = time.monotonic()
    try:
        result = launcher.start_service(tmp_path, command, "127.0.0.1", port, wait_seconds=wait_seconds)
        elapsed = time.monotonic() - started
        assert result == expected
        assert elapsed < wait_seconds + 1.0
        assert launcher.pid_is_live(pid_file)
        log = (tmp_path / "output/runtime/mssystem_v16_service.log").read_text()
        assert log.count("health-request") == 1
    finally:
        if pid_file.exists():
            pid = int(pid_file.read_text())
            try:
                os.kill(pid, signal.SIGTERM)
                os.waitpid(pid, 0)
            except (ProcessLookupError, ChildProcessError):
                pass
