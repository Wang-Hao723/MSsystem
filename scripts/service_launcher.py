#!/usr/bin/env python3
"""Start this checkout independently of the invoking terminal, with a health gate.

This is not a login service or an automatic restart loop. A timed-out startup
keeps its PID and log for diagnosis; it is never reported as a healthy service.
"""
from __future__ import annotations

import argparse
import errno
import fcntl
import json
import os
import subprocess
import sys
import tempfile
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def atomic_text(path: Path, text: str) -> None:
    """Publish small process-control files without exposing truncated contents."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix=path.name + ".", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def pid_is_live(pid_file: Path) -> bool:
    try:
        pid = int(pid_file.read_text(encoding="utf-8").strip())
    except (FileNotFoundError, ValueError):
        return False
    if pid <= 1:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # Never replace another user's live process-control record.
    return True


def check_health(project_root: Path, host: str, port: int, mode: str,
                 timeout: float = 30.0) -> dict:
    address = {"0.0.0.0": "127.0.0.1", "::": "::1"}.get(host, host)
    if ":" in address and not address.startswith("["):
        address = "[" + address + "]"
    endpoint = "/_stcore/health" if mode == "streamlit" else "/api/health"
    # Local service probes must not go through a workstation HTTP proxy.
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    with opener.open(f"http://{address}:{port}{endpoint}", timeout=timeout) as response:
        body = response.read(1024 * 1024)
    if mode == "streamlit":
        if body.strip() != b"ok":
            raise ValueError("The Streamlit health response is not ready")
        return {"status": "ok", "mode": mode}
    payload = json.loads(body)
    if (not isinstance(payload, dict) or payload.get("status") != "ok"
            or payload.get("project_root") != str(project_root.resolve())
            or payload.get("active_port") != port
            or not str(payload.get("service_version", "")).startswith("16.")):
        raise ValueError("Health response does not identify this V16 checkout and port")
    return payload


def start_service(project_root: Path, command: list[str], host: str, port: int,
                  mode: str = "fastapi", wait_seconds: float = 30) -> int:
    project_root = project_root.resolve()
    state = project_root / "output" / "runtime"
    state.mkdir(parents=True, exist_ok=True)
    pid_file = state / "mssystem_v16.pid"
    log_path = state / "mssystem_v16_service.log"
    with (state / "mssystem_v16.start.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            if exc.errno not in (errno.EACCES, errno.EAGAIN):
                raise
            print("Another V16 startup is in progress; no second process was started.", file=sys.stderr)
            return 2
        # Recheck under the lock: the shell's earlier ownership check can race.
        if pid_is_live(pid_file):
            print("A live managed PID already exists; inspect bash run.sh --status.", file=sys.stderr)
            return 2
        with log_path.open("ab", buffering=0) as log:
            stamp = datetime.now(timezone.utc).isoformat()
            log.write(f"\n[{stamp}] Starting {mode} at {host}:{port}\n".encode())
            process = subprocess.Popen(command, cwd=project_root, stdin=subprocess.DEVNULL,
                                       stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True, close_fds=True)
        try:
            atomic_text(pid_file, f"{process.pid}\n")
            atomic_text(state / "mssystem_v16.endpoint", f"{host} {port}\n")
        except OSError:
            # Only this just-created child may be stopped on publication failure.
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)
            raise
        deadline = time.monotonic() + wait_seconds
        last_error = "The health endpoint has not responded"
        while time.monotonic() < deadline:
            result = process.poll()
            if result is not None:
                print(f"V16 startup failed (exit {result}). Log: {log_path}", file=sys.stderr)
                return 3
            try:
                # A cold integrity check may take minutes. Keep one request
                # open within the startup deadline instead of abandoning it
                # every second while the server continues the same work.
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    break
                check_health(project_root, host, port, mode, timeout=remaining)
            except (OSError, ValueError) as exc:
                last_error = str(exc)
                time.sleep(min(.1, max(0.0, deadline - time.monotonic())))
                continue
            if process.poll() is None:
                print(f"MSsystem V16 is ready in the background (PID {process.pid}).")
                print(f"Log: {log_path}")
                return 0
        print(f"V16 PID {process.pid} is not yet healthy: {last_error}. "
              f"It remains available for diagnosis. Log: {log_path}", file=sys.stderr)
        return 4


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8516)
    parser.add_argument("--mode", choices=("fastapi", "streamlit"), default="fastapi")
    parser.add_argument("--wait-seconds", type=float, default=30)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.port <= 65535 or not 0 < args.wait_seconds <= 3600:
        parser.error("Use a valid port and a startup wait greater than 0 and at most 3600 seconds")
    if args.check:
        try:
            payload = check_health(args.project_root, args.host, args.port, args.mode)
        except (OSError, ValueError) as exc:
            print(f"V16 is not healthy: {exc}", file=sys.stderr)
            return 4
        print(json.dumps(payload, ensure_ascii=False))
        return 0
    if args.mode == "streamlit":
        command = [args.python, "-u", "-m", "streamlit", "run", "web/app.py",
                   "--server.port", str(args.port), "--server.address", args.host]
    else:
        command = [args.python, "-u", "-m", "uvicorn", "web.server:app",
                   "--host", args.host, "--port", str(args.port)]
    try:
        return start_service(args.project_root, command, args.host, args.port,
                             args.mode, args.wait_seconds)
    except OSError as exc:
        print(f"Could not start V16: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
