"""Bounded, full-content verification of registered external source files.

An OS read can remain blocked after an HTTP request times out. Keep all source
I/O in an owned subprocess, so one stuck descriptor cannot hold a server worker
indefinitely. A retry opens a fresh descriptor; neither attempt skips the hash.
"""
from __future__ import annotations

import hashlib
import json
import logging
import math
import os
from pathlib import Path
import stat
import subprocess
import sys
import time

_LOG = logging.getLogger(__name__)
_WORKER = Path(__file__).resolve()
_STAT_FIELDS = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")


def _identity(value):
    return tuple(getattr(value, field) for field in _STAT_FIELDS)


def _read_source(path):
    # This function runs only in the disposable worker, including stat/open.
    path = Path(path)
    before = path.stat()
    if not stat.S_ISREG(before.st_mode):
        raise ValueError("Registered source must be a regular file")
    digest = hashlib.sha256()
    size = 0
    with path.open("rb", buffering=0) as handle:
        opened = os.fstat(handle.fileno())
        if _identity(before) != _identity(opened):
            raise ValueError("Source file changed while opening")
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
            size += len(chunk)
        after = os.fstat(handle.fileno())
    if (_identity(before) != _identity(after)
            or _identity(before) != _identity(path.stat())
            or size != before.st_size):
        raise ValueError("Source file changed while hashing")
    return {"sha256": digest.hexdigest(), "size_bytes": size}


def verify_source_file(path, expected, *, timeout_seconds=30.0, attempts=2):
    """Verify every byte against the registry or raise within a bounded wait.

Only a timeout is retried. Missing files, altered bytes and changed file
identities fail immediately. The source and its registry are never modified.
"""
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("Source verification timeout must be positive and finite")
    if not isinstance(attempts, int) or isinstance(attempts, bool) or not 1 <= attempts <= 2:
        raise ValueError("Source verification permits one or two attempts")
    expected_size = expected.get("size_bytes")
    expected_hash = expected.get("sha256")
    if (not isinstance(expected_size, int) or isinstance(expected_size, bool)
            or expected_size < 0 or not isinstance(expected_hash, str)
            or len(expected_hash) != 64
            or any(c not in "0123456789abcdef" for c in expected_hash)):
        raise ValueError("Invalid registered source content identity")
    name = str(path)
    started = time.monotonic()
    for attempt in range(1, attempts + 1):
        _LOG.info("Verifying source %s (attempt %d/%d)", name, attempt, attempts)
        process = subprocess.Popen(
            [sys.executable, str(_WORKER), "--worker"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, close_fds=True,
        )
        try:
            stdout, stderr = process.communicate(
                json.dumps({"path": name}), timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            process.kill()  # Only the verifier subprocess created above.
            try:
                process.communicate(timeout=1.0)
            except subprocess.TimeoutExpired:
                # An uninterruptible OS read must not block the server again.
                _LOG.error("Source verifier PID %d still awaits OS cleanup", process.pid)
            if attempt < attempts:
                _LOG.warning("Source read timed out; retrying a fresh descriptor: %s", name)
                continue
            raise TimeoutError(
                f"Source verification timed out after {attempts} attempts: {name}") from None
        if process.returncode != 0:
            raise ValueError(f"Source verification failed: {name}: {stderr.strip()[:500]}")
        try:
            observed = json.loads(stdout)
        except (TypeError, ValueError) as exc:
            raise ValueError("Invalid source verification worker result: " + name) from exc
        if observed != {"size_bytes": expected_size, "sha256": expected_hash}:
            raise ValueError("Source file content changed: " + name)
        return {**observed, "path": name, "attempts": attempt,
                "elapsed_seconds": time.monotonic() - started}


def verify_registered_sources(source_files):
    return [verify_source_file(name, expected)
            for name, expected in source_files.items()]


if __name__ == "__main__":
    if sys.argv[1:] != ["--worker"]:
        raise SystemExit("Use --worker with a JSON path on standard input")
    try:
        print(json.dumps(_read_source(json.load(sys.stdin)["path"])))
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
