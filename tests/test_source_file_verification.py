import hashlib
import time

import pytest

from runtime import source_file_verification as verifier


def identity(data):
    return {"size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def test_full_content_not_only_file_size_is_verified(tmp_path):
    path = tmp_path / "source.mgf"
    data = b"BEGIN IONS\n" + bytes(range(256)) * 8192 + b"END IONS\n"
    path.write_bytes(data)
    result = verifier.verify_source_file(path, identity(data))
    assert result["sha256"] == identity(data)["sha256"]
    assert result["size_bytes"] == len(data)
    changed = bytearray(data)
    changed[len(changed) // 2] ^= 1
    path.write_bytes(changed)
    with pytest.raises(ValueError, match="content changed"):
        verifier.verify_source_file(path, identity(data))


def test_missing_source_is_reported_without_hanging(tmp_path):
    with pytest.raises(ValueError, match="FileNotFoundError"):
        verifier.verify_source_file(tmp_path / "missing.mgf", identity(b""))


def test_wrong_size_fails(tmp_path):
    path = tmp_path / "source.mgf"
    path.write_bytes(b"x")
    expected = identity(b"x")
    expected["size_bytes"] += 1
    with pytest.raises(ValueError, match="content changed"):
        verifier.verify_source_file(path, expected)


def test_timeout_kills_owned_worker_and_does_not_block_server(tmp_path, monkeypatch):
    worker = tmp_path / "stuck_worker.py"
    worker.write_text("import time\ntime.sleep(60)\n")
    monkeypatch.setattr(verifier, "_WORKER", worker)
    started = time.monotonic()
    with pytest.raises(TimeoutError, match="after 2 attempts"):
        verifier.verify_source_file("unreadable", identity(b""),
                                    timeout_seconds=0.15)
    assert time.monotonic() - started < 3


def test_timeout_retry_uses_new_worker_then_checks_hash(tmp_path, monkeypatch):
    marker = tmp_path / "first_attempt"
    worker = tmp_path / "retry_worker.py"
    worker.write_text(
        "import time,json\nfrom pathlib import Path\n"
        f"marker=Path({str(marker)!r})\n"
        "if not marker.exists():\n marker.touch()\n time.sleep(60)\n"
        f"print(json.dumps({identity(b'actual')!r}))\n")
    monkeypatch.setattr(verifier, "_WORKER", worker)
    result = verifier.verify_source_file("source", identity(b"actual"),
                                         timeout_seconds=0.3)
    assert result["attempts"] == 2


@pytest.mark.parametrize("timeout", [0, -1, float("nan"), float("inf")])
def test_invalid_deadline_rejected(timeout):
    with pytest.raises(ValueError, match="timeout"):
        verifier.verify_source_file("source", identity(b""), timeout_seconds=timeout)


def test_nonregular_file_rejected(tmp_path):
    with pytest.raises(ValueError, match="regular file"):
        verifier.verify_source_file(tmp_path, identity(b""))
