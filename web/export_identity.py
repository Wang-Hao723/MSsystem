"""Content-bound identities for persisted report exports.

The export cache deliberately hashes the small, task-local artifacts that are
actually read by report generation.  Large shared library stores remain bound
through their already-verified runtime/library snapshot identities.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time
from typing import Any, Iterable, Mapping


def file_sha256(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    """Return a content identity; missing optional inputs are explicit."""
    path = Path(path)
    started = time.perf_counter()
    try:
        stat = path.stat()
    except FileNotFoundError:
        return {"path": str(path), "status": "missing"}
    if not path.is_file():
        return {"path": str(path), "status": "not_file"}
    return {
        "path": str(path),
        "status": "file",
        "size_bytes": stat.st_size,
        "sha256": file_sha256(path),
        "hash_wall_seconds": time.perf_counter() - started,
    }


def build_export_identity(
    artifacts: Mapping[str, Path] | Iterable[tuple[str, Path]],
    *,
    bindings: Mapping[str, Any],
) -> dict[str, Any]:
    """Hash named task artifacts and stable, already-verified bindings."""
    items = artifacts.items() if isinstance(artifacts, Mapping) else artifacts
    records = {
        str(name): file_identity(Path(path))
        for name, path in sorted(items, key=lambda item: str(item[0]))
    }
    identity_records = {
        name: {key: value for key, value in record.items()
               if key != "hash_wall_seconds"}
        for name, record in records.items()
    }
    payload = {"schema": "mssystem-export-identity/1",
               "artifacts": identity_records, "bindings": dict(bindings)}
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
        default=str,
    ).encode("utf-8")
    return {
        "sha256": hashlib.sha256(canonical).hexdigest(),
        "artifacts": identity_records,
        "bindings": dict(bindings),
        "identity_wall_seconds": sum(
            float(record.get("hash_wall_seconds") or 0.0)
            for record in records.values()
        ),
    }
