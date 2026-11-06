"""Atomic V12 EviRecall-MV candidate-expansion sidecar contract."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from matching.evirank_retrieval import (
    build_candidate_recall_bundle,
)
from runtime.verified_sql_retrieval import create_verified_sql_retriever


EVIRECALL_ARTIFACT_SUFFIX = "_v12_candidate_recall.json"
_SAMPLE_BYTES = 1024 * 1024


def evirecall_artifact_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}{EVIRECALL_ARTIFACT_SUFFIX}"


def _sampled_file_sha256(path: Path) -> str:
    size = path.stat().st_size
    offsets = sorted(
        {
            0,
            max(0, (size - _SAMPLE_BYTES) // 2),
            max(0, size - _SAMPLE_BYTES),
        }
    )
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for offset in offsets:
            handle.seek(offset)
            digest.update(str(offset).encode("ascii"))
            digest.update(b"\0")
            digest.update(handle.read(_SAMPLE_BYTES))
    return digest.hexdigest()


def write_evirecall_artifact(
    task_id: str,
    v11_artifact: Path,
    legacy_results: Path,
    database: Path,
    output_path: Path,
    *,
    library_snapshot_id: str,
    library_sampled_sha256: str,
    detection_method_id: str = "unspecified",
    runtime_snapshot_id: str | None = None,
) -> dict[str, Any]:
    """Build EviRecall-MV without mutating V11, raw results or the library."""

    v11_path = Path(v11_artifact).resolve()
    source_path = Path(legacy_results).resolve()
    db_path = Path(database).resolve()
    for path in (v11_path, source_path, db_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    v11_before = v11_path.read_bytes()
    source_before = source_path.read_bytes()
    v11_payload = json.loads(v11_before)
    rows = json.loads(source_before)
    if not isinstance(v11_payload, dict):
        raise ValueError("V11 artifact must contain a JSON object")
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise ValueError("Legacy result artifact must contain a JSON array of objects")
    declared_task = str((v11_payload.get("metadata") or {}).get("task_id") or "")
    if declared_task and declared_task != task_id:
        raise ValueError("Task id does not match the V11 artifact")
    db_before = {
        "size_bytes": db_path.stat().st_size,
        "mtime_ns": db_path.stat().st_mtime_ns,
        "sampled_sha256": _sampled_file_sha256(db_path),
    }
    expected_sample = str(library_sampled_sha256 or "").strip().lower()
    if len(expected_sample) != 64 or db_before["sampled_sha256"] != expected_sample:
        raise ValueError("Shared library sampled hash does not match the frozen binding")
    retriever = create_verified_sql_retriever(
        db_path,
        library_snapshot_id=library_snapshot_id,
    )
    payload = build_candidate_recall_bundle(
        v11_payload,
        rows,
        retriever=retriever,
        detection_method_id=detection_method_id,
    )
    metadata = payload["metadata"]
    metadata.update(
        {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "runtime_snapshot_id": runtime_snapshot_id,
            "legacy_result_file_sha256": hashlib.sha256(source_before).hexdigest(),
            "artifact_role": "shadow_v12_candidate_recall_expansion",
            "production_artifact_role": "v11_flagship_unchanged",
            "database_open_mode": "sqlite_mode_ro_query_only",
            "library_sampled_sha256": expected_sample,
        }
    )
    db_after = {
        "size_bytes": db_path.stat().st_size,
        "mtime_ns": db_path.stat().st_mtime_ns,
        "sampled_sha256": _sampled_file_sha256(db_path),
    }
    if db_after != db_before:
        raise RuntimeError("Shared spectral library changed during EviRecall materialization")
    if v11_path.read_bytes() != v11_before:
        raise RuntimeError("V11 artifact changed during EviRecall materialization")
    if source_path.read_bytes() != source_before:
        raise RuntimeError("Legacy results changed during EviRecall materialization")
    metadata["database_sampled_fingerprint_unchanged"] = True

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, separators=(",", ":"), default=str)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, output_path)
    return {
        **metadata,
        "status": "available",
        "artifact": str(output_path.resolve()),
        "size_bytes": output_path.stat().st_size,
    }


def read_evirecall_artifact(
    path: Path,
    *,
    offset: int = 0,
    limit: int = 500,
    entity_id: str | None = None,
) -> dict[str, Any]:
    if offset < 0 or limit < 1:
        raise ValueError("offset must be >=0 and limit must be >=1")
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    entities = list(payload.get("entities") or [])
    wanted = str(entity_id or "").strip()
    if wanted:
        entities = [
            entity
            for entity in entities
            if str(entity.get("entity_id") or "") == wanted
        ]
    return {
        "metadata": payload.get("metadata") or {},
        "invariants": payload.get("invariants") or {},
        "entities": entities[offset : offset + limit],
        "entities_total": len(entities),
        "offset": offset,
        "limit": limit,
        "has_more": offset + limit < len(entities),
    }


__all__ = [
    "EVIRECALL_ARTIFACT_SUFFIX",
    "evirecall_artifact_path",
    "read_evirecall_artifact",
    "write_evirecall_artifact",
]
