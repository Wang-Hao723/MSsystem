"""Versioned sidecar contract for the V11 flagship entity decisions."""

from __future__ import annotations

import hashlib
import json
from runtime.json_encoding import write_json_encoded
import os
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from matching.flagship_ranker import build_flagship_result


FLAGSHIP_ARTIFACT_SUFFIX = "_v11_flagship.json"


def flagship_artifact_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}{FLAGSHIP_ARTIFACT_SUFFIX}"


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def write_flagship_artifact(
    task_id: str,
    legacy_results: Iterable[Mapping[str, Any]],
    output_path: Path,
    *,
    runtime_snapshot_id: str | None = None,
) -> dict[str, Any]:
    """Build and atomically write the primary V11 decision artifact."""

    payload = build_flagship_result(task_id, legacy_results)
    entity_digest = hashlib.sha256()
    for entity in payload.get("entities") or []:
        entity_digest.update(_canonical(entity))
        entity_digest.update(b"\n")
    metadata = payload["metadata"]
    metadata.update(
        {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "runtime_snapshot_id": runtime_snapshot_id,
            "entity_semantic_sha256": entity_digest.hexdigest(),
            "artifact_role": "primary_v11_entity_identification",
            "legacy_scan_result_role": "candidate_generation_and_raw_evidence",
        }
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        write_json_encoded(payload, handle, ensure_ascii=False, separators=(",", ":"), default=str)
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


def read_flagship_artifact(
    path: Path,
    *,
    offset: int = 0,
    limit: int = 500,
    decision_status: str | None = None,
    entity_id: str | None = None,
) -> dict[str, Any]:
    """Read a paginated entity projection without loading a large sidecar."""

    if offset < 0 or limit < 1:
        raise ValueError("offset must be >=0 and limit must be >=1")
    wanted_status = str(decision_status or "").strip()
    wanted_entity = str(entity_id or "").strip()
    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        entities = list(payload.get("entities") or [])
        if wanted_status:
            entities = [
                entity for entity in entities
                if str((entity.get("decision") or {}).get("status") or "") == wanted_status
            ]
        if wanted_entity:
            entities = [entity for entity in entities if str(entity.get("entity_id") or "") == wanted_entity]
        return {
            "metadata": payload.get("metadata") or {},
            "invariants": payload.get("invariants") or {},
            "entities": entities[offset : offset + limit],
            "entities_total": len(entities),
            "offset": offset,
            "limit": limit,
            "has_more": offset + limit < len(entities),
        }

    with path.open("rb") as handle:
        metadata = next(ijson.items(handle, "metadata", use_float=True), {})
    with path.open("rb") as handle:
        invariants = next(ijson.items(handle, "invariants", use_float=True), {})
    selected = []
    matched_total = 0
    unfiltered_total = (
        int(metadata.get("entity_count") or 0)
        if not wanted_status and not wanted_entity
        else None
    )
    with path.open("rb") as handle:
        for entity in ijson.items(handle, "entities.item", use_float=True):
            if wanted_status and str((entity.get("decision") or {}).get("status") or "") != wanted_status:
                continue
            if wanted_entity and str(entity.get("entity_id") or "") != wanted_entity:
                continue
            if matched_total >= offset and len(selected) < limit:
                selected.append(entity)
            matched_total += 1
            if unfiltered_total is not None and len(selected) >= limit:
                break
    if unfiltered_total is not None:
        matched_total = unfiltered_total
    return {
        "metadata": metadata,
        "invariants": invariants,
        "entities": selected,
        "entities_total": matched_total,
        "offset": offset,
        "limit": limit,
        "has_more": offset + len(selected) < matched_total,
    }


__all__ = [
    "FLAGSHIP_ARTIFACT_SUFFIX",
    "flagship_artifact_path",
    "read_flagship_artifact",
    "write_flagship_artifact",
]
