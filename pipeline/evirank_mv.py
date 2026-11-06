"""Atomic sidecar contract for V12-Flagship / EviRank-MV."""

from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from matching.evirank_mv import build_evirank_mv_result


EVIRANK_ARTIFACT_SUFFIX = "_v12_evirank_mv.json"


def evirank_artifact_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}{EVIRANK_ARTIFACT_SUFFIX}"


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def _read_optional_json(value: Path | Mapping[str, Any] | None) -> Mapping[str, Any] | None:
    if value is None:
        return None
    if isinstance(value, Mapping):
        return value
    payload = json.loads(Path(value).read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise ValueError(f"Expected a JSON object in {value}")
    return payload


def write_evirank_artifact(
    task_id: str,
    v11_artifact: Path | Mapping[str, Any],
    output_path: Path,
    *,
    detection_method_id: str = "unspecified",
    view_bundle: Path | Mapping[str, Any] | None = None,
    candidate_recall_bundle: Path | Mapping[str, Any] | None = None,
    residual_model: Path | Mapping[str, Any] | None = None,
    runtime_snapshot_id: str | None = None,
) -> dict[str, Any]:
    """Build V12 without mutating the V11 source artifact."""

    source_path = Path(v11_artifact) if not isinstance(v11_artifact, Mapping) else None
    before_bytes = source_path.read_bytes() if source_path is not None else None
    v11_payload = (
        json.loads(before_bytes)
        if before_bytes is not None
        else dict(v11_artifact)
    )
    if not isinstance(v11_payload, Mapping):
        raise ValueError("V11 artifact must be a JSON object")
    declared_task = str((v11_payload.get("metadata") or {}).get("task_id") or "")
    if declared_task and declared_task != task_id:
        raise ValueError("Task id does not match the V11 artifact")
    payload = build_evirank_mv_result(
        v11_payload,
        detection_method_id=detection_method_id,
        view_bundle=_read_optional_json(view_bundle),
        candidate_recall_bundle=_read_optional_json(candidate_recall_bundle),
        residual_model=_read_optional_json(residual_model),
    )
    entity_digest = hashlib.sha256()
    for entity in payload.get("entities") or []:
        entity_digest.update(_canonical(entity))
        entity_digest.update(b"\n")
    metadata = payload["metadata"]
    metadata.update({
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "runtime_snapshot_id": runtime_snapshot_id,
        "entity_semantic_sha256": entity_digest.hexdigest(),
        "artifact_role": "shadow_v12_evirank_mv_entity_ranking",
        "production_artifact_role": "v11_flagship_unchanged",
    })
    if source_path is not None and source_path.read_bytes() != before_bytes:
        raise RuntimeError("V11 artifact changed while V12 was materialized")
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


def read_evirank_artifact(
    path: Path,
    *,
    offset: int = 0,
    limit: int = 500,
    entity_id: str | None = None,
) -> dict[str, Any]:
    if offset < 0 or limit < 1:
        raise ValueError("offset must be >=0 and limit must be >=1")
    wanted = str(entity_id or "").strip()
    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        entities = list(payload.get("entities") or [])
        if wanted:
            entities = [entity for entity in entities if str(entity.get("entity_id") or "") == wanted]
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
    total = 0
    with path.open("rb") as handle:
        for entity in ijson.items(handle, "entities.item", use_float=True):
            if wanted and str(entity.get("entity_id") or "") != wanted:
                continue
            if total >= offset and len(selected) < limit:
                selected.append(entity)
            total += 1
    return {
        "metadata": metadata,
        "invariants": invariants,
        "entities": selected,
        "entities_total": total,
        "offset": offset,
        "limit": limit,
        "has_more": offset + len(selected) < total,
    }


__all__ = [
    "EVIRANK_ARTIFACT_SUFFIX",
    "evirank_artifact_path",
    "read_evirank_artifact",
    "write_evirank_artifact",
]
