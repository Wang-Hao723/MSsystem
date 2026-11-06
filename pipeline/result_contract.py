"""Materialize and read the versioned IdentificationRecord sidecar.

The legacy result JSON remains the immutable production source during P1.  This
artifact is additive: it proves that every consumer can converge on one DTO
without changing search scores, Top-1 selection or existing report payloads.
"""

from __future__ import annotations

import hashlib
import json
from runtime.json_encoding import write_json_encoded
import os
from collections import Counter
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from schemas.identification import IDENTIFICATION_SCHEMA_VERSION
from schemas.legacy import legacy_row_to_identification


CONTRACT_SCHEMA_VERSION = "1.0.0"


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def identification_contract_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}_identification_records.json"


def runtime_snapshot_id_from_manifest(project_root: Path) -> str | None:
    path = project_root / "config" / "runtime_manifest.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    value = str(payload.get("runtime_snapshot_id") or "").strip()
    return value or None


def iter_legacy_result_artifact(path: Path) -> Iterator[Mapping[str, Any]]:
    """Stream legacy result rows without materializing raw peak arrays in RAM."""
    try:
        import ijson

        with path.open("rb") as handle:
            for record in ijson.items(handle, "item", use_float=True):
                if isinstance(record, Mapping):
                    yield record
        return
    except ImportError:
        pass
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise ValueError(f"Expected a result array in {path}")
    for record in payload:
        if isinstance(record, Mapping):
            yield record


def write_identification_contract(
    task_id: str,
    legacy_results: Iterable[Mapping[str, Any]],
    output_path: Path,
    *,
    runtime_snapshot_id: str | None = None,
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    scientific_digest = hashlib.sha256()
    claim_scopes: Counter[str] = Counter()
    claim_statuses: Counter[str] = Counter()
    module_states: dict[str, Counter[str]] = {}
    evidence_nodes = 0
    ranking_changed = 0

    for row_index, row in enumerate(legacy_results):
        if not isinstance(row, Mapping):
            continue
        record = legacy_row_to_identification(
            task_id,
            row,
            row_index=row_index,
            runtime_snapshot_id=runtime_snapshot_id,
        )
        if record is None:
            continue
        payload = record.model_dump(mode="json", exclude_none=True)
        records.append(payload)
        scientific_digest.update(_canonical(payload))
        scientific_digest.update(b"\n")
        claim_scopes[record.decision.scope.value] += 1
        claim_statuses[record.decision.status] += 1
        evidence_nodes += len(record.evidence)
        ranking_changed += int(record.scores.ranking_changed)
        for module, state in record.module_states.items():
            module_states.setdefault(module, Counter())[state.value] += 1

    metadata = {
        "contract_schema_version": CONTRACT_SCHEMA_VERSION,
        "identification_schema_version": IDENTIFICATION_SCHEMA_VERSION,
        "task_id": task_id,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "runtime_snapshot_id": runtime_snapshot_id,
        "record_count": len(records),
        "evidence_node_count": evidence_nodes,
        "claim_scope_counts": dict(sorted(claim_scopes.items())),
        "claim_status_counts": dict(sorted(claim_statuses.items())),
        "module_state_counts": {
            module: dict(sorted(counts.items()))
            for module, counts in sorted(module_states.items())
        },
        "production_ranking_changed_records": ranking_changed,
        "record_semantic_sha256": scientific_digest.hexdigest(),
        "ordinary_ms2_identifications_only": True,
        "expected_compounds_counted": False,
        "legacy_result_mutated": False,
    }
    output = {"metadata": metadata, "records": records}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        write_json_encoded(output, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, output_path)
    metadata["artifact"] = str(output_path.resolve())
    metadata["size_bytes"] = output_path.stat().st_size
    return metadata


def read_identification_contract(
    path: Path,
    *,
    offset: int = 0,
    limit: int = 500,
) -> dict[str, Any]:
    if offset < 0 or limit < 1:
        raise ValueError("offset must be >=0 and limit must be >=1")
    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        records = payload.get("records") or []
        metadata = payload.get("metadata") or {}
        return {
            "metadata": metadata,
            "records": records[offset : offset + limit],
            "offset": offset,
            "limit": limit,
            "has_more": offset + limit < len(records),
        }

    with path.open("rb") as handle:
        metadata = next(ijson.items(handle, "metadata", use_float=True), {})
    records = []
    seen = 0
    with path.open("rb") as handle:
        for seen, record in enumerate(
            ijson.items(handle, "records.item", use_float=True), start=1
        ):
            index = seen - 1
            if index < offset:
                continue
            if len(records) >= limit:
                break
            records.append(record)
    total = int(metadata.get("record_count") or 0)
    return {
        "metadata": metadata,
        "records": records,
        "offset": offset,
        "limit": limit,
        "has_more": offset + len(records) < total,
    }
