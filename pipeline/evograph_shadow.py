"""Atomic sidecar materialization for the EvoGraph-F shadow prototype."""

from __future__ import annotations

import hashlib
import json
import os
import threading
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from matching.evograph_fragmentation import (
    SCHEMA_VERSION as GRAPH_SCHEMA_VERSION,
    build_evograph_graph,
    default_policy_path,
    load_evograph_policy,
)
from pipeline.result_contract import iter_legacy_result_artifact


SIDECAR_SCHEMA_VERSION = "1.0.0"
_lock_guard = threading.Lock()
_artifact_locks: dict[str, threading.Lock] = {}


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def evograph_sidecar_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}_evograph_fragmentation.json"


def _artifact_lock(path: Path) -> threading.Lock:
    key = str(path.resolve())
    with _lock_guard:
        return _artifact_locks.setdefault(key, threading.Lock())


def _read_metadata(path: Path) -> dict[str, Any]:
    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        return dict(payload.get("metadata") or {})
    with path.open("rb") as handle:
        return dict(next(ijson.items(handle, "metadata", use_float=True), {}) or {})


def materialize_evograph_sidecar(
    task_id: str,
    legacy_result_path: Path,
    output_path: Path,
    *,
    policy_path: Path | None = None,
    runtime_snapshot_id: str | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Read an immutable legacy artifact and atomically write a shadow graph."""

    if not legacy_result_path.is_file():
        raise FileNotFoundError(f"Legacy result artifact not found: {legacy_result_path}")
    policy_path = policy_path or default_policy_path()
    policy = load_evograph_policy(policy_path)
    source_sha256 = _sha256_file(legacy_result_path)
    policy_sha256 = _sha256_file(policy_path)

    with _artifact_lock(output_path):
        if output_path.is_file() and not force:
            existing = _read_metadata(output_path)
            if (
                existing.get("source_results_sha256") == source_sha256
                and existing.get("policy_sha256") == policy_sha256
            ):
                return {**existing, "artifact": output_path.name, "cache_hit": True}

        records: list[dict[str, Any]] = []
        semantic_digest = hashlib.sha256()
        status_counts: Counter[str] = Counter()
        formula_context_counts: Counter[str] = Counter()
        edge_status_counts: Counter[str] = Counter()
        graph_totals: Counter[str] = Counter()
        maximum = int(
            (policy.get("input_policy") or {}).get("maximum_spectra_per_task") or 5000
        )
        truncated = False

        for result_index, row in enumerate(iter_legacy_result_artifact(legacy_result_path)):
            if result_index >= maximum:
                truncated = True
                break
            graph = build_evograph_graph(row, result_index=result_index, policy=policy)
            records.append(graph)
            semantic_digest.update(_canonical(graph))
            semantic_digest.update(b"\n")
            status_counts[str(graph.get("status") or "unknown")] += 1
            context = graph.get("context") or {}
            formula_context_counts[
                str(context.get("formula_assignment_status") or "unknown")
            ] += 1
            edge_status_counts.update(graph.get("edge_status_counts") or {})
            audit = graph.get("graph_audit") or {}
            for field in (
                "node_count",
                "fragment_node_count",
                "candidate_edge_count",
                "selected_edge_count",
                "cycle_count",
                "direction_violations",
                "formula_conservation_violations",
                "formula_enumeration_truncated_nodes",
            ):
                graph_totals[field] += int(audit.get(field) or 0)

        source_sha256_after = _sha256_file(legacy_result_path)
        if source_sha256_after != source_sha256:
            raise RuntimeError("Legacy result artifact changed during EvoGraph materialization")

        generated_at = datetime.now(timezone.utc).isoformat()
        safety_passed = not (
            graph_totals["cycle_count"]
            or graph_totals["direction_violations"]
            or graph_totals["formula_conservation_violations"]
        )
        metadata = {
            "sidecar_schema_version": SIDECAR_SCHEMA_VERSION,
            "graph_schema_version": GRAPH_SCHEMA_VERSION,
            "task_id": task_id,
            "module_id": "evograph_fragmentation",
            "module_state": "shadow",
            "policy_id": policy.get("policy_id"),
            "policy_sha256": policy_sha256,
            "runtime_snapshot_id": runtime_snapshot_id,
            "generated_at": generated_at,
            "record_count": len(records),
            "truncated_by_task_limit": truncated,
            "maximum_spectra_per_task": maximum,
            "status_counts": dict(sorted(status_counts.items())),
            "formula_context_counts": dict(sorted(formula_context_counts.items())),
            "edge_status_counts": dict(sorted(edge_status_counts.items())),
            "graph_totals": dict(sorted(graph_totals.items())),
            "graph_semantic_sha256": semantic_digest.hexdigest(),
            "source_results_sha256": source_sha256,
            "source_results_sha256_after": source_sha256_after,
            "source_results_unchanged": True,
            "legacy_result_mutated": False,
            "production_ranking_changed_records": 0,
            "legacy_confidence_changed_records": 0,
            "high_promotions": 0,
            "log_bayes_factor_available": False,
            "heuristic_marginals_are_probabilities": False,
            "external_calibration_status": "unavailable",
            "safety_invariants_passed": safety_passed,
        }
        output = {"metadata": metadata, "records": records}
        output_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = output_path.with_name(f".{output_path.name}.tmp-{os.getpid()}-{threading.get_ident()}")
        try:
            with temporary.open("w", encoding="utf-8") as handle:
                json.dump(output, handle, ensure_ascii=False, separators=(",", ":"))
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, output_path)
        finally:
            temporary.unlink(missing_ok=True)
        return {
            **metadata,
            "artifact": output_path.name,
            "size_bytes": output_path.stat().st_size,
            "cache_hit": False,
        }


def read_evograph_sidecar(
    path: Path,
    *,
    offset: int = 0,
    limit: int = 100,
    result_index: int | None = None,
) -> dict[str, Any]:
    if offset < 0 or limit < 1:
        raise ValueError("offset must be >= 0 and limit must be >= 1")
    if result_index is not None and result_index < 0:
        raise ValueError("result_index must be >= 0")
    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        records = list(payload.get("records") or [])
        metadata = dict(payload.get("metadata") or {})
        if result_index is not None:
            selected = [record for record in records if int(record.get("result_index", -1)) == result_index]
            return {
                "metadata": metadata,
                "records": selected,
                "offset": result_index,
                "limit": 1,
                "has_more": False,
            }
        return {
            "metadata": metadata,
            "records": records[offset : offset + limit],
            "offset": offset,
            "limit": limit,
            "has_more": offset + limit < len(records),
        }

    metadata = _read_metadata(path)
    records: list[dict[str, Any]] = []
    with path.open("rb") as handle:
        for index, record in enumerate(ijson.items(handle, "records.item", use_float=True)):
            if result_index is not None:
                if int(record.get("result_index", -1)) == result_index:
                    records.append(record)
                    break
                continue
            if index < offset:
                continue
            if len(records) >= limit:
                break
            records.append(record)
    total = int(metadata.get("record_count") or 0)
    return {
        "metadata": metadata,
        "records": records,
        "offset": result_index if result_index is not None else offset,
        "limit": 1 if result_index is not None else limit,
        "has_more": False if result_index is not None else offset + len(records) < total,
    }


__all__ = [
    "SIDECAR_SCHEMA_VERSION",
    "evograph_sidecar_path",
    "materialize_evograph_sidecar",
    "read_evograph_sidecar",
]
