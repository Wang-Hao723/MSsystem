"""Bounded, read-only, on-demand candidate comparison for saved V16 entities."""
from __future__ import annotations

import asyncio
import copy
import hashlib
import re
import threading
from collections import OrderedDict
from pathlib import Path

import ijson
from fastapi import HTTPException, Query
from fastapi.responses import JSONResponse
from pydantic import BaseModel, ConfigDict, Field

from matching.library_record_store import semantic_sha256
from pipeline.v16_candidate_comparison import ALGORITHM_ID
from pipeline.v16_complete_candidate_panel import project_query, rescore_complete_panel


class CandidateComparisonRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    entity_id: str = Field(min_length=1, max_length=256)


def _connectivity(item):
    if not isinstance(item, dict):
        return None
    key = (item.get("connectivity_identity") or item.get("identity_key")
           or (item.get("candidate") or {}).get("inchikey") or item.get("inchikey") or "")
    key = str(key).split(":")[-1].split("-")[0].strip().upper()
    return key if re.fullmatch(r"[A-Z]{14}", key) else None


def _stamp(path):
    path = Path(path)
    stat = path.stat()
    return (str(path.resolve()), stat.st_ino, stat.st_size, stat.st_mtime_ns)


class CandidateComparisonService:
    """No production artifact writes, scan conversion, or full-library preload.

    One request at a time; at most 200 saved queries and 40 candidate keys.
    A small immutable-result cache avoids repeating the expensive comparison.
    All dependencies are explicit so tests need no production database.
    """
    def __init__(self, *, task_directory, task_lookup, artifact_path, artifact_current,
                 read_entity, retriever_factory, runtime_snapshot):
        self.task_directory = task_directory
        self.task_lookup = task_lookup
        self.artifact_path = artifact_path
        self.artifact_current = artifact_current
        self.read_entity = read_entity
        self.retriever_factory = retriever_factory
        self.runtime_snapshot = runtime_snapshot
        self.lock = threading.Lock()
        self.cache = OrderedDict()
        self.maximum_cache_bytes = 32 * 1024 * 1024

    def compute(self, task_id, entity_id):
        if not re.fullmatch(r"[A-Za-z0-9_-]{1,96}", task_id):
            raise HTTPException(404, "Task not found")
        task = self.task_lookup(task_id)
        if not task:
            raise HTTPException(404, "Task not found")
        if task.get("status") != "completed":
            raise HTTPException(409, "Only completed tasks can be compared; running analysis is not interrupted")
        if not self.lock.acquire(blocking=False):
            raise HTTPException(429, "One candidate comparison is already running; retry after it finishes")
        try:
            return self._compute(task_id, entity_id)
        finally:
            self.lock.release()

    def _compute(self, task_id, entity_id):
        directory = Path(self.task_directory())
        artifact = Path(self.artifact_path(directory, task_id))
        source = directory / f"{task_id}_results.json"
        if not source.is_file() or not artifact.is_file() or not self.artifact_current(artifact):
            raise HTTPException(409, "Saved query spectra and a current V16 flagship artifact are required; this endpoint does not rebuild them")
        before = {"query_artifact": _stamp(source), "flagship_artifact": _stamp(artifact)}
        payload = self.read_entity(artifact, entity_id)
        entity = next((x for x in payload.get("entities", []) if x.get("entity_id") == entity_id), None)
        if entity is None:
            raise HTTPException(404, "Entity is not in this task")
        raw_indices = entity.get("member_result_indices") or []
        if any(isinstance(i, bool) or not isinstance(i, int) or i < 0 for i in raw_indices):
            raise HTTPException(422, "Saved entity contains invalid scan indices")
        indices = sorted(set(raw_indices))
        if not indices or len(indices) > 200:
            raise HTTPException(422, "Comparison requires 1 to 200 saved entity spectra; no spectra are silently truncated")
        selected = {}
        with source.open("rb") as handle:
            for index, row in enumerate(ijson.items(handle, "item", use_float=True)):
                if index in indices:
                    selected[index] = row
                if index >= indices[-1]:
                    break
        if set(selected) != set(indices):
            raise HTTPException(409, "One or more saved query spectra are missing")
        queries = [project_query(selected[i]) for i in indices]
        native = entity.get("candidates") or []
        vetoed = {_connectivity(x) for x in native if x.get("vetoed")}
        initial = sorted({_connectivity(x) for x in native if not x.get("vetoed")} - {None})
        if len(initial) > 40:
            raise HTTPException(422, "Saved native pool exceeds the 40-candidate bound")
        winner = entity.get("winner")
        baseline = _connectivity(winner) if isinstance(winner, dict) else None
        retriever = self.retriever_factory()
        database_path = Path(retriever.db_path)
        before["database"] = _stamp(database_path)
        wal = Path(str(database_path) + "-wal")
        before["database_wal"] = _stamp(wal) if wal.exists() else None
        source_modules = (Path(__file__), Path(__file__).parents[1] / "pipeline/v16_candidate_comparison.py",
                          Path(__file__).parents[1] / "pipeline/v16_complete_candidate_panel.py")
        code_hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_modules}
        cache_key = semantic_sha256({"algorithm": ALGORITHM_ID, "inputs": before,
            "query_projection": queries, "entity": entity, "scoring_code": code_hashes,
            "library_snapshot": retriever.library_snapshot_id, "recall_policy": retriever.policy})
        if cache_key in self.cache:
            self.cache.move_to_end(cache_key)
            result = copy.deepcopy(self.cache[cache_key][0])
            result["cache_hit"] = True
            return result
        # Keep the native candidates and generate the same production Top40
        # expansion. The pool is fixed before any paired evidence is computed.
        recall = retriever.retrieve_entity(entity_id, queries,
            existing_candidates=[{"identity_key": "std2d:" + key} for key in initial])
        pool = sorted((set(initial) | {_connectivity(x) for x in recall.get("expanded_candidates", [])}) - {None} - vetoed)
        if not pool or len(pool) > 40:
            raise HTTPException(422, "No bounded candidate pool is available; candidates are not fabricated or truncated")
        panel = rescore_complete_panel(queries, pool, retriever, include_comparisons=True,
                                      baseline_connectivity=baseline)
        display_order = []
        for item in [*native, *recall.get("expanded_candidates", [])]:
            key = _connectivity(item)
            if key in pool and key not in display_order:
                display_order.append(key)
        after = {"query_artifact": _stamp(source), "flagship_artifact": _stamp(artifact),
                 "database": _stamp(database_path), "database_wal": _stamp(wal) if wal.exists() else None}
        if before != after or any(hashlib.sha256(p.read_bytes()).hexdigest() != code_hashes[str(p)] for p in source_modules):
            raise HTTPException(409, "Evidence, database, or comparison code changed during computation; result discarded")
        result = {"algorithm_id": ALGORITHM_ID, "task_id": task_id, "entity_id": entity_id,
            "comparison_id": cache_key,
            "runtime_snapshot_id": self.runtime_snapshot(), "panel": panel,
            "input_binding": before, "scoring_code_sha256": code_hashes,
            "frozen_candidate_pool_sha256": semantic_sha256(pool),
            "library_recall_policy_sha256": semantic_sha256(retriever.policy),
            "excluded_vetoed_connectivities": sorted(vetoed - {None}),
            "candidate_display_order": display_order,
            "display_order_rule": "saved_native_order_then_production_recall_order",
            "fragment_only_sorting_accepted_for_production": False,
            "source_result_indices": indices, "cache_hit": False,
            "production_results_modified": False, "primary_ranking_changed": False,
            "comparison_order_is_identification": False, "library_write_enabled": False}
        import json
        size = len(json.dumps(result, ensure_ascii=False, separators=(",", ":")).encode())
        if size <= self.maximum_cache_bytes // 2:
            self.cache[cache_key] = (copy.deepcopy(result), size)
            while len(self.cache) > 3 or sum(x[1] for x in self.cache.values()) > self.maximum_cache_bytes:
                self.cache.popitem(last=False)
        return result


def register_candidate_comparison_routes(app, service):
    async def compute(task_id, entity_id):
        try:
            return await asyncio.to_thread(service.compute, task_id, entity_id)
        except (ValueError, TypeError, KeyError) as exc:
            raise HTTPException(422, f"Saved evidence cannot be compared: {exc}") from exc
        except OSError as exc:
            raise HTTPException(409, "Saved evidence or the registered database is unavailable") from exc

    @app.post("/api/tasks/{task_id}/candidate-comparison")
    async def compare_task_candidates(task_id: str, request: CandidateComparisonRequest):
        return await compute(task_id, request.entity_id)

    @app.get("/api/tasks/{task_id}/candidate-comparison/export")
    async def export_task_comparison(task_id: str,
            entity_id: str = Query(min_length=1, max_length=256),
            comparison_id: str = Query(pattern=r"^[a-f0-9]{64}$")):
        result = await compute(task_id, entity_id)
        if result["comparison_id"] != comparison_id:
            raise HTTPException(409, "Inputs changed since the displayed comparison; recompute before exporting")
        return JSONResponse(result, headers={"Cache-Control": "no-store, max-age=0",
            "Content-Disposition": f'attachment; filename="{task_id}_v16_candidate_comparison.json"'})
