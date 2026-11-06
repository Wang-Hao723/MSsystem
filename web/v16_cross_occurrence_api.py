"""Explicit-cohort API for the controlled V16 cross-occurrence ranker."""

from __future__ import annotations

import asyncio
import hashlib
import json
import os
import re
import shutil
import tempfile
import threading
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ijson
from fastapi import HTTPException
from pydantic import BaseModel, ConfigDict, Field

from pipeline.v16_cross_occurrence_ranker import (
    ALGORITHM_ID,
    apply_consensus_to_payload,
    build_consensus_plan,
    cluster_events,
    load_policy,
    project_entity,
    score_event_direct_neutral,
    semantic_sha256,
)


class CrossOccurrenceRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    task_ids: list[str]
    same_lc_method_confirmed: bool
    activate: bool = True
    cohort_label: str | None = Field(default=None, max_length=120)


class CrossOccurrenceRollbackRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    cohort_id: str = Field(pattern=r"^XCR-[0-9]{8}T[0-9]{6}Z-[a-f0-9]{16}$")
    confirm_rollback: bool


def _stamp(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "inode": stat.st_ino,
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(
                payload,
                handle,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
                default=str,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _read_metadata(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        return dict(next(ijson.items(handle, "metadata", use_float=True), {}) or {})


def _select_source_rows(path: Path, indices: set[int]) -> dict[int, dict[str, Any]]:
    if not indices:
        return {}
    if min(indices) < 0:
        raise ValueError("Saved result indices must be non-negative")
    selected: dict[int, dict[str, Any]] = {}
    last = max(indices)
    with path.open("rb") as handle:
        for index, row in enumerate(ijson.items(handle, "item", use_float=True)):
            if index in indices:
                selected[index] = row
            if index >= last:
                break
    if set(selected) != indices:
        raise ValueError("One or more saved query spectra are missing")
    return selected


class CrossOccurrenceService:
    """Build, atomically activate, and roll back one explicit cohort ranking."""

    def __init__(
        self,
        *,
        task_directory,
        task_lookup,
        artifact_path,
        artifact_current,
        retriever_factory,
        runtime_snapshot,
        runtime_idle,
        ensure_artifact=None,
    ):
        self.task_directory = task_directory
        self.task_lookup = task_lookup
        self.artifact_path = artifact_path
        self.artifact_current = artifact_current
        self.retriever_factory = retriever_factory
        self.runtime_snapshot = runtime_snapshot
        self.runtime_idle = runtime_idle
        self.ensure_artifact = ensure_artifact
        self.lock = threading.Lock()

    def _validated_tasks(self, request: CrossOccurrenceRequest, policy: dict[str, Any]):
        activation = policy["activation"]
        task_ids = [str(value or "").strip() for value in request.task_ids]
        if len(task_ids) != len(set(task_ids)):
            raise HTTPException(422, "Each task may appear only once in a cohort")
        if not int(activation["minimum_task_count"]) <= len(task_ids) <= int(activation["maximum_task_count"]):
            raise HTTPException(422, "Select 3 to 100 completed tasks from one LC method")
        if any(re.fullmatch(r"[A-Za-z0-9_-]{1,96}", task_id) is None for task_id in task_ids):
            raise HTTPException(422, "A task ID is invalid")
        if request.same_lc_method_confirmed is not True:
            raise HTTPException(422, "Explicit same-LC-method confirmation is required")
        tasks = {}
        for task_id in sorted(task_ids):
            task = self.task_lookup(task_id)
            if not task:
                raise HTTPException(404, f"Task not found: {task_id}")
            if task.get("status") != "completed":
                raise HTTPException(409, f"Task {task_id} is not completed")
            tasks[task_id] = task
        method_ids = {
            str((task.get("params") or {}).get("rt_method_id") or "").strip()
            for task in tasks.values()
        } - {""}
        if len(method_ids) > 1:
            raise HTTPException(422, "Selected tasks declare different LC method IDs")
        return tasks, sorted(method_ids)

    def compute(self, request: CrossOccurrenceRequest) -> dict[str, Any]:
        if not self.lock.acquire(blocking=False):
            raise HTTPException(429, "A cross-occurrence operation is already running")
        try:
            return self._compute(request)
        finally:
            self.lock.release()

    def _compute(self, request: CrossOccurrenceRequest) -> dict[str, Any]:
        policy = load_policy()
        tasks, method_ids = self._validated_tasks(request, policy)
        if not self.runtime_idle():
            raise HTTPException(409, "Wait until the analysis queue is idle before cohort reranking")
        directory = Path(self.task_directory())
        artifacts: dict[str, Path] = {}
        sources: dict[str, Path] = {}
        input_stamps: dict[str, dict[str, Any]] = {}
        input_metadata: dict[str, dict[str, Any]] = {}
        events: list[dict[str, Any]] = []
        for task_id in sorted(tasks):
            artifact = Path(self.artifact_path(directory, task_id))
            source = directory / f"{task_id}_results.json"
            if not source.is_file():
                raise HTTPException(409, f"Task {task_id} lacks current V16 ranking and saved spectra")
            if not artifact.is_file() or not self.artifact_current(artifact):
                if self.ensure_artifact is None:
                    raise HTTPException(409, f"Task {task_id} lacks current V16 ranking and saved spectra")
                if not self.runtime_idle():
                    raise HTTPException(409, "Analysis queue became active before V16 artifact materialization")
                self.ensure_artifact(task_id, tasks[task_id], source, artifact)
            if not artifact.is_file() or not self.artifact_current(artifact):
                raise HTTPException(409, f"Task {task_id} V16 ranking could not be materialized")
            metadata = _read_metadata(artifact)
            if metadata.get("cross_occurrence_consensus_status") == "activated":
                raise HTTPException(409, f"Task {task_id} already has an active cohort ranking; roll it back or reanalyse first")
            artifacts[task_id], sources[task_id] = artifact, source
            input_stamps[task_id] = {
                "flagship": _stamp(artifact),
                "source": _stamp(source),
            }
            input_metadata[task_id] = {
                "source_semantic_sha256": metadata.get("source_semantic_sha256"),
                "entity_semantic_sha256": metadata.get("entity_semantic_sha256"),
                "runtime_snapshot_id": metadata.get("runtime_snapshot_id"),
            }
            with artifact.open("rb") as handle:
                for entity in ijson.items(handle, "entities.item", use_float=True):
                    event = project_entity(task_id, entity)
                    if event is not None:
                        events.append(event)
        raw_clusters = cluster_events(events, policy)
        potential_ids: set[str] = set()
        events_by_id = {str(event["event_id"]): event for event in events}
        for cluster in raw_clusters:
            members = [events_by_id[str(row["feature_id"])] for row in cluster]
            formula_counts = Counter(str(member["baseline_formula"]) for member in members)
            potential_ids.update(
                str(member["event_id"])
                for member in members
                if formula_counts[str(member["baseline_formula"])] >= 3
            )
        wanted_by_task: dict[str, set[int]] = defaultdict(set)
        for event in events:
            if event["event_id"] not in potential_ids:
                event["direct_top1_connectivity"] = event["baseline_top1_connectivity"]
                event["candidate_scores"] = {}
                continue
            indices = list(event.get("member_result_indices") or [])
            if not indices or len(indices) > 200:
                raise HTTPException(422, "An aligned entity has an invalid or unbounded saved-spectrum set")
            wanted_by_task[str(event["task_id"])].update(indices)
        retriever = self.retriever_factory()
        database = Path(retriever.db_path)
        database_stamp = _stamp(database)
        scored_pairs = 0
        scored_events = 0
        for task_id, indices in sorted(wanted_by_task.items()):
            selected = _select_source_rows(sources[task_id], indices)
            for event in events:
                if event["task_id"] != task_id or event["event_id"] not in potential_ids:
                    continue
                rows = [selected[index] for index in event["member_result_indices"]]
                known_families = {
                    str(row.get("instrument_type") or "").strip().casefold()
                    for row in rows
                } - {""}
                if any("qtof" not in family for family in known_families):
                    raise HTTPException(422, "The validated ranker currently accepts QTOF acquisitions only")
                result = score_event_direct_neutral(event, rows, retriever, policy)
                event.update(result)
                scored_pairs += int(result["scored_pair_count"])
                scored_events += 1
        plan = build_consensus_plan(events, policy)
        plan.update(
            {
                "task_ids": sorted(tasks),
                "same_lc_method_confirmed": True,
                "declared_rt_method_ids": method_ids,
                "cohort_label": request.cohort_label,
                "potentially_vote_eligible_event_count": len(potential_ids),
                "spectrally_scored_event_count": scored_events,
                "scored_query_reference_pair_count": scored_pairs,
                "runtime_snapshot_id": self.runtime_snapshot(),
                "library_snapshot_id": retriever.library_snapshot_id,
                "policy_sha256": semantic_sha256(policy),
                "input_metadata": input_metadata,
                "activation_requested": bool(request.activate),
            }
        )
        current_stamps = {
            task_id: {"flagship": _stamp(artifacts[task_id]), "source": _stamp(sources[task_id])}
            for task_id in sorted(tasks)
        }
        if current_stamps != input_stamps or _stamp(database) != database_stamp:
            raise HTTPException(409, "An input artifact or the spectral library changed during scoring")
        if not request.activate:
            plan["status"] = "dry_run"
            plan["runtime_effect"] = False
            return plan
        if int(plan["applied_count"]) == 0:
            plan["status"] = "no_change"
            plan["runtime_effect"] = False
            return plan
        if not self.runtime_idle():
            raise HTTPException(409, "Analysis queue became active before cohort activation")
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        cohort_digest = semantic_sha256(
            {
                "tasks": sorted(tasks),
                "input": input_metadata,
                "plan": plan["proposals"],
                "policy": plan["policy_sha256"],
            }
        )[:16]
        cohort_id = f"XCR-{timestamp}-{cohort_digest}"
        ranking_root = directory / "cross_occurrence_rankings"
        backup_root = directory / "cross_occurrence_backups" / cohort_id
        if backup_root.exists() or (ranking_root / f"{cohort_id}.json").exists():
            raise HTTPException(409, "A cohort receipt with this identity already exists")
        required = 2 * sum(path.stat().st_size for path in artifacts.values()) + 100 * 1024 * 1024
        free = shutil.disk_usage(directory).free
        if free < required:
            raise HTTPException(507, "Insufficient free space for atomic backups and replacements")
        backup_root.mkdir(parents=True, exist_ok=False)
        replacements: dict[str, dict[str, Any]] = {}
        temporary_outputs: dict[str, Path] = {}
        replaced: list[str] = []
        try:
            for task_id in sorted(tasks):
                artifact = artifacts[task_id]
                original_sha = _file_sha256(artifact)
                backup = backup_root / artifact.name
                shutil.copy2(artifact, backup)
                if _file_sha256(backup) != original_sha:
                    raise OSError("Backup verification failed")
                payload = json.loads(artifact.read_text(encoding="utf-8"))
                updated, changed = apply_consensus_to_payload(
                    payload,
                    plan["proposals"],
                    cohort_id=cohort_id,
                    policy=policy,
                )
                temporary = artifact.with_name(f".{artifact.name}.{cohort_id}.tmp")
                _write_json_atomic(temporary, updated)
                activated_sha = _file_sha256(temporary)
                temporary_outputs[task_id] = temporary
                replacements[task_id] = {
                    "artifact": str(artifact.resolve()),
                    "backup": str(backup.resolve()),
                    "original_sha256": original_sha,
                    "activated_sha256": activated_sha,
                    "changed_entity_count": changed,
                }
            if {
                task_id: {"flagship": _stamp(artifacts[task_id]), "source": _stamp(sources[task_id])}
                for task_id in sorted(tasks)
            } != input_stamps:
                raise RuntimeError("Task artifacts changed before activation")
            for task_id in sorted(tasks):
                os.replace(temporary_outputs[task_id], artifacts[task_id])
                replaced.append(task_id)
            receipt = {
                **plan,
                "schema_version": "mssystem.v16.cross_occurrence_consensus.receipt/1.0.0",
                "status": "activated",
                "runtime_effect": True,
                "cohort_id": cohort_id,
                "activated_at": datetime.now(timezone.utc).isoformat(),
                "artifacts": replacements,
                "rollback_available": True,
                "compound_conclusion_sidecars_require_freshness_rebuild": True,
            }
            receipt["receipt_semantic_sha256"] = semantic_sha256(receipt)
            _write_json_atomic(ranking_root / f"{cohort_id}.json", receipt)
            return receipt
        except Exception:
            for task_id in replaced:
                backup = Path(replacements[task_id]["backup"])
                artifact = artifacts[task_id]
                restore = artifact.with_name(f".{artifact.name}.{cohort_id}.restore.tmp")
                shutil.copy2(backup, restore)
                os.replace(restore, artifact)
            raise
        finally:
            for temporary in temporary_outputs.values():
                temporary.unlink(missing_ok=True)

    def get_receipt(self, cohort_id: str) -> dict[str, Any]:
        if re.fullmatch(r"XCR-[0-9]{8}T[0-9]{6}Z-[a-f0-9]{16}", cohort_id) is None:
            raise HTTPException(404, "Cohort receipt not found")
        path = Path(self.task_directory()) / "cross_occurrence_rankings" / f"{cohort_id}.json"
        if not path.is_file():
            raise HTTPException(404, "Cohort receipt not found")
        return json.loads(path.read_text(encoding="utf-8"))

    def rollback(self, request: CrossOccurrenceRollbackRequest) -> dict[str, Any]:
        if request.confirm_rollback is not True:
            raise HTTPException(422, "Explicit rollback confirmation is required")
        if not self.lock.acquire(blocking=False):
            raise HTTPException(429, "A cross-occurrence operation is already running")
        try:
            if not self.runtime_idle():
                raise HTTPException(409, "Wait until the analysis queue is idle before rollback")
            receipt = self.get_receipt(request.cohort_id)
            if receipt.get("status") != "activated":
                raise HTTPException(409, "This cohort ranking is not active")
            restored = []
            for task_id, record in sorted((receipt.get("artifacts") or {}).items()):
                artifact = Path(str(record.get("artifact") or ""))
                backup = Path(str(record.get("backup") or ""))
                if not artifact.is_file() or not backup.is_file():
                    raise HTTPException(409, f"Rollback files are missing for task {task_id}")
                if _file_sha256(artifact) != record.get("activated_sha256"):
                    raise HTTPException(409, f"Task {task_id} changed after activation; rollback refused")
                if _file_sha256(backup) != record.get("original_sha256"):
                    raise HTTPException(409, f"Backup verification failed for task {task_id}")
            for task_id, record in sorted((receipt.get("artifacts") or {}).items()):
                artifact = Path(record["artifact"])
                restore = artifact.with_name(f".{artifact.name}.{request.cohort_id}.restore.tmp")
                shutil.copy2(Path(record["backup"]), restore)
                os.replace(restore, artifact)
                restored.append(task_id)
            receipt["status"] = "rolled_back"
            receipt["runtime_effect"] = False
            receipt["rollback_available"] = False
            receipt["rolled_back_at"] = datetime.now(timezone.utc).isoformat()
            receipt["restored_task_ids"] = restored
            receipt.pop("receipt_semantic_sha256", None)
            receipt["receipt_semantic_sha256"] = semantic_sha256(receipt)
            path = Path(self.task_directory()) / "cross_occurrence_rankings" / f"{request.cohort_id}.json"
            _write_json_atomic(path, receipt)
            return {
                "status": "rolled_back",
                "cohort_id": request.cohort_id,
                "restored_task_ids": restored,
            }
        finally:
            self.lock.release()


def register_cross_occurrence_routes(app, service: CrossOccurrenceService) -> None:
    @app.post("/api/workflows/cross-occurrence-ranking")
    async def run_cross_occurrence_ranking(request: CrossOccurrenceRequest):
        try:
            return await asyncio.to_thread(service.compute, request)
        except HTTPException:
            raise
        except (ValueError, TypeError, KeyError) as exc:
            raise HTTPException(422, f"Cross-occurrence evidence is invalid: {exc}") from exc
        except OSError as exc:
            raise HTTPException(409, f"Cross-occurrence artifact operation failed: {exc}") from exc

    @app.get("/api/workflows/cross-occurrence-ranking/{cohort_id}")
    async def get_cross_occurrence_ranking(cohort_id: str):
        return await asyncio.to_thread(service.get_receipt, cohort_id)

    @app.post("/api/workflows/cross-occurrence-ranking/rollback")
    async def rollback_cross_occurrence_ranking(request: CrossOccurrenceRollbackRequest):
        try:
            return await asyncio.to_thread(service.rollback, request)
        except HTTPException:
            raise
        except OSError as exc:
            raise HTTPException(409, f"Cross-occurrence rollback failed: {exc}") from exc


__all__ = [
    "CrossOccurrenceRequest",
    "CrossOccurrenceRollbackRequest",
    "CrossOccurrenceService",
    "register_cross_occurrence_routes",
]
