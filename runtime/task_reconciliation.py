"""Restore a corrupted status projection only from an intact completed run."""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

import ijson


def task_reconciliation_fingerprint(task: dict[str, Any]) -> str:
    payload = {key: value for key, value in task.items() if key != "results"}
    raw = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(raw.encode()).hexdigest()


def _completed_ledger(task: dict[str, Any]) -> dict[str, Any]:
    history = task.get("state_history") or []
    if task.get("status") != "failed" or not history:
        raise ValueError("Only a failed status with a completed ledger can be reconciled")
    last = history[-1]
    run_id = task.get("active_run_id")
    if not (
        last.get("to") == "completed"
        and last.get("reason") == "v11_flagship_and_supporting_artifacts_committed"
        and last.get("version") == task.get("state_version")
        and run_id and last.get("run_id") == run_id
        and task.get("result_run_id") == run_id
    ):
        raise ValueError("Completed ledger, result run and current generation do not agree")
    return last


def _stamp(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "inode": stat.st_ino}


def _content_binding(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
    return {**_stamp(path), "sha256": sha256}


def _metadata(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        result = next(ijson.items(handle, "metadata", use_float=True), None)
    if not isinstance(result, dict):
        raise ValueError(f"Missing completed artifact metadata: {path.name}")
    return result


def verify_completed_task_artifacts(task: dict[str, Any], task_directory: Path) -> dict[str, Any]:
    """Read explicit committed artifacts; never recreate results or alter a task."""
    last = _completed_ledger(task)
    directory = Path(task_directory).resolve()
    task_id = str(task["id"])
    source = directory / f"{task_id}_results.json"
    compound_summary = task.get("compound_conclusion_summary") or {}
    flagship_summary = task.get("flagship_summary") or {}
    paths = [source]
    for summary in (compound_summary, flagship_summary):
        value = summary.get("artifact")
        if not value:
            raise ValueError("Completed artifact path is not registered")
        path = Path(value).resolve()
        if not path.is_relative_to(directory) or not path.name.startswith(task_id + "_"):
            raise ValueError("Completed artifact is outside this task directory")
        paths.append(path)
    before = [_stamp(path) for path in paths]
    compound = _metadata(paths[1])
    flagship = _metadata(paths[2])
    if any(meta.get("task_id") != task_id for meta in (compound, flagship)):
        raise ValueError("Completed artifact belongs to another task")
    if Path(compound.get("source_artifact") or "").resolve() != source:
        raise ValueError("Compound source binding differs from the saved scan file")
    for field, actual in (("source_artifact_size_bytes", before[0]["size"]),
                          ("source_artifact_mtime_ns", before[0]["mtime_ns"]),
                          ("primary_identity_artifact_size_bytes", before[2]["size"]),
                          ("primary_identity_artifact_mtime_ns", before[2]["mtime_ns"])):
        if compound.get(field) != actual:
            raise ValueError(f"Completed artifact binding changed: {field}")
    if Path(compound.get("primary_identity_artifact") or "").resolve() != paths[2]:
        raise ValueError("Completed primary identity binding differs")
    for key in ("source_row_count", "source_semantic_sha256"):
        if compound.get(key) != compound_summary.get(key):
            raise ValueError(f"Compound completion receipt differs: {key}")
        if flagship.get(key) != flagship_summary.get(key):
            raise ValueError(f"Flagship completion receipt differs: {key}")
    count = 0
    with source.open("rb") as handle:
        # A saved analysis is a JSON array, not an arbitrary object with an
        # `item` member. Check the root before streaming its potentially large rows.
        while True:
            first = handle.read(1)
            if not first or not first.isspace():
                break
        if first != b"[":
            raise ValueError("Saved scan results are not an array")
        handle.seek(0)
        for row in ijson.items(handle, "item", use_float=True):
            if not isinstance(row, dict) or row.get("spectrum_id") is None:
                raise ValueError("Saved scan result is malformed")
            count += 1
    expected = task.get("total_spectra")
    if not count or any(value != count for value in
                        (expected, compound.get("source_row_count"), flagship.get("source_row_count"))):
        raise ValueError("Completed scan count does not match the retained artifacts")
    if before != [_stamp(path) for path in paths]:
        raise ValueError("Completed artifacts changed during verification")
    bindings = [_content_binding(path) for path in paths]
    if before != [
        {key: value for key, value in binding.items() if key != "sha256"}
        for binding in bindings
    ]:
        raise ValueError("Completed artifacts changed during content binding")
    return {"task_id": task_id, "run_id": last["run_id"], "state_version": last["version"],
            "task_fingerprint": task_reconciliation_fingerprint(task),
            "source_row_count": count, "artifact_stamps": bindings,
            "verification": "completed_run_and_retained_source_bindings",
            "historical_source_content_hash_available": False,
            "current_artifact_content_sha256_bound": True,
            "score_rank_or_evidence_changes": False}


def reconciled_completed_task(task: dict[str, Any], evidence: dict[str, Any]) -> dict[str, Any]:
    last = _completed_ledger(task)
    if evidence.get("task_fingerprint") != task_reconciliation_fingerprint(task):
        raise ValueError("Task changed during verification")
    if not evidence.get("artifact_stamps"):
        raise ValueError("Missing verified artifact bindings")
    for expected in evidence["artifact_stamps"]:
        if _content_binding(Path(expected["path"])) != expected:
            raise ValueError("Artifact changed before status reconciliation")
    result = deepcopy({key: value for key, value in task.items() if key != "results"})
    result["results"] = task.get("results", [])
    now = datetime.now(timezone.utc).isoformat()
    version = int(task["state_version"]) + 1
    event = {"version": version, "from": task["status"], "to": "completed", "at": now,
             "reason": "reconciled_completed_ledger_after_unrecorded_status_override",
             "run_id": last["run_id"], "prior_error": task.get("error"),
             "prior_message": task.get("message"), "verification": evidence}
    result.update(status="completed", state_version=version, progress=100, error=None,
                  end_time=last["at"], message="Completed results restored from the verified run ledger")
    result["state_history"] = list(task["state_history"]) + [event]
    result["state_reconciliations"] = list(task.get("state_reconciliations") or []) + [event]
    return result
