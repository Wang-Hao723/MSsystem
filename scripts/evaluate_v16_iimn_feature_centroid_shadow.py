#!/usr/bin/env python3
"""Evaluate the V16 feature-centred IIMN shadow for explicit task IDs.

The evaluator never derives polarity, sample identity, or pairing from a file
name.  Every task must be named explicitly, completed in the supplied task
index, and bound to both a saved result JSON and readable source MS data.
Output is written only after every requested task passes its source-hash and
shadow invariants.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import sys
import tempfile
import time
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.v16_iimn_feature_centroid_shadow import (  # noqa: E402
    build_v16_iimn_feature_centroid_shadow,
    load_v16_iimn_feature_centroid_shadow_policy,
    semantic_sha256,
)


EVALUATION_SCHEMA_VERSION = (
    "mssystem.v16.iimn-feature-centroid-shadow.evaluation/1.0.0"
)
EVALUATOR_ID = "mssystem-v16-iimn-feature-centroid-shadow-evaluator-1.0.0"


def _read_json(path: str | Path) -> Any:
    selected = Path(path).expanduser().resolve()
    try:
        return json.loads(selected.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Unreadable JSON input: {selected}: {exc}") from exc


def _explicit_task_ids(
    task_ids: Sequence[str] | None,
    task_id_file: str | Path | None,
) -> list[str]:
    selected = [str(value or "").strip() for value in task_ids or []]
    if task_id_file is not None:
        payload = _read_json(task_id_file)
        if isinstance(payload, Mapping):
            payload = payload.get("task_ids")
        if not isinstance(payload, list):
            raise ValueError("Task-ID file must be a JSON list or {task_ids: [...]}")
        selected.extend(str(value or "").strip() for value in payload)
    if not selected or any(not value for value in selected):
        raise ValueError("At least one non-empty explicit task ID is required")
    duplicates = sorted(
        task_id for task_id, count in Counter(selected).items() if count > 1
    )
    if duplicates:
        raise ValueError(f"Duplicate explicit task ID: {duplicates[0]}")
    return selected


def _stream_candidate_independent_rows(
    path: str | Path,
    *,
    allowed_fields: Sequence[str],
) -> list[dict[str, Any]]:
    """Stream only locked top-level scalar fields from a result array."""

    selected = Path(path).expanduser().resolve()
    if not selected.is_file():
        raise FileNotFoundError(f"Saved result JSON is unavailable: {selected}")
    try:
        import ijson
    except ImportError as exc:
        raise RuntimeError(
            "ijson is required for bounded-memory V16 shadow evaluation"
        ) from exc

    allowed = set(str(field) for field in allowed_fields)
    rows = []
    current: dict[str, Any] | None = None
    try:
        with selected.open("rb") as handle:
            for prefix, event, value in ijson.parse(handle, use_float=True):
                if prefix == "item" and event == "start_map":
                    current = {}
                    continue
                if prefix == "item" and event == "end_map":
                    if current is None:
                        raise RuntimeError("Malformed top-level result row")
                    rows.append(current)
                    current = None
                    continue
                if current is None or event not in {
                    "string",
                    "number",
                    "boolean",
                    "null",
                }:
                    continue
                if not prefix.startswith("item."):
                    continue
                field = prefix[5:]
                if "." not in field and field in allowed:
                    current[field] = value
    except (OSError, UnicodeError, ValueError) as exc:
        raise RuntimeError(f"Unreadable saved result JSON: {selected}: {exc}") from exc
    if current is not None:
        raise RuntimeError(f"Incomplete result row in {selected}")
    if not rows:
        raise RuntimeError(f"Saved result JSON contains no rows: {selected}")
    return rows


def _resolve_source_data(task_id: str, task: Mapping[str, Any]) -> Path:
    source = str(
        task.get("analysis_file_path")
        or task.get("raw_input_file_path")
        or task.get("file_path")
        or ""
    ).strip()
    if not source:
        raise RuntimeError(f"Task {task_id} has no source MS-data binding")
    path = Path(source).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Task {task_id} source MS data is unavailable: {path}")
    return path


def _task_record(
    task_id: str,
    task: Mapping[str, Any],
    result_path: Path,
    source_path: Path,
    rows: Sequence[Mapping[str, Any]],
    *,
    builder: Callable[..., dict[str, Any]],
    repeat_check: bool,
) -> dict[str, Any]:
    started = time.perf_counter()
    payload = builder(
        task,
        rows,
        source_results_path=result_path,
        source_data_path=source_path,
    )
    first_semantic = str(
        (payload.get("metadata") or {}).get("sidecar_semantic_sha256") or ""
    )
    if not first_semantic:
        raise RuntimeError(f"Task {task_id} shadow sidecar lacks a semantic hash")
    repeat_semantic = None
    if repeat_check:
        repeated = builder(
            task,
            rows,
            source_results_path=result_path,
            source_data_path=source_path,
        )
        repeat_semantic = str(
            (repeated.get("metadata") or {}).get("sidecar_semantic_sha256") or ""
        )
        if first_semantic != repeat_semantic:
            raise RuntimeError(f"Task {task_id} shadow repeat semantic hash changed")

    comparison = payload.get("comparison") or {}
    current = comparison.get("current_precursor_center") or {}
    shadow = comparison.get("feature_centroid_center") or {}
    delta = comparison.get("delta") or {}
    trace_audit = payload.get("trace_center_audit") or {}
    source_binding = payload.get("source_binding") or {}
    record = {
        "task_id": task_id,
        "polarity": str((payload.get("metadata") or {}).get("polarity") or ""),
        "candidate_independent_result_rows": len(rows),
        "trace_center_counts": copy.deepcopy(trace_audit.get("counts") or {}),
        "current_precursor_center": {
            "summary": copy.deepcopy(current.get("summary") or {}),
            "network_semantic_sha256": str(
                current.get("network_semantic_sha256") or ""
            ),
        },
        "linked_feature_centroid_center": {
            "summary": copy.deepcopy(shadow.get("summary") or {}),
            "network_semantic_sha256": str(
                shadow.get("network_semantic_sha256") or ""
            ),
        },
        "delta": {
            "network_count": int(delta.get("network_count") or 0),
            "networked_feature_count": int(
                delta.get("networked_feature_count") or 0
            ),
            "edge_count": int(delta.get("edge_count") or 0),
            "added_network_count": len(delta.get("added_network_fingerprints") or []),
            "removed_network_count": len(
                delta.get("removed_network_fingerprints") or []
            ),
        },
        "source_binding": {
            "task_metadata_sha256": str(
                (payload.get("metadata") or {}).get("task_metadata_sha256") or ""
            ),
            "candidate_independent_projection_sha256": str(
                (payload.get("metadata") or {}).get(
                    "candidate_independent_projection_sha256"
                )
                or ""
            ),
            "results_file": copy.deepcopy(source_binding.get("results_file")),
            "source_data_file": copy.deepcopy(
                source_binding.get("source_data_file")
            ),
            "current_trace_semantic_sha256": str(
                source_binding.get("current_trace_semantic_sha256") or ""
            ),
            "feature_centroid_trace_semantic_sha256": str(
                source_binding.get("feature_centroid_trace_semantic_sha256") or ""
            ),
            "source_files_unchanged": source_binding.get("source_files_unchanged")
            is True,
        },
        "sidecar_semantic_sha256": first_semantic,
        "repeatability": {
            "checked_in_process": bool(repeat_check),
            "first_semantic_sha256": first_semantic,
            "repeat_semantic_sha256": repeat_semantic,
            "consistent": (first_semantic == repeat_semantic)
            if repeat_check
            else None,
        },
        "runtime_seconds": round(time.perf_counter() - started, 6),
        "invariants": copy.deepcopy(payload.get("invariants") or {}),
    }
    if record["source_binding"]["source_files_unchanged"] is not True:
        raise RuntimeError(f"Task {task_id} source-file invariant failed")
    if record["polarity"] not in {"POS", "NEG"}:
        raise RuntimeError(f"Task {task_id} lacks an exact evaluated polarity")
    record["task_semantic_sha256"] = semantic_sha256(
        {key: value for key, value in record.items() if key != "runtime_seconds"}
    )
    return record


def _aggregate(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for polarity in ("ALL", "POS", "NEG"):
        selected = [
            record
            for record in records
            if polarity == "ALL" or record.get("polarity") == polarity
        ]
        if not selected:
            continue

        def summary(record: Mapping[str, Any], key: str) -> Mapping[str, Any]:
            return (record.get(key) or {}).get("summary") or {}

        counts = [record.get("trace_center_counts") or {} for record in selected]
        deltas = [record.get("delta") or {} for record in selected]
        output[polarity] = {
            "task_count": len(selected),
            "candidate_independent_result_rows": sum(
                int(record.get("candidate_independent_result_rows") or 0)
                for record in selected
            ),
            "rows_using_feature_centroid": sum(
                int(item.get("rows_using_feature_centroid") or 0) for item in counts
            ),
            "rows_using_precursor_fallback": sum(
                int(item.get("rows_using_precursor_fallback") or 0) for item in counts
            ),
            "feature_centroid_outside_current_xic_window": sum(
                int(item.get("feature_centroid_outside_current_xic_window") or 0)
                for item in counts
            ),
            "current_network_count": sum(
                int(summary(record, "current_precursor_center").get("network_count") or 0)
                for record in selected
            ),
            "feature_centroid_network_count": sum(
                int(
                    summary(record, "linked_feature_centroid_center").get(
                        "network_count"
                    )
                    or 0
                )
                for record in selected
            ),
            "network_count_delta": sum(
                int(delta.get("network_count") or 0) for delta in deltas
            ),
            "current_networked_feature_count": sum(
                int(
                    summary(record, "current_precursor_center").get(
                        "networked_features"
                    )
                    or 0
                )
                for record in selected
            ),
            "feature_centroid_networked_feature_count": sum(
                int(
                    summary(record, "linked_feature_centroid_center").get(
                        "networked_features"
                    )
                    or 0
                )
                for record in selected
            ),
            "current_edge_count": sum(
                int(summary(record, "current_precursor_center").get("edge_count") or 0)
                for record in selected
            ),
            "feature_centroid_edge_count": sum(
                int(
                    summary(record, "linked_feature_centroid_center").get(
                        "edge_count"
                    )
                    or 0
                )
                for record in selected
            ),
            "tasks_with_network_gain": sum(
                int(delta.get("network_count") or 0) > 0 for delta in deltas
            ),
            "tasks_with_no_network_change": sum(
                int(delta.get("network_count") or 0) == 0 for delta in deltas
            ),
            "tasks_with_network_loss": sum(
                int(delta.get("network_count") or 0) < 0 for delta in deltas
            ),
            "added_network_count": sum(
                int(delta.get("added_network_count") or 0) for delta in deltas
            ),
            "removed_network_count": sum(
                int(delta.get("removed_network_count") or 0) for delta in deltas
            ),
            "runtime_seconds": round(
                sum(float(record.get("runtime_seconds") or 0.0) for record in selected),
                6,
            ),
            "all_source_files_unchanged": all(
                (record.get("source_binding") or {}).get("source_files_unchanged")
                is True
                for record in selected
            ),
        }
    return output


def _semantic_material(payload: Mapping[str, Any]) -> dict[str, Any]:
    records = []
    for record in payload.get("tasks") or []:
        records.append(
            {key: value for key, value in record.items() if key != "runtime_seconds"}
        )
    aggregate = copy.deepcopy(payload.get("aggregate") or {})
    for summary in aggregate.values():
        if isinstance(summary, dict):
            summary.pop("runtime_seconds", None)
    return {
        "schema_version": payload.get("schema_version"),
        "evaluator_id": payload.get("evaluator_id"),
        "policy_sha256": payload.get("policy_sha256"),
        "requested_task_ids": payload.get("requested_task_ids"),
        "tasks": records,
        "aggregate": aggregate,
        "candidate_identity_used": False,
        "filename_used_for_polarity_or_pairing": False,
    }


def evaluate_tasks(
    *,
    task_ids: Sequence[str],
    tasks_index_path: str | Path,
    tasks_dir: str | Path,
    repeat_check: bool = False,
    builder: Callable[..., dict[str, Any]] = build_v16_iimn_feature_centroid_shadow,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Evaluate explicit task IDs and return one deterministic JSON payload."""

    explicit_ids = _explicit_task_ids(task_ids, None)
    index = _read_json(tasks_index_path)
    if not isinstance(index, Mapping):
        raise RuntimeError("Task index must be a JSON object keyed by task ID")
    task_root = Path(tasks_dir).expanduser().resolve()
    if not task_root.is_dir():
        raise FileNotFoundError(f"Task artifact directory is unavailable: {task_root}")
    policy = load_v16_iimn_feature_centroid_shadow_policy()
    allowed = policy["candidate_independent_projection"]["allowed_result_fields"]

    records = []
    started = time.perf_counter()
    for ordinal, task_id in enumerate(explicit_ids, start=1):
        task = index.get(task_id)
        if not isinstance(task, Mapping):
            raise KeyError(f"Explicit task ID is absent from the task index: {task_id}")
        if str(task.get("status") or "") != "completed":
            raise RuntimeError(f"Explicit task {task_id} is not completed")
        result_path = (task_root / f"{task_id}_results.json").resolve()
        if result_path.parent != task_root:
            raise RuntimeError(f"Task result path escaped the task directory: {task_id}")
        rows = _stream_candidate_independent_rows(
            result_path, allowed_fields=allowed
        )
        source_path = _resolve_source_data(task_id, task)
        records.append(
            _task_record(
                task_id,
                task,
                result_path,
                source_path,
                rows,
                builder=builder,
                repeat_check=repeat_check,
            )
        )
        if progress is not None:
            progress(
                f"evaluated {ordinal}/{len(explicit_ids)} task_id={task_id}"
            )

    payload = {
        "schema_version": EVALUATION_SCHEMA_VERSION,
        "evaluator_id": EVALUATOR_ID,
        "policy_sha256": semantic_sha256(policy),
        "requested_task_ids": list(explicit_ids),
        "task_count": len(records),
        "tasks": records,
        "aggregate": _aggregate(records),
        "runtime_seconds": round(time.perf_counter() - started, 6),
        "candidate_identity_used": False,
        "filename_used_for_polarity_or_pairing": False,
        "pairing_performed": False,
        "all_source_files_unchanged": all(
            (record.get("source_binding") or {}).get("source_files_unchanged")
            is True
            for record in records
        ),
    }
    payload["evaluation_semantic_sha256"] = semantic_sha256(
        _semantic_material(payload)
    )
    return payload


def _atomic_write_json(path: str | Path, payload: Mapping[str, Any]) -> Path:
    selected = Path(path).expanduser().resolve()
    selected.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{selected.name}.", suffix=".tmp", dir=selected.parent
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, selected)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
    return selected


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task-id",
        action="append",
        default=[],
        help="Explicit completed task ID; repeat for multiple tasks.",
    )
    parser.add_argument(
        "--task-id-file",
        help="JSON list, or object with task_ids, containing explicit task IDs.",
    )
    parser.add_argument(
        "--tasks-index",
        default=str(ROOT / "output" / "tasks_index.json"),
    )
    parser.add_argument(
        "--tasks-dir",
        default=str(ROOT / "output" / "tasks"),
    )
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--repeat-check",
        action="store_true",
        help="Run each task twice and require identical sidecar semantic hashes.",
    )
    parser.add_argument(
        "--compare-semantic-to",
        help="Require the completed evaluation hash to match an earlier JSON output.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    task_ids = _explicit_task_ids(args.task_id, args.task_id_file)
    payload = evaluate_tasks(
        task_ids=task_ids,
        tasks_index_path=args.tasks_index,
        tasks_dir=args.tasks_dir,
        repeat_check=bool(args.repeat_check),
        progress=lambda message: print(message, file=sys.stderr, flush=True),
    )
    if args.compare_semantic_to:
        previous = _read_json(args.compare_semantic_to)
        previous_hash = str(previous.get("evaluation_semantic_sha256") or "")
        current_hash = str(payload.get("evaluation_semantic_sha256") or "")
        if not previous_hash or previous_hash != current_hash:
            raise RuntimeError(
                "Evaluation semantic hash does not match the comparison artifact"
            )
        payload["repeatability_against_artifact"] = {
            "path": str(Path(args.compare_semantic_to).expanduser().resolve()),
            "semantic_sha256": current_hash,
            "consistent": True,
        }
    output = _atomic_write_json(args.output, payload)
    print(
        json.dumps(
            {
                "output": str(output),
                "task_count": payload["task_count"],
                "evaluation_semantic_sha256": payload[
                    "evaluation_semantic_sha256"
                ],
                "aggregate": payload["aggregate"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
