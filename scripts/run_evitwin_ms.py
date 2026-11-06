#!/usr/bin/env python3
"""Materialize EviTwin-MS authentic-reference ranking sidecars.

The runner never reads the validation workbook.  Tasks can be selected by
explicit ID or filename regular expression.  A standard injection cannot
match the authentic entry built from that same task.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from governance.authentic_anchor_registry import (  # noqa: E402
    load_active_authentic_anchor_library,
    resolve_active_library_path,
)
from matching.evitwin_ms import (  # noqa: E402
    build_evitwin_result,
    write_evitwin_result,
)
from matching.evitwin_runtime import (  # noqa: E402
    resolve_evitwin_query_instrument_family,
)
from pipeline.v15_flagship import v15_flagship_artifact_path  # noqa: E402


DEFAULT_REGISTRY = PROJECT_ROOT / "config" / "authentic_anchor_registry.json"
DEFAULT_TASK_DIR = PROJECT_ROOT / "output" / "tasks"
DEFAULT_TASK_INDEX = PROJECT_ROOT / "output" / "tasks_index.json"


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _active_library(registry_path: Path) -> tuple[Path, dict[str, Any]]:
    _, registry = load_active_authentic_anchor_library(
        registry_path,
        project_root=PROJECT_ROOT,
    )
    return resolve_active_library_path(
        registry,
        project_root=PROJECT_ROOT,
    ), registry


def _task_polarity(task: dict[str, Any]) -> str:
    filename = str(task.get("filename") or "")
    match = re.search(r"(?:^|_)(POS|NEG)(?:_|\.|$)", filename, re.I)
    if match:
        return match.group(1).upper()
    declared = str(task.get("ion_mode") or "").upper()
    if declared in {"POS", "NEG"}:
        return declared
    counts = task.get("ion_mode_counts") or {}
    nonzero = [
        mode for mode in ("POS", "NEG") if int(counts.get(mode) or 0) > 0
    ]
    return nonzero[0] if len(nonzero) == 1 else ""


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-id", action="append", default=[])
    parser.add_argument(
        "--filename-regex",
        help="Select completed tasks whose saved filename matches this regex",
    )
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--task-dir", type=Path, default=DEFAULT_TASK_DIR)
    parser.add_argument(
        "--primary-dir",
        "--v11-dir",
        dest="primary_dir",
        type=Path,
        help=(
            "Optional directory containing active-version primary flagship "
            "artifacts; --v11-dir is retained as a deprecated alias."
        ),
    )
    parser.add_argument(
        "--primary-map",
        "--v11-map",
        dest="primary_map",
        type=Path,
        help=(
            "Optional JSON object mapping task IDs to explicit primary result paths; "
            "useful when replay artifacts have cohort-specific names."
        ),
    )
    parser.add_argument("--task-index", type=Path, default=DEFAULT_TASK_INDEX)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_TASK_DIR,
        help="Sidecar destination; defaults to the live task artifact directory",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    library_path, registry = _active_library(args.registry.resolve())
    library, _ = load_active_authentic_anchor_library(
        args.registry.resolve(),
        project_root=PROJECT_ROOT,
    )
    tasks = _load_json(args.task_index)
    if not isinstance(tasks, dict):
        raise ValueError("task index must be a JSON object")
    wanted = {str(value).strip() for value in args.task_id if str(value).strip()}
    pattern = re.compile(args.filename_regex, re.I) if args.filename_regex else None
    if not wanted and pattern is None:
        raise ValueError("select tasks with --task-id or --filename-regex")
    selected = []
    for task_id, task in tasks.items():
        if not isinstance(task, dict) or task.get("status") != "completed":
            continue
        if task_id in wanted or (
            pattern is not None and pattern.search(str(task.get("filename") or ""))
        ):
            selected.append((task_id, task))
    unknown = wanted - {task_id for task_id, _ in selected}
    if unknown:
        raise ValueError(f"requested completed task IDs were not found: {sorted(unknown)}")
    records = []
    primary_dir = args.primary_dir or args.task_dir
    primary_map = (
        _load_json(args.primary_map) if args.primary_map is not None else {}
    )
    if not isinstance(primary_map, dict):
        raise ValueError("--primary-map must contain a JSON object")
    for task_id, task in sorted(selected):
        polarity = _task_polarity(task)
        if polarity not in {"POS", "NEG"}:
            records.append(
                {
                    "task_id": task_id,
                    "filename": task.get("filename"),
                    "status": "skipped_unknown_polarity",
                }
            )
            continue
        result_path = args.task_dir / f"{task_id}_results.json"
        declared_primary = str(primary_map.get(task_id) or "").strip()
        if declared_primary:
            declared_path = Path(declared_primary).expanduser()
            primary_path = (
                declared_path.resolve()
                if declared_path.is_absolute()
                else (PROJECT_ROOT / declared_path).resolve()
            )
        else:
            primary_path = v15_flagship_artifact_path(primary_dir, task_id)
        destination = args.output_dir / f"{task_id}_evitwin_ms.json"
        if destination.is_file() and not args.overwrite:
            records.append(
                {
                    "task_id": task_id,
                    "filename": task.get("filename"),
                    "status": "retained_existing",
                    "artifact": str(destination.resolve()),
                }
            )
            continue
        if not result_path.is_file() or not primary_path.is_file():
            records.append(
                {
                    "task_id": task_id,
                    "filename": task.get("filename"),
                    "status": "skipped_missing_source_artifact",
                }
            )
            continue
        rows = _load_json(result_path)
        primary = _load_json(primary_path)
        if not isinstance(rows, list) or not isinstance(primary, dict):
            raise ValueError(f"{task_id}: source artifacts have invalid shapes")
        instrument_audit = resolve_evitwin_query_instrument_family(
            task,
            reference_instrument_family=str(
                registry.get("instrument_family") or "agilent-qtof-6546"
            ),
        )
        payload = build_evitwin_result(
            primary,
            rows,
            library,
            task_id=task_id,
            polarity=polarity,
            query_instrument_family=str(
                instrument_audit["scoring_instrument_family"]
            ),
            evirank_payload=None,
        )
        payload["metadata"]["instrument_applicability"] = instrument_audit
        write_evitwin_result(destination, payload)
        records.append(
            {
                "task_id": task_id,
                "filename": task.get("filename"),
                "status": "materialized",
                "artifact": str(destination.resolve()),
                **payload["summary"],
            }
        )
    summary = {
        "selected_task_count": len(selected),
        "materialized_task_count": sum(
            record["status"] == "materialized" for record in records
        ),
        "priority_compound_count": sum(
            int(record.get("priority_compound_count") or 0) for record in records
        ),
        "display_rank_change_count": sum(
            int(record.get("display_rank_change_count") or 0) for record in records
        ),
        "library": str(library_path),
        "records": records,
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
