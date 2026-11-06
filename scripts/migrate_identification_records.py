#!/usr/bin/env python3
"""Backfill unified IdentificationRecord artifacts from completed legacy tasks."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pipeline.result_contract import (
    identification_contract_path,
    iter_legacy_result_artifact,
    runtime_snapshot_id_from_manifest,
    write_identification_contract,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_ids", nargs="*")
    parser.add_argument("--all-completed", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    task_dir = PROJECT_ROOT / "output" / "tasks"
    tasks_index_path = PROJECT_ROOT / "output" / "tasks_index.json"
    tasks_index = json.loads(tasks_index_path.read_text(encoding="utf-8"))
    selected = list(args.task_ids)
    if args.all_completed:
        selected.extend(
            task_id
            for task_id, task in tasks_index.items()
            if task.get("status") == "completed"
        )
    selected = list(dict.fromkeys(selected))
    if not selected:
        parser.error("provide task_ids or --all-completed")
    runtime_snapshot_id = runtime_snapshot_id_from_manifest(PROJECT_ROOT)
    failed = 0
    for task_id in selected:
        result_path = task_dir / f"{task_id}_results.json"
        output_path = identification_contract_path(task_dir, task_id)
        if output_path.exists() and not args.force:
            print(f"{task_id}: exists ({output_path.name})")
            continue
        if not result_path.is_file():
            failed += 1
            print(f"{task_id}: missing {result_path}", file=sys.stderr)
            continue
        metadata = write_identification_contract(
            task_id,
            iter_legacy_result_artifact(result_path),
            output_path,
            runtime_snapshot_id=runtime_snapshot_id,
        )
        print(
            f"{task_id}: {metadata['record_count']} records, "
            f"ranking_changed={metadata['production_ranking_changed_records']}, "
            f"sha256={metadata['record_semantic_sha256']}"
        )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
