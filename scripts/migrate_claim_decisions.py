#!/usr/bin/env python3
"""Backfill compact EviGate 2.0 claim decisions for completed V9 tasks."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from governance.claim_compiler import default_task_claim_context
from pipeline.claim_contract import (
    claim_contract_path,
    iter_identification_records,
    write_claim_contract,
)
from pipeline.result_contract import identification_contract_path


def _existing_summary(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    metadata = dict(payload.get("metadata") or {})
    metadata.update(
        {
            "status": "available",
            "artifact": str(path.resolve()),
            "size_bytes": path.stat().st_size,
        }
    )
    return metadata


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
    failures = 0
    changed = False
    for task_id in selected:
        source = identification_contract_path(task_dir, task_id)
        output = claim_contract_path(task_dir, task_id)
        if output.is_file() and not args.force:
            summary = _existing_summary(output)
            print(f"{task_id}: exists ({output.name}); indexed")
        elif not source.is_file():
            failures += 1
            print(f"{task_id}: missing {source}", file=sys.stderr)
            continue
        else:
            summary = write_claim_contract(
                task_id,
                iter_identification_records(source),
                output,
                context=default_task_claim_context(task_id),
            )
            summary["status"] = "available"
            print(
                f"{task_id}: records={summary['record_count']}, "
                f"released={summary['released_count']}, "
                f"calibrated={summary['calibrated_count']}, "
                f"ranking_changed={summary['production_ranking_changed_records']}"
            )
        task = tasks_index.get(task_id)
        if isinstance(task, dict) and task.get("claim_contract_summary") != summary:
            task["claim_contract_summary"] = summary
            changed = True
    if changed:
        temporary = tasks_index_path.with_name(
            f".{tasks_index_path.name}.tmp-{os.getpid()}"
        )
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(tasks_index, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, tasks_index_path)
        print(f"updated {tasks_index_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
