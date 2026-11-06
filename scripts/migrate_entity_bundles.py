#!/usr/bin/env python3
"""Backfill candidate-independent entity bundles for completed V9 tasks."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pipeline.entity_builder import entity_bundle_path, write_entity_bundle
from pipeline.result_contract import iter_legacy_result_artifact


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
    index_changed = False
    for task_id in selected:
        result_path = task_dir / f"{task_id}_results.json"
        output_path = entity_bundle_path(task_dir, task_id)
        if output_path.exists() and not args.force:
            payload = json.loads(output_path.read_text(encoding="utf-8"))
            summary = dict(payload.get("summary") or {})
            summary.update(
                {
                    "status": "available",
                    "schema_version": payload.get("schema_version"),
                    "graph_sha256": (payload.get("entity_graph") or {}).get(
                        "graph_sha256"
                    ),
                    "artifact": str(output_path.resolve()),
                    "size_bytes": output_path.stat().st_size,
                }
            )
            print(f"{task_id}: exists ({output_path.name}); indexed")
        else:
            if not result_path.is_file():
                failures += 1
                print(f"{task_id}: missing {result_path}", file=sys.stderr)
                continue
            summary = write_entity_bundle(
                task_id,
                iter_legacy_result_artifact(result_path),
                output_path,
            )
            print(
                f"{task_id}: scans={summary['scan_count']}, "
                f"events={summary['precursor_event_count']}, "
                f"features={summary['chromatographic_feature_count']}, "
                f"ion_groups={summary['ion_identity_group_count']}, "
                f"ranking_changed={summary['feature_consensus_ranking_changed']}"
            )
        task = tasks_index.get(task_id)
        if isinstance(task, dict) and task.get("entity_bundle_summary") != summary:
            task["entity_bundle_summary"] = summary
            index_changed = True
    if index_changed:
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
