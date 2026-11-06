#!/usr/bin/env python3
"""Isolated worker for one task workbook; keeps heavy rendering out of the API GIL."""

from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    from web import server

    task_id = str(args.task_id)
    if not task_id or any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for character in task_id):
        raise ValueError("Invalid task identifier")
    index = json.loads(server.TASKS_FILE.read_text(encoding="utf-8"))
    task = index.get(task_id)
    if not isinstance(task, dict):
        raise KeyError(f"Task not found: {task_id}")
    result_path = server.TASK_DIR / f"{task_id}_results.json"
    if not result_path.is_file():
        raise FileNotFoundError(result_path)
    started = time.perf_counter()
    cpu_started = time.process_time()
    results = json.loads(result_path.read_text(encoding="utf-8"))
    read_finished = time.perf_counter()
    if not isinstance(results, list):
        raise ValueError("Saved task result is not a JSON array")
    prepared = server._project_claims_for_task(task_id, results)
    prepared_finished = time.perf_counter()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.tmp-{os.getpid()}")
    projection_started = time.perf_counter()
    payload = server._build_analysis_excel_bytes(
        task_id,
        task,
        prepared,
        language="zh",
        output_path=temporary,
    )
    os.replace(temporary, args.output)
    finished = time.perf_counter()
    usage = resource.getrusage(resource.RUSAGE_SELF)
    print(
        json.dumps(
            {
                "task_id": task_id,
                "output": str(args.output),
                "size_bytes": args.output.stat().st_size,
                "runtime_snapshot_id": server._current_runtime_snapshot_id(),
                "timings": {
                    "read_wall_seconds": read_finished - started,
                    "projection_wall_seconds": prepared_finished - read_finished,
                    "projection_and_render_wall_seconds": finished - projection_started,
                    "total_wall_seconds": finished - started,
                    "process_cpu_seconds": time.process_time() - cpu_started,
                    "max_rss_kib": int(usage.ru_maxrss // 1024
                                       if sys.platform == "darwin"
                                       else usage.ru_maxrss),
                },
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
