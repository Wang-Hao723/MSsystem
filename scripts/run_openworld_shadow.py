#!/usr/bin/env python3
"""Run fail-closed open-world model sidecars for an existing V10 task."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.openworld_shadow import (
    FioraRuntime,
    build_fiora_forward_shadow,
    write_shadow_artifact,
)


DEFAULT_FIORA_EXECUTABLE = Path("/opt/anaconda3/envs/fiora/bin/fiora-predict")
DEFAULT_FIORA_MODEL = Path(
    "/Users/zizi/Desktop/fiora/fiora/resources/models/fiora_OS_v1.0.0.pt"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run candidate-conditioned FIORA forward prediction as a separate "
            "shadow artifact; production results are never modified."
        )
    )
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fiora-executable", type=Path, default=DEFAULT_FIORA_EXECUTABLE)
    parser.add_argument("--fiora-model", type=Path, default=DEFAULT_FIORA_MODEL)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--max-candidates", type=int, default=5)
    parser.add_argument("--max-spectra", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_path = args.results or (
        ROOT / "output" / "tasks" / f"{args.task_id}_results.json"
    )
    output_path = args.output or (
        ROOT / "output" / "tasks" / f"{args.task_id}_openworld_shadow.json"
    )
    results = json.loads(results_path.read_text(encoding="utf-8"))
    if not isinstance(results, list):
        raise SystemExit("Results artifact must contain a JSON array")
    payload = build_fiora_forward_shadow(
        results,
        task_id=args.task_id,
        runtime=FioraRuntime(
            executable=args.fiora_executable,
            model_path=args.fiora_model,
            timeout_seconds=args.timeout,
            device=args.device,
        ),
        max_candidates_per_spectrum=args.max_candidates,
        max_spectra=args.max_spectra,
    )
    write_shadow_artifact(payload, output_path)
    print(json.dumps({
        "status": payload["status"],
        "output": str(output_path.resolve()),
        "summary": payload["summary"],
        "production_results_modified": False,
    }, ensure_ascii=False, indent=2))
    return 0 if payload["status"] in {"completed_shadow", "not_applicable"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
