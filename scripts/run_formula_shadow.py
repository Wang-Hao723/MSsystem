#!/usr/bin/env python3
"""Run msFIDDLE QTOF formula prediction for an existing V10 task."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.openworld_shadow import (
    FiddleRuntime,
    build_fiddle_formula_shadow,
    write_shadow_artifact,
)


DEFAULT_FIDDLE_PYTHON = Path("/opt/anaconda3/envs/fiddle/bin/python")
DEFAULT_FIDDLE_WORKER = ROOT / "scripts" / "fiddle_worker.py"
DEFAULT_FIDDLE_MODEL = Path(
    "/Users/zizi/Desktop/Mssystem_model_assets/fiddle/fiddle_tcn_qtof.pt"
)
DEFAULT_FIDDLE_RESCORE_MODEL = Path(
    "/Users/zizi/Desktop/Mssystem_model_assets/fiddle/fiddle_rescore_qtof.pt"
)
DEFAULT_FIDDLE_CONFIG = Path(
    "/opt/anaconda3/envs/fiddle/lib/python3.10/site-packages/"
    "msfiddle/config/fiddle_tcn_qtof.yml"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run QTOF formula prediction as a separate shadow artifact; "
            "production results and High assignments are never modified."
        )
    )
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fiddle-python", type=Path, default=DEFAULT_FIDDLE_PYTHON)
    parser.add_argument("--fiddle-worker", type=Path, default=DEFAULT_FIDDLE_WORKER)
    parser.add_argument("--fiddle-model", type=Path, default=DEFAULT_FIDDLE_MODEL)
    parser.add_argument(
        "--fiddle-rescore-model", type=Path, default=DEFAULT_FIDDLE_RESCORE_MODEL
    )
    parser.add_argument("--fiddle-config", type=Path, default=DEFAULT_FIDDLE_CONFIG)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--max-spectra", type=int)
    parser.add_argument(
        "--lockmass-shadow",
        type=Path,
        help=(
            "Optional qualified candidate-independent lock-mass artifact used "
            "only for formula-model precursor input."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_path = args.results or (
        ROOT / "output" / "tasks" / f"{args.task_id}_results.json"
    )
    output_path = args.output or (
        ROOT / "output" / "tasks" / f"{args.task_id}_formula_shadow.json"
    )
    results = json.loads(results_path.read_text(encoding="utf-8"))
    if not isinstance(results, list):
        raise SystemExit("Results artifact must contain a JSON array")
    precursor_correction = None
    if args.lockmass_shadow:
        precursor_correction = json.loads(
            args.lockmass_shadow.read_text(encoding="utf-8")
        )
        if not isinstance(precursor_correction, dict):
            raise SystemExit("Lock-mass shadow artifact must contain a JSON object")
    payload = build_fiddle_formula_shadow(
        results,
        task_id=args.task_id,
        runtime=FiddleRuntime(
            python_executable=args.fiddle_python,
            worker_path=args.fiddle_worker,
            model_path=args.fiddle_model,
            rescore_model_path=args.fiddle_rescore_model,
            config_path=args.fiddle_config,
            timeout_seconds=args.timeout,
            device=args.device,
            batch_size=args.batch_size,
        ),
        top_k=args.top_k,
        max_spectra=args.max_spectra,
        precursor_correction=precursor_correction,
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
