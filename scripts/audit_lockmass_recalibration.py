#!/usr/bin/env python3
"""Create a candidate-independent lock-mass shadow artifact from mzML."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.lockmass_recalibration import (
    AGILENT_NEGATIVE_REFERENCE_MZ,
    audit_mzml_lockmass,
    build_selected_match_counterfactual,
)
from matching.openworld_shadow import write_shadow_artifact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit an explicit vendor lock mass without using library matches; "
            "the resulting correction remains shadow-only."
        )
    )
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--mzml", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--results", type=Path)
    parser.add_argument(
        "--reference-mz", type=float, default=AGILENT_NEGATIVE_REFERENCE_MZ
    )
    parser.add_argument("--search-window-da", type=float, default=0.006)
    parser.add_argument("--minimum-intensity", type=float, default=1000.0)
    parser.add_argument("--minimum-observations", type=int, default=50)
    parser.add_argument("--minimum-presence-fraction", type=float, default=0.25)
    parser.add_argument("--maximum-mad-ppm", type=float, default=5.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = args.output or (
        ROOT / "output" / "tasks" / f"{args.task_id}_lockmass_shadow.json"
    )
    payload = audit_mzml_lockmass(
        args.mzml,
        task_id=args.task_id,
        reference_mz=args.reference_mz,
        search_window_da=args.search_window_da,
        minimum_intensity=args.minimum_intensity,
        minimum_observations=args.minimum_observations,
        minimum_presence_fraction=args.minimum_presence_fraction,
        maximum_mad_ppm=args.maximum_mad_ppm,
    )
    results_path = args.results or (
        ROOT / "output" / "tasks" / f"{args.task_id}_results.json"
    )
    if results_path.is_file():
        results = json.loads(results_path.read_text(encoding="utf-8"))
        if isinstance(results, list):
            payload["selected_match_counterfactual"] = (
                build_selected_match_counterfactual(results, payload)
            )
    write_shadow_artifact(payload, output)
    print(json.dumps({
        "status": payload["status"],
        "output": str(output.resolve()),
        "summary": payload.get("summary") or {},
        "correction": payload.get("correction"),
        "production_results_modified": False,
    }, ensure_ascii=False, indent=2))
    return 0 if payload["status"] in {"qualified_shadow", "not_qualified"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
