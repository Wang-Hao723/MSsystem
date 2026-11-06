#!/usr/bin/env python3
"""Map msFIDDLE formulas to read-only local-library structure hypotheses."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config.settings import DB_PATH
from matching.formula_candidate_shadow import (
    build_formula_conditioned_library_shadow,
)
from matching.openworld_shadow import write_shadow_artifact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Enumerate formula-conditioned structure hypotheses in a separate "
            "shadow artifact without changing production ranking."
        )
    )
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--formula-shadow", type=Path)
    parser.add_argument("--database", type=Path, default=DB_PATH)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-candidates-per-formula", type=int, default=10)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    formula_path = args.formula_shadow or (
        ROOT / "output" / "tasks" / f"{args.task_id}_formula_shadow.json"
    )
    output_path = args.output or (
        ROOT / "output" / "tasks" / f"{args.task_id}_formula_structure_shadow.json"
    )
    formula_artifact = json.loads(formula_path.read_text(encoding="utf-8"))
    if not isinstance(formula_artifact, dict):
        raise SystemExit("Formula shadow artifact must contain a JSON object")
    payload = build_formula_conditioned_library_shadow(
        formula_artifact,
        db_path=args.database,
        max_candidates_per_formula=args.max_candidates_per_formula,
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
