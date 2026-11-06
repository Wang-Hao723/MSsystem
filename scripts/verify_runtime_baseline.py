#!/usr/bin/env python3
"""Verify frozen V10 task summaries and semantic result projections."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from runtime.golden import verify_golden_baseline


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        type=Path,
        default=PROJECT_ROOT / "config" / "golden_baseline.json",
    )
    parser.add_argument(
        "--tasks-index",
        type=Path,
        default=PROJECT_ROOT / "output" / "tasks_index.json",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=PROJECT_ROOT / "output" / "tasks",
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Verify task-level values without reading result artifacts.",
    )
    args = parser.parse_args()
    report = verify_golden_baseline(
        args.baseline,
        args.tasks_index,
        args.results_root,
        verify_semantic_digests=not args.summary_only,
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
