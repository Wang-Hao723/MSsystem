#!/usr/bin/env python3
"""Run the preregistered EvoGraph Phase 0 anti-leakage audit."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from eval.evograph.phase0_audit import (  # noqa: E402
    audit_phase0_manifests,
    load_manifest,
    write_audit_report,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Audit identity/scaffold/spectrum leakage and peak-blind negative "
            "controls before EvoGraph training or evaluation."
        )
    )
    parser.add_argument("--train", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--allow-scaffold-overlap",
        action="store_true",
        help="Report scaffold overlap without making it a gate failure.",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Always exit zero after writing the report (gate status is unchanged).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    train = load_manifest(args.train)
    test = load_manifest(args.test)
    report = audit_phase0_manifests(
        train,
        test,
        train_path=args.train,
        test_path=args.test,
        strict_scaffold=not args.allow_scaffold_overlap,
    )
    write_audit_report(report, args.output)
    print(json.dumps({
        "output": str(args.output.resolve()),
        "passed": report["gate"]["passed"],
        "blocking_violations": report["gate"]["blocking_violations"],
        "train_rows": report["train"]["row_count"],
        "test_rows": report["test"]["row_count"],
    }, ensure_ascii=False))
    if args.report_only or report["gate"]["passed"]:
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
