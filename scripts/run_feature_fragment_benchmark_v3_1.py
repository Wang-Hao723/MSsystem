#!/usr/bin/env python3
"""Run the read-only V3.1 amendment over the frozen V3 corpus."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark_v3_1.workflow import run_benchmark_v3_1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--corpus",
        type=Path,
        default=ROOT
        / ".codex_tmp/feature_fragment_benchmark_v3/full_corpus/spectra.jsonl",
        help="Frozen V3 spectra JSONL; it is opened read-only.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "docs/audits/feature_fragment_benchmark_v3_1",
    )
    args = parser.parse_args()
    receipt = run_benchmark_v3_1(
        root=ROOT,
        frozen_corpus_path=args.corpus,
        output_dir=args.output_dir,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
