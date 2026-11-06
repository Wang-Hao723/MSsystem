#!/usr/bin/env python3
"""Run discovery, calibration, and sealed testing for the V15 benchmark."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark.corpus import file_sha256, iter_jsonl
from matching.feature_fragment_benchmark.workflow import run_benchmark


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus-jsonl", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--seed", default="v15-feature-fragment-benchmark-v2"
    )
    parser.add_argument("--minimum-target-identities", type=int, default=20)
    parser.add_argument("--minimum-background-identities", type=int, default=20)
    parser.add_argument("--minimum-target-present", type=int, default=5)
    parser.add_argument("--maximum-q-value", type=float, default=0.05)
    parser.add_argument("--minimum-lr-positive", type=float, default=2.0)
    parser.add_argument("--maximum-rules-per-class-stratum", type=int, default=12)
    args = parser.parse_args()
    records = list(iter_jsonl(args.corpus_jsonl))
    receipt = run_benchmark(
        records,
        output_dir=args.output_dir,
        corpus_jsonl_sha256=file_sha256(args.corpus_jsonl),
        seed=args.seed,
        discovery_parameters={
            "minimum_target_identities": args.minimum_target_identities,
            "minimum_background_identities": args.minimum_background_identities,
            "minimum_target_present": args.minimum_target_present,
            "maximum_q_value": args.maximum_q_value,
            "minimum_lr_positive": args.minimum_lr_positive,
            "maximum_rules_per_class_stratum": args.maximum_rules_per_class_stratum,
        },
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
