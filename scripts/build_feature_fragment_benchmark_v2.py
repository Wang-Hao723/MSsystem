#!/usr/bin/env python3
"""Build the V15 public MassBank feature-fragment corpus."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark.corpus import build_massbank_corpus


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT
        / "config"
        / "feature_fragment_benchmark_v2"
        / "public_massbank_sources.json",
    )
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--download-missing", action="store_true")
    parser.add_argument("--max-records", type=int)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    receipt = build_massbank_corpus(
        args.manifest,
        cache_dir=args.cache_dir,
        output_jsonl=args.output_dir / "spectra.jsonl",
        rejects_jsonl=args.output_dir / "rejects.jsonl",
        download_missing=args.download_missing,
        max_records=args.max_records,
    )
    receipt_path = args.output_dir / "corpus_receipt.json"
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
