#!/usr/bin/env python3
"""Run the V3.1 full-manifest readiness and intake audit."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark_v31.common import iter_gzip_jsonl
from matching.feature_fragment_benchmark_v31.workflow import run_v31_readiness


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--v3-manifest",
        type=Path,
        default=ROOT
        / "docs/audits/feature_fragment_benchmark_v3/record_manifest.jsonl.gz",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "docs/audits/feature_fragment_benchmark_v31",
    )
    args = parser.parse_args()
    config = ROOT / "config/feature_fragment_benchmark_v31"
    receipt = run_v31_readiness(
        list(iter_gzip_jsonl(args.v3_manifest)),
        source_manifest_path=args.v3_manifest,
        output_dir=args.output_dir,
        preregistration_path=config / "preregistration.json",
        metadata_contract_path=config / "metadata_derivation_contract.json",
        source_screening_contract_path=config / "source_screening_contract.json",
        chemistry_contract_path=config / "chemistry_hypothesis_contract.json",
        class_contract_path=config / "class_hierarchy_contract.json",
        rule_schema_path=config / "rule_schema.json",
        runtime_gate_path=ROOT / "config/feature_fragments/runtime_gate_v2.json",
        v2_artifact_dir=ROOT / "docs/audits/feature_fragment_benchmark_v2",
        v3_artifact_dir=ROOT / "docs/audits/feature_fragment_benchmark_v3",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
