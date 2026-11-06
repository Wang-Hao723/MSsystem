#!/usr/bin/env python3
"""Run the V3 preregistered research benchmark and write audit artifacts."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark.corpus import iter_jsonl
from matching.feature_fragment_benchmark_v3.workflow import run_benchmark_v3


def _v2_ids(path: Path) -> set[str]:
    values = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                values.add(str(json.loads(line)["spectrum_id"]))
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus-dir", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "docs/audits/feature_fragment_benchmark_v3",
    )
    parser.add_argument(
        "--preregistration",
        type=Path,
        default=ROOT / "config/feature_fragment_benchmark_v3/preregistration.json",
    )
    parser.add_argument(
        "--public-sources",
        type=Path,
        default=ROOT / "config/feature_fragment_benchmark_v3/public_sources.json",
    )
    parser.add_argument(
        "--instrument-domains",
        type=Path,
        default=ROOT / "config/feature_fragment_benchmark_v3/instrument_domains.json",
    )
    parser.add_argument(
        "--license-policy",
        type=Path,
        default=ROOT / "config/feature_fragment_benchmark_v3/license_policy.json",
    )
    parser.add_argument(
        "--evidence-sources",
        type=Path,
        default=ROOT / "config/feature_fragment_benchmark_v3/evidence_sources.json",
    )
    parser.add_argument(
        "--metadata-schema",
        type=Path,
        default=ROOT
        / "config/feature_fragment_benchmark_v3/spectrum_metadata_schema.json",
    )
    args = parser.parse_args()
    corpus_receipt = json.loads(
        (args.corpus_dir / "corpus_receipt.json").read_text(encoding="utf-8")
    )
    records = list(iter_jsonl(args.corpus_dir / "spectra.jsonl"))
    v2_dir = ROOT / "docs/audits/feature_fragment_benchmark_v2"
    receipt = run_benchmark_v3(
        records,
        output_dir=args.output_dir,
        corpus_receipt=corpus_receipt,
        rejects_jsonl=args.corpus_dir / "rejects.jsonl",
        preregistration_path=args.preregistration,
        public_sources_path=args.public_sources,
        instrument_domains_path=args.instrument_domains,
        license_policy_path=args.license_policy,
        evidence_sources_path=args.evidence_sources,
        metadata_schema_path=args.metadata_schema,
        runtime_gate_path=ROOT / "config/feature_fragments/runtime_gate_v2.json",
        v2_artifact_dir=v2_dir,
        v2_record_ids=_v2_ids(v2_dir / "record_manifest.jsonl.gz"),
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
