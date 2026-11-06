#!/usr/bin/env python3
"""Validate V15 feature-fragment catalogs against the explicit public fixture."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragments.catalog import DEFAULT_CONFIG_DIR, load_catalog
from matching.feature_fragments.validation import file_sha256, validate_public_pilot


DEFAULT_FIXTURE = (
    ROOT
    / "tests"
    / "fixtures"
    / "feature_fragments"
    / "massbank_antwerp_phospholipids_v1.json"
)
DEFAULT_OUTPUT = (
    ROOT / "docs" / "audits" / "V15_FEATURE_FRAGMENT_PUBLIC_PILOT_RECEIPT_20260826.json"
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config-dir", type=Path, default=DEFAULT_CONFIG_DIR)
    parser.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    catalog = load_catalog(args.config_dir)
    with args.fixture.open("r", encoding="utf-8") as handle:
        fixture = json.load(handle)
    receipt = validate_public_pilot(catalog, fixture)
    receipt["generated_on"] = date.today().isoformat()
    receipt["inputs"] = {
        "source_catalog": {
            "path": str(args.config_dir / "source_catalog.json"),
            "sha256": file_sha256(args.config_dir / "source_catalog.json"),
        },
        "instrument_domains": {
            "path": str(args.config_dir / "instrument_domains.json"),
            "sha256": file_sha256(args.config_dir / "instrument_domains.json"),
        },
        "rules": {
            "path": str(args.config_dir / "rules_v1.json"),
            "sha256": file_sha256(args.config_dir / "rules_v1.json"),
        },
        "fixture": {
            "path": str(args.fixture),
            "sha256": file_sha256(args.fixture),
            "embedded_content_sha256": fixture["content_sha256"],
        },
    }
    receipt["implementation"] = {
        path: file_sha256(ROOT / path)
        for path in (
            "matching/feature_fragments/schema.py",
            "matching/feature_fragments/evidence.py",
            "matching/feature_fragments/validation.py",
        )
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    summary = {
        "output": str(args.output),
        "record_count": receipt["fixture_summary"]["record_count"],
        "identity_count": receipt["fixture_summary"]["independent_2d_identity_count"],
        "formula_family_count": receipt["fixture_summary"][
            "normalized_formula_family_count"
        ],
        "rule_count": receipt["catalog_summary"]["rule_count"],
        "instrument_source_audit_passed": receipt["instrument_family_source_audit"][
            "passed"
        ],
        "runtime_influence": receipt["runtime_influence"],
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
