#!/usr/bin/env python3
"""Recompute V3 invariants and emit the phase-1 structured receipt."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark.corpus import file_sha256
from matching.feature_fragment_benchmark_v3.workflow import V2_ARTIFACT_NAMES


DB_IDENTIFIER_NAME = re.compile(
    r"^(?:NCGC[-_A-Z0-9]*|CCMSLIB[A-Z0-9_-]*|HMDB\d+|CHEBI:\d+|"
    r"CID:\d+|DTXSID\d+|MSBNK[-_A-Z0-9]+|LM[A-Z0-9]+)$",
    re.IGNORECASE,
)


def _json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _gzip_rows(path: Path) -> Iterable[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def audit(artifact_dir: Path) -> dict[str, Any]:
    benchmark = _json(artifact_dir / "benchmark_receipt.json")
    runtime_export = _json(artifact_dir / "runtime_export_gate.json")
    review = _json(artifact_dir / "review_candidates.json")
    tested = _json(artifact_dir / "sealed_test_rules.json")
    records = list(_gzip_rows(artifact_dir / "record_manifest.jsonl.gz"))
    runtime_gate_path = ROOT / "config/feature_fragments/runtime_gate_v2.json"
    runtime_gate = _json(runtime_gate_path)
    v2_dir = ROOT / "docs/audits/feature_fragment_benchmark_v2"
    current_v2 = {
        name: file_sha256(v2_dir / name)
        for name in V2_ARTIFACT_NAMES
        if (v2_dir / name).exists()
    }
    artifact_hash_checks = {
        name: {
            "expected": expected,
            "observed": file_sha256(artifact_dir / name),
            "matches": file_sha256(artifact_dir / name) == expected,
        }
        for name, expected in benchmark["artifact_sha256"].items()
    }
    database_identifier_names = [
        row["spectrum_id"]
        for row in records
        if DB_IDENTIFIER_NAME.fullmatch(str(row["compound"]["display_name"]))
    ]
    group_peak_hash_violations = sum(
        row["independence"].get("peak_hash_used_in_acquisition_replication_group")
        is not False
        for row in records
    )
    stage_violations = sum(
        row.get("evidence_stage") != "raw_msms_spectrum"
        or row.get("identification_scope", {}).get("standard_confirmed_by_mssystem")
        is not False
        for row in records
    )
    cross_unit_conversions = sum(
        row["acquisition"]["collision_energy"].get("conversion_applied") is not False
        for row in records
    )
    scan_probe = next(
        (
            row
            for row in records
            if row["source"]["accession"] == "MSBNK-Antwerp_Univ-AN111301"
        ),
        None,
    )
    scan_probe_ok = bool(
        scan_probe
        and scan_probe["acquisition"]["scan_range"]["minimum_mz"] == 55.0
        and scan_probe["acquisition"]["scan_range"]["maximum_mz"] == 1254.0
        and scan_probe["quantitative_eligibility"]["range_contains_precursor"]
        is True
    )
    rule_upgrade_violations = sum(
        bool(rule.get(field))
        for rule in tested
        for field in (
            "identity_upgrade_allowed",
            "formula_upgrade_allowed",
            "isomer_resolution_allowed",
            "standard_confirmation_allowed",
        )
    )
    runtime_fail_closed = bool(
        runtime_gate.get("deployment_state") == "fail_closed"
        and runtime_gate.get("runtime_evaluation_enabled") is False
        and runtime_gate.get("runtime_adapter_available") is False
        and runtime_gate.get("exported_rule_count") == 0
        and runtime_gate.get("runtime_influence") is False
        and runtime_gate.get("identity_rank_influence") is False
        and float(runtime_gate.get("score_contribution", 0.0)) == 0.0
    )
    config_paths = {
        "preregistration": ROOT
        / "config/feature_fragment_benchmark_v3/preregistration.json",
        "public_sources.json": ROOT
        / "config/feature_fragment_benchmark_v3/public_sources.json",
        "instrument_domains.json": ROOT
        / "config/feature_fragment_benchmark_v3/instrument_domains.json",
        "license_policy.json": ROOT
        / "config/feature_fragment_benchmark_v3/license_policy.json",
        "evidence_sources.json": ROOT
        / "config/feature_fragment_benchmark_v3/evidence_sources.json",
        "spectrum_metadata_schema.json": ROOT
        / "config/feature_fragment_benchmark_v3/spectrum_metadata_schema.json",
    }
    configuration_hash_checks = {
        name: file_sha256(path) for name, path in config_paths.items()
    }
    checks = {
        "artifact_hashes_match": all(
            row["matches"] for row in artifact_hash_checks.values()
        ),
        "v2_artifacts_unchanged": current_v2
        == benchmark["immutability"]["v2_artifact_hashes_before"]
        == benchmark["immutability"]["v2_artifact_hashes_after"],
        "runtime_gate_remains_fail_closed": runtime_fail_closed,
        "runtime_export_is_zero": runtime_export["runtime_exported_rule_count"]
        == review["runtime_exported_rule_count"]
        == benchmark["runtime"]["runtime_exported_rule_count"]
        == 0,
        "runtime_ranking_and_score_remain_off": runtime_export[
            "runtime_ranking_enabled"
        ]
        is False
        and float(runtime_export["score_contribution"]) == 0.0,
        "database_identifiers_not_used_as_compound_names": not database_identifier_names,
        "evidence_stage_and_standard_claims_valid": stage_violations == 0,
        "acquisition_group_excludes_peak_hash": group_peak_hash_violations == 0,
        "cross_ce_unit_conversion_is_zero": cross_unit_conversions == 0,
        "real_hyphenated_scan_range_is_fixed": scan_probe_ok,
        "rules_cannot_upgrade_identity_formula_isomer_or_standard": rule_upgrade_violations
        == 0,
        "next_state_waits_for_user_review": benchmark["next_state"]
        == review["next_state"]
        == runtime_export["next_state"]
        == "waiting_user_review",
        "config_hashes_match_benchmark": (
            configuration_hash_checks["preregistration"]
            == benchmark["preregistration"]["sha256"]
            and configuration_hash_checks["public_sources.json"]
            == benchmark["configuration_sha256"]["public_sources.json"]
            and configuration_hash_checks["instrument_domains.json"]
            == benchmark["configuration_sha256"]["instrument_domains.json"]
            and configuration_hash_checks["license_policy.json"]
            == benchmark["configuration_sha256"]["license_policy.json"]
            and configuration_hash_checks["evidence_sources.json"]
            == benchmark["configuration_sha256"]["evidence_sources.json"]
            and configuration_hash_checks["spectrum_metadata_schema.json"]
            == benchmark["configuration_sha256"]["spectrum_metadata_schema.json"]
        ),
    }
    failures = sorted(name for name, passed in checks.items() if not passed)
    grade_counts = Counter(rule.get("evidence_grade", "unavailable") for rule in tested)
    code_paths = [
        ROOT / "matching/feature_fragment_benchmark/corpus.py",
        *sorted((ROOT / "matching/feature_fragment_benchmark_v3").glob("*.py")),
        ROOT / "scripts/build_feature_fragment_benchmark_v3.py",
        ROOT / "scripts/run_feature_fragment_benchmark_v3.py",
        ROOT / "scripts/audit_feature_fragment_benchmark_v3.py",
    ]
    return {
        "schema_version": "3.0.0",
        "receipt_id": "v15-feature-fragment-benchmark-v3-phase1-20260827",
        "status": "pass" if not failures else "fail",
        "benchmark_state": benchmark["status"],
        "runtime_state": "fail_closed",
        "next_state": "waiting_user_review",
        "checks": checks,
        "failures": failures,
        "code_sha256": {
            str(path.relative_to(ROOT)): file_sha256(path) for path in code_paths
        },
        "configuration_sha256": configuration_hash_checks,
        "artifact_hash_checks": artifact_hash_checks,
        "v2_artifact_sha256": current_v2,
        "corpus": benchmark["inputs"],
        "public_spectrum_additions": benchmark["public_spectrum_additions"],
        "metadata": benchmark["metadata"],
        "mirrors": benchmark["mirror_audit"],
        "split": benchmark["split"],
        "collision_energy_models": benchmark["collision_energy_models"],
        "rules": {
            **benchmark["rules"],
            "recomputed_evidence_grade_counts": dict(sorted(grade_counts.items())),
        },
        "runtime": benchmark["runtime"],
        "deviations": benchmark["deviations"],
        "limitations": benchmark["limitations"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=ROOT / "docs/audits/feature_fragment_benchmark_v3",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_RECEIPT_20260827.json",
    )
    args = parser.parse_args()
    receipt = audit(args.artifact_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if receipt["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
