#!/usr/bin/env python3
"""Recompute V3.1 full-manifest invariants and emit a structured receipt."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark_v31.common import (
    file_sha256,
    iter_gzip_jsonl,
    read_json,
    write_json,
)
from scripts.frozen_artifact_paths import resolve_frozen_v15_artifact


DB_IDENTIFIER_NAME = re.compile(
    r"^(?:NCGC[-_A-Z0-9]*|CCMSLIB[A-Z0-9_-]*|HMDB\d+|CHEBI:\d+|"
    r"CID:\d+|DTXSID\d+|MSBNK[-_A-Z0-9]+|LM[A-Z0-9]+)$",
    re.IGNORECASE,
)


def _snapshot_directory(path: Path) -> dict[str, str]:
    return {
        str(item.relative_to(path)): file_sha256(item)
        for item in sorted(path.rglob("*"))
        if item.is_file()
    }


def audit(artifact_dir: Path) -> dict[str, Any]:
    benchmark = read_json(artifact_dir / "benchmark_receipt.json")
    eligibility = read_json(artifact_dir / "eligibility_summary.json")
    metadata = read_json(artifact_dir / "metadata_derivation_audit.json")
    source_screening = read_json(artifact_dir / "source_screening_receipt.json")
    ce = read_json(artifact_dir / "collision_energy_readiness.json")
    split = read_json(artifact_dir / "split_isolation_readiness.json")
    near_mirror = read_json(artifact_dir / "near_mirror_calibration_gate.json")
    chemistry = read_json(artifact_dir / "chemistry_hypothesis_receipt.json")
    rule_contract = read_json(artifact_dir / "rule_contract_receipt.json")
    runtime_export = read_json(artifact_dir / "runtime_export_gate.json")
    runtime_gate_path = ROOT / "config/feature_fragments/runtime_gate_v2.json"
    runtime_gate = read_json(runtime_gate_path)
    manifest_path = resolve_frozen_v15_artifact(
        benchmark["input"]["v3_record_manifest_path"],
        benchmark["input"]["v3_record_manifest_sha256"], ROOT,
    )

    artifact_hash_checks = {
        name: {
            "expected": expected,
            "observed": file_sha256(artifact_dir / name),
            "matches": file_sha256(artifact_dir / name) == expected,
        }
        for name, expected in benchmark["artifact_sha256"].items()
    }
    config_dir = ROOT / "config/feature_fragment_benchmark_v31"
    config_paths = {
        "preregistration.json": config_dir / "preregistration.json",
        "metadata_derivation_contract.json": config_dir
        / "metadata_derivation_contract.json",
        "source_screening_contract.json": config_dir
        / "source_screening_contract.json",
        "chemistry_hypothesis_contract.json": config_dir
        / "chemistry_hypothesis_contract.json",
        "class_hierarchy_contract.json": config_dir
        / "class_hierarchy_contract.json",
        "rule_schema.json": config_dir / "rule_schema.json",
    }
    config_hashes = {name: file_sha256(path) for name, path in config_paths.items()}
    gap_count = 0
    gap_passes = Counter()
    gap_runtime_violations = 0
    for row in iter_gzip_jsonl(artifact_dir / "eligibility_gap_matrix.jsonl.gz"):
        gap_count += 1
        gap_runtime_violations += row.get("runtime_effect") is not False
        for lane, result in row["evidence_lanes"].items():
            gap_passes[lane] += result["gap"][
                "passes_count_and_precision_planning_gate"
            ]
    database_identifier_names = []
    stage_or_claim_violations = 0
    manifest_record_count = 0
    for row in iter_gzip_jsonl(manifest_path):
        manifest_record_count += 1
        name = str((row.get("compound") or {}).get("display_name") or "")
        if DB_IDENTIFIER_NAME.fullmatch(name):
            database_identifier_names.append(str(row.get("spectrum_id")))
        scope = row.get("identification_scope") or {}
        if (
            row.get("evidence_stage") != "raw_msms_spectrum"
            or scope.get("standard_confirmed_by_mssystem") is not False
        ):
            stage_or_claim_violations += 1
    current_v2 = _snapshot_directory(
        ROOT / "docs/audits/feature_fragment_benchmark_v2"
    )
    current_v3 = _snapshot_directory(
        ROOT / "docs/audits/feature_fragment_benchmark_v3"
    )
    runtime_fail_closed = bool(
        runtime_gate.get("deployment_state") == "fail_closed"
        and runtime_gate.get("runtime_evaluation_enabled") is False
        and runtime_gate.get("runtime_adapter_available") is False
        and runtime_gate.get("exported_rule_count") == 0
        and runtime_gate.get("identity_rank_influence") is False
        and runtime_gate.get("runtime_influence") is False
        and float(runtime_gate.get("score_contribution", 0.0)) == 0.0
    )
    checks = {
        "artifact_hashes_match": all(
            row["matches"] for row in artifact_hash_checks.values()
        ),
        "configuration_hashes_match": config_hashes
        == benchmark["configuration_sha256"],
        "input_manifest_hash_and_count_match": (
            file_sha256(manifest_path)
            == benchmark["input"]["v3_record_manifest_sha256"]
            and manifest_record_count == benchmark["input"]["record_count"]
        ),
        "v2_artifacts_unchanged": current_v2
        == benchmark["immutability"]["v2_artifact_hashes_before"]
        == benchmark["immutability"]["v2_artifact_hashes_after"],
        "v3_artifacts_unchanged": current_v3
        == benchmark["immutability"]["v3_artifact_hashes_before"]
        == benchmark["immutability"]["v3_artifact_hashes_after"],
        "gap_matrix_recomputes_summary": (
            gap_count == eligibility["cell_count"]
            and dict(sorted(gap_passes.items()))
            == eligibility["cells_passing_count_and_precision_planning_gate"]
        ),
        "gap_matrix_has_no_runtime_effect": gap_runtime_violations == 0,
        "current_and_v31_quantitative_lanes_are_empty": all(
            eligibility["lane_record_counts"][lane] == 0
            for lane in (
                "current_positive",
                "current_negative",
                "v31_positive",
                "v31_negative",
            )
        ),
        "all_source_datasets_remain_audit_only": source_screening[
            "status_counts"
        ]
        == {"audit_only": source_screening["source_count"]},
        "static_derivation_changed_only_permitted_field": (
            set(metadata["derived_field_counts"])
            <= set(metadata["permitted_vendor_manual_derived_fields"])
            and metadata["protected_field_violation_count"] == 0
            and metadata["source_quantitative_eligibility_flags_modified"] is False
        ),
        "near_mirror_gate_uses_confidence_bound_and_stays_off": (
            near_mirror["minimum_independent_negative_pairs_if_zero_false_links"]
            == 299
            and near_mirror["near_mirror_edges_enabled"] is False
        ),
        "three_isolation_views_fail_without_positive_records": (
            split["positive_lr_eligible_record_count"] == 0
            and all(
                view[
                    "three_nonempty_zero_overlap_partitions_structurally_possible"
                ]
                is False
                for view in split["views"].values()
            )
        ),
        "ce_models_and_cross_unit_conversion_are_zero": (
            ce["ce_response_model_count"] == 0
            and ce["cross_unit_conversion_applied"] is False
        ),
        "chemistry_and_rule_contract_add_no_fixed_score": (
            chemistry["formula_compatibility_is_a_fixed_score"] is False
            and chemistry["runtime_effect"] is False
            and rule_contract["validated_rule_count"] == 0
            and rule_contract["runtime_exported_rule_count"] == 0
            and bool(rule_contract["forbidden_fixed_score_fields"])
        ),
        "shared_runtime_and_v31_export_stay_fail_closed": (
            runtime_fail_closed
            and runtime_export["runtime_exported_rule_count"] == 0
            and runtime_export["runtime_ranking_enabled"] is False
            and runtime_export["identity_rank_influence"] is False
            and float(runtime_export["score_contribution"]) == 0.0
        ),
        "database_identifiers_not_used_as_names": not database_identifier_names,
        "evidence_stage_and_standard_claims_valid": stage_or_claim_violations == 0,
    }
    failures = sorted(name for name, passed in checks.items() if not passed)
    code_paths = [
        *sorted((ROOT / "matching/feature_fragment_benchmark_v31").glob("*.py")),
        ROOT / "scripts/run_feature_fragment_benchmark_v31.py",
        ROOT / "scripts/audit_feature_fragment_benchmark_v31.py",
        ROOT / "scripts/frozen_artifact_paths.py",
    ]
    return {
        "schema_version": "3.1.0",
        "receipt_id": "v15-feature-fragment-benchmark-v31-audit-20260827",
        "status": "pass" if not failures else "fail",
        "benchmark_state": benchmark["status"],
        "runtime_state": "fail_closed",
        "checks": checks,
        "failures": failures,
        "code_sha256": {
            str(path.relative_to(ROOT)): file_sha256(path) for path in code_paths
        },
        "configuration_sha256": config_hashes,
        "artifact_hash_checks": artifact_hash_checks,
        "corpus": benchmark["corpus"],
        "eligibility": {
            key: eligibility[key]
            for key in (
                "cell_count",
                "frozen_leaf_class_count",
                "identity_class_conflict_count",
                "lane_record_counts",
                "lane_identity_counts",
                "cells_passing_count_and_precision_planning_gate",
                "planning_required_identities_per_target_and_background_arm",
                "metadata_derivation_event_counts",
            )
        },
        "metadata_derivation": metadata,
        "source_screening": {
            "source_count": source_screening["source_count"],
            "status_counts": source_screening["status_counts"],
        },
        "split_isolation": split,
        "collision_energy": {
            "group_count": ce["group_count"],
            "ce_response_model_count": ce["ce_response_model_count"],
            "cross_unit_conversion_applied": ce["cross_unit_conversion_applied"],
        },
        "near_mirror": near_mirror,
        "rules": benchmark["rules"],
        "runtime": benchmark["runtime"],
        "next_state": benchmark["next_state"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=ROOT / "docs/audits/feature_fragment_benchmark_v31",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V31_RECEIPT_20260827.json",
    )
    args = parser.parse_args()
    receipt = audit(args.artifact_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, receipt)
    print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if receipt["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
