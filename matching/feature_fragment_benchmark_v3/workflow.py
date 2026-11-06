"""End-to-end V3 research benchmark with immutable audit artifacts."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from copy import deepcopy
from pathlib import Path
from typing import Any

from matching.feature_fragment_benchmark.corpus import canonical_bytes, file_sha256

from .export import export_review_candidates
from .independence import cluster_mirrors
from .rules import (
    build_identity_observations,
    calibrate_rules,
    discover_rules,
    evaluate_sealed_rules,
)
from .split import build_origin_isolated_split


V2_ARTIFACT_NAMES = (
    "benchmark_receipt.json",
    "calibration_model.json",
    "ce_dependence.json",
    "corpus_receipt.json",
    "discovery_rules.json",
    "mirror_clusters.jsonl.gz",
    "record_manifest.jsonl.gz",
    "rejects.jsonl.gz",
    "runtime_export_gate.json",
    "sealed_test_rules.json",
    "split_manifest.jsonl.gz",
)


def _json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_gzip_jsonl(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    with gzip.GzipFile(filename=str(path), mode="wb", mtime=0) as compressed:
        with io.TextIOWrapper(compressed, encoding="utf-8") as handle:
            for row in rows:
                handle.write(canonical_bytes(row).decode("utf-8") + "\n")


def _copy_jsonl_to_gzip(source: Path, destination: Path) -> None:
    with source.open("rt", encoding="utf-8") as input_handle, gzip.GzipFile(
        filename=str(destination), mode="wb", mtime=0
    ) as compressed:
        with io.TextIOWrapper(compressed, encoding="utf-8") as output_handle:
            for line in input_handle:
                if line.strip():
                    output_handle.write(line)


def _snapshot_v2(v2_dir: Path) -> dict[str, str]:
    return {
        name: file_sha256(v2_dir / name)
        for name in V2_ARTIFACT_NAMES
        if (v2_dir / name).exists()
    }


def _manifest_record(
    record: Mapping[str, Any], split_assignment: Mapping[str, str]
) -> dict[str, Any]:
    row = deepcopy(record)
    row.pop("peaks", None)
    spectrum_id = str(record["spectrum_id"])
    row["split"] = split_assignment.get(spectrum_id, "not_quantitatively_eligible")
    return row


def run_benchmark_v3(
    records: Sequence[Mapping[str, Any]],
    *,
    output_dir: str | Path,
    corpus_receipt: Mapping[str, Any],
    rejects_jsonl: str | Path,
    preregistration_path: str | Path,
    public_sources_path: str | Path,
    instrument_domains_path: str | Path,
    license_policy_path: str | Path,
    evidence_sources_path: str | Path,
    metadata_schema_path: str | Path,
    runtime_gate_path: str | Path,
    v2_artifact_dir: str | Path,
    v2_record_ids: set[str] | None = None,
) -> dict[str, Any]:
    """Run V3 once; passing rules remain review-only by construction."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    preregistration_path = Path(preregistration_path)
    public_sources_path = Path(public_sources_path)
    instrument_domains_path = Path(instrument_domains_path)
    license_policy_path = Path(license_policy_path)
    evidence_sources_path = Path(evidence_sources_path)
    metadata_schema_path = Path(metadata_schema_path)
    runtime_gate_path = Path(runtime_gate_path)
    v2_artifact_dir = Path(v2_artifact_dir)
    rejects_jsonl = Path(rejects_jsonl)
    preregistration = _json(preregistration_path)
    public_sources = _json(public_sources_path)
    runtime_gate = _json(runtime_gate_path)
    v2_before = _snapshot_v2(v2_artifact_dir)

    working_records = [deepcopy(row) for row in records]
    mirror = cluster_mirrors(working_records)
    for row in working_records:
        row["independence"]["mirror_cluster_id"] = mirror["assignment"][
            str(row["spectrum_id"])
        ]
    quantitative = [
        row
        for row in working_records
        if row["quantitative_eligibility"]["positive_lr_eligible"] is True
    ]
    split = build_origin_isolated_split(
        quantitative,
        seed=preregistration["splitting"]["seed"],
        test_fraction=float(preregistration["splitting"]["sealed_test_fraction"]),
        calibration_fraction_within_development=float(
            preregistration["splitting"][
                "calibration_fraction_within_development"
            ]
        ),
    )
    discovery_rules: list[dict[str, Any]] = []
    calibrated: list[dict[str, Any]] = []
    tested: list[dict[str, Any]] = []
    observation_counts: dict[str, int] = {}
    if split["split_feasible"]:
        records_by_split = {
            name: [
                row
                for row in quantitative
                if split["assignment"][str(row["spectrum_id"])] == name
            ]
            for name in ("discovery", "calibration", "test")
        }
        observations = {
            name: build_identity_observations(rows)
            for name, rows in records_by_split.items()
        }
        observation_counts = {
            name: len(rows) for name, rows in observations.items()
        }
        discovery_rules = discover_rules(observations["discovery"])
        calibrated = calibrate_rules(discovery_rules, observations["calibration"])
        tested = evaluate_sealed_rules(calibrated, observations["test"])

    review_export = export_review_candidates(tested)
    ce_models: list[dict[str, Any]] = []
    v2_ids = v2_record_ids or set()
    new_audit_ids = {
        str(row["spectrum_id"])
        for row in working_records
        if str(row["spectrum_id"]) not in v2_ids
    }
    new_quantitative_ids = {
        str(row["spectrum_id"])
        for row in quantitative
        if str(row["spectrum_id"]) not in v2_ids
    }
    eligibility_reasons = Counter(
        reason
        for row in working_records
        for reason in row["quantitative_eligibility"]["ineligible_reason_codes"]
    )
    evidence_grades = Counter(row.get("evidence_grade", "unavailable") for row in tested)
    rule_failures = Counter(
        failure
        for row in tested
        for failure in row.get("sealed_test_gate", {}).get("failures", [])
    )

    _write_json(output_dir / "discovery_rules.json", discovery_rules)
    _write_json(output_dir / "calibration_model.json", calibrated)
    _write_json(output_dir / "sealed_test_rules.json", tested)
    _write_json(output_dir / "ce_models.json", ce_models)
    _write_json(output_dir / "review_candidates.json", review_export)
    _write_json(
        output_dir / "runtime_export_gate.json",
        {
            "schema_version": "3.0.0",
            "status": "fail_closed_phase_1",
            "runtime_exported_rule_count": 0,
            "runtime_catalog_modified": False,
            "runtime_gate_modified": False,
            "runtime_ranking_enabled": False,
            "score_contribution": 0.0,
            "review_candidate_count": review_export["review_candidate_count"],
            "next_state": "waiting_user_review",
        },
    )
    _write_gzip_jsonl(
        output_dir / "record_manifest.jsonl.gz",
        (
            _manifest_record(row, split.get("assignment", {}))
            for row in sorted(working_records, key=lambda item: item["spectrum_id"])
        ),
    )
    _write_gzip_jsonl(
        output_dir / "mirror_clusters.jsonl.gz", mirror["clusters"]
    )
    _write_gzip_jsonl(
        output_dir / "split_manifest.jsonl.gz",
        (
            {
                "spectrum_id": spectrum_id,
                "split": split_name,
            }
            for spectrum_id, split_name in sorted(split.get("assignment", {}).items())
        ),
    )
    _copy_jsonl_to_gzip(rejects_jsonl, output_dir / "rejects.jsonl.gz")

    artifact_names = (
        "calibration_model.json",
        "ce_models.json",
        "discovery_rules.json",
        "mirror_clusters.jsonl.gz",
        "record_manifest.jsonl.gz",
        "rejects.jsonl.gz",
        "review_candidates.json",
        "runtime_export_gate.json",
        "sealed_test_rules.json",
        "split_manifest.jsonl.gz",
    )
    artifact_hashes = {
        name: file_sha256(output_dir / name) for name in artifact_names
    }
    v2_after = _snapshot_v2(v2_artifact_dir)
    v2_unchanged = v2_before == v2_after
    runtime_fail_closed = bool(
        runtime_gate.get("deployment_state") == "fail_closed"
        and runtime_gate.get("runtime_evaluation_enabled") is False
        and runtime_gate.get("runtime_adapter_available") is False
        and runtime_gate.get("exported_rule_count") == 0
        and runtime_gate.get("identity_rank_influence") is False
        and runtime_gate.get("runtime_influence") is False
        and float(runtime_gate.get("score_contribution", 0.0)) == 0.0
    )
    passed_rules = review_export["review_candidate_count"]
    status = "candidates_waiting_review" if passed_rules else "failed_closed"
    receipt = {
        "schema_version": "3.0.0",
        "receipt_id": "v15-feature-fragment-benchmark-v3-20260827",
        "status": status,
        "generated_on": "2026-08-27",
        "next_state": "waiting_user_review",
        "runtime": {
            "enabled": False,
            "fail_closed_verified": runtime_fail_closed,
            "runtime_gate_path": str(runtime_gate_path),
            "runtime_gate_sha256": file_sha256(runtime_gate_path),
            "runtime_exported_rule_count": 0,
            "runtime_ranking_enabled": False,
            "score_contribution": 0.0,
        },
        "preregistration": {
            "path": str(preregistration_path),
            "sha256": file_sha256(preregistration_path),
            "id": preregistration["preregistration_id"],
        },
        "configuration_sha256": {
            "public_sources.json": file_sha256(public_sources_path),
            "instrument_domains.json": file_sha256(instrument_domains_path),
            "license_policy.json": file_sha256(license_policy_path),
            "evidence_sources.json": file_sha256(evidence_sources_path),
            "spectrum_metadata_schema.json": file_sha256(metadata_schema_path),
        },
        "inputs": {
            "repository_release": corpus_receipt.get("repository_release"),
            "raw_record_count": corpus_receipt.get("raw_record_count"),
            "accepted_audit_record_count": len(working_records),
            "rejected_record_count": corpus_receipt.get("rejected_record_count"),
            "source_file_sha256": corpus_receipt.get("source_file_sha256", {}),
            "license_tier_counts": corpus_receipt.get("license_tier_counts", {}),
            "source_count": len(
                {row["source"]["source_id"] for row in working_records}
            ),
            "laboratory_count": len(
                {row["source"]["laboratory_id"] for row in working_records}
            ),
            "two_dimensional_identity_count": len(
                {row["compound"]["connectivity"] for row in working_records}
            ),
            "formula_family_count": len(
                {row["compound"].get("formula") for row in working_records}
            ),
            "exact_instrument_domain_count": len(
                {
                    row["acquisition"]["instrument_domain_id"]
                    for row in working_records
                    if row["acquisition"]["instrument_domain_id"]
                }
            ),
        },
        "public_spectrum_additions": {
            "audit_only_records_absent_from_v2_manifest": len(new_audit_ids),
            "quantitative_records_absent_from_v2_manifest": len(
                new_quantitative_ids
            ),
            "screening_register": public_sources.get("screening_register", []),
        },
        "metadata": {
            "scan_range_status_counts": corpus_receipt.get(
                "scan_range_status_counts", {}
            ),
            "scan_range_contains_precursor_counts": corpus_receipt.get(
                "scan_range_contains_precursor_counts", {}
            ),
            "precursor_isolation_status_counts": corpus_receipt.get(
                "precursor_isolation_status_counts", {}
            ),
            "collision_energy_unit_mode_counts": corpus_receipt.get(
                "collision_energy_unit_mode_counts", {}
            ),
            "exact_condition_metadata_complete_counts": corpus_receipt.get(
                "exact_condition_metadata_complete_counts", {}
            ),
            "acquisition_replication_group_count": corpus_receipt.get(
                "acquisition_replication_group_count", 0
            ),
            "positive_lr_eligible_record_count": len(quantitative),
            "negative_lr_eligible_record_count": sum(
                row["quantitative_eligibility"]["negative_lr_eligible"]
                for row in working_records
            ),
            "quantitative_ineligibility_reason_counts": dict(
                sorted(eligibility_reasons.items())
            ),
            "missing_field_counts": corpus_receipt.get(
                "metadata_missing_field_counts", {}
            ),
        },
        "mirror_audit": {
            key: mirror[key]
            for key in (
                "component_count",
                "exact_duplicate_edge_count",
                "near_duplicate_edge_count",
                "component_size_counts",
                "near_calibration",
                "near_candidate_status",
            )
            if key in mirror
        },
        "split": split,
        "identity_observation_counts": observation_counts,
        "collision_energy_models": {
            "model_count": len(ce_models),
            "cross_unit_conversion_applied": False,
        },
        "rules": {
            "discovered": len(discovery_rules),
            "calibrated_total": len(calibrated),
            "selected_for_sealed_test": sum(
                bool(row.get("selected_for_sealed_test")) for row in calibrated
            ),
            "sealed_test_evaluated": len(tested),
            "evidence_grade_counts": dict(sorted(evidence_grades.items())),
            "failure_reason_counts": dict(sorted(rule_failures.items())),
            "review_candidate_count": passed_rules,
            "runtime_exported_rule_count": 0,
            "identity_formula_isomer_standard_upgrade_allowed": False,
        },
        "immutability": {
            "v2_artifact_hashes_before": v2_before,
            "v2_artifact_hashes_after": v2_after,
            "v2_artifacts_unchanged": v2_unchanged,
        },
        "artifact_sha256": artifact_hashes,
        "deviations": [],
        "limitations": [
            "No public record is promoted to Mssystem standard confirmation.",
            "A repository identifier is never used as a compound name.",
            "Missing exact collision-cell, scan, isolation, source, laboratory, or license metadata makes quantitative evidence unavailable rather than zero.",
            "Uncalibrated near-mirror thresholds are not applied.",
            "Passing research candidates, if any, still require a separate user-reviewed runtime adapter phase."
        ],
    }
    _write_json(output_dir / "benchmark_receipt.json", receipt)
    return receipt
