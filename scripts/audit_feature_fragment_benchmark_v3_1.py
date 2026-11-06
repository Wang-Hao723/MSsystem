#!/usr/bin/env python3
"""Recompute V3.1 invariants and write the structured phase-1 receipt."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark.corpus import file_sha256
from matching.feature_fragment_benchmark_v3_1.independence import cluster_mirrors_v31
from matching.feature_fragment_benchmark_v3_1.raw_metadata import (
    load_psi_ms_mapping,
    validate_psi_ms_mapping,
)
from matching.feature_fragment_benchmark_v3_1.workflow import (
    ARTIFACT_NAMES,
    V31_CONFIG_NAMES,
    _dependency_hashes,
)


DB_IDENTIFIER_NAME = re.compile(
    r"^(?:NCGC[-_A-Z0-9]*|CCMSLIB[A-Z0-9_-]*|HMDB\d+|CHEBI:\d+|"
    r"CID:\d+|DTXSID\d+|MSBNK[-_A-Z0-9]+|LM[A-Z0-9]+)$",
    re.IGNORECASE,
)

POSITIVE_ONLY_MODES = {
    "stepped",
    "CES",
    "ramp",
    "rolling",
    "dynamic",
    "MSE",
    "bbCID",
    "all_ions",
}


def _json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _gzip_rows(path: Path) -> Iterable[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(f"invalid_gzip_jsonl_line:{line_number}") from error


def _minimal_record(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "spectrum_id": record["spectrum_id"],
        "source": record["source"],
        "compound": record["compound"],
        "acquisition": record["acquisition"],
        "raw_acquisition": record["raw_acquisition"],
        "independence": record["independence"],
        "quantitative_eligibility": record["quantitative_eligibility"],
        "v3_1_amendment": record["v3_1_amendment"],
    }


def _runtime_fail_closed(runtime_gate: dict[str, Any]) -> bool:
    return bool(
        runtime_gate.get("deployment_state") == "fail_closed"
        and runtime_gate.get("runtime_evaluation_enabled") is False
        and runtime_gate.get("runtime_adapter_available") is False
        and runtime_gate.get("exported_rule_count") == 0
        and runtime_gate.get("runtime_influence") is False
        and runtime_gate.get("identity_rank_influence") is False
        and float(runtime_gate.get("score_contribution", 0.0)) == 0.0
    )


def audit(artifact_dir: Path) -> dict[str, Any]:
    benchmark = _json(artifact_dir / "benchmark_receipt.json")
    amendment_path = ROOT / "config/feature_fragment_benchmark_v3_1/amendment.json"
    amendment = _json(amendment_path)
    qc_policy = _json(
        ROOT / "config/feature_fragment_benchmark_v3_1/spectral_qc_policy.json"
    )
    source_amendments = _json(
        ROOT
        / "config/feature_fragment_benchmark_v3_1/public_source_amendments.json"
    )
    mirror_policy = _json(
        ROOT / "config/feature_fragment_benchmark_v3_1/mirror_policy.json"
    )
    receipt_schema = _json(
        ROOT
        / "config/feature_fragment_benchmark_v3_1/benchmark_receipt_schema.json"
    )
    record_schema = _json(
        ROOT
        / "config/feature_fragment_benchmark_v3_1/spectrum_metadata_amendment_schema.json"
    )
    receipt_schema_errors = list(Draft202012Validator(receipt_schema).iter_errors(benchmark))
    artifact_hash_checks = {
        name: {
            "expected": benchmark["artifact_sha256"][name],
            "observed": file_sha256(artifact_dir / name),
            "matches": file_sha256(artifact_dir / name)
            == benchmark["artifact_sha256"][name],
        }
        for name in ARTIFACT_NAMES
    }
    configuration_sha256 = {
        name: file_sha256(
            ROOT / "config/feature_fragment_benchmark_v3_1" / name
        )
        for name in V31_CONFIG_NAMES
    }
    runtime_gate_path = ROOT / "config/feature_fragments/runtime_gate_v2.json"
    runtime_gate = _json(runtime_gate_path)
    formal_v3 = _json(
        ROOT / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_RECEIPT_20260827.json"
    )
    base_v3_benchmark = _json(
        ROOT / "docs/audits/feature_fragment_benchmark_v3/benchmark_receipt.json"
    )
    base_paths = {
        "machine_preregistration": ROOT
        / amendment["base_release"]["machine_preregistration_path"],
        "human_preregistration": ROOT
        / amendment["base_release"]["human_preregistration_path"],
        "phase1_receipt": ROOT / amendment["base_release"]["phase1_receipt_path"],
        "benchmark_receipt": ROOT
        / amendment["base_release"]["benchmark_receipt_path"],
    }
    base_expected = {
        "machine_preregistration": amendment["base_release"]
        ["machine_preregistration_sha256"],
        "human_preregistration": amendment["base_release"]
        ["human_preregistration_sha256"],
        "phase1_receipt": amendment["base_release"]["phase1_receipt_sha256"],
        "benchmark_receipt": amendment["base_release"]["benchmark_receipt_sha256"],
    }
    base_observed = {name: file_sha256(path) for name, path in base_paths.items()}
    base_artifact_checks = {
        name: {
            "expected": row["expected"],
            "observed": file_sha256(
                ROOT / "docs/audits/feature_fragment_benchmark_v3" / name
            ),
            "matches": row["expected"]
            == file_sha256(
                ROOT / "docs/audits/feature_fragment_benchmark_v3" / name
            ),
        }
        for name, row in formal_v3["artifact_hash_checks"].items()
    }

    record_validator = Draft202012Validator(record_schema)
    record_count = 0
    record_schema_error_count = 0
    duplicate_ids = 0
    seen_ids: set[str] = set()
    monotonic_positive_violations = 0
    monotonic_negative_violations = 0
    multi_energy_negative_violations = 0
    conversion_violations = 0
    derived_ev_violations = 0
    database_identifier_names = 0
    evidence_stage_violations = 0
    standard_claim_violations = 0
    mirror_assignment_violations = 0
    peak_hash_acquisition_group_violations = 0
    minimum_peak_gate_active_records = 0
    reason_counts: Counter[str] = Counter()
    minimal_records: list[dict[str, Any]] = []
    for record in _gzip_rows(artifact_dir / "record_manifest.jsonl.gz"):
        record_count += 1
        record_schema_error_count += sum(1 for _ in record_validator.iter_errors(record))
        spectrum_id = str(record["spectrum_id"])
        duplicate_ids += spectrum_id in seen_ids
        seen_ids.add(spectrum_id)
        amendment_row = record["v3_1_amendment"]
        eligibility = record["quantitative_eligibility"]
        monotonic_positive_violations += bool(
            eligibility["positive_lr_eligible"]
            and not amendment_row["base_positive_lr_eligible"]
        )
        monotonic_negative_violations += bool(
            eligibility["negative_lr_eligible"]
            and not amendment_row["base_negative_lr_eligible"]
        )
        mode = str((record["acquisition"].get("collision_energy") or {}).get("mode"))
        multi_energy_negative_violations += bool(
            eligibility["negative_lr_eligible"]
            and (
                mode in POSITIVE_ONLY_MODES
                or record["acquisition"].get("activation_location") == "in_source"
            )
        )
        conversion_violations += amendment_row["nce_to_ev_conversion_applied"] is not False
        derived_ev_violations += amendment_row["nce_to_ev_derived_value"] is not None
        database_identifier_names += bool(
            DB_IDENTIFIER_NAME.fullmatch(str(record["compound"]["display_name"]))
        )
        evidence_stage_violations += record.get("evidence_stage") != "raw_msms_spectrum"
        standard_claim_violations += (
            record.get("identification_scope", {}).get(
                "standard_confirmed_by_mssystem"
            )
            is not False
        )
        peak_hash_acquisition_group_violations += (
            record["independence"].get(
                "peak_hash_used_in_acquisition_replication_group"
            )
            is not False
        )
        minimum_peak_gate_active_records += (
            amendment_row["spectral_qc"]["minimum_peak_gate_active"] is True
        )
        reason_counts.update(eligibility.get("ineligible_reason_codes") or [])
        minimal_records.append(_minimal_record(record))

    recomputed_mirror = cluster_mirrors_v31(minimal_records)
    for record in minimal_records:
        spectrum_id = str(record["spectrum_id"])
        mirror_assignment_violations += (
            record["independence"].get("mirror_cluster_id")
            != recomputed_mirror["assignment"][spectrum_id]
        )
    policy_summary = _json(artifact_dir / "policy_summary.json")
    feasibility = _json(artifact_dir / "domain_feasibility_report.json")
    runtime_export = _json(artifact_dir / "runtime_export_gate.json")
    review = _json(artifact_dir / "review_candidates.json")
    tested = _json(artifact_dir / "sealed_test_rules.json")
    feasibility_text = json.dumps(feasibility, ensure_ascii=False)
    mapping = load_psi_ms_mapping()
    psi_mapping_valid = True
    try:
        validate_psi_ms_mapping(mapping)
    except ValueError:
        psi_mapping_valid = False
    dependency_sha256 = _dependency_hashes(ROOT)
    imported_v3_dependency_checks = {
        name: {
            "expected": formal_v3["code_sha256"].get(name),
            "observed": digest,
            "matches": formal_v3["code_sha256"].get(name) in {None, digest},
        }
        for name, digest in dependency_sha256.items()
        if name != "matching/feature_fragment_benchmark/mirror.py"
    }
    recomputed_mirror_summary = {
        key: recomputed_mirror[key]
        for key in (
            "exact_edge_count_by_type",
            "unique_exact_edge_count",
            "near_duplicate_edge_count",
            "component_count",
            "component_size_counts",
            "near_candidate_status",
        )
    }
    recomputed_mirror_summary["component_size_counts"] = {
        str(key): value
        for key, value in recomputed_mirror_summary["component_size_counts"].items()
    }
    checks = {
        "artifact_hashes_match": all(
            row["matches"] for row in artifact_hash_checks.values()
        ),
        "receipt_schema_valid": not receipt_schema_errors,
        "record_schema_valid": record_schema_error_count == 0,
        "frozen_record_count_matches": record_count
        == amendment["base_release"]["accepted_audit_record_count"]
        == benchmark["input"]["record_count"],
        "spectrum_ids_unique": duplicate_ids == 0,
        "v3_1_never_promotes_v3_ineligible_records": monotonic_positive_violations
        == monotonic_negative_violations
        == 0,
        "multi_energy_and_in_source_lr_negative_forbidden": multi_energy_negative_violations
        == 0,
        "nce_to_ev_conversion_absent": conversion_violations
        == derived_ev_violations
        == policy_summary["nce_to_ev_conversion_count"]
        == 0,
        "database_identifiers_not_used_as_compound_names": database_identifier_names
        == 0,
        "evidence_stage_remains_raw_reference_spectrum": evidence_stage_violations
        == standard_claim_violations
        == 0,
        "acquisition_replication_groups_exclude_peak_hash": peak_hash_acquisition_group_violations
        == 0,
        "mirror_assignments_recompute": mirror_assignment_violations == 0,
        "mirror_summary_recomputes": recomputed_mirror_summary
        == {
            key: benchmark["mirrors"][key]
            for key in (
                "exact_edge_count_by_type",
                "unique_exact_edge_count",
                "near_duplicate_edge_count",
                "component_count",
                "component_size_counts",
                "near_candidate_status",
            )
        },
        "near_mirrors_remain_uncalibrated_and_unlinked": benchmark["mirrors"]
        ["near_duplicate_edge_count"]
        == 0
        and benchmark["mirrors"]["near_calibration"].get(
            "near_mirror_links_allowed"
        )
        is False,
        "minimum_peak_gate_remains_inactive": qc_policy["peak_presence"]
        ["minimum_deduplicated_peak_count_gate"]["active"]
        is False
        and minimum_peak_gate_active_records == 0,
        "frozen_qc_values_are_explicit": qc_policy["peak_presence"]
        ["minimum_relative_intensity"]
        == 0.01
        and qc_policy["peak_presence"]["mass_tolerance"]
        == {
            "rule": "maximum_of_absolute_da_and_ppm",
            "absolute_da": 0.01,
            "ppm": 20.0,
        },
        "precursor_purity_unknown_blocks_lr_negative": policy_summary[
            "v3_1_negative_lr_eligible_record_count"
        ]
        == 0
        and qc_policy["precursor_purity"]["threshold_status"]
        == "not_calibrated",
        "psi_ms_mapping_is_semantically_valid": psi_mapping_valid,
        "psi_ms_1002482_is_not_activation": "MS:1002482"
        not in mapping["mzml"]["activation"]
        and "MS:1002482" in mapping["mzml"]["terms_not_activation"],
        "mzml_and_mzxml_origins_are_distinct": mapping["mzml"]["value_origin"]
        == "mzml_cv_param"
        and mapping["mzxml"]["value_origin"] == "mzxml_attribute",
        "all_exact_domains_infeasible_only_for_frozen_release": feasibility[
            "domain_count"
        ]
        == 12
        and feasibility["status_counts"]
        == {"infeasible_for_frozen_release": 12}
        and "infeasible_forever" not in feasibility_text,
        "partial_pooling_cannot_rescue_failed_domain": feasibility[
            "partial_pooling"
        ]["may_rescue_failed_exact_domain"]
        is False,
        "mona_is_screening_only_not_independent_replication": source_amendments[
            "sources"
        ][0]["integration_status"]
        == "screened_not_integrated"
        and source_amendments["sources"][0][
            "aggregator_counts_as_new_source_or_laboratory"
        ]
        is False
        and benchmark["input"]["new_public_spectra_integrated"] == 0,
        "base_v3_named_hashes_unchanged": base_observed == base_expected,
        "base_v3_artifacts_unchanged": all(
            row["matches"] for row in base_artifact_checks.values()
        ),
        "imported_v3_dependencies_match_frozen_receipt": all(
            row["matches"] for row in imported_v3_dependency_checks.values()
        ),
        "v3_1_config_hashes_match_benchmark": configuration_sha256
        == benchmark["configuration_sha256"],
        "imported_dependency_hashes_match_benchmark": dependency_sha256
        == benchmark["imported_dependency_sha256"],
        "runtime_gate_remains_fail_closed": _runtime_fail_closed(runtime_gate),
        "runtime_gate_hash_unchanged": file_sha256(runtime_gate_path)
        == formal_v3["runtime"]["runtime_gate_sha256"]
        == benchmark["runtime"]["runtime_gate_sha256"],
        "runtime_export_ranking_and_score_remain_zero": runtime_export[
            "runtime_exported_rule_count"
        ]
        == review["runtime_exported_rule_count"]
        == benchmark["rules"]["runtime_exported_rule_count"]
        == 0
        and runtime_export["runtime_ranking_enabled"] is False
        and float(runtime_export["score_contribution"]) == 0.0,
        "phase_2_remains_unapproved": runtime_export["phase_2_allowed"] is False,
        "no_identity_formula_isomer_or_standard_upgrade": not any(
            bool(rule.get(field))
            for rule in tested
            for field in (
                "identity_upgrade_allowed",
                "formula_upgrade_allowed",
                "isomer_resolution_allowed",
                "standard_confirmation_allowed",
            )
        ),
        "benchmark_fails_closed_at_frozen_release_capacity_gate": benchmark[
            "status"
        ]
        == "infeasible_for_frozen_release"
        and policy_summary["v3_1_positive_lr_eligible_record_count"] == 0,
        "next_state_waits_for_user_review": benchmark["next_state"]
        == runtime_export["next_state"]
        == review["next_state"]
        == "waiting_user_review",
        "formal_v3_was_pass_and_failed_closed": formal_v3["status"] == "pass"
        and base_v3_benchmark["status"] == "failed_closed",
        "mirror_policy_has_three_separate_exact_edges": [
            row["edge_type"] for row in mirror_policy["exact_edges"]
        ]
        == ["canonical_usi_core", "origin_raw_locator", "exact_peak_hash"]
        and mirror_policy["exclusive_first_key"] is False,
    }
    failures = sorted(name for name, passed in checks.items() if not passed)
    code_paths = [
        *sorted((ROOT / "matching/feature_fragment_benchmark_v3_1").glob("*.py")),
        ROOT / "scripts/run_feature_fragment_benchmark_v3_1.py",
        ROOT / "scripts/audit_feature_fragment_benchmark_v3_1.py",
    ]
    return {
        "schema_version": "3.1.0",
        "receipt_id": "v15-feature-fragment-benchmark-v3.1-phase1-20260827",
        "status": "pass" if not failures else "fail",
        "benchmark_state": benchmark["status"],
        "runtime_state": "fail_closed",
        "next_state": "waiting_user_review",
        "checks": checks,
        "failures": failures,
        "code_sha256": {
            str(path.relative_to(ROOT)): file_sha256(path) for path in code_paths
        },
        "configuration_sha256": configuration_sha256,
        "artifact_hash_checks": artifact_hash_checks,
        "base_v3_named_hash_checks": {
            name: {
                "expected": base_expected[name],
                "observed": base_observed[name],
                "matches": base_expected[name] == base_observed[name],
            }
            for name in base_expected
        },
        "base_v3_artifact_hash_checks": base_artifact_checks,
        "imported_v3_dependency_hash_checks": imported_v3_dependency_checks,
        "input": benchmark["input"],
        "policy": {
            **policy_summary,
            "recomputed_ineligibility_reason_counts": dict(
                sorted(reason_counts.items())
            ),
        },
        "mirrors": benchmark["mirrors"],
        "feasibility": benchmark["feasibility"],
        "split": benchmark["split"],
        "rules": benchmark["rules"],
        "runtime": benchmark["runtime"],
        "license_and_source": benchmark["license_and_source"],
        "limitations": benchmark["limitations"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=ROOT / "docs/audits/feature_fragment_benchmark_v3_1",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_1_RECEIPT_20260827.json",
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
