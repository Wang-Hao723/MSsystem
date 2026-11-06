"""Deterministic V3.1 amendment workflow over the frozen V3 public corpus."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from collections import Counter
from collections.abc import Iterable, Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

from matching.feature_fragment_benchmark.corpus import canonical_bytes, file_sha256
from matching.feature_fragment_benchmark_v3.export import export_review_candidates
from matching.feature_fragment_benchmark_v3.rules import build_identity_observations
from matching.feature_fragment_benchmark_v3.split import build_origin_isolated_split

from .feasibility import build_domain_feasibility_report
from .independence import cluster_mirrors_v31
from .policy import apply_v31_policy
from .rules import (
    calibrate_rules,
    discover_rules_with_policy,
    evaluate_sealed_rules,
)


ARTIFACT_NAMES = (
    "calibration_model.json",
    "ce_models.json",
    "discovery_rules.json",
    "domain_feasibility_report.json",
    "mirror_clusters.jsonl.gz",
    "policy_summary.json",
    "record_manifest.jsonl.gz",
    "review_candidates.json",
    "runtime_export_gate.json",
    "sealed_test_rules.json",
    "split_manifest.jsonl.gz",
)

V31_CONFIG_NAMES = (
    "amendment.json",
    "benchmark_receipt_schema.json",
    "evidence_source_amendments.json",
    "mirror_policy.json",
    "psi_ms_mapping.json",
    "public_source_amendments.json",
    "spectral_qc_policy.json",
    "spectrum_metadata_amendment_schema.json",
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


def _iter_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(f"invalid_jsonl_line:{line_number}") from error


def _hash_paths(paths: Mapping[str, Path]) -> dict[str, str]:
    return {name: file_sha256(path) for name, path in paths.items()}


def _base_paths(root: Path) -> dict[str, Path]:
    return {
        "config/feature_fragment_benchmark_v3/preregistration.json": root
        / "config/feature_fragment_benchmark_v3/preregistration.json",
        "docs/V15_FEATURE_FRAGMENT_BENCHMARK_V3_PREREGISTRATION.md": root
        / "docs/V15_FEATURE_FRAGMENT_BENCHMARK_V3_PREREGISTRATION.md",
        "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_RECEIPT_20260827.json": root
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_RECEIPT_20260827.json",
        "docs/audits/feature_fragment_benchmark_v3/benchmark_receipt.json": root
        / "docs/audits/feature_fragment_benchmark_v3/benchmark_receipt.json",
        "config/feature_fragments/runtime_gate_v2.json": root
        / "config/feature_fragments/runtime_gate_v2.json",
    }


def _validate_base_release(
    root: Path,
    amendment: Mapping[str, Any],
    formal_receipt: Mapping[str, Any],
    benchmark_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    base = amendment["base_release"]
    expected = {
        "config/feature_fragment_benchmark_v3/preregistration.json": base[
            "machine_preregistration_sha256"
        ],
        "docs/V15_FEATURE_FRAGMENT_BENCHMARK_V3_PREREGISTRATION.md": base[
            "human_preregistration_sha256"
        ],
        "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_RECEIPT_20260827.json": base[
            "phase1_receipt_sha256"
        ],
        "docs/audits/feature_fragment_benchmark_v3/benchmark_receipt.json": base[
            "benchmark_receipt_sha256"
        ],
    }
    observed = _hash_paths(_base_paths(root))
    failures = [
        name for name, expected_hash in expected.items() if observed[name] != expected_hash
    ]
    runtime_expected = formal_receipt["runtime"]["runtime_gate_sha256"]
    if observed["config/feature_fragments/runtime_gate_v2.json"] != runtime_expected:
        failures.append("config/feature_fragments/runtime_gate_v2.json")
    if formal_receipt.get("status") != "pass":
        failures.append("formal_v3_receipt_not_pass")
    if benchmark_receipt.get("status") != "failed_closed":
        failures.append("base_v3_benchmark_not_failed_closed")
    if failures:
        raise ValueError("base_v3_freeze_mismatch:" + ",".join(sorted(failures)))
    return {
        "expected": expected,
        "observed": observed,
        "runtime_expected": runtime_expected,
        "validated": True,
    }


def _minimal_record(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "spectrum_id": record["spectrum_id"],
        "source": deepcopy(record.get("source") or {}),
        "compound": deepcopy(record.get("compound") or {}),
        "acquisition": deepcopy(record.get("acquisition") or {}),
        "raw_acquisition": deepcopy(record.get("raw_acquisition") or {}),
        "independence": deepcopy(record.get("independence") or {}),
        "quantitative_eligibility": deepcopy(
            record.get("quantitative_eligibility") or {}
        ),
        "v3_1_amendment": deepcopy(record.get("v3_1_amendment") or {}),
    }


def _runtime_fail_closed(runtime_gate: Mapping[str, Any]) -> bool:
    return bool(
        runtime_gate.get("deployment_state") == "fail_closed"
        and runtime_gate.get("runtime_evaluation_enabled") is False
        and runtime_gate.get("runtime_adapter_available") is False
        and runtime_gate.get("exported_rule_count") == 0
        and runtime_gate.get("runtime_influence") is False
        and runtime_gate.get("identity_rank_influence") is False
        and float(runtime_gate.get("score_contribution", 0.0)) == 0.0
    )


def _dependency_hashes(root: Path) -> dict[str, str]:
    paths = [
        root / "matching/feature_fragment_benchmark/corpus.py",
        root / "matching/feature_fragment_benchmark/mirror.py",
        root / "matching/feature_fragment_benchmark_v3/export.py",
        root / "matching/feature_fragment_benchmark_v3/rules.py",
        root / "matching/feature_fragment_benchmark_v3/split.py",
        root / "matching/feature_fragment_benchmark_v3/statistics.py",
    ]
    return {str(path.relative_to(root)): file_sha256(path) for path in paths}


def run_benchmark_v3_1(
    *,
    root: str | Path,
    frozen_corpus_path: str | Path,
    output_dir: str | Path,
) -> dict[str, Any]:
    """Derive V3.1 audit artifacts without modifying V3 or runtime files."""
    root = Path(root).resolve()
    frozen_corpus_path = Path(frozen_corpus_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config_dir = root / "config/feature_fragment_benchmark_v3_1"
    v3_config_dir = root / "config/feature_fragment_benchmark_v3"
    amendment = _json(config_dir / "amendment.json")
    qc_policy = _json(config_dir / "spectral_qc_policy.json")
    mirror_policy = _json(config_dir / "mirror_policy.json")
    source_amendments = _json(config_dir / "public_source_amendments.json")
    domains_config = _json(v3_config_dir / "instrument_domains.json")
    base_preregistration = _json(v3_config_dir / "preregistration.json")
    formal_receipt_path = (
        root / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_RECEIPT_20260827.json"
    )
    benchmark_receipt_path = (
        root / "docs/audits/feature_fragment_benchmark_v3/benchmark_receipt.json"
    )
    formal_receipt = _json(formal_receipt_path)
    benchmark_receipt = _json(benchmark_receipt_path)
    base_before = _validate_base_release(
        root, amendment, formal_receipt, benchmark_receipt
    )
    corpus_sha256 = file_sha256(frozen_corpus_path)
    if corpus_sha256 != amendment["base_release"]["frozen_corpus_sha256"]:
        raise ValueError("frozen_v3_corpus_sha256_mismatch")

    minimal_records: list[dict[str, Any]] = []
    quantitative_records: list[dict[str, Any]] = []
    reason_counts: Counter[str] = Counter()
    ce_mode_counts: Counter[str] = Counter()
    base_positive = 0
    base_negative = 0
    amended_positive = 0
    amended_negative = 0
    nce_to_ev_conversions = 0
    line_count = 0
    for base_record in _iter_jsonl(frozen_corpus_path):
        line_count += 1
        revised = apply_v31_policy(base_record, qc_policy=qc_policy)
        amendment_row = revised["v3_1_amendment"]
        eligibility = revised["quantitative_eligibility"]
        base_positive += amendment_row["base_positive_lr_eligible"] is True
        base_negative += amendment_row["base_negative_lr_eligible"] is True
        amended_positive += eligibility["positive_lr_eligible"] is True
        amended_negative += eligibility["negative_lr_eligible"] is True
        nce_to_ev_conversions += amendment_row["nce_to_ev_conversion_applied"] is True
        ce_mode_counts[str((revised["acquisition"].get("collision_energy") or {}).get("mode") or "unresolved")] += 1
        reason_counts.update(eligibility.get("ineligible_reason_codes") or [])
        minimal_records.append(_minimal_record(revised))
        if eligibility["positive_lr_eligible"] is True:
            quantitative_records.append(revised)
    if line_count != int(amendment["base_release"]["accepted_audit_record_count"]):
        raise ValueError("frozen_v3_corpus_record_count_mismatch")

    mirror = cluster_mirrors_v31(minimal_records)
    for row in minimal_records:
        row["independence"]["mirror_cluster_id"] = mirror["assignment"][
            str(row["spectrum_id"])
        ]
    for row in quantitative_records:
        row["independence"]["mirror_cluster_id"] = mirror["assignment"][
            str(row["spectrum_id"])
        ]

    feasibility = build_domain_feasibility_report(
        minimal_records,
        instrument_domains=domains_config["instrument_domains"],
        thresholds=amendment["unchanged_primary_thresholds"],
        frozen_release_id=str(benchmark_receipt["inputs"]["repository_release"]),
        frozen_corpus_sha256=corpus_sha256,
    )
    split = build_origin_isolated_split(
        quantitative_records,
        seed=base_preregistration["splitting"]["seed"],
        test_fraction=float(
            base_preregistration["splitting"]["sealed_test_fraction"]
        ),
        calibration_fraction_within_development=float(
            base_preregistration["splitting"]
            ["calibration_fraction_within_development"]
        ),
    )
    discovery_rules: list[dict[str, Any]] = []
    calibrated_rules: list[dict[str, Any]] = []
    tested_rules: list[dict[str, Any]] = []
    observation_counts: dict[str, int] = {}
    if split["split_feasible"]:
        by_split = {
            name: [
                row
                for row in quantitative_records
                if split["assignment"][str(row["spectrum_id"])] == name
            ]
            for name in ("discovery", "calibration", "test")
        }
        observations = {
            name: build_identity_observations(rows) for name, rows in by_split.items()
        }
        observation_counts = {
            name: len(rows) for name, rows in observations.items()
        }
        discovery_rules = discover_rules_with_policy(
            observations["discovery"],
            qc_policy=qc_policy,
            maximum_q_value=float(
                base_preregistration["discovery"]["benjamini_hochberg_q_maximum"]
            ),
            minimum_lr_positive=float(
                base_preregistration["discovery"]["minimum_exploratory_lr_positive"]
            ),
        )
        calibrated_rules = calibrate_rules(
            discovery_rules, observations["calibration"]
        )
        tested_rules = evaluate_sealed_rules(
            calibrated_rules, observations["test"]
        )
    review = export_review_candidates(tested_rules)
    ce_models: list[dict[str, Any]] = []
    policy_summary = {
        "schema_version": "3.1.0",
        "frozen_corpus_record_count": line_count,
        "base_positive_lr_eligible_record_count": base_positive,
        "base_negative_lr_eligible_record_count": base_negative,
        "v3_1_positive_lr_eligible_record_count": amended_positive,
        "v3_1_negative_lr_eligible_record_count": amended_negative,
        "collision_energy_mode_counts": dict(sorted(ce_mode_counts.items())),
        "ineligibility_reason_counts": dict(sorted(reason_counts.items())),
        "nce_to_ev_conversion_count": nce_to_ev_conversions,
        "spectral_qc_policy_id": qc_policy["policy_id"],
        "minimum_peak_gate_active": qc_policy["peak_presence"]
        ["minimum_deduplicated_peak_count_gate"]["active"],
        "precursor_purity_threshold_status": qc_policy["precursor_purity"]
        ["threshold_status"],
    }

    _write_json(output_dir / "discovery_rules.json", discovery_rules)
    _write_json(output_dir / "calibration_model.json", calibrated_rules)
    _write_json(output_dir / "sealed_test_rules.json", tested_rules)
    _write_json(output_dir / "ce_models.json", ce_models)
    _write_json(output_dir / "review_candidates.json", review)
    _write_json(output_dir / "domain_feasibility_report.json", feasibility)
    _write_json(output_dir / "policy_summary.json", policy_summary)
    _write_json(
        output_dir / "runtime_export_gate.json",
        {
            "schema_version": "3.1.0",
            "status": "fail_closed_phase_1_amendment",
            "runtime_exported_rule_count": 0,
            "runtime_catalog_modified": False,
            "runtime_gate_modified": False,
            "runtime_ranking_enabled": False,
            "score_contribution": 0.0,
            "review_candidate_count": review["review_candidate_count"],
            "phase_2_allowed": False,
            "next_state": "waiting_user_review",
        },
    )
    _write_gzip_jsonl(
        output_dir / "mirror_clusters.jsonl.gz", mirror["clusters"]
    )
    _write_gzip_jsonl(
        output_dir / "split_manifest.jsonl.gz",
        (
            {"spectrum_id": spectrum_id, "split": split_name}
            for spectrum_id, split_name in sorted(split.get("assignment", {}).items())
        ),
    )

    def manifest_rows() -> Iterable[dict[str, Any]]:
        for base_record in _iter_jsonl(frozen_corpus_path):
            revised = apply_v31_policy(base_record, qc_policy=qc_policy)
            spectrum_id = str(revised["spectrum_id"])
            revised["independence"]["mirror_cluster_id"] = mirror["assignment"][
                spectrum_id
            ]
            revised.pop("peaks", None)
            revised["split"] = split.get("assignment", {}).get(
                spectrum_id, "not_quantitatively_eligible"
            )
            yield revised

    _write_gzip_jsonl(output_dir / "record_manifest.jsonl.gz", manifest_rows())
    artifact_sha256 = {
        name: file_sha256(output_dir / name) for name in ARTIFACT_NAMES
    }
    runtime_gate_path = root / "config/feature_fragments/runtime_gate_v2.json"
    runtime_gate = _json(runtime_gate_path)
    config_sha256 = {
        name: file_sha256(config_dir / name) for name in V31_CONFIG_NAMES
    }
    base_after = _validate_base_release(
        root, amendment, formal_receipt, benchmark_receipt
    )
    runtime_fail_closed = _runtime_fail_closed(runtime_gate)
    benchmark_state = (
        "infeasible_for_frozen_release"
        if amended_positive == 0
        else "candidates_waiting_review"
        if review["review_candidate_count"]
        else "failed_closed"
    )
    receipt = {
        "schema_version": "3.1.0",
        "receipt_id": "v15-feature-fragment-benchmark-v3.1-20260827",
        "generated_on": "2026-08-27",
        "status": benchmark_state,
        "runtime_state": "fail_closed",
        "next_state": "waiting_user_review",
        "input": {
            "frozen_corpus_path": str(frozen_corpus_path),
            "frozen_corpus_sha256": corpus_sha256,
            "record_count": line_count,
            "repository_release": benchmark_receipt["inputs"][
                "repository_release"
            ],
            "new_public_spectra_integrated": source_amendments[
                "new_public_spectra_integrated_in_this_amendment"
            ],
        },
        "base_v3_freeze": {
            "before": base_before,
            "after": base_after,
            "unchanged": base_before == base_after,
        },
        "configuration_sha256": config_sha256,
        "imported_dependency_sha256": _dependency_hashes(root),
        "policy": policy_summary,
        "mirrors": {
            key: mirror[key]
            for key in (
                "exact_edge_count_by_type",
                "unique_exact_edge_count",
                "near_duplicate_edge_count",
                "component_count",
                "component_size_counts",
                "near_calibration",
                "near_candidate_status",
            )
        },
        "feasibility": {
            "status_counts": feasibility["status_counts"],
            "domain_count": feasibility["domain_count"],
            "unresolved_domain_record_count": feasibility[
                "unresolved_domain_audit_bucket"
            ]["record_count"],
            "permanent_infeasibility_claim_allowed": False,
            "partial_pooling_may_rescue_failed_exact_domain": False,
        },
        "split": split,
        "identity_observation_counts": observation_counts,
        "collision_energy_models": {
            "model_count": len(ce_models),
            "native_units_retained": True,
            "nce_to_ev_conversion_applied": False,
            "cross_unit_conversion_applied": False,
        },
        "rules": {
            "discovered": len(discovery_rules),
            "calibrated": len(calibrated_rules),
            "sealed_test_evaluated": len(tested_rules),
            "review_candidate_count": review["review_candidate_count"],
            "runtime_exported_rule_count": 0,
            "identity_upgrade_allowed": False,
            "formula_upgrade_allowed": False,
            "isomer_resolution_allowed": False,
            "standard_confirmation_allowed": False,
        },
        "runtime": {
            "fail_closed_verified": runtime_fail_closed,
            "runtime_gate_path": str(runtime_gate_path),
            "runtime_gate_sha256": file_sha256(runtime_gate_path),
            "runtime_exported_rule_count": 0,
            "runtime_ranking_enabled": False,
            "score_contribution": 0.0,
            "runtime_files_modified": False,
        },
        "license_and_source": {
            "mona_status": source_amendments["sources"][0][
                "integration_status"
            ],
            "mona_aggregator_counts_as_source_or_laboratory": False,
            "per_spectrum_origin_and_license_required": True,
            "new_public_spectra_integrated": 0,
        },
        "artifact_sha256": artifact_sha256,
        "limitations": [
            "The feasibility result applies only to the frozen public release.",
            "No record is promoted from a raw reference spectrum to a Mssystem-confirmed standard.",
            "No database identifier is used as a compound name.",
            "No NCE-to-eV conversion is produced.",
            "No uncalibrated near-mirror edge is applied.",
            "Any later public-spectrum integration requires a separately frozen source release, origin trace, and per-spectrum license audit."
        ],
    }
    _write_json(output_dir / "benchmark_receipt.json", receipt)
    return receipt
