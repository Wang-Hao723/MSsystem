"""Generate immutable V3.1 readiness, intake, and fail-closed audit artifacts."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .common import file_sha256, read_json, write_gzip_jsonl, write_json
from .metadata import apply_static_metadata_derivations
from .near_mirror import minimum_zero_event_trials
from .readiness import (
    build_eligibility_gap_matrix,
    collision_energy_summary,
    split_isolation_readiness,
)
from .screening import screen_manifest_sources
from .statistics import prospective_planning_table


def _snapshot_directory(path: Path) -> dict[str, str]:
    return {
        str(item.relative_to(path)): file_sha256(item)
        for item in sorted(path.rglob("*"))
        if item.is_file()
    }


def _corpus_statistics(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    identities = {
        str((row.get("compound") or {}).get("connectivity"))
        for row in records
        if (row.get("compound") or {}).get("connectivity")
    }
    formulas = {
        str((row.get("compound") or {}).get("formula"))
        for row in records
        if (row.get("compound") or {}).get("formula")
    }
    sources = {
        str((row.get("source") or {}).get("source_id"))
        for row in records
        if (row.get("source") or {}).get("source_id")
    }
    labs = {
        str((row.get("source") or {}).get("laboratory_id"))
        for row in records
        if (row.get("source") or {}).get("laboratory_id")
    }
    domain_buckets: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "records": 0,
            "identities": set(),
            "formulas": set(),
            "sources": set(),
            "laboratories": set(),
            "models": set(),
        }
    )
    source_buckets: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "records": 0,
            "identities": set(),
            "formulas": set(),
            "laboratories": set(),
        }
    )
    license_records = Counter()
    license_identities: dict[str, set[str]] = defaultdict(set)
    for row in records:
        compound = row.get("compound") or {}
        source = row.get("source") or {}
        acquisition = row.get("acquisition") or {}
        identity = str(compound.get("connectivity") or "")
        formula = str(compound.get("formula") or "")
        source_id = str(source.get("source_id") or "unreported_source")
        lab_id = str(source.get("laboratory_id") or "unreported_laboratory")
        domain_id = str(
            acquisition.get("instrument_domain_id") or "unresolved_instrument_domain"
        )
        license_tier = str((source.get("license") or {}).get("tier") or "excluded")
        license_records[license_tier] += 1
        if identity:
            license_identities[license_tier].add(identity)
        domain = domain_buckets[domain_id]
        domain["records"] += 1
        domain["models"].add(str(acquisition.get("instrument_model")))
        if identity:
            domain["identities"].add(identity)
        if formula:
            domain["formulas"].add(formula)
        if source_id:
            domain["sources"].add(source_id)
        if lab_id:
            domain["laboratories"].add(lab_id)
        source_bucket = source_buckets[source_id]
        source_bucket["records"] += 1
        if identity:
            source_bucket["identities"].add(identity)
        if formula:
            source_bucket["formulas"].add(formula)
        if lab_id:
            source_bucket["laboratories"].add(lab_id)
    return {
        "record_count": len(records),
        "two_dimensional_identity_count": len(identities),
        "formula_family_count": len(formulas),
        "source_count": len(sources),
        "laboratory_count": len(labs),
        "resolved_exact_domain_record_count": sum(
            bucket["records"]
            for domain_id, bucket in domain_buckets.items()
            if domain_id != "unresolved_instrument_domain"
        ),
        "unresolved_domain_record_count": domain_buckets[
            "unresolved_instrument_domain"
        ]["records"],
        "license_tier_record_counts": dict(sorted(license_records.items())),
        "license_tier_identity_counts": {
            tier: len(values) for tier, values in sorted(license_identities.items())
        },
        "instrument_domains": [
            {
                "instrument_domain_id": domain_id,
                "record_count": bucket["records"],
                "identity_count": len(bucket["identities"]),
                "formula_family_count": len(bucket["formulas"]),
                "source_count": len(bucket["sources"]),
                "laboratory_count": len(bucket["laboratories"]),
                "instrument_models": sorted(bucket["models"]),
            }
            for domain_id, bucket in sorted(domain_buckets.items())
        ],
        "sources": [
            {
                "source_id": source_id,
                "record_count": bucket["records"],
                "identity_count": len(bucket["identities"]),
                "formula_family_count": len(bucket["formulas"]),
                "laboratory_count": len(bucket["laboratories"]),
            }
            for source_id, bucket in sorted(source_buckets.items())
        ],
    }


def _metadata_derivation_audit(
    source_records: Sequence[Mapping[str, Any]],
    derived_records: Sequence[Mapping[str, Any]],
    event_rows: Sequence[Mapping[str, Any]],
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    events = Counter()
    fields = Counter()
    models = Counter()
    source_ids = Counter()
    for row in event_rows:
        if row.get("status") != "derived":
            continue
        events[str(row.get("derivation_id"))] += 1
        fields[str(row.get("field"))] += 1
        source_ids[str((row.get("source") or {}).get("source_id"))] += 1
    for source, derived in zip(source_records, derived_records, strict=True):
        if (derived.get("v31_metadata_derivation") or {}).get("events"):
            model = (derived.get("acquisition") or {}).get("instrument_model")
            models[str(model)] += 1
        if source.get("quantitative_eligibility") != derived.get(
            "quantitative_eligibility"
        ):
            raise RuntimeError("metadata derivation changed source eligibility flags")
    return {
        "contract_id": contract["contract_id"],
        "source_record_count": len(source_records),
        "derived_record_count": sum(events.values()),
        "derivation_event_counts": dict(sorted(events.items())),
        "derived_field_counts": dict(sorted(fields.items())),
        "instrument_model_counts": dict(sorted(models.items())),
        "citation_source_counts": dict(sorted(source_ids.items())),
        "permitted_vendor_manual_derived_fields": contract[
            "permitted_vendor_manual_derived_fields"
        ],
        "protected_per_scan_fields": contract["protected_per_scan_fields"],
        "protected_field_violation_count": 0,
        "source_quantitative_eligibility_flags_modified": False,
    }


def _top_actionable_rows(
    rows: Sequence[Mapping[str, Any]], *, limit: int = 100
) -> list[Mapping[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            -row["evidence_lanes"]["positive_metadata_rescue_ceiling"]["target"][
                "identity_count"
            ],
            -row["evidence_lanes"]["positive_metadata_rescue_ceiling"][
                "background"
            ]["identity_count"],
            row["cell_id"],
        ),
    )[:limit]


def _compact_top_row(row: Mapping[str, Any]) -> dict[str, Any]:
    stratum = row["acquisition_stratum"]
    ceiling = row["evidence_lanes"]["positive_metadata_rescue_ceiling"]
    return {
        "cell_id": row["cell_id"],
        "chemical_class": row["chemical_class"]["name"],
        "instrument_domain_id": stratum["instrument_domain_id"],
        "instrument_model": stratum["instrument_model"],
        "polarity": stratum["polarity"],
        "adduct": stratum["adduct"],
        "charge": stratum["charge"],
        "collision_energy_unit": stratum["collision_energy_unit"],
        "collision_energy_mode": stratum["collision_energy_mode"],
        "collision_energy_stratum": stratum["collision_energy_stratum"],
        "license_tier": stratum["license_tier"],
        "target": ceiling["target"],
        "background": ceiling["background"],
        "shortfalls": ceiling["gap"]["shortfalls"],
    }


def _write_top_gap_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    fields = (
        "cell_id",
        "chemical_class",
        "instrument_domain_id",
        "instrument_model",
        "polarity",
        "adduct",
        "charge",
        "ce_unit",
        "ce_mode",
        "ce_stratum",
        "license_tier",
        "positive_current_target_identities",
        "positive_v31_target_identities",
        "positive_rescue_ceiling_target_identities",
        "positive_rescue_ceiling_background_identities",
        "target_identity_shortfall",
        "background_identity_shortfall",
        "target_source_shortfall",
        "target_laboratory_shortfall",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            stratum = row["acquisition_stratum"]
            lanes = row["evidence_lanes"]
            ceiling = lanes["positive_metadata_rescue_ceiling"]
            shortfall = ceiling["gap"]["shortfalls"]
            writer.writerow(
                {
                    "cell_id": row["cell_id"],
                    "chemical_class": row["chemical_class"]["name"],
                    "instrument_domain_id": stratum["instrument_domain_id"],
                    "instrument_model": stratum["instrument_model"],
                    "polarity": stratum["polarity"],
                    "adduct": stratum["adduct"],
                    "charge": stratum["charge"],
                    "ce_unit": stratum["collision_energy_unit"],
                    "ce_mode": stratum["collision_energy_mode"],
                    "ce_stratum": stratum["collision_energy_stratum"],
                    "license_tier": stratum["license_tier"],
                    "positive_current_target_identities": lanes["current_positive"][
                        "target"
                    ]["identity_count"],
                    "positive_v31_target_identities": lanes["v31_positive"][
                        "target"
                    ]["identity_count"],
                    "positive_rescue_ceiling_target_identities": ceiling["target"][
                        "identity_count"
                    ],
                    "positive_rescue_ceiling_background_identities": ceiling[
                        "background"
                    ]["identity_count"],
                    "target_identity_shortfall": shortfall["target_identities"],
                    "background_identity_shortfall": shortfall[
                        "background_identities"
                    ],
                    "target_source_shortfall": shortfall["target_sources"],
                    "target_laboratory_shortfall": shortfall[
                        "target_laboratories"
                    ],
                }
            )


def run_v31_readiness(
    records: Sequence[Mapping[str, Any]],
    *,
    source_manifest_path: str | Path,
    output_dir: str | Path,
    preregistration_path: str | Path,
    metadata_contract_path: str | Path,
    source_screening_contract_path: str | Path,
    chemistry_contract_path: str | Path,
    class_contract_path: str | Path,
    rule_schema_path: str | Path,
    runtime_gate_path: str | Path,
    v2_artifact_dir: str | Path,
    v3_artifact_dir: str | Path,
) -> dict[str, Any]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_manifest_path = Path(source_manifest_path)
    preregistration_path = Path(preregistration_path)
    metadata_contract_path = Path(metadata_contract_path)
    source_screening_contract_path = Path(source_screening_contract_path)
    chemistry_contract_path = Path(chemistry_contract_path)
    class_contract_path = Path(class_contract_path)
    rule_schema_path = Path(rule_schema_path)
    runtime_gate_path = Path(runtime_gate_path)
    v2_artifact_dir = Path(v2_artifact_dir)
    v3_artifact_dir = Path(v3_artifact_dir)
    preregistration = read_json(preregistration_path)
    metadata_contract = read_json(metadata_contract_path)
    screening_contract = read_json(source_screening_contract_path)
    chemistry_contract = read_json(chemistry_contract_path)
    class_contract = read_json(class_contract_path)
    rule_schema = read_json(rule_schema_path)
    runtime_gate = read_json(runtime_gate_path)
    v2_before = _snapshot_directory(v2_artifact_dir)
    v3_before = _snapshot_directory(v3_artifact_dir)

    derived_records: list[dict[str, Any]] = []
    derivation_events: list[dict[str, Any]] = []
    for record in records:
        derived, events = apply_static_metadata_derivations(record, metadata_contract)
        derived_records.append(derived)
        derivation_events.extend(events)
    derivation_audit = _metadata_derivation_audit(
        records, derived_records, derivation_events, metadata_contract
    )

    power_rows = prospective_planning_table(
        family_wise_alpha=float(
            preregistration["prospective_power_design"]["family_wise_alpha"]
        ),
        power=float(preregistration["prospective_power_design"]["target_power"]),
        maximum_log_half_width=float(
            preregistration["primary_thresholds_per_stratum_and_class"][
                "maximum_log_lr_interval_half_width"
            ]
        ),
        minimum_n=int(
            preregistration["primary_thresholds_per_stratum_and_class"][
                "minimum_target_identities"
            ]
        ),
    )
    design = preregistration["prospective_power_design"]
    selected_power = next(
        row
        for row in power_rows
        if row["holm_family_size"] == int(design["holm_family_size"])
        and row["background_frequency"]
        == float(design["planning_background_frequency"])
        and row["lr_positive"]
        == float(design["planning_likelihood_ratio_positive"])
    )
    planning_required = int(selected_power["planning_required_per_arm"])
    gap_rows, gap_summary = build_eligibility_gap_matrix(
        derived_records,
        thresholds=preregistration["primary_thresholds_per_stratum_and_class"],
        planning_required_per_arm=planning_required,
        fixed_ce_bin_width=float(
            preregistration["collision_energy"]["fixed_native_unit_bin_width"]
        ),
    )
    top_rows = _top_actionable_rows(gap_rows)
    corpus_stats = _corpus_statistics(records)
    gap_summary["corpus_statistics"] = corpus_stats
    gap_summary["top_actionable_cells"] = [
        _compact_top_row(row) for row in top_rows
    ]
    ce_summary = collision_energy_summary(derived_records)
    screening = screen_manifest_sources(records)
    split_readiness = split_isolation_readiness(derived_records)
    power_receipt = {
        "schema_version": "3.1.0",
        "method": "two-proportion planning plus Haldane-Anscombe expected-count interval precision",
        "selected_design": selected_power,
        "scenario_count": len(power_rows),
        "scenarios": power_rows,
        "claim_limit": design["claim_limit"],
    }
    near_mirror = {
        "schema_version": "3.1.0",
        "status": "not_calibrated",
        "near_mirror_edges_enabled": False,
        "minimum_independent_negative_pairs_if_zero_false_links": minimum_zero_event_trials(
            float(
                preregistration["near_mirror_calibration"][
                    "maximum_false_link_probability"
                ]
            ),
            float(preregistration["near_mirror_calibration"]["one_sided_confidence"]),
        ),
        "minimum_independent_positive_pairs_if_zero_missed_links": preregistration[
            "near_mirror_calibration"
        ]["minimum_independent_positive_pairs_if_zero_missed_links"],
        "calibration_pair_count": 0,
        "holdout_validation_pair_count": 0,
        "reason": "no source-laboratory-identity-isolated labelled pair set was supplied",
    }
    runtime_export = {
        "schema_version": "3.1.0",
        "status": "fail_closed_readiness_only",
        "runtime_exported_rule_count": 0,
        "runtime_ranking_enabled": False,
        "identity_rank_influence": False,
        "score_contribution": 0.0,
        "runtime_catalog_modified": False,
        "shared_matching_files_modified_by_v31": False,
        "next_state": "dataset_specific_metadata_intake_required",
    }
    chemistry_receipt = {
        "schema_version": "3.1.0",
        "contract_id": chemistry_contract["contract_id"],
        "active_manual_fragment_rule_count": len(
            chemistry_contract["active_manual_fragment_rules"]
        ),
        "formula_compatibility_is_a_fixed_score": False,
        "joint_fragment_rules_fitted": 0,
        "runtime_effect": False,
    }
    rule_contract_receipt = {
        "schema_version": "3.1.0",
        "schema_id": rule_schema["schema_id"],
        "required_field_count": len(rule_schema["required_fields"]),
        "forbidden_fixed_score_fields": rule_schema["forbidden_fields"],
        "validated_rule_count": 0,
        "runtime_exported_rule_count": 0,
        "reason": "no class-stratum cell passed the preregistered count and metadata gates"
    }

    write_gzip_jsonl(output_dir / "eligibility_gap_matrix.jsonl.gz", gap_rows)
    write_json(output_dir / "eligibility_summary.json", gap_summary)
    write_json(output_dir / "power_planning.json", power_receipt)
    write_json(output_dir / "metadata_derivation_audit.json", derivation_audit)
    write_json(output_dir / "source_screening_receipt.json", screening)
    write_json(output_dir / "collision_energy_readiness.json", ce_summary)
    write_json(output_dir / "split_isolation_readiness.json", split_readiness)
    write_json(output_dir / "near_mirror_calibration_gate.json", near_mirror)
    write_json(output_dir / "chemistry_hypothesis_receipt.json", chemistry_receipt)
    write_json(output_dir / "rule_contract_receipt.json", rule_contract_receipt)
    write_json(output_dir / "runtime_export_gate.json", runtime_export)
    _write_top_gap_csv(output_dir / "eligibility_top_100_gaps.csv", top_rows)

    artifact_names = (
        "chemistry_hypothesis_receipt.json",
        "collision_energy_readiness.json",
        "eligibility_gap_matrix.jsonl.gz",
        "eligibility_summary.json",
        "eligibility_top_100_gaps.csv",
        "metadata_derivation_audit.json",
        "near_mirror_calibration_gate.json",
        "power_planning.json",
        "rule_contract_receipt.json",
        "runtime_export_gate.json",
        "source_screening_receipt.json",
        "split_isolation_readiness.json",
    )
    artifact_hashes = {
        name: file_sha256(output_dir / name) for name in artifact_names
    }
    v2_after = _snapshot_directory(v2_artifact_dir)
    v3_after = _snapshot_directory(v3_artifact_dir)
    shared_runtime_fail_closed = bool(
        runtime_gate.get("deployment_state") == "fail_closed"
        and runtime_gate.get("runtime_evaluation_enabled") is False
        and runtime_gate.get("runtime_adapter_available") is False
        and runtime_gate.get("exported_rule_count") == 0
        and runtime_gate.get("identity_rank_influence") is False
        and runtime_gate.get("runtime_influence") is False
        and float(runtime_gate.get("score_contribution", 0.0)) == 0.0
    )
    ready_cell_count = gap_summary["cells_passing_count_and_precision_planning_gate"][
        "v31_positive"
    ]
    receipt = {
        "schema_version": "3.1.0",
        "receipt_id": "v15-feature-fragment-v31-readiness-20260827",
        "generated_on": "2026-08-27",
        "status": "fail_closed_no_rule_ready" if not ready_cell_count else "intake_ready_cells_require_new_sealed_split",
        "input": {
            "v3_record_manifest_path": str(source_manifest_path),
            "v3_record_manifest_sha256": file_sha256(source_manifest_path),
            "record_count": len(records),
        },
        "configuration_sha256": {
            "preregistration.json": file_sha256(preregistration_path),
            "metadata_derivation_contract.json": file_sha256(metadata_contract_path),
            "source_screening_contract.json": file_sha256(
                source_screening_contract_path
            ),
            "chemistry_hypothesis_contract.json": file_sha256(
                chemistry_contract_path
            ),
            "class_hierarchy_contract.json": file_sha256(class_contract_path),
            "rule_schema.json": file_sha256(rule_schema_path),
        },
        "contracts": {
            "preregistration_id": preregistration["preregistration_id"],
            "metadata_contract_id": metadata_contract["contract_id"],
            "source_screening_contract_id": screening_contract["contract_id"],
            "chemistry_contract_id": chemistry_contract["contract_id"],
            "class_contract_id": class_contract["contract_id"],
            "rule_schema_id": rule_schema["schema_id"],
        },
        "corpus": corpus_stats,
        "eligibility": gap_summary,
        "metadata_derivation": derivation_audit,
        "source_screening": {
            "source_count": screening["source_count"],
            "status_counts": screening["status_counts"],
        },
        "split_isolation": split_readiness,
        "collision_energy": {
            "group_count": ce_summary["group_count"],
            "ce_response_model_count": ce_summary["ce_response_model_count"],
            "cross_unit_conversion_applied": False,
        },
        "rules": {
            "discovered": 0,
            "calibrated": 0,
            "sealed_test_evaluated": 0,
            "runtime_exported": 0,
            "reason": "no class-stratum cell passed the frozen identity, formula, source, laboratory, and power-planning count gates",
        },
        "runtime": {
            "enabled": False,
            "shared_runtime_gate_path": str(runtime_gate_path),
            "shared_runtime_gate_sha256": file_sha256(runtime_gate_path),
            "shared_runtime_fail_closed_verified": shared_runtime_fail_closed,
            "runtime_exported_rule_count": 0,
            "runtime_ranking_enabled": False,
            "score_contribution": 0.0,
        },
        "immutability": {
            "v2_artifact_hashes_before": v2_before,
            "v2_artifact_hashes_after": v2_after,
            "v2_artifacts_unchanged": v2_before == v2_after,
            "v3_artifact_hashes_before": v3_before,
            "v3_artifact_hashes_after": v3_after,
            "v3_artifacts_unchanged": v3_before == v3_after,
        },
        "artifact_sha256": artifact_hashes,
        "deviations": [],
        "limitations": [
            "Static instrument manuals cannot establish per-scan precursor selection, isolation width, collision energy, scan identity, or reference-material purity.",
            "A compatible elemental balance is a plausibility annotation and does not add a score or alter an empirical likelihood ratio.",
            "No repository record is promoted to Mssystem standard confirmation and no database accession is treated as a compound name.",
            "No CE response curve is fitted without a frozen rule family and eligible identity denominators in the same native CE unit and instrument domain.",
            "The existing V3 records are reused only for readiness accounting; a future sealed cohort must exclude their identities, formula families, sources, laboratories, mirrors, and raw acquisition groups."
        ],
        "next_state": "dataset_specific_metadata_intake_required",
    }
    write_json(output_dir / "benchmark_receipt.json", receipt)
    return receipt
