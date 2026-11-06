"""End-to-end discovery/calibration/sealed-test benchmark workflow."""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .corpus import canonical_bytes, file_sha256
from .mirror import cluster_mirror_records
from .split import build_isolated_split
from .statistics import (
    build_identity_observations,
    calibrate_rules,
    discover_rules,
    evaluate_ce_dependence,
    test_frozen_rules,
)


def _write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_gzip_jsonl(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            for row in rows:
                compressed.write(canonical_bytes(row) + b"\n")


def _rules_hash(rules: Sequence[Mapping[str, Any]]) -> str:
    frozen = [
        {
            "rule_id": rule["rule_id"],
            "feature_id": rule["feature_id"],
            "stratum_id": rule["stratum_id"],
            "mass_tolerance": rule["mass_tolerance"],
            "minimum_relative_intensity": rule["minimum_relative_intensity"],
            "independent_fragmentation_group": rule["independent_fragmentation_group"],
            "selected_for_sealed_test": rule.get("selected_for_sealed_test"),
            "frozen_parameter_sha256": rule.get("frozen_parameter_sha256"),
        }
        for rule in sorted(rules, key=lambda row: row["rule_id"])
    ]
    return hashlib.sha256(canonical_bytes(frozen)).hexdigest()


def run_benchmark(
    records: Sequence[dict[str, Any]],
    *,
    output_dir: str | Path,
    corpus_jsonl_sha256: str,
    seed: str = "v15-feature-fragment-benchmark-v2",
    discovery_parameters: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Run the three-way benchmark and write auditable derived artifacts."""
    if not records:
        raise ValueError("benchmark requires at least one normalized spectrum")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    discovery_parameters = dict(discovery_parameters or {})

    mirror = cluster_mirror_records(records, allow_near_candidates=False)
    for record in records:
        record["independence"]["mirror_cluster_id"] = mirror["mapping"][
            record["spectrum_id"]
        ]
    split = build_isolated_split(records, seed=seed)
    if not split["split_feasible"]:
        raise ValueError("identity/formula/mirror-isolated three-way split is infeasible")
    assignment = split["assignment"]

    split_rows = [
        {"spectrum_id": spectrum_id, **row}
        for spectrum_id, row in sorted(assignment.items())
    ]
    _write_gzip_jsonl(output_dir / "split_manifest.jsonl.gz", split_rows)
    cluster_members: dict[str, list[str]] = defaultdict(list)
    for spectrum_id, cluster_id in mirror["mapping"].items():
        cluster_members[cluster_id].append(spectrum_id)
    mirror_rows = [
        {
            "mirror_cluster_id": cluster_id,
            "member_count": len(members),
            "member_spectrum_ids": sorted(members),
            "link_basis": "exact_order_and_intensity_scale_invariant_peak_hash",
        }
        for cluster_id, members in sorted(cluster_members.items())
        if len(members) > 1
    ]
    _write_gzip_jsonl(output_dir / "mirror_clusters.jsonl.gz", mirror_rows)

    record_manifest_rows = []
    for record in sorted(records, key=lambda row: row["spectrum_id"]):
        record_manifest_rows.append(
            {
                "spectrum_id": record["spectrum_id"],
                "evidence_stage": record["evidence_stage"],
                "query_confirmation_status": record["query_confirmation_status"],
                "source": record["source"],
                "compound": {
                    "display_name": record["compound"]["display_name"],
                    "name_provenance": record["compound"]["name_provenance"],
                    "inchikey": record["compound"]["inchikey"],
                    "connectivity": record["compound"]["connectivity"],
                    "formula": record["compound"]["formula"],
                    "chemical_class": record["compound"]["chemical_class"],
                },
                "acquisition": record["acquisition"],
                "independence": record["independence"],
                "quality": record["quality"],
                "split": assignment[record["spectrum_id"]],
            }
        )
    _write_gzip_jsonl(output_dir / "record_manifest.jsonl.gz", record_manifest_rows)

    records_by_split: dict[str, list[dict[str, Any]]] = {
        "discovery": [],
        "calibration": [],
        "test": [],
    }
    for record in records:
        records_by_split[assignment[record["spectrum_id"]]["split"]].append(record)
    observations = {
        split_name: build_identity_observations(rows)
        for split_name, rows in records_by_split.items()
    }
    discovered = discover_rules(observations["discovery"], **discovery_parameters)
    calibrated = calibrate_rules(discovered, observations["calibration"])
    frozen_rule_hash = _rules_hash(calibrated)
    tested = test_frozen_rules(calibrated, observations["test"])
    ce_dependence = evaluate_ce_dependence(tested, observations["test"])

    _write_json(output_dir / "discovery_rules.json", discovered)
    _write_json(output_dir / "calibration_model.json", calibrated)
    _write_json(output_dir / "sealed_test_rules.json", tested)
    _write_json(output_dir / "ce_dependence.json", ce_dependence)

    zero_upgrade_violations = sum(
        bool(rule.get("identity_upgrade_allowed"))
        or bool(rule.get("formula_upgrade_allowed"))
        or bool(rule.get("isomer_resolution_allowed"))
        or bool(rule.get("standard_confirmation_allowed"))
        for rule in tested
    )
    output_files = [
        "split_manifest.jsonl.gz",
        "mirror_clusters.jsonl.gz",
        "record_manifest.jsonl.gz",
        "discovery_rules.json",
        "calibration_model.json",
        "sealed_test_rules.json",
        "ce_dependence.json",
    ]
    receipt = {
        "schema_version": "2.0.0",
        "status": "sealed_test_completed",
        "scientific_scope": "public_repository_experimental_msms_spectra",
        "runtime_effect": "none_shadow_benchmark_only",
        "input": {
            "normalized_spectrum_count": len(records),
            "corpus_jsonl_sha256": corpus_jsonl_sha256,
            "source_count": len({record["source"]["source_id"] for record in records}),
            "laboratory_count": len(
                {record["source"]["laboratory_id"] for record in records}
            ),
            "instrument_domain_count": len(
                {
                    record["acquisition"]["instrument_domain_id"]
                    for record in records
                }
            ),
            "two_dimensional_identity_count": len(
                {record["compound"]["connectivity"] for record in records}
            ),
            "formula_family_count": len(
                {record["compound"]["formula"] for record in records}
            ),
        },
        "mirror_audit": {
            key: value for key, value in mirror.items() if key != "mapping"
        },
        "split": {
            key: value
            for key, value in split.items()
            if key not in {"assignment", "components"}
        },
        "identity_observation_counts": {
            key: len(value) for key, value in observations.items()
        },
        "rule_counts": {
            "discovered": len(discovered),
            "calibrated_total": len(calibrated),
            "selected_for_sealed_test": sum(
                bool(rule["selected_for_sealed_test"]) for rule in calibrated
            ),
            "sealed_test_evaluated": len(tested),
            "sealed_test_evidence_grades": dict(
                sorted(Counter(rule["evidence_grade"] for rule in tested).items())
            ),
        },
        "ce_dependence": {
            "series_count": len(ce_dependence),
            "status_counts": dict(
                sorted(Counter(row["status"] for row in ce_dependence).items())
            ),
            "minimum_same_unit_fixed_ce_settings": 3,
            "minimum_shared_test_identities": 20,
            "cross_unit_conversion_applied": False,
        },
        "sealed_test_guard": {
            "sealed_test_spectrum_ids_sha256": split[
                "sealed_test_spectrum_ids_sha256"
            ],
            "frozen_rule_set_sha256_before_test": frozen_rule_hash,
            "test_did_not_create_rules": len(tested)
            <= sum(bool(rule["selected_for_sealed_test"]) for rule in calibrated),
        },
        "evidence_stage_guard": {
            "raw_input_stage": "raw_msms_spectrum",
            "output_stage": "class_level_msms_support_only",
            "identity_formula_isomer_standard_upgrade_violation_count": zero_upgrade_violations,
        },
        "limitations": [
            "source contributors overlap discovery, calibration, and test; source holdout is not claimed",
            "near-mirror thresholds lack a labelled cross-repository duplicate calibration set and were not applied",
            "records without a reported scan mass range contribute only to explicitly labelled exploratory absence denominators",
            "no tested rule may distinguish an exact identity or an isomer",
        ],
        "artifact_sha256": {
            name: file_sha256(output_dir / name) for name in output_files
        },
    }
    _write_json(output_dir / "benchmark_receipt.json", receipt)
    return receipt
