from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

from matching.feature_fragment_benchmark.export import export_empirically_validated_rules
from matching.feature_fragment_benchmark.mirror import exact_peak_hash
from matching.feature_fragment_benchmark.split import build_isolated_split
from matching.feature_fragment_benchmark.workflow import run_benchmark


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_export_gate_rejects_provisional_rule_and_accepts_only_exact_domain_grade() -> None:
    base = {
        "rule_id": "FFV2-TEST",
        "status": "sealed_test_evaluated_shadow",
        "evidence_grade": "provisional_shadow",
        "empirical_validation_gate_unmet": ["test_sources_lt_3"],
        "evidence_stage_output": "class_level_msms_support_only",
        "identity_upgrade_allowed": False,
        "formula_upgrade_allowed": False,
        "isomer_resolution_allowed": False,
        "standard_confirmation_allowed": False,
        "frozen_parameter_sha256": "a" * 64,
    }
    rejected = export_empirically_validated_rules([base])
    assert rejected["exported_rule_count"] == 0
    assert rejected["integration_status"] == "no_rule_eligible_for_runtime_adapter"

    accepted_rule = {
        **base,
        "evidence_grade": "empirically_validated_exact_domain",
        "empirical_validation_gate_unmet": [],
    }
    accepted = export_empirically_validated_rules([accepted_rule])
    assert accepted["exported_rule_count"] == 1
    assert accepted["accepted_rules"][0]["rule_id"] == "FFV2-TEST"


def test_changing_only_sealed_test_peaks_leaves_discovery_and_calibration_frozen(
    tmp_path: Path, benchmark_record_factory
) -> None:
    records = [
        benchmark_record_factory(
            f"workflow-{index:03d}",
            identity=f"workflow-identity-{index:03d}",
            formula=f"C{index + 10}H{2 * index + 20}O2",
            mirror_cluster_id=f"workflow-mirror-{index:03d}",
            peaks=[[50.0, 1.0], [180.0, 0.5], [250.0 + index / 1000.0, 0.1]],
        )
        for index in range(100)
    ]
    initial_split = build_isolated_split(records, seed="workflow-test-seal")
    by_split: dict[str, list[dict]] = {
        "discovery": [],
        "calibration": [],
        "test": [],
    }
    for record in records:
        split_name = initial_split["assignment"][record["spectrum_id"]]["split"]
        by_split[split_name].append(record)

    for split_records in by_split.values():
        split_records.sort(key=lambda row: row["spectrum_id"])
        target_count = len(split_records) // 2
        for index, record in enumerate(split_records):
            is_target = index < target_count
            record["compound"]["chemical_class"]["leaf_class"] = (
                "Target" if is_target else "Background"
            )
            unique_mz = 250.0 + int(record["spectrum_id"].rsplit("-", 1)[1]) / 1000.0
            record["peaks"] = [
                [50.0, 1.0],
                [100.0 if is_target else 120.0, 0.8],
                [unique_mz, 0.1],
            ]
            record["independence"]["exact_peak_hash"] = exact_peak_hash(
                record["peaks"]
            )

    first_records = deepcopy(records)
    second_records = deepcopy(records)
    changed_test_id = by_split["test"][0]["spectrum_id"]
    changed = next(
        record for record in second_records if record["spectrum_id"] == changed_test_id
    )
    changed["peaks"] = [[50.0, 1.0], [120.0, 0.8], [333.333, 0.1]]
    changed["independence"]["exact_peak_hash"] = exact_peak_hash(changed["peaks"])

    common = {
        "seed": "workflow-test-seal",
        "discovery_parameters": {
            "minimum_target_identities": 10,
            "minimum_background_identities": 10,
            "minimum_target_present": 10,
            "maximum_q_value": 1.0,
            "minimum_lr_positive": 1.1,
        },
    }
    first_dir = tmp_path / "first"
    second_dir = tmp_path / "second"
    first = run_benchmark(
        first_records,
        output_dir=first_dir,
        corpus_jsonl_sha256="a" * 64,
        **common,
    )
    second = run_benchmark(
        second_records,
        output_dir=second_dir,
        corpus_jsonl_sha256="b" * 64,
        **common,
    )

    assert first["sealed_test_guard"]["frozen_rule_set_sha256_before_test"] == second[
        "sealed_test_guard"
    ]["frozen_rule_set_sha256_before_test"]
    assert _read_json(first_dir / "discovery_rules.json") == _read_json(
        second_dir / "discovery_rules.json"
    )
    assert _read_json(first_dir / "calibration_model.json") == _read_json(
        second_dir / "calibration_model.json"
    )
    assert _read_json(first_dir / "sealed_test_rules.json") != _read_json(
        second_dir / "sealed_test_rules.json"
    )
