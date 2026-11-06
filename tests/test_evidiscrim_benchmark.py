from __future__ import annotations

import copy
from pathlib import Path

import pytest

from eval.evidiscrim_benchmark import (
    audit_identity_disjoint_split,
    build_whole_spectrum_fiora_prediction,
    evaluate_locked_blind_predictions,
    freeze_identity_disjoint_blind_set,
    load_benchmark_policy,
    lock_method_prediction,
    lock_prediction_bundle,
)
from matching.evidiscrim import build_evidiscrim_shadow, load_evidiscrim_policy
from matching.qtof_calibration import (
    apply_qtof_intensity_calibrator,
    fit_qtof_intensity_calibrator,
)


ROOT = Path(__file__).resolve().parents[1]


def _benchmark_policy() -> dict:
    policy = load_benchmark_policy(
        ROOT / "config" / "evidiscrim_benchmark_policy.json"
    )
    policy = copy.deepcopy(policy)
    policy["split"]["minimum_total_connectivity_identities"] = 3
    policy["split"]["minimum_blind_test_connectivity_identities"] = 1
    policy["statistics"]["paired_cluster_bootstrap_iterations"] = 300
    policy["qtof_calibration"]["minimum_connectivity_identities_per_domain"] = 2
    policy["qtof_calibration"]["minimum_predicted_peaks_per_domain"] = 4
    return policy


PAIRS = (
    ("CCO", "COC", "C2H6O"),
    ("CCCO", "CC(O)C", "C3H8O"),
    ("CCC=O", "CC(=O)C", "C3H6O"),
    ("CCCC", "CC(C)C", "C4H10"),
    ("CCCCC", "CC(C)CC", "C5H12"),
    ("CCCCCC", "CC(C)CCC", "C6H14"),
)


def _standards() -> list[dict]:
    rows = []
    for index, (truth, competitor, formula) in enumerate(PAIRS, start=1):
        rows.append({
            "standard_id": f"standard-{index}",
            "blind_code": f"blind-{index}",
            "replicate_id": f"replicate-{index}",
            "spectrum_id": f"scan={index}",
            "source_dataset": "independent-authentic-study",
            "source_record_id": f"source-{index}",
            "truth_provenance": "authenticated_reference_material",
            "authentic_standard_lot": f"lot-{index}",
            "true_candidate_id": f"truth-{index}",
            "true_smiles": truth,
            "true_inchikey": "",
            "formula": formula,
            "adduct": "[M+H]+",
            "polarity": "positive",
            "precursor_mz": 100.0 + index,
            "collision_energy": 20.0,
            "collision_energy_unit": "eV",
            "instrument_vendor": "Agilent",
            "instrument_model": "6545 Q-TOF",
            "analyzer_type": "Q-TOF",
            "activation_type": "CID",
            "acquisition_batch": f"batch-{1 + index % 2}",
            "run_order": index,
            "peaks": [
                [31.0 + index, 100.0],
                [45.0 + index, 40.0],
                [60.0 + index, 5.0],
            ],
            "candidates": [
                {
                    "candidate_id": f"truth-{index}",
                    "smiles": truth,
                    "formula": formula,
                    "inchikey": "",
                    "source": "frozen-candidate-generator",
                },
                {
                    "candidate_id": f"competitor-{index}",
                    "smiles": competitor,
                    "formula": formula,
                    "inchikey": "",
                    "source": "frozen-candidate-generator",
                },
            ],
        })
    return rows


def test_identity_disjoint_freeze_is_deterministic_and_seals_truth():
    policy = _benchmark_policy()
    first = freeze_identity_disjoint_blind_set(
        _standards(), policy=policy, sealing_salt="test-salt"
    )
    second = freeze_identity_disjoint_blind_set(
        _standards(), policy=policy, sealing_salt="test-salt"
    )

    assert first["status"] == "ready"
    assert first["leakage_audit"]["status"] == "passed"
    assert first["component_manifest"] == second["component_manifest"]
    assert (
        first["blind_test_public"]["blind_test_cases"]
        == second["blind_test_public"]["blind_test_cases"]
    )
    public_text = str(first["blind_test_public"])
    assert "true_candidate_id" not in public_text
    assert "standard_id" not in public_text
    assert first["sealed_labels"]["labels"]


def test_split_audit_detects_candidate_identity_leakage():
    policy = _benchmark_policy()
    bundle = freeze_identity_disjoint_blind_set(
        _standards(), policy=policy, sealing_salt="test-salt"
    )
    train = copy.deepcopy(bundle["train_records"])
    blind = copy.deepcopy(bundle["calibration_records"][:1])
    blind[0]["candidate_connectivity_keys"] = train[0][
        "candidate_connectivity_keys"
    ]
    audit = audit_identity_disjoint_split(
        {"train": train, "calibration": [], "blind_test": blind}, policy
    )
    assert audit["status"] == "failed"
    assert any(
        item["key_type"] == "candidate_connectivity"
        for item in audit["blocking_violations"]
    )


def test_stereoisomer_only_candidate_duplicate_is_rejected():
    policy = _benchmark_policy()
    row = _standards()[0]
    row = copy.deepcopy(row)
    row.update({
        "true_smiles": "F[C@H](Cl)Br",
        "formula": "CHBrClF",
    })
    row["candidates"] = [
        {
            "candidate_id": "left",
            "smiles": "F[C@H](Cl)Br",
            "formula": "CHBrClF",
            "inchikey": "",
            "source": "frozen",
        },
        {
            "candidate_id": "right",
            "smiles": "F[C@@H](Cl)Br",
            "formula": "CHBrClF",
            "inchikey": "",
            "source": "frozen",
        },
    ]
    row["true_candidate_id"] = "left"
    with pytest.raises(ValueError, match="duplicate-connectivity"):
        freeze_identity_disjoint_blind_set([row], policy=policy)


def _evidiscrim_input():
    left = {
        "compound_name": "ethanol",
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "smiles": "CCO",
        "formula": "C2H6O",
        "adduct": "[M+H]+",
        "rank_score": 0.8,
    }
    right = {
        "compound_name": "dimethyl ether",
        "inchikey": "LCGLNKUTAGEVQW-UHFFFAOYSA-N",
        "smiles": "COC",
        "formula": "C2H6O",
        "adduct": "[M+H]+",
        "rank_score": 0.7,
    }
    left["top_candidates"] = [dict(left), right]
    results = [{
        "spectrum_id": "scan=1",
        "precursor_mz": 47.0,
        "collision_energy": 20.0,
        "collision_energy_unit": "eV",
        "raw_peaks": [[31.0, 100.0], [29.0, 10.0]],
        "match": left,
    }]
    return results, left, right


def _prediction(identity: str, peaks: list[list[float]], *, two_models: bool):
    record = {
        "row_index": 0,
        "identity_key": f"ik::{identity}",
        "collision_energy": 20.0,
        "collision_energy_unit": "eV",
        "energy_transfer_qualified": True,
        "predicted_peaks": peaks,
        "provider_id": "model-a",
        "architecture_family": "architecture-a",
    }
    if two_models:
        record["prediction_members"] = [
            {
                "provider_id": "model-a",
                "architecture_family": "architecture-a",
                "predicted_peaks": peaks,
                "energy_transfer_qualified": True,
            },
            {
                "provider_id": "model-b",
                "architecture_family": "architecture-b",
                "predicted_peaks": peaks,
                "energy_transfer_qualified": True,
            },
        ]
    return record


def test_single_model_cannot_resolve_but_two_architectures_can():
    results, left, right = _evidiscrim_input()
    policy = load_evidiscrim_policy(ROOT / "config" / "evidiscrim_policy.json")
    singleton = build_evidiscrim_shadow(
        results,
        task_id="singleton",
        prediction_records=[
            _prediction(left["inchikey"], [[31.0, 1.0], [29.0, 0.1]], two_models=False),
            _prediction(right["inchikey"], [[45.0, 1.0], [29.0, 0.1]], two_models=False),
        ],
        policy=policy,
    )["spectra"][0]["pairwise_comparisons"][0]
    ensemble = build_evidiscrim_shadow(
        results,
        task_id="ensemble",
        prediction_records=[
            _prediction(left["inchikey"], [[31.0, 1.0], [29.0, 0.1]], two_models=True),
            _prediction(right["inchikey"], [[45.0, 1.0], [29.0, 0.1]], two_models=True),
        ],
        policy=policy,
    )["spectra"][0]["pairwise_comparisons"][0]

    assert singleton["status"] == "abstained"
    assert singleton["abstention_reason"] == (
        "insufficient_independent_architecture_families"
    )
    assert ensemble["status"] == "resolved_counterfactual"
    assert ensemble["ensemble_qualified"] is True
    assert ensemble["common_architecture_family_count"] == 2


def test_repeated_runs_of_one_architecture_do_not_fake_an_ensemble():
    results, left, right = _evidiscrim_input()
    policy = load_evidiscrim_policy(ROOT / "config" / "evidiscrim_policy.json")

    def repeated(identity: str, peaks: list[list[float]]) -> dict:
        record = _prediction(identity, peaks, two_models=True)
        for index, member in enumerate(record["prediction_members"]):
            member["provider_id"] = f"checkpoint-{index}"
            member["architecture_family"] = "same-architecture"
        return record

    pair = build_evidiscrim_shadow(
        results,
        task_id="same-family-repeats",
        prediction_records=[
            repeated(left["inchikey"], [[31.0, 1.0], [29.0, 0.1]]),
            repeated(right["inchikey"], [[45.0, 1.0], [29.0, 0.1]]),
        ],
        policy=policy,
    )["spectra"][0]["pairwise_comparisons"][0]

    assert pair["status"] == "abstained"
    assert pair["common_architecture_family_count"] == 1
    assert pair["abstention_reason"] == (
        "insufficient_independent_architecture_families"
    )


def test_architecture_disagreement_is_measured_and_penalizes_separation():
    results, left, right = _evidiscrim_input()
    policy = load_evidiscrim_policy(ROOT / "config" / "evidiscrim_policy.json")

    def record(identity: str, first: list[list[float]], second: list[list[float]]):
        return {
            "row_index": 0,
            "identity_key": f"ik::{identity}",
            "collision_energy": 20.0,
            "collision_energy_unit": "eV",
            "prediction_members": [
                {
                    "provider_id": "model-a",
                    "architecture_family": "architecture-a",
                    "predicted_peaks": first,
                    "energy_transfer_qualified": True,
                },
                {
                    "provider_id": "model-b",
                    "architecture_family": "architecture-b",
                    "predicted_peaks": second,
                    "energy_transfer_qualified": True,
                },
            ],
        }

    agreed = build_evidiscrim_shadow(
        results,
        task_id="agreed",
        prediction_records=[
            record(
                left["inchikey"],
                [[31.0, 1.0], [29.0, 0.1]],
                [[31.0, 1.0], [29.0, 0.1]],
            ),
            record(
                right["inchikey"],
                [[45.0, 1.0], [29.0, 0.1]],
                [[45.0, 1.0], [29.0, 0.1]],
            ),
        ],
        policy=policy,
    )["spectra"][0]["pairwise_comparisons"][0]["energy_separations"][0]
    disagreed = build_evidiscrim_shadow(
        results,
        task_id="disagreed",
        prediction_records=[
            record(
                left["inchikey"],
                [[31.0, 1.0], [29.0, 0.1]],
                [[31.0, 0.2], [45.0, 1.0]],
            ),
            record(
                right["inchikey"],
                [[45.0, 1.0], [29.0, 0.1]],
                [[31.0, 1.0], [45.0, 0.2]],
            ),
        ],
        policy=policy,
    )["spectra"][0]["pairwise_comparisons"][0]["energy_separations"][0]

    assert agreed["between_architecture_uncertainty"] == 0
    assert disagreed["between_architecture_uncertainty"] > 0
    assert disagreed["uncertainty_reliability_factor"] < agreed[
        "uncertainty_reliability_factor"
    ]
    assert disagreed["spectral_distribution_divergence"] < agreed[
        "spectral_distribution_divergence"
    ]


def test_qtof_calibration_preserves_base_architecture_and_exact_domain():
    policy = _benchmark_policy()
    examples = []
    for index, smiles in enumerate(("CCO", "CCCO"), start=1):
        examples.append({
            "split": "calibration",
            "true_smiles": smiles,
            "formula": "C2H6O" if index == 1 else "C3H8O",
            "true_inchikey": "",
            "spectrum_id": f"cal-{index}",
            "observed_peaks": [[31.0, 1.0], [45.0, 0.2]],
            "predicted_peaks": [[31.0, 1.0], [45.0, 0.8]],
            "precursor_mz": 100.0 + index,
            "instrument_vendor": "Agilent",
            "instrument_model": "6545 Q-TOF",
            "analyzer_type": "Q-TOF",
            "activation_type": "CID",
            "polarity": "positive",
            "adduct": "[M+H]+",
            "collision_energy": 20.0,
            "collision_energy_unit": "eV",
            "provider_id": "fiora",
            "base_architecture_family": "fiora-family",
        })
    calibrator = fit_qtof_intensity_calibrator(examples, policy=policy)
    assert calibrator["status"] == "fitted"
    assert calibrator["counts_as_independent_architecture"] is False
    calibrated = apply_qtof_intensity_calibrator(
        {
            "predicted_peaks": [[31.0, 1.0], [45.0, 0.8]],
            "provider_id": "fiora",
            "architecture_family": "fiora-family",
        },
        calibrator=calibrator,
        instrument_context=examples[0],
    )
    assert calibrated["architecture_family"] == "fiora-family"
    assert calibrated["counts_as_independent_architecture"] is False
    mismatch = dict(examples[0], instrument_model="6546 Q-TOF")
    with pytest.raises(ValueError, match="no exact QTOF calibration curve"):
        apply_qtof_intensity_calibrator(
            calibrated, calibrator=calibrator, instrument_context=mismatch
        )
    with pytest.raises(ValueError, match="base architecture does not match"):
        apply_qtof_intensity_calibrator(
            {
                "predicted_peaks": [[31.0, 1.0], [45.0, 0.8]],
                "provider_id": "different-model",
                "architecture_family": "different-family",
            },
            calibrator=calibrator,
            instrument_context=examples[0],
        )


def test_qtof_calibration_rejects_blind_identity_overlap():
    policy = _benchmark_policy()
    example = {
        "split": "calibration",
        "true_smiles": "CCO",
        "formula": "C2H6O",
        "true_inchikey": "",
        "spectrum_id": "cal-overlap",
        "observed_peaks": [[31.0, 1.0], [45.0, 0.2]],
        "predicted_peaks": [[31.0, 1.0], [45.0, 0.8]],
        "precursor_mz": 100.0,
        "instrument_vendor": "Agilent",
        "instrument_model": "6545 Q-TOF",
        "analyzer_type": "quadrupole time-of-flight",
        "activation_type": "CID",
        "polarity": "positive",
        "adduct": "[M+H]+",
        "collision_energy": 20.0,
        "collision_energy_unit": "eV",
        "provider_id": "fiora",
        "base_architecture_family": "fiora-family",
    }
    with pytest.raises(ValueError, match="overlaps the sealed blind test"):
        fit_qtof_intensity_calibrator(
            [example],
            policy=policy,
            blind_connectivity_keys=["LFQSCWFLJHTTHZ"],
        )


def _five_way_inputs():
    policy = _benchmark_policy()
    bundle = freeze_identity_disjoint_blind_set(
        _standards(), policy=policy, sealing_salt="test-salt"
    )
    predictions = []
    for case in bundle["blind_test_public"]["blind_test_cases"]:
        label = next(
            item for item in bundle["sealed_labels"]["labels"]
            if item["blind_code"] == case["blind_code"]
        )
        truth = label["true_candidate_id"]
        wrong = next(
            candidate["candidate_id"] for candidate in case["candidates"]
            if candidate["candidate_id"] != truth
        )
        methods = {
            "production_baseline": (wrong, ()),
            "whole_spectrum_fiora": (
                wrong, ("fiora_local_neighborhood_fragmentation_gnn",)
            ),
            "evidiscrim_full": (
                truth,
                (
                    "fiora_local_neighborhood_fragmentation_gnn",
                    "3dmolms_3d_conformation_gnn",
                ),
            ),
            "evidiscrim_no_uncertainty": (
                truth,
                (
                    "fiora_local_neighborhood_fragmentation_gnn",
                    "3dmolms_3d_conformation_gnn",
                ),
            ),
            "evidiscrim_no_abstention": (
                truth,
                (
                    "fiora_local_neighborhood_fragmentation_gnn",
                    "3dmolms_3d_conformation_gnn",
                ),
            ),
        }
        for method, (top1, families) in methods.items():
            predictions.append(lock_method_prediction({
                "blind_code": case["blind_code"],
                "method": method,
                "ranked_candidate_ids": [top1],
                "abstained": False,
                "released": False,
                "release_score": 0.8,
                "architecture_families": families,
                "architecture_model_artifact_sha256": {
                    family: ("f" if "fiora" in family else "3") * 64
                    for family in families
                },
            }))
    fdr = []
    for index in range(20):
        for method, false in (
            ("production_baseline", True),
            ("evidiscrim_full", False),
        ):
            fdr.append({
                "evaluation_unit_id": f"project-unit-{index}",
                "method": method,
                "released": True,
                "false_discovery": false,
                "source_artifact_sha256": "a" * 64,
            })
    prediction_bundle = lock_prediction_bundle(
        bundle["blind_test_public"], predictions
    )
    return policy, bundle, prediction_bundle, fdr


def test_five_way_gate_requires_top1_ci_and_nonworse_project_fdr():
    policy, bundle, predictions, fdr = _five_way_inputs()
    report = evaluate_locked_blind_predictions(
        public_blind_artifact=bundle["blind_test_public"],
        sealed_labels=bundle["sealed_labels"],
        predictions=predictions,
        project_fdr_observations=fdr,
        split_leakage_audit=bundle["leakage_audit"],
        policy=policy,
    )
    assert report["promotion_gate"]["production_discussion_eligible"] is True
    assert report["promotion_gate"]["automatic_deployment_authorized"] is False
    assert (
        report["paired_top1_comparisons"][
            "evidiscrim_full_minus_production_baseline"
        ]["ci_lower"] > 0
    )
    assert report["project_fdr_comparison"]["not_worse"] is True


def test_prediction_lock_tampering_is_rejected():
    policy, bundle, predictions, fdr = _five_way_inputs()
    predictions = copy.deepcopy(predictions)
    predictions["records"][0]["release_score"] = 0.81
    with pytest.raises(ValueError, match="prediction bundle lock"):
        evaluate_locked_blind_predictions(
            public_blind_artifact=bundle["blind_test_public"],
            sealed_labels=bundle["sealed_labels"],
            predictions=predictions,
            project_fdr_observations=fdr,
            split_leakage_audit=bundle["leakage_audit"],
            policy=policy,
        )


def test_worse_project_fdr_blocks_production_discussion():
    policy, bundle, predictions, fdr = _five_way_inputs()
    for observation in fdr:
        observation["false_discovery"] = (
            observation["method"] == "evidiscrim_full"
        )
    report = evaluate_locked_blind_predictions(
        public_blind_artifact=bundle["blind_test_public"],
        sealed_labels=bundle["sealed_labels"],
        predictions=predictions,
        project_fdr_observations=fdr,
        split_leakage_audit=bundle["leakage_audit"],
        policy=policy,
    )
    assert report["project_fdr_comparison"]["not_worse"] is False
    assert report["promotion_gate"]["production_discussion_eligible"] is False
    assert report["promotion_gate"]["decision"] == "remain_shadow"


def test_duplicate_project_fdr_units_are_rejected():
    policy, bundle, predictions, fdr = _five_way_inputs()
    fdr.append(copy.deepcopy(fdr[0]))
    with pytest.raises(ValueError, match="duplicate project FDR observation"):
        evaluate_locked_blind_predictions(
            public_blind_artifact=bundle["blind_test_public"],
            sealed_labels=bundle["sealed_labels"],
            predictions=predictions,
            project_fdr_observations=fdr,
            split_leakage_audit=bundle["leakage_audit"],
            policy=policy,
        )


def test_mismatched_project_fdr_source_cohort_blocks_gate():
    policy, bundle, predictions, fdr = _five_way_inputs()
    for observation in fdr:
        if observation["method"] == "evidiscrim_full":
            observation["source_artifact_sha256"] = "b" * 64
    report = evaluate_locked_blind_predictions(
        public_blind_artifact=bundle["blind_test_public"],
        sealed_labels=bundle["sealed_labels"],
        predictions=predictions,
        project_fdr_observations=fdr,
        split_leakage_audit=bundle["leakage_audit"],
        policy=policy,
    )
    assert report["project_fdr_comparison"]["status"] == (
        "unavailable_mismatched_project_fdr_source_cohort"
    )
    assert report["promotion_gate"]["production_discussion_eligible"] is False


def test_whole_spectrum_fiora_comparator_uses_matchms_hungarian():
    policy = _benchmark_policy()
    case = {
        "blind_code": "blind-x",
        "precursor_mz": 100.0,
        "peaks": [[31.0, 1.0], [45.0, 0.2]],
        "candidates": [
            {"candidate_id": "correct"},
            {"candidate_id": "wrong"},
        ],
    }
    prediction = build_whole_spectrum_fiora_prediction(
        case,
        {
            "correct": [[31.005, 1.0], [45.005, 0.2]],
            "wrong": [[70.0, 1.0], [80.0, 0.2]],
        },
        policy=policy,
        model_artifact_sha256="f" * 64,
    )
    assert prediction["ranked_candidate_ids"][0] == "correct"
    assert prediction["method"] == "whole_spectrum_fiora"
