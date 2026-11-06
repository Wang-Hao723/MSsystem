from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from eval import v11_external_cohort as cohort
from eval import v11_external_validation as external
from matching.chemical_identity import rdkit_identity_profile


SALT = "v11-external-cohort-unit-test-sealing-salt-20260808"


def _policy(*, minimum_identities: int = 8) -> dict:
    policy = copy.deepcopy(external.load_external_validation_policy())
    policy["statistics"]["minimum_independent_identities"] = minimum_identities
    policy["statistics"]["paired_bootstrap_iterations"] = 500
    return policy


def _boundaries() -> dict:
    return {
        "connectivity": set(),
        "exact_structures": set(),
        "scaffolds": set(),
        "registry_name_formula": set(),
        "development_truth_manifest_sha256": "a" * 64,
        "identity_registry_sha256": "b" * 64,
    }


def _patch_external_contract(monkeypatch: pytest.MonkeyPatch, policy: dict) -> None:
    contract = {
        field: policy["flagship_contract"][field]
        for field in (
            "algorithm_id",
            "policy_semantic_sha256",
            "identity_registry_sha256",
        )
    }
    monkeypatch.setattr(
        external,
        "validate_frozen_contract",
        lambda selected_policy, *, project_root: contract,
    )
    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: _boundaries(),
    )


def _standards(count: int) -> list[dict]:
    rows = []
    for index in range(1, count + 1):
        smiles = "C" * (index + 1) + "O"
        profile = rdkit_identity_profile(smiles)
        assert profile["valid"]
        rows.append(
            {
                "standard_id": f"external-standard-{index:03d}",
                "validation_track": "novel_identity",
                "detection_method_id": (
                    "agilent-qtof-cid-pos-v1"
                    if index % 2
                    else "agilent-qtof-cid-neg-v1"
                ),
                "compound_name": f"External homolog {index}",
                "formula": profile["structure_formula"],
                "inchikey": profile["inchikey"],
                "smiles": smiles,
                "truth_status": "authentic_standard",
                "truth_provenance": "certified_reference_material",
                "standard_lot": f"lot-{index:03d}",
                "chemical_class": ("class-a", "class-b", "class-c")[index % 3],
                "instrument_domain": "Agilent-6546-QTOF-CID",
                "polarity": "POS" if index % 2 else "NEG",
                "adduct": "[M+H]+" if index % 2 else "[M-H]-",
                "acquisition_batch": "external-batch-01",
            }
        )
    return rows


def _schedule(
    monkeypatch: pytest.MonkeyPatch,
    *,
    count: int = 12,
    seed: int = 20260808,
    repeat_fraction: float = 0.2,
) -> tuple[dict, dict]:
    policy = _policy(minimum_identities=min(8, count))
    _patch_external_contract(monkeypatch, policy)
    schedule = cohort.create_randomized_schedule(
        _standards(count),
        cohort_id="external-confirmatory-01",
        seed=seed,
        policy=policy,
        technical_repeat_fraction=repeat_fraction,
    )
    return policy, schedule


def _write_sidecars(
    tmp_path: Path,
    schedule: dict,
    policy: dict,
) -> list[dict]:
    method_cases = [
        row
        for row in schedule["records"]
        if row["analysis_role"] != external.REPEAT_ROLE
    ]
    truth_keys = {
        row["method_case_id"]: rdkit_identity_profile(row["smiles"])["connectivity_key"]
        for row in method_cases
    }
    ordered_truths = [truth_keys[row["method_case_id"]] for row in method_cases]
    contract = policy["flagship_contract"]
    bindings = []
    for index, row in enumerate(schedule["records"]):
        task_id = f"external-task-{index:03d}"
        entity_id = f"external-entity-{index:03d}"
        truth = truth_keys[row["method_case_id"]]
        wrong = ordered_truths[(ordered_truths.index(truth) + 1) % len(ordered_truths)]
        payload = {
            "metadata": {
                "task_id": task_id,
                "algorithm_id": contract["algorithm_id"],
                "policy_sha256": contract["policy_semantic_sha256"],
                "identity_registry_sha256": contract["identity_registry_sha256"],
                "source_semantic_sha256": f"{index + 1:064x}",
            },
            "entities": [
                {
                    "entity_id": entity_id,
                    "scan_top1_audit_baseline": {"identity_key": f"std2d:{wrong}"},
                    "winner": {"identity_key": f"std2d:{truth}"},
                    "decision": {
                        "status": "selected",
                        "scope": "2d_connectivity_candidate",
                    },
                }
            ],
        }
        path = tmp_path / f"{task_id}_v11_flagship.json"
        path.write_text(
            json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )
        bindings.append(
            {
                "record_id": row["record_id"],
                "source_task_id": task_id,
                "source_entity_id": entity_id,
                "flagship_artifact_path": str(path),
            }
        )
    return bindings


def test_seeded_schedule_is_reproducible_and_seed_sensitive(monkeypatch):
    policy = _policy()
    _patch_external_contract(monkeypatch, policy)
    standards = _standards(12)
    first = cohort.create_randomized_schedule(
        standards,
        cohort_id="external-confirmatory-01",
        seed=11,
        policy=policy,
        technical_repeat_fraction=0.2,
    )
    second = cohort.create_randomized_schedule(
        standards,
        cohort_id="external-confirmatory-01",
        seed=11,
        policy=policy,
        technical_repeat_fraction=0.2,
    )
    changed = cohort.create_randomized_schedule(
        standards,
        cohort_id="external-confirmatory-01",
        seed=12,
        policy=policy,
        technical_repeat_fraction=0.2,
    )

    assert first["records"] == second["records"]
    assert (
        first["randomization_schedule_sha256"]
        == second["randomization_schedule_sha256"]
    )
    assert first["records"] != changed["records"]
    assert (
        first["randomization_schedule_sha256"]
        != changed["randomization_schedule_sha256"]
    )


def test_schedule_has_one_primary_per_identity_and_spaced_secondary_repeats(
    monkeypatch,
):
    policy, schedule = _schedule(monkeypatch)
    records = schedule["records"]
    primaries = [
        row for row in records if row["analysis_role"] == external.PRIMARY_ROLE
    ]
    repeats = [row for row in records if row["analysis_role"] == external.REPEAT_ROLE]

    assert len(primaries) == 12
    assert len({row["standard_id"] for row in primaries}) == 12
    assert len(repeats) == 3
    assert schedule["design_summary"]["independent_primary_identities"] == 12
    assert schedule["design_summary"]["technical_repeat_records"] == 3
    assert (
        schedule["design_summary"]["minimum_observed_technical_repeat_run_separation"]
        >= 5
    )
    assert schedule["design_summary"]["maximum_consecutive_same_chemical_class"] < len(
        records
    )
    assert (
        cohort.validate_schedule_artifact(schedule, policy=policy)["status"] == "passed"
    )


def test_schedule_hash_rejects_truth_role_and_run_order_tampering(monkeypatch):
    policy, schedule = _schedule(monkeypatch)
    primary_index = next(
        index
        for index, row in enumerate(schedule["records"])
        if row["analysis_role"] == external.PRIMARY_ROLE
    )
    for index, field, replacement in (
        (0, "chemical_class", "post-hoc-class"),
        (0, "detection_method_id", "post-hoc-method"),
        (primary_index, "analysis_role", "technical_repeat"),
        (0, "run_order", 999),
    ):
        tampered = copy.deepcopy(schedule)
        tampered["records"][index][field] = replacement
        with pytest.raises(ValueError, match="schedule SHA-256 mismatch"):
            cohort.validate_schedule_artifact(tampered, policy=policy)

    relabeled_track = copy.deepcopy(schedule)
    relabeled_track["records"][primary_index]["validation_track"] = "method_transfer"
    relabeled_track["records"][primary_index]["analysis_role"] = "method_variant"
    with pytest.raises(ValueError, match="schedule SHA-256 mismatch"):
        cohort.validate_schedule_artifact(relabeled_track, policy=policy)

    relabeled_fraction = copy.deepcopy(schedule)
    for row in relabeled_fraction["records"]:
        row["technical_repeat_fraction_requested"] = 0.0
    relabeled_fraction["technical_repeat_fraction_requested"] = 0.0
    replacement_hash = external.compute_randomization_schedule_sha256(
        relabeled_fraction["records"],
        cohort_id=relabeled_fraction["cohort_id"],
        policy=policy,
    )
    relabeled_fraction["randomization_schedule_sha256"] = replacement_hash
    for row in relabeled_fraction["records"]:
        row["randomization_schedule_sha256"] = replacement_hash
    with pytest.raises(ValueError, match="technical-repeat count differs"):
        cohort.validate_schedule_artifact(relabeled_fraction, policy=policy)


def test_schedule_rejects_duplicate_identity_method_cell_and_novel_primary_overlap(
    monkeypatch,
):
    policy = _policy()
    _patch_external_contract(monkeypatch, policy)
    duplicated = _standards(8)
    duplicated[1].update(
        {field: duplicated[0][field] for field in ("formula", "inchikey", "smiles")}
    )
    duplicated[1]["detection_method_id"] = duplicated[0]["detection_method_id"]
    with pytest.raises(ValueError, match="duplicate identity-detection-method cells"):
        cohort.create_randomized_schedule(
            duplicated,
            cohort_id="duplicate-connectivity",
            seed=1,
            policy=policy,
            technical_repeat_fraction=0,
        )

    overlap = rdkit_identity_profile(_standards(1)[0]["smiles"])["connectivity_key"]
    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: {
            **_boundaries(),
            "connectivity": {overlap},
        },
    )
    with pytest.raises(ValueError, match="novel-identity primary overlaps"):
        cohort.create_randomized_schedule(
            _standards(8),
            cohort_id="development-overlap",
            seed=1,
            policy=policy,
            technical_repeat_fraction=0,
        )


def test_same_compound_different_methods_are_retained_and_known_identity_is_allowed(
    monkeypatch,
):
    policy = _policy()
    _patch_external_contract(monkeypatch, policy)
    standards = _standards(8)
    second_method = copy.deepcopy(standards[0])
    second_method.update(
        {
            "validation_track": "method_transfer",
            "detection_method_id": "orbitrap-hcd-pos-v1",
            "instrument_domain": "Thermo-Orbitrap-HCD",
            "acquisition_batch": "external-batch-02",
        }
    )
    schedule = cohort.create_randomized_schedule(
        [*standards, second_method],
        cohort_id="same-compound-cross-method",
        seed=1,
        policy=policy,
        technical_repeat_fraction=0,
    )

    assert schedule["design_summary"]["independent_primary_identities"] == 8
    assert schedule["design_summary"]["method_variant_records"] == 1
    assert schedule["design_summary"]["planned_identity_method_cases"] == 9
    assert schedule["design_summary"]["identities_with_multiple_detection_methods"] == 1

    known = copy.deepcopy(standards[0])
    known.update(
        {
            "validation_track": "method_transfer",
            "detection_method_id": "orbitrap-hcd-pos-known-identity",
            "instrument_domain": "Thermo-Orbitrap-HCD",
        }
    )
    overlap = rdkit_identity_profile(known["smiles"])["connectivity_key"]
    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: {
            **_boundaries(),
            "connectivity": {overlap},
        },
    )
    known_schedule = cohort.create_randomized_schedule(
        [known],
        cohort_id="known-identity-new-method",
        seed=2,
        policy=policy,
        technical_repeat_fraction=0,
    )
    assert known_schedule["design_summary"]["independent_primary_identities"] == 0
    assert known_schedule["design_summary"]["method_variant_records"] == 1


def test_binding_is_explicit_complete_and_rejects_sidecar_reuse(
    monkeypatch,
    tmp_path,
):
    policy, schedule = _schedule(
        monkeypatch,
        count=8,
        repeat_fraction=0,
    )
    bindings = _write_sidecars(tmp_path, schedule, policy)

    with pytest.raises(ValueError, match="bindings are incomplete"):
        cohort.bind_flagship_artifacts(
            schedule,
            bindings[:-1],
            policy=policy,
        )

    empty_path = copy.deepcopy(bindings)
    empty_path[0]["flagship_artifact_path"] = ""
    with pytest.raises(ValueError, match="flagship_artifact_path is required"):
        cohort.bind_flagship_artifacts(
            schedule,
            empty_path,
            policy=policy,
        )

    reused = copy.deepcopy(bindings)
    reused[1].update(
        {
            "source_task_id": reused[0]["source_task_id"],
            "source_entity_id": reused[0]["source_entity_id"],
            "flagship_artifact_path": reused[0]["flagship_artifact_path"],
        }
    )
    with pytest.raises(ValueError, match="bound to multiple records"):
        cohort.bind_flagship_artifacts(
            schedule,
            reused,
            policy=policy,
        )


def test_plan_bind_freeze_lock_and_evaluate_is_confirmatory_end_to_end(
    monkeypatch,
    tmp_path,
):
    policy, schedule = _schedule(
        monkeypatch,
        count=8,
        repeat_fraction=0.125,
    )
    bindings = _write_sidecars(tmp_path, schedule, policy)
    prepared = cohort.prepare_bound_external_cohort(
        schedule,
        bindings,
        policy=policy,
        sealing_salt=SALT,
    )

    assert (
        prepared["public"]["cohort_summary"]["independent_connectivity_identities"] == 8
    )
    assert prepared["public"]["cohort_summary"]["novel_identity_primary_records"] == 8
    assert prepared["public"]["cohort_summary"]["technical_repeat_records"] == 1
    rendered_public = json.dumps(prepared["public"], ensure_ascii=False)
    assert "compound_name" not in rendered_public
    assert "standard_id" not in rendered_public
    assert set(prepared["artifact_map"]["records"][0]) == {
        "blind_code",
        "flagship_artifact_path",
    }

    locked = external.lock_predictions_from_flagship_artifacts(
        prepared["public"],
        prepared["artifact_map"]["records"],
        policy=policy,
    )
    report = external.evaluate_external_validation(
        public=prepared["public"],
        sealed_truth=prepared["sealed_truth"],
        prediction_bundle=locked,
        policy=policy,
    )

    assert report["primary"]["independent_identities"] == 8
    assert report["primary"]["technical_repeat_records_excluded"] == 1
    assert report["primary"]["paired_outcomes"]["gained"] == 8
    assert report["primary"]["paired_outcomes"]["lost"] == 0
    assert report["confirmatory_accuracy_claim_allowed"] is True
    assert report["automatic_deployment_authorized"] is False
