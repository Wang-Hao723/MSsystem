from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from eval import v11_external_validation as external
from matching.chemical_identity import rdkit_identity_profile


ROOT = Path(__file__).resolve().parents[1]
SALT = "v11-external-validation-unit-test-salt-20260808"
SMILES = (
    "CCO",
    "CCCO",
    "CCCCO",
    "CCCCCO",
    "CCCCCCO",
    "CCCCCCCO",
    "CCCCCCCCO",
    "CCCCCCCCCO",
)


def _policy(*, minimum_identities: int = 3) -> dict:
    policy = copy.deepcopy(external.load_external_validation_policy())
    policy["statistics"]["minimum_independent_identities"] = minimum_identities
    policy["statistics"]["paired_bootstrap_iterations"] = 1000
    return policy


def _boundaries(
    *,
    connectivity: set[str] | None = None,
    exact_structures: set[str] | None = None,
    scaffolds: set[str] | None = None,
    registry_name_formula: set[tuple[str, str]] | None = None,
) -> dict:
    return {
        "connectivity": connectivity or set(),
        "exact_structures": exact_structures or set(),
        "scaffolds": scaffolds or set(),
        "registry_name_formula": registry_name_formula or set(),
        "development_truth_manifest_sha256": "a" * 64,
        "identity_registry_sha256": "b" * 64,
    }


def _manifest(count: int, *, add_repeat: bool = False) -> list[dict]:
    rows = []
    for index, smiles in enumerate(SMILES[:count], start=1):
        profile = rdkit_identity_profile(smiles)
        assert profile["valid"]
        rows.append(
            {
                "record_id": f"record-{index}-primary",
                "standard_id": f"standard-{index}",
                "method_case_id": f"methodcase-{index}",
                "detection_method_id": "agilent-qtof-cid-pos-v1",
                "validation_track": "novel_identity",
                "analysis_role": "primary",
                "replicate_id": "injection-1",
                "compound_name": f"External alcohol {index}",
                "formula": profile["structure_formula"],
                "inchikey": profile["inchikey"],
                "smiles": smiles,
                "truth_status": "authentic_standard",
                "truth_provenance": "certified_reference_material",
                "standard_lot": f"lot-{index}",
                "chemical_class": "aliphatic alcohol",
                "instrument_domain": "Agilent-6546-QTOF-CID",
                "polarity": "POS",
                "adduct": "[M+H]+",
                "acquisition_batch": f"batch-{1 + index % 2}",
                "run_order": index,
                "randomization_schedule_sha256": "d" * 64,
                "source_task_id": f"task-{index}",
                "source_entity_id": f"entity-{index}",
                "flagship_artifact_sha256": f"{index:064x}",
                "source_semantic_sha256": f"{index + 100:064x}",
            }
        )
    if add_repeat:
        repeat = copy.deepcopy(rows[0])
        repeat.update(
            {
                "record_id": "record-1-repeat",
                "analysis_role": "technical_repeat",
                "replicate_id": "injection-2",
                "source_task_id": "task-1-repeat",
                "source_entity_id": "entity-1-repeat",
                "run_order": 100,
                "flagship_artifact_sha256": "f" * 64,
                "source_semantic_sha256": "e" * 64,
            }
        )
        rows.append(repeat)
    return rows


def _method_variant(source: dict, *, suffix: str = "orbitrap") -> dict:
    row = copy.deepcopy(source)
    row.update(
        {
            "record_id": f"{source['record_id']}-{suffix}",
            "method_case_id": f"{source['method_case_id']}-{suffix}",
            "detection_method_id": f"thermo-orbitrap-hcd-pos-{suffix}",
            "validation_track": "method_transfer",
            "analysis_role": "method_variant",
            "instrument_domain": "Thermo-Orbitrap-HCD",
            "acquisition_batch": f"batch-{suffix}",
            "run_order": 1,
            "source_task_id": f"{source['source_task_id']}-{suffix}",
            "source_entity_id": f"{source['source_entity_id']}-{suffix}",
            "flagship_artifact_sha256": "a" * 64,
            "source_semantic_sha256": "b" * 64,
        }
    )
    return row


def _add_schedule_contract(
    rows: list[dict],
    policy: dict,
    *,
    cohort_id: str,
) -> None:
    repeat_count = sum(row["analysis_role"] == "technical_repeat" for row in rows)
    for row in rows:
        row.update(
            {
                "cohort_id": cohort_id,
                "randomization_algorithm_id": policy["acquisition_controls"][
                    "randomization_algorithm_id"
                ],
                "randomization_seed": 20260808,
                "technical_repeat_fraction_requested": 0.2 if repeat_count else 0.0,
            }
        )
    schedule_hash = external.compute_randomization_schedule_sha256(
        rows,
        cohort_id=cohort_id,
        policy=policy,
    )
    for row in rows:
        row["randomization_schedule_sha256"] = schedule_hash


def _freeze(monkeypatch: pytest.MonkeyPatch, rows: list[dict], policy: dict):
    _add_schedule_contract(rows, policy, cohort_id="external-unit-test")
    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: _boundaries(),
    )
    return external.freeze_external_cohort(
        rows,
        cohort_id="external-unit-test",
        policy=policy,
        sealing_salt=SALT,
    )


def _prediction_rows(
    public: dict,
    sealed: dict,
    policy: dict,
    outcomes: dict[str, tuple[str, str, str, str]],
) -> list[dict]:
    contract = policy["flagship_contract"]
    rows = []
    for case in public["cases"]:
        blind_code = case["blind_code"]
        record_id = sealed["cases"][blind_code]["record_id"]
        baseline, flagship, status, scope = outcomes[record_id]
        rows.append(
            {
                "blind_code": blind_code,
                "public_cohort_sha256": public["public_cohort_sha256"],
                "flagship_artifact_sha256": case["flagship_artifact_sha256"],
                "source_semantic_sha256": case["source_semantic_sha256"],
                "baseline_candidate_connectivity_key": baseline,
                "flagship_winner_connectivity_key": flagship,
                "flagship_decision_status": status,
                "flagship_decision_scope": scope,
                "flagship_algorithm_id": contract["algorithm_id"],
                "flagship_policy_sha256": contract["policy_semantic_sha256"],
                "identity_registry_sha256": contract["identity_registry_sha256"],
            }
        )
    return rows


def test_frozen_contract_matches_the_current_v11_policy_and_registry():
    policy = external.load_external_validation_policy()
    actual = external.validate_frozen_contract(policy)

    assert actual == {
        field: policy["flagship_contract"][field]
        for field in (
            "algorithm_id",
            "policy_semantic_sha256",
            "identity_registry_sha256",
        )
    }


def test_frozen_sample_size_sensitivity_table_is_exactly_reproducible():
    audit = external.audit_sample_size_plan(external.load_external_validation_policy())

    assert audit["passed"] is True
    assert all(check["passed"] for check in audit["checks"])
    primary = audit["checks"][1]
    assert primary["calculated_n_80"] == 185
    assert primary["calculated_n_90"] == 239


def test_freeze_hides_truth_and_does_not_count_technical_repeats(monkeypatch):
    policy = _policy(minimum_identities=3)
    public, sealed = _freeze(monkeypatch, _manifest(3, add_repeat=True), policy)

    assert public["cohort_summary"]["records"] == 4
    assert public["cohort_summary"]["independent_connectivity_identities"] == 3
    assert public["cohort_summary"]["novel_identity_primary_records"] == 3
    assert public["cohort_summary"]["technical_repeat_records"] == 1
    assert public["cohort_summary"]["confirmatory_minimum_reached"] is True
    rendered_public = json.dumps(public, ensure_ascii=False)
    for forbidden in (
        "compound_name",
        "standard_id",
        "canonical_isomeric_smiles",
        '"inchikey"',
    ):
        assert forbidden not in rendered_public
    assert sealed["sealed_truth_sha256"] == public["sealed_truth_sha256"]
    assert len(sealed["cases"]) == 4


def test_freeze_rejects_overlap_only_for_novel_identity_primary(monkeypatch):
    policy = _policy()
    rows = _manifest(1)
    profile = rdkit_identity_profile(rows[0]["smiles"])

    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: _boundaries(
            connectivity={profile["connectivity_key"]}
        ),
    )
    _add_schedule_contract(rows, policy, cohort_id="overlap")
    with pytest.raises(ValueError, match="novel-identity primary overlaps"):
        external.freeze_external_cohort(
            rows,
            cohort_id="overlap",
            policy=policy,
            sealing_salt=SALT,
        )

    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: _boundaries(
            registry_name_formula={("external alcohol 1", rows[0]["formula"])}
        ),
    )
    _add_schedule_contract(rows, policy, cohort_id="registry-overlap")
    with pytest.raises(ValueError, match="novel-identity primary is covered"):
        external.freeze_external_cohort(
            rows,
            cohort_id="registry-overlap",
            policy=policy,
            sealing_salt=SALT,
        )


def test_known_identity_and_same_compound_cross_method_records_are_retained(
    monkeypatch,
):
    policy = _policy(minimum_identities=1)
    primary = _manifest(1)[0]
    variant = _method_variant(primary)
    profile = rdkit_identity_profile(primary["smiles"])

    known_only = [copy.deepcopy(variant)]
    _add_schedule_contract(known_only, policy, cohort_id="known-method-transfer")
    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: _boundaries(
            connectivity={profile["connectivity_key"]},
            registry_name_formula={("external alcohol 1", primary["formula"])},
        ),
    )
    public_known, _ = external.freeze_external_cohort(
        known_only,
        cohort_id="known-method-transfer",
        policy=policy,
        sealing_salt=SALT,
    )
    assert public_known["cohort_summary"]["novel_identity_primary_records"] == 0
    assert public_known["cohort_summary"]["method_variant_records"] == 1
    assert (
        public_known["externality_audit"]["method_variant_development_overlap_records"]
        == 1
    )
    assert (
        public_known["externality_audit"][
            "method_variant_identity_registry_overlap_records"
        ]
        == 1
    )

    rows = [primary, variant]
    _add_schedule_contract(rows, policy, cohort_id="same-compound-two-methods")
    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: _boundaries(),
    )
    public, sealed = external.freeze_external_cohort(
        rows,
        cohort_id="same-compound-two-methods",
        policy=policy,
        sealing_salt=SALT,
    )
    truths = {
        truth["record_id"]: truth["truth_connectivity_key"]
        for truth in sealed["cases"].values()
    }
    wrong = rdkit_identity_profile(SMILES[1])["connectivity_key"]
    outcomes = {
        record_id: (
            wrong,
            truth,
            "selected",
            "2d_connectivity_candidate",
        )
        for record_id, truth in truths.items()
    }
    locked = external.lock_predictions(
        public,
        _prediction_rows(public, sealed, policy, outcomes),
        policy=policy,
        derivation="verified_projection_from_pinned_v11_flagship_artifacts",
    )
    report = external.evaluate_external_validation(
        public=public,
        sealed_truth=sealed,
        prediction_bundle=locked,
        policy=policy,
    )
    assert report["primary"]["independent_identities"] == 1
    assert (
        report["primary"]["method_variant_records_excluded_from_novel_identity_n"] == 1
    )
    assert report["secondary"]["method_transfer"]["records"] == 1
    assert report["secondary"]["method_transfer"]["identity_clusters"] == 1
    assert (
        report["secondary"]["same_identity_cross_method_agreement"]["identities"] == 1
    )


def test_direct_freeze_rejects_duplicate_identity_detection_method_cell(
    monkeypatch,
):
    policy = _policy(minimum_identities=1)
    primary = _manifest(1)[0]
    duplicate = _method_variant(primary, suffix="duplicate-method")
    duplicate["detection_method_id"] = primary["detection_method_id"]
    rows = [primary, duplicate]
    _add_schedule_contract(rows, policy, cohort_id="duplicate-method-cell")
    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: _boundaries(),
    )

    with pytest.raises(ValueError, match="duplicate identity-detection-method cells"):
        external.freeze_external_cohort(
            rows,
            cohort_id="duplicate-method-cell",
            policy=policy,
            sealing_salt=SALT,
        )


def test_freeze_rejects_duplicate_run_order_inside_an_acquisition_block(monkeypatch):
    policy = _policy()
    rows = _manifest(2)
    rows[1]["acquisition_batch"] = rows[0]["acquisition_batch"]
    rows[1]["run_order"] = rows[0]["run_order"]
    monkeypatch.setattr(
        external,
        "_development_boundaries",
        lambda selected_policy, *, project_root: _boundaries(),
    )
    _add_schedule_contract(rows, policy, cohort_id="duplicate-run-order")

    with pytest.raises(ValueError, match="run_order must be unique"):
        external.freeze_external_cohort(
            rows,
            cohort_id="duplicate-run-order",
            policy=policy,
            sealing_salt=SALT,
        )


def test_prediction_lock_rejects_truth_fields_and_hash_mismatch(monkeypatch):
    policy = _policy(minimum_identities=3)
    public, sealed = _freeze(monkeypatch, _manifest(3), policy)
    truths = {
        truth["record_id"]: truth["truth_connectivity_key"]
        for truth in sealed["cases"].values()
    }
    outcomes = {
        record_id: (identity, identity, "selected", "2d_connectivity_candidate")
        for record_id, identity in truths.items()
    }
    rows = _prediction_rows(public, sealed, policy, outcomes)
    poisoned = copy.deepcopy(rows)
    poisoned[0]["truth_connectivity_key"] = truths["record-1-primary"]
    with pytest.raises(ValueError, match="prediction schema mismatch"):
        external.lock_predictions(public, poisoned, policy=policy)

    wrong_hash = copy.deepcopy(rows)
    wrong_hash[0]["source_semantic_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="source_semantic_sha256 differs"):
        external.lock_predictions(public, wrong_hash, policy=policy)

    locked = external.lock_predictions(public, rows, policy=policy)
    assert len(locked["predictions"]) == 3
    assert len(locked["prediction_bundle_sha256"]) == 64


def test_artifact_projection_derives_the_frozen_baseline_and_flagship(
    monkeypatch, tmp_path
):
    policy = _policy(minimum_identities=1)
    manifest = _manifest(1)
    truth = rdkit_identity_profile(manifest[0]["smiles"])["connectivity_key"]
    wrong = rdkit_identity_profile(SMILES[1])["connectivity_key"]
    source_semantic = "c" * 64
    contract = policy["flagship_contract"]
    artifact = {
        "metadata": {
            "task_id": manifest[0]["source_task_id"],
            "algorithm_id": contract["algorithm_id"],
            "policy_sha256": contract["policy_semantic_sha256"],
            "identity_registry_sha256": contract["identity_registry_sha256"],
            "source_semantic_sha256": source_semantic,
        },
        "entities": [
            {
                "entity_id": manifest[0]["source_entity_id"],
                "scan_top1_audit_baseline": {"identity_key": f"std2d:{wrong}"},
                "winner": {"identity_key": f"std2d:{truth}"},
                "decision": {
                    "status": "selected",
                    "scope": "2d_connectivity_candidate",
                },
            }
        ],
    }
    artifact_path = tmp_path / "task-1_v11_flagship.json"
    artifact_path.write_text(
        json.dumps(artifact, ensure_ascii=False, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    manifest[0]["flagship_artifact_sha256"] = external.file_sha256(artifact_path)
    manifest[0]["source_semantic_sha256"] = source_semantic
    public, sealed = _freeze(monkeypatch, manifest, policy)
    locked = external.lock_predictions_from_flagship_artifacts(
        public,
        [
            {
                "blind_code": public["cases"][0]["blind_code"],
                "flagship_artifact_path": str(artifact_path),
            }
        ],
        policy=policy,
    )

    prediction = locked["predictions"][0]
    assert locked["prediction_derivation"] == (
        "verified_projection_from_pinned_v11_flagship_artifacts"
    )
    assert prediction["baseline_candidate_connectivity_key"] == wrong
    assert prediction["flagship_winner_connectivity_key"] == truth
    report = external.evaluate_external_validation(
        public=public,
        sealed_truth=sealed,
        prediction_bundle=locked,
        policy=policy,
    )
    assert (
        report["gate_checks"][
            "paired_predictions_derived_from_pinned_flagship_artifacts"
        ]
        is True
    )
    assert report["primary"]["paired_outcomes"]["gained"] == 1


def test_strict_endpoint_counts_ambiguous_winner_as_incorrect_and_excludes_repeat(
    monkeypatch,
):
    policy = _policy(minimum_identities=4)
    public, sealed = _freeze(monkeypatch, _manifest(4, add_repeat=True), policy)
    truths = {
        truth["record_id"]: truth["truth_connectivity_key"]
        for truth in sealed["cases"].values()
    }
    wrong = truths["record-4-primary"]
    outcomes = {
        "record-1-primary": (
            truths["record-1-primary"],
            truths["record-1-primary"],
            "selected",
            "2d_connectivity_candidate",
        ),
        "record-1-repeat": (
            truths["record-1-repeat"],
            truths["record-1-repeat"],
            "selected",
            "2d_connectivity_candidate",
        ),
        "record-2-primary": (
            truths["record-2-primary"],
            truths["record-2-primary"],
            "selected_candidate_only",
            "ambiguous_connectivity_candidate",
        ),
        "record-3-primary": (
            wrong,
            truths["record-3-primary"],
            "selected",
            "2d_connectivity_candidate",
        ),
        "record-4-primary": (
            "",
            truths["record-4-primary"],
            "selected",
            "2d_connectivity_candidate",
        ),
    }
    locked = external.lock_predictions(
        public,
        _prediction_rows(public, sealed, policy, outcomes),
        policy=policy,
    )
    report = external.evaluate_external_validation(
        public=public,
        sealed_truth=sealed,
        prediction_bundle=locked,
        policy=policy,
    )

    assert report["primary"]["independent_identities"] == 4
    assert report["primary"]["technical_repeat_records_excluded"] == 1
    assert report["primary"]["baseline_scan_top1"]["correct"] == 2
    assert report["primary"]["v11_raw_winner_secondary"]["correct"] == 4
    assert report["primary"]["v11_strict_connectivity_release"]["correct"] == 3
    assert report["primary"]["v11_strict_connectivity_release"]["released"] == 3
    assert report["primary"]["paired_outcomes"] == {
        "gained": 2,
        "lost": 1,
        "both_correct": 1,
        "both_incorrect": 0,
    }
    assert (
        report["secondary"]["repeatability"]["technical_repeats_added_to_primary_n"]
        == 0
    )
    assert report["confirmatory_accuracy_claim_allowed"] is False


def test_all_confirmatory_gates_can_pass_only_after_lock_and_minimum_n(monkeypatch):
    policy = _policy(minimum_identities=8)
    public, sealed = _freeze(monkeypatch, _manifest(8), policy)
    truths = {
        truth["record_id"]: truth["truth_connectivity_key"]
        for truth in sealed["cases"].values()
    }
    wrong = truths["record-8-primary"]
    outcomes = {}
    for index in range(1, 9):
        record_id = f"record-{index}-primary"
        truth = truths[record_id]
        if index <= 6:
            outcomes[record_id] = (
                wrong,
                truth,
                "selected",
                "2d_connectivity_candidate",
            )
        else:
            outcomes[record_id] = (
                wrong,
                wrong,
                "selected",
                "2d_connectivity_candidate",
            )
    locked = external.lock_predictions(
        public,
        _prediction_rows(public, sealed, policy, outcomes),
        policy=policy,
        derivation="verified_projection_from_pinned_v11_flagship_artifacts",
    )
    report = external.evaluate_external_validation(
        public=public,
        sealed_truth=sealed,
        prediction_bundle=locked,
        policy=policy,
    )

    assert report["primary"]["paired_outcomes"]["gained"] == 6
    assert report["primary"]["paired_outcomes"]["lost"] == 0
    assert report["primary"]["two_sided_exact_mcnemar_p_value"] == pytest.approx(
        0.03125
    )
    assert report["primary"]["accuracy_delta"]["ci_lower"] > 0
    assert all(report["gate_checks"].values())
    assert report["confirmatory_accuracy_claim_allowed"] is True
    assert report["automatic_deployment_authorized"] is False
