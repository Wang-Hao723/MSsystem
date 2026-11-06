from __future__ import annotations

import copy
import hashlib
import json

import pytest

from matching.evirank_mv import (
    FEATURE_ORDER,
    build_evirank_mv_result,
    build_method_case_index,
    load_evirank_policy,
)
from matching.evirank_retrieval import (
    evirecall_policy_sha256,
    load_evirecall_policy,
)
from matching.flagship_ranker import build_flagship_result, load_flagship_policy
from pipeline.evirank_mv import read_evirank_artifact, write_evirank_artifact


IK_A = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
IK_B = "DDDDDDDDDDDDDD-EEEEEEEEEE-F"


def _candidate(name: str, inchikey: str, strength: float, source: str) -> dict:
    return {
        "compound_name": name,
        "inchikey": inchikey,
        "formula": "C8H10N4O2",
        "smiles": "",
        "cosine": strength,
        "entropy": strength,
        "reverse_match": strength,
        "fragment_score": strength,
        "intensity_score": strength,
        "mass_ppm": 0.0,
        "dataset_source": source,
        "expert_score": 1.0,
        "classification": "test-class",
    }


def _v11_payload() -> dict:
    policy = copy.deepcopy(load_flagship_policy())
    policy["scan_score"]["reciprocal_rank_fraction"] = 0.01
    policy["scan_score"]["spectral_fraction"] = 0.99
    policy["entity_score"]["mean_collision_energy_support"] = 0.5499
    policy["entity_score"]["collision_energy_top1_fraction"] = 0.0001
    candidates = [
        _candidate("baseline-a", IK_A, 0.82, "MassBank"),
        _candidate("challenger-b", IK_B, 0.82, "GNPS"),
    ]
    rows = []
    for ce in (10.0, 20.0, 30.0):
        match = copy.deepcopy(candidates[0])
        match["top_candidates"] = copy.deepcopy(candidates)
        rows.append({
            "spectrum_id": f"scan-{ce:g}",
            "precursor_feature_id": "PF-1",
            "precursor_mz": 300.0,
            "retention_time": 1.0,
            "collision_energy": ce,
            "collision_energy_unit": "eV",
            "ion_mode": "POS",
            "tic": 1000.0 + ce,
            "raw_peaks": [[100.0, 100.0], [150.0, 50.0]],
            "query_peaks": [[100.0, 100.0], [150.0, 50.0]],
            "match": match,
        })
    return build_flagship_result("v12-test", rows, policy=policy)


def _model(*, coefficient_feature: str = "forward_spectrum_ensemble", coefficient: float = 4.0) -> dict:
    means = {feature: 0.0 for feature in FEATURE_ORDER}
    scales = {feature: 1.0 for feature in FEATURE_ORDER}
    coefficients = {feature: 0.0 for feature in FEATURE_ORDER}
    means[coefficient_feature] = 0.5
    scales[coefficient_feature] = 0.5
    coefficients[coefficient_feature] = coefficient
    return {
        "schema_version": "1.0.0",
        "model_id": "synthetic-qualified-shadow",
        "algorithm_id": "mssystem-v12-flagship-evirank-mv-0.2.0",
        "state": "qualified_shadow",
        "feature_order": list(FEATURE_ORDER),
        "feature_mean": means,
        "feature_scale": scales,
        "coefficients": coefficients,
        "intercept": 0.0,
        "maximum_absolute_delta": 0.08,
        "temperature": 1.0,
        "training_contract": {
            "connectivity_identity_disjoint": True,
            "same_identity_different_methods_retained": True,
            "technical_repeats_nonindependent": True,
            "ground_truth_unavailable_at_inference": True,
        },
    }


def _forward_bundle(payload: dict, *, domain_match: bool = True) -> dict:
    metadata = payload["metadata"]
    entity = payload["entities"][0]
    candidates = entity["candidates"]
    scores = {
        candidates[0]["identity_key"]: 0.0,
        candidates[1]["identity_key"]: 1.0,
    }
    return {
        "metadata": {"source_semantic_sha256": metadata["source_semantic_sha256"]},
        "entities": {
            entity["entity_id"]: {
                "candidates": {
                    identity: {
                        "forward_spectrum_ensemble": {
                            "score": score,
                            "eligible": True,
                            "candidate_conditioned": True,
                            "instrument_domain_match": domain_match,
                            "independent_architecture_families": 2,
                            "source_semantic_sha256": metadata["source_semantic_sha256"],
                        }
                    }
                    for identity, score in scores.items()
                }
            }
        },
    }


def _canonical_sha256(value: dict) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode("utf-8")
    ).hexdigest()


def _candidate_recall_bundle(
    payload: dict,
    method_id: str,
    *,
    qualified_switch: bool = False,
) -> dict:
    policy = load_evirecall_policy()
    consensus_policy = policy["source_paired_consensus"]
    thresholds = {
        "maximum_v11_top1_margin": consensus_policy["maximum_v11_top1_margin"],
        "minimum_same_formula_connectivities": consensus_policy[
            "minimum_same_formula_connectivities"
        ],
        "maximum_same_formula_connectivities": consensus_policy[
            "maximum_same_formula_connectivities"
        ],
        "maximum_expanded_same_formula_challengers": consensus_policy[
            "maximum_expanded_same_formula_challengers"
        ],
        "minimum_common_provenance_strata": consensus_policy[
            "minimum_common_provenance_strata"
        ],
        "minimum_unique_provenance_evidence_profiles": consensus_policy[
            "minimum_unique_provenance_evidence_profiles"
        ],
        "minimum_distinct_collision_energy_conditions": consensus_policy[
            "minimum_distinct_collision_energy_conditions"
        ],
        "minimum_spectral_views": consensus_policy["minimum_spectral_views"],
        "minimum_paired_observations": consensus_policy[
            "minimum_paired_observations"
        ],
        "minimum_median_score_delta": consensus_policy[
            "minimum_median_score_delta"
        ],
        "minimum_observation_win_fraction": consensus_policy[
            "minimum_observation_win_fraction"
        ],
        "minimum_provenance_stratum_win_fraction": consensus_policy[
            "minimum_provenance_stratum_win_fraction"
        ],
        "minimum_unique_provenance_profile_win_fraction": consensus_policy[
            "minimum_unique_provenance_profile_win_fraction"
        ],
    }
    expanded_identity = "std2d:GGGGGGGGGGGGGG"
    entities = []
    for entity in payload["entities"]:
        base_candidates = [
            candidate
            for candidate in entity["candidates"]
            if not candidate.get("vetoed")
        ]
        baseline_identity = base_candidates[0]["identity_key"]
        baseline_connectivity = baseline_identity.split(":", 1)[1]
        baseline_margin = round(
            float(base_candidates[0]["entity_score"])
            - float(base_candidates[1]["entity_score"]),
            8,
        )
        comparisons = []
        if qualified_switch:
            comparisons = [
                {
                    "challenger_identity_key": expanded_identity,
                    "challenger_connectivity_identity": "GGGGGGGGGGGGGG",
                    "paired_observation_count": 8,
                    "technical_repeat_observation_count": 8,
                    "common_provenance_stratum_count": 2,
                    "common_provenance_strata": ["source-a", "source-b"],
                    "unique_provenance_evidence_profile_count": 2,
                    "unique_provenance_evidence_profile_sizes": [1, 1],
                    "distinct_collision_energy_condition_count": 2,
                    "spectral_view_count": 2,
                    "spectral_views": ["library_match", "sensitive"],
                    "median_score_delta": 0.05,
                    "observation_win_count": 7,
                    "observation_win_fraction": 0.875,
                    "provenance_stratum_win_count": 2,
                    "provenance_stratum_win_fraction": 1.0,
                    "unique_provenance_profile_win_count": 2,
                    "unique_provenance_profile_win_fraction": 1.0,
                    "failed_gates": [],
                    "qualified": True,
                    "score_is_probability": False,
                }
            ]
        consensus = {
            "schema_version": "1.0.0",
            "module_id": consensus_policy["module_id"],
            "algorithm_id": consensus_policy["algorithm_id"],
            "deployment_state": "shadow",
            "status": "qualified_switch" if qualified_switch else "baseline_retained",
            "reason": (
                "qualified_source_paired_same_formula_challenger"
                if qualified_switch
                else "no_qualified_source_paired_challenger"
            ),
            "evaluated": qualified_switch,
            "baseline_identity_key": baseline_identity,
            "baseline_connectivity_identity": baseline_connectivity,
            "baseline_margin": baseline_margin,
            "formula": "C8H10N4O2",
            "same_formula_connectivity_count": 2,
            "thresholds": thresholds,
            "comparisons": comparisons,
            "selected_identity_key": (
                expanded_identity if qualified_switch else baseline_identity
            ),
            "selected_connectivity_identity": (
                "GGGGGGGGGGGGGG" if qualified_switch else baseline_connectivity
            ),
            "ranking_changed_from_v11": qualified_switch,
            "technical_repeats_collapsed_by_median": True,
            "reference_copies_collapsed_within_provenance_stratum": True,
            "mirrored_reference_content_collapsed_across_provenance_strata": True,
            "provenance_strata_are_independent_sources": False,
            "ground_truth_used_at_inference": False,
            "candidate_names_affect_scoring": False,
            "score_is_probability": False,
            "production_effect": False,
        }
        entities.append(
            {
                "schema_version": "2.0.0",
                "module_id": "v12_evirank_candidate_recall",
                "algorithm_id": policy["algorithm_id"],
                "algorithm_name": policy["display_name"],
                "deployment_state": "shadow",
                "entity_id": entity["entity_id"],
                "detection_method_id": method_id,
                "library_snapshot_id": "synthetic-snapshot",
                "policy_sha256": evirecall_policy_sha256(policy),
                "existing_connectivity_count": len(entity["candidates"]),
                "target_total_connectivities": 20,
                "expanded_candidate_count": 1,
                "expanded_candidates": [
                    {
                        "identity_key": expanded_identity,
                        "connectivity_identity": "GGGGGGGGGGGGGG",
                        "candidate": {"formula": "C8H10N4O2"},
                        "candidate_origin": "evirecall_expansion",
                        "enters_shadow_candidate_union": True,
                        "enters_residual_model_scoring": False,
                        "source_paired_selection_qualified": qualified_switch,
                        "score_is_probability": False,
                        "production_effect": False,
                    }
                ],
                "source_paired_consensus": consensus,
                "technical_repeats_increase_condition_count": False,
                "candidate_names_affect_retrieval_or_scoring": False,
                "expanded_candidates_enter_shadow_candidate_union": True,
                "expanded_candidates_enter_residual_model_scoring": False,
                "blind_window_widening": False,
                "blind_adduct_enumeration": False,
                "score_is_probability": False,
                "production_effect": False,
            }
        )
    return {
        "metadata": {
            "schema_version": "2.0.0",
            "module_id": "v12_evirank_candidate_recall",
            "algorithm_id": policy["algorithm_id"],
            "algorithm_name": policy["display_name"],
            "deployment_state": "shadow",
            "task_id": payload["metadata"]["task_id"],
            "detection_method_id": method_id,
            "base_algorithm_id": payload["metadata"]["algorithm_id"],
            "base_policy_sha256": payload["metadata"]["policy_sha256"],
            "base_source_semantic_sha256": payload["metadata"]["source_semantic_sha256"],
            "base_artifact_sha256": _canonical_sha256(payload),
            "library_snapshot_id": "synthetic-snapshot",
            "policy_sha256": evirecall_policy_sha256(policy),
            "source_row_count": payload["metadata"]["source_row_count"],
            "entity_count": len(entities),
            "entities_with_expansion": len(entities),
            "expanded_candidate_count": len(entities),
            "entities_with_qualified_source_paired_switch": (
                len(entities) if qualified_switch else 0
            ),
            "source_membership_partition_validated": True,
            "score_is_probability": False,
            "production_effect": False,
        },
        "policy": policy,
        "invariants": {
            "shadow_only": True,
            "ground_truth_used_at_inference": False,
            "candidate_names_affect_retrieval_or_scoring": False,
            "expanded_candidates_enter_shadow_candidate_union": True,
            "expanded_candidates_enter_residual_model_scoring": False,
            "source_paired_consensus_ground_truth_used_at_inference": False,
            "provenance_strata_are_independent_sources": False,
            "mirrored_reference_content_profiles_create_extra_votes": False,
            "same_connectivity_different_detection_methods_preexcluded": False,
            "existing_candidate_exclusion_scope": "within_one_entity_and_detection_method_only",
            "technical_repeats_increase_condition_count": False,
            "blind_window_widening": False,
            "blind_adduct_enumeration": False,
            "score_is_probability": False,
            "production_effect": False,
        },
        "entities": entities,
    }
def test_default_v12_is_named_shadow_and_preserves_v11_order():
    v11 = _v11_payload()
    result = build_evirank_mv_result(v11, detection_method_id="qtof-neg-v1")
    entity = result["entities"][0]

    assert result["metadata"]["algorithm_name"] == "V12-Flagship / EviRank-MV"
    assert result["metadata"]["deployment_state"] == "shadow"
    assert result["metadata"]["residual_model_state"] == "unavailable_no_frozen_model"
    assert result["metadata"]["production_identity_selector"] is False
    assert result["metadata"]["production_ranking_changed"] is False
    assert entity["winner"]["identity_key"] == v11["entities"][0]["winner"]["identity_key"]
    assert entity["ranking_changed_from_v11"] is False
    assert all(candidate["evirank_mv_residual"] == 0.0 for candidate in entity["candidates"])


def test_qualified_multiview_model_can_only_apply_a_bounded_shadow_residual():
    v11 = _v11_payload()
    result = build_evirank_mv_result(
        v11,
        view_bundle=_forward_bundle(v11),
        residual_model=_model(),
    )
    entity = result["entities"][0]

    assert entity["ranking_changed_from_v11"] is True
    assert entity["winner"]["identity_key"] == "ik14:DDDDDDDDDDDDDD"
    assert max(abs(candidate["evirank_mv_residual"]) for candidate in entity["candidates"]) <= 0.08
    assert all(candidate["production_effect"] is False for candidate in entity["candidates"])
    assert result["metadata"]["entities_reranked_from_v11"] == 1
    assert result["metadata"]["production_ranking_changed"] is False


def test_domain_mismatched_optional_view_is_neutral_not_negative():
    v11 = _v11_payload()
    baseline = build_evirank_mv_result(v11, residual_model=_model())
    mismatched = build_evirank_mv_result(
        v11,
        view_bundle=_forward_bundle(v11, domain_match=False),
        residual_model=_model(),
    )

    assert [c["identity_key"] for c in mismatched["entities"][0]["candidates"]] == [
        c["identity_key"] for c in baseline["entities"][0]["candidates"]
    ]
    assert all(c["evirank_mv_residual"] == 0.0 for c in mismatched["entities"][0]["candidates"])
    reasons = {
        audit["reason"]
        for candidate in mismatched["entities"][0]["candidates"]
        for audit in candidate["evirank_mv_view_audit"]
        if audit["view"] == "forward_spectrum_ensemble"
    }
    assert reasons == {"instrument_domain_mismatch"}


def test_names_and_expert_fields_cannot_enter_v12_feature_vector():
    original = _v11_payload()
    perturbed = copy.deepcopy(original)
    for candidate in perturbed["entities"][0]["candidates"]:
        candidate["candidate"]["compound_name"] = "ground-truth-looking-name"
        candidate["candidate"]["classification"] = "changed"
        candidate["candidate"]["expert_score"] = 999.0
    first = build_evirank_mv_result(original, residual_model=_model(coefficient_feature="ce_support_stability"))
    second = build_evirank_mv_result(perturbed, residual_model=_model(coefficient_feature="ce_support_stability"))

    assert [c["identity_key"] for c in first["entities"][0]["candidates"]] == [
        c["identity_key"] for c in second["entities"][0]["candidates"]
    ]
    assert [c["evirank_mv_score"] for c in first["entities"][0]["candidates"]] == [
        c["evirank_mv_score"] for c in second["entities"][0]["candidates"]
    ]


def test_residual_model_rejects_method_exclusion_contract():
    model = _model()
    model["training_contract"]["same_identity_different_methods_retained"] = False
    with pytest.raises(ValueError, match="training contract"):
        build_evirank_mv_result(_v11_payload(), residual_model=model)


def test_same_compound_different_methods_are_retained_and_repeats_are_nested():
    index = build_method_case_index([
        {
            "record_id": "qtof-primary",
            "connectivity_identity": "IK14-A",
            "detection_method_id": "agilent-qtof-neg-cid",
            "role": "primary",
        },
        {
            "record_id": "qtof-repeat",
            "connectivity_identity": "IK14-A",
            "detection_method_id": "agilent-qtof-neg-cid",
            "role": "technical_repeat",
        },
        {
            "record_id": "orbitrap-primary",
            "connectivity_identity": "IK14-A",
            "detection_method_id": "orbitrap-pos-hcd",
            "role": "primary",
        },
    ])

    assert index["records_retained"] == 3
    assert index["connectivity_identity_count"] == 1
    assert index["method_case_count"] == 2
    assert index["technical_repeat_count"] == 1
    assert {case["detection_method_id"] for case in index["method_cases"]} == {
        "agilent-qtof-neg-cid",
        "orbitrap-pos-hcd",
    }
    assert index["technical_repeats_increase_independent_identity_n"] is False
    assert index["method_cases_increase_independent_identity_n"] is False


def test_materialization_is_atomic_and_does_not_modify_v11(tmp_path):
    source = tmp_path / "v12-test_v11_flagship.json"
    source.write_text(json.dumps(_v11_payload(), ensure_ascii=False), encoding="utf-8")
    before = source.read_bytes()
    output = tmp_path / "v12-test_v12_evirank_mv.json"

    summary = write_evirank_artifact(
        "v12-test",
        source,
        output,
        detection_method_id="qtof-neg-v1",
    )
    page = read_evirank_artifact(output, limit=10)

    assert source.read_bytes() == before
    assert summary["v11_artifact_mutated"] is False
    assert summary["production_ranking_changed"] is False
    assert page["entities_total"] == 1
    assert page["metadata"]["algorithm_name"] == "V12-Flagship / EviRank-MV"
    assert not list(tmp_path.glob(".*.tmp-*"))


def test_policy_hard_codes_method_retention_and_non_probability_semantics():
    policy = load_evirank_policy()
    method = policy["method_case_contract"]

    assert method["same_identity_different_methods_are_retained"] is True
    assert method["technical_repeats_increase_independent_identity_n"] is False
    assert policy["residual_model"]["score_is_probability"] is False

    unsafe = copy.deepcopy(policy)
    unsafe["candidate_recall"]["expanded_candidates_enter_residual_model_scoring"] = True
    with pytest.raises(ValueError, match="residual-model scoring"):
        build_evirank_mv_result(_v11_payload(), policy=unsafe)


def test_unqualified_recall_candidate_enters_shadow_union_but_cannot_be_top1():
    v11 = _v11_payload()
    method_id = "qtof-pos-cid-v1"
    bundle = _candidate_recall_bundle(v11, method_id)

    result = build_evirank_mv_result(
        v11,
        detection_method_id=method_id,
        candidate_recall_bundle=bundle,
    )
    entity = result["entities"][0]

    assert result["metadata"]["candidate_recall_state"] == (
        "attached_shadow_union_source_paired_gated"
    )
    assert result["metadata"]["entities_with_candidate_recall_expansion"] == 1
    assert result["metadata"]["expanded_candidates_entered_shadow_candidate_union"] is True
    assert result["metadata"]["expanded_candidates_entered_residual_model_scoring"] is False
    assert entity["candidate_recall_expanded"] is True
    assert entity["candidate_recall_expanded_candidate_count"] == 1
    assert entity["candidate_recall_expansion"]["expanded_candidates"][0][
        "identity_key"
    ] == "std2d:GGGGGGGGGGGGGG"
    assert entity["expanded_candidates_entered_shadow_candidate_union"] is True
    assert entity["expanded_candidates_entered_residual_model_scoring"] is False
    assert entity["candidates"][-1]["identity_key"] == "std2d:GGGGGGGGGGGGGG"
    assert entity["candidates"][-1]["shadow_selection_eligible"] is False
    assert entity["candidates"][-1]["evirank_mv_score"] is None
    assert entity["winner"]["identity_key"] == v11["entities"][0]["winner"][
        "identity_key"
    ]


def test_qualified_source_paired_consensus_can_select_expanded_shadow_candidate():
    v11 = _v11_payload()
    method_id = "qtof-pos-cid-v1"
    bundle = _candidate_recall_bundle(v11, method_id, qualified_switch=True)

    result = build_evirank_mv_result(
        v11,
        detection_method_id=method_id,
        candidate_recall_bundle=bundle,
    )
    entity = result["entities"][0]

    assert entity["winner"]["identity_key"] == "std2d:GGGGGGGGGGGGGG"
    assert entity["winner"]["candidate_origin"] == "evirecall_expansion"
    assert entity["winner"]["evirank_mv_score"] is None
    assert entity["winner"]["source_paired_selection_qualified"] is True
    assert entity["ranking_changed_from_v11"] is True
    assert entity["source_paired_switch_applied"] is True
    assert entity["decision"]["status"] == "selected_candidate_only"
    assert entity["decision"]["scope"] == (
        "source_paired_same_formula_connectivity_candidate"
    )
    assert "same_formula_connectivity_requires_orthogonal_confirmation" in entity[
        "decision"
    ]["blockers"]
    assert result["metadata"]["entities_with_qualified_source_paired_switch"] == 1
    assert result["metadata"]["production_ranking_changed"] is False


def test_forged_source_paired_threshold_or_qualification_is_rejected():
    v11 = _v11_payload()
    method_id = "method-a"
    bundle = _candidate_recall_bundle(v11, method_id, qualified_switch=True)
    forged = copy.deepcopy(bundle)
    forged["entities"][0]["source_paired_consensus"]["thresholds"][
        "minimum_common_provenance_strata"
    ] = 1

    with pytest.raises(ValueError, match="thresholds differ"):
        build_evirank_mv_result(
            v11,
            detection_method_id=method_id,
            candidate_recall_bundle=forged,
        )

    weakened_policy = copy.deepcopy(bundle)
    weakened_policy["policy"]["source_paired_consensus"][
        "minimum_unique_provenance_evidence_profiles"
    ] = 1
    weakened_policy["entities"][0]["source_paired_consensus"]["thresholds"][
        "minimum_unique_provenance_evidence_profiles"
    ] = 1
    weakened_policy["metadata"]["policy_sha256"] = _canonical_sha256(
        weakened_policy["policy"]
    )

    with pytest.raises(ValueError, match="frozen policy hash"):
        build_evirank_mv_result(
            v11,
            detection_method_id=method_id,
            candidate_recall_bundle=weakened_policy,
        )


def test_candidate_recall_bundle_refuses_cross_method_rebinding():
    v11 = _v11_payload()
    bundle = _candidate_recall_bundle(v11, "method-a")

    with pytest.raises(ValueError, match="detection-method binding"):
        build_evirank_mv_result(
            v11,
            detection_method_id="method-b",
            candidate_recall_bundle=bundle,
        )
