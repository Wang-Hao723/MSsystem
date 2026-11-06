from __future__ import annotations

import copy
import json
import sqlite3

import pytest

from eval.v15_top40_restricted_rerank.build_candidate_matrix import (
    _query_reference_overlap_audit,
)
from eval.v15_top40_restricted_rerank.evaluate_phase1_f002_f008_generalization import (
    _cluster_effect,
    _load_policy as load_generalization_policy,
    _mark_same_injection_truth_collisions,
    _stream_selected_entities,
)
from eval.v15_top40_restricted_rerank.evaluate_applicability_gate_exploration import (
    _load_policy as load_applicability_policy,
    applicability_gate_decisions,
)

from eval.v15_top40_restricted_rerank.phase0_evidence import (
    assert_truth_free_matrix,
    build_candidate_evidence,
    collapse_technical_repeat_metrics,
    fuse_top40_candidates,
    load_phase0_policy,
    normalize_adduct,
    normalize_charge,
    project_authentic_query,
    query_intensity_coverage,
    query_metadata_gate,
)
from eval.v15_top40_restricted_rerank.phase1_pairwise_ranker import (
    FEATURE_FIELDS,
    candidate_connectivity,
    feature_vector,
    fit_model,
    load_policy as load_phase1_policy,
    rank_rows,
    validate_model,
)
from eval.v15_top40_restricted_rerank.phase1_governance import (
    governed_shadow_decision,
    load_governance_policy,
)


def _entry() -> dict:
    return {
        "reference_id": "Q-001:NEG",
        "source_task_id": "task-1",
        "source_standard_lot": "lab-a",
        "instrument_family": "agilent-qtof-6546",
        "polarity": "NEG",
        "compound_name": "must not enter matrix",
        "formula": "C10H12O2",
        "inchikey": "ABCDEFGHIJKLMN-AAAAAAAAAA-N",
        "authentic_standard": True,
        "spectra": [
            {
                "source_spectrum_id": "scan-10",
                "collision_energy_ev": 10,
                "precursor_mz": 163.0754,
                "adduct": "[M-H]-",
                "peaks": [[50.0, 10.0], [80.0, 50.0], [120.0, 100.0]],
            },
            {
                "source_spectrum_id": "scan-20",
                "collision_energy_ev": 20,
                "precursor_mz": 163.0754,
                "adduct": "[M-H]-",
                "peaks": [[50.0, 20.0], [90.0, 60.0], [120.0, 100.0]],
            },
        ],
    }


def _scan_result(query_id: str, condition: str, first: str, second: str) -> dict:
    return {
        "status": "available",
        "query_id": query_id,
        "condition_id": condition,
        "rankings": {
            "library_match": [first, second],
            "sensitive": [second, first],
        },
        "candidates": {
            first: {
                "connectivity_identity": first,
                "candidate": {"formula": "C10H12O2", "compound_name": "ignored"},
            },
            second: {
                "connectivity_identity": second,
                "candidate": {"formula": "C10H12O2", "compound_name": "ignored"},
            },
        },
    }


def test_policy_is_phase0_only_and_runtime_inert() -> None:
    policy = load_phase0_policy()
    assert policy["approval"] == {
        "phase_0_read_only_matrix": True,
        "phase_1_model_development": False,
        "phase_2_sealed_validation": False,
        "phase_3_runtime_integration": False,
    }
    assert policy["execution"]["runtime_effect"] is False
    assert policy["execution"]["reference_library_write_allowed"] is False
    assert policy["execution"]["primary_payload_write_allowed"] is False


def test_query_projection_is_invariant_to_name_truth_and_standard_fields() -> None:
    original = _entry()
    changed = copy.deepcopy(original)
    changed.update(
        {
            "compound_name": "shuffled display name",
            "formula": "C999H999",
            "inchikey": "ZZZZZZZZZZZZZZ-ZZZZZZZZZZ-Z",
            "authentic_standard": False,
            "truth_connectivity": "SHOULDNOTENTER",
            "truth_rank": 99,
        }
    )
    assert project_authentic_query(original) == project_authentic_query(changed)


def test_forbidden_identity_or_label_fields_are_rejected() -> None:
    policy = load_phase0_policy()
    with pytest.raises(ValueError, match="Forbidden identity/label/display"):
        assert_truth_free_matrix(
            {"candidate": {"truth_connectivity": "ABCDEFGHIJKLMN"}},
            policy,
        )
    with pytest.raises(ValueError, match="Forbidden identity/label/display"):
        assert_truth_free_matrix(
            {"candidate": {"compound_name": "not a feature"}},
            policy,
        )


def test_rrf_union_is_invariant_to_scan_result_input_order() -> None:
    policy = load_phase0_policy()
    a = "AAAAAAAAAAAAAA"
    b = "BBBBBBBBBBBBBB"
    c = "CCCCCCCCCCCCCC"
    native = [{"connectivity": a, "formula": "C10H12O2", "entity_score": 0.7}]
    scans = [
        _scan_result("scan-10", "10:ev", b, c),
        _scan_result("scan-20", "20:ev", c, b),
    ]
    forward = fuse_top40_candidates(native, scans, policy)
    reverse = fuse_top40_candidates(native, list(reversed(scans)), policy)
    assert forward == reverse
    assert [row["candidate_connectivity"] for row in forward] == [a, b, c]


def test_missing_activation_and_scan_range_fail_closed() -> None:
    policy = load_phase0_policy()
    gate = query_metadata_gate(project_authentic_query(_entry()), policy)
    assert gate["runtime_switch_eligible"] is False
    assert "activation_method" in gate["missing_required_fields"]
    assert "scan_range_mz" in gate["missing_required_fields"]
    assert gate["distinct_collision_energy_count"] == 2


def test_query_intensity_coverage_uses_direct_and_shifted_matches() -> None:
    coverage = query_intensity_coverage(
        [[50.0, 10.0], [80.0, 30.0], [100.0, 60.0]],
        [[50.0, 5.0], [90.0, 20.0]],
        tolerance_da=0.02,
        query_precursor_mz=110.0,
        reference_precursor_mz=100.0,
    )
    # 50 is a direct match and 100 matches shifted 90; 80 is unmatched.
    assert coverage == pytest.approx(0.7)


def test_adduct_charge_and_technical_repeats_are_normalized_without_extra_votes() -> None:
    assert normalize_adduct(" [M-H]1- ") == "[M-H]-"
    assert normalize_charge("1-") == -1
    assert normalize_charge(None, adduct="[M+2H]2+") == 2
    repeated = [
        {
            "condition_id": "20:ev",
            "spectral_view": "library_match",
            "reference_content_fingerprint": "fp-1",
            "reference_spectrum_id": "r1",
            "dataset_source": "source-a",
            "modified_cosine": 0.8,
            "spectral_entropy": 0.7,
            "reverse_match": 0.6,
            "matched_peaks": 5,
            "query_intensity_coverage": 0.5,
            "polarity_compatible": True,
            "charge_compatible": True,
            "adduct_compatible": True,
        },
        {
            "condition_id": "20:ev",
            "spectral_view": "library_match",
            "reference_content_fingerprint": "fp-1",
            "reference_spectrum_id": "r2",
            "dataset_source": "source-a",
            "modified_cosine": 0.6,
            "spectral_entropy": 0.5,
            "reverse_match": 0.4,
            "matched_peaks": 3,
            "query_intensity_coverage": 0.3,
            "polarity_compatible": True,
            "charge_compatible": True,
            "adduct_compatible": True,
        },
    ]
    collapsed = collapse_technical_repeat_metrics(repeated)
    assert len(collapsed) == 1
    assert collapsed[0]["technical_repeat_count"] == 2
    assert collapsed[0]["modified_cosine"] == pytest.approx(0.7)
    assert collapsed[0]["matched_peaks"] == pytest.approx(4.0)


def test_candidate_origin_is_audit_only_and_mirror_profiles_do_not_add_votes() -> None:
    policy = load_phase0_policy()
    identity = "ABCDEFGHIJKLMN"
    scan_results = [
        {
            "query_id": "scan-10",
            "condition_id": "10:ev",
            "rankings": {},
            "candidates": {
                identity: {
                    "source_reference_records": {
                        "library_match": {
                            "source-a": {
                                "spectrum_id": "r1",
                                "content_fingerprint": "fp-1",
                                "score": 0.8,
                            },
                            "source-b": {
                                "spectrum_id": "r2",
                                "content_fingerprint": "fp-1",
                                "score": 0.7,
                            },
                            "source-c": {
                                "spectrum_id": "r3",
                                "content_fingerprint": "fp-2",
                                "score": 0.6,
                            },
                        }
                    }
                }
            },
        }
    ]
    query_gate = {
        "runtime_switch_eligible": False,
        "failure_reasons": ["missing:activation_method"],
    }
    base_candidate = {
        "candidate_connectivity": identity,
        "candidate_formula": "C10H12O2",
        "formula_matches_current_winner": True,
        "rrf_utility": 0.05,
        "candidate_origin": "native",
        "baseline_rank": 1,
    }
    native = build_candidate_evidence(
        base_candidate,
        scan_results,
        [],
        {},
        query_gate,
        policy,
        query_source_group="query-lab",
    )
    expanded_candidate = {**base_candidate, "candidate_origin": "expanded", "baseline_rank": None}
    expanded = build_candidate_evidence(
        expanded_candidate,
        scan_results,
        [],
        {},
        query_gate,
        policy,
        query_source_group="query-lab",
    )
    assert native["feature_vector"] == expanded["feature_vector"]
    assert native["audit_only"]["candidate_origin"] == "native"
    assert expanded["audit_only"]["candidate_origin"] == "expanded"
    assert (
        native["source_stability_audit"][
            "distinct_reference_content_fingerprint_count"
        ]
        == 2
    )
    assert "candidate_origin" not in native["feature_vector"]


def test_query_reference_overlap_audit_detects_scaled_content_mirror(tmp_path) -> None:
    database = tmp_path / "library.db"
    connection = sqlite3.connect(database)
    connection.execute(
        "CREATE TABLE Spectra ("
        "spectrum_id TEXT, precursor_mz REAL, peaks_json TEXT, dataset_source TEXT)"
    )
    connection.execute(
        "INSERT INTO Spectra VALUES (?,?,?,?)",
        (
            "ref-1",
            163.0754,
            json.dumps([[50.0, 20.0], [80.0, 100.0], [120.0, 200.0]]),
            "library-a",
        ),
    )
    connection.commit()
    connection.close()
    audit = _query_reference_overlap_audit(
        database,
        [
            {
                "spectrum_id": "query-1",
                "precursor_mz": 163.0754,
                "peaks": [[50.0, 10.0], [80.0, 50.0], [120.0, 100.0]],
            }
        ],
        "query-lab",
        {"mass_retrieval": {"direct_precursor_da": 0.02, "direct_precursor_ppm": 10}},
    )
    assert audit["development_matrix_eligible"] is False
    assert audit["query_content_mirror_count"] == 1
    assert audit["query_content_mirror_pairs"] == [
        {"query_spectrum_id": "query-1", "reference_spectrum_id": "ref-1"}
    ]


def _phase1_row(identity: str, offset: float) -> dict:
    return {
        "query_id": "q1",
        "candidate": {
            "candidate_connectivity": f"std2d:{identity}",
            "feature_vector": {
                field: offset + index / 100.0
                for index, field in enumerate(FEATURE_FIELDS)
            },
            "audit_only": {"candidate_origin": "native", "baseline_rank": 1},
        },
    }


def test_phase1_policy_and_model_remain_shadow_only_with_stable_solver() -> None:
    policy = load_phase1_policy()
    assert policy["model"]["solver"] == "liblinear"
    assert policy["execution"]["runtime_effect"] is False
    assert policy["approval"]["phase_2_sealed_validation"] is False
    assert policy["approval"]["phase_3_runtime_integration"] is False


def test_phase1_rank_is_invariant_to_name_origin_and_input_order() -> None:
    truth = _phase1_row("AAAAAAAAAAAAAA", 0.9)
    other = _phase1_row("BBBBBBBBBBBBBB", 0.1)
    rows_by_query = {"q1": [truth, other]}
    labels = {"q1": {"truth_connectivity": "AAAAAAAAAAAAAA"}}
    model = fit_model(rows_by_query, labels, ["q1"])
    changed = copy.deepcopy(other)
    changed["candidate"].update(
        {
            "compound_name": "must not affect score",
            "candidate_origin": "expanded",
            "audit_only": {"candidate_origin": "expanded", "baseline_rank": 999},
        }
    )
    assert feature_vector(other) == feature_vector(changed)
    forward = [candidate_connectivity(row) for row, _ in rank_rows([truth, changed], model)]
    reverse = [candidate_connectivity(row) for row, _ in rank_rows([changed, truth], model)]
    assert forward == reverse == ["AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB"]


def test_phase1_model_is_deterministic_and_hash_tampering_is_rejected() -> None:
    truth = _phase1_row("AAAAAAAAAAAAAA", 0.9)
    other = _phase1_row("BBBBBBBBBBBBBB", 0.1)
    rows_by_query = {"q1": [truth, other]}
    labels = {"q1": {"truth_connectivity": "AAAAAAAAAAAAAA"}}
    first = fit_model(rows_by_query, labels, ["q1"])
    second = fit_model(rows_by_query, labels, ["q1"])
    assert first == second
    validate_model(first)
    tampered = copy.deepcopy(first)
    tampered["coefficient"][0] += 0.1
    with pytest.raises(ValueError, match="semantic hash mismatch"):
        validate_model(tampered)


def test_phase1_governance_never_releases_same_formula_switch_as_unique_compound() -> None:
    baseline = _phase1_row("AAAAAAAAAAAAAA", 0.1)
    challenger = _phase1_row("BBBBBBBBBBBBBB", 0.9)
    baseline["candidate"]["candidate_formula"] = "C10H12O2"
    challenger["candidate"]["candidate_formula"] = "C10H12O2"
    rows_by_query = {"q1": [challenger, baseline]}
    labels = {"q1": {"truth_connectivity": "BBBBBBBBBBBBBB"}}
    model = fit_model(rows_by_query, labels, ["q1"])
    decision = governed_shadow_decision([baseline, challenger], model)
    assert decision["status"] == "unresolved_same_formula_isomer_group"
    assert decision["reportable_scope"] == "molecular_formula_isomer_group"
    assert decision["candidate_group_connectivities"] == [
        "AAAAAAAAAAAAAA",
        "BBBBBBBBBBBBBB",
    ]
    assert decision["reportable_single_connectivity"] is None
    assert decision["compound_count_increment"] == 0


def test_phase1_governance_never_counts_cross_formula_switch_as_compound() -> None:
    baseline = _phase1_row("AAAAAAAAAAAAAA", 0.1)
    challenger = _phase1_row("BBBBBBBBBBBBBB", 0.9)
    baseline["candidate"]["candidate_formula"] = "C10H12O2"
    challenger["candidate"]["candidate_formula"] = "C9H10O3"
    model = fit_model(
        {"q1": [challenger, baseline]},
        {"q1": {"truth_connectivity": "BBBBBBBBBBBBBB"}},
        ["q1"],
    )
    decision = governed_shadow_decision([baseline, challenger], model)
    assert decision["status"] == "unresolved_cross_formula_event"
    assert decision["reportable_scope"] == "unresolved_precursor_event"
    assert decision["reportable_single_connectivity"] is None
    assert decision["compound_count_increment"] == 0


def test_phase1_governance_policy_forbids_names_translation_and_runtime_effect() -> None:
    policy = load_governance_policy()
    assert policy["execution"]["runtime_effect"] is False
    assert policy["identity_rules"]["database_name_used_for_scoring"] is False
    assert policy["identity_rules"]["database_name_translation_allowed"] is False


def test_generalization_policy_is_frozen_shadow_only_and_excludes_automs() -> None:
    policy = load_generalization_policy()
    assert policy["execution"]["runtime_effect"] is False
    assert policy["execution"]["model_refit_allowed"] is False
    assert policy["cohort"]["automs_samples_allowed"] is False
    assert (
        policy["analysis"]["primary_applicability_population"]
        == "mass_selected_noncollision_events_with_a_frozen_candidate_matrix"
    )
    assert policy["analysis"]["zero_candidate_events_in_primary_reranker_endpoint"] is False


def test_applicability_exploration_is_explicitly_post_hoc_and_runtime_inert() -> None:
    policy = load_applicability_policy()
    assert policy["execution"]["runtime_effect"] is False
    assert policy["execution"]["production_integration_allowed"] is False
    assert policy["execution"]["model_refit_allowed"] is False
    assert policy["cohort"]["automs_samples_allowed"] is False
    assert (
        policy["selection_timing"][
            "gate_hypotheses_selected_after_f002_f008_outcomes_were_examined"
        ]
        is True
    )
    assert policy["selection_timing"]["nominal_intervals_are_confirmatory"] is False


def test_applicability_gate_uses_only_current_formula_and_polarity() -> None:
    development = {"C10H12O2"}
    assert applicability_gate_decisions("C10H12O2", "NEG", development) == {
        "formula_seen": True,
        "formula_seen_negative": True,
    }
    assert applicability_gate_decisions("C10H12O2", "POS", development) == {
        "formula_seen": True,
        "formula_seen_negative": False,
    }
    assert applicability_gate_decisions("C9H10O3", "NEG", development) == {
        "formula_seen": False,
        "formula_seen_negative": False,
    }


def test_generalization_anchor_selection_ignores_candidate_availability(tmp_path) -> None:
    artifact = tmp_path / "flagship.json"
    artifact.write_text(
        json.dumps(
            {
                "entities": [
                    {
                        "entity_id": "lower-tic-with-candidate",
                        "representative_search_precursor_mz": 100.0001,
                        "representative_tic": 100.0,
                        "candidates": [{"connectivity": "SHOULD-NOT-MATTER"}],
                    },
                    {
                        "entity_id": "higher-tic-without-candidate",
                        "representative_search_precursor_mz": 100.0008,
                        "representative_tic": 200.0,
                        "candidates": [],
                    },
                ]
            }
        ),
        encoding="utf-8",
    )
    selected = _stream_selected_entities(
        artifact,
        [{"ion_hypotheses": [{"mz": 100.0, "label": "[M-H]-"}]}],
        tolerance_ppm=10.0,
    )
    assert selected[0]["entity"]["entity_id"] == "higher-tic-without-candidate"


def test_identity_cluster_effect_does_not_count_repeated_injections_as_votes() -> None:
    rows = [
        {
            "transfer_metric_eligible": True,
            "truth_connectivity": "A",
            "baseline_top1_connectivity_correct": False,
            "variant_top1_connectivity_correct": True,
        }
        for _ in range(3)
    ]
    rows.append(
        {
            "transfer_metric_eligible": True,
            "truth_connectivity": "B",
            "baseline_top1_connectivity_correct": True,
            "variant_top1_connectivity_correct": False,
        }
    )
    effect = _cluster_effect(rows, "truth_connectivity", replicates=1000, seed=7)
    assert effect["cluster_count"] == 2
    assert effect["equal_weight_percentage_point_gain"] == pytest.approx(0.0)
    assert effect["positive_cluster_count"] == 1
    assert effect["negative_cluster_count"] == 1


def test_truth_collision_is_detected_even_when_entity_has_no_candidates() -> None:
    labels = [
        {
            "run": "f005-neg1",
            "selected_entity_id": "shared-event",
            "query_spectrum_sha256": None,
        },
        {
            "run": "f005-neg1",
            "selected_entity_id": "shared-event",
            "query_spectrum_sha256": None,
        },
        {
            "run": "f005-neg1",
            "selected_entity_id": "unique-event",
            "query_spectrum_sha256": None,
        },
    ]
    _mark_same_injection_truth_collisions(labels)
    assert labels[0]["same_injection_duplicate_query_bundle"] is True
    assert labels[1]["same_injection_duplicate_query_bundle"] is True
    assert labels[2]["same_injection_duplicate_query_bundle"] is False
