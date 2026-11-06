from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from eval.v15_top40_restricted_rerank.minimal_consensus_shadow import (
    ALGORITHM_ID,
    build_shadow_decision,
    load_policy,
)
from eval.v15_top40_restricted_rerank.freeze_minimal_consensus_predictions import (
    freeze_predictions,
    validate_query_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = (
    ROOT
    / "docs/audits/v15_top40_restricted_rerank_phase1/pairwise_ranker_model.json"
)


def _model() -> dict:
    import json

    return json.loads(MODEL_PATH.read_text(encoding="utf-8"))


def _row(
    connectivity: str,
    formula: str,
    rank: int | None,
    values: list[float | None],
    *,
    name: str,
) -> dict:
    fields = load_policy()["score_fields"]
    return {
        "query_id": "query-1",
        "candidate": {
            "candidate_connectivity": f"IK14:{connectivity}",
            "candidate_formula": formula,
            "feature_vector": dict(zip(fields, values)),
            "audit_only": {
                "baseline_rank": rank,
                "candidate_origin": "native" if rank is not None else "expanded",
            },
            "compound_name": name,
            "compound_name_cn": f"translated-{name}",
            "candidate_origin": "ignored",
        },
    }


def test_policy_is_frozen_shadow_and_excludes_automs() -> None:
    policy = load_policy()
    assert policy["algorithm_id"] == ALGORITHM_ID
    assert policy["execution"]["runtime_effect"] is False
    assert policy["execution"]["production_ranking_changed"] is False
    assert policy["execution"]["production_compound_count_changed"] is False
    assert policy["validation"]["automs_samples_allowed"] is False
    assert policy["selection_timing"]["existing_transfer_results_are_confirmatory"] is False


def test_three_score_field_wins_propose_only_an_unresolved_shadow_order() -> None:
    baseline = _row(
        "AAAAAAAAAAAAAA",
        "C15H10O5",
        1,
        [0.5, 0.5, 0.5, 0.9, 0.9, 0.9, 0.9, 0.9],
        name="Database baseline",
    )
    challenger = _row(
        "BBBBBBBBBBBBBB",
        "C15H10O5",
        2,
        [0.8, 0.8, 0.8, 0.1, 0.1, 0.1, 0.1, 0.1],
        name="Database challenger",
    )
    decision = build_shadow_decision(
        [baseline, challenger], _model(), polarity="NEG"
    )
    assert decision["gate_applied"] is True
    assert decision["score_field_vote"]["strictly_better_count"] == 3
    assert decision["production_top1_connectivity"] == "AAAAAAAAAAAAAA"
    assert decision["shadow_display_lead_connectivity"] == "BBBBBBBBBBBBBB"
    assert decision["production_top5_connectivities"][:2] == [
        "AAAAAAAAAAAAAA",
        "BBBBBBBBBBBBBB",
    ]
    assert decision["shadow_top5_connectivities"][:2] == [
        "BBBBBBBBBBBBBB",
        "AAAAAAAAAAAAAA",
    ]
    assert decision["reportable_scope"] == "molecular_formula_isomer_group"
    assert decision["reportable_single_connectivity"] is None
    assert decision["compound_count_increment"] == 0
    assert decision["production_ranking_changed"] is False
    assert decision["future_runtime_switch_eligible"] is False


def test_two_score_field_wins_retain_the_current_v15_winner() -> None:
    baseline = _row(
        "AAAAAAAAAAAAAA",
        "C15H10O5",
        1,
        [0.5, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9, 0.9],
        name="Baseline",
    )
    challenger = _row(
        "BBBBBBBBBBBBBB",
        "C15H10O5",
        2,
        [0.8, 0.8, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1],
        name="Challenger",
    )
    decision = build_shadow_decision(
        [baseline, challenger], _model(), polarity="NEG"
    )
    assert decision["gate_applied"] is False
    assert "challenger_wins_fewer_than_three_score_fields" in decision[
        "gate_failure_reasons"
    ]
    assert decision["shadow_display_lead_connectivity"] == "AAAAAAAAAAAAAA"


def test_support_only_candidate_never_enters_production_top5() -> None:
    baseline = _row(
        "AAAAAAAAAAAAAA",
        "C15H10O5",
        1,
        [0.2] * 8,
        name="Native baseline",
    )
    support_only = _row(
        "BBBBBBBBBBBBBB",
        "C15H10O5",
        None,
        [0.9] * 8,
        name="Expanded support candidate",
    )
    decision = build_shadow_decision(
        [baseline, support_only], _model(), polarity="NEG"
    )
    assert decision["gate_applied"] is True
    assert decision["production_top5_connectivities"] == ["AAAAAAAAAAAAAA"]
    assert decision["shadow_top5_connectivities"] == [
        "BBBBBBBBBBBBBB",
        "AAAAAAAAAAAAAA",
    ]
    assert decision["support_only_candidate_connectivities"] == [
        "BBBBBBBBBBBBBB"
    ]
    assert decision["production_ranking_changed"] is False


def test_positive_mode_unseen_formula_and_cross_formula_challengers_cannot_switch() -> None:
    baseline = _row(
        "AAAAAAAAAAAAAA",
        "C15H10O5",
        1,
        [0.1] * 8,
        name="Baseline",
    )
    cross_formula = _row(
        "BBBBBBBBBBBBBB",
        "C99H99O99",
        2,
        [1.0] * 8,
        name="Cross formula",
    )
    positive = build_shadow_decision(
        [baseline, cross_formula], _model(), polarity="POS"
    )
    assert positive["gate_applied"] is False
    assert positive["shadow_display_lead_connectivity"] == "AAAAAAAAAAAAAA"
    assert "polarity_outside_frozen_NEG_domain" in positive["gate_failure_reasons"]

    unseen_baseline = deepcopy(baseline)
    unseen_baseline["candidate"]["candidate_formula"] = "C99H99O99"
    unseen = build_shadow_decision(
        [unseen_baseline, cross_formula], _model(), polarity="NEG"
    )
    assert unseen["gate_applied"] is False
    assert "current_v15_winner_formula_outside_frozen_domain" in unseen[
        "gate_failure_reasons"
    ]


def test_candidate_names_and_translations_do_not_change_the_decision() -> None:
    baseline = _row(
        "AAAAAAAAAAAAAA",
        "C15H10O5",
        1,
        [0.4] * 8,
        name="Database name A",
    )
    challenger = _row(
        "BBBBBBBBBBBBBB",
        "C15H10O5",
        2,
        [0.8] * 8,
        name="Database name B",
    )
    original = build_shadow_decision(
        [baseline, challenger], _model(), polarity="NEG"
    )
    changed = deepcopy([baseline, challenger])
    changed[0]["candidate"]["compound_name"] = "完全不同的名称"
    changed[0]["candidate"]["compound_name_cn"] = "乱翻译一"
    changed[1]["candidate"]["compound_name"] = "another arbitrary name"
    changed[1]["candidate"]["compound_name_cn"] = "乱翻译二"
    renamed = build_shadow_decision(changed, _model(), polarity="NEG")
    comparable_keys = (
        "gate_applied",
        "production_top1_connectivity",
        "shadow_display_lead_connectivity",
        "shadow_candidate_order_connectivities",
        "score_field_vote",
        "compound_count_increment",
    )
    assert {key: original[key] for key in comparable_keys} == {
        key: renamed[key] for key in comparable_keys
    }
    assert renamed["candidate_name_used_for_scoring"] is False
    assert renamed["database_name_translation_allowed"] is False


def test_sealed_query_manifest_rejects_truth_identity_and_automs() -> None:
    clean = {
        "schema_version": "mssystem.v15.minimal_consensus_shadow.query_manifest/1.0.0",
        "sealed_truth_attached": False,
        "automs_sample_count": 0,
        "sample_type": "blind_authentic_standard",
        "cohort_id": "opaque-cohort-1",
        "queries": [
            {
                "query_id": "opaque-query-1",
                "record_id": "opaque-record-1",
                "polarity": "NEG",
                "task_id": "opaque-task-1",
                "detection_method_id": "method-1",
            }
        ],
    }
    assert validate_query_manifest(clean)[0]["polarity"] == "NEG"

    with_truth = deepcopy(clean)
    with_truth["queries"][0]["truth_connectivity"] = "AAAAAAAAAAAAAA"
    import pytest

    with pytest.raises(ValueError, match="Truth or identity field is forbidden"):
        validate_query_manifest(with_truth)

    automs = deepcopy(clean)
    automs["queries"][0]["task_id"] = "AutoMS-S3"
    with pytest.raises(ValueError, match="AutoMS"):
        validate_query_manifest(automs)


def test_prediction_freeze_retains_query_without_candidate_rows(tmp_path) -> None:
    import gzip
    import json

    matrix = tmp_path / "empty_matrix.jsonl.gz"
    with gzip.open(matrix, "wt", encoding="utf-8"):
        pass
    manifest = tmp_path / "query_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": (
                    "mssystem.v15.minimal_consensus_shadow.query_manifest/1.0.0"
                ),
                "sealed_truth_attached": False,
                "automs_sample_count": 0,
                "sample_type": "blind_authentic_standard",
                "cohort_id": "opaque-empty-cohort",
                "queries": [
                    {
                        "query_id": "opaque-query-without-candidate",
                        "record_id": "opaque-record-without-candidate",
                        "polarity": "NEG",
                        "task_id": "opaque-task",
                        "detection_method_id": "method",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "predictions.json"
    result = freeze_predictions(matrix, manifest, output)
    assert result["query_count"] == 1
    assert result["query_without_candidate_count"] == 1
    assert result["gate_applied_query_count"] == 0
    assert result["predictions"][0]["matrix_candidate_count"] == 0
    assert result["predictions"][0]["decision"]["status"] == (
        "unmatched_precursor_event"
    )
