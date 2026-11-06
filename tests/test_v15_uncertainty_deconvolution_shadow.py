import json

import pytest

from eval.v15_uncertainty_deconvolution.audit_selective_deconvolution import (
    _count_binary_tokens,
    _result_statuses,
)
from eval.v15_uncertainty_deconvolution.candidate_sets import (
    build_candidate_set,
    candidate_set_scope,
    conformal_higher_quantile,
    shadow_claim_action,
)


def test_split_conformal_quantile_uses_finite_sample_higher_order_statistic():
    result = conformal_higher_quantile([0.0, 0.1, 0.2, 0.3], alpha=0.25)
    assert result["requested_order_statistic"] == 4
    assert result["used_order_statistic"] == 4
    assert result["threshold"] == pytest.approx(0.3)
    assert result["order_statistic_capped"] is False


def test_candidate_set_ignores_names_and_retains_current_v15_winner():
    candidates = [
        {
            "connectivity": "std2d:AAAAAAAAAAAAAA",
            "formula": "C10H10O2",
            "score": 2.0,
            "compound_name": "must not be read",
            "translated_compound_name": "不得读取",
        },
        {
            "connectivity": "std2d:BBBBBBBBBBBBBB",
            "formula": "C10H10O2",
            "score": 1.0,
            "compound_name": "also ignored",
        },
    ]
    selected = build_candidate_set(
        candidates,
        score_gap_threshold=0.2,
        current_v15_winner="std2d:BBBBBBBBBBBBBB",
    )
    assert [row["connectivity"] for row in selected] == [
        "AAAAAAAAAAAAAA",
        "BBBBBBBBBBBBBB",
    ]
    assert all("compound_name" not in row for row in selected)
    assert candidate_set_scope(
        selected, current_v15_winner="BBBBBBBBBBBBBB"
    ) == "unresolved_isomer_candidate_group"


def test_different_formula_candidates_become_one_unresolved_event():
    selected = [
        {"connectivity": "AAAAAAAAAAAAAA", "formula": "C10H10O2", "score": 2.0},
        {"connectivity": "BBBBBBBBBBBBBB", "formula": "C9H8O3", "score": 1.9},
    ]
    scope = candidate_set_scope(
        selected, current_v15_winner="AAAAAAAAAAAAAA"
    )
    assert scope == "unresolved_precursor_event"
    action = shadow_claim_action("molecular_formula_isomer_group", scope)
    assert action["shadow_action"] == "shadow_demote_to_unresolved_precursor_event"
    assert action["unique_connectivity_release_created"] is False
    assert action["compound_count_increment"] == 0
    assert action["production_effect"] is False


def test_singleton_never_promotes_existing_candidate_scope():
    selected = [
        {"connectivity": "AAAAAAAAAAAAAA", "formula": "C10H10O2", "score": 2.0}
    ]
    scope = candidate_set_scope(
        selected, current_v15_winner="AAAAAAAAAAAAAA"
    )
    assert scope == "singleton_current_candidate_scope_retained_without_promotion"
    action = shadow_claim_action("single_collision_energy_candidate", scope)
    assert action["shadow_action"] == "retain_existing_v15_scope_without_promotion"
    assert action["unique_connectivity_release_created"] is False


def test_raw_metadata_counter_does_not_infer_missing_offsets(tmp_path):
    raw = tmp_path / "example.mzML"
    raw.write_bytes(
        b'<cvParam accession="MS:1000827"/><precursor spectrumRef="scan=1"/>'
        b'<cvParam accession="MS:1000744"/>'
    )
    counts = _count_binary_tokens(raw)
    assert counts["isolation_target_accession"] == 1
    assert counts["isolation_lower_offset_accession"] == 0
    assert counts["isolation_upper_offset_accession"] == 0
    assert counts["precursor_spectrum_reference"] == 1


def test_result_audit_separates_unknown_isolation_from_clean(tmp_path):
    result = tmp_path / "result.json"
    result.write_text(
        json.dumps(
            [
                {
                    "precursor_isolation_evidence": {
                        "status": "isolation_window_unavailable"
                    },
                    "match": {
                        "chimeric_deconvolution": {
                            "status": "isolation_window_unavailable",
                            "accepted": False,
                        }
                    },
                },
                {
                    "precursor_isolation_evidence": {"status": "mixed"},
                    "match": {
                        "chimeric_deconvolution": {
                            "status": "accepted_score_applied",
                            "accepted": True,
                        }
                    },
                },
            ]
        ),
        encoding="utf-8",
    )
    audit = _result_statuses(result)
    assert audit["result_count"] == 2
    assert audit["mixed_or_severe_isolation_count"] == 1
    assert audit["accepted_deconvolution_count"] == 1
    assert audit["precursor_isolation_status_counts"] == {
        "isolation_window_unavailable": 1,
        "mixed": 1,
    }
