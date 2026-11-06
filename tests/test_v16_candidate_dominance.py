from __future__ import annotations

import copy

import pytest

from eval.v16_candidate_ranking.dominance import FIELDS, compare, propose


def row(identity, score, rank=None, formula="C10H10O4"):
    return {"candidate": {"candidate_connectivity": identity, "candidate_formula": formula,
        "feature_vector": dict.fromkeys(FIELDS, score), "audit_only": {"baseline_rank": rank,
        "condition_view_contributions": [{"condition_id": "ce:10|ev|fixed", "spectral_view": "library_match"}]}}}


def test_unique_all_field_dominator_can_propose_same_formula_candidate():
    a, b = row("a", .5, 1), row("b", .6, 2)
    original = copy.deepcopy([a, b])
    result = propose([a, b])
    assert result["proposed_top1"] == "b"
    assert result["runtime_effect"] is False
    assert result["compound_count_increment"] == 0
    assert result["reportable_scope"] == "unresolved_same_formula_isomer_candidates"
    assert [a, b] == original


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), -1.0])
def test_missing_or_invalid_score_cannot_vote(value):
    a, b = row("a", .5, 1), row("b", .9, 2)
    b["candidate"]["feature_vector"][FIELDS[0]] = value
    assert propose([a, b])["changed"] is False


def test_real_zero_is_observed_but_explicit_missing_flag_is_not():
    a, b = row("a", 0, 1), row("b", .1, 2)
    assert compare(a, b)["dominates"] is True
    a["candidate"]["feature_missing"] = {FIELDS[0]: True}
    assert compare(a, b)["comparable"] is False


def test_one_bad_metric_is_not_compensated_by_six_good_metrics():
    a, b = row("a", .5, 1), row("b", .9, 2)
    b["candidate"]["feature_vector"][FIELDS[0]] = .4
    assert propose([a, b])["changed"] is False


def test_missing_query_condition_is_not_accepted_as_equal_comparison():
    a, b = row("a", .5, 1), row("b", .9, 2)
    b["candidate"]["audit_only"]["condition_view_contributions"] = []
    assert propose([a, b])["changed"] is False


def test_cross_formula_candidate_cannot_displace_baseline():
    assert propose([row("a", .5, 1), row("b", .9, 2, "C11H12O4")])["changed"] is False


def test_ties_and_incomparable_alternatives_preserve_baseline():
    assert propose([row("a", .5, 1), row("b", .5, 2)])["changed"] is False
    a, b, c = row("a", .4, 1), row("b", .6, 2), row("c", .6, 3)
    b["candidate"]["feature_vector"][FIELDS[0]] = .7
    c["candidate"]["feature_vector"][FIELDS[1]] = .7
    result = propose([a, b, c])
    assert result["proposed_top1"] == "a"
    assert result["reason"] == "multiple_incomparable_dominators"


def test_ranking_is_invariant_to_input_order_and_display_name():
    a, b = row("a", .5, 1), row("b", .6, 2)
    expected = propose([a, b])
    a["candidate"]["name"] = "known standard"
    b["candidate"]["name"] = "wrong label"
    assert propose([b, a]) == expected


def test_expansion_reorders_top5_without_deleting_other_candidates():
    rows = [row(str(index), .5, index + 1) for index in range(5)] + [row("new", .8)]
    result = propose(rows)
    assert result["proposed_top5"] == ["new", "0", "1", "2", "3"]
    assert len(rows) == 6


def test_duplicate_candidate_identity_is_rejected():
    with pytest.raises(ValueError, match="one matrix row"):
        propose([row("a", .5, 1), row("a", .6, 2)])
