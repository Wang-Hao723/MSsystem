from copy import deepcopy

import numpy as np
import pytest

from eval.v16_fragment_discrimination_20260907.run import (
    compare, content_unique, idf, kernel, rank_candidate_scores, score_event,
)


def fixture_event():
    def ref(key, peaks):
        return {"id": key, "fingerprint": key, "peaks": peaks}
    return {"query_id": "anonymous", "polarity": "POS", "current_top1": "A",
        "current_order": ["A", "B"], "current_formula": "C10H20O",
        "query_peaks": [[100., 10.], [150., 9.], [200., 8.]],
        "candidates": [
            {"connectivity": "A", "formula": "C10H20O", "references": [ref("a", [[100.02, 10.], [151., 9.], [202., 8.]])]},
            {"connectivity": "B", "formula": "C10H20O", "references": [ref("b", [[100., 10.], [150., 9.], [200., 8.]])]},
        ]}


def test_mass_kernel_tapers_continuously():
    result = kernel(np.asarray([100.]), np.asarray([100., 100.005, 100.015, 100.025, 100.1]))
    assert np.allclose(result, [[1., 1., .5, 0., 0.]], atol=1e-10)


def test_identical_spectra_score_one_and_disjoint_score_zero():
    q = [[100., 1.], [200., 3.], [300., 2.]]
    assert compare(q, q)["score"] == pytest.approx(1.)
    assert compare(q, [[110., 1.], [210., 3.], [310., 2.]])["score"] == 0.


def test_one_reference_peak_cannot_explain_two_query_peaks():
    result = compare([[100., 1.], [100.001, 1.]], [[100., 1.]])
    assert len(result["pairs"]) == 1
    assert sum(x > 0 for x in result["query_peak_support"]) == 1


def test_idf_counts_candidates_not_references():
    p = {"A": np.asarray([100., 100., 100., 150.]), "B": np.asarray([100.])}
    assert idf(np.asarray([100., 150.]), p)[1] > idf(np.asarray([100., 150.]), p)[0]
    assert idf(np.asarray([100.]), p)[0] == 1.


def test_missing_candidate_is_not_deleted_or_assigned_zero():
    e = fixture_event()
    e["candidates"][0]["references"] = []
    result = score_event(e)
    assert result["status"] == "incomplete_candidate_reference_panel"
    assert set(result["variants"].values()) == {"A"}
    assert result["scores"] == []


def test_spectrum_signal_can_change_rank_without_truth():
    e = fixture_event()
    before = deepcopy(e)
    result = score_event(e)
    assert result["variants"]["three_method_agreement"] == "B"
    assert e == before


def test_candidate_and_reference_order_and_duplicates_are_invariant():
    e = fixture_event()
    expected = score_event(e)
    e["candidates"].reverse()
    for c in e["candidates"]:
        c["references"] = list(reversed(c["references"] * 10))
    assert score_event(e) == expected


def test_mass_only_candidate_presence_does_not_change_rank():
    e = fixture_event()
    e["query_peaks"] = [[11., 1.], [17., 1.], [23., 1.]]
    result = score_event(e)
    assert result["variants"]["three_method_agreement"] == "A"
    assert result["variants"]["mass_kernel_cosine_all_native"] == "A"


def test_current_formula_is_not_replaced_by_agreement():
    e = fixture_event()
    e["candidates"][1]["formula"] = "C9H16O2"
    result = score_event(e)
    assert result["variants"]["mass_kernel_cosine_all_native"] == "B"
    assert result["variants"]["three_method_agreement"] == "A"


def test_truth_fields_rejected():
    e = fixture_event()
    e["truth_connectivity"] = "B"
    with pytest.raises(ValueError, match="Truth"):
        score_event(e)


def test_tied_candidates_cannot_be_selected_by_alphabetical_identity():
    records = [{"connectivity": "B", "score": 1.}, {"connectivity": "C", "score": 1.}]
    assert rank_candidate_scores(records, "score", "A") == ("A", 0.)
