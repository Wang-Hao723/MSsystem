from __future__ import annotations

import copy

import pytest

from eval.v16_candidate_ranking.complete_panel import rescore_complete_panel
from matching.evirank_retrieval import load_evirecall_policy

A, B, C = "A" * 14, "B" * 14, "C" * 14


def query(spectrum_id="q1"):
    return {"spectrum_id": spectrum_id, "precursor_mz": 200., "ion_mode": "NEG",
        "raw_peaks": [[80., 100.], [100., 50.], [120., 60.], [160., 20.]],
        "raw_peaks_complete": True, "collision_energy": 10., "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed", "activation_method": "CID", "instrument_type": "QTOF",
        "instrument": "Agilent 6546", "query_adduct": "[M-H]-", "precursor_charge": -1}


def reference(identity, peaks, spectrum_id="ref"):
    result = query(spectrum_id)
    result.update({"connectivity_identity": identity, "peaks": peaks, "charge": -1})
    result.pop("raw_peaks")
    return result


class Retriever:
    library_snapshot_id = "synthetic-frozen-library"
    def __init__(self, records):
        self.records = records
        self.policy = load_evirecall_policy()
    def _fetch_window(self, precursor, polarity):
        assert precursor == 200.
        assert polarity == "NEG"
        return copy.deepcopy(self.records), {"status": "available"}


def refs():
    return [reference(A, [[80., 100.], [100., 55.], [120., 58.], [160., 21.]], "a"),
            reference(B, [[83., 100.], [104., 50.], [124., 60.], [165., 20.]], "b")]


def test_every_candidate_is_scored_on_every_query_and_view_without_topk_filtering():
    rows = [query("q1"), query("q2")]
    rows[1]["match"] = {"top_candidates": [{"inchikey": A}]}
    original = copy.deepcopy(rows)
    panel = rescore_complete_panel(rows, [A, B], Retriever(refs()))
    assert len(panel["cells"]) == 2 * 2 * 2
    assert panel["cell_status_counts"] == {"measured": 8}
    assert panel["measured_zero_modified_cosine_cells"] == 4
    assert panel["condition_matched_cells"] == 8
    assert panel["new_winner_selected"] is False
    assert rows == original


def test_missing_reference_is_null_not_measured_zero():
    panel = rescore_complete_panel([query()], [A, B, C], Retriever(refs()))
    missing = [cell for cell in panel["cells"] if cell["candidate_connectivity"] == C]
    assert len(missing) == 2
    assert all(cell["all_condition_reference_metrics"] is None for cell in missing)
    assert all(cell["status"] != "measured" for cell in missing)


def test_unknown_conditions_do_not_become_condition_matched_support():
    row = query()
    row.pop("activation_method")
    panel = rescore_complete_panel([row], [A], Retriever(refs()))
    assert panel["cell_status_counts"] == {"measured": 2}
    assert panel["condition_matched_cells"] == 0
    assert all(cell["condition_matched_reference_metrics"] is None for cell in panel["cells"])


def test_reference_mirrors_are_excluded_and_duplicates_do_not_create_votes():
    duplicate = copy.deepcopy(refs()[0])
    duplicate["spectrum_id"] = "duplicate"
    mirrored = reference(C, query()["raw_peaks"], "query-copy")
    panel = rescore_complete_panel([query()], [A, C], Retriever([*refs(), duplicate, mirrored]))
    assert all(cell["reference_profile_count"] == 1 for cell in panel["cells"] if cell["candidate_connectivity"] == A)
    assert all(cell["all_condition_reference_metrics"] is None for cell in panel["cells"] if cell["candidate_connectivity"] == C)


def test_candidate_order_names_and_truth_fields_do_not_enter_scores():
    row = query()
    expected = rescore_complete_panel([row], [A, B], Retriever(refs()))
    row.update(compound_name="label", ground_truth=A, expected_compound=B, expert_score=1.)
    assert rescore_complete_panel([row], [B, A], Retriever(refs())) == expected


def test_full_40_candidate_panel_retains_all_cells_even_when_references_are_absent():
    import string
    ids = ["Z" * 12 + a + b for a in string.ascii_uppercase for b in string.ascii_uppercase][:40]
    panel = rescore_complete_panel([query()], ids, Retriever([]))
    assert panel["expected_cells"] == 80
    assert len(panel["cells"]) == 80
    assert all(cell["all_condition_reference_metrics"] is None for cell in panel["cells"])


def test_duplicate_query_ids_are_rejected():
    with pytest.raises(ValueError, match="Unique explicit"):
        rescore_complete_panel([query(), query()], [A], Retriever(refs()))
