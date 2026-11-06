from __future__ import annotations

import copy

import pytest

from matching.evirank_retrieval import load_evirecall_policy
from pipeline.v16_candidate_comparison import _align, _common_tracks
from pipeline.v16_complete_candidate_panel import rescore_complete_panel

A, B, C = "A" * 14, "B" * 14, "C" * 14


def query(energy=10., name=None):
    return {"spectrum_id": name or f"query-{energy}", "precursor_mz": 200.,
        "ion_mode": "NEG", "precursor_charge": -1, "query_adduct": "[M-H]-",
        "collision_energy": energy, "collision_energy_unit": "eV", "collision_energy_mode": "fixed",
        "instrument": "Agilent 6546", "instrument_type": "QTOF", "activation_method": "CID",
        "raw_peaks_complete": True,
        "raw_peaks": ([[70., 100.], [90., 80.], [110., 50.], [130., 30.], [150., 20.]] if energy == 10.
                      else [[70., 20.], [90., 90.], [110., 100.], [130., 40.], [150., 20.]])}


def reference(identity, energy=10., name=None):
    peaks = {
        (A, 10.): [[70., 100.], [90., 60.], [110., 60.], [130., 30.], [150., 20.]],
        (A, 20.): [[70., 22.], [90., 85.], [110., 95.], [130., 45.], [150., 22.]],
        (B, 10.): [[70., 100.], [90., 60.], [115., 60.], [135., 30.], [155., 20.]],
        (B, 20.): [[70., 95.], [90., 25.], [115., 60.], [135., 40.], [155., 25.]],
    }
    row = query(energy)
    row.pop("raw_peaks")
    row.update(spectrum_id=name or f"{identity}-{energy}", connectivity_identity=identity,
               peaks=peaks[identity, energy], charge=-1, adduct="[M-H]-", dataset_source="MassBank",
               molecular_formula="C9H10O4")
    return row


class Retriever:
    library_snapshot_id = "synthetic-test-not-validation"

    def __init__(self, records):
        self.records = copy.deepcopy(records)
        self.policy = load_evirecall_policy()

    def _fetch_window(self, precursor, polarity):
        assert precursor == 200. and polarity == "NEG"
        return copy.deepcopy(self.records), {"status": "available"}


def compare(rows=None, refs=None, pool=None):
    rows = [query(10.), query(20.)] if rows is None else rows
    refs = [reference(k, e) for k in (A, B) for e in (10., 20.)] if refs is None else refs
    return rescore_complete_panel(rows, pool or [A, B], Retriever(refs),
        include_comparisons=True, baseline_connectivity=A)["comparison_evidence"]


def test_complete_comparison_finds_observed_difference_ions_and_ce_trajectory():
    result = compare()
    pair = result["pairs"][0]
    assert pair["status"] == "compared"
    assert result["exact_condition_pair_count"] == result["trajectory_pair_count"] == 1
    assert pair["strict_evidence_preferred_candidate"] == A
    assert result["candidates"][0]["connectivity"] == A
    assert result["primary_selected_connectivity"] == A
    assert result["primary_ranking_changed"] is False
    assert all(t["energies"] == [10., 20.] for t in pair["ce_trajectories"])


def test_shared_ions_and_derived_neutral_losses_do_not_get_extra_votes():
    pair = compare()["pairs"][0]
    observation = pair["observations"][0]
    assert observation["shared_observed_product_ions_excluded"] == 2
    assert {x["query_mz"] for x in observation["a_observed_contrast_ions"]} == {110., 130., 150.}
    assert not observation["b_observed_contrast_ions"]
    assert observation["fragment_delta_a_minus_b"] < 1
    assert all(x["neutral_loss_da"] == 200 - x["query_mz"] for x in observation["a_observed_contrast_ions"])
    assert all(x["other_candidate_absence_proven"] is False for x in observation["a_observed_contrast_ions"])


def test_no_query_peak_absence_penalty_even_with_truncated_raw_declaration():
    rows = [query()]
    rows[0]["raw_peaks_complete"] = False
    result = compare(rows=rows)
    assert result["query_absence_penalties_applied"] is False
    assert result["trajectory_pair_count"] == 0
    assert all(x["query_mz"] in {p[0] for p in rows[0]["raw_peaks"]}
               for x in result["pairs"][0]["observations"][0]["a_observed_contrast_ions"])


@pytest.mark.parametrize("field", ["activation_method", "collision_energy_mode", "collision_energy_unit", "instrument", "query_adduct", "precursor_charge"])
def test_missing_metadata_stays_exploratory_and_does_not_become_same_condition(field):
    rows = [query(10.), query(20.)]
    for row in rows:
        row.pop(field)
    result = compare(rows=rows)
    assert result["compared_pair_count"] == 1
    assert result["exact_condition_pair_count"] == result["trajectory_pair_count"] == 0
    assert result["pairs"][0]["strict_evidence_preferred_candidate"] is None


@pytest.mark.parametrize("field,value", [("collision_energy_unit", "NCE"), ("collision_energy_mode", "stepped"),
    ("instrument", "Sciex TripleTOF 6600"), ("activation_method", "HCD"), ("adduct", "[M+Cl]-"),
    ("charge", -2), ("dataset_source", "")])
def test_incompatible_reference_strata_never_produce_ce_trajectory(field, value):
    refs = [reference(k, e) for k in (A, B) for e in (10., 20.)]
    for ref in refs:
        if ref["connectivity_identity"] == B:
            ref[field] = value
    result = compare(refs=refs)
    assert result["trajectory_pair_count"] == 0
    assert result["strict_preference_pair_count"] == 0


def test_near_but_unequal_reference_energy_is_not_exact_ce():
    refs = [reference(k, e) for k in (A, B) for e in (10., 20.)]
    for ref in refs:
        ref["collision_energy"] += 1
    assert compare(refs=refs)["exact_condition_pair_count"] == 0


def test_replicate_and_reference_mirror_duplicates_do_not_increase_scores():
    before = compare()
    refs = [reference(k, e) for k in (A, B) for e in (10., 20.)]
    duplicate = copy.deepcopy(refs[0])
    duplicate["spectrum_id"] = "z-duplicate"
    repeated_query = query(10., "z-repeat")
    after = compare(rows=[query(10.), query(20.), repeated_query], refs=refs + [duplicate] * 8)
    assert after["candidates"] == before["candidates"]
    assert after["pairs"][0]["exact_condition_paired_bins"] == before["pairs"][0]["exact_condition_paired_bins"]
    assert after["pairs"][0]["ce_trajectories"] == before["pairs"][0]["ce_trajectories"]


def test_names_truth_and_input_order_cannot_change_comparison():
    before = compare()
    rows = [query(20.), query(10.)]
    refs = [reference(k, e) for k in (B, A) for e in (20., 10.)]
    for row in rows + refs:
        row.update(compound_name="untrusted label", truth_key=B, target=A, expected_compound=B)
    assert compare(rows=rows, refs=refs, pool=[B, A]) == before


def test_absent_candidates_retain_null_scores_and_cannot_win():
    result = compare(pool=[A, B, C])
    candidate = next(x for x in result["candidates"] if x["connectivity"] == C)
    assert candidate["median_positive_fragment_contrast"] is None
    assert candidate["strict_pair_wins"] == 0
    assert result["pair_count"] == 3


def test_formula_conflicts_do_not_create_an_isomer_comparison():
    refs = [reference(A), reference(B)]
    refs[1]["molecular_formula"] = "C8H10O4"
    result = compare(rows=[query()], refs=refs)
    assert result["compared_pair_count"] == 0
    assert result["pairs"][0]["reason"] == "formula_missing_conflicting_or_different"


def test_single_energy_or_view_duplicates_do_not_create_trajectory():
    result = compare(rows=[query(10.), query(10., "repeat")])
    assert result["trajectory_pair_count"] == 0
    assert result["spectral_views_are_independent_replicates"] is False


def test_alignment_and_trajectory_intersection_are_one_to_one():
    alignment = _align([(100., 1.), (100.015, 1.)], [(100.007, 1.)], .02)
    assert len(alignment) == 1
    tracks = _common_tracks([
        {"tracks": [{"mz": 100., "query": .3, "a": .4, "b": .2},
                    {"mz": 100.015, "query": .6, "a": .5, "b": .7}]},
        {"tracks": [{"mz": 100.007, "query": .3, "a": .4, "b": .2}]}], .02)
    assert len(tracks) == 1


def test_missing_observation_at_one_energy_is_not_filled_with_zero_in_trajectory():
    rows = [query(10.), query(20.)]
    rows[1]["raw_peaks"] = [p for p in rows[1]["raw_peaks"] if p[0] != 90.]
    result = compare(rows=rows)
    assert result["trajectory_pair_count"] == 0  # only one shared observed track remains


def test_complete_40_candidate_pool_is_not_silently_truncated():
    ids = [chr(65 + i // 26) + chr(65 + i % 26) + "Z" * 12 for i in range(40)]
    result = rescore_complete_panel([query()], ids, Retriever([]), include_comparisons=True)["comparison_evidence"]
    assert result["candidate_count"] == 40
    assert result["pair_count"] == 780
    assert len(result["candidates"]) == 40
    assert result["compared_pair_count"] == 0


def test_baseline_outside_pool_is_rejected():
    with pytest.raises(ValueError, match="outside"):
        rescore_complete_panel([query()], [B], Retriever([reference(B)]),
            include_comparisons=True, baseline_connectivity=A)
