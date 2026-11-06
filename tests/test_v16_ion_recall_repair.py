from copy import deepcopy
import json
from pathlib import Path

import pytest

from eval.v16_ion_recall_repair.core import IONS, adduct, resolve_query_ion, neutral_recall


def network_row():
    ion = IONS["[M+NH4]+"]
    return {"precursor_mz": 300 + ion.shift, "ion_mode": "POS", "precursor_charge": 1,
        "lcms_adduct": "[M+H]+", "ion_identity": ion.label,
        "ion_identity_status": "networked", "ion_identity_network_id": "network-1",
        "ion_identity_neutral_mass": 300,
        "ion_identity_evidence": {"status": "networked", "network_id": "network-1", "ion_identity": ion.label,
            "runtime_influence": True, "identity_conflict": False, "edge_count": 1,
            "best_shape_correlation": .95, "best_overlap_fraction": .8}}


def test_ms1_network_wins_over_unvalidated_lcms_default():
    row = network_row()
    original = deepcopy(row)
    result = resolve_query_ion(row)
    assert result["status"] == "supported"
    assert result["adduct"] == "[M+NH4]+"
    assert result["neutral_mass"] == pytest.approx(300)
    assert row == original


@pytest.mark.parametrize("edit, expected", [
    ({"precursor_charge": 2}, "charge_conflict"),
    ({"ion_mode": "NEG"}, "polarity_conflict"),
    ({"ion_identity_neutral_mass": 310}, "network_neutral_mass_conflict"),
    ({"query_adduct": "[M+Na]+", "query_adduct_reliable": True}, "unsupported_or_conflicting_independent_ions"),
    ({"query_adduct": "Positive", "query_adduct_reliable": True}, "unsupported_or_conflicting_independent_ions"),
    ({"ion_identity": "[M+Na]+"}, "conflicting_network_fields"),
])
def test_conflicts_fail_closed(edit, expected):
    row = network_row()
    row.update(edit)
    assert resolve_query_ion(row)["status"] == expected


@pytest.mark.parametrize("key, value", [("identity_conflict", True), ("runtime_influence", False),
    ("best_shape_correlation", .84), ("best_overlap_fraction", .49), ("edge_count", 0), ("network_id", "other")])
def test_bad_network_not_evidence(key, value):
    row = network_row()
    row["ion_identity_evidence"][key] = value
    assert resolve_query_ion(row)["adduct"] is None


def test_candidate_or_polarity_cannot_choose_ion():
    row = {"precursor_mz": 301, "ion_mode": "POS", "lcms_adduct": "[M+H]+",
           "match": {"adduct": "[M+Na]+", "formula": "C10H20O10"}}
    assert resolve_query_ion(row)["adduct"] is None


@pytest.mark.parametrize("key", ["truth_formula", "truth_connectivity", "acceptable_connectivities", "ground_truth"])
def test_truth_fields_rejected(key):
    row = network_row()
    row[key] = "anything"
    with pytest.raises(ValueError, match="Truth-bearing"):
        resolve_query_ion(row)


@pytest.mark.parametrize("value", ["[M+Hac-H]-", "[M+Ac-H]-", "[M+CH3COO]-"])
def test_acetate_aliases(value):
    assert adduct(value) == "[M+Ac-H]-"


def test_reference_windows_require_query_evidence():
    class NoCalls:
        def _fetch_window(self, *args):
            raise AssertionError("No query ion means no neutral retrieval")
    result = neutral_recall({"precursor_mz": 301, "ion_mode": "POS"}, NoCalls(), {})
    assert result["candidates"] == []


def test_different_reference_ion_never_grants_winner_permission():
    from matching.library_quality import formula_exact_mass
    mass = formula_exact_mass("C6H12O6")
    row = {"precursor_mz": mass + IONS["[M+NH4]+"].shift, "ion_mode": "POS",
           "query_adduct": "[M+NH4]+", "query_adduct_reliable": True}
    class Library:
        def _fetch_window(self, mz, mode):
            return [{"precursor_mz": mass + IONS["[M+H]+"].shift,
                "adduct": "[M+H]+", "library_polarity": "POS", "molecular_formula": "C6H12O6",
                "connectivity_identity": "AAAAAAAAAAAAAA", "spectrum_id": "experimental-1"}], {"status": "available"}
    policy = json.loads((Path(__file__).parents[1] / "eval/v16_ion_recall_repair/protocol.json").read_text())
    result = neutral_recall(row, Library(), policy)
    assert len(result["candidates"]) == 1
    assert result["candidates"][0]["same_query_ion_reference_ids"] == []
    assert result["can_replace_winner"] is False


def integration_fixture():
    from matching.library_quality import formula_exact_mass
    mass = formula_exact_mass("C6H12O6")
    mz = mass + IONS["[M+H]+"].shift
    native = {"inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-N", "formula": "C6H12O6",
        "adduct": "[M+H]+", "raw_cosine": .2, "raw_entropy": .2,
        "raw_reverse_match": .2, "raw_fragment_score": .2, "raw_intensity_score": .2,
        "mass_ppm": 0, "top_candidates": []}
    row = {"spectrum_id": "measured-1", "precursor_mz": mz, "ion_mode": "POS",
        "query_adduct": "[M+H]+", "query_adduct_reliable": True,
        "precursor_charge": 1, "collision_energy": 20, "collision_energy_unit": "eV",
        "raw_peaks": [[60., 100.], [90., 30.], [120., 50.]], "raw_peaks_complete": True,
        "match": native}
    candidates = [{"connectivity": "AAAAAAAAAAAAAA", "origin": "native"},
                  {"connectivity": "CCCCCCCCCCCCCC", "origin": "expanded"}]
    reference = {"inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-N", "molecular_formula": "C6H12O6",
        "connectivity_identity": "CCCCCCCCCCCCCC", "precursor_mz": mz, "adduct": "[M+H]+",
        "library_polarity": "POS", "spectrum_id": "measured-ref-1", "source_tier": "experimental",
        "dataset_source": "massbank", "peaks": [[60., 90.], [90., 32.], [120., 49.]]}
    class Library:
        def _fetch_window(self, *args):
            return [reference], {"status": "available"}
    return row, candidates, reference, Library()


def test_expanded_reference_enters_native_entity_scoring_without_mutation():
    from eval.v16_ion_recall_repair.core import inject_expanded_native
    row, candidates, _, retriever = integration_fixture()
    before = deepcopy(row)
    result = inject_expanded_native(row, candidates, retriever)
    assert result["expanded_scored_count"] == 1
    assert result["winner"] == "CCCCCCCCCCCCCC"
    assert result["winner_is_expanded"] is True
    assert row == before


@pytest.mark.parametrize("edit, reason", [
    ({"library_polarity": "NEG"}, "reference_polarity"),
    ({"adduct": "[M+Na]+"}, "different_or_unspecified_reference_ion"),
    ({"precursor_mz": 182.}, "direct_precursor_outside_10ppm"),
    ({"peaks": [[60., 90.]]}, "insufficient_fragment_peaks"),
])
def test_ineligible_expanded_references_cannot_win(edit, reason):
    from eval.v16_ion_recall_repair.core import inject_expanded_native
    row, candidates, reference, retriever = integration_fixture()
    reference.update(edit)
    result = inject_expanded_native(row, candidates, retriever)
    assert result["expanded_scored_count"] == 0
    assert not result["winner_is_expanded"]
    assert result["rejected"][reason] == 1


def test_query_mirror_is_not_new_reference_evidence():
    from eval.v16_ion_recall_repair.core import inject_expanded_native
    row, candidates, reference, retriever = integration_fixture()
    reference["peaks"] = deepcopy(row["raw_peaks"])
    result = inject_expanded_native(row, candidates, retriever)
    assert result["expanded_scored_count"] == 0
    assert result["rejected"]["query_content_mirror"] == 1


def test_targeted_projection_cannot_overwrite_validated_other_ion():
    from eval.v16_ion_recall_repair.targeted_ms1 import project_network_evidence
    row = network_row()
    before = deepcopy(row)
    pair = {"neutral_mass": 300., "shape": {"correlation": .99, "overlap_fraction": .99}}
    assert project_network_evidence(row, "query", "[M+H]+", pair, feature_id="peer") is None
    assert row == before


def test_targeted_projection_cannot_erase_network_conflict():
    from eval.v16_ion_recall_repair.targeted_ms1 import project_network_evidence
    row = network_row(); row["ion_identity_evidence"]["identity_conflict"] = True
    pair = {"neutral_mass": 300., "shape": {"correlation": .99, "overlap_fraction": .99}}
    assert project_network_evidence(row, "query", "[M+NH4]+", pair, feature_id="peer") is None


def test_targeted_pair_requires_actual_mass_not_target_value():
    from eval.v16_ion_recall_repair.targeted_ms1 import assess_pairs
    shape = [1., 3., 5., 10., 5., 3., 1.]
    traces = {"q": {"rt": list(range(7)), "intensity": shape,
                      "observed_mz": [301.007276]*7, "target_mz": 301.007276},
              "p": {"rt": list(range(7)), "intensity": shape,
                      "observed_mz": [319.033823]*7, "target_mz": 318.033823}}
    pairs = [{"partner_target": "p", "query_ion": "[M+H]+", "partner_ion": "[M+NH4]+"}]
    result = assess_pairs("q", pairs, traces)
    assert result["status"] == "no_supported_pair"
    assert result["rejected"]["observed_neutral_mass_mismatch"] == 1


def test_native_signal_floor_blocks_mass_only_winner():
    from eval.v16_ion_recall_repair.network_peer_signal_guard import score_same_ion_row
    from matching.v15_flagship_ranker import load_v15_flagship_policy
    row, _, reference, retriever = integration_fixture()
    row["raw_peaks"] = [[11.12345,100.], [22.23456,40.], [33.34567,20.]]
    scored, audit = score_same_ion_row(row, {reference["connectivity_identity"]}, retriever, load_v15_flagship_policy())
    assert scored["match"] is None
    assert audit["eligible_reference_count"] == 0
    assert audit["rejected"]["below_existing_native_search_cosine_threshold"] == 1


def test_complete_panel_keeps_low_cells_but_disallows_signal_free_promotion():
    from eval.v16_ion_recall_repair.network_peer_complete_signal_gate import score_same_ion_row, winner_has_measured_signal
    from matching.v15_flagship_ranker import load_v15_flagship_policy
    row, _, reference, retriever = integration_fixture()
    row["raw_peaks"] = [[11.12345,100.], [22.23456,40.], [33.34567,20.]]
    scored, audit = score_same_ion_row(row, {reference["connectivity_identity"]}, retriever, load_v15_flagship_policy())
    assert scored["match"] is not None  # measured zero is not a missing observation
    assert scored["match"]["raw_cosine"] == 0
    assert not winner_has_measured_signal({"winner": {"collision_energy_evidence": [
        {"conditioned_multilibrary_evidence": {"reference_records": [{"score_components": {"cosine": 0.}}]}}]}})
    assert winner_has_measured_signal({"winner": {"collision_energy_evidence": [
        {"conditioned_multilibrary_evidence": {"reference_records": [{"score_components": {"cosine": .4}}]}}]}})
