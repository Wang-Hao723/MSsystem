from __future__ import annotations

from copy import deepcopy

from matching.evograph_fragmentation import (
    build_evograph_graph,
    canonical_explicit_adduct,
)


def _glucose_row() -> dict:
    return {
        "spectrum_id": "glucose-water-loss",
        "precursor_mz": 181.070664,
        "ion_mode": "POS",
        "precursor_charge": 1,
        "collision_energy": 20,
        "raw_peaks": [[163.060100, 100.0], [135.065200, 30.0]],
        "match": {
            "formula": "C6H12O6",
            "adduct": "[M+H]+",
            "final_score": 0.91,
            "confidence": "High",
        },
    }


def test_evograph_builds_conserving_water_loss_without_runtime_influence():
    row = _glucose_row()
    before = deepcopy(row)

    graph = build_evograph_graph(row)

    selected = [edge for edge in graph["edges"] if edge["selected"]]
    water = next(edge for edge in selected if edge["loss_formula"] == "H2O")
    assert water["parent_formula"] == "C6H12O6"
    assert water["child_formula"] == "C6H10O5"
    assert water["formula_conservation"] is True
    assert abs(water["mass_error_da"]) < 0.001
    assert water["calibrated"] is False
    assert graph["graph_audit"]["cycle_count"] == 0
    assert graph["production_ranking_changed"] is False
    assert graph["legacy_confidence_changed"] is False
    assert graph["high_promotions"] == 0
    assert graph["calibration"]["log_bayes_factor"] is None
    assert row == before


def test_evograph_does_not_guess_adduct_from_polarity_or_generic_label():
    row = _glucose_row()
    row["match"]["adduct"] = "Positive"

    graph = build_evograph_graph(row)

    assert canonical_explicit_adduct("Positive") is None
    assert graph["context"]["adduct"] is None
    assert graph["context"]["formula_assignment_status"] == "unavailable"
    assert graph["formula_graph_available"] is False


def test_evograph_abstains_from_formula_assignment_for_multicharge_precursor():
    row = _glucose_row()
    row["precursor_charge"] = 2

    graph = build_evograph_graph(row)

    assert graph["context"]["formula_assignment_status"] == "unavailable"
    assert "unsupported_precursor_charge" in graph["context"]["reasons"]
    assert graph["formula_graph_available"] is False


def test_evograph_rejects_ion_mode_and_explicit_adduct_polarity_conflict():
    row = _glucose_row()
    row["ion_mode"] = "NEG"

    graph = build_evograph_graph(row)

    assert graph["context"]["formula_assignment_status"] == "contradicted"
    assert "ion_mode_adduct_polarity_inconsistent" in graph["context"]["reasons"]
    assert graph["formula_graph_available"] is False


def test_evograph_is_deterministic_and_marks_library_formula_candidate_conditioned():
    first = build_evograph_graph(_glucose_row())
    second = build_evograph_graph(_glucose_row())

    assert first == second
    assert first["context"]["formula_source"] == "selected_library_candidate"
    assert first["context"]["formula_candidate_independent"] is False
    assert all(
        edge["posterior_semantics"]
        == "normalized_local_heuristic_weight_not_probability"
        for edge in first["edges"]
    )


def test_evograph_uses_measured_ms1_search_mass_and_preserves_selected_ion():
    row = _glucose_row()
    row["precursor_mz"] = 181.0689
    row["search_precursor_mz"] = 181.070664
    row["search_precursor_mz_source"] = (
        "measured_continuous_ms1_feature_centroid"
    )
    row["lcms_feature_used_for_search"] = True

    graph = build_evograph_graph(row)

    assert graph["precursor_mz"] == 181.070664
    assert graph["search_precursor_mz"] == 181.070664
    assert graph["selected_ion_mz"] == 181.0689
    assert graph["lcms_feature_used_for_search"] is True
    assert abs(graph["context"]["mass_consistency_error_ppm"]) < 1.0
    assert graph["nodes"][0]["mz"] == 181.070664
