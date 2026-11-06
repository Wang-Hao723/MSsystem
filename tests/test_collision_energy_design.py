import json
from pathlib import Path

from acquisition.action_space import candidate_actions
from acquisition.ambiguity_graph import build_acquisition_needs
from acquisition.energy_design import build_collision_energy_recommendation


def _policy():
    return json.loads(Path("config/acquisition_policy.json").read_text())


def _row(*, unit="eV", energy=15.0, instrument="Agilent Q-TOF"):
    return {
        "spectrum_id": "scan-1",
        "precursor_mz": 325.0677,
        "retention_time": 4.2,
        "ion_mode": "NEG",
        "instrument_type": "tof",
        "instrument": instrument,
        "collision_energy": energy,
        "collision_energy_unit": unit,
        "precursor_feature_collision_energies": [energy],
        "query_peaks": [[125.0, 100.0], [179.0, 40.0]],
        "match": {
            "compound_name": "A",
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "formula": "C15H18O8",
            "rank_score": 0.80,
            "lib_peaks": [[125.0, 100.0], [151.0, 25.0]],
            "top_candidates": [{
                "compound_name": "B",
                "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                "formula": "C15H18O8",
                "rank_score": 0.77,
                "lib_peaks": [[125.0, 100.0], [179.0, 35.0]],
            }],
        },
    }


def test_agilent_qtof_single_energy_uses_sequential_wide_bracket():
    recommendation = build_collision_energy_recommendation(_row(), _policy())

    assert recommendation["initial_energies"] == (10.0, 20.0, 40.0)
    design = recommendation["design"]
    assert design["evidence_status"] == "single_energy_only_no_response_curve"
    assert design["observed_collision_energies"] == (15.0,)
    assert design["automatic_ceiling"] == 60.0
    stages = {stage["stage_id"]: stage for stage in design["stages"]}
    assert stages["labile_low_energy_refinement"]["collision_energies"] == (5.0,)
    assert stages["transition_midpoint_refinement"]["collision_energies"] == (30.0,)
    assert stages["resistant_high_energy_extension"]["collision_energies"] == (60.0,)
    assert design["universal_complete_fragmentation_assumed"] is False
    assert design["cross_instrument_conversion_applied"] is False


def test_current_15_ev_is_not_expanded_to_unvalidated_5_15_25_triplet():
    row = _row()
    need = build_acquisition_needs("task", [row])[0]
    action = next(
        proposal for proposal in candidate_actions(need, row, _policy())
        if proposal["action_kind"] == "targeted_ms2_multi_ce"
    )

    assert action["collision_energies"] == (10.0, 20.0, 40.0)
    assert action["collision_energies"] != (5.0, 15.0, 25.0)


def test_nce_domain_never_receives_agilent_ev_numbers_or_conversion():
    recommendation = build_collision_energy_recommendation(
        _row(unit="NCE", energy=30.0), _policy()
    )

    assert recommendation["initial_energies"] == ()
    assert recommendation["collision_energy_unit"] == "NCE"
    assert recommendation["design"]["numeric_schedule_available"] is False
    assert recommendation["design"]["cross_instrument_conversion_applied"] is False


def test_opposite_polarity_does_not_reuse_observed_energy_curve():
    row = _row()
    need = build_acquisition_needs("task", [row])[0]
    action = next(
        proposal for proposal in candidate_actions(need, row, _policy())
        if proposal["action_kind"] == "opposite_polarity_targeted_ms2"
    )

    assert action["polarity_override"] == "POS"
    assert action["collision_energy_design"]["observed_collision_energies"] == ()
    assert action["collision_energy_design"]["existing_collision_energy"] is None
    assert action["collision_energy_design"]["evidence_status"] == (
        "no_observed_energy_response_curve"
    )


def test_agilent_qtof_uses_native_medium_and_narrow_isolation_settings():
    row = _row()
    row["match"]["chimeric_deconvolution"] = {"isolation_status": "mixed"}
    need = build_acquisition_needs("task", [row])[0]
    actions = candidate_actions(need, row, _policy())
    targeted = next(
        action for action in actions
        if action["action_kind"] == "targeted_ms2_multi_ce"
    )
    narrow = next(
        action for action in actions
        if action["action_kind"] == "narrow_isolation_ms2"
    )

    assert targeted["isolation_width_da"] == 4.0
    assert narrow["isolation_width_da"] == 1.3
    assert narrow["isolation_width_da"] != 0.4
    assert narrow["isolation_width_profile_id"] == (
        "agilent_qtof_native_isolation_v1.0.0"
    )


def test_unknown_instrument_domain_does_not_get_numeric_narrow_width():
    row = _row(unit="NCE", energy=30.0, instrument="unqualified instrument")
    row["instrument_type"] = "unknown"
    row["match"]["chimeric_deconvolution"] = {"isolation_status": "mixed"}
    need = build_acquisition_needs("task", [row])[0]
    narrow = next(
        action for action in candidate_actions(need, row, _policy())
        if action["action_kind"] == "narrow_isolation_ms2"
    )

    assert narrow["isolation_width_da"] is None
    assert narrow["isolation_width_evidence_status"] == (
        "instrument_domain_unqualified_no_numeric_width"
    )
