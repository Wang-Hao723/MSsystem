"""Synthetic checks for missing evidence and physical numeric inputs."""

import math

import pytest

from acquisition.energy_design import build_collision_energy_recommendation
from acquisition.observable_catalog import (
    candidate_observables,
    pair_fragment_separation,
)


@pytest.mark.parametrize("left,right", [([], [100.0, 200.0]), ([100.0], []), ([], [])])
def test_missing_candidate_fragments_never_become_discriminating_evidence(left, right):
    catalog = {
        "left": {"diagnostic_fragments": left},
        "right": {"diagnostic_fragments": right},
    }

    assert pair_fragment_separation("left", "right", catalog) == (0.0, ())


def test_missing_candidate_record_has_no_separation_credit():
    catalog = {"right": {"diagnostic_fragments": [100.0, 200.0]}}

    assert pair_fragment_separation("absent", "right", catalog) == (0.0, ())


def test_bilateral_fragments_retain_the_existing_set_difference_score():
    catalog = {
        "left": {"diagnostic_fragments": [100.0, 150.0]},
        "right": {"diagnostic_fragments": [100.0, 175.0]},
    }

    score, ions = pair_fragment_separation("left", "right", catalog)

    assert score == pytest.approx(2.0 / 3.0)
    assert ions == (150.0, 175.0)


def test_flat_mass_lists_keep_every_mass_and_nested_peaks_exclude_intensities():
    observable = candidate_observables({
        "diagnostic_ions": [100.0, 200.0, 300.0],
        "lib_peaks": [[120.0, 999.0], [150.0, 888.0]],
        "neutral_losses": [18.0, 44.0],
        "isomer_fragment_detail": {
            "candidate_specific_product_ions": [{"mz": 250.0}, {"product_mz": 275.0}],
        },
    })

    assert observable["diagnostic_fragments"] == [100.0, 120.0, 150.0, 200.0, 250.0, 275.0, 300.0]
    assert observable["neutral_losses"] == [18.0, 44.0]
    assert observable["confirmatory"] is False


def test_nonfinite_fragment_masses_are_ignored_in_all_supported_encodings():
    observable = candidate_observables({
        "diagnostic_ions": [100.0, math.inf, math.nan, -math.inf],
        "lib_peaks": [[math.inf, 10.0], [150.0, 20.0]],
        "product_ions": [{"mz": math.nan}, {"mz": math.inf}, {"mz": 175.0}],
        "neutral_losses": [math.inf, 18.0, math.nan],
    })

    assert observable["diagnostic_fragments"] == [100.0, 150.0, 175.0]
    assert observable["neutral_losses"] == [18.0]
    assert pair_fragment_separation("left", "right", {
        "left": {"diagnostic_fragments": [math.inf, math.nan]},
        "right": observable,
    }) == (0.0, ())


def test_nonfinite_observed_and_configured_collision_energies_are_not_exported():
    policy = {"collision_energy_design": {"agilent_qtof_ev": {
        "initial_scout_energies_eV": [10.0, math.inf, math.nan, 20.0],
        "conditional_low_energy_eV": [-math.inf, 5.0],
        "conditional_midpoint_energy_eV": [30.0, math.inf],
        "conditional_high_energy_eV": [60.0, math.nan],
        "automatic_ceiling_eV": math.inf,
    }}}
    recommendation = build_collision_energy_recommendation({
        "instrument": "Agilent Q-TOF",
        "collision_energy_unit": "eV",
        "collision_energy": math.inf,
        "precursor_feature_collision_energies": [math.inf, math.nan, 15.0],
    }, policy)

    design = recommendation["design"]
    assert recommendation["initial_energies"] == (10.0, 20.0)
    assert design["observed_collision_energies"] == (15.0,)
    assert design["existing_collision_energy"] == 15.0
    assert design["automatic_ceiling"] is None
    assert all(
        math.isfinite(energy)
        for stage in design["stages"]
        for energy in stage["collision_energies"]
    )


def test_qe_nce_with_nonfinite_energy_still_has_no_qtof_schedule():
    recommendation = build_collision_energy_recommendation({
        "instrument": "Thermo Q Exactive Orbitrap",
        "collision_energy_unit": "NCE",
        "collision_energy": math.inf,
    }, {})

    assert recommendation["initial_energies"] == ()
    assert recommendation["collision_energy_unit"] == "NCE"
    assert recommendation["design"]["existing_collision_energy"] is None
    assert recommendation["design"]["observed_collision_energies"] == ()
    assert recommendation["design"]["numeric_schedule_available"] is False
    assert recommendation["design"]["cross_instrument_conversion_applied"] is False
