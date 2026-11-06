from __future__ import annotations

import copy

import pytest

from matching.conditioned_multilibrary import aggregate_references


POLICY = {
    "best_reference_fraction": 0.5,
    "consensus_fraction": 0.5,
    "minimum_provenance_strata": 2,
    "minimum_unique_content_profiles": 2,
    "near_collision_energy_delta": {"eV": 5.0, "NCE": 5.0},
}

QUERY = {
    "instrument_type": "qtof",
    "instrument": "Agilent 6546 Q-TOF",
    "collision_energy": 20,
    "collision_energy_unit": "eV",
    "collision_energy_mode": "fixed",
    "activation_method": "CID",
}

PEAKS_A = [[50.0, 10.0], [75.0, 100.0], [90.0, 25.0]]
PEAKS_B = [[50.0, 20.0], [75.0, 60.0], [91.0, 100.0]]
PEAKS_C = [[52.0, 25.0], [78.0, 100.0], [93.0, 30.0]]


def _reference(
    score: float,
    peaks: list[list[float]],
    source: str,
    spectrum_id: str,
) -> dict:
    return {
        "spectral_score": score,
        "dataset_source": source,
        "source_tier": "experimental",
        "lib_spectrum_id": spectrum_id,
        "lib_source_accession": spectrum_id,
        "lib_peaks": peaks,
        "lib_instrument_type": "LC-ESI-QTOF",
        "lib_instrument": "Agilent 6546 Q-TOF",
        "lib_collision_energy": 20,
        "lib_collision_energy_unit": "eV",
        "lib_collision_energy_mode": "fixed",
        "lib_activation_method": "CID",
    }


def _mirror(record: dict, source: str, spectrum_id: str) -> dict:
    mirrored = copy.deepcopy(record)
    mirrored["dataset_source"] = source
    mirrored["lib_spectrum_id"] = spectrum_id
    mirrored["lib_source_accession"] = spectrum_id
    return mirrored


def _summary(result: dict) -> tuple:
    return (
        result["status"],
        result["conditioned_multilibrary_score"],
        result["provenance_consensus_score"],
        result["provenance_stratum_count"],
        result["conditioned_provenance_stratum_count"],
        result["rank_influence_allowed"],
        result["selected_stratum_scores"],
        result["display_only_all_condition_stratum_scores"],
    )


def test_multiple_profiles_mirrored_from_one_source_do_not_create_extra_votes():
    massbank_a = _reference(0.90, PEAKS_A, "MassBank_NIST", "mb-a")
    massbank_b = _reference(0.70, PEAKS_B, "MassBank_NIST", "mb-b")
    before = aggregate_references(QUERY, [massbank_a, massbank_b], POLICY)
    after = aggregate_references(
        QUERY,
        [
            massbank_a,
            massbank_b,
            _mirror(massbank_a, "ALL_GNPS", "gnps-a"),
            _mirror(massbank_b, "ALL_GNPS", "gnps-b"),
        ],
        POLICY,
    )

    assert before["provenance_stratum_count"] == 1
    assert after["provenance_stratum_count"] == 1
    assert after["conditioned_provenance_stratum_count"] == 1
    assert after["conditioned_multilibrary_score"] == pytest.approx(0.90)
    assert after["rank_influence_allowed"] is False


def test_true_independent_source_still_contributes_one_additional_vote():
    massbank_a = _reference(0.92, PEAKS_A, "MassBank_NIST", "mb-a")
    references = [
        massbank_a,
        _mirror(massbank_a, "ALL_GNPS", "gnps-a"),
        _reference(0.70, PEAKS_B, "MoNA-export", "mona-b"),
    ]

    result = aggregate_references(QUERY, references, POLICY)

    assert result["unique_content_profile_count"] == 2
    assert result["provenance_stratum_count"] == 2
    assert result["conditioned_provenance_stratum_count"] == 2
    assert result["conditioned_multilibrary_score"] == pytest.approx(0.865)
    assert result["rank_influence_allowed"] is True


def test_reference_input_order_does_not_change_grouping_or_score():
    massbank_a = _reference(0.90, PEAKS_A, "MassBank_NIST", "mb-a")
    massbank_b = _reference(0.70, PEAKS_B, "MassBank_NIST", "mb-b")
    references = [
        massbank_a,
        _mirror(massbank_a, "ALL_GNPS", "gnps-a"),
        massbank_b,
        _mirror(massbank_b, "ALL_GNPS", "gnps-b"),
        _reference(0.60, PEAKS_C, "MoNA-export", "mona-c"),
    ]

    forward = aggregate_references(QUERY, references, POLICY)
    reverse = aggregate_references(QUERY, list(reversed(references)), POLICY)

    assert _summary(forward) == _summary(reverse)


def test_partial_cross_source_mirror_cannot_increase_linked_group_vote_count():
    massbank_a = _reference(0.90, PEAKS_A, "MassBank_NIST", "mb-a")
    massbank_b = _reference(0.70, PEAKS_B, "MassBank_NIST", "mb-b")
    gnps_c = _reference(0.60, PEAKS_C, "ALL_GNPS", "gnps-c")
    before = aggregate_references(
        QUERY, [massbank_a, massbank_b, gnps_c], POLICY
    )
    after = aggregate_references(
        QUERY,
        [
            massbank_a,
            massbank_b,
            gnps_c,
            _mirror(massbank_a, "ALL_GNPS", "gnps-mirror-a"),
        ],
        POLICY,
    )

    assert before["provenance_stratum_count"] == 2
    assert before["rank_influence_allowed"] is True
    assert after["provenance_stratum_count"] == 1
    assert after["conditioned_provenance_stratum_count"] == 1
    assert after["rank_influence_allowed"] is False
    assert after["conditioned_multilibrary_score"] == pytest.approx(0.90)
