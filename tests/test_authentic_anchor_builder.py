from __future__ import annotations

from scripts.build_authentic_anchor_library import _select_references


def _row(
    spectrum_id: str,
    *,
    rt: float,
    ce: float,
    signal: float,
    adduct: str = "[M-H]-",
) -> dict:
    return {
        "spectrum_id": spectrum_id,
        "retention_time": rt,
        "collision_energy": ce,
        "tic": signal,
        "precursor_ms1_intensity": signal / 2,
        "search_precursor_mz": 285.04,
        "anchor_adduct": adduct,
        "anchor_mass_error_abs_ppm": 1.0,
        "query_peaks": [[80.0, signal], [100.0, signal / 2], [120.0, signal / 3]],
    }


def test_reference_selection_never_splices_an_intense_single_ce_background_scan():
    rows = [
        _row("target-10", rt=4.000, ce=10, signal=1_000),
        _row("target-20", rt=4.006, ce=20, signal=1_200),
        _row("target-40", rt=4.012, ce=40, signal=900),
        _row("background-40", rt=4.500, ce=40, signal=1_000_000),
    ]

    spectra, audit = _select_references(
        rows,
        1,
        rt_cluster_tolerance_min=0.08,
    )

    assert {row["source_spectrum_id"] for row in spectra} == {
        "target-10",
        "target-20",
        "target-40",
    }
    assert audit["selected_collision_energies_ev"] == [10.0, 20.0, 40.0]
    assert audit["discarded_mass_compatible_spectrum_count"] == 1


def test_high_signal_two_ce_cluster_can_beat_low_signal_three_ce_background():
    rows = [
        _row("low-10", rt=4.000, ce=10, signal=100),
        _row("low-20", rt=4.006, ce=20, signal=100),
        _row("low-40", rt=4.012, ce=40, signal=100),
        _row("target-20", rt=5.000, ce=20, signal=100_000),
        _row("target-40", rt=5.006, ce=40, signal=120_000),
    ]

    spectra, audit = _select_references(
        rows,
        1,
        rt_cluster_tolerance_min=0.08,
    )

    assert {row["source_spectrum_id"] for row in spectra} == {
        "target-20",
        "target-40",
    }
    assert audit["selected_rt_cluster_center_min"] == 5.003


def test_reference_clusters_do_not_mix_adducts_at_the_same_retention_time():
    rows = [
        _row("deprotonated-10", rt=4.000, ce=10, signal=1_000),
        _row("formate-20", rt=4.005, ce=20, signal=2_000, adduct="[M+HCOO]-"),
    ]

    spectra, audit = _select_references(
        rows,
        1,
        rt_cluster_tolerance_min=0.08,
    )

    assert len(spectra) == 1
    assert {row["adduct"] for row in spectra} == {audit["selected_adduct"]}
    assert audit["candidate_cluster_count"] == 2
