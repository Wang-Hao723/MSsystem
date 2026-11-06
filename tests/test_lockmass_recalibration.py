from __future__ import annotations

import pytest

from matching.lockmass_recalibration import (
    AGILENT_NEGATIVE_REFERENCE_MZ,
    build_selected_match_counterfactual,
    corrected_mz,
    summarize_lockmass_observations,
)


def test_stable_candidate_independent_lockmass_qualifies_for_shadow(tmp_path):
    reference = AGILENT_NEGATIVE_REFERENCE_MZ
    observations = []
    for index in range(60):
        ppm = 15.0 + ((index % 5) - 2) * 0.4
        observations.append({
            "observed_mz": reference * (1.0 + ppm * 1e-6),
            "intensity": 2000.0 + index,
            "ppm_error": ppm,
            "retention_time_min": index / 10,
        })

    artifact = summarize_lockmass_observations(
        observations,
        task_id="task",
        source_mzml=tmp_path / "sample.mzML",
        reference_mz=reference,
        ms1_scan_count=100,
        search_window_da=0.006,
        minimum_intensity=1000.0,
    )

    assert artifact["status"] == "qualified_shadow"
    assert artifact["candidate_independent"] is True
    assert artifact["chemical_identity_used_for_fit"] is False
    assert artifact["runtime_influence"] is False
    assert artifact["production_authorized"] is False
    observed = 500.0 * (1.0 + 15.0e-6)
    assert corrected_mz(observed, artifact) == pytest.approx(500.0, abs=1e-6)
    results = [{
        "spectrum_id": "scan=1",
        "match": {
            "compound_name": "A",
            "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "confidence": "High",
            "mass_query_value": observed,
            "mass_library_value": 500.0,
        },
    }]
    counterfactual = build_selected_match_counterfactual(results, artifact)
    assert counterfactual["summary"]["before_within_10ppm"] == 0
    assert counterfactual["summary"]["after_within_10ppm"] == 1
    assert counterfactual["ranking_changed"] is False


def test_unstable_or_missing_lockmass_fails_closed(tmp_path):
    artifact = summarize_lockmass_observations(
        [],
        task_id="task",
        source_mzml=tmp_path / "sample.mzML",
        reference_mz=AGILENT_NEGATIVE_REFERENCE_MZ,
        ms1_scan_count=100,
        search_window_da=0.006,
        minimum_intensity=1000.0,
    )

    assert artifact["status"] == "not_qualified"
    assert artifact["correction"] is None
    assert corrected_mz(500.0, artifact) is None
