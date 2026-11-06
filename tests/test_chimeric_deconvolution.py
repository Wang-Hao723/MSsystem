import math
import importlib

import numpy as np
import pytest

from matching.chimeric_deconvolution import (
    deconvolve_chimeric_spectrum,
    nonnegative_lasso,
)
from matching.honest_calibration import annotate as honest_annotate


search_module = importlib.import_module("matching.search")


def _mix(left, right, left_weight=1.0, right_weight=1.0):
    values = {}
    for mz, intensity in left:
        values[mz] = values.get(mz, 0.0) + intensity * left_weight
    for mz, intensity in right:
        values[mz] = values.get(mz, 0.0) + intensity * right_weight
    return sorted([[mz, value] for mz, value in values.items()])


TARGET = [[80.0, 15.0], [100.0, 100.0], [120.0, 50.0], [140.0, 25.0]]
CONTAMINANT = [[80.0, 10.0], [90.0, 80.0], [110.0, 100.0], [130.0, 40.0]]
TARGET_PRECURSOR_MZ = 493.57067832056  # C35H72 [M+H]+
COISOLATED_PRECURSOR_MZ = TARGET_PRECURSOR_MZ + 0.5


def test_nonnegative_lasso_requires_positive_l1_penalty():
    with pytest.raises(ValueError):
        nonnegative_lasso(np.eye(2), np.ones(2), penalty=0.0)


def test_extreme_finite_intensities_are_scaled_before_bin_accumulation():
    huge_target = [
        [100.0000, 1.0e308],
        [100.0010, 1.0e308],
        [120.0000, 5.0e307],
        [140.0000, 2.5e307],
    ]
    huge_contaminant = [
        [90.0000, 1.0e308],
        [90.0010, 1.0e308],
        [110.0000, 8.0e307],
        [130.0000, 4.0e307],
    ]
    query = [
        *[[mz, intensity * 0.70] for mz, intensity in huge_target],
        *[[mz, intensity * 0.45] for mz, intensity in huge_contaminant],
    ]
    with np.errstate(all="raise"):
        result = deconvolve_chimeric_spectrum(
            query,
            [
                {"id": "target", "peaks": huge_target},
                {"id": "contaminant", "peaks": huge_contaminant},
            ],
            "target",
        )
    assert result["status"] != "numerical_instability"
    assert np.isfinite(result["reconstruction_cosine"])
    assert all(
        np.isfinite(component["coefficient"])
        for component in result["components"]
    )


def test_nonnegative_lasso_rejects_nonfinite_input():
    with pytest.raises(ValueError, match="finite"):
        nonnegative_lasso(
            np.asarray([[1.0, np.inf], [0.0, 1.0]]),
            np.ones(2),
            penalty=0.02,
        )


def test_sparse_mixture_is_deconvolved_and_returns_observed_allocation():
    query = _mix(TARGET, CONTAMINANT, 0.70, 0.45)
    result = deconvolve_chimeric_spectrum(
        query,
        [
            {"id": "target", "peaks": TARGET},
            {"id": "contaminant", "peaks": CONTAMINANT},
            {"id": "irrelevant", "peaks": [[70, 100], [150, 90], [160, 50]]},
        ],
        "target",
    )
    assert result["accepted"] is True
    assert result["active_component_count"] == 2
    assert result["reconstruction_cosine"] > 0.98
    assert 0.40 < result["target_fraction"] < 0.80
    purified_mz = {round(peak[0], 1) for peak in result["purified_target_peaks"]}
    assert {100.0, 120.0, 140.0}.issubset(purified_mz)
    assert 90.0 not in purified_mz
    assert result["runtime_influence"] is True


def test_pure_spectrum_does_not_claim_chimeric_improvement():
    result = deconvolve_chimeric_spectrum(
        TARGET,
        [
            {"id": "target", "peaks": TARGET},
            {"id": "contaminant", "peaks": CONTAMINANT},
        ],
        "target",
    )
    assert result["accepted"] is False
    assert result["status"] == "single_component_solution"
    assert result["runtime_influence"] is False


def test_null_noise_does_not_pass_evidence_gate():
    noise = [[61.0, 100.0], [73.0, 55.0], [151.0, 80.0], [177.0, 35.0]]
    result = deconvolve_chimeric_spectrum(
        noise,
        [
            {"id": "target", "peaks": TARGET},
            {"id": "contaminant", "peaks": CONTAMINANT},
        ],
        "target",
    )
    assert result["accepted"] is False
    assert result["reconstruction_cosine"] < 0.75


def test_redundant_component_is_flagged_non_unique():
    duplicate = [[mz + 0.001, intensity * 2] for mz, intensity in TARGET]
    query = _mix(TARGET, CONTAMINANT, 0.70, 0.45)
    result = deconvolve_chimeric_spectrum(
        query,
        [
            {"id": "target", "peaks": TARGET},
            {"id": "duplicate", "peaks": duplicate},
            {"id": "contaminant", "peaks": CONTAMINANT},
        ],
        "target",
    )
    assert result["accepted"] is False
    assert result["status"] == "non_unique_components"
    assert math.isclose(result["target_max_component_similarity"], 1.0, abs_tol=1e-6)


def _search_record(index, precursor_mz, peaks, inchikey):
    return {
        "precursor_mz": precursor_mz,
        "neutral_mass": None,
        "neutral_mass_source": "",
        "peaks": peaks,
        "clean_peaks": peaks,
        "spectrum_id": f"spec-{index}",
        "compound_name": f"compound-{index}",
        "source_compound_name": "",
        "name_resolution_status": "",
        # The production search now rejects a full structure whose formula /
        # adduct cannot physically produce the reported precursor. Keep the
        # synthetic deconvolution target chemically self-consistent so this
        # test exercises deconvolution rather than the precursor safety gate.
        "smiles": "C" * 35 if index == 0 else "CCCO",
        "has_structure": True,
        "name_only_lead": False,
        "inchikey": inchikey,
        "reported_inchikey": inchikey,
        "formula": "C35H72" if index == 0 else "C3H8O",
        "adduct": "[M+H]+",
        "charge": 1,
        "collision_energy": None,
        "retention_time": None,
        "classification": "other",
        "classification_source": "test",
        "classification_confidence": None,
        "dataset_source": "test",
        "source_tier": "experimental",
        "intensity_origin": "measured",
        "quality_flags": [],
        "pathway_json": "",
        "polarity": "POS",
    }


def _disable_optional_search_layers(monkeypatch):
    monkeypatch.setattr(search_module, "EXPERT_SYSTEM_ENABLED", False)
    monkeypatch.setattr(search_module, "USE_EXPERT_V3", False)
    monkeypatch.setattr(search_module, "USE_PATHWAY_ALIGN", False)
    monkeypatch.setattr(search_module, "USE_CONSISTENCY_FILTER", False)
    monkeypatch.setattr(search_module, "USE_ADAPTIVE_PEAK_SEARCH", False)
    monkeypatch.setattr(search_module, "USE_ISOMER_FRAGMENT_RERANK", False)
    monkeypatch.setattr(search_module, "MODEL_PREDICTOR_AVAILABLE", False)
    monkeypatch.setattr(search_module, "_ensure_clf", lambda: False)


def test_search_applies_guarded_deconvolution_and_can_lift_medium_to_high(monkeypatch):
    _disable_optional_search_layers(monkeypatch)
    records = [
        _search_record(
            0, TARGET_PRECURSOR_MZ, TARGET, "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
        ),
        _search_record(
            1, COISOLATED_PRECURSOR_MZ, CONTAMINANT, "CCCCCCCCCCCCCC-DDDDDDDDDD-E"
        ),
    ]
    mz_index = {}
    for index, record in enumerate(records):
        bucket = int(record["precursor_mz"] / search_module.MZ_BUCKET_DALTON)
        mz_index.setdefault(bucket, []).append(index)
    query = _mix(TARGET, CONTAMINANT, 1.0, 0.5)
    isolation = {
        "status": "severe",
        "purity": 0.42,
        "isolation_window_source": "test",
        "strongest_coisolated_peaks": [[COISOLATED_PRECURSOR_MZ, 800.0]],
    }
    result = search_module.search(
        query,
        TARGET_PRECURSOR_MZ,
        records=records,
        mz_index=mz_index,
        query_polarity="POS",
        query_charge=1,
        query_isolation_evidence=isolation,
        enable_chimeric_deconvolution=True,
        cos_thresh=0.0,
        final_thresh=0.0,
        strict_evidence_gate=False,
    )
    assert result is not None
    evidence = result["chimeric_deconvolution"]
    assert evidence["status"] == "accepted_score_applied"
    assert result["cosine"] > result["raw_cosine"] + 0.03
    assert result["final_score"] > result["raw_final_score"]
    assert result["rank_score"] == result["raw_rank_score"]
    assert result["top_candidates"][0]["cosine"] == pytest.approx(
        result["raw_cosine"], abs=1e-4
    )
    assert (
        result["top_candidates"][0]["candidate_comparison_source"]
        == "raw_measured_spectrum"
    )
    assert evidence["top1_identity_preserved"] is True
    assert evidence["deconvolved_candidate_margin_used"] is False
    assert evidence["candidate_comparison_policy"] == "raw_measured_spectrum_only"
    assert evidence["all_score_increase_gates_passed"] is True
    for gate in (
        "reconstruction_quality",
        "regularization_stability",
        "target_fraction",
        "non_redundancy",
        "spectral_gain",
    ):
        assert evidence["score_increase_gates"][gate] is True
    assert honest_annotate(result, result.get("top_candidates"))[
        "compound_confidence"
    ] == "High"


def test_failed_deconvolution_keeps_raw_high_without_score_gain_and_warns():
    best = {
        "cosine": 0.97,
        "final_score": 0.80,
        "rank_score": 0.80,
        "mass_ppm": 1.0,
        "reverse_match": 0.90,
        "compound_name": "Example",
        "smiles": "CCO",
        "has_structure": True,
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "formula": "C2H6O",
        "classification": "other",
        "source_tier": "experimental",
        "top_candidates": [],
        "chimeric_deconvolution": {
            "status": "poor_reconstruction",
            "accepted": False,
            "runtime_influence": False,
            "isolation_status": "severe",
        },
    }
    annotation = honest_annotate(best, [])
    assert annotation["compound_confidence"] == "High"
    assert annotation["class_confidence"] == "Low"
    assert annotation["deconvolution_score_influence"] is False
    assert annotation["deconvolution_evidence_state"] == "no_score_change"
    assert "poor_reconstruction" in annotation["chimeric_isolation_warning"]
    assert "no deconvolution score gain" in annotation["chimeric_isolation_warning"]
    assert annotation["evidence_tier"] == "High (raw MS2; chimeric warning)"


def test_accepted_deconvolution_is_the_only_score_influencing_state():
    best = {
        "cosine": 0.97,
        "raw_cosine": 0.91,
        "final_score": 0.80,
        "raw_final_score": 0.74,
        "rank_score": 0.76,
        "raw_rank_score": 0.76,
        "mass_ppm": 1.0,
        "reverse_match": 0.90,
        "compound_name": "Example",
        "smiles": "CCO",
        "has_structure": True,
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "formula": "C2H6O",
        "classification": "other",
        "source_tier": "experimental",
        "top_candidates": [],
        "chimeric_deconvolution": {
            "status": "accepted_score_applied",
            "accepted": True,
            "runtime_influence": True,
            "score_contribution": 0.06,
            "isolation_status": "severe",
        },
    }
    annotation = honest_annotate(best, [])
    assert annotation["compound_confidence"] == "High"
    assert annotation["deconvolution_score_influence"] is True
    assert annotation["deconvolution_evidence_state"] == "score_applied"
    assert annotation["chimeric_isolation_warning"] is None


def test_same_run_clean_unknown_component_handles_library_gap(monkeypatch):
    _disable_optional_search_layers(monkeypatch)
    target_record = _search_record(
        0, TARGET_PRECURSOR_MZ, TARGET, "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    )
    records = [target_record]
    mz_index = {
        int(TARGET_PRECURSOR_MZ / search_module.MZ_BUCKET_DALTON): [0]
    }
    run_spectra = [
        {
            "spectrum_id": "clean-contaminant",
            "precursor_mz": COISOLATED_PRECURSOR_MZ,
            "precursor_charge": 1,
            "polarity": "POS",
            "peaks": [[mz, intensity * 1000] for mz, intensity in CONTAMINANT],
            "precursor_isolation_evidence": {"status": "clean", "purity": 0.95},
        },
        {
            "spectrum_id": "not-pure",
            "precursor_mz": COISOLATED_PRECURSOR_MZ + 0.1,
            "polarity": "POS",
            "peaks": [[mz, intensity * 1000] for mz, intensity in CONTAMINANT],
            "precursor_isolation_evidence": {"status": "mixed", "purity": 0.70},
        },
    ]
    run_pool = search_module.build_deconvolution_run_pool(run_spectra)
    assert len(run_pool[0]) == 1
    query = _mix(TARGET, CONTAMINANT, 1.0, 0.5)
    result = search_module.search(
        query,
        TARGET_PRECURSOR_MZ,
        records=records,
        mz_index=mz_index,
        query_polarity="POS",
        query_charge=1,
        query_isolation_evidence={
            "status": "severe",
            "purity": 0.42,
            "strongest_coisolated_peaks": [[COISOLATED_PRECURSOR_MZ, 800.0]],
        },
        enable_chimeric_deconvolution=True,
        deconvolution_run_pool=run_pool,
        cos_thresh=0.0,
        final_thresh=0.0,
    )
    assert result["chimeric_deconvolution"]["status"] == "accepted_score_applied"
    roles = {
        item["role"]
        for item in result["chimeric_deconvolution"]["component_candidates"]
    }
    assert "same_run_clean_spectrum" in roles
