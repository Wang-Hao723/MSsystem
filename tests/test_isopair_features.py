"""Unit tests for the locked IsoPair-MS feature contract."""

import numpy as np

from matching.isopair_features import (
    FEATURE_NAMES,
    spectral_pair_feature_dict,
    spectral_pair_features,
)


def _features(left, right, left_pmz=200.0, right_pmz=200.0):
    return spectral_pair_feature_dict(
        left,
        right,
        query_precursor_mz=left_pmz,
        reference_precursor_mz=right_pmz,
    )


def test_feature_schema_is_locked_and_contains_no_precursor_mass_shortcut():
    assert len(FEATURE_NAMES) == 62
    assert len(set(FEATURE_NAMES)) == len(FEATURE_NAMES)
    assert all("precursor_delta" not in name for name in FEATURE_NAMES)
    assert all("ppm" not in name for name in FEATURE_NAMES)


def test_identical_spectra_have_strong_strict_and_intensity_agreement():
    peaks = [(55.01, 100.0), (91.03, 45.0), (123.04, 20.0)]
    features = _features(peaks, peaks)

    assert features["product_cosine_0p01"] > 0.999
    assert features["neutral_loss_cosine_0p01"] > 0.999
    assert features["product_matched_intensity_corr_0p02"] > 0.999
    assert features["product_log_intensity_mae_0p02"] < 1e-6


def test_strict_and_broad_views_separate_small_mass_shifts():
    query = [(50.000, 100.0), (80.000, 60.0), (110.000, 30.0)]
    shifted = [(50.035, 100.0), (80.035, 60.0), (110.035, 30.0)]
    features = _features(query, shifted)

    assert features["product_cosine_0p02"] == 0.0
    assert features["product_cosine_0p05"] > 0.999
    assert features["strict_product_fraction"] < 0.01


def test_relative_intensity_mismatch_is_explicitly_measured():
    query = [(50.0, 100.0), (80.0, 50.0), (110.0, 10.0)]
    reordered = [(50.0, 10.0), (80.0, 50.0), (110.0, 100.0)]
    features = _features(query, reordered)

    assert features["product_matched_intensity_corr_0p02"] < 0.0
    assert features["product_log_intensity_mae_0p02"] > 1.0


def test_neutral_loss_view_can_agree_when_product_ions_do_not():
    query = [(100.0, 100.0), (130.0, 60.0), (160.0, 30.0)]
    reference = [(120.0, 100.0), (150.0, 60.0), (180.0, 30.0)]
    features = _features(query, reference, 200.0, 220.0)

    assert features["neutral_loss_cosine_0p02"] > 0.999
    assert features["product_cosine_0p02"] == 0.0


def test_vector_output_is_finite_and_matches_schema_order():
    vector = spectral_pair_features(
        [(50.0, 10.0), (70.0, 5.0), (90.0, 1.0)],
        [(50.0, 9.0), (70.01, 4.0), (100.0, 2.0)],
        query_precursor_mz=150.0,
        reference_precursor_mz=150.0,
    )

    assert vector.shape == (len(FEATURE_NAMES),)
    assert vector.dtype == np.float32
    assert np.isfinite(vector).all()
