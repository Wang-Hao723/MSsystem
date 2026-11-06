from utils.spectrum import process_peaks, resolve_peak_processing_config


def test_dominant_excluded_precursor_does_not_erase_measured_products():
    config = resolve_peak_processing_config(
        "library_match",
        relative_intensity_threshold=0.05,
        precursor_exclusion_da=0.5,
        top_n=15,
    )
    peaks, audit = process_peaks(
        [
            [50.0, 100.0],
            [75.0, 80.0],
            [125.0, 60.0],
            [199.99, 10_000.0],
        ],
        200.0,
        config=config,
    )

    assert peaks == [[50.0, 100.0], [75.0, 80.0], [125.0, 60.0]]
    assert audit["schema_version"] == 2
    assert audit["precursor_dominance_rescue_triggered"] is True


def test_rescue_does_not_invent_minimum_fragment_evidence():
    config = resolve_peak_processing_config(
        "library_match",
        relative_intensity_threshold=0.05,
        precursor_exclusion_da=0.5,
        top_n=15,
    )
    peaks, audit = process_peaks(
        [[50.0, 100.0], [75.0, 80.0], [199.99, 10_000.0]],
        200.0,
        config=config,
    )

    assert peaks == []
    assert audit["precursor_dominance_rescue_triggered"] is False
