from matching.precursor_features import annotate_precursor_features


def _row(scan, mz, rt, energy, polarity="POS", confidence=None):
    match = None
    if confidence:
        match = {
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "compound_name": "Example",
            "smiles": "CCO",
            "confidence": confidence,
            "confidence_class": confidence,
        }
    return {
        "spectrum_id": scan,
        "precursor_mz": mz,
        "retention_time": rt,
        "ion_mode": polarity,
        "collision_energy": energy,
        "match": match,
    }


def test_repeated_stepped_energy_scans_form_one_precursor_feature():
    rows = [
        _row("1", 500.2000, 4.000, 20.0, confidence="Medium"),
        _row("2", 500.2008, 4.010, 40.0, confidence="High"),
        _row("3", 500.1997, 4.020, 60.0, confidence="Medium"),
        _row("4", 500.2000, 4.200, 20.0, confidence="High"),
    ]

    summary = annotate_precursor_features(rows)

    assert rows[0]["precursor_feature_id"] == rows[1]["precursor_feature_id"]
    assert rows[1]["precursor_feature_id"] == rows[2]["precursor_feature_id"]
    assert rows[3]["precursor_feature_id"] != rows[0]["precursor_feature_id"]
    assert rows[0]["precursor_feature_spectra"] == 3
    assert rows[0]["precursor_feature_multi_energy"] is True
    assert rows[0]["precursor_feature_collision_energies"] == [20.0, 40.0, 60.0]
    assert rows[0]["precursor_feature_identity_consistent"] is True
    assert summary["total_spectra"] == 4
    assert summary["precursor_features"] == 2
    assert summary["matched_precursor_features"] == 2
    assert summary["high_precursor_features"] == 2


def test_polarity_and_missing_retention_time_are_never_merged():
    rows = [
        _row("1", 300.0, 1.0, 20.0, polarity="POS"),
        _row("2", 300.0, 1.01, 20.0, polarity="NEG"),
        _row("3", 300.0, None, 20.0, polarity="POS"),
        _row("4", 300.0, None, 40.0, polarity="POS"),
    ]

    summary = annotate_precursor_features(rows)

    assert len({row["precursor_feature_id"] for row in rows}) == 4
    assert summary["precursor_features"] == 4
    assert summary["approximate"] is True


def test_openms_feature_identity_overrides_fallback_event_for_summary_counts():
    rows = [
        _row("1", 500.2000, 4.000, 20.0, confidence="Medium"),
        _row("2", 500.2000, 4.200, 40.0, confidence="High"),
        _row("3", 300.1000, 2.000, 20.0, confidence="Medium"),
    ]
    # The first two scans are too far apart for the MS2-only fallback cluster,
    # but OpenMS attached both to the same measured chromatographic feature.
    rows[0]["lcms_feature_id"] = "LF000001"
    rows[1]["lcms_feature_id"] = "LF000001"

    summary = annotate_precursor_features(rows)

    assert rows[0]["precursor_feature_id"] != rows[1]["precursor_feature_id"]
    assert summary["precursor_features"] == 2
    assert summary["lcms_features_with_ms2"] == 1
    assert summary["lcms_repeated_ms2_features"] == 1
    assert summary["fallback_precursor_features"] == 1
    assert summary["matched_precursor_features"] == 2
    assert summary["high_precursor_features"] == 1
    assert summary["multi_energy_precursor_features"] == 1
