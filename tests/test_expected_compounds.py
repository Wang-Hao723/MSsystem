from __future__ import annotations

import json
from pathlib import Path
import sqlite3

import numpy as np

from matching.expected_compounds import (
    COMPOUND_DISCOVERER_METHOD_REFERENCE,
    FORMULA_EXPECTED_DATASET_MAP,
    analyse_expected_compound_spectra,
    build_expected_ion_hypotheses,
    detect_sparse_xic_peaks,
    link_expected_features_to_openms,
    list_expected_datasets,
    load_expected_compound_catalog,
    resolve_expected_dataset_selection,
    run_expected_compounds_workflow,
    write_expected_compounds_artifact,
    _target_runtime_activation,
)


INCHIKEY = "WQZGKKKJIJFFOK-GASJEMHNSA-N"


def test_empty_target_selection_fails_closed_without_opening_database(tmp_path):
    missing = tmp_path / "does-not-exist.db"
    assert load_expected_compound_catalog(missing, dataset_sources=[]) == []
    assert not missing.exists()


def test_formula_filename_does_not_silently_bind_a_target_catalogue():
    available = [
        {"dataset_source": value}
        for value in FORMULA_EXPECTED_DATASET_MAP.values()
    ]
    resolved = resolve_expected_dataset_selection(
        "9a7c_F001_NEG_1.mzML", [], available
    )
    assert resolved["selected_dataset_sources"] == []
    assert resolved["selection_source"] == "none"
    assert resolved["selection_status"] == (
        "unresolved_explicit_selection_required"
    )
    assert resolved["mapped_dataset_source"] == (
        FORMULA_EXPECTED_DATASET_MAP["F001"]
    )
    assert resolved["filename_inference_applied"] is False


def test_filename_autoselect_is_available_only_when_explicitly_opted_in():
    available = [
        {"dataset_source": value}
        for value in FORMULA_EXPECTED_DATASET_MAP.values()
    ]
    resolved = resolve_expected_dataset_selection(
        "9a7c_F001_NEG_1.mzML",
        [],
        available,
        allow_filename_autoselect=True,
    )
    assert resolved["selected_dataset_sources"] == [
        FORMULA_EXPECTED_DATASET_MAP["F001"]
    ]
    assert resolved["selection_source"] == "filename_autoselect"
    assert resolved["selection_status"] == "resolved_from_filename"
    assert resolved["filename_inference_applied"] is True


def test_target_dataset_resolution_never_uses_loose_filename_substrings():
    available = list(FORMULA_EXPECTED_DATASET_MAP.values())
    resolved = resolve_expected_dataset_selection(
        "XF0012_NEG.mzML", [], available
    )
    assert resolved["selected_dataset_sources"] == []
    assert resolved["selection_status"] == (
        "unresolved_explicit_selection_required"
    )


def test_unknown_explicit_target_dataset_is_reported_for_rejection():
    resolved = resolve_expected_dataset_selection(
        "F001_NEG_1.mzML",
        ["not-a-real-catalogue"],
        list(FORMULA_EXPECTED_DATASET_MAP.values()),
    )
    assert resolved["selected_dataset_sources"] == []
    assert resolved["selection_status"] == (
        "invalid_unknown_requested_sources"
    )
    assert resolved["unknown_requested_sources"] == [
        "not-a-real-catalogue"
    ]


def test_unmapped_filename_runs_no_target_layer_instead_of_all_catalogues(
    tmp_path,
):
    path = tmp_path / "unknown_NEG.mzML"
    path.write_text("", encoding="utf-8")
    payload = run_expected_compounds_workflow(
        path,
        db_path=tmp_path / "missing.db",
        polarity="NEG",
        dataset_sources=[],
    )
    assert payload["status"] == "unavailable_no_selected_target_dataset"
    assert "never searches all catalogues" in payload["summary"]["warning"]


def test_expected_compounds_artifact_is_compact_and_round_trips(tmp_path):
    path = tmp_path / "expected.json"
    payload = {
        "summary": {"detected_feature_count": 2},
        "features": [
            {"feature_id": "F1", "candidate_names": ["化合物 A"]},
            {"feature_id": "F2", "candidate_names": ["Compound B"]},
        ],
    }

    write_expected_compounds_artifact(path, payload)

    serialized = path.read_text(encoding="utf-8")
    assert serialized.endswith("\n")
    assert "\n  \"" not in serialized
    assert json.loads(serialized) == payload


def test_compound_discoverer_reference_is_explicitly_non_equivalent():
    reference = COMPOUND_DISCOVERER_METHOD_REFERENCE

    assert reference["reference_software"].endswith("3.5")
    assert "not asserted to reproduce" in reference["relationship"]
    assert (
        reference["compound_discoverer_reference_defaults"][
            "expected_mass_tolerance_ppm"
        ]
        == 5.0
    )
    assert (
        reference["compound_discoverer_reference_defaults"][
            "minimum_total_isotope_peaks"
        ]
        == 2
    )


def _target(
    *,
    target_id: str,
    name: str,
    inchikey: str = INCHIKEY,
    diagnostic_ions=(89.024, 119.035),
):
    return {
        "catalog_row_id": int(target_id.lstrip("t") or 1),
        "target_id": target_id,
        "compound_name": name,
        "formula": "C6H12O6",
        "exact_mass": 180.06338810418,
        "precursor_mz": 181.070664,
        "adduct": "[M+H]+",
        "inchikey": inchikey,
        "dataset_source": "user_contrib:test",
        "diagnostic_ions": list(diagnostic_ions),
    }


def _synthetic_ms1(ion, *, one_scan_only=False):
    if one_scan_only:
        intensities = [0, 0, 6000, 0, 0, 0, 0]
    else:
        intensities = [0, 1200, 3000, 6000, 3000, 1200, 0]
    spectra = []
    for index, abundance in enumerate(intensities):
        mzs = [100.0]
        values = [10.0]
        if abundance:
            mzs.extend(
                [
                    ion["target_mz"],
                    ion["target_mz"] + 1.0033548,
                    ion["target_mz"] + 2.0067096,
                ]
            )
            values.extend(
                [abundance, abundance * 0.065, abundance * 0.012]
            )
        order = np.argsort(mzs)
        spectra.append(
            {
                "ms level": 1,
                "id": f"scan={index}",
                "scanList": {
                    "scan": [{"scan start time": index * 0.05}]
                },
                "m/z array": np.asarray(mzs, dtype=float)[order],
                "intensity array": np.asarray(values, dtype=float)[order],
            }
        )
    return spectra


def _protonated_ion(targets):
    hypotheses = build_expected_ion_hypotheses(targets, "POS")
    return next(
        item for item in hypotheses if item["adduct"] == "[M+H]+"
    )


def test_formula_ion_is_primary_and_indistinguishable_names_collapse():
    targets = [
        _target(target_id="t1", name="Glucose"),
        _target(target_id="t2", name="Dextrose", inchikey=""),
    ]
    ion = _protonated_ion(targets)

    assert abs(ion["target_mz"] - 181.070664570801) < 1e-6
    assert ion["formula_theoretical_primary"] is True
    assert ion["reported_precursor_priority"] is False
    assert ion["candidate_names"] == ["Dextrose", "Glucose"]
    assert ion["catalog_row_count"] == 2
    assert ion["isomer_ambiguous"] is True

    adducts = {
        item["adduct"]
        for item in build_expected_ion_hypotheses(targets, "POS")
    }
    assert {"[M+H]+", "[M+Na]+", "[M+NH4]+"} <= adducts


def test_duplicate_external_target_ids_do_not_erase_catalogue_rows():
    first = _target(target_id="t1", name="Glucose")
    second = _target(target_id="t2", name="Dextrose")
    second["target_id"] = "t1"
    second["catalog_row_id"] = 2

    ion = _protonated_ion([first, second])

    assert ion["catalog_row_count"] == 2
    assert {member["catalog_row_id"] for member in ion["members"]} == {1, 2}


def test_wrong_reported_precursor_is_suppressed_from_production_hypotheses():
    target = _target(target_id="t1", name="Glucose")
    target["precursor_mz"] = 181.0
    protonated = [
        item
        for item in build_expected_ion_hypotheses([target], "POS")
        if item["adduct"] == "[M+H]+"
    ]
    assert len(protonated) == 1
    theoretical = protonated[0]
    assert abs(theoretical["target_mz"] - 181.070664570801) < 1e-6
    assert theoretical["generation_sources"] == ["formula_exact_mass"]
    assert "reported_precursor_conflicts_with_theoretical_ion" in (
        theoretical["ion_quality_flags"]
    )

    # The bad source value remains available only in an explicit audit/recall
    # mode; production XIC extraction does not execute it.
    audited = [
        item
        for item in build_expected_ion_hypotheses(
            [target],
            "POS",
            include_reported_precursor_conflicts=True,
        )
        if item["adduct"] == "[M+H]+"
    ]
    assert len(audited) == 2
    alternative = next(
        item for item in audited
        if item["reported_precursor_conflict_candidate"]
    )
    assert alternative["target_mz"] == 181.0
    assert alternative["requires_orthogonal_evidence_for_upgrade"] is True
    assert "reported_precursor_conflict_candidate" in alternative[
        "ion_origins"
    ]
    assert "reported_precursor_retained_as_alternative_hypothesis" in (
        alternative["ion_quality_flags"]
    )


def test_reported_precursor_within_twenty_ppm_does_not_duplicate_formula_ion():
    target = _target(target_id="t1", name="Glucose")
    theoretical = 181.070664570801
    target["precursor_mz"] = theoretical * (1 + 10e-6)
    protonated = [
        item
        for item in build_expected_ion_hypotheses(
            [target], "POS", include_safe_common_adducts=False
        )
        if item["adduct"] == "[M+H]+"
    ]
    assert len(protonated) == 1
    assert abs(protonated[0]["target_mz"] - theoretical) < 1e-6
    assert protonated[0]["reported_precursor_conflict_candidate"] is False


def test_conflict_candidate_never_collapses_into_formula_at_wide_xic_ppm():
    target = _target(target_id="t1", name="Glucose")
    theoretical = 181.070664570801
    target["precursor_mz"] = theoretical * (1 + 30e-6)
    protonated = [
        item
        for item in build_expected_ion_hypotheses(
            [target],
            "POS",
            include_safe_common_adducts=False,
            include_reported_precursor_conflicts=True,
            ambiguity_ppm=50,
        )
        if item["adduct"] == "[M+H]+"
    ]
    assert len(protonated) == 2
    assert {
        item["reported_precursor_conflict_candidate"]
        for item in protonated
    } == {False, True}


def test_whole_ms1_xic_detects_target_without_any_ms2():
    ion = _protonated_ion([_target(target_id="t1", name="Glucose")])
    result = analyse_expected_compound_spectra(
        _synthetic_ms1(ion),
        [ion],
        ms2_rows=[],
    )

    assert result["summary"]["detected_feature_count"] == 1
    feature = result["features"][0]
    assert feature["scan_points"] == 5
    assert type(feature["xic_nonzero_points"]) is int
    assert feature["apex_rt_minutes"] == 0.15
    assert feature["apex_intensity"] == 6000
    assert feature["formula_isotope_evidence"]["status"] == "supporting"
    assert feature["evidence_level"] == "MS1_ISOTOPE_SUPPORTED"
    assert feature["counted_in_standard_spectral_library_matches"] is False
    assert result["summary"]["standard_spectral_library_match_increment"] == 0


def test_one_scan_noise_does_not_become_expected_compound_feature():
    ion = _protonated_ion([_target(target_id="t1", name="Glucose")])
    result = analyse_expected_compound_spectra(
        _synthetic_ms1(ion, one_scan_only=True),
        [ion],
    )
    assert result["summary"]["detected_feature_count"] == 0


def test_two_reported_diagnostic_ions_raise_probable_evidence():
    ion = _protonated_ion([_target(target_id="t1", name="Glucose")])
    result = analyse_expected_compound_spectra(
        _synthetic_ms1(ion),
        [ion],
        ms2_rows=[
            {
                "spectrum_id": "ms2-1",
                "precursor_mz": ion["target_mz"],
                "retention_time": 0.15,
                "retention_time_unit": "minutes",
                "ion_mode": "POS",
                "query_peaks": [[89.024, 100.0], [119.035, 50.0]],
                "match": None,
            }
        ],
    )
    feature = result["features"][0]
    assert feature["evidence_level"] == "MS1_MS2_DIAGNOSTIC_PROBABLE"
    assert feature["evidence_display_label"] == (
        "MS1 + reported diagnostic-ion screening support"
    )
    assert feature["standard_confirmed"] is False
    assert feature["exact_identity_confirmed"] is False
    assert feature["feature_fragment_calibration_used"] is False
    assert feature["diagnostic_rule_validation_status"] == (
        "legacy_fixed_threshold_not_v15_conditionally_calibrated"
    )
    assert len(feature["ms2_diagnostic_support"][0][
        "matched_diagnostic_ions"
    ]) == 2


def test_confirmed_requires_high_exact_identity_and_measured_library_tier():
    target = _target(
        target_id="t1",
        name="Glucose",
        diagnostic_ions=(),
    )
    ion = _protonated_ion([target])
    base_row = {
        "spectrum_id": "ms2-1",
        "precursor_mz": ion["target_mz"],
        "retention_time": 0.15,
        "retention_time_unit": "minutes",
        "ion_mode": "POS",
        "query_peaks": [
            [89.0241, 100.0],
            [119.0352, 50.0],
            [149.0453, 25.0],
        ],
    }
    predicted = analyse_expected_compound_spectra(
        _synthetic_ms1(ion),
        [ion],
        ms2_rows=[
            {
                **base_row,
                "match": {
                    "inchikey": INCHIKEY,
                    "source_tier": "predicted",
                    "confidence": "High",
                    "has_structure": True,
                    "isomer_ambiguous": False,
                },
            }
        ],
    )
    assert predicted["features"][0]["evidence_level"] == (
        "MS1_ISOTOPE_SUPPORTED"
    )

    measured = analyse_expected_compound_spectra(
        _synthetic_ms1(ion),
        [ion],
        ms2_rows=[
            {
                **base_row,
                "match": {
                    "inchikey": INCHIKEY,
                    "source_tier": "experimental",
                    "confidence": "High",
                    "has_structure": True,
                    "isomer_ambiguous": False,
                    "lib_spectrum_id": "MSBNK-test",
                    "dataset_source": "MassBank",
                    "cosine": 0.95,
                    "final_score": 0.9,
                    "lib_peaks": [
                        [89.0240, 100.0],
                        [119.0350, 50.0],
                        [149.0450, 25.0],
                    ],
                },
            }
        ],
    )
    feature = measured["features"][0]
    assert feature["evidence_level"] == "EXPERIMENTAL_LIBRARY_CONFIRMED"
    assert feature["evidence_display_label"] == (
        "MS/MS experimental-library candidate support"
    )
    assert feature["standard_confirmed"] is False
    assert feature["exact_identity_confirmed"] is False
    assert feature["experimental_library_support"][0][
        "library_source_tier"
    ] == "experimental"
    assert feature["experimental_library_support"][0][
        "high_accuracy_fragment_support"
    ]["matched_pair_count"] == 3


def test_broad_library_similarity_does_not_upgrade_without_high_accuracy_pairs():
    target = _target(
        target_id="t1", name="Glucose", diagnostic_ions=()
    )
    ion = _protonated_ion([target])
    result = analyse_expected_compound_spectra(
        _synthetic_ms1(ion),
        [ion],
        ms2_rows=[{
            "spectrum_id": "ms2-broad-only",
            "precursor_mz": ion["target_mz"],
            "retention_time": 0.15,
            "retention_time_unit": "minutes",
            "ion_mode": "POS",
            "query_peaks": [
                [89.034, 100.0],
                [119.045, 50.0],
                [149.055, 25.0],
            ],
            "match": {
                "inchikey": INCHIKEY,
                "source_tier": "experimental",
                "confidence": "High",
                "has_structure": True,
                "isomer_ambiguous": False,
                "lib_spectrum_id": "MSBNK-broad-only",
                "dataset_source": "MassBank",
                "cosine": 0.95,
                "final_score": 0.9,
                "lib_peaks": [
                    [89.024, 100.0],
                    [119.035, 50.0],
                    [149.045, 25.0],
                ],
            },
        }],
    )

    feature = result["features"][0]
    assert feature["evidence_level"] == "MS1_ISOTOPE_SUPPORTED"
    assert feature["experimental_library_support"] == []
    rejection = feature["experimental_library_rejections"][0]
    assert rejection["reason"] == (
        "insufficient_high_accuracy_fragment_support"
    )
    assert rejection["high_accuracy_fragment_support"][
        "matched_pair_count"
    ] == 0


def test_strong_measured_different_structure_blocks_diagnostic_promotion():
    ion = _protonated_ion([_target(target_id="t1", name="Glucose")])
    result = analyse_expected_compound_spectra(
        _synthetic_ms1(ion),
        [ion],
        ms2_rows=[{
            "spectrum_id": "ms2-conflict",
            "precursor_mz": ion["target_mz"],
            "retention_time": 0.15,
            "retention_time_unit": "minutes",
            "ion_mode": "POS",
            "query_peaks": [[89.024, 100.0], [119.035, 50.0]],
            "match": {
                "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                "compound_name": "Different compound",
                "source_tier": "experimental",
                "confidence": "Low",
                "final_score": 0.75,
                "cosine": 0.79,
                "has_structure": True,
                "isomer_ambiguous": False,
            },
        }],
    )
    feature = result["features"][0]
    assert feature["evidence_level"] == "MS1_ISOTOPE_SUPPORTED"
    assert feature["ms2_diagnostic_support"] == []
    assert feature["ms2_diagnostic_conflicts"][0]["reason"] == (
        "orthogonal_library_identity_conflict"
    )


def test_connectivity_match_with_undefined_stereo_does_not_veto_diagnostics():
    ion = _protonated_ion([_target(target_id="t1", name="Glucose")])
    result = analyse_expected_compound_spectra(
        _synthetic_ms1(ion),
        [ion],
        ms2_rows=[{
            "spectrum_id": "ms2-connectivity",
            "precursor_mz": ion["target_mz"],
            "retention_time": 0.15,
            "retention_time_unit": "minutes",
            "ion_mode": "POS",
            "query_peaks": [[89.024, 100.0], [119.035, 50.0]],
            "match": {
                "inchikey": "WQZGKKKJIJFFOK-UHFFFAOYSA-N",
                "compound_name": "Glucose (stereo unspecified)",
                "source_tier": "experimental",
                "confidence": "Low",
                "final_score": 0.75,
                "cosine": 0.79,
                "has_structure": True,
                "isomer_ambiguous": False,
            },
        }],
    )
    feature = result["features"][0]
    assert feature["evidence_level"] == "MS1_MS2_DIAGNOSTIC_PROBABLE"
    assert feature["ms2_diagnostic_conflicts"] == []


def test_possible_isotope_subpattern_does_not_upgrade_mass_candidate():
    ion = _protonated_ion([_target(target_id="t1", name="Glucose")])
    spectra = _synthetic_ms1(ion)
    spacing = 1.00335483507
    for spectrum in spectra:
        mzs = list(spectrum["m/z array"])
        intensities = list(spectrum["intensity array"])
        target_index = next(
            (
                index for index, mz in enumerate(mzs)
                if abs(float(mz) - float(ion["target_mz"])) < 1e-6
            ),
            None,
        )
        if target_index is None:
            continue
        mzs.append(float(ion["target_mz"]) - spacing)
        intensities.append(float(intensities[target_index]) * 0.5)
        order = np.argsort(mzs)
        spectrum["m/z array"] = np.asarray(mzs)[order]
        spectrum["intensity array"] = np.asarray(intensities)[order]

    feature = analyse_expected_compound_spectra(
        spectra, [ion], ms2_rows=[]
    )["features"][0]
    assert feature["isotope_envelope"]["pattern_integrity"] == (
        "possible_subpattern"
    )
    assert feature["formula_isotope_evidence"]["status"] == "supporting"
    assert feature["evidence_level"] == "MS1_EXACT_MASS_CANDIDATE"


def test_expected_compound_isotope_envelope_rejects_nearby_nonisotope_peak():
    ion = _protonated_ion([_target(target_id="t1", name="Glucose")])
    spectra = _synthetic_ms1(ion)
    spacing = 1.00335483507
    for spectrum in spectra:
        mzs = list(spectrum["m/z array"])
        intensities = list(spectrum["intensity array"])
        target_index = next(
            (
                index for index, mz in enumerate(mzs)
                if abs(float(mz) - float(ion["target_mz"])) < 1e-6
            ),
            None,
        )
        if target_index is None:
            continue
        target_intensity = float(intensities[target_index])
        retained = [
            (mz, intensity)
            for mz, intensity in zip(mzs, intensities)
            if abs(float(mz) - (float(ion["target_mz"]) + spacing))
            > 0.001
            and abs(
                float(mz) - (float(ion["target_mz"]) + 2 * spacing)
            )
            > 0.001
        ]
        mzs = [item[0] for item in retained]
        intensities = [item[1] for item in retained]
        mzs.append(float(ion["target_mz"]) + spacing + 0.006)
        intensities.append(target_intensity * 0.5)
        order = np.argsort(mzs)
        spectrum["m/z array"] = np.asarray(mzs)[order]
        spectrum["intensity array"] = np.asarray(intensities)[order]

    feature = analyse_expected_compound_spectra(
        spectra, [ion], ms2_rows=[]
    )["features"][0]

    assert feature["isotope_envelope"]["mz_tolerance_ppm"] == 5.0
    assert feature["isotope_envelope"]["n_isotopes_observed"] == 0
    assert feature["formula_isotope_evidence"] is None
    assert feature["evidence_level"] == "MS1_EXACT_MASS_CANDIDATE"


def test_peak_gap_and_fwhm_quality_gates_are_enforced():
    ion = _protonated_ion([_target(target_id="t1", name="Glucose")])
    scan_rts = [index * 0.5 for index in range(7)]
    scan_ids = [f"scan={index}" for index in range(7)]
    hits = {
        index: {
            "intensity": intensity,
            "observed_mz": ion["target_mz"],
            "isotope_envelope": None,
        }
        for index, intensity in zip(
            [1, 2, 3, 4, 5],
            [1200, 3000, 6000, 3000, 1200],
        )
    }
    assert detect_sparse_xic_peaks(
        ion,
        scan_rts,
        scan_ids,
        hits,
        maximum_fwhm_seconds=30,
    ) == []

    gapped_hits = {
        index: {
            "intensity": intensity,
            "observed_mz": ion["target_mz"],
            "isotope_envelope": None,
        }
        for index, intensity in zip(
            [0, 1, 4, 5, 6],
            [1200, 3000, 6000, 3000, 1200],
        )
    }
    assert detect_sparse_xic_peaks(
        ion,
        [index * 0.01 for index in range(7)],
        scan_ids,
        gapped_hits,
        maximum_missing_scans=3,
        maximum_gap_ratio=0.20,
    ) == []


def test_openms_cross_check_is_chromatographic_not_identity_evidence():
    features = [{
        "feature_id": "ECI_1_P001",
        "target_mz": 181.070664,
        "observed_apex_mz": 181.070700,
        "apex_rt_minutes": 1.25,
        "boundary_start_rt_minutes": 1.20,
        "boundary_end_rt_minutes": 1.30,
    }]
    summary = link_expected_features_to_openms(
        features,
        [{
            "feature_id": "LF0001",
            "mz": 181.07069,
            "rt_minutes": 1.251,
            "rt_start_minutes": 1.22,
            "rt_end_minutes": 1.29,
            "intensity": 9000,
        }],
    )
    assert features[0]["openms_feature_supported"] is True
    assert features[0]["openms_feature_id"] == "LF0001"
    assert summary["linked_targeted_features"] == 1
    assert summary["identity_confidence_increment"] == 0


def test_catalog_loader_filters_dataset_and_preserves_diagnostic_ions(tmp_path):
    database = tmp_path / "targets.db"
    connection = sqlite3.connect(database)
    connection.executescript(
        """
        CREATE TABLE TargetCompounds (
            id INTEGER PRIMARY KEY,
            target_id TEXT,
            dataset_source TEXT,
            compound_name TEXT,
            molecular_formula TEXT,
            smiles TEXT,
            inchi TEXT,
            inchikey TEXT,
            reported_inchikey TEXT,
            structure_identity_verified INTEGER,
            structure_identity_source TEXT,
            exact_mass REAL,
            precursor_mz REAL,
            adduct TEXT,
            charge INTEGER,
            plant_source TEXT,
            source_file TEXT,
            source_content_sha256 TEXT,
            dataset_version TEXT,
            source_row_index INTEGER,
            provenance TEXT,
            doi TEXT,
            dataset_accession TEXT,
            spectrum_accession TEXT,
            quality_flags_json TEXT,
            source_tier TEXT
        );
        CREATE TABLE DiagnosticIons (
            id INTEGER PRIMARY KEY,
            target_id INTEGER,
            fragment_mz REAL,
            ordinal INTEGER
        );
        """
    )
    for identifier, dataset in ((1, "dataset-a"), (2, "dataset-b")):
        connection.execute(
            """
            INSERT INTO TargetCompounds (
                id, target_id, dataset_source, compound_name,
                molecular_formula, exact_mass, precursor_mz, adduct,
                smiles, structure_identity_verified,
                source_tier, quality_flags_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                identifier,
                f"t{identifier}",
                dataset,
                "Glucose",
                "C6H12O6",
                180.063388,
                181.070664,
                "[M+H]+",
                "C1(" if identifier == 2 else "",
                1 if identifier == 2 else 0,
                "literature_presence_only",
                '["presence_only"]',
            ),
        )
        connection.executemany(
            "INSERT INTO DiagnosticIons (target_id, fragment_mz, ordinal) "
            "VALUES (?, ?, ?)",
            [(identifier, 89.024, 1), (identifier, 119.035, 2)],
        )
    connection.execute(
        """
        INSERT INTO TargetCompounds (
            id, target_id, dataset_source, compound_name,
            molecular_formula, exact_mass, precursor_mz, adduct,
            smiles, structure_identity_verified,
            source_tier, quality_flags_json
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            3,
            "t3",
            "dataset-b",
            "Corrupted formula row",
            "383",
            None,
            383.0,
            "[M-H]-",
            "",
            0,
            "literature_presence_only",
            '["identity_status:rejected_identity_formula_mismatch"]',
        ),
    )
    connection.commit()
    connection.close()

    catalog = load_expected_compound_catalog(
        database, dataset_sources=["dataset-b"]
    )
    assert len(catalog) == 2
    valid = next(item for item in catalog if item["target_id"] == "t2")
    rejected = next(item for item in catalog if item["target_id"] == "t3")
    assert valid["dataset_source"] == "dataset-b"
    assert valid["diagnostic_ions"] == [89.024, 119.035]
    assert valid["smiles"] == ""
    assert valid["reported_smiles"] == "C1("
    assert valid["structure_validation_status"] == (
        "invalid_reported_structure_quarantined"
    )
    assert valid["structure_identity_verified"] is False
    assert valid["reported_structure_identity_verified"] is True
    assert valid["quality_flags"] == [
        "invalid_reported_structure_quarantined", "presence_only"
    ]
    assert rejected["runtime_active"] is False
    assert rejected["reported_formula"] == "383"
    assert set(rejected["runtime_quarantine_reasons"]) == {
        "identity_status:rejected_identity_formula_mismatch",
        "invalid_nonempty_molecular_formula",
    }
    assert all(
        "Corrupted formula row" not in item["candidate_names"]
        for item in build_expected_ion_hypotheses(catalog, "NEG")
    )
    datasets = list_expected_datasets(database)
    assert [item["dataset_source"] for item in datasets] == [
        "dataset-a", "dataset-b"
    ]


def test_malformed_reported_adduct_is_not_an_executable_ion_hypothesis():
    target = _target(target_id="t1", name="Glucose")
    target["adduct"] = "577.1354[M-H]-"
    target["reported_adduct"] = "577.1354[M-H]-"
    hypotheses = build_expected_ion_hypotheses([target], "POS")

    assert hypotheses
    assert {item["adduct"] for item in hypotheses} == {
        "[M+H]+", "[M+Na]+", "[M+NH4]+"
    }
    assert all(
        "reported_precursor_without_theoretical_ion"
        not in item["ion_origins"]
        for item in hypotheses
    )


def test_bare_reported_precursor_without_ion_identity_is_quarantined():
    target = {
        "target_id": "bare-mz",
        "reported_formula": "",
        "formula": "",
        "exact_mass": None,
        "reported_precursor_mz": 383.1234,
        "precursor_mz": 383.1234,
        "adduct": "",
        "quality_flags": [],
    }

    activation = _target_runtime_activation(target)

    assert activation["runtime_active"] is False
    assert activation["runtime_quarantine_reasons"] == [
        "reported_precursor_without_executable_adduct"
    ]


def test_reported_precursor_with_explicit_ion_identity_remains_executable():
    target = {
        "target_id": "explicit-mz",
        "compound_name": "Reported ion only",
        "reported_formula": "",
        "formula": "",
        "exact_mass": None,
        "reported_precursor_mz": 383.1234,
        "precursor_mz": 383.1234,
        "adduct": "[M-H]-",
        "reported_adduct": "[M-H]-",
        "quality_flags": [],
    }
    target.update(_target_runtime_activation(target))

    hypotheses = build_expected_ion_hypotheses([target], "NEG")

    assert target["runtime_active"] is True
    assert len(hypotheses) == 1
    assert hypotheses[0]["target_mz"] == 383.1234
    assert hypotheses[0]["ion_origins"] == [
        "reported_precursor_without_theoretical_ion"
    ]


def test_mgf_is_explicitly_unavailable_for_expected_compounds(tmp_path):
    input_path = tmp_path / "sample.mgf"
    input_path.write_text("BEGIN IONS\nEND IONS\n", encoding="utf-8")
    payload = run_expected_compounds_workflow(
        input_path,
        db_path=tmp_path / "missing.db",
        polarity="POS",
    )
    assert payload["status"] == "unavailable_requires_mzml_with_ms1"
    assert "MGF/HTML" in payload["summary"]["warning"]
    assert payload["features"] == []


def test_expected_compounds_endpoint_reads_separate_artifact(
    tmp_path, monkeypatch
):
    import web.server as server

    task_id = "expected-compounds-test"
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    server._tasks[task_id] = {
        "id": task_id,
        "status": "completed",
        "expected_compounds_summary": {"status": "ok"},
    }
    payload = {
        "task_id": task_id,
        "status": "ok",
        "summary": {
            "detected_feature_count": 1,
            "counted_in_standard_spectral_library_matches": False,
        },
        "features": [{"feature_id": "ECI_1_P001"}],
    }
    try:
        write_expected_compounds_artifact(
            tmp_path / f"{task_id}_expected_compounds.json",
            payload,
        )
        returned = server.get_task_expected_compounds(
            task_id,
            supported_only=False,
            evidence_level=None,
            limit=500,
            offset=0,
        )
        assert returned["features"][0]["feature_id"] == "ECI_1_P001"
        assert returned["summary"][
            "counted_in_standard_spectral_library_matches"
        ] is False
    finally:
        server._tasks.pop(task_id, None)


def test_server_layer_persists_explicit_mgf_unavailable_status(
    tmp_path, monkeypatch
):
    import web.server as server

    input_path = tmp_path / "sample.mgf"
    input_path.write_text("BEGIN IONS\nEND IONS\n", encoding="utf-8")
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    payload, summary = server._run_expected_compounds_layer(
        "mgf-task",
        file_path=str(input_path),
        detected_ion_mode="NEG",
        params={"enable_expected_compounds": True},
        serialized_ms2_rows=[],
    )
    assert payload["status"] == "unavailable_requires_mzml_with_ms1"
    assert summary["available"] is False
    assert summary["counted_in_standard_spectral_library_matches"] is False
    assert (tmp_path / "mgf-task_expected_compounds.json").is_file()


def test_analysis_params_use_cd_expected_compound_defaults():
    from web.server import AnalysisParams

    params = AnalysisParams()
    assert params.enable_expected_compounds is False
    assert params.expected_mass_tolerance_ppm == 5.0
    assert params.expected_min_peak_intensity == 1000.0
    assert params.expected_min_snr == 1.5
    assert params.expected_min_scan_points == 5
    assert params.expected_max_gap_ratio == 0.35
    assert params.expected_max_fwhm_seconds == 60.0
    assert params.expected_min_relative_valley_depth == 0.10
    assert params.expected_openms_rt_tolerance_minutes == 0.05
    assert params.expected_background_sample_blank_ratio == 5.0
    assert params.expected_generate_common_adducts is True


def test_results_ui_groups_expected_features_and_collapses_dense_details():
    results_js = (
        Path(__file__).resolve().parents[1]
        / "web"
        / "static"
        / "js"
        / "results.js"
    ).read_text(encoding="utf-8")

    assert "function _groupExpectedFeatures(features)" in results_js
    assert "const groups = _groupExpectedFeatures(features);" in results_js
    assert "候选组（已去重）" in results_js
    assert "展开 ${group.rawFeatureCount} 个原始 XIC 峰" in results_js
    assert "主表每个候选组只显示一次" in results_js

    mass_only = results_js.index("仅精确质量的去重候选组")
    mass_only_container = results_js[mass_only - 220 : mass_only]
    assert "<details " in mass_only_container
    assert "<details open" not in mass_only_container

    warning = results_js.index("MS1 模块适用性提示")
    warning_container = results_js[warning - 260 : warning]
    assert "<details " in warning_container
    assert "<details open" not in warning_container
