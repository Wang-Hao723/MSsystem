"""Regression tests for HRMS precursor screening and cross-adduct matching."""

from __future__ import annotations

import importlib

import pytest

from config.settings import MZ_TOLERANCE_DA
from csv_excel_to_mgf import (
    _ion_mode_from_adduct,
    normalize_adduct,
    normalize_charge,
    write_mgf,
)
from utils.mgf_io import parse_mgf_file
from matching.scores import mass_score, reverse_match_ratio


search_module = importlib.import_module("matching.search")


def _record(index, precursor_mz=500.0, *, polarity="POS",
            source_tier="experimental", adduct="[M+H]+",
            neutral_mass=None):
    peaks = [[100.0, 1.0], [150.0, 0.8], [200.0, 0.6]]
    return {
        "precursor_mz": precursor_mz,
        "neutral_mass": neutral_mass,
        "neutral_mass_source": "compound_exact_mass" if neutral_mass else "",
        "peaks": peaks,
        "clean_peaks": peaks,
        "spectrum_id": f"spec-{index}",
        "compound_name": f"compound-{index}",
        "source_compound_name": "",
        "name_resolution_status": "",
        "smiles": "",
        "has_structure": False,
        "name_only_lead": False,
        "inchikey": "",
        "formula": "",
        "adduct": adduct,
        "retention_time": None,
        "classification": "other",
        "dataset_source": "test",
        "source_tier": source_tier,
        "intensity_origin": "measured",
        "quality_flags": [],
        "pathway_json": "",
        "polarity": polarity,
    }


def test_mgf_preserves_specific_adduct_and_uses_it_for_polarity(tmp_path):
    mgf = tmp_path / "adducts.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=sodium\nPEPMASS=322.989218\n"
        "ADDUCT=[M+Na]+\n100 100\n150 50\nEND IONS\n"
        "BEGIN IONS\nTITLE=deprotonated\nPEPMASS=298.992724\n"
        "IONTYPE=[M-H]-\n100 100\n150 50\nEND IONS\n",
        encoding="utf-8",
    )

    records = parse_mgf_file(mgf)

    assert [record["adduct"] for record in records] == ["[M+Na]+", "[M-H]-"]
    assert [record["polarity"] for record in records] == ["POS", "NEG"]
    assert all(record["adduct_reliable"] for record in records)
    assert all(record["polarity_source"] == "ADDUCT" for record in records)


def test_mgf_charge_adduct_conflict_disables_neutral_matching(tmp_path):
    mgf = tmp_path / "conflict.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=conflict\nPEPMASS=322.989218\n"
        "CHARGE=1-\nADDUCT=[M+Na]+\n100 100\n150 50\nEND IONS\n",
        encoding="utf-8",
    )

    record = parse_mgf_file(mgf)[0]

    assert record["adduct"] == "[M+Na]+"
    assert record["polarity_conflict"] is True
    assert record["polarity"] is None
    assert record["adduct_reliable"] is False


def test_mgf_charge_magnitude_must_agree_with_adduct(tmp_path):
    mgf = tmp_path / "charge_magnitude_conflict.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=wrong-charge\nPEPMASS=322.989218\n"
        "CHARGE=2+\nADDUCT=[M+Na]+\n100 100\n150 50\nEND IONS\n",
        encoding="utf-8",
    )

    record = parse_mgf_file(mgf)[0]

    assert record["polarity"] == "POS"
    assert record["polarity_conflict"] is False
    assert record["adduct_reliable"] is False


def test_mgf_exposes_charge_and_collision_energy_as_expert_evidence(tmp_path):
    mgf = tmp_path / "expert_metadata.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=metadata\nPEPMASS=322.989218 4567\n"
        "CHARGE=2+\nCOLLISIONENERGY=35 eV\n"
        "100 100\n150 50\nEND IONS\n",
        encoding="utf-8",
    )

    record = parse_mgf_file(mgf)[0]

    assert record["precursor_charge"] == 2
    assert record["collision_energy"] == pytest.approx(35.0)
    assert record["collision_energy_unit"] == "eV"
    assert record["isotope_evidence_status"] == "ms1_not_available_in_mgf"
    assert record["precursor_isotope_evidence"] is None


def test_mgf_preserves_normalized_collision_energy_unit(tmp_path):
    mgf = tmp_path / "nce_metadata.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=nce\nPEPMASS=300\nCHARGE=1+\n"
        "NCE=28%\n100 100\n150 50\nEND IONS\n",
        encoding="utf-8",
    )

    record = parse_mgf_file(mgf)[0]

    assert record["collision_energy"] == pytest.approx(28.0)
    assert record["collision_energy_unit"] == "NCE"


def test_only_specific_polarity_consistent_adduct_creates_neutral_mass():
    neutral, polarity = search_module._trusted_query_neutral_mass(
        322.989218, "[M+Na]+", "POS"
    )
    assert neutral == pytest.approx(300.0, abs=1e-9)
    assert polarity == "POS"

    assert search_module._trusted_query_neutral_mass(
        322.989218, "Positive", "POS"
    ) == (None, None)


def test_theoretical_ion_mz_does_not_guess_from_polarity_only_label():
    from utils.adduct import ion_mz

    assert ion_mz(300.0, "Positive") is None
    assert ion_mz(300.0, "[M+H]+") == pytest.approx(301.007276)
    assert search_module._trusted_query_neutral_mass(
        322.989218, "[M+Na]+", "NEG"
    ) == (None, None)


def test_single_precursor_adduct_detection_is_explicitly_ambiguous():
    from utils.adduct import detect_adduct

    result = detect_adduct(301.007276, mode="positive")

    assert result["status"] == "ambiguous"
    assert result["adduct"] is None
    assert result["confidence"] == "unknown"
    assert {item["adduct"] for item in result["hypotheses"]} >= {
        "[M+H]+",
        "[M+Na]+",
    }


def test_independent_neutral_mass_or_coobserved_pair_can_assign_adducts():
    from utils.adduct import batch_detect_adducts, detect_adduct

    assigned = detect_adduct(
        322.989218,
        mode="positive",
        neutral_mass_reference=300.0,
    )
    assert assigned["status"] == "assigned"
    assert assigned["adduct"] == "[M+Na]+"
    assert assigned["confidence"] == "high"

    paired = batch_detect_adducts(
        [301.007276, 322.989218, 450.0],
        mode="positive",
    )
    assert paired[0]["adduct"] == "[M+H]+"
    assert paired[1]["adduct"] == "[M+Na]+"
    assert paired[0]["group_id"] == paired[1]["group_id"]
    assert paired[2]["status"] == "ambiguous"
    assert paired[2]["adduct"] is None


def test_cross_adduct_ppm_is_calculated_on_neutral_mass():
    record = _record(
        0, precursor_mz=301.007276, adduct="[M+H]+", neutral_mass=300.001
    )

    _, ppm, basis, query_value, library_value = search_module._candidate_mass_metrics(
        322.989218,
        record,
        query_neutral_mass=300.0,
        cross_adduct=True,
        ppm_tolerance=10.0,
    )

    assert basis == "neutral_mass"
    assert query_value == pytest.approx(300.0)
    assert library_value == pytest.approx(300.001)
    assert ppm == pytest.approx(3.333333, rel=1e-5)
    assert ppm < 10.0


def test_mass_score_uses_the_per_request_absolute_tolerance():
    record = _record(0, precursor_mz=500.1)

    narrow = search_module._candidate_mass_metrics(
        500.0, record, ppm_tolerance=10.0, da_tolerance=0.02
    )[0]
    broad = search_module._candidate_mass_metrics(
        500.0, record, ppm_tolerance=10.0, da_tolerance=0.5
    )[0]

    assert broad > narrow


def test_mass_score_is_monotonic_across_ppm_to_da_boundary():
    query_mz = 500.0
    errors_ppm = [0, 1, 5, 9.9, 10.0, 10.1, 20, 30, 39.9]
    scores = [
        mass_score(
            query_mz,
            query_mz + query_mz * error_ppm / 1e6,
            ppm_tol=10.0,
            da_tol=0.02,
        )
        for error_ppm in errors_ppm
    ]

    assert all(left >= right for left, right in zip(scores, scores[1:]))
    assert scores[4] >= scores[5]


def test_reverse_match_accepts_modified_cosine_precursor_shift_pairs():
    query = [[50.0, 1.0], [70.0, 0.5]]
    library = [[40.0, 1.0], [60.0, 0.5]]

    direct_only = reverse_match_ratio(query, library, mz_tol=0.1)
    shifted = reverse_match_ratio(
        query,
        library,
        mz_tol=0.1,
        query_pmz=110.0,
        lib_pmz=100.0,
    )

    assert direct_only == 0.0
    assert shifted == pytest.approx(1.0)


def test_direct_match_prefers_theoretical_ion_mz_over_rounded_library_value():
    record = _record(0, precursor_mz=301.01, adduct="[M+H]+", neutral_mass=300.0)
    record["neutral_mass_source"] = "compound_exact_mass"

    _, ppm, basis, query_value, library_value = search_module._candidate_mass_metrics(
        301.007276, record, ppm_tolerance=10.0
    )

    assert basis == "theoretical_ion_mz"
    assert query_value == pytest.approx(301.007276)
    assert library_value == pytest.approx(301.007276)
    assert ppm == pytest.approx(0.0, abs=1e-9)


def test_direct_match_exposes_incompatible_exact_mass_and_adduct_metadata():
    record = _record(0, precursor_mz=301.007, adduct="[M+Na]+", neutral_mass=300.0)
    record["neutral_mass_source"] = "compound_exact_mass"

    _, ppm, basis, _, library_value = search_module._candidate_mass_metrics(
        301.007, record, ppm_tolerance=10.0
    )

    assert basis == "precursor_mz_metadata_conflict"
    assert library_value == pytest.approx(301.007)
    assert ppm == pytest.approx(0.0, abs=1e-9)


def test_gross_structure_adduct_precursor_conflict_is_a_hard_search_gate():
    audit = search_module._structure_adduct_precursor_audit({
        "precursor_mz": 134.898,
        "reported_precursor_mz": 134.898,
        "adduct": "[M-H]-",
        "smiles": (
            "CCC(/C=C/C(C)C1CCC2(C)C3CCC45CC(O)CCC4(C5)"
            "C3=CCC12C)C(C)C"
        ),
        "inchikey": "QECYIHDHCKWWPU-CMDGGOBGSA-N",
        "formula": "C30H48O",
        "has_structure": True,
    })

    assert audit["status"] == (
        "gross_structure_adduct_precursor_conflict"
    )
    assert audit["hard_reject"] is True
    assert audit["theoretical_precursor_mz"] == pytest.approx(
        423.363240, abs=1e-6
    )


def test_valid_structure_adduct_precursor_metadata_remains_searchable():
    audit = search_module._structure_adduct_precursor_audit({
        "precursor_mz": 47.049141,
        "reported_precursor_mz": 47.049141,
        "adduct": "[M+H]+",
        "smiles": "CCO",
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "formula": "C2H6O",
        "has_structure": True,
    })

    assert audit["status"] == "consistent"
    assert audit["hard_reject"] is False


def test_missing_library_adduct_requires_a_structure_mass_explanation():
    valid = search_module._structure_adduct_precursor_audit({
        "precursor_mz": 301.03540,
        "reported_precursor_mz": 301.03540,
        "adduct": "",
        "polarity": "NEG",
        "smiles": "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "inchikey": "REFJWTPEDVJJIY-UHFFFAOYSA-N",
        "formula": "C15H10O7",
        "has_structure": True,
    }, query_polarity="NEG")

    assert valid["hard_reject"] is False
    assert valid["status"] == (
        "consistent_via_structure_formula_adduct_hypothesis"
    )
    assert valid["inferred_adduct"] == "[M-H]-"

    impossible = search_module._structure_adduct_precursor_audit({
        "precursor_mz": 253.14189,
        "reported_precursor_mz": 253.14189,
        "adduct": "",
        "polarity": "NEG",
        "smiles": "CCCCC(=O)OC",
        "inchikey": "XEKOWRVHYACXOJ-UHFFFAOYSA-N",
        "formula": "C6H12O2",
        "has_structure": True,
    }, query_polarity="NEG")

    assert impossible["hard_reject"] is True
    assert impossible["status"] == (
        "gross_structure_formula_precursor_conflict_without_adduct"
    )


def test_hrms_query_uses_ppm_as_hard_boundary_not_point_zero_two_da(
    monkeypatch,
):
    # 0.0165 Da is inside the legacy 0.02-Da retrieval floor but is 33 ppm at
    # m/z 500. A QE/Orbitrap query configured for 10 ppm must reject it before
    # spectral similarity can promote it.
    record = _record(0, precursor_mz=500.0165)
    records = [record]
    mz_index = {
        int(record["precursor_mz"] / search_module.MZ_BUCKET_DALTON): [0]
    }
    monkeypatch.setattr(search_module, "EXPERT_SYSTEM_ENABLED", False)
    monkeypatch.setattr(search_module, "USE_EXPERT_V3", False)
    monkeypatch.setattr(search_module, "USE_CONSISTENCY_FILTER", False)
    monkeypatch.setattr(search_module, "USE_ADAPTIVE_PEAK_SEARCH", False)
    monkeypatch.setattr(search_module, "USE_ISOMER_FRAGMENT_RERANK", False)
    monkeypatch.setattr(search_module, "MODEL_PREDICTOR_AVAILABLE", False)
    monkeypatch.setattr(search_module, "_ensure_clf", lambda: False)
    query = [[100.0, 1.0], [150.0, 0.8], [200.0, 0.6]]
    kwargs = dict(
        records=records,
        mz_index=mz_index,
        query_polarity="POS",
        mz_da=0.02,
        mz_ppm=10.0,
        cos_thresh=0.0,
        final_thresh=0.0,
        strict_evidence_gate=False,
    )

    assert search_module.search(
        query,
        500.0,
        query_instrument_type="Orbitrap",
        **kwargs,
    ) is None
    assert search_module.search(
        query,
        500.0,
        query_instrument_type="ion trap",
        **kwargs,
    ) is not None


@pytest.mark.parametrize(
    ("raw", "normalized", "mode"),
    [
        ("M+HCOO-", "[M+HCOO]-", "Negative"),
        ("M+Cl-", "[M+Cl]-", "Negative"),
        ("M+H+", "[M+H]+", "Positive"),
        ("M-H-", "[M-H]-", "Negative"),
        ("M+2H2+", "[M+2H]2+", "Positive"),
    ],
)
def test_adduct_normalization_respects_terminal_charge(raw, normalized, mode):
    assert normalize_adduct(raw) == normalized
    assert _ion_mode_from_adduct(normalized) == mode


@pytest.mark.parametrize(
    ("raw", "expected"),
    [("1-", -1), ("2+", 2), ("+2", 2), ("-2", -2), ("+", 1), ("-", -1)],
)
def test_charge_normalization_accepts_suffix_and_prefix_signs(raw, expected):
    assert normalize_charge(raw) == expected


@pytest.mark.parametrize(
    ("raw", "expected"),
    [("1-", "NEG"), ("-2", "NEG"), (-1, "NEG"),
     ("2+", "POS"), ("+2", "POS"), (1, "POS"), (0, None)],
)
def test_library_charge_is_a_polarity_fallback(raw, expected):
    assert search_module._charge_polarity(raw) == expected


def test_mgf_writer_never_labels_unitless_rt_as_seconds(tmp_path):
    base = {
        "name": "RT example",
        "precursor_mz": 301.1,
        "fragments": [100.0, 150.0],
        "fragment_intensities": [100.0, 50.0],
        "retention_time": 2.0,
    }
    unknown = tmp_path / "unknown_rt.mgf"
    minutes = tmp_path / "minutes_rt.mgf"

    write_mgf([base], unknown)
    write_mgf([{**base, "retention_time_unit": "minutes"}], minutes)

    assert "RTINSECONDS=" not in unknown.read_text()
    assert "COMMENT=RT_UNIT: unknown" in unknown.read_text()
    assert "RTINSECONDS=120.00" in minutes.read_text()


def test_default_precursor_window_is_not_half_a_dalton():
    assert MZ_TOLERANCE_DA == pytest.approx(0.02)
    records = [
        _record(0, 500.0),
        _record(1, 500.019),
        _record(2, 500.021),
        _record(3, 500.4),
    ]
    mz_index = {}
    for index, record in enumerate(records):
        bucket = int(record["precursor_mz"] / search_module.MZ_BUCKET_DALTON)
        mz_index.setdefault(bucket, []).append(index)

    candidates = search_module._query_candidates(500.0, records, mz_index)

    assert candidates == [0, 1]


def test_source_and_polarity_are_filtered_before_dense_cap(monkeypatch):
    records = [
        _record(0, source_tier="predicted", polarity="POS"),
        _record(1, source_tier="experimental", polarity="NEG"),
        _record(2, source_tier="experimental", polarity="POS"),
        _record(3, source_tier="curated", polarity="POS"),
    ]
    mz_index = {5000: [0, 1, 2, 3]}
    seen = {}

    class StopAfterPrescreen(Exception):
        pass

    def spy(candidates, *_args, **_kwargs):
        seen["candidates"] = list(candidates)
        raise StopAfterPrescreen

    monkeypatch.setattr(search_module, "MAX_CANDIDATES_PER_QUERY", 1)
    monkeypatch.setattr(search_module, "_prescreen_dense_candidates", spy)

    with pytest.raises(StopAfterPrescreen):
        search_module.search(
            [[100.0, 1.0], [150.0, 0.8], [200.0, 0.6]],
            500.0,
            records=records,
            mz_index=mz_index,
            query_polarity="POS",
        )

    assert seen["candidates"] == [2, 3]


def test_search_cross_adduct_requires_metadata_and_reports_neutral_ppm(monkeypatch):
    record = _record(
        0, precursor_mz=301.007276, adduct="[M+H]+", neutral_mass=300.0
    )
    record.update({"has_structure": True, "smiles": "CC", "name_only_lead": False})
    records = [record]
    monkeypatch.setattr(search_module, "_lib_records", records)
    monkeypatch.setattr(search_module, "_mz_index", {})
    monkeypatch.setattr(search_module, "_nm_index", {3000: [0]})
    monkeypatch.setattr(search_module, "EXPERT_SYSTEM_ENABLED", False)
    monkeypatch.setattr(search_module, "USE_EXPERT_V3", False)
    monkeypatch.setattr(search_module, "USE_CONSISTENCY_FILTER", False)
    monkeypatch.setattr(search_module, "USE_ADAPTIVE_PEAK_SEARCH", False)
    monkeypatch.setattr(search_module, "USE_ISOMER_FRAGMENT_RERANK", False)
    monkeypatch.setattr(search_module, "MODEL_PREDICTOR_AVAILABLE", False)
    monkeypatch.setattr(search_module, "_ensure_clf", lambda: False)

    query = [[100.0, 1.0], [150.0, 0.8], [200.0, 0.6]]
    result = search_module.search(
        query,
        322.989218,
        query_adduct="[M+Na]+",
        query_adduct_reliable=True,
        query_polarity="POS",
        cos_thresh=0.0,
        final_thresh=0.0,
        strict_evidence_gate=False,
    )

    assert result is not None
    assert result["cross_adduct_match"] is True
    assert result["mass_match_basis"] == "neutral_mass"
    assert result["mass_ppm"] == pytest.approx(0.0, abs=1e-9)
    assert result["query_neutral_mass"] == pytest.approx(300.0, abs=1e-9)
    assert result["adduct_detected"] == "[M+Na]+"
    assert result["adduct_detection_source"] == "query_metadata"

    # A naked neutral mass (or an explicitly conflicted adduct) must not open
    # the neutral index and therefore cannot retrieve this distant precursor.
    assert search_module.search(
        query,
        322.989218,
        query_neutral_mass=300.0,
        query_polarity="POS",
        cos_thresh=0.0,
        final_thresh=0.0,
    ) is None
    assert search_module.search(
        query,
        322.989218,
        query_adduct="[M+Na]+",
        query_adduct_reliable=False,
        query_polarity="POS",
        cos_thresh=0.0,
        final_thresh=0.0,
    ) is None


def test_per_request_weight_override_changes_the_final_score(monkeypatch):
    record = _record(0)
    record.update({"has_structure": True, "smiles": "CC", "name_only_lead": False})
    records = [record]
    mz_index = {int(500.0 / search_module.MZ_BUCKET_DALTON): [0]}
    monkeypatch.setattr(search_module, "EXPERT_SYSTEM_ENABLED", False)
    monkeypatch.setattr(search_module, "USE_EXPERT_V3", False)
    monkeypatch.setattr(search_module, "USE_PATHWAY_ALIGN", False)
    monkeypatch.setattr(search_module, "USE_CONSISTENCY_FILTER", False)
    monkeypatch.setattr(search_module, "USE_ADAPTIVE_PEAK_SEARCH", False)
    monkeypatch.setattr(search_module, "USE_ISOMER_FRAGMENT_RERANK", False)
    monkeypatch.setattr(search_module, "MODEL_PREDICTOR_AVAILABLE", False)
    monkeypatch.setattr(search_module, "_ensure_clf", lambda: False)
    monkeypatch.setattr(search_module, "enhanced_composite_score", lambda *_args: 0.5)

    kwargs = dict(
        records=records,
        mz_index=mz_index,
        query_polarity="POS",
        cos_thresh=0.0,
        final_thresh=0.0,
        strict_evidence_gate=False,
    )
    query = [[100.0, 1.0], [150.0, 0.8], [200.0, 0.6]]
    default = search_module.search(query, 500.0, **kwargs)
    without_cosine_weight = search_module.search(query, 500.0, w_cos=0.0, **kwargs)

    assert default is not None and without_cosine_weight is not None
    assert without_cosine_weight["final_score"] < default["final_score"]


def test_runtime_structure_conflict_is_subject_to_name_only_evidence_gate(monkeypatch):
    record = _record(0)
    record.update({
        "has_structure": True,
        "name_only_lead": False,
        "smiles": "not-a-valid-smiles",
        "formula": "C2H6",
    })
    records = [record]
    mz_index = {int(500.0 / search_module.MZ_BUCKET_DALTON): [0]}
    monkeypatch.setattr(search_module, "EXPERT_SYSTEM_ENABLED", False)
    monkeypatch.setattr(search_module, "USE_EXPERT_V3", False)
    monkeypatch.setattr(search_module, "USE_CONSISTENCY_FILTER", False)
    monkeypatch.setattr(search_module, "USE_ADAPTIVE_PEAK_SEARCH", False)
    monkeypatch.setattr(search_module, "USE_ISOMER_FRAGMENT_RERANK", False)
    monkeypatch.setattr(search_module, "MODEL_PREDICTOR_AVAILABLE", False)
    monkeypatch.setattr(search_module, "_ensure_clf", lambda: False)

    result = search_module.search(
        [[100.0, 1.0], [150.0, 0.8]],
        500.0,
        records=records,
        mz_index=mz_index,
        query_polarity="POS",
        cos_thresh=0.0,
        final_thresh=0.0,
        strict_evidence_gate=False,
    )

    assert result is None


def test_mgf_rtinminutes_preserves_value_and_unit(tmp_path):
    mgf = tmp_path / "minutes.mgf"
    mgf.write_text(
        "BEGIN IONS\nPEPMASS=200\nRTINMINUTES=5.25\nIONMODE=POSITIVE\n50 10\nEND IONS\n",
        encoding="utf-8",
    )

    record = parse_mgf_file(mgf)[0]

    assert record["retention_time"] == pytest.approx(5.25)
    assert record["retention_time_unit"] == "minutes"
