from pathlib import Path

import numpy as np
import pytest

from matching.blank_subtraction import (
    BlankProfile,
    assess_blank_features,
    assess_blank_replicates,
    assess_blank_role_groups,
    blank_polarity_coverage,
    create_blank_profile_from_list,
    load_blank_mgf,
    subtract_blank,
)


SAMPLE_PEAKS = [[100.0, 5.0], [150.0, 50.0]]


def _context_profile(*, precursor=500.0, rt=5.0):
    return create_blank_profile_from_list(
        [(100.0, 100.0)],
        precursor_mz=precursor,
        retention_time=rt,
        retention_time_unit="minutes",
    )


def test_context_free_legacy_profile_is_inert():
    profile = create_blank_profile_from_list([(100.0, 1000.0)])

    assert subtract_blank(SAMPLE_PEAKS, profile, precursor_mz=500.0) == SAMPLE_PEAKS


def test_missing_sample_precursor_is_inert():
    profile = _context_profile()

    assert subtract_blank(
        SAMPLE_PEAKS,
        profile,
        retention_time=5.0,
        retention_time_unit="minutes",
    ) == SAMPLE_PEAKS


def test_same_fragment_at_unrelated_precursor_is_not_removed():
    profile = _context_profile(precursor=700.0, rt=5.0)

    assert subtract_blank(
        SAMPLE_PEAKS,
        profile,
        precursor_mz=500.0,
        retention_time=5.0,
        retention_time_unit="minutes",
    ) == SAMPLE_PEAKS


def test_same_precursor_but_different_rt_is_not_removed():
    profile = _context_profile(precursor=500.0, rt=8.0)

    assert subtract_blank(
        SAMPLE_PEAKS,
        profile,
        precursor_mz=500.0,
        retention_time=5.0,
        retention_time_unit="minutes",
    ) == SAMPLE_PEAKS


def test_blank_rt_requires_sample_rt_and_known_unit():
    profile = _context_profile()

    assert subtract_blank(SAMPLE_PEAKS, profile, precursor_mz=500.0) == SAMPLE_PEAKS
    assert subtract_blank(
        SAMPLE_PEAKS,
        profile,
        precursor_mz=500.0,
        retention_time=5.0,
        retention_time_unit="cycles",
    ) == SAMPLE_PEAKS


def test_blank_polarity_coverage_rejects_known_opposite_polarity():
    neg_profile = create_blank_profile_from_list([(100.0, 1000.0)])
    neg_profile.available_polarities = frozenset({"NEG"})
    coverage = blank_polarity_coverage(
        [{"precursor_mz": 200.0, "polarity": "POS"}],
        [neg_profile],
    )
    assert coverage == {
        "sample_polarities": ["POS"],
        "blank_polarities": ["NEG"],
        "unknown_blank_profiles": 0,
        "missing_sample_polarities": ["POS"],
        "applicable": False,
    }


def test_blank_polarity_coverage_keeps_unknown_mgf_profile_compatible():
    unknown_profile = create_blank_profile_from_list([(100.0, 1000.0)])
    coverage = blank_polarity_coverage(
        [{"precursor_mz": 200.0, "polarity": "POS"}],
        [unknown_profile],
    )
    assert coverage["applicable"] is True
    assert coverage["unknown_blank_profiles"] == 1


def test_aligned_precursor_and_rt_allow_subtraction():
    profile = _context_profile()

    result = subtract_blank(
        SAMPLE_PEAKS,
        profile,
        precursor_mz=500.002,
        retention_time=306.0,
        retention_time_unit="seconds",
    )

    assert result == [[150.0, 50.0]]


def test_precursor_only_is_allowed_when_blank_has_no_rt():
    profile = create_blank_profile_from_list(
        [(100.0, 100.0)],
        precursor_mz=500.0,
    )

    assert subtract_blank(
        SAMPLE_PEAKS,
        profile,
        precursor_mz=500.0,
    ) == [[150.0, 50.0]]


def test_mgf_loader_keeps_spectrum_context_and_does_not_use_global_max(tmp_path: Path):
    mgf = tmp_path / "blank.mgf"
    mgf.write_text(
        """BEGIN IONS
PEPMASS=500.0
RTINSECONDS=300
100.0 10
END IONS
BEGIN IONS
PEPMASS=700.0
RTINSECONDS=300
100.0 1000
END IONS
""",
        encoding="utf-8",
    )
    profile = load_blank_mgf(mgf)

    result = subtract_blank(
        [[100.0, 50.0]],
        profile,
        precursor_mz=500.0,
        retention_time=5.0,
        retention_time_unit="minutes",
    )

    assert len(profile.context_spectra) == 2
    assert result[0][0] == pytest.approx(100.0)
    assert result[0][1] == pytest.approx(40.0)


def test_feature_filter_mgf_fallback_excludes_below_ratio_and_respects_polarity(tmp_path: Path):
    mgf = tmp_path / "blank.mgf"
    mgf.write_text(
        """BEGIN IONS
PEPMASS=500.0
RTINMINUTES=5
IONMODE=POS
PRECURSORINTENSITY=100
100 10
150 20
END IONS
""",
        encoding="utf-8",
    )
    profile = load_blank_mgf(mgf)
    samples = [
        {"precursor_mz": 500.002, "retention_time": 5.05,
         "retention_time_unit": "minutes", "polarity": "POS",
         "precursor_intensity": 200, "peaks": [[100, 20], [150, 30]]},
        {"precursor_mz": 500.002, "retention_time": 5.05,
         "retention_time_unit": "minutes", "polarity": "POS",
         "precursor_intensity": 400, "peaks": [[100, 20], [150, 30]]},
        {"precursor_mz": 500.002, "retention_time": 5.05,
         "retention_time_unit": "minutes", "polarity": "NEG",
         "precursor_intensity": 100, "peaks": [[100, 20], [150, 30]]},
    ]

    assessed = assess_blank_features(samples, profile, sample_to_blank_ratio=3.0)

    assert assessed[0]["blank_status"] == "excluded_below_ratio"
    assert assessed[0]["blank_excluded"] is True
    assert assessed[0]["blank_sample_ratio"] == pytest.approx(2.0)
    assert assessed[1]["blank_status"] == "retained_above_ratio"
    assert assessed[1]["blank_excluded"] is False
    assert assessed[2]["blank_status"] == "blank_polarity_not_applicable"


def test_feature_filter_mgf_without_precursor_intensity_is_audited_but_not_excluded(
    tmp_path: Path,
):
    mgf = tmp_path / "blank.mgf"
    mgf.write_text(
        "BEGIN IONS\nPEPMASS=500\nRTINMINUTES=5\nIONMODE=POS\n100 1000\nEND IONS\n",
        encoding="utf-8",
    )
    profile = load_blank_mgf(mgf)

    assessed = assess_blank_features([{
        "precursor_mz": 500.0,
        "retention_time": 5.0,
        "retention_time_unit": "minutes",
        "polarity": "POS",
        "peaks": [[100, 10]],
    }], profile)

    assert assessed[0]["blank_feature_matched"] is True
    assert assessed[0]["blank_status"] == "sample_intensity_missing"
    assert assessed[0]["blank_excluded"] is False
    assert assessed[0]["blank_evidence_source"] == "ms2_precursor_intensity"


def test_feature_filter_uses_blank_ms1_even_without_blank_ms2(monkeypatch, tmp_path: Path):
    blank_path = tmp_path / "blank.mzML"
    blank_path.write_text("placeholder")
    profile = BlankProfile(
        np.asarray([], dtype=float), np.asarray([], dtype=float),
        np.asarray([], dtype=float), np.asarray([], dtype=float),
        np.asarray([], dtype=int), source_path=blank_path, ms1_scan_count=1,
    )
    scan = {
        "ms level": 1,
        "positive scan": True,
        "scanList": {"scan": [{"scan start time": 5.0}]},
        "m/z array": np.asarray([499.0, 500.001, 501.0]),
        "intensity array": np.asarray([1.0, 100.0, 1.0]),
    }

    class Reader:
        def __enter__(self): return iter([scan])
        def __exit__(self, *args): return False

    from pyteomics import mzml
    monkeypatch.setattr(mzml, "read", lambda _: Reader())
    monkeypatch.setattr("utils.mzml_io._file_polarity", lambda _: "POS")
    samples = [
        {"precursor_mz": 500.0, "retention_time": 5.05,
         "retention_time_unit": "minutes", "polarity": "POS",
         "precursor_ms1_intensity": 250.0, "peaks": [[100, 10]]},
        {"precursor_mz": 500.0, "retention_time": 5.05,
         "retention_time_unit": "minutes", "polarity": "NEG",
         "precursor_ms1_intensity": 100.0, "peaks": [[100, 10]]},
    ]

    assessed = assess_blank_features(samples, profile, sample_to_blank_ratio=3.0)

    assert assessed[0]["blank_evidence_source"] == "blank_ms1"
    assert assessed[0]["blank_sample_ratio"] == pytest.approx(2.5)
    assert assessed[0]["blank_excluded"] is True
    assert assessed[1]["blank_status"] == "blank_polarity_not_applicable"


def test_feature_filter_does_not_exclude_when_sample_ms1_intensity_is_missing(monkeypatch, tmp_path: Path):
    blank_path = tmp_path / "blank.mzML"
    blank_path.write_text("placeholder")
    profile = BlankProfile(
        np.asarray([], dtype=float), np.asarray([], dtype=float),
        np.asarray([], dtype=float), np.asarray([], dtype=float),
        np.asarray([], dtype=int), source_path=blank_path, ms1_scan_count=1,
    )
    scan = {
        "ms level": 1,
        "positive scan": True,
        "scanList": {"scan": [{"scan start time": 5.0}]},
        "m/z array": np.asarray([500.0]),
        "intensity array": np.asarray([100.0]),
    }

    class Reader:
        def __enter__(self): return iter([scan])
        def __exit__(self, *args): return False

    from pyteomics import mzml
    monkeypatch.setattr(mzml, "read", lambda _: Reader())
    monkeypatch.setattr("utils.mzml_io._file_polarity", lambda _: "POS")
    assessed = assess_blank_features([
        {"precursor_mz": 500.0, "retention_time": 5.0,
         "retention_time_unit": "minutes", "polarity": "POS", "peaks": [[100, 10]]},
    ], profile)

    assert assessed[0]["blank_feature_matched"] is True
    assert assessed[0]["blank_status"] == "sample_intensity_missing"
    assert assessed[0]["blank_excluded"] is False


def test_symmetric_xic_uses_same_apex_statistic_for_sample_and_blank(
    monkeypatch, tmp_path: Path,
):
    sample_path = tmp_path / "sample.mzML"
    blank_path = tmp_path / "blank.mzML"
    sample_path.write_text("sample")
    blank_path.write_text("blank")
    profile = BlankProfile(
        np.asarray([], dtype=float), np.asarray([], dtype=float),
        np.asarray([], dtype=float), np.asarray([], dtype=float),
        np.asarray([], dtype=int), source_path=blank_path, ms1_scan_count=2,
        available_polarities=["POS"],
    )

    def scan(rt, intensity):
        return {
            "ms level": 1,
            "positive scan": True,
            "scanList": {"scan": [{"scan start time": rt}]},
            "m/z array": np.asarray([500.0]),
            "intensity array": np.asarray([intensity]),
        }

    class Reader:
        def __init__(self, scans): self.scans = scans
        def __enter__(self): return iter(self.scans)
        def __exit__(self, *args): return False

    from pyteomics import mzml

    def read(path):
        if Path(path) == sample_path:
            return Reader([scan(4.9, 40.0), scan(5.0, 400.0), scan(5.1, 40.0)])
        return Reader([scan(4.9, 10.0), scan(5.0, 100.0), scan(5.1, 10.0)])

    monkeypatch.setattr(mzml, "read", read)
    monkeypatch.setattr("utils.mzml_io._file_polarity", lambda _: "POS")
    assessed = assess_blank_features(
        [{
            "precursor_mz": 500.0,
            "retention_time": 5.0,
            "retention_time_unit": "minutes",
            "polarity": "POS",
            # This deliberately disagrees with the XIC apex and must not be used.
            "precursor_ms1_intensity": 5.0,
            "peaks": [[100, 10]],
        }],
        profile,
        sample_path=sample_path,
        sample_to_blank_ratio=3.0,
    )

    assert assessed[0]["blank_evidence_source"] == "paired_ms1_xic_apex"
    assert assessed[0]["blank_quantitation_mode"] == "symmetric_xic_apex"
    assert assessed[0]["blank_sample_intensity"] == pytest.approx(400.0)
    assert assessed[0]["blank_intensity"] == pytest.approx(100.0)
    assert assessed[0]["blank_sample_ratio"] == pytest.approx(4.0)
    assert assessed[0]["blank_status"] == "retained_above_ratio"
    assert assessed[0]["blank_sample_xic_area"] == pytest.approx(44.0)
    assert assessed[0]["blank_xic_area"] == pytest.approx(11.0)


def test_replicate_consensus_retains_one_off_blank_spike(tmp_path: Path):
    sample = [{
        "precursor_mz": 500.0,
        "retention_time": 5.0,
        "retention_time_unit": "minutes",
        "polarity": "POS",
        "precursor_intensity": 200.0,
        "peaks": [[100, 10]],
    }]

    def profile(name, precursor):
        path = tmp_path / name
        path.write_text(
            "BEGIN IONS\nPEPMASS=%s\nRTINMINUTES=5\nIONMODE=POS\n"
            "PRECURSORINTENSITY=100\n100 10\nEND IONS\n" % precursor,
            encoding="utf-8",
        )
        return load_blank_mgf(path)

    assessed = assess_blank_replicates(
        sample,
        [profile("blank1.mgf", 500.0), profile("blank2.mgf", 700.0)],
        sample_to_blank_ratio=3.0,
        min_detection_fraction=0.5,
    )

    assert assessed[0]["blank_feature_matched"] is True
    assert assessed[0]["blank_replicates_applicable"] == 2
    assert assessed[0]["blank_replicates_detected"] == 1
    assert assessed[0]["blank_required_detection_count"] == 2
    assert assessed[0]["blank_status"] == "retained_low_blank_prevalence"
    assert assessed[0]["blank_excluded"] is False


def test_replicate_consensus_uses_median_blank_intensity(tmp_path: Path):
    sample = [{
        "precursor_mz": 500.0,
        "retention_time": 5.0,
        "retention_time_unit": "minutes",
        "polarity": "POS",
        "precursor_intensity": 250.0,
        "peaks": [[100, 10]],
    }]

    profiles = []
    for index, intensity in enumerate((80, 100, 1000), start=1):
        path = tmp_path / f"blank{index}.mgf"
        path.write_text(
            "BEGIN IONS\nPEPMASS=500\nRTINMINUTES=5\nIONMODE=POS\n"
            f"PRECURSORINTENSITY={intensity}\n100 10\nEND IONS\n",
            encoding="utf-8",
        )
        profiles.append(load_blank_mgf(path))

    assessed = assess_blank_replicates(
        sample,
        profiles,
        sample_to_blank_ratio=3.0,
        min_detection_fraction=0.5,
    )

    assert assessed[0]["blank_replicates_detected"] == 3
    assert assessed[0]["blank_intensity"] == pytest.approx(100.0)
    assert assessed[0]["blank_sample_ratio"] == pytest.approx(2.5)
    assert assessed[0]["blank_status"] == "excluded_below_ratio"


def test_soft_blank_policy_keeps_sample_above_blank_but_below_threefold(
    tmp_path: Path,
):
    blank = tmp_path / "method_blank.mgf"
    blank.write_text(
        "BEGIN IONS\nPEPMASS=500\nRTINMINUTES=5\nIONMODE=POS\n"
        "PRECURSORINTENSITY=100\n100 10\nEND IONS\n",
        encoding="utf-8",
    )
    sample = [{
        "precursor_mz": 500.0,
        "retention_time": 5.0,
        "retention_time_unit": "minutes",
        "polarity": "POS",
        "precursor_intensity": 200.0,
        "peaks": [[100, 10]],
    }]

    assessed = assess_blank_features(
        sample,
        load_blank_mgf(blank),
        sample_to_blank_ratio=3.0,
        filter_policy="soft_evidence",
        blank_role="method_blank",
    )

    assert assessed[0]["blank_sample_ratio"] == pytest.approx(2.0)
    assert assessed[0]["blank_evidence_category"] == "blank_associated"
    assert assessed[0]["blank_status"] == "retained_blank_associated"
    assert assessed[0]["blank_flagged"] is True
    assert assessed[0]["blank_excluded"] is False
    assert assessed[0]["blank_clean_view_excluded"] is False


def test_soft_clean_view_requires_replicated_typed_background_blank(
    tmp_path: Path,
):
    profiles = []
    for index in (1, 2):
        path = tmp_path / f"method_blank_{index}.mgf"
        path.write_text(
            "BEGIN IONS\nPEPMASS=500\nRTINMINUTES=5\nIONMODE=POS\n"
            "PRECURSORINTENSITY=100\n100 10\nEND IONS\n",
            encoding="utf-8",
        )
        profiles.append(load_blank_mgf(path))
    sample = [{
        "precursor_mz": 500.0,
        "retention_time": 5.0,
        "retention_time_unit": "minutes",
        "polarity": "POS",
        "precursor_intensity": 80.0,
        "peaks": [[100, 10]],
    }]

    one_blank = assess_blank_replicates(
        sample,
        profiles[:1],
        filter_policy="soft_evidence",
        blank_role="method_blank",
    )[0]
    repeated = assess_blank_replicates(
        sample,
        profiles,
        filter_policy="soft_evidence",
        blank_role="method_blank",
    )[0]

    assert one_blank["blank_evidence_category"] == "blank_dominant"
    assert one_blank["blank_clean_view_excluded"] is False
    assert repeated["blank_evidence_category"] == "blank_dominant"
    assert repeated["blank_replicates_detected"] == 2
    assert repeated["blank_clean_view_excluded"] is True
    assert repeated["blank_excluded"] is False


def test_carryover_blanks_never_hide_identification(tmp_path: Path):
    profiles = []
    for index in (1, 2):
        path = tmp_path / f"carryover_{index}.mgf"
        path.write_text(
            "BEGIN IONS\nPEPMASS=500\nRTINMINUTES=5\nIONMODE=POS\n"
            "PRECURSORINTENSITY=100\n100 10\nEND IONS\n",
            encoding="utf-8",
        )
        profiles.append(load_blank_mgf(path))
    sample = [{
        "precursor_mz": 500.0,
        "retention_time": 5.0,
        "retention_time_unit": "minutes",
        "polarity": "POS",
        "precursor_intensity": 50.0,
        "peaks": [[100, 10]],
    }]

    assessed = assess_blank_role_groups(
        sample,
        profiles,
        ["carryover_blank", "carryover_blank"],
        filter_policy="soft_evidence",
    )[0]

    assert assessed["blank_evidence_role"] == "carryover_blank"
    assert assessed["blank_carryover_suspected"] is True
    assert assessed["blank_flagged"] is True
    assert assessed["blank_status"] == "retained_carryover_evidence"
    assert assessed["blank_clean_view_excluded"] is False
    assert assessed["blank_excluded"] is False


def test_different_blank_roles_are_not_pooled_as_replicates(tmp_path: Path):
    profiles = []
    for name in ("method", "solvent"):
        path = tmp_path / f"{name}.mgf"
        path.write_text(
            "BEGIN IONS\nPEPMASS=500\nRTINMINUTES=5\nIONMODE=POS\n"
            "PRECURSORINTENSITY=100\n100 10\nEND IONS\n",
            encoding="utf-8",
        )
        profiles.append(load_blank_mgf(path))
    sample = [{
        "precursor_mz": 500.0,
        "retention_time": 5.0,
        "retention_time_unit": "minutes",
        "polarity": "POS",
        "precursor_intensity": 50.0,
        "peaks": [[100, 10]],
    }]

    assessed = assess_blank_role_groups(
        sample,
        profiles,
        ["method_blank", "solvent_blank"],
        filter_policy="soft_evidence",
    )[0]

    assert assessed["blank_replicates_detected"] == 2
    assert assessed["blank_clean_view_excluded"] is False
    assert set(assessed["blank_role_evidence"]) == {
        "method_blank", "solvent_blank",
    }
