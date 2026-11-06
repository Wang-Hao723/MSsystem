"""Regression tests for precursor-level expert evidence."""

import json

import pytest
import numpy as np

from utils.isotope import (
    detect_precursor_isotope_envelope,
    evaluate_formula_isotope_consistency,
    theoretical_formula_isotope_ratios,
)
from web.server import _apply_sample_consensus
from matching.search import _apply_formula_isotope_rerank
from utils.mzml_io import (
    _extract_mzml_precursor_evidence,
    evaluate_precursor_isolation_purity,
    extract_tic_overview,
    parse_mzml_file,
)


class _UnitFloat(float):
    def __new__(cls, value, unit):
        instance = super().__new__(cls, value)
        instance.unit_info = unit
        return instance


def test_ms1_precursor_envelope_respects_charge_spacing():
    precursor = 500.2
    spacing = 1.00335483507 / 2
    peaks = [
        [precursor, 1000.0],
        [precursor + spacing, 220.0],
        [precursor + spacing * 2, 45.0],
    ]

    evidence = detect_precursor_isotope_envelope(peaks, precursor, charge=2)

    assert evidence["source"] == "adjacent_ms1_scan"
    assert evidence["charge"] == 2
    assert evidence["n_isotopes_observed"] == 2
    assert evidence["peaks"][1]["relative_to_m"] == 0.22
    assert evidence["mz_tolerance_ppm"] == 5.0
    assert evidence["mz_tolerance_da"] == pytest.approx(0.002501)


def test_ms1_precursor_envelope_rejects_peak_outside_qe_ppm_window():
    precursor = 230.24877
    spacing = 1.00335483507
    peaks = [
        [precursor, 1000.0],
        [precursor + spacing + 0.006, 500.0],
    ]

    evidence = detect_precursor_isotope_envelope(peaks, precursor, charge=1)

    assert evidence["mz_tolerance_da"] == pytest.approx(0.001151)
    assert evidence["n_isotopes_observed"] == 0
    assert evidence["peaks"][1]["mz"] is None


def test_ms1_precursor_envelope_allows_explicit_absolute_tolerance():
    precursor = 230.24877
    spacing = 1.00335483507
    peaks = [
        [precursor, 1000.0],
        [precursor + spacing + 0.006, 500.0],
    ]

    evidence = detect_precursor_isotope_envelope(
        peaks, precursor, charge=1, tol_da=0.01
    )

    assert evidence["mz_tolerance_source"] == "explicit_da"
    assert evidence["mz_tolerance_ppm"] is None
    assert evidence["n_isotopes_observed"] == 1


def test_formula_isotope_evidence_supports_compatible_carbon_envelope():
    theory = theoretical_formula_isotope_ratios("C20H20O5")
    envelope = {
        "source": "adjacent_ms1_scan",
        "peaks": [
            {"order": 0, "relative_to_m": 1.0},
            {"order": 1, "relative_to_m": theory["m1_over_m"]},
            {"order": 2, "relative_to_m": theory["m2_over_m"]},
        ],
    }

    result = evaluate_formula_isotope_consistency("C20H20O5", envelope)

    assert result["status"] == "supporting"
    assert result["score"] > 0.95


def test_formula_isotope_evidence_reranks_without_inflating_final_score():
    theory = theoretical_formula_isotope_ratios("C30H40O5")
    envelope = {
        "source": "adjacent_ms1_scan",
        "status": "observed",
        "n_isotopes_observed": 2,
        "peaks": [
            {"order": 0, "relative_to_m": 1.0},
            {"order": 1, "relative_to_m": theory["m1_over_m"]},
            {"order": 2, "relative_to_m": theory["m2_over_m"]},
        ],
    }
    results = [
        {"formula": "C30H40O5", "final_score": 0.70, "rank_score": 0.70},
        {"formula": "C5H10O5", "final_score": 0.70, "rank_score": 0.70},
    ]

    _apply_formula_isotope_rerank(results, envelope)

    assert results[0]["rank_score"] > results[1]["rank_score"]
    assert results[0]["final_score"] == results[1]["final_score"] == 0.70
    assert abs(results[0]["isotope_rank_adjustment"]) <= 0.01


def test_precursor_isolation_purity_includes_target_isotope_cluster():
    evidence = evaluate_precursor_isolation_purity(
        [
            (499.6, 20.0),
            (500.0, 1000.0),
            (501.0033548, 120.0),
            (501.2, 40.0),
        ],
        500.0,
        0.5,
        1.5,
        charge=1,
    )

    assert evidence["status"] == "clean"
    assert evidence["purity"] == pytest.approx(1120 / 1180, rel=1e-5)
    assert evidence["coisolated_peak_count"] == 2


def test_precursor_isolation_purity_flags_chimeric_window():
    evidence = evaluate_precursor_isolation_purity(
        [(500.0, 100.0), (500.2, 900.0)],
        500.0,
        0.5,
        0.5,
    )

    assert evidence["status"] == "severe"
    assert evidence["purity"] == pytest.approx(0.1)
    assert "chimeric" in evidence["interpretation"].lower()


def test_explicit_full_isolation_width_fills_missing_mzml_offsets(monkeypatch):
    from pyteomics import mzml

    spectra = [
        {
            "id": "scan=1",
            "ms level": 1,
            "m/z array": np.asarray([499.8, 500.0, 500.2]),
            "intensity array": np.asarray([10.0, 100.0, 90.0]),
        },
        {
            "id": "scan=2",
            "ms level": 2,
            "scanList": {"scan": [{
                "scan start time": _UnitFloat(90.0, "second"),
            }]},
            "precursorList": {"precursor": [{
                "spectrumRef": "scan=1",
                "selectedIonList": {"selectedIon": [{
                    "selected ion m/z": 500.0,
                    "charge state": 1,
                }]},
                "isolationWindow": {"isolation window target m/z": 500.0},
                "activation": {},
            }]},
        },
    ]

    class Reader:
        def __enter__(self):
            return iter(spectra)

        def __exit__(self, *_args):
            return False

    monkeypatch.setattr(mzml, "read", lambda _path: Reader())
    evidence_index = _extract_mzml_precursor_evidence(
        "missing-offsets.mzML",
        isolation_window_width_da=1.2,
        include_ordinal_index=True,
    )
    evidence = evidence_index[2]

    assert evidence["isolation_window_lower_offset"] == pytest.approx(0.6)
    assert evidence["isolation_window_upper_offset"] == pytest.approx(0.6)
    assert evidence["retention_time"] == pytest.approx(1.5)
    assert evidence["retention_time_unit"] == "minutes"
    assert evidence["native_scan_id"] == 2
    assert evidence_index[("ms2_ordinal", 0)] is evidence
    assert evidence["precursor_isolation_evidence"][
        "isolation_window_source"
    ] == "user_full_width_override"
    assert evidence["precursor_isolation_evidence"]["status"] == "mixed"


def test_titleless_matchms_spectrum_joins_native_ms1_evidence(monkeypatch):
    import utils.mzml_io as mzml_io

    class FakePeaks:
        mz = np.asarray([100.0, 150.0])
        intensities = np.asarray([10.0, 20.0])

    class FakeSpectrum:
        metadata = {
            "precursor_mz": 500.0,
            "title": None,
            "scan_number": None,
            "scan_start_time": [1.25],
        }
        peaks = FakePeaks()

        def get(self, key):
            return self.metadata.get(key)

    evidence = {
        "native_scan_id": 42,
        "native_id": "scanId=42",
        "precursor_charge": 1,
        "spectrum_representation": "centroid",
        "spectrum_representation_source": "mzml_cv",
        "isotope_evidence_status": "observed_from_ms1",
        "precursor_isolation_evidence": {"status": "clean", "purity": 0.95},
    }
    evidence_index = {
        42: evidence,
        ("ms2_ordinal", 0): evidence,
    }
    monkeypatch.setattr(mzml_io, "_file_polarity", lambda _path: "NEG")
    monkeypatch.setattr(
        mzml_io,
        "_extract_mzml_precursor_evidence",
        lambda *_args, **_kwargs: evidence_index,
    )
    monkeypatch.setattr(
        mzml_io,
        "load_from_mzml",
        lambda _path: iter([FakeSpectrum()]),
    )

    records = parse_mzml_file("titleless.mzML")

    assert len(records) == 1
    assert records[0]["spectrum_id"] == 42
    assert records[0]["scan_number"] == 42
    assert records[0]["title"] == "scanId=42"
    assert records[0]["precursor_charge"] == 1
    assert records[0]["isotope_evidence_status"] == "observed_from_ms1"


def test_matchms_scan_start_time_seconds_are_normalized_to_minutes(monkeypatch):
    import utils.mzml_io as mzml_io

    class FakePeaks:
        mz = np.asarray([100.0, 150.0])
        intensities = np.asarray([10.0, 20.0])

    class FakeSpectrum:
        metadata = {
            "precursor_mz": 500.0,
            "title": "scanId=42",
            "scan_start_time": [_UnitFloat(90.0, "second")],
        }
        peaks = FakePeaks()

        def get(self, key):
            return self.metadata.get(key)

    evidence = {
        "native_scan_id": 42,
        "native_id": "scanId=42",
        "spectrum_representation": "centroid",
        "spectrum_representation_source": "mzml_cv",
    }
    monkeypatch.setattr(mzml_io, "_file_polarity", lambda _path: "NEG")
    monkeypatch.setattr(
        mzml_io,
        "_extract_mzml_precursor_evidence",
        lambda *_args, **_kwargs: {
            42: evidence,
            ("ms2_ordinal", 0): evidence,
        },
    )
    monkeypatch.setattr(
        mzml_io,
        "load_from_mzml",
        lambda _path: iter([FakeSpectrum()]),
    )

    records = parse_mzml_file("openms-seconds.mzML")

    assert records[0]["retention_time"] == pytest.approx(1.5)
    assert records[0]["retention_time_unit"] == "minutes"


def test_tic_overview_normalizes_openms_second_cv_unit(monkeypatch):
    from pyteomics import mzml

    spectra = [{
        "id": "scan=1",
        "ms level": 1,
        "scanList": {"scan": [{
            "scan start time": _UnitFloat(120.0, "second"),
        }]},
        "total ion current": 42.0,
        "base peak intensity": 18.0,
        "base peak m/z": 321.1234,
    }]

    class Reader:
        def __enter__(self):
            return iter(spectra)

        def __exit__(self, *_args):
            return False

    read_options = {}

    def fake_read(_path, **kwargs):
        read_options.update(kwargs)
        return Reader()

    monkeypatch.setattr(mzml, "read", fake_read)

    overview = extract_tic_overview("openms-seconds.mzML")

    assert overview["rt_unit"] == "minute"
    assert overview["ms1"][0]["rt"] == pytest.approx(2.0)
    assert overview["ms1"][0]["bpc"] == pytest.approx(18.0)
    assert overview["ms1"][0]["base_peak_mz"] == pytest.approx(321.1234)
    assert overview["bpc_available"] is True
    assert overview["chromatogram_schema_version"] == 3
    assert read_options["decode_binary"] is False


def test_legacy_tic_overview_is_lazily_backfilled_with_real_ms1_bpc(
    tmp_path,
    monkeypatch,
):
    from web import server

    source_path = tmp_path / "sample.mzML"
    source_path.write_text("", encoding="utf-8")
    artifact_path = tmp_path / "legacy_tic_overview.json"
    artifact_path.write_text(json.dumps({
        "ms1": [{"rt": 1.0, "tic": 42.0, "scan_id": 1}],
        "ms2": [],
        "rt_unit": "minute",
    }), encoding="utf-8")
    enriched = {
        "chromatogram_schema_version": 2,
        "ms1": [{
            "rt": 1.0,
            "tic": 42.0,
            "bpc": 18.0,
            "base_peak_mz": 321.1234,
            "scan_id": 1,
        }],
        "ms2": [],
        "rt_unit": "minute",
        "bpc_available": True,
        "bpc_scope": "MS1 base peak intensity per survey scan",
    }
    monkeypatch.setattr(
        "utils.mzml_io.extract_tic_overview",
        lambda path: enriched if path == source_path else None,
    )

    result = server._load_tic_overview_with_bpc_backfill(
        "task-bpc",
        {"file_path": str(source_path)},
        artifact_path,
    )

    assert result["bpc_available"] is True
    assert result["bpc_backfilled"] is True
    assert result["ms1"][0]["bpc"] == pytest.approx(18.0)
    persisted = json.loads(artifact_path.read_text(encoding="utf-8"))
    assert persisted["bpc_backfill_task_id"] == "task-bpc"


@pytest.mark.parametrize(
    ("representation", "message"),
    [
        ("profile", "unpicked profile MS2 spectrum"),
        ("unknown", "omits centroid/profile CV terms"),
    ],
)
def test_mzml_matching_fails_closed_for_unpicked_or_unknown_ms2_representation(
    monkeypatch,
    representation,
    message,
):
    import utils.mzml_io as mzml_io

    class FakePeaks:
        mz = np.asarray([100.0, 150.0])
        intensities = np.asarray([10.0, 20.0])

    class FakeSpectrum:
        metadata = {"precursor_mz": 500.0, "title": "scanId=42"}
        peaks = FakePeaks()

        def get(self, key):
            return self.metadata.get(key)

    evidence = {
        "native_scan_id": 42,
        "native_id": "scanId=42",
        "spectrum_representation": representation,
        "spectrum_representation_source": (
            "mzml_cv" if representation == "profile" else "missing_cv"
        ),
    }
    monkeypatch.setattr(mzml_io, "_file_polarity", lambda _path: "NEG")
    monkeypatch.setattr(
        mzml_io,
        "_extract_mzml_precursor_evidence",
        lambda *_args, **_kwargs: {
            42: evidence,
            ("ms2_ordinal", 0): evidence,
        },
    )
    monkeypatch.setattr(
        mzml_io,
        "load_from_mzml",
        lambda _path: iter([FakeSpectrum()]),
    )

    with pytest.raises(ValueError, match=message):
        parse_mzml_file("representation.mzML")


def test_audited_assume_centroid_policy_can_accept_missing_cv(monkeypatch):
    import utils.mzml_io as mzml_io

    class FakePeaks:
        mz = np.asarray([100.0, 150.0])
        intensities = np.asarray([10.0, 20.0])

    class FakeSpectrum:
        metadata = {"precursor_mz": 500.0, "title": "scanId=42"}
        peaks = FakePeaks()

        def get(self, key):
            return self.metadata.get(key)

    evidence = {
        "native_scan_id": 42,
        "native_id": "scanId=42",
        "spectrum_representation": "unknown",
        "spectrum_representation_source": "missing_cv",
    }
    monkeypatch.setattr(mzml_io, "_file_polarity", lambda _path: "NEG")
    monkeypatch.setattr(
        mzml_io,
        "_extract_mzml_precursor_evidence",
        lambda *_args, **_kwargs: {
            42: evidence,
            ("ms2_ordinal", 0): evidence,
        },
    )
    monkeypatch.setattr(
        mzml_io,
        "load_from_mzml",
        lambda _path: iter([FakeSpectrum()]),
    )

    records = parse_mzml_file(
        "representation.mzML",
        allow_unknown_representation=True,
    )

    assert records[0]["spectrum_representation"] == "unknown"
    assert records[0]["spectrum_representation_source"] == "missing_cv"


def _consensus_row(rt, collision_energy):
    return {
        "retention_time": rt,
        "collision_energy": collision_energy,
        "match": {
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "has_structure": True,
            "confidence": "Medium",
            "evidence_tier": "Medium",
            "cosine": 0.85,
            "final_score": 0.62,
            "mass_ppm": 3.0,
            "class_agreement": 1.0,
            "isomer_ambiguous": False,
            "compound_name": "Example",
        },
    }


def test_sample_consensus_uses_local_rt_clusters_and_multi_energy_support():
    rows = [
        _consensus_row(1.00, 20.0),
        _consensus_row(1.02, 40.0),
        _consensus_row(5.00, 20.0),
    ]

    updated = _apply_sample_consensus(rows)

    assert updated[0]["match"]["evidence_tier"] == "Strong Candidate"
    assert updated[1]["match"]["multi_energy_consensus"] is True
    assert updated[1]["match"]["sample_consensus_support"]["n_spectra"] == 2
    assert updated[2]["match"]["evidence_tier"] == "Medium"


def test_sample_consensus_does_not_mix_evidence_across_failing_scans():
    strong_cos = _consensus_row(1.00, 20.0)
    strong_cos["match"].update(cosine=0.90, final_score=0.40, mass_ppm=3.0)
    strong_final = _consensus_row(1.02, 40.0)
    strong_final["match"].update(cosine=0.60, final_score=0.70, mass_ppm=3.0)

    updated = _apply_sample_consensus([strong_cos, strong_final])

    assert all(row["match"]["evidence_tier"] == "Medium" for row in updated)
    assert all("sample_consensus_support" not in row["match"] for row in updated)
