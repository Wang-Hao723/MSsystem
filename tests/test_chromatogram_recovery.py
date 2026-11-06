"""Regression coverage for mzML TIC/BPC recovery and honest display rules."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from utils import mzml_io


class _SpectrumReader:
    def __init__(self, spectra):
        self.spectra = spectra

    def __enter__(self):
        return iter(self.spectra)

    def __exit__(self, *_args):
        return False


class _ChromatogramReader:
    def __init__(self, chromatograms):
        self.chromatograms = chromatograms

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def iterfind(self, name):
        assert name == "chromatogram"
        return iter(self.chromatograms)


def _ms2_spectrum():
    return {
        "id": "scan=2",
        "ms level": 2,
        "scanList": {"scan": [{"scan start time": 1.1}]},
        "total ion current": 25.0,
        "base peak intensity": 12.0,
    }


def _tic_chromatogram():
    return {
        "id": "TIC",
        "total ion current chromatogram": "",
        "time array": np.asarray([1.0, 1.1, 1.2, 1.3]),
        "intensity array": np.asarray([100.0, 25.0, 300.0, 40.0]),
        "ms level": np.asarray([1, 2, 1, 2]),
    }


def test_source_tic_recovers_ms1_trace_without_inventing_bpc(monkeypatch):
    from pyteomics import mzml

    monkeypatch.setattr(
        mzml,
        "read",
        lambda _path, **_kwargs: _SpectrumReader([_ms2_spectrum()]),
    )
    monkeypatch.setattr(
        mzml,
        "MzML",
        lambda _path, **_kwargs: _ChromatogramReader([_tic_chromatogram()]),
    )
    monkeypatch.setattr(
        mzml_io, "_mzml_chromatogram_time_unit", lambda _path: "minute",
    )

    overview = mzml_io.extract_tic_overview("ms2-spectra-with-source-tic.mzML")

    assert overview["chromatogram_schema_version"] == 3
    assert [point["rt"] for point in overview["ms1"]] == [1.0, 1.2]
    assert [point["tic"] for point in overview["ms1"]] == [100.0, 300.0]
    assert all(point["bpc"] is None for point in overview["ms1"])
    assert overview["ms1_spectra_available"] is False
    assert overview["ms1_spectrum_count"] == 0
    assert overview["ms1_tic_point_count"] == 2
    assert overview["tic_source"] == "mzml_tic_chromatogram_ms_level_array"
    assert overview["bpc_available"] is False
    assert overview["bpc_source"] == "unavailable"


def test_explicit_source_bpc_chromatogram_is_aligned_to_ms1_points(monkeypatch):
    from pyteomics import mzml

    bpc = {
        "id": "BPC",
        "basepeak chromatogram": "",
        "time array": np.asarray([1.0, 1.1, 1.2, 1.3]),
        "intensity array": np.asarray([55.0, 12.0, 180.0, 20.0]),
    }
    monkeypatch.setattr(
        mzml,
        "read",
        lambda _path, **_kwargs: _SpectrumReader([_ms2_spectrum()]),
    )
    monkeypatch.setattr(
        mzml,
        "MzML",
        lambda _path, **_kwargs: _ChromatogramReader([_tic_chromatogram(), bpc]),
    )
    monkeypatch.setattr(
        mzml_io, "_mzml_chromatogram_time_unit", lambda _path: "minute",
    )

    overview = mzml_io.extract_tic_overview("source-tic-and-bpc.mzML")

    assert [point["bpc"] for point in overview["ms1"]] == [55.0, 180.0]
    assert overview["bpc_available"] is True
    assert overview["bpc_source"] == "mzml_basepeak_chromatogram"


def test_schema2_failed_backfill_retries_and_persists_source_tic(
    tmp_path,
    monkeypatch,
):
    from web import server

    source_path = tmp_path / "sample.mzML"
    source_path.write_text("", encoding="utf-8")
    artifact_path = tmp_path / "tic_overview.json"
    artifact_path.write_text(json.dumps({
        "chromatogram_schema_version": 2,
        "ms1": [],
        "ms2": [{"rt": 1.1, "tic": 25.0, "scan_id": 2}],
        "bpc_available": False,
        "bpc_backfill_attempted": True,
    }), encoding="utf-8")
    enriched = {
        "chromatogram_schema_version": 3,
        "ms1": [
            {"rt": 1.0, "tic": 100.0, "bpc": None, "scan_id": 0},
            {"rt": 1.2, "tic": 300.0, "bpc": None, "scan_id": 2},
        ],
        "ms2": [{"rt": 1.1, "tic": 25.0, "scan_id": 2}],
        "tic_source": "mzml_tic_chromatogram_ms_level_array",
        "ms1_spectra_available": False,
        "ms1_spectrum_count": 0,
        "ms1_tic_point_count": 2,
        "bpc_available": False,
        "bpc_source": "unavailable",
        "bpc_scope": (
            "unavailable_tic_only_without_ms1_peak_arrays_or_bpc_chromatogram"
        ),
    }
    monkeypatch.setattr(
        "utils.mzml_io.extract_tic_overview",
        lambda path: enriched if path == source_path else None,
    )

    result = server._load_tic_overview_with_bpc_backfill(
        "task-source-tic",
        {"file_path": str(source_path)},
        artifact_path,
    )

    assert result["chromatogram_schema_version"] == 3
    assert result["tic_chromatogram_backfilled"] is True
    assert result["bpc_available"] is False
    assert result["bpc_backfill_attempted"] is True
    assert len(result["ms1"]) == 2
    persisted = json.loads(artifact_path.read_text(encoding="utf-8"))
    assert persisted["tic_source"] == "mzml_tic_chromatogram_ms_level_array"


def test_frontend_rejects_null_bpc_instead_of_coercing_it_to_zero():
    source = Path("web/static/js/report-generator.js").read_text(encoding="utf-8")

    assert "point?.[valueKey] !== null" in source
    assert "TIC 不能反推出真实 BPC" in source
    assert "系统不会用 MS2 基峰冒充真实 BPC" in source
