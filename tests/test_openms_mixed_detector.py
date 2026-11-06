from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

from utils import mzml_io, openms_backend


class _FakeMzMLReader:
    def __init__(self, spectra):
        self.spectra = spectra

    def __enter__(self):
        return iter(self.spectra)

    def __exit__(self, exc_type, exc, traceback):
        return False


def test_representation_audit_keeps_openms_mislabeled_dad_out_of_ms1(
    monkeypatch,
    tmp_path,
):
    spectra = [
        {
            "id": "merged=2 row=0",
            "ms level": 1,
            "profile spectrum": "",
            "electromagnetic radiation spectrum": "",
        },
        {
            "id": "merged=3 row=1",
            "ms level": 1,
            "profile spectrum": "",
            "wavelength array": [190.0, 192.0],
        },
        {"id": "scan=1", "ms level": 1, "centroid spectrum": ""},
        {"id": "scan=2", "ms level": 2, "centroid spectrum": ""},
    ]
    monkeypatch.setattr(
        "pyteomics.mzml.read",
        lambda _path: _FakeMzMLReader(spectra),
    )
    path = tmp_path / "mixed-detector.mzML"

    audit = openms_backend.inspect_mzml_representations(path)

    assert audit["counts"]["ms1"] == {
        "centroid": 1,
        "profile": 0,
        "unknown": 0,
        "spectra": 1,
    }
    assert audit["counts"]["ms2"]["centroid"] == 1
    assert audit["counts"]["auxiliary"] == {
        "centroid": 0,
        "profile": 2,
        "unknown": 0,
        "spectra": 2,
    }
    assert audit["mass_spectra"] == 2
    assert audit["profile_mass_spectra"] == 0
    assert audit["unknown_mass_spectra"] == 0
    assert audit["examples"] == {"profile": [], "unknown": []}


def test_tic_overview_excludes_optical_spectra_with_synthesized_ms1(
    monkeypatch,
):
    spectra = [
        {
            "id": "scan=1",
            "ms level": 1,
            "scanList": {"scan": [{"scan start time": 1.0}]},
            "total ion current": 100.0,
            "base peak intensity": 50.0,
            "base peak m/z": 123.4,
        },
        {
            "id": "scan=99",
            "ms level": 1,
            "profile spectrum": "",
            "electromagnetic radiation spectrum": "",
            "wavelength array": [190.0, 192.0],
            "scanList": {"scan": [{"scan start time": 1.1}]},
            "total ion current": 999999.0,
            "base peak intensity": 999999.0,
        },
        {
            "id": "scan=2",
            "ms level": 2,
            "scanList": {"scan": [{"scan start time": 1.2}]},
            "total ion current": 25.0,
            "base peak intensity": 12.0,
        },
    ]
    monkeypatch.setattr(
        "pyteomics.mzml.read",
        lambda _path, **_kwargs: _FakeMzMLReader(spectra),
    )

    overview = mzml_io.extract_tic_overview("mixed-detector.mzML")

    assert [point["scan_id"] for point in overview["ms1"]] == [1]
    assert [point["scan_id"] for point in overview["ms2"]] == [2]
    assert overview["ms1"][0]["tic"] == 100.0
    assert all(point["scan_id"] != 99 for point in overview["ms1"])


def test_all_centroid_branch_reports_auxiliary_spectra_at_top_level(
    monkeypatch,
    tmp_path,
):
    audit = {
        "counts": {
            "ms1": {"centroid": 1, "profile": 0, "unknown": 0, "spectra": 1},
            "ms2": {"centroid": 1, "profile": 0, "unknown": 0, "spectra": 1},
            "auxiliary": {"centroid": 0, "profile": 1, "unknown": 0, "spectra": 1},
        },
        "mass_spectra": 2,
        "profile_mass_spectra": 0,
        "unknown_mass_spectra": 0,
    }
    monkeypatch.setattr(openms_backend, "inspect_mzml_representations", lambda _path: audit)

    source, payload = openms_backend.prepare_mzml_for_analysis(
        tmp_path / "already-centroided.mzML",
        tmp_path / "artifacts",
    )

    assert source == (tmp_path / "already-centroided.mzML").resolve()
    assert payload["status"] == "not_required_all_mass_spectra_centroid"
    assert payload["auxiliary_spectra_excluded_from_output"] == 1


def _load_openms_runner_with_fake_binding(monkeypatch):
    class _SpectrumType:
        PROFILE = 1
        CENTROID = 2
        UNKNOWN = 3

    fake_openms = types.SimpleNamespace(
        __version__="test",
        SpectrumSettings=types.SimpleNamespace(SpectrumType=_SpectrumType),
    )
    monkeypatch.setitem(sys.modules, "pyopenms", fake_openms)
    path = Path(__file__).resolve().parents[1] / "scripts" / "openms_untargeted.py"
    spec = importlib.util.spec_from_file_location(
        "_test_openms_untargeted_mixed_detector", path
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_peak_pick_output_retains_only_source_explicit_ms1_ms2(
    monkeypatch,
    tmp_path,
):
    runner = _load_openms_runner_with_fake_binding(monkeypatch)

    class FakeSpectrum:
        def __init__(self, native_id, ms_level):
            self.native_id = native_id
            self.ms_level = ms_level

        def getNativeID(self):
            return self.native_id

        def getMSLevel(self):
            return self.ms_level

    class FakeExperiment:
        def __init__(self):
            self.spectra = [
                FakeSpectrum("scan=1", 1),
                FakeSpectrum("scan=2", 2),
                # OpenMS synthesized MS1 for an original DAD/UV spectrum.
                FakeSpectrum("merged=2 row=0", 1),
            ]

        def __iter__(self):
            return iter(self.spectra)

        def setSpectra(self, spectra):
            self.spectra = list(spectra)

    class FakeMzMLFile:
        stored_ids = None

        def load(self, _path, _experiment):
            return None

        def store(self, _path, experiment):
            type(self).stored_ids = [
                spectrum.getNativeID() for spectrum in experiment
            ]

    runner.oms.MSExperiment = FakeExperiment
    runner.oms.MzMLFile = FakeMzMLFile
    runner._explicit_mzml_ms_level_ids = lambda _path: (
        {1: {"scan=1"}, 2: {"scan=2"}},
        1,
    )

    def fake_pick(experiment, **_kwargs):
        assert [spectrum.getNativeID() for spectrum in experiment] == [
            "scan=1",
            "scan=2",
        ]
        return {
            "1": {"unknown": 0, "profile_picked": 1},
            "2": {"unknown": 0, "profile_picked": 1},
        }

    runner._pick_profile_spectra = fake_pick
    args = types.SimpleNamespace(
        input=str(tmp_path / "mixed.mzML"),
        output=str(tmp_path / "centroided.mzML"),
        strategy="hires",
        signal_to_noise=1.0,
    )

    payload = runner.command_peak_pick(args)

    assert FakeMzMLFile.stored_ids == ["scan=1", "scan=2"]
    assert payload["auxiliary_spectra_excluded_from_output"] == 1
    assert payload["analysis_output_filter"] == {
        "policy": "retain_only_source_explicit_ms1_ms2_native_ids",
        "retained_spectra": 2,
        "retained_by_ms_level": {"1": 1, "2": 1},
        "excluded_auxiliary_spectra": 1,
        "missing_explicit_mass_spectra": 0,
    }


def test_source_contract_rejects_optical_spectrum_even_with_synthesized_ms1(
    monkeypatch,
    tmp_path,
):
    runner = _load_openms_runner_with_fake_binding(monkeypatch)
    path = tmp_path / "openms-roundtrip.mzML"
    path.write_text(
        """<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run><spectrumList count="3">
    <spectrum id="scan=1" index="0" defaultArrayLength="0">
      <cvParam accession="MS:1000511" name="ms level" value="1"/>
    </spectrum>
    <spectrum id="scan=2" index="1" defaultArrayLength="0">
      <cvParam accession="MS:1000511" name="ms level" value="2"/>
    </spectrum>
    <spectrum id="dad=1" index="2" defaultArrayLength="0">
      <cvParam accession="MS:1000511" name="ms level" value="1"/>
      <cvParam accession="MS:1000804"
               name="electromagnetic radiation spectrum"/>
    </spectrum>
  </spectrumList></run>
</mzML>
""",
        encoding="utf-8",
    )

    levels, auxiliary = runner._explicit_mzml_ms_level_ids(path)

    assert levels == {1: {"scan=1"}, 2: {"scan=2"}}
    assert auxiliary == 1
