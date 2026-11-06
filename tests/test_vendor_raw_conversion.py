from __future__ import annotations

import asyncio
import hashlib
import json
import platform
import threading
from copy import deepcopy
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from utils.vendor_raw_conversion import (
    ConverterUnavailableError,
    VendorLicenseNotAcceptedError,
    VendorRawValidationError,
    build_docker_msconvert_command,
    build_local_msconvert_command,
    inspect_agilent_d,
    inspect_thermo_raw,
    validate_mzml,
)


def _section(identifier: str, value: str, name: str = "") -> str:
    return (
        "<SectionInfo>"
        f"<Name>{name}</Name>"
        f"<ID>{identifier}</ID>"
        f"<Value>{value}</Value>"
        "</SectionInfo>"
    )


def _make_agilent_d(
    root: Path,
    name: str = "sample_NEG.d",
    storage_mode: str = "Centroid",
) -> Path:
    raw = root / name
    acq = raw / "AcqData"
    acq.mkdir(parents=True)
    (acq / "MSScan.bin").write_bytes(b"scan")
    (acq / "MSPeak.bin").write_bytes(b"peak")
    rc_report = (
        "<DSRdlReport><Device><Name>Q-TOF</Name></Device></DSRdlReport>"
    )
    scic_report = (
        "<SectionReport>"
        + _section("ionPolarity", "Negative", "Ion Polarity")
        + _section("MS1MinRange", "100", "MS Min Range (m/z)")
        + _section("MS1MaxRange", "1700", "MS Max Range (m/z)")
        + _section("MS1ScanRate", "3.00", "MS Scan Rate")
        + _section("MS2MinRange", "100", "MS/MS Min Range (m/z)")
        + _section("MS2MaxRange", "1700", "MS/MS Max Range (m/z)")
        + _section("MS2ScanRate", "3.00", "MS/MS Scan Rate")
        + _section("IsolationWidthMS2", "Medium (~4 amu)", "Isolation Width MS/MS")
        + _section(
            "TS_1_AcqAutoMS2_CollisionEnergyFixed",
            "30.00,40.00",
            "Use Fixed Collision Energies",
        )
        + "<SectionInfo><Name>Acquisition Mode AutoMS2</Name>"
        "<ID>TS_1_AcqAutoMS2</ID></SectionInfo>"
        + _section("StorageMode", storage_mode, "Storage Mode")
        + "</SectionReport>"
    )
    outer = (
        "<DataFileMethodReport><MethodReport>"
        "<MethodName>synthetic_neg.m</MethodName>"
        f"<RCDevicesXml>{rc_report.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')}</RCDevicesXml>"
        f"<SCICDevicesXml>{scic_report.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')}</SCICDevicesXml>"
        "</MethodReport></DataFileMethodReport>"
    )
    (acq / "AcqMethod.xml").write_text(outer, encoding="utf-8")
    return raw


def _make_thermo_raw(
    root: Path,
    name: str = "Standard-MIX.raw",
    payload: bytes = b"synthetic thermo raw",
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    raw = root / name
    raw.write_bytes(payload)
    return raw


def _patch_local_conversion_runtime(monkeypatch, tmp_path, server):
    status = {
        "local_msconvert": {"available": False, "path": None},
        "docker": {
            "path": "/usr/local/bin/docker",
            "cli_available": True,
            "daemon_available": True,
            "image": "proteowizard/test:pinned",
            "image_available": True,
            "image_identity": "sha256:test-image",
            "wine_runner": "wine64_anyuser",
            "apple_silicon_pinned_image": True,
            "wine_runner_verified": True,
        },
        "recommended_backend": "docker",
    }
    monkeypatch.setattr(
        server,
        "VENDOR_CONVERSION_CACHE_ROOT",
        tmp_path / "vendor_conversion_cache",
    )
    monkeypatch.setattr(
        server,
        "VENDOR_CONVERSION_LEDGER_FILE",
        tmp_path / "vendor_conversion_ledger.json",
    )
    monkeypatch.setattr(
        server,
        "_cached_or_prepare_vendor_converter_status",
        lambda **kwargs: deepcopy(status),
    )
    return status


def _write_mzml(
    path: Path,
    *,
    include_ms1: bool = True,
    include_auxiliary: bool = False,
    include_offsets: bool = False,
    source_sha1: str = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
) -> None:
    spectra = []
    if include_auxiliary:
        spectra.append(
            """
            <spectrum index="0" id="dad=1" defaultArrayLength="2">
              <cvParam accession="MS:1000128" name="profile spectrum" value=""/>
            </spectrum>
            """
        )
    if include_ms1:
        spectra.append(
            """
            <spectrum index="0" id="scanId=1" defaultArrayLength="2">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <cvParam accession="MS:1000127" name="centroid spectrum" value=""/>
              <cvParam accession="MS:1000129" name="negative scan" value=""/>
            </spectrum>
            """
        )
    offset_xml = (
        '<cvParam accession="MS:1000828" name="isolation window lower offset" value="2"/>'
        '<cvParam accession="MS:1000829" name="isolation window upper offset" value="2"/>'
        if include_offsets else ""
    )
    spectra.append(
        f"""
        <spectrum index="{len(spectra)}" id="scanId=2" defaultArrayLength="3">
          <cvParam accession="MS:1000511" name="ms level" value="2"/>
          <cvParam accession="MS:1000127" name="centroid spectrum" value=""/>
          <cvParam accession="MS:1000129" name="negative scan" value=""/>
          <precursor spectrumRef="scanId=1">
            <isolationWindow>
              <cvParam accession="MS:1000827" name="isolation window target m/z" value="200"/>
              {offset_xml}
            </isolationWindow>
          </precursor>
        </spectrum>
        """
    )
    xml = f"""<?xml version="1.0" encoding="utf-8"?>
    <indexedmzML xmlns="http://psi.hupo.org/ms/mzml">
      <mzML>
        <fileDescription><sourceFileList count="1">
          <sourceFile id="MSScan.bin" name="MSScan.bin">
            <cvParam accession="MS:1000569" name="SHA-1" value="{source_sha1}"/>
          </sourceFile>
        </sourceFileList></fileDescription>
        <softwareList count="1"><software id="pwiz" version="test"/></softwareList>
        <run><spectrumList count="{len(spectra)}">{''.join(spectra)}</spectrumList></run>
      </mzML>
    </indexedmzML>
    """
    path.write_text(xml, encoding="utf-8")


def test_inspect_agilent_d_recovers_ms1_and_isolation_width(tmp_path):
    report = inspect_agilent_d(_make_agilent_d(tmp_path))
    acquisition = report["acquisition"]
    assert report["format"] == "Agilent MassHunter .d"
    assert acquisition["instrument"] == "Q-TOF"
    assert acquisition["acquisition_mode"] == "AutoMS2"
    assert acquisition["polarities"] == ["NEG"]
    assert acquisition["ms1_expected"] is True
    assert acquisition["ms2_expected"] is True
    assert acquisition["ms1_scan_rates_hz"] == [3.0]
    assert acquisition["isolation_width_da"] == 4.0
    assert acquisition["collision_energies"] == [30.0, 40.0]
    assert acquisition["centroid_acquisition_present"] is True
    assert acquisition["profile_acquisition_present"] is False
    scan_signature = report["required_file_signatures"]["AcqData/MSScan.bin"]
    assert scan_signature["sha256"] == hashlib.sha256(b"scan").hexdigest()
    assert scan_signature["sha1_for_mzml_provenance"] == hashlib.sha1(
        b"scan", usedforsecurity=False
    ).hexdigest()
    assert scan_signature["size_bytes"] == 4
    assert isinstance(scan_signature["mtime_ns"], int)


def test_inspect_agilent_both_storage_mode_reports_mixed_acquisition(tmp_path):
    report = inspect_agilent_d(
        _make_agilent_d(tmp_path, storage_mode="Both")
    )
    acquisition = report["acquisition"]

    assert acquisition["storage_modes"] == ["Both"]
    assert acquisition["centroid_acquisition_present"] is True
    assert acquisition["profile_acquisition_present"] is True


def test_inspect_thermo_raw_fingerprints_file_without_guessing_binary_metadata(
    tmp_path,
):
    raw = _make_thermo_raw(tmp_path / "neg", payload=b"thermo-qe")

    report = inspect_thermo_raw(raw, expected_polarities=["NEG"])

    assert report["source_kind"] == "thermo_raw"
    assert report["format"] == "Thermo Fisher RAW"
    assert report["file_count"] == 1
    assert report["total_bytes"] == len(b"thermo-qe")
    assert report["source_file_signature"]["sha256"] == hashlib.sha256(
        b"thermo-qe"
    ).hexdigest()
    assert report["source_file_signature"][
        "sha1_for_mzml_provenance"
    ] == hashlib.sha1(b"thermo-qe", usedforsecurity=False).hexdigest()
    assert report["acquisition"]["polarities"] == ["NEG"]
    assert report["acquisition"]["instrument"] == ""
    assert report["acquisition"]["ms1_expected"] is None


def test_inspect_thermo_raw_rejects_directory_or_wrong_suffix(tmp_path):
    directory = tmp_path / "not-a-file.raw"
    directory.mkdir()
    with pytest.raises(VendorRawValidationError, match="must be one .raw file"):
        inspect_thermo_raw(directory)

    wrong = tmp_path / "sample.bin"
    wrong.write_bytes(b"raw")
    with pytest.raises(VendorRawValidationError, match=".raw file suffix"):
        inspect_thermo_raw(wrong)


def test_validate_mzml_counts_ms1_ms2_and_audits_method_fallback(tmp_path):
    path = tmp_path / "converted.mzML"
    _write_mzml(path)
    report = validate_mzml(
        path,
        expected_ms1=True,
        expected_polarities=["NEG"],
        isolation_width_fallback_da=4.0,
    )
    assert report["status"] == "valid"
    assert report["actual_spectrum_count"] == 2
    assert report["mass_spectrum_count"] == 2
    assert report["auxiliary_spectrum_count"] == 0
    assert report["ms1_count"] == 1
    assert report["ms2_count"] == 1
    assert report["precursor_spectrum_ref_count"] == 1
    assert report["isolation_offsets_complete"] is False
    assert report["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert any("AcqMethod.xml (4 Da)" in item for item in report["warnings"])


def test_validate_mzml_separates_auxiliary_detector_spectra(tmp_path):
    path = tmp_path / "with_dad.mzML"
    _write_mzml(path, include_auxiliary=True)

    report = validate_mzml(path, expected_ms1=True)

    assert report["actual_spectrum_count"] == 3
    assert report["mass_spectrum_count"] == 2
    assert report["auxiliary_spectrum_count"] == 1
    assert report["nonempty_auxiliary_spectrum_count"] == 1
    assert any("auxiliary detector spectra" in item for item in report["warnings"])


def test_validate_mzml_extracts_qe_instrument_and_acquisition_cv_terms(tmp_path):
    path = tmp_path / "qe.mzML"
    source_sha1 = "1234567890abcdef1234567890abcdef12345678"
    path.write_text(
        f'''<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <fileDescription><sourceFileList count="1"><sourceFile id="RAW1" name="sample.raw">
    <cvParam accession="MS:1000563" name="Thermo RAW format" value=""/>
    <cvParam accession="MS:1000569" name="SHA-1" value="{source_sha1}"/>
  </sourceFile></sourceFileList></fileDescription>
  <referenceableParamGroupList count="1"><referenceableParamGroup id="CommonInstrumentParams">
    <cvParam accession="MS:1001911" name="Q Exactive" value=""/>
    <cvParam accession="MS:1000529" name="instrument serial number" value="QE-490"/>
  </referenceableParamGroup></referenceableParamGroupList>
  <softwareList count="1"><software id="pwiz" version="3.0.25142"/></softwareList>
  <instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">
    <referenceableParamGroupRef ref="CommonInstrumentParams"/>
    <componentList count="4">
      <source order="1"><cvParam accession="MS:1000073" name="electrospray ionization" value=""/></source>
      <analyzer order="2"><cvParam accession="MS:1000081" name="quadrupole" value=""/></analyzer>
      <analyzer order="3"><cvParam accession="MS:1000484" name="orbitrap" value=""/></analyzer>
      <detector order="4"><cvParam accession="MS:1000624" name="inductive detector" value=""/></detector>
    </componentList>
  </instrumentConfiguration></instrumentConfigurationList>
  <run id="qe-run" startTimeStamp="2023-03-27T10:14:49Z" defaultInstrumentConfigurationRef="IC1">
    <spectrumList count="2">
      <spectrum index="0" id="scan=1" defaultArrayLength="2">
        <cvParam accession="MS:1000511" name="ms level" value="1"/>
        <cvParam accession="MS:1000128" name="profile spectrum" value=""/>
        <cvParam accession="MS:1000129" name="negative scan" value=""/>
        <scanList count="1"><scan>
          <cvParam accession="MS:1000800" name="mass resolving power" value="70000"/>
          <cvParam accession="MS:1000512" name="filter string" value="FTMS - p ESI Full ms"/>
        </scan></scanList>
      </spectrum>
      <spectrum index="1" id="scan=2" defaultArrayLength="2">
        <cvParam accession="MS:1000511" name="ms level" value="2"/>
        <cvParam accession="MS:1000128" name="profile spectrum" value=""/>
        <cvParam accession="MS:1000129" name="negative scan" value=""/>
        <scanList count="1"><scan><cvParam accession="MS:1000800" name="mass resolving power" value="17500"/></scan></scanList>
        <precursor spectrumRef="scan=1"><isolationWindow>
          <cvParam accession="MS:1000827" name="isolation window target m/z" value="200"/>
          <cvParam accession="MS:1000828" name="isolation window lower offset" value="0.75"/>
          <cvParam accession="MS:1000829" name="isolation window upper offset" value="0.75"/>
        </isolationWindow><activation><cvParam accession="MS:1000045" name="collision energy" value="30"/></activation></precursor>
      </spectrum>
    </spectrumList>
  </run>
</mzML>''',
        encoding="utf-8",
    )

    report = validate_mzml(
        path,
        expected_ms1=True,
        expected_polarities=["NEG"],
        expected_source_sha1={"sample.raw": source_sha1},
    )

    assert report["instrument"]["models"] == [{
        "accession": "MS:1001911",
        "name": "Q Exactive",
        "value": "",
    }]
    assert report["instrument"]["serial_numbers"] == ["QE-490"]
    assert [row["name"] for row in report["instrument"]["mass_analyzers"]] == [
        "quadrupole",
        "orbitrap",
    ]
    assert report["instrument"]["run"]["startTimeStamp"] == (
        "2023-03-27T10:14:49Z"
    )
    assert report["acquisition"]["resolving_powers_by_ms_level"] == {
        "1": [70000.0],
        "2": [17500.0],
    }
    assert report["acquisition"]["collision_energies_by_ms_level"] == {
        "2": [30.0]
    }
    assert report["acquisition"]["isolation_widths_da"] == [1.5]
    assert report["source_checksum_matches"] == {"sample.raw": True}


def test_validate_mzml_rejects_conversion_that_drops_expected_ms1(tmp_path):
    path = tmp_path / "ms2_only.mzML"
    _write_mzml(path, include_ms1=False)
    with pytest.raises(VendorRawValidationError, match="produced 0 MS1"):
        validate_mzml(path, expected_ms1=True)


def test_validate_mzml_rejects_stale_output_from_another_raw_directory(tmp_path):
    path = tmp_path / "stale.mzML"
    _write_mzml(path)
    with pytest.raises(VendorRawValidationError, match="checksum mismatch"):
        validate_mzml(
            path,
            expected_source_sha1={
                "MSScan.bin": "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
            },
        )


def test_cli_pairs_existing_mzml_with_raw_method_and_source(capsys, tmp_path):
    from scripts.convert_vendor_raw import main

    raw = _make_agilent_d(tmp_path)
    path = tmp_path / "paired.mzML"
    _write_mzml(path, source_sha1=hashlib.sha1(b"scan").hexdigest())

    exit_code = main([
        "--validate-mzml",
        str(path),
        "--raw-dir",
        str(raw),
    ])

    assert exit_code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["ms1_count"] == 1
    assert payload["ms2_count"] == 1
    assert payload["source_checksum_matches"]["MSScan.bin"] is True
    assert payload["isolation_width_fallback_da"] == 4.0
    assert payload["paired_raw"]["acquisition"]["ms1_expected"] is True
    assert any(
        "AcqMethod.xml (4 Da)" in warning
        for warning in payload["warnings"]
    )


def test_msconvert_commands_never_add_mslevel_filter(tmp_path):
    raw = tmp_path / "sample.d"
    raw.mkdir()
    output = tmp_path / "out"
    output.mkdir()
    local = build_local_msconvert_command(
        "/opt/pwiz/msconvert", raw, output, "sample.mzML"
    )
    docker = build_docker_msconvert_command(
        "docker", "proteowizard/test:latest", raw, output, "sample.mzML"
    )
    assert "--mzML" in local and "--64" in local and "--zlib" in local
    assert "--mzML" in docker and "--64" in docker and "--zlib" in docker
    assert not any("msLevel" in token for token in local + docker)
    assert "--filter" not in local
    assert "--filter" not in docker
    expected_runner = (
        "wine64_anyuser"
        if platform.machine().lower() in {"arm64", "aarch64"}
        else "wine"
    )
    assert docker[docker.index("proteowizard/test:latest") + 1] == expected_runner


def test_agilent_auto_peak_picking_centroids_both_and_preserves_pure_centroid(
    monkeypatch, tmp_path
):
    from utils import vendor_raw_conversion

    status = {
        "local_msconvert": {
            "path": "/opt/pwiz/msconvert",
            "available": True,
        },
        "docker": {
            "path": "",
            "cli_available": False,
            "daemon_available": False,
            "image_available": False,
        },
    }
    commands = []
    monkeypatch.setattr(
        vendor_raw_conversion,
        "converter_status",
        lambda **kwargs: status,
    )
    monkeypatch.setattr(
        vendor_raw_conversion,
        "_run_probe",
        lambda *args, **kwargs: type(
            "Completed",
            (),
            {"returncode": 0, "stdout": "msconvert test", "stderr": ""},
        )(),
    )

    def fake_run(command, **kwargs):
        commands.append(command)
        output_dir = Path(command[command.index("--outdir") + 1])
        output_name = command[command.index("--outfile") + 1]
        _write_mzml(
            output_dir / output_name,
            source_sha1=hashlib.sha1(
                b"scan", usedforsecurity=False
            ).hexdigest(),
        )
        return type(
            "Completed",
            (),
            {"returncode": 0, "stdout": "converted", "stderr": ""},
        )()

    monkeypatch.setattr(vendor_raw_conversion.subprocess, "run", fake_run)

    both = _make_agilent_d(
        tmp_path / "both", name="both.d", storage_mode="Both"
    )
    both_manifest = vendor_raw_conversion.convert_agilent_d(
        both,
        tmp_path / "both.mzML",
        backend="local",
        peak_picking_strategy="auto",
    )
    both_command = commands[-1]
    assert both_manifest["conversion_policy"][
        "peak_picking_strategy_resolved"
    ] == "vendor"
    assert both_command[both_command.index("--filter") + 1] == (
        "peakPicking true 1-"
    )

    centroid = _make_agilent_d(
        tmp_path / "centroid",
        name="centroid.d",
        storage_mode="Centroid",
    )
    centroid_manifest = vendor_raw_conversion.convert_agilent_d(
        centroid,
        tmp_path / "centroid.mzML",
        backend="local",
        peak_picking_strategy="auto",
    )
    centroid_command = commands[-1]
    assert centroid_manifest["conversion_policy"][
        "peak_picking_strategy_resolved"
    ] == "preserve"
    assert "--filter" not in centroid_command


def test_prepare_converter_requires_explicit_license_and_pull_permission(
    monkeypatch,
):
    from utils import vendor_raw_conversion

    status = {
        "local_msconvert": {"available": False},
        "docker": {
            "path": "/usr/local/bin/docker",
            "cli_available": True,
            "daemon_available": True,
            "image_available": False,
        },
    }
    monkeypatch.setattr(
        vendor_raw_conversion,
        "converter_status",
        lambda **kwargs: status,
    )

    with pytest.raises(VendorLicenseNotAcceptedError):
        vendor_raw_conversion.prepare_converter(
            accept_vendor_licenses=False,
            allow_container_pull=True,
        )
    with pytest.raises(ConverterUnavailableError, match="one-time"):
        vendor_raw_conversion.prepare_converter(
            accept_vendor_licenses=True,
            allow_container_pull=False,
        )


def test_prepare_converter_pulls_pinned_image_once(monkeypatch):
    from utils import vendor_raw_conversion

    before = {
        "local_msconvert": {"available": False},
        "docker": {
            "path": "/usr/local/bin/docker",
            "cli_available": True,
            "daemon_available": True,
            "image_available": False,
        },
    }
    after = {
        "local_msconvert": {"available": False},
        "docker": {
            "path": "/usr/local/bin/docker",
            "cli_available": True,
            "daemon_available": True,
            "image_available": True,
        },
    }
    statuses = iter([before, after])
    commands = []
    monkeypatch.setattr(
        vendor_raw_conversion,
        "converter_status",
        lambda **kwargs: next(statuses),
    )
    monkeypatch.setattr(vendor_raw_conversion, "_arm64_host", lambda: True)

    def fake_run(command, **kwargs):
        commands.append(command)
        return type(
            "Completed",
            (),
            {"returncode": 0, "stdout": "downloaded", "stderr": ""},
        )()

    monkeypatch.setattr(vendor_raw_conversion.subprocess, "run", fake_run)
    image = "proteowizard/test:pinned"

    result = vendor_raw_conversion.prepare_converter(
        accept_vendor_licenses=True,
        allow_container_pull=True,
        docker_image=image,
    )

    assert result["ready"] is True
    assert result["downloaded"] is True
    assert commands[0] == [
        "/usr/local/bin/docker",
        "pull",
        "--platform",
        "linux/amd64",
        image,
    ]
    assert commands[1][:6] == [
        "/usr/local/bin/docker",
        "run",
        "--rm",
        "--platform",
        "linux/amd64",
        image,
    ]
    assert "command -v wine64_anyuser" in commands[1][-1]


def test_conversion_never_pulls_a_missing_image(monkeypatch, tmp_path):
    from utils import vendor_raw_conversion

    raw = _make_agilent_d(tmp_path)
    status = {
        "local_msconvert": {"path": "", "available": False},
        "docker": {
            "path": "/usr/local/bin/docker",
            "cli_available": True,
            "daemon_available": True,
            "image_available": False,
        },
    }
    monkeypatch.setattr(
        vendor_raw_conversion,
        "converter_status",
        lambda **kwargs: status,
    )

    def unexpected_run(*args, **kwargs):
        pytest.fail("conversion attempted a subprocess before local image check")

    monkeypatch.setattr(vendor_raw_conversion.subprocess, "run", unexpected_run)
    with pytest.raises(ConverterUnavailableError, match="never downloads images"):
        vendor_raw_conversion.convert_agilent_d(
            raw,
            tmp_path / "out.mzML",
            backend="docker",
            accept_vendor_licenses=True,
            allow_container_pull=True,
        )


def test_docker_command_can_force_audited_wine64_runner(tmp_path):
    raw = tmp_path / "sample.d"
    raw.mkdir()
    output = tmp_path / "out"
    output.mkdir()

    command = build_docker_msconvert_command(
        "docker",
        "proteowizard/test:pinned",
        raw,
        output,
        "sample.mzML",
        wine_runner="wine64_anyuser",
    )

    assert command[command.index("proteowizard/test:pinned") + 1] == (
        "wine64_anyuser"
    )


def test_docker_command_supports_thermo_raw_and_vendor_filter_is_first(tmp_path):
    raw = _make_thermo_raw(tmp_path / "neg")
    output = tmp_path / "out"
    output.mkdir()

    command = build_docker_msconvert_command(
        "docker",
        "proteowizard/test:pinned",
        raw,
        output,
        "Standard-MIX_NEG.mzML",
        wine_runner="wine64_anyuser",
        vendor_peak_picking=True,
    )

    assert "/input/Standard-MIX.raw" in command
    assert f"{raw.parent}:/input:ro" in command
    filter_index = command.index("--filter")
    assert command[filter_index + 1] == "peakPicking true 1-"
    assert not any("msLevel" in token for token in command)
    assert command.index("--filter") < command.index("--outdir")


def test_raw_method_width_is_inherited_but_user_override_wins():
    from web.server import AnalysisParams, _analysis_params_with_task_provenance

    task = {
        "recommended_isolation_window_width_da": 4.0,
        "recommended_isolation_window_width_source": (
            "Agilent AcqMethod.xml / IsolationWidthMS2"
        ),
    }
    inherited = _analysis_params_with_task_provenance(task, AnalysisParams())
    assert inherited["isolation_window_width_da"] == 4.0
    assert inherited["isolation_window_width_source"].startswith("raw_method:")

    overridden = _analysis_params_with_task_provenance(
        task, AnalysisParams(isolation_window_width_da=1.2)
    )
    assert overridden["isolation_window_width_da"] == 1.2
    assert overridden["isolation_window_width_source"] == "user_full_width_override"


def test_vendor_folder_api_reconstructs_one_d_and_creates_task(
    monkeypatch, tmp_path
):
    from utils import vendor_raw_conversion
    from web import server

    upload_dir = tmp_path / "uploads"
    vendor_dir = tmp_path / "vendor"
    upload_dir.mkdir()
    vendor_dir.mkdir()
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "VENDOR_RAW_DIR", vendor_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks.json")
    monkeypatch.setattr(server, "_tasks", {})

    def fake_convert(raw_dir, output_path, **kwargs):
        raw_dir = Path(raw_dir)
        assert (raw_dir / "AcqData" / "MSScan.bin").read_bytes() == b"scan"
        assert (raw_dir / "AcqData" / "MSPeak.bin").read_bytes() == b"peak"
        assert kwargs["allow_container_pull"] is False
        output_path = Path(output_path)
        _write_mzml(output_path)
        output = validate_mzml(
            output_path,
            expected_ms1=True,
            expected_polarities=["NEG"],
            isolation_width_fallback_da=4.0,
        )
        return {
            "status": "completed",
            "source": {
                "source_path": str(raw_dir),
                "raw_name": raw_dir.name,
                "acquisition": {
                    "isolation_width_da": 4.0,
                    "isolation_width_source": (
                        "Agilent AcqMethod.xml / IsolationWidthMS2"
                    ),
                },
            },
            "output": output,
            "converter": {"backend": "test"},
        }

    monkeypatch.setattr(vendor_raw_conversion, "convert_agilent_d", fake_convert)
    relative = [
        "sample.d/AcqData/MSScan.bin",
        "sample.d/AcqData/MSPeak.bin",
        "sample.d/AcqData/AcqMethod.xml",
    ]
    files = [
        ("files", ("MSScan.bin", b"scan", "application/octet-stream")),
        ("files", ("MSPeak.bin", b"peak", "application/octet-stream")),
        ("files", ("AcqMethod.xml", b"<method/>", "application/xml")),
    ]
    with TestClient(server.app) as client:
        response = client.post(
            "/api/vendor/convert",
            files=files,
            data={
                "manifest": json.dumps(relative),
                "accept_vendor_licenses": "true",
                "allow_container_pull": "false",
                "backend": "auto",
            },
        )
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["uploaded"][0]["filename"] == "sample.mzML"
    assert body["uploaded"][0]["ms1_count"] == 1
    task_id = body["uploaded"][0]["task_id"]
    assert server._tasks[task_id]["recommended_isolation_window_width_da"] == 4.0
    assert Path(server._tasks[task_id]["file_path"]).is_file()
    assert list(vendor_dir.iterdir()) == []


def test_native_vendor_folder_api_selects_and_converts_without_browser_upload(
    monkeypatch, tmp_path
):
    from utils import vendor_raw_conversion
    from web import server

    raw = _make_agilent_d(tmp_path)
    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_vendor_local_selections", {})
    converter_status = _patch_local_conversion_runtime(
        monkeypatch, tmp_path, server
    )
    monkeypatch.setattr(
        server,
        "_choose_macos_vendor_directory",
        lambda: raw,
    )

    def fake_convert(raw_dir, output_path, **kwargs):
        raw_dir = Path(raw_dir)
        assert raw_dir == raw
        assert kwargs["accept_vendor_licenses"] is True
        assert kwargs["allow_container_pull"] is False
        output_path = Path(output_path)
        _write_mzml(output_path)
        output = validate_mzml(
            output_path,
            expected_ms1=True,
            expected_polarities=["NEG"],
            isolation_width_fallback_da=4.0,
        )
        source = deepcopy(kwargs["raw_report"])
        source["acquisition"].update({
            "isolation_width_da": 4.0,
            "isolation_width_source": (
                "Agilent AcqMethod.xml / IsolationWidthMS2"
            ),
        })
        _, converter_identity = server._vendor_converter_identity(
            converter_status,
            requested_backend=kwargs["backend"],
            accept_vendor_licenses=True,
        )
        return {
            "status": "completed",
            "source": source,
            "output": output,
            "converter": converter_identity,
            "conversion_policy": {
                "peak_picking_strategy_resolved": "preserve",
                "preserves_ms1_and_ms2": True,
            },
        }

    monkeypatch.setattr(vendor_raw_conversion, "convert_agilent_d", fake_convert)
    with TestClient(server.app) as client:
        selected = client.post("/api/vendor/select-local-folder", json={})
        assert selected.status_code == 200, selected.text
        selection = selected.json()
        assert selection["raw_name"] == raw.name
        assert selection["selection_count"] == 1
        assert selection["selections"][0]["raw_name"] == raw.name
        assert selection["browser_upload_required"] is False
        assert selection["spectrum_representation"] == "centroid"
        assert selection["recommended_peak_picking_strategy"] == "preserve"
        assert selection["acquisition"]["centroid_acquisition_present"] is True
        assert selection["acquisition"]["profile_acquisition_present"] is False

        converted = client.post(
            "/api/vendor/convert-local",
            json={
                "selection_token": selection["selection_token"],
                "accept_vendor_licenses": True,
                "allow_container_pull": False,
                "backend": "auto",
                "peak_picking_strategy": "auto",
            },
        )

    assert converted.status_code == 200, converted.text
    body = converted.json()
    assert body["selection_mode"] == "macos_native_folder_picker"
    assert body["uploaded"][0]["filename"] == "sample_NEG.mzML"
    task_id = body["uploaded"][0]["task_id"]
    task = server._tasks[task_id]
    assert task["source_raw_filename"] == raw.name
    assert task["vendor_conversion"]["source"]["source_path"] == raw.name
    assert task["vendor_local_selection"]["mode"] == (
        "macos_native_folder_picker"
    )
    assert selection["selection_token"] not in server._vendor_local_selections


def test_native_vendor_folder_api_accepts_parent_with_multiple_d_directories(
    monkeypatch,
    tmp_path,
):
    from web import server

    parent = tmp_path / "batch"
    first = _make_agilent_d(parent, "sample_01.d")
    second = _make_agilent_d(parent, "sample_02.d")
    monkeypatch.setattr(server, "_vendor_local_selections", {})
    monkeypatch.setattr(
        server,
        "_choose_macos_vendor_directory",
        lambda: parent,
    )

    with TestClient(server.app) as client:
        selected = client.post("/api/vendor/select-local-folder", json={})

    assert selected.status_code == 200, selected.text
    payload = selected.json()
    assert payload["selection_mode"] == (
        "macos_native_parent_folder_picker"
    )
    assert payload["selected_root_name"] == parent.name
    assert payload["selection_count"] == 2
    assert payload["invalid_directories"] == []
    assert [row["raw_name"] for row in payload["selections"]] == [
        first.name,
        second.name,
    ]
    tokens = {
        row["selection_token"] for row in payload["selections"]
    }
    assert len(tokens) == 2
    assert tokens == set(server._vendor_local_selections)


def test_native_thermo_batch_discovers_pos_neg_one_level_and_keeps_context(
    monkeypatch,
    tmp_path,
):
    from web import server

    root = tmp_path / "qe_blind"
    neg = _make_thermo_raw(root / "neg")
    pos = _make_thermo_raw(root / "pos")
    monkeypatch.setattr(server, "_vendor_local_selections", {})
    monkeypatch.setattr(server, "_choose_macos_thermo_directory", lambda: root)

    with TestClient(server.app) as client:
        selected = client.post("/api/vendor/select-local-thermo-folder", json={})

    assert selected.status_code == 200, selected.text
    payload = selected.json()
    assert payload["selection_mode"] == (
        "macos_native_thermo_parent_folder_picker"
    )
    assert payload["selection_count"] == 2
    rows = {row["relative_source"]: row for row in payload["selections"]}
    assert rows[f"neg/{neg.name}"]["expected_polarities"] == ["NEG"]
    assert rows[f"pos/{pos.name}"]["expected_polarities"] == ["POS"]
    assert all(row["source_kind"] == "thermo_raw" for row in rows.values())
    assert len(server._vendor_local_selections) == 2


def test_native_thermo_conversion_creates_polarity_safe_qe_task(
    monkeypatch,
    tmp_path,
):
    from utils import vendor_raw_conversion
    from web import server

    root = tmp_path / "qe_blind"
    raw = _make_thermo_raw(root / "neg")
    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_vendor_local_selections", {})
    converter_status = _patch_local_conversion_runtime(
        monkeypatch, tmp_path, server
    )
    monkeypatch.setattr(server, "_choose_macos_thermo_directory", lambda: root)

    def fake_convert(raw_file, output_path, **kwargs):
        assert Path(raw_file) == raw
        assert kwargs["expected_polarities"] == ["NEG"]
        assert kwargs["peak_picking_strategy"] == "auto"
        assert kwargs["allow_container_pull"] is False
        output_path = Path(output_path)
        _write_mzml(output_path)
        output = validate_mzml(
            output_path,
            expected_ms1=True,
            expected_polarities=["NEG"],
        )
        output["source_checksum_matches"] = {raw.name: True}
        output["instrument"] = {
            "models": [{"accession": "MS:1001911", "name": "Q Exactive", "value": ""}],
            "serial_numbers": ["QE-490"],
            "mass_analyzers": [{"accession": "MS:1000484", "name": "orbitrap", "value": ""}],
        }
        output["acquisition"] = {"isolation_widths_da": [1.5]}
        source = deepcopy(kwargs["raw_report"])
        _, converter_identity = server._vendor_converter_identity(
            converter_status,
            requested_backend=kwargs["backend"],
            accept_vendor_licenses=True,
        )
        return {
            "status": "completed",
            "source_kind": "thermo_raw",
            "format": "Thermo Fisher RAW to mzML",
            "source": source,
            "output": output,
            "converter": converter_identity,
            "conversion_policy": {
                "peak_picking_strategy_resolved": "vendor",
                "preserves_ms1_and_ms2": True,
            },
        }

    monkeypatch.setattr(vendor_raw_conversion, "convert_thermo_raw", fake_convert)
    with TestClient(server.app) as client:
        selection = client.post(
            "/api/vendor/select-local-thermo-folder", json={}
        ).json()
        converted = client.post(
            "/api/vendor/convert-local",
            json={
                "selection_token": selection["selection_token"],
                "accept_vendor_licenses": True,
                "allow_container_pull": False,
                "backend": "auto",
                "peak_picking_strategy": "auto",
            },
        )
        assert converted.status_code == 200, converted.text
        task_id = converted.json()["uploaded"][0]["task_id"]
        audit = client.get(f"/api/tasks/{task_id}/input-audit")

    task = server._tasks[task_id]
    assert task["filename"] == "Standard-MIX_NEG.mzML"
    assert task["source_raw_filename"] == raw.name
    assert task["vendor_local_selection"]["source_kind"] == "thermo_raw"
    assert task["vendor_local_selection"]["expected_polarities"] == ["NEG"]
    assert task["recommended_isolation_window_width_da"] == 1.5
    assert task["vendor_conversion"]["source"]["source_path"] == raw.name
    assert audit.status_code == 200, audit.text
    assert audit.json()["verdict"] == "verified"
    assert audit.json()["input_kind"] == "Thermo Fisher RAW"
    assert audit.json()["acquisition"]["instrument"] == "Q Exactive"


def test_native_thermo_conversion_can_register_explicit_blank_without_task(
    monkeypatch,
    tmp_path,
):
    from utils import vendor_raw_conversion
    from web import server

    root = tmp_path / "qe_blind"
    raw = _make_thermo_raw(root / "neg", name="Blank-meoh-2.raw")
    blank_dir = tmp_path / "blanks"
    upload_dir = tmp_path / "uploads"
    blank_dir.mkdir()
    upload_dir.mkdir()
    monkeypatch.setattr(server, "BLANK_DIR", blank_dir)
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_vendor_local_selections", {})
    converter_status = _patch_local_conversion_runtime(
        monkeypatch, tmp_path, server
    )
    monkeypatch.setattr(server, "_choose_macos_thermo_directory", lambda: root)

    def fake_convert(raw_file, output_path, **kwargs):
        assert Path(raw_file) == raw
        output_path = Path(output_path)
        assert output_path.parent == blank_dir
        output_path.write_bytes(b"converted blank")
        source = deepcopy(kwargs["raw_report"])
        _, converter_identity = server._vendor_converter_identity(
            converter_status,
            requested_backend=kwargs["backend"],
            accept_vendor_licenses=True,
        )
        return {
            "status": "completed",
            "source_kind": "thermo_raw",
            "source": source,
            "output": {
                "status": "valid",
                "size_bytes": output_path.stat().st_size,
                "sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
                "ms1_count": 10,
                "ms2_count": 5,
                "polarities": {"POS": 0, "NEG": 15},
                "source_checksum_matches": {raw.name: True},
                "instrument": {"models": [{"name": "Q Exactive"}]},
                "acquisition": {"isolation_widths_da": [1.5]},
                "warnings": [],
            },
            "converter": converter_identity,
            "conversion_policy": {
                "peak_picking_strategy_resolved": "vendor"
            },
        }

    def fake_register(path, display_name, **kwargs):
        assert Path(path).is_file()
        assert display_name == "Blank-meoh-2_NEG.mzML"
        assert kwargs["blank_role"] == "solvent_blank"
        assert kwargs["extra_metadata"]["vendor_local_selection"][
            "intake_role"
        ] == "blank"
        return {
            "blank_path": Path(path).name,
            "blank_filename": display_name,
            "polarities": ["NEG"],
            "blank_role": "solvent_blank",
            "n_ms1_scans": 10,
            "n_ms2_contexts": 5,
        }

    monkeypatch.setattr(vendor_raw_conversion, "convert_thermo_raw", fake_convert)
    monkeypatch.setattr(server, "_register_managed_blank_file", fake_register)
    with TestClient(server.app) as client:
        selection = client.post(
            "/api/vendor/select-local-thermo-folder", json={}
        ).json()
        converted = client.post(
            "/api/vendor/convert-local",
            json={
                "selection_token": selection["selection_token"],
                "accept_vendor_licenses": True,
                "backend": "auto",
                "intake_role": "blank",
                "blank_role": "solvent_blank",
                "peak_picking_strategy": "auto",
            },
        )

    assert converted.status_code == 200, converted.text
    payload = converted.json()
    assert payload["uploaded"] == []
    assert payload["blank_count"] == 1
    assert payload["blank_uploads"][0]["polarities"] == ["NEG"]
    assert server._tasks == {}
    assert Path(blank_dir / payload["blank_uploads"][0]["blank_path"]).is_file()


def test_native_vendor_operations_reject_cross_site_browser_requests():
    from web import server

    with TestClient(server.app) as client:
        response = client.post(
            "/api/vendor/converter/prepare",
            json={
                "accept_vendor_licenses": True,
                "allow_container_pull": True,
            },
            headers={"Origin": "https://untrusted.example"},
        )

    assert response.status_code == 403
    assert "Cross-site" in response.json()["detail"]


def test_native_vendor_conversion_reuses_verified_cache_with_independent_tasks(
    monkeypatch,
    tmp_path,
):
    from utils import vendor_raw_conversion
    from web import server

    raw = _make_agilent_d(tmp_path)
    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_vendor_local_selections", {})
    converter_status = _patch_local_conversion_runtime(
        monkeypatch, tmp_path, server
    )
    monkeypatch.setattr(server, "_choose_macos_vendor_directory", lambda: raw)
    converter_calls = 0

    def fake_convert(raw_dir, output_path, **kwargs):
        nonlocal converter_calls
        converter_calls += 1
        assert Path(raw_dir) == raw
        output_path = Path(output_path)
        _write_mzml(output_path)
        output = validate_mzml(
            output_path,
            expected_ms1=True,
            expected_polarities=["NEG"],
            isolation_width_fallback_da=4.0,
        )
        source = deepcopy(kwargs["raw_report"])
        _, converter_identity = server._vendor_converter_identity(
            converter_status,
            requested_backend=kwargs["backend"],
            accept_vendor_licenses=True,
        )
        return {
            "status": "completed",
            "source_kind": "agilent_d",
            "format": "Agilent MassHunter .d to mzML",
            "source": source,
            "output": output,
            "converter": converter_identity,
            "conversion_policy": {
                "peak_picking_strategy_requested": "auto",
                "peak_picking_strategy_resolved": "preserve",
                "preserves_ms1_and_ms2": True,
            },
        }

    monkeypatch.setattr(vendor_raw_conversion, "convert_agilent_d", fake_convert)
    responses = []
    with TestClient(server.app) as client:
        for _ in range(2):
            selection = client.post(
                "/api/vendor/select-local-folder", json={}
            ).json()
            response = client.post(
                "/api/vendor/convert-local",
                json={
                    "selection_token": selection["selection_token"],
                    "accept_vendor_licenses": True,
                    "allow_container_pull": False,
                    "backend": "auto",
                    "peak_picking_strategy": "auto",
                },
            )
            assert response.status_code == 200, response.text
            responses.append(response.json())

    assert converter_calls == 1
    assert responses[0]["conversion"]["cache"]["hit"] is False
    assert responses[0]["conversion"]["cache"]["status"] == "stored"
    assert responses[1]["conversion"]["cache"]["hit"] is True
    assert responses[1]["conversion"]["cache"]["status"] == "hit"
    assert (
        responses[0]["conversion"]["cache"]["cache_key"]
        == responses[1]["conversion"]["cache"]["cache_key"]
    )
    task_ids = [response["uploaded"][0]["task_id"] for response in responses]
    assert task_ids[0] != task_ids[1]
    task_paths = [Path(server._tasks[task_id]["file_path"]) for task_id in task_ids]
    assert all(path.is_file() and not path.is_symlink() for path in task_paths)
    assert task_paths[0].read_bytes() == task_paths[1].read_bytes()
    assert task_paths[0].stat().st_ino != task_paths[1].stat().st_ino
    assert responses[1]["conversion"]["cache"][
        "task_output_has_independent_inode"
    ] is True


def test_validated_vendor_registration_rejects_symlink_and_wrong_hash(
    monkeypatch,
    tmp_path,
):
    from web import server

    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    external = tmp_path / "outside.mzML"
    external.write_bytes(b"outside payload")
    linked = upload_dir / "linked.mzML"
    linked.symlink_to(external)

    def details(path: Path, digest: str) -> dict:
        return {
            "destination": str(path),
            "display_name": "linked.mzML",
            "raw_name": "source.d",
            "source_kind": "agilent_d",
            "intake_role": "sample",
            "blank_role": "background_blank",
            "selection_mode": "test",
            "selection_audit": {},
            "conversion": {
                "source": {"acquisition": {}},
                "output": {
                    "status": "valid",
                    "size_bytes": external.stat().st_size,
                    "sha256": digest,
                },
            },
        }

    with pytest.raises(ValueError, match="symlink"):
        server._register_validated_vendor_conversion(
            details(linked, hashlib.sha256(external.read_bytes()).hexdigest()),
            registration_token="registration-token",
        )

    regular = upload_dir / "regular.mzML"
    regular.write_bytes(external.read_bytes())
    with pytest.raises(ValueError, match="SHA-256 changed"):
        server._register_validated_vendor_conversion(
            details(regular, "0" * 64),
            registration_token="registration-token",
        )


def test_queued_ledger_selection_survives_picker_ttl_in_same_process(
    monkeypatch,
    tmp_path,
):
    from web import server

    raw = _make_agilent_d(tmp_path)
    token = "a" * 32
    monkeypatch.setattr(
        server,
        "VENDOR_CONVERSION_LEDGER_FILE",
        tmp_path / "vendor_conversion_ledger.json",
    )
    monkeypatch.setattr(
        server,
        "_vendor_local_selections",
        {
            token: {
                "source_path": str(raw),
                "raw_name": raw.name,
                "source_kind": "agilent_d",
                "expected_polarities": [],
                "selection_mode": "test-picker",
                "created_at": 1.0,
            }
        },
    )
    server._vendor_conversion_ledger().queue(
        token,
        metadata={
            "source_path": str(raw),
            "raw_name": raw.name,
            "source_kind": "agilent_d",
            "expected_polarities": [],
            "selection_mode": "test-picker",
            "selection_created_at_epoch": 1.0,
            "display_name": "sample_NEG.mzML",
            "destination": str(tmp_path / "sample_NEG.mzML"),
            "request_contract": {},
        },
    )

    resolved = server._resolve_local_vendor_selection(token)

    assert resolved == (token, raw.resolve(), "agilent_d", [], "test-picker")


def test_same_vendor_selection_token_has_one_atomic_conversion_owner(
    monkeypatch,
    tmp_path,
):
    from utils import vendor_raw_conversion
    from web import server

    raw = _make_agilent_d(tmp_path)
    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_vendor_local_selections", {})
    converter_status = _patch_local_conversion_runtime(
        monkeypatch, tmp_path, server
    )
    monkeypatch.setattr(server, "_choose_macos_vendor_directory", lambda: raw)
    converter_started = threading.Event()
    allow_converter_finish = threading.Event()
    converter_calls = 0

    def fake_convert(raw_dir, output_path, **kwargs):
        nonlocal converter_calls
        converter_calls += 1
        converter_started.set()
        assert allow_converter_finish.wait(timeout=5)
        output_path = Path(output_path)
        _write_mzml(output_path)
        output = validate_mzml(
            output_path,
            expected_ms1=True,
            expected_polarities=["NEG"],
            isolation_width_fallback_da=4.0,
        )
        _, converter_identity = server._vendor_converter_identity(
            converter_status,
            requested_backend=kwargs["backend"],
            accept_vendor_licenses=True,
        )
        return {
            "status": "completed",
            "source_kind": "agilent_d",
            "format": "Agilent MassHunter .d to mzML",
            "source": deepcopy(kwargs["raw_report"]),
            "output": output,
            "converter": converter_identity,
            "conversion_policy": {
                "peak_picking_strategy_requested": "auto",
                "peak_picking_strategy_resolved": "preserve",
                "preserves_ms1_and_ms2": True,
            },
        }

    monkeypatch.setattr(vendor_raw_conversion, "convert_agilent_d", fake_convert)
    with TestClient(server.app) as client:
        selection = client.post(
            "/api/vendor/select-local-folder", json={}
        ).json()
        request_payload = {
            "selection_token": selection["selection_token"],
            "accept_vendor_licenses": True,
            "allow_container_pull": False,
            "backend": "auto",
            "peak_picking_strategy": "auto",
        }
        with ThreadPoolExecutor(max_workers=2) as executor:
            first = executor.submit(
                client.post,
                "/api/vendor/convert-local",
                json=request_payload,
            )
            assert converter_started.wait(timeout=5)
            second = executor.submit(
                client.post,
                "/api/vendor/convert-local",
                json=request_payload,
            )
            second_response = second.result(timeout=5)
            allow_converter_finish.set()
            first_response = first.result(timeout=10)

    assert first_response.status_code == 200, first_response.text
    assert second_response.status_code == 409, second_response.text
    assert converter_calls == 1
    task_id = first_response.json()["uploaded"][0]["task_id"]
    task_path = Path(server._tasks[task_id]["file_path"])
    assert task_path.is_file()
    assert task_path.name.startswith(f"{selection['selection_token']}_")


def test_same_eight_character_token_prefix_keeps_distinct_vendor_outputs(
    monkeypatch,
    tmp_path,
):
    from utils import vendor_raw_conversion
    from web import server

    first_raw = _make_agilent_d(tmp_path / "first-source")
    second_raw = _make_agilent_d(tmp_path / "second-source")
    first_token = "deadbeef" + "0" * 24
    second_token = "deadbeef" + "1" * 24
    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(
        server,
        "_vendor_local_selections",
        {
            first_token: {
                "source_path": str(first_raw),
                "raw_name": first_raw.name,
                "source_kind": "agilent_d",
                "expected_polarities": ["NEG"],
                "selection_mode": "test-picker",
                "created_at": server.time.time(),
            },
            second_token: {
                "source_path": str(second_raw),
                "raw_name": second_raw.name,
                "source_kind": "agilent_d",
                "expected_polarities": ["NEG"],
                "selection_mode": "test-picker",
                "created_at": server.time.time(),
            },
        },
    )
    converter_status = _patch_local_conversion_runtime(
        monkeypatch, tmp_path, server
    )

    def fake_convert(raw_dir, output_path, **kwargs):
        output_path = Path(output_path)
        _write_mzml(output_path)
        output = validate_mzml(
            output_path,
            expected_ms1=True,
            expected_polarities=["NEG"],
            isolation_width_fallback_da=4.0,
        )
        _, converter_identity = server._vendor_converter_identity(
            converter_status,
            requested_backend=kwargs["backend"],
            accept_vendor_licenses=True,
        )
        return {
            "status": "completed",
            "source_kind": "agilent_d",
            "format": "Agilent MassHunter .d to mzML",
            "source": deepcopy(kwargs["raw_report"]),
            "output": output,
            "converter": converter_identity,
            "conversion_policy": {
                "peak_picking_strategy_requested": "auto",
                "peak_picking_strategy_resolved": "preserve",
                "preserves_ms1_and_ms2": True,
            },
        }

    monkeypatch.setattr(vendor_raw_conversion, "convert_agilent_d", fake_convert)
    responses = []
    with TestClient(server.app) as client:
        for token in (first_token, second_token):
            response = client.post(
                "/api/vendor/convert-local",
                json={
                    "selection_token": token,
                    "accept_vendor_licenses": True,
                    "allow_container_pull": False,
                    "backend": "auto",
                    "peak_picking_strategy": "auto",
                },
            )
            assert response.status_code == 200, response.text
            responses.append(response.json())

    task_ids = [response["uploaded"][0]["task_id"] for response in responses]
    task_paths = [Path(server._tasks[task_id]["file_path"]) for task_id in task_ids]
    assert task_paths[0] != task_paths[1]
    assert task_paths[0].name.startswith(f"{first_token}_")
    assert task_paths[1].name.startswith(f"{second_token}_")
    assert all(path.is_file() for path in task_paths)
    assert task_paths[0].read_bytes() == task_paths[1].read_bytes()


def test_legacy_eight_character_queued_destination_is_rejected_unchanged(
    monkeypatch,
    tmp_path,
):
    from web import server

    raw = _make_agilent_d(tmp_path)
    token = "cafebabe" + "2" * 24
    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    legacy_destination = upload_dir / f"{token[:8]}_{raw.stem}.mzML"
    original_bytes = b"legacy output must remain untouched"
    legacy_destination.write_bytes(original_bytes)
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_vendor_local_selections", {})
    _patch_local_conversion_runtime(monkeypatch, tmp_path, server)
    server._vendor_conversion_ledger().queue(
        token,
        metadata={
            "source_path": str(raw.resolve()),
            "raw_name": raw.name,
            "source_kind": "agilent_d",
            "expected_polarities": ["NEG"],
            "selection_mode": "test-picker",
            "selection_created_at_epoch": 1.0,
            "display_name": f"{raw.stem}.mzML",
            "destination": str(legacy_destination),
            "request_contract": {
                "backend": "auto",
                "intake_role": "sample",
                "blank_role": "background_blank",
                "peak_picking_strategy": "auto",
            },
        },
    )

    with TestClient(server.app) as client:
        response = client.post(
            "/api/vendor/convert-local",
            json={
                "selection_token": token,
                "accept_vendor_licenses": True,
                "allow_container_pull": False,
                "backend": "auto",
                "peak_picking_strategy": "auto",
            },
        )

    assert response.status_code == 409, response.text
    assert "legacy non-unique path" in response.json()["detail"]
    assert legacy_destination.read_bytes() == original_bytes
    assert server._vendor_conversion_ledger().get(token)["state"] == "queued"


def test_lost_vendor_claim_never_publishes_attempt_staging(
    monkeypatch,
    tmp_path,
):
    from utils import vendor_raw_conversion
    from utils.vendor_conversion_cache import ConversionLedger, LedgerStateError
    from web import server

    raw = _make_agilent_d(tmp_path)
    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(server, "TASKS_FILE", tmp_path / "tasks.json")
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_vendor_local_selections", {})
    converter_status = _patch_local_conversion_runtime(
        monkeypatch, tmp_path, server
    )
    monkeypatch.setattr(server, "_choose_macos_vendor_directory", lambda: raw)

    def fake_convert(raw_dir, output_path, **kwargs):
        output_path = Path(output_path)
        assert ".part.mzML" in output_path.name
        _write_mzml(output_path)
        output = validate_mzml(
            output_path,
            expected_ms1=True,
            expected_polarities=["NEG"],
            isolation_width_fallback_da=4.0,
        )
        _, converter_identity = server._vendor_converter_identity(
            converter_status,
            requested_backend=kwargs["backend"],
            accept_vendor_licenses=True,
        )
        return {
            "status": "completed",
            "source_kind": "agilent_d",
            "format": "Agilent MassHunter .d to mzML",
            "source": deepcopy(kwargs["raw_report"]),
            "output": output,
            "converter": converter_identity,
            "conversion_policy": {
                "peak_picking_strategy_requested": "auto",
                "peak_picking_strategy_resolved": "preserve",
                "preserves_ms1_and_ms2": True,
            },
        }

    monkeypatch.setattr(vendor_raw_conversion, "convert_agilent_d", fake_convert)
    monkeypatch.setattr(
        ConversionLedger,
        "renew_lease",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            LedgerStateError("synthetic stale claim")
        ),
    )

    with TestClient(server.app) as client:
        selection = client.post(
            "/api/vendor/select-local-folder", json={}
        ).json()
        converted = client.post(
            "/api/vendor/convert-local",
            json={
                "selection_token": selection["selection_token"],
                "accept_vendor_licenses": True,
                "allow_container_pull": False,
                "backend": "auto",
                "peak_picking_strategy": "auto",
            },
        )

    assert converted.status_code == 409, converted.text
    assert "stale claim" in converted.json()["detail"]
    assert list(upload_dir.iterdir()) == []
    ledger_item = server._vendor_conversion_ledger().get(
        selection["selection_token"]
    )
    assert ledger_item is not None
    assert ledger_item["state"] == "failed"


def test_cancelled_validated_registration_keeps_operation_fenced(
    monkeypatch,
    tmp_path,
):
    from fastapi import Request
    from web import server

    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    destination = upload_dir / "validated.mzML"
    destination.write_bytes(b"validated output")
    monkeypatch.setattr(server, "UPLOAD_DIR", upload_dir)
    monkeypatch.setattr(
        server,
        "VENDOR_CONVERSION_LEDGER_FILE",
        tmp_path / "vendor_conversion_ledger.json",
    )
    monkeypatch.setattr(server, "_vendor_conversion_active_claims", set())
    monkeypatch.setattr(server, "_vendor_conversion_accepting_claims", True)
    token = "b" * 32
    ledger = server._vendor_conversion_ledger()
    ledger.queue(token, metadata={"source_path": str(tmp_path / "source.d")})
    claimed = ledger.mark_converting(token)
    details = {
        "destination": str(destination),
        "display_name": "validated.mzML",
        "raw_name": "source.d",
        "source_kind": "agilent_d",
        "intake_role": "sample",
        "blank_role": "background_blank",
        "selection_mode": "test",
        "selection_audit": {},
        "conversion": {
            "source": {"acquisition": {}},
            "output": {
                "status": "valid",
                "size_bytes": destination.stat().st_size,
                "sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
            },
        },
    }
    validated = ledger.mark_validated(
        token,
        claimed["claim_token"],
        details=details,
    )
    registration_started = threading.Event()
    allow_registration_finish = threading.Event()

    def gated_registration(*_args, **_kwargs):
        registration_started.set()
        assert allow_registration_finish.wait(timeout=5)
        return {
            "uploaded": [],
            "blank_uploads": [],
            "count": 0,
            "failed": [],
            "failed_count": 0,
            "conversion": details["conversion"],
            "selection_mode": "test",
        }

    monkeypatch.setattr(
        server,
        "_register_validated_vendor_conversion",
        gated_registration,
    )
    request = Request({
        "type": "http",
        "http_version": "1.1",
        "method": "POST",
        "scheme": "http",
        "path": "/api/vendor/convert-local",
        "raw_path": b"/api/vendor/convert-local",
        "query_string": b"",
        "headers": [],
        "client": ("127.0.0.1", 12345),
        "server": ("127.0.0.1", 8516),
    })
    payload = server.VendorLocalConvertRequest(
        selection_token=token,
        accept_vendor_licenses=True,
    )

    async def scenario():
        operation = asyncio.create_task(
            server.convert_local_vendor_folder(payload, request)
        )
        assert await asyncio.to_thread(registration_started.wait, 5)
        operation.cancel()
        await asyncio.sleep(0.02)
        marker = (
            f"registration:{token}:{validated['registration_token']}"
        )
        assert marker in server._vendor_conversion_active_claims
        assert not operation.done()
        await server.shutdown()
        assert marker in server._vendor_conversion_active_claims
        allow_registration_finish.set()
        with pytest.raises(asyncio.CancelledError):
            await operation
        assert marker not in server._vendor_conversion_active_claims

    asyncio.run(scenario())
    registered = server._vendor_conversion_ledger().get(token)
    assert registered is not None
    assert registered["state"] == "registered"


def test_thread_worker_finishes_before_cancellation_is_propagated():
    from web import server

    started = threading.Event()
    finish = threading.Event()

    def gated_worker():
        started.set()
        assert finish.wait(timeout=5)
        return "finished"

    async def scenario():
        operation = asyncio.create_task(
            server._await_thread_completion(gated_worker)
        )
        assert await asyncio.to_thread(started.wait, 5)
        operation.cancel()
        await asyncio.sleep(0.02)
        assert not operation.done()
        finish.set()
        result, cancellation = await operation
        assert result == "finished"
        assert isinstance(cancellation, asyncio.CancelledError)

    asyncio.run(scenario())


def test_vendor_raw_ui_explains_whole_directory_selection():
    root = Path(__file__).resolve().parents[1]
    html = (root / "web" / "static" / "index.html").read_text(encoding="utf-8")
    upload_js = (
        root / "web" / "static" / "js" / "upload.js"
    ).read_text(encoding="utf-8")

    assert "选择 Agilent .d / 父目录" in html
    assert "选择 Thermo QE .raw / 批次目录" in html
    assert "macOS 本机选择" not in html
    assert "浏览器选择 Agilent .d（Chrome）" in html
    assert "内嵌窗口" in html
    assert "不要进入" in html
    assert "MSScan.bin" in html
    assert "MSPeak.bin" in html
    assert "AcqMethod.xml" in html
    assert 'id="vendorRawInput" webkitdirectory directory multiple' in html
    assert 'id="vendorConverterSetup"' in html
    assert 'id="vendorLicenseAccepted"' in html
    assert 'id="vendorAllowPull"' in html
    assert 'id="btnPrepareVendorConverter"' in html
    assert 'class="intake-step-rail"' in html
    assert "参考谱库" in html
    assert "选择样本" in html
    assert "可选质控" in html
    assert "参数与运行" in html
    assert "centroid 数据不二次质心化" in html
    assert "不做二次质心化" in html
    assert "镜像缓存后无需重复下载" in html
    assert "原始数据不会上传到云端" in html
    assert "本机已有可用的" in html
    assert "确认一次后仅记在本机浏览器" in html
    assert "分析和厂商数据转换过程绝不联网拉取镜像" in html
    assert "开始分析只使用本机已有转换器，绝不联网" in html
    assert 'id="advancedAnalysisSettings"' in html
    assert "selectLocalVendorFolder" in upload_js
    assert "selectLocalThermoFolder" in upload_js
    assert "/vendor/select-local-thermo-folder" in upload_js
    assert "updateVendorSourceRole" in upload_js
    assert "作为样本分析" in upload_js
    assert "作为空白背景" in upload_js
    assert "intake_role" in upload_js
    assert "prepareVendorConverter" in upload_js
    assert "recommended_peak_picking_strategy" in upload_js
    assert "spectrum_representation" in upload_js
    assert "VENDOR_LICENSE_STORAGE_KEY" in upload_js
    assert "window.localStorage.setItem" in upload_js
    assert "_vendorLocalSelections" in upload_js
    assert "/vendor/convert-local" in upload_js
