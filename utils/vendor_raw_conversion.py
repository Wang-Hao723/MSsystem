"""Validated conversion of vendor raw mass-spectrometry data to mzML.

Agilent MassHunter ``.d`` datasets and Thermo Fisher ``.raw`` files are
proprietary inputs.  They cannot be decoded safely by renaming a file or by a
pure-Python parser.  This module therefore orchestrates the vendor readers
distributed with ProteoWizard and validates the resulting mzML before
MSsystem can analyse it.

Two execution backends are supported:

* a local ``msconvert`` executable (normally Windows);
* the ProteoWizard Docker/Wine image for macOS/Linux.

The conversion deliberately applies no ``msLevel`` filter.  In particular,
``msLevel 2-`` is forbidden because it silently removes the survey MS1 scans
needed for XIC, ion-identity networking and auditable precursor isolation.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Iterable, Mapping


DEFAULT_PWIZ_DOCKER_IMAGE = (
    "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:latest"
)
# The July 2025 Wine 10.6 base used by current tags invokes an i386 Unix
# helper even for the 64-bit msconvert executable.  Docker Desktop's Rosetta
# path translates x86_64 but not i386, while QEMU-i386 aborts in Wine's fixed
# mmap setup on Apple Silicon.  This last official pre-Wine-10 build exposes a
# true ``wine64`` entry point and is pinned for reproducibility on ARM hosts.
APPLE_SILICON_PWIZ_DOCKER_IMAGE = (
    "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:"
    "3.0.25142-eb0d7d9"
)
AGILENT_REQUIRED_FILES = (
    "AcqData/MSScan.bin",
    "AcqData/MSPeak.bin",
    "AcqData/AcqMethod.xml",
)
THERMO_RAW_SUFFIX = ".raw"


class VendorRawConversionError(RuntimeError):
    """Base class for conversion and validation failures."""


class VendorRawValidationError(VendorRawConversionError):
    """The supplied raw directory or produced mzML is invalid."""


class ConverterUnavailableError(VendorRawConversionError):
    """No usable ProteoWizard execution backend is available."""


class VendorLicenseNotAcceptedError(VendorRawConversionError):
    """The Docker backend was requested without explicit licence acceptance."""


def _arm64_host() -> bool:
    return platform.machine().lower() in {"arm64", "aarch64"}


def _default_docker_image() -> str:
    return (
        APPLE_SILICON_PWIZ_DOCKER_IMAGE
        if _arm64_host()
        else DEFAULT_PWIZ_DOCKER_IMAGE
    )


def _docker_wine_runner() -> str:
    # A 64-bit Windows msconvert must never be routed through Wine's i386
    # helper on ARM.  Historical official images expose this wrapper.
    return "wine64_anyuser" if _arm64_host() else "wine"


def _verify_docker_wine_runner(docker: str, image: str) -> None:
    """Fail before a batch when an ARM image lacks the true 64-bit runner."""

    if not _arm64_host():
        return
    probe_command = [
        docker,
        "run",
        "--rm",
        "--platform",
        "linux/amd64",
        image,
        "/bin/sh",
        "-lc",
        "command -v wine64 >/dev/null && command -v wine64_anyuser >/dev/null",
    ]
    runner_probe = subprocess.run(
        probe_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=60,
        check=False,
    )
    if runner_probe.returncode != 0:
        raise ConverterUnavailableError(
            "The selected ProteoWizard image has no true 64-bit Wine "
            "runner and cannot be used safely on ARM. Use the pinned "
            f"official image {APPLE_SILICON_PWIZ_DOCKER_IMAGE}."
        )


def _local_name(tag: str) -> str:
    return str(tag).rsplit("}", 1)[-1]


def _sha256_file(path: Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _sha1_file(path: Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    """Return SHA-1 only for matching ProteoWizard's source-file provenance."""
    digest = hashlib.sha1(usedforsecurity=False)
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _fingerprint_file(
    path: Path,
    chunk_size: int = 4 * 1024 * 1024,
) -> dict[str, Any]:
    """Hash a large source once and reject files changing during inspection."""

    before = path.stat()
    sha256 = hashlib.sha256()
    sha1 = hashlib.sha1(usedforsecurity=False)
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            sha256.update(chunk)
            sha1.update(chunk)
    after = path.stat()
    if (
        before.st_dev != after.st_dev
        or before.st_ino != after.st_ino
        or before.st_size != after.st_size
        or before.st_mtime_ns != after.st_mtime_ns
        or before.st_ctime_ns != after.st_ctime_ns
    ):
        raise VendorRawValidationError(
            f"Vendor raw source changed while it was being fingerprinted: {path}"
        )
    return {
        "size_bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "sha256": sha256.hexdigest(),
        "sha1_for_mzml_provenance": sha1.hexdigest(),
        # Excluded from the cache key: this binds the just-computed hashes to
        # the current filesystem object so key construction can avoid a second
        # full read of a multi-GB vendor payload.
        "inspection_stat": {
            "device": after.st_dev,
            "inode": after.st_ino,
            "ctime_ns": after.st_ctime_ns,
        },
    }


class _DigestingReader:
    """Update a digest while a streaming parser consumes a binary file."""

    def __init__(self, handle: Any, digest: Any):
        self._handle = handle
        self._digest = digest

    def read(self, size: int = -1) -> bytes:
        data = self._handle.read(size)
        if data:
            self._digest.update(data)
        return data


def _safe_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else number


def _unique_numbers(values: Iterable[Any]) -> list[float]:
    numbers = {
        round(number, 9)
        for value in values
        if (number := _safe_float(value)) is not None
    }
    return sorted(numbers)


def _find_first_local(root: ET.Element, name: str) -> ET.Element | None:
    for element in root.iter():
        if _local_name(element.tag) == name:
            return element
    return None


def _child_texts(element: ET.Element) -> dict[str, str]:
    out: dict[str, str] = {}
    for child in list(element):
        text = str(child.text or "").strip()
        if text:
            out[_local_name(child.tag)] = text
    return out


def _parse_isolation_width(value: str) -> float | None:
    text = str(value or "")
    match = re.search(
        r"(?:~|approximately\s*)?(\d+(?:\.\d+)?)\s*(?:amu|da)\b",
        text,
        re.IGNORECASE,
    )
    if match:
        width = _safe_float(match.group(1))
        return width if width and 0 < width <= 50 else None
    number = _safe_float(text)
    return number if number and 0 < number <= 50 else None


def parse_agilent_acquisition_method(raw_dir: str | Path) -> dict[str, Any]:
    """Extract auditable acquisition settings from ``AcqMethod.xml``.

    Agilent stores the detailed report as escaped XML inside a second XML
    document.  ElementTree decodes that text for us; no vendor binary parsing
    is attempted here.
    """
    raw_path = Path(raw_dir).expanduser().resolve()
    method_path = raw_path / "AcqData" / "AcqMethod.xml"
    if not method_path.is_file():
        return {
            "status": "missing",
            "method_file": "AcqData/AcqMethod.xml",
            "warnings": ["AcqMethod.xml is missing; acquisition settings unavailable."],
        }

    try:
        outer = ET.parse(method_path).getroot()
    except (OSError, ET.ParseError) as exc:
        return {
            "status": "unreadable",
            "method_file": "AcqData/AcqMethod.xml",
            "warnings": [f"AcqMethod.xml could not be parsed: {exc}"],
        }

    method_name_node = _find_first_local(outer, "MethodName")
    method_name = str(method_name_node.text or "").strip() if method_name_node is not None else ""
    embedded_reports = [
        element
        for element in outer.iter()
        if _local_name(element.tag) in {"RCDevicesXml", "SCICDevicesXml"}
        and str(element.text or "").strip()
    ]
    if not embedded_reports:
        return {
            "status": "partial",
            "method_file": "AcqData/AcqMethod.xml",
            "method_name": method_name,
            "warnings": ["Detailed Q-TOF acquisition report is absent."],
        }

    reports: list[ET.Element] = []
    report_errors: list[str] = []
    for embedded in embedded_reports:
        try:
            reports.append(ET.fromstring(str(embedded.text)))
        except ET.ParseError as exc:
            report_errors.append(f"{_local_name(embedded.tag)}: {exc}")
    if not reports:
        return {
            "status": "partial",
            "method_file": "AcqData/AcqMethod.xml",
            "method_name": method_name,
            "warnings": [
                "Embedded acquisition reports could not be parsed: "
                + "; ".join(report_errors)
            ],
        }

    sections: list[dict[str, str]] = []
    devices: list[str] = []
    for report in reports:
        for element in report.iter():
            local = _local_name(element.tag)
            if local == "SectionInfo":
                sections.append(_child_texts(element))
            elif local == "Device":
                values = _child_texts(element)
                if values.get("Name"):
                    devices.append(values["Name"])

    def values_for_id(identifier: str) -> list[str]:
        return [
            section["Value"]
            for section in sections
            if section.get("ID") == identifier and section.get("Value")
        ]

    names = [section.get("Name", "") for section in sections]
    polarities = []
    for value in values_for_id("ionPolarity"):
        normalized = str(value).strip().upper()
        if normalized.startswith("NEG"):
            polarities.append("NEG")
        elif normalized.startswith("POS"):
            polarities.append("POS")
    polarities = sorted(set(polarities))

    isolation_raw_values = list(dict.fromkeys(values_for_id("IsolationWidthMS2")))
    isolation_widths = [
        width
        for value in isolation_raw_values
        if (width := _parse_isolation_width(value)) is not None
    ]
    isolation_widths = sorted(set(isolation_widths))

    collision_energies: list[float] = []
    for section in sections:
        if "CollisionEnergyFixed" not in section.get("ID", ""):
            continue
        for token in re.split(r"[,;\s]+", section.get("Value", "")):
            number = _safe_float(token)
            if number is not None:
                collision_energies.append(number)

    storage_modes = list(dict.fromkeys(values_for_id("StorageMode")))
    normalized_storage_modes = {
        str(value).strip().casefold() for value in storage_modes
    }
    ms1_rates = _unique_numbers(values_for_id("MS1ScanRate"))
    ms2_rates = _unique_numbers(values_for_id("MS2ScanRate"))
    acquisition_mode = "AutoMS2" if any("Acquisition Mode AutoMS2" in name for name in names) else ""
    instrument = "Q-TOF" if any("Q-TOF" in device.upper() for device in devices) else ""

    warnings: list[str] = [
        f"One embedded acquisition report was unreadable: {error}"
        for error in report_errors
    ]
    if len(isolation_widths) > 1:
        warnings.append(
            "Multiple MS/MS isolation widths occur in the acquisition method; "
            "a single automatic fallback will not be applied."
        )
    if not polarities:
        warnings.append("Ion polarity was not found in the acquisition method.")

    return {
        "status": "ok",
        "method_file": "AcqData/AcqMethod.xml",
        "method_name": method_name,
        "instrument": instrument,
        "devices": list(dict.fromkeys(devices)),
        "acquisition_mode": acquisition_mode,
        "polarities": polarities,
        "ms1_expected": bool(ms1_rates and max(ms1_rates) > 0),
        "ms2_expected": bool(ms2_rates and max(ms2_rates) > 0),
        "ms1_range_mz": {
            "min": _unique_numbers(values_for_id("MS1MinRange")),
            "max": _unique_numbers(values_for_id("MS1MaxRange")),
        },
        "ms2_range_mz": {
            "min": _unique_numbers(values_for_id("MS2MinRange")),
            "max": _unique_numbers(values_for_id("MS2MaxRange")),
        },
        "ms1_scan_rates_hz": ms1_rates,
        "ms2_scan_rates_hz": ms2_rates,
        "storage_modes": storage_modes,
        # Agilent reports ``Both`` when the acquisition stores profile and
        # centroid representations together.  Treating it as an unknown mode
        # made automatic conversion skip vendor peak picking and emitted very
        # large profile mzML files.
        "centroid_acquisition_present": bool(
            normalized_storage_modes & {"centroid", "both"}
        ),
        "profile_acquisition_present": bool(
            normalized_storage_modes & {"profile", "both"}
        ),
        "isolation_width_raw": isolation_raw_values,
        "isolation_width_da": (
            isolation_widths[0] if len(isolation_widths) == 1 else None
        ),
        "isolation_width_source": (
            "Agilent AcqMethod.xml / IsolationWidthMS2"
            if len(isolation_widths) == 1 else None
        ),
        "collision_energies": sorted(set(collision_energies)),
        "warnings": warnings,
    }


def inspect_agilent_d(raw_dir: str | Path) -> dict[str, Any]:
    """Validate an Agilent MassHunter ``.d`` directory and fingerprint it."""
    raw_path = Path(raw_dir).expanduser().resolve()
    errors: list[str] = []
    if not raw_path.exists():
        raise VendorRawValidationError(f"Raw directory does not exist: {raw_path}")
    if not raw_path.is_dir():
        raise VendorRawValidationError(f"Agilent input must be a .d directory: {raw_path}")
    if raw_path.suffix.lower() != ".d":
        errors.append("Agilent MassHunter raw input must have a .d directory suffix.")

    signatures: dict[str, dict[str, Any]] = {}
    for relative in AGILENT_REQUIRED_FILES:
        path = raw_path / relative
        if not path.is_file():
            errors.append(f"Missing required Agilent file: {relative}")
            continue
        signatures[relative] = _fingerprint_file(path)
    if errors:
        raise VendorRawValidationError("; ".join(errors))

    file_count = 0
    total_bytes = 0
    for path in raw_path.rglob("*"):
        if path.is_file():
            file_count += 1
            total_bytes += path.stat().st_size

    acquisition = parse_agilent_acquisition_method(raw_path)
    return {
        "status": "ok",
        "source_kind": "agilent_d",
        "format": "Agilent MassHunter .d",
        "raw_name": raw_path.name,
        "source_path": str(raw_path),
        "file_count": file_count,
        "total_bytes": total_bytes,
        "required_file_signatures": signatures,
        "acquisition": acquisition,
    }


def inspect_thermo_raw(
    raw_file: str | Path,
    *,
    expected_polarities: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Validate and fingerprint one Thermo Fisher ``.raw`` file.

    Instrument model, stored spectrum representation and acquisition settings
    are intentionally not guessed from the binary container or filename.  They
    are read through the vendor reader and audited from the converted mzML.
    ``expected_polarities`` is optional contextual evidence (for example, a
    user-selected POS/NEG batch) and is checked against every converted scan.
    """

    raw_path = Path(raw_file).expanduser().resolve()
    if not raw_path.exists():
        raise VendorRawValidationError(
            f"Thermo RAW file does not exist: {raw_path}"
        )
    if not raw_path.is_file():
        raise VendorRawValidationError(
            f"Thermo input must be one .raw file: {raw_path}"
        )
    if raw_path.suffix.lower() != THERMO_RAW_SUFFIX:
        raise VendorRawValidationError(
            "Thermo Fisher raw input must have a .raw file suffix."
        )
    if raw_path.stat().st_size <= 0:
        raise VendorRawValidationError(f"Thermo RAW file is empty: {raw_path}")

    polarities = sorted({
        str(mode).strip().upper()
        for mode in (expected_polarities or [])
        if str(mode).strip().upper() in {"POS", "NEG"}
    })
    signature = _fingerprint_file(raw_path)
    return {
        "status": "ok",
        "source_kind": "thermo_raw",
        "format": "Thermo Fisher RAW",
        "raw_name": raw_path.name,
        "source_path": str(raw_path),
        "file_count": 1,
        "total_bytes": signature["size_bytes"],
        "source_file_signature": signature,
        "acquisition": {
            "status": "verified_after_vendor_conversion",
            "instrument": "",
            "instrument_verification": "converted_mzML_cv_terms",
            "acquisition_mode": "",
            "polarities": polarities,
            "polarity_expectation_source": (
                "user_selected_batch_context" if polarities else None
            ),
            "ms1_expected": None,
            "ms2_expected": True,
            "centroid_acquisition_present": None,
            "profile_acquisition_present": None,
            "isolation_width_da": None,
            "isolation_width_source": None,
            "warnings": ([] if polarities else [
                "No polarity expectation was supplied; converted per-spectrum "
                "polarity remains authoritative."
            ]),
        },
    }


def inspect_vendor_raw(
    source: str | Path,
    *,
    expected_polarities: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Dispatch validation without interpreting proprietary binary payloads."""

    path = Path(source).expanduser().resolve()
    if path.suffix.lower() == ".d":
        return inspect_agilent_d(path)
    if path.suffix.lower() == THERMO_RAW_SUFFIX:
        return inspect_thermo_raw(
            path,
            expected_polarities=expected_polarities,
        )
    raise VendorRawValidationError(
        "Unsupported vendor raw source; expected Agilent .d or Thermo .raw."
    )


def resolve_peak_picking_strategy(
    raw_report: Mapping[str, Any],
    peak_picking_strategy: str = "auto",
) -> str:
    """Resolve the converter policy before a cache lookup or conversion.

    This function is shared by the content-addressed cache and the converter so
    ``auto`` cannot resolve differently on the hit and miss paths.
    """

    requested = str(peak_picking_strategy or "auto").strip().lower()
    if requested not in {"auto", "preserve", "vendor", "openms"}:
        raise ValueError(
            "peak_picking_strategy must be auto, preserve, vendor, or openms"
        )
    if requested != "auto":
        return requested
    source_kind = str(raw_report.get("source_kind") or "agilent_d")
    if source_kind == "thermo_raw":
        return "vendor"
    acquisition = raw_report.get("acquisition") or {}
    centroid_present = acquisition.get("centroid_acquisition_present") is True
    profile_present = acquisition.get("profile_acquisition_present") is True
    return "preserve" if centroid_present and not profile_present else "vendor"


def _reuse_inspected_raw_report(
    raw_source: str | Path,
    raw_report: Mapping[str, Any],
    *,
    expected_source_kind: str,
) -> dict[str, Any]:
    """Validate a caller-supplied inspection before avoiding a second read."""

    if not isinstance(raw_report, Mapping) or raw_report.get("status") != "ok":
        raise VendorRawValidationError(
            "Reused vendor raw inspection must have status=ok."
        )
    report_kind = str(raw_report.get("source_kind") or "")
    if report_kind != expected_source_kind:
        raise VendorRawValidationError(
            "Reused vendor raw inspection has the wrong source kind."
        )
    source = Path(raw_source).expanduser().resolve()
    report_source = Path(
        str(raw_report.get("source_path") or "")
    ).expanduser().resolve()
    if report_source != source:
        raise VendorRawValidationError(
            "Reused vendor raw inspection belongs to a different source path."
        )
    try:
        return json.loads(json.dumps(raw_report, allow_nan=False))
    except (TypeError, ValueError) as exc:
        raise VendorRawValidationError(
            "Reused vendor raw inspection is not valid JSON provenance."
        ) from exc


def _find_local_msconvert(explicit: str | Path | None = None) -> str | None:
    candidates: list[str] = []
    if explicit:
        candidates.append(str(explicit))
    configured = os.environ.get("MSSYSTEM_MSCONVERT_BIN", "").strip()
    if configured:
        candidates.append(configured)
    discovered = shutil.which("msconvert") or shutil.which("msconvert.exe")
    if discovered:
        candidates.append(discovered)
    if os.name == "nt":
        for root_name in ("ProgramFiles", "ProgramFiles(x86)"):
            root = os.environ.get(root_name)
            if not root:
                continue
            for path in sorted(Path(root).glob("ProteoWizard*/msconvert.exe"), reverse=True):
                candidates.append(str(path))
    for candidate in candidates:
        path = Path(candidate).expanduser()
        if path.is_file():
            return str(path.resolve())
    return None


def _run_probe(command: list[str], timeout: float = 8.0) -> subprocess.CompletedProcess:
    return subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=timeout,
        check=False,
    )


def converter_status(
    *,
    msconvert_bin: str | Path | None = None,
    docker_image: str | None = None,
) -> dict[str, Any]:
    """Return non-mutating backend diagnostics for the UI and CLI."""
    image = (
        docker_image
        or os.environ.get("MSSYSTEM_PWIZ_DOCKER_IMAGE", "").strip()
        or _default_docker_image()
    )
    local = _find_local_msconvert(msconvert_bin)
    docker = shutil.which("docker")
    daemon_available = False
    daemon_error = ""
    image_available = False
    image_identity = ""
    if docker:
        try:
            probe = _run_probe([docker, "info", "--format", "{{.ServerVersion}}"])
            daemon_available = probe.returncode == 0
            if not daemon_available:
                daemon_error = (probe.stderr or probe.stdout).strip()
        except (OSError, subprocess.TimeoutExpired) as exc:
            daemon_error = str(exc)
        if daemon_available:
            try:
                inspected = _run_probe(
                    [
                        docker,
                        "image",
                        "inspect",
                        image,
                        "--format",
                        "{{.Id}}",
                    ]
                )
                image_available = inspected.returncode == 0
                if image_available:
                    image_identity = inspected.stdout.strip()
            except (OSError, subprocess.TimeoutExpired):
                pass

    recommended = (
        "local" if local else "docker"
        if docker and daemon_available and image_available else None
    )
    return {
        "platform": {
            "system": platform.system(),
            "machine": platform.machine(),
        },
        "local_msconvert": {
            "available": bool(local),
            "path": local,
        },
        "docker": {
            "cli_available": bool(docker),
            "path": docker,
            "daemon_available": daemon_available,
            "daemon_error": daemon_error,
            "image": image,
            "image_available": image_available,
            "image_identity": image_identity,
            "requires_vendor_license_acceptance": True,
            "requires_amd64_emulation": _arm64_host(),
            "wine_runner": _docker_wine_runner(),
            "apple_silicon_pinned_image": (
                image == APPLE_SILICON_PWIZ_DOCKER_IMAGE
            ),
        },
        "recommended_backend": recommended,
    }


def build_local_msconvert_command(
    msconvert_bin: str,
    raw_dir: str | Path,
    output_dir: str | Path,
    output_name: str,
    *,
    vendor_peak_picking: bool = False,
) -> list[str]:
    """Build a loss-avoiding ProteoWizard command.

    No ``--filter msLevel`` or generic peak-picking filter is added.  Agilent
    centroid storage is read through the native vendor reader and both MS1 and
    MS2 are retained.
    """
    command = [
        str(msconvert_bin),
        str(Path(raw_dir).resolve()),
        "--mzML",
        "--64",
        "--zlib",
    ]
    # ProteoWizard requires vendor peak picking to be the first filter.  No
    # msLevel filter is ever added, so both survey MS1 and product-ion MS2 are
    # retained.
    if vendor_peak_picking:
        command.extend(["--filter", "peakPicking true 1-"])
    command.extend([
        "--outdir",
        str(Path(output_dir).resolve()),
        "--outfile",
        output_name,
    ])
    return command


def build_docker_msconvert_command(
    docker_bin: str,
    image: str,
    raw_dir: str | Path,
    output_dir: str | Path,
    output_name: str,
    *,
    wine_runner: str | None = None,
    vendor_peak_picking: bool = False,
) -> list[str]:
    raw_path = Path(raw_dir).resolve()
    output_path = Path(output_dir).resolve()
    command = [
        docker_bin,
        "run",
        "--rm",
    ]
    if _arm64_host():
        command.extend(["--platform", "linux/amd64"])
    command.extend(
        [
            "-e",
            "WINEDEBUG=-all",
            "-v",
            f"{raw_path.parent}:/input:ro",
            "-v",
            f"{output_path}:/output",
            image,
            wine_runner or _docker_wine_runner(),
            "msconvert",
            f"/input/{raw_path.name}",
            "--mzML",
            "--64",
            "--zlib",
        ]
    )
    if vendor_peak_picking:
        command.extend(["--filter", "peakPicking true 1-"])
    command.extend(
        [
            "--outdir",
            "/output",
            "--outfile",
            output_name,
        ]
    )
    return command


def _docker_execution_image(
    image_tag: str,
    docker_status: Mapping[str, Any],
    *,
    verified_status_reused: bool,
) -> str:
    """Bind a verified preflight to the immutable image it inspected.

    Direct library callers that do not supply a preflight retain the historical
    tag-based command for compatibility.  The server path supplies one verified
    status report, so execution must use its Docker image ID and cannot drift if
    the mutable tag is moved between cache-key construction and ``docker run``.
    """

    if not verified_status_reused:
        return image_tag
    image_identity = str(docker_status.get("image_identity") or "").strip()
    if not image_identity.startswith("sha256:") or not image_identity[7:]:
        raise ConverterUnavailableError(
            "The verified ProteoWizard preflight has no immutable sha256 image "
            "identity; refusing a tag-based Docker execution."
        )
    return image_identity


def _cv_value(element: ET.Element, accession: str) -> str | None:
    for child in element.iter():
        if _local_name(child.tag) != "cvParam":
            continue
        if child.attrib.get("accession") == accession:
            return child.attrib.get("value")
    return None


def _cv_terms(element: ET.Element) -> list[dict[str, str]]:
    return [
        {
            "accession": str(child.attrib.get("accession") or ""),
            "name": str(child.attrib.get("name") or ""),
            "value": str(child.attrib.get("value") or ""),
        }
        for child in element.iter()
        if _local_name(child.tag) == "cvParam"
    ]


def validate_mzml(
    mzml_path: str | Path,
    *,
    expected_ms1: bool | None = None,
    expected_ms2: bool = True,
    expected_polarities: Iterable[str] | None = None,
    isolation_width_fallback_da: float | None = None,
    expected_source_sha1: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Stream-validate an mzML conversion without decoding every peak array."""
    path = Path(mzml_path).expanduser().resolve()
    if not path.is_file():
        raise VendorRawValidationError(f"Converted mzML does not exist: {path}")
    if path.stat().st_size <= 0:
        raise VendorRawValidationError(f"Converted mzML is empty: {path}")

    declared_count: int | None = None
    actual_count = 0
    ms_levels: dict[int, int] = {}
    nonempty_levels: dict[int, int] = {}
    centroid_count = 0
    profile_count = 0
    representation_by_ms_level: dict[int, dict[str, int]] = {}
    polarities = {"POS": 0, "NEG": 0}
    precursor_count = 0
    precursor_ref_count = 0
    isolation_target_count = 0
    isolation_lower_count = 0
    isolation_upper_count = 0
    software: list[dict[str, str]] = []
    source_file_sha1: dict[str, str] = {}
    source_file_formats: list[dict[str, str]] = []
    referenceable_groups: dict[str, list[dict[str, str]]] = {}
    instrument_models: list[dict[str, str]] = []
    instrument_serial_numbers: list[str] = []
    ion_sources: list[dict[str, str]] = []
    mass_analyzers: list[dict[str, str]] = []
    detectors: list[dict[str, str]] = []
    run_metadata: dict[str, str] = {}
    collision_energies: dict[int, set[float]] = {}
    resolving_powers: dict[int, set[float]] = {}
    isolation_widths: set[float] = set()
    filter_string_examples: dict[int, list[str]] = {}
    root_names: set[str] = set()

    mzml_sha256 = hashlib.sha256()
    before_parse = path.stat()
    mzml_handle = None
    try:
        mzml_handle = path.open("rb")
        parser_input = _DigestingReader(mzml_handle, mzml_sha256)
        events = ET.iterparse(parser_input, events=("start", "end"))
        for event, element in events:
            local = _local_name(element.tag)
            if event == "start":
                if not root_names:
                    root_names.add(local)
                if local == "spectrumList" and declared_count is None:
                    try:
                        declared_count = int(element.attrib.get("count", ""))
                    except ValueError:
                        declared_count = None
                elif local == "run" and not run_metadata:
                    run_metadata = {
                        key: str(element.attrib.get(key) or "")
                        for key in (
                            "id",
                            "startTimeStamp",
                            "defaultInstrumentConfigurationRef",
                            "defaultSourceFileRef",
                        )
                        if str(element.attrib.get(key) or "")
                    }
                continue

            if local == "software":
                software.append(
                    {
                        "id": str(element.attrib.get("id") or ""),
                        "version": str(element.attrib.get("version") or ""),
                    }
                )
                continue
            if local == "sourceFile":
                source_name = str(
                    element.attrib.get("name")
                    or element.attrib.get("id")
                    or ""
                ).strip()
                checksum = _cv_value(element, "MS:1000569")
                if source_name and checksum:
                    source_file_sha1[source_name] = str(checksum).strip().lower()
                for term in _cv_terms(element):
                    name = term["name"].lower()
                    if "format" in name or "nativeid" in name:
                        source_file_formats.append({
                            "source_name": source_name,
                            **term,
                        })
                continue
            if local == "referenceableParamGroup":
                group_id = str(element.attrib.get("id") or "")
                if group_id:
                    referenceable_groups[group_id] = _cv_terms(element)
                continue
            if local == "instrumentConfiguration":
                configuration_terms = [
                    {
                        "accession": str(child.attrib.get("accession") or ""),
                        "name": str(child.attrib.get("name") or ""),
                        "value": str(child.attrib.get("value") or ""),
                    }
                    for child in list(element)
                    if _local_name(child.tag) == "cvParam"
                ]
                for child in element.iter():
                    if _local_name(child.tag) != "referenceableParamGroupRef":
                        continue
                    configuration_terms.extend(
                        referenceable_groups.get(
                            str(child.attrib.get("ref") or ""), []
                        )
                    )
                for term in configuration_terms:
                    if term["accession"] == "MS:1000529":
                        if term["value"]:
                            instrument_serial_numbers.append(term["value"])
                    elif term["name"] and term not in instrument_models:
                        # Component terms are reported separately below.  A
                        # referenced common-instrument group normally contains
                        # the exact vendor model CV term (for example Q Exactive).
                        instrument_models.append(term)
                for component_name, destination in (
                    ("source", ion_sources),
                    ("analyzer", mass_analyzers),
                    ("detector", detectors),
                ):
                    for component in element.iter():
                        if _local_name(component.tag) != component_name:
                            continue
                        for term in _cv_terms(component):
                            if term not in destination:
                                destination.append(term)
                continue
            if local != "spectrum":
                continue

            actual_count += 1
            level = _safe_float(_cv_value(element, "MS:1000511"))
            level_int = int(level) if level is not None and level.is_integer() else 0
            ms_levels[level_int] = ms_levels.get(level_int, 0) + 1
            try:
                array_length = int(element.attrib.get("defaultArrayLength", "0"))
            except ValueError:
                array_length = 0
            if array_length > 0:
                nonempty_levels[level_int] = nonempty_levels.get(level_int, 0) + 1

            accessions = {
                child.attrib.get("accession")
                for child in element.iter()
                if _local_name(child.tag) == "cvParam"
            }
            if "MS:1000127" in accessions:
                centroid_count += 1
                representation_by_ms_level.setdefault(
                    level_int, {"centroid": 0, "profile": 0, "unknown": 0}
                )["centroid"] += 1
            if "MS:1000128" in accessions or "MS:1000126" in accessions:
                profile_count += 1
                representation_by_ms_level.setdefault(
                    level_int, {"centroid": 0, "profile": 0, "unknown": 0}
                )["profile"] += 1
            if not (
                "MS:1000127" in accessions
                or "MS:1000128" in accessions
                or "MS:1000126" in accessions
            ):
                representation_by_ms_level.setdefault(
                    level_int, {"centroid": 0, "profile": 0, "unknown": 0}
                )["unknown"] += 1
            if "MS:1000130" in accessions:
                polarities["POS"] += 1
            if "MS:1000129" in accessions:
                polarities["NEG"] += 1

            collision_energy = _safe_float(_cv_value(element, "MS:1000045"))
            if collision_energy is not None:
                collision_energies.setdefault(level_int, set()).add(
                    round(collision_energy, 9)
                )
            resolving_power = _safe_float(_cv_value(element, "MS:1000800"))
            if resolving_power is not None:
                resolving_powers.setdefault(level_int, set()).add(
                    round(resolving_power, 9)
                )
            filter_string = str(
                _cv_value(element, "MS:1000512") or ""
            ).strip()
            examples = filter_string_examples.setdefault(level_int, [])
            if filter_string and filter_string not in examples and len(examples) < 5:
                examples.append(filter_string)

            for child in element.iter():
                child_local = _local_name(child.tag)
                if child_local == "precursor":
                    precursor_count += 1
                    if str(child.attrib.get("spectrumRef") or "").strip():
                        precursor_ref_count += 1
                elif child_local == "cvParam":
                    accession = child.attrib.get("accession")
                    if accession == "MS:1000827":
                        isolation_target_count += 1
                    elif accession == "MS:1000828":
                        isolation_lower_count += 1
                    elif accession == "MS:1000829":
                        isolation_upper_count += 1
            lower = _safe_float(_cv_value(element, "MS:1000828"))
            upper = _safe_float(_cv_value(element, "MS:1000829"))
            if lower is not None and upper is not None:
                isolation_widths.add(round(lower + upper, 9))
            element.clear()
    except (OSError, ET.ParseError) as exc:
        raise VendorRawValidationError(f"Converted mzML is not valid XML: {exc}") from exc
    finally:
        if mzml_handle is not None:
            mzml_handle.close()
    after_parse = path.stat()
    if (
        before_parse.st_size != after_parse.st_size
        or before_parse.st_mtime_ns != after_parse.st_mtime_ns
    ):
        raise VendorRawValidationError(
            f"Converted mzML changed while it was being validated: {path}"
        )

    errors: list[str] = []
    warnings: list[str] = []
    if actual_count == 0:
        errors.append("mzML contains no spectra.")
    if declared_count is not None and declared_count != actual_count:
        errors.append(
            f"spectrumList declares {declared_count} spectra but {actual_count} were read."
        )

    ms1_count = ms_levels.get(1, 0)
    ms2_count = ms_levels.get(2, 0)
    mass_spectrum_count = sum(
        count for level, count in ms_levels.items() if level >= 1
    )
    auxiliary_spectrum_count = ms_levels.get(0, 0)
    nonempty_ms1 = nonempty_levels.get(1, 0)
    nonempty_ms2 = nonempty_levels.get(2, 0)
    nonempty_auxiliary = nonempty_levels.get(0, 0)
    if expected_ms1 is True and ms1_count == 0:
        errors.append(
            "Acquisition method expects survey MS1, but conversion produced 0 MS1 spectra. "
            "Do not use an msLevel 2 filter."
        )
    if expected_ms1 is True and ms1_count > 0 and nonempty_ms1 == 0:
        errors.append("All converted MS1 spectra are empty.")
    if expected_ms2 and ms2_count == 0:
        errors.append("Conversion produced 0 MS2 spectra.")
    if expected_ms2 and ms2_count > 0 and nonempty_ms2 == 0:
        errors.append("All converted MS2 spectra are empty.")

    expected_modes = {
        str(mode).strip().upper()
        for mode in (expected_polarities or [])
        if str(mode).strip().upper() in {"POS", "NEG"}
    }
    observed_modes = {mode for mode, count in polarities.items() if count}
    if expected_modes and observed_modes and not observed_modes.issubset(expected_modes):
        errors.append(
            "Converted polarity does not agree with AcqMethod.xml: "
            f"expected {sorted(expected_modes)}, observed {sorted(observed_modes)}."
        )
    if not observed_modes:
        warnings.append("No per-spectrum polarity CV terms were found.")
    if auxiliary_spectrum_count:
        warnings.append(
            f"mzML retains {auxiliary_spectrum_count} auxiliary detector spectra "
            "without an MS level (for example DAD); they are not counted as "
            "MS1/MS2 and are excluded from compound matching."
        )

    if ms2_count and precursor_ref_count < ms2_count:
        warnings.append(
            f"{ms2_count - precursor_ref_count} MS2 spectra lack an MS1 spectrumRef."
        )
    if ms2_count and isolation_target_count < ms2_count:
        warnings.append(
            f"{ms2_count - isolation_target_count} MS2 spectra lack an isolation target m/z."
        )
    offsets_complete = (
        ms2_count > 0
        and isolation_lower_count >= ms2_count
        and isolation_upper_count >= ms2_count
    )
    if ms2_count and not offsets_complete:
        if isolation_width_fallback_da is not None:
            warnings.append(
                "mzML omits isolation lower/upper offsets; the full width from "
                f"Agilent AcqMethod.xml ({isolation_width_fallback_da:g} Da) is available "
                "as an audited analysis fallback."
            )
        else:
            warnings.append(
                "mzML omits isolation lower/upper offsets and no acquisition-method "
                "fallback is available; precursor deconvolution must remain disabled."
            )
    if centroid_count == 0 and profile_count == 0:
        warnings.append("Spectrum representation (centroid/profile) is not declared.")
    expected_hashes = {
        Path(str(name)).name: str(value).strip().lower()
        for name, value in (expected_source_sha1 or {}).items()
        if str(name).strip() and str(value).strip()
    }
    source_checksum_matches: dict[str, bool | None] = {}
    for name, expected_hash in expected_hashes.items():
        observed_hash = source_file_sha1.get(name)
        if observed_hash and observed_hash != expected_hash:
            source_checksum_matches[name] = False
            errors.append(
                f"mzML source checksum mismatch for {name}; output does not belong "
                "to the supplied raw directory."
            )
        elif observed_hash is None:
            source_checksum_matches[name] = None
            warnings.append(
                f"mzML does not expose a SHA-1 provenance checksum for {name}."
            )
        else:
            source_checksum_matches[name] = True

    report = {
        "status": "valid" if not errors else "invalid",
        "path": str(path),
        "size_bytes": after_parse.st_size,
        "sha256": mzml_sha256.hexdigest(),
        "declared_spectrum_count": declared_count,
        "actual_spectrum_count": actual_count,
        "ms_levels": {str(key): value for key, value in sorted(ms_levels.items())},
        "mass_spectrum_count": mass_spectrum_count,
        "auxiliary_spectrum_count": auxiliary_spectrum_count,
        "nonempty_auxiliary_spectrum_count": nonempty_auxiliary,
        "ms1_count": ms1_count,
        "ms2_count": ms2_count,
        "nonempty_ms1_count": nonempty_ms1,
        "nonempty_ms2_count": nonempty_ms2,
        "centroid_spectrum_count": centroid_count,
        "profile_spectrum_count": profile_count,
        "representation_by_ms_level": {
            str(key): value
            for key, value in sorted(representation_by_ms_level.items())
        },
        "polarities": polarities,
        "precursor_count": precursor_count,
        "precursor_spectrum_ref_count": precursor_ref_count,
        "isolation_target_count": isolation_target_count,
        "isolation_lower_offset_count": isolation_lower_count,
        "isolation_upper_offset_count": isolation_upper_count,
        "isolation_offsets_complete": offsets_complete,
        "isolation_width_fallback_da": isolation_width_fallback_da,
        "software": software,
        "source_file_formats": source_file_formats,
        "source_file_sha1": source_file_sha1,
        "source_checksum_matches": source_checksum_matches,
        "instrument": {
            "models": instrument_models,
            "serial_numbers": sorted(set(instrument_serial_numbers)),
            "ion_sources": ion_sources,
            "mass_analyzers": mass_analyzers,
            "detectors": detectors,
            "run": run_metadata,
        },
        "acquisition": {
            "collision_energies_by_ms_level": {
                str(level): sorted(values)
                for level, values in sorted(collision_energies.items())
            },
            "resolving_powers_by_ms_level": {
                str(level): sorted(values)
                for level, values in sorted(resolving_powers.items())
            },
            "isolation_widths_da": sorted(isolation_widths),
            "filter_string_examples_by_ms_level": {
                str(level): values
                for level, values in sorted(filter_string_examples.items())
                if values
            },
        },
        "errors": errors,
        "warnings": warnings,
    }
    if errors:
        raise VendorRawValidationError(" ".join(errors))
    return report


def _subprocess_error(command: list[str], result: subprocess.CompletedProcess) -> str:
    output = "\n".join(
        value.strip()
        for value in (result.stdout or "", result.stderr or "")
        if value.strip()
    )
    if len(output) > 6000:
        output = output[-6000:]
    executable = Path(command[0]).name
    return (
        f"{executable} conversion failed with exit code {result.returncode}."
        + (f"\n{output}" if output else "")
    )


def prepare_converter(
    *,
    accept_vendor_licenses: bool = False,
    allow_container_pull: bool = False,
    docker_image: str | None = None,
    timeout_seconds: float = 4 * 60 * 60,
) -> dict[str, Any]:
    """Prepare the converter independently of a raw-data upload.

    This makes the one-time Docker image pull explicit in the UI.  It never
    downloads the licence-bearing image unless the user has both accepted the
    vendor licences and allowed the pull.
    """

    image = (
        docker_image
        or os.environ.get("MSSYSTEM_PWIZ_DOCKER_IMAGE", "").strip()
        or _default_docker_image()
    )
    status = converter_status(docker_image=image)
    if status["local_msconvert"]["available"]:
        return {
            "ready": True,
            "downloaded": False,
            "backend": "local",
            "status": status,
        }
    docker = status["docker"]["path"]
    if not docker:
        raise ConverterUnavailableError(
            "Docker CLI is not installed. Install and start Docker Desktop."
        )
    if not status["docker"]["daemon_available"]:
        raise ConverterUnavailableError(
            "Docker Desktop/daemon is not running. Start Docker Desktop and retry."
        )
    if not accept_vendor_licenses:
        raise VendorLicenseNotAcceptedError(
            "Preparing the ProteoWizard Docker image requires explicit "
            "acceptance of the bundled vendor licences."
        )
    if status["docker"]["image_available"]:
        _verify_docker_wine_runner(docker, image)
        status["docker"]["wine_runner_verified"] = True
        return {
            "ready": True,
            "downloaded": False,
            "backend": "docker",
            "status": status,
        }
    if not allow_container_pull:
        raise ConverterUnavailableError(
            "The ProteoWizard image is not installed. Check the one-time "
            "download permission and retry."
        )

    command = [docker, "pull"]
    if _arm64_host():
        command.extend(["--platform", "linux/amd64"])
    command.append(image)
    pulled = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=timeout_seconds,
        check=False,
    )
    if pulled.returncode != 0:
        raise ConverterUnavailableError(_subprocess_error(command, pulled))
    refreshed = converter_status(docker_image=image)
    if not refreshed["docker"]["image_available"]:
        raise ConverterUnavailableError(
            "Docker reported a successful pull, but the ProteoWizard image "
            "cannot be inspected locally."
        )
    _verify_docker_wine_runner(docker, image)
    refreshed["docker"]["wine_runner_verified"] = True
    return {
        "ready": True,
        "downloaded": True,
        "backend": "docker",
        "status": refreshed,
    }


def _convert_vendor_source(
    raw_source: str | Path,
    output_path: str | Path,
    *,
    raw_report: dict[str, Any],
    backend: str = "auto",
    msconvert_bin: str | Path | None = None,
    docker_image: str | None = None,
    accept_vendor_licenses: bool = False,
    allow_container_pull: bool = False,
    overwrite: bool = False,
    timeout_seconds: float = 4 * 60 * 60,
    peak_picking_strategy: str = "auto",
    converter_status_report: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Convert one already-inspected vendor source and return a manifest.

    Conversion is deliberately network-free.  ``allow_container_pull`` remains
    in the public signature for compatibility with older callers, but image
    installation is only performed by :func:`prepare_converter` after an
    explicit user action.
    """
    started = time.time()
    raw_path = Path(raw_source).expanduser().resolve()
    source_kind = str(raw_report.get("source_kind") or "agilent_d")
    destination = Path(output_path).expanduser().resolve()
    if destination.suffix.lower() != ".mzml":
        raise VendorRawValidationError("Output path must end in .mzML")
    if destination.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    requested_peak_picking = str(
        peak_picking_strategy or "auto"
    ).strip().lower()
    acquisition = raw_report.get("acquisition") or {}
    resolved_peak_picking = resolve_peak_picking_strategy(
        raw_report,
        requested_peak_picking,
    )
    vendor_peak_picking = resolved_peak_picking == "vendor"

    selected = str(backend or "auto").strip().lower()
    if selected not in {"auto", "local", "docker"}:
        raise ValueError("backend must be auto, local, or docker")
    image = (
        docker_image
        or os.environ.get("MSSYSTEM_PWIZ_DOCKER_IMAGE", "").strip()
        or _default_docker_image()
    )
    preflight_reused = isinstance(converter_status_report, Mapping)
    status = (
        json.loads(json.dumps(converter_status_report))
        if preflight_reused
        else converter_status(msconvert_bin=msconvert_bin, docker_image=image)
    )
    local = status["local_msconvert"]["path"]
    if selected == "auto":
        if local:
            selected = "local"
        elif status["docker"]["cli_available"]:
            selected = "docker"
        else:
            raise ConverterUnavailableError(
                "No local msconvert or Docker CLI is available. Install ProteoWizard "
                "on Windows, or Docker Desktop on macOS/Linux."
            )

    temp_dir = Path(
        tempfile.mkdtemp(
            prefix=f".{source_kind}_convert_",
            dir=str(destination.parent),
        )
    )
    output_name = destination.name
    converter_identity: dict[str, Any] = {"backend": selected}
    command: list[str]
    try:
        if selected == "local":
            if not local:
                raise ConverterUnavailableError(
                    "Local msconvert was not found. Set MSSYSTEM_MSCONVERT_BIN or "
                    "use the Docker backend."
                )
            command = build_local_msconvert_command(
                local,
                raw_path,
                temp_dir,
                output_name,
                vendor_peak_picking=vendor_peak_picking,
            )
            version = _run_probe([local, "--version"], timeout=15)
            converter_identity.update(
                {
                    "executable": local,
                    "version": (version.stdout or version.stderr).strip(),
                }
            )
        else:
            if not accept_vendor_licenses:
                raise VendorLicenseNotAcceptedError(
                    "ProteoWizard's Docker image requires explicit acceptance of "
                    "the bundled vendor licences."
                )
            docker = status["docker"]["path"]
            if not docker:
                raise ConverterUnavailableError("Docker CLI is not installed.")
            if not status["docker"]["daemon_available"]:
                raise ConverterUnavailableError(
                    "Docker Desktop/daemon is not running. Start Docker Desktop and retry."
                )
            if not status["docker"]["image_available"]:
                raise ConverterUnavailableError(
                    f"ProteoWizard image is not installed locally: {image}. "
                    "Analysis never downloads images. Use the separate "
                    "Prepare / Download converter action once while online, "
                    "or import a previously saved image with docker load."
                )
            execution_image = _docker_execution_image(
                image,
                status["docker"],
                verified_status_reused=preflight_reused,
            )
            wine_runner = _docker_wine_runner()
            if not status["docker"].get("wine_runner_verified"):
                _verify_docker_wine_runner(docker, execution_image)
            command = build_docker_msconvert_command(
                docker,
                execution_image,
                raw_path,
                temp_dir,
                output_name,
                wine_runner=wine_runner,
                vendor_peak_picking=vendor_peak_picking,
            )
            converter_identity.update(
                {
                    "image": image,
                    "image_identity": status["docker"].get("image_identity", ""),
                    "execution_image": execution_image,
                    "platform": "linux/amd64",
                    "wine_runner": wine_runner,
                    "apple_silicon_compatibility_pinned": (
                        image == APPLE_SILICON_PWIZ_DOCKER_IMAGE
                    ),
                }
            )

        converter_started = time.perf_counter()
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
        if result.returncode != 0:
            raise VendorRawConversionError(_subprocess_error(command, result))
        converter_seconds = time.perf_counter() - converter_started

        produced = temp_dir / output_name
        if not produced.is_file():
            alternatives = sorted(temp_dir.glob("*.mzML")) + sorted(
                temp_dir.glob("*.mzml")
            )
            alternatives = list(dict.fromkeys(alternatives))
            if len(alternatives) == 1:
                produced = alternatives[0]
            else:
                raise VendorRawConversionError(
                    "msconvert reported success but did not produce exactly one mzML file."
                )

        expected_source_sha1 = {
            Path(relative).name: details.get("sha1_for_mzml_provenance", "")
            for relative, details in (
                raw_report.get("required_file_signatures") or {}
            ).items()
        }
        source_file_signature = raw_report.get("source_file_signature") or {}
        if source_file_signature:
            expected_source_sha1[raw_path.name] = source_file_signature.get(
                "sha1_for_mzml_provenance", ""
            )
        validation_started = time.perf_counter()
        validation = validate_mzml(
            produced,
            expected_ms1=acquisition.get("ms1_expected"),
            expected_ms2=bool(acquisition.get("ms2_expected", True)),
            expected_polarities=acquisition.get("polarities") or [],
            isolation_width_fallback_da=acquisition.get("isolation_width_da"),
            expected_source_sha1=expected_source_sha1,
        )
        validation_seconds = time.perf_counter() - validation_started
        finalise_started = time.perf_counter()
        if destination.exists() and overwrite:
            destination.unlink()
        produced.replace(destination)
        finalise_seconds = time.perf_counter() - finalise_started
        validation["path"] = str(destination)

        command_audit = []
        for token in command:
            audited = str(token)
            audited = audited.replace(str(temp_dir), "<CONVERSION_OUTPUT>")
            audited = audited.replace(str(raw_path.parent), "<RAW_PARENT>")
            audited = audited.replace(
                str(raw_path), f"<VENDOR_SOURCE>/{raw_path.name}"
            )
            audited = audited.replace(
                f"/input/{raw_path.name}", "<VENDOR_SOURCE>"
            )
            command_audit.append(audited)
        manifest = {
            "status": "completed",
            "format": f"{raw_report.get('format') or 'Vendor raw'} to mzML",
            "source_kind": source_kind,
            "source": raw_report,
            "output": validation,
            "converter": converter_identity,
            "conversion_policy": {
                "output_format": "mzML",
                "precision_bits": 64,
                "zlib": True,
                "ms_level_filter": None,
                "peak_picking_strategy_requested": requested_peak_picking,
                "peak_picking_strategy_resolved": resolved_peak_picking,
                "peak_picking_filter": (
                    "peakPicking true 1-" if vendor_peak_picking else None
                ),
                "peak_picking_filter_order": (
                    "first_filter" if vendor_peak_picking else None
                ),
                "downstream_openms_peak_picking": (
                    resolved_peak_picking == "openms"
                ),
                "preserves_ms1_and_ms2": True,
                "converter_preflight_reused": preflight_reused,
                "command_audit_without_source_path": command_audit,
            },
            "duration_seconds": round(time.time() - started, 3),
            "timing": {
                "converter_seconds": round(converter_seconds, 3),
                "stream_validation_and_sha256_seconds": round(
                    validation_seconds, 3
                ),
                "atomic_output_move_seconds": round(finalise_seconds, 3),
                "source_fingerprinting_passes_per_required_file": 1,
                "mzml_validation_passes": 1,
            },
        }
        return manifest
    except subprocess.TimeoutExpired as exc:
        raise VendorRawConversionError(
            f"Conversion exceeded {timeout_seconds:g} seconds and was stopped."
        ) from exc
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def convert_agilent_d(
    raw_dir: str | Path,
    output_path: str | Path,
    *,
    backend: str = "auto",
    msconvert_bin: str | Path | None = None,
    docker_image: str | None = None,
    accept_vendor_licenses: bool = False,
    allow_container_pull: bool = False,
    overwrite: bool = False,
    timeout_seconds: float = 4 * 60 * 60,
    peak_picking_strategy: str = "auto",
    converter_status_report: Mapping[str, Any] | None = None,
    raw_report: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Convert one Agilent MassHunter ``.d`` directory."""

    inspected = (
        _reuse_inspected_raw_report(
            raw_dir,
            raw_report,
            expected_source_kind="agilent_d",
        )
        if raw_report is not None
        else inspect_agilent_d(raw_dir)
    )
    return _convert_vendor_source(
        raw_dir,
        output_path,
        raw_report=inspected,
        backend=backend,
        msconvert_bin=msconvert_bin,
        docker_image=docker_image,
        accept_vendor_licenses=accept_vendor_licenses,
        allow_container_pull=allow_container_pull,
        overwrite=overwrite,
        timeout_seconds=timeout_seconds,
        peak_picking_strategy=peak_picking_strategy,
        converter_status_report=converter_status_report,
    )


def convert_thermo_raw(
    raw_file: str | Path,
    output_path: str | Path,
    *,
    expected_polarities: Iterable[str] | None = None,
    backend: str = "auto",
    msconvert_bin: str | Path | None = None,
    docker_image: str | None = None,
    accept_vendor_licenses: bool = False,
    allow_container_pull: bool = False,
    overwrite: bool = False,
    timeout_seconds: float = 4 * 60 * 60,
    peak_picking_strategy: str = "auto",
    converter_status_report: Mapping[str, Any] | None = None,
    raw_report: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Convert one Thermo Fisher ``.raw`` file with paired provenance."""

    inspected = (
        _reuse_inspected_raw_report(
            raw_file,
            raw_report,
            expected_source_kind="thermo_raw",
        )
        if raw_report is not None
        else inspect_thermo_raw(
            raw_file,
            expected_polarities=expected_polarities,
        )
    )
    return _convert_vendor_source(
        raw_file,
        output_path,
        raw_report=inspected,
        backend=backend,
        msconvert_bin=msconvert_bin,
        docker_image=docker_image,
        accept_vendor_licenses=accept_vendor_licenses,
        allow_container_pull=allow_container_pull,
        overwrite=overwrite,
        timeout_seconds=timeout_seconds,
        peak_picking_strategy=peak_picking_strategy,
        converter_status_report=converter_status_report,
    )


def convert_vendor_raw(
    source: str | Path,
    output_path: str | Path,
    *,
    expected_polarities: Iterable[str] | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Unified public entry point for supported proprietary raw formats."""

    path = Path(source).expanduser().resolve()
    if path.suffix.lower() == ".d":
        return convert_agilent_d(path, output_path, **kwargs)
    if path.suffix.lower() == THERMO_RAW_SUFFIX:
        return convert_thermo_raw(
            path,
            output_path,
            expected_polarities=expected_polarities,
            **kwargs,
        )
    raise VendorRawValidationError(
        "Unsupported vendor raw source; expected Agilent .d or Thermo .raw."
    )


def write_conversion_manifest(
    manifest: dict[str, Any], path: str | Path
) -> Path:
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    temporary.replace(destination)
    return destination
