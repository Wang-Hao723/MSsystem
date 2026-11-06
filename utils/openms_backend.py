"""Controlled subprocess access to the V7 PyOpenMS backend."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from functools import lru_cache
from pathlib import Path
from typing import Any

from utils.spectrum_contract import (
    is_auxiliary_detector_spectrum as _is_auxiliary_detector_spectrum,
)


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "openms_untargeted.py"


class OpenMSBackendError(RuntimeError):
    """OpenMS is unavailable or rejected an input."""


def _runtime_environment() -> dict[str, str]:
    env = dict(os.environ)
    conda_lib = Path(sys.executable).resolve().parent.parent / "lib"
    if conda_lib.is_dir():
        existing = env.get("DYLD_LIBRARY_PATH", "")
        values = [str(conda_lib)]
        if existing:
            values.append(existing)
        env["DYLD_LIBRARY_PATH"] = os.pathsep.join(values)
    return env


def _run_openms(
    command: list[str],
    *,
    manifest_path: Path,
    timeout_seconds: int = 3600,
) -> dict[str, Any]:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    argv = [
        sys.executable,
        str(RUNNER),
        "--json-output",
        str(manifest_path),
        *command,
    ]
    try:
        completed = subprocess.run(
            argv,
            cwd=str(ROOT),
            env=_runtime_environment(),
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise OpenMSBackendError(
            "PyOpenMS optional workflow exceeded its bounded runtime "
            f"({int(timeout_seconds)} s) and was terminated so library "
            "matching and the remaining queue can continue."
        ) from exc
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise OpenMSBackendError(
            "PyOpenMS backend produced no readable manifest. "
            f"exit={completed.returncode}; stderr={completed.stderr[-2000:]}"
        ) from exc
    payload["process"] = {
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-2000:],
        "stderr_tail": completed.stderr[-2000:],
    }
    diagnostics = []
    if "Feature's charge is 0" in completed.stdout:
        diagnostics.append(
            {
                "code": "openms_unassigned_feature_charge",
                "severity": "information",
                "interpretation": (
                    "OpenMS left some features without an ion-identity charge; "
                    "they remain measured LC-MS features but receive no adduct "
                    "or neutral-mass claim."
                ),
            }
        )
    if "gapped, even-numbered charge ladders" in completed.stdout:
        diagnostics.append(
            {
                "code": "openms_charge_ladder_diagnostic",
                "severity": "warning",
                "interpretation": (
                    "OpenMS emitted its generic charge-ladder diagnostic. "
                    "Adduct groups are retained only as ion-identity hypotheses "
                    "and never upgraded to compound identifications."
                ),
            }
        )
    if diagnostics:
        payload["algorithm_diagnostics"] = diagnostics
    if completed.returncode != 0 or payload.get("status") == "error":
        raise OpenMSBackendError(
            str(payload.get("error") or "PyOpenMS backend failed")
        )
    return payload


@lru_cache(maxsize=1)
def probe_openms() -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="mssystem_openms_probe_") as directory:
        manifest = Path(directory) / "probe.json"
        try:
            return _run_openms(["probe"], manifest_path=manifest, timeout_seconds=60)
        except Exception as exc:
            return {
                "status": "unavailable",
                "backend": "PyOpenMS",
                "error": str(exc),
                "python": sys.executable,
            }


def inspect_mzml_representations(path: str | Path) -> dict[str, Any]:
    """Count representation CV terms per MS level, excluding detector traces."""

    from pyteomics import mzml

    target = Path(path).expanduser().resolve()
    counts = {
        "ms1": {"centroid": 0, "profile": 0, "unknown": 0, "spectra": 0},
        "ms2": {"centroid": 0, "profile": 0, "unknown": 0, "spectra": 0},
        "auxiliary": {"centroid": 0, "profile": 0, "unknown": 0, "spectra": 0},
    }
    examples = {"profile": [], "unknown": []}
    with mzml.read(str(target)) as reader:
        for index, spectrum in enumerate(reader):
            level = spectrum.get("ms level")
            if _is_auxiliary_detector_spectrum(spectrum):
                key = "auxiliary"
            else:
                key = (
                    "ms1" if level == 1 else "ms2" if level == 2 else "auxiliary"
                )
            bucket = counts[key]
            bucket["spectra"] += 1
            if "centroid spectrum" in spectrum:
                representation = "centroid"
            elif "profile spectrum" in spectrum:
                representation = "profile"
            else:
                representation = "unknown"
            bucket[representation] += 1
            if (
                representation in examples
                and len(examples[representation]) < 10
                and key in {"ms1", "ms2"}
            ):
                examples[representation].append(
                    str(spectrum.get("id") or f"index={index}")
                )
    mass_spectra = counts["ms1"]["spectra"] + counts["ms2"]["spectra"]
    return {
        "schema_version": 1,
        "input": str(target),
        "decision_scope": (
            "per-spectrum MS1/MS2 CV terms; optical detector CV terms take "
            "precedence over synthesized MS levels"
        ),
        "counts": counts,
        "mass_spectra": mass_spectra,
        "profile_mass_spectra": (
            counts["ms1"]["profile"] + counts["ms2"]["profile"]
        ),
        "unknown_mass_spectra": (
            counts["ms1"]["unknown"] + counts["ms2"]["unknown"]
        ),
        "examples": examples,
    }


def prepare_mzml_for_analysis(
    path: str | Path,
    artifact_dir: str | Path,
    *,
    policy: str = "auto",
    peak_picker: str = "hires",
    signal_to_noise: float = 1.0,
) -> tuple[Path, dict[str, Any]]:
    """Peak-pick only profile MS1/MS2 spectra and preserve centroid spectra."""

    source = Path(path).expanduser().resolve()
    is_mzxml = source.suffix.lower() == ".mzxml"
    if is_mzxml:
        # mzXML does not carry the mzML representation CV contract in the same
        # form.  OpenMS performs the per-spectrum signal-shape decision while
        # normalizing the file to mzML.
        audit = {
            "schema_version": 1,
            "input": str(source),
            "decision_scope": "OpenMS per-spectrum signal-shape detection during mzXML normalization",
            "counts": {},
            "mass_spectra": None,
            "profile_mass_spectra": None,
            "unknown_mass_spectra": None,
        }
    else:
        audit = inspect_mzml_representations(source)
    profile_count = int(audit["profile_mass_spectra"] or 0)
    unknown_count = int(audit["unknown_mass_spectra"] or 0)
    normalized_policy = str(policy or "auto").strip().lower()
    if normalized_policy not in {"auto", "openms", "reject", "assume_centroid"}:
        raise ValueError(
            "profile_peak_picking must be auto, openms, reject, or assume_centroid"
        )

    if not is_mzxml and profile_count == 0 and unknown_count == 0:
        auxiliary_count = int(
            ((audit.get("counts") or {}).get("auxiliary") or {}).get("spectra")
            or 0
        )
        return source, {
            "status": "not_required_all_mass_spectra_centroid",
            "policy": normalized_policy,
            "representation_audit": audit,
            "auxiliary_spectra_excluded_from_output": auxiliary_count,
            "backend": None,
            "input": str(source),
            "output": str(source),
        }
    if normalized_policy == "reject":
        raise OpenMSBackendError(
            f"Input contains {profile_count} profile and {unknown_count} "
            "representation-unknown MS1/MS2 spectra; peak picking is required."
        )
    if normalized_policy == "assume_centroid":
        if is_mzxml:
            raise OpenMSBackendError(
                "mzXML must be normalized to mzML by OpenMS before analysis; "
                "assume_centroid is available only for mzML peak lists missing CV terms."
            )
        if profile_count:
            raise OpenMSBackendError(
                "assume_centroid cannot override explicit profile-spectrum CV terms"
            )
        auxiliary_count = int(
            ((audit.get("counts") or {}).get("auxiliary") or {}).get("spectra")
            or 0
        )
        return source, {
            "status": "accepted_unknown_as_centroid_by_explicit_user_policy",
            "policy": normalized_policy,
            "representation_audit": audit,
            "auxiliary_spectra_excluded_from_output": auxiliary_count,
            "backend": None,
            "input": str(source),
            "output": str(source),
            "warning": (
                "Representation CV terms were missing; the user explicitly chose "
                "to treat the peak lists as centroid data."
            ),
        }

    probe = probe_openms()
    if probe.get("status") != "ok":
        raise OpenMSBackendError(
            "Profile/unknown spectra require OpenMS peak picking, but PyOpenMS "
            f"is unavailable: {probe.get('error')}"
        )
    directory = Path(artifact_dir).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    output = directory / f"{source.stem}.openms-centroided.mzML"
    manifest_path = directory / f"{source.stem}.openms-peak-picking.json"
    payload = _run_openms(
        [
            "peak-pick",
            "--input",
            str(source),
            "--output",
            str(output),
            "--strategy",
            peak_picker,
            "--signal-to-noise",
            str(float(signal_to_noise)),
        ],
        manifest_path=manifest_path,
    )
    payload.update(
        {
            "policy": normalized_policy,
            "representation_audit_before": audit,
            "normalized_from_mzxml": is_mzxml,
            "manifest": str(manifest_path),
        }
    )
    payload["representation_audit_after"] = inspect_mzml_representations(output)
    if (
        payload["representation_audit_after"]["profile_mass_spectra"]
        or payload["representation_audit_after"]["unknown_mass_spectra"]
    ):
        raise OpenMSBackendError(
            "OpenMS output still contains profile/unknown MS1/MS2 spectra"
        )
    return output, payload


def detect_untargeted_features(
    path: str | Path,
    artifact_dir: str | Path,
    *,
    ion_mode: str,
    peak_picker: str = "hires",
    signal_to_noise: float = 1.0,
    mass_error_ppm: float = 10.0,
    noise_threshold: float = 1000.0,
    chrom_peak_snr: float = 3.0,
    min_trace_length: float = 5.0,
    min_fwhm: float = 1.0,
    max_fwhm: float = 60.0,
    remove_single_traces: bool = False,
    charge_low: int = 1,
    charge_high: int = 1,
    group_adducts: bool = True,
    max_adduct_grouping_features: int = 6_000,
    timeout_seconds: int = 300,
) -> dict[str, Any]:
    directory = Path(artifact_dir).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    source = Path(path).expanduser().resolve()
    feature_xml = directory / f"{source.stem}.featureXML"
    groups_xml = directory / f"{source.stem}.adduct-groups.consensusXML"
    manifest = directory / f"{source.stem}.untargeted-features.json"
    command = [
        "features",
        "--input",
        str(source),
        "--output-features",
        str(feature_xml),
        "--output-adduct-groups",
        str(groups_xml),
        "--peak-picker",
        peak_picker,
        "--signal-to-noise",
        str(float(signal_to_noise)),
        "--mass-error-ppm",
        str(float(mass_error_ppm)),
        "--noise-threshold",
        str(float(noise_threshold)),
        "--chrom-peak-snr",
        str(float(chrom_peak_snr)),
        "--min-trace-length",
        str(float(min_trace_length)),
        "--min-fwhm",
        str(float(min_fwhm)),
        "--max-fwhm",
        str(float(max_fwhm)),
        "--charge-low",
        str(int(charge_low)),
        "--charge-high",
        str(int(charge_high)),
        "--ion-mode",
        str(ion_mode).upper(),
        "--max-adduct-grouping-features",
        str(max(0, int(max_adduct_grouping_features))),
    ]
    if remove_single_traces:
        command.append("--remove-single-traces")
    if group_adducts:
        command.append("--group-adducts")
    payload = _run_openms(
        command,
        manifest_path=manifest,
        timeout_seconds=max(30, min(3600, int(timeout_seconds))),
    )
    payload["manifest"] = str(manifest)
    payload["feature_xml"] = str(feature_xml)
    payload["adduct_groups_xml"] = (
        str(groups_xml) if groups_xml.exists() else None
    )
    return payload


def align_feature_maps(
    feature_xml_paths: list[str | Path],
    output_dir: str | Path,
    *,
    rt_tolerance_seconds: float = 20.0,
    mz_tolerance: float = 10.0,
    mz_unit: str = "ppm",
) -> dict[str, Any]:
    directory = Path(output_dir).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    consensus = directory / "aligned.consensusXML"
    manifest = directory / "rt-alignment.json"
    payload = _run_openms(
        [
            "align",
            "--inputs",
            *[str(Path(path).expanduser().resolve()) for path in feature_xml_paths],
            "--output-consensus",
            str(consensus),
            "--rt-tolerance-seconds",
            str(float(rt_tolerance_seconds)),
            "--mz-tolerance",
            str(float(mz_tolerance)),
            "--mz-unit",
            mz_unit,
        ],
        manifest_path=manifest,
    )
    payload["manifest"] = str(manifest)
    return payload
