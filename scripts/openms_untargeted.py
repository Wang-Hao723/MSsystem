#!/usr/bin/env python3
"""Isolated PyOpenMS backend for representation-aware LC-MS processing.

This process boundary is deliberate.  On the current macOS/Conda installation
PyOpenMS needs Conda's shared libraries ahead of Homebrew's libraries.  The web
server invokes this script with a controlled environment, records the exact
OpenMS version, and never silently substitutes the legacy heuristic picker.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import traceback
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any

import pyopenms as oms


PROFILE = int(oms.SpectrumSettings.SpectrumType.PROFILE)
CENTROID = int(oms.SpectrumSettings.SpectrumType.CENTROID)
UNKNOWN = int(oms.SpectrumSettings.SpectrumType.UNKNOWN)

_AUXILIARY_SPECTRUM_ACCESSIONS = {"MS:1000804"}
_AUXILIARY_SPECTRUM_NAMES = {
    "electromagnetic radiation spectrum",
    "absorption spectrum",
    "emission spectrum",
    "wavelength array",
}

DEFAULT_POS_ADDUCTS = "H:+:0.4,Na:+:0.25,NH4:+:0.25,K:+:0.1,H-2O-1:0:0.05"
DEFAULT_NEG_ADDUCTS = "H-1:-:0.6,Cl:-:0.2,H-3O-1:-:0.1,CH2O2H-1:-:0.1"
# MetaboliteFeatureDeconvolution is a graph/ILP problem, not a linear pass over
# features.  A real 7,656-feature POS map exceeded 300 s while the same map's
# MassTraceDetection -> FeatureFindingMetabo path completed in 2.24 s.  A
# 5,959-feature production map completed the ILP quickly, so 6,000 is the
# measured high-throughput boundary for this runtime.  Treat it as a hard
# safety ceiling even when an older queued task requests the former 10,000
# limit (or 0, which historically disabled the guard).
MAX_SAFE_ADDUCT_GROUPING_FEATURES = 6_000
DEFAULT_MAX_ADDUCT_GROUPING_FEATURES = MAX_SAFE_ADDUCT_GROUPING_FEATURES


def _effective_adduct_grouping_limit(configured_max_features: int) -> int:
    """Return the enforced high-throughput ceiling for the optional ILP."""

    configured = max(0, int(configured_max_features))
    if configured == 0:
        return MAX_SAFE_ADDUCT_GROUPING_FEATURES
    return min(configured, MAX_SAFE_ADDUCT_GROUPING_FEATURES)


def _local_name(tag: str) -> str:
    return str(tag).rsplit("}", 1)[-1]


def _explicit_mzml_ms_level_ids(
    path: str | Path,
) -> tuple[dict[int, set[str]], int]:
    """Return native IDs carrying an explicit mzML MS-level CV term.

    OpenMS maps spectra without ``MS:1000511`` to MS level 1.  That behavior is
    unsafe for multi-detector files because DAD/UV spectra have no MS level but
    can otherwise enter metabolomics feature detection as survey scans.  Read
    the source contract directly and allow only explicit MS1/MS2 native IDs.
    """

    levels: dict[int, set[str]] = {1: set(), 2: set()}
    reference_levels: dict[str, int] = {}
    reference_auxiliary_groups: set[str] = set()
    auxiliary = 0
    for event, element in ET.iterparse(str(path), events=("end",)):
        local = _local_name(element.tag)
        if local == "referenceableParamGroup":
            level = None
            is_auxiliary_detector = False
            for child in element.iter():
                if _local_name(child.tag) != "cvParam":
                    continue
                accession = str(child.attrib.get("accession") or "")
                name = str(child.attrib.get("name") or "").strip().lower()
                if (
                    accession in _AUXILIARY_SPECTRUM_ACCESSIONS
                    or name in _AUXILIARY_SPECTRUM_NAMES
                ):
                    is_auxiliary_detector = True
                if accession == "MS:1000511":
                    try:
                        level = int(float(child.attrib.get("value", "")))
                    except (TypeError, ValueError):
                        level = None
            group_id = str(element.attrib.get("id") or "")
            if group_id and level in {1, 2}:
                reference_levels[group_id] = level
            if group_id and is_auxiliary_detector:
                reference_auxiliary_groups.add(group_id)
            element.clear()
            continue
        if local != "spectrum":
            continue
        level = None
        is_auxiliary_detector = False
        for child in element.iter():
            child_local = _local_name(child.tag)
            if child_local == "cvParam":
                accession = str(child.attrib.get("accession") or "")
                name = str(child.attrib.get("name") or "").strip().lower()
                if (
                    accession in _AUXILIARY_SPECTRUM_ACCESSIONS
                    or name in _AUXILIARY_SPECTRUM_NAMES
                ):
                    is_auxiliary_detector = True
                if accession == "MS:1000511":
                    try:
                        level = int(float(child.attrib.get("value", "")))
                    except (TypeError, ValueError):
                        level = None
            if child_local == "referenceableParamGroupRef":
                group_ref = str(child.attrib.get("ref") or "")
                referenced = reference_levels.get(group_ref)
                if referenced in {1, 2}:
                    level = referenced
                if group_ref in reference_auxiliary_groups:
                    is_auxiliary_detector = True
        native_id = str(element.attrib.get("id") or "")
        if not is_auxiliary_detector and level in {1, 2} and native_id:
            levels[level].add(native_id)
        else:
            auxiliary += 1
        element.clear()
    return levels, auxiliary


def _retain_explicit_mass_spectra(
    experiment: oms.MSExperiment,
    allowed_native_ids: dict[int, set[str]],
) -> dict[str, Any]:
    """Remove auxiliary detector spectra before an analysis mzML is stored.

    PyOpenMS maps spectra without an explicit ``MS:1000511`` term to MS level
    1 when it loads and stores an experiment.  In LC-MS files containing DAD or
    UV spectra, that turns optical spectra into apparent profile MS1 spectra.
    The source mzML contract parsed by :func:`_explicit_mzml_ms_level_ids` is
    authoritative, so retain only its explicit MS1/MS2 native IDs and verify
    that OpenMS neither changed nor dropped any of them.
    """

    expected_level_by_id: dict[str, int] = {}
    for level in (1, 2):
        for native_id in allowed_native_ids.get(level, set()):
            previous = expected_level_by_id.setdefault(str(native_id), level)
            if previous != level:
                raise RuntimeError(
                    f"Duplicate mzML native ID {native_id!r} has conflicting "
                    f"MS levels {previous} and {level}"
                )

    retained = []
    retained_by_level = {1: 0, 2: 0}
    observed_native_ids: set[str] = set()
    excluded = 0
    for spectrum in experiment:
        native_id = str(spectrum.getNativeID())
        expected_level = expected_level_by_id.get(native_id)
        if expected_level is None:
            excluded += 1
            continue
        if native_id in observed_native_ids:
            raise RuntimeError(
                f"OpenMS loaded duplicate explicit mass-spectrum native ID "
                f"{native_id!r}"
            )
        actual_level = int(spectrum.getMSLevel())
        if actual_level != expected_level:
            raise RuntimeError(
                f"OpenMS changed explicit mzML native ID {native_id!r} from "
                f"MS{expected_level} to MS{actual_level}"
            )
        retained.append(spectrum)
        observed_native_ids.add(native_id)
        retained_by_level[expected_level] += 1

    missing = sorted(set(expected_level_by_id) - observed_native_ids)
    if missing:
        examples = ", ".join(repr(item) for item in missing[:5])
        raise RuntimeError(
            f"OpenMS dropped {len(missing)} explicit MS1/MS2 spectra while "
            f"loading mzML; examples: {examples}"
        )

    experiment.setSpectra(retained)
    return {
        "policy": "retain_only_source_explicit_ms1_ms2_native_ids",
        "retained_spectra": len(retained),
        "retained_by_ms_level": {
            str(level): retained_by_level[level] for level in (1, 2)
        },
        "excluded_auxiliary_spectra": excluded,
        "missing_explicit_mass_spectra": 0,
    }


def _write_json(path: str | Path, payload: dict[str, Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str) + "\n",
        encoding="utf-8",
    )


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _jsonable(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    try:
        return value.toString().decode("utf-8", errors="replace")
    except Exception:
        return str(value)


def _single_text(value: Any) -> str:
    """Return one stable text value from OpenMS scalar/list metadata."""

    normalized = _jsonable(value)
    if isinstance(normalized, list):
        normalized = normalized[0] if len(normalized) == 1 else ""
    return str(normalized or "").strip()


def _adduct_grouping_decision(
    raw_feature_count: int,
    requested: bool,
    max_features: int,
) -> tuple[bool, str]:
    """Bound the optional ILP before it can block an analysis queue.

    ``MetaboliteFeatureDeconvolution`` solves an integer program whose runtime
    depends on both feature count and local graph density.  It is useful
    ion-identity evidence, but it does not participate in MS/MS library
    matching or High-confidence calibration.  Large maps therefore retain the
    measured FeatureFindingMetabo map and skip only the optional ILP.
    """

    count = max(0, int(raw_feature_count))
    limit = _effective_adduct_grouping_limit(max_features)
    if count == 0:
        return False, "not_run_no_features"
    if not requested:
        return False, "disabled_preserved_feature_map"
    if count > limit:
        return False, "skipped_feature_limit_preserved_feature_map"
    return True, "running"


def _canonical_adduct_label(meta: dict[str, Any], charge: int) -> tuple[str | None, str | None]:
    """Translate OpenMS formula-style adduct metadata to common MS notation.

    ``MetaboliteFeatureDeconvolution`` records the machine-readable formula in
    ``dc_charge_adducts``.  Its additional humanized label can render formate
    as ``[M+C+H+2O]-`` and a water loss as ``[M-3H-O]-``.  Preserve the raw
    value for audit, but present the standard analytical notation.
    """

    raw = _single_text(
        meta.get("dc_charge_adducts")
        or meta.get("dc_adduct")
        or meta.get("adduct")
    )
    compact = raw.replace(" ", "")
    polarity_sign = "-" if charge < 0 else "+"
    canonical = {
        # Negative-mode defaults.
        "H-1": "[M-H]-",
        "Cl": "[M+Cl]-",
        "Cl1": "[M+Cl]-",
        "H-3O-1": "[M-H2O-H]-",
        "C1H1O2": "[M+HCOO]-",
        "CH1O2": "[M+HCOO]-",
        "CHO2": "[M+HCOO]-",
        "CH2O2H-1": "[M+HCOO]-",
        # Positive-mode defaults.
        "H": "[M+H]+",
        "H1": "[M+H]+",
        "Na": "[M+Na]+",
        "Na1": "[M+Na]+",
        "NH4": "[M+NH4]+",
        "N1H4": "[M+NH4]+",
        "K": "[M+K]+",
        "K1": "[M+K]+",
        "H-2O-1": "[M-H2O]+",
    }.get(compact)
    if canonical:
        return canonical, raw or None

    displayed = meta.get("adducts")
    if isinstance(displayed, list) and len(displayed) == 1:
        displayed = displayed[0]
    fallback = _single_text(displayed)
    if fallback:
        return fallback, raw or fallback
    if compact.startswith("[M") and compact.endswith(("+", "-")):
        return compact, raw or compact
    if raw:
        return f"{raw} ({polarity_sign} mode; OpenMS formula)", raw
    return None, None


def _spectrum_representation(spectrum: oms.MSSpectrum) -> tuple[str, str]:
    """Return ``(representation, evidence)`` for one spectrum."""

    explicit = int(spectrum.getType())
    if explicit == CENTROID:
        return "centroid", "mzml_cv"
    if explicit == PROFILE:
        return "profile", "mzml_cv"
    # OpenMS' own per-spectrum signal-shape detector is the only fallback for
    # an mzML that omitted both CV terms.  It is recorded as inferred evidence.
    # The Python binding exposes only getType(), unlike the C++ overload with
    # an estimate flag. Use the same OpenMS estimator through its bound API.
    inferred = int(oms.PeakTypeEstimator().estimateType(spectrum))
    if inferred == CENTROID:
        return "centroid", "openms_signal_shape"
    if inferred == PROFILE:
        return "profile", "openms_signal_shape"
    return "unknown", "openms_signal_shape"


def _picker(strategy: str, signal_to_noise: float):
    if strategy == "iterative":
        picker = oms.PeakPickerIterative()
        params = picker.getParameters()
        params.setValue("signal_to_noise_", float(signal_to_noise))
    else:
        picker = oms.PeakPickerHiRes()
        params = picker.getParameters()
        params.setValue("signal_to_noise", float(signal_to_noise))
    picker.setParameters(params)
    return picker


def _pick_profile_spectra(
    experiment: oms.MSExperiment,
    *,
    strategy: str,
    signal_to_noise: float,
    ms_levels: set[int] = {1, 2},
    allowed_native_ids: dict[int, set[str]] | None = None,
) -> dict[str, Any]:
    picker = _picker(strategy, signal_to_noise)
    counts: dict[str, dict[str, int]] = {
        str(level): {
            "spectra": 0,
            "explicit_centroid": 0,
            "explicit_profile": 0,
            "inferred_centroid": 0,
            "inferred_profile": 0,
            "unknown": 0,
            "profile_picked": 0,
            "input_points": 0,
            "output_points": 0,
        }
        for level in sorted(ms_levels)
    }
    for index in range(experiment.getNrSpectra()):
        spectrum = experiment.getSpectrum(index)
        level = int(spectrum.getMSLevel())
        if level not in ms_levels:
            continue
        if (
            allowed_native_ids is not None
            and str(spectrum.getNativeID())
            not in allowed_native_ids.get(level, set())
        ):
            continue
        bucket = counts[str(level)]
        bucket["spectra"] += 1
        bucket["input_points"] += int(spectrum.size())
        representation, evidence = _spectrum_representation(spectrum)
        bucket[f"{'explicit' if evidence == 'mzml_cv' else 'inferred'}_{representation}"] = (
            bucket.get(
                f"{'explicit' if evidence == 'mzml_cv' else 'inferred'}_{representation}",
                0,
            )
            + 1
        )
        if representation == "profile":
            picked = oms.MSSpectrum()
            picker.pick(spectrum, picked)
            experiment[index] = picked
            bucket["profile_picked"] += 1
            bucket["output_points"] += int(picked.size())
        elif representation == "centroid":
            if evidence == "openms_signal_shape":
                # Persist the estimator's decision: otherwise an inferred
                # centroid is written back as UNKNOWN and fails the mzML
                # representation check despite its unchanged discrete peaks.
                spectrum.setType(CENTROID)
                experiment[index] = spectrum
            bucket["output_points"] += int(spectrum.size())
        else:
            bucket["unknown"] += 1
            bucket["output_points"] += int(spectrum.size())
    return counts


def command_probe(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "status": "ok",
        "backend": "PyOpenMS",
        "openms_version": str(oms.__version__),
        "python": sys.executable,
        "algorithms": [
            "PeakPickerHiRes",
            "PeakPickerIterative",
            "MassTraceDetection",
            "ElutionPeakDetection",
            "FeatureFindingMetabo",
            "MetaboliteFeatureDeconvolution",
            "MapAlignmentAlgorithmPoseClustering",
            "FeatureGroupingAlgorithmQT",
        ],
    }


def command_peak_pick(args: argparse.Namespace) -> dict[str, Any]:
    experiment = oms.MSExperiment()
    input_path = Path(args.input)
    allowed_native_ids = None
    auxiliary_spectra = None
    output_filter = None
    if input_path.suffix.lower() == ".mzxml":
        oms.MzXMLFile().load(str(input_path), experiment)
    else:
        allowed_native_ids, auxiliary_spectra = _explicit_mzml_ms_level_ids(
            input_path
        )
        oms.MzMLFile().load(str(input_path), experiment)
        output_filter = _retain_explicit_mass_spectra(
            experiment, allowed_native_ids
        )
    counts = _pick_profile_spectra(
        experiment,
        strategy=args.strategy,
        signal_to_noise=args.signal_to_noise,
        allowed_native_ids=allowed_native_ids,
    )
    unresolved = sum(level["unknown"] for level in counts.values())
    if unresolved:
        raise RuntimeError(
            f"OpenMS could not determine representation for {unresolved} MS1/MS2 spectra"
        )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    oms.MzMLFile().store(str(output), experiment)
    return {
        "status": "ok",
        "backend": f"PyOpenMS {oms.__version__}",
        "algorithm": (
            "PeakPickerIterative" if args.strategy == "iterative"
            else "PeakPickerHiRes"
        ),
        "signal_to_noise": float(args.signal_to_noise),
        "input": str(input_path.resolve()),
        "output": str(output.resolve()),
        "counts_by_ms_level": counts,
        "auxiliary_spectra_excluded_from_peak_picking": auxiliary_spectra,
        "analysis_output_filter": output_filter,
        "auxiliary_spectra_excluded_from_output": (
            output_filter.get("excluded_auxiliary_spectra")
            if output_filter is not None
            else None
        ),
        "profile_spectra_picked": sum(
            level["profile_picked"] for level in counts.values()
        ),
    }


def _feature_meta(feature: oms.Feature) -> dict[str, Any]:
    keys: list[Any] = []
    try:
        feature.getKeys(keys)
    except Exception:
        keys = []
    output = {}
    for raw_key in keys:
        key = _jsonable(raw_key)
        if not key:
            continue
        try:
            output[str(key)] = _jsonable(feature.getMetaValue(raw_key))
        except Exception:
            continue
    return output


def _serialize_features(
    feature_map: oms.FeatureMap,
    adduct_groups: oms.ConsensusMap | None = None,
) -> list[dict[str, Any]]:
    group_by_uid: dict[int, str] = {}
    if adduct_groups is not None:
        for index, consensus_feature in enumerate(adduct_groups):
            group_id = f"AG{index + 1:06d}"
            for handle in consensus_feature.getFeatureList():
                group_by_uid[int(handle.getUniqueId())] = group_id

    feature_meta = [_feature_meta(feature) for feature in feature_map]
    raw_group_counts = Counter(
        str(meta.get("Group") or "")
        for meta in feature_meta
        if str(meta.get("Group") or "")
    )
    multi_group_ids = {
        raw_group: f"AG{index + 1:06d}"
        for index, raw_group in enumerate(
            sorted(
                raw_group
                for raw_group, count in raw_group_counts.items()
                if count > 1
            )
        )
    }

    rows = []
    for index, (feature, meta) in enumerate(zip(feature_map, feature_meta)):
        uid = int(feature.getUniqueId())
        rt_seconds = _finite(feature.getRT())
        width_seconds = max(0.0, _finite(feature.getWidth()) or 0.0)
        rt_start = (
            rt_seconds - width_seconds / 2.0 if rt_seconds is not None else None
        )
        rt_end = (
            rt_seconds + width_seconds / 2.0 if rt_seconds is not None else None
        )
        charge = int(feature.getCharge())
        adduct, adduct_raw = _canonical_adduct_label(meta, charge)
        neutral_mass = (
            _finite(meta.get("dc_neutral_mass"))
            or _finite(meta.get("neutral_mass"))
        )
        if neutral_mass is None and charge:
            adduct_mass = _finite(meta.get("dc_charge_adduct_mass"))
            if adduct_mass is not None:
                neutral_mass = (
                    abs(charge) * float(feature.getMZ()) - adduct_mass
                )
        raw_group = str(meta.get("Group") or "")
        isotope_trace_count = max(
            1,
            int(_finite(meta.get("num_of_masstraces")) or 1),
            1 + len(feature.getSubordinates()),
        )
        rows.append(
            {
                "feature_id": f"LF{index + 1:07d}",
                "openms_unique_id": uid,
                "mz": _finite(feature.getMZ()),
                "rt_seconds": rt_seconds,
                "rt_minutes": (
                    rt_seconds / 60.0 if rt_seconds is not None else None
                ),
                "rt_start_seconds": rt_start,
                "rt_end_seconds": rt_end,
                "rt_start_minutes": (
                    rt_start / 60.0 if rt_start is not None else None
                ),
                "rt_end_minutes": (
                    rt_end / 60.0 if rt_end is not None else None
                ),
                "width_seconds": width_seconds,
                "intensity": _finite(feature.getIntensity()),
                "charge": charge,
                "quality": _finite(feature.getOverallQuality()),
                "isotope_trace_count": isotope_trace_count,
                "adduct": _jsonable(adduct),
                "adduct_raw_openms": _jsonable(adduct_raw),
                "neutral_mass": neutral_mass,
                "adduct_group_id": (
                    multi_group_ids.get(raw_group)
                    if raw_group
                    else group_by_uid.get(uid)
                ),
                "openms_meta": meta,
            }
        )
    return rows


def command_features(args: argparse.Namespace) -> dict[str, Any]:
    experiment = oms.MSExperiment()
    input_path = Path(args.input)
    allowed_native_ids, auxiliary_spectra = _explicit_mzml_ms_level_ids(
        input_path
    )
    oms.MzMLFile().load(str(input_path), experiment)
    representation = _pick_profile_spectra(
        experiment,
        strategy=args.peak_picker,
        signal_to_noise=args.signal_to_noise,
        ms_levels={1},
        allowed_native_ids=allowed_native_ids,
    )

    ms1 = oms.MSExperiment()
    for spectrum in experiment:
        if (
            int(spectrum.getMSLevel()) == 1
            and str(spectrum.getNativeID()) in allowed_native_ids[1]
        ):
            ms1.addSpectrum(spectrum)
    ms1.sortSpectra(True)
    ms1_count = int(ms1.getNrSpectra())

    feature_map = oms.FeatureMap()
    mass_traces = []
    split_traces = []
    if ms1_count:
        detector = oms.MassTraceDetection()
        params = detector.getDefaults()
        params.setValue("mass_error_ppm", float(args.mass_error_ppm))
        params.setValue("noise_threshold_int", float(args.noise_threshold))
        params.setValue("chrom_peak_snr", float(args.chrom_peak_snr))
        params.setValue("min_trace_length", float(args.min_trace_length))
        detector.setParameters(params)
        detector.run(ms1, mass_traces, 0)

        elution = oms.ElutionPeakDetection()
        params = elution.getDefaults()
        params.setValue("width_filtering", args.width_filtering)
        params.setValue("chrom_peak_snr", float(args.chrom_peak_snr))
        params.setValue("min_fwhm", float(args.min_fwhm))
        params.setValue("max_fwhm", float(args.max_fwhm))
        elution.setParameters(params)
        elution.detectPeaks(mass_traces, split_traces)

        finder = oms.FeatureFindingMetabo()
        params = finder.getDefaults()
        params.setValue("isotope_filtering_model", args.isotope_model)
        params.setValue(
            "remove_single_traces", "true" if args.remove_single_traces else "false"
        )
        params.setValue("charge_lower_bound", int(args.charge_low))
        params.setValue("charge_upper_bound", int(args.charge_high))
        params.setValue("report_convex_hulls", "true")
        finder.setParameters(params)
        chromatograms = []
        finder.run(split_traces, feature_map, chromatograms)
        feature_map.setUniqueIds()

    raw_feature_count = int(feature_map.size())
    feature_finder_charge_counts = dict(
        sorted(Counter(int(feature.getCharge()) for feature in feature_map).items())
    )
    configured_deconvolution_limit = int(args.max_adduct_grouping_features)
    effective_deconvolution_limit = _effective_adduct_grouping_limit(
        configured_deconvolution_limit
    )
    run_deconvolution, deconvolution_status = _adduct_grouping_decision(
        raw_feature_count,
        bool(args.group_adducts),
        configured_deconvolution_limit,
    )
    decharged = feature_map
    groups = oms.ConsensusMap()
    edges = oms.ConsensusMap()
    deconvolution_error = (
        (
            "Optional MetaboliteFeatureDeconvolution skipped: "
            f"{raw_feature_count} features exceed the effective high-throughput "
            f"safety limit of {effective_deconvolution_limit} "
            f"(requested {configured_deconvolution_limit}). The measured "
            "FeatureFindingMetabo map is preserved for MS2 attachment."
        )
        if deconvolution_status == "skipped_feature_limit_preserved_feature_map"
        else ""
    )
    if run_deconvolution:
        try:
            deconvolver = oms.MetaboliteFeatureDeconvolution()
            params = deconvolver.getDefaults()
            adducts = args.adducts or (
                DEFAULT_NEG_ADDUCTS
                if args.ion_mode == "NEG"
                else DEFAULT_POS_ADDUCTS
            )
            params.setValue(
                "potential_adducts",
                [item.strip().encode() for item in adducts.split(",") if item.strip()],
            )
            # FeatureFindingMetabo reports charge magnitudes even for negative
            # acquisitions, but MetaboliteFeatureDeconvolution requires signed
            # search bounds.  Positive bounds with ``negative_mode=true`` make
            # the implicit [M-H]- count negative and abort the solver.
            if args.ion_mode == "NEG":
                deconv_charge_min = -int(args.charge_high)
                deconv_charge_max = -int(args.charge_low)
            else:
                deconv_charge_min = int(args.charge_low)
                deconv_charge_max = int(args.charge_high)
            params.setValue("charge_min", deconv_charge_min)
            params.setValue("charge_max", deconv_charge_max)
            params.setValue(
                "charge_span_max",
                int(args.charge_high) - int(args.charge_low) + 1,
            )
            params.setValue("q_try", "heuristic")
            params.setValue("mass_max_diff", float(args.adduct_mass_max_diff))
            params.setValue("retention_max_diff", float(args.adduct_rt_max_diff))
            params.setValue(
                "retention_max_diff_local", float(args.adduct_rt_max_diff)
            )
            params.setValue(
                "negative_mode", "true" if args.ion_mode == "NEG" else "false"
            )
            deconvolver.setParameters(params)
            candidate = oms.FeatureMap()
            deconvolver.compute(feature_map, candidate, groups, edges)
            candidate.setUniqueIds()
            decharged = candidate
            deconvolution_status = "ok"
        except Exception as exc:
            deconvolution_status = "failed_preserved_feature_map"
            deconvolution_error = str(exc)

    output = Path(args.output_features)
    output.parent.mkdir(parents=True, exist_ok=True)
    oms.FeatureXMLFile().store(str(output), decharged)
    if args.output_adduct_groups and groups.size():
        oms.ConsensusXMLFile().store(str(args.output_adduct_groups), groups)

    features = _serialize_features(
        decharged, groups if deconvolution_status == "ok" else None
    )
    multi_member_group_ids = {
        str(feature.get("adduct_group_id"))
        for feature in features
        if feature.get("adduct_group_id")
    }
    grouped_features = sum(
        bool(feature.get("adduct_group_id")) for feature in features
    )
    assigned_adducts = sum(bool(feature.get("adduct")) for feature in features)
    unassigned_adducts = len(features) - assigned_adducts
    return {
        "status": "ok" if ms1_count else "unavailable_no_ms1_scans",
        "backend": f"PyOpenMS {oms.__version__}",
        "pipeline": [
            "per-spectrum PeakPickerHiRes/PeakPickerIterative for profile MS1 only",
            "MassTraceDetection",
            "ElutionPeakDetection",
            "FeatureFindingMetabo",
            *(
                ["MetaboliteFeatureDeconvolution"]
                if args.group_adducts
                else []
            ),
        ],
        "input": str(Path(args.input).resolve()),
        "output_features": str(output.resolve()),
        "representation_audit": representation,
        "auxiliary_spectra_excluded": auxiliary_spectra,
        "ms1_spectra": ms1_count,
        "mass_traces": len(mass_traces),
        "split_mass_traces": len(split_traces),
        "raw_features": raw_feature_count,
        "features": features,
        "feature_count": len(features),
        "adduct_group_count": len(multi_member_group_ids),
        "adduct_grouped_feature_count": grouped_features,
        "molecular_entity_group_count_including_singletons": int(groups.size()),
        "adduct_edge_count": int(edges.size()),
        "feature_finder_charge_counts": feature_finder_charge_counts,
        "ion_identity_unassigned_feature_count": unassigned_adducts,
        "isotope_grouped_feature_count": sum(
            int(feature.get("isotope_trace_count") or 1) > 1
            for feature in features
        ),
        "adduct_assignment_count": assigned_adducts,
        "adduct_assignment_coverage_pct": (
            round(100.0 * assigned_adducts / len(features), 3)
            if features
            else 0.0
        ),
        "adduct_deconvolution_status": deconvolution_status,
        "adduct_deconvolution_error": deconvolution_error,
        "adduct_grouping_interpretation": (
            "inferred ion-identity hypothesis; not a compound identification"
        ),
        "parameters": {
            "peak_picker": args.peak_picker,
            "peak_picker_signal_to_noise": args.signal_to_noise,
            "mass_error_ppm": args.mass_error_ppm,
            "noise_threshold_int": args.noise_threshold,
            "chrom_peak_snr": args.chrom_peak_snr,
            "min_trace_length_seconds": args.min_trace_length,
            "width_filtering": args.width_filtering,
            "min_fwhm_seconds": args.min_fwhm,
            "max_fwhm_seconds": args.max_fwhm,
            "isotope_model": args.isotope_model,
            "remove_single_traces": args.remove_single_traces,
            "charge_low": args.charge_low,
            "charge_high": args.charge_high,
            "deconvolution_signed_charge_bounds": (
                [-int(args.charge_high), -int(args.charge_low)]
                if args.ion_mode == "NEG"
                else [int(args.charge_low), int(args.charge_high)]
            ),
            "deconvolution_charge_span_max": (
                int(args.charge_high) - int(args.charge_low) + 1
            ),
            "ion_mode": args.ion_mode,
            "group_adducts": args.group_adducts,
            "max_adduct_grouping_features": int(
                args.max_adduct_grouping_features
            ),
            "effective_max_adduct_grouping_features": (
                effective_deconvolution_limit
            ),
            "adduct_mass_max_diff_da": float(args.adduct_mass_max_diff),
            "adduct_rt_max_diff_seconds": float(args.adduct_rt_max_diff),
        },
    }


def command_align(args: argparse.Namespace) -> dict[str, Any]:
    feature_maps: list[oms.FeatureMap] = []
    labels = []
    for path_text in args.inputs:
        path = Path(path_text)
        feature_map = oms.FeatureMap()
        oms.FeatureXMLFile().load(str(path), feature_map)
        feature_map.setUniqueIds()
        feature_maps.append(feature_map)
        labels.append(path.stem)

    if len(feature_maps) < 2:
        raise ValueError("RT alignment requires at least two featureXML inputs")

    reference_index = max(
        range(len(feature_maps)), key=lambda index: feature_maps[index].size()
    )
    aligner = oms.MapAlignmentAlgorithmPoseClustering()
    aligner.setReference(feature_maps[reference_index])
    transformer = oms.MapAlignmentTransformer()
    alignment_status = []
    for index, feature_map in enumerate(feature_maps):
        if index == reference_index:
            alignment_status.append({"index": index, "status": "reference"})
            continue
        transformation = oms.TransformationDescription()
        try:
            aligner.align(feature_map, transformation)
            transformer.transformRetentionTimes(feature_map, transformation, True)
            alignment_status.append({"index": index, "status": "aligned"})
        except Exception as exc:
            alignment_status.append(
                {"index": index, "status": "failed_unaligned", "error": str(exc)}
            )

    grouper = oms.FeatureGroupingAlgorithmQT()
    params = grouper.getParameters()
    params.setValue("distance_RT:max_difference", float(args.rt_tolerance_seconds))
    params.setValue("distance_MZ:max_difference", float(args.mz_tolerance))
    params.setValue("distance_MZ:unit", args.mz_unit)
    grouper.setParameters(params)

    consensus = oms.ConsensusMap()
    headers = consensus.getColumnHeaders()
    for index, feature_map in enumerate(feature_maps):
        header = oms.ColumnHeader()
        header.filename = labels[index]
        header.size = feature_map.size()
        header.unique_id = feature_map.getUniqueId()
        headers[index] = header
    consensus.setColumnHeaders(headers)
    grouper.group(feature_maps, consensus)
    consensus.setUniqueIds()

    output = Path(args.output_consensus)
    output.parent.mkdir(parents=True, exist_ok=True)
    oms.ConsensusXMLFile().store(str(output), consensus)
    rows = []
    for index, consensus_feature in enumerate(consensus):
        members = [
            {
                "sample_index": int(handle.getMapIndex()),
                "sample": labels[int(handle.getMapIndex())],
                "feature_unique_id": int(handle.getUniqueId()),
                "mz": _finite(handle.getMZ()),
                "rt_seconds": _finite(handle.getRT()),
                "intensity": _finite(handle.getIntensity()),
                "charge": int(handle.getCharge()),
            }
            for handle in consensus_feature.getFeatureList()
        ]
        rows.append(
            {
                "consensus_id": f"CF{index + 1:07d}",
                "mz": _finite(consensus_feature.getMZ()),
                "rt_seconds": _finite(consensus_feature.getRT()),
                "intensity": _finite(consensus_feature.getIntensity()),
                "charge": int(consensus_feature.getCharge()),
                "sample_count": len(members),
                "members": members,
            }
        )
    sample_count_distribution = dict(
        sorted(Counter(row["sample_count"] for row in rows).items())
    )
    return {
        "status": "ok",
        "backend": f"PyOpenMS {oms.__version__}",
        "pipeline": [
            "MapAlignmentAlgorithmPoseClustering",
            "FeatureGroupingAlgorithmQT",
        ],
        "inputs": [str(Path(path).resolve()) for path in args.inputs],
        "reference_index": reference_index,
        "reference": labels[reference_index],
        "alignment_status": alignment_status,
        "rt_tolerance_seconds": float(args.rt_tolerance_seconds),
        "mz_tolerance": float(args.mz_tolerance),
        "mz_unit": args.mz_unit,
        "consensus_count": len(rows),
        "multi_sample_consensus_count": sum(
            row["sample_count"] > 1 for row in rows
        ),
        "single_sample_consensus_count": sum(
            row["sample_count"] == 1 for row in rows
        ),
        "consensus_sample_count_distribution": sample_count_distribution,
        "consensus_features": rows,
        "output_consensus": str(output.resolve()),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-output", required=True)
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("probe")

    peak = subparsers.add_parser("peak-pick")
    peak.add_argument("--input", required=True)
    peak.add_argument("--output", required=True)
    peak.add_argument("--strategy", choices=["hires", "iterative"], default="hires")
    peak.add_argument("--signal-to-noise", type=float, default=1.0)

    features = subparsers.add_parser("features")
    features.add_argument("--input", required=True)
    features.add_argument("--output-features", required=True)
    features.add_argument("--output-adduct-groups")
    features.add_argument("--peak-picker", choices=["hires", "iterative"], default="hires")
    features.add_argument("--signal-to-noise", type=float, default=1.0)
    features.add_argument("--mass-error-ppm", type=float, default=10.0)
    features.add_argument("--noise-threshold", type=float, default=1000.0)
    features.add_argument("--chrom-peak-snr", type=float, default=3.0)
    features.add_argument("--min-trace-length", type=float, default=5.0)
    features.add_argument(
        "--width-filtering", choices=["fixed", "auto"], default="fixed"
    )
    features.add_argument("--min-fwhm", type=float, default=1.0)
    features.add_argument("--max-fwhm", type=float, default=60.0)
    features.add_argument(
        "--isotope-model", default="metabolites (5% RMS)"
    )
    features.add_argument("--remove-single-traces", action="store_true")
    features.add_argument("--charge-low", type=int, default=1)
    features.add_argument("--charge-high", type=int, default=1)
    features.add_argument("--ion-mode", choices=["POS", "NEG"], required=True)
    features.add_argument("--group-adducts", action="store_true")
    features.add_argument(
        "--max-adduct-grouping-features",
        type=int,
        default=DEFAULT_MAX_ADDUCT_GROUPING_FEATURES,
        help=(
            "Skip only the optional adduct-grouping ILP above this raw "
            "FeatureFindingMetabo feature count. Values above the production "
            f"safety ceiling ({MAX_SAFE_ADDUCT_GROUPING_FEATURES}) are clamped; "
            "0 uses that ceiling."
        ),
    )
    features.add_argument("--adducts")
    features.add_argument("--adduct-mass-max-diff", type=float, default=0.05)
    features.add_argument("--adduct-rt-max-diff", type=float, default=10.0)

    align = subparsers.add_parser("align")
    align.add_argument("--inputs", nargs="+", required=True)
    align.add_argument("--output-consensus", required=True)
    align.add_argument("--rt-tolerance-seconds", type=float, default=20.0)
    align.add_argument("--mz-tolerance", type=float, default=10.0)
    align.add_argument("--mz-unit", choices=["ppm", "Da"], default="ppm")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    try:
        if args.command == "probe":
            payload = command_probe(args)
        elif args.command == "peak-pick":
            payload = command_peak_pick(args)
        elif args.command == "features":
            payload = command_features(args)
        else:
            payload = command_align(args)
    except Exception as exc:
        payload = {
            "status": "error",
            "backend": f"PyOpenMS {getattr(oms, '__version__', 'unknown')}",
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        _write_json(args.json_output, payload)
        raise
    _write_json(args.json_output, payload)


if __name__ == "__main__":
    main()
