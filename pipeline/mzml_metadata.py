"""Preserve query acquisition quantities without converting instrument domains.

The accepted reference-library normalizer is deliberately unchanged. mzML CV
quantity names and their units, rather than a string search over the entire
activation block, determine the collision-energy domain.
"""
from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from pathlib import Path

from utils.collision_energy import canonical_collision_energy_unit
from utils.instrument_metadata import normalize_instrument_family


def _local(tag):
    return str(tag).rsplit("}", 1)[-1]


def _quantity(value):
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) and number >= 0 else None


def query_acquisition_projection(record):
    """Keep the acquisition conditions through analysis and result serialization."""
    keys = ("collision_energy_mode", "activation_method", "activation_raw",
            "instrument_domain_id", "collision_cell_region",
            "collision_energy_quantities", "instrument_configuration_id",
            "acquisition_metadata_source", "ion_injection_time",
            "mass_resolving_power", "preset_scan_configuration", "scan_filter_string",
            "reported_precursor_charge", "precursor_charge_status",
            "adduct_charge_evidence")
    return {key: record[key] for key in keys if key in record}


def read_instrument_header(path):
    """Resolve parameter-group and per-scan configuration references in the header.

    Stop before the first spectrum; binary arrays are neither scanned nor read.
    Model terms are direct configuration CV terms (including group references),
    never component detector names or an inferred model from a file name.
    """
    groups, configs, default = {}, {}, ""
    with Path(path).open("rb") as handle:
        for event, elem in ET.iterparse(handle, events=("start", "end")):
            name = _local(elem.tag)
            if event == "start":
                if name == "run":
                    default = elem.get("defaultInstrumentConfigurationRef", "")
                if name in {"spectrumList", "chromatogramList"}:
                    break
                continue
            if name == "referenceableParamGroup":
                groups[elem.get("id", "")] = [dict(c.attrib) for c in elem if _local(c.tag) == "cvParam"]
            elif name == "instrumentConfiguration":
                params = [dict(c.attrib) for c in elem if _local(c.tag) == "cvParam"]
                for child in elem:
                    if _local(child.tag) == "referenceableParamGroupRef":
                        params.extend(groups.get(child.get("ref", ""), []))
                model_names = [p.get("name", "") for p in params
                               if p.get("name", "").lower() not in {
                                   "instrument model", "instrument serial number", "instrument attribute"
                               } and normalize_instrument_family(p.get("name", ""))]
                component_names = [c.get("name", "") for c in elem.iter() if _local(c.tag) == "cvParam"]
                model = model_names[0] if len(set(model_names)) == 1 else ""
                configs[elem.get("id", "")] = {
                    "instrument": model,
                    "instrument_type": normalize_instrument_family(model or " ".join(component_names)),
                }
    if not default and len(configs) == 1:
        default = next(iter(configs))
    return {"default": default, "configurations": configs}


def acquisition_metadata(raw, header):
    """Read only explicit scan/activation metadata; retain separate CE quantities."""
    scans = (raw.get("scanList") or {}).get("scan") or [{}]
    scan = scans[0]
    config_id = scan.get("instrumentConfigurationRef") or header["default"]
    instrument = header["configurations"].get(config_id, {})
    precursors = (raw.get("precursorList") or {}).get("precursor") or [{}]
    if len(precursors) > 1 or len((precursors[0].get("selectedIonList") or {}).get("selectedIon") or []) > 1:
        raise ValueError("mzML contains a multi-precursor MS2 spectrum; single-precursor library matching requires deconvolution or a precursor-specific export")
    activation = precursors[0].get("activation") or {}
    terms = []
    for key, forced_unit in (("collision energy", ""), ("normalized collision energy", "NCE")):
        value = activation.get(key)
        if value is None:
            continue
        values = value if isinstance(value, (list, tuple)) else [value]
        for item in values:
            terms.append({"term": key, "value": _quantity(item),
                          "unit": forced_unit or canonical_collision_energy_unit(getattr(item, "unit_info", "")),
                          "raw": str(item)})
    valid = [t for t in terms if t["value"] is not None]
    # Actual eV/voltage and normalized settings are distinct quantities. Prefer
    # the explicit absolute quantity when supplied; never attach NCE to it.
    chosen = next((t for t in valid if t["term"] == "collision energy" and t["unit"] == "eV"), None)
    chosen = chosen or next((t for t in valid if t["term"] == "normalized collision energy"), None)
    chosen = chosen or next(iter(valid), None)
    activation_names = " ".join(str(k).lower() for k in activation)
    repeated = any(sum(t["term"] == term for t in valid) > 1 for term in {t["term"] for t in valid})
    mode = "stepped" if repeated or "step" in activation_names else "ramp" if "ramp" in activation_names else "fixed" if chosen else "missing"
    ce = chosen["value"] if chosen and mode == "fixed" else None
    result = {
        **instrument,
        "instrument_configuration_id": config_id,
        "collision_energy": ce,
        "collision_energy_unit": chosen["unit"] if chosen else "",
        "collision_energy_mode": mode,
        "collision_energy_quantities": terms,
        "activation_method": " ".join(str(k) for k in activation
                                       if "dissociation" in str(k).lower() or str(k).lower() in {"hcd", "cid"}),
        "acquisition_metadata_source": "mzml_explicit_scan_cv_v1",
    }
    for key, output in (("ion injection time", "ion_injection_time"),
                        ("mass resolving power", "mass_resolving_power"),
                        ("preset scan configuration", "preset_scan_configuration")):
        value = scan.get(key, raw.get(key))
        if value is not None:
            result[output] = {"value": _quantity(value), "unit": str(getattr(value, "unit_info", "") or ""), "raw": str(value)}
    if "filter string" in scan:
        result["scan_filter_string"] = str(scan["filter string"])
    if "positive scan" in raw and "negative scan" not in raw:
        result["scan_polarity"] = "POS"
    elif "negative scan" in raw and "positive scan" not in raw:
        result["scan_polarity"] = "NEG"
    return result
