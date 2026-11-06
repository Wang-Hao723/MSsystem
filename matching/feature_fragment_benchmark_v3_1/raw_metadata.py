"""Read-only mzML and mzXML metadata extraction for the V3.1 amendment.

mzML values retain PSI-MS CV accessions and use ``mzml_cv_param`` as their
origin.  mzXML values come from literal XML attributes and use
``mzxml_attribute``.  The interfaces never convert collision-energy units and
never fill an absent field from an instrument-family assumption.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET


DEFAULT_MAPPING_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/feature_fragment_benchmark_v3_1/psi_ms_mapping.json"
)


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def load_psi_ms_mapping(path: str | Path = DEFAULT_MAPPING_PATH) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def validate_psi_ms_mapping(mapping: Mapping[str, Any]) -> None:
    mzml = mapping.get("mzml") or {}
    activation = mzml.get("activation") or {}
    excluded = mzml.get("terms_not_activation") or {}
    expected = {
        "MS:1000422",
        "MS:1002481",
        "MS:1002472",
        "MS:1000133",
    }
    if not expected.issubset(activation):
        raise ValueError("psi_ms_activation_mapping_incomplete")
    if "MS:1002482" in activation or "MS:1002482" not in excluded:
        raise ValueError("psi_ms_1002482_must_be_non_activation")
    if (mzml.get("collision_energy") or {}).get("absolute") != "MS:1000045":
        raise ValueError("psi_ms_collision_energy_mapping_invalid")
    if (mzml.get("scan_window") or {}).get("lower_limit") != "MS:1000501":
        raise ValueError("psi_ms_scan_lower_mapping_invalid")
    if (mzml.get("scan_window") or {}).get("upper_limit") != "MS:1000500":
        raise ValueError("psi_ms_scan_upper_mapping_invalid")


def _cv_params(element: ET.Element) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for node in element.iter():
        if _local_name(node.tag) != "cvParam":
            continue
        rows.append(
            {
                "accession": node.attrib.get("accession"),
                "name": node.attrib.get("name"),
                "value": node.attrib.get("value"),
                "unit_accession": node.attrib.get("unitAccession"),
                "unit_name": node.attrib.get("unitName"),
                "value_origin": "mzml_cv_param",
            }
        )
    return rows


def _first_cv_value(
    params: Sequence[Mapping[str, Any]], accession: str
) -> dict[str, Any] | None:
    return next((dict(row) for row in params if row.get("accession") == accession), None)


def _all_cv_values(
    params: Sequence[Mapping[str, Any]], accession: str
) -> list[dict[str, Any]]:
    return [dict(row) for row in params if row.get("accession") == accession]


def _native_unit(param: Mapping[str, Any], *, normalized: bool = False) -> str | None:
    if normalized:
        return "NCE"
    unit_accession = str(param.get("unit_accession") or "")
    unit_name = str(param.get("unit_name") or "").strip().lower()
    if unit_accession == "UO:0000266" or unit_name in {
        "electronvolt",
        "electron volt",
        "electronvolts",
        "ev",
    }:
        return "eV"
    if unit_name in {"volt", "volts", "v"}:
        return "V"
    if unit_name in {"percent", "%"}:
        return "percent_nominal"
    return None


def _mzml_collision_energy(
    activation_params: Sequence[Mapping[str, Any]], mapping: Mapping[str, Any]
) -> dict[str, Any]:
    terms = mapping["mzml"]["collision_energy"]
    absolute = _all_cv_values(activation_params, terms["absolute"])
    normalized = _all_cv_values(activation_params, terms["normalized"])
    ramp_ev_start = _first_cv_value(activation_params, terms["ramp_start_ev"])
    ramp_ev_end = _first_cv_value(activation_params, terms["ramp_end_ev"])
    ramp_pct_start = _first_cv_value(activation_params, terms["ramp_start_percent"])
    ramp_pct_end = _first_cv_value(activation_params, terms["ramp_end_percent"])
    base = {
        "value_origin": "mzml_cv_param",
        "conversion_applied": False,
        "cross_unit_equivalence_assumed": False,
        "cross_vendor_equivalence_assumed": False,
        "equation": None,
        "weights": [],
    }
    if ramp_ev_start and ramp_ev_end:
        values = [_float(ramp_ev_start.get("value")), _float(ramp_ev_end.get("value"))]
        if all(value is not None for value in values):
            return {
                **base,
                "status": "resolved",
                "mode": "ramp",
                "unit": "eV",
                "values": values,
                "weights": [None, None],
                "minimum": min(values),
                "maximum": max(values),
                "source_cv_params": [ramp_ev_start, ramp_ev_end],
            }
    if ramp_pct_start and ramp_pct_end:
        values = [_float(ramp_pct_start.get("value")), _float(ramp_pct_end.get("value"))]
        if all(value is not None for value in values):
            return {
                **base,
                "status": "resolved",
                "mode": "ramp",
                "unit": "percent_nominal",
                "values": values,
                "weights": [None, None],
                "minimum": min(values),
                "maximum": max(values),
                "source_cv_params": [ramp_pct_start, ramp_pct_end],
            }
    candidates = [(row, False) for row in absolute] + [
        (row, True) for row in normalized
    ]
    numeric = [
        (row, is_normalized, _float(row.get("value")))
        for row, is_normalized in candidates
    ]
    numeric = [row for row in numeric if row[2] is not None]
    if numeric:
        units = {
            _native_unit(row, normalized=is_normalized)
            for row, is_normalized, _ in numeric
        }
        if len(units) != 1 or None in units:
            return {
                **base,
                "status": "unresolved",
                "mode": "unresolved",
                "unit": None,
                "values": [value for _, _, value in numeric],
                "reason": "collision_energy_unit_missing_or_mixed",
                "source_cv_params": [row for row, _, _ in numeric],
            }
        values = [value for _, _, value in numeric]
        return {
            **base,
            "status": "resolved",
            "mode": "fixed" if len(values) == 1 else "stepped",
            "unit": next(iter(units)),
            "values": values,
            "weights": [1.0] if len(values) == 1 else [None] * len(values),
            "center": values[0] if len(values) == 1 else None,
            "minimum": min(values),
            "maximum": max(values),
            "source_cv_params": [row for row, _, _ in numeric],
        }
    return {
        **base,
        "status": "unresolved",
        "mode": "unresolved",
        "unit": None,
        "values": [],
        "reason": "collision_energy_cv_param_missing",
        "source_cv_params": [],
    }


def _mzml_isolation(
    spectrum: ET.Element, mapping: Mapping[str, Any]
) -> dict[str, Any]:
    terms = mapping["mzml"]["isolation_window"]
    windows = [
        element
        for element in spectrum.iter()
        if _local_name(element.tag) == "isolationWindow"
    ]
    if not windows:
        return {
            "status": "unreported",
            "target_mz": None,
            "lower_offset_mz": None,
            "upper_offset_mz": None,
            "lower_mz": None,
            "upper_mz": None,
            "value_origin": "mzml_cv_param",
            "source_cv_params": [],
        }
    params = _cv_params(windows[0])
    target_row = _first_cv_value(params, terms["target_mz"])
    lower_row = _first_cv_value(params, terms["lower_offset"])
    upper_row = _first_cv_value(params, terms["upper_offset"])
    target = _float(target_row.get("value")) if target_row else None
    lower_offset = _float(lower_row.get("value")) if lower_row else None
    upper_offset = _float(upper_row.get("value")) if upper_row else None
    verified = bool(
        target is not None
        and lower_offset is not None
        and upper_offset is not None
        and lower_offset >= 0
        and upper_offset >= 0
    )
    return {
        "status": "verified" if verified else "incomplete",
        "target_mz": target,
        "lower_offset_mz": lower_offset,
        "upper_offset_mz": upper_offset,
        "lower_mz": target - lower_offset if verified else None,
        "upper_mz": target + upper_offset if verified else None,
        "width_mz": lower_offset + upper_offset if verified else None,
        "value_origin": "mzml_cv_param",
        "source_cv_params": [
            row for row in (target_row, lower_row, upper_row) if row is not None
        ],
    }


def _mzml_scan_range(
    spectrum: ET.Element, mapping: Mapping[str, Any]
) -> dict[str, Any]:
    terms = mapping["mzml"]["scan_window"]
    windows = [
        element
        for element in spectrum.iter()
        if _local_name(element.tag) == "scanWindow"
    ]
    if not windows:
        return {
            "status": "unresolved",
            "minimum_mz": None,
            "maximum_mz": None,
            "value_origin": "mzml_cv_param",
            "source_cv_params": [],
        }
    params = _cv_params(windows[0])
    lower_row = _first_cv_value(params, terms["lower_limit"])
    upper_row = _first_cv_value(params, terms["upper_limit"])
    lower = _float(lower_row.get("value")) if lower_row else None
    upper = _float(upper_row.get("value")) if upper_row else None
    valid = bool(lower is not None and upper is not None and 0 <= lower < upper)
    return {
        "status": "resolved" if valid else "unresolved",
        "minimum_mz": lower if valid else None,
        "maximum_mz": upper if valid else None,
        "value_origin": "mzml_cv_param",
        "source_cv_params": [
            row for row in (lower_row, upper_row) if row is not None
        ],
    }


def _mzml_selected_ion(spectrum: ET.Element) -> dict[str, Any]:
    params = _cv_params(spectrum)
    precursor = _first_cv_value(params, "MS:1000744")
    charge = _first_cv_value(params, "MS:1000041")
    return {
        "precursor_mz": _float(precursor.get("value")) if precursor else None,
        "charge": _int(charge.get("value")) if charge else None,
        "value_origin": "mzml_cv_param",
        "source_cv_params": [row for row in (precursor, charge) if row is not None],
    }


def _mzml_activation(
    spectrum: ET.Element, mapping: Mapping[str, Any]
) -> dict[str, Any]:
    elements = [
        element
        for element in spectrum.iter()
        if _local_name(element.tag) == "activation"
    ]
    params = _cv_params(elements[0]) if elements else []
    activation_map = mapping["mzml"]["activation"]
    methods = [
        {
            "accession": row["accession"],
            **activation_map[row["accession"]],
            "source_cv_param": row,
        }
        for row in params
        if row.get("accession") in activation_map
    ]
    forbidden = [
        row
        for row in params
        if row.get("accession")
        in mapping["mzml"].get("terms_not_activation", {})
    ]
    return {
        "status": "resolved" if methods else "unresolved",
        "methods": methods,
        "value_origin": "mzml_cv_param",
        "non_activation_terms_seen": forbidden,
        "collision_energy": _mzml_collision_energy(params, mapping),
    }


def extract_mzml_spectrum_metadata(
    path: str | Path,
    *,
    spectrum_id: str | None = None,
    index: int | None = None,
    mapping_path: str | Path = DEFAULT_MAPPING_PATH,
) -> dict[str, Any]:
    """Extract one mzML spectrum's literal acquisition metadata.

    Exactly one of ``spectrum_id`` or ``index`` must identify the spectrum.
    The returned file hash covers the raw mzML bytes.
    """
    if (spectrum_id is None) == (index is None):
        raise ValueError("provide_exactly_one_of_spectrum_id_or_index")
    path = Path(path)
    mapping = load_psi_ms_mapping(mapping_path)
    validate_psi_ms_mapping(mapping)
    instrument_configurations: dict[str, list[dict[str, Any]]] = {}
    selected: ET.Element | None = None
    selected_instrument_ref: str | None = None
    for _, element in ET.iterparse(path, events=("end",)):
        name = _local_name(element.tag)
        if name == "instrumentConfiguration":
            config_id = element.attrib.get("id")
            if config_id:
                instrument_configurations[config_id] = _cv_params(element)
            element.clear()
            continue
        if name != "spectrum":
            continue
        matches = (
            spectrum_id is not None
            and element.attrib.get("id") == spectrum_id
        ) or (
            index is not None and _int(element.attrib.get("index")) == index
        )
        if matches:
            selected = element
            selected_instrument_ref = next(
                (
                    node.attrib.get("instrumentConfigurationRef")
                    for node in element.iter()
                    if _local_name(node.tag) == "scan"
                    and node.attrib.get("instrumentConfigurationRef")
                ),
                None,
            )
            break
        element.clear()
    if selected is None:
        raise KeyError("mzml_spectrum_not_found")
    params = _cv_params(selected)
    polarity = (
        "positive"
        if _first_cv_value(params, "MS:1000130")
        else "negative"
        if _first_cv_value(params, "MS:1000129")
        else None
    )
    ms_level_row = _first_cv_value(params, "MS:1000511")
    result = {
        "schema_version": "3.1.0",
        "file_format": "mzML",
        "source_file": str(path.resolve()),
        "source_file_sha256": _sha256(path),
        "spectrum_locator": {
            "id": selected.attrib.get("id"),
            "index": _int(selected.attrib.get("index")),
        },
        "value_origin": "mzml_cv_param",
        "ms_level": _int(ms_level_row.get("value")) if ms_level_row else None,
        "polarity": polarity,
        "selected_ion": _mzml_selected_ion(selected),
        "precursor_isolation": _mzml_isolation(selected, mapping),
        "activation": _mzml_activation(selected, mapping),
        "scan_range": _mzml_scan_range(selected, mapping),
        "instrument_configuration_ref": selected_instrument_ref,
        "instrument_configuration_cv_params": instrument_configurations.get(
            selected_instrument_ref or "", []
        ),
        "conversion_applied": False,
    }
    result["blocking_missing_fields"] = sorted(
        name
        for name, missing in {
            "instrument_configuration": not result["instrument_configuration_cv_params"],
            "activation": result["activation"]["status"] != "resolved",
            "collision_energy": result["activation"]["collision_energy"]["status"]
            != "resolved",
            "scan_range": result["scan_range"]["status"] != "resolved",
            "precursor_isolation": result["precursor_isolation"]["status"]
            != "verified",
            "polarity": polarity is None,
            "precursor_mz": result["selected_ion"]["precursor_mz"] is None,
            "charge": result["selected_ion"]["charge"] is None,
        }.items()
        if missing
    )
    return result


def _attribute_value(
    attributes: Mapping[str, str], *names: str
) -> tuple[str | None, str | None]:
    for name in names:
        if name in attributes:
            return attributes[name], name
    return None, None


def _reported_attributes(
    *pairs: tuple[str | None, Any]
) -> dict[str, Any]:
    return {name: value for name, value in pairs if name is not None}


def extract_mzxml_scan_metadata(
    path: str | Path, *, scan_number: int | str
) -> dict[str, Any]:
    """Extract one mzXML scan without pretending attributes are mzML CV terms."""
    path = Path(path)
    target = str(scan_number)
    selected: ET.Element | None = None
    for _, element in ET.iterparse(path, events=("end",)):
        if _local_name(element.tag) != "scan":
            continue
        if str(element.attrib.get("num")) == target:
            selected = element
            break
        element.clear()
    if selected is None:
        raise KeyError("mzxml_scan_not_found")
    precursor_node = next(
        (
            node
            for node in selected.iter()
            if _local_name(node.tag) == "precursorMz"
        ),
        None,
    )
    precursor_attributes = precursor_node.attrib if precursor_node is not None else {}
    precursor_mz = _float(precursor_node.text) if precursor_node is not None else None
    width_raw, width_name = _attribute_value(
        precursor_attributes,
        "windowWideness",
        "isolationWindowWidth",
    )
    width = _float(width_raw)
    lower_raw, lower_name = _attribute_value(
        precursor_attributes, "isolationWindowLowerOffset", "lowerOffset"
    )
    upper_raw, upper_name = _attribute_value(
        precursor_attributes, "isolationWindowUpperOffset", "upperOffset"
    )
    lower_offset = _float(lower_raw)
    upper_offset = _float(upper_raw)
    if width is not None and width >= 0 and lower_offset is None and upper_offset is None:
        lower_offset = width / 2.0
        upper_offset = width / 2.0
    isolation_verified = bool(
        precursor_mz is not None
        and lower_offset is not None
        and upper_offset is not None
        and lower_offset >= 0
        and upper_offset >= 0
    )
    collision_raw, collision_attribute = _attribute_value(
        precursor_attributes, "collisionEnergy"
    )
    if collision_raw is None:
        collision_raw, collision_attribute = _attribute_value(
            selected.attrib, "collisionEnergy"
        )
    unit_raw, unit_attribute = _attribute_value(
        precursor_attributes,
        "collisionEnergyUnit",
        "collisionEnergyUnits",
    )
    if unit_raw is None:
        unit_raw, unit_attribute = _attribute_value(
            selected.attrib, "collisionEnergyUnit", "collisionEnergyUnits"
        )
    collision_value = _float(collision_raw)
    collision_status = (
        "resolved"
        if collision_value is not None and unit_raw
        else "unit_unresolved"
        if collision_value is not None
        else "unresolved"
    )
    activation_raw, activation_attribute = _attribute_value(
        precursor_attributes, "activationMethod"
    )
    if activation_raw is None:
        activation_raw, activation_attribute = _attribute_value(
            selected.attrib, "activationMethod"
        )
    lower_mz = _float(selected.attrib.get("lowMz"))
    upper_mz = _float(selected.attrib.get("highMz"))
    scan_valid = bool(
        lower_mz is not None and upper_mz is not None and 0 <= lower_mz < upper_mz
    )
    polarity_raw = selected.attrib.get("polarity")
    polarity = (
        "positive" if polarity_raw == "+" else "negative" if polarity_raw == "-" else None
    )
    result = {
        "schema_version": "3.1.0",
        "file_format": "mzXML",
        "source_file": str(path.resolve()),
        "source_file_sha256": _sha256(path),
        "spectrum_locator": {"scan_number": target},
        "value_origin": "mzxml_attribute",
        "ms_level": _int(selected.attrib.get("msLevel")),
        "polarity": polarity,
        "selected_ion": {
            "precursor_mz": precursor_mz,
            "charge": _int(precursor_attributes.get("precursorCharge")),
            "value_origin": "mzxml_attribute",
            "source_attributes": {
                "precursorMz_text": precursor_node.text if precursor_node is not None else None,
                "precursorCharge": precursor_attributes.get("precursorCharge"),
            },
        },
        "precursor_isolation": {
            "status": "verified" if isolation_verified else "unreported",
            "target_mz": precursor_mz,
            "lower_offset_mz": lower_offset,
            "upper_offset_mz": upper_offset,
            "lower_mz": precursor_mz - lower_offset if isolation_verified else None,
            "upper_mz": precursor_mz + upper_offset if isolation_verified else None,
            "width_mz": lower_offset + upper_offset if isolation_verified else None,
            "value_origin": "mzxml_attribute",
            "source_attributes": _reported_attributes(
                (width_name, width_raw),
                (lower_name, lower_raw),
                (upper_name, upper_raw),
            ),
        },
        "activation": {
            "status": "resolved" if activation_raw else "unresolved",
            "literal": activation_raw,
            "value_origin": "mzxml_attribute",
            "source_attribute": activation_attribute,
            "collision_energy": {
                "status": collision_status,
                "mode": "fixed" if collision_value is not None else "unresolved",
                "unit": unit_raw,
                "values": [collision_value] if collision_value is not None else [],
                "weights": [1.0] if collision_value is not None else [],
                "value_origin": "mzxml_attribute",
                "source_attributes": _reported_attributes(
                    (collision_attribute, collision_raw),
                    (unit_attribute, unit_raw),
                ),
                "conversion_applied": False,
            },
        },
        "scan_range": {
            "status": "resolved" if scan_valid else "unresolved",
            "minimum_mz": lower_mz if scan_valid else None,
            "maximum_mz": upper_mz if scan_valid else None,
            "value_origin": "mzxml_attribute",
            "source_attributes": {
                "lowMz": selected.attrib.get("lowMz"),
                "highMz": selected.attrib.get("highMz"),
            },
        },
        "conversion_applied": False,
    }
    result["blocking_missing_fields"] = sorted(
        name
        for name, missing in {
            "activation": result["activation"]["status"] != "resolved",
            "collision_energy_unit": collision_status != "resolved",
            "scan_range": result["scan_range"]["status"] != "resolved",
            "precursor_isolation": result["precursor_isolation"]["status"]
            != "verified",
            "polarity": polarity is None,
            "precursor_mz": precursor_mz is None,
            "charge": result["selected_ion"]["charge"] is None,
        }.items()
        if missing
    )
    return result
