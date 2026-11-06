"""Strict parsers for scan ranges and collision-energy programs."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from typing import Any

from matching.feature_fragment_benchmark.corpus import parse_scan_range as _parse_scan_range


_NUMBER = r"(?:\d+(?:\.\d*)?|\.\d+)"
_UNIT = r"(?:eV|V|NCE|%\s*\(nominal\)|%|percent(?:_nominal)?)"


def parse_scan_range(raw: Any) -> dict[str, Any]:
    """Return an auditable scan-range object or an explicit unresolved state."""
    parsed = _parse_scan_range(raw)
    if parsed is None:
        return {
            "raw": None if raw is None else str(raw).strip(),
            "minimum_mz": None,
            "maximum_mz": None,
            "status": "unresolved",
            "parser": "explicit_positive_range_separator_v3",
        }
    return {
        "raw": str(raw).strip(),
        "minimum_mz": parsed[0],
        "maximum_mz": parsed[1],
        "status": "resolved",
        "parser": "explicit_positive_range_separator_v3",
    }


def _unit(raw: str) -> str:
    compact = re.sub(r"\s+", "", raw).lower()
    if compact == "ev":
        return "eV"
    if compact == "v":
        return "V"
    if compact == "nce":
        return "NCE"
    if compact in {"%", "%(nominal)", "percent", "percent_nominal"}:
        return "percent_nominal"
    raise ValueError(f"unsupported collision-energy unit: {raw}")


def _base(raw: str) -> dict[str, Any]:
    return {
        "raw": raw,
        "status": "resolved",
        "conversion_applied": False,
        "cross_unit_equivalence_assumed": False,
        "cross_vendor_equivalence_assumed": False,
        "equation": None,
        "values": [],
        "weights": [],
        "center": None,
        "minimum": None,
        "maximum": None,
    }


def parse_collision_energy_program(raw_value: Any) -> dict[str, Any]:
    """Parse fixed, ramped, stepped, CES, dynamic, and DIA CE programs.

    The returned unit is the literal native unit.  No NCE/eV/V conversion is
    performed.  Composite spectra retain their complete value vector rather
    than being reduced to a midpoint.
    """
    raw = re.sub(r"\s+", " ", str(raw_value or "")).strip()
    if not raw:
        return {
            **_base(raw),
            "status": "unresolved",
            "mode": "unresolved",
            "unit": None,
            "reason": "collision_energy_missing",
        }

    lowered = raw.lower()
    acquisition_mode = None
    if re.search(r"\bms\s*e\b|\bmse\b", lowered):
        acquisition_mode = "MSE"
    elif "bbcid" in lowered:
        acquisition_mode = "bbCID"
    elif "all ions" in lowered or "all-ion" in lowered:
        acquisition_mode = "all_ions"

    ces = re.search(
        rf"(?:ce\s*[:=]?\s*)?({_NUMBER})\s*(V|eV)\b.*?\bCES\s*(?:[:=]|±|\+/-)?\s*({_NUMBER})\s*(V|eV)?",
        raw,
        flags=re.IGNORECASE,
    )
    if ces:
        center = float(ces.group(1))
        spread = float(ces.group(3))
        unit = _unit(ces.group(2))
        spread_unit = _unit(ces.group(4) or ces.group(2))
        if unit != spread_unit or spread <= 0 or center - spread < 0:
            return {
                **_base(raw),
                "status": "unresolved",
                "mode": "CES",
                "unit": None,
                "reason": "ces_unit_or_range_invalid",
            }
        values = [center - spread, center, center + spread]
        return {
            **_base(raw),
            "mode": "CES",
            "unit": unit,
            "values": values,
            "weights": [None, None, None],
            "center": center,
            "minimum": values[0],
            "maximum": values[-1],
            "spread": spread,
            "acquisition_mode": acquisition_mode,
        }

    equation = re.search(
        rf"(?:CE\s*=|collision energy\s*=)\s*([^;]+(?:m/z|mz)[^;]*|[^;]+)",
        raw,
        flags=re.IGNORECASE,
    )
    if equation and any(token in lowered for token in ("rolling", "dynamic", "slope", "offset", "m/z", "mz")):
        unit_match = re.search(rf"({_UNIT})\b", raw, flags=re.IGNORECASE)
        return {
            **_base(raw),
            "mode": "rolling" if "rolling" in lowered else "dynamic",
            "unit": _unit(unit_match.group(1)) if unit_match else None,
            "equation": equation.group(1).strip(),
            "acquisition_mode": acquisition_mode,
        }

    ramp = re.search(
        rf"(?:ramp(?:ed)?\s*)?({_NUMBER})\s*(?:-|–|—|\bto\b)\s*({_NUMBER})\s*({_UNIT})",
        raw,
        flags=re.IGNORECASE,
    )
    explicit_ramp = any(token in lowered for token in ("ramp", "ramped", "mse", "bbcid", "all ions", "all-ion"))
    if ramp and explicit_ramp:
        minimum, maximum = sorted((float(ramp.group(1)), float(ramp.group(2))))
        if minimum == maximum:
            explicit_ramp = False
        else:
            return {
                **_base(raw),
                "mode": acquisition_mode or "ramp",
                "unit": _unit(ramp.group(3)),
                "values": [minimum, maximum],
                "weights": [None, None],
                "center": None,
                "minimum": minimum,
                "maximum": maximum,
                "acquisition_mode": acquisition_mode,
                "precursor_selected": acquisition_mode not in {"MSE", "bbCID", "all_ions"},
            }

    stepped_prefix = bool(re.search(r"\b(?:stepped?|steps?)\b", lowered))
    unit_match = re.search(rf"({_UNIT})\s*$", raw, flags=re.IGNORECASE)
    if stepped_prefix and unit_match:
        unit = _unit(unit_match.group(1))
        numeric_part = raw[: unit_match.start()]
        values = [float(value) for value in re.findall(_NUMBER, numeric_part)]
        if len(values) >= 2 and all(math.isfinite(value) and value >= 0 for value in values):
            return {
                **_base(raw),
                "mode": "stepped",
                "unit": unit,
                "values": values,
                "weights": [None] * len(values),
                "center": None,
                "minimum": min(values),
                "maximum": max(values),
                "acquisition_mode": acquisition_mode,
            }

    fixed = re.fullmatch(
        rf"(?:CE\s*[:=]?\s*)?({_NUMBER})\s*({_UNIT})",
        raw,
        flags=re.IGNORECASE,
    )
    if fixed:
        value = float(fixed.group(1))
        if math.isfinite(value) and value >= 0:
            return {
                **_base(raw),
                "mode": "fixed",
                "unit": _unit(fixed.group(2)),
                "values": [value],
                "weights": [1.0],
                "center": value,
                "minimum": value,
                "maximum": value,
                "acquisition_mode": acquisition_mode,
                "precursor_selected": True,
            }

    nominal = re.fullmatch(rf"({_NUMBER})\s*%\s*\(nominal\)", raw, flags=re.I)
    if nominal:
        value = float(nominal.group(1))
        return {
            **_base(raw),
            "mode": "fixed",
            "unit": "percent_nominal",
            "values": [value],
            "weights": [1.0],
            "center": value,
            "minimum": value,
            "maximum": value,
            "acquisition_mode": acquisition_mode,
            "precursor_selected": True,
        }

    return {
        **_base(raw),
        "status": "unresolved",
        "mode": acquisition_mode or "unresolved",
        "unit": None,
        "reason": "collision_energy_program_unresolved",
        "acquisition_mode": acquisition_mode,
    }


def ce_program_key(program: Mapping[str, Any]) -> tuple[Any, ...]:
    """Return a lossless hash key for a resolved CE program."""
    return (
        program.get("unit"),
        program.get("mode"),
        tuple(program.get("values") or ()),
        tuple(program.get("weights") or ()),
        program.get("equation"),
        program.get("acquisition_mode"),
    )


def ce_program_has_complete_vector(program: Mapping[str, Any]) -> bool:
    """Check that composite programs were not reduced to a midpoint."""
    mode = program.get("mode")
    values: Sequence[Any] = program.get("values") or []
    if mode == "fixed":
        return len(values) == 1
    if mode in {"stepped", "CES", "ramp", "MSE", "bbCID", "all_ions"}:
        return len(values) >= 2
    if mode in {"rolling", "dynamic"}:
        return bool(program.get("equation"))
    return False
