"""Collision-energy metadata helpers.

Units are preserved explicitly because a bare value cannot be assumed to be
electronvolts or normalized collision energy across instruments.
"""

from __future__ import annotations

import math
import re
from typing import Any, Mapping


# Collision energy is non-negative.  Deliberately omit a leading sign so the
# hyphen in ``19.9-29.9 eV`` is interpreted as a range separator, not ``-29.9``.
_NUMBER_RE = re.compile(r"(?:\d+(?:\.\d*)?|\.\d+)")


def canonical_collision_energy_unit(value: Any) -> str:
    """Return a conservative CE unit without converting unlike quantities.

    ``V`` is collision voltage, not collision energy in eV, and a generic
    percentage is not assumed to be vendor-normalized collision energy.  They
    remain separate domains so downstream code cannot compare them by accident.
    """
    text = str(value or "").strip().lower()
    compact = re.sub(r"[\s_-]+", "", text)
    if not compact:
        return ""
    if compact in {"ev", "electronvolt", "electronvolts"}:
        return "eV"
    if compact in {
        "nce",
        "normalized",
        "normalizedcollisionenergy",
        "normalisedcollisionenergy",
    }:
        return "NCE"
    if compact in {"v", "volt", "volts"}:
        return "V"
    if compact in {"%", "percent", "percentage", "nominalpercent"}:
        return "percent"
    return ""


def parse_collision_energy(value: Any, unit_hint: Any = None) -> dict[str, Any]:
    """Parse heterogeneous MassBank/vendor CE text while preserving provenance.

    Examples include ``10(NCE)``, ``20 eV``, ``6V``, ``Ramp 19.9-29.9 eV``
    and ``30 % (nominal)``.  The result is deliberately descriptive: only a
    single value with eV/NCE is marked directly comparable to a query CE.
    """
    raw = "" if value is None else str(value).strip()
    hinted_unit = canonical_collision_energy_unit(unit_hint)
    if not raw:
        return {
            "raw": "",
            "value": None,
            "unit": hinted_unit,
            "mode": "missing",
            "minimum": None,
            "maximum": None,
            "values": [],
            "comparable": False,
        }

    lower = raw.lower()
    detected_unit = ""
    if re.search(r"\bnce\b|normali[sz]ed(?:\s+collision\s+energy)?", lower):
        detected_unit = "NCE"
    elif re.search(
        r"(?:^|[^A-Za-z])e\s*v\b|electron\s*volts?",
        lower,
        flags=re.IGNORECASE,
    ):
        detected_unit = "eV"
    elif "%" in raw or re.search(r"\bpercent(?:age)?\b", lower):
        detected_unit = "percent"
    elif re.search(r"(?<![A-Za-z])v\b|\bvolts?\b", lower, flags=re.IGNORECASE):
        detected_unit = "V"
    unit = detected_unit or hinted_unit

    numbers: list[float] = []
    for token in _NUMBER_RE.findall(raw):
        try:
            number = float(token)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number) and 0 <= number < 10000:
            numbers.append(number)

    if not numbers:
        mode = "unparsed"
    elif re.search(r"\bramp(?:ed)?\b", lower):
        mode = "ramp"
    elif re.search(r"\bstep(?:ped)?\b", lower):
        mode = "stepped"
    elif len(numbers) >= 2 and re.search(
        r"\d(?:\.\d+)?\s*(?:-|–|—|to)\s*\d", lower
    ):
        mode = "range"
    elif len(numbers) > 1:
        mode = "stepped"
    else:
        mode = "single"

    minimum = min(numbers) if numbers else None
    maximum = max(numbers) if numbers else None
    single_value = numbers[0] if mode == "single" and numbers else None
    return {
        "raw": raw,
        "value": single_value,
        "unit": unit,
        "mode": mode,
        "minimum": minimum,
        "maximum": maximum,
        "values": numbers,
        "comparable": bool(
            single_value is not None and unit in {"eV", "NCE"}
        ),
    }


def collision_energy_evidence(record: Mapping[str, Any]) -> Any:
    """Return structured CE when the unit is known, otherwise the raw value."""
    value = record.get("collision_energy")
    if value is None:
        return None
    unit = str(record.get("collision_energy_unit") or "").strip()
    if unit in {"eV", "NCE"}:
        return {"value": value, "unit": unit}
    return value
