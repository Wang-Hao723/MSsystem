"""Normalize instrument metadata used by condition-aware MS/MS comparisons."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


_FAMILY_PATTERNS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("qtof", (
        "qtof", "q-tof", "quadrupole time-of-flight",
        "quadrupole time of flight",
    )),
    ("orbitrap", (
        "orbitrap", "q exactive", "q-exactive", "exploris",
        "fusion lumos", "orbitrap eclipse", "orbitrap ascend",
    )),
    ("fticr", (
        "ft-icr", "fticr", "fourier transform ion cyclotron",
    )),
    ("triple_quadrupole", (
        "triple quadrupole", "triplequad", "qqq",
    )),
    ("ion_trap", (
        "ion trap", "linear trap", "ltq",
    )),
    ("tof", ("time-of-flight", "time of flight", " tof")),
    ("quadrupole", ("quadrupole",)),
    ("sector", ("magnetic sector", "sector")),
)

_GENERIC_COMPONENT_NAMES = {
    "instrument model",
    "mass spectrometer file format",
    "electrospray ionization",
    "nanoelectrospray",
    "atmospheric pressure chemical ionization",
    "quadrupole",
    "time-of-flight",
    "time of flight",
    "orbitrap",
    "ion trap",
    "electron multiplier",
    "microchannel plate detector",
}


def normalize_instrument_family(*values: Any) -> str:
    """Map heterogeneous vendor/model labels to a conservative analyzer family."""
    text = " ".join(str(value or "") for value in values).strip().lower()
    compact = re.sub(r"[^a-z0-9]+", "", text)
    if not compact:
        return ""
    for family, patterns in _FAMILY_PATTERNS:
        for pattern in patterns:
            normalized = re.sub(r"[^a-z0-9]+", "", pattern.lower())
            if normalized and normalized in compact:
                return family
    return ""


def _local_name(tag: str) -> str:
    return str(tag or "").rsplit("}", 1)[-1]


def extract_mzml_instrument_metadata(file_path: str | Path) -> dict[str, str]:
    """Read the first mzML instrumentConfiguration without scanning spectra.

    mzML instrument metadata precedes the spectrum list. ``iterparse`` stops at
    the first configuration, so binary peak arrays are not decoded.
    """
    path = Path(file_path)
    if not path.exists():
        return {"instrument_type": "", "instrument": ""}
    try:
        context = ET.iterparse(str(path), events=("end",))
        for _, elem in context:
            if _local_name(elem.tag) != "instrumentConfiguration":
                continue
            names: list[str] = []
            for child in elem.iter():
                if _local_name(child.tag) != "cvParam":
                    continue
                name = str(child.attrib.get("name") or "").strip()
                if name:
                    names.append(name)
            family = normalize_instrument_family(*names)
            model_candidates = [
                name for name in names
                if name.strip().lower() not in _GENERIC_COMPONENT_NAMES
                and (
                    any(char.isdigit() for char in name)
                    or any(
                        brand in name.lower()
                        for brand in (
                            "agilent", "bruker", "thermo", "waters",
                            "sciex", "shimadzu", "jeol",
                        )
                    )
                )
            ]
            instrument = (
                max(model_candidates, key=len)
                if model_candidates else ""
            )
            elem.clear()
            return {
                "instrument_type": family,
                "instrument": instrument,
            }
    except (ET.ParseError, OSError, ValueError):
        pass
    return {"instrument_type": "", "instrument": ""}
