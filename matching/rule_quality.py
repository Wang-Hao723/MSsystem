"""Shared quality gates for expert-rule data.

The expert libraries contain charged fragment formulae such as ``C5H15NO4P+``.
This module intentionally performs only a conservative monoisotopic-mass check;
it does not attempt formula inference or silently repair curated records.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Any, Optional


_MONOISOTOPIC_MASS = {
    "H": 1.00782503223,
    "B": 11.00930536,
    "C": 12.0,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "F": 18.99840316273,
    "Na": 22.9897692820,
    "Mg": 23.985041697,
    "Al": 26.98153853,
    "Si": 27.97692653465,
    "P": 30.97376199842,
    "S": 31.9720711744,
    "Cl": 34.968852682,
    "K": 38.9637064864,
    "Ca": 39.962590863,
    "Fe": 55.93493633,
    "Br": 78.9183376,
    "I": 126.904473,
}

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def normalize_polarity(value: object) -> Optional[str]:
    """Return ``positive``/``negative`` or ``None`` for unknown values."""
    text = str(value or "").strip().lower()
    if text in {"pos", "positive", "+", "+1", "1+"}:
        return "positive"
    if text in {"neg", "negative", "-", "-1", "1-"}:
        return "negative"
    return None


def formula_exact_mass(formula: str) -> Optional[float]:
    """Calculate a simple elemental monoisotopic mass.

    Electron-mass corrections (<0.001 Da for singly charged ions) are omitted;
    the production quality gate uses a 0.01 Da tolerance.  Isotopically labelled,
    parenthesised, or otherwise non-elemental formulae return ``None`` rather than
    being partially parsed.
    """
    if not isinstance(formula, str) or not formula.strip():
        return None
    text = formula.strip().replace(" ", "")
    if text.startswith("[") and "]" in text:
        closing = text.rfind("]")
        text = text[1:closing] + text[closing + 1:]
    # Strip only the charge suffix.  Do not consume the final element count
    # (e.g. the ``7`` in C7H7+).
    text = re.sub(r"(?:\^\d*[+-]|[+-]\d*)$", "", text)
    if not text:
        return None

    mass = 0.0
    cursor = 0
    for match in _TOKEN.finditer(text):
        if match.start() != cursor:
            return None
        element, count_text = match.groups()
        if element not in _MONOISOTOPIC_MASS:
            return None
        count = int(count_text) if count_text else 1
        mass += _MONOISOTOPIC_MASS[element] * count
        cursor = match.end()
    if cursor != len(text):
        return None
    return mass


def formula_mz_error(formula: str, mz: float) -> Optional[float]:
    """Return observed minus formula mass, or ``None`` if not calculable."""
    exact_mass = formula_exact_mass(formula)
    if exact_mass is None:
        return None
    return float(mz) - exact_mass


def formula_mz_consistent(formula: str, mz: float, tolerance: float = 0.01) -> bool:
    """Whether a singly charged fragment formula agrees with its curated m/z."""
    error = formula_mz_error(formula, mz)
    return error is not None and abs(error) <= tolerance


@lru_cache(maxsize=20000)
def smiles_matches_any(smiles: str, smarts_patterns: tuple[str, ...]) -> bool:
    """Whether a valid candidate structure matches at least one SMARTS guard."""
    if not smiles or not smarts_patterns:
        return False
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        for smarts in smarts_patterns:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None and mol.HasSubstructMatch(pattern):
                return True
    except Exception:
        return False
    return False


@lru_cache(maxsize=20000)
def smiles_matches_all(smiles: str, smarts_patterns: tuple[str, ...]) -> bool:
    """Whether a valid candidate structure matches every SMARTS guard."""
    if not smiles or not smarts_patterns:
        return False
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        for smarts in smarts_patterns:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None or not mol.HasSubstructMatch(pattern):
                return False
        return True
    except Exception:
        return False


def normalize_adduct(value: object) -> str:
    """Return the canonical spelling of a recognized adduct, else ``""``.

    Rule applicability must never infer a specific adduct from polarity alone.
    The neutral-mass table is reused only as an alias catalogue; the original
    annotation must identify a concrete ion species.
    """
    text = str(value or "").strip()
    if not text:
        return ""
    compact = text.upper().replace(" ", "")
    aliases = {
        "[M]+": "[M]+", "M+": "[M]+",
        "[M+H]+": "[M+H]+", "M+H": "[M+H]+",
        "[M+NA]+": "[M+Na]+", "M+NA": "[M+Na]+",
        "[M+K]+": "[M+K]+", "M+K": "[M+K]+",
        "[M+NH4]+": "[M+NH4]+", "M+NH4": "[M+NH4]+",
        "[M+H-H2O]+": "[M+H-H2O]+",
        "[M-H2O+H]+": "[M+H-H2O]+",
        "[M]-": "[M]-", "M-": "[M]-",
        "[M-H]-": "[M-H]-", "M-H": "[M-H]-",
        "[M+CL]-": "[M+Cl]-", "M+CL": "[M+Cl]-",
        "[M+FA-H]-": "[M+FA-H]-",
        "[M+HCOO]-": "[M+FA-H]-",
        "[M+HCOOH-H]-": "[M+FA-H]-",
        "[M+HAC-H]-": "[M+Hac-H]-",
        "[M+CH3COO]-": "[M+Hac-H]-",
        "[M+CH3COOH-H]-": "[M+Hac-H]-",
        "[M-H-H2O]-": "[M-H-H2O]-",
    }
    if compact in aliases:
        return aliases[compact]
    return text


def collision_energy_value_unit(value: Any) -> tuple[Optional[float], str]:
    """Parse collision energy without pretending that unitless values are eV.

    Supported inputs are ``{"value": 30, "unit": "eV"}``, ``"30 eV"`` and
    ``"25 NCE"``.  A bare number is returned with an empty unit and therefore
    cannot activate a unit-conditioned rule.
    """
    unit = ""
    raw = value
    if isinstance(value, dict):
        raw = value.get("value", value.get("collision_energy"))
        unit = str(value.get("unit") or "").strip()
    elif isinstance(value, str):
        match = re.search(
            r"([-+]?\d+(?:\.\d+)?)\s*(nce|normalized(?:\s+collision\s+energy)?|ev|%)?",
            value,
            flags=re.IGNORECASE,
        )
        if not match:
            return None, ""
        raw = match.group(1)
        unit = str(match.group(2) or "").strip()
    try:
        number = float(raw)
    except (TypeError, ValueError):
        return None, ""
    if not (number >= 0.0 and number < 10000.0):
        return None, ""
    normalized_unit = unit.lower().replace(" ", "_")
    if normalized_unit in {"ev"}:
        normalized_unit = "eV"
    elif normalized_unit in {"nce", "normalized", "normalized_collision_energy", "%"}:
        normalized_unit = "NCE"
    else:
        normalized_unit = ""
    return number, normalized_unit
