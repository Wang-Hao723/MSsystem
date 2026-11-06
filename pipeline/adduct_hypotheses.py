"""Formula-conditioned report hypotheses using explicit monoisotopic ion masses.

These hypotheses do not assign a sample adduct or modify the accepted library
mass index. Atom masses: NIST isotopic compositions; proton/electron: CODATA
2022. Aliases with identical stoichiometry represent one ion hypothesis.
"""
from __future__ import annotations

import math

H = 1.00782503223
O = 15.99491461957
N = 14.00307400443
ELECTRON = 0.0005485799090441
PROTON = 1.0072764665789
WATER = 2 * H + O

POSITIVE = (
    ("[M+H]+", PROTON, 1),
    ("[M+Na]+", 22.9897692820 - ELECTRON, 1),
    ("[M+K]+", 38.9637064864 - ELECTRON, 1),
    ("[M+NH4]+", N + 4 * H - ELECTRON, 1),
    ("[M-H2O+H]+", PROTON - WATER, 1),
    ("[M+2H]2+", 2 * PROTON, 2),
)
NEGATIVE = (
    ("[M-H]-", -PROTON, 1),
    ("[M+Cl]-", 34.968852682 + ELECTRON, 1),
    ("[M+FA-H]-", 12 + H + 2 * O + ELECTRON, 1),
    ("[M+Hac-H]-", 24 + 3 * H + 2 * O + ELECTRON, 1),
    ("[M-H-H2O]-", -PROTON - WATER, 1),
)


def detect_adduct(precursor_mz, peaks=None, tol_da=.05, mode="positive", neutral_mass_reference=None):
    """Return formula/mass-compatible hypotheses; absence is not structural evidence."""
    try:
        mz, tolerance = float(precursor_mz), float(tol_da)
        reference = None if neutral_mass_reference is None else float(neutral_mass_reference)
    except (TypeError, ValueError, OverflowError):
        return None
    if not math.isfinite(mz) or mz <= 0 or not math.isfinite(tolerance) or tolerance < 0:
        return None
    if reference is not None and (not math.isfinite(reference) or reference <= 0):
        return None
    forms = POSITIVE if mode in {"positive", "POS"} else NEGATIVE if mode in {"negative", "NEG"} else ()
    hypotheses = []
    for name, shift, charge in forms:
        neutral = mz * charge - shift
        if not 0 < neutral <= 2000:
            continue
        error = None if reference is None else abs(neutral - reference)
        # Compare before rounding, otherwise a boundary error can change whether
        # a report claims a unique mass-compatible hypothesis.
        if error is None or error <= tolerance:
            hypotheses.append({"adduct": name, "neutral_mass": round(neutral, 9),
                               "mass_error_da": error, "charge": charge})
    assigned = reference is not None and len(hypotheses) == 1
    chosen = hypotheses[0] if assigned else {}
    return {"status": "assigned" if assigned else "ambiguous" if hypotheses else "no_match",
            "adduct": chosen.get("adduct"), "neutral_mass": chosen.get("neutral_mass"),
            "mass_error_da": chosen.get("mass_error_da"), "charge": chosen.get("charge"),
            "confidence": "hypothesis", "hypotheses": hypotheses,
            "source": "monoisotopic_formula_precursor_hypothesis_v1",
            "candidate_independent": False}
