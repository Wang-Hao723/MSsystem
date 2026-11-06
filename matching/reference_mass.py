"""Derive a neutral monoisotopic mass from mutually consistent source fields.

The source's exact_mass value remains untouched. A repaired neutral axis is
eligible only when a neutral, single-component, unlabeled molecular structure
agrees with the submitted formula and the explicitly specified precursor ion.
"""
from functools import lru_cache
import math
import re

from matching.library_quality import formula_exact_mass, MONOISOTOPIC_MASS
from utils.adduct import NEUTRAL_SHIFTS, ion_mz

DERIVED_NEUTRAL_MASS_SOURCE = 'structure_formula_verified_exact_mass'
_SPECIFIC_ADDUCTS = {
    key.upper().replace(' ', '') for key in NEUTRAL_SHIFTS
} - {'[M]+', '[M]-', 'M+', 'M-'}


def _positive_number(value):
    return (not isinstance(value, bool) and isinstance(value, (int, float))
            and math.isfinite(value) and value > 0)


# Formula identity is an elemental-composition test. A second atomic-mass
# table cannot independently establish identity and its rounding must not veto
# an otherwise equal formula (for example, several 127I atoms).
_FORMULA_ELEMENT = re.compile(r"([A-Z][a-z]?)([0-9]*)", re.ASCII)
_SUBSCRIPTS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def _complete_formula_counts(formula):
    """Parse a whole neutral formula; never discard hydrate/salt suffixes."""
    text = re.sub(r"\s+", "", str(formula or "").translate(_SUBSCRIPTS))
    counts, end = {}, 0
    for token in _FORMULA_ELEMENT.finditer(text):
        if token.start() != end:
            return None
        element, digits = token.groups()
        count = int(digits) if digits else 1
        if element not in MONOISOTOPIC_MASS or count <= 0:
            return None
        counts[element] = counts.get(element, 0) + count
        end = token.end()
    return counts if counts and end == len(text) else None



@lru_cache(maxsize=65536)
def _formula_mass(formula):
    return formula_exact_mass(formula)


@lru_cache(maxsize=131072)
def _verified_structure_mass(smiles, formula, reported_mass):
    counts = _complete_formula_counts(formula)
    if counts is None:
        return None
    formula_mass = _formula_mass(formula)
    if not _positive_number(formula_mass):
        return None
    if abs(reported_mass-formula_mass) <= max(.0001, formula_mass*10e-6):
        return None
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    mol = Chem.MolFromSmiles(smiles)
    if (mol is None or Chem.GetFormalCharge(mol) != 0
            or len(Chem.GetMolFrags(mol)) != 1
            or any(atom.GetIsotope() for atom in mol.GetAtoms())
            or any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())):
        return None
    if _complete_formula_counts(rdMolDescriptors.CalcMolFormula(mol)) != counts:
        return None
    exact = rdMolDescriptors.CalcExactMolWt(mol)
    if not _positive_number(exact):
        return None
    return exact


def verified_neutral_mass_correction(record):
    """Return a justified correction, or None to preserve the original axis.

    The 50 ppm / 0.02 Da test checks consistency within the reference record;
    it does not widen the query's identity-search tolerance. It cannot promote
    a candidate to a higher evidence level by itself.
    """
    reported_mass = record.get('exact_mass')
    smiles = str(record.get('smiles') or '').strip()
    formula = str(record.get('formula') or '').strip()
    if not _positive_number(reported_mass) or not smiles or not formula:
        return None
    adduct = str(record.get('adduct') or '').strip()
    if adduct.upper().replace(' ', '') not in _SPECIFIC_ADDUCTS:
        return None
    observed = record.get('reported_precursor_mz')
    if observed is None:
        observed = record.get('precursor_mz')
    if not _positive_number(observed):
        return None
    exact = _verified_structure_mass(smiles, formula, reported_mass)
    if exact is None:
        return None
    ion = ion_mz(exact, adduct)
    if not _positive_number(ion) or abs(observed-ion) > max(.02, ion*50e-6):
        return None
    return exact
