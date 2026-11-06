"""Regression cases for whole-formula identity and multiply iodinated sources."""
from copy import deepcopy

import pytest

from matching.reference_mass import verified_neutral_mass_correction


@pytest.mark.parametrize("formula", [
    "C3H10N2.H2O", "C3H10N2·H2O", "C3H10N2.garbage", "C3H10N2..",
    "C3H10N0", "C3H10N2+", "C3H10N2^0", "C3H10N2Xx", "C3H10N2.0",
    "C3H10N2(-)", "[13C]3H10N2",
])
def test_formula_suffix_and_different_composition_cannot_authorize_repair(formula):
    record = dict(smiles="NCCCN", formula=formula, exact_mass=74.125,
                  adduct="[M+H]+", precursor_mz=75.09)
    assert verified_neutral_mass_correction(record) is None


@pytest.mark.parametrize("formula", [
    "C3H10N2", "H10N2C3", "C₃H₁₀N₂", " C3 H10 N2 ", "CH4CH4CH2N2",
])
def test_equivalent_whole_formula_is_an_element_count_comparison(formula):
    record = dict(smiles="NCCCN", formula=formula, exact_mass=74.125,
                  adduct="[M+H]+", precursor_mz=75.09167477)
    assert verified_neutral_mass_correction(record) == pytest.approx(74.08439832, abs=1e-9)


@pytest.mark.parametrize("smiles,formula,reported,adduct,precursor,expected", [
    ("O=C(NCC(O)CO)c1c(I)c(C(=O)NCC(O)CO)c(I)c(N(CCO)C(=O)CO)c1I",
     "C18H24I3N3O9", 6.864673, "[M+H]+", 807.87195, 806.864673348),
    ("O=C1OC2(c3ccccc31)c1cc(I)c(O)c(I)c1Oc1c2cc(I)c(O)c1I",
     "C20H8I4O5", 35.655065, "[M+H]+", 836.66234, 835.655065356),
    ("O=C(O)c1c(I)cc(I)c(O)c1I", "C7H3I3O3", 15.721638,
     "[M-H]-", 514.71436, 515.721637956),
    ("CC(=O)Nc1c(I)cc(I)c(C(=O)O)c1I", "C9H6I3NO3", 56.748187,
     "[M+H]+", 557.75546, 556.748187052),
])
def test_iodine_table_rounding_does_not_veto_equal_composition(
        smiles, formula, reported, adduct, precursor, expected):
    record = dict(smiles=smiles, formula=formula, exact_mass=reported,
                  adduct=adduct, precursor_mz=precursor)
    original = deepcopy(record)
    assert verified_neutral_mass_correction(record) == pytest.approx(expected, abs=1e-8)
    assert record == original
    assert verified_neutral_mass_correction({**record, "formula": ""}) is None
    assert verified_neutral_mass_correction({**record, "precursor_mz": precursor + 1}) is None


def test_neutral_radical_is_not_a_closed_shell_repair_source():
    record = dict(smiles="[CH3]", formula="CH3", exact_mass=15.5,
                  adduct="[M+H]+", precursor_mz=16.0307515)
    assert verified_neutral_mass_correction(record) is None
