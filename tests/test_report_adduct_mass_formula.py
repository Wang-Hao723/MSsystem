"""Check ion stoichiometry and the report consumer, without sample data."""
import pytest
from pipeline import adduct_hypotheses as ion
from web.report_common import _resolved_adduct


def test_water_loss_is_seventeen_da_not_ten_millida():
    correct=300+ion.PROTON-ion.WATER
    answer=ion.detect_adduct(correct,neutral_mass_reference=300,tol_da=1e-6)
    assert answer['adduct'] == '[M-H2O+H]+'
    assert answer['status'] == 'assigned'
    assert len(answer['hypotheses']) == 1
    assert ion.detect_adduct(299.989476,neutral_mass_reference=300,tol_da=1e-5)['status'] == 'no_match'


@pytest.mark.parametrize('adduct,shift', [('[M+Cl]-',34.969401261909),('[M-H-H2O]-',-19.017841150609)])
def test_negative_ion_electron_and_water_mass(adduct,shift):
    result=ion.detect_adduct(300+shift,mode='negative',neutral_mass_reference=300,tol_da=1e-7)
    assert result['adduct'] == adduct
    assert result['candidate_independent'] is False


def test_mass_boundary_is_not_rounded_into_acceptance():
    result=ion.detect_adduct(300+ion.PROTON+.0010004,neutral_mass_reference=300,tol_da=.001)
    assert result['status'] == 'no_match'


def test_report_consumer_exposes_one_alias_as_candidate_only():
    from matching.rule_quality import formula_exact_mass
    mass=formula_exact_mass('C6H12O6')
    row={'precursor_mz':mass+ion.PROTON-ion.WATER,'ion_mode':'POS','precursor_charge':1}
    name,source,_=_resolved_adduct({'formula':'C6H12O6'},row)
    assert name == '[M-H2O+H]+'
    assert source == 'candidate_formula_precursor_unique_10ppm'


def test_invalid_numbers_are_unavailable():
    assert ion.detect_adduct(float('nan')) is None
    assert ion.detect_adduct(300,neutral_mass_reference=float('inf')) is None
