import pytest

from utils.adduct import NEGATIVE_ADDUCTS, POSITIVE_ADDUCTS, neutral_mass


@pytest.mark.parametrize('mass', [100., 200., 500., 1000.])
@pytest.mark.parametrize('adduct,shift', [
    ('[M+Cl]-', 34.968852682 + 0.000548579909),
    ('[M-H-H2O]-', -(1.007276466621 + 2 * 1.00782503223 + 15.99491461957)),
])
def test_ionic_atomic_mass_roundtrip(mass, adduct, shift):
    assert neutral_mass(mass + shift, adduct) == pytest.approx(mass, abs=1e-9)


def test_public_dehydration_candidate_uses_proton_minus_water():
    shifts = {name: shift for name, shift, charge in POSITIVE_ADDUCTS if charge == 1}
    water = 2 * 1.00782503223 + 15.99491461957
    expected = 1.007276466621 - water
    assert shifts['[M+H-H2O]+'] == pytest.approx(expected, abs=1e-9)


def test_chloride_candidate_table_includes_electron_mass():
    shifts = {name: shift for name, shift, charge in NEGATIVE_ADDUCTS if charge == 1}
    assert shifts['[M+Cl]-'] == pytest.approx(34.968852682 + .000548579909, abs=1e-12)
