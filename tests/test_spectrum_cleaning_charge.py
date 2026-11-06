"""Charge-domain regression tests for MS/MS peak cleaning."""

from utils.spectrum import clean_peaks


def test_single_charge_removes_products_above_precursor():
    cleaned = clean_peaks(
        [[80.0, 20.0], [150.0, 100.0], [240.0, 80.0]],
        150.0,
        precursor_charge=1,
    )

    assert [peak[0] for peak in cleaned] == [80.0]


def test_double_charge_keeps_legal_higher_mz_product():
    cleaned = clean_peaks(
        [[80.0, 20.0], [150.0, 100.0], [240.0, 80.0], [310.0, 70.0]],
        150.0,
        precursor_charge="2+",
    )

    assert [peak[0] for peak in cleaned] == [80.0, 240.0]
