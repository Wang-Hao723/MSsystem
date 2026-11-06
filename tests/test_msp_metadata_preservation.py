"""MSP acquisition metadata must survive import without inventing energy units."""
import pytest
import build_from_mgf_msp as builder


def load(tmp_path, monkeypatch, extra):
    monkeypatch.setattr(builder, '_resolve_classification_safely', lambda **kw: 'other')
    p = tmp_path / 'reference.msp'
    p.write_text('NAME: reference\nPRECURSORMZ: 289.0718\n'
                 'PRECURSORTYPE: [M-H]-\nIONMODE: Negative\nFORMULA: C15H14O6\n'
                 'INCHIKEY: FSYDWKPCKNCRDI-UHFFFAOYSA-N\n'
                 'Comment: DB#=MSBNK-Eawag-EQ01200358; origin=MassBank\n'
                 + extra + '\nNum Peaks: 3\n55.019 100\n108.0216 90\n109.0296 80\n')
    records, count = builder.parse_msp_file(p)
    assert count == 1
    rec = records[0]
    assert rec['inchikey'] == 'FSYDWKPCKNCRDI-UHFFFAOYSA-N'
    assert rec['adduct'] == '[M-H]-'
    assert rec['source_accession'] == 'MSBNK-Eawag-EQ01200358'
    assert len(rec['peaks']) == 3
    return rec


@pytest.mark.parametrize('raw,unit,value,mode,lo,hi', [
    ('150 % (nominal)', 'percent', 150, 'single', 150, 150),
    ('35', '', 35, 'single', 35, 35),
    ('6V', 'V', 6, 'single', 6, 6),
    ('Ramp 19.9-29.9 eV', 'eV', None, 'ramp', 19.9, 29.9),
    ('10 [NCE]', 'NCE', 10, 'single', 10, 10),
])
def test_energy_and_separate_instrument_fields(tmp_path, monkeypatch, raw, unit, value, mode, lo, hi):
    r = load(tmp_path, monkeypatch, 'INSTRUMENTTYPE: LC-ESI-QFT\n'
             'INSTRUMENT: Exploris 240 Orbitrap Thermo Scientific\nCOLLISIONENERGY: ' + raw)
    assert r['instrument_type'] == 'LC-ESI-QFT'
    assert r['instrument'] == 'Exploris 240 Orbitrap Thermo Scientific'
    assert r['collision_energy_raw'] == raw
    assert r['collision_energy_unit'] == unit
    assert r['collision_energy'] == value
    assert r['collision_energy_mode'] == mode
    assert (r['collision_energy_min'], r['collision_energy_max']) == (lo, hi)


def test_explicit_energy_unit(tmp_path, monkeypatch):
    r = load(tmp_path, monkeypatch, 'Collision_energy: 20\nCollision_energy_unit: eV')
    assert (r['collision_energy'], r['collision_energy_unit']) == (20, 'eV')


def test_missing_energy_stays_missing(tmp_path, monkeypatch):
    r = load(tmp_path, monkeypatch, 'INSTRUMENTTYPE: LC-ESI-QFT')
    assert r['collision_energy'] is None
    assert r['collision_energy_unit'] == ''
    assert r['collision_energy_mode'] == 'missing'
    assert r['instrument'] is None


def test_model_field_cannot_hide_predicted_origin(tmp_path, monkeypatch):
    r = load(tmp_path, monkeypatch, 'INSTRUMENTTYPE: LC-ESI-QTOF\nINSTRUMENT: CFM-ID predicted')
    assert r['source_tier'] == 'predicted'
