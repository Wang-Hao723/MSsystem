"""Chemical eligibility, provenance and retrieval regressions for mass repair."""
import copy
import hashlib
import importlib
import json
import sqlite3

import pytest

from matching.reference_mass import (
    DERIVED_NEUTRAL_MASS_SOURCE, verified_neutral_mass_correction,
)
from utils.adduct import ion_mz


def reference(**updates):
    row = dict(smiles='NCCCN', formula='C3H10N2', exact_mass=74.125,
               adduct='[M+H]+', precursor_mz=75.09,
               neutral_mass=74.125, neutral_mass_source='compound_exact_mass')
    return {**row, **updates}


def test_formula_structure_and_reported_ion_repair_average_mass_without_mutation():
    row = reference()
    before = copy.deepcopy(row)
    assert verified_neutral_mass_correction(row) == pytest.approx(74.08439832, abs=1e-8)
    assert row == before


@pytest.mark.parametrize('updates', [
    {'exact_mass': None}, {'exact_mass': float('nan')}, {'exact_mass': float('inf')},
    {'exact_mass': True}, {'exact_mass': -1.}, {'exact_mass': '74.125'},
    {'exact_mass': 74.0844}, {'smiles': ''}, {'smiles': 'invalid smiles'},
    {'formula': ''}, {'formula': 'C4H12N2'},
    {'smiles': '[NH3+]CCCN', 'formula': 'C3H11N2'},
    {'smiles': 'NCCCN.Cl', 'formula': 'C3H11ClN2'},
    {'smiles': 'N[13CH2]CCN'},
    {'adduct': 'Positive'}, {'adduct': 'unknown'}, {'adduct': '[M]+'},
    {'adduct': '[M]-'}, {'adduct': '[M+2H]2+'}, {'adduct': '[2M+H]+'}, {'adduct': '[M+Na]+', 'precursor_mz': 75.09},
    {'precursor_mz': 100.}, {'precursor_mz': float('inf')}, {'precursor_mz': True},
])
def test_incomplete_inconsistent_or_unsupported_reference_keeps_original(updates):
    assert verified_neutral_mass_correction(reference(**updates)) is None


@pytest.mark.parametrize('adduct', ['[M+H]+', '[M+Na]+', '[M-H]-', '[M+Cl]-'])
def test_explicit_supported_ion_is_checked_without_guessing_query_adduct(adduct):
    mass = 74.08439832
    row = reference(adduct=adduct, precursor_mz=ion_mz(mass, adduct))
    assert verified_neutral_mass_correction(row) == pytest.approx(mass, abs=1e-8)


def test_reported_precursor_conflict_cannot_be_hidden_by_search_axis():
    row = reference(reported_precursor_mz=100., search_precursor_mz=75.09167432)
    assert verified_neutral_mass_correction(row) is None


def test_corrected_theoretical_mass_recovers_rounded_reference_inside_query_tolerance():
    search = importlib.import_module('matching.search')
    row = reference()
    mass, source = search._library_record_neutral_mass(row)
    assert source == DERIVED_NEUTRAL_MASS_SOURCE
    row.update(neutral_mass=mass, neutral_mass_source=source)
    score, ppm, basis, query_mass, library_mass = search._candidate_mass_metrics(
        75.09153344922562, row, ppm_tolerance=10., da_tolerance=.01)
    assert basis == 'theoretical_ion_mz'
    assert library_mass == pytest.approx(75.09167432, abs=1e-8)
    assert 1.87 < ppm < 1.88 and score > 0
    _, outside_ppm, _, _, _ = search._candidate_mass_metrics(
        75.1, row, ppm_tolerance=10., da_tolerance=.01)
    assert outside_ppm > 100.
    assert row['exact_mass'] == 74.125


def test_neutral_axis_and_explicit_audited_precursor_keep_their_semantics():
    search = importlib.import_module('matching.search')
    row = reference(neutral_mass=74.08439832, neutral_mass_source=DERIVED_NEUTRAL_MASS_SOURCE)
    _, ppm, basis, _, mass = search._candidate_mass_metrics(
        97.07361632, row, query_neutral_mass=74.08439832, cross_adduct=True)
    assert ppm == 0. and mass == 74.08439832 and basis == 'neutral_mass'
    row.update(search_precursor_mz=75.092, precursor_mz_source='verified_metadata')
    _, ppm, basis, _, mass = search._candidate_mass_metrics(75.092, row)
    assert ppm == 0. and mass == 75.092 and basis == 'verified_metadata'


def test_streaming_and_legacy_loader_share_repaired_mass_and_preserve_database(tmp_path, monkeypatch):
    search = importlib.import_module('matching.search')
    path = tmp_path/'mass.db'
    with sqlite3.connect(path) as db:
        db.executescript('''
            CREATE TABLE Spectra (id INTEGER PRIMARY KEY, spectrum_id TEXT, precursor_mz REAL,
                peaks_json TEXT, adduct TEXT, retention_time REAL, charge INTEGER, source_tier TEXT);
            CREATE TABLE Compounds (id INTEGER PRIMARY KEY, compound_name TEXT, smiles TEXT,
                inchikey TEXT, molecular_formula TEXT, exact_mass REAL, classification TEXT);
            CREATE TABLE Mapping (spectrum_id INTEGER, compound_id INTEGER);
        ''')
        for i, tier in enumerate(['experimental', 'predicted'], 1):
            db.execute('INSERT INTO Spectra VALUES (?,?,?,?,?,?,?,?)',
                (i, str(i), 75.09, json.dumps([[20., 3.], [35., 12.5], [45., 4.]]), '[M+H]+', 1., 1, tier))
            db.execute('INSERT INTO Compounds VALUES (?,?,?,?,?,?,?)',
                (i, '1,3-Diaminopropane', 'NCCCN', 'XFNJVJPLKCPIBV-UHFFFAOYSA-N', 'C3H10N2', 74.125, 'Amine'))
            db.execute('INSERT INTO Mapping VALUES (?,?)', (i, i))
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    monkeypatch.setenv('MSSYSTEM_LIBRARY_RECORD_STORE', 'off')
    for key, value in [('_lib_records', None), ('_lib_db_path', None), ('_mz_index', None),
                       ('_nm_index', None), ('_lib_source_families', ())]:
        monkeypatch.setattr(search, key, value)
    streamed = list(search.iter_normalized_library_records(path))
    loaded, _ = search.load_library(path)
    assert streamed == loaded and len(loaded) == 1
    assert loaded[0]['neutral_mass_source'] == DERIVED_NEUTRAL_MASS_SOURCE
    assert loaded[0]['neutral_mass'] == pytest.approx(74.08439832, abs=1e-8)
    assert loaded[0]['exact_mass'] == 74.125
    assert sum(len(indices) for indices in search._nm_index.values()) == 1
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
