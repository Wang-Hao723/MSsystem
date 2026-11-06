"""Source-binding failures and actual API loader wiring must remain observable."""
import copy
import hashlib
import importlib
import json
from pathlib import Path
import pytest
from runtime import library_metadata_repairs as M

class Base(list):
    def read_uncached(self, index):
        return self[index]

@pytest.fixture
def bound(tmp_path):
    source = tmp_path / 'source.mgf'
    source.write_text('CHARGE=1-\n')
    base = Base([
        {'spectrum_id': 'control', 'polarity': 'POS', 'charge': 1, 'peaks': [[100, 1]]},
        {'spectrum_id': 'repair', 'polarity': 'POS', 'charge': 1, 'peaks': [[101, 1]], 'precursor_mz': 202},
    ])
    base.manifest = {'normalized_records_sha256': 'base', 'record_count': 2, 'source_db': {'sha256': 'db'}}
    registry = {'schema_version': 'mssystem.source-metadata-repairs/1.0.0', 'registry_id': 'test',
        'base_normalized_records_sha256': 'base', 'base_record_count': 2, 'source_database_sha256': 'db',
        'source_files': {str(source): {'size_bytes': source.stat().st_size, 'sha256': hashlib.sha256(source.read_bytes()).hexdigest()}},
        'rows': [{'record_index': 1, 'spectrum_id': 'repair', 'original_record_sha256': M.digest(base[1]),
                  'updates': {'charge': -1, 'source_ion_mode': 'NEGATIVE', 'adduct': '[M-H]-'}, 'source_polarity': 'NEG',
                  'source_binding': {'source_file': str(source), 'block_sha256': 'block', 'structure_key_equal': True,
                                     'raw_peaks_equal': True, 'existing_importer_peak_transform_equal': False}}]}
    seal(registry)
    return base, registry, source

def seal(registry):
    registry['semantic_sha256'] = M.digest({k:v for k,v in registry.items() if k != 'semantic_sha256'})

def test_restore_only_bound_fields_and_keep_membership(bound):
    base, registry, _ = bound
    before = copy.deepcopy(list(base))
    fixed = M.VerifiedMetadataRecords(base, registry)
    assert len(fixed) == len(base)
    assert list(base) == before
    assert fixed[0] is base[0]
    assert fixed[1]['charge'] == -1 and fixed[1]['polarity'] == 'NEG'
    assert fixed[-1] == fixed[1] == fixed.read_uncached(1)
    assert fixed[:][1] == fixed[1]
    assert fixed[1]['peaks'] == before[1]['peaks']
    assert fixed[1]['precursor_mz'] == before[1]['precursor_mz']
    with pytest.raises(IndexError): fixed[-3]
    assert M.library_metadata_repair_status(fixed, True)['applied_reference_count'] == 1
    assert M.library_metadata_repair_status(fixed, False)['applied_reference_count'] == 0

@pytest.mark.parametrize('damage', ['registry', 'base', 'source', 'original', 'duplicate', 'peaks_update', 'wrong_polarity'])
def test_reject_tampered_or_incompatible_bindings(bound, damage):
    base, reg, source = bound
    if damage == 'registry': reg['registry_id'] = 'tampered'
    if damage == 'base': base.manifest['normalized_records_sha256'] = 'other'
    if damage == 'source': source.write_text('CHARGE=1+\n')
    if damage == 'original': base[1]['precursor_mz'] += 1
    if damage == 'duplicate': reg['rows'].append(copy.deepcopy(reg['rows'][0])); seal(reg)
    if damage == 'peaks_update': reg['rows'][0]['updates']['peaks'] = []; seal(reg)
    if damage == 'wrong_polarity': reg['rows'][0]['source_polarity'] = 'POS'; seal(reg)
    with pytest.raises(ValueError): M.VerifiedMetadataRecords(base, reg)

def test_full_loader_shared_identity_and_small_unchanged(bound, tmp_path, monkeypatch):
    # This fixture registers only spectral metadata. The combined mass+metadata
    # path has its own source-bound atomic-publication tests.
    from runtime import source_compound_mass_repairs as mass_repairs
    monkeypatch.setattr(mass_repairs, '_STATES', {})
    monkeypatch.setattr(mass_repairs, '_read_config', lambda root: None)
    import config.settings as settings
    S = importlib.import_module('matching.search')
    base, reg, _ = bound
    full = tmp_path / 'full.db'; small = tmp_path / 'small.db'; index = {'same': True}
    path = tmp_path / 'registry.json'; path.write_text(json.dumps(reg))
    monkeypatch.setattr(M, 'REGISTRY_PATH', path)
    monkeypatch.setattr(settings, 'DB_PATH', full)
    monkeypatch.setattr(S, '_lib_records', base)
    monkeypatch.setattr(S, 'load_library', lambda _: (S._lib_records, index))
    fixed, outindex = M.load_library_with_verified_repairs(full)
    assert outindex is index and S._lib_records is fixed
    assert M.load_library_with_verified_repairs(full)[0] is fixed
    monkeypatch.setattr(S, '_lib_records', base)
    assert M.load_library_with_verified_repairs(small)[0] is base
    base.manifest['record_count'] = 3
    with pytest.raises(ValueError): M.load_library_with_verified_repairs(full)
    assert S._lib_records is None

def test_web_uses_verified_loader_and_reports_actual_application():
    from web import server
    assert server.load_library is M.load_library_with_verified_repairs
    import inspect
    assert 'library_metadata_repair_status(_lib_records, _lib_loaded)' in inspect.getsource(server.health_check)

def test_preflight_detects_source_change_before_load(bound, tmp_path):
    base, reg, source = bound
    config = tmp_path / 'config'; config.mkdir()
    cache = tmp_path / 'cache'; cache.mkdir()
    (config/'library_metadata_repairs.json').write_text(json.dumps(reg))
    (config/'library_record_store.json').write_text(json.dumps({'enabled': True, 'cache_dir': str(cache)}))
    (cache/'manifest.json').write_text(json.dumps(base.manifest))
    assert M.metadata_repair_preflight(tmp_path)['compatible'] is True
    source.write_text('CHARGE=1+\n')
    assert M.metadata_repair_preflight(tmp_path)['compatible'] is False
