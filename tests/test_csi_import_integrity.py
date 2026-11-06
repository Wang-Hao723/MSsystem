"""Reject incomplete and stale CSI evidence before native service calls."""
import importlib.util
import json
from pathlib import Path
import pytest

MODULE = Path(__file__).resolve().parents[1] / 'eval/msystem_method_breakthrough_20260908/csi_import_integrity.py'
spec = importlib.util.spec_from_file_location('tested_csi_integrity', MODULE)
integrity = importlib.util.module_from_spec(spec)
spec.loader.exec_module(integrity)


@pytest.fixture
def bundle(tmp_path, monkeypatch):
    inp, db, protocol, jar = [tmp_path / n for n in ['input.tsv', 'native.siriusdb', 'protocol.json', 'native.jar']]
    inp.write_text('CCO\tAAAAAAAAAAAAAA\tname\nCCN\tBBBBBBBBBBBBBB\tname\n')
    db.write_bytes(b'native database snapshot')
    protocol.write_text('{}')
    jar.write_bytes(b'native program snapshot')
    monkeypatch.setattr(integrity, 'native_jars', lambda: [jar])
    qa = tmp_path / 'IMPORT_QA.json'
    fp = tmp_path / 'IMPORT_QA_FINGERPRINTS.jsonl'
    fp.write_text(''.join(json.dumps({'raw14': k, 'set_bits': 0, 'bits_sha256': 'a' * 64}) + '\n'
                          for k in ['AAAAAAAAAAAAAA', 'BBBBBBBBBBBBBB']))
    receipt = {'schema': integrity.SCHEMA, 'passed': True, 'database_name': 'mssystem_negfull',
               'input_rows': 2, 'expected_raw14': 2, 'observed_raw14': 2, 'fingerprints_checked': 2,
               'missing': [], 'extra': [], 'fingerprint_manifest': integrity.file_binding(fp),
               'native_jars': [integrity.file_binding(jar)],
               **{name: integrity.file_binding(p) for name, p in [('input', inp), ('database', db), ('protocol', protocol)]}}
    qa.write_text(json.dumps(receipt))
    return {'qa': qa, 'input': inp, 'database': db, 'protocol': protocol, 'jar': jar,
            'fingerprints': fp, 'receipt': receipt}


def check(b):
    return integrity.require_import_receipt(b['qa'], b['input'], b['database'], b['protocol'], database_name='mssystem_negfull')


def test_complete_bound_receipt_passes(bundle):
    assert check(bundle)['passed'] is True


@pytest.mark.parametrize('passed', [False, 'true', 1, None])
def test_existence_or_truthy_values_do_not_pass(bundle, passed):
    bundle['receipt']['passed'] = passed
    bundle['qa'].write_text(json.dumps(bundle['receipt']))
    with pytest.raises(RuntimeError, match='did not pass'):
        check(bundle)


@pytest.mark.parametrize('field', ['input', 'database', 'protocol', 'jar', 'fingerprints'])
def test_changed_input_db_program_or_manifest_is_rejected(bundle, field):
    with bundle[field].open('ab') as stream:
        stream.write(b'changed')
    with pytest.raises(RuntimeError, match='Content changed'):
        check(bundle)


def test_same_bytes_at_different_database_path_rejected(bundle, tmp_path):
    other = tmp_path / 'other.siriusdb'
    other.write_bytes(bundle['database'].read_bytes())
    bundle['database'] = other
    with pytest.raises(RuntimeError, match='wrong path'):
        check(bundle)


def test_old_identity_only_receipt_rejected(bundle):
    bundle['qa'].write_text(json.dumps({'passed': True, 'expected_raw14': 2, 'observed_raw14': 2}))
    with pytest.raises(RuntimeError, match='lacks native fingerprint'):
        check(bundle)


def test_manifest_cannot_substitute_a_duplicate_for_missing_structure(bundle):
    row = json.loads(bundle['fingerprints'].read_text().splitlines()[0])
    bundle['fingerprints'].write_text((json.dumps(row) + '\n') * 2)
    bundle['receipt']['fingerprint_manifest'] = integrity.file_binding(bundle['fingerprints'])
    bundle['qa'].write_text(json.dumps(bundle['receipt']))
    with pytest.raises(RuntimeError, match='Incomplete or duplicate'):
        check(bundle)


def test_all_zero_fingerprint_is_readable_not_missing():
    assert len(integrity.validate_sparse_bits([], {0, 2, 5})) == 64


@pytest.mark.parametrize('bits', [[-1], [6], [2, 2], [5, 2], [1.5], [True]])
def test_invalid_native_bit_layout_rejected(bits):
    with pytest.raises(RuntimeError):
        integrity.validate_sparse_bits(bits, {0, 2, 5})


def test_native_mask_indices_are_not_assumed_contiguous():
    assert integrity.validate_sparse_bits([0, 5], {0, 2, 5})


def test_isolated_native_target_stays_in_its_experiment_store(tmp_path):
    out, store = tmp_path / 'receipts', tmp_path / 'store'
    out.mkdir();store.mkdir()
    record = {'database_name': 'mssystem_negfull_integrity', 'path': str(store / 'new.siriusdb')}
    (out / 'DATABASE_TARGET.json').write_text(json.dumps(record))
    assert integrity.database_target(out, store) == record
    record['path'] = str(tmp_path / 'unrelated.siriusdb')
    (out / 'DATABASE_TARGET.json').write_text(json.dumps(record))
    with pytest.raises(RuntimeError, match='escaped'):
        integrity.database_target(out, store)


def test_live_stage_lock_rejects_duplicate_and_releases_on_close(tmp_path):
    spec = importlib.util.spec_from_file_location('tested_stage_lock', MODULE.with_name('bottleneck_execution_lock.py'))
    locking = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(locking)
    path = tmp_path / 'stage.lock'
    first = locking.acquire_execution_lock(path)
    try:
        original = path.read_bytes()
        with pytest.raises(RuntimeError, match='already running'):
            locking.acquire_execution_lock(path)
        assert path.read_bytes() == original
    finally:
        first.close()
    locking.acquire_execution_lock(path).close()
