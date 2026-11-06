import hashlib
import json
import pytest
from scripts.frozen_artifact_paths import (
    resolve_frozen_v15_artifact,
    resolve_identity_scope_migrated_v15_artifact,
)


def test_relocated_copy_requires_exact_bytes(tmp_path):
    root = tmp_path / 'mssystem_V16'
    p = root / 'output/tasks/frozen.json'
    p.parent.mkdir(parents=True)
    p.write_bytes(b'original frozen artifact')
    old = tmp_path / 'mssystem_V15/output/tasks/frozen.json'
    digest = hashlib.sha256(p.read_bytes()).hexdigest()
    assert resolve_frozen_v15_artifact(old, digest, root) == p
    p.write_bytes(b'new analysis with the same filename')
    with pytest.raises(ValueError, match='SHA-256 mismatch'):
        resolve_frozen_v15_artifact(old, digest, root)


def test_existing_mismatched_original_is_not_hidden_by_good_copy(tmp_path):
    root = tmp_path / 'mssystem_V16'
    root.mkdir()
    original = tmp_path / 'mssystem_V15/file.json'
    original.parent.mkdir()
    original.write_bytes(b'changed')
    (root/'file.json').write_bytes(b'frozen')
    with pytest.raises(ValueError, match='SHA-256 mismatch'):
        resolve_frozen_v15_artifact(original, hashlib.sha256(b'frozen').hexdigest(), root)


@pytest.mark.parametrize('recorded', ['elsewhere/file.json', 'mssystem_V15/../outside.json'])
def test_unrelated_or_escaping_path_is_not_relocated(tmp_path, recorded):
    root = tmp_path / 'mssystem_V16'
    root.mkdir()
    with pytest.raises(FileNotFoundError):
        resolve_frozen_v15_artifact(tmp_path/recorded, '0'*64, root)


def test_original_file_also_requires_hash(tmp_path):
    p=tmp_path/'frozen.json';p.write_bytes(b'original')
    assert resolve_frozen_v15_artifact(p, hashlib.sha256(b'original').hexdigest(), tmp_path)==p
    with pytest.raises(ValueError, match='requires'):
        resolve_frozen_v15_artifact(p, None, tmp_path)


def _migration_receipt(tmp_path, *, before, after, ranking_after=None):
    ranking_before = hashlib.sha256(b'ranking').hexdigest()
    identity = hashlib.sha256(b'identity').hexdigest()
    payload = {
        'schema_version': 'mssystem.v15.production-identity-scope-migration/1.0.0',
        'library_retrieval_repeated': False,
        'candidate_ranking_changed': False,
        'candidate_scores_changed': False,
        'database_scientific_identity_fields_changed': False,
        'migration_failure_count': 0,
        'tasks': [{
            'task_id': 'task-1',
            'status': 'completed',
            'before_flagship_sha256': before,
            'after_flagship_sha256': after,
            'ranking_unchanged': True,
            'before_ranking_sha256': ranking_before,
            'after_ranking_sha256': ranking_after or ranking_before,
            'database_scientific_identity_fields_unchanged': True,
            'before_scientific_identity_sha256': identity,
            'after_scientific_identity_sha256': identity,
        }],
    }
    path = tmp_path / 'receipt.json'
    path.write_text(json.dumps(payload), encoding='utf-8')
    return path, hashlib.sha256(path.read_bytes()).hexdigest()


def test_bound_identity_scope_migration_accepts_unchanged_ranking(tmp_path):
    artifact = tmp_path / 'frozen.json'
    artifact.write_bytes(b'after identity-scope migration')
    before = hashlib.sha256(b'before identity-scope migration').hexdigest()
    after = hashlib.sha256(artifact.read_bytes()).hexdigest()
    receipt, receipt_sha = _migration_receipt(tmp_path, before=before, after=after)

    assert resolve_identity_scope_migrated_v15_artifact(
        artifact,
        before,
        tmp_path,
        task_id='task-1',
        migration_receipt_path=receipt,
        migration_receipt_sha256=receipt_sha,
    ) == artifact


def test_bound_identity_scope_migration_rejects_changed_ranking(tmp_path):
    artifact = tmp_path / 'frozen.json'
    artifact.write_bytes(b'after identity-scope migration')
    before = hashlib.sha256(b'before identity-scope migration').hexdigest()
    after = hashlib.sha256(artifact.read_bytes()).hexdigest()
    receipt, receipt_sha = _migration_receipt(
        tmp_path,
        before=before,
        after=after,
        ranking_after=hashlib.sha256(b'changed ranking').hexdigest(),
    )

    with pytest.raises(ValueError, match='invariant failed'):
        resolve_identity_scope_migrated_v15_artifact(
            artifact,
            before,
            tmp_path,
            task_id='task-1',
            migration_receipt_path=receipt,
            migration_receipt_sha256=receipt_sha,
        )
