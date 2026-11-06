import json
from pathlib import Path

import pytest

from runtime import verified_sql_artifacts as binding


def write_artifact(path, recall, ce):
    path.write_text(json.dumps({'metadata': {
        'candidate_recall_support_library_snapshot_id': recall,
        'source_paired_ce_gate_library_snapshot_id': ce,
    }, 'entities': []}))


@pytest.mark.parametrize('recall,ce,expected', [
    ('current', 'current', True), ('old', 'current', False),
    ('current', 'old', False), (None, None, False),
])
def test_both_sql_evidence_generations_must_match(tmp_path, monkeypatch, recall, ce, expected):
    path = tmp_path / 'entity.json'
    write_artifact(path, recall, ce)
    monkeypatch.setattr(binding, 'expected_sql_snapshot', lambda root: 'current')
    before = path.read_bytes()
    assert binding.sql_artifact_is_current(path, project_root=tmp_path) is expected
    assert path.read_bytes() == before


def test_fixture_without_production_configuration_is_unchanged(tmp_path):
    assert binding.sql_artifact_is_current(tmp_path / 'unused', project_root=tmp_path)


def test_old_sql_evidence_cannot_take_no_search_upgrade(tmp_path, monkeypatch):
    import web.server as server
    source = tmp_path / 'task_results.json'
    source.write_text('[]')
    artifact = tmp_path / 'task_v16_flagship.json'
    write_artifact(artifact, 'old', 'old')
    baseline_bytes = artifact.read_bytes()
    monkeypatch.setattr(server, 'sql_artifact_is_current', lambda *a, **k: False)
    monkeypatch.setattr(server, '_upgrade_flagship_artifact',
                        lambda *a, **k: pytest.fail('Old SQL evidence was upgraded without retrieval'))
    monkeypatch.setattr(server, '_v11_flagship_artifact_path', lambda *a: source)
    retriever = object()
    monkeypatch.setattr(server, '_v15_ce_evidence_retriever', lambda: retriever)
    calls = []
    def rebuild(*args, **kwargs):
        calls.append(kwargs)
        return {'status': 'rebuilt'}
    monkeypatch.setattr(server, '_write_flagship_artifact', rebuild)
    assert server._materialize_current_v15_flagship('task', {}, source, artifact) == {'status': 'rebuilt'}
    assert calls[0]['ce_evidence_retriever'] is retriever
    assert calls[0]['candidate_recall_retriever'] is retriever
    assert artifact.read_bytes() == baseline_bytes


def test_server_factory_uses_verified_adapter(tmp_path, monkeypatch):
    import web.server as server
    (tmp_path / 'config').mkdir()
    db = tmp_path / 'library.db'
    db.touch()
    (tmp_path / 'config/runtime_manifest.json').write_text(json.dumps({
        'library': {'production_full': {'snapshot_id': 'bound-db', 'path': str(db)}}
    }))
    monkeypatch.setattr(server, 'PROJECT_ROOT', tmp_path)
    monkeypatch.setattr(server, 'DB_PATH', db)
    sentinel = object()
    captured = []
    def factory(*args, **kwargs):
        captured.append((args, kwargs))
        return sentinel
    monkeypatch.setattr(server, 'create_verified_sql_retriever', factory)
    assert server._v15_ce_evidence_retriever() is sentinel
    assert captured[0][1]['library_snapshot_id'] == 'bound-db'
    assert captured[0][1]['project_root'] == tmp_path


def test_main_repairs_cannot_silently_lose_sql_configuration(tmp_path):
    (tmp_path / 'config').mkdir()
    (tmp_path / 'config/library_metadata_repairs.json').write_text('{}')
    with pytest.raises(ValueError, match='SQL metadata configuration is missing'):
        binding.expected_sql_snapshot(tmp_path)


def test_matching_declared_generation_cannot_bypass_failed_dependency_validation(tmp_path, monkeypatch):
    path = tmp_path / 'entity.json'
    write_artifact(path, 'current', 'current')
    monkeypatch.setattr(binding, 'expected_sql_snapshot', lambda root: 'current')
    monkeypatch.setattr(binding, 'sql_metadata_repair_preflight',
                        lambda root: {'compatible': False, 'error': 'changed source'})
    with pytest.raises(ValueError, match='changed source'):
        binding.sql_artifact_is_current(path, project_root=tmp_path)
