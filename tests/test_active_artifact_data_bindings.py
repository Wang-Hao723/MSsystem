"""A changed library cannot overwrite an earlier version's report path."""
import copy
import hashlib

import pytest

from matching.v15_flagship_ranker import load_v15_flagship_policy
from pipeline import active_artifact_binding as binding
from pipeline.derived_artifact_paths import derived_artifact_path


def test_new_library_changes_destination_and_stale_loaded_binding_is_rejected(tmp_path, monkeypatch):
    data = tmp_path/'library.json'
    data.write_text('{"accepted_cache":"old"}')
    key = 'data:test-library'
    monkeypatch.setattr(binding, '_DATA_PATHS', {key: data})
    monkeypatch.setattr(binding, '_LOADED_DATA_DIGESTS', {key: hashlib.sha256(data.read_bytes()).hexdigest()})
    policy = copy.deepcopy(load_v15_flagship_policy())
    original = binding.active_implementation_sha256(policy)
    policy['artifact_namespace']['implementation_sha256'] = original
    old_path = derived_artifact_path(tmp_path, 'sample', '_v16_flagship.json', policy=policy,
                                     actual_implementation_sha256=original)
    old_path.write_bytes(b'original report')
    data.write_text('{"accepted_cache":"new"}')
    with pytest.raises(ValueError, match='library binding changed'):
        binding.active_implementation_sha256(policy)
    # Model a reviewed process restart: load the newly verified data identity.
    monkeypatch.setattr(binding, '_LOADED_DATA_DIGESTS', {key: hashlib.sha256(data.read_bytes()).hexdigest()})
    revised = binding.active_implementation_sha256(policy)
    assert revised != original
    policy['artifact_namespace']['implementation_sha256'] = revised
    new_path = derived_artifact_path(tmp_path, 'sample', '_v16_flagship.json', policy=policy,
                                     actual_implementation_sha256=revised)
    assert new_path != old_path
    new_path.write_bytes(b'new report')
    assert old_path.read_bytes() == b'original report'


def test_retrieval_and_mass_implementation_are_report_dependencies():
    assert {'matching.search', 'matching.reference_mass'} <= set(binding._SOURCE_PATHS)
    assert {'data:config/library_record_store.json', 'data:config/sql_metadata_repairs.json'} <= set(binding._DATA_PATHS)


def test_source_stage_and_sql_projection_are_report_dependencies():
    assert {'matching.source_ms_stage', 'matching.evirank_retrieval'} <= set(binding._SOURCE_PATHS)


def test_source_stage_change_cannot_overwrite_report(tmp_path, monkeypatch):
    source = tmp_path / 'source_ms_stage.py'
    source.write_text('stage_rule = 1\n')
    key = 'matching.source_ms_stage'
    monkeypatch.setattr(binding, '_SOURCE_PATHS', {key: source})
    monkeypatch.setattr(binding, '_LOADED_SOURCE_DIGESTS', {key: hashlib.sha256(source.read_bytes()).hexdigest()})
    policy = copy.deepcopy(load_v15_flagship_policy())
    original = binding.active_implementation_sha256(policy)
    policy['artifact_namespace']['implementation_sha256'] = original
    old_path = derived_artifact_path(tmp_path, 'source-stage-test', '_v16_flagship.json', policy=policy,
                                    actual_implementation_sha256=original)
    old_path.write_bytes(b'original source-stage report')
    source.write_text('stage_rule = 2\n')
    with pytest.raises(ValueError, match='implementation changed'):
        binding.active_implementation_sha256(policy)
    monkeypatch.setattr(binding, '_LOADED_SOURCE_DIGESTS', {key: hashlib.sha256(source.read_bytes()).hexdigest()})
    revised = binding.active_implementation_sha256(policy)
    policy['artifact_namespace']['implementation_sha256'] = revised
    new_path = derived_artifact_path(tmp_path, 'source-stage-test', '_v16_flagship.json', policy=policy,
                                    actual_implementation_sha256=revised)
    assert new_path != old_path and old_path.read_bytes() == b'original source-stage report'


def test_full_source_verifier_is_report_dependency():
    assert 'runtime.source_file_verification' in binding._SOURCE_PATHS


def test_source_verifier_same_size_rewrite_requires_new_report(tmp_path, monkeypatch):
    import os
    source = tmp_path / 'source_file_verification.py'
    source.write_text('verify_mode = 1\n')
    key = 'runtime.source_file_verification'
    original_stat = source.stat()
    monkeypatch.setattr(binding, '_SOURCE_PATHS', {key: source})
    monkeypatch.setattr(binding, '_LOADED_SOURCE_DIGESTS', {key: hashlib.sha256(source.read_bytes()).hexdigest()})
    policy = copy.deepcopy(load_v15_flagship_policy())
    old = binding.active_implementation_sha256(policy)
    source.write_text('verify_mode = 2\n')
    os.utime(source, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
    with pytest.raises(ValueError, match='implementation changed'):
        binding.active_implementation_sha256(policy)
    monkeypatch.setattr(binding, '_LOADED_SOURCE_DIGESTS', {key: hashlib.sha256(source.read_bytes()).hexdigest()})
    assert binding.active_implementation_sha256(policy) != old
