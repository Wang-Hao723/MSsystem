"""Versioned writes must preserve old evidence through actual pipeline calls."""
import copy
import importlib
import json
import os
from pathlib import Path

import pytest

from pipeline import active_artifact_binding as binding
from pipeline import compound_conclusions as compounds
from pipeline import v15_flagship as flagship
from tests.test_spectral_tie_query_pipeline import row


@pytest.fixture
def versioned_policy(monkeypatch):
    loader = importlib.import_module('matching.v15_flagship_ranker')
    policy = copy.deepcopy(loader.load_v15_flagship_policy())
    policy['algorithm_id'] = 'mssystem-v16-shared-spectral-ties-multilibrary-1.0.0'
    policy['artifact_namespace'] = {'name': 'exact-spectral-ties-v1',
                                    'implementation_sha256': '0' * 64}
    policy['artifact_namespace']['implementation_sha256'] = binding.active_implementation_sha256(policy)
    for module in (loader, flagship, compounds):
        monkeypatch.setattr(module, 'load_v15_flagship_policy', lambda: copy.deepcopy(policy))
    return policy


def materialize(tmp_path):
    rows = [row('ABCDEF')]
    source = tmp_path / 'synthetic_results.json'
    source.write_text(json.dumps(rows))
    primary = flagship.v15_flagship_artifact_path(tmp_path, 'synthetic')
    output = compounds.compound_conclusion_artifact_path(tmp_path, 'synthetic')
    flagship.write_v15_flagship_artifact('synthetic', rows, primary)
    compounds.write_compound_conclusion_artifact('synthetic', rows, output,
        source_artifact=source, primary_identity_artifact=primary)
    return rows, source, primary, output


def test_actual_writers_preserve_legacy_and_do_not_promote_hidden_ties(versioned_policy, tmp_path):
    old = {tmp_path / ('synthetic' + suffix): b'original artifact bytes\n'
           for suffix in ('_v16_flagship.json', '_v11_compounds.json')}
    for path, data in old.items():
        path.write_bytes(data)
    rows, source, primary, output = materialize(tmp_path)
    assert flagship.artifact_is_current(primary)
    assert compounds.artifact_is_current(output, source, primary)
    assert all(path.read_bytes() == data for path, data in old.items())
    p, c = json.loads(primary.read_text()), json.loads(output.read_text())
    assert p['metadata']['derived_artifact_binding'] == c['metadata']['derived_artifact_binding']
    assert p['entities'][0]['decision']['scope'] == 'unresolved_precursor_event'
    assert c['metadata']['high_confidence_2d_connectivity_count'] == 0
    assert c['metadata']['high_confidence_isomer_group_count'] == 0
    assert c['metadata']['spectral_high_confidence_scope_group_count'] == 0
    assert c['metadata']['derived_primary_hash_basis'] == 'file_bytes'
    before = primary.read_bytes()
    assert flagship.upgrade_v15_flagship_identity_scope_artifact(primary)['upgrade_status'] == 'already_current'
    assert primary.read_bytes() == before


def test_wrong_destination_rejected_before_scoring(versioned_policy, monkeypatch, tmp_path):
    old = tmp_path / 'synthetic_v16_flagship.json'
    old.write_bytes(b'old')
    def forbidden(*args, **kwargs):
        raise AssertionError('Scoring must not start for an old destination')
    monkeypatch.setattr(flagship, 'build_v15_flagship_result', forbidden)
    with pytest.raises(ValueError, match='content-bound destination'):
        flagship.write_v15_flagship_artifact('synthetic', [row('AB')], old)
    assert old.read_bytes() == b'old'


def test_policy_change_creates_new_path_and_cannot_upgrade_old(versioned_policy, tmp_path):
    _, _, primary, output = materialize(tmp_path)
    before = primary.read_bytes(), output.read_bytes()
    versioned_policy['decision']['minimum_entity_score'] += .01
    assert flagship.v15_flagship_artifact_path(tmp_path, 'synthetic') != primary
    assert compounds.compound_conclusion_artifact_path(tmp_path, 'synthetic') != output
    assert not flagship.artifact_is_current(primary)
    with pytest.raises(ValueError, match='new bound materialization'):
        flagship.upgrade_v15_flagship_identity_scope_artifact(primary)
    assert (primary.read_bytes(), output.read_bytes()) == before


def test_primary_bytes_are_checked_even_with_preserved_stat(versioned_policy, tmp_path):
    _, source, primary, output = materialize(tmp_path)
    stat = primary.stat()
    original = primary.read_bytes()
    token = b'"entity_score":0.'
    at = original.index(token) + len(token)
    changed = original[:at] + (b'8' if original[at:at+1] != b'8' else b'7') + original[at+1:]
    primary.write_bytes(changed)
    os.utime(primary, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert primary.stat().st_size == stat.st_size
    assert primary.stat().st_mtime_ns == stat.st_mtime_ns
    assert not compounds.artifact_is_current(output, source, primary)


def test_wrong_declared_implementation_never_falls_back(versioned_policy, tmp_path):
    versioned_policy['artifact_namespace']['implementation_sha256'] = 'f' * 64
    with pytest.raises(ValueError, match='differs from loaded'):
        flagship.v15_flagship_artifact_path(tmp_path, 'synthetic')
    assert not list(tmp_path.iterdir())


def test_rank_free_policy_cannot_use_exact_tie_provider(versioned_policy):
    versioned_policy['algorithm_id'] = 'mssystem-rank-free-spectral-entity-research-1.0.0'
    with pytest.raises(ValueError, match='No reviewed'):
        binding.active_implementation_sha256(versioned_policy)


def test_changed_loaded_projection_dependency_requires_restart(versioned_policy, monkeypatch):
    monkeypatch.setitem(binding._LOADED_SOURCE_DIGESTS, 'pipeline.compound_conclusions', 'f' * 64)
    with pytest.raises(ValueError, match='restart and rebind'):
        binding.active_implementation_sha256(versioned_policy)
