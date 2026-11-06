"""Integration contracts using synthetic candidates, never performance labels."""
import copy
import importlib
import json
from pathlib import Path

import pytest

from pipeline import spectral_tie_ranker as ranker
from pipeline import spectral_tie_scope as scope_guard
from pipeline import v15_flagship as pipeline
from pipeline.entity_confidence import single_energy_high_evidence
from pipeline.compound_conclusions import _project_primary_identity
from matching.v15_identity_resolution import apply_v15_identity_resolution


def candidate(letter):
    return {'compound_name': 'synthetic-contract-only', 'inchikey': letter * 14 + '-BBBBBBBBBB-C',
            'smiles': '', 'formula': 'C6H12O6', 'cosine': .9, 'entropy': .9,
            'reverse_match': .9, 'mass_score': .9, 'fragment_score': .9,
            'intensity_score': .9, 'final_score': .9, 'rank_score': .9,
            'has_structure': True, 'lib_spectrum_id': 'synthetic-' + letter}


def row(letters, index=0):
    candidates = [candidate(letter) for letter in letters]
    return {'spectrum_id': str(index), 'precursor_feature_id': 'fixed-synthetic-feature',
            'retention_time': 1., 'retention_time_unit': 'min', 'ion_mode': 'POS',
            'collision_energy': 10., 'collision_energy_unit': 'eV',
            'raw_peaks': [[10., 1.], [20., 2.]], 'tic': 3.,
            'match': {**copy.deepcopy(candidates[0]), 'top_candidates': candidates,
                      'confidence': 'High', 'two_dimensional_connectivity_confidence': 'High',
                      'identification_scope': 'exact_2d_connectivity'}}


@pytest.fixture
def active_policy(monkeypatch):
    loader = importlib.import_module('matching.v15_flagship_ranker')
    original = loader.load_v15_flagship_policy()
    policy = copy.deepcopy(original)
    policy['algorithm_id'] = 'mssystem-v16-shared-spectral-ties-multilibrary-1.0.0'
    policy['spectral_tie_handling'] = {
        'enabled': True, 'algorithm_id': ranker.ALGORITHM_ID,
        'scan_score_decimal_places': 8, 'ce_support_decimal_places': 6,
        'rank_support': 'mean_over_tied_ordinal_positions',
        'top1_support': 'equal_share_among_tied_leaders',
        'identity_key_affects_evidence_support': False, 'library_records_changed': False,
    }
    monkeypatch.setattr(loader, 'load_v15_flagship_policy', lambda: copy.deepcopy(policy))
    monkeypatch.setattr(pipeline, 'load_v15_flagship_policy', lambda: copy.deepcopy(policy))
    return policy


@pytest.mark.parametrize('letters,limit,incomplete', [('ABCDEF', 20, True), ('AB', 1, True), ('AB', 20, False)])
def test_complete_pipeline_high_scope(active_policy, tmp_path, letters, limit, incomplete):
    active_policy['candidate_pool']['maximum_candidates_per_scan'] = limit
    rows = [row(letters)]
    output = pipeline.v15_flagship_artifact_path(tmp_path, 'synthetic')
    pipeline.write_v15_flagship_artifact('synthetic', rows, output)
    payload = json.loads(output.read_text())
    entity = payload['entities'][0]
    high = single_energy_high_evidence(entity, rows)
    projection = _project_primary_identity(rows[0], entity, high)['match']
    assert pipeline.artifact_is_current(output)
    if incomplete:
        assert entity['decision']['scope'] == 'unresolved_precursor_event'
        assert entity['identity_resolution']['reportable_scope'] == 'unresolved_precursor_event'
        assert high is None
        assert projection['confidence'] == projection['two_dimensional_connectivity_confidence'] == projection['isomer_group_confidence'] == 'Unknown'
        assert entity['decision']['reportable_identity_key'] is None
        assert entity['decision']['ranking_lead_is_reportable_single_connectivity'] is False
        assert entity['decision_assurance']['decision_scope'] == 'unresolved_precursor_event'
    else:
        assert entity['decision']['scope'] == 'molecular_formula_isomer_group'
        assert high['scope'] == 'isomer_group'
        assert projection['isomer_group_confidence'] == 'High'


def test_missingness_changes_only_replicate_rule(active_policy, tmp_path):
    output = pipeline.v15_flagship_artifact_path(tmp_path, 'repeated')
    pipeline.write_v15_flagship_artifact('repeated', [row('AB', 0), row('AB', 1)], output)
    payload = json.loads(output.read_text())
    assert payload['metadata']['replicate_missingness_tie_algorithm_id'] == ranker.ALGORITHM_ID
    assert payload['metadata']['replicate_missingness_affects_primary_ranking'] is False
    assert payload['entities'][0]['candidates'][0]['entity_score'] == payload['entities'][0]['candidates'][1]['entity_score']


@pytest.mark.parametrize('field', ['spectral_tie_implementation_sha256', 'spectral_tie_base_source_sha256',
                                  'spectral_tie_handling_algorithm_id', 'spectral_tie_scope_guard_source_sha256'])
def test_wrong_implementation_cannot_be_reused_or_identity_upgraded(active_policy, tmp_path, field):
    output = pipeline.v15_flagship_artifact_path(tmp_path, 'synthetic')
    pipeline.write_v15_flagship_artifact('synthetic', [row('AB')], output)
    payload = json.loads(output.read_text())
    payload['metadata'][field] = 'wrong-implementation-binding'
    output.write_text(json.dumps(payload))
    assert not pipeline.artifact_is_current(output)
    with pytest.raises(ValueError, match='implementation changed'):
        pipeline.upgrade_v15_flagship_identity_scope_artifact(output, task_id='synthetic')


def test_current_identity_upgrade_is_noop(active_policy, tmp_path):
    output = pipeline.v15_flagship_artifact_path(tmp_path, 'synthetic')
    pipeline.write_v15_flagship_artifact('synthetic', [row('ABCDEF')], output)
    before = output.read_bytes()
    receipt = pipeline.upgrade_v15_flagship_identity_scope_artifact(output, task_id='synthetic')
    assert receipt['upgrade_status'] == 'already_current'
    assert output.read_bytes() == before


def test_scope_guard_idempotent_and_preserves_ranking_and_upstream(active_policy):
    base = ranker.build_flagship_result('synthetic', [row('ABCDEF')], policy=active_policy)
    resolved = apply_v15_identity_resolution(base)
    snapshot = scope_guard._ranking_snapshot(resolved)
    guarded = scope_guard.apply_spectral_tie_scope_guard(resolved)
    assert scope_guard.apply_spectral_tie_scope_guard(guarded) == guarded
    assert scope_guard._ranking_snapshot(guarded) == snapshot
    assert guarded['entities'][0]['identity_resolution']['upstream_decision_status'] == resolved['entities'][0]['identity_resolution']['upstream_decision_status']


def test_loaded_code_must_match_source_even_when_disk_hash_is_expected(active_policy, monkeypatch):
    base = importlib.import_module('matching.flagship_ranker')
    monkeypatch.setattr(base, '_candidate_rows', lambda row, policy: [])
    with pytest.raises((ValueError, KeyError)):
        ranker.build_flagship_result('synthetic', [row('AB')], policy=active_policy)


def test_source_change_cannot_be_labelled_as_loaded_implementation(active_policy, monkeypatch):
    monkeypatch.setattr(ranker, 'LOADED_IMPLEMENTATION_SHA256', 'old-loaded-source')
    with pytest.raises(ValueError, match='source changed'):
        ranker.build_flagship_result('synthetic', [row('AB')], policy=active_policy)
