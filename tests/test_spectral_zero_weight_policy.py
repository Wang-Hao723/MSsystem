"""Explicit zero disables evidence; default scores retain the verified baseline.

Frozen default digests were produced by the independently imported pre-repair
1.0.1 file in research/20260914_system_scientific_audit/before/pipeline/.
They contain every numeric entity leaf plus identity ordering, excluding only
artifact metadata. The repository test does not require that external file.
"""
import copy
import hashlib
import json

import pytest

from matching.v15_flagship_ranker import load_v15_flagship_policy
from pipeline import spectral_tie_ranker as ranker


ENTITY_KEYS = ('mean_collision_energy_support', 'median_collision_energy_support',
               'collision_energy_coverage', 'collision_energy_top1_fraction')


def _candidate(letter, score):
    return dict(compound_name='synthetic-' + letter,
                inchikey=letter * 14 + '-BBBBBBBBBB-C', smiles='', formula='C6H12O6',
                cosine=score, entropy=score, reverse_match=score, mass_score=score,
                fragment_score=score, intensity_score=score, final_score=score,
                rank_score=score, has_structure=True, lib_spectrum_id='synthetic-' + letter)


def _row(index, energy, scores):
    candidates = [_candidate(letter, score) for letter, score in zip('ABC', scores)]
    return dict(spectrum_id=str(index), precursor_feature_id='fixed-synthetic-feature',
                retention_time=1., retention_time_unit='min', ion_mode='POS',
                collision_energy=energy, collision_energy_unit='eV',
                raw_peaks=[[10., 1.], [20., 2.]], tic=3.,
                match={**copy.deepcopy(candidates[0]), 'top_candidates': candidates})


def _policy(module=ranker):
    policy = copy.deepcopy(load_v15_flagship_policy())
    policy['spectral_tie_handling']['algorithm_id'] = module.ALGORITHM_ID
    return policy


def _cases():
    return {
        'single_distinct': [_row(0, 10., [.90, .89, .60])],
        'single_tied': [_row(0, 10., [.90, .90, .60])],
        'multiple_ce': [_row(0, 10., [.90, .70, .60]), _row(1, 20., [.65, .88, .60])],
        'replicates_and_ce': [_row(0, 10., [.90, .70]), _row(1, 10., [.80, .75]),
                              _row(2, 20., [.85, .85])],
    }


def _numeric_entity_projection(payload):
    values = {}
    def visit(value, path):
        if isinstance(value, dict):
            for key, child in value.items():
                visit(child, (*path, str(key)))
        elif isinstance(value, list):
            for index, child in enumerate(value):
                visit(child, (*path, str(index)))
        elif isinstance(value, (int, float)):
            values['/'.join(path)] = value
        elif path and path[-1] in {'identity_key', 'entity_id'}:
            values['/'.join(path)] = value
    visit(payload['entities'], ('entities',))
    return values


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


BASELINE_NUMERIC_SHA256 = {
    'single_distinct': '1f9c84dc5011ea95813eb234150c870b9938ca9f759fec0a92baef5a852eef96',
    'single_tied': 'fb91f224f05dfa7c56711c48cdc75d7298fe303e61354b81185d574fb0eedbe6',
    'multiple_ce': 'e5a636a6da4a3d24dd19716132a6290b2b66e327295c96347ccfdb31f12061ff',
    'replicates_and_ce': '4ec4ce28e324e2fc5aff4924988739d27a10c852fb420fc440c95b24fdef95f3',
}


@pytest.mark.parametrize('case', list(BASELINE_NUMERIC_SHA256))
def test_nonzero_default_scores_match_frozen_pre_repair_implementation(case):
    output = ranker.build_flagship_result('synthetic-zero-audit', _cases()[case], policy=_policy())
    projection = _numeric_entity_projection(output)
    assert len(projection) > 100
    assert _digest(projection) == BASELINE_NUMERIC_SHA256[case]


@pytest.mark.parametrize('spectral,reciprocal', [(1., 0.), (0., 1.)])
def test_explicit_zero_scan_weight_removes_its_contribution(spectral, reciprocal):
    policy = _policy()
    policy['conditioned_multilibrary']['enabled'] = False
    policy['scan_score']['spectral_fraction'] = spectral
    policy['scan_score']['reciprocal_rank_fraction'] = reciprocal
    candidates = ranker._ranker_namespace()['_candidate_rows'](_row(0, 10., [.90, .70]), policy)
    assert len(candidates) == 2
    for item in candidates:
        expected = round(spectral * item['spectral_score'] + reciprocal * item['reciprocal_rank_support'], 6)
        assert item['scan_support'] == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize('disabled', ENTITY_KEYS)
def test_each_entity_component_can_be_explicitly_disabled(disabled):
    policy = _policy()
    policy['entity_score'][disabled] = 0.
    output = ranker.build_flagship_result('synthetic-zero-audit', _cases()['multiple_ce'], policy=policy)
    for candidate in output['entities'][0]['candidates']:
        components = candidate['score_components']
        expected = sum(policy['entity_score'][key] * components[key] for key in ENTITY_KEYS)
        assert candidate['entity_score'] == pytest.approx(expected, abs=1e-6)


def test_all_zero_entity_weights_yield_zero_and_keep_original_legacy_audit():
    policy = _policy()
    policy['entity_score'] = {key: 0. for key in ENTITY_KEYS}
    output = ranker.build_flagship_result('synthetic-zero-audit', _cases()['single_distinct'], policy=policy)
    candidates = output['entities'][0]['candidates']
    assert candidates and all(candidate['entity_score'] == 0. for candidate in candidates)
    assert all(candidate['legacy_v14_entity_score'] > 0. for candidate in candidates)
    assert output['entities'][0]['winner'] is None
