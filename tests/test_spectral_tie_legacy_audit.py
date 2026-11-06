"""Production regression: original audit must not inherit recovered tie coverage."""
import copy

import matching.flagship_ranker as base
from matching.v15_flagship_ranker import load_v15_flagship_policy
from pipeline import spectral_tie_ranker as tie


def test_hidden_tie_majority_does_not_rewrite_original_legacy_coverage():
    policy = copy.deepcopy(load_v15_flagship_policy())
    policy['candidate_pool'].update(maximum_candidates_per_scan=1, shortlist_size=5,
                                    same_energy_replicate_support='legacy_zero_fill')
    policy['identity_registry'].update(enabled=False, required=False)
    policy['spectral_tie_handling'] = {
        'enabled': True, 'algorithm_id': tie.ALGORITHM_ID,
        'scan_score_decimal_places': 8, 'ce_support_decimal_places': 6,
        'rank_support': 'mean_over_tied_ordinal_positions',
        'top1_support': 'equal_share_among_tied_leaders',
        'identity_key_affects_evidence_support': False, 'library_records_changed': False,
    }

    def candidate(letter):
        return dict(compound_name='synthetic-' + letter,
                    inchikey=letter * 14 + '-BBBBBBBBBB-C', smiles='', formula='',
                    cosine=.8, entropy=.8, reverse_match=.8, fragment_score=.8,
                    intensity_score=.8, mass_ppm=0., final_score=.8, rank_score=.8,
                    dataset_source='synthetic-test')

    a, b = candidate('A'), candidate('B')
    rows = []
    for i, candidates in enumerate(([b], [a, b], [a, b])):
        match = copy.deepcopy(candidates[0])
        match['top_candidates'] = copy.deepcopy(candidates)
        rows.append(dict(spectrum_id=str(i), precursor_feature_id='one-entity',
                         precursor_mz=300., precursor_charge=1, retention_time=1.,
                         collision_energy=20., collision_energy_unit='eV', ion_mode='POS',
                         raw_peaks=[[100., 100.], [150., 50.]],
                         query_peaks=[[100., 100.], [150., 50.]], match=match))
    before_rows, before_policy = copy.deepcopy(rows), copy.deepcopy(policy)
    functions = base._candidate_rows, base._finalize_entity, base.build_flagship_result
    original = base.build_flagship_result('synthetic', rows, policy=policy)['entities'][0]
    fixed = tie.build_flagship_result('synthetic', rows, policy=policy)['entities'][0]
    bm = {c['identity_key']: c for c in original['candidates']}
    fm = {c['identity_key']: c for c in fixed['candidates']}
    identity = 'ik14:BBBBBBBBBBBBBB'
    assert bm[identity]['legacy_v14_entity_score'] == 0.
    old_coverage = bm[identity]['score_components']['collision_energy_coverage']
    new_coverage = fm[identity]['score_components']['collision_energy_coverage']
    assert (old_coverage, new_coverage) == (0., 1.)
    # The former bug reused the repaired coverage in the historical formula.
    assert bm[identity]['legacy_v14_entity_score'] + .15 * (new_coverage - old_coverage) == .15
    for key in bm:
        assert fm[key]['legacy_v14_entity_score'] == bm[key]['legacy_v14_entity_score']
        assert fm[key]['legacy_v14_entity_rank'] == bm[key]['legacy_v14_entity_rank']
    assert fixed['v14_audit_baseline'] == original['v14_audit_baseline']
    assert len({c['entity_score'] for c in fixed['candidates']}) == 1
    assert rows == before_rows and policy == before_policy
    assert functions == (base._candidate_rows, base._finalize_entity, base.build_flagship_result)
