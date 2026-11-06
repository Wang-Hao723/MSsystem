from copy import deepcopy

from pipeline.entity_confidence import single_energy_high_evidence
from matching.flagship_ranker import connectivity_identity
from pipeline.compound_conclusions import _project_primary_identity
from pipeline.compound_role_projection import (
    apply_report_compound_roles, annotate_match_compound_roles,
    spectral_contaminant_veto,
)


def fixture():
    match = {
        'compound_name': 'caffeine', 'smiles': 'Cn1c(=O)c2c(ncn2C)n(C)c1=O',
        'formula': 'C8H10N4O2', 'confidence': 'High',
        'two_dimensional_connectivity_confidence': 'High',
        'identification_scope': 'exact_2d_connectivity', 'cosine': .95,
        'lib_spectrum_id': 'test-measured-reference',
    }
    identity = connectivity_identity(match)[0]
    candidate = {'identity_key': identity, 'candidate': {k: match[k] for k in ['compound_name', 'smiles', 'formula']},
                 'entity_score': .91, 'entity_rank': 1, 'vetoed': False}
    entity = {'entity_id': 'test', 'distinct_collision_energy_bins': 1, 'member_result_indices': [0],
              'winner': candidate, 'candidates': [candidate],
              'decision': {'status': 'selected_candidate_only', 'scope': 'single_collision_energy_candidate',
                           'blockers': ['insufficient_distinct_collision_energies']},
              'identity_resolution': {'upstream_decision_status': 'selected_candidate_only',
                                      'competitive_candidate_identity_keys': [identity]}}
    return entity, [{'spectrum_id': 'scan=1', 'match': match}]


def test_same_identity_measured_high_survives_single_energy_projection():
    entity, rows = fixture()
    before = deepcopy((entity, rows))
    evidence = single_energy_high_evidence(entity, rows)
    projected = _project_primary_identity(rows[0], entity, evidence)['match']
    assert projected['two_dimensional_connectivity_confidence'] == 'High'
    assert projected['exact_structure_confidence'] == 'Unknown'
    assert projected['entity_high_evidence']['standard_confirmed'] is False
    assert projected['rank_score'] == entity['winner']['entity_score']
    assert (entity, rows) == before


def test_high_cannot_transfer_to_a_different_ranked_structure():
    entity, rows = fixture()
    entity['winner']['identity_key'] = 'std2d:AAAAAAAAAAAAAA'
    assert single_energy_high_evidence(entity, rows) is None


def test_cross_formula_conflict_and_abstention_stay_unresolved():
    entity, rows = fixture()
    entity['decision']['scope'] = 'unresolved_precursor_event'
    assert single_energy_high_evidence(entity, rows) is None
    entity['decision']['scope'] = 'none'
    assert single_energy_high_evidence(entity, rows) is None


def test_group_evidence_cannot_be_converted_to_one_structure():
    entity, rows = fixture()
    rows[0]['match'].update(identification_scope='isomer_group', isomer_group_confidence='High')
    assert single_energy_high_evidence(entity, rows) is None


def test_group_high_retains_all_competing_structures():
    entity, rows = fixture()
    other = deepcopy(entity['winner'])
    other['identity_key'] = 'ik14:AAAAAAAAAAAAAA'
    other['entity_rank'] = 2
    entity['candidates'].append(other)
    entity['decision'].update(scope='molecular_formula_isomer_group', status='selected_isomer_group')
    rows[0]['match'].update(identification_scope='isomer_group', isomer_group_confidence='High',
                            isomer_set=[{'inchikey': 'AAAAAAAAAAAAAA-BBBBBBBBBB-C', 'spectrally_competitive': True}])
    evidence = single_energy_high_evidence(entity, rows)
    assert evidence['scope'] == 'isomer_group'
    projected = _project_primary_identity(rows[0], entity, evidence)['match']
    assert projected['isomer_group_confidence'] == 'High'
    assert projected['two_dimensional_connectivity_confidence'] == 'Unknown'
    entity['candidates'].pop()
    assert single_energy_high_evidence(entity, rows) is None


def test_low_tier_contaminant_and_transformation_conflict_cannot_be_promoted():
    for change in ('low', 'contaminant', 'transformation'):
        entity, rows = fixture()
        if change == 'low':
            rows[0]['match'].update(confidence='Medium', two_dimensional_connectivity_confidence='Medium')
        elif change == 'contaminant':
            rows[0]['match']['potential_contaminant'] = True
        else:
            entity['winner']['chemistry_gate'] = {'precursor_transformation_ambiguity': True}
        assert single_energy_high_evidence(entity, rows) is None


def test_repeated_single_energy_scans_do_not_require_high_on_every_member():
    entity, rows = fixture()
    rows.insert(0, deepcopy(rows[0]))
    rows[0]['match'].update(confidence='Medium', two_dimensional_connectivity_confidence='Medium')
    entity['member_result_indices'] = [0, 1]
    evidence = single_energy_high_evidence(entity, rows)
    assert [s['result_index'] for s in evidence['supporting_spectra']] == [1]


def _report_role_fixture():
    entity, rows = fixture()
    match = rows[0]['match']
    match.update(compound_name='Laureth-3', smiles='CCCCCCCCCCCCOCCOCCOCCO',
                 inchikey='FKMHSNTVILORFA-UHFFFAOYSA-N', formula='C18H38O4',
                 potential_contaminant=False)
    identity = connectivity_identity(match)[0]
    entity['winner']['identity_key'] = identity
    entity['winner']['candidate'] = {k: match[k] for k in ('compound_name','smiles','inchikey','formula')}
    entity['identity_resolution']['competitive_candidate_identity_keys'] = [identity]
    return entity, rows


def test_new_report_role_preserves_high_and_repeated_projection_is_idempotent():
    entity, rows = _report_role_fixture()
    for fn in (annotate_match_compound_roles, apply_report_compound_roles,
               annotate_match_compound_roles, apply_report_compound_roles):
        fn(rows[0]['match'])
        assert rows[0]['match']['potential_contaminant'] is True
        assert spectral_contaminant_veto(rows[0]['match']) is False
        assert single_energy_high_evidence(entity, rows)['tier'] == 'High'
    assert rows[0]['match']['report_role_projection']['pre_report_contaminant_veto'] is False
    once = deepcopy(rows[0]['match'])
    annotate_match_compound_roles(rows[0]['match'])
    assert rows[0]['match'] == once
    assert 'source_compound_role' not in rows[0]['match']


def test_original_contaminant_veto_survives_new_role_projection_and_repetition():
    entity, rows = _report_role_fixture()
    rows[0]['match']['potential_contaminant'] = True
    for fn in (annotate_match_compound_roles, apply_report_compound_roles,
               annotate_match_compound_roles):
        fn(rows[0]['match'])
        assert spectral_contaminant_veto(rows[0]['match']) is True
        assert single_energy_high_evidence(entity, rows) is None


def test_original_veto_is_not_cleared_for_unlisted_or_conflicting_structures():
    for fn in (annotate_match_compound_roles, apply_report_compound_roles):
        for conflict in (False, True):
            entity, rows = fixture()
            rows[0]['match']['potential_contaminant'] = True
            if conflict:
                rows[0]['match']['inchikey'] = 'FKMHSNTVILORFA-UHFFFAOYSA-N'
            fn(rows[0]['match'])
            assert rows[0]['match']['potential_contaminant'] is True
            assert spectral_contaminant_veto(rows[0]['match']) is True
            assert single_energy_high_evidence(entity, rows) is None


def test_incomplete_legacy_or_identity_changed_marker_cannot_bypass_veto():
    for change in ('missing_pre_veto', 'old_version', 'wrong_key', 'changed_structure'):
        entity, rows = _report_role_fixture()
        annotate_match_compound_roles(rows[0]['match'])
        marker = rows[0]['match']['report_role_projection']
        if change == 'missing_pre_veto':
            marker.pop('pre_report_contaminant_veto')
        elif change == 'old_version':
            marker['catalog_version'] = 'old-unverified'
        elif change == 'wrong_key':
            marker['matched_inchikey'] = 'AAAAAAAAAAAAAA-UHFFFAOYSA-N'
        else:
            rows[0]['match']['smiles'] = 'CCO'
        assert spectral_contaminant_veto(rows[0]['match']) is True
        assert single_energy_high_evidence(entity, rows) is None
