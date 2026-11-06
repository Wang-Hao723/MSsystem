import copy, importlib
import pytest
S = importlib.import_module('matching.search')
from matching.source_ms_stage import source_ms_stage


@pytest.mark.parametrize('text,status,anomaly,explicit,ineligible', [
    ('MS2','explicit_ms2_only',None,True,False),
    ('MS2-MS2 Composite','explicit_ms2_only',None,True,False),
    ('MS2-MS3 Composite','explicit_non_ms2','mixed_ms2_msn_composite',False,True),
    ('MS2-MS4 Composite','explicit_non_ms2','mixed_ms2_msn_composite',False,True),
    ('MS2-MS5 Composite','explicit_non_ms2','mixed_ms2_msn_composite',False,True),
    ('MS3','explicit_non_ms2','higher_order_spectrum',False,True),
    ('MS4','explicit_non_ms2','higher_order_spectrum',False,True),
    ('MS1','explicit_non_ms2','source_stage_below_ms2',False,True),
    ('MS3-MS5 Composite','explicit_non_ms2','higher_order_composite',False,True),
    ('', 'missing', None, False, False),
    (None, 'missing', None, False, False),
    ('MS/MS','unrecognized','unrecognized_source_stage_text',False,False),
    ('2','unrecognized','unrecognized_source_stage_text',False,False),
    ('MS4-MS2 Composite','unrecognized','invalid_stage_range',False,False),
])
def test_explicit_source_stage(text,status,anomaly,explicit,ineligible):
    r=source_ms_stage(text)
    assert (r['status'],r['source_anomaly_type'],r['explicit_ms2_only'],r['ineligible_for_ordinary_ms2_contrast']) == (status,anomaly,explicit,ineligible)


def rec(stage):
    return {'adduct':'[M+H]+','query_adduct':'[M+H]+','query_adduct_reliable':True,
            'lib_collision_energy':30,'lib_collision_energy_unit':'eV',
            'lib_instrument':'Thermo Q Exactive HF','lib_source_ms_level':stage}


def test_unknown_level_cannot_complete_conditions_but_is_not_guessed_or_removed():
    a,b=rec(None),rec('MS2')
    check=S._candidate_pair_comparability(a,b)
    assert check['comparable'] and not S._condition_metadata_complete(check)
    q=S._query_candidate_comparability(a,{'value':30,'unit':'eV'},query_instrument='Thermo Q Exactive HF')
    assert q['source_ms_stage']['status']=='missing' and not S._query_metadata_complete(q)


def test_explicit_ms2_conditions_can_still_be_complete():
    a,b=rec('MS2'),rec('MS2-MS2 Composite')
    assert S._condition_metadata_complete(S._candidate_pair_comparability(a,b))
    q=S._query_candidate_comparability(b,{'value':30,'unit':'eV'},query_instrument='Thermo Q Exactive HF')
    assert S._query_metadata_complete(q)


@pytest.mark.parametrize('other',['MS2',None,'MS3','MS2-MS5 Composite'])
def test_explicit_high_order_is_not_a_ms2_contrast_even_against_same_high_order(other):
    a,b=rec('MS2-MS3 Composite'),rec(other)
    before=copy.deepcopy((a,b));check=S._candidate_pair_comparability(a,b)
    assert not check['comparable'] and not S._condition_metadata_complete(check)
    assert (a,b)==before
    q=S._query_candidate_comparability(a,{'value':30,'unit':'eV'},query_instrument='Thermo Q Exactive HF')
    assert q['source_ms_stage']['ineligible_for_ordinary_ms2_contrast'] and not S._query_metadata_complete(q)


def test_instrument_and_peaks_cannot_prove_stage():
    for model in ['Thermo Finnigan LTQ','Thermo Q Exactive HF','SCIEX ZenoTOF 7600']:
        a=rec(None);a['lib_instrument']=model;a['lib_peaks']=[[100,100],[118,20]]
        check=S._candidate_pair_comparability(a,a)
        assert check['comparable']
        assert check['source_ms_stage_comparability']['status']=='source_ms_level_incomplete'


def test_existing_non_stage_conditions_remain_enforced():
    a,b=rec('MS2'),rec('MS2');b['lib_instrument']='Agilent 6530 Q-TOF'
    assert not S._candidate_pair_comparability(a,b)['comparable']
    b=rec('MS2');b['lib_collision_energy_unit']='NCE'
    assert not S._candidate_pair_comparability(a,b)['comparable']

