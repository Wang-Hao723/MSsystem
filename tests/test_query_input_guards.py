import copy
import importlib

import pytest

from pipeline.query_input_guards import (annotate_precursor_features,
    resolve_isopair_adduct_context, apply_isopair_ranker)
from pipeline.entity_builder import build_entity_graph


def context(**kwargs):
    return resolve_isopair_adduct_context(formula='C6H12O6',query_precursor_mz=181.070664,
        **({'query_adduct':'[M+H]+','query_polarity':'POS','query_charge':1}|kwargs))


@pytest.mark.parametrize('kwargs,reason',[
    ({'query_adduct':'BOGUS'},'unsupported_explicit_query_adduct'),
    ({'query_charge':2},'explicit_adduct_charge_conflict'),
    ({'query_charge':-1},'query_polarity_charge_conflict'),
    ({'query_charge':float('nan')},'invalid_query_charge'),
    ({'query_charge':1.5},'invalid_query_charge'),
])
def test_invalid_explicit_annotation_cannot_qualify(kwargs,reason):
    result=context(**kwargs)
    assert result['reason']==reason
    assert not result['ranking_eligible'] and not result['high_precision_gate_eligible']


def test_valid_explicit_input_keeps_its_evidence():
    assert context()['ranking_eligible'] is True
    assert context()['independent_query_evidence'] is True


@pytest.mark.parametrize('charge', [1, 1.0, '1', -1, '-1', '1-'])
def test_negative_mode_accepts_legacy_magnitude_and_matching_signed_charge(charge):
    result = context(query_adduct='[M-H]-', query_polarity='NEG', query_charge=charge)
    assert result['reason'] == 'explicit_query_adduct'
    assert result['query_polarity'] == 'NEG'
    assert result['ranking_eligible'] is True


@pytest.mark.parametrize('charge', ['+1', '1+'])
def test_explicit_positive_charge_still_conflicts_with_negative_mode(charge):
    result = context(query_adduct='[M-H]-', query_polarity='NEG', query_charge=charge)
    assert result['reason'] == 'query_polarity_charge_conflict'
    assert result['ranking_eligible'] is False


def test_unsigned_magnitude_does_not_invent_query_polarity():
    result = context(query_adduct=None, query_polarity=None, query_charge=1)
    assert result['reason'] == 'query_polarity_unavailable'
    assert result['ranking_eligible'] is False


def test_explicit_signed_charge_can_supply_missing_polarity():
    result = context(query_adduct=None, query_polarity=None, query_charge='1+')
    assert result['reason'] == 'unique_formula_mass_consistent_adduct'
    assert result['query_polarity'] == 'POS'
    assert result['independent_query_evidence'] is False


def test_actual_ranker_rejects_input_before_attempting_model_load():
    best={'formula':'C6H12O6','adduct':'[M+H]+','inchikey':'AAAAAAAAAAAAAA-UHFFFAOYSA-N'}
    result,evidence=apply_isopair_ranker(best,[best],[[100,20]],query_precursor_mz=181.070664,
        query_adduct='BOGUS',query_charge=1,query_polarity='POS')
    assert result is best and evidence['reason']=='unsupported_explicit_query_adduct'
    assert evidence['evaluated'] is False


def rows():
    return [{'spectrum_id':str(i),'precursor_mz':300.,'retention_time':5.+i*.005,
             'ion_mode':'POS','precursor_charge':charge,'collision_energy':20.,
             'collision_energy_unit':unit} for i,(charge,unit) in enumerate([(1,'eV'),(2,'eV'),(1,'NCE')])]


def test_charge_partitions_events_and_energy_units_are_not_collapsed():
    data=rows();summary=annotate_precursor_features(data)
    assert data[0]['precursor_feature_id']==data[2]['precursor_feature_id']
    assert data[0]['precursor_feature_id']!=data[1]['precursor_feature_id']
    assert data[0]['precursor_feature_condition_count']==2
    assert summary['precursor_features']==2
    assert summary['multiple_acquisition_condition_features']==1
    assert summary['multi_energy_precursor_features']==0  # units differ; no calibrated comparison


def test_same_unit_different_settings_still_count_as_multiple_energies():
    data=[rows()[0],dict(rows()[0],spectrum_id='b',collision_energy=30.)]
    assert annotate_precursor_features(data)['multi_energy_precursor_features']==1


def test_unknown_charge_does_not_bridge_incompatible_known_charges():
    data=rows();data[2]['precursor_charge']=None
    assert annotate_precursor_features(data)['precursor_features']==3


def test_same_ms1_reference_does_not_merge_different_precursor_ions():
    data=rows()
    for row in data:row['precursor_scan_ref']='scan=100'
    data[2]['precursor_mz']=500.
    graph=build_entity_graph('synthetic',data)
    assert graph.summary['precursor_event_count']==3
    assert build_entity_graph('synthetic',list(reversed(data))).graph_sha256==graph.graph_sha256


def test_library_module_globals_are_not_patched():
    legacy=importlib.import_module('matching.recap_adduct')
    assert legacy.resolve_isopair_adduct_context is not resolve_isopair_adduct_context
    assert legacy.canonical_adduct_key('BOGUS')=='BOGUS'
