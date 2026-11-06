"""Adversarial leakage, full-content, and split boundary tests for public training."""
from copy import deepcopy
import numpy as np
import pytest
from eval.v16_public_qtof_20260907.curate import fingerprint,identity,full_peaks
from eval.v16_public_qtof_20260907.prepare import check_splits,split_for_formula
from eval.v16_public_qtof_20260907.comparator import panel_features


def example():
    query={'peaks':[[50,1],[70,.4],[90,.2],[110,.1]],'precursor_mz':200}
    candidates=[{'connectivity':'a','references':[{'peaks':[[50,1],[70,.3],[90,.4],[120,.1]],'precursor_mz':200}]},
                {'connectivity':'b','references':[{'peaks':[[55,1],[70,.2],[100,.3],[120,.1]],'precursor_mz':200}]}]
    return query,candidates


def test_full_content_dedup_invariant_to_scale_order_and_exact_mass_duplicates():
    p=[[50,2],[70,3],[70,1],[90,.00001]]
    assert fingerprint(p)==fingerprint([[90,.00002],[70,8],[50,4]])
    assert fingerprint(p)!=fingerprint(p[:-1])


@pytest.mark.parametrize('field',['truth','ground_truth','truth_connectivity','acceptable_connectivities','label','compound_name'])
def test_labels_rejected_recursively(field):
    q,c=example();c[1]['references'][0][field]='injected'
    with pytest.raises(ValueError,match='forbidden'):
        panel_features(q,c)


def test_candidate_reversal_and_duplicate_reference_do_not_change_features():
    q,c=example();a,status=panel_features(q,c)
    assert status=='scored'
    d=deepcopy(c[::-1]);d[0]['references']*=3
    b,_=panel_features(q,d)
    np.testing.assert_allclose(a,b[::-1],rtol=0,atol=1e-14)


def test_identity_names_do_not_change_features():
    q,c=example();a,_=panel_features(q,c)
    c[0]['connectivity']='different-label';c[1]['connectivity']='anything'
    b,_=panel_features(q,c)
    np.testing.assert_array_equal(a,b)


def test_processed_mirror_not_accepted_as_independent_reference():
    q,c=example()
    c[0]['references']=[{'peaks':q['peaks']+[[200,2]],'precursor_mz':200}]
    x,status=panel_features(q,c)
    assert x is None and 'nonmirror' in status


@pytest.mark.parametrize('field',['connectivity','formula','fingerprint'])
def test_cross_split_leak_rejected(field):
    a={'connectivity':'a','formula':'A','fingerprint':'1','split':'train','role':'reference'}
    b={'connectivity':'b','formula':'B','fingerprint':'2','split':'evaluation','role':'query'}
    b[field]=a[field]
    with pytest.raises(ValueError,match='crosses splits'):
        check_splits([a,b])


def test_query_reference_leak_rejected_within_split():
    a={'connectivity':'a','formula':'A','fingerprint':'1','split':'train','role':'reference'}
    with pytest.raises(ValueError,match='query/reference'):
        check_splits([a,{**a,'role':'query'}])


def test_structural_tautomer_and_stereo_exclusion_keys():
    a=identity('CC(=O)C');b=identity('C=C(O)C')
    assert a['standardized_connectivity']==b['standardized_connectivity']
    a=identity('C[C@H](O)C(=O)O');b=identity('C[C@@H](O)C(=O)O')
    assert a['standardized_connectivity']==b['standardized_connectivity']
    assert not identity('CCO','C3H8O')['valid']


def test_split_never_uses_polarity_or_query_labels():
    assert split_for_formula('C10H12O2')==split_for_formula('C10H12O2')


def test_invalid_peak_rejected():
    with pytest.raises(ValueError):
        full_peaks([[50,float('nan')]])


def test_public_prediction_does_not_open_any_label_file(monkeypatch):
    import eval.v16_public_qtof_20260907.train_predict as tp
    def forbidden_read(path):
        raise AssertionError(f'Prediction attempted file read: {path}')
    monkeypatch.setattr(tp,'read',forbidden_read)
    class FixedModel:
        def predict_proba(self,x):
            return np.array([[.3,.7],[.8,.2]])
    q,c=example();x,_=panel_features(q,c)
    records=[{'query_id':'arbitrary','polarity':'POS','split':'evaluation',
              'keys':['a','b'],'x':x.tolist(),'status':'scored'}]
    result=tp.predict_public(records,{'POS':FixedModel()})
    assert result[0]['model_order']==['a','b']


def test_fit_opens_only_training_labels(monkeypatch,tmp_path):
    import eval.v16_public_qtof_20260907.train_predict as tp
    rows=[];labels=[]
    for pol in ('POS','NEG'):
        for i in range(30):
            key=f'{pol}-{i}'
            rows.append({'query_id':key,'polarity':pol,'split':'train',
                'keys':[key,'negative'],'x':[[1.],[0.]],'status':'scored'})
            labels.append({'query_id':key,'connectivity':key})
    # A malicious holdout row must never enter fit, even with absurd features.
    rows.append({'query_id':'evaluation-secret','polarity':'POS','split':'evaluation',
                 'keys':['secret','negative'],'x':[[999.],[0.]],'status':'scored'})
    opened=[]
    def guarded_read(path):
        opened.append(path.name)
        assert path.name=='train_labels.json'
        return labels
    class FitSpy:
        def __init__(self,**kwargs):
            pass
        def fit(self,x,y,sample_weight):
            assert len(x)==60 and x.max()==1
            assert sum(y)==30
    monkeypatch.setattr(tp,'read',guarded_read)
    monkeypatch.setattr(tp,'HistGradientBoostingClassifier',FitSpy)
    monkeypatch.setattr(tp,'ASSETS',tmp_path)
    monkeypatch.setattr(tp,'OUT',tmp_path)
    monkeypatch.setattr(tp.joblib,'dump',lambda m,p:p.write_text('test model'))
    monkeypatch.setattr(tp,'write',lambda p,v:None)
    (tmp_path/'protocol.json').write_text('{}')
    tp.fit(rows,{'model':{'type':'test'}})
    assert opened==['train_labels.json']


@pytest.mark.parametrize('status',['network_neutral_mass_conflict','polarity_conflict',
                                  'charge_conflict','invalid_neutral_mass',
                                  'unsupported_or_conflicting_independent_ions'])
def test_all_explicit_ion_conflicts_abstain(status):
    from eval.v16_public_qtof_20260907.train_predict import ion_conflict
    assert ion_conflict(status)
    assert not ion_conflict('supported')
    assert not ion_conflict('no_independent_ion_evidence')
