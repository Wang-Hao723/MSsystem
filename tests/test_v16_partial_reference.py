import numpy as np
import pytest
from eval.v16_partial_reference_20260907.experiment import vector,replace,enter_top5,signal,BASELINE_INDEX
from eval.v16_partial_reference_20260907.admission import make_pool

def feat(score):
    x=np.zeros(56);x[0]=np.log1p(6);x[14]=.6;x[BASELINE_INDEX]=score
    return x

def test_missing_incumbent_is_unknown_not_zero():
    before=['a','b'];assert replace(before,{'a':None,'b':feat(.99)})==before

def test_missing_competitor_keeps_its_slot():
    before=['a','unknown','b'];after=replace(before,{'a':feat(.4),'unknown':None,'b':feat(.8)})
    assert after==['b','unknown','a']

def test_expansion_cannot_displace_unknown_boundary():
    before=['a','b','c','d','unknown'];features={k:feat(.4) for k in before};features['unknown']=None;features['new']=feat(.99)
    assert enter_top5(before,before+['new'],features)==before

def test_cross_formula_pair_is_not_scored():
    assert replace(['a','b'],{'a':feat(.4),'b':feat(.9)},formulas={'a':'C2H6','b':'C3H8'})==['a','b']

def test_no_match_never_promotes():
    x=feat(0);x[14]=0
    assert replace(['a','b'],{'a':x,'b':x})==['a','b']
    assert enter_top5(['a'],['a','b'],{'a':x,'b':x})==['a']

def test_large_reference_peak_count_cannot_inflate_signal():
    x=feat(.9);x[0]=np.log1p(3);x[1]=np.log1p(100);x[14]=.1
    assert not signal(x)

def test_signal_needs_three_observed_matches():
    x=feat(.9);x[14]=2/6;assert not signal(x)
    x[14]=3/6;assert signal(x)

def test_processed_mirror_rejected_even_with_different_precursor_peak():
    q={'peaks':[[10,1],[20,.5],[30,.25],[100,.1]],'precursor_mz':100}
    ref={'peaks':[[10,1],[20,.5],[30,.25],[100,.9]],'precursor_mz':100}
    assert vector(q,[ref]) is None

def test_quality_filter_precedes_three_reference_cap():
    q={'peaks':[[10,1],[20,.5],[30,.25],[40,.1]],'precursor_mz':100}
    good={'peaks':[[10,1],[20,.4],[30,.2]],'precursor_mz':100}
    bad=[{'peaks':[[50+i,1]],'precursor_mz':100} for i in range(5)]
    expected=vector(q,[good]);assert expected is not None
    np.testing.assert_allclose(vector(q,bad+[good]),expected)

def test_reference_duplicate_and_order_invariance():
    q={'peaks':[[10,1],[20,.5],[30,.25],[40,.1]],'precursor_mz':100}
    refs=[{'peaks':[[10,1],[20,.4],[30,.2]],'precursor_mz':100},
          {'peaks':[[10,1],[20,.3],[40,.2]],'precursor_mz':100}]
    np.testing.assert_array_equal(vector(q,refs),vector(q,refs[::-1]*3))

@pytest.mark.parametrize('field',['truth','ground_truth','truth_connectivity','acceptable_connectivities','label','compound_name'])
def test_forbidden_fields_rejected(field):
    q={'peaks':[[10,1],[20,.5],[30,.2]],'precursor_mz':100,field:'secret'}
    with pytest.raises(ValueError):vector(q,[])

def test_fixed_pool_budget_and_candidate_permutation():
    cs=[{'connectivity':f'k{i}','origin':'native' if i<5 else 'extended','rank':i+1,'union_rank':i+1,'formula':'C2H6','smiles':'CC'} for i in range(40)]
    q={'current_connectivity':'k0','candidates':cs}
    recall={'neutral_recall':{'candidates':[{'connectivity':f'n{i}','neutral_error_ppm':i,'formula':'C2H6'} for i in range(29)]}}
    before,pool,_=make_pool(q,recall)
    assert len(pool)==40 and len(set(pool))==40 and pool[:5]==before
    assert make_pool(dict(q,candidates=cs[::-1]),recall)[:2]==(before,pool)

def test_expanded_gain_requires_real_evidence():
    before=['a','b','c','d','e'];features={k:feat(.4) for k in before};features['new']=feat(.9)
    after=enter_top5(before,before+['new'],features)
    assert len(after)==5 and after[0]=='new' and 'e' not in after

def test_learned_challenger_can_correct_maximum_cosine_incumbent():
    class Model:
        def predict_proba(self,x):return np.array([[.01,.99]])
    features={'a':feat(.9),'b':feat(.6)}
    assert replace(['a','b'],features)==['a','b']
    assert replace(['a','b'],features,Model(),.9)==['b','a']
