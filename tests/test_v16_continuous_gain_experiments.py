import numpy as np
import pytest
from eval.v16_incumbent_task_20260907.run import candidate_features,pair as model_pair,choose
from eval.v16_full_pool_consensus_20260907.run import consensus
from eval.v16_full_pool_reference_20260907.run import enter_fixed_head
from eval.v16_peak_jackknife_20260907.run import pair,references
from eval.v16_joint_jackknife_20260907.run import two_features,predict
from eval.v16_public_qtof_20260907.comparator import pair_features

def feature(score):
    x=np.zeros(56);x[0]=np.log1p(4);x[13]=score;x[14]=1.;return x

def peer(qid,run,eligible=True):
    return {'query_id':qid,'run':run,'polarity':'POS','rt':1.,'mz':200.,'before':['A','B'],
      'proposal':['B','A'],'vote_eligible':eligible,'forms':{'A':'C2H6','B':'C2H6'}}

def test_primary_peak_filter_below_three_is_unknown():
    assert candidate_features({'peaks':[[50,1],[70,.01],[90,.001]],'precursor_mz':200},[{'peaks':[[50,1],[70,.5],[90,.2]],'precursor_mz':200}]) is None

def test_incumbent_pair_dimensions():
    a=np.arange(63);b=a+1
    x=model_pair(a,b);assert len(x)==189;np.testing.assert_array_equal(x[-63:],np.ones(63))

def test_model_cannot_promote_unmatched_challenger():
    class AlwaysYes:
        def predict_proba(self,x):return np.array([[0,1]])
    assert choose(['A','B'],{'A':np.r_[feature(.5),np.zeros(7)],'B':np.zeros(63)},AlwaysYes(),.5)==['A','B']

def test_two_independent_supported_runs_required():
    q=peer('q','r1');p=peer('p','r2');dup=peer('p2','r2')
    assert consensus([q,p,dup])[0]['after']==q['before']
    assert consensus([q,p,peer('z','r3')])[0]['after']==q['proposal']

def test_unsupported_incumbent_cannot_cast_vote():
    q=peer('q','r1');r=consensus([q,peer('p','r2'),peer('z','r3',False)])
    assert r[0]['after']==q['before']

def test_conflicting_votes_within_run_abstain():
    q=peer('q','r1');p=peer('p','r2');bad=peer('p2','r2');bad['proposal']=['A','B']
    r=consensus([q,p,bad,peer('z','r3')])[0]
    assert r['after']==q['before'] and r['ambiguous_runs']==['r2']

def test_fixed_head_remains_while_tail_receives_candidate():
    f={k:feature(s) for k,s in zip('ABCDEF',[.1,.2,.3,.4,.5,.9])}
    out=enter_fixed_head(list('ABCDE'),list('ABCDEF'),f,{k:'C2H6' for k in f})
    assert out[0]=='A' and len(out)==5 and len(set(out))==5 and 'F' in out

@pytest.mark.parametrize('seed',range(6))
def test_fast_direct_matches_original_feature_implementation(seed):
    rng=np.random.default_rng(seed);m=np.arange(50.,70.);q=np.c_[m,rng.uniform(.01,1,20)];r=np.c_[m+rng.uniform(-.03,.03,20),rng.uniform(.01,1,20)]
    direct,ratio=pair(q,r);full=pair_features(q,r)
    assert direct==pytest.approx(full[13],abs=1e-12);assert ratio==full[14]

def test_original_mirror_never_reenters_after_perturbation():
    peaks=[[50,1],[70,.5],[90,.2],[110,.1]]
    q={'query':{'peaks':peaks,'precursor_mz':200},'panel':{'A':[{'peaks':peaks,'precursor_mz':200}]}}
    assert references(q)=={'A':[]}

def test_reference_deletion_requires_two_distinct_contents():
    q={'query_id':'q','query':{'peaks':[[50,1],[70,.5],[90,.2],[110,.1]],'precursor_mz':200},
      'before':['A','B'],'pool':['A','B'],'forms':{'A':'C2H6','B':'C2H6'},
      'panel':{k:[{'peaks':[[50,1],[70,.5],[90,.2],[120,.1]],'precursor_mz':200}] for k in ('A','B')}}
    bank=references(q);f,d=two_features(q['query']['peaks'],200,bank)
    assert f['B'] is not None and d['B'] is None
    out=predict(q,q['query']['peaks'],bank)
    assert all(out[f'jackknife_{t}']==q['before'] for t in (.6,.8,1.))
