from eval.v16_selective_evidence_20260907.run import accept,evidence
from eval.v16_selective_evidence_20260907.evaluate import summary

def test_requires_highest_reference_three_matches():
 assert not accept({'score':.99,'margin':.5,'matches':2,'references':5},(.3,.1,2))

def test_two_references_and_margin_required():
 assert not accept({'score':.99,'margin':.5,'matches':3,'references':1},(.3,.1,2))
 assert not accept({'score':.99,'margin':.09,'matches':3,'references':2},(.3,.1,2))
 assert accept({'score':.3,'margin':.1,'matches':3,'references':2},(.3,.1,2))

def test_empty_ranking_is_not_evidence():
 e=evidence({},[])
 assert not accept(e,(.3,.1,2))

def test_withdrawing_error_is_not_correction():
 rows=[{'query_id':'a','ranking':['wrong'],'accepted':False},{'query_id':'b','ranking':['right'],'accepted':False},{'query_id':'c','ranking':[],'accepted':False}]
 r=summary(rows,{k:['right'] for k in ('a','b','c')})
 assert r['queries']==3 and r['wrong_withdrawn']==1 and r['correct_withdrawn']==1 and r['empty_original']==1
 assert r['correct_accepted']==0 and r['accepted_precision'] is None

def test_acceptable_isomer_group_scored_as_given():
 r=summary([{'query_id':'a','ranking':['b'],'accepted':True}],{'a':['a','b']})
 assert r['correct_before']==1 and r['correct_accepted']==1
