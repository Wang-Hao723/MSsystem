import copy
import pytest
from eval.v16_runtime_candidate_bridge_20260907.predict_jsonl import validate_request,score_request
from eval.v16_runtime_candidate_bridge_20260907.run import build_input,identity

def request():
 p=identity('COc1ccccc1')
 return {'query_id':'test','precursor_mz':107.049,'polarity':'NEG','peaks':[[50.,1.],[60.,.5],[70.,.3]],'ion_status':'no_independent_ion_evidence','adduct':None,'current_order':[p['standardized_connectivity']],'current_candidates':[{'connectivity':p['standardized_connectivity'],'smiles':'COc1ccccc1','formula':'C7H8O'}]}

def test_truth_field_rejected():
 r=request();r['truth_connectivity']='anything'
 with pytest.raises(ValueError,match='unexpected_request_fields'):validate_request(r)

def test_historical_candidate_union_not_accepted():
 r=request();r['current_candidates'].append(dict(r['current_candidates'][0],connectivity='extra'))
 with pytest.raises(ValueError,match='candidate_count_mismatch'):validate_request(r)

@pytest.mark.parametrize('status',['network_neutral_mass_conflict','polarity_conflict','charge_conflict','conflicting_network_fields','unsupported_or_conflicting_independent_ions','invalid_neutral_mass'])
def test_all_existing_conflict_statuses_abstain_without_library_access(status):
 r=request();r['ion_status']=status
 out=score_request(r,None,None)
 assert out['after']==r['current_order'] and out['status']=='abstained_ion_conflict'

class Retriever:
 def __init__(self,records):self.records=records
 def _fetch_window(self,mz,pol):return self.records,{'read_only':True}

def extra_record():
 smiles='Cc1ccccc1O';p=identity(smiles)
 assert p['valid'] and p['normalized_formula']=='C7H8O'
 return {'connectivity_identity':p['standardized_connectivity'],'smiles':smiles,'molecular_formula':'C7H8O','inchikey':'','library_polarity':'NEG','precursor_mz':107.049,'adduct':'[M-H]-','peaks':[[50.,.9],[60.,.4],[70.,.1]],'spectrum_id':'r1'}

def test_pool_uses_indexed_candidates_and_rejects_wrong_polarity():
 q=request();r=extra_record();bad=dict(r,library_polarity='POS');ret=Retriever([bad,r])
 out=build_input(q,q['current_candidates'],q['current_order'],ret)
 assert out['pool']==q['current_order']+[r['connectivity_identity']]
 assert out['retrieval_audit']['rejections']['polarity_or_mass']==1

def test_missing_current_formula_does_not_borrow_formula_from_library():
 q=request();q['current_candidates'][0]['formula']=None
 out=build_input(q,q['current_candidates'],q['current_order'],Retriever([extra_record()]))
 assert out['pool']==q['current_order']

def test_query_mirror_cannot_admit_extra_candidate():
 q=request();r=extra_record();r['peaks']=q['peaks']
 out=build_input(q,q['current_candidates'],q['current_order'],Retriever([r]))
 assert out['pool']==q['current_order']

def test_duplicate_reference_does_not_change_pool():
 q=request();r=extra_record()
 a=build_input(q,q['current_candidates'],q['current_order'],Retriever([r]))
 b=build_input(q,q['current_candidates'],q['current_order'],Retriever([r,r]))
 assert a['pool']==b['pool']

@pytest.mark.parametrize('peaks',[[],[[50.,1.],[60.,.5]],[[50.,0.],[60.,0.],[70.,0.]]])
def test_empty_or_sparse_signal_abstains_without_library_access(peaks):
 r=request();r['peaks']=peaks
 out=score_request(r,None,None)
 assert out['after']==r['current_order'] and out['status']=='abstained_insufficient_query_signal'
