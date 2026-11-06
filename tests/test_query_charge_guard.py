"""Query-only charge/adduct safety; synthetic files and no production library."""
import base64
import struct
import importlib
import inspect
from unittest.mock import patch
import pytest
from pipeline import query_charge_guard as guard
from pipeline.query_search import search
from pipeline.mzml_input import parse_mzml_file

@pytest.mark.parametrize('mode,adduct,charge',[
 ('POS','[M+H]+',1),('POS','[M+H]+','1'),('POS','[M+H]+',1.0),
 ('POS','[M+H]+','+1'),('POS','[M+H]+','1+'),
 ('NEG','[M-H]-',1),('NEG','[M-H]-','1'),('NEG','[M-H]-',1.0),
 ('NEG','[M-H]-',-1),('NEG','[M-H]-','-1'),('NEG','[M-H]-','1-'),
 ('POS','[M+Na]+',None),('NEG','[M+Cl]-',None),
 ('POS','[M+H]+',''),('NEG','[M-H]-',''),
 ('POS','[M+H]+',0),('NEG','[M-H]-',0),('NEG','[M-H]-','0'),
 ('POS','[M+H]+',0.0),
 (None,'[M-H]-',1),(None,'[M+H]+','1+'),
 ('POS','[M+2H]2+',2),('NEG','[M-2H]2-',2),
 ('NEG','[M-2H]2-','2-'),('POS','',2),('NEG',None,2),
])
def test_valid_signed_or_magnitude_input_is_preserved(mode,adduct,charge):
 row={'adduct':adduct,'precursor_charge':charge,'polarity':mode,'adduct_reliable':True,
      'raw_peaks':[[151.,99.]],'extra':'keep'}
 out=guard.guard_query_record(row)
 assert not out['adduct_charge_evidence']['conflict']
 assert all(out[key]==value for key,value in row.items())
 assert row=={k:v for k,v in out.items() if k!='adduct_charge_evidence'}
 if charge in (1,1.0,'1',2):assert out['adduct_charge_evidence']['charge_polarity'] is None

@pytest.mark.parametrize('mode,adduct,charge,reason',[
 ('POS','[M+H]+',2,'explicit_adduct_charge_magnitude_conflict'),
 ('NEG','[M-H]-',2,'explicit_adduct_charge_magnitude_conflict'),
 ('POS','[M+H]+','2+','explicit_adduct_charge_magnitude_conflict'),
 ('NEG','[M-H]-','2-','explicit_adduct_charge_magnitude_conflict'),
 ('NEG','[M-H]-','+1','query_ion_polarity_conflict'),
 ('NEG','[M-H]-','1+','query_ion_polarity_conflict'),
 ('POS','[M+H]+',-1,'query_ion_polarity_conflict'),
 ('POS','[M+H]+','1-','query_ion_polarity_conflict'),
 ('NEG','[M+H]+',1,'query_ion_polarity_conflict'),
 ('NEG','[M-H]-',1.5,'invalid_query_charge'),
 ('POS','[M+H]+',float('inf'),'invalid_query_charge'),
 ('POS','[M+H]+','NaN','invalid_query_charge'),
 ('POS','[M+H]+',True,'invalid_query_charge'),
 ('POS','[M+H]+',False,'invalid_query_charge'),
 ('POS','[M+H]+','2 or 3','invalid_query_charge'),
])
def test_conflict_disables_reliability_but_preserves_original_metadata(mode,adduct,charge,reason):
 row={'adduct':adduct,'precursor_charge':charge,'polarity':mode,'adduct_reliable':True,'peaks':[[100.,4.]]}
 out=guard.guard_query_record(row)
 assert out['adduct_reliable'] is False
 assert reason in out['adduct_charge_evidence']['reasons']
 assert out['adduct']==adduct and out['precursor_charge']==charge and out['polarity']==mode
 assert out['peaks']==row['peaks'] and row['adduct_reliable'] is True

@pytest.mark.parametrize('old_flag',[False,None])
def test_existing_parser_rejection_never_promoted(old_flag):
 row={'adduct':'[M+H]+','precursor_charge':1,'polarity':'POS','adduct_reliable':old_flag}
 assert guard.guard_query_record(row)['adduct_reliable'] is old_flag



def cv(accession,name,value=''):
 return f'<cvParam cvRef="MS" accession="{accession}" name="{name}" value="{value}"/>'

def write_mzml(path,charge=2,adduct='[M+H]+',polarity='POS'):
 mz=base64.b64encode(struct.pack('<3d',50.,75.,100.)).decode();iv=base64.b64encode(struct.pack('<3d',100.,70.,50.)).decode()
 ch='' if charge is None else cv('MS:1000041','charge state',charge)
 path.write_text(f'''<?xml version="1.0"?><mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
 <instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">{cv('MS:1001911','Q Exactive')}</instrumentConfiguration></instrumentConfigurationList>
 <run id="test" defaultInstrumentConfigurationRef="IC1"><spectrumList count="1"><spectrum id="scan=7" index="0" defaultArrayLength="3">
 {cv('MS:1000511','ms level',2)}{cv('MS:1000127','centroid spectrum')}{cv('MS:1000130' if polarity=='POS' else 'MS:1000129','positive scan' if polarity=='POS' else 'negative scan')}
 <userParam name="adduct" value="{adduct}"/><userParam name="precursor_type" value="{adduct}"/>
 <scanList count="1"><scan>{cv('MS:1000016','scan start time',2.)}</scan></scanList>
 <precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>{cv('MS:1000744','selected ion m/z',151.007276)}{ch}<userParam name="adduct" value="{adduct}"/></selectedIon></selectedIonList><activation>{cv('MS:1000422','beam-type collision-induced dissociation')}</activation></precursor></precursorList>
 <binaryDataArrayList count="2"><binaryDataArray encodedLength="32">{cv('MS:1000523','64-bit float')}{cv('MS:1000576','no compression')}{cv('MS:1000514','m/z array')}<binary>{mz}</binary></binaryDataArray><binaryDataArray encodedLength="32">{cv('MS:1000523','64-bit float')}{cv('MS:1000576','no compression')}{cv('MS:1000515','intensity array')}<binary>{iv}</binary></binaryDataArray></binaryDataArrayList>
 </spectrum></spectrumList></run></mzML>''')

def search_probe(search, *, adduct='[M+H]+', charge=2, polarity='POS', reliable=True):
 """Run actual production query wrapper through lookup, intercept pre-score."""
 base=importlib.import_module('matching.search')
 record={'precursor_mz':172.989218,'neutral_mass':150.,'neutral_mass_source':'compound_exact_mass',
         'polarity':'POS','adduct':'[M+Na]+','source_tier':'experimental','dataset_source':'synthetic'}
 records=[record];mz_index={int(record['precursor_mz']/base.MZ_BUCKET_DALTON):[0]}
 nm_index={int(150./base.MZ_BUCKET_DALTON):[0]}
 seen={'neutral_lookup_calls':[], 'entered_scoring':False}
 actual=base._neutral_candidates
 def neutral_spy(nm,*a,**k):
  result=actual(nm,*a,**k);seen['neutral_lookup_calls'].append({'neutral_mass':nm,'reference_indices':list(result)});return result
 class ReachedScoring(Exception):pass
 def stop(*a,**k):seen['entered_scoring']=True;raise ReachedScoring
 with patch.multiple(base,_lib_records=records,_mz_index=mz_index,_nm_index=nm_index,
                     _lib_source_families=(),USE_NEUTRAL_MASS_MATCHING=True,
                     _neutral_candidates=neutral_spy,process_peaks=stop):
  try:
   result=search([[50.,100.],[75.,70.],[100.,50.]],151.007276,
                 query_adduct=adduct,query_charge=charge,query_polarity=polarity,
                 query_adduct_reliable=reliable)
   seen['result_is_none']=result is None
  except ReachedScoring:seen['result_is_none']=None
 return seen


@pytest.mark.parametrize('charge',[None,0,1,'1+','+1'])
def test_valid_cross_adduct_retrieval_keeps_expected_candidate(charge):
    result=search_probe(search,charge=charge)
    assert result['neutral_lookup_calls']==[{'neutral_mass':150.,'reference_indices':[0]}]
    assert result['entered_scoring']

@pytest.mark.parametrize('charge',[2,'2+','+2',-1,'1-','-1'])
@pytest.mark.parametrize('flag',[True,None])
def test_conflict_is_blocked_before_neutral_lookup(charge,flag):
    result=search_probe(search,charge=charge,reliable=flag)
    assert result['neutral_lookup_calls']==[] and result['result_is_none']

@pytest.mark.parametrize('mode,charge',[('POS',1),('NEG',1),('POS',2),('NEG',2),('POS',None),('NEG',None),('POS',0),('NEG',0)])
def test_real_mzml_keeps_charge_polarity_and_all_measured_peaks(tmp_path,mode,charge):
    path=tmp_path/'synthetic.mzML'
    write_mzml(path,charge=charge,adduct='[M+H]+' if mode=='POS' else '[M-H]-',polarity=mode)
    row,=parse_mzml_file(path)
    assert row['precursor_charge']==charge and row['polarity']==mode
    assert row['peaks']==[(50.,100.),(75.,70.),(100.,50.)]
    assert 'adduct_charge_evidence' in row
    assert not row['adduct_charge_evidence']['conflict']


def test_positional_search_metadata_cannot_bypass_guard():
    base=importlib.import_module('matching.search').search
    bound=inspect.signature(base).bind_partial([[100.,1.]],151.007276,query_adduct='[M+H]+',query_charge=2,query_adduct_reliable=True)
    bound.apply_defaults()
    args,kw,evidence=guard.guard_search_arguments(base,tuple(bound.arguments.values()),{})
    rebuilt=inspect.signature(base).bind_partial(*args,**kw)
    assert rebuilt.arguments['query_adduct_reliable'] is False
    assert rebuilt.arguments['query_adduct']=='[M+H]+' and rebuilt.arguments['query_charge']==2
    assert evidence['conflict']


def test_normal_multicharge_keeps_direct_precursor_retrieval(monkeypatch):
    base=importlib.import_module('matching.search')
    record={'precursor_mz':151.007276,'polarity':'NEG','source_tier':'experimental','dataset_source':'synthetic'}
    records=[record];index={int(151.007276/base.MZ_BUCKET_DALTON):[0]}
    class ReachedScoring(Exception):pass
    def stop(*a,**kw):
        assert kw['precursor_charge']==2
        raise ReachedScoring
    monkeypatch.setattr(base,'process_peaks',stop)
    with pytest.raises(ReachedScoring):
        search([[50.,100.]],151.007276,records=records,mz_index=index,
          query_adduct='[M-2H]2-',query_charge=2,query_polarity='NEG',query_adduct_reliable=True)

@pytest.mark.parametrize('mode,charge,conflict',[('POS',1,False),('NEG',1,False),('POS',-1,True),('NEG',-1,False),('POS',2,False),('NEG',2,False),('POS',-2,True),('NEG',-2,False),('POS','+1',False),('NEG','+1',True),('POS','+2',False),('NEG','+2',True),('POS',0,False),('NEG',0,False),('POS','+0',False),('NEG','+0',False),('POS','-0',False),('NEG','-0',False),('POS',None,False),('NEG',None,False),('POS','',False),('NEG','',False)])
def test_actual_xml_charge_sign_and_unknown_are_retained(tmp_path,mode,charge,conflict):
    path=tmp_path/'charge_values.mzML';write_mzml(path,charge=charge,polarity=mode)
    row,=parse_mzml_file(path)
    assert row['polarity']==mode
    assert row['precursor_charge']==(None if charge in (None,'') else abs(int(charge)))
    if charge is not None:
        assert row['reported_precursor_charge']==str(charge)
    else:
        assert 'reported_precursor_charge' not in row
    assert row['adduct_charge_evidence']['conflict'] is conflict
    if charge in (0,''):
        assert row['precursor_charge_status']==('unknown_zero' if charge==0 else 'missing_value')
    if conflict:
        assert row['adduct_reliable'] is False

@pytest.mark.parametrize('mode,adduct',[('POS','[M+H]+'),('NEG','[M-H]-')])
@pytest.mark.parametrize('charge',['+0','-0','0+','0-'])
def test_signed_zero_supplies_no_polarity_evidence(mode,adduct,charge):
    row={'adduct':adduct,'precursor_charge':charge,'polarity':mode,'adduct_reliable':True}
    out=guard.guard_query_record(row)
    assert out['adduct_reliable'] is True
    assert out['precursor_charge']==charge
    assert not out['adduct_charge_evidence']['conflict']
    assert out['adduct_charge_evidence']['charge_polarity'] is None
    assert out['adduct_charge_evidence']['charge_status']=='unknown_zero'

@pytest.mark.parametrize('mode',['POS','NEG'])
@pytest.mark.parametrize('charge',['1.5','-1.5','1e0','1.0'])
def test_xml_noninteger_charge_literal_is_rejected_before_truncation(tmp_path,mode,charge):
    path=tmp_path/'noninteger_charge.mzML';write_mzml(path,charge=charge,polarity=mode)
    with pytest.raises(ValueError,match='integer literal'):
        parse_mzml_file(path)
