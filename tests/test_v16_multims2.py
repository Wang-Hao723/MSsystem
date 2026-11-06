import json
from eval.v16_multims2_20260907.score import condition_panel
from eval.v16_multims2_20260907.prepare import file_role
from matchms.importing import load_from_mgf
from eval.v16_multims2_20260907.multispectrum import sample_id,companions

def test_condition_missing_candidate_falls_back_whole_panel():
 q={'activation':'CID','panel':{'a':[{'activation':'CID'}],'b':[{'activation':'EAD'}]}}
 assert condition_panel(q)==(q['panel'],False)

def test_condition_keeps_every_candidate():
 q={'activation':'EAD','panel':{'a':[{'activation':'CID'},{'activation':'EAD'}],'b':[{'activation':'EAD'}]}}
 p,applied=condition_panel(q)
 assert applied and set(p)==set(q['panel']) and all(r['activation']=='EAD' for rr in p.values() for r in rr)

def test_file_split_is_record_independent():
 assert file_role('injection.mzML')==file_role('injection.mzML')
 assert {file_role(f'injection{i}.mzML') for i in range(20)}=={'reference','query'}

def test_mgf_field_aliases_even_without_harmonization(tmp_path):
 p=tmp_path/'test.mgf';p.write_text('BEGIN IONS\nFILENAME=injection.mzML\nINCHIAUX=ABC\nPRECURSOR_MZ=100\nCOLLISION_ENERGY=[12.0]\n50 1\n60 2\n70 3\nEND IONS\n')
 m=next(load_from_mgf(p,metadata_harmonization=False)).metadata
 assert m['file_name']=='injection.mzML' and m['inchi_aux']=='ABC'
 assert json.loads(m['collision_energy'])==[12.0]

def test_sample_coordinates_preserved():
 assert sample_id('20230912_plate_1_pos_H10_CID_40ev.mzML')==sample_id('20230912_plate_1_pos_H10_EAD_12ev.mzML')
 assert sample_id('20230912_plate_1_pos_H10_CID_40ev.mzML')!=sample_id('20230912_plate_1_pos_H11_CID_40ev.mzML')
 assert sample_id('unknown.mzML') is None

def test_multispectrum_rejects_mass_or_candidate_mismatch():
 q={'filename':'plate_H10_CID_40ev.mzML','polarity':'POS','panel':{'a':[],'b':[]},'query':{'precursor_mz':100.}}
 same=q|{'filename':'plate_H10_EAD_12ev.mzML'}
 mass=same|{'query':{'precursor_mz':101.}}
 candidate=same|{'panel':{'a':[]}}
 assert companions(q,[q,same,mass,candidate])==[q,same]
