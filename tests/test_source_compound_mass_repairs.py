"""Portable runtime tests: immutable DB, bound literals and atomic publication."""
import copy,hashlib,importlib,json,sqlite3
from pathlib import Path
from types import SimpleNamespace
import pytest
from runtime import source_compound_mass_repairs as mass
from runtime import source_compound_mass_core as core
from runtime import library_metadata_repairs as metadata

def _sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

@pytest.fixture
def fixture(tmp_path):
    source=tmp_path/'one.mgf';smi='CCO';key='LFQSCWFLJHTTHZ-UHFFFAOYSA-N';formula='C2H6O';correct=46.041864812;old=6.041864812
    source.write_text(f'BEGIN IONS\nPEPMASS=47.049140812\nSMILES={smi}\nFORMULA={formula}\nEXACTMASS={correct}\nADDUCT=[M+H]+\n10 1\n20 2\nEND IONS\n')
    db=tmp_path/'library.db';conn=sqlite3.connect(db)
    conn.executescript('CREATE TABLE Spectra (id INTEGER PRIMARY KEY,spectrum_id TEXT,source_file TEXT,precursor_mz REAL,charge INTEGER,adduct TEXT,ion_mode TEXT,source_accession TEXT,peaks_json TEXT);CREATE TABLE Compounds(id INTEGER PRIMARY KEY,smiles TEXT,inchikey TEXT,molecular_formula TEXT,exact_mass REAL);CREATE TABLE Mapping(spectrum_id INTEGER,compound_id INTEGER);')
    raw={'spectrum_id':'one','source_file':str(source),'precursor_mz':47.049140812,'charge':1,'adduct':'[M+H]+','ion_mode':'POS','source_accession':None,'inchikey':key,'smiles':smi,'peaks_json':'[[10,1],[20,2]]','molecular_formula':None,'exact_mass':old,'database_spectrum_id':1}
    conn.execute('INSERT INTO Spectra VALUES(?,?,?,?,?,?,?,?,?)',(1,'one',str(source),raw['precursor_mz'],1,'[M+H]+','POS',None,raw['peaks_json']))
    conn.execute('INSERT INTO Compounds VALUES(?,?,?,?,?)',(1,smi,key,None,old));conn.execute('INSERT INTO Mapping VALUES(1,1)');conn.commit();conn.close()
    record={'spectrum_id':'one','inchikey':key,'smiles':smi,'exact_mass':old,'formula':'','neutral_mass':old,'neutral_mass_source':'compound_exact_mass','precursor_mz':raw['precursor_mz'],'charge':1,'adduct':'[M+H]+','peaks':[[10,1],[20,2]]}
    manifest={'record_count':1,'normalized_records_sha256':'a'*64,'record_alignment_sha256':'b'*64,'source_db':{'sha256':_sha(db)}}
    class Records(list):
        def read_uncached(self,i):return self[i]
    records=Records([record]);records.manifest=manifest
    ev={'source_file':str(source),'byte_offset':0,'byte_length':source.stat().st_size,'block_sha256':_sha(source),'formula':formula,'source_exact_mass':correct,'full_inchikey':key}
    row={'record_index':0,'spectrum_id':'one','sql_compound_id':1,'full_inchikey':key,'original_record_sha256':core.digest(record),'original_values':{k:record[k] for k in core.FIELDS},'old_neutral_mass':old,'updates':{'exact_mass':correct,'formula':formula,'neutral_mass':correct,'neutral_mass_source':'compound_exact_mass'}}
    registry={'schema_version':core.SCHEMA,'registry_id':'fixture','base_record_count':1,'base_normalized_records_sha256':manifest['normalized_records_sha256'],'base_record_alignment_sha256':manifest['record_alignment_sha256'],'source_database_sha256':_sha(db),'source_files':{str(source):{'sha256':_sha(source),'size_bytes':source.stat().st_size}},'excluded_quarantined_spectrum_ids':[],'compounds':[{'full_inchikey':key,'smiles':smi,'formula':formula,'source_exact_mass':correct,'original_source_evidence':[ev]}],'rows':[row]};registry['semantic_sha256']=core.digest(registry)
    binding={'record_index':0,'database_compound_id':1,'original_record_sha256':core.digest(record),'original_sql_values':raw,'original_sql_values_sha256':core.digest(raw),'updates':{'exact_mass':correct,'molecular_formula':formula}}
    config={'registry_semantic_sha256':registry['semantic_sha256'],'registered_reference_count':1,'base_normalized_records_sha256':manifest['normalized_records_sha256'],'base_record_alignment_sha256':manifest['record_alignment_sha256'],'source_database_sha256':_sha(db),'base_library_snapshot_id':'fixture','source_database_path':str(db),'sql_bindings_sha256':'c'*64,'semantic_sha256':'d'*64}
    class Raw:
        db_path=db;library_snapshot_id='fixture';_sqlite_query_strategy={'from_clause':'Spectra s JOIN Mapping m ON s.id=m.spectrum_id JOIN Compounds c ON c.id=m.compound_id'}
        def _connect(self):
            c=sqlite3.connect(db.as_uri()+'?mode=ro',uri=True);c.row_factory=sqlite3.Row;return c
    return SimpleNamespace(source=source,db=db,raw=raw,record=record,records=records,manifest=manifest,registry=registry,bindings=[binding],config=config,retriever=Raw(),correct=correct)

def _state(f,**kw):return mass.MassRepairState(f.config,f.registry,f.bindings,f.records,f.retriever,verify_sources=False,**kw)

def test_full_literal_and_sql_binding(fixture):
    f=fixture;before=_sha(f.db);state=mass.MassRepairState(f.config,f.registry,f.bindings,f.records,f.retriever,verify_sources=True)
    result=state.apply_sql(f.raw);assert result['exact_mass']==f.correct;assert result['molecular_formula']=='C2H6O';assert f.raw['exact_mass']!=f.correct;assert _sha(f.db)==before
    assert not ({'exact_mass','molecular_formula','formula','neutral_mass'}&metadata.FIELDS)

@pytest.mark.parametrize('field,value',[('exact_mass',99.),('molecular_formula','H2O'),('source_file','other'),('spectrum_id','other'),('inchikey','other'),('smiles','CC'),('precursor_mz',99.),('charge',2),('adduct','[M+Na]+'),('peaks_json','[]')])
def test_changed_original_sql_rejected(fixture,field,value):
    f=fixture;state=_state(f)
    with pytest.raises(ValueError):state.apply_sql({**f.raw,field:value})

def test_nonrepair_sql_exact(fixture):
    state=_state(fixture);row={**fixture.raw,'database_spectrum_id':999};assert state.apply_sql(row)==row

def test_base_alignment_rejected(fixture):
    fixture.records.manifest={**fixture.manifest,'record_alignment_sha256':'x'*64}
    with pytest.raises(ValueError,match='base binding'):_state(fixture)

def test_metadata_overlap_rejected(fixture):
    with pytest.raises(ValueError,match='Overlapping'):_state(fixture,metadata_indices=[0])

def test_duplicate_mapping_rejected(fixture):
    conn=sqlite3.connect(fixture.db);conn.execute('INSERT INTO Mapping VALUES(1,1)');conn.commit();conn.close()
    with pytest.raises(ValueError,match='one-to-one'):_state(fixture)

def test_raw_source_changed_rejected(fixture):
    f=fixture;f.source.write_text(f.source.read_text().replace('C2H6O','CH4O'))
    with pytest.raises(ValueError):mass.MassRepairState(f.config,f.registry,f.bindings,f.records,f.retriever,verify_sources=True)

def test_dependency_change_rejected(fixture,tmp_path):
    dep=tmp_path/'cfg.json';dep.write_text('{}');state=_state(fixture,dependencies=[dep]);dep.write_text('{"changed":1}')
    with pytest.raises(ValueError,match='dependencies changed'):state.assert_fresh()

def _loader_fixture(f,monkeypatch,tmp_path):
    state=_state(f);monkeypatch.setitem(mass._STATES,str(mass.ROOT.resolve()),state)
    empty={'schema_version':'mssystem.source-metadata-repairs/1.0.0','registry_id':'fixture-meta','rows':[],'source_files':{},'base_record_count':1,'base_normalized_records_sha256':f.manifest['normalized_records_sha256'],'source_database_sha256':f.manifest['source_db']['sha256']};empty['semantic_sha256']=core.digest(empty)
    rp=tmp_path/'metadata.json';rp.write_text(json.dumps(empty));monkeypatch.setattr(metadata,'REGISTRY_PATH',rp)
    search=importlib.import_module('matching.search');index={int(f.raw['precursor_mz']/search.MZ_BUCKET_DALTON):[0]};neutral={int(f.raw['exact_mass']/search.MZ_BUCKET_DALTON):[0]}
    for name,value in [('_lib_records',None),('_mz_index',index),('_nm_index',neutral),('_lib_db_path',None),('_lib_source_families',())]:monkeypatch.setattr(search,name,value)
    monkeypatch.setattr(search,'load_library',lambda *a:(search._lib_records if search._lib_records is not None else f.records,index))
    return state,search,index,neutral

def test_atomic_normalized_records_and_neutral_index(fixture,monkeypatch,tmp_path):
    f=fixture;state,search,index,oldneutral=_loader_fixture(f,monkeypatch,tmp_path)
    records,mz=metadata.load_library_with_verified_repairs();assert mz is index and search._mz_index is index;assert records[0]['neutral_mass']==f.correct
    assert 0 not in search._neutral_candidates(f.raw['exact_mass'],records,search._nm_index,tol_da=.001,tol_ppm=0)
    assert search._neutral_candidates(f.correct,records,search._nm_index,tol_da=.001,tol_ppm=0)==[0]
    assert f.records[0]['neutral_mass']==f.raw['exact_mass'];assert mass.source_compound_mass_status(records,True)['normalized_overlay_applied']
    again,again_index=metadata.load_library_with_verified_repairs();assert again is records and again_index is index

def test_failed_joint_publication_is_not_reported_applied(fixture,monkeypatch,tmp_path):
    f=fixture;state,search,index,oldneutral=_loader_fixture(f,monkeypatch,tmp_path)
    original=metadata.VerifiedMetadataRecords
    class FailMetadata(original):
        def __init__(self,*a,**k):raise ValueError('injected later metadata failure')
    monkeypatch.setattr(metadata,'VerifiedMetadataRecords',FailMetadata)
    with pytest.raises(ValueError,match='injected'):metadata.load_library_with_verified_repairs()
    assert search._lib_records is search._mz_index is search._nm_index is None
    assert mass.source_compound_mass_status(None,False)['normalized_overlay_applied'] is False
