import gzip
import json
import sqlite3
import pytest
from scripts.consolidate_jiawei_sanwujiao_targets import (
    MERGED_DATASET, consolidate, verify_merge_scope,
)

def database(tmp_path):
    p=tmp_path/'targets.db'
    with sqlite3.connect(p) as c:
        c.executescript('''
        CREATE TABLE TargetCompounds(id INTEGER PRIMARY KEY,target_id TEXT,dataset_source TEXT,compound_name TEXT,molecular_formula TEXT,source_tier TEXT,provenance TEXT,reported_collision_energy_raw TEXT);
        CREATE TABLE DiagnosticIons(id INTEGER PRIMARY KEY,target_id INTEGER REFERENCES TargetCompounds(id),fragment_mz REAL,ordinal INTEGER);
        CREATE TABLE Spectra(id INTEGER PRIMARY KEY,marker TEXT);
        CREATE TABLE Compounds(id INTEGER PRIMARY KEY,marker TEXT);
        CREATE TABLE Mapping(id INTEGER PRIMARY KEY,marker TEXT);
        INSERT INTO Spectra VALUES(1,'preserve');
        ''')
        for i,source in enumerate(['JWSWJ_加味三乌胶_targeted_20261020','JWSWJ_加味三乌胶_全文增补_targeted_20261021_v10','another-prescription'],1):
            c.execute('INSERT INTO TargetCompounds VALUES(?,?,?,?,?,?,?,?)',(i,f't{i}',source,'Glucose','C6H12O6','literature_presence_only',json.dumps({'original_msms':{'products':[89.024,119.035],'units':'eV','source_row':str(i)}}),'40 eV'))
            c.execute('INSERT INTO DiagnosticIons VALUES(?,?,?,?)',(i,i,89.024,1))
    return p

def test_preview_preserves_all_rows_and_writes_no_preimage(tmp_path):
    db=database(tmp_path);before=db.read_bytes();pre=tmp_path/'pre.json.gz'
    r=consolidate(db,pre,tmp_path/'receipt.json')
    assert r['merged_target_rows']==2 and r['source_dataset_count']==2
    assert r['named_formula_entities']==1
    assert not pre.exists() and db.read_bytes()==before

def test_physical_merge_keeps_methods_provenance_and_fragment_foreign_keys(tmp_path,monkeypatch):
    db=database(tmp_path);pre=tmp_path/'pre.json.gz';receipt=tmp_path/'receipt.json'
    monkeypatch.setenv('MSSYSTEM_LIBRARY_WRITE_ENABLED','1')
    r=consolidate(db,pre,receipt,apply=True)
    info,rows,scope=verify_merge_scope(db,pre,r)
    assert info['dataset_source']==MERGED_DATASET and len(rows)==2
    assert r['inserted_target_rows']==r['deleted_target_rows']==0
    with gzip.open(pre,'rt') as f:old=json.load(f)['rows']
    with sqlite3.connect(db) as c:
        assert c.execute('SELECT DISTINCT dataset_source FROM TargetCompounds WHERE id IN(1,2)').fetchall()==[(MERGED_DATASET,)]
        assert c.execute('SELECT target_id,fragment_mz FROM DiagnosticIons ORDER BY id').fetchall()==[(1,89.024),(2,89.024),(3,89.024)]
        assert c.execute('SELECT dataset_source FROM TargetCompounds WHERE id=3').fetchone()[0]=='another-prescription'
        assert c.execute('SELECT marker FROM Spectra').fetchone()[0]=='preserve'
    for row,previous in zip(rows,old):
        p=json.loads(row['provenance']);assert p.pop('original_dataset_source')==previous['dataset_source']
        assert p==json.loads(previous['provenance'])
        assert all(row[k]==v for k,v in previous.items() if k not in {'dataset_source','provenance'})
    assert scope['diagnostics_sha256']==r['pre_database_state']['diagnostics_sha256']

def test_other_table_trigger_is_blocked_and_entire_merge_rolls_back(tmp_path,monkeypatch):
    db=database(tmp_path);monkeypatch.setenv('MSSYSTEM_LIBRARY_WRITE_ENABLED','1')
    with sqlite3.connect(db) as c:
        c.executescript("CREATE TRIGGER bad AFTER UPDATE ON TargetCompounds BEGIN UPDATE Spectra SET marker='changed'; END;")
    with pytest.raises(sqlite3.DatabaseError,match='not authorized'):
        consolidate(db,tmp_path/'pre.json.gz',tmp_path/'receipt.json',apply=True)
    with sqlite3.connect(db) as c:
        assert c.execute('SELECT COUNT(*) FROM TargetCompounds WHERE dataset_source=?',(MERGED_DATASET,)).fetchone()[0]==0
        assert c.execute('SELECT marker FROM Spectra').fetchone()[0]=='preserve'

def test_post_merge_evidence_tampering_is_rejected(tmp_path,monkeypatch):
    db=database(tmp_path);pre=tmp_path/'pre.json.gz';monkeypatch.setenv('MSSYSTEM_LIBRARY_WRITE_ENABLED','1')
    r=consolidate(db,pre,tmp_path/'receipt.json',apply=True)
    with sqlite3.connect(db) as c:c.execute("UPDATE TargetCompounds SET reported_collision_energy_raw='999 eV' WHERE id=1")
    with pytest.raises(RuntimeError,match='differ from their source rows'):verify_merge_scope(db,pre,r)

def test_production_write_requires_explicit_enable(tmp_path):
    db=database(tmp_path)
    with pytest.raises(RuntimeError,match='requires MSSYSTEM_LIBRARY_WRITE_ENABLED'):
        consolidate(db,tmp_path/'pre.json.gz',tmp_path/'receipt.json',apply=True)
