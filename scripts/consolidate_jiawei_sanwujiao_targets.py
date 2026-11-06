#!/usr/bin/env python3
"""Physically merge every Jiawei target dataset while preserving row evidence.

Only dataset_source and the original-dataset provenance key may change. Target
IDs, diagnostic-ion foreign keys, masses, structures and method fields remain
byte-for-byte unchanged. A compressed preimage permits an exact restoration.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json, os, sqlite3, sys
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.install_jiawei_sanwujiao_targets import DEFAULT_DB, _canonical_json, _write_json_atomic

MERGED_DATASET='加味三乌胶靶向库'
OLD_PREFIX='JWSWJ_加味三乌胶_'
SCHEMA='mssystem-jiawei-target-consolidation/1.0'

def _hash_rows(rows):
    h=hashlib.sha256()
    for row in rows:h.update(_canonical_json(dict(row) if isinstance(row,sqlite3.Row) else row).encode());h.update(b'\n')
    return h.hexdigest()

def _targets(c):return [dict(r) for r in c.execute('SELECT * FROM TargetCompounds ORDER BY id')]

def _scope(c,db,selected_ids):
    stat=db.stat()
    schema=[list(r) for r in c.execute("SELECT type,name,tbl_name,sql FROM sqlite_master WHERE name NOT LIKE 'sqlite_%' ORDER BY type,name")]
    tables={r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    targets=_targets(c)
    return dict(file_stat=dict(size_bytes=stat.st_size,mtime_ns=stat.st_mtime_ns),
        schema_sha256=_hash_rows(schema),schema_version=c.execute('PRAGMA schema_version').fetchone()[0],
        table_row_counts={t:c.execute('SELECT COUNT(*) FROM "'+t+'"').fetchone()[0] for t in ['Spectra','Compounds','Mapping','TargetCompounds','DiagnosticIons'] if t in tables},
        non_dataset_target_rows=sum(r['id'] not in selected_ids for r in targets),
        non_dataset_targets_sha256=_hash_rows([r for r in targets if r['id'] not in selected_ids]),
        diagnostics_sha256=_hash_rows(c.execute('SELECT * FROM DiagnosticIons ORDER BY id')))

def _expected(row):
    row=dict(row);old=row['dataset_source'];provenance=json.loads(row['provenance'])
    if not isinstance(provenance,dict):raise ValueError('target provenance must be an object')
    if 'original_dataset_source' in provenance and provenance['original_dataset_source']!=old:
        raise ValueError('conflicting original dataset provenance')
    provenance['original_dataset_source']=old
    row['dataset_source']=MERGED_DATASET;row['provenance']=_canonical_json(provenance)
    return row

def _authorizer(action,table,column,database,trigger):
    if action==sqlite3.SQLITE_UPDATE:
        return sqlite3.SQLITE_OK if table=='TargetCompounds' and column in {'dataset_source','provenance'} and trigger is None else sqlite3.SQLITE_DENY
    if action in {sqlite3.SQLITE_INSERT,sqlite3.SQLITE_DELETE,sqlite3.SQLITE_CREATE_TABLE,sqlite3.SQLITE_DROP_TABLE,sqlite3.SQLITE_ALTER_TABLE,sqlite3.SQLITE_CREATE_TRIGGER,sqlite3.SQLITE_DROP_TRIGGER,sqlite3.SQLITE_ATTACH,sqlite3.SQLITE_DETACH}:
        return sqlite3.SQLITE_DENY
    return sqlite3.SQLITE_OK

def consolidate(db:Path,preimage:Path,receipt_path:Path,*,apply=False):
    if apply and os.environ.get('MSSYSTEM_LIBRARY_WRITE_ENABLED')!='1':raise RuntimeError('apply requires MSSYSTEM_LIBRARY_WRITE_ENABLED=1')
    c=sqlite3.connect('file:'+str(db.resolve())+('?mode=rw' if apply else '?mode=ro'),uri=True,timeout=30);c.row_factory=sqlite3.Row
    try:
        if apply:c.execute('BEGIN IMMEDIATE')
        rows=[r for r in _targets(c) if str(r['dataset_source']).startswith(OLD_PREFIX)]
        if not rows:raise RuntimeError('no separate Jiawei datasets remain')
        if any(r['source_tier']!='literature_presence_only' for r in rows):raise RuntimeError('unexpected target evidence tier')
        ids={r['id'] for r in rows};before=_scope(c,db,ids);expected=[_expected(r) for r in rows]
        if c.execute('SELECT COUNT(*) FROM TargetCompounds WHERE dataset_source=?',(MERGED_DATASET,)).fetchone()[0]:raise RuntimeError('canonical dataset already exists; review before another consolidation')
        sources=sorted({r['dataset_source'] for r in rows})
        result=dict(schema_version=SCHEMA,status='preview',database_path=str(db.resolve()),dataset_source=MERGED_DATASET,
            source_datasets=sources,source_dataset_count=len(sources),merged_target_rows=len(rows),
            named_formula_entities=len({(r['compound_name'].strip().lower(),r['molecular_formula'] or '') for r in rows}),
            source_row_ids=[r['id'] for r in rows],pre_database_state=before,target_semantic_sha256=_hash_rows(expected),
            pre_target_semantic_sha256=_hash_rows(rows),unchanged_target_columns=[k for k in rows[0] if k not in {'dataset_source','provenance'}],
            inserted_target_rows=0,inserted_diagnostic_rows=0,deleted_target_rows=0,generated_at=datetime.now(timezone.utc).isoformat())
        if not apply:_write_json_atomic(receipt_path,result);return result
        preimage.parent.mkdir(parents=True,exist_ok=True)
        with preimage.open('xb') as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',mtime=0) as z:
                z.write(_canonical_json(dict(schema_version=SCHEMA,rows=rows,pre_database_state=before)).encode())
        result['preimage_path']=str(preimage.resolve());result['preimage_sha256']=hashlib.sha256(preimage.read_bytes()).hexdigest()
        c.set_authorizer(_authorizer)
        for row in expected:
            changed=c.execute('UPDATE TargetCompounds SET dataset_source=?,provenance=? WHERE id=? AND dataset_source=?',
                (MERGED_DATASET,row['provenance'],row['id'],json.loads(row['provenance'])['original_dataset_source'])).rowcount
            if changed!=1:raise RuntimeError('source target changed during consolidation')
        current=[r for r in _targets(c) if r['id'] in ids]
        if current!=expected:raise RuntimeError('target fields differ from the allowed consolidation')
        after=_scope(c,db,ids)
        for k in before:
            if k!='file_stat' and before[k]!=after[k]:raise RuntimeError('protected database content changed: '+k)
        c.commit();result.update(status='consolidated',post_database_state=after,updated_target_rows=len(rows))
        _write_json_atomic(receipt_path,result);return result
    except Exception:c.rollback();raise
    finally:c.close()

def verify_merge_scope(db:Path,preimage:Path,receipt):
    if receipt.get('schema_version')!=SCHEMA or receipt.get('status')!='consolidated':raise RuntimeError('consolidation receipt required')
    if str(db.resolve())!=receipt.get('database_path'):raise RuntimeError('consolidation database mismatch')
    if hashlib.sha256(preimage.read_bytes()).hexdigest()!=receipt.get('preimage_sha256'):raise RuntimeError('preimage hash mismatch')
    with gzip.open(preimage,'rt',encoding='utf-8') as f:original=json.load(f)
    rows=original['rows'];expected=[_expected(r) for r in rows]
    if original['pre_database_state']!=receipt['pre_database_state']:raise RuntimeError('pre-database proof mismatch')
    if _hash_rows(rows)!=receipt['pre_target_semantic_sha256'] or _hash_rows(expected)!=receipt['target_semantic_sha256']:raise RuntimeError('target consolidation proof mismatch')
    if len(rows)!=receipt['merged_target_rows'] or [r['id'] for r in rows]!=receipt['source_row_ids']:raise RuntimeError('row set mismatch')
    ids={r['id'] for r in rows}
    c=sqlite3.connect('file:'+str(db.resolve())+'?mode=ro',uri=True);c.row_factory=sqlite3.Row
    try:
        actual=[r for r in _targets(c) if r['dataset_source']==MERGED_DATASET]
        if actual!=expected:raise RuntimeError('merged data differ from their source rows')
        if any(r['dataset_source'].startswith(OLD_PREFIX) for r in _targets(c)):raise RuntimeError('unmerged Jiawei dataset remains')
        current=_scope(c,db,ids)
        for k,v in receipt['pre_database_state'].items():
            if k!='file_stat' and current[k]!=v:raise RuntimeError('protected database proof differs: '+k)
    finally:c.close()
    info=dict(dataset_source=MERGED_DATASET,dataset_version='consolidated-20261023',source_content_sha256=receipt['preimage_sha256'])
    return info,expected,current

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--db',type=Path,default=DEFAULT_DB)
    p.add_argument('--preimage',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);p.add_argument('--apply',action='store_true');a=p.parse_args()
    r=consolidate(a.db,a.preimage,a.receipt,apply=a.apply)
    print(json.dumps({k:r[k] for k in ['status','dataset_source','source_dataset_count','merged_target_rows','named_formula_entities']},ensure_ascii=False))
if __name__=='__main__':main()
