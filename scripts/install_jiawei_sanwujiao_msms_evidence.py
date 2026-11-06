#!/usr/bin/env python3
"""Install reviewed source-specific MS/MS targets and diagnostic ions only.

No Spectra, Compounds, Mapping, schema, or existing target rows are writable.
Selected fragment lists retain their original MSn/adduct/source context in
TargetCompounds.provenance; missing spectral intensities are never invented.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, re, sqlite3, sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

PROJECT_ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT_ROOT))
from scripts.install_jiawei_sanwujiao_targets import (
    TARGET_COLUMNS, DEFAULT_DB, _canonical_json, _sha256_json, _write_json_atomic,
    _database_scope_snapshot, _existing_dataset, _semantic_projection,
    _structure_fields, _target_only_authorizer,
)
from matching.expected_compounds import _normalise_adduct
from matching.library_quality import formula_exact_mass, normalize_formula_text
from matching.library_record_store import file_sha256
from utils.adduct import ion_mz

CATALOG_SCHEMA='mssystem-jiawei-sanwujiao-msms-catalog/1.0'
RECEIPT_SCHEMA='mssystem-jiawei-msms-target-install/1.0'
DEFAULT_ROOT=PROJECT_ROOT.parent/'research/20261020_jiawei_sanwujiao_targeted/msms_v3'
ELECTRON=0.000548579909065

def _original_precursor_ms2(event,observed):
    if int(event['ms_level'])!=2:return False
    isolated=float(event['isolated_precursor_mz'])
    path=event.get('source_precursor_path')
    if path is not None and path!=[isolated]:return False
    if abs(isolated-observed)<=.55:return True
    # Table S9 prints truncated integer isolation labels alongside an exact-mass
    # precursor in the same row. Preserve both source values and their binding.
    if event.get('precursor_label_binding')!='nominal_integer_truncation_same_table_row':return False
    return (path==[isolated] and isolated.is_integer() and isolated<=observed<isolated+1
            and event.get('source_row_precursor_mz')==observed)

def load_msms_catalog(path:Path):
    payload=json.loads(path.read_text())
    if payload.get('schema_version')!=CATALOG_SCHEMA:raise ValueError('invalid MSMS catalogue schema')
    dataset=str(payload.get('dataset_source') or '')
    version_match=re.fullmatch(r'JWSWJ_加味三乌胶_.+_v([3-9]|10)',dataset)
    if not version_match:
        raise ValueError('catalogue is not a supported reviewed Jiawei v3/v4/v5/v6/v7/v8/v9/v10 dataset')
    major=version_match[1]
    if payload.get('dataset_version')!=major+'.0.0':
        raise ValueError('catalogue dataset/version mismatch')
    if payload.get('plant_source')!='加味三乌胶':raise ValueError('wrong formula source')
    files={str(Path(x['path']).resolve()):x['sha256'] for x in payload.get('source_files',[])}
    input_files={str(Path(x['path']).resolve()):x['sha256'] for x in payload.get('source_inputs',[])}
    for name,digest in {**files,**input_files}.items():
        if not re.fullmatch('[0-9a-f]{64}',digest) or file_sha256(name)!=digest:
            raise ValueError('source file hash mismatch: '+name)
    records=payload.get('rows') or []
    if not records:raise ValueError('empty catalogue')
    seen=set();rows=[];diagnostics={}
    for item in records:
        target=dict(item['target']);tid=target.get('target_id')
        if not tid or not str(tid).startswith('JWSWJ'+major+'-') or tid in seen:raise ValueError('duplicate or invalid target_id')
        seen.add(tid)
        if set(target)!=set(TARGET_COLUMNS):raise ValueError('unexpected target fields')
        if target['dataset_source']!=dataset or target['dataset_version']!=payload['dataset_version']:
            raise ValueError('target dataset differs')
        if target['source_tier']!='literature_presence_only':raise ValueError('unreviewed evidence tier')
        if not target['compound_name'] or not target['plant_source']:raise ValueError('missing compound/material name')
        formula=normalize_formula_text(target['molecular_formula']);mass=formula_exact_mass(formula)
        if not mass or abs(float(target['exact_mass'])-mass)/mass*1e6>5:raise ValueError('formula mass mismatch')
        flags=json.loads(target['quality_flags_json'])
        zs={int(f.split(':')[1]) for f in flags if str(f).startswith('formal_charge:')}
        if len(zs)!=1:raise ValueError('explicit intrinsic formal charge required')
        z=next(iter(zs))
        if z not in [-1,0,1]:raise ValueError('unsupported intrinsic charge')
        _structure_fields({**target,'structure_identity_verified':bool(target['structure_identity_verified'])},formula,z,where=tid)
        ad=_normalise_adduct(target['adduct'])
        if not ad or ad!=target['adduct']:raise ValueError('invalid adduct')
        if int(target['charge'])!=(1 if ad.endswith('+') else -1):raise ValueError('adduct/charge mismatch')
        if z and ad!=('[M]+' if z>0 else '[M]-'):raise ValueError('intrinsic ion adduct mismatch')
        expected=mass-z*ELECTRON if z else ion_mz(mass,ad)
        if expected is None or abs(float(target['precursor_mz'])-expected)>1e-6:raise ValueError('calculated precursor mismatch')
        ions=[float(v) for v in item.get('diagnostic_ions',[])]
        if ions!=sorted(set(ions)) or any(not math.isfinite(v) or v<=10 or v>=expected-.025 for v in ions):
            raise ValueError('invalid product/precursor separation')
        provenance=json.loads(target['provenance']);e=provenance.get('msms_evidence')
        if ions:
            if not e or e.get('runtime_diagnostic_mz')!=ions:raise ValueError('diagnostics lack source-specific evidence')
            if e.get('source_fragment_stage_disclosed') is False or e.get('source_precursor_isolation_disclosed') is False:
                raise ValueError('unreported source MS stage or precursor isolation cannot supply runtime diagnostics')
            if e.get('missing_intensities_filled') is not False:raise ValueError('invented intensities forbidden')
            if e.get('product_adduct_unambiguous') is not True or e.get('adduct')!=ad:
                raise ValueError('ambiguous product-ion/adduct assignment')
            observed=float(e['reported_precursor_mz'])
            if abs(observed-expected)/expected*1e6>20:raise ValueError('source precursor mass conflict')
            allowed=set()
            for event in e.get('reported_events',[]):
                if not _original_precursor_ms2(event,observed):continue
                for p in event['products']:
                    text=str(p['mz_text'])
                    if '.' in text and len(text.split('.')[1])>=2:allowed.add(float(p['mz']))
            if not set(ions)<=allowed:raise ValueError('MS3, wrong-precursor, or nominal-only ion leaked into runtime diagnostics')
        if e:
            p=str(Path(e['source_path']).resolve())
            if files.get(p)!=e['source_sha256'] or target['source_content_sha256']!=e['source_sha256']:
                raise ValueError('source-specific target is not hash-bound')
            source_adduct=e.get('adduct_as_reported',ad)
            if target['reported_precursor_mz']!=e['reported_precursor_mz'] or target['reported_adduct']!=source_adduct:
                raise ValueError('reported precursor/adduct erased or changed')
        diagnostics[tid]=ions;rows.append(target)
    info={'dataset_source':dataset,'dataset_version':payload['dataset_version'],
          'source_content_sha256':file_sha256(path),'target_row_count':len(rows),
          'compound_count':len({(r['compound_name'],r['molecular_formula']) for r in rows}),
          'diagnostic_ion_rows':sum(map(len,diagnostics.values()))}
    return info,rows,diagnostics

def _diagnostic_projection(c,dataset):
    return [list(r) for r in c.execute(
        'SELECT t.target_id,d.fragment_mz,d.ordinal FROM DiagnosticIons d JOIN TargetCompounds t ON t.id=d.target_id '
        'WHERE t.dataset_source=? ORDER BY t.target_id,d.ordinal,d.id',(dataset,))]

def _other_diagnostics(c,dataset):
    h=hashlib.sha256();count=0
    for row in c.execute('SELECT d.id,d.target_id,d.fragment_mz,d.ordinal FROM DiagnosticIons d '
        'LEFT JOIN TargetCompounds t ON t.id=d.target_id WHERE t.dataset_source IS NULL OR t.dataset_source<>? ORDER BY d.id',(dataset,)):
        h.update(_canonical_json(list(row)).encode());h.update(b'\n');count+=1
    return {'rows':count,'sha256':h.hexdigest()}

def scope_snapshot(c,path,dataset):
    r=_database_scope_snapshot(c,path,dataset)
    r['non_dataset_diagnostics']=_other_diagnostics(c,dataset)
    return r

def _expected_diagnostic_projection(diagnostics):
    return [[tid,mz,n] for tid,ions in sorted(diagnostics.items()) for n,mz in enumerate(ions,1)]

def _check_scope(before,after,rows,diagnostics):
    for k in ['schema_sha256','schema_version','non_dataset_target_rows','non_dataset_targets_sha256','non_dataset_diagnostics']:
        if before[k]!=after[k]:raise RuntimeError('unrelated database content changed: '+k)
    b=before['table_row_counts'];a=after['table_row_counts']
    for k in ['Spectra','Compounds','Mapping']:
        if b[k]!=a[k]:raise RuntimeError('spectral table changed: '+k)
    if a['TargetCompounds']!=b['TargetCompounds']+len(rows):raise RuntimeError('unexpected target increment')
    if a['DiagnosticIons']!=b['DiagnosticIons']+sum(map(len,diagnostics.values())):raise RuntimeError('unexpected diagnostic increment')

def msms_authorizer(action,arg1,arg2,db,trigger):
    if action==sqlite3.SQLITE_INSERT and arg1=='DiagnosticIons':return sqlite3.SQLITE_OK
    # The inherited authorizer disallows schema mutations and every other table.
    if action in [sqlite3.SQLITE_UPDATE,sqlite3.SQLITE_DELETE]:return sqlite3.SQLITE_DENY
    return _target_only_authorizer(action,arg1,arg2,db,trigger)

def install_catalog(catalog_path:Path,db_path:Path,receipt_path:Path,*,apply=False):
    if apply and os.environ.get('MSSYSTEM_LIBRARY_WRITE_ENABLED','').lower() not in ['1','true','yes','on']:
        raise RuntimeError('apply requires MSSYSTEM_LIBRARY_WRITE_ENABLED=1')
    info,rows,diagnostics=load_msms_catalog(catalog_path)
    c=sqlite3.connect(str(db_path) if apply else 'file:'+str(db_path.resolve())+'?mode=ro',uri=not apply,timeout=20)
    c.row_factory=sqlite3.Row
    try:
        before=scope_snapshot(c,db_path,info['dataset_source'])
        existing=_existing_dataset(c,info['dataset_source'])
        if existing:
            if _semantic_projection(existing)!=_semantic_projection(rows) or _diagnostic_projection(c,info['dataset_source'])!=_expected_diagnostic_projection(diagnostics):
                raise RuntimeError('existing dataset has different content')
            status='already_installed_exact';inserted=0;inserted_d=0;after=before
        elif not apply:
            status='preview';inserted=0;inserted_d=0;after=before
        else:
            c.set_authorizer(msms_authorizer);c.execute('BEGIN IMMEDIATE')
            for row in rows:
                values=dict(row);values['created_at']=datetime.now(timezone.utc).isoformat()
                cur=c.execute('INSERT INTO TargetCompounds ('+','.join(TARGET_COLUMNS)+') VALUES ('+','.join('?' for _ in TARGET_COLUMNS)+')',tuple(values[k] for k in TARGET_COLUMNS))
                for n,mz in enumerate(diagnostics[row['target_id']],1):
                    c.execute('INSERT INTO DiagnosticIons(target_id,fragment_mz,ordinal) VALUES(?,?,?)',(cur.lastrowid,mz,n))
            after=scope_snapshot(c,db_path,info['dataset_source'])
            _check_scope(before,after,rows,diagnostics)
            if _semantic_projection(_existing_dataset(c,info['dataset_source']))!=_semantic_projection(rows):raise RuntimeError('inserted targets differ')
            if _diagnostic_projection(c,info['dataset_source'])!=_expected_diagnostic_projection(diagnostics):raise RuntimeError('inserted diagnostics differ')
            c.commit();status='installed';inserted=len(rows);inserted_d=sum(map(len,diagnostics.values()))
        result={'schema_version':RECEIPT_SCHEMA,'status':status,'generated_at':datetime.now(timezone.utc).isoformat(),
            **info,'catalog_path':str(catalog_path.resolve()),'database_path':str(db_path.resolve()),
            'target_semantic_sha256':_sha256_json(_semantic_projection(rows)),
            'diagnostic_semantic_sha256':_sha256_json(_expected_diagnostic_projection(diagnostics)),
            'inserted_target_rows':inserted,'inserted_diagnostic_rows':inserted_d,
            'pre_database_state':before,'post_database_state':after,
            'allowed_write_tables':['TargetCompounds','DiagnosticIons'],'spectra_compounds_mapping_written':False}
        _write_json_atomic(receipt_path,result)
        return result
    except Exception:
        c.rollback();raise
    finally:c.close()

def verify_install_scope(db_path:Path,catalog_path:Path,receipt:Mapping):
    info,rows,diagnostics=load_msms_catalog(catalog_path)
    if receipt.get('schema_version')!=RECEIPT_SCHEMA or receipt.get('status')!='installed':raise RuntimeError('MSMS installer receipt required')
    for key in ['dataset_source','dataset_version','source_content_sha256']:
        if receipt.get(key)!=info[key]:raise RuntimeError('receipt/catalogue mismatch: '+key)
    if receipt['target_semantic_sha256']!=_sha256_json(_semantic_projection(rows)):raise RuntimeError('target proof mismatch')
    if receipt['diagnostic_semantic_sha256']!=_sha256_json(_expected_diagnostic_projection(diagnostics)):raise RuntimeError('diagnostic proof mismatch')
    c=sqlite3.connect('file:'+str(db_path.resolve())+'?mode=ro',uri=True);c.row_factory=sqlite3.Row
    try:
        current=scope_snapshot(c,db_path,info['dataset_source'])
        _check_scope(receipt['pre_database_state'],current,rows,diagnostics)
        if _semantic_projection(_existing_dataset(c,info['dataset_source']))!=_semantic_projection(rows):raise RuntimeError('installed targets changed')
        if _diagnostic_projection(c,info['dataset_source'])!=_expected_diagnostic_projection(diagnostics):raise RuntimeError('installed diagnostics changed')
    finally:c.close()
    return info,rows,current

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog',type=Path,default=DEFAULT_ROOT/'catalog_v3.json')
    p.add_argument('--db',type=Path,default=DEFAULT_DB)
    p.add_argument('--receipt',type=Path,default=DEFAULT_ROOT/'INSTALL_RECEIPT_V3.json')
    p.add_argument('--apply',action='store_true');a=p.parse_args()
    r=install_catalog(a.catalog.resolve(),a.db.resolve(),a.receipt.resolve(),apply=a.apply)
    print(json.dumps({k:r[k] for k in ['status','compound_count','target_row_count','diagnostic_ion_rows','inserted_target_rows','inserted_diagnostic_rows']},ensure_ascii=False))

if __name__=='__main__':main()
