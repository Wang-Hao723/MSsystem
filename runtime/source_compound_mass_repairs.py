"""Source-literal mass/formula repairs shared by normalized and SQL retrieval.

The accepted base cache and source database are immutable. Only an explicitly
configured, fully bound generation can be applied; verification failures abort.
"""
from __future__ import annotations
from pathlib import Path
import hashlib,json,os,sqlite3,threading
from runtime.source_compound_mass_core import (
    apply_verified_mass_repairs,validate_registry,verify_source_evidence,digest,
)
ROOT=Path(__file__).resolve().parents[1]
SCHEMA='mssystem.source-compound-mass-runtime/1.0.0'
_LOCK=threading.RLock()
_STATES={}
SQL_FIELDS=('spectrum_id','source_file','precursor_mz','charge','adduct','ion_mode','source_accession','inchikey','smiles','peaks_json','molecular_formula','exact_mass','database_spectrum_id')

def _sha(path):
    with Path(path).open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()
def _path(root,value):
    p=Path(value).expanduser();return (p if p.is_absolute() else Path(root)/p).resolve()
def _signature(paths):
    result=[]
    for p in paths:
        try:
            s=Path(p).stat();result.append((str(p),s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,None if s.st_size==0 and str(p).endswith(('-wal','-journal')) else s.st_ctime_ns))
        except FileNotFoundError:result.append((str(p),'missing'))
    return tuple(result)
def _read_config(root):
    p=Path(root)/'config/source_compound_mass_repairs.json'
    if not p.exists():return None
    c=json.loads(p.read_text())
    if c.get('schema_version')!=SCHEMA or c.get('semantic_sha256')!=digest({k:v for k,v in c.items() if k!='semantic_sha256'}) or c.get('enabled') is not True or c.get('database_writes_allowed') is not False:
        raise ValueError('Invalid source compound mass configuration')
    return c

def source_compound_mass_dependency_paths(project_root=ROOT):
    root=Path(project_root).resolve();previous=_STATES.get(str(root))
    if previous is not None:return list(previous.dependencies)
    paths=[root/'config/source_compound_mass_repairs.json']
    try:
        c=_read_config(root)
        if c is not None:
            rp=_path(root,c['registry_path']);paths.extend([rp,_path(root,c['sql_bindings_path']),root/'config/library_record_store.json',root/'config/library_metadata_repairs.json',root/'runtime/source_compound_mass_core.py',root/'runtime/source_compound_mass_repairs.py'])
            r=json.loads(rp.read_text());paths.extend(Path(x) for x in r['source_files']);db=_path(root,c['source_database_path']);paths.extend([db,Path(str(db)+'-wal'),Path(str(db)+'-journal')])
            sc=json.loads((root/'config/library_record_store.json').read_text());cache=_path(root,os.environ.get('MSSYSTEM_LIBRARY_RECORD_STORE','').strip() or sc['cache_dir'])
            from matching.library_record_store import FILE_NAMES
            paths.extend(cache/n for n in (*FILE_NAMES,'manifest.json','acceptance.json'))
    except (OSError,ValueError,TypeError,KeyError):pass
    return list(dict.fromkeys(paths))
def source_compound_mass_dependency_signature(project_root=ROOT):
    return (os.environ.get('MSSYSTEM_LIBRARY_RECORD_STORE',''),_signature(source_compound_mass_dependency_paths(project_root)))
def _sql_values(raw):return {k:raw.get(k) for k in SQL_FIELDS}

class MassRepairState:
    def __init__(self,config,registry,sql_bindings,backing,raw_retriever,*,dependencies=(),verify_sources=True,metadata_indices=()):
        self.config=dict(config);self.registry=registry;self.dependencies=tuple(dependencies);self.signature=_signature(self.dependencies);self.override=os.environ.get('MSSYSTEM_LIBRARY_RECORD_STORE','');self.sql_rows={};self.normalized_applied=False;self.sql_applications=0
        compounds=validate_registry(registry);manifest=backing.manifest
        if registry['semantic_sha256']!=config['registry_semantic_sha256'] or len(registry['rows'])!=config['registered_reference_count'] or len(sql_bindings)!=len(registry['rows']):raise ValueError('Compound mass registry membership differs')
        if (manifest['normalized_records_sha256']!=registry['base_normalized_records_sha256'] or manifest['normalized_records_sha256']!=config['base_normalized_records_sha256'] or manifest['record_alignment_sha256']!=registry['base_record_alignment_sha256'] or manifest['record_alignment_sha256']!=config['base_record_alignment_sha256'] or manifest['record_count']!=registry['base_record_count'] or len(backing)!=registry['base_record_count'] or manifest['source_db']['sha256']!=registry['source_database_sha256'] or manifest['source_db']['sha256']!=config['source_database_sha256']):raise ValueError('Source compound mass base binding differs')
        if raw_retriever.library_snapshot_id!=config['base_library_snapshot_id'] or Path(raw_retriever.db_path).resolve()!=Path(config['source_database_path']).resolve():raise ValueError('Source compound mass SQL database differs')
        by_index={x['record_index']:x for x in registry['rows']}
        if set(by_index)&set(metadata_indices):raise ValueError('Overlapping spectral metadata and compound mass repairs require explicit joint review')
        by_pk={};seen=set()
        for b in sql_bindings:
            i=b['record_index'];r=by_index.get(i);raw=b['original_sql_values'];pk=raw['database_spectrum_id']
            if r is None or i in seen or isinstance(pk,bool) or not isinstance(pk,int) or pk<1 or pk in by_pk:raise ValueError('SQL compound mass membership differs')
            if (set(raw)!=set(SQL_FIELDS) or digest(raw)!=b['original_sql_values_sha256'] or raw['spectrum_id']!=r['spectrum_id'] or raw['inchikey']!=r['full_inchikey'] or b['database_compound_id']!=r['sql_compound_id'] or raw['source_file'] not in registry['source_files'] or raw['exact_mass']!=r['original_values']['exact_mass'] or b['original_record_sha256']!=r['original_record_sha256'] or b['updates']!={'exact_mass':r['updates']['exact_mass'],'molecular_formula':r['updates']['formula']}):raise ValueError('SQL compound mass binding differs from normalized registry')
            original=backing.read_uncached(i)
            if digest(original)!=r['original_record_sha256']:raise ValueError('Original compound mass normalized record differs')
            seen.add(i);by_pk[pk]=b
        if seen!=set(by_index):raise ValueError('SQL compound mass binding incomplete')
        if verify_sources:verify_source_evidence(registry,compounds)
        # Use the original verified Mapping join, removing only the precursor
        # index hint for this complete exact-primary-key validation.
        strategy=raw_retriever._sqlite_query_strategy;from_clause=strategy['from_clause']
        if strategy.get('precursor_index'):from_clause=from_clause.replace(' INDEXED BY '+raw_retriever._quoted_identifier(strategy['precursor_index']),'')
        select='SELECT s.id AS database_spectrum_id,s.spectrum_id,s.source_file,s.precursor_mz,s.charge,s.adduct,s.ion_mode,s.source_accession,s.peaks_json,c.inchikey,c.smiles,c.molecular_formula,c.exact_mass,c.id AS database_compound_id,m.spectrum_id AS mapping_spectrum_id FROM '+from_clause+' WHERE s.id IN (SELECT value FROM json_each(?))'
        conn=raw_retriever._connect()
        try:
            observed=conn.execute(select,(json.dumps(sorted(by_pk)),)).fetchall()
        finally:conn.close()
        if len(observed)!=len(by_pk):raise ValueError('SQL compound mass join is not one-to-one')
        from rdkit import Chem
        for item in observed:
            raw=dict(item);pk=raw['database_spectrum_id'];b=by_pk.get(pk)
            if b is None or pk in self.sql_rows or raw['mapping_spectrum_id']!=pk or raw['database_compound_id']!=b['database_compound_id'] or _sql_values(raw)!=b['original_sql_values']:raise ValueError('Original SQL compound mass record differs')
            mol=Chem.MolFromSmiles(raw.get('smiles') or '')
            if mol is None or Chem.MolToInchiKey(mol)!=raw['inchikey']:raise ValueError('Original SQL full structure identity differs')
            self.sql_rows[pk]=b
        self.assert_fresh()
    def assert_fresh(self):
        if self.override!=os.environ.get('MSSYSTEM_LIBRARY_RECORD_STORE','') or self.signature!=_signature(self.dependencies):raise ValueError('Source compound mass dependencies changed; restart and reverify')
    def apply_sql(self,row):
        """Check a standalone row before and after applying its bound values."""
        self.assert_fresh()
        result=self._apply_sql_checked_row(row)
        self.assert_fresh()
        return result
    def _apply_sql_checked_row(self,row):
        # Only called inside an entry/exit generation-checked operation.
        # Full original values remain checked for every repaired row.
        raw=dict(row);binding=self.sql_rows.get(raw.get('database_spectrum_id'))
        if binding is None:return raw
        if _sql_values(raw)!=binding['original_sql_values']:raise ValueError('Original SQL candidate mass fields changed after verification')
        raw.update(binding['updates']);raw['source_compound_mass_repair']={'registry_semantic_sha256':self.registry['semantic_sha256'],'sql_bindings_sha256':self.config['sql_bindings_sha256'],'record_index':binding['record_index'],'original_exact_mass':binding['original_sql_values']['exact_mass'],'original_molecular_formula':binding['original_sql_values']['molecular_formula'],'source_database_written':False}
        self.sql_applications+=1;return raw
    def status(self):
        return {'status':'applied' if self.normalized_applied else 'bindings_verified','registered_reference_count':len(self.registry['rows']),'complete_sql_and_original_cache_verified_count':len(self.sql_rows),'normalized_overlay_applied':self.normalized_applied,'sql_candidate_applications':self.sql_applications,'registry_semantic_sha256':self.registry['semantic_sha256'],'config_semantic_sha256':self.config['semantic_sha256'],'source_database_written':False,'precursor_index_changed':False}

def get_source_compound_mass_state(project_root=ROOT,*,backing=None,raw_retriever=None):
    root=Path(project_root).resolve()
    with _LOCK:
        old=_STATES.get(str(root))
        if old is not None:old.assert_fresh();return old
        config=_read_config(root)
        if config is None:return None
        paths=source_compound_mass_dependency_paths(root);signature=_signature(paths);registry_path=_path(root,config['registry_path']);sql_path=_path(root,config['sql_bindings_path'])
        if _sha(registry_path)!=config['registry_sha256'] or _sha(sql_path)!=config['sql_bindings_sha256'] or sql_path.stat().st_size!=config['sql_bindings_size_bytes']:raise ValueError('Source compound mass registry files differ')
        registry=json.loads(registry_path.read_text());bindings=list(map(json.loads,sql_path.open()));owned=False
        try:
            if backing is None:
                from matching.library_record_store import processing_identity
                from runtime.library_cache_identity import load_compatible_record_store
                sc=json.loads((root/'config/library_record_store.json').read_text());override=os.environ.get('MSSYSTEM_LIBRARY_RECORD_STORE','').strip()
                if override.lower() in {'off','0','false'} or not override and sc.get('enabled') is not True:raise ValueError('Source compound mass repairs require accepted base cache')
                backing,mz,nm=load_compatible_record_store(_path(root,override or sc['cache_dir']),source_db=_path(root,config['source_database_path']),current_identity=processing_identity(),cache_bytes=0,project_root=root);owned=True
            if raw_retriever is None:
                from matching.evirank_retrieval import SQLiteCandidateRecallRetriever
                raw_retriever=SQLiteCandidateRecallRetriever(_path(root,config['source_database_path']),library_snapshot_id=config['base_library_snapshot_id'])
            metadata_path=root/'config/library_metadata_repairs.json';indices=[]
            if metadata_path.is_file():
                import ijson
                with metadata_path.open('rb') as handle:indices=[r['record_index'] for r in ijson.items(handle,'rows.item')]
            state=MassRepairState(config,registry,bindings,backing,raw_retriever,dependencies=paths,metadata_indices=indices)
        finally:
            if owned:backing.close();del mz,nm
        if signature!=_signature(paths):raise ValueError('Source compound mass dependencies changed during verification')
        _STATES[str(root)]=state;return state

def apply_source_compound_mass_records(records,neutral_index,*,project_root=ROOT,bucket_size=None):
    state=get_source_compound_mass_state(project_root,backing=records)
    if state is None:return records,neutral_index,None
    state.assert_fresh()
    wrapped,new,identity=apply_verified_mass_repairs(records,neutral_index,state.registry,bucket_size=bucket_size,verify_sources=False)
    return wrapped,new,identity

def mark_source_compound_mass_published(records,*,project_root=ROOT):
    state=get_source_compound_mass_state(project_root)
    if state is None:return
    state.assert_fresh()
    if not source_compound_mass_status(records,True,project_root=project_root).get('normalized_overlay_applied'):
        raise ValueError('Compound mass records were not prepared for publication')
    state.normalized_applied=True

def source_compound_mass_snapshot_id(base_id,project_root=ROOT):
    config=_read_config(Path(project_root).resolve())
    if config is None:return str(base_id)
    state=_STATES.get(str(Path(project_root).resolve()))
    if state is not None:state.assert_fresh()
    return str(base_id)+':compound-mass:'+config['registry_semantic_sha256']+':sql:'+config['sql_bindings_sha256']

def source_compound_mass_preflight(project_root=ROOT):
    try:
        state=get_source_compound_mass_state(project_root)
        return {'compatible':True,'status':'not_configured'} if state is None else {'compatible':True,**state.status()}
    except (OSError,ValueError,KeyError,TypeError,RuntimeError,sqlite3.Error) as exc:return {'compatible':False,'status':'invalid_binding','error':str(exc)}

def source_compound_mass_status(records=None,loaded=False,*,project_root=ROOT):
    state=_STATES.get(str(Path(project_root).resolve()))
    if state is None:return {'status':'not_initialized','normalized_overlay_applied':False}
    try:state.assert_fresh()
    except ValueError as exc:return {'status':'invalid_binding','error':str(exc)}
    result=state.status()
    if not loaded:
        result['normalized_overlay_applied']=False
        result['status']='bindings_verified'
    if loaded:
        cursor=records;found=False
        for _ in range(4):
            if cursor is None:break
            if getattr(cursor,'overlay_identity',{}).get('registry_semantic_sha256')==state.registry['semantic_sha256']:found=True;break
            cursor=getattr(cursor,'backing',None)
        result['normalized_overlay_applied']=found
        if not found:result['status']='not_applied'
    return result
