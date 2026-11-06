from __future__ import annotations
import hashlib, importlib.util, json, sqlite3
from pathlib import Path
import pytest
from matching.expected_compounds import _associate_ms2
from matching.library_quality import formula_exact_mass
from utils.adduct import ion_mz
from scripts import install_jiawei_sanwujiao_msms_evidence as installer
from scripts.install_jiawei_sanwujiao_targets import TARGET_COLUMNS

def _database(p):
    text={'target_id','source_file','source_content_sha256','dataset_version','dataset_source','compound_name',
          'molecular_formula','smiles','inchi','inchikey','reported_inchikey','structure_identity_source',
          'adduct','reported_adduct','reported_retention_time_unit','reported_collision_energy_raw','instrument',
          'target_searchability','plant_source','provenance','doi','dataset_accession','spectrum_accession',
          'source_tier','quality_flags_json','created_at'}
    integers={'charge','structure_identity_verified','source_row_index'}
    defs=['id INTEGER PRIMARY KEY']+[k+(' TEXT' if k in text else ' INTEGER' if k in integers else ' REAL') for k in TARGET_COLUMNS]
    with sqlite3.connect(p) as c:
        c.execute('CREATE TABLE TargetCompounds('+','.join(defs)+')')
        c.execute('CREATE TABLE DiagnosticIons(id INTEGER PRIMARY KEY,target_id INTEGER,fragment_mz REAL,ordinal INTEGER)')
        for n in ['Spectra','Compounds','Mapping']:
            c.execute('CREATE TABLE '+n+'(id INTEGER PRIMARY KEY,marker TEXT)')
            c.execute('INSERT INTO '+n+' VALUES(1,?)',('preserved',))

def _catalog(tmp):
    source=tmp/'source.txt';source.write_text('Isolated precursor 195.0877; MS2 product ions 138.0650, 164.0450. No intensities reported.')
    digest=hashlib.sha256(source.read_bytes()).hexdigest();mass=formula_exact_mass('C8H10N4O2');mz=ion_mz(mass,'[M+H]+')
    evidence={'adduct':'[M+H]+','reported_precursor_mz':195.0877,'source_path':str(source),'source_sha256':digest,
        'runtime_diagnostic_mz':[138.065,164.045],'product_adduct_unambiguous':True,'missing_intensities_filled':False,
        'reported_events':[{'ms_level':2,'isolated_precursor_mz':195.0877,'products':[
            {'mz':138.065,'mz_text':'138.0650','relative_intensity_reported':None},
            {'mz':164.045,'mz_text':'164.0450','relative_intensity_reported':None}]}]}
    t={k:None for k in TARGET_COLUMNS}
    t.update(target_id='JWSWJ3-test-1',dataset_source='JWSWJ_加味三乌胶_test_v3',dataset_version='3.0.0',
        compound_name='Caffeine',plant_source='test material',molecular_formula='C8H10N4O2',exact_mass=mass,
        precursor_mz=mz,adduct='[M+H]+',charge=1,structure_identity_verified=0,
        source_file=str(source),source_content_sha256=digest,source_tier='literature_presence_only',
        quality_flags_json=json.dumps(['formal_charge:0','literature_presence_only','calculated_ion_hypothesis']),
        reported_adduct='[M+H]+',reported_precursor_mz=195.0877,provenance=json.dumps({'msms_evidence':evidence}))
    return {'schema_version':installer.CATALOG_SCHEMA,'dataset_source':t['dataset_source'],'dataset_version':'3.0.0',
        'plant_source':'加味三乌胶','source_files':[{'path':str(source),'sha256':digest}],
        'rows':[{'target':t,'diagnostic_ions':[138.065,164.045]}]}

def test_msms_install_preserves_all_spectral_tables_and_is_idempotent(tmp_path,monkeypatch):
    db=tmp_path/'library.db';_database(db);payload=_catalog(tmp_path)
    catalog=tmp_path/'catalog.json';catalog.write_text(json.dumps(payload));receipt=tmp_path/'receipt.json'
    before=hashlib.sha256(db.read_bytes()).hexdigest()
    preview=installer.install_catalog(catalog,db,receipt,apply=False)
    assert preview['status']=='preview' and hashlib.sha256(db.read_bytes()).hexdigest()==before
    monkeypatch.setenv('MSSYSTEM_LIBRARY_WRITE_ENABLED','1')
    result=installer.install_catalog(catalog,db,receipt,apply=True)
    assert result['inserted_target_rows']==1 and result['inserted_diagnostic_rows']==2
    installer.verify_install_scope(db,catalog,result)
    assert installer.install_catalog(catalog,db,receipt,apply=True)['status']=='already_installed_exact'
    with sqlite3.connect(db) as c:
        assert c.execute('SELECT COUNT(*) FROM DiagnosticIons').fetchone()[0]==2
        for n in ['Spectra','Compounds','Mapping']:
            assert c.execute('SELECT marker FROM '+n).fetchall()==[('preserved',)]

def test_unknown_source_content_and_forbidden_ms3_products_fail_before_write(tmp_path):
    payload=_catalog(tmp_path);catalog=tmp_path/'catalog.json'
    evidence=json.loads(payload['rows'][0]['target']['provenance'])['msms_evidence']
    evidence['reported_events'][0]['ms_level']=3
    payload['rows'][0]['target']['provenance']=json.dumps({'msms_evidence':evidence})
    catalog.write_text(json.dumps(payload))
    with pytest.raises(ValueError,match='MS3'):
        installer.load_msms_catalog(catalog)
    payload=_catalog(tmp_path);catalog.write_text(json.dumps(payload))
    Path(payload['source_files'][0]['path']).write_text('changed source')
    with pytest.raises(ValueError,match='hash mismatch'):
        installer.load_msms_catalog(catalog)

@pytest.mark.parametrize('major',[4,5,6,7,8,9,10])
def test_reviewed_increments_are_version_bound_and_future_versions_rejected(tmp_path,major):
    payload=_catalog(tmp_path)
    payload.update(dataset_source=f'JWSWJ_加味三乌胶_原文补充_test_v{major}',dataset_version=f'{major}.0.0')
    row=payload['rows'][0]['target']
    row.update(target_id=f'JWSWJ{major}-test-1',dataset_source=payload['dataset_source'],dataset_version=f'{major}.0.0')
    catalog=tmp_path/f'v{major}.json';catalog.write_text(json.dumps(payload))
    assert installer.load_msms_catalog(catalog)[0]['target_row_count']==1
    payload['dataset_version']='3.0.0';catalog.write_text(json.dumps(payload))
    with pytest.raises(ValueError,match='version mismatch'):installer.load_msms_catalog(catalog)
    payload['dataset_source']=payload['dataset_source'].removesuffix(f'_v{major}')+'_v11'
    catalog.write_text(json.dumps(payload))
    with pytest.raises(ValueError,match='supported reviewed'):installer.load_msms_catalog(catalog)

@pytest.mark.parametrize('field',['source_fragment_stage_disclosed','source_precursor_isolation_disclosed'])
def test_unassigned_source_stage_or_isolation_cannot_be_promoted_to_runtime_diagnostics(tmp_path,field):
    payload=_catalog(tmp_path);target=payload['rows'][0]['target']
    provenance=json.loads(target['provenance']);provenance['msms_evidence'][field]=False
    target['provenance']=json.dumps(provenance);catalog=tmp_path/'unstaged.json';catalog.write_text(json.dumps(payload))
    with pytest.raises(ValueError,match='unreported source MS stage'):installer.load_msms_catalog(catalog)
    # The exact source list remains inspectable when no runtime ion is emitted.
    payload['rows'][0]['diagnostic_ions']=[]
    provenance['msms_evidence']['runtime_diagnostic_mz']=[];target['provenance']=json.dumps(provenance)
    catalog.write_text(json.dumps(payload));assert installer.load_msms_catalog(catalog)[0]['diagnostic_ion_rows']==0

def test_msms_authorizer_rejects_spectral_and_schema_mutations(tmp_path):
    db=tmp_path/'db.sqlite';_database(db)
    with sqlite3.connect(db) as c:
        c.set_authorizer(installer.msms_authorizer)
        for sql in ['DELETE FROM TargetCompounds','UPDATE Spectra SET marker="changed"',
                    'INSERT INTO Spectra(marker) VALUES("bad")','CREATE TABLE bad(id INT)']:
            with pytest.raises(sqlite3.DatabaseError):c.execute(sql)

def _associate(adduct,members):
    mass=formula_exact_mass('C6H12O6');mz=ion_mz(mass,adduct)
    feature={'expected_ion_id':'g1','boundary_start_rt_minutes':0.0,'boundary_end_rt_minutes':1.0,
             'apex_rt_minutes':.5,'fwhm_seconds':10,'target_mz':mz}
    ion={'polarity':'POS','adduct':adduct,'members':members}
    row={'spectrum_id':'synthetic','precursor_mz':mz,'retention_time':.5,'retention_time_unit':'minutes',
         'ion_mode':'POS','query_peaks':[[89.024,100],[119.035,50]],'match':None}
    _associate_ms2([feature],{'g1':ion},[row]);return feature

def test_reported_fragments_do_not_transfer_to_generated_adducts():
    member={'target_id':'source1','compound_name':'isomer candidates','active_adduct':'[M+H]+',
            'diagnostic_ions':[89.024,119.035]}
    assert len(_associate('[M+H]+',[member])['ms2_diagnostic_support'])==1
    assert _associate('[M+Na]+',[member])['ms2_diagnostic_support']==[]

def test_one_mrm_or_separate_one_ion_sources_do_not_form_two_ion_rule():
    members=[{'target_id':'article-A','active_adduct':'[M+H]+','diagnostic_ions':[89.024]},
             {'target_id':'article-B','active_adduct':'[M+H]+','diagnostic_ions':[119.035]}]
    assert _associate('[M+H]+',members)['ms2_diagnostic_support']==[]

@pytest.fixture(scope='module')
def parser():
    p=Path(__file__).resolve().parents[2]/'research/20261020_jiawei_sanwujiao_targeted/msms_v3/build_msms_catalog.py'
    spec=importlib.util.spec_from_file_location('jiawei_msms_source_builder',p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.parse_product_events

def test_source_parser_retains_ms3_as_separate_events(parser):
    events=parser('MS2[632]: 572.28, 540.22 MS3[572]: 512.20, 480.19',632.3065)
    assert [x['is_original_precursor_ms2'] for x in events]==[True,False]
    assert [p['mz'] for p in events[0]['products']]==[572.28,540.22]
    assert [p['mz'] for p in events[1]['products']]==[512.20,480.19]

def test_ppm_annotations_are_not_product_masses_or_intensities(parser):
    event=parser('115.0033 [M-H-H2O]− (−3.4)',133.0140)[0]
    assert event['products']==[{'mz':115.0033,'mz_text':'115.0033','relative_intensity_reported':None}]

def test_explicit_major_peak_intensities_are_retained_without_weak_peak_invention(parser):
    event=parser('336.4 (21), 188.3 (100)',354.1,relative_intensity=True)[0]
    assert [(p['mz'],p['relative_intensity_reported']) for p in event['products']]==[(336.4,21.0),(188.3,100.0)]

@pytest.fixture(scope='module')
def v4_parser():
    p=Path(__file__).resolve().parents[2]/'research/20261020_jiawei_sanwujiao_targeted/fulltext_v4/build_incremental_catalog_v4.py'
    spec=importlib.util.spec_from_file_location('jiawei_v4_incremental_source_builder',p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.parse_events

def test_ms3_and_ms4_isolation_paths_cannot_supply_ms2_diagnostics(v4_parser):
    events=v4_parser('MS2 [410]: 392.2435 MS3 [410→392]: 374.2330 MS4 [410→392→374]: 356.2225',410.2519)
    assert [e['source_precursor_path'] for e in events]==[[410.0],[410.0,392.0],[410.0,392.0,374.0]]
    assert [e['isolated_precursor_mz'] for e in events]==[410.0,392.0,374.0]
    assert [e['is_original_precursor_ms2'] for e in events]==[True,False,False]
    assert [p['mz'] for e in events if installer._original_precursor_ms2(e,410.2519) for p in e['products']]==[392.2435]

def test_truncated_nominal_precursor_requires_same_row_binding(v4_parser):
    ms2,ms3=v4_parser('MS2 [902]: 586.3010 MS3 [902→586]: 568.2908',902.5591)
    assert ms2['isolated_precursor_mz']==902.0
    assert ms2['source_row_precursor_mz']==902.5591
    assert installer._original_precursor_ms2(ms2,902.5591)
    assert not installer._original_precursor_ms2(ms3,902.5591)
    assert not installer._original_precursor_ms2({k:v for k,v in ms2.items() if k!='source_row_precursor_mz'},902.5591)
    assert not installer._original_precursor_ms2(ms2,903.5591)
    assert not installer._original_precursor_ms2({**ms2,'source_precursor_path':[902.0,586.0]},902.5591)
