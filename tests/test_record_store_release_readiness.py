"""A fresh source manifest must not hide a cache that the loader will reject."""
import json
from pathlib import Path
import pytest

import matching.library_record_store as store
import runtime.manifest as runtime


@pytest.fixture
def accepted_cache(tmp_path, monkeypatch):
    monkeypatch.delenv('MSSYSTEM_LIBRARY_RECORD_STORE', raising=False)
    identity = {'source_files_sha256': 'accepted-source', 'python': [3, 13],
        'numpy': 'test', 'load_name_only_leads': False,
        'library_top_n_peaks': 15, 'library_relative_threshold': .05}
    monkeypatch.setattr(store, 'processing_identity', lambda: dict(identity))
    cache=tmp_path/'cache'; cache.mkdir(); (tmp_path/'config').mkdir()
    (tmp_path/'config/library_record_store.json').write_text(json.dumps(
        {'enabled': True, 'cache_dir': str(cache)}))
    manifest={'processing_identity': dict(identity), 'record_count': 7,
        'normalized_records_sha256': 'serialized-records'}
    manifest['semantic_sha256']=store.semantic_sha256(manifest)
    (cache/'manifest.json').write_text(json.dumps(manifest))
    acceptance={'schema_version':'mssystem.library-record-store-acceptance/1.0.0',
        'manifest_semantic_sha256':manifest['semantic_sha256'],
        'full_record_roundtrip_sha256':'serialized-records','verified_record_count':7,
        'candidate_order_cases':10000,'candidate_order_mismatches':0,
        'paired_scoring_queries':100,'score_mismatches':0,'shared_database_unchanged':True}
    acceptance['semantic_sha256']=store.semantic_sha256(acceptance)
    (cache/'acceptance.json').write_text(json.dumps(acceptance))
    return tmp_path, cache, identity


def test_accepted_cache_identity_is_ready(accepted_cache):
    root,_,_=accepted_cache
    result=runtime.library_record_store_compatibility(root)
    assert result['compatible'] is True and result['record_count']==7


def test_new_matching_module_invalidates_cached_processing_identity(accepted_cache):
    root,_,identity=accepted_cache
    identity['source_files_sha256']='new-matching-module'
    result=runtime.library_record_store_compatibility(root)
    assert result['compatible'] is False
    assert result['status']=='processing_identity_mismatch'


def test_changed_preprocessing_options_are_also_incompatible(accepted_cache):
    root,_,identity=accepted_cache
    identity['library_top_n_peaks']=40
    assert runtime.library_record_store_compatibility(root)['compatible'] is False


def test_missing_or_tampered_acceptance_is_not_ready(accepted_cache):
    root,cache,_=accepted_cache
    (cache/'acceptance.json').write_text('{}')
    assert runtime.library_record_store_compatibility(root)['compatible'] is False
    (cache/'acceptance.json').unlink()
    assert runtime.library_record_store_compatibility(root)['compatible'] is False


def test_disabled_cache_does_not_block_sqlite_loader(accepted_cache, monkeypatch):
    root,cache,_=accepted_cache
    (cache/'manifest.json').unlink()
    monkeypatch.setenv('MSSYSTEM_LIBRARY_RECORD_STORE','off')
    assert runtime.library_record_store_compatibility(root)=={
        'status':'disabled','compatible':True,'required':False}


def test_cache_metadata_change_invalidates_health_cache(accepted_cache):
    root,cache,_=accepted_cache
    before=runtime.runtime_manifest_dependency_signature(root)
    (cache/'acceptance.json').write_text('{}')
    assert runtime.runtime_manifest_dependency_signature(root)!=before


def test_fresh_manifest_is_not_release_ready_with_incompatible_cache(accepted_cache, monkeypatch):
    root,_,identity=accepted_cache
    source={'source_tree_sha256':'fresh-source','file_count':1,'files':[]}
    catalogs={'compound_roles':{'sha256':'roles'},'screening':{
        'sha256':'screening','library_coverage':{'status':'current'}}}
    manifest={'source':source,'role_catalogs':catalogs,'library':{}}
    path=root/'config/runtime_manifest.json';path.write_text(json.dumps(manifest))
    monkeypatch.setattr(runtime,'source_tree_snapshot',lambda _:source)
    monkeypatch.setattr(runtime,'role_catalog_snapshot',lambda _:catalogs)
    monkeypatch.setattr(runtime,'library_snapshot_freshness',lambda _: {'matches':True})
    import governance.authentic_anchor_registry as anchors
    monkeypatch.setattr(anchors,'audit_active_authentic_anchor_registry',lambda *a,**k:{'passed':True})
    assert runtime.runtime_manifest_freshness(root,path)['release_ready'] is True
    identity['source_files_sha256']='incompatible-current-source'
    status=runtime.runtime_manifest_freshness(root,path)
    assert status['fresh'] is True
    assert status['release_ready'] is False
    assert status['library_record_store']['status']=='processing_identity_mismatch'
