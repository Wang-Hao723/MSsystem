"""Measured CE representatives remain explicit and saved-order independent."""
from copy import deepcopy
import csv
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from matching.evidence_network import prepare_feature_nodes
from matching.molecular_network import node_public, to_graphml, to_nodes_csv
from matching.untargeted_workflow import build_lcms_consensus_spectra


def rows():
    return [dict(lcms_feature_id='feature1', spectrum_id='scan-low', precursor_mz=301.,
                 retention_time=5., lcms_feature_rt_minutes=5.,
                 collision_energy=10.,collision_energy_unit='eV',collision_energy_mode='fixed',
                 activation_method='HCD',instrument_configuration_id='IC-low',ion_mode='POS',
                 precursor_charge=1, raw_peaks=[[70.,100.],[90.,50.],[110.,30.]],
                 query_peaks=[[70.,100.],[90.,50.],[110.,30.]]),
            dict(lcms_feature_id='feature1', spectrum_id='scan-high', precursor_mz=301.0004,
                 retention_time=5.01, lcms_feature_rt_minutes=5.,
                 collision_energy=30.,collision_energy_unit='eV',collision_energy_mode='fixed',
                 activation_method='CID',instrument_configuration_id='IC-high',ion_mode='POS',
                 precursor_charge=1,raw_peaks=[[80.,100.],[100.,50.],[120.,30.],[140.,10.]],
                 query_peaks=[[80.,100.],[100.,50.],[120.,30.],[140.,10.]])]


def test_network_uses_saved_primary_representative_in_either_result_order():
    result=rows();artifacts,summary=build_lcms_consensus_spectra(result)
    expected=next(a for a in artifacts if a['primary_for_feature'])
    before=deepcopy(result)
    for ordered in [result,list(reversed(result))]:
        nodes,_=prepare_feature_nodes(ordered)
        node=nodes[0]
        assert node['peaks']==[tuple(p)for p in expected['peaks']]
        assert node['peak_source']=='lcms_feature_representative'
        assert node['usable_member_count']==1
        assert node['peak_evidence']['spectrum_id']==expected['representative_spectrum_id']
        assert node['peak_evidence']['collision_energy']==30.
        assert node['peak_evidence']['artifact_id']==expected['consensus_spectrum_id']
    assert result==before
    assert summary['actual_consensus_spectra']==0


def test_legacy_feature_consensus_keeps_existing_label():
    result=rows()[:1];result[0]['lcms_consensus_peaks']=result[0]['query_peaks'];result[0]['lcms_consensus_spectrum_id']='CS:feature1'
    node=prepare_feature_nodes(result)[0][0]
    assert node['peak_source']=='lcms_feature_consensus'
    assert node['peak_evidence'] is None


def test_public_network_exports_preserve_peak_and_annotation_provenance():
    result=rows();artifacts,_=build_lcms_consensus_spectra(result)
    expected=next(a for a in artifacts if a['primary_for_feature'])
    result[0]['match']={'compound_name':'annotation','final_score':0.9,'cosine':0.9}
    node=prepare_feature_nodes(result)[0][0]
    network={'nodes':[node],'edges':[]}

    public=node_public(node)
    assert public['peak_evidence']['spectrum_id']==expected['representative_spectrum_id']
    assert public['peak_evidence']['result_index']==1
    assert public['spectrum_id']=='scan-high'
    assert public['peak_spectrum_id']=='scan-high'
    assert public['annotation_spectrum_id']=='scan-low'
    assert public['annotation_result_index']==0
    assert public['peak_evidence']['activation_method']=='CID'
    assert public['collision_energy']==30.
    assert public['activation_method']=='CID'
    assert public['instrument_configuration_id']=='IC-high'
    assert public['precursor_mz']==301.0004
    assert public['retention_time']==5.01

    csv_row=next(csv.DictReader(io.StringIO(to_nodes_csv(network))))
    csv_evidence=json.loads(csv_row['peak_evidence'])
    assert csv_evidence==public['peak_evidence']
    assert csv_row['spectrum_id']=='scan-high'
    assert csv_row['annotation_spectrum_id']=='scan-low'
    assert csv_row['annotation_result_index']=='0'

    graphml=to_graphml(network)
    assert graphml is not None
    root=ET.fromstring(graphml)
    keys={item.attrib['id']:item.attrib.get('attr.name') for item in root.findall('{http://graphml.graphdrawing.org/xmlns}key')}
    data={keys[item.attrib['key']]:item.text for item in root.find('.//{http://graphml.graphdrawing.org/xmlns}node').findall('{http://graphml.graphdrawing.org/xmlns}data')}
    assert json.loads(data['peak_evidence'])==public['peak_evidence']
    assert data['spectrum_id']=='scan-high'
    assert data['annotation_spectrum_id']=='scan-low'
    assert data['annotation_result_index']=='0'


def test_network_node_endpoint_returns_peak_provenance_separately(monkeypatch):
    from web import server

    result=rows();build_lcms_consensus_spectra(result)
    result[0]['match']={'compound_name':'annotation','final_score':0.9,'cosine':0.9}
    monkeypatch.setitem(server._tasks,'peak-evidence-task',{'status':'completed'})
    monkeypatch.setattr(server,'_get_task_results',lambda task_id:result)
    payload=server.get_task_network_node(
        'peak-evidence-task',0,cos_thresh=0.7,min_peaks=6,topk=10,
        topk_policy='mutual',mz_tol=0.1,max_component_size=0,
    )
    assert payload['peak_source']=='lcms_feature_representative'
    assert payload['peak_evidence']['spectrum_id']=='scan-high'
    assert payload['peak_evidence']['result_index']==1
    assert payload['spectrum_id']=='scan-high'
    assert payload['collision_energy']==30.
    assert payload['activation_method']=='CID'
    assert payload['instrument_configuration_id']=='IC-high'
    assert payload['precursor_mz']==301.0004
    assert payload['retention_time']==5.01
    assert payload['annotation_spectrum_id']=='scan-low'
    assert payload['annotation_result_index']==0
    assert payload['match']['compound_name']=='annotation'


def test_network_detail_ui_separates_measured_peak_and_annotation_rows():
    source=(Path(__file__).parents[1]/'web/static/js/network.js').read_text()
    assert "['峰谱图编号', peakMeta.spectrum_id" in source
    assert "['峰谱碰撞能', ce]" in source
    assert "['峰谱活化方式', peakMeta.activation_method" in source
    assert "['候选注释谱图', detail.annotation_spectrum_id" in source
    assert "['候选注释结果索引', detail.annotation_result_index" in source
    assert "['Spectrum ID', detail.spectrum_id" not in source


def test_cache_replacement_accepts_exact_reviewed_bytes_and_rejects_unknown_change(tmp_path,monkeypatch):
    import runtime.library_cache_identity as resolver
    import hashlib
    root=tmp_path; (root/'matching').mkdir();p=root/'matching/network.py';p.write_bytes(b'NEW=2\n')
    def record(data):return {'path':'matching/network.py','size_bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
    before=record(b'OLD=1\n');after=record(p.read_bytes());replacement={'before':{k:v for k,v in before.items()if k!='path'},'after':{k:v for k,v in after.items()if k!='path'}}
    monkeypatch.setattr(resolver,'_ISOLATED_ADDITIONS',{})
    monkeypatch.setattr(resolver,'_NON_NORMALIZING_REPLACEMENTS',{'matching/network.py':replacement})
    monkeypatch.setattr(resolver,'_processing_files',lambda root:[after])
    accepted={'source_files_sha256':resolver._semantic_sha256([before]),'options':1}
    current={'source_files_sha256':resolver._semantic_sha256([after]),'options':1}
    assert resolver.processing_identity_compatibility(accepted,current,project_root=root)['compatible']
    tampered={**after,'sha256':'0'*64}
    monkeypatch.setattr(resolver,'_processing_files',lambda root:[tampered])
    unknown={'source_files_sha256':resolver._semantic_sha256([tampered]),'options':1}
    assert not resolver.processing_identity_compatibility(accepted,unknown,project_root=root)['compatible']
    assert not resolver.processing_identity_compatibility(accepted,{**current,'options':2},project_root=root)['compatible']
