from html.parser import HTMLParser
from types import SimpleNamespace
import pytest
from web.frozen_batch_report import render_frozen_sample_report


class TableRows(HTMLParser):
    def __init__(self):
        super().__init__(); self.table=None; self.body=False; self.rows={}
    def handle_starttag(self, tag, attrs):
        if tag=='table': self.table=dict(attrs).get('id')
        if tag=='tbody': self.body=True
        if tag=='tr' and self.body: self.rows[self.table]=self.rows.get(self.table,0)+1
    def handle_endtag(self, tag):
        if tag=='tbody': self.body=False
        if tag=='table': self.table=None


def bridge(*, truncate=False):
    summary={'parsed_searched_ms2_spectra':1310,'matched_ms2_spectra':187,
             'ms1_feature_records':149,'blank_flagged_ms2':880,'blank_clean_view_excluded_ms2':506}
    rows=[{'frozen_compound_index':i,'record_type':'unannotated_ion_feature',
           'conclusion':{'key':'unannotated'},'representative':{'spectrum_id':f'scan-{i}',
           'blank_flagged':False}, 'plant_report':{}} for i in range(1310)]
    rows[168].update(record_type='unresolved_precursor_event',conclusion={'key':'unresolved_event'})
    rows[0].update(lead_candidate={'compound_name':'<script>alert(1)</script>'})
    def read(task, operation, **options):
        if operation=='expected':return {'features':[{'feature_id':f'ms1-{i}',
                    'evidence_level':'MS1_EXACT_MASS_CANDIDATE','exact_identity_confirmed':False,
                    'associated_ms2':[{'spectrum_id':'scan-0'}] if i==0 else []} for i in range(149)]}
        offset=options['offset']; limit=options['limit']
        assert options['result_view']=='native' and options['report_scope']=='all'
        return {'compounds':[] if truncate and offset else rows[offset:offset+limit],
                'native_total':1310,'xbt_linked_total':1,'has_more':offset+limit<1310}
    return SimpleNamespace(tasks={'t':{'filename':'sample_POS3.mzML','completed_batch':{'summary':summary}}},
        entries={'t':{'source_sha256':'s'*64,'results':{'sha256':'n'*64},'conclusions':{'sha256':'c'*64}}},read=read)


def test_report_keeps_full_native_and_all_independent_ms1_rows():
    html=render_frozen_sample_report(bridge(),'t')
    parser=TableRows(); parser.feed(html)
    assert parser.rows=={'all-native-conclusions':1310,'all-ms1-features':149}
    assert '原生结论 1310 条' in html and '独立MS1靶向特征 149 个' in html
    assert '附加筛选有 1 条' in html and 'unresolved_precursor_event' in html
    assert 'ms1-148' in html and 'scan-1309' in html
    assert '来源未核实；谱库假设QC' in html
    assert '关联MS1证据与配对空白' in html and 'XBT_targeted_20260817' in html
    assert '&lt;script&gt;' in html and '<script>alert' not in html
    assert '/api/tasks/t/frozen-tables/native_compounds' in html
    assert '/api/tasks/t/download/native-conclusions.json.xz' in html


def test_report_rejects_truncated_pages():
    with pytest.raises(RuntimeError,match='分页中断'):
        render_frozen_sample_report(bridge(truncate=True),'t')
