const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { test } = require('node:test');

const root = path.resolve(__dirname, '..');
const source = fs.readFileSync(path.join(root, 'web/static/js/completed-batches.js'), 'utf8');
const index = fs.readFileSync(path.join(root, 'web/static/index.html'), 'utf8');
const batchId = 'xbt-20261021-ms1-msms';
const base = `/api/completed-batches/${batchId}`;
const flush = () => new Promise(resolve => setImmediate(resolve));

class Element {
  constructor(dataset = {}) {
    this.dataset = dataset; this.innerHTML = ''; this.textContent = ''; this.value = '';
    this.hidden = false; this.attrs = {}; this.events = {};
    this.classes = new Set();
    this.classList = { toggle: (name, enabled) => enabled ? this.classes.add(name) : this.classes.delete(name) };
  }
  addEventListener(type, callback) { this.events[type] = callback; }
  setAttribute(name, value) { this.attrs[name] = String(value); }
  removeAttribute(name) { delete this.attrs[name]; if (name === 'src') this.src = ''; }
  getAttribute(name) { return this.attrs[name] ?? this[name] ?? null; }
  scrollIntoView() {}
  async fire(type) { await this.events[type]?.({ target: this }); await flush(); }
}

function setup() {
  const nodes = new Map(Array.from(index.matchAll(/\bid="([^"]+)"/g), match => [match[1], new Element()]));
  const entries = [new Element(), new Element(), new Element()];
  const tabs = ['overview', 'report', 'exports'].map(cbTab => new Element({ cbTab }));
  const reportButton = new Element({ cbTab: 'report' });
  const panels = ['overview', 'report', 'exports'].map(cbPanel => new Element({ cbPanel }));
  const samples = Array.from({ length: 149 }, (_, n) => ({ task_id: `task-${n}`, sample_key: `sample-${n}`, sample_id: `F005-${n}`, polarity: n % 2 ? 'NEG' : 'POS', replicate: 1, searched_ms2_spectra: 10, native_conclusion_records: 10, matched_ms2_spectra: 2, ms1_features: 5, xbt_named_candidate_rows: n ? 0 : 2, blank_flagged_ms2: 3, blank_clean_view_excluded_ms2: 1 }));
  const keys = ['file_manifest', 'xbt_csv', 'xbt_jsonl', 'qc', 'ms1', 'msms_table', 'ms1_table', 'qc_table', 'summary'];
  const batch = { id: batchId, title: '宣痹汤 / XBT', completed_samples: 149, total_samples: 149, paired_blank_files: 102, paired_groups: 50, report_url: `${base}/report`, index_url: `${base}/index`, download_all_url: `${base}/download-all`, xbt_dataset: 'XBT_targeted_20260817', ms1_evidence_counts: { MS1_EXACT_MASS_CANDIDATE: 103161, MS1_ISOTOPE_SUPPORTED: 31448, MS1_MS2_DIAGNOSTIC_PROBABLE: 11 }, totals: { ms1_features: 134620, searched_ms2_spectra: 253234, matched_ms2_spectra: 48361, xbt_named_candidate_rows: 1872, ambiguity_group_rows: 1182, standard_confirmed: 0, background_qc_rows: 686, external_qc_rows: 1563, name_qc_rows: 96, unlinked_qc_rows: 250673 }, samples, exports: keys.map(name => ({ name, path: `/exports/${name}.csv.xz`, size: 1024, rows: 149, download_url: `${base}/exports/${name}`, csv_url: ['file_manifest', 'xbt_csv', 'msms_table', 'ms1_table', 'qc_table'].includes(name) ? `${base}/exports/${name}?decompress=true` : null })) };
  const calls = [];
  const native = {native_conclusion_rows:1490, record_type_counts:{compound_or_isomer_group:100,unannotated_ion_feature:1390}, artifacts:{native_compounds:{path:'/exports/all_native_conclusions.csv.xz',size:4096}},
    report_url:`${base}/native-report`,index_url:`${base}/native-index`,download_all_url:`${base}/native-download-all`,native_compounds_csv_url:`${base}/native-exports/native_compounds?decompress=true`,native_compounds_xz_url:`${base}/native-exports/native_compounds`,native_table_url:`${base}/native-table`};
  let delayed = null;
  const context = {
    console, URLSearchParams, setTimeout,
    window: { location: { hash: `#completed-batches?batch=${batchId}` } },
    document: {
      getElementById: id => nodes.get(id),
      querySelectorAll: selector => selector === '[data-completed-batch-entry]' ? entries : selector === '[data-cb-panel]' ? panels : selector === '.cb-tabs [data-cb-tab]' ? tabs : selector === '[data-cb-tab]' ? [reportButton, ...tabs] : [],
    },
    fetch: async url => {
      calls.push(url);
      if (url === '/api/completed-batches') return { ok: true, json: async () => ({ batches: [batch] }) };
      if (url === base) return { ok: true, json: async () => batch };
      if (url === `${base}/native-delivery`) return {ok:true,json:async()=>native};
      if (delayed) { const wait = delayed; delayed = null; return wait; }
      const filtered = url.includes('sample_key=sample-1');
      const payload = { columns: ['sample_key', 'display_name', 'formula'], rows: filtered ? [] : [{ sample_key: 'sample-0', display_name: '<img onerror=bad>', formula: 'C6H12O6' }, { sample_key: 'sample-0', display_name: 'second saved candidate', formula: 'C7H8O' }], total: filtered ? 0 : 2, has_more: false };
      return { ok: true, json: async () => payload };
    },
  };
  vm.createContext(context); vm.runInContext(source, context);
  return { nodes, entries, tabs, panels, reportButton, calls, batch, module: vm.runInContext('CompletedBatches', context), delay: value => { delayed = value; } };
}

test('completed-batch entry is discoverable in results, exports and dashboard', async () => {
  const fixture = setup(); fixture.module.init(); await fixture.module.enterPage();
  assert.equal(fixture.entries.length, 3);
  for (const entry of fixture.entries) {
    assert.match(entry.innerHTML, /149 \/ 149/); assert.match(entry.innerHTML, /102/);
    assert.match(entry.innerHTML, /#completed-batches\?batch=xbt-20261021-ms1-msms/);
  }
  assert.equal((fixture.nodes.get('completedBatchSamples').innerHTML.match(/<tr/g) || []).length, 149);
  assert.equal((fixture.nodes.get('completedBatchExports').innerHTML.match(/<tr/g) || []).length, 10);
  assert.equal(fixture.nodes.get('completedBatchDownloadAll').href, `${base}/native-download-all`);
  assert.ok(fixture.calls.some(url=>url.startsWith(`${base}/native-table?`)));
  assert.match(fixture.nodes.get('completedBatchSamples').innerHTML,/\/api\/tasks\/task-148\/report/);
  assert.match(fixture.nodes.get('completedBatchMetrics').innerHTML,/完整原生结论记录/);
  assert.doesNotMatch(fixture.nodes.get('completedBatchMetrics').innerHTML,/XBT命名候选记录/);
  assert.ok(fixture.calls.every(url => url.startsWith('/api/completed-batches')));
});

test('preview preserves saved row order, escapes names and clears stale rows when filters change', async () => {
  const fixture = setup(); fixture.module.init(); await fixture.module.enterPage();
  const table = fixture.nodes.get('completedBatchPreviewRows');
  assert.match(table.innerHTML, /&lt;img onerror=bad&gt;/); assert.doesNotMatch(table.innerHTML, /<img/);
  assert.ok(table.innerHTML.indexOf('&lt;img') < table.innerHTML.indexOf('second saved candidate'));
  let release;
  fixture.delay(new Promise(resolve => { release = resolve; }));
  const filter = fixture.nodes.get('completedBatchSampleFilter'); filter.value = 'sample-1';
  const pending = filter.fire('change'); await flush();
  assert.equal(table.innerHTML, '');
  release({ ok: true, json: async () => ({ columns: ['display_name'], rows: [], total: 0, has_more: false }) });
  await pending; await flush();
  assert.equal(table.innerHTML, ''); assert.match(fixture.nodes.get('completedBatchPreviewStatus').textContent, /没有记录/);
  assert.equal(fixture.nodes.get('completedBatchNext').disabled, true);
  assert.match(fixture.nodes.get('completedBatchTableCsv').href, /decompress=true/);
});

test('report opens the frozen report inside V16 and retains download links', async () => {
  const fixture = setup(); fixture.module.init(); await fixture.module.enterPage();
  await fixture.reportButton.fire('click');
  assert.equal(fixture.nodes.get('completedBatchReportFrame').src, `${base}/native-report`);
  assert.equal(fixture.panels.find(panel => panel.dataset.cbPanel === 'report').hidden, false);
  assert.equal(fixture.panels.find(panel => panel.dataset.cbPanel === 'overview').hidden, true);
  assert.equal(fixture.nodes.get('completedBatchDownloadReport').href, `${base}/native-report?download=true`);
  assert.equal(fixture.nodes.get('completedBatchDownloadIndex').href, `${base}/native-index`);
});
