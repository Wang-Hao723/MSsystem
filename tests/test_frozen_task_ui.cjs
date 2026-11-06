const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const {test} = require('node:test');

const root = path.resolve(__dirname, '..');
const task = {id: 'abcdef123456', filename: '<img src=x onerror=bad>.mzML', completed_batch: {
  id: 'xbt-20261021-ms1-msms', sample_key: 'F005_POS1__0011223344556677',
  report_url: '/api/completed-batches/xbt-20261021-ms1-msms/report',
  download_all_url: '/api/completed-batches/xbt-20261021-ms1-msms/download-all',
  summary: {parsed_searched_ms2_spectra: 20, matched_ms2_spectra: 3,
    ms1_feature_records: 11, xbt_named_candidate_rows: 2, blank_flagged_ms2: 5,
    blank_clean_view_excluded_ms2: 4}
}};
const escapeHtml = value => String(value).replace(/[&<>"']/g, ch => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[ch]));

function setup() {
  const calls = [];
  const area = {innerHTML: '', textContent: ''};
  const toolbar = {classList: {add() {}, remove() {}}};
  const context = {console, URLSearchParams, TextEncoder, setTimeout,
    Tasks: {escapeHtml}, window: {},
    document: {getElementById: id => id === 'sampleReportToolbar' ? toolbar : area},
    API: {BASE: '/api', get: async url => {
      calls.push(url);
      if (url.endsWith('/status')) return task;
      if (url.includes('/compounds?')) return {compounds: [{frozen_compound_index: 0, record_type: 'unresolved_precursor_event', conclusion: {key: 'unresolved_event'}, representative: {blank_flagged: false, precursor_mz: 0}}], filtered_total: 1310, has_more: true};
      return {columns: ['name'], rows: [{name: '<script>bad</script>'}], total: 51, has_more: true};
    }},
    SampleResultData: {load: () => {throw new Error('Must not build a new report bundle');}},
  };
  vm.createContext(context);
  vm.runInContext(fs.readFileSync(path.join(root, 'web/static/js/frozen-task-ui.js'), 'utf8'), context);
  return {context, calls, area, ui: vm.runInContext('FrozenTaskUI', context)};
}

test('native result panel keeps MS1/MSMS/blank counts, all tables and safe names', () => {
  const {ui} = setup();
  const html = ui.panel(task);
  assert.match(html, /20 个 MS\/MS/);
  assert.match(html, /11 个 MS1/);
  assert.match(html, /标准品确证 0/);
  assert.match(html, /&lt;img/);
  assert.doesNotMatch(html, /<img src=x/);
  for (const key of ['native_compounds', 'xbt_csv', 'msms_table', 'ms1_table', 'qc_table']) {
    assert.ok(html.includes(`/api/tasks/${task.id}/frozen-tables/${key}`));
  }
  assert.ok(html.includes(task.completed_batch.report_url));
  assert.ok(html.includes(`/api/completed-batches/${task.completed_batch.id}/native-download-all`));
  assert.ok(html.includes(`/api/tasks/${task.id}/report`));
  assert.ok(html.includes(`/api/tasks/${task.id}/download/native-conclusions.json.xz`));
  assert.match(html, /附加视图/);
  assert.equal(ui.panel({id: 'old-task'}), '');
});

test('native preview sends exact sample key and paginates original saved table', async () => {
  const {ui, calls, area} = setup();
  const button = {closest: () => ({querySelector: () => area})};
  await ui.preview(button, task.id, 'ms1_table', 50);
  assert.equal(calls.length, 2);
  const url = new URL(calls[1], 'http://localhost');
  assert.equal(url.searchParams.get('sample_key'), task.completed_batch.sample_key);
  assert.equal(url.searchParams.get('offset'), '50');
  assert.match(area.innerHTML, /&lt;script&gt;/);
  assert.match(area.innerHTML, /共 51 条/);
});

test('each frozen sample report loads its own full evidence HTML', async () => {
  const {context, calls, area} = setup();
  vm.runInContext(fs.readFileSync(path.join(root, 'web/static/js/report-generator.js'), 'utf8'), context);
  await vm.runInContext('SampleReport', context).load(task.id);
  assert.deepEqual(calls, [`/tasks/${task.id}/status`]);
  assert.match(area.innerHTML, /<iframe/);
  assert.ok(area.innerHTML.includes(`/api/tasks/${task.id}/report`));
  assert.match(area.innerHTML, /当前样本全量/);
  assert.doesNotMatch(area.innerHTML, /src="\/api\/completed-batches/);
});

test('native preview reads full source and preserves zero and false evidence values', async () => {
  const {ui, calls, area} = setup();
  await ui.preview({closest: () => ({querySelector: () => area})}, task.id, 'native_compounds', 0);
  assert.match(calls[1], /result_view=native/);
  assert.match(area.innerHTML, /共 1,310 条/);
  assert.match(area.innerHTML, /unresolved_precursor_event/);
  assert.match(area.innerHTML, /<td>0<\/td>/);
  assert.match(area.innerHTML, /<td>false<\/td>/);
});

test('task queue uses source MS1/native/blank counts and a per-file report', () => {
  const context = {console, window: {}, document: {createElement: () => {
    let value=''; return {set textContent(text) {value=escapeHtml(text);}, get innerHTML() {return value;}};
  }}};
  vm.createContext(context);
  const source = fs.readFileSync(path.join(root, 'web/static/js/tasks.js'), 'utf8')
    .replace('    init, loadTasks,', '    renderTaskRow, init, loadTasks,');
  vm.runInContext(source, context);
  const saved = {...task,status:'completed',total_spectra:1310,matched_spectra:187,
    match_rate:14.3,progress:100,unique_compounds:null,params:{},completed_batch:{...task.completed_batch,
    summary:{compound_result_records:1310,ms1_feature_records:149,blank_flagged_ms2:880,blank_clean_view_excluded_ms2:506}},
    blank_summary:{requested:true,applicable:true,filter_policy:'soft_evidence'}};
  const row=vm.runInContext('Tasks',context).renderTaskRow(saved);
  assert.match(row,/187\/1310/); assert.match(row,/1,310 条原生结论/); assert.match(row,/149 个MS1特征/);
  assert.match(row,/标记 880 · 净化视图排除 506/);
  assert.doesNotMatch(row,/0 条被扣除/);
  assert.ok(row.includes(`/api/tasks/${task.id}/report`));
});
