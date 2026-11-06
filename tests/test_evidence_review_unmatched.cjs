const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

function createWorkspace(get, english = false) {
  const pane = { innerHTML: '' };
  const graph = { innerHTML: '' };
  const timers = [];
  const plots = [];
  const context = {
    window: { MSSLanguage: { getLanguage: () => english ? 'en' : 'zh' } },
    document: {
      getElementById(id) {
        if (id === 'unmatchedEvidencePane') return pane;
        if (id === 'reviewUnmatchedSpectrum' && pane.innerHTML.includes('id="reviewUnmatchedSpectrum"')) return graph;
        return null;
      },
      querySelectorAll: () => [],
    },
    Tasks: { escapeHtml: value => String(value).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])) },
    API: { BASE: '/api', get },
    SampleResultData: {
      load: async () => ({ groups: [], matched: [], task: { filename: 'fixture.mgf' }, flagship: { entities: [] } }),
      retentionTimeMinutes: row => row.retention_time == null ? null
        : Number(row.retention_time) / (row.retention_time_unit === 'seconds' ? 60 : 1),
    },
    SpectrumComparator: class { render(...args) { plots.push(args); } },
    setTimeout: callback => timers.push(callback),
    URLSearchParams, console,
  };
  vm.createContext(context);
  vm.runInContext(fs.readFileSync(path.resolve(__dirname, '../web/static/js/evidence-review.js'), 'utf8')
    + '\nthis.workspace = EvidenceReview;', context);
  return { workspace: context.workspace, pane, graph, plots,
    flush: () => { while (timers.length) timers.shift()(); } };
}

const row = id => ({ spectrum_id: `scan-${id}`, precursor_mz: 300.12345,
  retention_time: 120, retention_time_unit: 'seconds',
  query_peaks: [[100 + id, 10], [200 + id, 0]],
  dl_confidence: null, dl_top3: [['Class A', 0.8], ['Class B', null]] });

async function testSpectrumFirstAndMissingClassValue() {
  const fixture = row(1);
  const original = JSON.stringify(fixture);
  const w = createWorkspace(async () => ({ result: fixture }));
  await w.workspace.loadTask('task a/b');
  await w.workspace.selectUnmatched(1, false);
  const html = w.pane.innerHTML;
  assert.ok(html.indexOf('id="reviewUnmatchedSpectrum"') < html.indexOf('<details'));
  assert.ok(!html.includes('<details open'));
  assert.ok(html.includes('RT 2.000 min'));
  assert.ok(html.includes('80.00%'));
  assert.ok(html.includes('Class B</span><b>—</b>'));
  assert.ok(!html.includes('Class B</span><b>0.00%'));
  assert.ok(html.includes('/api/tasks/task%20a%2Fb/results/item/1'));
  w.flush();
  assert.equal(w.plots.length, 1);
  assert.equal(w.plots[0][2], null);
  assert.equal(JSON.stringify(w.plots[0][1].peaks), JSON.stringify(fixture.query_peaks));
  assert.equal(JSON.stringify(fixture), original);
}

async function testEnglishAndBoundedTable() {
  const fixture = row(2);
  fixture.query_peaks = Array.from({ length: 305 }, (_, i) => [50 + i, i]);
  const w = createWorkspace(async () => ({ result: fixture }), true);
  await w.workspace.loadTask('test');
  await w.workspace.selectUnmatched(2, false);
  assert.ok(w.pane.innerHTML.includes('Measured MS/MS spectrum'));
  assert.ok(w.pane.innerHTML.includes('Showing 300 / 305'));
  assert.equal((w.pane.innerHTML.match(/<tr>/g) || []).length, 301);
  w.flush();
  assert.equal(w.plots[0][1].peaks.length, 305);
}

async function testLateResponseCannotReplaceSelectedSpectrum() {
  const pending = new Map();
  const w = createWorkspace(url => new Promise(resolve => pending.set(url.split('/').pop(), resolve)));
  await w.workspace.loadTask('test');
  const first = w.workspace.selectUnmatched(1, false);
  const second = w.workspace.selectUnmatched(2, false);
  pending.get('2')({ result: row(2) });
  await second;
  pending.get('1')({ result: row(1) });
  await first;
  assert.ok(w.pane.innerHTML.includes('scan-2'));
  assert.ok(!w.pane.innerHTML.includes('scan-1'));
  w.flush();
  assert.equal(w.plots.length, 1);
  assert.equal(w.plots[0][1].peaks[0][0], 102);
}

async function testDelayedPlotCannotDrawOldSpectrum() {
  let resolveSecond;
  const w = createWorkspace(url => url.endsWith('/1') ? Promise.resolve({ result: row(1) })
    : new Promise(resolve => { resolveSecond = resolve; }));
  await w.workspace.loadTask('test');
  await w.workspace.selectUnmatched(1, false);
  const second = w.workspace.selectUnmatched(2, false);
  w.flush();
  assert.equal(w.plots.length, 0);
  resolveSecond({ result: row(2) });
  await second;
  w.flush();
  assert.equal(w.plots.length, 1);
  assert.equal(w.plots[0][1].peaks[0][0], 102);
}

async function testInvalidValuesStayMissingAndSourceIndicesRemain() {
  const fixture = row(3);
  fixture.precursor_mz = ' ';
  fixture.query_peaks = [[null, 1], ['', 2], [true, 3], [100, ' '], [101, 0], [102, 5]];
  fixture.dl_top3 = [['<script>bad</script>', ' '], ['Other', false]];
  const w = createWorkspace(async () => ({ result: fixture }));
  await w.workspace.loadTask('test');
  await w.workspace.selectUnmatched(3, false);
  assert.ok(w.pane.innerHTML.includes('m/z 未记录'));
  assert.ok(w.pane.innerHTML.includes('2 / 6'));
  assert.ok(w.pane.innerHTML.includes('<td>5</td>'));
  assert.ok(!w.pane.innerHTML.includes('<script>'));
  assert.ok(!w.pane.innerHTML.includes('0.00%'));
  w.flush();
  assert.equal(JSON.stringify(w.plots[0][1].peaks), '[[101,0],[102,5]]');
}

(async () => {
  await testSpectrumFirstAndMissingClassValue();
  await testEnglishAndBoundedTable();
  await testLateResponseCannotReplaceSelectedSpectrum();
  await testDelayedPlotCannotDrawOldSpectrum();
  await testInvalidValuesStayMissingAndSourceIndicesRemain();
  console.log('evidence-review-unmatched: 5/5 behavior groups passed');
})().catch(error => { console.error(error); process.exitCode = 1; });
