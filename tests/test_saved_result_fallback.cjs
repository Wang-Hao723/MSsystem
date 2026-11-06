const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');

const sourcePath = require('path').resolve(__dirname, '../web/static/js/saved-result-fallback.js');
const source = fs.readFileSync(sourcePath, 'utf8');

function loadModule(get) {
  const window = {};
  const context = vm.createContext({ window, API: { get }, console, encodeURIComponent });
  vm.runInContext(source, context, { filename: sourcePath });
  return window.SavedResultFallback;
}

async function testPaginationAndGrouping() {
  const requested = [];
  const pages = [
    {
      has_more: true,
      results: [
        { result_index: 8, spectrum_id: 'scan-8', match: { identity_key: 'saved-A', compound_name: 'A', smiles: 'CC', final_score: 0.4, confidence: 'High' } },
        { result_index: 9, spectrum_id: 'scan-9', query_peaks: [[10, 2]] },
        { result_index: 12, spectrum_id: 'scan-12', match: { identity_key: 'saved-A', compound_name: 'A', smiles: 'CC', final_score: 0.8, confidence: 'High' } },
      ],
    },
    {
      has_more: false,
      results: [
        { result_index: 41, spectrum_id: 'scan-41', match: { inchikey: 'ABCDEFGHIJKLMN-ABCDEFGHIJ-A', compound_name: 'B', smiles: 'CO', final_score: 0.6, confidence: 'Medium' } },
      ],
    },
  ];
  const fallback = loadModule(async url => {
    requested.push(url);
    return pages.shift();
  });
  const priorError = Object.assign(new Error('compound artifact unavailable'), { status: 500, code: 'SOURCE_MISSING' });
  const bundle = await fallback.load('task a/b', { id: 'task a/b', filename: 'saved.mgf' }, priorError);

  assert.equal(requested.length, 2);
  assert.match(requested[0], /tasks\/task%20a%2Fb\/results\?projection=summary&limit=5000&offset=0$/);
  assert.match(requested[1], /offset=3$/);
  assert.equal(bundle.savedOnly, true);
  assert.equal(bundle.groups.length, 2);
  assert.equal(bundle.matched.length, 3);
  assert.equal(bundle.unmatched.length, 1);
  assert.strictEqual(bundle.unannotatedRecords, bundle.unmatched);
  assert.equal(bundle.flagship.entities.length, 0);
  assert.equal(bundle.groups[0].rep.result_index, 12);
  assert.equal(bundle.groups[0].rep.match.final_score, 0.8);
  assert.equal(bundle.groups[0].rows.length, 2);
  assert.equal(bundle.groups[0].status.key, 'saved_candidate');
  assert.equal(bundle.groups[0].status.countedAsSampleCompound, false);
  assert.equal(bundle.groups[0].status.countedAsPlantCandidate, false);
  assert.equal(bundle.groups[0].rep.match.confidence, 'High');
  assert.equal(bundle.groups[0].rep.match.smiles, 'CC');
  assert.equal(bundle.compounds.metadata.saved_candidate_group_count, 2);
  assert.equal(bundle.compounds.metadata.main_result_count, 0);
  assert.equal(bundle.sourceFailure.status, 500);
  assert.equal(bundle.sourceFailure.code, 'SOURCE_MISSING');
}

async function testEmptyResults() {
  const fallback = loadModule(async () => ({ results: [], has_more: false }));
  const bundle = await fallback.load('empty-task', { id: 'empty-task' });
  assert.equal(bundle.groups.length, 0);
  assert.equal(bundle.matched.length, 0);
  assert.equal(bundle.unmatched.length, 0);
  assert.equal(bundle.compounds.metadata.saved_result_count, 0);
}

async function testInterfaceFailureStillThrows() {
  const expected = new Error('results endpoint failed');
  const fallback = loadModule(async () => { throw expected; });
  await assert.rejects(() => fallback.load('broken-task', {}, new Error('original failure')), expected);
}

async function testMissingIdentityAndIndexAreNotInvented() {
  const fallback = loadModule(async () => ({ has_more: false, results: [
    { result_index: null, match: { compound_name: '甲候选', inchikey: 'nan', smiles: 'nan', formula: 'C2H6O' } },
    { result_index: 9, match: { compound_name: '乙候选', inchikey: 'nan', smiles: 'nan', formula: 'C2H6O' } },
  ] }));
  const bundle = await fallback.load('missing-identity');
  assert.equal(bundle.groups.length, 2);
  assert.equal(bundle.groups.find(g => g.rep.match.compound_name === '甲候选').rep.result_index, null);
  assert.ok(bundle.groups.every(g => g.peakClusters === null));
}

(async () => {
  await testPaginationAndGrouping();
  await testEmptyResults();
  await testInterfaceFailureStillThrows();
  await testMissingIdentityAndIndexAreNotInvented();
  process.stdout.write('saved-result-fallback: 4/4 behavior checks passed\n');
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
