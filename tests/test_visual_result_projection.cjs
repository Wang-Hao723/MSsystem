const fs = require('fs');
const vm = require('vm');
const assert = require('assert/strict');
const root = require('path').resolve(__dirname, '..');
const context = { SampleResultData: { retentionTimeMinutes(row) {
  return Number(row.retention_time) / (row.retention_time_unit === 'seconds' ? 60 : 1);
} } };
vm.createContext(context);
vm.runInContext(fs.readFileSync(root + '/web/static/js/visual-results.js', 'utf8') + '\nthis.visual = VisualResults;', context);
const visual = context.visual;
function group(i) {
  return { index: i, status: { key: 'high' }, spectra: 2,
    rep: { result_index: 1000 - i, precursor_mz: 100 + i,
      retention_time: i * 60, retention_time_unit: 'seconds',
      match: { compound_name: 'Candidate ' + i, smiles: 'CCO', formula: 'C2H6O', cosine: .8, mass_ppm: null } },
    compound: {}, leadCandidate: {}, candidateMemberCount: 1 };
}
const before = [group(0), group(1), group(2), group(3)];
before[2].rep.retention_time = null;
before[3].rep.precursor_mz = '';
const html = visual.summary({ groups: before, task: {} }, {}, false);
assert.equal((html.match(/<circle /g) || []).length, 2);
assert.equal(visual.rt(before[0].rep), 0);
assert.equal(visual.rt(before[1].rep), 1);
assert.equal(visual.rt(before[2].rep), null);
assert.equal(visual.finite(''), null);
assert.equal(visual.finite(null), null);
assert.equal(visual.finite(0), 0);

const groups = Array.from({ length: 70 }, (_, i) => group(i)).reverse();
const frozen = JSON.stringify(groups);
const allIndices = [];
for (let page = 0; page < 3; page++) {
  const pageHtml = visual.catalog(groups, { taskId: 'snapshot', page });
  const links = [...pageHtml.matchAll(/SampleOverview.openEvidence\((\d+)\)/g)].map(m => Number(m[1]));
  allIndices.push(...links);
  links.forEach(i => assert.ok(pageHtml.includes(`/result/${1000 - i}?format=svg`)));
}
assert.equal(allIndices.length, 70);
assert.equal(new Set(allIndices).size, 70);
assert.deepEqual(allIndices, groups.map(g => g.index));
assert.equal(JSON.stringify(groups), frozen);

const noIndex = group(0);
noIndex.rep.result_index = null;
assert.ok(!visual.catalog([noIndex], { taskId: 'snapshot' }).includes('/result/0?'));
const noMass = visual.catalog([noIndex], { taskId: 'snapshot', mode: 'table' });
assert.ok(noMass.includes('— ppm'));
const malicious = group(0);
malicious.rep.match.compound_name = '<script>alert(1)</script>';
assert.ok(!visual.catalog([malicious]).includes('<script>'));
const saved = group(1); saved.status.key = 'saved_candidate';
const savedHtml = visual.summary({ groups: [saved], task: {}, savedOnly: true }, { high: 0, isomerHigh: 0, confirmed: 0 });
assert.equal((savedHtml.match(/<b>—<\/b>/g) || []).length, 3);
assert.ok(!savedHtml.includes('data-overview-filter="confirmed"'));
console.log('visual-result-projection: coordinates, immutable indices/pagination, missing values/escaping, saved-only evidence: 4/4 passed');
