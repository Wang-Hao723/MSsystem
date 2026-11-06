const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const sourcePath = path.resolve(__dirname, '../web/static/js/spectrum-compare.js');
const source = fs.readFileSync(sourcePath, 'utf8');

function render(query, library, pairs, options) {
  const calls = [];
  const container = { id: 'plot' };
  const window = {};
  const context = vm.createContext({
    window,
    document: { getElementById: id => id === 'plot' ? container : null },
    Plotly: {
      newPlot: (...args) => calls.push(args),
      purge() {},
    },
    console,
  });
  vm.runInContext(source, context, { filename: sourcePath });
  const comparator = new window.SpectrumComparator();
  comparator.render('plot', query, library, pairs, options);
  assert.equal(calls.length, 1);
  const [, traces, layout] = calls[0];
  return { comparator, traces, layout };
}

function traceByName(traces, fragment) {
  return traces.find(trace => trace.name.includes(fragment));
}

function finiteCoordinates(traces) {
  return traces.every(trace => trace.x.concat(trace.y)
    .filter(value => value !== null)
    .every(value => typeof value === 'number' && Number.isFinite(value)));
}

function testExactRecordedMatchesOnly() {
  const query = { name: 'Q', peaks: [[100.001, 20], [100.004, 10], [150, 5]] };
  const library = { name: 'L', peaks: [[200.001, 8], [200.004, 4]] };
  const frozen = JSON.stringify({ query, library });
  const result = render(query, library, [[100.001, 200.001]], {});
  const matched = result.traces.filter(trace => trace.name.includes('recorded match'));
  assert.ok(result.traces.every(trace => trace.type === 'scatter' && trace.mode === 'lines+markers'));
  assert.deepEqual(Array.from(matched[0].x).filter(x => x !== null), [100.001, 100.001]);
  assert.equal(matched[0].customdata[1][0], '100.001');
  assert.equal(matched[0].customdata[1][1], 100);
  assert.deepEqual(Array.from(traceByName(result.traces, 'other peak').x).filter(x => x !== null),
    [100.004, 100.004, 150, 150]);
  assert.equal(matched.length, 2);
  assert.equal(JSON.stringify({ query, library }), frozen);
}

function testPairShapesAreEquivalent() {
  const query = { peaks: [[100, 1], [101, 2]] };
  const library = { peaks: [[200, 3], [201, 4]] };
  const arrayResult = render(query, library, [[101, 201]], {});
  const dictResult = render(query, library, [{ query_mz: 101, lib_mz: 201 }], {});
  assert.equal(JSON.stringify(arrayResult.traces), JSON.stringify(dictResult.traces));
  assert.equal(JSON.stringify(arrayResult.layout.annotations), JSON.stringify(dictResult.layout.annotations));
}

function testMissingAndAmbiguousEndpointsDoNotHighlight() {
  const query = { peaks: [[100, 1], [100, 2], [110, 3]] };
  const library = { peaks: [[200, 2], [210, 3]] };
  const result = render(query, library, [
    [100, 200],
    [110, 999],
    { query_mz: '', lib_mz: 210 },
  ], {});
  assert.ok(!result.traces.some(trace => trace.name.includes('recorded match')));
  assert.ok(result.layout.annotations.some(annotation => annotation.text.endsWith('0/3')));
}

function testSavedScoreContract() {
  const spectrum = { peaks: [[100, 2]] };
  const zero = render(spectrum, spectrum, [], { score: 0 });
  assert.equal(zero.layout.title.text, 'MS/MS spectrum comparison — Saved cosine: 0');
  const stored = render(spectrum, spectrum, [[100, 100]], { score: 0.8123, scoreLabel: 'Library cosine' });
  assert.equal(stored.layout.title.text, 'MS/MS spectrum comparison — Library cosine: 0.8123');
  const absent = render(spectrum, spectrum, [[100, 100]], {});
  assert.equal(absent.layout.title.text, 'MS/MS spectrum comparison');
  assert.ok(!absent.layout.title.text.includes('Cosine'));
  for (const invalid of [null, '', Infinity, NaN, '0.9', false, {}]) {
    assert.equal(render(spectrum, spectrum, [], { score: invalid }).layout.title.text,
      'MS/MS spectrum comparison');
  }
}

function testMissingEmptyZeroAndBadPeaks() {
  const queryOnly = render({ peaks: [
    [50, 0], ['', 2], [70, ''], [80, -1], [90, Infinity], [100, false], [{}, 2],
  ] }, null, null, {});
  assert.deepEqual(Array.from(queryOnly.layout.yaxis.range), [0, 115]);
  assert.ok(queryOnly.layout.annotations.some(annotation => annotation.text === 'Reference spectrum unavailable'));
  assert.ok(finiteCoordinates(queryOnly.traces));
  assert.deepEqual(Array.from(queryOnly.traces).flatMap(trace => Array.from(trace.y))
    .filter(value => value !== null), [0, 0]);

  const libraryOnly = render(null, { peaks: [[75, 4]] }, [], {});
  assert.deepEqual(Array.from(libraryOnly.layout.yaxis.range), [-115, 115]);
  assert.deepEqual(Array.from(libraryOnly.traces[0].y), [0, -100, null]);
  assert.ok(finiteCoordinates(libraryOnly.traces));

  const bothEmpty = render(undefined, { peaks: [['', ''], [NaN, 2], [20, -3]] }, undefined, {});
  assert.equal(bothEmpty.traces.length, 0);
  assert.ok(bothEmpty.layout.annotations.some(annotation => annotation.text === 'No valid spectrum peaks available'));
  assert.ok(Number.isFinite(bothEmpty.layout.xaxis.range[1]));
}

function testToleranceAndEscaping() {
  const withinSerializationError = render(
    { name: '<img src=x onerror=1>', peaks: [[100.0000000005, 1]] },
    { name: '<script>alert(1)</script>', peaks: [[200, 1]] },
    [{ query_mz: 100, lib_mz: 200 }],
    { title: '<b>unsafe</b>', score: 0, text: { locatedPairs: '<svg onload=1>' } },
  );
  assert.equal(withinSerializationError.traces.filter(trace => trace.name.includes('recorded match')).length, 2);
  const serialized = JSON.stringify(withinSerializationError.layout);
  assert.ok(!serialized.includes('<script>'));
  assert.ok(!serialized.includes('<img'));
  assert.ok(!serialized.includes('<svg'));
  assert.ok(withinSerializationError.layout.title.text.includes('&lt;b&gt;unsafe&lt;/b&gt;'));
}

testExactRecordedMatchesOnly();
testPairShapesAreEquivalent();
testMissingAndAmbiguousEndpointsDoNotHighlight();
testSavedScoreContract();
testMissingEmptyZeroAndBadPeaks();
testToleranceAndEscaping();
const indexed = render({peaks:[[100,1],[100,2]]}, {peaks:[[100,1],[100,2]]}, [
  {query_peak_index:1,library_peak_index:0,query_mz:100,lib_mz:100},
], {});
assert.equal(indexed.traces.filter(t=>t.name.includes('recorded match')).length,2);
assert.deepEqual(Array.from(indexed.traces.find(t=>t.name==='Reference — recorded match').y),[0,-50,null]);
const wrongIndex = render({peaks:[[100,1]]}, {peaks:[[100,1]]}, [
  {query_peak_index:5,library_peak_index:0,query_mz:100,lib_mz:100},
], {});
assert.ok(!wrongIndex.traces.some(t=>t.name.includes('recorded match')));
assert.equal(indexed.comparator.constructor.savedEvidence({spectral_evidence:{schema:'wrong'}}), null);
const preciseScore = render({peaks:[[100,1]]}, {peaks:[[100,1]]}, [], {score:0.12345678,scorePrecision:4});
assert.ok(preciseScore.layout.title.text.endsWith('0.1235'));
console.log('spectrum-compare: exact matching, pair formats, ambiguity, saved score, missing/bad peaks, escaping, indexed evidence: 7/7 passed');
