const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const {test} = require('node:test');
const source = fs.readFileSync(path.join(__dirname, '../web/static/js/sample-overview.js'), 'utf8');

function setup(total, {duplicate = false, emptyPage = false, readError = false} = {}) {
  const calls = [];
  const metadata = {main_result_count: 168, unresolved_precursor_event_count: 19, original_field: 'preserved'};
  const rows = Array.from({length: total}, (_, index) => ({
    compound_id: `event:${index}`, identity_key: `identity:${index}`,
    frozen_compound_index: duplicate && index === 500 ? 499 : index,
    record_type: index === 168 ? 'unresolved_precursor_event' : 'unannotated_ion_feature',
    conclusion: {key: index === 168 ? 'unresolved_event' : 'unannotated'},
    representative: {result_index: index, spectrum_id: `scan-${index}`},
  }));
  const context = {console, window: {}, API: {get: async address => {
    calls.push(address);
    const url = new URL(address, 'http://localhost');
    if (url.pathname.endsWith('/status')) return {id: 'saved', completed_batch: {summary: {}}};
    if (url.pathname.endsWith('/flagship')) {
      const offset = Number(url.searchParams.get('offset') || 0);
      return {entities: offset ? [{entity_id: 'last', member_result_indices: [total-1]}] : [{entity_id: 'first', member_result_indices: [0]}], entities_total: 2, has_more: !offset};
    }
    if (url.pathname.endsWith('/compounds')) {
      if (readError) throw new Error('saved SHA mismatch');
      assert.equal(url.searchParams.get('result_view'), 'native');
      const limit = Number(url.searchParams.get('limit'));
      assert.equal(limit, 500);
      const offset = Number(url.searchParams.get('offset') || 0);
      return {metadata, native_total: total, xbt_linked_total: total ? 1 : 0,
        compounds: emptyPage && offset ? [] : rows.slice(offset, offset + limit), has_more: offset + limit < total};
    }
    if (url.pathname.endsWith('/replicate-consistency')) return {status: 'unavailable'};
    return {records: []};
  }}, SavedResultFallback: {load: () => {throw new Error('frozen reads must not fall back to empty local results');}}};
  vm.createContext(context); vm.runInContext(source, context);
  return {data: vm.runInContext('SampleResultData', context), calls, metadata};
}

for (const count of [0, 1310, 2506]) test(`frozen sample loads all ${count} native rows and paged entity evidence`, async () => {
  const {data, calls, metadata} = setup(count);
  const bundle = await data.load('saved');
  assert.equal(bundle.featureGroups.length, count);
  assert.equal(bundle.compounds.native_total, count);
  assert.equal(bundle.compounds.metadata, metadata);
  assert.equal(bundle.flagship.entities.length, 2);
  assert.equal(calls.filter(x => x.includes('/compounds?')).length, Math.max(1, Math.ceil(count/500)));
  if (count) {
    assert.equal(bundle.featureGroups[168].status.key, 'unresolved_event');
    assert.equal(bundle.featureGroups[count-1].index, count-1);
    assert.equal(bundle.featureGroups[count-1].entity.entity_id, 'last');
  }
});
test('duplicate native indices fail instead of claiming a complete result', async () => {
  await assert.rejects(setup(1310, {duplicate:true}).data.load('saved'), /索引缺失或重复/);
});
test('missing page fails instead of showing the first 500 as all results', async () => {
  await assert.rejects(setup(1310, {emptyPage:true}).data.load('saved'), /分页不完整/);
});
test('source verification error remains visible instead of an empty local fallback', async () => {
  await assert.rejects(setup(1310, {readError:true}).data.load('saved'), /saved SHA mismatch/);
});
