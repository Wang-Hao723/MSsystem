const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const {test} = require('node:test');

const root = path.resolve(__dirname, '..');

function dataLayer() {
  const context = {
    console,
    window: {MSSLanguage: {getLanguage: () => 'en'}},
    document: {},
  };
  vm.createContext(context);
  vm.runInContext(
    fs.readFileSync(path.join(root, 'web/static/js/sample-overview.js'), 'utf8'),
    context,
  );
  return vm.runInContext('SampleResultData', context);
}

test('saved candidate stays distinct from evaluated entity Medium', () => {
  const data = dataLayer();
  const saved = data.statusFromConclusion({
    key: 'candidate', label: '中等证据候选', counted_as_sample_compound: true,
  });
  assert.equal(saved.key, 'candidate');
  assert.doesNotMatch(saved.label, /中等|medium/i);
  assert.equal(data.statusFromConclusion({
    key: 'candidate', entity_confidence: 'Medium', entity_evaluated: true,
  }).key, 'medium');
  assert.equal(data.statusFromConclusion({key: 'future_tier'}).key, 'ungraded');
});

function compound({leadKey = 'AAAAAAAAAAAAAA-BBBBBBBBBB-C', leadIndex = 7} = {}) {
  return {
    compound_id: 'compound:one',
    identity_key: 'std2d:AAAAAAAAAAAAAA|formula:C6H6',
    identity_scope: '2d_connectivity',
    formula: 'C6H6',
    conclusion: {key: 'candidate'},
    lead_candidate: {
      compound_name: 'Bound candidate', formula: 'C6H6', inchikey: leadKey,
      result_index: leadIndex, cosine: 0.934, mass_ppm: -1.7,
    },
    representative: {
      result_index: 7, spectrum_id: 'scan-7',
      match: {compound_name: 'Bound candidate', formula: 'C6H6',
        inchikey: 'AAAAAAAAAAAAAA-BBBBBBBBBB-C'},
    },
  };
}

test('lead scores project only for the same candidate and representative scan', () => {
  const data = dataLayer();
  const accepted = data.groupCompounds({compounds: [compound()]})[0].rep.match;
  assert.equal(accepted.cosine, 0.934);
  assert.equal(accepted.mass_ppm, -1.7);
  assert.deepEqual(Array.from(accepted.saved_lead_score_projection.fields), ['cosine', 'mass_ppm']);

  const wrongCandidate = data.groupCompounds({compounds: [compound({
    leadKey: 'CCCCCCCCCCCCCC-DDDDDDDDDD-E',
  })]})[0].rep.match;
  assert.equal(wrongCandidate.cosine, undefined);
  assert.equal(wrongCandidate.mass_ppm, undefined);

  const wrongScan = data.groupCompounds({compounds: [compound({leadIndex: 8})]})[0].rep.match;
  assert.equal(wrongScan.cosine, undefined);
  assert.equal(wrongScan.mass_ppm, undefined);
});

test('explicit spectral binding requires the same scan, full candidate and reference', () => {
  const data = dataLayer();
  const row = compound();
  row.lead_candidate.lib_spectrum_id = 'LIB-7';
  row.representative.match.lib_spectrum_id = 'LIB-7';
  row.representative.match.spectral_evidence_binding = {
    schema: 'same-scan-selected-reference/1', status: 'bound_saved_reference',
    spectrum_id: 'scan-7', result_index: 7,
    candidate_inchikey: 'AAAAAAAAAAAAAA-BBBBBBBBBB-C',
    lib_spectrum_id: 'LIB-7', saved_confidence: 'Low',
    metrics: {cosine: 0.812, mass_ppm: 2.4},
  };
  const accepted = data.groupCompounds({compounds: [row]})[0].rep.match;
  assert.equal(accepted.cosine, 0.812);
  assert.equal(accepted.mass_ppm, 2.4);
  assert.equal(accepted.native_spectrum_confidence, 'Low');

  row.representative.match.cosine = 0.999;
  row.representative.match.mass_ppm = 0.1;
  row.representative.match.spectral_evidence_binding.lib_spectrum_id = 'OTHER-LIB';
  const rejected = data.groupCompounds({compounds: [row]})[0].rep.match;
  assert.equal(rejected.cosine, undefined);
  assert.equal(rejected.mass_ppm, undefined);
  assert.equal(rejected.unbound_saved_display_values.cosine, 0.999);
});

test('explicit external roles are routed to review without guessing from names', () => {
  const data = dataLayer();
  assert.equal(data.explicitExternalReview({
    leadCandidate: {compound_role: 'external_contaminant'},
    rep: {match: {}}, compound: {plant_report: {}},
  }), true);
  assert.equal(data.explicitExternalReview({
    leadCandidate: {compound_role: 'pesticide_screening_target'},
    rep: {match: {}}, compound: {plant_report: {}},
  }), true);
  assert.equal(data.explicitExternalReview({
    leadCandidate: {compound_name: 'Unfamiliar compound', compound_role: 'unclassified'},
    rep: {match: {}}, compound: {plant_report: {}},
  }), false);
});
