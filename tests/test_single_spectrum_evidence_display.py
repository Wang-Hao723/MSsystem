"""Library evidence must remain visible when an entity lacks multi-energy data."""
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'web/static/js/sample-overview.js'


def run_js(body):
    setup = f"""
const fs = require('fs'); const vm = require('vm');
const window = {{MSSLanguage: {{getLanguage: () => 'zh'}}}};
vm.runInThisContext(fs.readFileSync({json.dumps(str(SOURCE))}, 'utf8'));
"""
    result = subprocess.run(['node', '-e', setup + body], cwd=ROOT,
                            capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def test_high_library_evidence_survives_single_energy_entity_medium():
    result = run_js("""
const entity = {candidate_count: 2, distinct_collision_energy_bins: 1,
  decision: {status: 'selected_isomer_group', blockers: ['insufficient_distinct_collision_energies', 'competitive_same_formula_connectivities']}};
const bundle = {task: {compound_confidence_summary: {
  representative_high_reportable_scope_unique_compounds: 66,
  representative_2d_connectivity_high_unique_compounds: 36,
  representative_isomer_group_high_unique_compounds: 30}},
  compounds: {metadata: {high_confidence_scope_group_count: 0}},
  flagship: {entities_total: 1, has_more: false, entities: [entity]}};
const before = JSON.stringify(bundle);
console.log(JSON.stringify({axes: SampleResultData.confidenceAxes(bundle),
  zh: SampleResultData.confidenceAxesHtml(bundle, 'zh'),
  en: SampleResultData.confidenceAxesHtml(bundle, 'en'),
  reason: SampleResultData.entityDecisionNote(entity, 'zh'),
  unchanged: before === JSON.stringify(bundle)}));
""")
    assert result['axes']['libraryHigh'] == 66
    assert result['axes']['entityHigh'] == 0
    assert result['axes']['singleEnergyEntities'] == 1
    assert result['axes']['energyBlockedEntities'] == 1
    assert result['unchanged']
    assert '单谱谱库 High：66' in result['zh']
    assert '36 个单一连接候选、30 个异构体组候选' in result['zh']
    assert 'Single-spectrum library High' in result['en']
    assert '多能量条件' in result['reason']
    assert '多个同分子式连接结构候选' in result['reason']


def test_missing_and_partial_data_do_not_become_zero_or_full_cohort_claims():
    result = run_js("""
const empty = {task: {}, compounds: {metadata: {}}};
const partial = {task: {compound_confidence_summary: {high_reportable_scope_unique_compounds: 7}},
 flagship: {entities_total: 100, has_more: true,
 entities: [{candidate_count: 1, distinct_collision_energy_bins: 1,
 decision: {blockers: ['insufficient_distinct_collision_energies']}}]}};
console.log(JSON.stringify({empty: SampleResultData.confidenceAxes(empty),
 emptyHtml: SampleResultData.confidenceAxesHtml(empty), partial: SampleResultData.confidenceAxes(partial),
 partialHtml: SampleResultData.confidenceAxesHtml(partial)}));
""")
    assert result['empty']['libraryHigh'] is None
    assert result['empty']['entityHigh'] is None
    assert result['emptyHtml'] == ''
    assert result['partial']['libraryHigh'] == 7
    assert result['partial']['matchedEntities'] is None
    assert result['partial']['energyBlockedEntities'] is None
    assert '个有候选前体实体' not in result['partialHtml']


def test_validated_same_candidate_high_explains_the_actual_support():
    result = run_js("""
const entity = {decision: {blockers: ['insufficient_distinct_collision_energies']}};
const support = {tier: 'High', scope: 'isomer_group', supporting_spectra: [{result_index: 3}]};
console.log(JSON.stringify({zh: SampleResultData.entityDecisionNote(entity, 'zh', support),
 en: SampleResultData.entityDecisionNote(entity, 'en', support)}));
""")
    assert '异构体组 High' in result['zh']
    assert '组内单体仍未区分' in result['zh']
    assert 'individual isomers remain unresolved' in result['en']
    assert '未通过' not in result['zh']


def test_inconsistent_summary_subtypes_are_not_presented_as_a_valid_breakdown():
    result = run_js("""
const bundle = {task: {compound_confidence_summary: {
 high_reportable_scope_unique_compounds: 20,
 representative_2d_connectivity_high_unique_compounds: 3,
 representative_isomer_group_high_unique_compounds: 2}}};
console.log(JSON.stringify({axes: SampleResultData.confidenceAxes(bundle), html: SampleResultData.confidenceAxesHtml(bundle)}));
""")
    assert result['axes']['libraryHigh'] == 20
    assert not result['axes']['breakdownConsistent']
    assert '3 个单一连接候选' not in result['html']


def test_representative_selection_preserves_score_order_and_separates_scope():
    result = run_js("""
function row(id, structure, score, fields) { return {spectrum_id: id,
 match: {smiles: structure, rank_score: score, final_score: score, ...fields}}; }
const rows = [
 row('a-low-score-high', 'CCO', 0.8, {confidence: 'High'}),
 row('a-representative-medium', 'CCO', 0.9, {confidence: 'Medium'}),
 row('b-connectivity-high', 'CCN', 0.8, {confidence: 'Medium', two_dimensional_connectivity_confidence: 'High'}),
 row('c-isomer-high', 'CCCO', 0.7, {confidence: 'Medium', identification_scope: 'isomer_group', isomer_group_confidence: 'High'}),
 row('d-class-only-high', 'CCCC', 0.7, {confidence: 'Medium', confidence_class: 'High'}),
 row('e-invalid', '', 0.9, {confidence: 'High'}),
 row('f-retired-fallback', 'C', 0.9, {confidence: 'High', dataset_source: 'iceberg_fallback'})];
const before = JSON.stringify(rows);
console.log(JSON.stringify({ids: SampleResultData.singleSpectrumHighRepresentatives(rows).map(r => r.spectrum_id),
 unchanged: before === JSON.stringify(rows)}));
""")
    assert result['ids'] == ['b-connectivity-high', 'c-isomer-high']
    assert result['unchanged']
