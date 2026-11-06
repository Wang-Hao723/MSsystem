import json
from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
OVERVIEW = ROOT / "web/static/js/sample-overview.js"


def test_result_table_merges_repeated_features_by_candidate_identity():
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required for the result-view regression test")
    program = f"""
const fs = require('fs'); const vm = require('vm');
global.window = {{MSSLanguage: {{getLanguage: () => 'zh'}}}};
vm.runInThisContext(fs.readFileSync({json.dumps(str(OVERVIEW))}, 'utf8'));
const make = (name, formula, inchikey, feature, score) => ({{
  compound_id: `compound:${{inchikey.slice(0, 14)}}`,
  result_id: `result:${{feature}}`, feature_id: `lcms:${{feature}}`,
  identity_key: `std2d:${{inchikey.slice(0, 14)}}|formula:${{formula}}`,
  display_name: name, formula,
  lead_candidate: {{compound_name:name, formula, inchikey, rank_score:score}},
  candidate_members: [{{compound_name:name, formula, inchikey, rank_score:score}}],
  candidate_member_count: 1,
  conclusion: {{key:'candidate', label:'中等证据候选', counted_as_sample_compound:true, counted_as_plant_candidate:true}},
  representative: {{spectrum_id:`scan-${{feature}}`, retention_time:Number(feature), retention_time_unit:'minute', match:{{compound_name:name, formula, inchikey, rank_score:score, cosine:score}}}},
  spectra_count: 1, peak_cluster_count: 1,
  occurrences: [{{spectrum_id:`scan-${{feature}}`, retention_time:Number(feature), retention_time_unit:'minute'}}],
  sample_evidence: {{sample_supported_spectra:1, category_counts:{{sample_supported:1}}}},
}});
const makeIsomer = (feature, score) => ({{
  compound_id: `compound:isomer-${{feature}}`,
  result_id: `result:isomer-${{feature}}`, feature_id: `lcms:${{feature}}`,
  identity_key: `isomer-event:lcms:${{feature}}|formula:C19H23NO4`,
  display_name: 'Reticuline', formula: 'C19H23NO4',
  lead_candidate: {{compound_name:'Reticuline', formula:'C19H23NO4', inchikey:'BHLYRWXGMIUIHG-OAHLLOKOSA-N', rank_score:score}},
  candidate_members: [
    {{compound_name:'Reticuline', formula:'C19H23NO4', inchikey:'BHLYRWXGMIUIHG-OAHLLOKOSA-N', rank_score:score}},
    {{compound_name:'Sinomenine', formula:'C19H23NO4', inchikey:'INYYVPJSBIVGPH-QHRIQVFBSA-N', rank_score:score - 0.01}},
  ],
  candidate_member_count: 2, candidate_names: ['Reticuline', 'Sinomenine'],
  conclusion: {{key:'isomer', label:'高可信异构体组检出', counted_as_sample_compound:true, counted_as_plant_candidate:true, high_confidence_group_detected:true}},
  representative: {{spectrum_id:`scan-isomer-${{feature}}`, retention_time:Number(feature), retention_time_unit:'minute', match:{{compound_name:'Reticuline', formula:'C19H23NO4', rank_score:score, cosine:score}}}},
  spectra_count: 1, peak_cluster_count: 1,
  occurrences: [{{spectrum_id:`scan-isomer-${{feature}}`, retention_time:Number(feature), retention_time_unit:'minute'}}],
  sample_evidence: {{sample_supported_spectra:1, category_counts:{{sample_supported:1}}}},
}});
const payload = {{compounds:[
  make('acanthoside B','C28H36O13','WEKCEGQSIIQPAQ-IRBNZIFYSA-N','1',0.80),
  make('acanthoside B','C28H36O13','WEKCEGQSIIQPAQ-IRBNZIFYSA-N','2',0.90),
  make('Sinomenine','C19H23NO4','INYYVPJSBIVGPH-QHRIQVFBSA-N','3',0.60),
  make('Sinomenine','C19H23NO4','INYYVPJSBIVGPH-QHRIQVFBSA-N','4',0.70),
  make('Sinomenine','C19H23NO4','INYYVPJSBIVGPH-QHRIQVFBSA-N','5',0.65),
  makeIsomer('6',0.72),
  makeIsomer('7',0.74),
]}};
const raw = SampleResultData.groupCompounds(payload);
const merged = SampleResultData.deduplicateCompoundGroups(raw);
console.log(JSON.stringify(merged.map(group => ({{
  name: group.leadCandidate.compound_name,
  featureResults: group.featureResultCount,
  spectra: group.spectra,
  peakClusters: group.peakClusters,
  rows: group.rows.length,
}}))));
"""
    result = subprocess.run(
        [node, "-e", program], cwd=ROOT, check=True,
        capture_output=True, text=True,
    )
    assert json.loads(result.stdout) == [
        {
            "name": "Sinomenine",
            "featureResults": 3,
            "spectra": 3,
            "peakClusters": 3,
            "rows": 3,
        },
        {
            "name": "Reticuline",
            "featureResults": 2,
            "spectra": 2,
            "peakClusters": 2,
            "rows": 2,
        },
        {
            "name": "acanthoside B",
            "featureResults": 2,
            "spectra": 2,
            "peakClusters": 2,
            "rows": 2,
        },
    ]
