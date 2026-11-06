"""The interface must describe the actual ranking rule, not an intended rule."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "web/static/js/evidence-review.js"


def render(entity, language):
    script = f"""
      const fs = require('fs');
      const vm = require('vm');
      global.window = {{MSSLanguage: {{getLanguage: () => {json.dumps(language)}}}}};
      vm.runInThisContext(fs.readFileSync({json.dumps(str(SOURCE))}, 'utf8'));
      const entity = {json.dumps(entity)};
      const before = JSON.stringify(entity);
      const copy = EvidenceReview.replicateMissingnessSummary(entity);
      if (JSON.stringify(entity) !== before) throw Error('Display changed ranking');
      process.stdout.write(copy);
    """
    return subprocess.run(["node", "-e", script], check=True, capture_output=True, text=True).stdout


@pytest.mark.parametrize("language", ["zh-CN", "en", "en-US"])
def test_audit_explains_unchanged_production_ranking_and_counterfactual(language):
    entity = {"replicate_missingness_audit": {"status": "incomplete_truncated_topk_vectors",
        "counterfactual_winner_differs": True, "primary_ranking_changed": False}}
    output = render(entity, language)
    assert ("不用于替换正式排序" if language == "zh-CN" else "do not replace production ranking") in output
    assert ("首选候选发生变化" if language == "zh-CN" else "changes the base-ranking leader") in output


@pytest.mark.parametrize("imputed, phrase", [(True, "历史规则"), (False, "只使用已观测分数")])
def test_copy_matches_declared_rule(imputed, phrase):
    entity = {"candidates": [{"collision_energy_evidence": [{
        "truncated_topk_absence_imputed_as_zero": imputed}]}]}
    assert phrase in render(entity, "zh-CN")


def test_old_artifact_does_not_claim_a_known_missingness_rule():
    assert "未声明" in render({}, "zh-CN")
    assert "does not declare" in render({}, "en")


def test_rendered_evidence_uses_the_dynamic_rule_summary():
    source = SOURCE.read_text()
    assert "${esc(replicateMissingnessSummary(entity))}" in source
    assert "Top‑20 截断之外的缺失不按 0 分处罚" not in source
