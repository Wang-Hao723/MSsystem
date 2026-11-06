from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "web/static/js/candidate-comparison.js"
A, B = "A" * 14, "B" * 14


def render(payload, language="zh-CN", function="renderResult", extra=None):
    script = f"""
        const fs = require('fs'), vm = require('vm');
        global.window = {{MSSLanguage: {{getLanguage: () => {json.dumps(language)}}}}};
        vm.runInThisContext(fs.readFileSync({json.dumps(str(SOURCE))}, 'utf8'));
        const payload = {json.dumps(payload)}, extra = {json.dumps(extra)};
        const before = JSON.stringify(payload);
        const output = CandidateComparison.{function}(payload, extra);
        if (JSON.stringify(payload) !== before) throw Error('Rendering changed evidence');
        process.stdout.write(output);
    """
    return subprocess.run(["node", "-e", script], check=True, capture_output=True, text=True).stdout


@pytest.mark.parametrize("language", ["zh-CN", "en", "en-US"])
def test_default_display_preserves_native_order_not_failed_fragment_ranking(language):
    payload = {"candidate_display_order": [A, B], "panel": {"candidate_count": 2,
        "comparison_evidence": {"candidates": [{"connectivity": B, "review_rank": 1}, {"connectivity": A, "review_rank": 2}]}}}
    html = render(payload, language)
    assert html.index(A) < html.index(B)
    assert ("不替换正式 Top 1" if language == "zh-CN" else "does not replace production Top 1") in html


def test_null_is_not_displayed_as_measured_zero():
    html = render({"panel": {"comparison_evidence": {"candidates": [{"connectivity": A,
        "median_positive_fragment_contrast": None, "median_paired_cosine_delta": 0.0}]}}})
    assert "<td>—</td><td>0.0000</td>" in html


def test_all_40_candidates_render_and_injected_markup_is_escaped():
    candidates = [{"connectivity": f"candidate-{i}", "review_rank": i, "formulae": ["<script>alert(1)</script>"]} for i in range(40)]
    html = render({"panel": {"comparison_evidence": {"candidates": candidates}}})
    assert html.count('class="mono"') == 40
    assert "<script>alert" not in html
    assert "&lt;script&gt;alert" in html


def test_empty_ce_evidence_is_explained_not_fabricated():
    pair = {"a": A, "b": B, "observations": [], "ce_trajectories": []}
    html = render(pair, function="renderPair", extra={"cells": []})
    assert "暂无可比较轨迹" in html
    assert "缺字段、跨单位或单一 CE 不会被补成轨迹" in html


def test_ui_is_connected_to_real_entity_view_and_explicit_post_only():
    source = SOURCE.read_text()
    review = (ROOT / "web/static/js/evidence-review.js").read_text()
    index = (ROOT / "web/static/index.html").read_text()
    assert "CandidateComparison.mount" in review
    assert 'id="v16CandidateComparison"' in review
    assert index.index("candidate-comparison.js") < index.index("evidence-review.js")
    assert "API.post(`/tasks/${encodeURIComponent(taskId)}/candidate-comparison`" in source
    assert "button.addEventListener('click', async" in source
    assert "apply_rank_override" not in source
    assert "download.href =" in source and "/candidate-comparison/export?" in source
    assert "createObjectURL" not in source
