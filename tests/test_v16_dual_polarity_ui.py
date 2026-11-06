"""Source-level browser contract for the V16 dual-polarity shadow card."""

from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
OVERVIEW_PATH = ROOT / "web/static/js/sample-overview.js"
STYLE_PATH = ROOT / "web/static/css/workspace-v7.css"
OVERVIEW = OVERVIEW_PATH.read_text(encoding="utf-8")
STYLE = STYLE_PATH.read_text(encoding="utf-8")


def _section(start: str, end: str) -> str:
    start_index = OVERVIEW.index(start)
    end_index = OVERVIEW.index(end, start_index)
    return OVERVIEW[start_index:end_index]


def test_card_is_injected_by_the_sample_overview_without_static_html_dependency():
    render = _section("function renderOverview(bundle)", "function matchesAnalysisFilter")

    assert "${dualPolarityShadowSection(task)}" in render
    assert 'id="v16DualPolarityShadowCard"' in OVERVIEW
    assert "V16 双极性影子证据" in OVERVIEW
    assert "V16 dual-polarity shadow evidence" in OVERVIEW
    assert "SHADOW · READ ONLY" in OVERVIEW


def test_eligible_tasks_require_completed_exact_opposite_polarity_without_filename_pairing():
    exact_mode = _section(
        "function exactDualPolarity(task)",
        "function isCompletedDualPolarityTask",
    )
    completed = _section(
        "function isCompletedDualPolarityTask(task)",
        "function eligibleOppositePolarityTasks",
    )
    eligible = _section(
        "function eligibleOppositePolarityTasks(currentTask)",
        "function unresolvedEventKind",
    )

    assert "task?.ion_mode" in exact_mode
    assert "params" not in exact_mode
    assert "mode === 'POS' || mode === 'NEG'" in exact_mode
    assert "=== 'completed'" in completed
    assert "oppositeMode" in eligible
    assert "exactDualPolarity(task) === oppositeMode" in eligible
    assert "String(task?.id || '') !== String(currentTask?.id || '')" in eligible
    assert "filename" not in eligible.lower()
    assert "API.get('/tasks?status=completed')" in OVERVIEW


def test_pairing_starts_unselected_and_requires_two_manual_confirmations():
    card = _section(
        "function dualPolarityShadowSection(task)",
        "function updateDualPolarityReadiness",
    )
    readiness = _section(
        "function updateDualPolarityReadiness()",
        "async function runDualPolarityShadow",
    )

    assert '<option value="">${optionLabel}</option>' in card
    assert "Nothing is preselected" in card
    assert "文件名相似不构成配对证据" in card
    assert 'id="v16DualSameSampleConfirmed" type="checkbox"' in card
    assert 'id="v16DualSameLcMethodConfirmed" type="checkbox"' in card
    assert 'id="v16DualPolarityRun"' in card
    assert " disabled " in card
    assert "sameSample && sameLcMethod" in readiness
    assert "!!selectedTask" in readiness


def test_post_payload_maps_manual_selection_to_explicit_pos_and_neg_task_ids():
    submit = _section(
        "async function runDualPolarityShadow()",
        "function referenceStandardQCSection",
    )

    assert submit.index("if (!isCompletedDualPolarityTask") < submit.index(
        "await API.post"
    )
    assert submit.index("if (!sameSample || !sameLcMethod)") < submit.index(
        "await API.post"
    )
    assert "positive_task_id: currentMode === 'POS'" in submit
    assert "negative_task_id: currentMode === 'NEG'" in submit
    assert "same_sample_confirmed: true" in submit
    assert "same_lc_method_confirmed: true" in submit
    assert "pair_label:" in submit
    assert "'/workflows/dual-polarity-shadow'" in submit


def test_result_exposes_required_shadow_counts_and_production_invariants():
    result = _section(
        "function dualPolarityShadowResult(response)",
        "function dualPolarityShadowSection",
    )

    assert "accepted_neutral_entity_count" in result
    assert "ambiguous_anchor_count" in result
    assert "not_evaluable_anchor_count" in result
    assert "dual_polarity_supported_connectivity_count" in result
    for field in (
        "production_top1_changed",
        "production_top5_changed",
        "production_score_changed",
        "production_decision_changed",
        "production_confidence_changed",
        "production_compound_count_changed",
    ):
        assert field in result
    assert "不会改变 V15 Top1/Top5、分数、决策、置信度或化合物计数" in result
    assert "跨极性 MS2 相似度仅作辅助注释" in result


def test_shadow_card_has_scoped_responsive_styles():
    for selector in (
        ".v16-dual-card",
        ".v16-dual-form",
        ".v16-dual-confirmations",
        ".v16-dual-metrics",
        ".v16-dual-invariant.is-verified",
        ".v16-dual-invariant.is-incomplete",
    ):
        assert selector in STYLE


def test_sample_overview_javascript_parses_with_node():
    completed = subprocess.run(
        ["node", "--check", str(OVERVIEW_PATH)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
