import html as html_module
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INDEX_HTML = ROOT / "web" / "static" / "index.html"
UPLOAD_JS = ROOT / "web" / "static" / "js" / "upload.js"
APP_JS = ROOT / "web" / "static" / "js" / "app.js"
DASHBOARD_JS = ROOT / "web" / "static" / "js" / "dashboard.js"
LIBRARY_CENTER_JS = ROOT / "web" / "static" / "js" / "library-center.js"


def test_user_facing_branding_is_version_neutral_and_bilingual():
    markup = INDEX_HTML.read_text(encoding="utf-8")
    upload_js = UPLOAD_JS.read_text(encoding="utf-8")
    app_js = APP_JS.read_text(encoding="utf-8")
    dashboard_js = DASHBOARD_JS.read_text(encoding="utf-8")
    library_js = LIBRARY_CENTER_JS.read_text(encoding="utf-8")

    assert "<title>MSsystem · I.D.E.A platform</title>" in markup
    assert 'aria-label="I.D.E.A platform">I.D.E.A</div>' in markup
    assert 'data-label="数据与质量控制"' in markup
    assert 'data-label-en="Data &amp; Quality Control"' in markup
    assert '<span class="page-title-en">Data &amp; Quality Control</span>' in markup
    assert '<span class="page-title-en">Identification Results</span>' in markup
    assert '<span class="page-title-en">Reports &amp; Exports</span>' in markup
    assert '<span class="page-title-en">Library Management</span>' in markup
    assert "选择 Agilent .d / 父目录" in markup
    assert "选择 Thermo QE .raw / 批次目录" in markup

    visible_markup = re.sub(
        r"<(?:script|style)\b.*?</(?:script|style)>",
        "",
        markup,
        flags=re.IGNORECASE | re.DOTALL,
    )
    visible_text = html_module.unescape(
        re.sub(r"<[^>]+>", " ", visible_markup)
    )
    assert re.search(r"\bV7\b", visible_text) is None
    assert "macOS" not in visible_text

    assert "macOS" not in upload_js
    assert "MSsystem V7" not in app_js
    assert "MSsystem V7" not in dashboard_js
    assert "当前 V7 实例" not in library_js


def test_analysis_settings_are_visible_before_sample_selection():
    html = INDEX_HTML.read_text(encoding="utf-8")
    upload_js = UPLOAD_JS.read_text(encoding="utf-8")

    assert 'class="card intake-step-card" id="settingsCard"' in html
    assert 'class="card intake-step-card hidden" id="settingsCard"' not in html
    assert "settingsCard.classList.remove('hidden')" in upload_js
    assert "settingsCard.classList.toggle('hidden', !hasFiles)" not in upload_js
    assert 'id="analysisRunReadiness"' in html


def test_dashboard_task_open_uses_compound_first_sample_overview():
    dashboard_js = DASHBOARD_JS.read_text(encoding="utf-8")
    open_task = dashboard_js.split(
        "async function openTask(taskId, status)", 1
    )[1].split("function init()", 1)[0]

    assert "await SampleOverview.openTask(taskId)" in open_task
    assert "Results.loadForTask" not in open_task
    assert "setTimeout" not in open_task


def test_analysis_level_is_first_class_and_reports_effective_thresholds():
    html = INDEX_HTML.read_text(encoding="utf-8")

    assert 'id="analysisLevel"' in html
    assert 'value="routine" selected' in html
    assert 'value="sensitive"' in html
    assert 'value="report"' in html
    assert 'value="custom"' in html
    assert "const ANALYSIS_LEVELS" in html
    assert "cosTh: 0.30" in html
    assert "finalTh: 0.25" in html
    assert "peakProcessingPreset: 'library_match'" in html
    assert "strictEvidenceGate: true" in html

    for element_id in (
        "analysisSummaryLevel",
        "analysisSummaryCosine",
        "analysisSummaryFinal",
        "analysisSummaryMass",
        "analysisSummaryPeaks",
        "analysisSummaryStrict",
        "advancedSettingsState",
    ):
        assert f'id="{element_id}"' in html


def test_instrument_preset_does_not_silently_change_evidence_level():
    html = INDEX_HTML.read_text(encoding="utf-8")
    instrument_block = html.split("const INSTRUMENT_PRESETS = {", 1)[1].split(
        "let applyingAnalysisPreset", 1
    )[0]
    apply_instrument = html.split("function applyInstrumentPreset", 1)[1].split(
        "function markAnalysisLevelCustom", 1
    )[0]

    assert "mzDa" in instrument_block
    assert "mzPpm" in instrument_block
    assert "cosTh" not in instrument_block
    assert "finalTh" not in instrument_block
    assert "setAnalysisControlValue('mzDa'" in apply_instrument
    assert "setAnalysisControlValue('mzPpm'" in apply_instrument
    assert "setAnalysisControlValue('cosTh'" not in apply_instrument
    assert "setAnalysisControlValue('finalTh'" not in apply_instrument
