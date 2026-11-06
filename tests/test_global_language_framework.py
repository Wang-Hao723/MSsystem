"""Regression coverage for the V16 global interface-language framework."""

from __future__ import annotations

import json
import re
import subprocess
from html.parser import HTMLParser
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "web" / "static"
INDEX = STATIC / "index.html"
CATALOG = STATIC / "js" / "i18n-catalog.js"
RESULTS_CATALOG = STATIC / "js" / "i18n-results-catalog.js"
REPORT_CATALOG = STATIC / "js" / "i18n-report-catalog.js"
EVIDENCE_CATALOG = STATIC / "js" / "i18n-evidence-catalog.js"
RUNTIME = STATIC / "js" / "i18n.js"
RESULT_PRESENTERS = (
    STATIC / "js" / "results.js",
    STATIC / "js" / "sample-overview.js",
    STATIC / "js" / "evidence-review.js",
    STATIC / "js" / "report-generator.js",
)


def _node(script: str) -> str:
    completed = subprocess.run(
        ["node", "-e", script],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


class _InterfaceTextParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.skip_depth = 0
        self.text: list[str] = []
        self.attributes: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in {"script", "style"}:
            self.skip_depth += 1
        for name, value in attrs:
            if name in {"placeholder", "title", "aria-label", "alt", "label", "data-empty-text"} and value:
                self.attributes.append(" ".join(value.split()))

    def handle_endtag(self, tag: str) -> None:
        if tag in {"script", "style"} and self.skip_depth:
            self.skip_depth -= 1

    def handle_data(self, data: str) -> None:
        if not self.skip_depth:
            value = " ".join(data.split())
            if value:
                self.text.append(value)


def _catalog() -> dict[str, str]:
    script = f"""
      global.window = global;
      require({json.dumps(str(CATALOG))});
      process.stdout.write(JSON.stringify(
        global.MSS_I18N_CATALOG.legacy['en-US']
      ));
    """
    return json.loads(_node(script))


def test_language_assets_load_before_application_modules() -> None:
    page = INDEX.read_text(encoding="utf-8")
    catalog_position = page.index('/static/js/i18n-catalog.js?v=16.2')
    results_catalog_position = page.index('/static/js/i18n-results-catalog.js?v=15.5')
    report_catalog_position = page.index('/static/js/i18n-report-catalog.js?v=15.8')
    evidence_catalog_position = page.index('/static/js/i18n-evidence-catalog.js?v=15.7')
    runtime_position = page.index('/static/js/i18n.js?v=15.10')
    api_position = page.index('/static/js/api.js?v=3.4')
    assert (
        catalog_position
        < results_catalog_position
        < report_catalog_position
        < evidence_catalog_position
        < runtime_position
        < api_position
    )
    assert 'data-language-switch' in page
    assert 'data-language-option="zh-CN"' in page
    assert 'data-language-option="en-US"' in page


def test_all_static_chinese_interface_copy_has_a_curated_english_route() -> None:
    parser = _InterfaceTextParser()
    parser.feed(INDEX.read_text(encoding="utf-8"))
    han = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")
    source = {value for value in parser.text + parser.attributes if han.search(value)}
    translations = _catalog()

    # The group label uses a semantic key; the native language name is kept as
    # "中文" and excluded from translation with data-i18n-skip.
    semantic_control_copy = {"界面语言", "中文"}
    assert source - set(translations) == semantic_control_copy
    assert len(translations) >= 590
    assert not {
        source_text: english
        for source_text, english in translations.items()
        if han.search(english)
    }


def test_locale_resolution_translation_and_persistence_are_deterministic() -> None:
    script = f"""
      const fs = require('fs');
      const vm = require('vm');
      global.window = global;
      Object.defineProperty(global, 'navigator', {{
        value: {{ languages: ['zh-CN', 'en-US'], language: 'zh-CN' }},
        configurable: true,
      }});
      const values = new Map([['mssystem.locale.v1', 'zh-CN']]);
      global.localStorage = {{
        getItem: key => values.get(key) || null,
        setItem: (key, value) => values.set(key, value),
        removeItem: key => values.delete(key),
      }};
      global.location = {{ href: 'http://127.0.0.1:8516/app?lang=en#tasks' }};
      let replaced = '';
      global.history = {{
        state: null,
        replaceState: (_state, _title, url) => {{ replaced = url; }},
      }};
      vm.runInThisContext(fs.readFileSync({json.dumps(str(CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(RUNTIME))}, 'utf8'));

      const initial = {{
        locale: MSSLanguage.getLocale(),
        source: MSSLanguage.getSource(),
        title: MSSLanguage.t('分析工作台'),
        count: MSSLanguage.translateText('17 个样本'),
        userData: MSSLanguage.translateText('F001_独活寄生汤_targeted_20260817'),
      }};
      MSSLanguage.setLocale('zh-CN');
      const manual = {{
        locale: MSSLanguage.getLocale(),
        saved: values.get(MSSLanguage.STORAGE_KEY),
        url: replaced,
        title: MSSLanguage.t('分析工作台'),
      }};
      MSSLanguage.setLocale('auto');
      const automatic = {{
        locale: MSSLanguage.getLocale(),
        source: MSSLanguage.getSource(),
        saved: values.has(MSSLanguage.STORAGE_KEY),
        url: replaced,
      }};
      process.stdout.write(JSON.stringify({{ initial, manual, automatic }}));
    """
    result = json.loads(_node(script))
    assert result["initial"] == {
        "locale": "en-US",
        "source": "url",
        "title": "Analysis Workspace",
        "count": "17 samples",
        "userData": "F001_独活寄生汤_targeted_20260817",
    }
    assert result["manual"] == {
        "locale": "zh-CN",
        "saved": "zh-CN",
        "url": "/app?lang=zh#tasks",
        "title": "分析工作台",
    }
    assert result["automatic"] == {
        "locale": "zh-CN",
        "source": "browser",
        "saved": False,
        "url": "/app#tasks",
    }


def test_dynamic_translation_hooks_and_safe_hash_navigation_remain_wired() -> None:
    runtime = RUNTIME.read_text(encoding="utf-8")
    app = (STATIC / "js" / "app.js").read_text(encoding="utf-8")
    export = (STATIC / "js" / "export.js").read_text(encoding="utf-8")
    report = (STATIC / "js" / "report-generator.js").read_text(encoding="utf-8")
    page = INDEX.read_text(encoding="utf-8")
    assert "MutationObserver" in runtime
    assert "mssystem:languagechange" in runtime
    assert "mssystem:languagechange" in app
    assert "`${url.pathname}${url.search}${url.hash}`" in app
    assert "interfaceLanguage(language)" in export
    assert "downloadHTML(formatType, language = null)" in export
    assert "downloadPDF(language = null)" in export
    assert 'value="html_standard|auto"' in page
    assert "window.MSSLanguage?.getLanguage() || 'zh'" in report


def test_scientific_names_are_protected_from_interface_translation() -> None:
    runtime = RUNTIME.read_text(encoding="utf-8")
    assert runtime.count("'[data-scientific-value]'") >= 2

    presenters = {
        path.name: path.read_text(encoding="utf-8") for path in RESULT_PRESENTERS
    }
    assert "_nameWithCn" not in presenters["results.js"]
    assert "m.compound_name_cn" not in presenters["results.js"]
    assert "rawMatch.compound_name_cn" not in presenters["sample-overview.js"]
    assert "match.compound_name_cn || match.compound_name" not in presenters["sample-overview.js"]
    assert "if (match?.compound_name_cn)" not in presenters["evidence-review.js"]
    assert "'compound_name_cn'" not in presenters["evidence-review.js"]
    assert "match.compound_name_cn\n      ?" not in presenters["report-generator.js"]
    for name, source in presenters.items():
        assert "data-scientific-value" in source, name
    assert "ConfidenceGuide.load('zh')" not in presenters["sample-overview.js"]


def test_result_catalog_translates_dynamic_copy_without_rewriting_names() -> None:
    script = f"""
      const fs = require('fs');
      const vm = require('vm');
      global.window = global;
      Object.defineProperty(global, 'navigator', {{
        value: {{ languages: ['en-US'], language: 'en-US' }}, configurable: true,
      }});
      global.localStorage = {{
        getItem: () => null, setItem: () => {{}}, removeItem: () => {{}},
      }};
      global.location = {{ href: 'http://127.0.0.1:8516/?lang=en#results' }};
      global.history = {{ state: null, replaceState: () => {{}} }};
      vm.runInThisContext(fs.readFileSync({json.dumps(str(CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(RESULTS_CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(RUNTIME))}, 'utf8'));
      process.stdout.write(JSON.stringify({{
        heading: MSSLanguage.translateText('样本结果摘要'),
        marker: MSSLanguage.translateText(
          '高可信单一二维连接结构检出：Acanthoside B · RT 3.392 min · m/z 579.20947 · 1 张局部 MS2'
        ),
        fragment: MSSLanguage.translateText(
          '417.1556 → 402.1339，推测丢失／转化 甲基自由基'
        ),
        summaryResidual: MSSLanguage.translateText(
          '个结果记录。样品支持的主要结果为'
        ),
        tentativeClass: MSSLanguage.translateText(
          '类别线索（待核查）：lipid'
        ),
        unavailableFragmentation: MSSLanguage.translateText(
          '未能推理：当前证据不足以形成候选特异碎裂途径'
        ),
        isomerConfidence: MSSLanguage.translateText(
          'High（异构体组） / High（完整结构身份）'
        ),
        databaseName: MSSLanguage.translateText('十二烷基二甲基氧化胺'),
      }}));
    """
    result = json.loads(_node(script))
    assert result == {
        "heading": "Sample result summary",
        "marker": (
            "High-confidence single 2-D connectivity detection: Acanthoside B "
            "· RT 3.392 min · m/z 579.20947 · 1 local MS2 scan"
        ),
        "fragment": (
            "417.1556 → 402.1339, proposed loss/transformation: methyl radical"
        ),
        "summaryResidual": (
            " result records. The sample-supported primary results total "
        ),
        "tentativeClass": "Class lead requiring review: lipid",
        "unavailableFragmentation": (
            "Unavailable: current evidence is insufficient to propose a "
            "candidate-specific fragmentation pathway"
        ),
        "isomerConfidence": (
            "High (isomer group) / High (full-structure identity)"
        ),
        "databaseName": "十二烷基二甲基氧化胺",
    }


def test_dynamic_result_presenter_emits_english_copy_without_waiting_for_dom_rewrite() -> None:
    overview = (STATIC / "js" / "sample-overview.js").read_text(encoding="utf-8")
    report = (STATIC / "js" / "report-generator.js").read_text(encoding="utf-8")

    assert "The sample-supported primary results comprise" in overview
    assert "Class lead requiring review:" in overview
    assert (
        "Unavailable: current evidence is insufficient to propose a "
        "candidate-specific fragmentation pathway"
    ) in overview
    assert "(isomer group) /" in overview
    assert "primary-result acquisition locations" in report


def test_report_catalog_localizes_evidence_prose_and_preserves_captured_values() -> None:
    script = f"""
      const fs = require('fs');
      const vm = require('vm');
      global.window = global;
      Object.defineProperty(global, 'navigator', {{
        value: {{ languages: ['en-US'], language: 'en-US' }}, configurable: true,
      }});
      global.localStorage = {{
        getItem: () => null, setItem: () => {{}}, removeItem: () => {{}},
      }};
      global.location = {{ href: 'http://127.0.0.1:8516/?lang=en#report' }};
      global.history = {{ state: null, replaceState: () => {{}} }};
      vm.runInThisContext(fs.readFileSync({json.dumps(str(CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(RESULTS_CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(REPORT_CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(RUNTIME))}, 'utf8'));
      process.stdout.write(JSON.stringify({{
        heading: MSSLanguage.translateText('LC–MS/MS 样本综合分析报告'),
        candidate: MSSLanguage.translateText('本针 Top1 排序先导：Acanthoside B'),
        ion: MSSLanguage.translateText('C28H36O13 · 样品离子形式 [M-H]-'),
        evidence: MSSLanguage.translateText(
          '观测支持：中性丢失（4）、规则覆盖（1）；未观察 2 维；不适用 1 维'
        ),
        ranked: MSSLanguage.translateText(
          '综合排序分 0.820；规则证据、分类、生源和碎裂解释不参与身份排序。'
        ),
        arbitraryName: MSSLanguage.translateText('十二烷基二甲基氧化胺'),
      }}));
    """
    result = json.loads(_node(script))
    assert result == {
        "heading": "Integrated LC–MS/MS Sample Analysis Report",
        "candidate": "Top-ranked candidate in this injection: Acanthoside B",
        "ion": "C28H36O13 · Sample ion form [M-H]-",
        "evidence": (
            "Observed support: neutral losses (4), rule coverage (1); "
            "2 dimensions not observed; 1 dimension not applicable"
        ),
        "ranked": (
            "Composite ranking score 0.820; rule evidence, classification, "
            "biosource context and fragmentation interpretation do not "
            "contribute to identity ranking."
        ),
        "arbitraryName": "十二烷基二甲基氧化胺",
    }


def test_evidence_catalog_localizes_audit_copy_without_touching_scientific_values() -> None:
    script = f"""
      const fs = require('fs');
      const vm = require('vm');
      global.window = global;
      Object.defineProperty(global, 'navigator', {{
        value: {{ languages: ['en-US'], language: 'en-US' }}, configurable: true,
      }});
      global.localStorage = {{
        getItem: () => null, setItem: () => {{}}, removeItem: () => {{}},
      }};
      global.location = {{ href: 'http://127.0.0.1:8516/?lang=en#review' }};
      global.history = {{ state: null, replaceState: () => {{}} }};
      vm.runInThisContext(fs.readFileSync({json.dumps(str(CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(RESULTS_CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(REPORT_CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(EVIDENCE_CATALOG))}, 'utf8'));
      vm.runInThisContext(fs.readFileSync({json.dumps(str(RUNTIME))}, 'utf8'));
      process.stdout.write(JSON.stringify({{
        panel: MSSLanguage.translateText('V15 实体基础分与最终选择'),
        locator: MSSLanguage.translateText(
          'AutoMS-S2-1-3-neg.mzML · spectrum 203544 · result_index 144。此扫描只用于原始证据定位；它的旧 Top‑1、规则证据和经验通路不能覆盖旗舰实体决策。'
        ),
        release: MSSLanguage.translateText(
          'passed；失败项 无。该门只改变发布范围，不改变 Top 1、Top 5 或实体分。'
        ),
        runtime: MSSLanguage.translateText(
          '0 条；30 谱试运行规则 未启用；研究 catalog 未计算。'
        ),
        representative: MSSLanguage.translateText('代表谱 High 净增'),
        loading: MSSLanguage.translateText('正在核对各模块的执行状态与适用边界…'),
        sugarLoss: MSSLanguage.translateText('−162.05 (己糖)→'),
        name: MSSLanguage.translateText('Acanthoside B'),
      }}));
    """
    result = json.loads(_node(script))
    assert result == {
        "panel": "V15 entity base score and final selection",
        "locator": (
            "AutoMS-S2-1-3-neg.mzML · spectrum 203544 · result_index 144. "
            "This scan is used only for raw-evidence location; its legacy Top-1 "
            "result, rule evidence, and empirical pathway cannot override the "
            "flagship entity decision."
        ),
        "release": (
            "Passed; failed checks: none. This gate changes release scope only "
            "and does not alter the Top 1, Top 5, or entity score."
        ),
        "runtime": (
            "0 rules; 30-spectrum pilot rules: disabled; research catalog: "
            "not evaluated."
        ),
        "representative": "Net increase in representative-spectrum High results",
        "loading": "Checking module execution status and applicability…",
        "sugarLoss": "−162.05 (hexose)→",
        "name": "Acanthoside B",
    }


def test_language_javascript_passes_node_syntax_check() -> None:
    for path in (
        CATALOG,
        RESULTS_CATALOG,
        REPORT_CATALOG,
        EVIDENCE_CATALOG,
        RUNTIME,
        STATIC / "js" / "app.js",
        STATIC / "js" / "dashboard.js",
        STATIC / "js" / "export.js",
        STATIC / "js" / "report-generator.js",
        STATIC / "js" / "tasks.js",
        STATIC / "js" / "upload.js",
    ):
        subprocess.run(["node", "--check", path], cwd=ROOT, check=True)
