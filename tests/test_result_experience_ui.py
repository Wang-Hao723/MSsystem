from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INDEX = (ROOT / "web/static/index.html").read_text(encoding="utf-8")
OVERVIEW = (ROOT / "web/static/js/sample-overview.js").read_text(encoding="utf-8")
REVIEW = (ROOT / "web/static/js/evidence-review.js").read_text(encoding="utf-8")
RESULTS = (ROOT / "web/static/js/results.js").read_text(encoding="utf-8")
MODULE_EFFECTS = (ROOT / "web/static/js/module-effects.js").read_text(encoding="utf-8")
REPORT = (ROOT / "web/static/js/report-generator.js").read_text(encoding="utf-8")
STYLE = (ROOT / "web/static/css/style.css").read_text(encoding="utf-8")
ACQUIRE = (ROOT / "web/static/js/eviacquire.js").read_text(encoding="utf-8")
TASKS = (ROOT / "web/static/js/tasks.js").read_text(encoding="utf-8")


def _page(markup: str, page_id: str) -> str:
    start = markup.index(f'id="page-{page_id}"')
    next_page = markup.find('<section class="page"', start + 1)
    return markup[start:] if next_page < 0 else markup[start:next_page]


def test_identification_page_keeps_compounds_and_core_analysis_in_one_workflow():
    result_page = _page(INDEX, "results")
    review_page = _page(INDEX, "review")

    assert 'id="sampleOverviewArea"' in result_page
    assert 'id="sampleIntegratedAnalysis"' not in result_page
    assert "去重化合物结果" in OVERVIEW
    assert "鉴定结果概览" in OVERVIEW
    assert "分类 / 生源 / 碎裂 / 规则证据" in OVERVIEW
    assert "分类分析" in OVERVIEW
    assert "生源背景" in OVERVIEW
    assert "规则证据评估" in OVERVIEW
    assert "推测碎裂途径" in OVERVIEW
    assert "biosynthetic_pathways" in OVERVIEW
    assert "未匹配谱图没有被删除" in OVERVIEW
    assert "进入未匹配谱图工作台" in OVERVIEW
    assert "setAnalysisFilter" in OVERVIEW
    assert "谱库证据级别" in OVERVIEW
    assert "样本／空白" in OVERVIEW
    assert "High 级鉴定组" in OVERVIEW
    assert "High 单一二维连接" in OVERVIEW
    assert "High 异构体组" in OVERVIEW
    assert "已确证存在" in OVERVIEW
    assert "类别线索（待核查）：" in OVERVIEW
    assert "类别待核查（模型提示：" in RESULTS
    assert "查看证据详情" in OVERVIEW
    assert "Top1 排序先导" in OVERVIEW
    assert "查看全部" in OVERVIEW
    assert "conciseDisplayName" in OVERVIEW
    assert "其余名称不再挤在概览中" in OVERVIEW
    assert "同一 RT 簇另有" in REPORT
    assert "/replicate-consistency" in OVERVIEW
    assert "重复稳健 High" in OVERVIEW
    assert 'id="moduleEffectPanel"' in review_page
    assert 'id="eviAcquirePanel"' in review_page
    assert 'id="resultsArea"' in review_page
    assert "Results.showDetailFromSummary" in REVIEW
    assert "/compounds?limit=2000" in OVERVIEW
    assert "showDetailFromSummary" in RESULTS
    assert "_mountInlineDetail" in RESULTS
    assert "异构体组：" in RESULTS
    assert "Top1 完整结构：" in RESULTS
    assert "未命中专属规则" not in RESULTS
    assert "已运行，未命中类别专属推理规则" not in RESULTS
    assert "逐维规则评价已完成，本次未观察到正向支持" in RESULTS
    assert "样本结果摘要" in OVERVIEW
    assert "sample-result-summary-heading" in OVERVIEW
    assert "sample-result-analysis-tags" in OVERVIEW


def test_task_queue_shows_final_main_count_and_unresolved_exclusions():
    assert "compound_count_semantics" in TASKS
    assert "个主要化合物结果" in TASKS
    assert "个未解析事件（未计数）" in TASKS
    assert "个去重化合物" not in TASKS
    assert "blankPairingConfirmed" in TASKS
    assert "方法可比性与背景排除已确认" in TASKS


def test_english_result_overview_localizes_repeatability_and_warning_counts():
    assert "Technical-injection consistency" in OVERVIEW
    assert "High groups per injection" in OVERVIEW
    assert "Repeat-robust High" in OVERVIEW
    assert "payload.interpretation ||" in OVERVIEW
    assert "View counts by injection" in OVERVIEW
    assert "analysis notice${warningCount === 1 ? '' : 's'}" in OVERVIEW
    assert "data-scientific-value" in OVERVIEW


def test_isomer_competition_helper_is_exported_across_overview_scopes():
    assert "groupHasIsomerCompetition," in OVERVIEW
    assert (
        "SampleResultData.groupHasIsomerCompetition(group)"
        in OVERVIEW
    )
    assert "const hasIsomerCompetition = groupHasIsomerCompetition(group)" not in OVERVIEW


def test_v15_candidate_identity_is_atomic_and_loaded_before_compound_projection():
    assert "const candidate = { ...(selected.candidate || {}) };" in REVIEW
    assert "'compound_name', 'display_class', 'classification'" not in REVIEW
    assert "Scientific identity is atomic" in REVIEW
    assert "Unresolved same-formula isomer group" in REVIEW
    flagship_load = OVERVIEW.index("const flagship = await API.get(")
    compound_load = OVERVIEW.index("/compounds?limit=2000")
    assert flagship_load < compound_load
    assert "candidate_order_source" not in OVERVIEW


def test_unmatched_workspace_surfaces_retrieval_only_candidates_without_promoting_them():
    assert "function retrievalOnlyCandidateSupport" in REVIEW
    assert "flagshipEntityForResultIndex" in REVIEW
    assert "Spectral-retrieval candidate hypotheses — no primary result" in REVIEW
    assert "Database names below are displayed exactly as stored" in REVIEW
    assert "不构成唯一二维结构指认" in REVIEW
    assert "化合物计数增加 0" in REVIEW
    assert "data-scientific-value>${esc(databaseName)}" in REVIEW
    assert "item.candidate?.compound_name" in REVIEW
    assert "interfaceCopy(item.candidate?.compound_name" not in REVIEW


def test_integrated_analysis_has_purpose_tabs_and_review_page_redirects():
    result_page = _page(INDEX, "results")
    review_page = _page(INDEX, "review")

    assert 'id="reviewTaskSelect"' not in result_page
    assert 'id="reviewTaskSelect"' in review_page
    assert 'data-review-tab="identification"' in review_page
    assert 'data-review-tab="unmatched"' in review_page
    assert 'data-review-tab="process"' in review_page
    assert 'data-review-tab="followup"' in review_page
    assert 'id="resultsArea"' in review_page
    assert 'id="unmatchedWorkbenchPanel"' in review_page
    assert 'id="moduleEffectPanel"' in review_page
    assert 'id="eviAcquirePanel"' in review_page
    assert "classList.toggle('active'" in REVIEW
    assert "App.navigateTo('review')" in REVIEW
    assert "identification-workspace" in REVIEW
    assert "candidateEvidencePane" in REVIEW
    assert "result_index" in REVIEW
    assert "#review?" in REVIEW
    assert "async function enterPage" in REVIEW
    assert "preferredResultIndex == null" in REVIEW
    assert "rawResultIndex == null" in REVIEW
    assert "EvidenceReview.enterPage()" in (ROOT / "web/static/js/app.js").read_text(encoding="utf-8")
    assert "候选证据工作台" in review_page
    assert "未匹配谱图" in review_page
    assert "分析质量与适用性" in review_page
    assert "确认实验建议" in review_page


def test_quality_and_followup_pages_use_scannable_visual_hierarchy():
    assert "analysis-quality-dashboard" in MODULE_EFFECTS
    assert "quality-completion-ring" in MODULE_EFFECTS
    assert "quality-meter-row" in MODULE_EFFECTS
    assert "quality-module-timeline" in MODULE_EFFECTS
    assert "查看指标与适用边界" in MODULE_EFFECTS
    assert "grid-template-columns:repeat(auto-fit" not in MODULE_EFFECTS
    assert "acquisition-route" in ACQUIRE
    assert "优先行动队列" in ACQUIRE
    assert "acquisition-plan-list" in ACQUIRE
    assert "INITIAL_PLAN_LIMIT = 8" in ACQUIRE
    assert "查看依据、缺失证据与导出" in ACQUIRE
    assert "grid-template-columns:repeat(auto-fit" not in ACQUIRE
    for css_class in (
        ".quality-completion-ring",
        ".quality-meter-row",
        ".quality-module-timeline",
        ".acquisition-route",
        ".acquisition-plan-row",
    ):
        assert css_class in STYLE


def test_formal_report_language_is_plain_and_technical_audit_is_collapsed():
    user_facing_sources = "\n".join((INDEX, OVERVIEW, RESULTS, MODULE_EFFECTS, REPORT))

    assert "正式报告条件（不改变候选鉴定）" in RESULTS
    assert "报告条件未满足" in RESULTS
    assert "暂缓确认" not in "\n".join((INDEX, OVERVIEW, RESULTS, REVIEW))
    assert "待确认" not in "\n".join((INDEX, OVERVIEW, RESULTS, REVIEW))
    assert "正式报告允许范围" in RESULTS
    assert "主要缺失条件（最多两项）" in RESULTS
    assert "技术审计信息（默认收起）" in RESULTS
    assert '<details class="candidate-report-audit">' in RESULTS
    for internal_wording in (
        "EviGate-MS 证据护照",
        "证据治理阻断",
        "最终声称出口",
        "未释放",
        "证据复核与后续实验",
        "EviGate 2 最终声称编译",
        "影子证据与反事实模块",
    ):
        assert internal_wording not in user_facing_sources
    assert "鉴定结论条件评估" in MODULE_EFFECTS
    assert "补充证据与一致性检查" in MODULE_EFFECTS
    assert "heuristically_supported_edges: '启发式支持边'" in MODULE_EFFECTS
    assert "truncated_by_task_limit: '达到影子任务上限'" in MODULE_EFFECTS
    assert "high_promotions: 'High 提升'" in MODULE_EFFECTS


def test_integrated_report_contains_decision_relevant_sections():
    report_page = _page(INDEX, "report")

    assert 'id="sampleReportArea"' in report_page
    assert "样本总体结论" in REPORT
    assert "数据质量与色谱概览" in REPORT
    assert "总离子流色谱图（TIC）" in REPORT
    assert "基峰色谱图（BPC）" in REPORT
    assert "tic_overview" in REPORT
    assert "report-tic-svg" in REPORT
    assert "去重化合物结论构成" in REPORT
    assert "去重候选身份与检测证据" in REPORT
    assert "${esc(displayName(group))}" in REPORT
    assert "hasIsomerCompetition ? leadName : displayName(group)" not in REPORT
    assert "未注释离子特征" in REPORT
    assert "High 级鉴定组" in REPORT
    assert "High 单一二维连接" in REPORT
    assert "High 异构体组" in REPORT
    assert "异构体候选排序明细" in REPORT
    assert "技术重复针一致性" in REPORT
    assert "isomer-high" in REPORT
    assert "isomer-non-high" in REPORT
    assert "group.status.hasIsomerHigh" in REPORT
    assert "中等证据异构体组" in REPORT
    assert "证据不足线索（默认收起）" in REPORT
    assert "独立靶向化合物线索" in REPORT
    assert "分析方法与参数" in REPORT
    assert "证据边界与结果局限" in REPORT
    assert "建议的确认实验" in REPORT
    assert "报告范围与追溯" in REPORT
    assert "不重新评分" in REPORT
    assert "structuredRuleEvidenceLabel" in REPORT
    assert "reasoning[0]" not in REPORT
    assert "未保存可审计的逐维规则证据" in REPORT
    assert "下载样本综合报告（HTML）" in report_page
    assert "批量选择并下载（ZIP）" in report_page
    assert "SampleReport.openBatchExport()" in report_page


def test_integrated_report_supports_hierarchical_expansion_in_preview_and_export():
    assert 'data-report-action="expand"' in REPORT
    assert 'data-report-action="collapse"' in REPORT
    assert 'data-report-action="primary"' in REPORT
    assert 'data-report-kind="compound"' in REPORT
    assert 'data-report-kind="isomer-group"' in REPORT
    assert 'class="report-compound-card"' in REPORT
    assert 'class="report-isomer-group"' in REPORT
    assert "bindReportInteractions(area)" in REPORT
    assert "standaloneInteractionScript()" in REPORT
    assert "beforeprint" in REPORT
    assert "本异构体组的完整竞争候选" in REPORT
    assert "Integrated report v3: hierarchical, expandable reading experience." in STYLE
    assert ".sample-report-details > summary:not(.report-section-summary)" in STYLE
    assert "grid-template-columns: minmax(0,1fr) max-content" in STYLE
    assert "sample-report-details>summary:not(.report-section-summary)" in REPORT
    assert "grid-template-columns:minmax(0,1fr) max-content" in REPORT


def test_integrated_report_supports_multi_sample_zip_from_batch_export():
    export_page = _page(INDEX, "export")
    export_script = (ROOT / "web/static/js/export.js").read_text(encoding="utf-8")

    assert 'id="batchExportCard"' in export_page
    assert 'value="html_integrated|zh"' in export_page
    assert "新版样本综合报告（每样本 1 份，ZIP）" in export_page
    assert "async function generateBatch" in REPORT
    assert "00_批量报告清单.csv" in REPORT
    assert "buildZip(entries, onProgress)" in REPORT
    assert "format === 'html_integrated'" in export_script
    assert "SampleReport.generateBatch(taskIds" in export_script


def test_integrated_report_uses_readable_expandable_compound_cards_and_native_tic():
    assert "鉴定范围与结论边界" in REPORT
    assert "分类、生源、碎裂与规则证据" in REPORT
    assert "支持与代表 MS2 原始定位" in REPORT
    assert "离子身份与排序证据" in REPORT
    assert "Integrated sample report v2: readable, evidence-first layout" in STYLE
    assert "Integrated report v3: hierarchical, expandable reading experience" in STYLE
    assert ".report-compound-card > summary" in STYLE
    assert ".report-tic-svg" in STYLE
    assert "font-size: 0.86rem" in STYLE
    assert "grid-template-columns: repeat(5, minmax(0,1fr))" in STYLE


def test_main_result_table_separates_evidence_grade_from_claim_decision():
    assert "function finalClaimPresentation" not in OVERVIEW
    assert "<th>鉴定结论</th>" in OVERVIEW
    assert "<th>关键证据</th>" in OVERVIEW
    assert "group.status.label" in OVERVIEW
    assert "claim.label" not in OVERVIEW


def test_sample_summary_restores_native_tic_and_bpc_without_optical_spectra():
    assert "/tic_overview" in OVERVIEW
    assert "summaryChromatograms" in OVERVIEW
    assert "function summaryChromatograms" in REPORT
    assert "TIC / BPC 色谱摘要" in REPORT
    assert "DAD/UV 辅助光学谱不进入这两张图" in REPORT
    assert "sample-result-chromatograms" in STYLE
    assert "sample-chromatogram-grid" in STYLE


def test_result_and_detail_text_is_darkened_and_uses_readable_expert_labels():
    assert "font-size:9px" not in RESULTS
    assert "#page-results," in STYLE
    assert "#page-review," in STYLE
    assert ".candidate-inline-detail .score-lbl" in STYLE
    assert ".sample-candidate-detail-action" in STYLE
    assert ".sample-compound-table col.scope-col" in STYLE


def test_candidate_detail_explains_sparse_precursor_high_gate_in_chinese():
    assert "m.precursor_high_evidence_gate_required" in RESULTS
    assert "m.precursor_high_evidence_gate_reason_zh" in RESULTS
    assert "High 证据门控：该候选已降为 Medium" in RESULTS
    assert "碎片配对稀疏且缺少独立前体支持" in RESULTS


def test_integrated_report_cover_is_compact_and_neutral():
    assert ".sample-report-cover{padding:25px 31px;background:#fff" in REPORT
    assert "box-shadow: inset 0 4px #1d756f" in STYLE


def test_main_result_page_exposes_contaminant_count_filter_and_metadata():
    assert "contaminant_unique_compounds" in OVERVIEW
    assert "data-overview-filter=\"contaminant\"" in OVERVIEW
    assert "contaminant_category" in OVERVIEW
    assert "screening_domain" in OVERVIEW
    assert "农药筛查对象" in OVERVIEW


def test_integrated_report_reads_nested_auxiliary_detector_audit():
    assert "representation.representation_audit" in REPORT
    assert "audit?.counts?.auxiliary?.spectra" in REPORT
    assert "Math.max(" in REPORT


def test_chromatograms_mark_all_priority_evidence_tiers_and_uncovered_major_peaks():
    assert "prominentChromatogramPeaks" in REPORT
    assert "tic-high-marker" in REPORT
    assert "tic-isomer-marker" in REPORT
    assert "tic-medium-marker" in REPORT
    assert "tic-unassigned-marker" in REPORT
    assert "异构体组检出" in REPORT
    assert "暂无重点结果邻近的主要峰" in REPORT
    assert "未建立该基峰离子 MS2 结果的主要峰" in REPORT
    assert "不应解释为“没有化合物”" in REPORT
    assert "chromatogramChart(bundle.groups || [], 'tic')" in REPORT
    assert "chromatogramChart(bundle.groups || [], 'bpc')" in REPORT
    assert "valueKey = isBpc ? 'bpc' : 'tic'" in REPORT
    assert "basePeakMz: Number(point.base_peak_mz)" in REPORT
    assert "function withinPpm(" in REPORT
    assert "bpc_precursor_not_selected" in REPORT
    assert "withinPpm(peak.basePeakMz, marker.precursorMz" in REPORT
    assert "同一 RT 有 MS2 但该基峰离子未被选作隔离目标" in REPORT
    assert "ticChart(highGroups)" not in REPORT
    for css_class in (
        ".bpc-main-line",
        ".legend-bpc",
        ".tic-isomer-marker",
        ".tic-medium-marker",
        ".tic-unassigned-marker",
        ".report-tic-coverage",
    ):
        assert css_class in STYLE


def test_internal_action_wording_is_replaced_by_professional_terms():
    assert "查看全部动作" not in ACQUIRE
    assert "候选动作" not in ACQUIRE
    assert "查看完整采集方案" in ACQUIRE
    assert "当前鉴定排序未改变" in ACQUIRE
    assert "预测信息不作为确认性证据" in ACQUIRE
