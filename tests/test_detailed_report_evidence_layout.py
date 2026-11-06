from __future__ import annotations

from web.app import _gen_rich
from web.report_common import (
    build_diagnostic_ion_audit,
    render_diagnostic_ion_audit,
    render_product_structure_candidates,
    row_to_dict,
)
from web.server import _serialize_results


def _detailed_match() -> dict:
    return {
        "compound_name": "9,10-DiHOME; LC-ESI-QIT; MS2; CE",
        "source_compound_name": "9,10-DiHOME; LC-ESI-QIT; MS2; CE",
        "formula": "C18H34O4",
        "smiles": "",
        "inchikey": "",
        "classification": "Oxylipins",
        "dataset_source": "user_contrib:bzl",
        "source_tier": "user_contrib",
        "lib_spectrum_id": "BZL-001",
        "lib_source_accession": "BZL-A-001",
        "adduct": "[M-H]-",
        "cosine": 0.91,
        "final_score": 0.84,
        "mass_ppm": -2.2,
        "mass_score": 0.92,
        "fragment_score": 0.81,
        "intensity_score": 0.73,
        "expert_score": 0.60,
        "confidence": "High",
        "confidence_class": "High",
        "isomer_group_confidence": "Unknown",
        "lib_peaks": [[100.0, 20.0], [171.1386, 100.0], [250.0, 8.0]],
        "pathway_json": {
            "n_spectra": 4,
            "diagnostic_ions": [
                {"mz": 171.1386, "annotation": "stored hit", "freq": 0.8},
                {"mz": 199.1000, "annotation": "stored miss", "freq": 0.6},
            ],
            "neutral_losses": [
                {
                    "loss_da": 18.0106,
                    "fragment_mz": 281.9894,
                    "annotation": "water loss",
                }
            ],
        },
        "expert_details": {
            "diagnostic_ion_score": 0.5,
            "dimension_assessments": {
                "diagnostic_ion": {
                    "status": "scored",
                    "score": 0.5,
                    "reason": "matched_scope_guarded_diagnostic_ion",
                    "evidence_count": 1,
                }
            },
            "diagnostic_ion_analysis": {
                "matched_ions": [
                    {
                        "observed_mz": 171.1044,
                        "expected_mz": 171.1386,
                        "ion_name": "苯甲酸根离子",
                        "ion_formula": "C7H7O5-",
                        "structural_origin": "芳香酸结构",
                        "rule_id": "diagnostic-rule-test",
                    }
                ]
            },
            "neutral_loss_analysis": {
                "best_pathway": [
                    {
                        "from_mz": 300.0 - index,
                        "to_mz": 299.0 - index,
                        "mass_difference": 1.0,
                        "neutral_loss": f"step-{index}",
                        "formula": "H",
                    }
                    for index in range(1, 11)
                ]
            },
        },
    }


def _detailed_row() -> dict:
    return {
        "precursor_mz": 300.0,
        "spectrum_id": "controllerType=0 controllerNumber=1 scan=27",
        "precursor_scan_ref": "scan=26",
        "peaks_count": 4,
        "retention_time": 12.345,
        "retention_time_unit": "minute",
        "collision_energy": 30.0,
        "collision_energy_unit": "eV",
        "ion_mode": "negative",
        "precursor_charge": -1,
        "tic": 123456.0,
        "base_peak": 1000.0,
        # The precursor peak must remain visible in the saved raw-peak table.
        "raw_peaks": [
            [100.0, 20.0],
            [171.1044, 100.0],
            [300.0, 50.0],
            [900.0, 1.0],
        ],
        "query_peaks": [[100.0, 20.0], [171.1044, 100.0]],
        "matches": [_detailed_match()],
    }


def test_detailed_report_separates_scan_counts_and_uses_readable_evidence_layout():
    html = _gen_rich(
        [_detailed_row()],
        source_name="20251112142044-Z098011-6-27.mzML",
        total_spectra=2494,
        matched_spectra=356,
        high_supported_compounds=7,
        high_confidence_scope_groups=8,
        high_confidence_exact_structures=0,
        high_confidence_isomer_groups=8,
    )

    assert '<div class="val">2494</div><div class="lbl">原始 MS2 谱图</div>' in html
    assert '<div class="val">356</div><div class="lbl">有谱库候选的 MS2</div>' in html
    assert '<div class="val">8</div>\n    <div class="lbl">High 级鉴定组</div>' in html
    assert '单一二维连接：0 · 异构体组：8' in html
    assert "9,10-DiHOME</div>" in html
    assert "9,10-DiHOME; LC-ESI-QIT; MS2; CE" in html
    assert "<th>Compound</th>" not in html
    assert '<td colspan="6"' in html
    assert 'data-spectrum-layout="mirror-shared-axis"' in html
    assert 'viewBox="0 0 1120 300"' in html
    assert 'max-width:1080px;margin:0 auto' in html
    assert 'fill="#007A78"' in html
    assert 'fill="#6F4E9B"' in html
    assert 'fill="#C43D55"' in html
    assert "样品 MS/MS（Query，向上）· 质控后用于匹配 2 个峰；原始工件保存 4 个峰" in html
    assert "参考库 MS/MS（Library，向下）· 图中与完整峰表均为 3 个峰" in html
    assert "（绘图/保存）" not in html
    assert ".container{width:calc(100% - 56px);max-width:1520px;margin:28px auto 48px" in html
    assert "@media(max-width:800px){.container{width:calc(100% - 24px);margin:12px auto 24px" in html
    assert 'class="rc-main-content-grid"' in html
    assert 'class="rc-source-box"' in html
    assert 'class="rc-expert-panel"' in html
    assert '规则证据评估' in html
    assert '.rc-frag-viz{width:100%;text-align:center' in html
    assert '.rc-cand-row{display:flex;gap:10px;flex-wrap:wrap;justify-content:center' in html
    assert 'data-query-saved-peaks="4" data-query-plotted-peaks="2"' in html
    assert "未显示的样品原始峰 2 个" in html
    assert "MS2 谱图编号 / native ID" in html
    assert "controllerType=0 controllerNumber=1 scan=27" in html
    assert "离子模式" in html and "negative" in html
    assert "前体电荷" in html and "-1" in html
    assert "300.00000" in html
    assert "900.00000" in html
    assert "旧版归属的当前策略复核" in html
    assert "step-8" in html
    assert "candidate_structure_or_class_not_applicable" in html
    assert "全部有效数值峰已保存" in html
    assert "报告层未做 Top-N 截断" in html
    assert "口径" not in html
    assert "命中" not in html


def test_reportable_scope_high_does_not_reappear_as_zero_or_medium_in_exports():
    row = _detailed_row()
    match = row["matches"][0]
    match.update({
        "confidence": "Medium",
        "compound_confidence": "Medium",
        "exact_structure_confidence": "Medium",
        "two_dimensional_connectivity_confidence": "High",
        "isomer_group_confidence": "High",
        "identification_scope": "exact_2d_connectivity",
        "honesty_tier": "connectivity_high_stereo_unresolved",
    })

    html = _gen_rich([row], source_name="scope-high.mzML")
    report_row = row_to_dict({**row, "match": match}, report_mode=True)
    standard_row = row_to_dict({**row, "match": match}, report_mode=False)

    assert 'data-conf="High" data-full-conf="Medium"' in html
    assert "High（二维连接）" in html
    assert "完整身份 Medium" in html
    assert report_row["Confidence"] == "High"
    assert report_row["Reportable_Identity_Scope"] == "exact_2d_connectivity"
    assert report_row["Exact_Structure_Confidence"] == "Medium"
    assert report_row["Two_Dimensional_Connectivity_Confidence"] == "High"
    assert standard_row["Confidence"] == "High"
    assert standard_row["Reportable_Identity_Label"] == "High（二维连接）"

    match.update({
        "identification_scope": "isomer_group",
        "two_dimensional_connectivity_confidence": "Medium",
    })
    isomer_html = _gen_rich([row], source_name="scope-isomer-high.mzML")
    assert 'data-conf="High" data-full-conf="Medium"' in isomer_html
    assert "High（异构体组）" in isomer_html


def test_diagnostic_ions_show_exact_hit_rows_and_do_not_conflate_na_with_miss():
    match = _detailed_match()
    peaks = [[100.0, 20.0], [171.1044, 100.0]]
    audit = build_diagnostic_ion_audit(match, peaks, tolerance_da=0.05)
    html = render_diagnostic_ion_audit(match, peaks, tolerance_da=0.05)

    assert audit["runtime_status"] == "scored"
    assert audit["pathway_observed_count"] == 1
    assert "符合候选适用条件且已观察到的规则诊断离子（1 条）" in html
    assert "苯甲酸根离子" in html
    assert "171.10440" in html
    assert "171.13860" in html
    assert "0.03420" in html
    assert "参考库记录中该候选的常见碎片离子：样品谱观察到 1/2" in html

    match["expert_details"]["dimension_assessments"]["diagnostic_ion"] = {
        "status": "not_applicable",
        "reason": "classification_unconfirmed",
    }
    match["expert_details"]["diagnostic_ion_analysis"] = {"matched_ions": []}
    not_applicable_html = render_diagnostic_ion_audit(match, peaks)
    assert "未评价（该候选无适用的诊断规则）" in not_applicable_html
    assert "尚无可靠类别归属" in not_applicable_html
    assert "口径" not in not_applicable_html
    assert "命中" not in not_applicable_html

    match["expert_details"]["dimension_assessments"]["diagnostic_ion"] = {
        "status": "not_applicable",
        "reason": "no_curated_rule_for_ontology_class",
        "reason_label_cn": "三级胺氧化物已有可靠类别归属；该维度暂无发布规则",
    }
    ontology_gap_html = render_diagnostic_ion_audit(match, peaks)
    assert "三级胺氧化物已有可靠类别归属；该维度暂无发布规则" in ontology_gap_html
    assert "分类未确认" not in ontology_gap_html


def test_product_structure_cards_expose_auditable_smiles_and_assignment_type():
    html = render_product_structure_candidates([{
        "observed_mz": 213.242535,
        "expected_mz": 213.245101,
        "delta": 0.002566,
        "fragment_smiles": "CCCCCCCCCCCC[N+](C)C",
        "adduct": "[M]+",
        "assignment_type": "candidate_structure_hypothesis",
    }])

    assert "观测峰支持的候选产物结构" in html
    assert 'title="候选产物结构 SMILES"' in html
    assert "CCCCCCCCCCCC[N+](C)C" in html
    assert "candidate_structure_hypothesis" in html


def test_detailed_report_reads_scores_from_the_structured_dimension_ledger():
    row = _detailed_row()
    row["matches"][0]["expert_details"]["dimension_assessments"] = {
        "rule_coverage": {
            "status": "scored",
            "score": 0.67,
            "evidence_count": 2,
            "reason": "assessment_only_score",
        },
        "neutral_loss": {
            "status": "no_evidence",
            "score": 0.0,
            "evidence_count": 0,
            "reason": "evaluated_not_observed",
        },
        "bond_energy": {
            "status": "unavailable",
            "score": 0.0,
            "evidence_count": 0,
            "reason": "model_not_released",
        },
    }
    html = _gen_rich([row], source_name="dimension-ledger.mzML")

    assert 'title="assessment_only_score">0.67 · 2项</span>' in html
    assert "有支持证据 1/3；已评价但未观察 1；不适用/未建立 1" in html
    assert "Rule-evidence dimensions with support" in html


def test_detailed_report_explains_sparse_precursor_high_gate_in_chinese():
    row = _detailed_row()
    match = row["matches"][0]
    match.update({
        "confidence": "Medium",
        "exact_structure_confidence": "Medium",
        "precursor_high_evidence_gate_required": True,
        "precursor_high_evidence_gate_result": (
            "vetoed_to_medium_sparse_ms2_without_precursor_support"
        ),
        "precursor_high_evidence_gate_reason": (
            "only 3 directly matched product ions; precursor not observed"
        ),
        "precursor_high_evidence_gate_reason_zh": (
            "仅有 3 个直接匹配的产物离子；相邻 MS1 未观测到前体"
        ),
    })

    html = _gen_rich([row], source_name="sparse-gate.mzML")

    assert html.count("High 证据门控：") == 2
    assert "该候选已降为 Medium" in html
    assert "仅有 3 个直接匹配的产物离子；相邻 MS1 未观测到前体" in html
    assert "only 3 directly matched product ions" not in html


def test_new_task_artifacts_save_all_raw_and_search_view_peaks_without_top100():
    raw_peaks = [[50.0 + index, 1000.0 - index] for index in range(130)]
    query_peaks = [[60.0 + index, 500.0 - index] for index in range(110)]
    row = {
        "spectrum_id": "scan=full-peaks",
        "precursor_mz": 500.0,
        "peaks_count": len(raw_peaks),
        "raw_peaks": raw_peaks,
        "query_peaks": query_peaks,
        "match": None,
        # Avoid invoking the unmatched classifier in this serialization test.
        "blank_excluded": True,
    }

    serialized = _serialize_results([row])[0]

    assert len(serialized["raw_peaks"]) == 130
    assert len(serialized["query_peaks"]) == 110
    assert serialized["raw_peaks_saved_count"] == 130
    assert serialized["raw_peaks_complete"] is True
    assert serialized["raw_peaks_storage_policy"] == "all_valid_numeric_pairs"


def test_historical_peak_limit_is_disclosed_without_hiding_saved_peaks():
    row = _detailed_row()
    row["peaks_count"] = 125
    row["raw_peaks"] = [
        [50.0 + index, 100.0 - index / 2.0] for index in range(100)
    ]

    html = _gen_rich([row], source_name="historical.mzML")

    assert "历史工件存在保存上限" in html
    assert "源记录为 125 峰，当前历史任务工件仅保存 100 峰" in html
    assert "V11 新任务现按全部有效数值峰保存" in html
    assert "149.00000" in html
