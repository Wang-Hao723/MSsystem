"""Full evidence reports from saved native conclusions and independent MS1."""

from __future__ import annotations

from collections import Counter
from html import escape
import json


def _text(value):
    if value is None or value == "":
        return "—"
    if isinstance(value, (dict, list)):
        value = json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return escape(str(value))


def _table(columns, rows, *, table_id):
    headers = "".join("<th>" + escape(label) + "</th>" for label in columns)
    body = "".join("<tr>" + "".join("<td>" + _text(value) + "</td>" for value in row) + "</tr>" for row in rows)
    return f'<div class="scroll"><table id="{table_id}"><thead><tr>{headers}</tr></thead><tbody>{body}</tbody></table></div>'


def render_frozen_sample_report(bridge, task_id):
    task = bridge.tasks[task_id]
    rows = []
    page = None
    while True:
        page = bridge.read(task_id, "compound_page", report_scope="all", result_view="native",
                           conclusion=None, offset=len(rows), limit=500)
        current = page["compounds"]
        if not current and page["has_more"]:
            raise RuntimeError("全量报告的原生结论分页中断")
        rows.extend(current)
        if not page["has_more"]:
            break
    if len(rows) != page["native_total"]:
        raise RuntimeError("全量报告原生结论行数与冻结产物不一致")
    indices = [row["frozen_compound_index"] for row in rows]
    if len(indices) != len(set(indices)):
        raise RuntimeError("全量报告存在重复原生结论索引")
    expected = bridge.read(task_id, "expected")
    features = expected["features"]
    summary = task["completed_batch"]["summary"]
    if len(features) != summary["ms1_feature_records"]:
        raise RuntimeError("全量报告MS1特征行数与冻结产物不一致")
    statuses = Counter(row.get("conclusion", {}).get("key", "unknown") for row in rows)
    record_types = Counter(row.get("record_type", "unknown") for row in rows)
    prefix = f"/api/tasks/{task_id}"
    # Join only the spectrum IDs already saved in MS1 associations. No new
    # mass/RT matching, scoring, or identity decision is performed here.
    ms1_by_spectrum = {}
    for feature in features:
        for association in feature.get("associated_ms2") or []:
            if isinstance(association, dict) and association.get("spectrum_id") is not None:
                ms1_by_spectrum.setdefault(str(association["spectrum_id"]), []).append(feature)
    native_rows = []
    for row in rows:
        representative = row.get("representative") or {}
        lead = row.get("lead_candidate") or {}
        match = representative.get("match") or {}
        role = row.get("plant_report") or {}
        spectrum_ids = [str(representative.get("spectrum_id", "")),
                        *[str(value) for value in (row.get("lineage") or {}).get("previous_scan_record_ids", [])]]
        linked_features = {}
        for spectrum_id in spectrum_ids:
            for feature in ms1_by_spectrum.get(spectrum_id, []):
                linked_features[feature["feature_id"]] = feature
        if role.get("excluded_from_plant_candidate_table"):
            source_scope = "外源／身份冲突QC"
        elif not role.get("plant_origin_verified"):
            source_scope = "来源未核实；谱库假设QC"
        else:
            source_scope = "来源字段已有核验记录"
        native_rows.append([
            row["frozen_compound_index"], row.get("record_type"), row.get("feature_id"),
            row.get("conclusion", {}).get("key"),
            row.get("conclusion", {}).get("label"),
            lead.get("compound_name") or row.get("display_name"), row.get("formula"),
            row.get("identity_scope"), source_scope, representative.get("spectrum_id"),
            representative.get("result_index"), representative.get("precursor_mz"),
            representative.get("retention_time"), representative.get("retention_time_unit"),
            match.get("cosine"), lead.get("rank_score", lead.get("final_score", match.get("rank_score", match.get("final_score")))),
            representative.get("blank_flagged"), representative.get("blank_clean_view_excluded"),
            representative.get("blank_sample_ratio"),
            list(linked_features),
            [{"feature_id": feature["feature_id"], "evidence_level": feature.get("evidence_level"),
              "blank_category": feature.get("blank_evidence_category"),
              "blank_sample_ratio": feature.get("blank_sample_ratio"),
              "blank_clean_view_excluded": feature.get("blank_clean_view_excluded")}
             for feature in linked_features.values()],
        ])
    ms1_rows = [[
        row.get("feature_id"), row.get("candidate_names", row.get("candidate_name")),
        row.get("formula"), row.get("adduct"), row.get("evidence_level"),
        row.get("apex_rt_minutes"), row.get("apex_intensity"),
        row.get("baseline_corrected_area_intensity_seconds"), row.get("mass_error_ppm"),
        row.get("blank_evidence_category"), row.get("blank_sample_ratio"),
        len(row.get("associated_ms2") or []), row.get("exact_identity_confirmed"),
        row.get("associated_ms2"), row.get("blank_flagged"), row.get("blank_clean_view_excluded"),
    ] for row in features]
    native_table = _table([
        "冻结记录索引", "原生记录类型", "特征编号", "原生判定状态", "原生证据标签",
        "原始谱库假设标题（QC）", "分子式", "身份范围", "来源复核分类", "代表谱图编号",
        "谱图原始索引", "前体m/z", "保留时间", "时间单位", "余弦", "原排序分",
        "空白标记", "净化视图排除", "样本／空白比",
        "已保存MS1关联特征", "关联MS1证据与配对空白",
    ], native_rows, table_id="all-native-conclusions")
    ms1_table = _table([
        "MS1特征编号", "靶标候选名称／组", "分子式", "离子形式", "MS1证据等级",
        "峰顶RT（分钟）", "峰顶强度", "基线校正面积（强度×秒）", "质量误差（ppm）",
        "空白分类", "样本／空白比", "关联MS/MS数", "原确证字段",
        "已保存MS/MS关联", "MS1空白标记", "MS1净化视图排除",
    ], ms1_rows, table_id="all-ms1-features")
    title = escape(task["filename"])
    source_sha = escape(bridge.entries[task_id]["source_sha256"])
    native_sha = escape(bridge.entries[task_id]["results"]["sha256"])
    conclusion_sha = escape(bridge.entries[task_id]["conclusions"]["sha256"])
    return f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{title} · 全量MS/MS与MS1证据报告</title>
<style>body{{font:15px/1.65 system-ui,sans-serif;color:#24323c;margin:32px;max-width:1500px}}h1{{font-size:27px}}h2{{margin-top:32px}}.card{{background:#eff7f6;border-left:4px solid #246c71;padding:18px}}.scroll{{overflow:auto;max-height:75vh;border:1px solid #d9e0e3}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{border:1px solid #d9e0e3;padding:7px;vertical-align:top;min-width:85px}}th{{background:#eef3f5;position:sticky;top:0}}code{{overflow-wrap:anywhere}}a{{color:#126368}}.links{{display:flex;gap:18px;flex-wrap:wrap}}@media print{{.scroll{{max-height:none;overflow:visible}}th{{position:static}}body{{margin:10mm}}}}</style></head><body>
<h1>{title} · 全量MS/MS与MS1证据报告</h1>
<div class="card"><b>保存结果全量读取</b><p>已检索 {summary['parsed_searched_ms2_spectra']} 张MS/MS谱，其中 {summary['matched_ms2_spectra']} 张产生谱库候选；原生结论 {len(rows)} 条，独立MS1靶向特征 {len(features)} 个。标准品确证 0。</p>
<p>原生记录类型：{_text(dict(record_types))}。原生判定状态：{_text(dict(statuses))}。</p>
<p>MS1靶向来自系统 XBT_targeted_20260817。XBT与MS1关联的附加筛选有 {page['xbt_linked_total']} 条记录，完整原生结果按上列 {len(rows)} 条读取；候选、异构体、未解析事件、空白和未匹配记录分别保留。下表用已保存的谱图编号连接MS1特征和空白字段。</p>
<p>配对空白标记 {summary['blank_flagged_ms2']} 个MS/MS事件，净化视图排除 {summary['blank_clean_view_excluded_ms2']} 个。下表保留全部事件对应的原生记录及原始空白字段。</p></div>
<p class="links"><a href="{prefix}/frozen-tables/native_compounds">原生结论全表CSV</a><a href="{prefix}/frozen-tables/msms_table">MS/MS全部事件CSV</a><a href="{prefix}/frozen-tables/ms1_table">MS1全部特征CSV</a><a href="{prefix}/download/native-conclusions.json.xz">完整原生结论与证据JSON.XZ</a></p>
<h2>完整原生结论及QC记录 · {len(rows)}条</h2><p>谱库标题保留为候选假设或QC追溯字段。分子来源、原生谱库等级与标准品确证分别记录；来源未核实的记录不列作已核验植物内源成分。</p>{native_table}
<h2>独立MS1靶向特征 · {len(features)}条</h2><p>包含未关联MS/MS的MS1特征；质量匹配、同位素支持及诊断碎片支持按保存字段分别列出。</p>{ms1_table}
<h2>冻结文件绑定</h2><p>源文件SHA-256：<code>{source_sha}</code><br>原生MS/MS压缩结果SHA-256：<code>{native_sha}</code><br>联合结论压缩结果SHA-256：<code>{conclusion_sha}</code></p>
<p>本报告读取既有结果，不执行谱图检索、靶向峰提取、评分或结论重建。跨文件数量按记录出现次数计算。</p></body></html>'''
