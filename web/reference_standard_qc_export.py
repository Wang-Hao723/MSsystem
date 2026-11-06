"""Report projections for the explicit reference-standard QC sidecar."""

from __future__ import annotations

import html
import csv
import io
from typing import Any, Mapping


def _ratio_text(numerator: Any, denominator: Any) -> str:
    try:
        numerator_value = int(numerator or 0)
        denominator_value = int(denominator or 0)
    except (TypeError, ValueError):
        return f"{numerator} / {denominator}"
    percentage = (
        100.0 * numerator_value / denominator_value
        if denominator_value
        else 0.0
    )
    return f"{numerator_value}/{denominator_value} ({percentage:.1f}%)"


def reference_standard_qc_summary_rows(
    payload: Mapping[str, Any], *, language: str = "zh"
) -> list[dict[str, Any]]:
    summary = payload.get("summary") or {}
    if not payload.get("reference_standard_set") or not summary:
        return []
    declared = summary.get("declared_truth_identities")
    observable = summary.get(
        "single_injection_compatible_precursor_identities"
    )
    if language == "en":
        rows = [
            ("Declared reference standards", declared, "Declared set denominator"),
            (
                "Observable in this injection",
                _ratio_text(observable, declared),
                "Compatible precursor event present in this polarity/injection",
            ),
            (
                "Observable standards retrieved in saved Top-K",
                _ratio_text(summary.get("observable_saved_topk_connectivity"), observable),
                "Denominator: observable standards",
            ),
            (
                "Observable standards recovered at Top-1",
                _ratio_text(summary.get("observable_top1_connectivity"), observable),
                "Denominator: observable standards",
            ),
            (
                "Observable standards supported at reportable scope",
                _ratio_text(summary.get("observable_reported_scope_supported"), observable),
                "Denominator: observable standards",
            ),
            (
                "Observable standards reaching High reportable scope",
                _ratio_text(summary.get("observable_high_reportable_scope_supported"), observable),
                "Denominator: observable standards; High is not a recovery denominator",
            ),
        ]
        return [
            {
                "Item": label,
                "Value": value,
                "Note": f"{note}. Post-hoc QC only; excluded from candidate generation, scoring and ranking",
            }
            for label, value, note in rows
        ]
    rows = [
        ("混标声明标准品", declared, "声明真值集分母"),
        (
            "本次进样当前极性可观测标准品",
            _ratio_text(observable, declared),
            "本次极性/进样中存在兼容前体事件",
        ),
        (
            "可观测标准品候选层找回",
            _ratio_text(summary.get("observable_saved_topk_connectivity"), observable),
            "分母为本针可观测标准品",
        ),
        (
            "可观测标准品 Top1 找回",
            _ratio_text(summary.get("observable_top1_connectivity"), observable),
            "分母为本针可观测标准品",
        ),
        (
            "可观测标准品可报告范围找回",
            _ratio_text(summary.get("observable_reported_scope_supported"), observable),
            "分母为本针可观测标准品",
        ),
        (
            "可观测标准品达到 High 报告范围",
            _ratio_text(summary.get("observable_high_reportable_scope_supported"), observable),
            "分母为本针可观测标准品；High 不是召回率分母",
        ),
    ]
    return [
        {
            "项目": label,
            "结果": value,
            "说明": f"{note}；仅后验质控，不参与候选生成、评分或排序",
        }
        for label, value, note in rows
    ]


def reference_standard_qc_detail_rows(
    payload: Mapping[str, Any], *, language: str = "zh"
) -> list[dict[str, Any]]:
    rows = []
    for item in payload.get("per_standard") or []:
        best = item.get("best_correct_top1") or {}
        observable = bool(item.get("observable_in_injection"))
        if language == "en":
            rows.append(
                {
                    "Standard": item.get("name_en") or item.get("name_zh"),
                    "Chinese_Name": item.get("name_zh"),
                    "Formula": item.get("formula"),
                    "Observable_In_Injection": "yes" if observable else "no",
                    "Compatible_Precursor_Events": item.get("compatible_precursor_event_count"),
                    "Saved_TopK_Retrieved": (
                        "yes" if item.get("saved_topk_connectivity") else "no"
                    ) if observable else "not applicable",
                    "Best_Saved_Rank": (
                        item.get("best_saved_rank_connectivity")
                        if observable else "not applicable"
                    ),
                    "Top1_Retrieved": (
                        "yes" if item.get("top1_connectivity_any_scan") else "no"
                    ) if observable else "not applicable",
                    "Reportable_Scope_Supported": (
                        "yes" if item.get("reported_scope_supported_any_scan") else "no"
                    ) if observable else "not applicable",
                    "High_Reportable_Scope": (
                        "High" if item.get("high_reportable_scope_supported_any_scan") else "no"
                    ) if observable else "not applicable",
                    "Best_Correct_Candidate": best.get("candidate_name"),
                    "Cosine": best.get("cosine"),
                }
            )
        else:
            rows.append(
                {
                    "标准品": item.get("name_en") or item.get("name_zh"),
                    "中文名": item.get("name_zh"),
                    "分子式": item.get("formula"),
                    "本针可观测": "是" if observable else "否",
                    "兼容前体事件数": item.get("compatible_precursor_event_count"),
                    "候选层找回": (
                        "是" if item.get("saved_topk_connectivity") else "否"
                    ) if observable else "不适用",
                    "最佳保存排名": (
                        item.get("best_saved_rank_connectivity")
                        if observable else "不适用"
                    ),
                    "Top1 找回": (
                        "是" if item.get("top1_connectivity_any_scan") else "否"
                    ) if observable else "不适用",
                    "可报告范围支持": (
                        "是" if item.get("reported_scope_supported_any_scan") else "否"
                    ) if observable else "不适用",
                    "达到 High 报告范围": (
                        "High" if item.get("high_reportable_scope_supported_any_scan") else "未达 High"
                    ) if observable else "不适用",
                    "最佳正确候选": best.get("candidate_name"),
                    "Cosine": best.get("cosine"),
                }
            )
    return rows


def render_reference_standard_qc_html(
    payload: Mapping[str, Any], *, language: str = "zh"
) -> str:
    standard_set = payload.get("reference_standard_set") or {}
    summary = payload.get("summary") or {}
    if not standard_set or not summary:
        return ""
    english = language == "en"
    title = (
        f"{standard_set.get('label') or standard_set.get('id')} — reference-standard recovery audit"
        if english
        else f"{standard_set.get('label') or standard_set.get('id')} · 混标回收审计"
    )
    labels = (
        ["Declared", "Observable", "Top-K retrieved", "Top-1 recovered", "Reportable scope", "High scope"]
        if english
        else ["声明标准品", "本针可观测", "候选层找回", "Top1 找回", "可报告范围找回", "达到 High 范围"]
    )
    declared = summary.get("declared_truth_identities")
    observable = summary.get(
        "single_injection_compatible_precursor_identities"
    )
    values = [
        declared,
        _ratio_text(observable, declared),
        _ratio_text(summary.get("observable_saved_topk_connectivity"), observable),
        _ratio_text(summary.get("observable_top1_connectivity"), observable),
        _ratio_text(summary.get("observable_reported_scope_supported"), observable),
        _ratio_text(summary.get("observable_high_reportable_scope_supported"), observable),
    ]
    header = (
        "<tr><th>Standard</th><th>Observable</th><th>Top-K</th><th>Top-1</th><th>Reportable scope</th><th>High scope</th></tr>"
        if english
        else "<tr><th>标准品</th><th>本针可观测</th><th>候选层</th><th>Top1</th><th>可报告范围</th><th>High 范围</th></tr>"
    )
    detail_rows = []
    yes = "yes" if english else "是"
    no = "no" if english else "否"
    not_applicable = "n/a" if english else "不适用"
    for item in payload.get("per_standard") or []:
        observable_item = bool(item.get("observable_in_injection"))

        def evidence_cell(key: str, *, high: bool = False) -> str:
            if not observable_item:
                return not_applicable
            if item.get(key):
                return "High" if high else yes
            return no

        detail_rows.append(
            "<tr>"
            f"<td><b>{html.escape(str(item.get('name_en') or item.get('name_zh') or '—'))}</b>"
            f"<small style='display:block;color:#64748b'>{html.escape(str(item.get('name_zh') or ''))}</small></td>"
            f"<td>{yes if observable_item else no}</td>"
            f"<td>{evidence_cell('saved_topk_connectivity')}</td>"
            f"<td>{evidence_cell('top1_connectivity_any_scan')}</td>"
            f"<td>{evidence_cell('reported_scope_supported_any_scan')}</td>"
            f"<td>{evidence_cell('high_reportable_scope_supported_any_scan', high=True)}</td>"
            "</tr>"
        )
    cards = "".join(
        f"<div style='padding:10px;border:1px solid #99f6e4;background:white;border-radius:8px'><b style='font-size:20px;color:#0f766e'>{html.escape(str(value))}</b><br><small>{html.escape(label)}</small></div>"
        for label, value in zip(labels, values)
    )
    subtitle = (
        "Post-hoc QC only. Truth was excluded from candidate generation, scoring and ranking."
        if english
        else "仅作盲检索完成后的质量控制；真值未参与候选生成、评分或排序。"
    )
    denominator_warning = str(
        payload.get("denominator_warning")
        if english
        else payload.get("denominator_warning_zh")
        or payload.get("denominator_warning")
        or ""
    )
    return f"""
<section id="reference-standard-qc-report" style="margin:24px 0;padding:18px;border:1px solid #5eead4;border-left:5px solid #0f766e;border-radius:10px;background:#f0fdfa;color:#172033">
  <h2 style="margin:0;color:#115e59;font-size:20px">{html.escape(title)}</h2>
  <p style="margin:4px 0;color:#475569;font-size:13px">{html.escape(subtitle)}</p>
  <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(135px,1fr));gap:9px;margin-top:14px">{cards}</div>
  <p style="margin:12px 0;padding:9px 11px;background:white;border-left:4px solid #0f766e;border-radius:7px;font-size:12px;line-height:1.55">{html.escape(denominator_warning)}</p>
  <div style="overflow-x:auto"><table style="width:100%;border-collapse:collapse;font-size:12px"><thead>{header}</thead><tbody>{''.join(detail_rows)}</tbody></table></div>
</section>
"""


def append_reference_standard_qc_csv_section(
    writer: Any,
    payload: Mapping[str, Any],
    *,
    language: str = "zh",
) -> None:
    """Append a denominator-safe, post-hoc mixed-standard QC section."""

    summary_rows = reference_standard_qc_summary_rows(
        payload, language=language
    )
    detail_rows = reference_standard_qc_detail_rows(
        payload, language=language
    )
    if not summary_rows:
        return
    english = language == "en"
    writer.writerow([])
    writer.writerow([
        "=== Reference-standard recovery audit (post-hoc; excluded from search/ranking) ==="
        if english
        else "=== 混标回收审计（仅后验质控；不参与检索或排序）==="
    ])
    summary_columns = ["Item", "Value", "Note"] if english else ["项目", "结果", "说明"]
    writer.writerow(summary_columns)
    for row in summary_rows:
        writer.writerow([row.get(column, "") for column in summary_columns])
    if not detail_rows:
        return
    writer.writerow([])
    columns = list(detail_rows[0])
    writer.writerow(columns)
    for row in detail_rows:
        writer.writerow([row.get(column, "") for column in columns])


def reference_standard_qc_csv_text(
    payload: Mapping[str, Any], *, language: str = "zh"
) -> str:
    """Return a standalone reference-standard QC CSV projection."""

    buffer = io.StringIO()
    writer = csv.writer(buffer)
    append_reference_standard_qc_csv_section(
        writer, payload, language=language
    )
    return buffer.getvalue().lstrip("\r\n")


def reference_standard_qc_pdf_story(
    payload: Mapping[str, Any],
    *,
    language: str = "zh",
    max_rows: int = 100,
) -> list[Any]:
    """Build a compact PDF section with explicit, non-mixed denominators."""

    standard_set = payload.get("reference_standard_set") or {}
    summary = payload.get("summary") or {}
    if not standard_set or not summary:
        return []

    from reportlab.lib import colors
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.platypus import Paragraph, Spacer, Table, TableStyle
    from web.expected_compounds_export import reportlab_font_names

    english = language == "en"
    regular_font, bold_font = reportlab_font_names("en" if english else "zh")
    styles = getSampleStyleSheet()
    for style_name in ("Normal", "Heading2"):
        styles[style_name].fontName = regular_font
    note = ParagraphStyle(
        f"ReferenceStandardQCNote_{'en' if english else 'zh'}",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=7,
        leading=9,
        textColor=colors.HexColor("#334155"),
    )
    cell = ParagraphStyle(
        f"ReferenceStandardQCCell_{'en' if english else 'zh'}",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=6.2,
        leading=7.4,
        textColor=colors.HexColor("#172033"),
    )
    declared = summary.get("declared_truth_identities")
    observable = summary.get(
        "single_injection_compatible_precursor_identities"
    )
    summary_headers = (
        ["Declared", "Observable", "Top-K", "Top-1", "Reportable", "High scope"]
        if english
        else ["声明数", "本针可观测", "候选层", "Top1", "可报告范围", "High 范围"]
    )
    summary_values = [
        str(declared),
        _ratio_text(observable, declared),
        _ratio_text(summary.get("observable_saved_topk_connectivity"), observable),
        _ratio_text(summary.get("observable_top1_connectivity"), observable),
        _ratio_text(summary.get("observable_reported_scope_supported"), observable),
        _ratio_text(summary.get("observable_high_reportable_scope_supported"), observable),
    ]
    title = (
        f"{standard_set.get('label') or standard_set.get('id')} — reference-standard recovery audit"
        if english
        else f"{standard_set.get('label') or standard_set.get('id')} · 混标回收审计"
    )
    boundary = (
        "This is a post-hoc quality-control projection. The declared identities were excluded from candidate generation, scoring and ranking. Top-K, Top-1, reportable-scope and High-scope rates use only standards with a compatible precursor event in this injection."
        if english
        else "本节仅作盲检索完成后的质量控制；声明真值未参与候选生成、评分或排序。候选层、Top1、可报告范围和 High 范围均以本针存在兼容前体事件的标准品为分母。"
    )
    story: list[Any] = [
        Spacer(1, 0.4 * cm),
        Paragraph(html.escape(title), styles["Heading2"]),
        Paragraph(html.escape(boundary), note),
    ]
    summary_table = Table(
        [summary_headers, summary_values],
        colWidths=[2.65 * cm] * 6,
        repeatRows=1,
    )
    summary_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#0F766E")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), bold_font),
        ("FONTNAME", (0, 1), (-1, -1), regular_font),
        ("FONTSIZE", (0, 0), (-1, -1), 7),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#99F6E4")),
    ]))
    story.extend([summary_table, Spacer(1, 0.15 * cm)])
    warning = str(
        payload.get("denominator_warning")
        if english
        else payload.get("denominator_warning_zh")
        or payload.get("denominator_warning")
        or ""
    ).strip()
    if warning:
        story.extend([
            Paragraph(html.escape(warning), note),
            Spacer(1, 0.12 * cm),
        ])

    details = list(payload.get("per_standard") or [])[:max_rows]
    if details:
        yes = "yes" if english else "是"
        no = "no" if english else "否"
        na = "n/a" if english else "不适用"
        headers = (
            ["Standard", "Formula", "Observable", "Top-K", "Top-1", "Reportable", "High"]
            if english
            else ["标准品", "分子式", "可观测", "候选层", "Top1", "可报告范围", "High"]
        )
        table_rows: list[list[Any]] = [headers]
        for item in details:
            is_observable = bool(item.get("observable_in_injection"))

            def flag(key: str, *, high: bool = False) -> str:
                if not is_observable:
                    return na
                if item.get(key):
                    return "High" if high else yes
                return no

            table_rows.append([
                Paragraph(
                    html.escape(str(item.get("name_en") or item.get("name_zh") or "—")),
                    cell,
                ),
                str(item.get("formula") or "—"),
                yes if is_observable else no,
                flag("saved_topk_connectivity"),
                flag("top1_connectivity_any_scan"),
                flag("reported_scope_supported_any_scan"),
                flag("high_reportable_scope_supported_any_scan", high=True),
            ])
        detail_table = Table(
            table_rows,
            colWidths=[4.2 * cm, 2.0 * cm, 1.75 * cm, 1.75 * cm, 1.65 * cm, 2.15 * cm, 1.65 * cm],
            repeatRows=1,
        )
        detail_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#115E59")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), bold_font),
            ("FONTNAME", (0, 1), (-1, -1), regular_font),
            ("FONTSIZE", (0, 0), (-1, -1), 6.2),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("ALIGN", (1, 1), (-1, -1), "CENTER"),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F0FDFA")]),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#99F6E4")),
        ]))
        story.append(detail_table)
    return story


def inject_reference_standard_qc_html(
    html_content: str, payload: Mapping[str, Any], *, language: str = "zh"
) -> str:
    section = render_reference_standard_qc_html(payload, language=language)
    if not section:
        return html_content
    lower = html_content.lower()
    for marker in (
        '<section id="expected-compounds-report"',
        '<div id="ms2-identification-report"></div>',
        '<section class="table-box"',
        '<div class="table-container"',
        "</body>",
    ):
        position = lower.find(marker) if marker != "</body>" else lower.rfind(marker)
        if position >= 0:
            return html_content[:position] + section + html_content[position:]
    return html_content + section


__all__ = [
    "append_reference_standard_qc_csv_section",
    "inject_reference_standard_qc_html",
    "reference_standard_qc_detail_rows",
    "reference_standard_qc_csv_text",
    "reference_standard_qc_pdf_story",
    "reference_standard_qc_summary_rows",
    "render_reference_standard_qc_html",
]
