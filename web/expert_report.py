"""
Structured Rule-Evidence HTML Report Generator
==============================================

在保持原有 Standard/Detailed 报告不变的基础上，新增逐维规则证据报告，
展示规则适用性、观测支持、碎片特征库与结构解析推理。

不替代任何现有报告，作为独立的规则证据审计补充报告。
"""

import base64
import io
import json
from datetime import datetime
from pathlib import Path

from web.report_common import classification_evidence_status


def _classification_label(match):
    evidence = classification_evidence_status(match)
    if evidence["status"] == "tentative":
        return f"待核查类别线索：{evidence['label']}"
    if evidence["status"] == "unresolved":
        return "类别尚未可靠确定"
    return evidence["label"]


def _mol_to_base64(smiles, size=(260, 190)):
    """SMILES -> base64 PNG"""
    if not smiles or str(smiles).upper() in ("N/A", "NONE", "NULL", "", "NONE"):
        return ""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        from rdkit.Chem import Draw, AllChem
        mol = Chem.MolFromSmiles(str(smiles).strip())
        if mol is None:
            # Try canonicalization fallback
            try:
                from web.structure import _sanitize_smiles
                sanitized = _sanitize_smiles(str(smiles).strip())
                if sanitized:
                    mol = Chem.MolFromSmiles(sanitized)
            except Exception:
                pass
        if mol is None:
            return ""
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=size)
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return base64.b64encode(buf.getvalue()).decode("ascii")
    except Exception:
        return ""


def _fmt_mz(v):
    """安全格式化 m/z, 容忍 None/str"""
    try:
        return "{:.2f}".format(float(v))
    except (TypeError, ValueError):
        return "—"


def _render_sugar_cascade(sugar):
    """渲染糖苷糖基顺序丢失级联"""
    if not sugar or not sugar.get("cascade"):
        return ""
    chain_parts = []
    for s in sugar["cascade"]:
        chain_parts.append(
            f'<span class="node">{_fmt_mz(s["from_mz"])}</span>'
            f'<span class="arrow">−{_fmt_mz(s["mass_difference"])} ({s.get("sugar_name_cn","")})→</span>'
        )
    chain_parts.append(f'<span class="node aglycone">{_fmt_mz(sugar["aglycone_mz"])}<span class="node-tag">苷元</span></span>')
    comp = sugar.get("sugar_composition", {})
    comp_str = " + ".join(f"{cnt} {name}" for name, cnt in comp.items())
    return (
        '<div class="cascade-box cascade-sugar">'
        '<b>🍬 糖基顺序丢失级联</b> <span class="cascade-sub">糖苷/皂苷 · Domon-Costello</span>'
        f'<div class="cascade-chain">{"".join(chain_parts)}</div>'
        f'<div class="cascade-meta">糖组成: <b>{comp_str}</b> · '
        f'共 {sugar["n_sugars"]} 个糖 · 总糖质量 {sugar["total_sugar_mass"]:.1f} Da</div>'
        '</div>'
    )


def _render_ester_cascade(ester):
    """渲染乌头碱型酯基顺序丢失级联"""
    if not ester or not ester.get("cascade"):
        return ""
    chain_parts = []
    for s in ester["cascade"]:
        chain_parts.append(
            f'<span class="node">{_fmt_mz(s["from_mz"])}</span>'
            f'<span class="arrow">−{_fmt_mz(s["mass_difference"])} ({s.get("loss_name_cn","")})→</span>'
        )
    chain_parts.append(f'<span class="node aglycone">{_fmt_mz(ester["aglycone_mz"])}<span class="node-tag">乌头原碱</span></span>')
    comp = ester.get("loss_composition", {})
    comp_str = " + ".join(f"{cnt}×{name}" for name, cnt in comp.items())
    subtype = ester.get("aconitine_subtype", "")
    return (
        '<div class="cascade-box cascade-ester">'
        '<b>⚗️ 乌头碱型酯基丢失级联</b> '
        f'<span class="subtype-badge">{subtype}</span>'
        f'<div class="cascade-chain">{"".join(chain_parts)}</div>'
        f'<div class="cascade-meta">丢失组成: <b>{comp_str}</b> · '
        f'共 {ester["n_losses"]} 步 · 累计 {ester["total_loss_mass"]:.1f} Da</div>'
        '</div>'
    )


def _render_pathway_table(best_pathway):
    """渲染 fragment→fragment 碎裂步骤表"""
    if not best_pathway:
        return ""
    rows = []
    for s in best_pathway[:6]:
        rows.append(
            f'<tr><td>{_fmt_mz(s.get("from_mz"))}</td><td>{_fmt_mz(s.get("to_mz"))}</td>'
            f'<td>{_fmt_mz(s.get("mass_difference"))}</td>'
            f'<td>{s.get("neutral_loss", "—")}</td>'
            f'<td>{s.get("mechanism", "—")}</td>'
            f'<td>{s.get("priority", "—")}</td></tr>'
        )
    return (
        '<div class="pathway-section"><b>🧬 碎裂步骤</b> <span class="cascade-sub">fragment→fragment 连通链</span>'
        '<table class="fragment-table"><thead><tr>'
        '<th>From m/z</th><th>To m/z</th><th>Δm</th>'
        '<th>中性丢失</th><th>机理</th><th>优先级</th>'
        '</tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>'
    )


def _render_diagnostic_ions(matched_ions):
    """渲染诊断离子表"""
    if not matched_ions:
        return ""
    rows = []
    for ion in matched_ions[:6]:
        spec = ion.get("intensity_characteristic", "")
        spec_cls = "spec-high" if spec == "high" else ("spec-med" if spec == "medium" else "")
        rows.append(
            f'<tr><td>{_fmt_mz(ion.get("observed_mz"))}</td>'
            f'<td>{ion.get("ion_name", "—")}</td>'
            f'<td>{ion.get("structural_origin", "—")}</td>'
            f'<td><span class="spec {spec_cls}">{spec or "—"}</span></td></tr>'
        )
    return (
        '<div class="pathway-section"><b>🔬 诊断离子</b>'
        '<table class="fragment-table"><thead><tr>'
        '<th>m/z</th><th>名称</th><th>结构来源</th><th>特异性</th>'
        '</tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>'
    )


def _render_neutral_losses(nl_list):
    """渲染 precursor→fragment 中性丢失汇总表"""
    if not nl_list:
        return ""
    rows = []
    shown = 0
    for item in nl_list:
        if shown >= 8:
            break
        if item.get("neutral_loss"):
            rows.append(
                f'<tr><td>{_fmt_mz(item.get("mass_difference"))}</td>'
                f'<td class="nl-explained">{item.get("neutral_loss", "")}</td>'
                f'<td>{item.get("formula", "—")}</td>'
                f'<td>{item.get("mechanism", "—")}</td></tr>'
            )
            shown += 1
        elif item.get("possible_candidates"):
            cand = item["possible_candidates"]
            names = ", ".join(c.get("name", "") for c in cand[:2]) if cand else "无候选"
            rows.append(
                f'<tr><td>{_fmt_mz(item.get("mass_difference"))}</td>'
                f'<td colspan="3" class="nl-unexplained">未匹配 (候选: {names})</td></tr>'
            )
            shown += 1
    if not rows:
        return ""
    return (
        '<div class="pathway-section"><b>📉 中性丢失</b> <span class="cascade-sub">precursor→fragment</span>'
        '<table class="fragment-table"><thead><tr>'
        '<th>Δm</th><th>中性丢失</th><th>分子式</th><th>机理</th>'
        '</tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>'
    )


def _render_fragmentation_pathway(ed):
    """渲染完整结构化碎裂途径卡片

    汇总 expert_details 中的级联/路径/诊断离子/中性丢失, 此前这些结构化数据
    仅以压缩 chemical_reasoning 文本呈现, 现以表格+级联图形式完整展示。
    """
    if not ed:
        return ""
    nl = ed.get("neutral_loss_analysis") or {}
    sections = []
    sections.append(_render_sugar_cascade(nl.get("sugar_loss_cascade")))
    sections.append(_render_ester_cascade(nl.get("alkaloid_ester_cascade")))
    sections.append(_render_pathway_table(nl.get("best_pathway")))
    di = ed.get("diagnostic_ion_analysis") or {}
    sections.append(_render_diagnostic_ions(di.get("matched_ions")))
    sections.append(_render_neutral_losses(nl.get("neutral_losses")))
    sections = [s for s in sections if s]
    if not sections:
        return ""
    return '<div class="pathway-card"><div class="pathway-title">🧬 经验碎裂证据与中性丢失序列 (Empirical Fragment Evidence)</div>' + "".join(sections) + '</div>'


def generate_expert_report(results, output_path, mgf_name):
    """
    Generate a structured rule-evidence HTML report.

    在现有 Standard/Detailed 报告之外，提供结构化规则证据分析。
    包含：
    - 8 维规则证据的适用性、覆盖度与观测状态（纯 CSS）
    - 碎片特征库匹配详情
    - 裂解规则覆盖度
    - 诊断离子检测
    - 化学推理摘要
    """
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        pass

    matched = [r for r in results if r.get("matches")]
    matched.sort(key=lambda r: -r["matches"][0].get("final_score", 0))
    total = len(results)
    matched_count = len(matched)
    rate = matched_count / total * 100 if total else 0

    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Structured rule-evidence coverage.  The historical aggregate is retained
    # only inside each card as an explicitly uncalibrated audit field.
    structured_assessments = 0
    assessed_dimension_records = 0
    evidence_bearing_dimension_records = 0
    unavailable_dimension_records = 0
    for r in matched:
        m = r["matches"][0]
        details = m.get("expert_details") or {}
        assessments = details.get("dimension_assessments") or {}
        if not isinstance(assessments, dict) or not assessments:
            continue
        structured_assessments += 1
        for assessment in assessments.values():
            assessment = assessment if isinstance(assessment, dict) else {}
            status = str(assessment.get("status") or "unavailable")
            if status == "unavailable":
                unavailable_dimension_records += 1
                continue
            assessed_dimension_records += 1
            score = float(assessment.get("score") or 0)
            evidence_count = int(assessment.get("evidence_count") or 0)
            if status == "scored" and (evidence_count > 0 or score > 0):
                evidence_bearing_dimension_records += 1

    # Classification stats
    class_counts = {}
    for r in matched:
        cls = _classification_label(r["matches"][0])
        class_counts[cls] = class_counts.get(cls, 0) + 1

    css = """
:root{--teal:#2C7A7B;--teal-dark:#1A5253;--text:#1A202C;--text-secondary:#718096;--border:#E2E8F0;--bg:#F7FAFC;--card-bg:#FFFFFF;--radius:8px;--shadow:0 1px 3px rgba(0,0,0,0.06)}
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',system-ui,sans-serif;background:var(--bg);color:var(--text);line-height:1.6;font-size:14px;padding:1.5rem;max-width:1600px;margin:0 auto}
.rpt-header{background:linear-gradient(135deg,#553C9A,#6B46C1);color:white;padding:1.5rem 2rem;border-radius:10px;margin-bottom:1.25rem;box-shadow:0 4px 6px rgba(0,0,0,0.07);display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:0.75rem}
.rpt-header h1{font-size:1.3rem;font-weight:700}
.rpt-header .sub{font-size:0.85rem;opacity:0.88}
.stats-bar{display:grid;grid-template-columns:repeat(5,1fr);gap:0.75rem;margin-bottom:1.25rem}
@media(max-width:700px){.stats-bar{grid-template-columns:repeat(2,1fr)}}
.stat-card{background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:0.75rem 1rem;text-align:center;box-shadow:var(--shadow)}
.stat-value{font-size:1.4rem;font-weight:700;color:#553C9A;line-height:1.2}
.stat-label{font-size:0.72rem;color:var(--text-secondary);text-transform:uppercase;letter-spacing:0.04em;margin-top:0.15rem}
.section-title{font-size:1rem;font-weight:600;color:var(--text);border-left:3px solid #553C9A;padding-left:0.75rem;margin:1.5rem 0 0.75rem}
.match-card{background:white;border:1px solid var(--border);border-radius:8px;padding:1rem 1.25rem;margin-bottom:0.75rem;box-shadow:var(--shadow)}
.match-header{display:flex;justify-content:space-between;align-items:start;flex-wrap:wrap;gap:0.5rem;margin-bottom:0.75rem}
.match-title{font-weight:700;font-size:1.05rem;color:#2D3748}
.match-class{display:inline-block;padding:2px 10px;border-radius:9999px;font-size:0.7rem;font-weight:600;background:#E9D8FD;color:#553C9A}
.match-meta{font-size:0.78rem;color:var(--text-secondary)}
.dim-grid{display:grid;grid-template-columns:repeat(4,1fr);gap:0.5rem;margin:0.75rem 0}
@media(max-width:700px){.dim-grid{grid-template-columns:repeat(2,1fr)}}
.dim-bar{margin-bottom:0.3rem}
.dim-label{display:flex;justify-content:space-between;font-size:0.72rem;margin-bottom:2px}
.dim-label span:first-child{color:var(--text-secondary)}
.dim-label span:last-child{font-weight:600;color:#2D3748}
.dim-track{background:#EDF2F7;border-radius:4px;height:7px;overflow:hidden}
.dim-fill{background:linear-gradient(90deg,#553C9A,#805AD5);height:100%;border-radius:4px;transition:width 0.3s}
.dim-fill-high{background:linear-gradient(90deg,#276749,#38A169)}
.dim-fill-medium{background:linear-gradient(90deg,#744210,#D69E2E)}
.dim-fill-low{background:linear-gradient(90deg,#C53030,#E53E3E)}
.dim-fill-none{background:#CBD5E0}
.reasoning-box{background:#FFFBEB;border:1px solid #F6E05E;border-radius:6px;padding:0.6rem 0.8rem;font-size:0.8rem;color:#744210;margin-top:0.5rem;line-height:1.5}
.evidence-list{font-size:0.78rem;color:var(--text-secondary);margin-top:0.4rem}
.evidence-list li{margin-bottom:2px}
.fragment-table{width:100%;font-size:0.75rem;border-collapse:collapse;margin:0.5rem 0}
.fragment-table th{background:#553C9A;color:white;padding:5px 8px;text-align:left;font-weight:600;font-size:0.68rem}
.fragment-table td{padding:4px 8px;border-bottom:1px solid #EDF2F7}
.fragment-table tr:nth-child(even){background:#FAFBFC}
footer{text-align:center;padding:1.5rem;color:var(--text-secondary);font-size:0.82rem;margin-top:1.5rem;border-top:1px solid var(--border)}
.pathway-card{margin-top:0.6rem;padding:0.6rem 0.75rem;background:#FFFAF0;border:1px solid #EDD6B0;border-radius:6px}
.pathway-title{font-size:0.82rem;font-weight:700;color:#7B341E;margin-bottom:0.4rem}
.pathway-section{margin-top:0.5rem}
.pathway-section>b{font-size:0.76rem;color:#9C4221;display:inline-block;margin-bottom:3px}
.cascade-sub{font-size:0.68rem;color:#A0AEC0;font-weight:400;margin-left:4px}
.cascade-box{border-radius:5px;padding:0.45rem 0.6rem;margin:0.35rem 0;font-size:0.76rem}
.cascade-box>b{font-size:0.78rem;color:#276749}
.cascade-sugar{background:#F0FFF4;border:1px solid #9AE6B4}
.cascade-ester{background:#FFF5F5;border:1px solid #FEB2B2}
.cascade-ester>b{color:#9B2C2C}
.cascade-chain{margin:0.35rem 0;line-height:1.9;word-break:break-word}
.cascade-chain .node{display:inline-block;background:white;border:1px solid #CBD5E0;border-radius:4px;padding:1px 6px;font-weight:600;font-size:0.74rem;color:#2D3748}
.cascade-chain .node.aglycone{background:#FEFCBF;border-color:#ECC94B;color:#744210;position:relative}
.cascade-chain .node-tag{display:block;font-size:0.58rem;font-weight:400;color:#975A16;line-height:1}
.cascade-chain .arrow{color:#4A5568;font-size:0.7rem;margin:0 2px}
.cascade-meta{font-size:0.7rem;color:#4A5568;margin-top:2px}
.subtype-badge{display:inline-block;background:#9B2C2C;color:white;padding:1px 8px;border-radius:9999px;font-size:0.66rem;font-weight:600;margin-left:4px}
.nl-explained{color:#276749;font-weight:600}
.nl-unexplained{color:#A0AEC0;font-style:italic}
.spec{display:inline-block;padding:1px 6px;border-radius:3px;font-size:0.64rem;font-weight:600}
.spec-high{background:#C6F6D5;color:#22543D}
.spec-med{background:#FEEBC8;color:#7B341E}
"""

    html_parts = [
        '<!DOCTYPE html>',
        f'<html lang="zh-CN">',
        '<head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
        f'<title>{mgf_name} — Structured Rule-Evidence Report</title>',
        '<style>', css, '</style></head><body>',

        # Header
        '<header class="rpt-header">',
        f'<div><h1>🧪 {mgf_name} — Structured Rule-Evidence Analysis</h1>',
        f'<div class="sub">Generated: {now} | structured rule-evidence report | '
        f'Total: {total} spectra | Matched: {matched_count} ({rate:.1f}%) | '
        f'Structured assessments: {structured_assessments}/{matched_count}</div></div>',
        '</header>',

        # Stats
        '<div class="stats-bar">',
        f'<div class="stat-card"><div class="stat-value">{matched_count}</div><div class="stat-label">Matched</div></div>',
        f'<div class="stat-card"><div class="stat-value">{structured_assessments}</div><div class="stat-label">Structured assessments</div></div>',
        f'<div class="stat-card"><div class="stat-value">{assessed_dimension_records}</div><div class="stat-label">Dimensions assessed</div></div>',
        f'<div class="stat-card"><div class="stat-value">{evidence_bearing_dimension_records}</div><div class="stat-label">Evidence-bearing dimensions</div></div>',
        f'<div class="stat-card"><div class="stat-value">{len(class_counts)}</div><div class="stat-label">Compound Classes</div></div>',
        '</div>',
    ]

    # Class distribution
    if class_counts:
        html_parts.append('<div class="section-title">Compound Class Distribution</div>')
        html_parts.append('<div style="display:flex;flex-wrap:wrap;gap:0.5rem;margin-bottom:1rem">')
        for cls_name, cnt in sorted(class_counts.items(), key=lambda x: -x[1]):
            html_parts.append(
                f'<span style="background:#E9D8FD;color:#553C9A;padding:4px 12px;'
                f'border-radius:9999px;font-size:0.78rem;font-weight:600">'
                f'{cls_name}: {cnt}</span>'
            )
        html_parts.append('</div>')

    # Match details
    html_parts.append(f'<div class="section-title">Match Details — Rule-Evidence Audit ({matched_count} entries)</div>')

    for i, r in enumerate(matched[:100], 1):  # Cap at 100 for reasonable size
        m = r["matches"][0]
        name = m.get("compound_name", "Unknown")
        smiles = m.get("smiles", "")
        classification = _classification_label(m)
        inchikey = m.get("inchikey", "")
        formula = m.get("formula", "")
        cosine = m.get("cosine", 0)
        final = m.get("final_score", 0)
        expert = m.get("expert_score")
        mass_ppm = m.get("mass_ppm", None)
        structural_img = _mol_to_base64(smiles) if smiles else ""

        # Expert sub-scores — support both top-level keys and expert_details
        ed = m.get("expert_details") or {}
        _get = lambda key: m.get(key, 0) or ed.get(key, 0)
        dims = [
            ("neutral_loss", "Neutral Loss", _get("neutral_loss_score"), 0.22),
            ("rule_coverage", "Rule Coverage", _get("rule_coverage_score"), 0.22),
            ("fg_signature", "FG Signature", _get("fg_signature_score"), 0.18),
            ("feature_library", "Class Feature Library", _get("feature_library_score"), 0.15),
            ("structure_fragment", "Structure-Frag", _get("structure_fragment_score"), 0.10),
            ("diagnostic_ion", "Diagnostic Ion", _get("diagnostic_ion_score"), 0.08),
            ("bond_energy", "Bond Energy", _get("bond_energy_score"), 0.03),
            ("isomer_differentiation", "Isomer Diff", _get("isomer_differentiation_score"), 0.02),
        ]
        assessments = ed.get("dimension_assessments") or {}
        fusion = ed.get("score_fusion") or {}

        assessment_values = [
            value for value in assessments.values() if isinstance(value, dict)
        ] if isinstance(assessments, dict) else []
        assessed_here = sum(
            str(value.get("status") or "unavailable") != "unavailable"
            for value in assessment_values
        )
        evidence_here = sum(
            str(value.get("status") or "") == "scored"
            and (
                int(value.get("evidence_count") or 0) > 0
                or float(value.get("score") or 0) > 0
            )
            for value in assessment_values
        )

        def _dim_color(s, status="scored"):
            if status in ("not_applicable", "unavailable", "no_evidence"):
                return "dim-fill-none"
            if s >= 0.6: return "dim-fill-high"
            if s >= 0.3: return "dim-fill-medium"
            return "dim-fill-low"

        html_parts.append('<div class="match-card">')
        html_parts.append('<div class="match-header">')
        html_parts.append(f'<div><span class="match-title">#{i}: {name[:80]}</span> '
                         f'<span class="match-class">{classification}</span></div>')
        ppm_text = f'{float(mass_ppm):.1f}' if mass_ppm is not None else 'N/A'
        not_applicable_here = sum(
            str(value.get("status") or "") == "not_applicable"
            for value in assessment_values
        )
        unavailable_here = sum(
            str(value.get("status") or "unavailable") == "unavailable"
            for value in assessment_values
        )
        html_parts.append(
            f'<div class="match-meta">m/z {float(r.get("precursor_mz") or 0):.4f} | '
            f'Cosine: {float(cosine or 0):.3f} | Final rank score: {float(final or 0):.3f} | '
            f'Rule dimensions with supporting observations: <b>{evidence_here}/{assessed_here}</b> | '
            f'Not applicable: {not_applicable_here} | Unavailable: {unavailable_here} | PPM: {ppm_text}'
        )
        html_parts.append(f'</div></div>')
        if expert is not None:
            html_parts.append(
                '<details style="margin:-.2rem 0 .55rem;font-size:.68rem;color:#718096">'
                '<summary>Legacy numeric audit field (not a probability; ranking contribution = 0)</summary>'
                f'<div style="padding:.25rem 0 0 .6rem">Uncalibrated aggregate: {float(expert):.3f}</div>'
                '</details>'
            )

        # Structure image + SMILES
        if structural_img or smiles:
            html_parts.append('<div style="display:flex;gap:1rem;align-items:start;margin-bottom:0.5rem">')
            if structural_img:
                html_parts.append(f'<img src="data:image/png;base64,{structural_img}" '
                                 f'style="border:1px solid #ddd;border-radius:4px;max-width:140px">')
            html_parts.append(f'<div style="font-size:0.75rem;color:var(--text-secondary);max-width:500px;word-break:break-all">'
                             f'<b>SMILES:</b> {smiles[:200]}<br>'
                             f'<b>Formula:</b> {formula}<br>'
                             f'<b>InChIKey:</b> {inchikey[:30]}...</div>')
            html_parts.append('</div>')

        # Fused-source provenance + 8-dimension applicability-aware chart
        fusion_components = fusion.get("components") or {}
        fusion_badges = []
        for source_key in ("base_8d", "v3_kb", "compound_pathway"):
            component = fusion_components.get(source_key) or {}
            if not component.get("available"):
                continue
            selected = source_key == fusion.get("selected_source")
            fusion_badges.append(
                f'<span style="display:inline-block;padding:2px 7px;margin:0 4px 5px 0;'
                f'border-radius:999px;border:1px solid {"#805AD5" if selected else "#CBD5E0"};'
                f'background:{"#FAF5FF" if selected else "#F7FAFC"};font-size:0.7rem">'
                f'{component.get("label_en") or source_key}'
                f'{" (recorded source)" if selected else ""}</span>'
            )
        if fusion_badges:
            html_parts.append(
                '<div style="margin:0.3rem 0">'
                + "".join(fusion_badges)
                + '<div style="font-size:0.68rem;color:#718096">'
                  'Legacy component values are shown for audit only. They are uncalibrated, '
                  'are not probabilities, and do not establish compound identity. '
                  'The bars below report applicability and observations for the base dimensions.</div></div>'
            )
        html_parts.append('<div class="dim-grid">')
        for dim_key, dim_name, dim_val, dim_weight in dims:
            assessment = assessments.get(dim_key) or {}
            status = assessment.get("status") or (
                "scored" if dim_val > 0 else "no_evidence"
            )
            pct = int(dim_val * 100)
            display_value = (
                "N/A" if status == "not_applicable"
                else "Unavailable" if status == "unavailable"
                else "No supporting observation" if status == "no_evidence"
                else f"{dim_val:.3f}"
            )
            html_parts.append(
                f'<div class="dim-bar">'
                f'<div class="dim-label"><span>{dim_name} (w={dim_weight})</span>'
                f'<span title="{assessment.get("reason", "")}">{display_value}</span></div>'
                f'<div class="dim-track"><div class="dim-fill {_dim_color(dim_val, status)}" style="width:{pct}%"></div></div>'
                f'</div>'
            )
        html_parts.append('</div>')

        # Chemical reasoning — check both top-level and expert_details
        reasoning = m.get("chemical_reasoning", "") or ed.get("chemical_reasoning", "")
        if reasoning:
            html_parts.append(f'<div class="reasoning-box"><b>🧠 Chemical Reasoning:</b> {reasoning}</div>')

        # 结构化碎裂途径 (级联/步骤表/诊断离子/中性丢失) — 此前仅压缩文本
        pathway_html = _render_fragmentation_pathway(ed)
        if pathway_html:
            html_parts.append(pathway_html)

        html_parts.append('</div>')

    html_parts.append(f'<footer>MSsystem Structured Rule-Evidence Report | Generated: {now} | '
                      f'Covering {len(class_counts)} compound classes</footer>')
    html_parts.append('</body></html>')

    html = "\n".join(html_parts)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
