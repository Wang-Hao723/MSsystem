"""
MSsystem Standard HTML Report Generator
Output: Flat 14-column table with structure images, filter/sort, print support.
Template format: # | m/z | Name | SMILES | Formula | Confidence | Cosine | Mass | Fragment | Intensity | Final | Adduct | Class | Structure
"""
import base64
import io
from datetime import datetime
from pathlib import Path
from matching.honest_calibration import annotate as _honest_annotate
from web.report_common import (
    _standard_compound_name,
    classification_evidence_status,
    exact_structure_confidence,
    final_claim_badge_html,
    isomer_candidate_projection,
    reportable_identity_confidence,
    reportable_identity_label,
    render_isomer_candidate_competition_html,
)
from web.raw_spectrum_locator import raw_spectrum_locator

COSINE_THRESHOLD = 0.7
FINAL_SCORE_THRESHOLD = 0.5


def _attr(s):
    """转义 HTML 属性值 —— 防止化合物名含引号/尖括号破坏 data-* 属性。
    若不转义, data-name="..." 遇到内嵌引号会解析错位, 导致该行的
    dataset.name 读不到 → 搜索该化合物时无结果 (修复 v3.3.1)"""
    return (str(s).replace("&", "&amp;").replace('"', "&quot;")
            .replace("<", "&lt;").replace(">", "&gt;"))


def _mol_to_img(smiles, size=(300, 220)):
    """SMILES → inline SVG or base64 PNG structure image tag. Returns HTML string."""
    if not smiles or str(smiles).upper() in ("N/A", "NONE", "NULL", "", "NONE"):
        return ""
    smiles = str(smiles).strip()
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        from rdkit.Chem import Draw, AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            try:
                from web.structure import _sanitize_smiles
                sanitized = _sanitize_smiles(smiles)
                if sanitized:
                    mol = Chem.MolFromSmiles(sanitized)
            except Exception:
                pass
        if mol is None:
            return ""
        AllChem.Compute2DCoords(mol)

        # Try PNG first (needs Cairo)
        try:
            img = Draw.MolToImage(mol, size=size)
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            b64 = base64.b64encode(buf.getvalue()).decode("ascii")
            return f"data:image/png;base64,{b64}"
        except Exception:
            pass

        # Fall back to SVG (no Cairo needed)
        drawer = Draw.MolDraw2DSVG(size[0], size[1])
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        # Strip XML declaration for inline use
        svg = svg.replace("<?xml version='1.0' encoding='iso-8859-1'?>\n", "")
        # Return as base64 SVG data URI
        b64 = base64.b64encode(svg.encode("utf-8")).decode("ascii")
        return f"data:image/svg+xml;base64,{b64}"
    except Exception:
        return ""


def _badge_css(label):
    """P4 诚实性校准: 置信度标签 → CSS class (单一映射, 配合 honest_calibration)。"""
    return {"High": "badge-high", "RT-confirmed High": "badge-high",
            "Class High": "badge-high", "Isomer-group High": "badge-high",
            "Strong Candidate": "badge-medium",
            "Medium": "badge-medium", "Low": "badge-low",
            "Putative": "badge-low"}.get(label, "badge-unknown")


def _score_color(s):
    if s >= 0.6: return "#1a5f2a"
    if s >= 0.4: return "#8b6914"
    return "#666666"


def _class_color(cls):
    if not cls: return "#555"
    c = str(cls).lower()
    if "alkaloid" in c: return "#7b2d8b"
    if "flavonoid" in c: return "#c4881a"
    if any(x in c for x in ("phenyl", "terpen", "lipid")): return "#2a6b8b"
    if any(x in c for x in ("phenol", "acid", "amino")): return "#1a6b2a"
    return "#333"


def _msms_fragments(query_peaks, lib_peaks, tol=0.02, top_n=20):
    """Build semicolon-separated MS/MS fragment m/z string from matched query peaks."""
    if not query_peaks or not lib_peaks:
        return ""
    lib_mz_list = [lp[0] for lp in lib_peaks]
    matched = []
    for qp in query_peaks:
        for lmz in lib_mz_list:
            if abs(qp[0] - lmz) <= tol:
                matched.append(f"{qp[0]:.4f}")
                break
    return "; ".join(matched[:top_n]) if matched else ""


def generate_html(
    results,
    output_path,
    mgf_name,
    *,
    high_confidence_scope_groups=None,
    high_confidence_exact_structures=None,
    high_confidence_isomer_groups=None,
):
    """Generate standard HTML report — flat table matching template format with filter/sort/print."""

    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        pass

    total = len(results)
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    matched_all = [r for r in results if r.get("matches")]
    unmatched_list = [r for r in results if not r.get("matches")]
    matched_all.sort(key=lambda r: (
        -float(
            r["matches"][0].get(
                "rank_score", r["matches"][0].get("final_score", 0)
            ) or 0
        ),
        -float(r["matches"][0].get("final_score", 0) or 0),
        -float(r["matches"][0].get("cosine", 0) or 0),
    ))

    def identity_key(row):
        match = row["matches"][0]
        return (
            match.get("inchikey")
            or match.get("smiles")
            or match.get("compound_name")
            or f"mz_{row['precursor_mz']:.4f}"
        )

    # Dedup: keep the raw-rank representative per compound.
    deduped = []
    seen_keys = set()
    for r in matched_all:
        key = identity_key(r)
        if key and key not in seen_keys:
            seen_keys.add(key)
            deduped.append(r)
    matched_list = deduped

    matched_count = len(matched_list)
    unmatched_count = total - matched_count
    dup_count = len(matched_all) - len(matched_list)
    rate_str = f"{matched_count/total*100:.1f}%" if total else "N/A"

    # Count High at the identity scope actually reportable by each row.
    high_supported_conf = len({
        identity_key(r)
        for r in matched_all
        if reportable_identity_confidence(r["matches"][0]) == "High"
    })
    # A user-facing High result is an identification *group*: either one
    # supported 2D structure or one unresolved High-confidence isomer group.
    # Keep the legacy exact-candidate count only as a compatibility fallback
    # for callers that do not have the compound-conclusion projection.
    high_scope_count = (
        max(0, int(high_confidence_scope_groups))
        if high_confidence_scope_groups is not None
        else high_supported_conf
    )
    high_exact_count = (
        max(0, int(high_confidence_exact_structures))
        if high_confidence_exact_structures is not None
        else high_supported_conf
    )
    high_isomer_count = (
        max(0, int(high_confidence_isomer_groups))
        if high_confidence_isomer_groups is not None
        else max(0, high_scope_count - high_exact_count)
    )

    # Class filter options
    classes_set = sorted(set(
        str(m["matches"][0].get("classification", "Unknown"))
        for m in matched_list
        if m["matches"][0].get("classification")
    ))

    css = """
:root{--primary:#2C7A7B;--primary-dark:#1A5253;--primary-light:#E6FFFA;--text:#1A202C;--text-secondary:#718096;--border:#E2E8F0;--bg:#F7FAFC;--card-bg:#FFFFFF;--radius:8px;--shadow:0 1px 3px rgba(0,0,0,0.06)}
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',system-ui,sans-serif;background:var(--bg);color:var(--text);line-height:1.6;font-size:15px;padding:1.5rem;max-width:1700px;margin:0 auto}
@media(max-width:768px){body{padding:0.75rem;font-size:14px}}
.rpt-header{background:linear-gradient(135deg,#1a5276,#2e86ab);color:white;padding:1.5rem 2rem;border-radius:10px;margin-bottom:1.25rem;box-shadow:0 4px 6px rgba(0,0,0,0.07);display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:0.75rem}
.rpt-header h1{font-size:1.4rem;font-weight:700;letter-spacing:-0.02em}
.rpt-header .sub{font-size:0.9rem;opacity:0.88;line-height:1.4}
.btn-print{background:rgba(255,255,255,0.2);color:white;border:1px solid rgba(255,255,255,0.4);padding:0.5rem 1.2rem;border-radius:6px;cursor:pointer;font-size:0.9rem;font-weight:600;font-family:inherit;transition:all 0.2s}
.btn-print:hover{background:rgba(255,255,255,0.3)}
.stats-bar{display:grid;grid-template-columns:repeat(4,1fr);gap:0.75rem;margin-bottom:1.25rem}
@media(max-width:600px){.stats-bar{grid-template-columns:repeat(2,1fr)}}
.stat-card{background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:0.75rem 1rem;text-align:center;box-shadow:var(--shadow)}
.stat-value{font-size:1.5rem;font-weight:700;color:#1a5276;line-height:1.2}
.stat-label{font-size:0.75rem;color:var(--text-secondary);text-transform:uppercase;letter-spacing:0.05em;margin-top:0.15rem}
.stat-detail{font-size:0.68rem;color:var(--text-secondary);margin-top:0.16rem;line-height:1.35}
.toolbar{background:white;padding:0.7rem 1rem;border-radius:8px;margin-bottom:0.75rem;box-shadow:var(--shadow);border:1px solid var(--border);display:flex;gap:0.6rem;flex-wrap:wrap;align-items:center;font-size:0.88rem}
.toolbar label{font-weight:600;color:var(--text-secondary);white-space:nowrap;font-size:0.82rem}
.toolbar select,.toolbar input{padding:0.4rem 0.55rem;border:1px solid var(--border);border-radius:4px;font-size:0.85rem;font-family:inherit;background:white}
.toolbar select:focus,.toolbar input:focus{outline:none;border-color:#2e86ab}
.toolbar .result-count{margin-left:auto;color:var(--text-secondary);font-size:0.82rem}
.table-wrap{background:white;border:1px solid var(--border);border-radius:8px;overflow:hidden;box-shadow:var(--shadow);overflow-x:auto}
table{width:100%;border-collapse:collapse;font-size:0.82rem}
#matchTable{min-width:2240px;table-layout:fixed}
thead{background:linear-gradient(135deg,#1a5276,#2e86ab);position:sticky;top:0;z-index:10}
th{color:white;padding:9px 6px;text-align:left;font-weight:600;font-size:0.72rem;text-transform:uppercase;letter-spacing:0.04em;white-space:nowrap}
td{padding:7px 6px;border-bottom:1px solid #EDF2F7;vertical-align:middle}
tbody tr{transition:background 0.12s}
tbody tr:hover{background:#e3f2fd}
tbody tr:nth-child(even){background:#FAFBFC}
tbody tr:nth-child(even):hover{background:#e3f2fd}
.col-idx{color:var(--text-secondary);text-align:center;width:34px;font-size:0.76rem}
.col-num{text-align:center;font-family:'SF Mono','Fira Code',monospace;font-size:0.78rem}
.col-name{font-weight:600;width:190px;max-width:190px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.confidence-cell{white-space:normal;line-height:1.3}
.confidence-cell small{display:block;margin-top:3px;color:var(--text-secondary);font-size:0.62rem}
.col-smiles{font-family:'SF Mono','Fira Code',monospace;font-size:0.72rem;color:var(--text-secondary);max-width:130px;word-break:break-all}
.col-mol{text-align:center}
.col-mol img{max-width:105px;max-height:80px;border:1px solid var(--border);border-radius:4px;padding:2px;display:inline-block}
.col-msms{font-family:'SF Mono','Fira Code',monospace;font-size:0.7rem;color:#555;max-width:160px;word-break:break-all}
.badge{display:inline-block;padding:2px 9px;border-radius:9999px;font-size:0.68rem;font-weight:600;white-space:nowrap}
.badge-high{background:#C6F6D5;color:#276749}
.badge-medium{background:#FEFCBF;color:#744210}
.badge-low{background:#FEEBC8;color:#7B341E}
.badge-unknown{background:#EDF2F7;color:#4A5568}
.sc-green{color:#1a5f2a;font-weight:700}
.sc-amber{color:#8b6914;font-weight:700}
.sc-gray{color:#666}
.section-title{font-size:1rem;font-weight:600;color:var(--text);border-left:3px solid var(--primary);padding-left:0.75rem;margin:1.25rem 0 0.75rem}
.unmatched-section{margin-top:1.25rem}
.unmatched-section h2{font-size:0.95rem;color:var(--text);margin-bottom:0.75rem;border-left:3px solid var(--primary);padding-left:0.75rem}
.no-match-row td{color:#A0AEC0;font-style:italic}
@media print{
  body{padding:0.4in;background:white;font-size:10px}
  .rpt-header{background:#1a5276!important;-webkit-print-color-adjust:exact;print-color-adjust:exact}
  .toolbar{display:none}
  .btn-print{display:none}
  table{font-size:0.63rem}
  th{font-size:0.58rem}
  td{padding:3px 4px}
}
footer{text-align:center;padding:1.5rem;color:var(--text-secondary);font-size:0.82rem;margin-top:1.5rem;border-top:1px solid var(--border)}
"""

    # ── Build HTML ──────────────────────────────────────────────────
    html_parts = [
        '<!DOCTYPE html>',
        '<html lang="zh-CN">',
        '<head>',
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width,initial-scale=1.0">',
        f'<title>{mgf_name} — LC-MS/MS Analysis Report</title>',
        '<style>', css, '</style>',
        '</head>',
        '<body>',

        # Header
        '<header class="rpt-header">',
        f'<div><h1>{mgf_name} — Analysis Report</h1>'
        f'<div class="sub">Generated: {now} &nbsp;|&nbsp; '
        f'Total: {total} spectra &nbsp;|&nbsp; '
        f'Matched: {matched_count} ({rate_str}) &nbsp;|&nbsp; '
        f'{("Dedup from " + str(len(matched_all)) + " | ") if dup_count > 0 else ""}'
        f'Thresholds: Cosine &ge; {COSINE_THRESHOLD}, Final Score &ge; {FINAL_SCORE_THRESHOLD}</div></div>',
        '<button class="btn-print" onclick="window.print()">Print / Save PDF</button>',
        '</header>',

        # Stats
        '<div class="stats-bar">',
        f'<div class="stat-card"><div class="stat-value">{total}</div><div class="stat-label">Total Spectra</div></div>',
        f'<div class="stat-card"><div class="stat-value">{matched_count}</div><div class="stat-label">Matched ({rate_str})</div></div>',
        f'<div class="stat-card"><div class="stat-value">{high_scope_count}</div><div class="stat-label">High 级鉴定组</div>'
        f'<div class="stat-detail">单一二维连接 {high_exact_count} · 未解异构体组 {high_isomer_count}</div></div>',
        f'<div class="stat-card"><div class="stat-value">{unmatched_count}</div><div class="stat-label">Unmatched</div></div>',
        '</div>',
    ]

    # ── Toolbar ──
    if matched_list:
        html_parts.append('<div class="toolbar">')
        html_parts.append('<label>Sort:</label>')
        html_parts.append('<select id="sortSelect" onchange="applyFilters()">'
                          '<option value="final_desc">Final Score &darr;</option>'
                          '<option value="cosine_desc">Cosine &darr;</option>'
                          '<option value="name_asc">Name A-Z</option>'
                          '</select>')
        html_parts.append('<label style="margin-left:6px">Filter name:</label>')
        html_parts.append('<input type="text" id="nameFilter" placeholder="Compound..." oninput="applyFilters()" style="width:150px">')
        html_parts.append('<label style="margin-left:6px">Confidence:</label>')
        html_parts.append('<select id="confFilter" onchange="applyFilters()">'
                          '<option value="all">All</option>'
                          '<option value="High">High</option>'
                          '<option value="Medium">Medium</option>'
                          '<option value="Low">Low</option>'
                          '</select>')
        html_parts.append('<label style="margin-left:6px">Class:</label>')
        html_parts.append('<select id="classFilter" onchange="applyFilters()"><option value="all">All</option>')
        for cls in classes_set[:50]:
            html_parts.append(f'<option value="{cls}">{cls}</option>')
        html_parts.append('</select>')
        html_parts.append('<span class="result-count" id="resultCount"></span>')
        html_parts.append('</div>')

    # ── Matched Table (flat, no expand) ──
    if matched_list:
        source_file = raw_spectrum_locator({}, mgf_name)["source_file"] or mgf_name
        html_parts.append(
            '<div style="margin:0.75rem 0;padding:0.65rem 0.8rem;border:1px solid #9cc9c6;'
            'background:#f0f8f7;border-radius:7px;font-size:0.78rem;line-height:1.55">'
            '<b>原始 MS2 人工复核：</b>以下碎片来自查询谱，不是库谱。优先在源文件 '
            f'<span style="font-family:monospace">{_attr(source_file)}</span> 中按 '
            '<b>MS2 谱图编号/native ID</b> 定位；若查看软件不支持编号检索，再联合使用 '
            '<b>RT（min）+ 前体 m/z</b> 查找，并以碎片峰表复核。</div>'
        )
        html_parts.append(
            f'<div class="section-title">Matched Spectra ({matched_count})</div>'
        )
        html_parts.append('<div class="table-wrap"><table id="matchTable"><thead><tr>'
                          '<th style="width:34px">#</th>'
                          '<th style="width:90px">Query m/z</th>'
                          '<th style="width:190px">Compound Name</th>'
                          '<th style="width:130px">SMILES</th>'
                          '<th style="width:90px">Formula</th>'
                          '<th style="width:140px;text-align:center">可报告置信度</th>'
                          '<th style="width:102px;text-align:center">鉴定结论</th>'
                          '<th style="width:65px;text-align:center">Cosine</th>'
                          '<th style="width:55px;text-align:center">Mass</th>'
                          '<th style="width:60px;text-align:center">Fragment</th>'
                          '<th style="width:60px;text-align:center">Intensity</th>'
                          '<th style="width:60px;text-align:center">Final</th>'
                          '<th style="width:58px;text-align:center">RT (min)</th>'
                          '<th style="width:65px;text-align:center">Peak_Height</th>'
                          '<th style="width:65px;text-align:center">Peak_Area</th>'
                          '<th style="width:70px">Adduct</th>'
                          '<th style="width:90px">Class</th>'
                          '<th style="width:120px">Dataset_Source</th>'
                          '<th style="width:175px">MS/MS · 原始谱图定位</th>'
                          '<th style="width:165px">结构推理 / Top-K</th>'
                          '<th style="width:110px;text-align:center">Structure</th>'
                          '</tr></thead><tbody id="tableBody">')

        for i, r in enumerate(matched_list):
            m = r["matches"][0]
            query_mz = r["precursor_mz"]
            smiles = m.get("smiles", "") or ""
            mol_img = _mol_to_img(smiles, size=(200, 160))
            # Preserve the complete-identity tier while displaying the tier at
            # the row's explicitly reportable identity scope.
            _ann = _honest_annotate(m, m.get("top_candidates"))
            conf = reportable_identity_confidence(m)
            conf_label = reportable_identity_label(m, language="zh")
            full_identity_conf = exact_structure_confidence(m)
            gate_required = bool(
                _ann.get("precursor_high_evidence_gate_required")
            )
            gate_reason = str(
                _ann.get("precursor_high_evidence_gate_reason_zh")
                or _ann.get("precursor_high_evidence_gate_reason")
                or ""
            ).strip()
            badge_cls = _badge_css(conf)
            formula_str = m.get("formula", "") or "-"
            class_evidence = classification_evidence_status(m)
            classification = str(class_evidence["label"])
            if class_evidence["status"] == "tentative":
                classification = f"待核查：{classification}"
            dataset_source = str(m.get("dataset_source", "") or "")
            fs = m.get("final_score", 0)
            cs = m.get("cosine", 0)
            mass_s = m.get("mass_score", 0)
            frag_s = m.get("fragment_score", 0)
            inten_s = m.get("intensity_score", 0)
            tic = r.get("tic", 0) or 0
            base_peak = r.get("base_peak", 0) or 0
            locator = raw_spectrum_locator(r, mgf_name)
            rt_minutes = locator["retention_time_min"]
            rt_display = f"{rt_minutes:.3f}" if rt_minutes is not None else "-"
            spectrum_id = _attr(locator["spectrum_id"] or "未记录")
            adduct = str(m.get("adduct", "") or "-")
            name = (
                _standard_compound_name(m.get("compound_name", ""), m)
                or "Unnamed compound"
            )[:70]
            name_cell_text = _attr(name[:55])
            badge_cell_html = (
                f'<span class="badge {badge_cls}" '
                f'title="可报告身份置信度；完整身份置信度：{_attr(full_identity_conf)}">'
                f'{_attr(conf_label)}</span>'
                + (
                    f'<small>完整身份：{_attr(full_identity_conf)}</small>'
                    if conf != full_identity_conf else ""
                )
            )
            final_claim_cell_html = final_claim_badge_html(m, language="zh")

            # MS/MS fragment summary
            qpeaks = r.get("query_peaks", [])
            lpeaks = m.get("lib_peaks", [])
            msms_str = _msms_fragments(qpeaks, lpeaks) if lpeaks else ""
            msms_display = msms_str if msms_str else "-"

            # Score color classes
            cs_cls = "sc-green" if cs >= 0.9 else ("sc-amber" if cs >= 0.8 else "sc-gray")
            fs_cls = "sc-green" if fs >= 0.7 else ("sc-amber" if fs >= 0.5 else "sc-gray")

            # Structure image cell
            img_tag = (
                f"<img src='{mol_img}' alt='Structure' "
                f"style='max-width:100px;max-height:75px'>"
            ) if mol_img else '<span style="color:#ccc;font-size:0.7rem">-</span>'

            # SMILES truncated
            smiles_short = (smiles[:45] + "...") if len(smiles) > 45 else (smiles or "-")

            # v3 Phase 4/5: 结构推理 + Top-K 近邻 (从 expert_details.v3 + top_candidates)
            reason_lines = []
            v3 = (m.get("expert_details") or {}).get("v3") or {}
            for line in (v3.get("reasoning") or [])[:4]:
                line = _attr(str(line))[:90]
                if (
                    "异构" not in line
                    and "isomer" not in line.lower()
                    and (
                        "母核" in line
                        or "诊断离子" in line
                        or "推测" in line
                        or "生源" in line
                    )
                ):
                    reason_lines.append(line)
            isopair = m.get("isopair_evidence") or {}
            if isopair.get("evaluated"):
                isopair_action = (
                    "2-D connectivity 已重排"
                    if isopair.get("ranking_changed")
                    else "2-D connectivity 排名未改变"
                )
                gate_result = str(
                    _ann.get("isopair_high_precision_gate_result")
                    or "evaluated_not_high_candidate"
                )
                reason_lines.append(
                    _attr(
                        f"IsoPair-MS: {isopair_action}; High gate={gate_result}; "
                        "scope=InChIKey-14, no stereochemistry claim"
                    )
                )
            top_c = m.get("top_candidates") or []
            isomer_projection = isomer_candidate_projection(m)
            # Unresolved isomers use the same complete member set as the
            # result page and CSV.  Ordinary rows retain a compact Top-K hint.
            runup_parts = []
            alternative_candidates = (
                isomer_projection.get("competing_candidates") or []
                if isomer_projection.get("status") != "not_applicable"
                else top_c[1:5]
            )
            for c in alternative_candidates:
                cn = _standard_compound_name(
                    c.get("compound_name", ""), c
                )
                if cn and cn.lower() != name.lower():
                    runup_parts.append(
                        f"#{int(c.get('rank') or len(runup_parts) + 2)} {_attr(cn)} "
                        f"(cos {float(c.get('cosine') or 0):.3f}; "
                        f"final {float(c.get('final_score') or 0):.3f})"
                    )
            runup_label = (
                f"异构体竞争候选（其余 {isomer_projection.get('competitive_candidate_count', 0)} 个）"
                if isomer_projection.get("status") != "not_applicable"
                else "其他谱库候选"
            )
            runup_html = (
                f'<span style="color:#1a5276">↳ {runup_label}: '
                f'{" · ".join(runup_parts)}</span>'
                if runup_parts else ""
            )
            reason_cell = ("<br>".join(reason_lines) +
                           (f"<br>{runup_html}" if runup_html else "")) or "—"
            if gate_required:
                gate_note = (
                    '<span style="display:block;margin-top:5px;padding:4px 6px;'
                    'border-left:3px solid #c47b00;background:#fff8e8;color:#7a4b00">'
                    '<b>High 证据门控：</b>该候选已降为 Medium。'
                    f'{_attr(gate_reason)}</span>'
                )
                reason_cell = (
                    gate_note
                    if reason_cell == "—"
                    else f"{reason_cell}<br>{gate_note}"
                )

            html_parts.append(
                f'<tr data-idx="{i}" data-cosine="{cs:.4f}" data-final="{fs:.4f}" '
                f'data-name="{_attr(name.lower())}" data-class="{_attr(classification)}" data-conf="{_attr(conf)}">'
                f'<td class="col-idx">{i+1}</td>'
                f'<td class="col-num">{query_mz:.4f}</td>'
                f'<td class="col-name" title="{_attr(name)}">{name_cell_text}</td>'
                f'<td class="col-smiles" title="{smiles}">{smiles_short}</td>'
                f'<td style="font-family:monospace;font-size:0.8rem">{formula_str}</td>'
                f'<td class="confidence-cell" style="text-align:center">{badge_cell_html}</td>'
                f'<td style="text-align:center">{final_claim_cell_html}</td>'
                f'<td class="col-num {cs_cls}">{cs:.4f}</td>'
                f'<td class="col-num" style="color:{_score_color(mass_s)}">{mass_s:.3f}</td>'
                f'<td class="col-num" style="color:{_score_color(frag_s)}">{frag_s:.3f}</td>'
                f'<td class="col-num" style="color:{_score_color(inten_s)}">{inten_s:.3f}</td>'
                f'<td class="col-num" style="font-weight:700;color:#1a5276">{fs:.4f}</td>'
                f'<td class="col-num" style="font-size:0.78rem">{rt_display}</td>'
                f'<td class="col-num" style="font-size:0.78rem">{base_peak:,.0f}</td>'
                f'<td class="col-num" style="font-size:0.78rem">{tic:,.0f}</td>'
                f'<td style="font-size:0.78rem">{adduct[:15]}</td>'
                f'<td><span style="background:{_class_color(classification)};color:white;padding:2px 7px;border-radius:4px;font-size:0.7rem;white-space:nowrap">{classification[:22]}</span></td>'
                f'<td style="font-size:0.72rem;color:#555" title="{dataset_source}">{dataset_source[:20] or "—"}</td>'
                '<td class="col-msms">'
                f'<div style="font-size:0.68rem;color:#176b68;font-weight:700;margin-bottom:3px">'
                f'RT {rt_display} min · MS2 {spectrum_id}</div>'
                f'<div title="{msms_str}">{msms_display[:80]}</div></td>'
                f'<td style="font-size:0.7rem;color:#444;line-height:1.35">{reason_cell}</td>'
                f'<td class="col-mol">{img_tag}</td>'
                f'</tr>'
            )

        html_parts.append('</tbody></table></div>')

        # The complete candidate set stays in the offline file, while native
        # disclosures keep the on-screen report compact and readable.
        html_parts.append(
            render_isomer_candidate_competition_html(
                matched_list, language="zh"
            )
        )

    # ── Unmatched Section (含 MSClassify 谱图判类) ──
    if unmatched_list:
        # 统计未命中谱图的 MSClassify 判类 (display_class), 区分高/低置信
        from collections import Counter as _Ctr
        _cls_tot = _Ctr()
        _n_hi = _n_lo = 0
        for r in unmatched_list:
            c = r.get("display_class")
            if not c:
                continue
            _cls_tot[c] += 1
            _cf = float(r.get("dl_confidence") or 0)
            if _cf >= 0.30:
                _n_hi += 1
            else:
                _n_lo += 1
        _n_classified = sum(_cls_tot.values())
        _summary_html = ""
        if _n_classified:
            _summary_html = (
                '<div style="margin:0.5rem 0 0.75rem;font-size:0.82rem">'
                f'🧪 <b>谱图判类汇总</b> (MSClassify 从 MS2 推断, 非库匹配; '
                f'{_n_classified}/{unmatched_count} 有判类 · {_n_hi} 高置信 / '
                f'<span style="color:#b9770e">{_n_lo} 低置信</span>):<br>'
                + ''.join(
                    f'<span style="display:inline-block;margin:2px 4px;padding:2px 9px;'
                    f'border-radius:10px;background:#e8f0fe;color:#1a56c4;font-size:0.74rem">'
                    f'{_attr(c)} · {n}</span>' for c, n in _cls_tot.most_common()
                ) + '</div>'
            )
        html_parts.append('<div class="unmatched-section">')
        html_parts.append(
            f'<h2>Unmatched Spectra ({unmatched_count})'
            + (f' — 含 MSClassify 谱图判类 ({_n_classified} 条有判类)</h2>' if _n_classified else '</h2>')
        )
        if _summary_html:
            html_parts.append(_summary_html)
        html_parts.append('<div class="table-wrap"><table><thead><tr>'
                          '<th>#</th><th>Query m/z</th><th>RT (min)</th><th>MS2 spectrum/native ID</th>'
                          '<th>Best Cosine</th><th>预测分类 (MSClassify)</th><th>置信度</th>'
                          '</tr></thead><tbody>')
        for i, r in enumerate(unmatched_list):
            best = r.get("best_cosine", 0)
            _dc = r.get("display_class") or ""
            _cf = r.get("dl_confidence")
            if _dc:
                _cf_val = float(_cf) if _cf is not None else 0.0
                _cls_cell = (f'<span style="background:{_class_color(_dc)};color:white;'
                             f'padding:2px 7px;border-radius:4px;font-size:0.7rem">{_attr(_dc[:22])}</span>')
                _conf_cell = (f'{_cf_val:.2f} ' + ('<span style="color:#1a5f2a">高</span>'
                              if _cf_val >= 0.30 else '<span style="color:#b9770e">低</span>'))
            else:
                _cls_cell = '<span style="color:#bbb">—</span>'
                _conf_cell = '—'
            _locator = raw_spectrum_locator(r, mgf_name)
            _rt = _locator["retention_time_min"]
            _rt_cell = f"{_rt:.3f}" if _rt is not None else "-"
            html_parts.append(
                f'<tr class="no-match-row">'
                f'<td class="col-idx">{matched_count + i + 1}</td>'
                f'<td class="col-num">{r["precursor_mz"]:.4f}</td>'
                f'<td class="col-num">{_rt_cell}</td>'
                f'<td class="col-num">{_attr(_locator["spectrum_id"] or "未记录")}</td>'
                f'<td class="col-num">{best:.3f}</td>'
                f'<td>{_cls_cell}</td>'
                f'<td style="font-size:0.74rem">{_conf_cell}</td>'
                f'</tr>'
            )
        html_parts.append('</tbody></table></div></div>')

    # ── Footer ──
    html_parts.append(
        f'<footer>Generated by MSsystem &middot; LC-MS/MS Spectral Matching &middot; {now}</footer>'
    )

    # ── JavaScript: Filter & Sort ──
    if matched_list:
        html_parts.append('''<script>
// 缓存原始行数据 —— 脚本加载时执行一次, 作为稳定的过滤/排序数据源。
// 修复 bug: 旧版每次从当前 DOM 读取行, 搜索一次后 DOM 只剩过滤结果,
// 再搜就在残缺子集上操作 → 搜不到 + 无法连续搜索。
var ALL_ROWS = (function() {
  var tbody = document.getElementById("tableBody");
  if (!tbody) return [];
  return Array.from(tbody.querySelectorAll("tr")).map(function(r) {
    return {
      row: r,
      cosine: parseFloat(r.dataset.cosine),
      final: parseFloat(r.dataset.final),
      name: r.dataset.name || "",
      cls: r.dataset.class || "",
      conf: r.dataset.conf || ""
    };
  });
})();
var TOTAL_COUNT = ALL_ROWS.length;

function applyFilters() {
  var sortBy = document.getElementById("sortSelect").value;
  var nameFilter = (document.getElementById("nameFilter").value || "").toLowerCase().trim();
  var confFilter = document.getElementById("confFilter").value;
  var classFilter = document.getElementById("classFilter").value;
  var tbody = document.getElementById("tableBody");

  // 始终从缓存的完整数据集过滤, 绝不依赖当前 DOM 内容
  var filtered = ALL_ROWS.filter(function(i) {
    if (nameFilter && i.name.indexOf(nameFilter) < 0) return false;
    if (confFilter !== "all" && i.conf !== confFilter) return false;
    if (classFilter !== "all" && i.cls !== classFilter) return false;
    return true;
  });

  if (sortBy === "final_desc") filtered.sort(function(a,b) { return b.final - a.final; });
  else if (sortBy === "cosine_desc") filtered.sort(function(a,b) { return b.cosine - a.cosine; });
  else if (sortBy === "name_asc") filtered.sort(function(a,b) { return a.name.localeCompare(b.name); });

  // 重编号
  filtered.forEach(function(item, idx) {
    item.row.querySelector(".col-idx").textContent = idx + 1;
  });

  // 按过滤+排序结果重排 DOM
  while (tbody.firstChild) tbody.removeChild(tbody.firstChild);
  filtered.forEach(function(item) { tbody.appendChild(item.row); });

  document.getElementById("resultCount").textContent = "Showing " + filtered.length + " of " + TOTAL_COUNT + " matches";
}
applyFilters();
</script>''')

    html_parts.append('</body></html>')

    # ── Write ──
    html = "\n".join(html_parts)
    out = Path(output_path) if not isinstance(output_path, Path) else output_path
    out.write_text(html, encoding="utf-8")
