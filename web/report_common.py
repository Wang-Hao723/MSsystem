"""web/report_common.py — 所有 Web 输出(standard/detailed HTML · PDF · CSV · Excel)的共享真源。

把此前散落在 server.py / app.py / report.py 里各自重复、口径漂移的三类逻辑收口到此:
  1. 置信度: 统一走 matching.honest_calibration (化合物级 + 类级双口径), 取代各处裸 m["confidence"]
     与 report.py 旧 KPI 定义 (cos>=0.9 & final>=0.5)。
  2. 裂解途径: render_db_pathway 直接消费 Compounds.fragmentation_pathway 预存 JSON
     (scripts/mine_pathways.py 挖出的化合物馆藏途径), 而非每次临时算的 expert_details。
  3. CFM-ID 碎裂: render_cfm_cascade 解耦展示 —— CFM 跑出预测产物就画(前体→预测产物),
     不再以"预测碎片质量与观测 m/z <=0.5Da"为是否显示的门; 命中观测峰仅作绿标注释。

零排序风险: 全部为呈现层纯函数, 不触碰 final_score / 匹配结果 → 红线 11/14 不受影响。
"""
from __future__ import annotations

import json
import logging
import re
import subprocess
from functools import lru_cache
from html import escape as _html_escape
from typing import Any, Mapping, Optional

from web.governance_blocked_projection import (
    GOVERNANCE_BLOCKED_EXPORT_COLUMNS,
    append_governance_blocked_csv_section,
    governance_blocked_export_rows,
    inject_governance_blocked_html,
    project_governance_blocked_candidates,
    render_governance_blocked_html,
)

logger = logging.getLogger(__name__)

# ── 视觉常量 (与 detailed 报告模板 A 一致的深蓝主色 + 碎裂橙/前体绿) ──
COLOR_PRIMARY = "#1a5276"
COLOR_PRECURSOR = "#22543D"
COLOR_FRAGMENT = "#7B341E"
_BADGE_CSS = {"High": "badge-high", "RT-confirmed High": "badge-high",
              "Class High": "badge-high", "Isomer-group High": "badge-high",
              "Strong Candidate": "badge-medium",
              "Medium": "badge-medium", "Low": "badge-low",
              "Putative": "badge-low", "Unknown": "badge-unknown"}


# ══════════════════════════════════════════════════════════════════
# Docker 探测 (进程级缓存, 避免每次报告都 shell out)
# ══════════════════════════════════════════════════════════════════
@lru_cache(maxsize=1)
def docker_available() -> bool:
    """docker daemon 可用 + cfmid 镜像在位 → True (进程级缓存)。"""
    try:
        out = subprocess.run(["docker", "images", "kykeonanalytics/cfmid-arm64"],
                             capture_output=True, text=True, timeout=8)
        if out.returncode != 0:
            return False
        # 输出含镜像名行 = 镜像存在
        return "kykeonanalytics/cfmid-arm64" in out.stdout
    except Exception:
        return False


# ══════════════════════════════════════════════════════════════════
# 置信度 (honest_calibration 单一真源)
# ══════════════════════════════════════════════════════════════════
def honest_annotation(m: Optional[dict]) -> dict:
    """对一条 match 计算诚实性校准注解 (化合物级/类级置信 + 异构集 + 诚实标题)。

    幂等: 若 match 已含 P4 序列化字段(server.py:_serialize_results 已算), 直接复用,
    不重复计算; 否则现场调 honest_calibration.annotate。
    """
    if not m:
        return {}
    if m.get("compound_confidence") or m.get("honesty_tier"):
        # 已序列化过 → 复用 (字段见 server.py:_serialize_results)
        return {
            "compound_confidence": m.get("confidence") or m.get("compound_confidence"),
            "class_confidence": m.get("confidence_class") or m.get("class_confidence"),
            "exact_structure_confidence": m.get(
                "exact_structure_confidence"
            ) or m.get("confidence") or m.get("compound_confidence"),
            "two_dimensional_connectivity_confidence": m.get(
                "two_dimensional_connectivity_confidence"
            ) or m.get("exact_structure_confidence") or m.get(
                "confidence"
            ) or m.get("compound_confidence"),
            "isomer_group_confidence": m.get(
                "isomer_group_confidence", "Unknown"
            ),
            "identification_scope": m.get("identification_scope"),
            "isomer_group_formula": m.get("isomer_group_formula"),
            "isomer_group_size": m.get("isomer_group_size"),
            "isomer_group_evidence": m.get("isomer_group_evidence") or {},
            "isomer_ambiguity_type": m.get(
                "isomer_ambiguity_type", "none"
            ),
            "isomer_relationship_counts": m.get(
                "isomer_relationship_counts"
            ) or {},
            "isomer_variant_subtype_counts": m.get(
                "isomer_variant_subtype_counts"
            ) or {},
            "isomer_resolution_status": m.get(
                "isomer_resolution_status", "not_applicable"
            ),
            "isomer_literature_evidence": m.get(
                "isomer_literature_evidence"
            ) or (
                (m.get("isomer_group_evidence") or {}).get(
                    "literature_evidence"
                ) or {}
            ),
            "honesty_tier": m.get("honesty_tier"),
            "honest_headline": m.get("honest_headline"),
            "isomer_ambiguous": m.get("isomer_ambiguous"),
            "isomer_n": m.get("isomer_n"),
            "isomer_set": m.get("isomer_set"),
            "class_consensus": m.get("class_consensus"),
            "evidence_tier": m.get("evidence_tier"),
            "annotation_level": m.get("annotation_level"),
            "rt_confirmed": m.get("rt_confirmed"),
            "rt_warning": m.get("rt_warning"),
            "deconvolution_score_influence": m.get(
                "deconvolution_score_influence"
            ),
            "deconvolution_evidence_state": m.get(
                "deconvolution_evidence_state"
            ),
            "chimeric_isolation_warning": m.get(
                "chimeric_isolation_warning"
            ),
        }
    try:
        from matching.honest_calibration import annotate as _annotate
        return _annotate(m, m.get("top_candidates"))
    except Exception:
        return {"compound_confidence": m.get("confidence", "Unknown"),
                "class_confidence": m.get("confidence_class", "Unknown")}


def compound_confidence(m: Optional[dict]) -> str:
    return str((honest_annotation(m) or {}).get("compound_confidence") or "Unknown")


def class_confidence(m: Optional[dict]) -> str:
    return str((honest_annotation(m) or {}).get("class_confidence") or "Unknown")


def isomer_group_confidence(m: Optional[dict]) -> str:
    return str(
        (honest_annotation(m) or {}).get("isomer_group_confidence")
        or "Unknown"
    )


def two_dimensional_connectivity_confidence(m: Optional[dict]) -> str:
    """Confidence at one concrete 2-D molecular connectivity.

    This is deliberately separate from full-identity confidence: reference
    spectra can support one InChIKey connectivity block while leaving
    stereochemistry unresolved.
    """
    ann = honest_annotation(m) or {}
    return str(
        ann.get("two_dimensional_connectivity_confidence")
        or ann.get("exact_structure_confidence")
        or ann.get("compound_confidence")
        or "Unknown"
    )


def exact_structure_confidence(m: Optional[dict]) -> str:
    """Confidence for the complete molecular identity, including stereo."""
    ann = honest_annotation(m) or {}
    return str(
        ann.get("exact_structure_confidence")
        or ann.get("compound_confidence")
        or "Unknown"
    )


def identification_scope(m: Optional[dict]) -> str:
    return str((honest_annotation(m) or {}).get("identification_scope") or "")


def reportable_identity_confidence(m: Optional[dict]) -> str:
    """Return confidence at the identity scope the result actually reports.

    A row scoped to one 2-D connectivity must not be displayed as Medium just
    because the full stereochemical identity remains Medium.  Likewise, an
    unresolved isomer group uses its group confidence.  Lower-scope class and
    putative leads continue to use the complete-candidate confidence and are
    never promoted by this presentation helper.
    """
    scope = identification_scope(m)
    if scope == "isomer_group":
        return isomer_group_confidence(m)
    if scope in {"exact_2d_connectivity", "exact_structure"}:
        return two_dimensional_connectivity_confidence(m)
    return compound_confidence(m)


def reportable_identity_scope_label(
    m: Optional[dict], language: str = "zh"
) -> str:
    scope = identification_scope(m)
    english = str(language).lower().startswith("en")
    labels = {
        "exact_structure": ("full identity", "完整身份"),
        "exact_2d_connectivity": ("2-D connectivity", "二维连接"),
        "isomer_group": ("isomer group", "异构体组"),
        "compound_class": ("chemical class", "化学类别"),
        "literature_target_lead": ("literature lead", "文献线索"),
        "putative_candidate": ("library candidate", "谱库候选"),
        "spectral_lead": ("spectral lead", "谱图线索"),
    }
    return labels.get(scope, ("candidate", "候选"))[0 if english else 1]


def reportable_identity_label(m: Optional[dict], language: str = "zh") -> str:
    confidence = reportable_identity_confidence(m)
    scope_label = reportable_identity_scope_label(m, language=language)
    if str(language).lower().startswith("en"):
        return f"{confidence} ({scope_label})"
    return f"{confidence}（{scope_label}）"


def evidence_tier(m: Optional[dict]) -> str:
    return str((honest_annotation(m) or {}).get("evidence_tier") or "Putative")


def annotation_level(m: Optional[dict]) -> str:
    return str((honest_annotation(m) or {}).get("annotation_level") or "Putative")


def confidence_badge_css(label: str) -> str:
    return _BADGE_CSS.get(str(label), "badge-unknown")


def confidence_badge_html(m: Optional[dict]) -> str:
    """Render confidence at the result's explicitly reportable identity scope.

    The analysis artifact retains the full calibration and ambiguity audit.
    Normal reports deliberately show one compound name plus one confidence
    label so internal evidence scopes do not replace the identification result.
    """
    if not m:
        return '<span class="badge badge-unknown">—</span>'
    comp = reportable_identity_confidence(m)
    label = reportable_identity_label(m, language="zh")
    full_identity = exact_structure_confidence(m)
    title = _html_escape(
        f"可报告身份置信度；完整身份置信度：{full_identity}", quote=True
    )
    return (
        f'<span class="badge {confidence_badge_css(comp)}" '
        f'title="{title}">{_html_escape(label)}</span>'
    )


def final_claim_decision(m: Optional[dict]) -> dict:
    """Return the immutable EviGate 2 projection attached by the API layer."""
    if not isinstance(m, dict):
        return {}
    decision = m.get("final_claim_decision")
    return decision if isinstance(decision, dict) else {}


def final_claim_qvalue(m: Optional[dict]):
    """Return only the q-value calibrated for the decision's declared scope.

    A q-value is deliberately kept separate from probability/confidence fields.
    Missing or out-of-domain calibration therefore exports as blank, never zero.
    """
    decision = final_claim_decision(m)
    calibration = decision.get("calibration") or {}
    if not isinstance(calibration, dict) or calibration.get("status") != "calibrated":
        return None
    field = {
        "authentic_standard": "q_connectivity",
        "2d_connectivity": "q_connectivity",
        "isomer_group": "q_isomer_group",
        "molecular_formula": "q_formula",
        "chemical_class": "q_class",
    }.get(str(decision.get("scope") or ""))
    return calibration.get(field) if field else None


def final_claim_badge_html(m: Optional[dict], language: str = "zh") -> str:
    """Render the identification conclusion separately from evidence grade."""
    decision = final_claim_decision(m)
    english = str(language).lower().startswith("en")
    if not decision:
        label = "Not assessed" if english else "尚未评估"
        title = (
            "Identification-conclusion conditions have not been assessed"
            if english else "尚未评估鉴定结论条件"
        )
        return (
            '<span style="display:inline-block;padding:2px 7px;border-radius:9px;'
            'background:#f1f5f9;color:#64748b;font-size:.68rem;font-weight:700" '
            f'title="{title}">{label}</span>'
        )
    released = bool(decision.get("released"))
    label = (
        "Released" if english else "可作结论"
    ) if released else (
        "Criteria not met" if english else "未达到当前结论条件"
    )
    scope = str(decision.get("scope") or "unresolved")
    calibration = decision.get("calibration") or {}
    cal_status = str(calibration.get("status") or "unavailable")
    qvalue = final_claim_qvalue(m)
    scope_labels = {
        "authentic_standard": ("same-method standard", "同方法标准品确认"),
        "2d_connectivity": ("2-D connectivity candidate", "具体化合物候选（连接结构）"),
        "isomer_group": ("molecular-formula isomer group", "分子式异构体组"),
        "molecular_formula": ("molecular formula", "分子式"),
        "chemical_class": ("chemical class", "化学类别"),
        "unresolved": ("library candidate", "谱库候选"),
    }
    scope_label = scope_labels.get(scope, ("library candidate", "谱库候选"))[
        0 if english else 1
    ]
    calibration_label = (
        "available" if cal_status == "calibrated" else "not yet available"
    ) if english else (
        "已建立" if cal_status == "calibrated" else "尚未建立"
    )
    q_text = (
        f"; q={qvalue:.4g}" if isinstance(qvalue, (int, float)) else ""
    )
    title = (
        f"Identification conclusion: {label}; reportable level: {scope_label}; "
        f"project calibration: {calibration_label}{q_text}"
        if english
        else f"鉴定结论：{label}；当前结论层级：{scope_label}；"
        f"项目校准：{calibration_label}{q_text}"
    )
    background = "#dcfce7" if released else "#fff3d6"
    color = "#166534" if released else "#8a5707"
    return (
        f'<span style="display:inline-block;padding:2px 7px;border-radius:9px;'
        f'background:{background};color:{color};font-size:.68rem;font-weight:700" '
        f'title="{title}">{label}</span>'
    )


# ══════════════════════════════════════════════════════════════════
# 共享 RDKit SVG (前体/碎片 2D 结构)
# ══════════════════════════════════════════════════════════════════
def svg_for_smiles(smi: str, size=(118, 88)) -> str:
    """SMILES → 内联 SVG (无 xml 头, 自适应宽); 失败→''。"""
    if not smi:
        return ""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Draw
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            return ""
        AllChem.Compute2DCoords(mol)
        d = Draw.MolDraw2DSVG(size[0], size[1])
        d.DrawMolecule(mol)
        d.FinishDrawing()
        s = d.GetDrawingText().split("?>", 1)[-1]
        return s.replace("<svg ", '<svg style="width:100%;height:auto" ', 1)
    except Exception:
        return ""


def _observed_mzs(query_peaks) -> list:
    if not query_peaks:
        return []
    out = []
    for p in query_peaks:
        try:
            out.append(float(p[0]))
        except (TypeError, ValueError, IndexError):
            continue
    return out


def _matched_in_query(mz: float, observed: list, tol: float = 0.05) -> bool:
    return any(abs(mz - q) <= tol for q in observed)


# ══════════════════════════════════════════════════════════════════
# 预存裂解途径渲染器 (问题2 核心) — 直接消费 DB 的 fragmentation_pathway JSON
# ══════════════════════════════════════════════════════════════════
def _parse_pathway_json(pathway_json) -> Optional[dict]:
    """惰性解析预存途径 JSON (字符串或已解析 dict 均可)。"""
    if not pathway_json:
        return None
    if isinstance(pathway_json, dict):
        return pathway_json
    try:
        obj = json.loads(pathway_json)
        return obj if isinstance(obj, dict) else None
    except (ValueError, TypeError):
        return None


_NON_CLASS_LABELS = {"", "unknown", "other", "unclassified", "none", "nan"}


def classification_evidence_status(m: Optional[dict]) -> dict:
    """Return the display class together with its independent-evidence status.

    ``display_class`` can be a spectrum-model fallback.  Its mere presence must
    therefore never be interpreted as a reliable chemical classification.
    """
    if not isinstance(m, dict):
        return {
            "label": "Unknown",
            "source": "missing",
            "reliable": False,
            "status": "unresolved",
            "reason": "missing_match",
        }
    details = m.get("expert_details") or {}
    provenance = details.get("classification_provenance") or {}
    display = str(m.get("display_class") or "").strip()
    raw = str(m.get("classification") or "").strip()
    analysis_class = str(provenance.get("analysis_class") or "").strip()
    candidates = (
        (display, "display"),
        (analysis_class, "analysis"),
        (raw, "raw"),
    )
    label, origin = next(
        (
            (value, candidate_origin)
            for value, candidate_origin in candidates
            if value.lower() not in _NON_CLASS_LABELS
        ),
        ("Unknown", "none"),
    )
    if origin == "none":
        return {
            "label": label,
            "source": str(
                m.get("display_classification_source")
                or m.get("classification_source")
                or "missing"
            ),
            "reliable": False,
            "status": "unresolved",
            "reason": str(
                provenance.get("reliability_reason") or "no_nontrivial_class"
            ),
        }
    source = (
        m.get("display_classification_source")
        if origin == "display"
        else m.get("classification_source")
    )
    if not source:
        source = m.get("classification_source")
    from matching.classification_provenance import classification_source_trust

    source_reliable, source_reason = classification_source_trust(source)
    explicitly_untrusted = source_reason.startswith("untrusted_source:")
    provenance_matches = (
        bool(provenance.get("reliable"))
        and analysis_class.lower() == label.lower()
    )
    explicit_reliable = m.get("classification_reliable") is True
    reliable = bool(
        source_reliable
        or provenance_matches
        or (explicit_reliable and not explicitly_untrusted)
    )
    if reliable:
        reason = (
            str(provenance.get("reliability_reason"))
            if provenance_matches and provenance.get("reliability_reason")
            else source_reason if source_reliable else "explicit_reliable_flag"
        )
    else:
        reason = str(
            provenance.get("reliability_reason")
            or source_reason
            or "independent_classification_evidence_unavailable"
        )
    return {
        "label": label,
        "source": str(source or "missing"),
        "reliable": reliable,
        "status": "reliable" if reliable else "tentative",
        "reason": reason,
    }


def fragmentation_evidence_status(m: Optional[dict]) -> dict:
    """Describe fragmentation evidence without implying a proven mechanism."""
    if not isinstance(m, dict):
        return {
            "evidence_type": "none",
            "has_library_signature": False,
            "has_empirical_neutral_loss_sequence": False,
            "candidate_pathway_eligible_for_fusion": False,
            "candidate_pathway_runtime_influence": False,
            "mechanism_confirmed": False,
        }
    details = m.get("expert_details") or {}
    neutral_loss = details.get("neutral_loss_analysis") or {}
    has_library = bool(_parse_pathway_json(m.get("pathway_json")))
    has_empirical = bool(
        neutral_loss.get("best_pathway")
        or neutral_loss.get("sequential_losses")
        or neutral_loss.get("neutral_losses")
    )
    pathway = details.get("pathway_align") or m.get("pathway_align") or {}
    eligible = pathway.get("eligible_for_fusion") is True
    score_contribution = pathway.get("score_contribution")
    try:
        nonzero_contribution = float(score_contribution or 0.0) != 0.0
    except (TypeError, ValueError):
        nonzero_contribution = False
    runtime_influence = bool(
        pathway.get("runtime_influence")
        or pathway.get("ranking_changed")
        or nonzero_contribution
    )
    evidence_types = []
    if has_library:
        evidence_types.append("library_derived_fragmentation_signature")
    if has_empirical:
        evidence_types.append("empirical_neutral_loss_sequence")
    if not evidence_types:
        evidence_types.append("none")
    return {
        "evidence_type": ";".join(evidence_types),
        "has_library_signature": has_library,
        "has_empirical_neutral_loss_sequence": has_empirical,
        "candidate_pathway_eligible_for_fusion": eligible,
        "candidate_pathway_runtime_influence": runtime_influence,
        "mechanism_confirmed": False,
    }


def render_db_pathway(
    pathway_json,
    query_peaks=None,
    max_items: Optional[int] = 8,
    *,
    candidate_smiles: str = "",
    classification: str = "",
    adduct: str = "",
    polarity: str = "",
) -> str:
    """渲染化合物**馆藏预存**裂解途径 (来自该化合物自身多张库谱挖掘)。

    数据源: Compounds.fragmentation_pathway (scripts/mine_pathways.py 产出)。这是该化合物
    稳定的特征碎裂指纹, 独立于任何单次 query —— 即用户要求的"匹配到化合物后直接调预存途径"。
    schema: {n_spectra, precursor_mz, neutral_losses[], diagnostic_ions[], loss_chains[]}
      - neutral_losses[]: {loss_da, fragment_mz, freq, rel_intensity, annotation, nl_key, recurring}
      - diagnostic_ions[]: {mz, freq, rel_intensity, annotation, diag_key, recurring}
      - loss_chains[]: list[list[float]]  (前体逐次减去的策展中性丢失序列)
    query_peaks 可选 → 在表格里标注哪些特征碎片在本 query 被观测到 (绿色✓)。
    max_items 默认维持紧凑预览；Detailed 报告传 None 时展示全部已保存条目。
    """
    pw = _parse_pathway_json(pathway_json)
    if not pw:
        return ""

    n_spec = pw.get("n_spectra", "?")
    pmz = pw.get("precursor_mz")
    pmz_s = f"{pmz:.4f}" if isinstance(pmz, (int, float)) else "—"
    diags = pw.get("diagnostic_ions") or []
    losses = pw.get("neutral_losses") or []
    chains = pw.get("loss_chains") or []
    if not (diags or losses or chains):
        return ""

    observed = _observed_mzs(query_peaks)
    parts = []

    def _items(values):
        values = list(values or [])
        if max_items is None:
            return values
        return values[:max(0, int(max_items))]

    def _recur_badge(it):
        if it.get("recurring"):
            return '<span style="font-size:8px;color:#276749;background:#C6F6D5;padding:0 4px;border-radius:3px">复现</span>'
        if it.get("nl_key") or it.get("diag_key"):
            return '<span style="font-size:8px;color:#2D3748;background:#EDF2F7;padding:0 4px;border-radius:3px">策展</span>'
        return ""

    def _obs_badge(mz):
        if mz and _matched_in_query(float(mz), observed):
            return '<span style="color:#38A169;font-weight:700">✓本谱</span>'
        return ""

    def _loss_annotation(it):
        annotation = str(it.get("annotation") or "").strip()
        key = str(it.get("nl_key") or "").strip()
        observed_loss = _diagnostic_float(it.get("loss_da"))
        if not annotation or not key or observed_loss is None:
            return annotation or "—", "mass_difference_only"
        try:
            from matching.fragmentation_rules import get_neutral_loss_by_mass
            candidates = get_neutral_loss_by_mass(
                observed_loss,
                tolerance=0.01,
                classification=classification,
                strict_fg=bool(classification),
                polarity=polarity,
                candidate_smiles=candidate_smiles,
                query_adduct=adduct,
            )
            supported = any(
                str(item.formula or "").upper() == key.upper()
                or str(item.name_en or "").lower() == annotation.lower()
                for item in candidates
            )
        except Exception:
            supported = False
        if supported:
            return annotation, "current_policy_supported"
        return f"未归属（旧标签 {annotation} 已驳回）", "legacy_annotation_rejected"

    # 诊断离子表
    if diags:
        rows = ""
        for ion in _items(diags):
            mz = ion.get("mz")
            try:
                mz_s = f"{float(mz):.3f}" if mz is not None else "—"
            except (TypeError, ValueError):
                mz_s = "—"
            ann = ion.get("annotation") or "—"
            freq = ion.get("freq")
            freq_s = f"{float(freq)*100:.0f}%" if isinstance(freq, (int, float)) else "—"
            ri = ion.get("rel_intensity")
            ri_s = f"{float(ri):.0f}" if isinstance(ri, (int, float)) else "—"
            rows += (f'<tr><td class="mono">{mz_s}</td><td>{ann}</td><td>{freq_s}</td>'
                     f'<td>{ri_s}</td><td>{_recur_badge(ion)}</td><td>{_obs_badge(mz)}</td></tr>')
        parts.append(
            '<div class="pathway-section"><b>🔬 参考谱复现碎片离子</b> <span class="cascade-sub">候选相关统计，不自动等同于类特异诊断离子</span>'
            '<table class="fragment-table"><thead><tr>'
            '<th>m/z</th><th>名称</th><th>出现率</th><th>相对强度</th><th>来源</th><th>本谱</th>'
            f'</tr></thead><tbody>{rows}</tbody></table></div>'
        )

    # 中性丢失表
    if losses:
        rows = ""
        for it in _items(losses):
            ld = it.get("loss_da")
            fmz = it.get("fragment_mz")
            try:
                ld_s = f"{float(ld):.3f}" if ld is not None else "—"
            except (TypeError, ValueError):
                ld_s = "—"
            try:
                fmz_s = f"{float(fmz):.4f}" if fmz is not None else "—"
            except (TypeError, ValueError):
                fmz_s = "—"
            ann, annotation_status = _loss_annotation(it)
            freq = it.get("freq")
            freq_s = f"{float(freq)*100:.0f}%" if isinstance(freq, (int, float)) else "—"
            status_badge = (
                '<span style="color:#b91c1c;font-weight:700">旧归属已驳回</span>'
                if annotation_status == "legacy_annotation_rejected" else _recur_badge(it)
            )
            rows += (f'<tr><td class="mono">{ld_s}</td><td class="mono">{fmz_s}</td><td>{_html_escape(ann)}</td>'
                     f'<td>{freq_s}</td><td>{status_badge}</td><td>{_obs_badge(fmz)}</td></tr>')
        parts.append(
            '<div class="pathway-section"><b>📉 中性丢失</b> <span class="cascade-sub">前体→产物</span>'
            '<table class="fragment-table"><thead><tr>'
            '<th>Δm</th><th>产物 m/z</th><th>归属</th><th>出现率</th><th>来源</th><th>本谱</th>'
            f'</tr></thead><tbody>{rows}</tbody></table></div>'
        )

    # Stored chains are m/z node sequences, not a list of losses.  The old
    # renderer subtracted those node masses again and produced impossible
    # cascades.  Show the observed/reference topology and adjacent deltas.
    if chains:
        chain_cells = ""
        shown = 0
        for ch in _items(chains):
            if not isinstance(ch, list) or not ch:
                continue
            node_values = []
            for value in ch:
                try:
                    node_values.append(float(value))
                except (TypeError, ValueError):
                    continue
            if not node_values:
                continue
            steps = f'<span class="mono">{node_values[0]:.4f}</span>'
            for parent_mz, child_mz in zip(node_values, node_values[1:]):
                steps += (
                    f' <span style="color:#9C4221">→ Δ{parent_mz-child_mz:.4f}</span> '
                    f'<span class="mono">{child_mz:.4f}</span>'
                )
            chain_cells += f'<div style="min-width:120px">{steps}</div>'
            shown += 1
        if shown:
            parts.append(
                '<div class="pathway-section"><b>🔗 参考谱复现 m/z 链</b> '
                '<span class="cascade-sub">相邻节点质量差为候选解释，不是已确认机理</span>'
                f'<div style="display:flex;gap:8px;flex-wrap:wrap;margin-top:4px">{chain_cells}</div></div>'
            )

    return (
        '<div class="pathway-card">'
        f'<div class="pathway-title">🧬 候选相关馆藏碎片摘要 · 来自 {n_spec} 张库谱 · 前体 m/z {pmz_s}'
        '<span style="font-size:9px;color:#A0AEC0;font-weight:400;margin-left:6px">'
        '(候选依赖的参考谱统计；不作为独立身份分数或已确认机理)</span></div>'
        + "".join(parts) + '</div>'
    )


_DIAGNOSTIC_REASON_LABELS_ZH = {
    "classification_unconfirmed": "尚无可靠类别归属，类特异诊断规则未启用",
    "no_curated_rule_for_ontology_class": (
        "候选已有可靠类别归属，但当前发布的规则库尚无该类别的特异诊断离子规则"
    ),
    "no_scope_guarded_diagnostic_ion_match": "诊断规则已评价；当前样品谱未观察到符合规则适用条件的诊断离子",
    "matched_scope_guarded_diagnostic_ion": "当前样品谱观察到符合候选极性、类别和结构条件的诊断离子",
    "base_8d_not_evaluated": "该历史结果未保存基础规则证据维度评价",
}


def _diagnostic_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number and abs(number) != float("inf") else None


def build_diagnostic_ion_audit(
    match: Optional[dict],
    query_peaks=None,
    tolerance_da: float = 0.05,
) -> dict:
    """Build an interpretation-safe diagnostic-ion audit for reports.

    Two evidence sources stay separate:
      1. runtime class/scope-guarded rule hits already saved in expert_details;
      2. candidate-specific recurring ions saved in fragmentation_pathway and
         checked against the report's saved query peak list.

    The second comparison is report-only and never changes score or ranking.
    """
    match = match if isinstance(match, dict) else {}
    details = match.get("expert_details") or {}
    diagnostic = details.get("diagnostic_ion_analysis") or {}
    assessment = (
        (details.get("dimension_assessments") or {}).get("diagnostic_ion")
        or {}
    )
    runtime_hits = [
        dict(item) for item in (diagnostic.get("matched_ions") or [])
        if isinstance(item, dict)
    ]
    runtime_status = str(assessment.get("status") or "").strip()
    if not runtime_status:
        if runtime_hits:
            runtime_status = "scored"
        elif details:
            runtime_status = "unavailable"
        else:
            runtime_status = "not_recorded"
    runtime_reason = str(assessment.get("reason") or "").strip()

    peaks = []
    for peak in query_peaks or []:
        try:
            mz = _diagnostic_float(peak[0])
            intensity = _diagnostic_float(peak[1])
        except (TypeError, IndexError, KeyError):
            continue
        if mz is None:
            continue
        peaks.append((mz, max(0.0, intensity or 0.0)))
    max_intensity = max((intensity for _, intensity in peaks), default=0.0)

    pathway = _parse_pathway_json(match.get("pathway_json")) or {}
    pathway_rows = []
    for ion in pathway.get("diagnostic_ions") or []:
        if not isinstance(ion, dict):
            continue
        expected = _diagnostic_float(ion.get("mz"))
        nearest = None
        if expected is not None and peaks:
            nearest = min(peaks, key=lambda peak: abs(peak[0] - expected))
        delta = (
            abs(nearest[0] - expected)
            if expected is not None and nearest is not None else None
        )
        observed = delta is not None and delta <= float(tolerance_da)
        pathway_rows.append({
            "expected_mz": expected,
            "observed": observed,
            "nearest_mz": nearest[0] if nearest is not None else None,
            "delta_da": delta,
            "query_intensity": nearest[1] if nearest is not None else None,
            "query_relative_intensity_pct": (
                nearest[1] / max_intensity * 100.0
                if nearest is not None and max_intensity > 0 else None
            ),
            "annotation": ion.get("annotation") or "",
            "frequency": _diagnostic_float(ion.get("freq")),
            "library_relative_intensity": _diagnostic_float(
                ion.get("rel_intensity")
            ),
            "source": (
                "recurring" if ion.get("recurring")
                else "curated" if ion.get("diag_key")
                else "stored_library_feature"
            ),
        })

    return {
        "runtime_status": runtime_status,
        "runtime_reason": runtime_reason,
        "runtime_reason_zh": str(
            assessment.get("reason_label_cn")
            or _DIAGNOSTIC_REASON_LABELS_ZH.get(runtime_reason, runtime_reason)
        ),
        "runtime_score": _diagnostic_float(
            assessment.get("score", details.get("diagnostic_ion_score"))
        ) or 0.0,
        "runtime_hits": runtime_hits,
        "pathway_rows": pathway_rows,
        "pathway_observed_count": sum(
            1 for row in pathway_rows if row["observed"]
        ),
        "query_peak_count": len(peaks),
        "tolerance_da": float(tolerance_da),
        "report_review_only": True,
    }


def render_diagnostic_ion_audit(
    match: Optional[dict],
    query_peaks=None,
    tolerance_da: float = 0.05,
) -> str:
    """Render an explicit observed/not-observed/not-applicable diagnostic-ion section."""
    audit = build_diagnostic_ion_audit(
        match, query_peaks=query_peaks, tolerance_da=tolerance_da
    )
    status = audit["runtime_status"]
    status_labels = {
        "scored": "已观察到符合规则的诊断离子",
        "no_evidence": "已评价，未观察到规则诊断离子",
        "not_applicable": "未评价（该候选无适用的诊断规则）",
        "unavailable": "未建立或未保存可执行规则",
        "not_recorded": "历史结果未保存诊断离子状态",
    }
    status_label = status_labels.get(status, status or "未记录")
    status_class = (
        "diag-ok" if status == "scored"
        else "diag-miss" if status == "no_evidence"
        else "diag-na"
    )
    reason = audit.get("runtime_reason_zh") or ""
    parts = [
        '<div class="diagnostic-audit" '
        f'data-diagnostic-status="{_html_escape(status)}">',
        '<div class="diag-summary">'
        '<div><b>🔬 诊断离子核查</b>'
        '<span class="diag-sub">规则诊断离子与参考库常见碎片离子分开列示</span></div>'
        f'<span class="diag-status {status_class}">{_html_escape(status_label)}</span>'
        '</div>',
    ]
    if reason:
        parts.append(f'<div class="diag-note">{_html_escape(reason)}</div>')

    runtime_hits = audit["runtime_hits"]
    if runtime_hits:
        rows = []
        for item in runtime_hits:
            expected = _diagnostic_float(item.get("expected_mz"))
            observed = _diagnostic_float(item.get("observed_mz"))
            delta = _diagnostic_float(item.get("delta"))
            if delta is None and expected is not None and observed is not None:
                delta = abs(observed - expected)
            observed_cell = (
                f'<td class="mono">{observed:.5f}</td>'
                if observed is not None else '<td>—</td>'
            )
            expected_cell = (
                f'<td class="mono">{expected:.5f}</td>'
                if expected is not None else '<td>—</td>'
            )
            delta_cell = (
                f'<td class="mono">{delta:.5f}</td>'
                if delta is not None else '<td>—</td>'
            )
            rows.append(
                '<tr>'
                + observed_cell
                + expected_cell
                + delta_cell
                + f'<td>{_html_escape(str(item.get("ion_name") or "—"))}</td>'
                + f'<td class="mono">{_html_escape(str(item.get("ion_formula") or "—"))}</td>'
                + f'<td>{_html_escape(str(item.get("structural_origin") or "—"))}</td>'
                + f'<td>{_html_escape(str(item.get("rule_id") or item.get("evidence_level") or "—"))}</td>'
                + '</tr>'
            )
        parts.append(
            '<div class="diag-block"><div class="diag-block-title">'
            f'符合候选适用条件且已观察到的规则诊断离子（{len(rows)} 条）</div>'
            '<div class="diag-table-wrap"><table class="diag-table"><thead><tr>'
            '<th>实测 m/z</th><th>理论 m/z</th><th>ΔDa</th><th>诊断离子</th>'
            '<th>分子式</th><th>结构来源/意义</th><th>规则</th>'
            f'</tr></thead><tbody>{"".join(rows)}</tbody></table></div></div>'
        )
    elif status == "not_applicable":
        parts.append(
            '<div class="diag-empty">该候选没有可执行的诊断离子规则，因此本节未作是否存在的判断。</div>'
        )
    elif status == "no_evidence":
        parts.append(
            '<div class="diag-empty">诊断规则已执行，但当前保存结果只记录实际观察到的离子，'
            '没有保存每条未出现规则的完整清单。因此，这里只能说明当前结果没有记录到符合规则的诊断离子，'
            '不能反推所有可能的诊断离子都不存在。</div>'
        )
    else:
        parts.append(
            '<div class="diag-empty">当前任务没有保存可逐项核对的规则诊断离子结果。</div>'
        )

    pathway_rows = audit["pathway_rows"]
    if pathway_rows:
        rows = []
        for row in pathway_rows:
            expected = row["expected_mz"]
            nearest = row["nearest_mz"]
            delta = row["delta_da"]
            freq = row["frequency"]
            rel = row["query_relative_intensity_pct"]
            observed = row["observed"]
            expected_cell = (
                f'<td class="mono">{expected:.5f}</td>'
                if expected is not None else '<td>—</td>'
            )
            nearest_cell = (
                f'<td class="mono">{nearest:.5f}</td>'
                if nearest is not None else '<td>—</td>'
            )
            delta_cell = (
                f'<td class="mono">{delta:.5f}</td>'
                if delta is not None else '<td>—</td>'
            )
            relative_intensity_cell = (
                f'<td>{rel:.1f}%</td>' if rel is not None else '<td>—</td>'
            )
            frequency_cell = (
                f'<td>{freq * 100:.0f}%</td>'
                if freq is not None else '<td>—</td>'
            )
            rows.append(
                '<tr>'
                + expected_cell
                + nearest_cell
                + delta_cell
                + relative_intensity_cell
                + f'<td>{_html_escape(str(row["annotation"] or "—"))}</td>'
                + frequency_cell
                + f'<td>{_html_escape(str(row["source"]))}</td>'
                + (
                    '<td><span class="diag-hit">✓ 已观察</span></td>'
                    if observed else '<td><span class="diag-not-seen">○ 未观察</span></td>'
                )
                + '</tr>'
            )
        parts.append(
            '<div class="diag-block"><div class="diag-block-title">'
            f'参考库记录中该候选的常见碎片离子：样品谱观察到 {audit["pathway_observed_count"]}/{len(pathway_rows)}'
            f' <span class="diag-sub">报告复核容差 ±{audit["tolerance_da"]:.3f} Da</span></div>'
            '<div class="diag-table-wrap"><table class="diag-table"><thead><tr>'
            '<th>候选 m/z</th><th>最近实测 m/z</th><th>ΔDa</th><th>实测相对强度</th>'
            '<th>注释</th><th>馆藏出现率</th><th>来源</th><th>状态</th>'
            f'</tr></thead><tbody>{"".join(rows)}</tbody></table></div></div>'
        )
    else:
        parts.append(
            '<div class="diag-block"><div class="diag-block-title">参考库记录中该候选的常见碎片离子</div>'
            '<div class="diag-empty">参考库没有为该候选提供可逐项核对的常见碎片离子清单。'
            '因此本节无法评价，而不是检查后判定这些离子全部不存在。</div></div>'
        )

    parts.append(
        '<div class="diag-footnote">数据来源说明：规则诊断离子直接读取任务保存的 expert_details；'
        '参考库常见碎片离子使用报告保存的样品谱峰表进行只读核对，不重新计算得分、排序或置信度。'
        '诊断离子是支持或限定候选的证据，不单独等同于唯一结构确证。</div></div>'
    )
    return "".join(parts)


# ══════════════════════════════════════════════════════════════════
# CFM-ID 候选碎片渲染器 — 明确区分模型预测与观测支持
# ══════════════════════════════════════════════════════════════════
def render_cfm_cascade(smiles: str, ion: str = "+", query_peaks=None,
                       allow_docker: bool = True, max_products: int = 6,
                       cache_only: Optional[bool] = None) -> str:
    """Render CFM-ID candidate product structures with evidence boundaries.

    与旧逻辑的关键区别: **CFM 跑出碎片就展示其预测产物**, 不再以"预测质量与观测 m/z <=0.5Da"
    作为是否显示的门; 观测吻合仅作绿色 ✓ 注释。Docker 开启 + 缓存就绪即有 CFM 输出。

    Args:
        smiles: 前体 SMILES
        ion: "+" / "-"
        query_peaks: 本谱观测峰, 用于标注哪些 CFM 产物被观测到 (绿✓); 不传则不标注
        allow_docker: True=缓存 miss 时即时跑 docker (单化合物接口用); False=仅读缓存 (批量/报告用)
        cache_only: 显式覆盖 (None 时 = not allow_docker)
    Returns: HTML 字符串; CFM 无产物时返回 "" (调用方退 RDKit 兜底)。
    """
    if not smiles or not str(smiles).strip():
        return ""
    try:
        from utils.cfm_fragment import get_fragments
    except Exception:
        return ""
    co = (not allow_docker) if cache_only is None else cache_only
    try:
        frags = get_fragments(str(smiles).strip(), depth=1, ionization=ion, cache_only=co)
    except Exception:
        frags = []
    if not frags:
        return ""

    observed = _observed_mzs(query_peaks)

    # ── 按质量去重: fraggraph-gen 对**同一个碎片离子 m/z** 会枚举多个候选结构(同质量、
    # 不同连接性, 动辄 10-26 个)。旧代码直接取 frags[1:7] → 6 条同质量记录 → 级联画成
    # "前体 → 6 个 m/z 完全相同的结构"(即用户反复指出的乱象)。
    # 去重精度 = 显示精度(0.1 Da, 与下方 "m/z %.1f" 标签一致): 既合并同离子多候选,
    # 也合并 0.1 Da 内的同量异构碎片(显示上本就分不开), 保证不出现两个相同 m/z 标签。
    # 每个独特 m/z 只保留首个候选(fraggraph 按合理性排序, 首个最 sane)。
    seen_mass = {}
    for f in frags:
        try:
            m = round(float(f.get("mass", 0)), 1)
        except (TypeError, ValueError):
            continue
        if m not in seen_mass:
            seen_mass[m] = f
    uniq = sorted(seen_mass.values(), key=lambda x: -float(x.get("mass", 0)))

    # ── 真前体: 用输入 SMILES。fraggraph "fragonly" 模式不含前体, 旧代码误把最大碎片
    # 当前体(其质量在前体列表里重复 10-26 次, 本就不是唯一的前体)。质量取 RDKit 精确质量
    # 并按离子模式加/减质子 → 与谱图 [M+H]⁺/[M-H]⁻ 前体 m/z 一致, 避免新的"前体质量对不上"困惑。
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        _mol = Chem.MolFromSmiles(str(smiles).strip())
        _pre_neutral = float(Descriptors.ExactMolWt(_mol)) if _mol else None
    except Exception:
        _pre_neutral = None
    if _pre_neutral is None:
        _pre_neutral = float(uniq[0].get("mass", 0)) if uniq else 0.0
    if ion == "-":
        _pre_mass = _pre_neutral - 1.007276   # [M-H]⁻
    elif ion == "+":
        _pre_mass = _pre_neutral + 1.007276   # [M+H]⁺
    else:
        _pre_mass = _pre_neutral
    precursor = {"smiles": str(smiles).strip(), "mass": _pre_mass}
    pre_svg = svg_for_smiles(precursor["smiles"])
    # 碎片节点 = 去重后的独特碎片离子(排除与前体同 m/z 的, 若有)
    _pre_r = round(_pre_mass, 1)
    products = [f for f in uniq if round(float(f.get("mass", 0)), 1) != _pre_r][:max_products]

    def _cell(frag, is_precursor=False):
        smi = frag.get("smiles", "")
        svg = svg_for_smiles(smi, (115, 86)) if not is_precursor else pre_svg
        if not svg:
            return ""
        try:
            mz = float(frag.get("mass", 0))
        except (TypeError, ValueError):
            return ""
        if is_precursor:
            return (
                '<div style="display:flex;flex-direction:column;align-items:center;min-width:96px">'
                '<div style="font-size:8px;color:#22543D;font-weight:700;background:#F0FFF4;padding:1px 5px;border-radius:3px;margin-bottom:2px">前体</div>'
                f'<div style="width:115px;border:1px solid #9AE6B4;border-radius:4px;background:white;padding:2px">{svg}</div>'
                f'<div style="font-size:9px;color:#22543D;font-weight:600;margin-top:2px">m/z {mz:.1f}</div></div>'
            )
        # 预测产物: 命中观测峰→绿框 + ✓, 否则橙框
        hit = _matched_in_query(mz, observed, tol=0.5)
        border = "#38A169" if hit else "#FBD38D"
        badge = ('<div style="font-size:8px;color:#276749;font-weight:700;background:#C6F6D5;padding:1px 5px;'
                 'border-radius:3px;margin-bottom:2px;white-space:nowrap">✓ 观测吻合</div>') if hit else ""
        return (
            '<div style="display:flex;flex-direction:column;align-items:center;min-width:96px">'
            f'{badge}<div style="width:115px;border:1px solid {border};border-radius:4px;background:white;padding:2px">{svg}</div>'
            f'<div style="font-size:9px;color:#7B341E;font-weight:600;margin-top:2px">m/z {mz:.1f}</div></div>'
        )

    cells = [c for c in (_cell(p) for p in products) if c]
    pre_cell = _cell(precursor, is_precursor=True)
    if not pre_cell:
        return ""
    # 无匹配观测的产物也算 (解耦): 只要有前体+至少一个产物结构就画
    arrow = '<div style="display:flex;align-items:center;color:#A0AEC0;font-size:14px;padding:0 2px">→</div>'
    src_tag = "即时预测 · Docker" if (allow_docker and not co) else "读自预计算缓存"
    note = (f'<div style="font-size:8px;color:#A0AEC0;margin-top:6px;line-height:1.4">'
            f'ⓘ CFM-ID depth-1 模型给出的候选产物结构（{src_tag}）；这不是已确认的气相机理。'
            '同一显示 m/z 的多个结构假设仅保留模型排序首位用于排版，不代表唯一结构；'
            '绿框只表示预测质量附近存在观测峰。</div>')
    body = pre_cell + (arrow + "".join(cells) if cells else "")
    return (f'<div style="margin-top:4px"><div style="display:flex;gap:6px;flex-wrap:wrap;align-items:center">'
            f'{body}</div>{note}</div>')


def render_product_structure_candidates(candidates) -> str:
    """Render observation-linked RDKit product hypotheses from saved evidence."""

    cells = []
    for item in candidates or []:
        if not isinstance(item, dict):
            continue
        fragment_smiles = str(item.get("fragment_smiles") or "").strip()
        svg = svg_for_smiles(fragment_smiles, (200, 130)) if fragment_smiles else ""
        if not svg:
            continue
        observed = _diagnostic_float(item.get("observed_mz"))
        expected = _diagnostic_float(item.get("expected_mz") or item.get("fragment_mass"))
        delta = _diagnostic_float(item.get("delta"))
        observed_text = f"{observed:.5f}" if observed is not None else "—"
        expected_text = f"{expected:.5f}" if expected is not None else "—"
        delta_text = f"<br>Δ {delta:.5f} Da" if delta is not None else ""
        assignment_type = str(item.get("assignment_type") or "").strip()
        cells.append(
            '<div style="width:230px;padding:10px;border:1px solid #93c5fd;'
            'border-radius:7px;background:#fff">'
            f'<div style="height:136px">{svg}</div>'
            f'<div style="font-size:10.5px;line-height:1.5;color:#1e3a8a">'
            f'观测 m/z {observed_text}<br>理论 m/z {expected_text}{delta_text}'
            f'<br>{_html_escape(str(item.get("adduct") or ""))}</div>'
            '<code title="候选产物结构 SMILES" style="display:block;margin-top:5px;'
            'font-size:9.5px;line-height:1.4;color:#334155;overflow-wrap:anywhere;'
            'white-space:normal">'
            f'{_html_escape(fragment_smiles)}</code>'
            + (
                '<div style="margin-top:4px;font-size:9px;color:#64748b;'
                'overflow-wrap:anywhere">'
                f'{_html_escape(assignment_type)}</div>'
                if assignment_type else ""
            )
            + '</div>'
        )
        if len(cells) >= 8:
            break
    if not cells:
        return ""
    return (
        '<div style="margin-top:4px"><div style="font-size:12px;font-weight:700;'
        'color:#1e3a8a;margin-bottom:5px">观测峰支持的候选产物结构</div>'
        '<div style="display:flex;gap:8px;flex-wrap:wrap">'
        + "".join(cells)
        + '</div><div style="font-size:9.5px;color:#64748b;margin-top:7px">'
        '结构由候选前体的单键断裂枚举与精确质量匹配得到；属于候选结构解释，'
        '不等于已确认的碎裂机理或唯一产物结构。</div></div>'
    )


def render_fragmentation_viz(smiles: str, ion: str = "+", losses=None,
                             query_peaks=None, allow_docker: bool = True,
                             cache_only: Optional[bool] = None,
                             product_candidates=None) -> str:
    """碎裂产物结构可视化 (统一 CFM→RDKit 兜底): 取代 app.py/server.py 两处重复逻辑。

    优先 CFM-ID 预测产物级联 (render_cfm_cascade); 无产物→RDKit 质量匹配法断键高亮
    (utils.frag_viz.rdkit_bond_highlight_svg, 无需 Docker); 都没有→"".
    """
    observed_candidates = render_product_structure_candidates(product_candidates)
    if observed_candidates:
        return observed_candidates
    cfm = render_cfm_cascade(smiles, ion=ion, query_peaks=query_peaks,
                             allow_docker=allow_docker, cache_only=cache_only)
    if cfm:
        return cfm
    try:
        from utils.frag_viz import rdkit_bond_highlight_svg
        rsvg, nhl = rdkit_bond_highlight_svg(smiles, losses or [])
    except Exception:
        rsvg, nhl = "", 0
    if rsvg and nhl > 0:
        return (
            '<div style="margin-top:4px">'
            '<div style="font-size:9px;color:#9C4221;font-weight:700;background:#FFF5E6;padding:2px 7px;'
            f'border-radius:3px;display:inline-block;margin-bottom:4px">断键位置推断 · RDKit 质量匹配法 (标注 {nhl} 个可能断裂的键)</div>'
            f'<div style="max-width:400px;border:1px solid #FBD38D;border-radius:6px;background:white;padding:6px">{rsvg}</div>'
            '<div style="font-size:8px;color:#A0AEC0;margin-top:4px;line-height:1.4">'
            'ⓘ 当前仅标出与质量差相容的候选断键位置；不是唯一产物结构，也不构成已确认机理。'
            '</div></div>'
        )
    return ""


# ══════════════════════════════════════════════════════════════════
# 规范列定义 (CSV / Excel 统一) — 列名 & 顺序单一真源
# ══════════════════════════════════════════════════════════════════
def _fmt(v):
    if v is None or v == "":
        return ""
    if isinstance(v, float):
        return f"{v:.4f}" if v != int(v) else f"{v:.1f}"
    return str(v)


def _reason_text(m: dict) -> str:
    """结构推理文本 (与 detailed 报告一致: v3 推理链优先, 退 chemical_reasoning)。"""
    ed = m.get("expert_details") or {}
    v3 = ed.get("v3") or {}
    parts = [str(x)[:55] for x in (v3.get("reasoning") or [])[:4]]
    if not parts:
        cr = ed.get("chemical_reasoning")
        if cr:
            parts = [str(cr)[:120]]
    return "; ".join(parts)


def _rt_evidence_text(m: dict) -> str:
    ev = m.get("rt_evidence") or {}
    if not ev:
        return ""
    status = ev.get("status", "")
    delta = ev.get("rt_delta_sec")
    tol = ev.get("tolerance_sec")
    std = ev.get("standard_name") or ""
    src = ev.get("source") or ""
    bits = [status]
    if std:
        bits.append(str(std))
    if delta is not None:
        bits.append(f"DeltaRT={delta}s")
    if tol is not None:
        bits.append(f"tol={tol}s")
    if src:
        bits.append(str(src))
    return " | ".join(bits)


def _topk_text(m: dict) -> str:
    """Top-K 候选文本 (按化合物去重, 与 server.py 旧 _topk_str 一致)。"""
    seen, out = set(), []
    for c in (m.get("top_candidates") or [])[:6]:
        rn = (c.get("compound_name") or "").lower().split("(")[0].split("[")[0].strip()
        nn = "".join(ch for ch in rn if ch.isalnum())
        k = f"{nn}|{c.get('formula')}"
        if k in seen:
            continue
        seen.add(k)
        try:
            cos = f"{float(c.get('cosine', 0) or 0):.3f}"
        except (TypeError, ValueError):
            cos = "?"
        out.append(f"{(c.get('compound_name') or '')[:25]}({cos})")
    return "; ".join(out)


def isomer_candidate_projection(m: Optional[dict]) -> dict[str, Any]:
    """Return one complete, export-safe unresolved-isomer candidate set.

    The same canonical member builder is used by the compound conclusion page,
    CSV, Excel and report renderers.  The first member is explicitly a ranking
    lead, never a uniquely resolved structure; all later members are competing
    isomer candidates retained by the saved analysis.
    """
    if not isinstance(m, dict) or not m:
        return {
            "status": "not_applicable",
            "candidate_count": 0,
            "competitive_candidate_count": 0,
            "lead_candidate": {},
            "candidate_members": [],
            "competing_candidates": [],
            "missing_candidate_record_count": 0,
            "competing_candidate_text": "",
            "candidate_ranking_text": "",
            "candidate_details_json": "[]",
        }

    annotation = honest_annotation(m) or {}
    enriched = dict(m)
    for key in (
        "identification_scope",
        "isomer_group_formula",
        "isomer_group_size",
        "isomer_group_evidence",
        "isomer_ambiguous",
        "isomer_n",
        "isomer_set",
        "same_formula_alternatives",
        "isomer_ambiguity_type",
        "isomer_resolution_status",
        "isomer_group_confidence",
    ):
        if enriched.get(key) in (None, "", [], {}):
            value = annotation.get(key)
            if value not in (None, "", [], {}):
                enriched[key] = value

    from pipeline.compound_conclusions import ranked_candidate_members

    members = ranked_candidate_members(enriched)
    scope = str(
        enriched.get("identification_scope")
        or annotation.get("identification_scope")
        or ""
    ).lower()
    unresolved = bool(enriched.get("isomer_ambiguous")) or scope in {
        "isomer_group",
        "molecular_formula_isomer_group",
        "stereoisomer_group",
    }
    if not unresolved:
        members = []

    # ``isomer_n`` historically meant "number of alternatives", while a few
    # older annotations populated it from an ``isomer_set`` that also included
    # the lead.  Likewise, some ``candidate_names`` arrays repeated the lead.
    # Correct those two legacy over-counts before using numeric declarations to
    # expose genuinely missing candidate records.
    isomer_set = list(enriched.get("isomer_set") or [])
    lead_name = str(enriched.get("compound_name") or "").strip().casefold()

    def candidate_name(value: Any) -> str:
        if isinstance(value, dict):
            value = value.get("compound_name")
        return str(value or "").strip().casefold()

    isomer_set_names = [candidate_name(value) for value in isomer_set]
    isomer_set_names = [value for value in isomer_set_names if value]
    isomer_set_includes_lead = bool(
        lead_name and lead_name in isomer_set_names
    )
    evidence = enriched.get("isomer_group_evidence") or {}
    evidence_names = [
        str(value or "").strip().casefold()
        for value in (
            evidence.get("candidate_names")
            if isinstance(evidence, dict)
            else []
        ) or []
        if str(value or "").strip()
    ]

    declared_candidate_count = 0
    try:
        legacy_isomer_n = int(enriched.get("isomer_n") or 0)
        isomer_n_total = (
            legacy_isomer_n
            if isomer_set_includes_lead
            else 1 + legacy_isomer_n
        )
        group_size_total = int(enriched.get("isomer_group_size") or 0)
        if evidence_names and group_size_total == len(evidence_names):
            group_size_total = len(set(evidence_names))
        if (
            isomer_set_includes_lead
            and group_size_total == len(isomer_set_names) + 1
        ):
            group_size_total -= 1
        declared_candidate_count = max(
            group_size_total,
            isomer_n_total if legacy_isomer_n else 0,
        )
    except (TypeError, ValueError):
        declared_candidate_count = 0
    candidate_count = max(
        len(members),
        declared_candidate_count if unresolved else 0,
    )
    missing_count = max(0, candidate_count - len(members))
    lead = dict(members[0]) if members else {}
    competitors = [dict(member) for member in members[1:]]

    def number(value: Any) -> str:
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            return ""
        if parsed != parsed or parsed in (float("inf"), float("-inf")):
            return ""
        return f"{parsed:.4f}".rstrip("0").rstrip(".")

    def member_text(member: dict[str, Any]) -> str:
        fields = [
            f"rank={int(member.get('rank') or 0)}",
            f"role={member.get('role') or ''}",
            f"name={member.get('compound_name') or 'Unnamed candidate'}",
        ]
        for label, key in (
            ("formula", "formula"),
            ("adduct", "adduct"),
            ("cosine", "cosine"),
            ("entropy", "entropy"),
            ("final", "final_score"),
            ("rank_score", "rank_score"),
            ("mass_ppm", "mass_ppm"),
            ("inchikey", "inchikey"),
            ("source", "dataset_source"),
            ("library_spectrum", "lib_spectrum_id"),
            ("relationship", "relationship"),
        ):
            raw = member.get(key)
            value = number(raw) if key in {
                "cosine", "entropy", "final_score", "rank_score", "mass_ppm"
            } else str(raw or "").strip()
            if value:
                fields.append(f"{label}={value}")
        return " | ".join(fields)

    ranking_text = "; ".join(member_text(member) for member in members)
    competitor_text = "; ".join(
        member_text(member) for member in competitors
    )
    if missing_count:
        missing_note = (
            f"missing_candidate_records={missing_count}"
        )
        ranking_text = "; ".join(filter(None, [ranking_text, missing_note]))
        competitor_text = "; ".join(
            filter(None, [competitor_text, missing_note])
        )
    details_payload = {
        "status": (
            "unresolved_competition"
            if unresolved and candidate_count > 1
            else "unresolved_candidate_records_incomplete"
            if unresolved
            else "not_applicable"
        ),
        "candidate_count": candidate_count,
        "competitive_candidate_count": max(0, candidate_count - 1),
        "missing_candidate_record_count": missing_count,
        "candidates": members,
    }
    return {
        **details_payload,
        "lead_candidate": lead,
        "candidate_members": members,
        "competing_candidates": competitors,
        "competing_candidate_text": competitor_text,
        "candidate_ranking_text": ranking_text,
        "candidate_details_json": json.dumps(
            details_payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ),
    }


ISOMER_CANDIDATE_EXPORT_COLUMNS = [
    "Result_Spectrum",
    "Precursor_mz",
    "Retention_Time",
    "Identification_Scope",
    "Isomer_Group_Confidence",
    "Isomer_Resolution_Status",
    "Isomer_Ambiguity_Type",
    "Candidate_Count",
    "Candidate_Record_Status",
    "Candidate_Rank",
    "Candidate_Role",
    "Is_Ranking_Lead",
    "Compound_Name",
    "Formula",
    "InChIKey",
    "SMILES",
    "Adduct",
    "Cosine",
    "Entropy",
    "Reverse_Match",
    "Mass_PPM",
    "Final_Score",
    "Rank_Score",
    "Confidence",
    "Classification",
    "Relationship",
    "Dataset_Source",
    "Library_Spectrum",
]


_ISOMER_CANDIDATE_EXPORT_LABELS_ZH = {
    "Result_Spectrum": "代表谱图编号",
    "Precursor_mz": "前体 m/z",
    "Retention_Time": "RT (min)",
    "Identification_Scope": "鉴定层级",
    "Isomer_Group_Confidence": "异构体组置信度",
    "Isomer_Resolution_Status": "异构体解析状态",
    "Isomer_Ambiguity_Type": "异构体歧义类型",
    "Candidate_Count": "组内候选总数",
    "Candidate_Record_Status": "候选记录状态",
    "Candidate_Rank": "候选排名",
    "Candidate_Role": "候选角色",
    "Is_Ranking_Lead": "是否当前最佳候选",
    "Compound_Name": "候选化合物名称",
    "Formula": "分子式",
    "InChIKey": "InChIKey",
    "SMILES": "SMILES",
    "Adduct": "候选谱库加合物",
    "Cosine": "MS/MS 相似度",
    "Entropy": "谱熵相似度",
    "Reverse_Match": "反向匹配覆盖",
    "Mass_PPM": "质量误差 (ppm)",
    "Final_Score": "综合匹配分",
    "Rank_Score": "候选排序分",
    "Confidence": "候选置信度",
    "Classification": "化学类别",
    "Relationship": "与最佳候选关系",
    "Dataset_Source": "谱库来源",
    "Library_Spectrum": "馆藏谱图编号",
}


def isomer_candidate_export_labels(language: str = "zh") -> list[str]:
    if language == "en":
        return list(ISOMER_CANDIDATE_EXPORT_COLUMNS)
    return [
        _ISOMER_CANDIDATE_EXPORT_LABELS_ZH[column]
        for column in ISOMER_CANDIDATE_EXPORT_COLUMNS
    ]


def isomer_candidate_export_rows(
    results: list[dict[str, Any]], *, language: str = "zh"
) -> list[dict[str, Any]]:
    """Return a normalised long table: one row per isomer candidate."""
    output: list[dict[str, Any]] = []
    seen_groups: set[str] = set()
    for result in results or []:
        match = _match_of(result)
        if not isinstance(match, dict):
            continue
        projection = isomer_candidate_projection(match)
        if projection["status"] == "not_applicable":
            continue
        members = list(projection.get("candidate_members") or [])
        signature = json.dumps(
            [
                member.get("inchikey")
                or [member.get("compound_name"), member.get("formula")]
                for member in members
            ],
            ensure_ascii=False,
            sort_keys=True,
        )
        group_key = "|".join([
            str(match.get("isomer_group_formula") or match.get("formula") or ""),
            signature,
        ])
        if group_key in seen_groups:
            continue
        seen_groups.add(group_key)
        annotation = honest_annotation(match) or {}
        common = {
            "Result_Spectrum": result.get("spectrum_id", ""),
            "Precursor_mz": _standard_number(result.get("precursor_mz")),
            "Retention_Time": _standard_number(result.get("retention_time")),
            "Identification_Scope": identification_scope(match),
            "Isomer_Group_Confidence": isomer_group_confidence(match),
            "Isomer_Resolution_Status": annotation.get(
                "isomer_resolution_status", ""
            ),
            "Isomer_Ambiguity_Type": annotation.get(
                "isomer_ambiguity_type", ""
            ),
            "Candidate_Count": projection.get("candidate_count", 0),
        }
        for member in members:
            role = str(member.get("role") or "competitive_isomer")
            output.append({
                **common,
                "Candidate_Record_Status": "available",
                "Candidate_Rank": member.get("rank", ""),
                "Candidate_Role": (
                    "当前最佳候选（排名首位，未代表唯一解析）"
                    if language == "zh" and role == "ranking_lead"
                    else "竞争异构体候选"
                    if language == "zh"
                    else "ranking lead (not uniquely resolved)"
                    if role == "ranking_lead"
                    else "competitive isomer candidate"
                ),
                "Is_Ranking_Lead": bool(member.get("is_ranking_lead")),
                "Compound_Name": member.get("compound_name", ""),
                "Formula": member.get("formula", ""),
                "InChIKey": member.get("inchikey", ""),
                "SMILES": member.get("smiles", ""),
                "Adduct": member.get("adduct", ""),
                "Cosine": member.get("cosine", ""),
                "Entropy": member.get("entropy", ""),
                "Reverse_Match": member.get("reverse_match", ""),
                "Mass_PPM": member.get("mass_ppm", ""),
                "Final_Score": member.get("final_score", ""),
                "Rank_Score": member.get("rank_score", ""),
                "Confidence": member.get("confidence", ""),
                "Classification": member.get("classification", ""),
                "Relationship": member.get("relationship", ""),
                "Dataset_Source": member.get("dataset_source", ""),
                "Library_Spectrum": member.get("lib_spectrum_id", ""),
            })
        for missing_index in range(
            int(projection.get("missing_candidate_record_count") or 0)
        ):
            output.append({
                **common,
                "Candidate_Record_Status": (
                    "旧任务未保存该竞争候选明细"
                    if language == "zh"
                    else "candidate detail absent from legacy artifact"
                ),
                "Candidate_Rank": len(members) + missing_index + 1,
                "Candidate_Role": (
                    "竞争异构体候选" if language == "zh"
                    else "competitive isomer candidate"
                ),
                "Is_Ranking_Lead": False,
            })
    return output


def isomer_candidate_export_rows_from_compounds(
    compounds: list[dict[str, Any]], *, language: str = "zh"
) -> list[dict[str, Any]]:
    """Project final compound conclusions to one row per retained candidate.

    Scan-level matches can contain broad Top-K provenance or stale ambiguity
    flags that are intentionally narrowed when the final compound conclusion is
    built.  CSV/Excel/report delivery must therefore prefer the conclusion's
    structured ``candidate_members`` over rebuilding a wider candidate set from
    every raw match.
    """
    output: list[dict[str, Any]] = []
    accepted_scopes = {
        "isomer_group",
        "molecular_formula_isomer_group",
        "stereoisomer_group",
    }
    for compound in compounds or []:
        if not isinstance(compound, dict):
            continue
        conclusion = compound.get("conclusion") or {}
        if not isinstance(conclusion, dict):
            conclusion = {}
        scope = str(compound.get("identity_scope") or "").lower()
        if scope not in accepted_scopes and not bool(
            conclusion.get("isomer_group_only")
        ):
            continue
        representative = compound.get("representative") or {}
        if not isinstance(representative, dict):
            representative = {}
        match = _match_of(representative) or {}
        annotation = honest_annotation(match) or {}
        supplied = [
            dict(member) for member in (
                compound.get("candidate_members") or []
            ) if isinstance(member, dict)
        ]
        if not supplied:
            seen_names: set[str] = set()
            for raw_name in compound.get("candidate_names") or []:
                name = str(raw_name or "").strip()
                key = name.casefold()
                if not name or key in seen_names:
                    continue
                seen_names.add(key)
                supplied.append({
                    "compound_name": name,
                    "formula": compound.get("formula", ""),
                    "candidate_record_status": (
                        "旧任务仅保存候选名称"
                        if language == "zh"
                        else "legacy artifact retained name only"
                    ),
                })
        members: list[dict[str, Any]] = []
        for index, raw_member in enumerate(supplied):
            member = dict(raw_member)
            member["rank"] = index + 1
            member["role"] = (
                "ranking_lead" if index == 0 else "competitive_isomer"
            )
            member["is_ranking_lead"] = index == 0
            members.append(member)
        try:
            declared_count = int(
                compound.get("candidate_member_count") or 0
            )
        except (TypeError, ValueError):
            declared_count = 0
        candidate_count = max(declared_count, len(members))
        common = {
            "Result_Spectrum": representative.get("spectrum_id", ""),
            "Precursor_mz": _standard_number(
                representative.get("precursor_mz")
            ),
            "Retention_Time": _standard_number(
                representative.get("retention_time")
            ),
            "Identification_Scope": scope or "isomer_group",
            "Isomer_Group_Confidence": (
                conclusion.get("isomer_group_confidence")
                or match.get("isomer_group_confidence")
                or annotation.get("isomer_group_confidence")
                or ""
            ),
            "Isomer_Resolution_Status": (
                match.get("isomer_resolution_status")
                or annotation.get("isomer_resolution_status")
                or ""
            ),
            "Isomer_Ambiguity_Type": (
                match.get("isomer_ambiguity_type")
                or annotation.get("isomer_ambiguity_type")
                or ""
            ),
            "Candidate_Count": candidate_count,
        }
        for member in members:
            lead = bool(member.get("is_ranking_lead"))
            output.append({
                **common,
                "Candidate_Record_Status": member.get(
                    "candidate_record_status", "available"
                ),
                "Candidate_Rank": member.get("rank", ""),
                "Candidate_Role": (
                    "当前最佳候选（排名首位，未代表唯一解析）"
                    if language == "zh" and lead
                    else "竞争异构体候选"
                    if language == "zh"
                    else "ranking lead (not uniquely resolved)"
                    if lead
                    else "competitive isomer candidate"
                ),
                "Is_Ranking_Lead": lead,
                "Compound_Name": member.get("compound_name", ""),
                "Formula": member.get("formula")
                or compound.get("formula", ""),
                "InChIKey": member.get("inchikey", ""),
                "SMILES": member.get("smiles", ""),
                "Adduct": member.get("adduct", ""),
                "Cosine": member.get("cosine", ""),
                "Entropy": member.get("entropy", ""),
                "Reverse_Match": member.get("reverse_match", ""),
                "Mass_PPM": member.get("mass_ppm", ""),
                "Final_Score": member.get("final_score", ""),
                "Rank_Score": member.get("rank_score", ""),
                "Confidence": member.get("confidence", ""),
                "Classification": member.get("classification", ""),
                "Relationship": member.get("relationship", ""),
                "Dataset_Source": member.get("dataset_source", ""),
                "Library_Spectrum": member.get("lib_spectrum_id", ""),
            })
        for missing_index in range(max(0, candidate_count - len(members))):
            output.append({
                **common,
                "Candidate_Record_Status": (
                    "旧任务未保存该竞争候选明细"
                    if language == "zh"
                    else "candidate detail absent from legacy artifact"
                ),
                "Candidate_Rank": len(members) + missing_index + 1,
                "Candidate_Role": (
                    "竞争异构体候选" if language == "zh"
                    else "competitive isomer candidate"
                ),
                "Is_Ranking_Lead": False,
            })
    return output


def render_isomer_candidate_competition_html(
    results: list[dict[str, Any]], *, language: str = "zh"
) -> str:
    """Render a two-level disclosure appendix for every offline HTML report."""
    rows = isomer_candidate_export_rows(results, language=language)
    if not rows:
        return ""

    grouped: list[list[dict[str, Any]]] = []
    group_index: dict[tuple[Any, ...], int] = {}
    for row in rows:
        key = (
            row.get("Result_Spectrum"),
            row.get("Precursor_mz"),
            row.get("Retention_Time"),
            row.get("Isomer_Ambiguity_Type"),
        )
        if key not in group_index:
            group_index[key] = len(grouped)
            grouped.append([])
        grouped[group_index[key]].append(row)

    english = language == "en"
    title = (
        "Unresolved isomer competition — complete offline candidate list"
        if english
        else "未解析异构体竞争候选——离线报告完整清单"
    )
    intro = (
        "Every group lists the current ranking lead first and then every retained "
        "competitive isomer. The ranking lead is not a uniquely resolved structure."
        if english
        else "每组第 1 名是当前排序最佳候选，随后列出全部保留的竞争异构体；第 1 名不等于已经唯一解析的结构。"
    )
    def cell(value: Any, missing: str = "—") -> str:
        if value in (None, ""):
            return missing
        return _html_escape(str(value))

    def compact_number(value: Any, decimals: int) -> str:
        if value in (None, ""):
            return "—"
        try:
            return f"{float(value):.{decimals}f}"
        except (TypeError, ValueError):
            return cell(value)

    def readable_code(value: Any, translations: dict[str, str]) -> str:
        raw = str(value or "").strip()
        if not raw:
            return "—"
        return translations.get(raw, raw.replace("_", " "))

    def field(
        label: str,
        value: Any,
        *,
        mono: bool = False,
        wide: bool = False,
    ) -> str:
        classes = ["icx-field"]
        if mono:
            classes.append("icx-field-mono")
        if wide:
            classes.append("icx-field-wide")
        return (
            f'<div class="{" ".join(classes)}">'
            f'<span>{_html_escape(label)}</span><b>{cell(value)}</b></div>'
        )

    labels = (
        {
            "groups": "unresolved groups",
            "candidates": "candidate records",
            "expand_group": "Show candidates",
            "collapse_group": "Hide candidates",
            "lead": "ranking lead",
            "group_confidence": "Isomer-group confidence",
            "resolution": "Resolution status",
            "show_evidence": "Show scores and evidence",
            "hide_evidence": "Hide scores and evidence",
            "rank": "Rank",
            "formula": "Formula",
            "status": "Record status",
            "classification": "Chemical class",
            "adduct": "Adduct",
            "confidence": "Candidate confidence",
            "cosine": "MS/MS cosine",
            "entropy": "Spectral entropy similarity",
            "reverse": "Reverse match",
            "final": "Final match score",
            "rank_score": "Candidate ranking score",
            "mass_ppm": "Mass error (ppm)",
            "relationship": "Relationship to ranking lead",
            "source": "Library dataset",
            "spectrum": "Library spectrum ID",
        }
        if english
        else {
            "groups": "个未解析异构体组",
            "candidates": "条候选记录",
            "expand_group": "展开组内候选",
            "collapse_group": "收起组内候选",
            "lead": "当前排序第 1",
            "group_confidence": "异构体组置信度",
            "resolution": "解析状态",
            "show_evidence": "查看分数与证据",
            "hide_evidence": "收起分数与证据",
            "rank": "组内排名",
            "formula": "分子式",
            "status": "候选记录状态",
            "classification": "候选化学类别",
            "adduct": "加合物",
            "confidence": "候选置信度",
            "cosine": "MS/MS 相似度",
            "entropy": "谱熵相似度",
            "reverse": "反向匹配分",
            "final": "综合匹配分",
            "rank_score": "候选排序分",
            "mass_ppm": "质量误差（ppm）",
            "relationship": "与当前第 1 名的关系",
            "source": "谱库数据集",
            "spectrum": "谱库谱图编号",
        }
    )
    resolution_labels = (
        {
            "competitor_features_observed_but_unreplicated": (
                "Competitive-candidate features observed; no independent-"
                "injection replicate"
            ),
            "conflicting_features_observed_but_unreplicated": (
                "Conflicting candidate features observed; no independent-"
                "injection replicate"
            ),
            "diagnostically_uninformative": (
                "Current MS/MS fragments do not distinguish group members"
            ),
            "mixed_isomer_group_requires_multi_axis_resolution": (
                "Multiple isomer relationships; current MS/MS does not "
                "uniquely resolve the structure"
            ),
            "same_connectivity_requires_chiral_or_orthogonal_standard": (
                "Same 2D connectivity; chiral separation or orthogonal "
                "authentic-standard evidence is required"
            ),
            "same_connectivity_structure_definition_specificity_ambiguity": (
                "Same 2D connectivity; library records use different levels "
                "of structure specificity"
            ),
            "same_connectivity_structure_definition_unresolved": (
                "Same 2D connectivity; current evidence does not determine "
                "the precise structure definition"
            ),
            "selected_candidate_supported_but_unreplicated": (
                "Current score selects the ranking lead; no independent-"
                "injection replicate"
            ),
            "single_selected_feature_insufficient": (
                "Only one discriminating feature was observed; the structure "
                "remains unresolved"
            ),
        }
        if english
        else {
            "competitor_features_observed_but_unreplicated": (
                "观察到竞争候选的判别特征，尚无独立进样重复"
            ),
            "conflicting_features_observed_but_unreplicated": (
                "观察到相互冲突的候选判别特征，尚无独立进样重复"
            ),
            "diagnostically_uninformative": (
                "当前 MS/MS 碎片不能区分组内异构体"
            ),
            "mixed_isomer_group_requires_multi_axis_resolution": (
                "组内含多种异构关系；当前 MS/MS 不能唯一解析结构"
            ),
            "same_connectivity_requires_chiral_or_orthogonal_standard": (
                "二维连接相同；需手性色谱或正交标准品证据区分"
            ),
            "same_connectivity_structure_definition_specificity_ambiguity": (
                "二维连接相同；谱库记录采用不同的结构定义精细程度"
            ),
            "same_connectivity_structure_definition_unresolved": (
                "二维连接相同；当前证据不能确定具体结构定义"
            ),
            "selected_candidate_supported_but_unreplicated": (
                "当前评分选择第 1 名，尚无独立进样重复"
            ),
            "single_selected_feature_insufficient": (
                "仅观察到一项判别特征，不能唯一解析结构"
            ),
        }
    )
    relationship_labels = (
        {
            "ranking_lead": "Ranking lead (isomer not uniquely resolved)",
            "constitutional_or_positional_isomer": (
                "Constitutional or positional isomer candidate"
            ),
            "same_connectivity_variant": "Same-2D-connectivity isomer candidate",
            "same_formula_competitor": (
                "Same-formula competitor; precise isomer relationship undetermined"
            ),
        }
        if english
        else {
            "ranking_lead": "当前排序第 1（异构体未唯一解析）",
            "constitutional_or_positional_isomer": "构造或位置异构体候选",
            "same_connectivity_variant": "二维连接相同的异构候选",
            "same_formula_competitor": (
                "同分子式竞争候选（具体异构关系未确定）"
            ),
        }
    )
    record_status_labels = (
        {"available": "Candidate detail saved in the task result"}
        if english
        else {"available": "候选明细已保存在任务结果中"}
    )
    label_separator = ": " if english else "："

    blocks = []
    for index, group in enumerate(grouped, 1):
        first = group[0]
        lead = next(
            (row for row in group if bool(row.get("Is_Ranking_Lead"))),
            first,
        )
        lead_name = cell(
            lead.get("Compound_Name"),
            "Unnamed candidate" if english else "未命名候选",
        )
        formula = cell(lead.get("Formula"))
        count = int(first.get("Candidate_Count") or len(group))
        locator = (
            f"MS2 {cell(first.get('Result_Spectrum'))}; "
            f"m/z {compact_number(first.get('Precursor_mz'), 5)}; "
            f"RT {compact_number(first.get('Retention_Time'), 3)} min"
        )
        resolution = readable_code(
            first.get("Isomer_Resolution_Status"), resolution_labels
        )
        group_number = f"Group {index}" if english else f"第 {index} 组"
        candidate_cards = []
        for row in group:
            rank = cell(row.get("Candidate_Rank"))
            name = cell(
                row.get("Compound_Name"),
                "Unnamed candidate" if english else "未命名候选",
            )
            role = cell(row.get("Candidate_Role"))
            candidate_formula = cell(row.get("Formula"))
            candidate_cards.append(
                '<details class="icx-candidate">'
                '<summary class="icx-candidate-summary">'
                f'<span class="icx-rank">{rank}</span>'
                '<span class="icx-candidate-identity">'
                f'<b>{name}</b><small>{role} · {candidate_formula}</small></span>'
                '<span class="icx-score-strip">'
                f'<span>MS/MS <b>{cell(row.get("Cosine"))}</b></span>'
                f'<span>{"Final" if english else "综合匹配分"} '
                f'<b>{cell(row.get("Final_Score"))}</b></span>'
                f'<span>ppm <b>{cell(row.get("Mass_PPM"))}</b></span>'
                '</span>'
                '<span class="icx-candidate-action">'
                f'<span class="icx-open-label">{labels["show_evidence"]}</span>'
                f'<span class="icx-close-label">{labels["hide_evidence"]}</span>'
                '<i aria-hidden="true"></i></span>'
                '</summary>'
                '<div class="icx-candidate-body">'
                '<div class="icx-score-grid">'
                + field(labels["cosine"], row.get("Cosine"))
                + field(labels["entropy"], row.get("Entropy"))
                + field(labels["reverse"], row.get("Reverse_Match"))
                + field(labels["final"], row.get("Final_Score"))
                + field(labels["rank_score"], row.get("Rank_Score"))
                + field(labels["mass_ppm"], row.get("Mass_PPM"))
                + '</div><div class="icx-evidence-grid">'
                + field(labels["rank"], row.get("Candidate_Rank"))
                + field(labels["formula"], row.get("Formula"), mono=True)
                + field(labels["adduct"], row.get("Adduct"), mono=True)
                + field(labels["confidence"], row.get("Confidence"))
                + field(labels["classification"], row.get("Classification"))
                + field(
                    labels["relationship"],
                    readable_code(row.get("Relationship"), relationship_labels),
                )
                + field("InChIKey", row.get("InChIKey"), mono=True, wide=True)
                + field("SMILES", row.get("SMILES"), mono=True, wide=True)
                + field(labels["source"], row.get("Dataset_Source"))
                + field(labels["spectrum"], row.get("Library_Spectrum"), mono=True)
                + field(
                    labels["status"],
                    readable_code(
                        row.get("Candidate_Record_Status"), record_status_labels
                    ),
                    wide=True,
                )
                + '</div>'
                + '</div></details>'
            )
        blocks.append(
            '<details class="icx-group">'
            '<summary class="icx-group-summary">'
            f'<span class="icx-group-number">{group_number}</span>'
            '<span class="icx-group-main">'
            f'<b>{lead_name}</b><small>{labels["lead"]} · {formula}</small></span>'
            '<span class="icx-group-meta">'
            f'<span class="icx-count">{count} {labels["candidates"]}</span>'
            f'<span>{locator}</span>'
            f'<span class="icx-confidence">{labels["group_confidence"]}{label_separator}'
            f'{cell(first.get("Isomer_Group_Confidence"))}</span>'
            f'<span class="icx-resolution">{labels["resolution"]}{label_separator}'
            f'{cell(resolution)}</span></span>'
            '<span class="icx-group-action">'
            f'<span class="icx-open-label">{labels["expand_group"]}</span>'
            f'<span class="icx-close-label">{labels["collapse_group"]}</span>'
            '<i aria-hidden="true"></i></span>'
            '</summary>'
            '<div class="icx-group-body">' + "".join(candidate_cards) + '</div>'
            '</details>'
        )
    return (
        '<style>'
        '.isomer-candidate-export{margin:26px 0;padding:20px;background:#f8fafc;'
        'border:1px solid #d8e2ea;border-radius:14px;color:#172033}'
        '.icx-heading{display:flex;align-items:flex-start;justify-content:space-between;'
        'gap:16px;margin-bottom:15px}'
        '.icx-heading h2{margin:0 0 6px;color:#123c69;font-size:21px;line-height:1.35}'
        '.icx-heading p{margin:0;max-width:920px;color:#536273;font-size:13px;line-height:1.65}'
        '.icx-total{flex:0 0 auto;display:flex;gap:7px;flex-wrap:wrap;justify-content:flex-end}'
        '.icx-total span{padding:5px 9px;border:1px solid #bfd5df;border-radius:999px;'
        'background:#edf6f8;color:#175c69;font-size:11px;font-weight:700;white-space:nowrap}'
        '.icx-group{margin-top:10px;border:1px solid #d6e0e8;border-radius:10px;'
        'background:#fff;overflow:hidden;box-shadow:0 2px 5px rgba(16,42,67,.04)}'
        '.icx-group-summary,.icx-candidate-summary{list-style:none;cursor:pointer}'
        '.icx-group-summary::-webkit-details-marker,.icx-candidate-summary::-webkit-details-marker{display:none}'
        '.icx-group-summary{display:grid;grid-template-columns:auto minmax(180px,1.1fr) '
        'minmax(280px,1.7fr) auto;align-items:center;gap:13px;padding:14px 15px;'
        'background:linear-gradient(90deg,#f4f9fb,#fff);transition:background .15s ease}'
        '.icx-group-summary:hover{background:#edf6f8}'
        '.icx-group-number{padding:5px 8px;border-radius:6px;background:#123c69;color:#fff;'
        'font-size:11px;font-weight:800;white-space:nowrap}'
        '.icx-group-main{min-width:0}.icx-group-main>b{display:block;color:#123c69;'
        'font-size:15px;overflow-wrap:anywhere}.icx-group-main>small{display:block;margin-top:2px;'
        'color:#667085;font-size:11px}'
        '.icx-group-meta{display:flex;align-items:center;gap:7px;flex-wrap:wrap;color:#59697a;'
        'font-size:11px}.icx-group-meta>span{padding:3px 7px;border-radius:5px;background:#f3f6f8}'
        '.icx-group-meta .icx-count{background:#e5f4ef;color:#17624a;font-weight:800}'
        '.icx-group-meta .icx-confidence{background:#e8f2f8;color:#174f70;font-weight:700}'
        '.icx-group-meta .icx-resolution{background:#fff7e8;color:#805a12}'
        '.icx-group-action,.icx-candidate-action{display:inline-flex;align-items:center;gap:6px;'
        'color:#176b68;font-size:11px;font-weight:800;white-space:nowrap}'
        '.icx-group-action i,.icx-candidate-action i{width:8px;height:8px;border-right:2px solid currentColor;'
        'border-bottom:2px solid currentColor;transform:rotate(45deg);transition:transform .18s ease}'
        '.icx-group[open] .icx-group-action i,.icx-candidate[open] .icx-candidate-action i{'
        'transform:rotate(225deg)}'
        '.icx-close-label{display:none}.icx-group[open]>.icx-group-summary .icx-open-label,'
        '.icx-candidate[open]>.icx-candidate-summary .icx-open-label{display:none}'
        '.icx-group[open]>.icx-group-summary .icx-close-label,'
        '.icx-candidate[open]>.icx-candidate-summary .icx-close-label{display:inline}'
        '.icx-group-body{padding:12px;background:#f5f8fa;border-top:1px solid #dde6ec}'
        '.icx-candidate{margin:8px 0;border:1px solid #d9e3ea;border-radius:8px;background:#fff;'
        'overflow:hidden}.icx-candidate:first-child{margin-top:0}.icx-candidate:last-child{margin-bottom:0}'
        '.icx-candidate-summary{display:grid;grid-template-columns:36px minmax(170px,1.25fr) '
        'minmax(260px,1fr) auto;align-items:center;gap:11px;padding:10px 12px}'
        '.icx-candidate-summary:hover{background:#f7fbfc}'
        '.icx-rank{display:flex;width:30px;height:30px;align-items:center;justify-content:center;'
        'border-radius:50%;background:#e8f2f6;color:#124e6b;font-size:12px;font-weight:900}'
        '.icx-candidate-identity{min-width:0}.icx-candidate-identity>b{display:block;color:#203a4b;'
        'font-size:13px;overflow-wrap:anywhere}.icx-candidate-identity>small{display:block;margin-top:2px;'
        'color:#718096;font-size:10.5px}'
        '.icx-score-strip{display:flex;gap:6px;flex-wrap:wrap}.icx-score-strip>span{padding:3px 7px;'
        'border-radius:5px;background:#f0f5f8;color:#607080;font-size:10px}.icx-score-strip b{color:#154e6b}'
        '.icx-candidate-body{padding:12px 14px 14px;border-top:1px solid #e2e8ee;background:#fbfcfd}'
        '.icx-score-grid{display:grid;grid-template-columns:repeat(6,minmax(90px,1fr));gap:7px;'
        'margin-bottom:9px}.icx-score-grid .icx-field{background:#eef5f7;border-color:#d5e5e9}'
        '.icx-score-grid .icx-field b{font-size:15px;color:#12546a}'
        '.icx-evidence-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:7px}'
        '.icx-field{min-width:0;padding:8px 9px;border:1px solid #e1e7ec;border-radius:6px;background:#fff}'
        '.icx-field>span{display:block;margin-bottom:2px;color:#778493;font-size:9.5px;'
        'letter-spacing:.02em}.icx-field>b{display:block;color:#263f4f;font-size:11px;font-weight:700;'
        'overflow-wrap:anywhere}.icx-field-mono>b{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;'
        'word-break:break-all}.icx-field-wide{grid-column:span 2}'
        '@media(max-width:900px){.icx-heading{display:block}.icx-total{justify-content:flex-start;margin-top:10px}'
        '.icx-group-summary{grid-template-columns:auto 1fr auto}.icx-group-meta{grid-column:1/-1}'
        '.icx-candidate-summary{grid-template-columns:36px 1fr auto}.icx-score-strip{grid-column:2/-1}'
        '.icx-score-grid{grid-template-columns:repeat(3,1fr)}.icx-evidence-grid{grid-template-columns:repeat(2,1fr)}}'
        '@media(max-width:560px){.isomer-candidate-export{padding:12px}.icx-heading h2{font-size:18px}'
        '.icx-group-summary{grid-template-columns:1fr;gap:8px}.icx-group-number{justify-self:start}'
        '.icx-group-action{justify-self:start}.icx-group-meta{grid-column:auto}'
        '.icx-candidate-summary{grid-template-columns:32px 1fr}.icx-candidate-action{grid-column:2}'
        '.icx-score-strip{grid-column:1/-1}.icx-score-grid,.icx-evidence-grid{grid-template-columns:1fr 1fr}'
        '.icx-field-wide{grid-column:1/-1}}'
        '@media print{.isomer-candidate-export{padding:0;border:0;background:#fff}'
        '.icx-group:not([open])>.icx-group-body,.icx-candidate:not([open])>.icx-candidate-body{display:block!important}'
        '.icx-group-action,.icx-candidate-action{display:none!important}.icx-group,.icx-candidate{box-shadow:none}'
        '.icx-group{break-inside:auto}.icx-candidate{break-inside:avoid}.icx-group-body{padding:6px;background:#fff}'
        '.icx-candidate-summary{padding:7px}.icx-candidate-body{padding:7px}.icx-field{padding:5px}}'
        '</style>'
        '<section class="isomer-candidate-export">'
        '<div class="icx-heading"><div>'
        f'<h2>{_html_escape(title)}</h2>'
        f'<p>{_html_escape(intro)}</p></div>'
        '<div class="icx-total">'
        f'<span>{len(grouped)} {labels["groups"]}</span>'
        f'<span>{len(rows)} {labels["candidates"]}</span>'
        '</div></div>'
        + "".join(blocks)
        + "</section>"
    )


def _sample_support_text(m: dict) -> str:
    s = m.get("sample_consensus_support") or {}
    if not s:
        return ""
    try:
        text = (
            f"{int(s.get('n_spectra') or 0)} MS2; "
            f"RT span {float(s.get('rt_span_min') or 0):.3f} min; "
            f"best cos {float(s.get('best_cosine') or 0):.3f}; "
            f"best final {float(s.get('best_final') or 0):.3f}"
        )
        energies = s.get("collision_energies") or []
        if energies:
            text += "; CE " + "/".join(str(x) for x in energies) + " eV"
        return text
    except Exception:
        return str(s)


def _isotope_evidence_text(m: dict, r: dict) -> str:
    ev = m.get("formula_isotope_evidence") or {}
    if ev:
        return (
            f"{ev.get('status', '')}; score={_fmt(ev.get('score'))}; "
            f"{ev.get('interpretation', '')}"
        )
    status = r.get("isotope_evidence_status") or "not_available"
    if status == "referenced_ms1_not_in_file":
        return "unavailable: referenced MS1 scan is not present in the input file"
    return str(status)


def _isolation_evidence_text(r: dict) -> str:
    evidence = r.get("precursor_isolation_evidence") or {}
    status = evidence.get("status") or "not_available"
    purity = evidence.get("purity")
    if purity is None:
        return str(status)
    return (
        f"{status}; purity={_fmt(purity)}; "
        f"coisolated_peaks={evidence.get('coisolated_peak_count', '')}; "
        f"{evidence.get('interpretation', '')}"
    )


def _contrastive_shadow_text(m: dict) -> str:
    evidence = m.get("contrastive_expert_shadow") or {}
    if not evidence:
        return ""
    return (
        f"{evidence.get('status', '')}; class={evidence.get('candidate_class', '')}; "
        f"ions={len(evidence.get('matched_ions') or [])}; "
        f"losses={len(evidence.get('matched_losses') or [])}; "
        f"weight={_fmt(evidence.get('descriptive_evidence_weight'))}; "
        "score_contribution=0"
    )


def _isomer_literature(m: dict) -> dict:
    annotation = honest_annotation(m) or {}
    return (
        annotation.get("isomer_literature_evidence")
        or (annotation.get("isomer_group_evidence") or {}).get(
            "literature_evidence"
        )
        or m.get("isomer_literature_evidence")
        or (m.get("isomer_group_evidence") or {}).get(
            "literature_evidence"
        )
        or {}
    )


def _spec_text(r: dict) -> str:
    """Return the processed MS/MS spectrum used for library comparison.

    Vendor-profile spectra can contain thousands of raw points and exceed the
    32,767-character spreadsheet-cell limit.  ``query_peaks`` is the cleaned,
    auditable spectrum actually used by the search.  Raw points remain in the
    detailed JSON/audit data and are only a fallback for legacy rows.
    """
    peaks = r.get("query_peaks") or r.get("raw_peaks") or []
    out = []
    for p in peaks:
        try:
            out.append(f"{float(p[0]):.5f}:{float(p[1]):.0f}")
        except (TypeError, ValueError, IndexError):
            continue
    return " ".join(out)


def _raw_spec_text(r: dict) -> str:
    """Return the saved, pre-processing MS/MS peak list.

    This is intentionally separate from :func:`_spec_text`: the latter is the
    spectrum actually used for library comparison, while this field lets a
    reader audit what was retained from the input before peak filtering.  Very
    large profile spectra are still protected by the CSV producer's explicit
    spreadsheet-cell limit and truncation marker.
    """
    out = []
    for peak in r.get("raw_peaks") or []:
        try:
            out.append(f"{float(peak[0]):.5f}:{float(peak[1]):.0f}")
        except (TypeError, ValueError, IndexError):
            continue
    return " ".join(out)


def _processed_peaks_count(r: dict):
    """Return the number of peaks represented by ``MS_MS_Spectrum``."""
    query = r.get("query_peaks")
    if isinstance(query, (list, tuple)) and query:
        return len(query)
    return r.get("processed_peaks_count", r.get("peaks_count", ""))


def _raw_peaks_count(r: dict):
    value = r.get("raw_peaks_input_count")
    if value not in (None, ""):
        return value
    raw = r.get("raw_peaks")
    if isinstance(raw, (list, tuple)):
        return len(raw)
    return r.get("peaks_count", "")


def _specific_adduct(value) -> str:
    """Return a chemically specific adduct, not a polarity placeholder."""
    text = str(value or "").strip()
    if text.lower() in {
        "", "nan", "none", "null", "n/a", "unknown",
        "pos", "positive", "+", "neg", "negative", "-",
    }:
        return ""
    return text


def _resolved_adduct(m: dict, r: dict) -> tuple[str, str, object]:
    """Return an auditable adduct value for presentation exports.

    Sample-side ion identity wins over candidate metadata.  A library adduct
    is retained as an explicitly candidate-conditioned hypothesis when no
    independent sample assignment exists.  Polarity alone never assigns an
    adduct.
    """
    lcms_adduct = _specific_adduct(r.get("lcms_adduct"))
    if lcms_adduct:
        return lcms_adduct, "lcms_ion_identity", ""

    query_adduct = _specific_adduct(
        m.get("query_adduct") or r.get("query_adduct")
    )
    if query_adduct and bool(
        m.get("query_adduct_reliable")
        or r.get("query_adduct_reliable")
    ):
        return query_adduct, "explicit_query_metadata", ""

    audited_candidate_hypothesis = _specific_adduct(
        m.get("candidate_adduct_hypothesis")
    )
    if audited_candidate_hypothesis:
        return (
            audited_candidate_hypothesis,
            "candidate_structure_formula_precursor_hypothesis",
            "",
        )

    try:
        charge = r.get("precursor_charge")
        if charge not in (None, "") and abs(int(float(charge))) != 1:
            reported = _specific_adduct(m.get("adduct"))
            return (
                (reported, "library_candidate_hypothesis", "") if reported
                else ("Unresolved", "unsupported_precursor_charge", "")
            )
        precursor_mz = float(r.get("precursor_mz"))
    except (TypeError, ValueError):
        reported = _specific_adduct(m.get("adduct"))
        return (
            (reported, "library_candidate_hypothesis", "") if reported
            else ("Unresolved", "missing_precursor_mz", "")
        )

    from matching.rule_quality import (
        formula_exact_mass,
        normalize_polarity,
    )
    from pipeline.adduct_hypotheses import detect_adduct

    neutral_mass = formula_exact_mass(str(m.get("formula") or ""))
    polarity = normalize_polarity(
        r.get("ion_mode") or m.get("query_ion_mode")
    )
    if neutral_mass is None or polarity is None:
        reported = _specific_adduct(m.get("adduct"))
        return (
            (reported, "library_candidate_hypothesis", "") if reported
            else ("Unresolved", "missing_formula_or_polarity", "")
        )

    tolerance_da = max(0.001, neutral_mass * 10e-6)
    detected = detect_adduct(
        precursor_mz,
        tol_da=tolerance_da,
        mode=polarity,
        neutral_mass_reference=neutral_mass,
    )
    if not detected or detected.get("status") != "assigned":
        reported = _specific_adduct(m.get("adduct"))
        return (
            (reported, "library_candidate_hypothesis", "") if reported
            else ("Unresolved", "no_unique_formula_mz_solution", "")
        )
    hypotheses = detected.get("hypotheses") or []
    unique_adducts = {
        str(item.get("adduct") or "").strip()
        for item in hypotheses
        if str(item.get("adduct") or "").strip()
    }
    if len(unique_adducts) != 1:
        reported = _specific_adduct(m.get("adduct"))
        return (
            (reported, "library_candidate_hypothesis", "") if reported
            else ("Unresolved", "no_unique_formula_mz_solution", "")
        )
    error_da = detected.get("mass_error_da")
    error_ppm = ""
    if error_da is not None and neutral_mass > 0:
        error_ppm = round(float(error_da) / neutral_mass * 1e6, 4)
    return (
        str(detected.get("adduct") or ""),
        "candidate_formula_precursor_unique_10ppm",
        error_ppm,
    )


def _adduct_candidate_independent(m: dict, r: dict) -> bool:
    return _resolved_adduct(m, r)[1] in {
        "lcms_ion_identity", "explicit_query_metadata"
    }


def _adduct_sample_library_consistency(m: dict, r: dict) -> str:
    sample = _specific_adduct(_resolved_adduct(m, r)[0])
    library = _specific_adduct(m.get("adduct"))
    if not sample or not library:
        return "not_tested"
    if _resolved_adduct(m, r)[1] in {
        "library_candidate_hypothesis",
        "candidate_structure_formula_precursor_hypothesis",
        "candidate_formula_precursor_unique_10ppm",
    }:
        return "not_independent"
    return (
        "consistent" if re.sub(r"\s+", "", sample).upper()
        == re.sub(r"\s+", "", library).upper()
        else "inconsistent"
    )


def _collision_energy_value(value):
    if isinstance(value, dict):
        return value.get("value", value.get("collision_energy", ""))
    return value


def _collision_energy_unit(value, explicit_unit="") -> str:
    if str(explicit_unit or "").strip():
        return str(explicit_unit).strip()
    if isinstance(value, dict):
        return str(value.get("unit") or "").strip()
    return ""


def _library_collision_energy(m: dict):
    raw = m.get("lib_collision_energy_raw")
    if raw not in (None, ""):
        return raw
    return m.get("lib_collision_energy", "")


def _library_ce_metadata_status(m: dict) -> str:
    value = _library_collision_energy(m)
    if value in (None, ""):
        return ""
    return (
        "reported"
        if str(m.get("lib_collision_energy_unit") or "").strip()
        else "unit_missing"
    )


def _matched_fragment_details(m: dict, r: dict) -> str:
    """Format matched query fragments as ``m/z:intensity`` pairs."""
    raw = m.get("ms_ms_fragments")
    targets = []
    if isinstance(raw, (list, tuple)):
        values = raw
    else:
        values = re.split(r"[;,\s]+", str(raw or "").strip())
    for value in values:
        if isinstance(value, dict):
            value = value.get("query_mz", value.get("mz"))
        try:
            targets.append(float(value))
        except (TypeError, ValueError):
            continue

    peaks = []
    for peak in r.get("query_peaks") or r.get("raw_peaks") or []:
        try:
            peaks.append((float(peak[0]), float(peak[1])))
        except (TypeError, ValueError, IndexError):
            continue

    output = []
    for target in targets:
        nearest = min(
            peaks,
            key=lambda peak: abs(peak[0] - target),
            default=None,
        )
        if nearest is not None and abs(nearest[0] - target) <= 0.02:
            output.append(f"{nearest[0]:.5f}:{nearest[1]:.6g}")
        else:
            output.append(f"{target:.5f}")
    return "; ".join(output)


# (列名, 取值函数(m, r)) — m=match dict, r=result row; 全部 Web CSV/Excel 共用此列名与顺序
RESULT_COLUMNS = [
    ("Spectrum", lambda m, r: r.get("spectrum_id", "")),
    ("Ion_Mode", lambda m, r: r.get("ion_mode", "")),
    ("Precursor_mz", lambda m, r: _fmt(r.get("precursor_mz"))),
    ("Search_Precursor_mz", lambda m, r: _fmt(
        r.get("search_precursor_mz", r.get("precursor_mz"))
    )),
    ("Search_Precursor_Source", lambda m, r: r.get(
        "search_precursor_mz_source", ""
    )),
    ("Selected_to_Search_Delta_PPM", lambda m, r: _fmt(
        r.get("search_precursor_delta_ppm")
    )),
    ("MS1_Feature_Used_For_Search", lambda m, r: bool(
        r.get("lcms_feature_used_for_search")
    )),
    ("Retention_Time", lambda m, r: _fmt(r.get("retention_time"))),
    ("Precursor_Charge", lambda m, r: _fmt(r.get("precursor_charge"))),
    ("Collision_Energy_eV", lambda m, r: _fmt(r.get("collision_energy"))),
    ("Query_Collision_Energy", lambda m, r: _fmt(
        _collision_energy_value(r.get("collision_energy"))
    )),
    ("Query_Collision_Energy_Unit", lambda m, r: _collision_energy_unit(
        r.get("collision_energy"), r.get("collision_energy_unit")
    )),
    ("Query_Instrument_Type", lambda m, r: (
        r.get("instrument_type")
        or m.get("query_instrument_type", "")
    )),
    ("Query_Instrument", lambda m, r: (
        r.get("instrument")
        or m.get("query_instrument", "")
    )),
    ("RT_Method_ID", lambda m, r: (
        r.get("rt_method_id")
        or m.get("query_rt_method_id", "")
    )),
    ("Precursor_Intensity", lambda m, r: _fmt(r.get("precursor_intensity"))),
    ("Peaks_Count", lambda m, r: _processed_peaks_count(r)),
    ("Raw_Peaks_Count", lambda m, r: _raw_peaks_count(r)),
    ("Height", lambda m, r: _fmt(r.get("base_peak"))),
    ("Compound_Name", lambda m, r: m.get("compound_name", "")),
    ("Formula", lambda m, r: m.get("formula", "")),
    ("Isomer_Competition_Status", lambda m, r: (
        isomer_candidate_projection(m)["status"]
    )),
    ("Isomer_Candidate_Count", lambda m, r: (
        isomer_candidate_projection(m)["candidate_count"]
    )),
    ("Isomer_Competitive_Candidate_Count", lambda m, r: (
        isomer_candidate_projection(m)["competitive_candidate_count"]
    )),
    ("Isomer_Best_Candidate", lambda m, r: (
        (isomer_candidate_projection(m).get("lead_candidate") or {}).get(
            "compound_name", ""
        )
    )),
    ("Isomer_Competing_Candidates", lambda m, r: (
        isomer_candidate_projection(m)["competing_candidate_text"]
    )),
    ("Isomer_Candidate_Ranking_Details", lambda m, r: (
        isomer_candidate_projection(m)["candidate_ranking_text"]
    )),
    ("Isomer_Candidate_Details_JSON", lambda m, r: (
        isomer_candidate_projection(m)["candidate_details_json"]
    )),
    ("SMILES", lambda m, r: m.get("smiles", "")),
    ("InChIKey", lambda m, r: m.get("inchikey", "")),
    ("Adduct", lambda m, r: _resolved_adduct(m, r)[0]),
    ("Adduct_Source", lambda m, r: _resolved_adduct(m, r)[1]),
    ("Adduct_Inference_Error_PPM", lambda m, r: _resolved_adduct(m, r)[2]),
    ("Adduct_Candidate_Independent", _adduct_candidate_independent),
    ("Adduct_Sample_Library_Consistency", _adduct_sample_library_consistency),
    ("LCMS_Neutral_Mass", lambda m, r: _fmt(
        r.get("lcms_neutral_mass") or r.get("ion_identity_neutral_mass")
    )),
    ("Library_Reported_Adduct", lambda m, r: m.get("adduct", "")),
    ("Candidate_Adduct_Hypothesis", lambda m, r: m.get(
        "candidate_adduct_hypothesis", ""
    )),
    ("Candidate_Adduct_Hypothesis_Source", lambda m, r: m.get(
        "candidate_adduct_hypothesis_source", ""
    )),
    ("Classification", lambda m, r: m.get("classification", "")),
    ("Display_Classification", lambda m, r: classification_evidence_status(m)[
        "label"
    ]),
    ("Classification_Source", lambda m, r: m.get(
        "classification_source", ""
    )),
    ("Display_Classification_Source", lambda m, r: classification_evidence_status(m)[
        "source"
    ]),
    ("Classification_Reliable", lambda m, r: classification_evidence_status(m)[
        "reliable"
    ]),
    ("Classification_Status", lambda m, r: classification_evidence_status(m)[
        "status"
    ]),
    ("Classification_Reliability_Reason", lambda m, r: classification_evidence_status(m)[
        "reason"
    ]),
    ("Classification_Confidence", lambda m, r: _fmt(
        m.get("classification_confidence")
    )),
    ("Fragmentation_Evidence_Type", lambda m, r: fragmentation_evidence_status(m)[
        "evidence_type"
    ]),
    ("Feature_Fragment_Status", lambda m, r: (
        (m.get("feature_fragment_evidence") or {}).get("status", "")
    )),
    ("Feature_Fragment_Independent_Count", lambda m, r: (
        (m.get("feature_fragment_evidence") or {}).get(
            "informative_independent_group_count", 0
        )
    )),
    ("Feature_Fragment_Present_Count", lambda m, r: (
        (m.get("feature_fragment_evidence") or {}).get(
            "present_independent_fragment_count", 0
        )
    )),
    ("Feature_Fragment_Score_Nats", lambda m, r: _fmt(
        (m.get("feature_fragment_evidence") or {}).get(
            "combined_score_nats",
            (m.get("feature_fragment_evidence") or {}).get(
                "combined_log_likelihood_ratio"
            ),
        )
    )),
    ("Feature_Fragment_Score_Bits", lambda m, r: _fmt(
        (m.get("feature_fragment_evidence") or {}).get("combined_score_bits")
    )),
    ("Feature_Fragment_Score_Decibans", lambda m, r: _fmt(
        (m.get("feature_fragment_evidence") or {}).get(
            "combined_score_decibans"
        )
    )),
    ("Feature_Fragment_Production_Score_Contribution", lambda m, r: _fmt(
        (m.get("feature_fragment_evidence") or {}).get(
            "score_contribution", 0.0
        )
    )),
    ("Feature_Fragment_Score_Standard_JSON", lambda m, r: json.dumps(
        (m.get("feature_fragment_evidence") or {}).get(
            "feature_fragment_score_standard", {}
        ),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )),
    ("Feature_Fragment_Rule_Evidence_JSON", lambda m, r: json.dumps(
        (m.get("feature_fragment_evidence") or {}).get("rule_evidence", []),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )),
    ("Candidate_Pathway_Eligible_For_Fusion", lambda m, r: fragmentation_evidence_status(m)[
        "candidate_pathway_eligible_for_fusion"
    ]),
    ("Candidate_Pathway_Runtime_Influence", lambda m, r: fragmentation_evidence_status(m)[
        "candidate_pathway_runtime_influence"
    ]),
    ("Dataset_Source", lambda m, r: m.get("dataset_source", "")),
    ("Source_Tier", lambda m, r: m.get("source_tier", "unknown")),
    ("Library_Spectrum", lambda m, r: m.get("lib_spectrum_id", "")),
    ("Library_Source_Accession", lambda m, r: m.get(
        "lib_source_accession", ""
    )),
    ("Library_Collision_Energy", lambda m, r: _fmt(
        _library_collision_energy(m)
    )),
    ("Library_Collision_Energy_Unit", lambda m, r: m.get(
        "lib_collision_energy_unit", ""
    )),
    ("Library_CE_Metadata_Status", lambda m, r: _library_ce_metadata_status(m)),
    ("Library_Instrument_Type", lambda m, r: m.get(
        "lib_instrument_type", ""
    )),
    ("Library_Instrument", lambda m, r: m.get("lib_instrument", "")),
    ("Intensity_Origin", lambda m, r: m.get("intensity_origin", "unknown")),
    ("Library_Quality_Flags", lambda m, r: "; ".join(m.get("quality_flags") or [])),
    ("Compound_Role", lambda m, r: m.get("compound_role", "unclassified")),
    ("Potential_Contaminant", lambda m, r: bool(m.get("potential_contaminant"))),
    ("Contaminant_Category", lambda m, r: m.get("contaminant_category", "")),
    ("Natural_Product_Candidate", lambda m, r: bool(
        m.get("natural_product_candidate", False)
    )),
    ("Cosine_Score", lambda m, r: _fmt(m.get("cosine"))),
    ("Mass_PPM", lambda m, r: _fmt(m.get("mass_ppm"))),
    ("Fragment_Score", lambda m, r: _fmt(m.get("fragment_score"))),
    ("Intensity_Score", lambda m, r: _fmt(m.get("intensity_score"))),
    ("Isomer_Fragment_Score", lambda m, r: _fmt(m.get("isomer_fragment_score"))),
    ("Isomer_Fragment_Bonus", lambda m, r: _fmt(m.get("isomer_fragment_bonus"))),
    ("Isomer_Diagnostic_Status", lambda m, r: (
        (m.get("isomer_fragment_detail") or {}).get("diagnostic_status", "")
    )),
    ("Isomer_Diagnostic_Product_Ions", lambda m, r: (
        f"{(m.get('isomer_fragment_detail') or {}).get('unique_supported', 0)}/"
        f"{(m.get('isomer_fragment_detail') or {}).get('unique_total', 0)}"
    )),
    ("Isomer_Diagnostic_Neutral_Losses", lambda m, r: (
        f"{(m.get('isomer_fragment_detail') or {}).get('neutral_loss_supported', 0)}/"
        f"{(m.get('isomer_fragment_detail') or {}).get('neutral_loss_total', 0)}"
    )),
    ("Observed_Isomer_Diagnostic_Product_mz", lambda m, r: "; ".join(
        str(item.get("query_mz"))
        for item in (
            (m.get("isomer_fragment_detail") or {}).get(
                "candidate_specific_product_ions"
            ) or []
        )
        if item.get("observed") and item.get("query_mz") is not None
    )),
    ("Observed_Isomer_Diagnostic_Neutral_Loss_Da", lambda m, r: "; ".join(
        str(item.get("query_loss_da"))
        for item in (
            (m.get("isomer_fragment_detail") or {}).get(
                "candidate_specific_neutral_losses"
            ) or []
        )
        if item.get("observed") and item.get("query_loss_da") is not None
    )),
    ("Isomer_Diagnostic_Metadata_Comparability", lambda m, r: (
        (m.get("isomer_fragment_detail") or {}).get(
            "metadata_comparability", ""
        )
    )),
    ("Isomer_Ambiguity_Type", lambda m, r: (
        (honest_annotation(m) or {}).get(
            "isomer_ambiguity_type", ""
        )
    )),
    ("Isomer_Resolution_Status", lambda m, r: (
        (honest_annotation(m) or {}).get(
            "isomer_resolution_status", ""
        )
    )),
    ("Isomer_Variant_Subtypes", lambda m, r: json.dumps(
        (honest_annotation(m) or {}).get(
            "isomer_variant_subtype_counts", {}
        ),
        ensure_ascii=False,
        sort_keys=True,
    )),
    ("Isomer_Spectrally_Competitive_Candidates", lambda m, r: (
        (m.get("isomer_fragment_detail") or {}).get(
            "spectrally_competitive_candidate_count", 0
        )
    )),
    ("Isomer_Opposing_Features_Supported", lambda m, r: (
        (m.get("isomer_fragment_detail") or {}).get(
            "opposing_supported_features", 0
        )
    )),
    ("Isomer_All_Same_Formula_Opposing_Features", lambda m, r: (
        (m.get("isomer_fragment_detail") or {}).get(
            "all_same_formula_opposing_supported_features", 0
        )
    )),
    ("Isomer_Evidence_Balance", lambda m, r: (
        (m.get("isomer_fragment_detail") or {}).get(
            "evidence_balance", 0
        )
    )),
    ("Isomer_Reference_Replication", lambda m, r: json.dumps(
        (m.get("isomer_fragment_detail") or {}).get(
            "reference_replication", {}
        ),
        ensure_ascii=False,
        sort_keys=True,
    )),
    ("Isomer_Diagnostic_Tolerance_Da", lambda m, r: _fmt(
        (m.get("isomer_fragment_detail") or {}).get(
            "diagnostic_mz_tolerance_da"
        )
    )),
    ("Isomer_Literature_Status", lambda m, r: (
        _isomer_literature(m).get("status", "")
    )),
    ("Isomer_Literature_PMIDs", lambda m, r: "; ".join(
        str(source.get("pmid"))
        for source in _isomer_literature(m).get("sources", [])
        if source.get("pmid")
    )),
    ("Isomer_Literature_DOIs", lambda m, r: "; ".join(
        str(source.get("doi"))
        for source in _isomer_literature(m).get("sources", [])
        if source.get("doi")
    )),
    ("Isomer_Literature_Missing_Experimental_Information", lambda m, r: (
        "; ".join(
            str(value)
            for value in _isomer_literature(m).get(
                "missing_experimental_information", []
            )
        )
    )),
    ("Isomer_Literature_PubMed_Query", lambda m, r: (
        _isomer_literature(m).get("pubmed_query_url", "")
    )),
    ("Isomer_Literature_Used_For_Scoring", lambda m, r: bool(
        _isomer_literature(m).get("used_for_scoring", False)
    )),
    ("Shadow_Direct_Cosine_Hungarian", lambda m, r: _fmt(
        ((m.get("multimetric_shadow") or {}).get(
            "direct_cosine_hungarian"
        ) or {}).get("score")
    )),
    ("Shadow_Modified_Cosine_Hungarian", lambda m, r: _fmt(
        ((m.get("multimetric_shadow") or {}).get(
            "modified_cosine_hungarian"
        ) or {}).get("score")
    )),
    ("Shadow_Strict_Modified_Cosine_Hungarian", lambda m, r: _fmt(
        ((m.get("multimetric_shadow") or {}).get(
            "strict_modified_cosine_hungarian"
        ) or {}).get("score")
    )),
    ("Shadow_Neutral_Loss_Cosine", lambda m, r: _fmt(
        ((m.get("multimetric_shadow") or {}).get(
            "neutral_loss_cosine"
        ) or {}).get("score")
    )),
    ("Expert_Score", lambda m, r: _fmt(m.get("expert_score"))),
    ("Contrastive_Expert_Shadow", lambda m, r: _contrastive_shadow_text(m)),
    ("Final_Score", lambda m, r: _fmt(m.get("final_score"))),
    ("Rank_Score", lambda m, r: _fmt(m.get("rank_score"))),
    ("Raw_Cosine_Score", lambda m, r: _fmt(m.get("raw_cosine"))),
    ("Raw_Final_Score", lambda m, r: _fmt(m.get("raw_final_score"))),
    ("Candidate_Comparison_Source", lambda m, r: m.get(
        "candidate_comparison_source", "measured_spectrum"
    )),
    ("Raw_Candidate_Cosine_Margin", lambda m, r: _fmt(
        m.get("candidate_margin")
    )),
    ("Chimeric_Deconvolution_Status", lambda m, r: (
        (m.get("chimeric_deconvolution") or {}).get("status", "")
    )),
    ("Deconvolution_Target_Fraction", lambda m, r: _fmt(
        (m.get("chimeric_deconvolution") or {}).get("target_fraction")
    )),
    ("Deconvolution_Score_Gain", lambda m, r: _fmt(
        (m.get("chimeric_deconvolution") or {}).get("score_contribution")
    )),
    ("Ion_Identity_Network", lambda m, r: r.get("ion_identity_network_id", "")),
    ("Ion_Identity", lambda m, r: r.get("ion_identity", "")),
    ("Ion_Identity_Neutral_Mass", lambda m, r: _fmt(
        r.get("ion_identity_neutral_mass")
    )),
    ("Ion_Identity_Representative", lambda m, r: bool(
        r.get("ion_identity_is_representative")
    )),
    ("Ion_Identity_Shape_Correlation", lambda m, r: _fmt(
        (r.get("ion_identity_evidence") or {}).get("best_shape_correlation")
    )),
    ("Ion_Identity_Conflict", lambda m, r: bool(
        (m or {}).get("ion_identity_conflict")
        or (m or {}).get("ion_identity_adduct_conflict")
    )),
    ("Ion_Identity_Linked_Annotation", lambda m, r: (
        (r.get("ion_identity_linked_annotation") or {}).get("compound_name", "")
    )),
    ("Search_Pass", lambda m, r: m.get("search_pass", "primary_top15")),
    # ``Confidence`` is the reader-facing confidence at the identity scope
    # actually reported by the row.  The three explicit scope columns below
    # preserve the complete calibration and prevent a 2-D/isomer-group High
    # result from being exported as a misleading Medium.
    ("Confidence", lambda m, r: reportable_identity_confidence(m)),
    ("Reportable_Identity_Scope", lambda m, r: identification_scope(m)),
    ("Reportable_Identity_Label", lambda m, r: reportable_identity_label(m)),
    ("Final_Claim_Status", lambda m, r: final_claim_decision(m).get("status", "")),
    ("Final_Claim_Scope", lambda m, r: final_claim_decision(m).get("scope", "")),
    ("Final_Claim_Released", lambda m, r: (
        final_claim_decision(m).get("released", "")
        if final_claim_decision(m) else ""
    )),
    ("Project_Q_Value", lambda m, r: _fmt(final_claim_qvalue(m))),
    ("Calibration_Status", lambda m, r: (
        (final_claim_decision(m).get("calibration") or {}).get("status", "")
    )),
    ("Independent_Injection_Count", lambda m, r: (
        final_claim_decision(m).get("independent_injection_count", "")
    )),
    ("Final_Claim_Reason_Codes", lambda m, r: "; ".join(
        final_claim_decision(m).get("reason_codes") or []
    )),
    ("EviGate2_Policy_ID", lambda m, r: final_claim_decision(m).get("policy_id", "")),
    ("Exact_Structure_Confidence", lambda m, r: exact_structure_confidence(m)),
    ("Two_Dimensional_Connectivity_Confidence", lambda m, r: (
        two_dimensional_connectivity_confidence(m)
    )),
    ("Isomer_Group_Confidence", lambda m, r: isomer_group_confidence(m)),
    ("Confidence_Class", lambda m, r: class_confidence(m)),
    ("Identification_Scope", lambda m, r: (
        (honest_annotation(m) or {}).get("identification_scope", "")
    )),
    ("Isomer_Group_Formula", lambda m, r: (
        (honest_annotation(m) or {}).get("isomer_group_formula", "")
    )),
    ("Isomer_Group_Size", lambda m, r: (
        (honest_annotation(m) or {}).get("isomer_group_size", 0)
    )),
    ("Evidence_Tier", lambda m, r: evidence_tier(m)),
    ("Annotation_Level", lambda m, r: annotation_level(m)),
    ("Sample_Consensus", lambda m, r: _sample_support_text(m)),
    ("MS1_Isotope_Evidence", lambda m, r: _isotope_evidence_text(m, r)),
    ("Precursor_Isolation_Evidence", lambda m, r: _isolation_evidence_text(r)),
    ("RT_Evidence", lambda m, r: _rt_evidence_text(m)),
    ("Library_RT", lambda m, r: _fmt(m.get("lib_retention_time"))),
    ("MS_MS_Fragments", lambda m, r: m.get("ms_ms_fragments", "")),
    ("Matched_Fragment_Details", lambda m, r: _matched_fragment_details(m, r)),
    ("Fragment_Count", lambda m, r: m.get("fragment_count", "")),
    ("MS_MS_Spectrum", lambda m, r: _spec_text(r)),
    ("Structure_Reasoning", lambda m, r: _reason_text(m)),
    ("TopK_Candidates", lambda m, r: _topk_text(m)),
    ("Blank_Status", lambda m, r: r.get("blank_status", "")),
    ("Blank_Feature_Matched", lambda m, r: bool(r.get("blank_feature_matched"))),
    ("Blank_Excluded", lambda m, r: bool(r.get("blank_excluded"))),
    ("Blank_Flagged", lambda m, r: bool(r.get("blank_flagged"))),
    ("Blank_Clean_View_Excluded", lambda m, r: bool(
        r.get("blank_clean_view_excluded")
    )),
    ("Blank_Carryover_Suspected", lambda m, r: bool(
        r.get("blank_carryover_suspected")
    )),
    ("Blank_Filter_Policy", lambda m, r: r.get("blank_filter_policy", "")),
    ("Blank_Evidence_Role", lambda m, r: r.get("blank_evidence_role", "")),
    ("Blank_Evidence_Category", lambda m, r: r.get(
        "blank_evidence_category", ""
    )),
    ("Blank_Evidence_Source", lambda m, r: r.get("blank_evidence_source", "")),
    ("Sample_Blank_Ratio", lambda m, r: _fmt(r.get("blank_sample_ratio"))),
    ("Blank_MZ_Error_PPM", lambda m, r: _fmt(r.get("blank_mz_error_ppm"))),
    ("Blank_RT_Delta_Min", lambda m, r: _fmt(r.get("blank_rt_delta_minutes"))),
    ("Blank_Quantitation_Mode", lambda m, r: r.get("blank_quantitation_mode", "")),
    ("Blank_Sample_XIC_Area", lambda m, r: _fmt(r.get("blank_sample_xic_area"))),
    ("Blank_XIC_Area", lambda m, r: _fmt(r.get("blank_xic_area"))),
    ("Blank_Replicates_Detected", lambda m, r: r.get("blank_replicates_detected", "")),
    ("Blank_Replicates_Applicable", lambda m, r: r.get("blank_replicates_applicable", "")),
    ("Blank_Detection_Fraction", lambda m, r: _fmt(r.get("blank_detection_fraction"))),
]

# 日常“MS2 鉴定表”保留读者实际复核所需的完整采集、匹配、碎片、
# 分类、谱库来源和异构体竞争候选身份。低层异构体诊断、shadow
# metrics、空白审计和内部布尔状态仍只放在 report 态，避免标准表
# 重新变成百余列的审计转储。
RESULTS_COLUMN_NAMES = [
    "Compound_Name",
    "Formula",
    "Isomer_Competition_Status",
    "Isomer_Candidate_Count",
    "Isomer_Competitive_Candidate_Count",
    "Isomer_Best_Candidate",
    "Isomer_Competing_Candidates",
    "Isomer_Candidate_Ranking_Details",
    "Isomer_Candidate_Details_JSON",
    "Isomer_Group_Confidence",
    "Isomer_Resolution_Status",
    "Adduct",
    "Adduct_Source",
    "Adduct_Inference_Error_PPM",
    "Ion_Mode",
    "Precursor_mz",
    "Search_Precursor_mz",
    "Search_Precursor_Source",
    "Selected_to_Search_Delta_PPM",
    "MS1_Feature_Used_For_Search",
    "Precursor_Charge",
    "Retention_Time",
    "Query_Collision_Energy",
    "Query_Collision_Energy_Unit",
    "Query_Instrument_Type",
    "Query_Instrument",
    "Confidence",
    "Reportable_Identity_Scope",
    "Reportable_Identity_Label",
    "Final_Claim_Status",
    "Final_Claim_Scope",
    "Final_Claim_Released",
    "Project_Q_Value",
    "Calibration_Status",
    "Independent_Injection_Count",
    "Final_Claim_Reason_Codes",
    "EviGate2_Policy_ID",
    "Final_Score",
    "Cosine_Score",
    "Mass_PPM",
    "Fragment_Count",
    "MS_MS_Fragments",
    "Matched_Fragment_Details",
    "Feature_Fragment_Status",
    "Feature_Fragment_Independent_Count",
    "Feature_Fragment_Present_Count",
    "Feature_Fragment_Score_Nats",
    "Feature_Fragment_Score_Bits",
    "Feature_Fragment_Score_Decibans",
    "Feature_Fragment_Production_Score_Contribution",
    "MS_MS_Spectrum",
    "Peaks_Count",
    "Classification",
    "Display_Classification",
    "Classification_Source",
    "Display_Classification_Source",
    "Classification_Reliable",
    "Classification_Status",
    "Classification_Reliability_Reason",
    "Classification_Confidence",
    "Compound_Role",
    "Potential_Contaminant",
    "Contaminant_Category",
    "Natural_Product_Candidate",
    "Fragmentation_Evidence_Type",
    "Candidate_Pathway_Eligible_For_Fusion",
    "Candidate_Pathway_Runtime_Influence",
    "Dataset_Source",
    "Source_Tier",
    "Library_Spectrum",
    "Library_Reported_Adduct",
    "Library_Collision_Energy",
    "Library_Collision_Energy_Unit",
    "Library_CE_Metadata_Status",
    "Library_Instrument_Type",
    "Library_Instrument",
    "InChIKey",
    "SMILES",
    "Spectrum",
]
REPORT_COLUMN_NAMES = [c[0] for c in RESULT_COLUMNS]
_RESULTS_DROP = set(REPORT_COLUMN_NAMES) - set(RESULTS_COLUMN_NAMES)

_STANDARD_NUMERIC_COLUMNS = {
    "Isomer_Candidate_Count",
    "Isomer_Competitive_Candidate_Count",
    "Precursor_mz",
    "Search_Precursor_mz",
    "Selected_to_Search_Delta_PPM",
    "Precursor_Charge",
    "Retention_Time",
    "Query_Collision_Energy",
    "Adduct_Inference_Error_PPM",
    "Cosine_Score",
    "Final_Score",
    "Mass_PPM",
    "Fragment_Count",
    "Peaks_Count",
    "Classification_Confidence",
    "Library_Collision_Energy",
    "Project_Q_Value",
    "Independent_Injection_Count",
    "Feature_Fragment_Independent_Count",
    "Feature_Fragment_Present_Count",
    "Feature_Fragment_Score_Nats",
    "Feature_Fragment_Score_Bits",
    "Feature_Fragment_Score_Decibans",
    "Feature_Fragment_Production_Score_Contribution",
}
_STANDARD_NAME_ISOMER_METADATA_RE = re.compile(
    r"\s*\(\s*known\s+isomers\s*:\s*\d+\s*;\s*"
    r"isobaric\s+peaks\s*:\s*\d+\s*\)\s*$",
    re.IGNORECASE,
)
_STANDARD_NAME_SOURCE_PREFIX_RE = re.compile(
    r"^(?:MassBank|GNPS|ReSpect)\s*:\s*\S+\s+",
    re.IGNORECASE,
)
_STANDARD_NAME_ACQUISITION_SUFFIX_RE = re.compile(
    r"\s*;\s*(?:LC|GC)-(?:ESI|APCI|EI)[^;|]*"
    r";\s*MS(?:2|/MS)?(?:\s*;\s*(?:CE|CID|HCD|NCE)(?:\s*[^;|]*)?)?\s*$",
    re.IGNORECASE,
)


def _standard_compound_name(value, match: Optional[dict] = None) -> str:
    """Keep the standard table focused on one readable compound name.

    Source accessions, alias lists and library-supplied isomer counters remain
    available in the explicit detailed-audit export.
    """
    text = str(value or "").strip()
    text = _STANDARD_NAME_SOURCE_PREFIX_RE.sub("", text)
    text = _STANDARD_NAME_ACQUISITION_SUFFIX_RE.sub("", text).strip()
    text = text.split("|", 1)[0].strip()
    text = _STANDARD_NAME_ISOMER_METADATA_RE.sub("", text).strip()
    if match:
        from matching.compound_names import display_compound_name

        text = display_compound_name(
            text,
            match.get("inchikey", ""),
        )
    return text


def _standard_number(value):
    """Return a real spreadsheet number instead of a formatted numeric string."""
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return value
    return int(parsed) if parsed.is_integer() else parsed


def _match_of(r: dict) -> Optional[dict]:
    """从 result row 取 match dict (兼容 serialized 'match' 与 report-style 'matches')。"""
    m = r.get("match")
    if m:
        return m
    ms = r.get("matches")
    if ms:
        return ms[0]
    return None


def row_to_dict(r: dict, report_mode: bool = True) -> dict:
    """把一条 result row → 规范列 dict (供 CSV/Excel 统一生成)。

    report_mode=True 返回全集列; False 返回精简列(去掉内部分析, 用于 results CSV)。
    match 缺失(未命中)的行只填 spectrum 级字段 + Confidence='未命中'。
    """
    m = _match_of(r)
    out = {}
    isomer_fields = {
        "Isomer_Competition_Status",
        "Isomer_Candidate_Count",
        "Isomer_Competitive_Candidate_Count",
        "Isomer_Best_Candidate",
        "Isomer_Competing_Candidates",
        "Isomer_Candidate_Ranking_Details",
        "Isomer_Candidate_Details_JSON",
    }
    for name, fn in RESULT_COLUMNS:
        if name in isomer_fields:
            out[name] = ""
            continue
        try:
            out[name] = (fn(m, r) if m else "")
        except Exception:
            out[name] = ""
    if m:
        projection = isomer_candidate_projection(m)
        lead = projection.get("lead_candidate") or {}
        out.update({
            "Isomer_Competition_Status": projection["status"],
            "Isomer_Candidate_Count": projection["candidate_count"],
            "Isomer_Competitive_Candidate_Count": projection[
                "competitive_candidate_count"
            ],
            "Isomer_Best_Candidate": lead.get("compound_name", ""),
            "Isomer_Competing_Candidates": projection[
                "competing_candidate_text"
            ],
            "Isomer_Candidate_Ranking_Details": projection[
                "candidate_ranking_text"
            ],
            "Isomer_Candidate_Details_JSON": projection[
                "candidate_details_json"
            ],
        })
    if not m:
        out["Confidence"] = "未命中"
        # IIMN can link an unmatched ion form to a network without turning it
        # into an independent library identification.  Preserve that explicit
        # evidence in exports even though ``match`` correctly remains null.
        out["Ion_Identity_Network"] = r.get("ion_identity_network_id", "")
        out["Ion_Identity"] = r.get("ion_identity", "")
        out["Ion_Identity_Neutral_Mass"] = _fmt(
            r.get("ion_identity_neutral_mass")
        )
        out["Ion_Identity_Representative"] = bool(
            r.get("ion_identity_is_representative")
        )
        out["Ion_Identity_Shape_Correlation"] = _fmt(
            (r.get("ion_identity_evidence") or {}).get("best_shape_correlation")
        )
        out["Ion_Identity_Linked_Annotation"] = (
            (r.get("ion_identity_linked_annotation") or {}).get("compound_name", "")
        )
    if report_mode:
        return out

    standard = {
        name: out.get(name, "")
        for name in RESULTS_COLUMN_NAMES
    }
    standard["Compound_Name"] = _standard_compound_name(
        standard.get("Compound_Name"),
        m,
    )
    for name in _STANDARD_NUMERIC_COLUMNS:
        standard[name] = _standard_number(standard.get(name))
    return standard


MS2_UNIFIED_CANDIDATE_COLUMNS = [
    "Record_Type",
    "Identification_Group_ID",
    "Feature_ID",
    "Candidate_Record_Status",
    "Candidate_Rank",
    "Candidate_Role",
    "Is_Ranking_Lead",
    "Candidate_Compound_Name",
    "Candidate_Formula",
    "Candidate_InChIKey",
    "Candidate_SMILES",
    "Candidate_Library_Adduct",
    "Candidate_Cosine",
    "Candidate_Entropy",
    "Candidate_Reverse_Match",
    "Candidate_Mass_PPM",
    "Candidate_Final_Score",
    "Candidate_Rank_Score",
    "Candidate_Confidence",
    "Candidate_Classification",
    "Candidate_Relationship",
    "Candidate_Dataset_Source",
    "Candidate_Library_Spectrum",
]

MS2_UNIFIED_COLUMN_NAMES = [
    *RESULTS_COLUMN_NAMES[:2],
    *MS2_UNIFIED_CANDIDATE_COLUMNS,
    *RESULTS_COLUMN_NAMES[2:],
]


# Reader-facing MS2 result table.  Candidate rows stay in the same CSV, but
# internal identifiers, JSON payloads, rule codes and duplicate booleans are
# deliberately excluded.  Scientific evidence is not removed: acquisition,
# adduct, similarity, fragment-ion, classification and provenance fields are
# retained in a stable reading order.
MS2_READER_COLUMN_NAMES = [
    "Group_Number",
    "Feature_ID",
    "Candidate_Compound_Name",
    "Candidate_Rank",
    "Candidate_Role",
    "Group_Candidate_Count",
    "Candidate_Formula",
    "Identification_Conclusion",
    "Group_Confidence",
    "Candidate_Confidence",
    "Retention_Time",
    "Precursor_mz",
    "Search_Precursor_mz",
    "Search_Precursor_Source",
    "Ion_Mode",
    "Query_Collision_Energy",
    "Query_Collision_Energy_Unit",
    "Query_Instrument_Type",
    "Adduct",
    "Adduct_Source",
    "Candidate_Library_Adduct",
    "Library_Collision_Energy",
    "Library_Collision_Energy_Unit",
    "Library_Instrument_Type",
    "Candidate_Cosine",
    "Candidate_Entropy",
    "Candidate_Final_Score",
    "Candidate_Mass_PPM",
    "Fragment_Count",
    "MS_MS_Fragments",
    "Matched_Fragment_Details",
    "MS_MS_Spectrum",
    "Peaks_Count",
    "Raw_MS_MS_Spectrum",
    "Raw_Peaks_Count",
    "Candidate_Classification",
    "Display_Classification",
    "Classification_Confidence",
    "Display_Classification_Source",
    "Classification_Status",
    "Classification_Reliability_Reason",
    "Candidate_Dataset_Source",
    "Candidate_Library_Spectrum",
    "Candidate_InChIKey",
    "Candidate_SMILES",
    "Spectrum",
    "Export_Note",
]


_MS2_READER_LABELS_ZH = {
    "Group_Number": "鉴定组",
    "Feature_ID": "离子特征编号",
    "Candidate_Compound_Name": "候选化合物",
    "Candidate_Rank": "组内排名",
    "Candidate_Role": "候选角色",
    "Group_Candidate_Count": "组内候选数",
    "Candidate_Formula": "分子式",
    "Identification_Conclusion": "鉴定结论",
    "Group_Confidence": "鉴定组置信度",
    "Candidate_Confidence": "候选置信度",
    "Retention_Time": "RT (min)",
    "Precursor_mz": "MS2 隔离目标 m/z",
    "Search_Precursor_mz": "检索用实测/回退 m/z",
    "Search_Precursor_Source": "检索质量轴来源",
    "Ion_Mode": "离子模式",
    "Query_Collision_Energy": "样本 CE",
    "Query_Collision_Energy_Unit": "样本 CE 单位",
    "Query_Instrument_Type": "样本仪器类型",
    "Adduct": "样品离子形式",
    "Adduct_Source": "样品离子形式依据",
    "Candidate_Library_Adduct": "候选谱库加合物",
    "Library_Collision_Energy": "代表馆藏 CE",
    "Library_Collision_Energy_Unit": "代表馆藏 CE 单位",
    "Library_Instrument_Type": "代表馆藏仪器类型",
    "Candidate_Cosine": "候选 MS/MS 相似度",
    "Candidate_Entropy": "候选谱熵相似度",
    "Candidate_Final_Score": "候选综合匹配分",
    "Candidate_Mass_PPM": "候选质量误差 (ppm)",
    "Fragment_Count": "代表谱匹配碎片数",
    "MS_MS_Fragments": "代表谱匹配碎片 m/z",
    "Matched_Fragment_Details": "代表谱匹配碎片 (m/z:强度)",
    "Feature_Fragment_Status": "特征碎片计分状态",
    "Feature_Fragment_Independent_Count": "独立特征碎裂证据数",
    "Feature_Fragment_Present_Count": "已观测独立特征碎片数",
    "Feature_Fragment_Score_Nats": "特征碎片证据分 (lnLR, nats)",
    "Feature_Fragment_Score_Bits": "特征碎片证据分 (bits)",
    "Feature_Fragment_Score_Decibans": "特征碎片证据分 (decibans)",
    "Feature_Fragment_Production_Score_Contribution": "特征碎片生产排序贡献",
    "MS_MS_Spectrum": "样本 MS2 碎片离子（实际检索谱，m/z:强度）",
    "Peaks_Count": "实际检索 MS2 峰数",
    "Raw_MS_MS_Spectrum": "原始 MS2 碎片离子（已保存输入峰，m/z:强度）",
    "Raw_Peaks_Count": "原始 MS2 输入峰数",
    "Candidate_Classification": "候选化学类别",
    "Display_Classification": "鉴定组化学类别",
    "Classification_Confidence": "分类置信度",
    "Display_Classification_Source": "分类依据",
    "Classification_Status": "分类证据状态",
    "Classification_Reliability_Reason": "分类可靠性说明",
    "Candidate_Dataset_Source": "谱库来源",
    "Candidate_Library_Spectrum": "馆藏谱图编号",
    "Candidate_InChIKey": "InChIKey",
    "Candidate_SMILES": "SMILES",
    "Spectrum": "样本谱图编号",
    "Export_Note": "备注",
}


_MS2_READER_LABELS_EN = {
    "Group_Number": "Identification group",
    "Feature_ID": "Ion feature ID",
    "Candidate_Compound_Name": "Candidate compound",
    "Candidate_Rank": "Rank within group",
    "Candidate_Role": "Candidate role",
    "Group_Candidate_Count": "Candidates in group",
    "Candidate_Formula": "Formula",
    "Identification_Conclusion": "Identification conclusion",
    "Group_Confidence": "Group confidence",
    "Candidate_Confidence": "Candidate confidence",
    "Retention_Time": "RT (min)",
    "Precursor_mz": "MS2 isolation target m/z",
    "Search_Precursor_mz": "Measured/fallback search m/z",
    "Search_Precursor_Source": "Search mass-axis source",
    "Ion_Mode": "Ion mode",
    "Query_Collision_Energy": "Query collision energy",
    "Query_Collision_Energy_Unit": "Query collision-energy unit",
    "Query_Instrument_Type": "Query instrument type",
    "Adduct": "Sample ion form",
    "Adduct_Source": "Sample ion-form basis",
    "Candidate_Library_Adduct": "Candidate library adduct",
    "Library_Collision_Energy": "Representative library collision energy",
    "Library_Collision_Energy_Unit": "Representative library CE unit",
    "Library_Instrument_Type": "Representative library instrument type",
    "Candidate_Cosine": "Candidate MS/MS cosine",
    "Candidate_Entropy": "Candidate spectral entropy",
    "Candidate_Final_Score": "Candidate final score",
    "Candidate_Mass_PPM": "Candidate mass error (ppm)",
    "Fragment_Count": "Representative matched-fragment count",
    "MS_MS_Fragments": "Representative matched-fragment m/z",
    "Matched_Fragment_Details": "Representative matched fragments (m/z:intensity)",
    "Feature_Fragment_Status": "Feature-fragment scoring status",
    "Feature_Fragment_Independent_Count": "Independent feature-fragment evidence count",
    "Feature_Fragment_Present_Count": "Observed independent feature-fragment count",
    "Feature_Fragment_Score_Nats": "Feature-fragment evidence score (lnLR, nats)",
    "Feature_Fragment_Score_Bits": "Feature-fragment evidence score (bits)",
    "Feature_Fragment_Score_Decibans": "Feature-fragment evidence score (decibans)",
    "Feature_Fragment_Production_Score_Contribution": "Feature-fragment production rank contribution",
    "MS_MS_Spectrum": "Sample MS2 ions used for search (m/z:intensity)",
    "Peaks_Count": "MS2 peaks used for search",
    "Raw_MS_MS_Spectrum": "Saved raw input MS2 ions (m/z:intensity)",
    "Raw_Peaks_Count": "Raw input MS2 peak count",
    "Candidate_Classification": "Candidate chemical class",
    "Display_Classification": "Group chemical class",
    "Classification_Confidence": "Classification confidence",
    "Display_Classification_Source": "Classification basis",
    "Classification_Status": "Classification evidence status",
    "Classification_Reliability_Reason": "Classification reliability note",
    "Candidate_Dataset_Source": "Library source",
    "Candidate_Library_Spectrum": "Library spectrum ID",
    "Candidate_InChIKey": "InChIKey",
    "Candidate_SMILES": "SMILES",
    "Spectrum": "Query spectrum ID",
    "Export_Note": "Note",
}


def ms2_reader_column_labels(language: str = "zh") -> list[str]:
    labels = _MS2_READER_LABELS_EN if language == "en" else _MS2_READER_LABELS_ZH
    return [labels[column] for column in MS2_READER_COLUMN_NAMES]


def _unified_candidate_fields(
    candidate: Mapping[str, Any],
    *,
    group_id: str,
    record_type: str,
) -> dict[str, Any]:
    return {
        "Record_Type": record_type,
        "Identification_Group_ID": group_id,
        "Candidate_Record_Status": candidate.get(
            "Candidate_Record_Status", "available"
        ),
        "Candidate_Rank": candidate.get("Candidate_Rank", ""),
        "Candidate_Role": candidate.get("Candidate_Role", ""),
        "Is_Ranking_Lead": bool(candidate.get("Is_Ranking_Lead")),
        "Candidate_Compound_Name": candidate.get("Compound_Name", ""),
        "Candidate_Formula": candidate.get("Formula", ""),
        "Candidate_InChIKey": candidate.get("InChIKey", ""),
        "Candidate_SMILES": candidate.get("SMILES", ""),
        "Candidate_Library_Adduct": candidate.get("Adduct", ""),
        "Candidate_Cosine": candidate.get("Cosine", ""),
        "Candidate_Entropy": candidate.get("Entropy", ""),
        "Candidate_Reverse_Match": candidate.get("Reverse_Match", ""),
        "Candidate_Mass_PPM": candidate.get("Mass_PPM", ""),
        "Candidate_Final_Score": candidate.get("Final_Score", ""),
        "Candidate_Rank_Score": candidate.get("Rank_Score", ""),
        "Candidate_Confidence": candidate.get("Confidence", ""),
        "Candidate_Classification": candidate.get("Classification", ""),
        "Candidate_Relationship": candidate.get("Relationship", ""),
        "Candidate_Dataset_Source": candidate.get("Dataset_Source", ""),
        "Candidate_Library_Spectrum": candidate.get(
            "Library_Spectrum", ""
        ),
    }


def _single_unified_candidate_fields(
    result: Mapping[str, Any],
    base: Mapping[str, Any],
    *,
    group_id: str,
    language: str,
) -> dict[str, Any]:
    match = _match_of(dict(result)) or {}
    return {
        "Record_Type": "single_identification",
        "Identification_Group_ID": group_id,
        "Candidate_Record_Status": "available",
        "Candidate_Rank": 1,
        "Candidate_Role": (
            "单一报告候选" if language == "zh" else "reported identity"
        ),
        "Is_Ranking_Lead": True,
        "Candidate_Compound_Name": base.get("Compound_Name", ""),
        "Candidate_Formula": base.get("Formula", ""),
        "Candidate_InChIKey": base.get("InChIKey", ""),
        "Candidate_SMILES": base.get("SMILES", ""),
        "Candidate_Library_Adduct": match.get("adduct", ""),
        "Candidate_Cosine": match.get("cosine", ""),
        "Candidate_Entropy": match.get("entropy", ""),
        "Candidate_Reverse_Match": match.get("reverse_match", ""),
        "Candidate_Mass_PPM": match.get("mass_ppm", ""),
        "Candidate_Final_Score": match.get("final_score", ""),
        "Candidate_Rank_Score": match.get(
            "rank_score", match.get("final_score", "")
        ),
        "Candidate_Confidence": match.get("confidence", ""),
        "Candidate_Classification": match.get("classification", ""),
        "Candidate_Relationship": (
            "single_reported_identity"
        ),
        "Candidate_Dataset_Source": match.get("dataset_source", ""),
        "Candidate_Library_Spectrum": match.get("lib_spectrum_id", ""),
    }


def _apply_final_isomer_summary(
    base: dict[str, Any], candidates: list[dict[str, Any]]
) -> None:
    names = [
        str(candidate.get("Compound_Name") or "").strip()
        for candidate in candidates
        if str(candidate.get("Compound_Name") or "").strip()
    ]
    candidate_count = max(
        [len(candidates)]
        + [int(candidate.get("Candidate_Count") or 0) for candidate in candidates]
    )
    details = [{
        "rank": candidate.get("Candidate_Rank"),
        "role": candidate.get("Candidate_Role"),
        "is_ranking_lead": bool(candidate.get("Is_Ranking_Lead")),
        "compound_name": candidate.get("Compound_Name", ""),
        "formula": candidate.get("Formula", ""),
        "inchikey": candidate.get("InChIKey", ""),
        "smiles": candidate.get("SMILES", ""),
        "adduct": candidate.get("Adduct", ""),
        "cosine": candidate.get("Cosine", ""),
        "entropy": candidate.get("Entropy", ""),
        "reverse_match": candidate.get("Reverse_Match", ""),
        "mass_ppm": candidate.get("Mass_PPM", ""),
        "final_score": candidate.get("Final_Score", ""),
        "rank_score": candidate.get("Rank_Score", ""),
        "relationship": candidate.get("Relationship", ""),
        "dataset_source": candidate.get("Dataset_Source", ""),
        "lib_spectrum_id": candidate.get("Library_Spectrum", ""),
    } for candidate in candidates]
    base.update({
        "Isomer_Competition_Status": "unresolved_competition",
        "Isomer_Candidate_Count": candidate_count,
        "Isomer_Competitive_Candidate_Count": max(0, candidate_count - 1),
        "Isomer_Best_Candidate": names[0] if names else "",
        "Isomer_Competing_Candidates": "; ".join(names[1:]),
        "Isomer_Candidate_Ranking_Details": "; ".join(
            f"rank={item['rank']} | name={item['compound_name']}"
            for item in details
        ),
        "Isomer_Candidate_Details_JSON": json.dumps(
            {
                "status": "unresolved_competition",
                "candidate_count": candidate_count,
                "competitive_candidate_count": max(0, candidate_count - 1),
                "candidates": details,
            },
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ),
        "Isomer_Group_Confidence": (
            candidates[0].get("Isomer_Group_Confidence", "")
            if candidates else ""
        ),
        "Isomer_Resolution_Status": (
            candidates[0].get("Isomer_Resolution_Status", "")
            if candidates else ""
        ),
    })


def _result_with_source_spectrum(
    representative: Mapping[str, Any],
    source_by_spectrum: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Restore scan-level peaks omitted from compact compound projections."""
    spectrum_id = str(representative.get("spectrum_id") or "")
    source = source_by_spectrum.get(spectrum_id) or {}
    merged = {**dict(source), **dict(representative)}
    source_match = _match_of(dict(source)) or {}
    representative_match = _match_of(dict(representative)) or {}
    if source_match or representative_match:
        merged["match"] = {**source_match, **representative_match}
    return merged


def _standard_ms2_base(result: Mapping[str, Any]) -> dict[str, Any]:
    prepared = dict(result)
    base = row_to_dict(prepared, report_mode=False)
    base.update({
        # These reader-only fields are intentionally not added to the legacy
        # 94-column machine schema.  They feed the compact delivery table.
        "Raw_MS_MS_Spectrum": _raw_spec_text(prepared),
        "Raw_Peaks_Count": _standard_number(_raw_peaks_count(prepared)),
        "_Raw_Peaks_Saved_Count": prepared.get(
            "raw_peaks_saved_count",
            len(prepared.get("raw_peaks") or []),
        ),
        "_Raw_Peaks_Complete": prepared.get("raw_peaks_complete"),
    })
    return base


def ms2_unified_candidate_rows(
    matched: list[dict[str, Any]],
    *,
    compounds: list[dict[str, Any]] | None = None,
    language: str = "zh",
) -> list[dict[str, Any]]:
    """Return one CSV row per final candidate in the ordinary MS2 table.

    Group-level MS2/query fields remain repeated on every row. Candidate-level
    identity and library evidence live in explicit ``Candidate_*`` columns, so
    unresolved isomers are complete without a second CSV and without implying
    that every candidate is an independently observed chromatographic signal.
    """
    output: list[dict[str, Any]] = []
    if compounds is not None:
        source_by_spectrum = {
            str(row.get("spectrum_id")): row
            for row in matched or []
            if isinstance(row, dict) and str(row.get("spectrum_id") or "")
        }
        for index, compound in enumerate(compounds):
            if not isinstance(compound, dict):
                continue
            representative = compound.get("representative") or {}
            if not isinstance(representative, dict):
                continue
            representative = _result_with_source_spectrum(
                representative, source_by_spectrum
            )
            base = _standard_ms2_base(representative)
            base["Feature_ID"] = compound.get("feature_id", "")
            group_id = str(
                compound.get("compound_id")
                or compound.get("identity_key")
                or f"compound-{index + 1}"
            )
            if not _match_of(representative):
                output.append({
                    **base,
                    "Record_Type": "unannotated_feature",
                    "Identification_Group_ID": group_id,
                    "Candidate_Record_Status": "no_library_candidate",
                    "Candidate_Rank": "",
                    "Candidate_Role": "",
                    "Is_Ranking_Lead": False,
                    "Candidate_Compound_Name": "",
                    "Candidate_Formula": "",
                    "Candidate_InChIKey": "",
                    "Candidate_SMILES": "",
                    "Candidate_Library_Adduct": "",
                    "Candidate_Cosine": "",
                    "Candidate_Entropy": "",
                    "Candidate_Reverse_Match": "",
                    "Candidate_Mass_PPM": "",
                    "Candidate_Final_Score": "",
                    "Candidate_Rank_Score": "",
                    "Candidate_Confidence": "",
                    "Candidate_Classification": "",
                    "Candidate_Relationship": "no_library_candidate",
                    "Candidate_Dataset_Source": "",
                    "Candidate_Library_Spectrum": "",
                })
                continue
            candidates = isomer_candidate_export_rows_from_compounds(
                [compound], language=language
            )
            if candidates:
                _apply_final_isomer_summary(base, candidates)
                for candidate in candidates:
                    output.append({
                        **base,
                        **_unified_candidate_fields(
                            candidate,
                            group_id=group_id,
                            record_type="isomer_candidate",
                        ),
                    })
            else:
                output.append({
                    **base,
                    **_single_unified_candidate_fields(
                        representative,
                        base,
                        group_id=group_id,
                        language=language,
                    ),
                })
        return output

    for index, result in enumerate(matched or []):
        if not isinstance(result, dict) or not _match_of(result):
            continue
        base = _standard_ms2_base(result)
        group_id = f"spectrum:{result.get('spectrum_id') or index + 1}"
        candidates = isomer_candidate_export_rows(
            [result], language=language
        )
        if candidates:
            _apply_final_isomer_summary(base, candidates)
            for candidate in candidates:
                output.append({
                    **base,
                    **_unified_candidate_fields(
                        candidate,
                        group_id=group_id,
                        record_type="isomer_candidate",
                    ),
                })
        else:
            output.append({
                **base,
                **_single_unified_candidate_fields(
                    result,
                    base,
                    group_id=group_id,
                    language=language,
                ),
            })
    return output


def _reader_conclusion_label(
    row: Mapping[str, Any],
    compound: Mapping[str, Any] | None,
    language: str,
) -> str:
    conclusion = (
        compound.get("conclusion")
        if isinstance(compound, Mapping)
        else {}
    ) or {}
    if language != "en":
        return str(
            conclusion.get("label")
            or row.get("Reportable_Identity_Label")
            or row.get("Confidence")
            or ""
        )
    key = str(conclusion.get("key") or "").lower()
    confidence = str(
        conclusion.get("isomer_group_confidence")
        or row.get("Isomer_Group_Confidence")
        or ""
    )
    labels = {
        "confirmed": "Confirmed present",
        "background": "Blank/background associated",
        "high": "High-confidence single 2-D connectivity detected",
        "candidate": "Medium-evidence candidate",
        "insufficient": "Insufficient evidence",
    }
    if key == "isomer":
        return (
            "High-confidence isomer group detected"
            if confidence == "High"
            else "Medium-evidence isomer group"
        )
    return labels.get(
        key,
        str(row.get("Reportable_Identity_Label") or row.get("Confidence") or ""),
    )


def _reader_classification_source(value: Any, language: str) -> str:
    raw = str(value or "").strip()
    if not raw:
        return ""
    key = raw.lower()
    labels = {
        "runtime_structure_rule:stored_class_agrees_with_smiles_structure_rule": (
            "SMILES structure rule agrees with the library class",
            "SMILES 结构规则与谱库类别一致",
        ),
        "runtime_structure_rule:stored_and_structure_classes_compatible": (
            "SMILES structure rule is compatible with the library class",
            "SMILES 结构规则与谱库类别相容",
        ),
        "runtime_structure_rule:trusted_source:curated_ciead_npclassifier": (
            "Structure rule plus curated NPClassifier class",
            "结构规则 + 策展 NPClassifier 类别",
        ),
        "runtime_structure_rule:smiles_structure_rule_overrode_provisional_stored_class": (
            "SMILES structure rule corrected a provisional library class",
            "SMILES 结构规则校正暂定谱库类别",
        ),
        "ml": (
            "MS/MS spectrum classification model",
            "MS/MS 谱图分类模型",
        ),
        "library": ("Library classification", "谱库类别"),
        "missing": ("No reliable classification basis", "未提供可靠分类依据"),
    }
    if key in labels:
        return labels[key][0 if language == "en" else 1]
    return (
        "See detailed audit for the classification source code"
        if language == "en"
        else "其他分类依据（详细代码见审计 CSV）"
    )


def _reader_classification_reason(value: Any, language: str) -> str:
    raw = str(value or "").strip()
    if not raw:
        return ""
    key = raw.lower()
    labels = {
        "stored_class_agrees_with_smiles_structure_rule": (
            "Library class and SMILES structure rule agree",
            "谱库类别与 SMILES 结构规则一致",
        ),
        "stored_and_structure_classes_compatible": (
            "Library class and SMILES structure rule are compatible",
            "谱库类别与 SMILES 结构规则相容",
        ),
        "trusted_source:curated_ciead_npclassifier": (
            "Supported by a curated NPClassifier source",
            "由策展 NPClassifier 来源支持",
        ),
        "smiles_structure_rule_overrode_provisional_stored_class": (
            "SMILES structure rule corrected the provisional library class",
            "SMILES 结构规则已校正暂定谱库类别",
        ),
        "untrusted_source:ml": (
            "Model-only classification; independent structural support is absent",
            "仅有模型分类，缺少独立结构证据",
        ),
        "untrusted_source:library": (
            "Library class label only; independent structural support is absent",
            "仅有谱库类别标签，缺少独立结构证据",
        ),
        "untrusted_source:missing": (
            "No reliable classification source was recorded",
            "未记录可靠的分类来源",
        ),
        "no_nontrivial_class": (
            "No interpretable chemical class was obtained",
            "未得到可解释的化学类别",
        ),
        "independent_classification_evidence_unavailable": (
            "Independent classification evidence is unavailable",
            "缺少独立分类证据",
        ),
        "missing_match": (
            "No library match is available for classification",
            "缺少可用于分类的谱库匹配",
        ),
    }
    if key in labels:
        return labels[key][0 if language == "en" else 1]
    if key.startswith("untrusted_source:"):
        source = key.split(":", 1)[1].replace("_", " ")
        return (
            f"Classification source is not independently trusted ({source})"
            if language == "en"
            else f"分类来源未达到独立可靠性门槛（{source}）"
        )
    return (
        "See detailed audit for the classification reliability code"
        if language == "en"
        else "其他分类可靠性原因（详细代码见审计 CSV）"
    )


def ms2_reader_rows(
    matched: list[dict[str, Any]],
    *,
    compounds: list[dict[str, Any]] | None = None,
    language: str = "zh",
) -> list[dict[str, Any]]:
    """Project the complete candidate table into a readable delivery CSV.

    The row unit remains one retained identity candidate.  Group-level MS2
    evidence is repeated so sorting/filtering never detaches a candidate from
    its spectrum, while internal IDs and machine-only payloads remain in the
    explicit detailed-audit export.
    """
    unified = ms2_unified_candidate_rows(
        matched,
        compounds=compounds,
        language=language,
    )
    compound_by_id: dict[str, Mapping[str, Any]] = {}
    for index, compound in enumerate(compounds or []):
        if not isinstance(compound, Mapping):
            continue
        group_id = str(
            compound.get("compound_id")
            or compound.get("identity_key")
            or f"compound-{index + 1}"
        )
        compound_by_id[group_id] = compound

    rows_by_group: dict[str, list[dict[str, Any]]] = {}
    ordered_group_ids: list[str] = []
    for row in unified:
        group_id = str(row.get("Identification_Group_ID") or "")
        if group_id not in rows_by_group:
            ordered_group_ids.append(group_id)
            rows_by_group[group_id] = []
        rows_by_group[group_id].append(row)
    group_numbers = {
        group_id: f"G{index:03d}"
        for index, group_id in enumerate(ordered_group_ids, 1)
    }

    output: list[dict[str, Any]] = []
    for row in unified:
        group_id = str(row.get("Identification_Group_ID") or "")
        group_rows = rows_by_group.get(group_id) or [row]
        declared_count = row.get("Isomer_Candidate_Count")
        try:
            declared_count = int(declared_count or 0)
        except (TypeError, ValueError):
            declared_count = 0
        candidate_count = max(1, declared_count, len(group_rows))
        lead = bool(row.get("Is_Ranking_Lead"))
        record_type = str(row.get("Record_Type") or "")
        single = record_type == "single_identification"
        unannotated = record_type == "unannotated_feature"
        if single:
            # A resolved final compound may still carry broad search Top-K
            # provenance in the representative match.  That provenance is not
            # the final group membership and must not make a single-candidate
            # row claim that it contains 3/7/etc. candidates.
            candidate_count = 1
        if unannotated:
            role = "未注释特征" if language == "zh" else "unannotated feature"
            candidate_count = 0
        elif language == "en":
            role = (
                "single candidate"
                if single
                else "ranking lead (isomer not uniquely resolved)"
                if lead
                else "competitive isomer candidate"
            )
        else:
            role = (
                "单一候选"
                if single
                else "最佳候选（异构体未唯一解析）"
                if lead
                else "竞争异构体候选"
            )

        notes = []
        record_status = str(row.get("Candidate_Record_Status") or "").strip()
        if record_status and record_status.lower() != "available":
            notes.append(record_status)
        if row.get("_Raw_Peaks_Complete") is False:
            notes.append(
                "saved raw peak list is incomplete"
                if language == "en"
                else "原始峰按保存策略未完整保留"
            )
        if len(str(row.get("Raw_MS_MS_Spectrum") or "")) > 32000:
            notes.append(
                "raw peak list exceeds one spreadsheet cell; the CSV cell is explicitly truncated"
                if language == "en"
                else "原始峰列表超过单个表格单元格容量；该单元格将明确截断"
            )

        compound = compound_by_id.get(group_id)
        conclusion = (
            compound.get("conclusion")
            if isinstance(compound, Mapping)
            else {}
        ) or {}
        group_confidence = (
            conclusion.get("isomer_group_confidence")
            or row.get("Isomer_Group_Confidence")
            or row.get("Confidence")
            or ""
        )
        projected = {
            "Group_Number": group_numbers.get(group_id, ""),
            "Feature_ID": row.get("Feature_ID", ""),
            "Candidate_Compound_Name": row.get("Candidate_Compound_Name", ""),
            "Candidate_Rank": row.get("Candidate_Rank", ""),
            "Candidate_Role": role,
            "Group_Candidate_Count": candidate_count,
            "Candidate_Formula": row.get("Candidate_Formula", ""),
            "Identification_Conclusion": _reader_conclusion_label(
                row, compound, language
            ),
            "Group_Confidence": group_confidence,
            "Candidate_Confidence": row.get("Candidate_Confidence", ""),
            "Retention_Time": row.get("Retention_Time", ""),
            "Precursor_mz": row.get("Precursor_mz", ""),
            "Search_Precursor_mz": row.get("Search_Precursor_mz", ""),
            "Search_Precursor_Source": row.get("Search_Precursor_Source", ""),
            "Ion_Mode": row.get("Ion_Mode", ""),
            "Query_Collision_Energy": row.get("Query_Collision_Energy", ""),
            "Query_Collision_Energy_Unit": row.get(
                "Query_Collision_Energy_Unit", ""
            ),
            "Query_Instrument_Type": row.get("Query_Instrument_Type", ""),
            "Adduct": row.get("Adduct", ""),
            "Adduct_Source": row.get("Adduct_Source", ""),
            "Candidate_Library_Adduct": row.get("Candidate_Library_Adduct", ""),
            "Library_Collision_Energy": row.get("Library_Collision_Energy", ""),
            "Library_Collision_Energy_Unit": row.get(
                "Library_Collision_Energy_Unit", ""
            ),
            "Library_Instrument_Type": row.get("Library_Instrument_Type", ""),
            "Candidate_Cosine": row.get("Candidate_Cosine", ""),
            "Candidate_Entropy": row.get("Candidate_Entropy", ""),
            "Candidate_Final_Score": row.get("Candidate_Final_Score", ""),
            "Candidate_Mass_PPM": row.get("Candidate_Mass_PPM", ""),
            "Fragment_Count": row.get("Fragment_Count", ""),
            "MS_MS_Fragments": row.get("MS_MS_Fragments", ""),
            "Matched_Fragment_Details": row.get("Matched_Fragment_Details", ""),
            "Feature_Fragment_Status": row.get("Feature_Fragment_Status", ""),
            "Feature_Fragment_Independent_Count": row.get(
                "Feature_Fragment_Independent_Count", ""
            ),
            "Feature_Fragment_Present_Count": row.get(
                "Feature_Fragment_Present_Count", ""
            ),
            "Feature_Fragment_Score_Nats": row.get(
                "Feature_Fragment_Score_Nats", ""
            ),
            "Feature_Fragment_Score_Bits": row.get(
                "Feature_Fragment_Score_Bits", ""
            ),
            "Feature_Fragment_Score_Decibans": row.get(
                "Feature_Fragment_Score_Decibans", ""
            ),
            "Feature_Fragment_Production_Score_Contribution": row.get(
                "Feature_Fragment_Production_Score_Contribution", ""
            ),
            "MS_MS_Spectrum": row.get("MS_MS_Spectrum", ""),
            "Peaks_Count": row.get("Peaks_Count", ""),
            "Raw_MS_MS_Spectrum": row.get("Raw_MS_MS_Spectrum", ""),
            "Raw_Peaks_Count": row.get("Raw_Peaks_Count", ""),
            "Candidate_Classification": row.get("Candidate_Classification", ""),
            "Display_Classification": row.get("Display_Classification", ""),
            "Classification_Confidence": row.get("Classification_Confidence", ""),
            "Display_Classification_Source": row.get(
                "Display_Classification_Source", ""
            ),
            "Classification_Status": row.get("Classification_Status", ""),
            "Classification_Reliability_Reason": row.get(
                "Classification_Reliability_Reason", ""
            ),
            "Candidate_Dataset_Source": row.get("Candidate_Dataset_Source", ""),
            "Candidate_Library_Spectrum": row.get(
                "Candidate_Library_Spectrum", ""
            ),
            "Candidate_InChIKey": row.get("Candidate_InChIKey", ""),
            "Candidate_SMILES": row.get("Candidate_SMILES", ""),
            "Spectrum": row.get("Spectrum", ""),
            "Export_Note": "; ".join(dict.fromkeys(notes)),
        }
        projected["Display_Classification_Source"] = _reader_classification_source(
            projected["Display_Classification_Source"], language
        )
        projected["Classification_Reliability_Reason"] = _reader_classification_reason(
            projected["Classification_Reliability_Reason"], language
        )
        output.append({
            column: projected.get(column, "")
            for column in MS2_READER_COLUMN_NAMES
        })
    return output


__all__ = [
    "docker_available", "honest_annotation", "compound_confidence",
    "exact_structure_confidence", "two_dimensional_connectivity_confidence",
    "isomer_group_confidence", "identification_scope",
    "reportable_identity_confidence", "reportable_identity_scope_label",
    "reportable_identity_label", "class_confidence",
    "evidence_tier", "annotation_level",
    "confidence_badge_css", "confidence_badge_html",
    "final_claim_decision", "final_claim_qvalue", "final_claim_badge_html",
    "classification_evidence_status", "fragmentation_evidence_status",
    "svg_for_smiles", "render_db_pathway", "build_diagnostic_ion_audit",
    "render_diagnostic_ion_audit", "render_cfm_cascade", "render_fragmentation_viz",
    "isomer_candidate_projection", "isomer_candidate_export_rows",
    "isomer_candidate_export_rows_from_compounds",
    "isomer_candidate_export_labels", "ISOMER_CANDIDATE_EXPORT_COLUMNS",
    "render_isomer_candidate_competition_html",
    "MS2_UNIFIED_CANDIDATE_COLUMNS", "MS2_UNIFIED_COLUMN_NAMES",
    "ms2_unified_candidate_rows", "MS2_READER_COLUMN_NAMES",
    "ms2_reader_column_labels", "ms2_reader_rows",
    "RESULT_COLUMNS", "RESULTS_COLUMN_NAMES", "REPORT_COLUMN_NAMES", "row_to_dict",
    "GOVERNANCE_BLOCKED_EXPORT_COLUMNS",
    "project_governance_blocked_candidates", "governance_blocked_export_rows",
    "append_governance_blocked_csv_section", "render_governance_blocked_html",
    "inject_governance_blocked_html",
]
