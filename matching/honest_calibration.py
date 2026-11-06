"""
mssystem v3 — 诚实性校准 (P4, Slice A 结构诚实)
==============================================
对 search() 返回的 best + top_candidates 计算"诚实呈现"字段, 解决系统过度声称:
  - 异构体未解析 (top-K 同分子式不同结构, 当前证据不足以点名) → 异构体组置信 + 具体结构歧义
  - 弱/类似物匹配 → putative 而非虚假自信的具体化合物

数据依据 (eval/bench_stratified_n100.json): cos≥0.99 在 clean_hit 占 78%(h1=100%) 但在
isomer_inseparable 占 66%、analogue_dominant 占 45%(h1=0%) → cos 无法区分"正确自信命中"与
"异构/类似物 miss", 只有结构信号(异构集+类共识)能挽回诚实度。故本模块纯结构规则, 无数值校准。

本模块只消费 search 输出并保持 pure-function / 无 I/O；检索端负责保证 Top-K
首项与最终所选 best 一致，呈现端仍对历史缓存做兼容重排。
集成: web/server.py (序列化) + web/report.py (HTML) + web/static/js/results.js (UI) 共用本单一真源,
取代原先分歧的 scores.confidence(API) 与 report.py local confidence_label(无 FP 降级)。

零排序风险: 只计算呈现字段, 不改 final_score / 不改匹配结果 → 红线 11/14 不受影响。
"""
from __future__ import annotations
from collections import Counter
from typing import Optional
from matching.library_quality import compound_identity_key, resolve_structure_formula
from matching.isomer_resolution import (
    connectivity_block,
    isomer_relationship,
    structure_definition_is_unspecified,
    summarize_isomer_relationships,
)
from matching.isomer_literature import build_isomer_literature_context
from matching.classification_provenance import trusted_analysis_class

# 复用既有类规范化 (不重造): 与 search.py:574 display_class 解析一致
from mssystem.classification import is_source_or_generic, CANONICAL_CLASSES
from matching.library_quality import source_tier_confidence_cap

# P-E: pathway-aware 置信度提升开关与阈值 (纯增量; 关=当前行为)
try:
    from config.settings import (
        USE_PATHWAY_CONFIDENCE, PATHWAY_CONF_HIGH_ALIGN, PATHWAY_CONF_HIGH_COS,
        PATHWAY_CONF_HIGH_FINAL, PATHWAY_CONF_HIGH_PPM,
        USE_MARGIN_DOWNGRADE, MARGIN_DOWNGRADE_THRESH,
        USE_MARGIN_UPGRADE, MARGIN_UPGRADE_THRESH,
        USE_RT_STANDARD_CONFIDENCE, RT_CONF_HIGH_COS, RT_CONF_HIGH_FINAL,
        RT_CONF_HIGH_PPM, USE_STRONG_CANDIDATE_TIER,
        STRONG_CANDIDATE_COS, STRONG_CANDIDATE_FINAL,
        STRONG_CANDIDATE_PPM, STRONG_CANDIDATE_CLASS_AGREE,
        PATHWAY_STRONG_ALIGN,
        USE_ISOPAIR_HIGH_PRECISION_VETO,
        ISOMER_AMBIGUITY_MAX_COSINE_GAP, ISOMER_AMBIGUITY_MAX_RANK_GAP,
    )
except Exception:  # 独立调用 / 单测时无 config → 关闭 P-E/margin, 走当前行为
    USE_PATHWAY_CONFIDENCE = False
    PATHWAY_CONF_HIGH_ALIGN = 0.5
    PATHWAY_CONF_HIGH_COS = 0.8
    PATHWAY_CONF_HIGH_FINAL = 0.55
    PATHWAY_CONF_HIGH_PPM = 50
    USE_MARGIN_DOWNGRADE = False
    MARGIN_DOWNGRADE_THRESH = 0.10
    USE_MARGIN_UPGRADE = False
    MARGIN_UPGRADE_THRESH = 0.10
    USE_RT_STANDARD_CONFIDENCE = False
    RT_CONF_HIGH_COS = 0.80
    RT_CONF_HIGH_FINAL = 0.50
    RT_CONF_HIGH_PPM = 50.0
    USE_STRONG_CANDIDATE_TIER = True
    STRONG_CANDIDATE_COS = 0.85
    STRONG_CANDIDATE_FINAL = 0.55
    STRONG_CANDIDATE_PPM = 50.0
    STRONG_CANDIDATE_CLASS_AGREE = 0.50
    PATHWAY_STRONG_ALIGN = 0.50
    USE_ISOPAIR_HIGH_PRECISION_VETO = True
    ISOMER_AMBIGUITY_MAX_COSINE_GAP = 0.05
    ISOMER_AMBIGUITY_MAX_RANK_GAP = 0.00

# ── 可调阈值 (FP 阈值逐字抄自 web/server.py:621-627 保红线稳定) ──
COS_PUTATIVE = 0.70          # cos < 此值 → putative (镜像 scores.confidence Low 下限)
CLASS_AGREE_DOWNGRADE = 0.50 # 类共识一致率 < 此值 → 降化合物置信
PPM_TO_LOW = 200             # |mass_ppm| > 此 → Low (质量明显不符)
PPM_HIGH_TO_MED = 50         # |mass_ppm| > 此且 High → Medium (与严格专家/RT门一致)
REV_HIGH_TO_MED = 0.6        # reverse_match < 此且 High → Medium (lib 特征峰 query 没有)

_RANK = {"Unknown": 0, "Low": 1, "Medium": 2, "High": 3}
_NA_SMILES = {"", "nan", "none", "null"}
_COMPARISON_SCORE_FIELDS = (
    "cosine",
    "entropy",
    "fragment_score",
    "intensity_score",
    "reverse_match",
    "final_score",
    "rank_score",
)


def candidate_comparison_view(candidate: dict) -> dict:
    """Return the immutable raw-spectrum view used for candidate comparisons.

    A successfully deconvolved Top-1 may carry stronger absolute evidence in
    ``cosine``/``final_score``.  Those values must never widen its Top-1 margin,
    reorder Top-K, or make a same-formula competitor disappear.  Search stores
    the pre-deconvolution values under ``raw_*``; this helper restores them for
    every relative/candidate-to-candidate decision while preserving metadata.
    """
    view = dict(candidate or {})
    raw_used = False
    for field in _COMPARISON_SCORE_FIELDS:
        raw_value = view.get(f"raw_{field}")
        if raw_value is not None:
            view[field] = raw_value
            raw_used = True
    already_frozen = (
        view.get("candidate_comparison_source") == "raw_measured_spectrum"
    )
    view["candidate_comparison_source"] = (
        "raw_measured_spectrum"
        if raw_used or already_frozen else "measured_spectrum"
    )
    return view


def _has_struct(best: dict) -> bool:
    """是否有可信结构 (SMILES 非空非占位)。与 server.py 无结构降级逻辑一致。"""
    if best.get("has_structure") is False:
        return False
    smi = best.get("smiles")
    return bool(smi and str(smi).strip().lower() not in _NA_SMILES)


def _canon(raw) -> str:
    """原始库 classification → 规范类 (CANONICAL_CLASSES 小写集合); 泛化/源名 → 'other'。"""
    s = (str(raw) if raw is not None else "").strip().lower()
    if not s or s in _NA_SMILES:
        return "other"
    if s in CANONICAL_CLASSES and not is_source_or_generic(raw):
        return s
    return "other"


def _trusted_class_vote(candidate: dict) -> str:
    """Return a class vote only when its provenance can influence decisions."""
    resolved, reliable, _reason = trusted_analysis_class(
        candidate,
        # Old unit-test/direct-call fixtures omitted this field entirely.  Live
        # search rows always populate it (usually as ``unknown``), so legacy
        # compatibility cannot turn a production missing source into trust.
        allow_legacy_missing=True,
    )
    return _canon(resolved) if reliable else "other"


def _eligible_pathway_score(candidate: dict) -> float:
    """Return only independently mined, explicitly production-eligible score."""
    pathway = ((candidate.get("expert_details") or {}).get("pathway_align") or {})
    if not isinstance(pathway, dict):
        return 0.0
    if not (
        pathway.get("candidate_conditioned") is False
        and pathway.get("independent_evidence") is True
        and pathway.get("eligible_for_fusion") is True
        and pathway.get("runtime_influence") is True
    ):
        return 0.0
    try:
        return max(0.0, min(1.0, float(pathway.get("score", 0.0) or 0.0)))
    except (TypeError, ValueError):
        return 0.0


def _base_confidence(cos: float, final: float) -> str:
    """scores.confidence 非专家分支 (与 matching/scores.py:276-282 一致)。"""
    if cos >= 0.9 and final >= 0.5:
        return "High"
    if cos >= 0.8 and final >= 0.4:
        return "Medium"
    if cos >= 0.7:
        return "Low"
    return "Unknown"


def _apply_confidence_cap(value: str, cap: Optional[str]) -> str:
    if not cap:
        return value
    return value if _RANK.get(value, 0) <= _RANK.get(cap, 0) else cap


def _fp_downgrade(conf: str, ppm, rev) -> str:
    """收编 server.py:621-627 假阳性降级 (ppm / reverse_match)。Unknown 不动。"""
    if conf == "Unknown":
        return conf
    if ppm is not None:
        ap = abs(ppm)
        if ap > PPM_TO_LOW:
            return "Low"
        if ap > PPM_HIGH_TO_MED and conf == "High":
            conf = "Medium"
    if rev is not None and rev < REV_HIGH_TO_MED and conf == "High":
        conf = "Medium"
    return conf


def _abs_ppm_ok(ppm, limit: float) -> bool:
    if ppm is None:
        return True
    try:
        return abs(float(ppm)) <= float(limit)
    except (TypeError, ValueError):
        return True


def _candidate_rank_value(candidate: dict) -> float:
    """Comparable presentation rank without importing search.py circularly."""
    candidate = candidate_comparison_view(candidate)
    for key in ("rank_score", "final_score", "cosine"):
        try:
            value = float(candidate.get(key, 0) or 0)
            if value == value:
                return value
        except (TypeError, ValueError):
            continue
    return 0.0


def _ordered_distinct_candidates(best: dict, top_candidates: list[dict]) -> list[dict]:
    """Return best once, followed by distinct alternative identities by rank.

    Production search normally stores best at Top-K position zero, but the
    final isomer tie-break can select a different near-equal candidate.  Old
    calibration then used the selected best itself as ``tc[1]`` and computed a
    zero margin.  This helper also prevents the best class from receiving two
    votes merely because Top-K already contains it.
    """
    best_key = compound_identity_key(best)
    seen = {best_key} if best_key else set()
    alternatives: list[dict] = []
    for candidate in top_candidates or []:
        if not isinstance(candidate, dict):
            continue
        key = compound_identity_key(candidate)
        if key and key in seen:
            continue
        if key:
            seen.add(key)
        alternatives.append(candidate)
    alternatives.sort(
        key=lambda item: (
            -_candidate_rank_value(item),
            -float(
                candidate_comparison_view(item).get("cosine", 0) or 0
            ),
            compound_identity_key(item),
        )
    )
    return [best, *alternatives]


def is_spectrally_competitive_isomer(
    best_cosine: float,
    best_rank: float,
    candidate: dict,
) -> bool:
    """Conservatively decide whether MS/MS leaves an isomer competitive.

    Cosine and the production ranking score carry overlapping but non-identical
    evidence.  Exact-structure separability requires a clear loss on both;
    absent/legacy scores therefore remain ambiguous rather than silently High.
    """
    candidate = candidate_comparison_view(candidate)
    try:
        candidate_cosine = float(candidate.get("cosine"))
    except (TypeError, ValueError):
        # Missing scores in a legacy cache cannot prove separability.
        return True
    if candidate_cosine != candidate_cosine:
        return True
    try:
        candidate_rank = _candidate_rank_value(candidate)
    except (TypeError, ValueError):
        return True
    cosine_gap = best_cosine - candidate_cosine
    rank_gap = best_rank - candidate_rank
    return (
        cosine_gap <= ISOMER_AMBIGUITY_MAX_COSINE_GAP
        or rank_gap <= ISOMER_AMBIGUITY_MAX_RANK_GAP
    )


# Backward-compatible private alias for existing callers and cached notebooks.
_is_spectrally_competitive_isomer = is_spectrally_competitive_isomer


def normalize_isomer_fragment_resolution_status(
    detail: dict | None,
) -> str:
    """Apply symmetric evidence gates to current and cached fragment details."""
    if not isinstance(detail, dict):
        return ""
    status = str(detail.get("resolution_status") or "")
    replication_complete = (
        (detail.get("reference_replication") or {}).get("status")
        == "complete"
    )
    metadata_complete = detail.get("metadata_comparability") == "complete"
    if status in {
        "diagnostic_evidence_favors_competitor",
        "competitor_supported_pending_validation",
    }:
        if not replication_complete:
            status = "competitor_features_observed_but_unreplicated"
        elif not metadata_complete:
            status = "competitor_features_observed_but_metadata_incomplete"
        else:
            status = "competitor_supported_pending_validation"
    elif status == "conflicting_diagnostic_evidence":
        if not replication_complete:
            status = "conflicting_features_observed_but_unreplicated"
        elif not metadata_complete:
            status = (
                "conflicting_features_observed_but_metadata_incomplete"
            )
    elif status == "selected_candidate_supported_pending_validation":
        if not replication_complete:
            status = "selected_candidate_supported_but_unreplicated"
        elif not metadata_complete:
            status = (
                "selected_candidate_supported_but_metadata_incomplete"
            )
    detail["resolution_status"] = status
    return status


def annotate(best: Optional[dict], top_candidates: Optional[list[dict]] = None) -> dict:
    """计算诚实呈现字段。

    Args:
        best: search() 返回的 best 匹配 dict (cosine/final_score/mass_ppm/reverse_match/
              smiles/formula/inchikey/compound_name/display_class/classification)。
        top_candidates: search() 返回的 top-K (None 则读 best['top_candidates'])。
    Returns:
        flat dict, 各 key 安全可并入 serialized match:
        isomer_set / isomer_n / isomer_ambiguous / class_consensus / class_agreement /
        compound_confidence / class_confidence / honesty_tier / honest_headline。
    """
    empty = dict(isomer_set=[], isomer_n=0, isomer_ambiguous=False,
                 same_formula_alternatives=[], same_formula_alternatives_n=0,
                 noncompetitive_same_formula_n=0,
                 isomer_ambiguity_type="none",
                 isomer_relationship_counts={},
                 isomer_variant_subtype_counts={},
                 isomer_resolution_status="not_applicable",
                 class_consensus="other", class_agreement=0.0,
                 compound_confidence="Unknown", class_confidence="Unknown",
                 confidence_semantics="rule_based_operational_evidence_tier",
                 confidence_is_probability=False,
                 confidence_is_project_qvalue=False,
                 legacy_unvalidated_expert_influence_detected=False,
                 margin_upgrade_counterfactual_eligible=False,
                 margin_upgrade_runtime_influence=False,
                 exact_structure_confidence="Unknown",
                 isomer_group_confidence="Unknown",
                 identification_scope="unmatched",
                 isomer_group_formula="", isomer_group_size=0,
                 isomer_group_evidence={},
                 isomer_literature_evidence={},
                 honesty_tier="putative", honest_headline="未匹配",
                 evidence_tier="Putative", annotation_level="Unmatched",
                 candidate_comparison_source="measured_spectrum",
                 candidate_comparison_cosine=None, candidate_margin=None,
                 rt_confirmed=False, rt_warning=None,
                 deconvolution_score_influence=False,
                 deconvolution_evidence_state="not_applicable",
                 chimeric_isolation_warning=None,
                 isopair_reranked_2d_connectivity=False,
                 isopair_high_precision_gate_required=False,
                 isopair_high_precision_support=False,
                 isopair_high_precision_gate_result="not_applicable",
                 isopair_high_eligible=False,
                 isopair_evidence={},
                 precursor_high_evidence_gate_required=False,
                 precursor_high_evidence_gate_result="not_applicable",
                 precursor_high_evidence_gate_reason=None,
                 precursor_high_evidence_gate_reason_zh=None,
                 direct_matched_product_ions=None,
                 two_dimensional_connectivity_confidence="Unknown")
    if not best or not isinstance(best, dict):
        return empty

    tc = top_candidates if top_candidates is not None else (best.get("top_candidates") or [])
    for candidate in [best, *tc]:
        if isinstance(candidate, dict):
            normalize_isomer_fragment_resolution_status(
                candidate.get("isomer_fragment_detail")
            )
    ordered_candidates = _ordered_distinct_candidates(best, tc)
    alternatives = ordered_candidates[1:]
    cos = float(best.get("cosine", 0) or 0)
    fin = float(best.get("final_score", 0) or 0)
    comparison_best = candidate_comparison_view(best)
    comparison_cos = float(comparison_best.get("cosine", 0) or 0)
    comparison_source = comparison_best["candidate_comparison_source"]
    best_rank = _candidate_rank_value(comparison_best)
    ppm = best.get("mass_ppm")
    rev = best.get("reverse_match")
    bf = resolve_structure_formula(
        str(best.get("smiles") or ""),
        str(best.get("formula") or ""),
        str(best.get("adduct") or ""),
    )
    bik = compound_identity_key(best)
    best_name = str(best.get("compound_name") or "") or "未命名"
    has_struct = _has_struct(best)
    rt_ev = best.get("rt_evidence") or {}
    rt_confirmed = False
    isopair_evidence = best.get("isopair_evidence") or {}
    isopair_selected_connectivity = str(
        isopair_evidence.get("selected_connectivity") or ""
    )
    isopair_model_evaluated = bool(
        isopair_evidence.get("evaluated")
        and isopair_evidence.get("identity_scope") == "2d_connectivity_only"
        and isopair_selected_connectivity
        and isopair_selected_connectivity == connectivity_block(best)
    )
    isopair_high_precision_gate_applicable = bool(
        USE_ISOPAIR_HIGH_PRECISION_VETO
        and isopair_model_evaluated
        and isopair_evidence.get("high_precision_gate_applicable")
    )
    isopair_reranked_2d_connectivity = bool(
        isopair_model_evaluated and isopair_evidence.get("ranking_changed")
    )
    isopair_high_precision_support = bool(
        isopair_high_precision_gate_applicable
        and best.get("isopair_high_precision_support")
        and isopair_evidence.get("high_precision_support")
        and not isopair_evidence.get("ranking_changed")
        and not isopair_evidence.get("high_promotion_allowed")
    )
    isopair_high_eligible = isopair_high_precision_support
    expert_details = best.get("expert_details") or {}
    expert_fusion = (
        expert_details.get("score_fusion")
        if isinstance(expert_details, dict) else {}
    ) or {}
    pathway_payload = (
        expert_details.get("pathway_align")
        if isinstance(expert_details, dict) else None
    )
    v3_payload = (
        expert_details.get("v3")
        if isinstance(expert_details, dict) else None
    )
    try:
        base_expert_score = float(
            expert_details.get("expert_total_score", 0.0) or 0.0
        )
    except (TypeError, ValueError):
        base_expert_score = 0.0
    try:
        untrusted_v3_score = (
            float(v3_payload.get("v3_score", 0.0) or 0.0)
            if isinstance(v3_payload, dict)
            and v3_payload.get("classification_reliable") is not True
            else 0.0
        )
    except (TypeError, ValueError):
        untrusted_v3_score = 0.0
    try:
        legacy_pathway_score = (
            float(pathway_payload.get("score", 0.0) or 0.0)
            if isinstance(pathway_payload, dict) else 0.0
        )
    except (TypeError, ValueError):
        legacy_pathway_score = 0.0
    legacy_unvalidated_expert_influence_detected = bool(
        (
            isinstance(pathway_payload, dict)
            and legacy_pathway_score > 0.0
            and expert_fusion.get("method")
            != "max_of_eligible_expert_sources"
        )
        or untrusted_v3_score > base_expert_score + 1e-9
    )
    rt_mismatch = bool(
        USE_RT_STANDARD_CONFIDENCE
        and rt_ev.get("method_comparable")
        and rt_ev.get("identity_scope") == "exact_structure"
        and rt_ev.get("chromatographic_status") == "mismatch"
    )

    # A high-similarity named library spectrum without a validated structure is
    # useful as a traceable lead, but it is neither a structure identification
    # nor a generic "putative compound".  Keep this distinction in the single
    # honesty source so serialization and legacy-cache upgrades cannot diverge.
    if best.get("name_only_lead") and not has_struct:
        return dict(
            isomer_set=[], isomer_n=0, isomer_ambiguous=False,
            same_formula_alternatives=[], same_formula_alternatives_n=0,
            noncompetitive_same_formula_n=0,
            isomer_ambiguity_type="none",
            isomer_relationship_counts={},
            isomer_variant_subtype_counts={},
            isomer_resolution_status="not_applicable",
            class_consensus="other", class_agreement=0.0,
            compound_confidence="Unknown", class_confidence="Unknown",
            confidence_semantics="rule_based_operational_evidence_tier",
            confidence_is_probability=False,
            confidence_is_project_qvalue=False,
            legacy_unvalidated_expert_influence_detected=False,
            margin_upgrade_counterfactual_eligible=False,
            margin_upgrade_runtime_influence=False,
            exact_structure_confidence="Unknown",
            isomer_group_confidence="Unknown",
            identification_scope="spectral_lead",
            isomer_group_formula="", isomer_group_size=0,
            isomer_group_evidence={},
            isomer_literature_evidence={},
            honesty_tier="name_only_lead",
            honest_headline=f"{best_name}（仅谱图线索，无可信结构）",
            evidence_tier="Spectral Lead",
            annotation_level="Library spectrum lead — no validated structure",
            candidate_comparison_source=comparison_source,
            candidate_comparison_cosine=round(comparison_cos, 6),
            candidate_margin=None,
            rt_confirmed=False, rt_warning=None,
            deconvolution_score_influence=False,
            deconvolution_evidence_state="not_applicable",
            chimeric_isolation_warning=None,
            precursor_high_evidence_gate_required=False,
            precursor_high_evidence_gate_result="not_applicable",
            precursor_high_evidence_gate_reason=None,
            precursor_high_evidence_gate_reason_zh=None,
            direct_matched_product_ions=None,
        )

    # ── (1) 异构集: top-K 中所有同 formula 且完整身份不同者 (全集, 非仅 top-2) ──
    same_formula_alternatives = []
    isomer_set = []
    if bf:
        for c in alternatives:
            cf = resolve_structure_formula(
                str(c.get("smiles") or ""),
                str(c.get("formula") or ""),
                str(c.get("adduct") or ""),
            )
            cik = compound_identity_key(c)
            if cf and cf == bf and cik and cik != bik:
                item = {
                    "compound_name": str(c.get("compound_name") or "")[:60],
                    "inchikey": c.get("inchikey") or "",
                    "formula": cf,
                    "cosine": round(
                        float(candidate_comparison_view(c).get("cosine", 0) or 0),
                        4,
                    ),
                    "rank_score": round(_candidate_rank_value(c), 4),
                    "spectrally_competitive": _is_spectrally_competitive_isomer(
                        comparison_cos, best_rank, c
                    ),
                    "relationship": isomer_relationship(best, c),
                }
                same_formula_alternatives.append(item)
                if item["spectrally_competitive"]:
                    isomer_set.append(item)
    isomer_n = len(isomer_set)
    isomer_ambiguous = isomer_n >= 1
    competitive_candidates = [
        candidate
        for candidate in alternatives
        if (
            resolve_structure_formula(
                str(candidate.get("smiles") or ""),
                str(candidate.get("formula") or ""),
                str(candidate.get("adduct") or ""),
            ) == bf
            and compound_identity_key(candidate)
            != compound_identity_key(best)
            and _is_spectrally_competitive_isomer(
                comparison_cos, best_rank, candidate
            )
        )
    ] if bf else []
    relationship_summary = summarize_isomer_relationships(
        best, competitive_candidates
    )
    ambiguity_type = relationship_summary["ambiguity_type"]
    # A full-InChIKey near tie is not automatically a 2-D connectivity tie.
    # When every competitive candidate shares the selected candidate's first
    # InChIKey block, only stereo/isotope/tautomer/definition layers remain
    # unresolved.  Preserve that full-identity boundary without discarding
    # evidence that is valid at two-dimensional connectivity scope.
    same_connectivity_only = bool(
        isomer_ambiguous
        and ambiguity_type == "same_connectivity_variant_group"
        and connectivity_block(best)
    )
    variant_subtype_counts = relationship_summary[
        "same_connectivity_variant_subtype_counts"
    ]
    definition_specificity_ambiguity = bool(
        variant_subtype_counts.get(
            "structure_definition_specificity_variant"
        )
    )
    protonation_layer_ambiguity = bool(
        variant_subtype_counts.get(
            "protonation_or_charge_layer_variant"
        )
    )
    stereo_or_definition_ambiguity = bool(
        variant_subtype_counts.get(
            "stereochemical_isotopic_or_tautomeric_variant"
        )
        or variant_subtype_counts.get(
            "same_connectivity_variant_unresolved"
        )
    )
    rt_exact_eligible = bool(
        USE_RT_STANDARD_CONFIDENCE
        and rt_ev.get("matched")
        and rt_ev.get("eligible_for_exact_identity")
    )
    same_connectivity_requires_stereo_separation = bool(
        isomer_ambiguous
        and stereo_or_definition_ambiguity
        and not rt_ev.get("stereochemistry_capable")
    )
    underspecified_selected_identity = bool(
        definition_specificity_ambiguity
        and structure_definition_is_unspecified(best)
    )
    rt_confirmed = bool(
        rt_exact_eligible
        and not same_connectivity_requires_stereo_separation
        and not underspecified_selected_identity
    )

    # ── (2) 类共识: top-K (含 best) 规范类多数投票 + 一致率 ──
    votes = [_trusted_class_vote(c) for c in ordered_candidates]
    candidate_class_hypotheses = dict(
        sorted(Counter(vote for vote in votes if vote != "other").items())
    )
    top1_class_vote = votes[0] if votes else "other"
    # A lower-ranked candidate's curated/structure class is evidence about
    # that alternative, not about the selected Top-1.  Without independent
    # support for Top-1 itself, transferring the alternative's class creates
    # exactly the LPA->glycoside anomaly seen in the current sample.
    pool = (
        [vote for vote in votes if vote != "other"]
        if top1_class_vote != "other"
        else []
    )
    if pool:
        cls_consensus, cls_cnt = Counter(pool).most_common(1)[0]
        cls_agreement = cls_cnt / len(pool)
        class_consensus_status = "top1_independently_supported"
    else:
        cls_consensus, cls_agreement = "other", 0.0
        class_consensus_status = "top1_unverified_no_cross_candidate_transfer"
    class_evidence_available = bool(pool)
    # 展示类: 共识优先, 退回 best display_class
    cls_disp = cls_consensus if cls_consensus != "other" else (
        _canon(best.get("display_class") or best.get("classification")))
    if cls_disp == "other":
        cls_disp = best.get("display_class") or "other"

    # ── (3) compound_confidence: base + FP + 诚实降级 ──
    comp = _base_confidence(cos, fin)
    margin_upgrade_counterfactual_eligible = False
    if not has_struct:
        comp = "Unknown"
    else:
        comp = _fp_downgrade(comp, ppm, rev)
        if comp != "Unknown":
            if best.get("mass_metadata_consistent") is False and comp == "High":
                comp = "Medium"
            if isomer_ambiguous and comp == "High":           # 异构 → 无法点名具体化合物
                comp = "Medium"
            if cos < COS_PUTATIVE:                             # 弱匹配
                comp = "Low"
            if (
                class_evidence_available
                and cls_agreement < CLASS_AGREE_DOWNGRADE
                and comp == "High"
            ):
                comp = "Medium"

        # ── (3b) P-E: pathway-aware 置信度提升 (纯增量; 关=当前行为) ──
        # 用 P-A 通路对齐开一条不依赖 cos≥0.9 的 High 通路 (直击"高置信度偏少": cos<0.9 的真匹配
        # 被卡 Medium)。守卫与 P4 一致避免过度声称: not isomer_ambiguous(当前证据未解析异构体) +
        # cls_agreement≥0.5(类共识一致, 对齐 P4 降级条件) + 质量符 + 结构存在。仅 Medium→High。
        if USE_PATHWAY_CONFIDENCE and comp == "Medium" and has_struct:
            _pa = _eligible_pathway_score(best)
            if (_pa >= PATHWAY_CONF_HIGH_ALIGN
                    and cos >= PATHWAY_CONF_HIGH_COS
                    and fin >= PATHWAY_CONF_HIGH_FINAL
                    and (ppm is None or abs(ppm) <= PATHWAY_CONF_HIGH_PPM)
                    and best.get("mass_metadata_consistent") is not False
                    and not isomer_ambiguous
                    and class_evidence_available
                    and cls_agreement >= CLASS_AGREE_DOWNGRADE):
                comp = "High"

        # ── (3c) margin-aware 置信度校准 (cosine top1−top2 分差; presentation-only) ──
        # 清晰胜出(大 margin)→提档; 模糊(小 margin)→降档防过度声称。依据
        # eval/struct_inference/validate_margin_calibration.py (443q): High 档 margin<0.10
        # 仅 0.636 正确, ≥0.10 者 0.843; Med 档 margin≥0.10 者 0.84(≥现 High)。升档守卫同 P4/P-E。
        if USE_MARGIN_DOWNGRADE or USE_MARGIN_UPGRADE:
            _margin = None
            if alternatives:
                try:
                    competitor = candidate_comparison_view(alternatives[0])
                    _margin = (
                        comparison_cos
                        - float(competitor.get("cosine", 0) or 0)
                    )
                except (TypeError, ValueError):
                    _margin = None
            if _margin is not None:
                if (USE_MARGIN_DOWNGRADE and comp == "High"
                        and _margin < MARGIN_DOWNGRADE_THRESH):
                    comp = "Medium"
                elif (USE_MARGIN_UPGRADE and comp == "Medium"
                        and _margin >= MARGIN_UPGRADE_THRESH
                        and not isomer_ambiguous
                        and class_evidence_available
                        and cls_agreement >= CLASS_AGREE_DOWNGRADE
                        and best.get("mass_metadata_consistent") is not False
                        and _abs_ppm_ok(ppm, PPM_HIGH_TO_MED)):
                    # The 443-query analysis selected this threshold on the
                    # same development cohort and did not calibrate a High
                    # operating point.  Preserve the counterfactual for audit,
                    # but never alter production confidence.
                    margin_upgrade_counterfactual_eligible = True

        # ── (3d) v5: local RT standard evidence ──
        # RT is orthogonal to MS/MS and can resolve otherwise inseparable isomers,
        # but only when it comes from same-method local standards.  Public-library
        # RT is intentionally not used here.
        if rt_mismatch and comp in ("High", "Medium"):
            comp = "Low"
        elif (rt_confirmed and has_struct
              and best.get("mass_metadata_consistent") is not False
              and cos >= RT_CONF_HIGH_COS
              and fin >= RT_CONF_HIGH_FINAL and _abs_ppm_ok(ppm, RT_CONF_HIGH_PPM)):
            comp = "High"

    # ── (3e) isomer_group_confidence: same-formula group, not exact identity ──
    # Use the immutable measured-spectrum view.  A same-formula competitive
    # alternative prevents exact-structure High, but it does not erase strong
    # evidence that the analyte belongs to this explicit molecular-formula
    # isomer group.  No diagnostic shadow score is allowed to promote this
    # field; it is derived only from the already calibrated base evidence.
    comparison_fin = float(comparison_best.get("final_score", 0) or 0)
    comparison_rev = comparison_best.get("reverse_match")
    if not has_struct or not bf:
        group_conf = "Unknown"
    else:
        group_conf = _fp_downgrade(
            _base_confidence(comparison_cos, comparison_fin),
            ppm,
            comparison_rev,
        )
        if (
            best.get("mass_metadata_consistent") is False
            and group_conf == "High"
        ):
            group_conf = "Medium"
        if comparison_cos < COS_PUTATIVE:
            group_conf = "Low"
        if rt_mismatch and group_conf in ("High", "Medium"):
            group_conf = "Low"
        elif (
            rt_confirmed
            and best.get("mass_metadata_consistent") is not False
            and cos >= RT_CONF_HIGH_COS
            and fin >= RT_CONF_HIGH_FINAL
            and _abs_ppm_ok(ppm, RT_CONF_HIGH_PPM)
        ):
            group_conf = "High"
    if not isomer_ambiguous:
        group_conf = comp

    # ── (4) class_confidence: base + FP, 异构不降级 (异构同类), 类共识约束 ──
    if not has_struct:
        cconf = "Unknown"
    else:
        cconf = _fp_downgrade(_base_confidence(cos, fin), ppm, rev)
        if cconf != "Unknown":
            if cls_consensus == "other":                      # 无可靠类信息 → 类级不可信
                cconf = "Low"
            elif cls_agreement < CLASS_AGREE_DOWNGRADE and cconf == "High":
                cconf = "Medium"
            # NOTE: 异构不降级 class (异构体同类, 类级稳健)

    # ── (4a) Source-quality hard cap ──
    # Must run after pathway/margin/RT upgrades.  A generated or presence-only
    # library entry cannot become an experimental Level 2a match merely because
    # cosine is high or the Top-1 margin is large.
    source_tier = str(best.get("source_tier") or "unknown").strip().lower()
    source_cap = source_tier_confidence_cap(source_tier)
    comp = _apply_confidence_cap(comp, source_cap)
    group_conf = _apply_confidence_cap(group_conf, source_cap)
    cconf = _apply_confidence_cap(cconf, source_cap)
    if legacy_unvalidated_expert_influence_detected:
        # Historical caches may have baked a self-mined pathway or an
        # unconfirmed ML-class v3 score into ``final_score``.  Without replaying
        # the raw spectrum that contribution cannot be removed exactly, so the
        # cached claim fails closed until reanalysis.
        comp = _apply_confidence_cap(comp, "Medium")
        group_conf = _apply_confidence_cap(group_conf, "Medium")
        cconf = _apply_confidence_cap(cconf, "Medium")
    if (
        str(best.get("lib_precursor_mz_basis") or "")
        == "structure_adduct_theoretical_index"
    ):
        # Theoretical retrieval rescues genuine low-precision/legacy spectra,
        # but a conflicting reported precursor remains a metadata limitation.
        # It can support a putative Level-2 match, never an unqualified High.
        comp = _apply_confidence_cap(comp, "Medium")
        group_conf = _apply_confidence_cap(group_conf, "Medium")
        cconf = _apply_confidence_cap(cconf, "Medium")
    quality_flags = {
        str(flag) for flag in (best.get("quality_flags") or [])
    }
    if "confidence_cap_medium" in quality_flags:
        # A real experimental spectrum can sometimes be recovered only from a
        # publication figure or another explicitly approximate representation
        # of relative intensities. Keep it searchable, but never let later
        # score or RT upgrades turn that curation limit into a High result.
        comp = _apply_confidence_cap(comp, "Medium")
        group_conf = _apply_confidence_cap(group_conf, "Medium")
        cconf = _apply_confidence_cap(cconf, "Medium")
    if "precursor_adduct_conflict_manual_review_required" in quality_flags:
        comp = _apply_confidence_cap(comp, "Low")
        group_conf = _apply_confidence_cap(group_conf, "Low")
        cconf = _apply_confidence_cap(cconf, "Medium")
    if source_tier == "literature_presence_only":
        cconf = "Unknown"

    # A high cosine from three sparse product-ion overlaps is not sufficient
    # for an unqualified High call when the selected precursor was explicitly
    # absent from adjacent continuous MS1 and search had to fall back to the
    # quadrupole selected-ion setpoint.  This gate is deliberately narrow:
    # it applies only when all limitations are recorded and a matchms audit
    # provides an inspectable matched-ion count.  It changes presentation
    # confidence only; candidate retrieval, scores and ranking remain intact.
    multimetric = best.get("multimetric_shadow") or {}
    direct_metric = multimetric.get("direct_cosine_hungarian") or {}
    try:
        direct_matched_product_ions = int(direct_metric.get("matched_peaks"))
    except (TypeError, ValueError):
        direct_matched_product_ions = None
    precursor_source = str(
        best.get("query_search_precursor_mz_source")
        or (best.get("query_precursor_mass_audit") or {}).get(
            "search_precursor_mz_source"
        )
        or ""
    ).strip()
    isotope_status = str(
        best.get("isotope_evidence_status") or ""
    ).strip().lower()
    explicit_query_adduct = str(best.get("query_adduct") or "").strip()
    independent_query_adduct = bool(
        explicit_query_adduct and best.get("query_adduct_reliable")
    )
    sparse_precursor_gate_context = bool(
        isotope_status == "precursor_not_observed_in_ms1"
        and precursor_source == "ms2_selected_ion_mz_fallback"
        and not independent_query_adduct
    )
    precursor_high_evidence_gate_required = bool(
        (comp == "High" or group_conf == "High")
        and sparse_precursor_gate_context
        and direct_matched_product_ions is not None
        and direct_matched_product_ions < 4
        and not rt_confirmed
    )
    precursor_high_evidence_gate_reason = None
    precursor_high_evidence_gate_reason_zh = None
    if precursor_high_evidence_gate_required:
        comp = _apply_confidence_cap(comp, "Medium")
        group_conf = _apply_confidence_cap(group_conf, "Medium")
        cconf = _apply_confidence_cap(cconf, "Medium")
        precursor_high_evidence_gate_result = (
            "vetoed_to_medium_sparse_ms2_without_precursor_support"
        )
        precursor_high_evidence_gate_reason = (
            f"only {direct_matched_product_ions} directly matched product ions; "
            "precursor not observed in adjacent MS1; selected-ion m/z fallback; "
            "no reliable candidate-independent query adduct"
        )
        precursor_high_evidence_gate_reason_zh = (
            f"仅有 {direct_matched_product_ions} 个直接匹配的产物离子；"
            "相邻 MS1 未观测到前体；检索使用选中离子 m/z 回退值；"
            "缺少与候选结构无关的可靠样品加合物证据"
        )
    elif (comp == "High" or group_conf == "High") and sparse_precursor_gate_context:
        precursor_high_evidence_gate_result = (
            "passed_minimum_direct_product_ion_support"
            if direct_matched_product_ions is not None
            else "not_evaluable_missing_matchms_peak_count"
        )
    elif comp == "High" or group_conf == "High":
        precursor_high_evidence_gate_result = "not_required_independent_precursor_context"
    else:
        precursor_high_evidence_gate_result = "not_high_candidate"

    # Deconvolution is an optional positive-evidence channel.  A failed or
    # inconclusive fit means "no score gain", not negative evidence against an
    # otherwise valid raw-spectrum match.  Keep the raw confidence unchanged and
    # expose the unresolved mixed/severe isolation as a separate warning.  This
    # preserves the immutable raw Top-1/margin contract while preventing an
    # unsuccessful optimisation attempt from manufacturing a confidence loss.
    deconvolution = best.get("chimeric_deconvolution") or {}
    isolation_status = str(
        deconvolution.get("isolation_status")
        or (best.get("precursor_isolation_evidence") or {}).get("status")
        or ""
    )
    deconvolution_accepted = bool(
        deconvolution.get("accepted")
        and deconvolution.get("runtime_influence")
        and deconvolution.get("status") == "accepted_score_applied"
    )
    deconvolution_unresolved = (
        isolation_status in {"mixed", "severe"} and not deconvolution_accepted
    )
    deconvolution_status = str(deconvolution.get("status") or "not_run")
    if deconvolution_accepted:
        deconvolution_evidence_state = "score_applied"
    elif deconvolution_unresolved:
        deconvolution_evidence_state = "no_score_change"
    else:
        deconvolution_evidence_state = "not_applicable"
    chimeric_isolation_warning = None
    if deconvolution_unresolved:
        chimeric_isolation_warning = (
            f"{isolation_status} precursor isolation; deconvolution status "
            f"'{deconvolution_status}' did not pass every score-increase gate. "
            "Raw-spectrum confidence retained; no deconvolution score gain applied."
        )

    # An ion-identity network is orthogonal positive evidence only when its
    # compound identity and query/library adduct checks agree.  Conflicts can
    # never coexist with an unqualified High label.
    ion_identity_conflict = bool(best.get("ion_identity_conflict"))
    ion_adduct_conflict = bool(best.get("ion_identity_adduct_conflict"))
    if ion_identity_conflict or ion_adduct_conflict:
        comp = _apply_confidence_cap(comp, "Medium")
        group_conf = _apply_confidence_cap(group_conf, "Medium")
        cconf = _apply_confidence_cap(cconf, "Medium")

    # Any identity selected by the near-tie model remains a Medium compound
    # call unless an authenticated same-method RT standard independently
    # confirms it.  This cap is separate from the selective High precision
    # gate.  Only the explicit-adduct branch may influence ranking; formula-
    # conditioned context remains shadow-only and is ineligible for that gate.
    if isopair_reranked_2d_connectivity and not rt_confirmed:
        comp = _apply_confidence_cap(comp, "Medium")

    # IsoPair-MS has two deliberately separate jobs.  Its explicit-adduct
    # near-tie branch may rerank 2-D connectivities, but never promotes that
    # difficult case to High.  In an already-clear MS/MS High case, the locked
    # selective score is instead used as a precision gate: lack of model
    # support caps only the exact candidate to Medium.  A same-method
    # authenticated RT standard is independent evidence and therefore
    # supersedes this model-only veto.
    isopair_high_precision_gate_required = bool(
        USE_ISOPAIR_HIGH_PRECISION_VETO
        and isopair_high_precision_gate_applicable
        and comp == "High"
        and not rt_confirmed
    )
    if rt_confirmed and isopair_high_precision_gate_applicable:
        isopair_high_precision_gate_result = "rt_standard_supersedes_model_gate"
    elif isopair_high_precision_gate_required:
        if isopair_high_precision_support:
            isopair_high_precision_gate_result = "supported"
        else:
            comp = "Medium"
            isopair_high_precision_gate_result = "vetoed_to_medium"
    elif isopair_high_precision_gate_applicable:
        isopair_high_precision_gate_result = "evaluated_not_high_candidate"
    elif isopair_model_evaluated:
        isopair_high_precision_gate_result = "ranking_only_not_applicable"
    else:
        isopair_high_precision_gate_result = "not_applicable"

    # ── (4b) v5: Strong Candidate presentation tier ──
    _pa = _eligible_pathway_score(best)
    strong_candidate = (
        USE_STRONG_CANDIDATE_TIER and has_struct and comp == "Medium"
        and best.get("mass_metadata_consistent") is not False
        and not isomer_ambiguous
        and class_evidence_available
        and cls_agreement >= STRONG_CANDIDATE_CLASS_AGREE
        and _abs_ppm_ok(ppm, STRONG_CANDIDATE_PPM)
        and (
            (cos >= STRONG_CANDIDATE_COS and fin >= STRONG_CANDIDATE_FINAL)
            or (_pa >= PATHWAY_STRONG_ALIGN and cos >= 0.80 and fin >= 0.50)
        )
    )

    # ── (5) honesty_tier ──
    if source_tier == "literature_presence_only":
        tier = "literature_target_lead"
    elif rt_confirmed and comp == "High":
        tier = "rt_confirmed"
    elif (
        comp == "High"
        and isopair_high_precision_gate_result == "supported"
    ):
        tier = "isopair_connectivity_high"
    elif (not has_struct) or cos < COS_PUTATIVE:
        tier = "putative"
    elif same_connectivity_only and group_conf == "High":
        tier = "connectivity_high_stereo_unresolved"
    elif same_connectivity_only:
        tier = "connectivity_candidate_stereo_unresolved"
    elif isomer_ambiguous:
        tier = "isomer_ambiguous"
    elif deconvolution_accepted and comp == "High":
        tier = "deconvolved_confident"
    elif comp == "High":
        tier = "confident"
    elif comp in ("Medium", "Low") and cconf in ("High", "Medium") and cls_consensus != "other":
        tier = "class_level"
    else:
        tier = "putative"

    # ── (6) honest_headline ──
    if tier == "literature_target_lead":
        headline = f"{best_name} · 文献靶标碎片线索（无实测强度，不能作为实验谱鉴定）"
    elif tier == "rt_confirmed":
        ref = rt_ev.get("standard_name") or "local RT standard"
        delta = rt_ev.get("rt_delta_sec")
        dstr = f", ΔRT {delta:+.1f}s" if isinstance(delta, (int, float)) else ""
        headline = f"{best_name} ({cls_disp}) · RT标准品确认: {ref}{dstr}"
    elif tier == "isopair_connectivity_high":
        headline = (
            f"{best_name} · IsoPair-MS 2-D connectivity High"
            "（不声明立体化学已解析）"
        )
    elif tier == "connectivity_high_stereo_unresolved":
        headline = (
            f"{best_name} · 二维连接结构 High；"
            "立体化学／完整结构定义层尚未解析"
        )
    elif tier == "connectivity_candidate_stereo_unresolved":
        headline = (
            f"{best_name} · 二维连接结构候选；"
            "立体化学／完整结构定义层尚未解析"
        )
    elif tier == "isomer_ambiguous":
        names = [c["compound_name"] for c in isomer_set][:4]
        who = " / ".join([best_name] + names)
        headline = (
            f"异构体组 {group_conf}（{bf}）· 可能为 {who} 之一；"
            f"具体结构 {comp}，仅凭当前 MS/MS 尚未解析"
        )
    elif tier == "deconvolved_confident":
        fraction = deconvolution.get("target_fraction")
        fraction_text = (
            f", target fraction {float(fraction):.2f}"
            if isinstance(fraction, (int, float)) else ""
        )
        headline = (
            f"{best_name} ({cls_disp}) · 嵌合 MS/MS 稀疏去卷积支持"
            f"{fraction_text}"
        )
    elif strong_candidate:
        headline = f"{best_name} ({cls_disp}) · Strong Candidate: 谱图/质量/Top-K类共识支持, 待RT或标品确认"
    elif tier == "class_level":
        headline = f"{cls_disp} (类级推测, cos {cos:.3f}) · 化合物待定 ({best_name})"
    elif tier == "putative":
        headline = f"{cls_disp or 'Unknown'} (putative, cos {cos:.3f})"
    else:  # confident
        headline = f"{best_name} ({cls_disp})"

    if tier == "literature_target_lead":
        evidence_tier = "Literature Target Lead"
        annotation_level = "Presence-only diagnostic-ion list; experimental spectrum/standard required"
    elif tier == "rt_confirmed":
        evidence_tier = "RT-confirmed High"
        annotation_level = "Level 1-like (MS/MS + local RT standard)"
    elif tier == "isopair_connectivity_high":
        evidence_tier = "IsoPair-MS 2-D connectivity High"
        annotation_level = (
            "Level 2a library match at 2-D connectivity scope; "
            "stereochemistry not claimed"
        )
    elif tier == "connectivity_high_stereo_unresolved":
        evidence_tier = "2-D connectivity High"
        annotation_level = (
            "Level 2a library match at 2-D connectivity scope; "
            "stereochemistry/full structure definition unresolved"
        )
    elif tier == "connectivity_candidate_stereo_unresolved":
        evidence_tier = "2-D connectivity candidate"
        annotation_level = (
            "Library candidate at 2-D connectivity scope; "
            "stereochemistry/full structure definition unresolved"
        )
    elif tier == "deconvolved_confident" and comp == "High":
        evidence_tier = "High (deconvolved chimeric MS2)"
        annotation_level = (
            "Level 2a library match — chimeric MS2 deconvolved with stable non-negative L1 fit"
        )
    elif comp == "High" and deconvolution_unresolved:
        evidence_tier = "High (raw MS2; chimeric warning)"
        annotation_level = (
            "Level 2a raw-spectrum library match; mixed/severe precursor isolation "
            "was unresolved and deconvolution contributed no score"
        )
    elif comp == "High":
        evidence_tier = "High"
        annotation_level = "Level 2a library match"
    elif isomer_ambiguous and group_conf == "High":
        evidence_tier = "Isomer-group High"
        annotation_level = (
            "Level 2/3 isomer-group library evidence; exact structure unresolved"
        )
    elif strong_candidate:
        evidence_tier = "Strong Candidate"
        annotation_level = "Level 2 candidate"
    elif cconf == "High":
        evidence_tier = "Class High"
        annotation_level = "Level 3 class-level"
    else:
        evidence_tier = "Putative" if comp in ("Unknown", "Low") else comp
        annotation_level = "Putative"

    if best.get("potential_contaminant"):
        role_label = "External contaminant"
        evidence_tier = f"{role_label} ({evidence_tier})"
        annotation_level = (
            f"{annotation_level}; excluded from natural-product counts"
        )
        headline = f"⚠ 潜在外源污染物：{best_name} · {headline}"

    rt_warning = None
    if rt_mismatch:
        rt_warning = (
            f"RT outside local standard window: ΔRT {rt_ev.get('rt_delta_sec')}s "
            f"(tol {rt_ev.get('tolerance_sec')}s)"
        )

    candidate_margin = None
    if alternatives:
        try:
            candidate_margin = (
                comparison_cos
                - float(
                    candidate_comparison_view(alternatives[0]).get(
                        "cosine", 0
                    ) or 0
                )
            )
        except (TypeError, ValueError):
            candidate_margin = None

    if source_tier == "literature_presence_only":
        identification_scope = "literature_target_lead"
    elif rt_confirmed and comp == "High":
        identification_scope = "exact_structure"
    elif tier == "isopair_connectivity_high":
        identification_scope = "exact_2d_connectivity"
    elif same_connectivity_only:
        identification_scope = "exact_2d_connectivity"
    elif isomer_ambiguous:
        identification_scope = "isomer_group"
    elif comp == "High":
        identification_scope = "exact_structure"
    elif cconf in ("High", "Medium") and cls_consensus != "other":
        identification_scope = "compound_class"
    else:
        identification_scope = "putative_candidate"

    isomer_detail = best.get("isomer_fragment_detail") or {}
    if tier == "isopair_connectivity_high":
        isomer_resolution_status = (
            "isopair_high_precision_supported_2d_connectivity_only"
        )
    elif not isomer_ambiguous:
        isomer_resolution_status = "not_ambiguous"
    elif rt_confirmed:
        isomer_resolution_status = (
            "resolved_by_same_method_exact_rt_standard"
        )
    elif ambiguity_type == "mixed_isomer_group":
        isomer_resolution_status = (
            "mixed_isomer_group_requires_multi_axis_resolution"
        )
    elif same_connectivity_requires_stereo_separation:
        isomer_resolution_status = (
            "same_connectivity_requires_chiral_or_orthogonal_standard"
        )
    elif underspecified_selected_identity:
        isomer_resolution_status = (
            "same_connectivity_structure_definition_unresolved"
        )
    elif definition_specificity_ambiguity:
        isomer_resolution_status = (
            "same_connectivity_structure_definition_specificity_ambiguity"
        )
    elif protonation_layer_ambiguity:
        isomer_resolution_status = (
            "same_connectivity_protonation_variant_requires_standardization"
        )
    else:
        isomer_resolution_status = isomer_detail.get(
            "resolution_status", "unresolved"
        )
    if ambiguity_type == "mixed_isomer_group":
        resolution_requirement = (
            "resolve every axis in the mixed group: complete underspecified "
            "structures, compare exact full-InChIKey standards under the same "
            "LC method, separate different-connectivity candidates by RT or "
            "validated CCS/diagnostic evidence, and demonstrate chiral or "
            "stereoselective capability when a true stereo layer is present"
        )
    elif same_connectivity_requires_stereo_separation:
        resolution_requirement = (
            "exact full-InChIKey reference under the same LC method plus "
            "demonstrated stereochemical separation (for example a chiral "
            "method/co-injection), or validated orthogonal evidence"
        )
    elif underspecified_selected_identity:
        resolution_requirement = (
            "complete the selected structure definition first; an "
            "underspecified InChIKey cannot be promoted to one exact "
            "stereochemical/isotopic identity by ordinary MS/MS"
        )
    elif definition_specificity_ambiguity:
        resolution_requirement = (
            "the competing library entry is structurally underspecified; "
            "confirm the selected complete full-InChIKey with an authenticated "
            "same-method standard. A chiral method is required only if a real "
            "stereochemical competitor is independently established"
        )
    elif protonation_layer_ambiguity:
        resolution_requirement = (
            "standardize the neutral parent, salt/counter-ion and ionization "
            "state, then compare an authenticated standard under the same LC "
            "method; a chiral method is not implied by protonation alone"
        )
    else:
        resolution_requirement = (
            "same-method exact-identity reference standard RT, validated "
            "ion mobility/CCS, or compound-split calibrated diagnostic "
            "fragment evidence"
        )
    literature_evidence = (
        build_isomer_literature_context(
            best,
            competitive_candidates,
            ambiguity_type=ambiguity_type,
            variant_subtype_counts=variant_subtype_counts,
            isomer_detail=isomer_detail,
        )
        if isomer_ambiguous else {}
    )
    isomer_group_evidence = {
        "formula": bf or "",
        "candidate_count": 1 + isomer_n if bf else 0,
        "candidate_names": [
            best_name,
            *[item["compound_name"] for item in isomer_set],
        ][:8] if bf else [],
        "confidence": group_conf,
        "exact_structure_resolved": bool(
            not isomer_ambiguous or (rt_confirmed and comp == "High")
        ),
        "two_dimensional_connectivity_resolved": bool(
            not isomer_ambiguous
            or same_connectivity_only
            or (rt_confirmed and comp == "High")
        ),
        "diagnostic_fragment_status": isomer_detail.get(
            "diagnostic_status", "not_available"
        ),
        "diagnostic_product_ions_supported": isomer_detail.get(
            "unique_supported", 0
        ),
        "diagnostic_neutral_losses_supported": isomer_detail.get(
            "neutral_loss_supported", 0
        ),
        "diagnostic_evidence_runtime_influence": False,
        "ambiguity_type": ambiguity_type,
        "relationship_counts": relationship_summary[
            "relationship_counts"
        ],
        "same_connectivity_variant_subtype_counts": (
            variant_subtype_counts
        ),
        "resolution_status": isomer_resolution_status,
        "opposing_diagnostic_features_supported": isomer_detail.get(
            "opposing_supported_features", 0
        ),
        "reference_replication": isomer_detail.get(
            "reference_replication", {}
        ),
        "metadata_comparability": isomer_detail.get(
            "metadata_comparability", "not_available"
        ),
        "rt_identity_scope": rt_ev.get("identity_scope"),
        "rt_method_status": rt_ev.get("method_status"),
        "rt_stereochemistry_capable": bool(
            rt_ev.get("stereochemistry_capable")
        ),
        "resolution_requirement": (
            resolution_requirement
            if isomer_ambiguous and not rt_confirmed else None
        ),
        "isopair_ms": {
            "evaluated": isopair_model_evaluated,
            "ranking_changed": isopair_reranked_2d_connectivity,
            "adduct_context": isopair_evidence.get("adduct_context") or {},
            "high_precision_gate_applicable": (
                isopair_high_precision_gate_applicable
            ),
            "high_precision_gate_required": (
                isopair_high_precision_gate_required
            ),
            "high_precision_support": isopair_high_precision_support,
            "high_precision_gate_result": (
                isopair_high_precision_gate_result
            ),
            "identity_scope": "2d_connectivity_only",
            "stereochemistry_claimed": False,
        },
        "literature_evidence": literature_evidence,
    }

    exact_structure_confidence = comp
    two_dimensional_connectivity_confidence = comp
    if tier == "isopair_connectivity_high" and not rt_confirmed:
        # The model's supervised label is the first InChIKey block.  Preserve
        # the High label at that validated scope without silently converting it
        # into a full stereochemical identity claim.
        exact_structure_confidence = _apply_confidence_cap(comp, "Medium")
        two_dimensional_connectivity_confidence = "High"
    elif same_connectivity_only and not rt_confirmed:
        # The group-level evidence applies to one shared 2-D connectivity.
        # Keep the unresolved full identity capped while preserving the honest
        # confidence at the connectivity scope.
        exact_structure_confidence = _apply_confidence_cap(comp, "Medium")
        two_dimensional_connectivity_confidence = group_conf

    return dict(
        isomer_set=isomer_set, isomer_n=isomer_n, isomer_ambiguous=isomer_ambiguous,
        same_formula_alternatives=same_formula_alternatives,
        same_formula_alternatives_n=len(same_formula_alternatives),
        noncompetitive_same_formula_n=len(same_formula_alternatives) - isomer_n,
        isomer_ambiguity_type=ambiguity_type,
        isomer_relationship_counts=relationship_summary[
            "relationship_counts"
        ],
        isomer_variant_subtype_counts=variant_subtype_counts,
        isomer_resolution_status=isomer_resolution_status,
        class_consensus=cls_consensus, class_agreement=round(cls_agreement, 3),
        class_consensus_status=class_consensus_status,
        class_candidate_hypotheses=candidate_class_hypotheses,
        compound_confidence=comp, class_confidence=cconf,
        confidence_semantics="rule_based_operational_evidence_tier",
        confidence_is_probability=False,
        confidence_is_project_qvalue=False,
        legacy_unvalidated_expert_influence_detected=(
            legacy_unvalidated_expert_influence_detected
        ),
        margin_upgrade_counterfactual_eligible=(
            margin_upgrade_counterfactual_eligible
        ),
        margin_upgrade_runtime_influence=False,
        exact_structure_confidence=exact_structure_confidence,
        two_dimensional_connectivity_confidence=(
            two_dimensional_connectivity_confidence
        ),
        isomer_group_confidence=group_conf,
        identification_scope=identification_scope,
        isomer_group_formula=bf or "",
        isomer_group_size=(1 + isomer_n if bf else 0),
        isomer_group_evidence=isomer_group_evidence,
        isomer_literature_evidence=literature_evidence,
        honesty_tier=tier, honest_headline=headline,
        evidence_tier=evidence_tier, annotation_level=annotation_level,
        candidate_comparison_source=comparison_source,
        candidate_comparison_cosine=round(comparison_cos, 6),
        candidate_margin=(
            round(candidate_margin, 6)
            if candidate_margin is not None else None
        ),
        rt_confirmed=rt_confirmed, rt_warning=rt_warning,
        deconvolution_score_influence=deconvolution_accepted,
        deconvolution_evidence_state=deconvolution_evidence_state,
        chimeric_isolation_warning=chimeric_isolation_warning,
        isopair_reranked_2d_connectivity=isopair_reranked_2d_connectivity,
        isopair_high_precision_gate_required=(
            isopair_high_precision_gate_required
        ),
        isopair_high_precision_support=isopair_high_precision_support,
        isopair_high_precision_gate_result=(
            isopair_high_precision_gate_result
        ),
        isopair_high_eligible=isopair_high_eligible,
        isopair_evidence=isopair_evidence,
        precursor_high_evidence_gate_required=(
            precursor_high_evidence_gate_required
        ),
        precursor_high_evidence_gate_result=(
            precursor_high_evidence_gate_result
        ),
        precursor_high_evidence_gate_reason=(
            precursor_high_evidence_gate_reason
        ),
        precursor_high_evidence_gate_reason_zh=(
            precursor_high_evidence_gate_reason_zh
        ),
        direct_matched_product_ions=direct_matched_product_ions,
    )
