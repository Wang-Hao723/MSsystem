"""
化学专家评分模块 — Chemical Expert Scorer (v3.0 Enhanced)
==========================================================

超越简单的数值相似度比较，引入质谱裂解原理的专家知识来评估匹配质量。
v3.0 整合了碎片特征标准库和结构解析引擎的知识。

核心评估维度 (v3.0):

1. **裂解规则符合度** — 碎片峰是否可以用已知裂解规则解释
2. **中性丢失化学合理性** — 质量差是否对应合理的化学基团丢失
3. **官能团特征覆盖** — 是否检测到化合物类型预期的特征碎片
4. **键能热力学合理性** — 推测的断键位点是否化学上合理
5. **结构-碎片一致性** — RDKit结构预测的碎片是否与观测一致
6. **碎片特征库验证 (v3.0新增)** — 碎片是否匹配该类化合物的已知特征碎片
7. **异构体区分证据 (v3.0新增)** — 是否有区分异构体的诊断离子对

知识来源 (v3.0):
  - 《质谱分析技术原理与应用》 台湾质谱学会 (2019)
  - 《中药质谱分析》 刘淑莹等 (2012)
  - Spectraverse (Gupta et al., Anal. Chem. 2026)

Author: Enhanced Expert System Module
Version: 3.0.0
"""

from __future__ import annotations

import logging
import math
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from matching.classification_provenance import classification_source_trust

import numpy as np

from .rule_quality import collision_energy_value_unit, normalize_polarity

from .fragmentation_rules import (
    FRAGMENTATION_RULES,
    COMMON_NEUTRAL_LOSSES,
    DIAGNOSTIC_IONS,
    FUNCTIONAL_GROUP_SIGNATURES,
    BOND_DISSOCIATION_ENERGIES,
    get_fg_signature,
    compute_possible_neutral_losses,
    compute_rule_coverage,
    match_diagnostic_ion,
    get_neutral_loss_by_mass,
    get_composite_neutral_loss_by_mass,
    get_rules_by_functional_group,
    estimate_bond_lability,
    normalize_classification,
    _is_nl_compatible,
    _is_diagnostic_ion_compatible,
    NeutralLoss,
    FragmentationRule,
    DiagnosticIon,
)

# v3.0: 整合碎片特征标准库
try:
    from .fragment_feature_library import (
        COMPOUND_CLASS_PROFILES,
        get_profile,
        get_fragments_by_class,
        match_fragment_to_class,
        classify_compound_from_fragments,
        get_diagnostic_ion_pairs,
        get_fragmentation_pathway,
        get_isomer_rules,
        FragmentFeature,
        CompoundClassProfile,
    )
    FEATURE_LIB_AVAILABLE = True
except ImportError:
    FEATURE_LIB_AVAILABLE = False

# v3.0: 整合结构解析引擎
try:
    from .structure_elucidation import (
        StructureElucidationEngine,
        explain_match,
    )
    ELUCIDATION_AVAILABLE = True
except ImportError:
    ELUCIDATION_AVAILABLE = False

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────
# 结构守卫 (v3.4.0) — 碎片判类不得违背 SMILES 结构
# ─────────────────────────────────────────────────────────────────
# 碎片相似度判类 (Layer 2) 易因质量重叠产生跨类误判 (例: 青藤碱碎片 137/153/165
# 与黄酮 B/A 环碎片重叠 → 被误判为黄酮)。当 SMILES 可用时, 对 *高特异类*
# (有可靠核心子结构的类) 做结构确认: 碎片推断的类必须在 SMILES 中找到对应核心,
# 否则拒绝 (退回 Unknown, 宁缺毋滥)。甾体/萜/酚酸/木脂素等无可靠通用 SMARTS 的类
# 不在此列, 仍信任标签/碎片。

# RDKit 懒加载 (避免无 RDKit 环境下 import 失败)
def _lazy_rdkit():
    try:
        from rdkit import Chem  # type: ignore
        return Chem
    except Exception:
        return None

# 高特异类 → 核心子结构 SMARTS (与 RuleBasedClassifier 同源, 仅保留经验证可靠者)
_CLASS_CORE_SMARTS: Dict[str, List[str]] = {
    "amine_oxide": ["[N+X4]([#6])([#6])([#6])[O-]"],
    "flavonoid": [
        "O=c1ccoc2ccccc12",                      # chromone (黄酮/黄酮醇/异黄酮)
        "O=C1CC(c2ccccc2)Oc2ccccc21",            # flavanone (黄烷酮)
        "O=C(/C=C/c1ccccc1)/C=C/c1ccccc1",       # chalcone (查耳酮)
    ],
    "alkaloid": ["[NX3;R]"],                       # 环内碱性氮 (兜底吗啡烷/哌啶/嘌呤等)
    "sugar": ["C1OC(O)C(O)C(O)C1O", "C1OC(O)C(O)C1O"],   # 吡喃糖/呋喃糖
    "amino_acid": ["N[CX4][CX3](=O)O"],           # alpha-氨基羧基
    "coumarin": ["c1ccc2oc(=O)ccc2c1"],           # 苯并吡喃-2-酮
    "quinone": ["O=C1c2ccccc2C(=O)c2ccccc21"],    # 9,10-蒽醌
    "nucleoside": ["[n;R]-[C;R]1[O;R][C;R][C;R][C;R]1"],
}


def _class_consistent_with_structure(classification: str, smiles: str) -> bool:
    """检查碎片推断的类别是否与 SMILES 结构一致。

    Returns:
        True  — 一致 (或无法判定: SMILES 缺失/非法/类无核心定义 → 宽松放行)
        False — SMILES 明确不含该类核心子结构 → 跨类矛盾, 应拒绝

    对 *高特异类* (在 _CLASS_CORE_SMARTS 中) 才做硬性确认; 其余类直接放行。
    """
    cores = _CLASS_CORE_SMARTS.get(classification)
    if not cores:
        return True  # 非高特异类: 无可靠 SMARTS, 信任碎片/标签
    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return True  # 无 SMILES: 无法判定, 宽松放行 (交碎片库负责)
    Chem = _lazy_rdkit()
    if Chem is None:
        return True  # 无 RDKit: 无法判定, 放行
    try:
        with __import__("warnings").catch_warnings():
            __import__("warnings").simplefilter("ignore")
            mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return False
    if mol is None:
        return False  # SMILES 非法: 不得作为结构守卫证据
    for patt in cores:
        p = Chem.MolFromSmarts(patt)
        if p is not None and mol.HasSubstructMatch(p):
            return True  # 找到该类核心 → 一致
    return False  # SMILES 合法但不含任何核心 → 跨类矛盾


_COMPATIBLE_CLASS_PAIRS = {
    frozenset({"glycoside", "flavonoid"}),
    frozenset({"glycoside", "phenolic_acid"}),
    frozenset({"glycoside", "coumarin"}),
    frozenset({"glycoside", "quinone"}),
    frozenset({"glycoside", "lignan"}),
    frozenset({"glycoside", "terpenoid"}),
    frozenset({"glycoside", "steroid"}),
    frozenset({"glycoside", "saponin"}),
    frozenset({"glycoside", "benzenoid"}),
    frozenset({"saponin", "terpenoid"}),
    frozenset({"saponin", "steroid"}),
    frozenset({"lipid", "fatty_acid"}),
    frozenset({"phenolic_acid", "benzenoid"}),
}

_EXPERT_DIMENSION_SPECS = (
    ("neutral_loss", "neutral_loss_score", "中性丢失", "Neutral Loss", 0.22),
    ("rule_coverage", "rule_coverage_score", "规则覆盖", "Rule Coverage", 0.22),
    ("fg_signature", "fg_signature_score", "官能团特征", "FG Signature", 0.18),
    ("feature_library", "feature_library_score", "大类特征库", "Class Feature Library", 0.15),
    ("structure_fragment", "structure_fragment_score", "结构-碎片", "Structure-Fragment", 0.10),
    ("diagnostic_ion", "diagnostic_ion_score", "诊断离子", "Diagnostic Ion", 0.08),
    ("bond_energy", "bond_energy_score", "键能支持", "Bond-Energy Support", 0.03),
    ("isomer_differentiation", "isomer_differentiation_score", "异构区分", "Isomer Differentiation", 0.02),
)


def _expert_float(value: Any, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def _classification_source_is_trusted(source: Optional[str]) -> tuple[bool, str]:
    """Return whether a stored class source can gate class-specific rules.

    ``None`` preserves compatibility for direct/legacy callers that supplied a
    class but had no provenance field.  Production search always supplies the
    actual database source, so ML/spectrum-only labels are not silently treated
    as confirmed structural classes.
    """
    return classification_source_trust(source, allow_legacy_missing=True)


def _classes_compatible(first: str, second: str) -> bool:
    if not first or not second or "Unknown" in (first, second):
        return True
    return first == second or frozenset({first, second}) in _COMPATIBLE_CLASS_PAIRS


def _dimension_entry(
    *,
    key: str,
    score_key: str,
    label: str,
    label_en: str,
    weight: float,
    score: float,
    status: str,
    reason: str,
    evidence_count: int = 0,
    expected_count: Optional[int] = None,
) -> Dict[str, Any]:
    entry: Dict[str, Any] = {
        "key": key,
        "score_key": score_key,
        "label": label,
        "label_en": label_en,
        "weight": weight,
        "score": round(max(0.0, min(1.0, score)), 4),
        "status": status,
        "reason": reason,
        "evidence_count": max(int(evidence_count or 0), 0),
    }
    if expected_count is not None:
        entry["expected_count"] = max(int(expected_count or 0), 0)
    return entry


def build_expert_dimension_assessments(details: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """Separate a real zero from an inapplicable/unavailable expert dimension."""
    base_evaluated = (
        "expert_total_score" in details
        or any(spec[1] in details for spec in _EXPERT_DIMENSION_SPECS)
    )
    evaluation_unavailable = details.get("evaluation_status") == "unavailable"
    if not base_evaluated or evaluation_unavailable:
        reason = (
            str(details.get("evaluation_reason") or "base_8d_not_evaluated")
            if evaluation_unavailable
            else "base_8d_not_evaluated"
        )
        return {
            spec[0]: _dimension_entry(
                key=spec[0],
                score_key=spec[1],
                label=spec[2],
                label_en=spec[3],
                weight=spec[4],
                score=_expert_float(details.get(spec[1])),
                status="unavailable",
                reason=reason,
            )
            for spec in _EXPERT_DIMENSION_SPECS
        }

    reliable_class = bool(details.get("classification_reliable"))
    knowledge = details.get("knowledge_applicability") or {}
    nl = details.get("neutral_loss_analysis") or {}
    di = details.get("diagnostic_ion_analysis") or {}
    rc = details.get("rule_coverage_analysis") or {}
    fg = details.get("fg_signature_analysis") or {}
    sf = details.get("structure_fragment_analysis") or {}
    be = details.get("bond_energy_analysis") or {}
    ff = details.get("feature_library_analysis") or {}
    iso = details.get("isomer_differentiation_analysis") or {}

    score_by_key = {
        spec[0]: _expert_float(details.get(spec[1]))
        for spec in _EXPERT_DIMENSION_SPECS
    }
    assessments: Dict[str, Dict[str, Any]] = {}

    def add(
        key: str,
        status: str,
        reason: str,
        evidence_count: int = 0,
        expected_count: Optional[int] = None,
    ) -> None:
        spec = next(item for item in _EXPERT_DIMENSION_SPECS if item[0] == key)
        assessments[key] = _dimension_entry(
            key=spec[0],
            score_key=spec[1],
            label=spec[2],
            label_en=spec[3],
            weight=spec[4],
            score=score_by_key[key],
            status=status,
            reason=reason,
            evidence_count=evidence_count,
            expected_count=expected_count,
        )

    nl_evidence = (
        int(nl.get("explained_count") or 0)
        + int(nl.get("sequential_explained") or 0)
        + int(nl.get("cascade_recovered") or 0)
    )
    if knowledge.get("neutral_loss_applicable") is False:
        add("neutral_loss", "not_applicable", "multiply_charged_precursor_without_product_charge")
    elif score_by_key["neutral_loss"] > 0:
        add("neutral_loss", "scored", "matched_curated_neutral_loss", nl_evidence)
    else:
        add("neutral_loss", "no_evidence", "no_curated_neutral_loss_match", nl_evidence)

    di_count = int(di.get("matched_count") or len(di.get("matched_ions") or []))
    if score_by_key["diagnostic_ion"] > 0:
        add("diagnostic_ion", "scored", "matched_scope_guarded_diagnostic_ion", di_count)
    elif not reliable_class:
        add("diagnostic_ion", "not_applicable", "classification_unconfirmed")
    else:
        add("diagnostic_ion", "no_evidence", "no_scope_guarded_diagnostic_ion_match")

    expected_rules = rc.get("expected_rules") or []
    excluded_rules = rc.get("excluded_rules") or {}
    if not reliable_class:
        add("rule_coverage", "not_applicable", "classification_unconfirmed")
    elif rc.get("note") == "unsupported_classification":
        add("rule_coverage", "not_applicable", "no_supported_rule_family_for_class")
    elif not expected_rules and excluded_rules:
        add(
            "rule_coverage",
            "not_applicable",
            "all_rules_outside_polarity_adduct_structure_or_energy_scope",
            expected_count=len(excluded_rules),
        )
    elif not expected_rules:
        add("rule_coverage", "unavailable", "no_executable_rules_for_class")
    elif score_by_key["rule_coverage"] > 0:
        add(
            "rule_coverage",
            "scored",
            "observed_losses_cover_executable_rules",
            len(rc.get("covered_rules") or []),
            len(expected_rules),
        )
    else:
        add(
            "rule_coverage",
            "no_evidence",
            "executable_rules_present_but_not_covered",
            0,
            len(expected_rules),
        )

    if not reliable_class:
        add("fg_signature", "not_applicable", "classification_unconfirmed")
    elif not fg.get("signature_found"):
        add("fg_signature", "unavailable", "no_functional_group_signature_for_class")
    elif score_by_key["fg_signature"] > 0:
        add(
            "fg_signature",
            "scored",
            "class_signature_partially_observed",
            len(fg.get("detected_losses") or []) + len(fg.get("detected_ions") or []),
        )
    else:
        add("fg_signature", "no_evidence", "class_signature_not_observed")

    sf_note = str(sf.get("consistency_note") or "")
    if sf_note in {"no_structure_or_no_fragments", "invalid_smiles", "rdkit_not_available"}:
        add("structure_fragment", "unavailable", sf_note or "structure_analysis_unavailable")
    elif score_by_key["structure_fragment"] > 0:
        add(
            "structure_fragment",
            "scored",
            "candidate_bond_fragments_match_observed_ions",
            int(sf.get("fragments_explained") or 0),
            int(sf.get("total_fragments_observed") or 0),
        )
    else:
        add(
            "structure_fragment",
            "no_evidence",
            "no_predicted_single_bond_fragment_match",
            0,
            int(sf.get("total_fragments_observed") or 0),
        )

    be_note = str(be.get("note") or "")
    if be_note == "no_structure":
        add("bond_energy", "unavailable", "candidate_structure_unavailable")
    elif score_by_key["bond_energy"] > 0:
        add(
            "bond_energy",
            "scored",
            "labile_bond_supported_by_observed_neutral_loss",
            int(be.get("explained_neutral_losses") or 0),
        )
    else:
        add("bond_energy", "no_evidence", be_note or "no_observed_loss_support_for_labile_bond")

    ff_note = str(ff.get("note") or "")
    ff_expected = int(ff.get("total_expected") or 0)
    ff_matches = int(ff.get("matched_count") or len(ff.get("matched_features") or []))
    if score_by_key["feature_library"] > 0:
        add("feature_library", "scored", "matched_class_feature_library", ff_matches, ff_expected)
    elif not reliable_class:
        add("feature_library", "not_applicable", "classification_unconfirmed")
    elif ff_note == "feature_library_unavailable":
        add("feature_library", "unavailable", ff_note)
    elif ff_note.startswith("no_profile_for_"):
        add("feature_library", "unavailable", "no_class_feature_profile")
    elif ff_expected <= 0:
        add("feature_library", "not_applicable", "no_scope_applicable_class_features")
    else:
        add("feature_library", "no_evidence", "class_features_present_but_not_observed", 0, ff_expected)

    iso_note = str(iso.get("note") or "")
    iso_evidence = int(iso.get("complete_pairs") or 0) + int(iso.get("partial_pairs") or 0)
    if score_by_key["isomer_differentiation"] > 0:
        add("isomer_differentiation", "scored", "diagnostic_isomer_ion_pair_observed", iso_evidence)
    elif not reliable_class:
        add("isomer_differentiation", "not_applicable", "classification_unconfirmed")
    elif iso_note in {"feature_library_unavailable", "no_isomer_rules_for_class"}:
        add("isomer_differentiation", "not_applicable", iso_note)
    else:
        add("isomer_differentiation", "no_evidence", "isomer_rules_present_but_no_diagnostic_pair")

    return assessments


def finalize_expert_evidence(
    details: Optional[Dict[str, Any]],
    fused_score: Optional[float] = None,
    classification_context: Optional[Dict[str, Any]] = None,
) -> Optional[Dict[str, Any]]:
    """Attach auditable 8-D applicability and max-fusion provenance.

    This helper is intentionally reusable for cached results: historical task
    JSON can be displayed truthfully without rerunning a large raw file.
    """
    if not isinstance(details, dict):
        return details

    if classification_context and not isinstance(details.get("classification_provenance"), dict):
        details["classification_provenance"] = {
            "input_classification": classification_context.get("classification") or "Unknown",
            "input_source": classification_context.get("classification_source") or "unknown",
            "input_confidence": classification_context.get("classification_confidence"),
            "effective_classification": details.get("effective_classification") or "Unknown",
            "analysis_class": details.get("analysis_class") or "Unknown",
            "reliable": bool(details.get("classification_reliable")),
            "reliability_reason": "legacy_result_recorded_at_analysis_time",
        }

    dimensions = build_expert_dimension_assessments(details)
    details["dimension_assessments"] = dimensions

    base_score = _expert_float(details.get("expert_total_score"))
    v3_payload = details.get("v3")
    pathway_payload = details.get("pathway_align")
    components: Dict[str, Dict[str, Any]] = {
        "base_8d": {
            "label": "基础 8 维",
            "label_en": "Base 8 dimensions",
            "score": round(base_score, 4),
            "available": "expert_total_score" in details,
            "independent_source": "curated_class_and_structure_rules",
            "independent_evidence": True,
            "eligible_for_fusion": "expert_total_score" in details,
        },
        "v3_kb": {
            "label": "v3 规则库",
            "label_en": "v3 knowledge base",
            "score": round(
                _expert_float(v3_payload.get("v3_score"))
                if isinstance(v3_payload, dict) else 0.0,
                4,
            ),
            "available": isinstance(v3_payload, dict),
            "independent_source": "numeric_fragmentation_rules_v3",
            "independent_evidence": True,
            "eligible_for_fusion": bool(
                isinstance(v3_payload, dict)
                and v3_payload.get("classification_reliable") is True
            ),
        },
        "compound_pathway": {
            "label": "候选专属通路",
            "label_en": "Compound-specific pathway",
            "score": round(
                _expert_float(
                    pathway_payload.get("score", pathway_payload.get("alignment_score"))
                )
                if isinstance(pathway_payload, dict) else 0.0,
                4,
            ),
            "available": isinstance(pathway_payload, dict),
            "independent_source": "candidate_conditioned_library_spectral_summary",
            "independent_evidence": False,
            "eligible_for_fusion": False,
            "suppression_reason": (
                "query_exclusion_and_source_independence_not_proven"
            ),
        },
    }
    eligible_components = [
        (key, value["score"])
        for key, value in components.items()
        if value["available"] and value["eligible_for_fusion"]
    ]
    component_max = max((score for _, score in eligible_components), default=0.0)
    # Never trust a caller-supplied fused score here: historical callers could
    # have max-merged candidate-conditioned pathway evidence before invoking
    # this helper.  Reconstruct the production score solely from explicitly
    # eligible components so cached and live results fail closed alike.
    final_score = component_max
    selected_sources = [
        key for key, score in eligible_components
        if abs(score - final_score) <= 1e-6
    ]
    suppressed_sources = [
        key for key, value in components.items()
        if value["available"] and not value["eligible_for_fusion"]
    ]

    details["fused_expert_score"] = round(final_score, 4)
    details["score_fusion"] = {
        "method": "max_of_eligible_expert_sources",
        "fused_score": round(final_score, 4),
        "selected_source": selected_sources[0] if selected_sources else "none",
        "selected_sources": selected_sources,
        "suppressed_sources": suppressed_sources,
        "caller_supplied_fused_score": (
            round(_expert_float(fused_score), 4)
            if fused_score is not None else None
        ),
        "components": components,
        "base_dimensions_all_zero": all(
            assessment["score"] == 0 for assessment in dimensions.values()
        ),
        "base_dimensions_not_applicable_or_unavailable": sum(
            assessment["status"] in {"not_applicable", "unavailable"}
            for assessment in dimensions.values()
        ),
        "version": "expert-evidence-provenance-v2",
    }
    return details


# ─────────────────────────────────────────────────────────────────
# 专家评分器核心类
# ─────────────────────────────────────────────────────────────────


class ChemicalExpertScorer:
    """
    化学专家评分器 — 从质谱裂解原理角度评估匹配质量

    Usage:
        scorer = ChemicalExpertScorer()
        expert_scores = scorer.evaluate(
            query_peaks, query_pmz, candidate_smiles, candidate_formula,
            candidate_classification
        )
    """

    def __init__(
        self,
        neutral_loss_tolerance: float = 0.02,
        diagnostic_ion_tolerance: float = 0.02,
        min_fragments_for_analysis: int = 3,
    ):
        """
        Args:
            neutral_loss_tolerance: 中性丢失质量匹配容差(Da)
            diagnostic_ion_tolerance: 诊断离子质量匹配容差(Da)
            min_fragments_for_analysis: 最少碎片数才进行详细分析
        """
        self.nl_tol = neutral_loss_tolerance
        self.di_tol = diagnostic_ion_tolerance
        self.min_frags = min_fragments_for_analysis

    # ── 主评估接口 ───────────────────────────────────────────────

    def evaluate(
        self,
        query_peaks: List[List[float]],
        query_pmz: float,
        candidate_smiles: str,
        candidate_formula: str = "",
        candidate_classification: str = "Unknown",
        candidate_inchikey: str = "",
        query_polarity: str = "",
        query_adduct: str = "",
        query_charge: Optional[int] = None,
        query_collision_energy: Any = None,
        candidate_adduct: str = "",
        candidate_classification_source: Optional[str] = None,
        candidate_classification_confidence: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        执行完整的专家评估

        Args:
            query_peaks: 查询峰列表 [[mz, intensity], ...]
            query_pmz: 查询前体离子 m/z
            candidate_smiles: 候选化合物 SMILES
            candidate_formula: 候选化合物分子式
            candidate_classification: 候选化合物分类
            candidate_inchikey: 候选化合物 InChIKey
            query_polarity: 查询谱图极性(POS/NEG/positive/negative)
            query_adduct: 查询谱图加合物
            query_charge: 查询前体电荷数；多电荷时禁用不安全的单电荷NL推断
            query_collision_energy: 碰撞能量元数据(用于适用性记录)
            candidate_adduct: 馆藏候选加合物；仅在查询加合物缺失且极性一致时作为
                明确标注的候选假设，不冒充独立实测加合物
            candidate_classification_source: 候选分类来源(curated/smarts/ml/...)
            candidate_classification_confidence: 候选分类来源给出的置信度

        Returns:
            评估结果字典，包含各维度得分和详细分析
        """
        if query_pmz <= 0:
            return self._empty_result("invalid_precursor_mz")

        # 碎片太少时仍进行分析，但标记为低精度模式
        low_precision_mode = len(query_peaks) < self.min_frags
        if not query_peaks:
            return self._empty_result("no_peaks")

        # 提取碎片 m/z 列表
        fragment_mzs = [p[0] for p in query_peaks if p[0] < query_pmz - 1.5]
        effective_polarity = normalize_polarity(query_polarity)
        if effective_polarity is None:
            adduct_text = str(query_adduct or "").strip()
            if adduct_text.endswith("+"):
                effective_polarity = "positive"
            elif adduct_text.endswith("-"):
                effective_polarity = "negative"
        if effective_polarity is None and query_charge is not None:
            try:
                effective_polarity = "positive" if int(query_charge) > 0 else (
                    "negative" if int(query_charge) < 0 else None
                )
            except (TypeError, ValueError):
                pass

        explicit_query_adduct = str(query_adduct or "").strip()
        candidate_adduct_text = str(candidate_adduct or "").strip()
        effective_adduct = explicit_query_adduct
        adduct_source = "query_metadata" if explicit_query_adduct else "unavailable"
        if not effective_adduct and candidate_adduct_text:
            candidate_adduct_polarity = (
                "positive" if candidate_adduct_text.endswith("+")
                else "negative" if candidate_adduct_text.endswith("-")
                else None
            )
            if (
                effective_polarity is None
                or candidate_adduct_polarity is None
                or candidate_adduct_polarity == effective_polarity
            ):
                effective_adduct = candidate_adduct_text
                adduct_source = "candidate_hypothesis"
                if effective_polarity is None:
                    effective_polarity = candidate_adduct_polarity

        try:
            effective_charge = max(abs(int(query_charge)), 1) if query_charge is not None else 1
        except (TypeError, ValueError):
            effective_charge = 1
        collision_energy_value, collision_energy_unit = collision_energy_value_unit(
            query_collision_energy
        )

        # ── 自动推断化合物分类 + 来源可靠性审计 ──
        # 不能再把所有非 Unknown 标签一概视为可靠：Full 库中相当一部分标签来自
        # 未校准 ML。只有策展/结构来源，或 ML 与独立 SMILES 规则一致时，才能启用
        # “按大类查规则”的维度。
        input_classification = normalize_classification(candidate_classification)
        effective_classification = input_classification
        source_trusted, source_reason = _classification_source_is_trusted(
            candidate_classification_source
        )
        structure_class_raw = "Unknown"
        structure_class = "Unknown"
        fragment_class = "Unknown"
        fragment_class_confidence = 0.0
        cls_reliable = False
        reliability_reason = source_reason

        if candidate_smiles:
            try:
                from mssystem.classification import RuleBasedClassifier

                structure_class_raw = RuleBasedClassifier().predict(
                    smiles=candidate_smiles
                )
                structure_class = normalize_classification(structure_class_raw)
            except Exception:
                structure_class_raw = "Unknown"
                structure_class = "Unknown"

        if structure_class != "Unknown":
            if input_classification == "Unknown":
                effective_classification = structure_class
                cls_reliable = True
                reliability_reason = "smiles_structure_rule"
            elif structure_class == input_classification:
                effective_classification = input_classification
                cls_reliable = True
                reliability_reason = (
                    "stored_class_agrees_with_smiles_structure_rule"
                )
            elif not source_trusted:
                effective_classification = structure_class
                cls_reliable = True
                reliability_reason = (
                    "smiles_structure_rule_overrode_provisional_stored_class"
                )
            elif candidate_classification_source is None:
                # Direct legacy callers intentionally provide a class without
                # source metadata. Preserve their historical contract.
                effective_classification = input_classification
                cls_reliable = True
                reliability_reason = "legacy_explicit_class_without_source_metadata"
            elif _classes_compatible(input_classification, structure_class):
                effective_classification = input_classification
                cls_reliable = True
                reliability_reason = "stored_and_structure_classes_compatible"
            else:
                # A trusted label and the supplied structure disagree.  Do not
                # choose a side for class-specific scoring; retain the label for
                # traceability and fail closed on class-gated rules.
                effective_classification = input_classification
                cls_reliable = False
                reliability_reason = "stored_class_conflicts_with_smiles_structure"
        elif input_classification != "Unknown":
            effective_classification = input_classification
            cls_reliable = source_trusted
            reliability_reason = source_reason

        # Final fallback: fragments may provide a descriptive class hypothesis,
        # but never unlock class-specific scoring on their own.
        if effective_classification == "Unknown" and FEATURE_LIB_AVAILABLE and fragment_mzs:
            try:
                cls_result = classify_compound_from_fragments(
                    fragment_mzs, polarity=effective_polarity or ""
                )
                if cls_result and cls_result.get("class_name"):
                    inferred = normalize_classification(
                        str(cls_result.get("class_name") or "")
                    )
                    inferred_confidence = _expert_float(
                        cls_result.get("confidence")
                    )
                    if inferred != "Unknown" and inferred_confidence >= 0.5:
                        if candidate_smiles and not _class_consistent_with_structure(
                            inferred, candidate_smiles
                        ):
                            logger.debug(
                                "Fragment class %s rejected by candidate structure",
                                inferred,
                            )
                        else:
                            fragment_class = inferred
                            fragment_class_confidence = inferred_confidence
                            effective_classification = inferred
                            cls_reliable = False
                            reliability_reason = (
                                "fragment_inference_only_class_rules_skipped"
                            )
            except Exception:
                pass

        analysis_class = effective_classification if cls_reliable else "Unknown"
        try:
            input_confidence = (
                float(candidate_classification_confidence)
                if candidate_classification_confidence is not None
                else None
            )
        except (TypeError, ValueError):
            input_confidence = None
        classification_provenance = {
            "input_classification": str(candidate_classification or "Unknown"),
            "normalized_input_classification": input_classification,
            "input_source": (
                str(candidate_classification_source)
                if candidate_classification_source is not None
                else "legacy_unspecified"
            ),
            "input_confidence": input_confidence,
            "structure_class_raw": structure_class_raw,
            "structure_class": structure_class,
            "fragment_class": fragment_class,
            "fragment_class_confidence": round(fragment_class_confidence, 4),
            "effective_classification": effective_classification,
            "analysis_class": analysis_class,
            "reliable": cls_reliable,
            "reliability_reason": reliability_reason,
        }

        # 1. 中性丢失分析
        nl_analysis = self._analyze_neutral_losses(
            query_pmz, fragment_mzs, classification=analysis_class,
            query_charge=effective_charge,
            polarity=effective_polarity or "",
            candidate_smiles=candidate_smiles,
            query_adduct=effective_adduct,
            query_collision_energy=query_collision_energy,
        )

        # 2. 诊断离子检测
        di_analysis = self._analyze_diagnostic_ions(
            fragment_mzs, classification=analysis_class,
            polarity=effective_polarity or "",
            candidate_smiles=candidate_smiles,
            raw_classification=candidate_classification,
            query_adduct=effective_adduct,
            query_collision_energy=query_collision_energy,
        )

        # 3. 裂解规则覆盖度(传入完整NL分析以包含顺序丢失)
        rc_analysis = self._analyze_rule_coverage(
            nl_analysis,
            analysis_class,
            polarity=effective_polarity or "",
            candidate_smiles=candidate_smiles,
            query_adduct=effective_adduct,
            query_collision_energy=query_collision_energy,
        )

        # 4. 官能团特征验证
        fg_analysis = self._analyze_functional_group_signature(
            nl_analysis["neutral_losses"],
            di_analysis["matched_ions"],
            analysis_class,
            polarity=effective_polarity or "",
        )

        # 5. 结构-碎片一致性(需要RDKit)
        sf_analysis = self._analyze_structure_fragment_consistency(
            candidate_smiles, fragment_mzs, query_pmz,
            classification=analysis_class,
            polarity=effective_polarity or "",
            query_adduct=effective_adduct,
            query_charge=effective_charge,
        )

        # 6. 键能合理性评估
        be_analysis = self._analyze_bond_energy_plausibility(
            candidate_smiles, nl_analysis["neutral_losses"]
        )

        # 7. v3.0: 碎片特征库验证 (Fragment Feature Library Verification)
        ff_analysis = self._analyze_feature_library(
            fragment_mzs, analysis_class, polarity=effective_polarity or "",
            candidate_smiles=candidate_smiles,
            raw_classification=candidate_classification,
            query_adduct=effective_adduct,
            query_collision_energy=query_collision_energy,
        )

        # 8. v3.0: 异构体诊断离子对检测 (Isomer Differentiation)
        id_analysis = self._analyze_isomer_differentiation(
            fragment_mzs, analysis_class, polarity=effective_polarity or ""
        )

        # 9. 计算综合专家得分 (v3.0 权重)
        expert_scores = self._compute_expert_scores(
            nl_analysis, di_analysis, rc_analysis,
            fg_analysis, sf_analysis, be_analysis,
            ff_analysis, id_analysis,
        )

        # 低精度模式：总分打折
        if low_precision_mode:
            expert_scores["total"] = round(expert_scores["total"] * 0.6, 4)

        result = {
            "expert_total_score": expert_scores["total"],
            "neutral_loss_score": expert_scores["neutral_loss"],
            "diagnostic_ion_score": expert_scores["diagnostic_ion"],
            "rule_coverage_score": expert_scores["rule_coverage"],
            "fg_signature_score": expert_scores["fg_signature"],
            "structure_fragment_score": expert_scores["structure_fragment"],
            "bond_energy_score": expert_scores["bond_energy"],
            "feature_library_score": expert_scores.get("feature_library", 0.0),
            "isomer_differentiation_score": expert_scores.get("isomer_differentiation", 0.0),
            "confidence_level": self._expert_confidence(expert_scores["total"]),
            "effective_classification": effective_classification,
            "analysis_class": analysis_class,
            "classification_reliable": cls_reliable,
            "classification_provenance": classification_provenance,
            "knowledge_applicability": {
                "polarity": effective_polarity or "unknown",
                "adduct": effective_adduct,
                "adduct_source": adduct_source,
                "query_adduct": explicit_query_adduct,
                "candidate_adduct": candidate_adduct_text,
                "charge": effective_charge,
                "collision_energy": query_collision_energy,
                "collision_energy_value": collision_energy_value,
                "collision_energy_unit": collision_energy_unit or "unknown",
                "neutral_loss_applicable": effective_charge == 1,
                "warnings": [
                    warning for warning in (
                        "polarity_unknown_rules_not_filtered" if effective_polarity is None else "",
                        "multiply_charged_precursor_neutral_loss_disabled" if effective_charge > 1 else "",
                        "collision_energy_unknown_energy_rules_not_activated"
                        if collision_energy_value is None else "",
                        "collision_energy_unit_unknown_conditioned_rules_not_activated"
                        if collision_energy_value is not None and not collision_energy_unit else "",
                        "candidate_adduct_used_as_hypothesis"
                        if adduct_source == "candidate_hypothesis" else "",
                        "classification_unconfirmed_class_specific_rules_skipped"
                        if not cls_reliable else "",
                    ) if warning
                ],
            },
            "neutral_loss_analysis": nl_analysis,
            "diagnostic_ion_analysis": di_analysis,
            "rule_coverage_analysis": rc_analysis,
            "fg_signature_analysis": fg_analysis,
            "structure_fragment_analysis": sf_analysis,
            "bond_energy_analysis": be_analysis,
            "feature_library_analysis": ff_analysis,
            "isomer_differentiation_analysis": id_analysis,
            "chemical_reasoning": self._generate_reasoning(
                nl_analysis, di_analysis, rc_analysis, fg_analysis, sf_analysis,
                ff_analysis, id_analysis,
                effective_classification=analysis_class,
                cls_reliable=cls_reliable,
            ),
        }
        return finalize_expert_evidence(
            result,
            fused_score=expert_scores["total"],
        )

    # ── 1. 中性丢失分析 ──────────────────────────────────────────

    def _analyze_neutral_losses(
        self, precursor_mz: float, fragment_mzs: List[float],
        classification: str = "",
        query_charge: int = 1,
        polarity: str = "",
        candidate_smiles: str = "",
        query_adduct: str = "",
        query_collision_energy: Any = None,
    ) -> Dict[str, Any]:
        """分析前体→碎片 AND 碎片→碎片之间的质量差

        两阶段分析:
        1. precursor→fragment: 传统中性丢失分析
        2. fragment→fragment: 顺序裂解分析(如 CO→HCN 序列)

        分类过滤: 使用 classification 参数过滤与化合物类型不兼容的中性丢失。
        """
        if not fragment_mzs or precursor_mz <= 0:
            return {"neutral_losses": [], "sequential_losses": [],
                    "explained_count": 0, "unexplained_count": 0,
                    "sequential_explained": 0, "explanation_rate": 0.0}
        if abs(int(query_charge or 1)) > 1:
            return {
                "neutral_losses": [], "sequential_losses": [],
                "connected_pathways": [], "best_pathway": [],
                "explained_count": 0, "unexplained_count": len(fragment_mzs),
                "sequential_explained": 0, "cascade_recovered": 0,
                "explanation_rate": 0.0, "priority_weighted_score": 0.0,
                "total_fragment_pairs": len(fragment_mzs),
                "note": "multiply_charged_precursor_neutral_loss_not_assigned_without_product_charge",
            }

        # ── 阶段1: precursor→fragment 中性丢失 ──
        nl_pairs = compute_possible_neutral_losses(
            precursor_mz, fragment_mzs, self.nl_tol,
            classification=classification, strict_fg=True,
            polarity=polarity, candidate_smiles=candidate_smiles,
            query_adduct=query_adduct,
            query_collision_energy=query_collision_energy,
        )

        explained = []
        unexplained = []
        for diff, nl in nl_pairs:
            if nl:
                explained.append({
                    "type": "precursor_to_fragment",
                    "mass_difference": round(diff, 4),
                    "neutral_loss": nl.name,
                    "neutral_loss_en": nl.name_en,
                    "formula": nl.formula,
                    "mechanism": nl.mechanism,
                    "priority": nl.priority,
                    "functional_groups": sorted(nl.functional_groups),
                })
            else:
                composites = get_composite_neutral_loss_by_mass(
                    diff,
                    tolerance=min(self.nl_tol, 0.05),
                    classification=classification,
                    strict_fg=True,
                )
                if composites:
                    composite = composites[0]
                    explained.append({
                        "type": "precursor_to_fragment",
                        "mass_difference": round(diff, 4),
                        "neutral_loss": f"累积丢失：{composite['name']}",
                        "neutral_loss_en": f"Cumulative loss: {composite['name_en']}",
                        "formula": composite["formula"],
                        "mechanism": composite["mechanism"],
                        "priority": composite["priority"],
                        "functional_groups": [],
                        "assignment_level": "tentative_composite",
                        "mass_delta": composite["delta"],
                    })
                    continue
                unexplained.append({
                    "type": "precursor_to_fragment",
                    "mass_difference": round(diff, 4),
                    "assignment_status": "unassigned_fragment_difference",
                    "possible_candidates": [
                        {"name": x.name, "mass": round(x.mass, 4), "formula": x.formula}
                        for x in get_neutral_loss_by_mass(diff, self.nl_tol * 2,
                                                          classification=classification,
                                                          strict_fg=True,
                                                          polarity=polarity,
                                                          candidate_smiles=candidate_smiles,
                                                          query_adduct=query_adduct,
                                                          query_collision_energy=query_collision_energy)[:3]
                    ],
                })

        # ── 阶段2: fragment→fragment 连通碎裂链分析 ──
        # 构建真实连通链 (A→B→C, 每步以上一步输出为起点)
        # 替代旧方案: 旧方案将所有无关碎片对拼接成伪路径
        sorted_frags = sorted(fragment_mzs, reverse=True)  # 从大到小排列
        pathways = self._build_connected_pathways(
            sorted_frags, classification,
            polarity=polarity,
            candidate_smiles=candidate_smiles,
            query_adduct=query_adduct,
            query_collision_energy=query_collision_energy,
        )

        # 向后兼容: 从最佳路径填充 sequential 列表
        best_pathway = pathways[0] if pathways else []
        sequential = list(best_pathway)  # 每个元素已是 dict
        seq_explained = len(best_pathway)

        # ── 阶段3: 类别专属级联 (从前体起链, 弥补阶段2从最大碎片起链的缺口) ──
        # 糖苷/皂苷: 糖基顺序丢失 [M+H]→逐级脱糖→苷元
        # 乌头碱型生物碱: 酯基顺序丢失 DDA→MDA→乌头原碱 (毒性构型判别)
        # 注: 这些丢失在阶段1/2中可能已被独立计为单步, 此处仅作"连通级联"
        # 的结构化补充, 不重复计入 explanation_rate, 避免双重计数。
        sugar_cascade = self._detect_sugar_loss_cascade(
            precursor_mz, fragment_mzs, classification
        )
        ester_cascade = self._detect_alkaloid_ester_loss_cascade(
            precursor_mz, fragment_mzs, classification
        )

        # ── 综合统计 ──
        # 阶段3级联可能解释了阶段1/2无法识别的丢失。现 AcOH/PhCOOH/CH3OH/H2O 均在
        # 通用库中, 阶段1对匹配类(生物碱/萜类/酯类)已可独立检出; 级联此处仍作
        # "连通级联"的结构化补充(从前体起链的完整路径)。把级联中"未被阶段1/2覆盖"
        # 的步数补回 explained_count, 修复盲点 ( provably 不双重计数)。
        explained_diffs = {round(e["mass_difference"], 2) for e in explained}
        seq_diffs = {round(s["mass_difference"], 2) for s in sequential if "mass_difference" in s}
        cascade_recovered = 0
        for cascade_obj in (sugar_cascade, ester_cascade):
            if not cascade_obj:
                continue
            for step in cascade_obj.get("cascade", []):
                d = round(step.get("mass_difference", 0), 2)
                if d > 0 and d not in explained_diffs and d not in seq_diffs:
                    cascade_recovered += 1

        total_pairs = max(len(nl_pairs) + len(sequential) + cascade_recovered, 1)
        total_explained = len(explained) + seq_explained + cascade_recovered
        explanation_rate = total_explained / total_pairs

        # 优先级加权得分(阶段1+阶段2)
        weighted_score = 0.0
        all_explained_priorities = [e["priority"] for e in explained]
        all_explained_priorities += [s["priority"] for s in sequential if "priority" in s]
        if all_explained_priorities:
            max_priority = max(all_explained_priorities)
            weighted_score = sum(p / max_priority for p in all_explained_priorities)
            weighted_score = min(1.0, weighted_score / max(len(all_explained_priorities), 1))

        return {
            "neutral_losses": explained + unexplained,      # 阶段1
            "sequential_losses": sequential,                 # 阶段2 (向后兼容)
            "connected_pathways": pathways,                  # 新增: 所有连通链
            "best_pathway": best_pathway,                    # 新增: 最优链
            "sugar_loss_cascade": sugar_cascade,             # 阶段3a: 糖苷糖基级联
            "alkaloid_ester_cascade": ester_cascade,         # 阶段3b: 乌头碱酯基级联
            "explained_count": len(explained) + cascade_recovered,
            "unexplained_count": len(unexplained),
            "sequential_explained": seq_explained,
            "cascade_recovered": cascade_recovered,
            "explanation_rate": explanation_rate,
            "priority_weighted_score": weighted_score,
            "total_fragment_pairs": total_pairs,
        }

    def _build_connected_pathways(
        self,
        sorted_frags: List[float],
        classification: str = "",
        max_depth: int = 5,
        max_pathways: int = 10,
        polarity: str = "",
        candidate_smiles: str = "",
        query_adduct: str = "",
        query_collision_energy: Any = None,
    ) -> List[List[Dict[str, Any]]]:
        """构建真实连通碎裂链

        与旧方案(所有碎片对扁平列表)不同，本方法构建的链满足:
        A → B → C 其中 B 是 A 丢失 NL1 后的碎片，C 是 B 丢失 NL2 后的碎片。

        算法:
        1. 构建邻接表: 对每个碎片 A，找到所有比它小的碎片 B，
           使得 mass(A) - mass(B) 匹配已知中性丢失
        2. DFS 搜索链: 从每个碎片出发深度优先搜索，收集长度 ≥ 2 的链
        3. 按链长度 + 优先级排序

        Args:
            sorted_frags: 从大到小排序的碎片 m/z 列表
            classification: 化合物分类(用于过滤中性丢失)
            max_depth: 最大搜索深度
            max_pathways: 最多返回的路径数

        Returns:
            连通路径列表，每条路径是步骤字典的列表
        """
        if len(sorted_frags) < 2:
            return []

        # ── 1. 构建邻接表 ──
        # children[parent_mz] = [(child_mz, nl_match)]
        children: Dict[float, List[Tuple[float, NeutralLoss]]] = {}

        for i, parent_mz in enumerate(sorted_frags):
            children[parent_mz] = []
            for j in range(i + 1, len(sorted_frags)):
                child_mz = sorted_frags[j]
                diff = parent_mz - child_mz
                if diff <= 0:
                    continue
                matches = get_neutral_loss_by_mass(
                    diff, self.nl_tol,
                    classification=classification,
                    strict_fg=True,
                    polarity=polarity,
                    candidate_smiles=candidate_smiles,
                    query_adduct=query_adduct,
                    query_collision_energy=query_collision_energy,
                )
                if matches:
                    children[parent_mz].append((child_mz, matches[0]))

        # ── 2. DFS 搜索连通链 ──
        pathways: List[List[Dict[str, Any]]] = []

        def dfs(node: float, path: List[Dict[str, Any]], visited: set):
            if len(path) >= 1:
                # 记录有效路径 (至少 1 步 = 2 个碎片)
                pathways.append(list(path))

            if len(path) >= max_depth:
                return

            for child_mz, nl_match in children.get(node, []):
                if child_mz in visited:
                    continue
                step = {
                    "type": "fragment_to_fragment",
                    "from_mz": round(node, 4),
                    "to_mz": round(child_mz, 4),
                    "mass_difference": round(node - child_mz, 4),
                    "neutral_loss": nl_match.name,
                    "neutral_loss_en": nl_match.name_en,
                    "formula": nl_match.formula,
                    "mechanism": nl_match.mechanism,
                    "priority": nl_match.priority,
                }
                path.append(step)
                visited.add(child_mz)
                dfs(child_mz, path, visited)
                path.pop()
                visited.discard(child_mz)

        for start_mz in sorted_frags:
            dfs(start_mz, [], {start_mz})

        # ── 3. 排序: 链长度优先，然后按优先级之和 ──
        pathways.sort(key=lambda p: (
            len(p),
            sum(s.get("priority", 0) for s in p),
        ), reverse=True)

        return pathways[:max_pathways]

    # ── 1b. 糖苷糖基顺序丢失级联 (从前体起链) ───────────────────

    def _detect_sugar_loss_cascade(
        self,
        precursor_mz: float,
        fragment_mzs: List[float],
        classification: str = "",
    ) -> Optional[Dict[str, Any]]:
        """糖苷/皂苷糖基顺序丢失级联检测

        从前体离子出发, 贪心链接连续的糖基中性丢失, 复现糖苷类的特征碎裂路径:
            [M+H]⁺ → [M+H−糖]⁺ → [M+H−2糖]⁺ → … → 苷元 (aglycone)

        与 _build_connected_pathways (fragment→fragment, 从最大碎片起链) 互补:
        本方法从前体起链, 从而捕获最关键的"首步脱糖", 给出完整糖链测序。

        化学依据:
          - Domon & Costello 糖链裂解命名 (Y₀/B₀ 离子)
          - 《中药质谱分析》刘淑莹等 (2012): 糖苷类糖链顺序丢失规律
          - 糖基丢失质量: 己糖162 / 脱氧己糖146 / 戊糖132 / 葡萄糖醛酸176

        糖基精确质量复用 COMMON_NEUTRAL_LOSSES (C6H10O5/C6H10O4/C5H8O4/C6H8O6),
        保留各自糖型语义, 不引入新的化学常数。

        Returns:
            级联描述 dict (cascade/n_sugars/aglycone_mz/sugar_composition),
            无级联或证据不足时返回 None。
        """
        if not fragment_mzs or precursor_mz <= 0:
            return None

        # ── 糖基丢失参考: 从 COMMON_NEUTRAL_LOSSES 取, 保留糖型语义 ──
        # (formula key → 糖型代号, 中文名, 缺省精确质量[库中无该 key 时回退])
        sugar_fallback = {
            "C6H10O5": ("hexose", "己糖", 162.0528),
            "C6H10O4": ("deoxyhexose", "脱氧己糖", 146.0579),
            "C5H8O4":  ("pentose", "戊糖", 132.0423),
            "C6H8O6":  ("glucuronic_acid", "葡萄糖醛酸", 176.0321),
        }
        sugar_refs = []
        for key, (sugar_type, sugar_cn, default_mass) in sugar_fallback.items():
            nl = COMMON_NEUTRAL_LOSSES.get(key)
            mass = nl.mass if nl else default_mass
            sugar_refs.append((sugar_type, sugar_cn, mass))

        # 级联连接允许 0.06 Da：输入文件常把高 m/z 碎片只保留 1 位小数，
        # 162.0/162.1 均应能连接到 162.0528；通用 HRMS 规则仍使用 0.05 Da。
        tol = max(self.nl_tol, 0.06)
        # 特征性含糖苷键的类别: 单步糖丢失即具诊断意义, 可直接报告;
        # 其余类 (酚酸/香豆素/醌/木脂素/生物碱/苯环类...) 单次 162 丢失可能为
        # 巧合或非糖中性丢失, 需 ≥2 步连续糖丢失作强证据, 避免误报。
        glyco_classes = {
            "saponin", "glycoside", "steroid", "flavonoid",
            "triterpenoid", "terpenoid",
        }
        cls = (classification or "").lower()
        is_glyco_class = any(g in cls for g in glyco_classes)

        # ── 从前体起链, 贪心匹配连续糖基丢失 ──
        cascade: List[Dict[str, Any]] = []
        used: set = set()
        current = float(precursor_mz)
        total_sugar_mass = 0.0
        frags_sorted = sorted(set(fragment_mzs), reverse=True)

        # 大皂苷(如人参皂苷 Rb1)含 4-5 个糖, 预留 8 步余量
        for _ in range(8):
            best = None  # (frag_mz, sugar_type, sugar_cn, mass_diff, err)
            for frag in frags_sorted:
                if frag in used or frag >= current - 1.0:
                    continue
                diff = current - frag
                if diff <= 0:
                    continue
                for sugar_type, sugar_cn, sugar_mass in sugar_refs:
                    err = abs(diff - sugar_mass)
                    if err <= tol and (best is None or err < best[4]):
                        best = (frag, sugar_type, sugar_cn, diff, err)
            if best is None:
                break
            frag, sugar_type, sugar_cn, diff, _ = best
            cascade.append({
                "from_mz": round(current, 4),
                "to_mz": round(frag, 4),
                "mass_difference": round(diff, 4),
                "sugar_type": sugar_type,
                "sugar_name_cn": sugar_cn,
            })
            total_sugar_mass += diff
            used.add(frag)
            current = frag

        n_sugars = len(cascade)

        # ── 报告门槛 (宁缺毋滥, 避免巧合的单次 162 丢失误报) ──
        if n_sugars == 0:
            return None
        if n_sugars < 2 and not is_glyco_class:
            return None

        # 糖组成统计, e.g. {"hexose": 2, "deoxyhexose": 1}
        composition: Dict[str, int] = {}
        for step in cascade:
            composition[step["sugar_type"]] = composition.get(step["sugar_type"], 0) + 1

        return {
            "cascade": cascade,
            "n_sugars": n_sugars,
            "total_sugar_mass": round(total_sugar_mass, 4),
            "aglycone_mz": round(current, 4),  # 末碎片 = 苷元/去糖骨架离子
            "sugar_composition": composition,
            "is_glycoside_evidence": True,
        }

    # ── 1c. 乌头碱型生物碱酯基顺序丢失级联 ─────────────────────

    def _detect_alkaloid_ester_loss_cascade(
        self,
        precursor_mz: float,
        fragment_mzs: List[float],
        classification: str = "",
    ) -> Optional[Dict[str, Any]]:
        """乌头碱型 C₁₉二萜生物碱 酯基顺序丢失级联检测

        乌头碱骨架(附子/川乌/草乌, 治类风湿寒湿痹证核心药)的酯取代基会产生
        特征中性丢失, 据此可区分毒性构型:

            双酯型 DDA (剧毒) ─脱乙酰基─▶ 单酯型 MDA (低毒) ─脱苯甲酰基─▶ 乌头原碱

        酯取代基 → 特征丢失:
            乙酰基   → −60   AcOH     (乙酸, C2H4O2)
            苯甲酰基 → −122  PhCOOH   (苯甲酸, C7H6O2)
            甲氧基   → −32   CH3OH    (甲醇)
            羟基     → −18   H2O      (水)

        化学依据: 《中药质谱分析》刘淑莹等 (2012) 第5章 C19二萜生物碱碎裂规律。
        四种酯基丢失 (AcOH/PhCOOH/CH3OH/H2O) 的精确质量统一取自
        COMMON_NEUTRAL_LOSSES (单一真值源, 与通用中性丢失检测对齐);
        此处保留乌头碱毒性构型判别用的叙述性名称 (脱乙酰/脱苯甲酰 = 去酯基)。
        注: PhCOOH (苯甲酸) 已提升进通用库 (含 alkaloid/terpenoid 官能团标签),
        故阶段1对生物碱类即可独立检出该丢失。

        Returns:
            级联描述 dict (cascade/n_losses/aglycone_mz/loss_composition/
            aconitine_subtype), 无级联或证据不足时返回 None。
        """
        if not fragment_mzs or precursor_mz <= 0:
            return None

        # 仅生物碱类触发 (乌头碱型属二萜生物碱)
        cls = (classification or "").lower()
        if "alkaloid" not in cls and "diterpene" not in cls:
            return None

        # ── 酯基丢失参考 (乌头碱型特征) ──
        # 精确质量统一从 COMMON_NEUTRAL_LOSSES 取 (单一真值源, 库无该 key 时回退);
        # 保留乌头碱毒性构型判别用的叙述性名称 (乙酰基/苯甲酰基 = 去酯基框架)。
        def _lib_mass(key: str, default: float) -> float:
            nl = COMMON_NEUTRAL_LOSSES.get(key)
            return nl.mass if nl else default

        ester_refs = [
            ("acetate", "乙酰基(−AcOH)", _lib_mass("CH3COOH", 60.0211)),
            ("benzoate", "苯甲酰基(−PhCOOH)", _lib_mass("C6H5COOH", 122.0368)),
            ("methoxy", "甲氧基(−CH3OH)", _lib_mass("CH3OH", 32.0262)),
            ("hydroxyl", "羟基(−H2O)", _lib_mass("H2O", 18.0106)),
        ]
        tol = self.nl_tol

        # ── 从前体起链, 贪心匹配连续酯基丢失 ──
        cascade: List[Dict[str, Any]] = []
        used: set = set()
        current = float(precursor_mz)
        total_loss_mass = 0.0
        frags_sorted = sorted(set(fragment_mzs), reverse=True)

        for _ in range(6):  # 乌头碱骨架酯基/取代基有限, 6 步余量足够
            best = None  # (frag_mz, loss_type, loss_cn, mass_diff, err)
            for frag in frags_sorted:
                if frag in used or frag >= current - 1.0:
                    continue
                diff = current - frag
                if diff <= 0:
                    continue
                for loss_type, loss_cn, loss_mass in ester_refs:
                    err = abs(diff - loss_mass)
                    if err <= tol and (best is None or err < best[4]):
                        best = (frag, loss_type, loss_cn, diff, err)
            if best is None:
                break
            frag, loss_type, loss_cn, diff, _ = best
            cascade.append({
                "from_mz": round(current, 4),
                "to_mz": round(frag, 4),
                "mass_difference": round(diff, 4),
                "loss_type": loss_type,
                "loss_name_cn": loss_cn,
            })
            total_loss_mass += diff
            used.add(frag)
            current = frag

        n_losses = len(cascade)
        if n_losses == 0:
            return None

        # ── 酯基丢失组合 → 推断乌头碱型亚类 (毒性分级) ──
        loss_types = {c["loss_type"] for c in cascade}
        if "acetate" in loss_types and "benzoate" in loss_types:
            # 观察到 C-8 乙酰基 + C-14 苯甲酰基两者丢失 → 双酯型母体
            aconitine_subtype = "双酯型(DDA, 毒性最高)"
        elif "benzoate" in loss_types:
            aconitine_subtype = "单酯型(MDA, 毒性降低)"
        elif "acetate" in loss_types:
            aconitine_subtype = "脂型/单酯型(LDA/MDA)"
        else:
            aconitine_subtype = "含酯基/甲氧基生物碱"

        composition: Dict[str, int] = {}
        for step in cascade:
            composition[step["loss_type"]] = composition.get(step["loss_type"], 0) + 1

        return {
            "cascade": cascade,
            "n_losses": n_losses,
            "total_loss_mass": round(total_loss_mass, 4),
            "aglycone_mz": round(current, 4),  # 末碎片 = 乌头原碱/去酯骨架离子
            "loss_composition": composition,
            "aconitine_subtype": aconitine_subtype,
            "is_alkaloid_evidence": True,
        }

    # ── 2. 诊断离子检测 ──────────────────────────────────────────

    def _analyze_diagnostic_ions(
        self, fragment_mzs: List[float],
        classification: str = "",
        polarity: str = "",
        candidate_smiles: str = "",
        raw_classification: str = "",
        query_adduct: str = "",
        query_collision_energy: Any = None,
    ) -> Dict[str, Any]:
        """检测谱图中是否含有已知的诊断离子

        分类过滤: 当 classification 提供时,仅匹配与该分类化学兼容的诊断离子,
        避免黄酮类诊断离子出现在酚酸类化合物的分析中。
        """
        matched = []
        unmatched = fragment_mzs.copy()

        for mz in fragment_mzs:
            ions = match_diagnostic_ion(
                mz, self.di_tol,
                classification=classification, strict_fg=True,
                polarity=polarity,
                candidate_smiles=candidate_smiles,
                raw_classification=raw_classification,
                query_adduct=query_adduct,
                query_collision_energy=query_collision_energy,
            )
            # 一个观测峰最多贡献一次诊断证据，避免同质量规则重复放大得分。
            if ions:
                intensity_rank = {"high": 2, "medium": 1, "low": 0}
                ions = sorted(
                    ions,
                    key=lambda ion: (
                        abs(mz - ion.mz),
                        -intensity_rank.get(ion.intensity_characteristic, 0),
                    ),
                )[:1]
            for ion in ions:
                matched.append({
                    "observed_mz": mz,
                    "expected_mz": ion.mz,
                    "delta": round(abs(mz - ion.mz), 5),
                    "ion_name": ion.name,
                    "ion_formula": ion.formula,
                    "functional_group": ion.functional_group,
                    "structural_origin": ion.structural_origin,
                    "intensity_characteristic": ion.intensity_characteristic,
                    "rule_id": ion.rule_id,
                    "polarity": ion.polarity,
                    "evidence_level": ion.evidence_level,
                    "references": list(ion.references),
                })
                if mz in unmatched:
                    unmatched.remove(mz)

        # 同一谱图可能含重复峰（或同一峰经多个上游步骤重复写入）。诊断证据的
        # 计数单位必须是“独立规则”，而不是峰记录数；否则同一核碱基/头基离子
        # 会重复抬高 expert score。每条规则只保留质量误差最小的观测。
        best_by_rule: Dict[Any, Dict[str, Any]] = {}
        for hit in matched:
            evidence_key = hit.get("rule_id") or (
                round(float(hit["expected_mz"]), 4),
                hit.get("ion_formula"),
                hit.get("functional_group"),
            )
            previous = best_by_rule.get(evidence_key)
            if previous is None or hit["delta"] < previous["delta"]:
                best_by_rule[evidence_key] = hit
        matched = list(best_by_rule.values())

        # 按官能团分类汇总
        fg_hits = defaultdict(int)
        for m in matched:
            fg_hits[m["functional_group"]] += 1

        return {
            "matched_ions": matched,
            "matched_count": len(matched),
            "unmatched_peaks": len(unmatched),
            "functional_group_hits": dict(fg_hits),
            "detection_rate": len(matched) / len(fragment_mzs) if fragment_mzs else 0.0,
        }

    # ── 3. 裂解规则覆盖度 ────────────────────────────────────────

    def _analyze_rule_coverage(
        self,
        nl_analysis: Dict[str, Any],
        classification: str,
        polarity: str = "",
        candidate_smiles: str = "",
        query_adduct: str = "",
        query_collision_energy: Any = None,
    ) -> Dict[str, Any]:
        """评估观测到的中性丢失覆盖了哪些已知裂解规则(含顺序丢失)"""
        classification_lower = classification.lower() if classification else ""

        # 映射分类到官能团
        classification_to_fg = {
            "amine_oxide": "amine_oxide",
            "amine oxide": "amine_oxide",
            "aminoxide": "amine_oxide",
            "flavonoid": "flavonoid",
            "alkaloid": "alkaloid",
            "terpenoid": "terpenoid",
            "terpene": "terpenoid",
            "diterpenoid": "terpenoid",
            "sesquiterpenoid": "terpenoid",
            "triterpenoid": "terpenoid",
            "phenolic_acid": "phenolic_acid",
            "phenolic": "phenolic_acid",
            "glycoside": "glycoside",
            "sugar": "sugar",
            "monosaccharide": "sugar",
            "carbohydrate": "sugar",
            "fatty_acid": "fatty_acid",
            "lipid": "lipid",
            "coumarin": "coumarin",
            "quinone": "quinone",
            "steroid": "steroid",
            "amino_acid": "amino_acid",
            "amino": "amino_acid",
            "benzenoid": "benzenoid",
            "sulfonate": "sulfonate",
            "sulfonic": "sulfonate",
            "lignan": "lignan",
            "saponin": "saponin",
            "nucleoside": "nucleoside",
            "nucleotide": "nucleoside",
        }

        fg = None
        for key, val in classification_to_fg.items():
            if key in classification_lower:
                fg = val
                break

        if not fg:
            return {
                "classification": classification,
                "matched_functional_group": None,
                "expected_rules": [],
                "covered_rules": [],
                "coverage_rate": 0.0,
                "note": "unsupported_classification",
            }

        # 从NL分析中提取用于规则匹配(包含阶段1和阶段2)
        nl_for_coverage = []
        neutral_losses = nl_analysis.get("neutral_losses", [])
        sequential_losses = nl_analysis.get("sequential_losses", [])

        # 阶段1: precursor→fragment
        for item in neutral_losses:
            if isinstance(item, dict):
                if "neutral_loss" in item:
                    nl_for_coverage.append((item["mass_difference"],
                        get_neutral_loss_by_mass(item["mass_difference"], self.nl_tol,
                                                 classification=classification,
                                                 strict_fg=True,
                                                 polarity=polarity,
                                                 candidate_smiles=candidate_smiles,
                                                 query_adduct=query_adduct,
                                                 query_collision_energy=query_collision_energy)[0]
                        if get_neutral_loss_by_mass(item["mass_difference"], self.nl_tol,
                                                    classification=classification,
                                                    strict_fg=True,
                                                    polarity=polarity,
                                                    candidate_smiles=candidate_smiles,
                                                    query_adduct=query_adduct,
                                                    query_collision_energy=query_collision_energy)
                        else None))
                elif "mass_difference" in item:
                    diff = item["mass_difference"]
                    matches = get_neutral_loss_by_mass(diff, self.nl_tol,
                                                       classification=classification,
                                                       strict_fg=True,
                                                       polarity=polarity,
                                                       candidate_smiles=candidate_smiles,
                                                       query_adduct=query_adduct,
                                                       query_collision_energy=query_collision_energy)
                    nl_for_coverage.append((diff, matches[0] if matches else None))

        # 阶段2: fragment→fragment 顺序丢失
        for item in sequential_losses:
            if isinstance(item, dict):
                diff = item.get("mass_difference", 0)
                if diff > 0:
                    matches = get_neutral_loss_by_mass(
                        diff, self.nl_tol,
                        classification=classification,
                        strict_fg=True,
                        polarity=polarity,
                        candidate_smiles=candidate_smiles,
                        query_adduct=query_adduct,
                        query_collision_energy=query_collision_energy,
                    )
                    nl_for_coverage.append((diff, matches[0] if matches else None))

        coverage = compute_rule_coverage(
            nl_for_coverage, fg, polarity=polarity,
            mass_tolerance=self.nl_tol,
            candidate_smiles=candidate_smiles,
            query_adduct=query_adduct,
            query_collision_energy=query_collision_energy,
        )
        sig = get_fg_signature(fg)

        # 获取预期规则详情
        expected_rule_ids = coverage.get("covered_rules", []) + coverage.get("uncovered_rules", [])
        expected_rules = []
        for rid in expected_rule_ids:
            rule = FRAGMENTATION_RULES.get(rid)
            if rule:
                expected_rules.append({
                    "rule_id": rid,
                    "name": rule.name,
                    "mechanism_type": rule.mechanism_type,
                    "confidence_score": rule.confidence_score,
                })

        covered_details = []
        for rid in coverage.get("covered_rules", []):
            rule = FRAGMENTATION_RULES.get(rid)
            if rule:
                covered_details.append({
                    "rule_id": rid,
                    "name": rule.name,
                    "mechanism_type": rule.mechanism_type,
                    "confidence_score": rule.confidence_score,
                })

        return {
            "classification": classification,
            "matched_functional_group": fg,
            "expected_rules": expected_rules,
            "covered_rules": covered_details,
            "coverage_rate": coverage.get("coverage_rate", 0.0),
            "uncovered_rules": coverage.get("uncovered_rules", []),
            "matched_neutral_losses": coverage.get("matched_neutral_losses", []),
            "excluded_rules": coverage.get("excluded_rules", {}),
            "fragmentation_pathway": sig["fragmentation_pathway"] if sig else "",
        }

    # ── 4. 官能团特征验证 ────────────────────────────────────────

    def _analyze_functional_group_signature(
        self,
        neutral_losses: List[Dict[str, Any]],
        diagnostic_ions: List[Dict[str, Any]],
        classification: str,
        polarity: str = "",
    ) -> Dict[str, Any]:
        """
        验证谱图是否显示目标化合物类型的特征裂解模式

        不同化合物类型有独特的"指纹"：黄酮有糖丢失+rDA，生物碱有α-裂解等。
        """
        sig = FUNCTIONAL_GROUP_SIGNATURES.get(
            self._classify_to_fg(classification)
        )

        if not sig:
            return {"matched": False, "signature_found": None,
                    "characteristic_loss_detection": 0.0,
                    "diagnostic_ion_detection": 0.0,
                    "overall_signature_score": 0.0}

        # 检测特征中性丢失
        expected_losses = sig["characteristic_losses"]
        detected_losses = []
        for item in neutral_losses:
            if "neutral_loss" in item or "mass_difference" in item:
                mass = item.get("mass_difference", item.get("mass", 0))
                for expected in expected_losses:
                    if abs(mass - expected) <= self.nl_tol:
                        detected_losses.append(expected)
                        break

        loss_detection = len(set(detected_losses)) / len(expected_losses) if expected_losses else 0.0

        # 检测特征诊断离子
        wanted_polarity = normalize_polarity(polarity)
        expected_ions = []
        for expected in sig["diagnostic_ions"]:
            catalog_entries = [
                ion for ion in DIAGNOSTIC_IONS.values()
                if abs(ion.mz - expected) <= self.di_tol
            ]
            if not catalog_entries:
                continue  # 无公式/极性元数据的旧数值不进入生产分母
            if wanted_polarity and not any(
                normalize_polarity(ion.polarity) in (None, wanted_polarity)
                for ion in catalog_entries
            ):
                continue
            expected_ions.append(expected)
        detected_ions = []
        for ion in diagnostic_ions:
            for expected in expected_ions:
                if abs(ion["observed_mz"] - expected) <= self.di_tol:
                    detected_ions.append(expected)
                    break

        ion_detection = len(set(detected_ions)) / len(expected_ions) if expected_ions else 0.0

        # 综合特征得分
        overall = (loss_detection * 0.5 + ion_detection * 0.5)

        return {
            "matched": overall > 0.3,
            "signature_found": sig,
            "expected_losses": expected_losses,
            "detected_losses": list(set(detected_losses)),
            "characteristic_loss_detection": loss_detection,
            "expected_ions": expected_ions,
            "detected_ions": list(set(detected_ions)),
            "diagnostic_ion_detection": ion_detection,
            "overall_signature_score": overall,
            "fragmentation_pathway": sig["fragmentation_pathway"],
        }

    # ── 5. 结构-碎片一致性 ──────────────────────────────────────

    def _analyze_structure_fragment_consistency(
        self,
        smiles: str,
        fragment_mzs: List[float],
        precursor_mz: float,
        classification: str = "",
        polarity: str = "",
        query_adduct: str = "",
        query_charge: int = 1,
    ) -> Dict[str, Any]:
        """
        评估候选化合物结构能否合理解释观测到的碎片

        使用 RDKit 计算可能的碎片质量，与观测比较。
        """
        if not smiles or not fragment_mzs:
            return {"score": 0.0, "plausible_fragments": [],
                    "consistency_note": "no_structure_or_no_fragments"}

        try:
            from rdkit import Chem
        except ImportError:
            return {"score": 0.0, "plausible_fragments": [],
                    "consistency_note": "rdkit_not_available"}

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"score": 0.0, "plausible_fragments": [],
                    "consistency_note": "invalid_smiles"}

        # 使用键断裂方法生成可能的碎片
        best_by_observed: Dict[float, Dict[str, Any]] = {}
        effective_polarity = normalize_polarity(polarity)
        adduct_text = str(query_adduct or "")
        if effective_polarity == "positive":
            adduct_shifts = [(1.007276, "[M+H]+"), (0.0, "[M]+")]
            if "Na" in adduct_text:
                adduct_shifts.insert(1, (22.989218, "[M+Na]+"))
        elif effective_polarity == "negative":
            adduct_shifts = [(-1.007276, "[M-H]-"), (0.0, "[M]-")]
        else:
            adduct_shifts = [
                (1.007276, "[M+H]+"),
                (-1.007276, "[M-H]-"),
                (22.989218, "[M+Na]+"),
                (0.0, "[M]"),
            ]
        fragment_tolerance = min(max(self.di_tol, 0.02), 0.05)
        try:
            frags = self._generate_possible_fragments(mol)
            for frag_smi, frag_mass in frags:
                for obs_mz in fragment_mzs:
                    for adduct_mass, adduct_label in adduct_shifts:
                        expected = frag_mass + adduct_mass
                        delta = abs(obs_mz - expected)
                        if delta <= fragment_tolerance:
                            candidate = {
                                "observed_mz": round(obs_mz, 4),
                                "fragment_smiles": frag_smi,
                                "fragment_mass": round(frag_mass, 4),
                                "expected_mz": round(expected, 4),
                                "delta": round(delta, 5),
                                "adduct": adduct_label,
                                "assignment_type": "rdkit_single_bond_fragment",
                            }
                            previous = best_by_observed.get(obs_mz)
                            if previous is None or candidate["delta"] < previous["delta"]:
                                best_by_observed[obs_mz] = candidate
        except Exception as e:
            logger.debug(f"Fragment generation failed: {e}")

        # 补充：基于中性丢失预测碎片(不依赖RDKit键断裂)
        # 如果前体质量已知，观测碎片 = 前体质量 - 中性丢失
        try:
            from .fragmentation_rules import get_neutral_loss_by_mass
            if abs(int(query_charge or 1)) == 1:
                for obs_mz in fragment_mzs:
                    if obs_mz in best_by_observed:
                        continue
                    # 对同为单电荷的前体/产物，m/z 差即中性丢失质量；不应先假设[M+H]+。
                    diff = precursor_mz - obs_mz
                    if diff > 5:
                        nl_matches = get_neutral_loss_by_mass(
                            diff, self.nl_tol,
                            classification=classification, strict_fg=True,
                        )
                        if not nl_matches:
                            continue
                        if effective_polarity == "negative":
                            nl_adduct = "[M-H]- (NL-based)"
                        elif effective_polarity == "positive":
                            nl_adduct = "[M+H]+ (NL-based)"
                        else:
                            nl_adduct = "single-charge (NL-based)"
                        best_by_observed[obs_mz] = {
                            "observed_mz": round(obs_mz, 4),
                            "fragment_smiles": "",
                            "fragment_mass": round(obs_mz, 4),
                            "expected_mz": round(obs_mz, 4),
                            "delta": 0.0,
                            "adduct": nl_adduct,
                            "neutral_loss": nl_matches[0].name,
                            "assignment_type": "neutral_loss_mass_difference",
                        }
        except Exception:
            pass

        plausible_fragments = [
            best_by_observed[mz] for mz in sorted(best_by_observed)
        ]

        if fragment_mzs:
            consistency_score = len(plausible_fragments) / len(fragment_mzs)
            consistency_score = min(1.0, consistency_score * 2)  # Amplify, ceiling at 1
        else:
            consistency_score = 0.0

        return {
            "score": consistency_score,
            "plausible_fragments": plausible_fragments,
            "total_fragments_observed": len(fragment_mzs),
            "fragments_explained": len(plausible_fragments),
            "consistency_note": f"{len(plausible_fragments)}/{len(fragment_mzs)} "
                               f"fragments structurally plausible",
        }

    def _generate_possible_fragments(
        self, mol
    ) -> List[Tuple[str, float]]:
        """
        通过断开单键生成可能的碎片

        这是一个计算昂贵的操作，只断开最不稳定的键。
        """
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        fragments = []
        bond_labels = estimate_bond_lability(Chem.MolToSmiles(mol))

        # 只考虑最不稳定的前5个键
        labile_bonds = bond_labels.get("labile_bonds", [])[:5]

        for bond_info in labile_bonds:
            bond_idx = bond_info["bond_idx"]
            bond = mol.GetBondWithIdx(bond_idx)
            if bond is None:
                continue

            # 创建可编辑分子并断键
            try:
                rw_mol = Chem.RWMol(mol)
                rw_mol.RemoveBond(
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                )
                # 获取碎片
                frag_mols = Chem.GetMolFrags(rw_mol.GetMol(), asMols=True,
                                              sanitizeFrags=True)
                for frag in frag_mols:
                    try:
                        frag_smi = Chem.MolToSmiles(frag)
                        frag_mass = Descriptors.ExactMolWt(frag)
                        fragments.append((frag_smi, frag_mass))
                    except Exception:
                        pass
            except Exception:
                continue

        return fragments

    # ── 6. 键能合理性评估 ────────────────────────────────────────

    def _analyze_bond_energy_plausibility(
        self,
        smiles: str,
        neutral_losses: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """评估推测的断键位点在能量上是否合理

        改进：结合观测到的中性丢失数量调节得分。
        - 如果结构有易断裂键 AND 观测到多种中性丢失 → 高分（合理裂解）
        - 如果结构有易断裂键但观测到的丢失很少 → 中等分（可能匹配不佳）
        - 如果结构稳定（无易断裂键）→ 低分
        """
        if not smiles:
            return {"score": 0.0, "analysis": [], "note": "no_structure"}

        lability = estimate_bond_lability(smiles)
        labile_bonds = lability.get("labile_bonds", [])

        if not labile_bonds:
            return {"score": 0.0, "analysis": [], "note": "no_labile_bonds_found"}

        # 统计已解释的中性丢失数量
        explained_nl = [nl for nl in neutral_losses
                        if isinstance(nl, dict) and nl.get("neutral_loss")]
        n_explained = len(explained_nl)
        if n_explained == 0:
            return {
                "score": 0.0,
                "analysis": labile_bonds[:5],
                "note": "labile_bonds_are_prior_only_without_observed_neutral_loss",
                "evidence_level": "heuristic_prior_not_scored",
            }

        # 检查是否有高倾向性键存在
        high_lability = [b for b in labile_bonds if b.get("lability_score", 0) >= 0.7]
        medium_lability = [b for b in labile_bonds if 0.5 <= b.get("lability_score", 0) < 0.7]

        # 结合键易断裂程度和实际观测到的中性丢失数量
        if high_lability:
            score = min(0.6, 0.2 + n_explained * 0.08)
        elif medium_lability:
            score = min(0.4, 0.1 + n_explained * 0.06)
        else:
            score = min(0.2, n_explained * 0.04)

        return {
            "score": round(score, 4),
            "high_lability_bonds": len(high_lability),
            "medium_lability_bonds": len(medium_lability),
            "total_labile_bonds": len(labile_bonds),
            "explained_neutral_losses": n_explained,
            "most_labile_region": lability.get("most_labile_region", "unknown"),
            "top_labile_bonds": labile_bonds[:5],
            "fragmentation_sites_estimated": lability.get("estimated_fragmentation_sites", 0),
            "evidence_level": "heuristic_supporting_only",
            "note": "High-lability bonds exist → fragmentation chemically plausible"
                    if high_lability else "No highly labile bonds → structure may be stable",
        }

    # ── 7. v3.0: 碎片特征库验证 ────────────────────────────────────

    def _analyze_feature_library(
        self,
        fragment_mzs: List[float],
        classification: str,
        polarity: str = "",
        candidate_smiles: str = "",
        raw_classification: str = "",
        query_adduct: str = "",
        query_collision_energy: Any = None,
    ) -> Dict[str, Any]:
        if not FEATURE_LIB_AVAILABLE or not fragment_mzs:
            return {"score": 0.0, "matched_features": [],
                    "total_expected": 0, "note": "feature_library_unavailable"}

        normalized_class = normalize_classification(classification)
        if normalized_class == "Unknown":
            return {
                "score": 0.0,
                "matched_features": [],
                "total_expected": 0,
                "note": "classification_unconfirmed",
            }

        # Class feature profiles are class-specific evidence.  Fragment-only
        # classification is retained as a descriptive hypothesis upstream, but
        # cannot select its own profile and then score the same fragments.
        profile = get_profile(normalized_class.lower())

        if not profile:
            return {"score": 0.0, "matched_features": [],
                    "total_expected": 0, "note": f"no_profile_for_{normalized_class}"}

        expected_fragments = get_fragments_by_class(
            profile.class_name, polarity=polarity, validated_only=True,
            candidate_smiles=candidate_smiles,
            raw_classification=raw_classification,
            enforce_scope=True,
            query_adduct=query_adduct,
            query_collision_energy=query_collision_energy,
        )
        matched = []
        for frag in expected_fragments:
            for obs_mz in fragment_mzs:
                if abs(obs_mz - frag.mz) <= self.di_tol:
                    matched.append({
                        "observed_mz": round(obs_mz, 4),
                        "expected_mz": frag.mz,
                        "fragment_name": frag.name,
                        "fragment_name_en": frag.name_en,
                        "intensity_rank": frag.intensity_rank,
                        "diagnostic_value": frag.diagnostic_value,
                        "structural_origin": frag.structural_origin,
                        "polarity": frag.polarity,
                        "formula_quality": "validated",
                    })
                    break

        weighted_score = 0.0
        for m in matched:
            if m["intensity_rank"] == "base":
                weighted_score += 0.4
            elif m["intensity_rank"] == "major":
                weighted_score += 0.3
            elif m["intensity_rank"] == "medium":
                weighted_score += 0.15
            else:
                weighted_score += 0.05

        total_expected = len(expected_fragments)
        raw_rate = len(matched) / max(total_expected, 1)
        score = min(1.0, raw_rate * 0.5 + weighted_score * 0.5 / max(total_expected * 0.2, 1))

        return {
            "score": round(score, 4),
            "matched_features": matched,
            "matched_count": len(matched),
            "total_expected": total_expected,
            "match_rate": raw_rate,
            "weighted_score": round(weighted_score, 4),
            "class_profile": profile.class_name_cn,
        }

    # ── 8. v3.0: 异构体区分分析 ────────────────────────────────────

    def _analyze_isomer_differentiation(
        self,
        fragment_mzs: List[float],
        classification: str,
        polarity: str = "",
    ) -> Dict[str, Any]:
        if not FEATURE_LIB_AVAILABLE:
            return {"score": 0.0, "detected_pairs": [],
                    "isomer_evidence": [], "note": "feature_library_unavailable"}

        ion_pairs = get_diagnostic_ion_pairs(
            classification.lower(), polarity=polarity
        )
        if not ion_pairs:
            return {"score": 0.0, "detected_pairs": [],
                    "isomer_evidence": [],
                    "note": "no_isomer_rules_for_class"}

        detected_pairs = []
        for pair in ion_pairs:
            mz1, mz2 = pair
            found1 = any(abs(mz - mz1) <= self.di_tol for mz in fragment_mzs)
            found2 = any(abs(mz - mz2) <= self.di_tol for mz in fragment_mzs)
            if found1 and found2:
                detected_pairs.append({
                    "ion_1_mz": mz1, "ion_2_mz": mz2,
                    "both_detected": True,
                })
            elif found1 or found2:
                detected_pairs.append({
                    "ion_1_mz": mz1, "ion_2_mz": mz2,
                    "both_detected": False, "partial_detection": True,
                })

        isomer_rules = get_isomer_rules(classification.lower())
        isomer_evidence = list(isomer_rules[:3]) if isomer_rules else []

        if detected_pairs:
            complete_pairs = [p for p in detected_pairs if p.get("both_detected")]
            partial_pairs = [p for p in detected_pairs if p.get("partial_detection")]
            score = min(1.0, len(complete_pairs) * 0.4 + len(partial_pairs) * 0.15)
        else:
            score = 0.0

        return {
            "score": round(score, 4),
            "detected_pairs": detected_pairs,
            "complete_pairs": len([p for p in detected_pairs if p.get("both_detected")]),
            "partial_pairs": len([p for p in detected_pairs if p.get("partial_detection")]),
            "isomer_evidence": isomer_evidence,
            "differentiation_power": "strong" if score >= 0.6 else (
                "moderate" if score >= 0.3 else "weak"),
        }

    # ── 综合评分计算 v3.0 ────────────────────────────────────────

    def _compute_expert_scores(
        self,
        nl_analysis: Dict[str, Any],
        di_analysis: Dict[str, Any],
        rc_analysis: Dict[str, Any],
        fg_analysis: Dict[str, Any],
        sf_analysis: Dict[str, Any],
        be_analysis: Dict[str, Any],
        ff_analysis: Optional[Dict[str, Any]] = None,
        id_analysis: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, float]:
        """v3.1 增强权重 — 聚合8个评估维度，动态权重重分配"""

        # 1. 中性丢失得分
        nl_score = nl_analysis.get("explanation_rate", 0.0) * 0.7
        nl_score += nl_analysis.get("priority_weighted_score", 0.0) * 0.3

        # 2. 诊断离子得分
        di_score = di_analysis.get("detection_rate", 0.0)
        has_high_specificity = any(
            ion.get("intensity_characteristic") == "high"
            for ion in di_analysis.get("matched_ions", [])
        )
        di_score = min(1.0, di_score + (0.2 if has_high_specificity else 0.0))

        # 3-6. 传统维度
        rc_score = rc_analysis.get("coverage_rate", 0.0)
        fg_score = fg_analysis.get("overall_signature_score", 0.0)
        sf_score = sf_analysis.get("score", 0.0)
        be_score = be_analysis.get("score", 0.0)

        # 7-8. v3.0 新增维度
        ff_score = ff_analysis.get("score", 0.0) if ff_analysis else 0.0
        id_score = id_analysis.get("score", 0.0) if id_analysis else 0.0

        # v3.1 基础权重组合
        weights = {
            "neutral_loss": 0.22,
            "rule_coverage": 0.22,
            "fg_signature": 0.18,
            "feature_library": 0.15,
            "structure_fragment": 0.10,
            "diagnostic_ion": 0.08,
            "bond_energy": 0.03,
            "isomer_differentiation": 0.02,
        }

        # ── 动态权重重分配 ──
        # 当分类未知时，rc/fg/ff 维度为0或不可靠，
        # 将其权重按比例分配到可靠的维度(nl, di, be, sf)
        classification = rc_analysis.get("classification", "Unknown")
        is_unclassified = classification and classification.lower() in (
            "n/a", "unknown", "", "other", "none", "unsupported_classification"
        )

        if is_unclassified:
            # 不可靠维度仅保留合计 0.04 的审计权重，其余再分配。旧实现先把
            # 全部 0.57 重分配、随后又保留 0.04，导致总权重变成 1.04。
            unreliable_weight = (
                weights["rule_coverage"]
                + weights["fg_signature"]
                + weights["feature_library"]
                + weights["isomer_differentiation"]
            ) - 0.04
            # 将不可靠权重按比例分配到可靠维度
            reliable_keys = ["neutral_loss", "diagnostic_ion",
                             "structure_fragment", "bond_energy"]
            reliable_base = sum(weights[k] for k in reliable_keys)
            for k in reliable_keys:
                weights[k] += unreliable_weight * (weights[k] / reliable_base)
            # 不可靠维度权重设为很小的值（不完全归零，保留微弱信号）
            weights["rule_coverage"] = 0.01
            weights["fg_signature"] = 0.01
            weights["feature_library"] = 0.01
            weights["isomer_differentiation"] = 0.01

        total = (
            nl_score * weights["neutral_loss"]
            + rc_score * weights["rule_coverage"]
            + fg_score * weights["fg_signature"]
            + ff_score * weights["feature_library"]
            + sf_score * weights["structure_fragment"]
            + di_score * weights["diagnostic_ion"]
            + be_score * weights["bond_energy"]
            + id_score * weights["isomer_differentiation"]
        )

        return {
            "total": round(total, 4),
            "neutral_loss": round(nl_score, 4),
            "diagnostic_ion": round(di_score, 4),
            "rule_coverage": round(rc_score, 4),
            "fg_signature": round(fg_score, 4),
            "structure_fragment": round(sf_score, 4),
            "bond_energy": round(be_score, 4),
            "feature_library": round(ff_score, 4),
            "isomer_differentiation": round(id_score, 4),
        }

    # ── 辅助函数 ─────────────────────────────────────────────────

    @staticmethod
    def _expert_confidence(score: float) -> str:
        """专家评分 → 置信度等级"""
        if score >= 0.70:
            return "High"
        elif score >= 0.45:
            return "Medium"
        elif score >= 0.25:
            return "Low"
        return "Insufficient"

    @staticmethod
    def _classify_to_fg(classification: str) -> str:
        """分类映射到官能团 (v3.3.0 expanded)"""
        c = classification.lower() if classification else ""
        mapping = {
            "amine_oxide": "amine_oxide", "amine oxide": "amine_oxide",
            "aminoxide": "amine_oxide",
            "flavonoid": "flavonoid", "alkaloid": "alkaloid",
            "terpenoid": "terpenoid", "terpene": "terpenoid",
            "phenolic_acid": "phenolic_acid", "phenolic": "phenolic_acid",
            "glycoside": "glycoside", "sugar": "sugar",
            "monosaccharide": "sugar", "carbohydrate": "sugar",
            "fatty_acid": "fatty_acid",
            "lipid": "lipid",
            "coumarin": "coumarin", "quinone": "quinone",
            "steroid": "steroid", "amino_acid": "amino_acid",
            "amino": "amino_acid", "benzenoid": "benzenoid",
            "sulfonate": "sulfonate", "sulfonic": "sulfonate",
            "lignan": "lignan", "saponin": "saponin",
        }
        for key, val in mapping.items():
            if key in c:
                return val
        return "unknown"

    @staticmethod
    def _empty_result(reason: str) -> Dict[str, Any]:
        result = {
            "expert_total_score": 0.0,
            "neutral_loss_score": 0.0,
            "diagnostic_ion_score": 0.0,
            "rule_coverage_score": 0.0,
            "fg_signature_score": 0.0,
            "structure_fragment_score": 0.0,
            "bond_energy_score": 0.0,
            "feature_library_score": 0.0,
            "isomer_differentiation_score": 0.0,
            "confidence_level": "Insufficient",
            "neutral_loss_analysis": {},
            "diagnostic_ion_analysis": {},
            "rule_coverage_analysis": {},
            "fg_signature_analysis": {},
            "structure_fragment_analysis": {},
            "bond_energy_analysis": {},
            "feature_library_analysis": {},
            "isomer_differentiation_analysis": {},
            "evaluation_status": "unavailable",
            "evaluation_reason": reason,
            "chemical_reasoning": f"Insufficient data: {reason}",
        }
        return finalize_expert_evidence(result, fused_score=0.0)

    def _generate_reasoning(
        self,
        nl_analysis: Dict[str, Any],
        di_analysis: Dict[str, Any],
        rc_analysis: Dict[str, Any],
        fg_analysis: Dict[str, Any],
        sf_analysis: Dict[str, Any],
        ff_analysis: Optional[Dict[str, Any]] = None,
        id_analysis: Optional[Dict[str, Any]] = None,
        effective_classification: str = "Unknown",
        cls_reliable: bool = True,
    ) -> str:
        """v3.0: 生成人类可读的化学推理摘要"""
        parts = []

        # 分类不可靠(仅碎片库推断)时, 类特异 diagnostic ion / rule coverage 已在上游
        # 用 Unknown 降级跳过 —— 明确告知用户分类待确认, 避免把生物碱误报为黄酮等
        if not cls_reliable:
            parts.append(
                "🏷️ Classification unconfirmed (inferred from fragments only; "
                "class-specific rules skipped to avoid cross-class misreporting)"
            )

        # 中性丢失
        nl_rate = nl_analysis.get("explanation_rate", 0)
        seq_explained = nl_analysis.get("sequential_explained", 0)
        nl_verb = "✅" if nl_rate > 0.5 else ("⚠️" if nl_rate > 0.2 else "❌")
        parts.append(
            f"{nl_verb} Neutral losses: {nl_rate:.0%} explained "
            f"(precursor→frag: {nl_analysis.get('explained_count', 0)}, "
            f"frag→frag sequential: {seq_explained})"
        )

        # 连通碎裂路径 (使用真实的连通链替代旧的扁平拼接)
        # v3.1 防御性验证: 即使上游过滤正确，也做最终一致性检查
        best_pathway = nl_analysis.get("best_pathway", [])
        if best_pathway:
            chain_parts = []
            incompatible_nls = []
            for step in best_pathway[:4]:  # 最多显示 4 步
                nl_name = step.get("neutral_loss", "")
                # 最终防御: 检查 NL 是否与有效分类兼容
                if effective_classification and effective_classification != "Unknown":
                    nl_obj = COMMON_NEUTRAL_LOSSES.get(
                        nl_name,
                        COMMON_NEUTRAL_LOSSES.get(
                            step.get("formula", ""), None
                        )
                    )
                    # 通过名称反向查找 NL 对象
                    if nl_obj is None:
                        for _key, _nl in COMMON_NEUTRAL_LOSSES.items():
                            if _nl.name == nl_name:
                                nl_obj = _nl
                                break
                    if nl_obj and not _is_nl_compatible(
                        nl_obj, effective_classification
                    ):
                        incompatible_nls.append(nl_name)
                        continue  # 跳过不兼容的 NL

                chain_parts.append(
                    f"{nl_name} "
                    f"({step.get('from_mz', 0):.1f}→{step.get('to_mz', 0):.1f})"
                )
            if chain_parts:
                parts.append(f"   Sequential pathway: {' → '.join(chain_parts)}")
            if incompatible_nls:
                parts.append(
                    f"   ⚠️ Filtered out incompatible NLs: {', '.join(incompatible_nls)} "
                    f"(not consistent with class '{effective_classification}')"
                )
        elif seq_explained == 0:
            # 无连通链时不显示虚假路径
            pass

        # 类别专属级联 (阶段3): 糖苷糖基 / 乌头碱酯基 — 从前体起链的完整级联
        sugar_cascade = nl_analysis.get("sugar_loss_cascade")
        if sugar_cascade and sugar_cascade.get("n_sugars", 0) >= 1:
            sc_steps = sugar_cascade["cascade"]
            start_mz = sc_steps[0]["from_mz"]
            chain_str = "→".join("{:.1f}".format(s["to_mz"]) for s in sc_steps)
            comp = sugar_cascade.get("sugar_composition", {})
            comp_str = ", ".join("{}×{}".format(cnt, name) for name, cnt in comp.items())
            parts.append(
                "🍬 Sugar-loss cascade: {:.1f}→{} ({} loss) → aglycone m/z {:.1f}, "
                "{} sugar(s)".format(
                    start_mz, chain_str, comp_str,
                    sugar_cascade["aglycone_mz"], sugar_cascade["n_sugars"],
                )
            )

        ester_cascade = nl_analysis.get("alkaloid_ester_cascade")
        if ester_cascade and ester_cascade.get("n_losses", 0) >= 1:
            ec_steps = ester_cascade["cascade"]
            start_mz = ec_steps[0]["from_mz"]
            chain_str = "→".join("{:.1f}".format(s["to_mz"]) for s in ec_steps)
            comp = ester_cascade.get("loss_composition", {})
            comp_str = ", ".join("{}×{}".format(cnt, name) for name, cnt in comp.items())
            parts.append(
                "⚗️ Alkaloid ester-loss cascade: {:.1f}→{} ({}) → {}".format(
                    start_mz, chain_str, comp_str,
                    ester_cascade["aconitine_subtype"],
                )
            )

        # 诊断离子 (含跨类一致性检查)
        di_count = di_analysis.get("matched_count", 0)
        if di_count > 0:
            ion_names = [ion["ion_name"] for ion in di_analysis.get("matched_ions", [])[:3]]
            parts.append(f"🔬 {di_count} diagnostic ion(s): {', '.join(ion_names)}")

            # ── 跨类矛盾检测 ──
            # 如果诊断离子属于与最终分类不同的化合物类，发出警告
            if effective_classification != "Unknown":
                from .fragmentation_rules import _is_diagnostic_ion_compatible
                for ion in di_analysis.get("matched_ions", []):
                    ion_fg = ion.get("functional_group", "").lower()
                    # 构造临时 DiagnosticIon 用于兼容性检查
                    tmp_ion = DiagnosticIon(
                        mz=0.0, formula="", name="",
                        structural_origin="", functional_group=ion_fg,
                        rule_id="",
                    )
                    if not _is_diagnostic_ion_compatible(tmp_ion, effective_classification):
                        parts.append(
                            f"   ⚠️ Cross-class match: '{ion.get('ion_name', '?')}' "
                            f"(fg={ion_fg}) not consistent with class '{effective_classification}'"
                        )

        # 裂解规则
        rc_rate = rc_analysis.get("coverage_rate", 0)
        covered = rc_analysis.get("covered_rules", [])
        if covered:
            rule_names = [r["name"] for r in covered[:3]]
            parts.append(f"📏 {rc_rate:.0%} rule coverage: {', '.join(rule_names)}")

        # v3.0: 碎片特征库验证 (类标签依赖可靠分类; 不可靠时跳过避免误导)
        if cls_reliable and ff_analysis and ff_analysis.get("matched_count", 0) > 0:
            ff_count = ff_analysis["matched_count"]
            ff_class = ff_analysis.get("class_profile", "")
            parts.append(f"📚 {ff_count} feature library matches ({ff_class})")

        # v3.0: 异构体区分
        if id_analysis and id_analysis.get("complete_pairs", 0) > 0:
            parts.append(
                f"🔍 Isomer differentiation: {id_analysis.get('differentiation_power', 'N/A')} "
                f"({id_analysis['complete_pairs']} complete ion pair(s))"
            )

        # 官能团特征
        fg_score = fg_analysis.get("overall_signature_score", 0)
        if fg_score > 0.5:
            path = fg_analysis.get("fragmentation_pathway", "")
            parts.append(f"🧩 FG signature: {fg_score:.0%}" + (f" (pathway: {path[:80]}...)" if path else ""))

        # 结构一致性
        sf_score = sf_analysis.get("score", 0)
        if sf_score > 0.4:
            parts.append(f"🏗️ Structure-fragment: {sf_score:.0%}")

        return " | ".join(parts) if parts else "Insufficient data for chemical reasoning"


# ─────────────────────────────────────────────────────────────────
# 便捷函数 — 无需实例化即可使用
# ─────────────────────────────────────────────────────────────────

# 全局单例
def enrich_neutral_loss_analysis(
    analysis: Dict[str, Any],
    classification: str = "",
    tolerance: float = 0.05,
) -> Dict[str, Any]:
    """Upgrade legacy neutral-loss details for display without re-searching.

    Old task JSON labels every precursor-to-fragment difference as a neutral
    loss.  This function assigns conservative two-step cumulative losses where
    supported and explicitly marks the remainder as unassigned fragment
    differences.  It mutates and returns ``analysis``.
    """
    if not isinstance(analysis, dict):
        return analysis
    items = analysis.get("neutral_losses") or []
    if not isinstance(items, list):
        return analysis
    for item in items:
        if not isinstance(item, dict) or item.get("neutral_loss"):
            continue
        diff = item.get("mass_difference")
        try:
            diff = float(diff)
        except (TypeError, ValueError):
            item["assignment_status"] = "unassigned_fragment_difference"
            continue
        composites = get_composite_neutral_loss_by_mass(
            diff, tolerance=tolerance,
            classification=classification, strict_fg=True,
        )
        if composites:
            composite = composites[0]
            item.update({
                "neutral_loss": f"累积丢失：{composite['name']}",
                "neutral_loss_en": f"Cumulative loss: {composite['name_en']}",
                "formula": composite["formula"],
                "mechanism": composite["mechanism"],
                "priority": composite["priority"],
                "assignment_level": "tentative_composite",
                "mass_delta": composite["delta"],
            })
            item.pop("assignment_status", None)
        else:
            item["assignment_status"] = "unassigned_fragment_difference"

    explained = sum(bool(item.get("neutral_loss")) for item in items if isinstance(item, dict))
    unexplained = sum(not item.get("neutral_loss") for item in items if isinstance(item, dict))
    sequential = int(analysis.get("sequential_explained") or 0)
    cascade = int(analysis.get("cascade_recovered") or 0)
    total = max(len(items) + sequential + cascade, 1)
    analysis["explained_count"] = explained + cascade
    analysis["unexplained_count"] = unexplained
    analysis["explanation_rate"] = (explained + sequential + cascade) / total
    return analysis


_scorer_instance: Optional[ChemicalExpertScorer] = None


def get_expert_scorer() -> ChemicalExpertScorer:
    """获取全局专家评分器实例"""
    global _scorer_instance
    if _scorer_instance is None:
        _scorer_instance = ChemicalExpertScorer()
    return _scorer_instance


def expert_evaluate(
    query_peaks: List[List[float]],
    query_pmz: float,
    candidate_smiles: str,
    candidate_formula: str = "",
    candidate_classification: str = "Unknown",
    candidate_inchikey: str = "",
    query_polarity: str = "",
    query_adduct: str = "",
    query_charge: Optional[int] = None,
    query_collision_energy: Any = None,
    candidate_adduct: str = "",
    candidate_classification_source: Optional[str] = None,
    candidate_classification_confidence: Optional[float] = None,
) -> Dict[str, Any]:
    """
    便捷函数：执行专家评估

    Usage:
        result = expert_evaluate(peaks, pmz, smiles, formula, classification)
        print(f"Expert Score: {result['expert_total_score']:.3f}")
    """
    scorer = get_expert_scorer()
    return scorer.evaluate(
        query_peaks=query_peaks,
        query_pmz=query_pmz,
        candidate_smiles=candidate_smiles,
        candidate_formula=candidate_formula,
        candidate_classification=candidate_classification,
        candidate_inchikey=candidate_inchikey,
        query_polarity=query_polarity,
        query_adduct=query_adduct,
        query_charge=query_charge,
        query_collision_energy=query_collision_energy,
        candidate_adduct=candidate_adduct,
        candidate_classification_source=candidate_classification_source,
        candidate_classification_confidence=candidate_classification_confidence,
    )
