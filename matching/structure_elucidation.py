"""
结构解析推理引擎 — Structure Elucidation Reasoning Engine
==========================================================

从MS/MS谱图出发，结合质谱裂解原理和碎片特征标准库，逐步推理化合物结构。
实现从"匹配"到"解析"的跨越 —— 不仅找出库中最相似的化合物，
更能解释为什么这个结构是合理的，以及如何从碎片反推结构。

核心推理链:
  1. 前体离子分析 → 确定分子量/元素组成/加合物类型
  2. 中性丢失分析 → 识别官能团和亚结构
  3. 诊断离子检测 → 确定化合物类别
  4. 裂解路径重建 → 重建从母离子到各碎片的完整路径
  5. 结构假设生成 → 基于证据生成候选亚结构
  6. 假设验证与排序 → 使用专家评分选择最佳匹配

知识来源:
  - 《质谱分析技术原理与应用》 台湾质谱学会 (2019)
  - 《中药质谱分析》 刘淑莹等 (2012)
  - Spectraverse (Gupta et al., Anal. Chem. 2026)

Author: Enhanced Expert System Module
Version: 3.0.0
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from .fragment_feature_library import (
    COMPOUND_CLASS_PROFILES,
    ADDUCT_PROFILES,
    STRUCTURE_ELUCIDATION_RULES,
    match_fragment_to_class,
    classify_compound_from_fragments,
    get_profile,
    get_diagnostic_ion_pairs,
    get_fragmentation_pathway,
    get_isomer_rules,
    get_structure_elucidation_rules,
    FragmentFeature,
    CompoundClassProfile,
)

from .fragmentation_rules import (
    FRAGMENTATION_RULES,
    COMMON_NEUTRAL_LOSSES,
    DIAGNOSTIC_IONS,
    get_neutral_loss_by_mass,
    match_diagnostic_ion,
    compute_possible_neutral_losses,
    compute_rule_coverage,
    get_fg_signature,
)

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────
# 数据结构
# ─────────────────────────────────────────────────────────────────


@dataclass
class ElucidationEvidence:
    """单条推理证据"""
    evidence_type: str           # evidence type
    description: str             # 人类可读描述
    supporting_mz: List[float]   # 支持此证据的碎片m/z
    confidence: float            # 置信度 (0-1)
    alternative_explanation: str = ""  # 其他可能的解释
    source: str = ""             # 知识来源


@dataclass
class FragmentationPath:
    """裂解路径节点"""
    mz: float                    # 当前节点 m/z
    parent_mz: Optional[float]   # 父节点 m/z
    neutral_loss: Optional[str]  # 中性丢失类型
    neutral_loss_mass: float = 0.0  # 中性丢失精确质量
    assigned_formula: str = ""   # 指定分子式
    confidence: float = 0.5      # 该路径步的置信度
    children: List["FragmentationPath"] = field(default_factory=list)


@dataclass
class StructureHypothesis:
    """结构假设"""
    smiles: str                  # 假设的SMILES
    compound_class: str          # 化合物类别
    substructure_evidence: List[ElucidationEvidence]
    fragmentation_path: FragmentationPath
    expert_score: float = 0.0
    mass_match_ppm: float = 0.0
    ranking: int = 0


@dataclass
class ElucidationReport:
    """完整结构解析报告"""
    precursor_mz: float
    adduct_type: str
    neutral_mass: float
    inferred_formula: str
    inferred_compound_class: str
    class_confidence: float
    diagnostic_evidence: List[ElucidationEvidence]
    fragmentation_tree: FragmentationPath
    structure_hypotheses: List[StructureHypothesis]
    best_hypothesis: Optional[StructureHypothesis]
    reasoning_summary: str


# ─────────────────────────────────────────────────────────────────
# 结构解析引擎
# ─────────────────────────────────────────────────────────────────


class StructureElucidationEngine:
    """
    结构解析推理引擎

    从MS/MS谱图逐步推理化合物结构，生成可解释的结构假设。

    Usage:
        engine = StructureElucidationEngine()
        report = engine.elucidate(peaks, precursor_mz, ionization_mode)
        print(report.reasoning_summary)
    """

    def __init__(
        self,
        fragment_tolerance: float = 0.02,
        neutral_loss_tolerance: float = 0.5,
        precursor_tolerance_ppm: float = 10.0,
    ):
        self.frag_tol = fragment_tolerance
        self.nl_tol = neutral_loss_tolerance
        self.ppm_tol = precursor_tolerance_ppm
        self.elucidation_rules = get_structure_elucidation_rules()

    # ── 主接口 ──────────────────────────────────────────────────

    def elucidate(
        self,
        peaks: List[List[float]],
        precursor_mz: float,
        ionization_mode: str = "positive",
        candidate_smiles: Optional[str] = None,
        candidate_formula: Optional[str] = None,
    ) -> ElucidationReport:
        """
        执行完整的结构解析推理

        Args:
            peaks: [[mz, intensity], ...]
            precursor_mz: 前体离子 m/z
            ionization_mode: "positive" 或 "negative"
            candidate_smiles: 候选SMILES(可选)
            candidate_formula: 候选分子式(可选)

        Returns:
            完整解析报告
        """
        # 步骤1: 前体离子分析
        adduct_analysis = self._analyze_precursor(precursor_mz, ionization_mode)

        # 步骤2: 中性丢失与官能团诊断
        nl_evidence = self._analyze_neutral_losses(peaks, precursor_mz)

        # 步骤3: 诊断离子与化合物分类
        di_evidence, compound_class, class_conf = self._classify_from_diagnostic_ions(
            peaks
        )

        # 步骤4: 构建裂解树
        frag_tree = self._build_fragmentation_tree(peaks, precursor_mz, compound_class)

        # 步骤5: 子结构证据汇总
        substructure_evidence = self._collect_substructure_evidence(
            peaks, precursor_mz, compound_class
        )

        # 步骤6: 生成候选假设
        hypotheses = self._generate_hypotheses(
            peaks, precursor_mz, compound_class,
            nl_evidence, di_evidence, substructure_evidence,
            candidate_smiles,
        )

        # 步骤7: 生成推理摘要
        best = hypotheses[0] if hypotheses else None
        summary = self._generate_summary(
            adduct_analysis, compound_class, class_conf,
            nl_evidence, di_evidence, substructure_evidence,
            best,
        )

        return ElucidationReport(
            precursor_mz=precursor_mz,
            adduct_type=adduct_analysis.get("most_likely_adduct", "unknown"),
            neutral_mass=adduct_analysis.get("neutral_mass", 0.0),
            inferred_formula=candidate_formula or "",
            inferred_compound_class=compound_class,
            class_confidence=class_conf,
            diagnostic_evidence=substructure_evidence,
            fragmentation_tree=frag_tree,
            structure_hypotheses=sorted(
                hypotheses, key=lambda h: h.expert_score, reverse=True
            ),
            best_hypothesis=best,
            reasoning_summary=summary,
        )

    # ── 步骤1: 前体离子分析 ─────────────────────────────────────

    def _analyze_precursor(
        self, precursor_mz: float, ionization_mode: str
    ) -> Dict[str, Any]:
        """分析前体离子 → 推断加合物类型和中性质量"""
        candidates = []

        for adduct_name, profile in ADDUCT_PROFILES.items():
            if profile.polarity == ionization_mode or profile.polarity == "both":
                neutral_mass = precursor_mz - profile.mass_shift
                candidates.append({
                    "adduct": adduct_name,
                    "mass_shift": profile.mass_shift,
                    "neutral_mass": round(neutral_mass, 4),
                    "prevalence": profile.prevalence,
                    "typical_compounds": profile.typical_compounds,
                })

        # 按普遍性排序
        candidates.sort(key=lambda x: x["prevalence"], reverse=True)

        return {
            "precursor_mz": precursor_mz,
            "ionization_mode": ionization_mode,
            "candidates": candidates[:5],
            "most_likely_adduct": candidates[0]["adduct"] if candidates else "unknown",
            "neutral_mass": candidates[0]["neutral_mass"] if candidates else 0.0,
        }

    # ── 步骤2: 中性丢失分析 ─────────────────────────────────────

    def _analyze_neutral_losses(
        self, peaks: List[List[float]], precursor_mz: float
    ) -> List[ElucidationEvidence]:
        """从中性丢失推断官能团"""
        evidence_list = []
        fragment_mzs = [p[0] for p in peaks if p[0] < precursor_mz - 1.5]

        # v3.1 修复: 不传递 classification 给中性丢失分析
        # 因为此时类别尚未确定，使用保守的未知分类过滤
        nl_pairs = compute_possible_neutral_losses(
            precursor_mz, fragment_mzs, self.nl_tol,
            classification="", strict_fg=True,
        )

        # 官能团-中性丢失映射
        fg_nl_map = {
            "carboxyl": [43.9898],        # CO2
            "alcohol": [18.0106],          # H2O
            "phenol": [18.0106, 43.9898],  # H2O, CO2
            "methoxy": [15.0235, 32.0262], # CH3, CH3OH
            "amine": [17.0265],            # NH3
            "carbonyl": [27.9949],         # CO
            "glycoside": [162.0528, 146.0579, 132.0423],  # 糖基
            "ester": [60.0211, 32.0262],   # CH3COOH, CH3OH
            "methylenedioxy": [30.0106],   # OCH2
        }

        detected_fgs = set()
        for diff, nl in nl_pairs:
            if nl is None:
                continue
            for fg, masses in fg_nl_map.items():
                for expected_mass in masses:
                    if abs(diff - expected_mass) <= self.nl_tol:
                        detected_fgs.add(fg)
                        break

        for fg in detected_fgs:
            evidence_list.append(ElucidationEvidence(
                evidence_type="neutral_loss_diagnosis",
                description=f"检测到 {fg} 官能团的特征中性丢失",
                supporting_mz=[
                    round(precursor_mz - diff, 4)
                    for diff, nl in nl_pairs
                    if (nl and any(
                        abs(diff - m) <= self.nl_tol
                        for m in fg_nl_map.get(fg, [])
                    ))
                ],
                confidence=0.85,
                source="中性丢失分析",
            ))

        return evidence_list

    # ── 步骤3: 诊断离子与分类 ──────────────────────────────────

    def _classify_from_diagnostic_ions(
        self, peaks: List[List[float]]
    ) -> Tuple[List[ElucidationEvidence], str, float]:
        """通过诊断离子对化合物进行分类

        v3.1 修复: 对诊断离子匹配结果进行跨类过滤。
        仅当推断出的化合物类别确定后，才保留与该类别兼容的证据。
        """
        fragment_mzs = [p[0] for p in peaks]
        class_scores = classify_compound_from_fragments(fragment_mzs, self.frag_tol)

        # 先确定最佳类别
        if class_scores:
            best_class = max(class_scores, key=class_scores.get)
            best_conf = class_scores[best_class]
        else:
            best_class = "unknown"
            best_conf = 0.0

        evidence_list = []
        for mz_i in fragment_mzs:
            matches = match_fragment_to_class(mz_i, self.frag_tol)
            for match in matches:
                # v3.1 防御性过滤: 仅保留与推断类别一致的证据
                match_class = match.get("class_name", "unknown")
                if best_class != "unknown" and match_class != best_class:
                    # 跨类匹配: 保留但降权并标注警告
                    if match["diagnostic_value"] >= 0.85:
                        evidence_list.append(ElucidationEvidence(
                            evidence_type="diagnostic_ion",
                            description=(
                                f"m/z {mz_i:.4f} 匹配 {match['class_name_cn']} "
                                f"特征碎片: {match['fragment'].name} "
                                f"(诊断值: {match['diagnostic_value']:.2f}) "
                                f"⚠️ 跨类匹配: 推断类别={best_class}, 碎片类别={match_class}"
                            ),
                            supporting_mz=[mz_i],
                            confidence=match["diagnostic_value"] * 0.3,  # 跨类降权
                            source=f"碎片特征标准库 — {match['class_name_cn']}档案(跨类)",
                            alternative_explanation=f"该碎片也可能属于{match_class}类化合物",
                        ))
                    continue

                if match["diagnostic_value"] >= 0.85:
                    evidence_list.append(ElucidationEvidence(
                        evidence_type="diagnostic_ion",
                        description=(
                            f"m/z {mz_i:.4f} 匹配 {match['class_name_cn']} "
                            f"特征碎片: {match['fragment'].name} "
                            f"(诊断值: {match['diagnostic_value']:.2f})"
                        ),
                        supporting_mz=[mz_i],
                        confidence=match["diagnostic_value"],
                        source=f"碎片特征标准库 — {match['class_name_cn']}档案",
                    ))

        # 最佳分类
        if class_scores:
            best_class = max(class_scores, key=class_scores.get)
            best_conf = class_scores[best_class]
        else:
            best_class = "unknown"
            best_conf = 0.0

        return evidence_list, best_class, best_conf

    # ── 步骤4: 构建裂解树 ──────────────────────────────────────

    def _build_fragmentation_tree(
        self,
        peaks: List[List[float]],
        precursor_mz: float,
        compound_class: str,
    ) -> FragmentationPath:
        """构建从母离子到各碎片的裂解路径树"""
        fragment_mzs = [p[0] for p in peaks if p[0] < precursor_mz - 1.5]

        # 根节点 = 母离子
        root = FragmentationPath(
            mz=precursor_mz, parent_mz=None,
            neutral_loss=None, confidence=1.0,
        )

        # 获取该化合物类别的已知裂解路径
        pathway_desc = get_fragmentation_pathway(compound_class)

        # 建立碎片与母离子的连接
        sorted_frags = sorted(fragment_mzs, reverse=True)
        assigned = set()
        nodes_by_mz = {precursor_mz: root}

        for frag_mz in sorted_frags:
            # 寻找最近的合理父离子
            best_parent = None
            best_nl = None
            best_confidence = 0.0

            for parent_mz in sorted(
                [precursor_mz] + [
                    f for f in sorted_frags if f > frag_mz + 1.5
                ], reverse=True
            ):
                diff = parent_mz - frag_mz
                # v3.1 修复: 传递 compound_class 进行化学兼容性过滤
                matches = get_neutral_loss_by_mass(
                    diff, self.nl_tol,
                    classification=compound_class, strict_fg=True,
                )
                if matches:
                    best_parent = parent_mz
                    best_nl = matches[0]
                    best_confidence = matches[0].priority / 10.0
                    break
                elif diff >= 1.0:
                    # 近匹配 (放宽容差但保持分类过滤)
                    candidates = get_neutral_loss_by_mass(
                        diff, self.nl_tol * 2,
                        classification=compound_class, strict_fg=True,
                    )
                    if candidates:
                        best_parent = parent_mz
                        best_nl = candidates[0]
                        best_confidence = 0.4
                        break

            if best_parent is not None:
                node = FragmentationPath(
                    mz=frag_mz, parent_mz=best_parent,
                    neutral_loss=best_nl.name if best_nl else None,
                    neutral_loss_mass=round(best_parent - frag_mz, 4),
                    confidence=best_confidence,
                )
                nodes_by_mz[frag_mz] = node

                if best_parent in nodes_by_mz:
                    nodes_by_mz[best_parent].children.append(node)

                assigned.add(frag_mz)

        return root

    # ── 步骤5: 子结构证据收集 ──────────────────────────────────

    def _collect_substructure_evidence(
        self,
        peaks: List[List[float]],
        precursor_mz: float,
        compound_class: str,
    ) -> List[ElucidationEvidence]:
        """从谱图中收集所有子结构证据"""
        evidence = []
        fragment_mzs = [p[0] for p in peaks]

        rules = self.elucidation_rules.get("substructure_diagnosis", {})

        for substructure, info in rules.items():
            matched_evidence = []
            for evidence_pattern in info.get("evidence", []):
                # 简单匹配: 检查是否有对应的m/z或中性丢失
                if "m/z" in evidence_pattern:
                    # 提取m/z值
                    import re
                    mz_values = re.findall(r'm/z (\d+\.?\d*)', evidence_pattern)
                    for expected_mz_str in mz_values:
                        expected_mz = float(expected_mz_str)
                        for obs_mz in fragment_mzs:
                            if abs(obs_mz - expected_mz) <= self.frag_tol:
                                matched_evidence.append(obs_mz)
                                break

                # 检查中性丢失
                if "脱" in evidence_pattern or "丢失" in evidence_pattern:
                    nl_values = re.findall(r'(\d+)Da', evidence_pattern)
                    for nl_str in nl_values:
                        nl_mass = float(nl_str)
                        for frag_mz in fragment_mzs:
                            diff = precursor_mz - frag_mz
                            if abs(diff - nl_mass) <= self.nl_tol:
                                matched_evidence.append(frag_mz)
                                break

            if matched_evidence:
                evidence.append(ElucidationEvidence(
                    evidence_type="substructure_diagnosis",
                    description=(
                        f"检测到{substructure}的证据: {info.get('evidence', [])[0] if info.get('evidence') else ''}"
                    ),
                    supporting_mz=matched_evidence,
                    confidence=info.get("confidence", 0.8),
                    source="结构解析推理规则库",
                ))

        return evidence

    # ── 步骤6: 生成结构假设 ────────────────────────────────────

    def _generate_hypotheses(
        self,
        peaks: List[List[float]],
        precursor_mz: float,
        compound_class: str,
        nl_evidence: List[ElucidationEvidence],
        di_evidence: List[ElucidationEvidence],
        substructure_evidence: List[ElucidationEvidence],
        candidate_smiles: Optional[str] = None,
    ) -> List[StructureHypothesis]:
        """综合所有证据生成结构假设"""
        hypotheses = []

        # 如果提供了候选SMILES，构建假设
        if candidate_smiles:
            all_evidence = nl_evidence + di_evidence + substructure_evidence
            confidence_sum = sum(e.confidence for e in all_evidence)

            # 基础专家评分
            base_score = min(0.9, confidence_sum / max(len(all_evidence), 1))

            # 裂解路径一致性加分
            profile = get_profile(compound_class)
            if profile:
                base_score += 0.1 * profile.confidence_score

            hypotheses.append(StructureHypothesis(
                smiles=candidate_smiles,
                compound_class=compound_class,
                substructure_evidence=all_evidence,
                fragmentation_path=FragmentationPath(
                    mz=precursor_mz, parent_mz=None,
                    neutral_loss=None, confidence=1.0,
                ),
                expert_score=min(1.0, base_score),
                ranking=1,
            ))

        # 即没有候选SMILES，也生成基于证据的类别假设
        if not hypotheses and compound_class != "unknown":
            profile = get_profile(compound_class)
            if profile and profile.reference_compounds:
                for ref_compound in profile.reference_compounds[:3]:
                    hypotheses.append(StructureHypothesis(
                        smiles="(需要从数据库中获取)",
                        compound_class=compound_class,
                        substructure_evidence=di_evidence,
                        fragmentation_path=FragmentationPath(
                            mz=precursor_mz, parent_mz=None,
                            neutral_loss=None, confidence=1.0,
                        ),
                        expert_score=0.5,
                        ranking=len(hypotheses) + 1,
                    ))

        return hypotheses

    # ── 步骤7: 推理摘要 ────────────────────────────────────────

    def _generate_summary(
        self,
        adduct_analysis: Dict[str, Any],
        compound_class: str,
        class_conf: float,
        nl_evidence: List[ElucidationEvidence],
        di_evidence: List[ElucidationEvidence],
        substructure_evidence: List[ElucidationEvidence],
        best_hypothesis: Optional[StructureHypothesis],
    ) -> str:
        """生成人类可读的推理摘要"""
        parts = []

        # 前体分析
        adduct = adduct_analysis.get("most_likely_adduct", "unknown")
        neutral_mass = adduct_analysis.get("neutral_mass", 0)
        parts.append(f"【前体离子分析】")
        parts.append(f"  观测前体m/z: {adduct_analysis.get('precursor_mz', 0):.4f}")
        parts.append(f"  最可能加合物: {adduct} → 中性质量: {neutral_mass:.4f} Da")

        # 化合物分类
        profile = get_profile(compound_class)
        class_cn = profile.class_name_cn if profile else compound_class
        parts.append(f"\n【化合物分类】")
        parts.append(f"  推断类别: {class_cn} (置信度: {class_conf:.0%})")
        if profile:
            parts.append(f"  预期裂解路径: {profile.fragmentation_pathway[:100]}...")

        # 官能团证据
        if nl_evidence:
            parts.append(f"\n【官能团证据】")
            for ev in nl_evidence[:5]:
                parts.append(f"  ✓ {ev.description} (置信度: {ev.confidence:.0%})")

        # 诊断离子
        if di_evidence:
            parts.append(f"\n【诊断离子】")
            for ev in di_evidence[:3]:
                parts.append(f"  ✓ {ev.description}")

        # 子结构
        if substructure_evidence:
            parts.append(f"\n【子结构诊断】")
            for ev in substructure_evidence[:5]:
                parts.append(f"  ✓ {ev.description}")

        # 最佳假设
        if best_hypothesis:
            parts.append(f"\n【最佳结构假设】")
            parts.append(f"  SMILES: {best_hypothesis.smiles}")
            parts.append(f"  类别: {best_hypothesis.compound_class}")
            parts.append(f"  专家评分: {best_hypothesis.expert_score:.3f}")
            parts.append(f"  证据数量: {len(best_hypothesis.substructure_evidence)}")

        # 异构体区分
        if compound_class != "unknown":
            isomer_rules = get_isomer_rules(compound_class)
            if isomer_rules:
                parts.append(f"\n【异构体区分规则】")
                for rule in isomer_rules[:3]:
                    parts.append(f"  → {rule}")

        return "\n".join(parts)


# ─────────────────────────────────────────────────────────────────
# 便捷函数
# ─────────────────────────────────────────────────────────────────

_engine_instance: Optional[StructureElucidationEngine] = None


def get_elucidation_engine() -> StructureElucidationEngine:
    """获取全局结构解析引擎"""
    global _engine_instance
    if _engine_instance is None:
        _engine_instance = StructureElucidationEngine()
    return _engine_instance


def elucidate_structure(
    peaks: List[List[float]],
    precursor_mz: float,
    ionization_mode: str = "positive",
    candidate_smiles: Optional[str] = None,
) -> ElucidationReport:
    """
    便捷函数: 从MS/MS谱图解析结构

    Usage:
        report = elucidate_structure(peaks, precursor_mz, "positive")
        print(report.reasoning_summary)
    """
    engine = get_elucidation_engine()
    return engine.elucidate(peaks, precursor_mz, ionization_mode, candidate_smiles)


def explain_match(
    query_peaks: List[List[float]],
    query_pmz: float,
    match_smiles: str,
    match_name: str = "",
    match_class: str = "Unknown",
) -> Dict[str, Any]:
    """
    解释一个匹配为何合理

    生成从化合物结构到观测碎片的一致性分析。
    """
    engine = get_elucidation_engine()
    report = engine.elucidate(
        query_peaks, query_pmz,
        candidate_smiles=match_smiles,
    )

    # 从特征库获取类别信息
    profile = get_profile(match_class.lower())

    explanation = {
        "match_compound": match_name,
        "match_smiles": match_smiles,
        "observed_precursor_mz": query_pmz,
        "inferred_adduct": report.adduct_type,
        "inferred_class": report.inferred_compound_class,
        "class_confidence": report.class_confidence,
        "expected_fragmentation_pathway": (
            profile.fragmentation_pathway if profile else ""
        ),
        "detected_characteristic_fragments": [
            {
                "mz": ev.supporting_mz[0] if ev.supporting_mz else 0,
                "description": ev.description,
                "confidence": ev.confidence,
            }
            for ev in report.diagnostic_evidence
        ],
        "substructure_evidence": [
            {
                "evidence": ev.description,
                "supporting_fragments": ev.supporting_mz,
            }
            for ev in report.structure_hypotheses[0].substructure_evidence
            if report.structure_hypotheses
        ],
        "isomer_considerations": (
            list(get_isomer_rules(match_class.lower()))
            if match_class != "Unknown" else []
        ),
        "reasoning_summary": report.reasoning_summary,
        "overall_plausibility": (
            report.structure_hypotheses[0].expert_score
            if report.structure_hypotheses else 0.0
        ),
    }

    return explanation
