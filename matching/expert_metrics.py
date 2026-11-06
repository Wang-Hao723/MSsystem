"""
专家评估指标模块 — Expert Evaluation Metrics
=============================================

设计一套基于专家规则的评估指标，集成到训练和测试流程中。
与传统的 Cosine/MSE 不同，这些指标衡量模型是否真正"理解"了质谱裂解原理。

核心指标:
1. **碎裂规则覆盖率 (FRC)** — 模型预测的碎片能用已知裂解规则解释的比例
2. **结构推理得分 (SRS)** — 从碎片逆向推理结构的合理性
3. **官能团诊断得分 (FDS)** — 特征碎片离子的检测与归属
4. **中性丢失化学合理性 (NLR)** — 碎片间质量差的可解释性
5. **键能热力学一致性 (BEC)** — 推测断键位点的能量合理性

Author: Expert System Module
Version: 2.0.0
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

import numpy as np

from .fragmentation_rules import (
    FRAGMENTATION_RULES,
    DIAGNOSTIC_IONS,
    FUNCTIONAL_GROUP_SIGNATURES,
    COMMON_NEUTRAL_LOSSES,
    get_neutral_loss_by_mass,
    compute_possible_neutral_losses,
    compute_rule_coverage,
    match_diagnostic_ion,
    get_fg_signature,
    estimate_bond_lability,
    get_rules_summary,
)

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────
# 指标计算函数
# ─────────────────────────────────────────────────────────────────


def fragmentation_rule_coverage(
    fragment_mzs: List[float],
    precursor_mz: float,
    compound_classification: str,
    neutral_loss_tol: float = 0.5,
) -> Dict[str, Any]:
    """
    碎裂规则覆盖率 (Fragmentation Rule Coverage, FRC)

    衡量观测到的碎片模式能用已知裂解规则解释的程度。
    高FRC意味着碎片模式符合已知的化学原理。

    Returns:
        {
            "frc_score": float (0-1),
            "covered_rules": [rule_ids],
            "uncovered_rules": [rule_ids],
            "explanation_details": [...],
        }
    """
    if not fragment_mzs or precursor_mz <= 0:
        return {"frc_score": 0.0, "covered_rules": [], "uncovered_rules": [],
                "note": "insufficient_data"}

    # 计算中性丢失
    nl_pairs = compute_possible_neutral_losses(precursor_mz, fragment_mzs,
                                                 neutral_loss_tol)

    # 映射分类到官能团
    classification_map = {
        "flavonoid": "flavonoid", "alkaloid": "alkaloid",
        "terpenoid": "terpenoid", "terpene": "terpenoid",
        "phenolic_acid": "phenolic_acid", "phenolic": "phenolic_acid",
        "glycoside": "glycoside", "sugar": "glycoside",
        "fatty_acid": "fatty_acid", "lipid": "fatty_acid",
        "coumarin": "coumarin", "quinone": "quinone",
        "steroid": "terpenoid",
    }

    c_lower = compound_classification.lower() if compound_classification else ""
    fg = "unknown"
    for key, val in classification_map.items():
        if key in c_lower:
            fg = val
            break

    # 获取相关规则
    relevant_rules = []
    if fg != "unknown":
        sig = get_fg_signature(fg)
        if sig:
            relevant_rules = sig["rules"]
    else:
        # 通用评估：检查所有规则中是否有任何匹配
        relevant_rules = list(FRAGMENTATION_RULES.keys())

    # 计算覆盖
    coverage = compute_rule_coverage(nl_pairs, fg)

    # FRC得分计算
    if relevant_rules:
        # 加权因素1: 规则覆盖率
        raw_coverage = coverage.get("coverage_rate", 0.0)

        # 加权因素2: 匹配质量(匹配到的NL优先级)
        matched_nls = coverage.get("matched_neutral_losses", [])
        nl_quality = 0.0
        if matched_nls:
            for nl_name in matched_nls:
                nl = COMMON_NEUTRAL_LOSSES.get(nl_name)
                if nl:
                    nl_quality += nl.priority / 10.0
            nl_quality = min(1.0, nl_quality / max(len(matched_nls), 1))

        # 组合
        frc_score = raw_coverage * 0.7 + nl_quality * 0.3
    else:
        frc_score = 0.5  # 无相关规则时中性评分

    return {
        "frc_score": round(frc_score, 4),
        "covered_rules": coverage.get("covered_rules", []),
        "uncovered_rules": coverage.get("uncovered_rules", []),
        "total_expected_rules": coverage.get("total_rules", 0),
        "functional_group": fg,
        "matched_neutral_losses": matched_nls,
        "explanation_quality": nl_quality if matched_nls else 0.0,
    }


def structural_reasoning_score(
    smiles: str,
    fragment_mzs: List[float],
    precursor_mz: float,
    diag_ion_tol: float = 0.02,
) -> Dict[str, Any]:
    """
    结构推理得分 (Structural Reasoning Score, SRS)

    评估从碎片谱图逆向推理出合理分子结构的能力。
    包含三个子维度：
    - 碎片质量与可能子结构的一致性
    - 诊断离子与官能团的对应关系
    - 中性丢失与结构特征的对应

    Returns:
        {
            "srs_score": float (0-1),
            "substructure_match_count": int,
            "diagnostic_assignments": [...],
            "reasoning_path": str,
        }
    """
    if not smiles or not fragment_mzs:
        return {"srs_score": 0.0, "note": "insufficient_data"}

    # 1. 诊断离子分析
    diagnostic_assignments = []
    for mz in fragment_mzs:
        ions = match_diagnostic_ion(mz, diag_ion_tol)
        for ion in ions:
            diagnostic_assignments.append({
                "mz": mz,
                "ion_name": ion.name,
                "functional_group": ion.functional_group,
                "structural_origin": ion.structural_origin,
                "rule_id": ion.rule_id,
            })

    # 2. 化学键倾向性分析
    lability = estimate_bond_lability(smiles)
    labile_bonds = lability.get("labile_bonds", [])

    # 3. 结构-碎片匹配
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        mol = Chem.MolFromSmiles(smiles)
        has_valid_structure = mol is not None
        if has_valid_structure:
            mol_mass = Descriptors.ExactMolWt(mol)
            mass_match = abs(precursor_mz - (mol_mass + 1.0078)) < 0.5
        else:
            mol_mass = 0.0
            mass_match = False
    except ImportError:
        has_valid_structure = False
        mol_mass = 0.0
        mass_match = False

    # 4. 计算SRS
    # 子维度1: 诊断离子赋值(0-0.4)
    diag_contribution = min(0.4, len(diagnostic_assignments) * 0.15)

    # 子维度2: 结构可碎裂性(0-0.3)
    lability_contribution = min(0.3, len(labile_bonds) * 0.06)

    # 子维度3: 质量一致性(0-0.3)
    mass_contribution = 0.3 if mass_match else (0.15 if has_valid_structure else 0.0)

    srs_score = diag_contribution + lability_contribution + mass_contribution

    return {
        "srs_score": round(srs_score, 4),
        "diagnostic_assignments": diagnostic_assignments,
        "num_diagnostic_assignments": len(diagnostic_assignments),
        "labile_bonds_found": len(labile_bonds),
        "structure_valid": has_valid_structure,
        "mass_consistent": mass_match,
        "reasoning_path": (
            f"{len(diagnostic_assignments)} diagnostic ions assigned; "
            f"{len(labile_bonds)} labile bonds identified; "
            f"mass {'consistent' if mass_match else 'inconsistent'}"
        ),
    }


def functional_group_diagnostic_score(
    fragment_mzs: List[float],
    compound_classification: str,
    diag_ion_tol: float = 0.02,
    nl_tol: float = 0.5,
) -> Dict[str, Any]:
    """
    官能团诊断得分 (Functional Group Diagnostic Score, FDS)

    衡量实验谱图中是否检测到特定化合物类型预期的特征碎片。
    这是化合物类别确认的核心指标。

    Returns:
        {
            "fds_score": float (0-1),
            "expected_signature": {...},
            "detected_features": [...],
            "missed_features": [...],
        }
    """
    c_lower = compound_classification.lower() if compound_classification else ""
    fg = "unknown"
    classification_map = {
        "flavonoid": "flavonoid", "alkaloid": "alkaloid",
        "terpenoid": "terpenoid", "terpene": "terpenoid",
        "phenolic_acid": "phenolic_acid", "phenolic": "phenolic_acid",
        "glycoside": "glycoside", "sugar": "glycoside",
        "fatty_acid": "fatty_acid", "lipid": "fatty_acid",
        "coumarin": "coumarin", "quinone": "quinone",
        "steroid": "terpenoid",
    }
    for key, val in classification_map.items():
        if key in c_lower:
            fg = val
            break

    sig = get_fg_signature(fg)
    if not sig:
        return {"fds_score": 0.0, "note": f"no_signature_for_classification: {fg}"}

    # 检测特征离子
    expected_ions = sig.get("diagnostic_ions", [])
    detected_ions = []
    for mz in fragment_mzs:
        for expected in expected_ions:
            if abs(mz - expected) <= diag_ion_tol:
                detected_ions.append({
                    "observed_mz": mz,
                    "expected_mz": expected,
                    "delta": round(abs(mz - expected), 5),
                })

    missed_ions = [e for e in expected_ions
                   if not any(abs(mz - e) <= diag_ion_tol for mz in fragment_mzs)]

    # 检测特征中性丢失
    expected_losses = sig.get("characteristic_losses", [])
    detected_losses = []
    for mz in fragment_mzs:
        for expected in expected_losses:
            diff = abs(mz - expected)
            if diff <= nl_tol:
                detected_losses.append(expected)

    # 综合得分
    ion_score = len(detected_ions) / len(expected_ions) if expected_ions else 0.0
    loss_score = len(set(detected_losses)) / len(expected_losses) if expected_losses else 0.0

    fds_score = ion_score * 0.6 + loss_score * 0.4

    return {
        "fds_score": round(fds_score, 4),
        "functional_group": fg,
        "expected_signature": sig.get("fragmentation_pathway", ""),
        "detected_ions": detected_ions,
        "missed_ions": missed_ions,
        "detected_losses": list(set(detected_losses)),
        "expected_losses": expected_losses,
        "ion_detection_rate": ion_score,
        "loss_detection_rate": loss_score,
    }


def neutral_loss_rationality(
    fragment_mzs: List[float],
    precursor_mz: float,
    tolerance: float = 0.5,
) -> Dict[str, Any]:
    """
    中性丢失化学合理性 (Neutral Loss Rationality, NLR)

    评估碎片之间的质量差是否能对应合理的化学基团丢失。
    不合理的质量差(如罕见的自由基丢失)会降低此得分。

    Returns:
        {
            "nlr_score": float (0-1),
            "rational_losses": [...],
            "irrational_losses": [...],
            "overall_rationality": float (0-1),
        }
    """
    if not fragment_mzs or precursor_mz <= 0:
        return {"nlr_score": 0.0, "note": "insufficient_data"}

    nl_pairs = compute_possible_neutral_losses(precursor_mz, fragment_mzs,
                                                 tolerance)

    rational = []
    irrational = []

    for diff, nl in nl_pairs:
        if nl:
            rational.append({
                "mass_difference": round(diff, 4),
                "neutral_loss": nl.name,
                "formula": nl.formula,
                "mechanism": nl.mechanism,
                "priority": nl.priority,
            })
        else:
            # 检查是否对应任何已知的元素组合
            matches = get_neutral_loss_by_mass(diff, tolerance * 2)
            irrational.append({
                "mass_difference": round(diff, 4),
                "closest_matches": [
                    {"name": m.name, "mass": round(m.mass, 4),
                     "delta": round(abs(m.mass - diff), 4)}
                    for m in matches[:3]
                ] if matches else [],
                "is_common_radical": diff < 30,  # 小分子自由基通常不合理
            })

    total = len(nl_pairs)
    if total == 0:
        return {"nlr_score": 0.0, "note": "no_fragment_pairs"}

    # NLR = 合理解释的比例 + 高优先级NL的加权
    rational_rate = len(rational) / total

    # 加权：高优先级的中性丢失更可信
    weighted_bonus = 0.0
    if rational:
        priorities = [r["priority"] for r in rational]
        max_p = max(priorities) if priorities else 1
        weighted_bonus = sum(p / max_p for p in priorities) / len(priorities) * 0.2

    nlr_score = min(1.0, rational_rate + weighted_bonus)

    return {
        "nlr_score": round(nlr_score, 4),
        "rational_losses": rational,
        "irrational_losses": irrational,
        "rational_count": len(rational),
        "irrational_count": len(irrational),
        "overall_rationality": rational_rate,
        "weighted_bonus": round(weighted_bonus, 4),
    }


def bond_energy_consistency(
    smiles: str,
    fragment_mzs: List[float],
    precursor_mz: float,
) -> Dict[str, Any]:
    """
    键能热力学一致性 (Bond Energy Consistency, BEC)

    评估推测的断键位点是否具有化学合理性。
    基于键解离能(BDE)数据库的启发式估算。
    """
    if not smiles:
        return {"bec_score": 0.5, "note": "no_structure_provided"}

    lability = estimate_bond_lability(smiles)

    labile_bonds = lability.get("labile_bonds", [])
    fragmentation_sites = lability.get("estimated_fragmentation_sites", 0)

    # BEC = 存在易断裂键 + 合理碎片数
    if not labile_bonds:
        return {"bec_score": 0.3, "note": "no_labile_bonds_identified"}

    # 有多少高倾向性键
    high_lability = [b for b in labile_bonds if b["lability_score"] >= 0.7]

    # 预期的碎片数量(基于可断裂键数)
    expected_max_fragments = len(labile_bonds)
    observed_fragments = len(fragment_mzs)

    # 碎片数合理性: 观测碎片数不应远大于可断裂键数
    if expected_max_fragments > 0:
        fragment_ratio = min(1.0, observed_fragments / expected_max_fragments)
    else:
        fragment_ratio = 0.5

    # 组合得分
    lability_score = min(1.0, len(high_lability) / 3.0) if high_lability else 0.3
    bec_score = lability_score * 0.6 + fragment_ratio * 0.4

    return {
        "bec_score": round(bec_score, 4),
        "total_labile_bonds": len(labile_bonds),
        "high_lability_bonds": len(high_lability),
        "estimated_fragmentation_sites": fragmentation_sites,
        "observed_fragments": observed_fragments,
        "fragment_to_site_ratio": fragment_ratio,
        "most_labile_region": lability.get("most_labile_region", "unknown"),
        "top_bonds": labile_bonds[:3],
    }


# ─────────────────────────────────────────────────────────────────
# 综合专家评估指标
# ─────────────────────────────────────────────────────────────────

@dataclass
class ExpertMetricSet:
    """一次专家评估的所有指标结果"""
    frc: Dict[str, Any]   # 碎裂规则覆盖率
    srs: Dict[str, Any]   # 结构推理得分
    fds: Dict[str, Any]   # 官能团诊断得分
    nlr: Dict[str, Any]   # 中性丢失合理性
    bec: Dict[str, Any]   # 键能一致性

    @property
    def composite_expert_score(self) -> float:
        """综合专家评分 (CES - Composite Expert Score)"""
        weights = {"frc": 0.25, "srs": 0.25, "fds": 0.20, "nlr": 0.20, "bec": 0.10}
        scores = {
            "frc": self.frc.get("frc_score", 0),
            "srs": self.srs.get("srs_score", 0),
            "fds": self.fds.get("fds_score", 0),
            "nlr": self.nlr.get("nlr_score", 0),
            "bec": self.bec.get("bec_score", 0.5),
        }
        return sum(scores[k] * weights[k] for k in weights)

    @property
    def summary(self) -> Dict[str, Any]:
        """生成指标摘要"""
        return {
            "composite_expert_score": round(self.composite_expert_score, 4),
            "frc": round(self.frc.get("frc_score", 0), 4),
            "srs": round(self.srs.get("srs_score", 0), 4),
            "fds": round(self.fds.get("fds_score", 0), 4),
            "nlr": round(self.nlr.get("nlr_score", 0), 4),
            "bec": round(self.bec.get("bec_score", 0.5), 4),
            "grade": self._grade(),
            "covered_rules": self.frc.get("covered_rules", []),
            "diagnostic_assignments": [
                d.get("ion_name") for d in self.srs.get("diagnostic_assignments", [])
            ],
            "rational_loss_count": self.nlr.get("rational_count", 0),
        }

    def _grade(self) -> str:
        """基于CES的质量等级"""
        ces = self.composite_expert_score
        if ces >= 0.70:
            return "A — Excellent chemical reasoning"
        elif ces >= 0.55:
            return "B — Good chemical reasoning"
        elif ces >= 0.40:
            return "C — Moderate chemical reasoning"
        elif ces >= 0.25:
            return "D — Limited chemical reasoning"
        return "F — Poor chemical reasoning"


def compute_all_expert_metrics(
    smiles: str,
    fragment_mzs: List[float],
    precursor_mz: float,
    compound_classification: str,
    nl_tolerance: float = 0.5,
    diag_tolerance: float = 0.02,
) -> ExpertMetricSet:
    """
    计算所有五项专家评估指标

    Args:
        smiles: 化合物SMILES
        fragment_mzs: 碎片离子m/z列表
        precursor_mz: 前体离子m/z
        compound_classification: 化合物分类(如 Flavonoid, Alkaloid)
        nl_tolerance: 中性丢失容差(Da)
        diag_tolerance: 诊断离子容差(Da)

    Returns:
        ExpertMetricSet 包含所有指标结果
    """
    frc = fragmentation_rule_coverage(
        fragment_mzs, precursor_mz, compound_classification, nl_tolerance
    )
    srs = structural_reasoning_score(
        smiles, fragment_mzs, precursor_mz, diag_tolerance
    )
    fds = functional_group_diagnostic_score(
        fragment_mzs, compound_classification, diag_tolerance, nl_tolerance
    )
    nlr = neutral_loss_rationality(
        fragment_mzs, precursor_mz, nl_tolerance
    )
    bec = bond_energy_consistency(
        smiles, fragment_mzs, precursor_mz
    )

    return ExpertMetricSet(frc=frc, srs=srs, fds=fds, nlr=nlr, bec=bec)


# ─────────────────────────────────────────────────────────────────
# 批量评估与报告
# ─────────────────────────────────────────────────────────────────

def evaluate_model_performance(
    predictions: List[Dict[str, Any]],
    ground_truth: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    评估模型的专家级性能

    对模型的一组预测结果进行专家评估，生成性能报告。
    可与传统指标(Cosine, MSE)并行使用。

    Args:
        predictions: 预测结果列表，每项含 {'smiles', 'peaks': [(mz,int)], 'pmz', ...}
        ground_truth: 真值列表，每项含 {'smiles', 'classification', ...}

    Returns:
        性能报告
    """
    all_metrics = []
    for pred, gt in zip(predictions, ground_truth):
        frag_mzs = [p[0] for p in pred.get("peaks", [])]
        pmz = pred.get("pmz", 0)
        smi = pred.get("smiles", "")
        cls = gt.get("classification", "Unknown")

        metrics = compute_all_expert_metrics(smi, frag_mzs, pmz, cls)
        all_metrics.append(metrics)

    if not all_metrics:
        return {"error": "no_valid_predictions"}

    # 聚合统计
    ces_scores = [m.composite_expert_score for m in all_metrics]
    frc_scores = [m.frc.get("frc_score", 0) for m in all_metrics]
    srs_scores = [m.srs.get("srs_score", 0) for m in all_metrics]
    fds_scores = [m.fds.get("fds_score", 0) for m in all_metrics]
    nlr_scores = [m.nlr.get("nlr_score", 0) for m in all_metrics]
    bec_scores = [m.bec.get("bec_score", 0.5) for m in all_metrics]

    # 等级分布
    grade_counts = defaultdict(int)
    for m in all_metrics:
        grade_counts[m._grade()[0]] += 1

    return {
        "num_evaluated": len(all_metrics),
        "composite_expert_score": {
            "mean": np.mean(ces_scores),
            "std": np.std(ces_scores),
            "median": np.median(ces_scores),
            "min": np.min(ces_scores),
            "max": np.max(ces_scores),
        },
        "frc": {"mean": np.mean(frc_scores), "std": np.std(frc_scores)},
        "srs": {"mean": np.mean(srs_scores), "std": np.std(srs_scores)},
        "fds": {"mean": np.mean(fds_scores), "std": np.std(fds_scores)},
        "nlr": {"mean": np.mean(nlr_scores), "std": np.std(nlr_scores)},
        "bec": {"mean": np.mean(bec_scores), "std": np.std(bec_scores)},
        "grade_distribution": dict(grade_counts),
        "knowledge_base_summary": get_rules_summary(),
    }


# ─────────────────────────────────────────────────────────────────
# 与训练流程集成的回调函数
# ─────────────────────────────────────────────────────────────────

class ExpertMetricCallback:
    """
    训练回调：在每个epoch后计算专家指标

    可以作为 PyTorch Lightning callback 或手动调用。
    """

    def __init__(self, eval_frequency: int = 1):
        self.eval_frequency = eval_frequency
        self.history: List[Dict[str, Any]] = []

    def on_validation_epoch_end(
        self,
        epoch: int,
        val_predictions: List[Dict[str, Any]],
        val_ground_truth: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """在每个验证epoch结束后调用"""
        if epoch % self.eval_frequency != 0:
            return {}

        report = evaluate_model_performance(val_predictions, val_ground_truth)
        report["epoch"] = epoch
        self.history.append(report)

        logger.info(
            "Epoch %d — Expert CES: %.4f (FRC=%.4f, SRS=%.4f, FDS=%.4f)",
            epoch,
            report["composite_expert_score"]["mean"],
            report["frc"]["mean"],
            report["srs"]["mean"],
            report["fds"]["mean"],
        )

        return report

    def get_best_epoch(self) -> Dict[str, Any]:
        """获取CES最高的epoch"""
        if not self.history:
            return {}
        return max(self.history, key=lambda x: x["composite_expert_score"]["mean"])
