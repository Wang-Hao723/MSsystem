"""
质谱数据质量控制模块 — Data Quality Control (Spectraverse-inspired)
====================================================================

基于 Spectraverse (Gupta et al., Anal. Chem. 2026) 的数据预处理标准，
以及《质谱分析技术原理与应用》中介绍的数据解析最佳实践，
提供全面的质谱数据质量控制评估。

核心功能:
  1. 谱图级质量控制 (Spectrum-level QC)
  2. 元数据完整性检查 (Metadata completeness check)
  3. 化学结构-谱图一致性验证 (Structure-spectrum coherence)
  4. 近重复谱图检测 (Near-duplicate detection)
  5. 数据库整体质量评估 (Global database quality assessment)

Author: Enhanced Expert System Module
Version: 3.0.0
"""

from __future__ import annotations

import hashlib
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from .fragment_feature_library import DATA_QUALITY_STANDARDS
from .library_quality import normalize_inchikey

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────
# 数据结构
# ─────────────────────────────────────────────────────────────────


@dataclass
class QualityReport:
    """单张谱图的质量报告"""
    spectrum_id: str
    passed: bool
    quality_score: float           # 0-1
    flags: List[str]               # 质量问题标记
    warnings: List[str]            # 警告(不影响通过/失败)
    recommendations: List[str]     # 改进建议
    metrics: Dict[str, Any]        # 详细指标


@dataclass
class BatchQualityReport:
    """批量谱图质量报告"""
    total_spectra: int
    passed_count: int
    failed_count: int
    flag_distribution: Dict[str, int]
    quality_score_distribution: Dict[str, float]
    recommendations: List[str]


# ─────────────────────────────────────────────────────────────────
# 质量控制引擎
# ─────────────────────────────────────────────────────────────────


class DataQualityController:
    """
    质谱数据质量控制引擎

    实现 Spectraverse 论文中描述的预处理检查点。

    Usage:
        dqc = DataQualityController()
        report = dqc.check_spectrum(peaks, precursor_mz, metadata)
        if report.passed:
            # 使用此谱图
            ...
    """

    def __init__(self):
        self.standards = DATA_QUALITY_STANDARDS
        self.spec_std = self.standards.get("spectrum_level", {})
        self.meta_std = self.standards.get("metadata_level", {})
        self.struct_std = self.standards.get("chemical_structure_level", {})
        self.dup_std = self.standards.get("duplicate_detection", {})
        self.quality_flags = self.standards.get("quality_flags", {})

    # ── 单谱图质检 ──────────────────────────────────────────────

    def check_spectrum(
        self,
        peaks: List[List[float]],
        precursor_mz: float,
        metadata: Optional[Dict[str, Any]] = None,
        smiles: Optional[str] = None,
        spectrum_id: str = "",
    ) -> QualityReport:
        """
        对单张MS/MS谱图进行全面质量检查

        Args:
            peaks: [[mz, intensity], ...]
            precursor_mz: 前体离子 m/z
            metadata: {"ionization_mode": ..., "adduct": ..., "collision_energy": ...}
            smiles: 标注的化合物SMILES
            spectrum_id: 谱图标识符

        Returns:
            QualityReport
        """
        flags = []
        warnings = []
        recommendations = []
        metrics = {}

        # ── 1. 碎片基础检查 ──
        fragment_mzs = [p[0] for p in peaks]
        fragment_intensities = [p[1] for p in peaks]
        num_fragments = len(fragment_mzs)

        # 1.1 最少碎片数
        min_frags = self.spec_std.get("min_fragment_ions", 3)
        metrics["num_fragments"] = num_fragments
        if num_fragments < min_frags:
            flags.append(f"INSUFFICIENT_FRAGMENTS: {num_fragments} < {min_frags}")
            recommendations.append("需要至少3个碎片离子进行可靠匹配")

        # 1.2 m/z范围检查
        for mz in fragment_mzs:
            if mz < self.spec_std.get("min_mz", 10):
                flags.append(f"FRAGMENT_BELOW_MIN_MZ: m/z {mz:.4f} < 10")
                break

        # 1.3 碎片 > 前体m/z
        above_precursor = [mz for mz in fragment_mzs if mz > precursor_mz + 0.5]
        if above_precursor:
            flags.append(
                f"FRAGMENTS_ABOVE_PRECURSOR: {len(above_precursor)} fragments "
                f"have m/z greater than precursor {precursor_mz:.4f}"
            )
            recommendations.append("去除m/z大于前体离子的碎片(可能来自共流出或污染)")

        # 1.4 均匀强度检查
        if num_fragments >= 2:
            unique_intensities = set(round(i, 4) for i in fragment_intensities)
            if len(unique_intensities) == 1:
                flags.append("UNIFORM_INTENSITIES: all fragments have identical intensity")
                recommendations.append("均匀强度的谱图通常是伪谱(processing artifact)")

        # 1.5 零强度碎片
        zero_intensity = [i for i in fragment_intensities if i <= 0]
        if zero_intensity:
            flags.append(f"ZERO_INTENSITY_FRAGMENTS: {len(zero_intensity)} fragments")
            recommendations.append("去除零强度碎片离子")

        # 1.6 单一碎片近前体检查
        if num_fragments == 1:
            diff = abs(precursor_mz - fragment_mzs[0])
            if diff <= 1.6:
                flags.append(f"SINGLE_FRAGMENT_NEAR_PRECURSOR: diff={diff:.2f}Da")
                recommendations.append("仅有一个碎片且接近前体→无结构信息量")

        # 1.7 低分辨检查
        decimal_places = set()
        for mz in fragment_mzs:
            s = f"{mz:.6f}"
            if "." in s:
                decimals = len(s.split(".")[1].rstrip("0"))
                decimal_places.add(decimals)
        avg_decimal_places = np.mean(list(decimal_places)) if decimal_places else 0
        metrics["avg_decimal_places"] = avg_decimal_places
        if avg_decimal_places <= 2:
            flags.append(f"LOW_RESOLUTION: average {avg_decimal_places:.1f} decimal places")
            recommendations.append("低分辨数据: 考虑仅用于数据库搜索而非ML训练")

        # ── 2. 前体质量精度检查 ──

        if smiles:
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    exact_mass = Descriptors.ExactMolWt(mol)
                    theoretical_mz = exact_mass + 1.0078
                    ppm_error = abs(precursor_mz - theoretical_mz) / theoretical_mz * 1e6
                    metrics["precursor_ppm_error"] = ppm_error

                    ppm_threshold = self.spec_std.get("precursor_ppm_tolerance", 10.0)
                    if ppm_error > ppm_threshold:
                        flags.append(f"PRECURSOR_PPM_TOO_HIGH: {ppm_error:.1f}ppm > {ppm_threshold}ppm")
                else:
                    metrics["precursor_ppm_error"] = None
            except Exception:
                metrics["precursor_ppm_error"] = None
        else:
            metrics["precursor_ppm_error"] = None

        # ── 3. 元数据完整性 ──
        if metadata:
            required_fields = self.meta_std.get("required_fields", [])
            missing_required = [f for f in required_fields if f not in metadata or not metadata[f]]
            if missing_required:
                warnings.append(f"MISSING_REQUIRED_METADATA: {', '.join(missing_required)}")
                recommendations.append(f"补充必要元数据字段: {missing_required}")

            recommended_fields = self.meta_std.get("recommended_fields", [])
            missing_recommended = [
                f for f in recommended_fields if f not in metadata or not metadata[f]
            ]
            if missing_recommended:
                warnings.append(f"MISSING_RECOMMENDED_METADATA: {', '.join(missing_recommended)[:100]}")

        # ── 4. 计算质量得分 ──
        base_score = 1.0
        for flag in flags:
            if "INSUFFICIENT_FRAGMENTS" in flag:
                base_score -= 0.4
            elif "FRAGMENTS_ABOVE_PRECURSOR" in flag:
                base_score -= 0.15
            elif "UNIFORM_INTENSITIES" in flag:
                base_score -= 0.4
            elif "ZERO_INTENSITY_FRAGMENTS" in flag:
                base_score -= 0.1
            elif "LOW_RESOLUTION" in flag:
                base_score -= 0.2
            elif "SINGLE_FRAGMENT" in flag:
                base_score -= 0.3
            elif "PRECURSOR_PPM" in flag:
                base_score -= 0.25
            else:
                base_score -= 0.1

        quality_score = max(0.0, min(1.0, base_score))

        return QualityReport(
            spectrum_id=spectrum_id,
            passed=len(flags) == 0,
            quality_score=quality_score,
            flags=flags,
            warnings=warnings,
            recommendations=recommendations,
            metrics=metrics,
        )

    # ── 近重复检测 ──────────────────────────────────────────────

    def detect_near_duplicates(
        self,
        spectra: List[Dict[str, Any]],
        cosine_threshold: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        检测近重复谱图

        策略 (来自Spectraverse):
        1. 按完整规范 InChIKey + 极性分组
        2. 组内计算余弦相似度
        3. 超过阈值的视为近重复
        4. 保留碎片数最多的为代表

        Args:
            spectra: [{"peaks": [[mz,int],...], "inchikey": ..., "polarity": ..., "id": ...}, ...]
            cosine_threshold: 余弦阈值(默认0.99)

        Returns:
            {
                "duplicate_groups": [...],
                "representatives": [...],
                "removed_count": int,
            }
        """
        threshold = cosine_threshold or self.dup_std.get("cosine_threshold", 0.99)

        # 分组
        groups = defaultdict(list)
        for i, spec in enumerate(spectra):
            inchikey = normalize_inchikey(spec.get("inchikey", ""))
            polarity = spec.get("polarity", "positive")
            # Invalid/partial keys must not merge unrelated or stereochemically
            # distinct records.  Fall back to a spectrum-specific identity.
            identity = inchikey or f"spectrum:{spec.get('id', i)}:{i}"
            group_key = f"{identity}_{polarity}"
            groups[group_key].append(i)

        duplicate_groups = []
        representatives = set(range(len(spectra)))
        removed_count = 0

        for group_key, indices in groups.items():
            if len(indices) <= 1:
                continue

            # 组内计算配对余弦相似度
            for i in range(len(indices)):
                for j in range(i + 1, len(indices)):
                    idx1, idx2 = indices[i], indices[j]
                    cos_sim = self._compute_cosine(
                        spectra[idx1].get("peaks", []),
                        spectra[idx2].get("peaks", []),
                    )

                    if cos_sim >= threshold:
                        # 保留碎片数多的
                        frags1 = len(spectra[idx1].get("peaks", []))
                        frags2 = len(spectra[idx2].get("peaks", []))

                        if frags1 >= frags2:
                            keep, remove = idx1, idx2
                        else:
                            keep, remove = idx2, idx1

                        duplicate_groups.append({
                            "kept_id": spectra[keep].get("id", ""),
                            "removed_id": spectra[remove].get("id", ""),
                            "cosine_similarity": cos_sim,
                            "kept_fragments": frags1,
                            "removed_fragments": frags2,
                        })

                        if remove in representatives:
                            representatives.remove(remove)
                            removed_count += 1

        return {
            "duplicate_groups": duplicate_groups,
            "representative_indices": sorted(representatives),
            "removed_count": removed_count,
            "kept_count": len(representatives),
            "duplicate_rate": removed_count / len(spectra) if spectra else 0.0,
        }

    @staticmethod
    def _compute_cosine(
        peaks1: List[List[float]],
        peaks2: List[List[float]],
        mz_tolerance: float = 0.1,
    ) -> float:
        """计算两张谱图的余弦相似度"""
        if not peaks1 or not peaks2:
            return 0.0

        # 简单实现: 将峰合并到公共m/z轴,计算余弦
        all_mz = sorted(set(
            [round(p[0], 1) for p in peaks1] +
            [round(p[0], 1) for p in peaks2]
        ))

        vec1 = np.zeros(len(all_mz))
        vec2 = np.zeros(len(all_mz))

        for p in peaks1:
            idx = all_mz.index(round(p[0], 1))
            vec1[idx] = max(vec1[idx], p[1] if len(p) > 1 else 1.0)

        for p in peaks2:
            idx = all_mz.index(round(p[0], 1))
            vec2[idx] = max(vec2[idx], p[1] if len(p) > 1 else 1.0)

        # 归一化
        n1 = np.linalg.norm(vec1)
        n2 = np.linalg.norm(vec2)
        if n1 == 0 or n2 == 0:
            return 0.0

        return float(np.dot(vec1, vec2) / (n1 * n2))

    # ── 批量质检报告 ────────────────────────────────────────────

    def check_batch(
        self,
        spectra: List[Dict[str, Any]],
    ) -> BatchQualityReport:
        """
        批量谱图质量评估

        Args:
            spectra: 谱图列表

        Returns:
            批量质量报告
        """
        all_flags = defaultdict(int)
        quality_scores = []
        passed = 0
        failed = 0

        for spec in spectra:
            report = self.check_spectrum(
                peaks=spec.get("peaks", []),
                precursor_mz=spec.get("precursor_mz", 0),
                metadata=spec.get("metadata"),
                smiles=spec.get("smiles"),
                spectrum_id=spec.get("id", ""),
            )

            quality_scores.append(report.quality_score)
            if report.passed:
                passed += 1
            else:
                failed += 1

            for flag in report.flags:
                flag_type = flag.split(":")[0]
                all_flags[flag_type] += 1

        return BatchQualityReport(
            total_spectra=len(spectra),
            passed_count=passed,
            failed_count=failed,
            flag_distribution=dict(all_flags),
            quality_score_distribution={
                "mean": np.mean(quality_scores) if quality_scores else 0,
                "std": np.std(quality_scores) if quality_scores else 0,
                "min": np.min(quality_scores) if quality_scores else 0,
                "max": np.max(quality_scores) if quality_scores else 0,
                "median": np.median(quality_scores) if quality_scores else 0,
            },
            recommendations=self._generate_batch_recommendations(all_flags),
        )

    def _generate_batch_recommendations(
        self, flag_distribution: Dict[str, int]
    ) -> List[str]:
        """生成批量级别的改进建议"""
        recommendations = []

        # 按问题严重性排序
        if "LOW_RESOLUTION" in flag_distribution:
            count = flag_distribution["LOW_RESOLUTION"]
            pct = count / sum(flag_distribution.values()) * 100
            recommendations.append(
                f"低分辨谱图占{pct:.0f}% — 优先使用高分辨(QTOF/Orbitrap)数据用于ML训练 (Spectraverse标准)"
            )

        if "UNIFORM_INTENSITIES" in flag_distribution:
            recommendations.append(
                "检测到均匀强度伪谱 — 建议检查MGF转换过程,确保保留原始峰强度 (Spectraverse §Quality Filters)"
            )

        if "INSUFFICIENT_FRAGMENTS" in flag_distribution:
            recommendations.append(
                "部分谱图碎片数不足 — 可能需要调整碰撞能量以获取更丰富的碎片信息 (台湾质谱学会 §谱图解析)"
            )

        if "FRAGMENTS_ABOVE_PRECURSOR" in flag_distribution:
            recommendations.append(
                "存在碎片m/z > 前体m/z — 建议进行空白扣除和同位素峰过滤 (参考MS-Dial/GNPS策略)"
            )

        return recommendations


# ─────────────────────────────────────────────────────────────────
# 便捷函数
# ─────────────────────────────────────────────────────────────────

_dqc_instance: Optional[DataQualityController] = None


def get_quality_controller() -> DataQualityController:
    """获取全局质量控制引擎"""
    global _dqc_instance
    if _dqc_instance is None:
        _dqc_instance = DataQualityController()
    return _dqc_instance


def quick_quality_check(
    peaks: List[List[float]],
    precursor_mz: float,
    smiles: Optional[str] = None,
) -> Tuple[bool, float, List[str]]:
    """
    快速质量检查

    Returns:
        (passed, quality_score, flags)
    """
    dqc = get_quality_controller()
    report = dqc.check_spectrum(peaks, precursor_mz, smiles=smiles)
    return report.passed, report.quality_score, report.flags
