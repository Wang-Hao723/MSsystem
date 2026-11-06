"""
专家评估日志系统 — Expert Evaluation Logger
=============================================

将专家评判结果详细记录到日志和结构化数据文件中，
方便后续回顾、审计和持续改进。

功能:
1. **结构化JSON日志** — 每次匹配的完整专家评估记录
2. **汇总统计日志** — 批量运行后的总体统计
3. **Markdown报告** — 人类可读的专家评估报告
4. **审计追踪** — 记录每次评估的时间、参数和结果

日志文件路径:
- logs/expert_evaluations.jsonl      — 逐条评估记录(JSONL格式)
- logs/expert_summary.json           — 批量汇总统计
- logs/expert_report.md              — Markdown专家报告
- logS/expert_audit.log              — 审计日志

Author: Expert System Module
Version: 2.0.0
"""

from __future__ import annotations

import json
import logging
import os
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────
# 配置
# ─────────────────────────────────────────────────────────────────

DEFAULT_LOG_DIR = Path(__file__).resolve().parent.parent / "logs"
EVAL_LOGFILE = "expert_evaluations.jsonl"
SUMMARY_LOGFILE = "expert_summary.json"
REPORT_LOGFILE = "expert_report.md"
AUDIT_LOGFILE = "expert_audit.log"


# ─────────────────────────────────────────────────────────────────
# 专家日志记录器
# ─────────────────────────────────────────────────────────────────


class ExpertLogger:
    """
    专家评估日志记录器

    Usage:
        elog = ExpertLogger(log_dir="logs")
        elog.log_evaluation(match_result, expert_scores, metrics)
        elog.write_summary()
        elog.generate_report()
    """

    def __init__(
        self,
        log_dir: Optional[str] = None,
        session_id: Optional[str] = None,
        enabled: bool = True,
    ):
        """
        Args:
            log_dir: 日志目录路径
            session_id: 会话ID(默认自动生成)
            enabled: 是否启用日志记录
        """
        self.log_dir = Path(log_dir) if log_dir else DEFAULT_LOG_DIR
        self.session_id = session_id or datetime.now().strftime("%Y%m%d_%H%M%S")
        self.enabled = enabled
        self.evaluations: List[Dict[str, Any]] = []
        self.start_time = time.time()

        if self.enabled:
            self.log_dir.mkdir(parents=True, exist_ok=True)
            self._setup_audit_log()
            self._log_session_start()

    def _setup_audit_log(self):
        """初始化审计日志"""
        audit_path = self.log_dir / AUDIT_LOGFILE
        audit_handler = logging.FileHandler(audit_path, encoding="utf-8")
        audit_handler.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        ))
        self.audit_logger = logging.getLogger(f"expert_audit_{self.session_id}")
        self.audit_logger.addHandler(audit_handler)
        self.audit_logger.setLevel(logging.DEBUG)
        self.audit_logger.propagate = False

    def _log_session_start(self):
        """记录会话开始"""
        self.audit_logger.info("=" * 60)
        self.audit_logger.info("Expert Evaluation Session Started")
        self.audit_logger.info(f"Session ID: {self.session_id}")
        self.audit_logger.info(f"Log Directory: {self.log_dir}")
        self.audit_logger.info(f"Timestamp: {datetime.now().isoformat()}")
        self.audit_logger.info("=" * 60)

    # ── 核心记录接口 ───────────────────────────────────────────────

    def log_evaluation(
        self,
        query_info: Dict[str, Any],
        match_result: Dict[str, Any],
        expert_scores: Dict[str, Any],
        expert_metrics: Optional[Dict[str, Any]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        """
        记录一次完整的专家评估

        Args:
            query_info: 查询信息 {spectrum_id, precursor_mz, num_peaks, ...}
            match_result: 匹配结果 {compound_name, smiles, cosine, final_score, ...}
            expert_scores: 专家评分 {expert_total_score, neutral_loss_score, ...}
            expert_metrics: 专家指标(可选)
            metadata: 额外元数据(可选)

        Returns:
            记录ID
        """
        record_id = f"eval_{self.session_id}_{len(self.evaluations):06d}"
        timestamp = datetime.now().isoformat()

        record = {
            "record_id": record_id,
            "session_id": self.session_id,
            "timestamp": timestamp,
            "query": {
                "spectrum_id": query_info.get("spectrum_id", "unknown"),
                "precursor_mz": query_info.get("precursor_mz"),
                "num_peaks": query_info.get("num_peaks"),
                "retention_time": query_info.get("retention_time"),
                "file_source": query_info.get("file_source", ""),
            },
            "match": {
                "compound_name": match_result.get("compound_name", ""),
                "smiles": match_result.get("smiles", ""),
                "inchikey": match_result.get("inchikey", ""),
                "formula": match_result.get("formula", ""),
                "classification": match_result.get("classification", "Unknown"),
                "cosine_score": match_result.get("cosine"),
                "final_score": match_result.get("final_score"),
                "mass_score": match_result.get("mass_score"),
                "fragment_score": match_result.get("fragment_score"),
                "intensity_score": match_result.get("intensity_score"),
                "mass_ppm": match_result.get("mass_ppm"),
            },
            "expert_evaluation": {
                "expert_total_score": expert_scores.get("expert_total_score"),
                "confidence_level": expert_scores.get("confidence_level"),
                "effective_classification": expert_scores.get(
                    "effective_classification", ""
                ),
                "sub_scores": {
                    "neutral_loss": expert_scores.get("neutral_loss_score"),
                    "diagnostic_ion": expert_scores.get("diagnostic_ion_score"),
                    "rule_coverage": expert_scores.get("rule_coverage_score"),
                    "fg_signature": expert_scores.get("fg_signature_score"),
                    "structure_fragment": expert_scores.get("structure_fragment_score"),
                    "bond_energy": expert_scores.get("bond_energy_score"),
                },
                "chemical_reasoning": expert_scores.get("chemical_reasoning", ""),
                "rule_coverage_details": expert_scores.get("rule_coverage_analysis", {}),
                "diagnostic_ions_found": [
                    ion.get("ion_name")
                    for ion in expert_scores.get("diagnostic_ion_analysis", {})
                                .get("matched_ions", [])
                ],
                "neutral_loss_analysis": expert_scores.get(
                    "neutral_loss_analysis", {}
                ),
            },
            "expert_metrics": expert_metrics or {},
            "metadata": metadata or {},
        }

        self.evaluations.append(record)

        # 写入JSONL
        if self.enabled:
            self._append_to_jsonl(record)
            self.audit_logger.info(
                f"Evaluation {record_id}: {match_result.get('compound_name', '?')[:50]} "
                f"— Expert={expert_scores.get('expert_total_score', 0):.3f} "
                f"(Cosine={match_result.get('cosine', 0):.3f})"
            )

        return record_id

    def _append_to_jsonl(self, record: Dict[str, Any]):
        """追加一条记录到JSONL文件"""
        eval_path = self.log_dir / EVAL_LOGFILE
        try:
            with open(eval_path, "a", encoding="utf-8") as f:
                f.write(json.dumps(record, ensure_ascii=False) + "\n")
        except Exception as e:
            self.audit_logger.error(f"Failed to write evaluation record: {e}")

    # ── 批量汇总 ───────────────────────────────────────────────────

    def write_summary(self) -> Dict[str, Any]:
        """
        生成并写入批量汇总统计

        Returns:
            汇总统计字典
        """
        if not self.evaluations:
            summary = {"status": "no_evaluations", "count": 0}
        else:
            summary = self._compute_summary()

        if self.enabled:
            summary_path = self.log_dir / SUMMARY_LOGFILE
            with open(summary_path, "w", encoding="utf-8") as f:
                json.dump(summary, f, ensure_ascii=False, indent=2)

            self.audit_logger.info(
                f"Summary written: {summary.get('count', 0)} evaluations "
                f"→ {summary_path}"
            )

        return summary

    def _compute_summary(self) -> Dict[str, Any]:
        """计算汇总统计"""
        n = len(self.evaluations)

        # 提取所有分数
        expert_totals = []
        cosine_scores = []
        final_scores = []
        confidence_levels = defaultdict(int)
        classification_dist = defaultdict(int)
        grade_dist = defaultdict(int)
        rule_coverage_rates = []
        diagnostic_ion_counts = []

        for e in self.evaluations:
            ee = e.get("expert_evaluation", {})
            mm = e.get("match", {})

            et = ee.get("expert_total_score", 0)
            if et is not None:
                expert_totals.append(et)

            cs = mm.get("cosine_score", 0)
            if cs is not None:
                cosine_scores.append(cs)

            fs = mm.get("final_score", 0)
            if fs is not None:
                final_scores.append(fs)

            # 置信度：从 expert_total_score 推算（防止 null）
            raw_conf = ee.get("confidence_level")
            if raw_conf and raw_conf not in ("None", "null", ""):
                confidence_levels[raw_conf] += 1
            else:
                # 从专家分数推算置信度
                et_val = et if et is not None else 0
                if et_val >= 0.70:
                    confidence_levels["High"] += 1
                elif et_val >= 0.45:
                    confidence_levels["Medium"] += 1
                elif et_val >= 0.25:
                    confidence_levels["Low"] += 1
                else:
                    confidence_levels["Insufficient"] += 1
            classification_dist[mm.get("classification", "Unknown")] += 1

            # 根据专家分数分级
            et_val = et if et is not None else 0
            if et_val >= 0.70:
                grade_dist["A"] += 1
            elif et_val >= 0.55:
                grade_dist["B"] += 1
            elif et_val >= 0.40:
                grade_dist["C"] += 1
            elif et_val >= 0.25:
                grade_dist["D"] += 1
            else:
                grade_dist["F"] += 1

            # 规则覆盖率
            rc = ee.get("rule_coverage_details", {})
            if rc:
                rule_coverage_rates.append(rc.get("coverage_rate", 0))

            di = len(ee.get("diagnostic_ions_found", []))
            diagnostic_ion_counts.append(di)

        import numpy as np

        return {
            "session_id": self.session_id,
            "timestamp": datetime.now().isoformat(),
            "count": n,
            "elapsed_seconds": round(time.time() - self.start_time, 1),
            "expert_scores": {
                "mean": round(float(np.mean(expert_totals)), 4) if expert_totals else None,
                "std": round(float(np.std(expert_totals)), 4) if expert_totals else None,
                "median": round(float(np.median(expert_totals)), 4) if expert_totals else None,
                "min": round(float(np.min(expert_totals)), 4) if expert_totals else None,
                "max": round(float(np.max(expert_totals)), 4) if expert_totals else None,
            },
            "cosine_scores": {
                "mean": round(float(np.mean(cosine_scores)), 4) if cosine_scores else None,
                "std": round(float(np.std(cosine_scores)), 4) if cosine_scores else None,
            },
            "final_scores": {
                "mean": round(float(np.mean(final_scores)), 4) if final_scores else None,
                "std": round(float(np.std(final_scores)), 4) if final_scores else None,
            },
            "expert_vs_cosine_correlation": round(
                float(np.corrcoef(expert_totals, cosine_scores)[0, 1]), 4
            ) if (len(expert_totals) >= 2
                   and len(cosine_scores) >= 2
                   and len(expert_totals) == len(cosine_scores)) else None,
            "confidence_distribution": dict(confidence_levels),
            "grade_distribution": dict(grade_dist),
            "classification_distribution": dict(classification_dist),
            "avg_rule_coverage_rate": round(
                float(np.mean(rule_coverage_rates)), 4
            ) if rule_coverage_rates else None,
            "avg_diagnostic_ions_per_match": round(
                float(np.mean(diagnostic_ion_counts)), 2
            ) if diagnostic_ion_counts else None,
            "top_performing_classifications": sorted(
                classification_dist.items(), key=lambda x: x[1], reverse=True
            )[:10],
        }

    # ── Markdown报告生成 ────────────────────────────────────────────

    def generate_report(self) -> str:
        """
        生成人类可读的Markdown专家评估报告

        Returns:
            报告内容(str)
        """
        summary = self.write_summary()

        report = self._build_report_content(summary)

        if self.enabled:
            report_path = self.log_dir / REPORT_LOGFILE
            with open(report_path, "w", encoding="utf-8") as f:
                f.write(report)
            self.audit_logger.info(f"Report generated: {report_path}")

        return report

    def _build_report_content(self, summary: Dict[str, Any]) -> str:
        """构建Markdown报告"""
        n = summary.get("count", 0)
        es = summary.get("expert_scores", {})
        cs = summary.get("cosine_scores", {})
        corr = summary.get("expert_vs_cosine_correlation")

        lines = [
            f"# 🧪 化学专家评估报告 — Expert Evaluation Report",
            f"",
            f"**会话ID**: `{self.session_id}`  ",
            f"**生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  ",
            f"**评估数量**: {n}  ",
            f"**总耗时**: {summary.get('elapsed_seconds', 0):.1f} 秒  ",
            f"",
            f"---",
            f"",
            f"## 📊 总体统计",
            f"",
            f"| 指标 | 均值 | 标准差 | 中位数 | 最小值 | 最大值 |",
            f"|------|------|--------|--------|--------|--------|",
            f"| **专家综合评分 (CES)** | {es.get('mean', 'N/A')} | {es.get('std', 'N/A')} | {es.get('median', 'N/A')} | {es.get('min', 'N/A')} | {es.get('max', 'N/A')} |",
            f"| **余弦相似度 (Cosine)** | {cs.get('mean', 'N/A')} | {cs.get('std', 'N/A')} | N/A | N/A | N/A |",
            f"| **综合得分 (Final)** | {summary.get('final_scores', {}).get('mean', 'N/A')} | {summary.get('final_scores', {}).get('std', 'N/A')} | N/A | N/A | N/A |",
            f"",
            f"### 专家评分 vs 余弦相似度",
            f"",
        ]

        if corr is not None:
            lines.append(f"- **相关系数 (Pearson r)**: {corr:.4f}")
            if abs(corr) < 0.5:
                lines.append(
                    f"- ⚠️ 专家评分与余弦相似度的相关性**较低** (|r| < 0.5)，"
                    f"说明专家评估提供了**独立于余弦相似度的额外信息**。"
                )
            elif abs(corr) < 0.8:
                lines.append(
                    f"- 📊 专家评分与余弦相似度**中等相关** (0.5 ≤ |r| < 0.8)，"
                    f"两者部分一致但仍提供了互补信息。"
                )
            else:
                lines.append(
                    f"- ✅ 专家评分与余弦相似度**高度相关** (|r| ≥ 0.8)，"
                    f"两者一致性高。"
                )

        lines.extend([
            f"",
            f"## 🏅 置信度分布",
            f"",
        ])
        for level in ["High", "Medium", "Low", "Insufficient"]:
            count = summary.get("confidence_distribution", {}).get(level, 0)
            pct = count / n * 100 if n > 0 else 0
            bar = "█" * int(pct / 5) + "░" * (20 - int(pct / 5))
            lines.append(f"- **{level}**: {count} ({pct:.1f}%) `{bar}`")

        lines.extend([
            f"",
            f"## 🎯 专家等级分布",
            f"",
        ])
        grade_desc = {
            "A": "优秀 — 化学推理极好",
            "B": "良好 — 化学推理较好",
            "C": "中等 — 化学推理一般",
            "D": "有限 — 化学推理较弱",
            "F": "差 — 化学推理不足",
        }
        for grade in ["A", "B", "C", "D", "F"]:
            count = summary.get("grade_distribution", {}).get(grade, 0)
            pct = count / n * 100 if n > 0 else 0
            bar = "█" * int(pct / 5) + "░" * (20 - int(pct / 5))
            lines.append(
                f"- **等级 {grade}** ({grade_desc[grade]}): "
                f"{count} ({pct:.1f}%) `{bar}`"
            )

        lines.extend([
            f"",
            f"## 🧬 化合物分类覆盖",
            f"",
        ])
        top_class = summary.get("top_performing_classifications", [])
        for cls, count in top_class:
            lines.append(f"- **{cls}**: {count} 次匹配")

        lines.extend([
            f"",
            f"## 📏 平均指标",
            f"",
            f"- **平均规则覆盖率**: {summary.get('avg_rule_coverage_rate', 'N/A')}",
            f"- **平均诊断离子数/匹配**: {summary.get('avg_diagnostic_ions_per_match', 'N/A')}",
            f"",
            f"---",
            f"",
            f"### 知识库状态",
            f"",
        ])

        try:
            from .fragmentation_rules import get_rules_summary
            kb = get_rules_summary()
            lines.extend([
                f"- **裂解规则**: {kb['total_rules']} 条",
                f"- **中性丢失条目**: {kb['total_neutral_losses']} 条",
                f"- **诊断离子**: {kb['total_diagnostic_ions']} 个",
                f"- **BDE条目**: {kb['total_bde_entries']} 条",
                f"- **官能团覆盖**: {kb['functional_groups_covered']} 类",
                f"- **规则平均置信度**: {kb['rule_confidence_avg']:.2f}",
                f"",
                f"**机理类型分布**:",
            ])
            for mech, count in kb.get("mechanism_types", {}).items():
                lines.append(f"  - {mech}: {count} 条规则")
        except Exception:
            pass

        lines.extend([
            f"",
            f"---",
            f"",
            f"*报告由 ExpertLogger 自动生成 | 版本 2.0.0*",
        ])

        return "\n".join(lines)

    # ── 对比分析 ───────────────────────────────────────────────────

    def log_comparison(
        self,
        method_a: str,
        scores_a: List[float],
        method_b: str,
        scores_b: List[float],
        metric_name: str = "score",
    ):
        """
        记录两种方法的对比分析

        用于比较专家评分与传统方法的差异
        """
        import numpy as np

        a_arr = np.array(scores_a)
        b_arr = np.array(scores_b)

        diff = a_arr - b_arr
        corr = np.corrcoef(a_arr, b_arr)[0, 1] if len(a_arr) > 1 else 0

        comparison = {
            "timestamp": datetime.now().isoformat(),
            "method_a": method_a,
            "method_b": method_b,
            "metric": metric_name,
            "method_a_stats": {
                "mean": float(np.mean(a_arr)),
                "std": float(np.std(a_arr)),
                "min": float(np.min(a_arr)),
                "max": float(np.max(a_arr)),
            },
            "method_b_stats": {
                "mean": float(np.mean(b_arr)),
                "std": float(np.std(b_arr)),
                "min": float(np.min(b_arr)),
                "max": float(np.max(b_arr)),
            },
            "difference": {
                "mean": float(np.mean(diff)),
                "std": float(np.std(diff)),
                "a_better_count": int(np.sum(diff > 0)),
                "b_better_count": int(np.sum(diff < 0)),
                "tie_count": int(np.sum(np.abs(diff) < 1e-6)),
            },
            "correlation": float(corr),
            "interpretation": (
                f"{method_a} outperforms {method_b} by "
                f"{float(np.mean(diff)):.4f} on average"
                if float(np.mean(diff)) > 0
                else f"{method_b} outperforms {method_a} by "
                     f"{-float(np.mean(diff)):.4f} on average"
            ),
        }

        comp_path = self.log_dir / f"comparison_{self.session_id}.json"
        with open(comp_path, "w", encoding="utf-8") as f:
            json.dump(comparison, f, ensure_ascii=False, indent=2)

        return comparison

    # ── 清理与归档 ─────────────────────────────────────────────────

    def archive(self, archive_name: Optional[str] = None):
        """归档当前会话的日志"""
        name = archive_name or f"archive_{self.session_id}"
        archive_dir = self.log_dir / name
        archive_dir.mkdir(parents=True, exist_ok=True)

        import shutil

        # 移动日志文件
        for fname in [EVAL_LOGFILE, SUMMARY_LOGFILE, REPORT_LOGFILE, AUDIT_LOGFILE]:
            src = self.log_dir / fname
            if src.exists():
                shutil.move(str(src), str(archive_dir / fname))

        self.audit_logger.info(f"Session archived to: {archive_dir}")
        return str(archive_dir)

    def get_session_stats(self) -> Dict[str, Any]:
        """获取当前会话的快速统计"""
        if not self.evaluations:
            return {"count": 0}

        n = len(self.evaluations)
        expert_scores = [
            e.get("expert_evaluation", {}).get("expert_total_score", 0)
            for e in self.evaluations
        ]
        cosine_scores = [
            e.get("match", {}).get("cosine_score", 0)
            for e in self.evaluations
        ]

        import numpy as np

        return {
            "count": n,
            "avg_expert_score": float(np.mean(expert_scores)),
            "avg_cosine_score": float(np.mean(cosine_scores)),
            "high_confidence_ratio": sum(
                1 for e in self.evaluations
                if e.get("expert_evaluation", {}).get("confidence_level") == "High"
            ) / n if n > 0 else 0,
            "session_duration_min": round((time.time() - self.start_time) / 60, 1),
        }


# ─────────────────────────────────────────────────────────────────
# 便捷函数
# ─────────────────────────────────────────────────────────────────

_default_logger: Optional[ExpertLogger] = None


def get_expert_logger(
    log_dir: Optional[str] = None,
    enabled: bool = True,
) -> ExpertLogger:
    """获取全局专家日志记录器"""
    global _default_logger
    if _default_logger is None:
        _default_logger = ExpertLogger(log_dir=log_dir, enabled=enabled)
    return _default_logger


def log_match(
    query_info: Dict[str, Any],
    match_result: Dict[str, Any],
    expert_scores: Dict[str, Any],
    expert_metrics: Optional[Dict[str, Any]] = None,
) -> str:
    """便捷函数：记录一次匹配的专家评估"""
    elog = get_expert_logger()
    return elog.log_evaluation(query_info, match_result, expert_scores, expert_metrics)
