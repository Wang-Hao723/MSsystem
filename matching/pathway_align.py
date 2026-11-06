"""P-A: 化合物专属裂解通路对齐 (per-compound pathway alignment)
================================================================
消费 Compounds.fragmentation_pathway 列 (scripts/mine_pathways.py 挖出的逐化合物
中性丢失/诊断离子, 来自该化合物自身多张库谱), 给候选打"该候选自身特征碎片 vs query
观测峰"的对齐分。这是 91 万行通路列首次被评分代码使用 (此前只有硬编码类级模板在用)。

与 expert/expert_v3 的区别: expert_v3 按 *类别名* 查通用诊断离子; 本模块用 *该候选
化合物自身* 的 mined 特征峰 → 化合物级特异性, 而非类级。

纯函数 + lru_cache 解析, 无副作用。该特征来自候选自身的馆藏谱，若挖矿时没有
严格排除本次 query/同源数据，会产生目标泄漏。因此生产只允许 shadow 计算；结果必须
明确标为 candidate-conditioned 且不得进入 expert_score、排序或置信度。

pathway_align_score 算法:
  1. 门控: pathway None 或特征峰全空 → 0 (覆盖 27.5% 空壳化合物)
  2. 建候选自身可信特征 m/z: diagnostic_ions.mz + neutral_losses.fragment_mz,
     经 recurring/key + freq>=min_support 过滤 (剔单/双谱噪声)
  3. 匹配: query 峰 abs(q - 特征mz) <= mz_tol (绝对 m/z, 不依赖 precursor→抗加合物差异)
  4. freq 加权召回: Σmatched_freq / Σall_freq
  5. 防稀疏通胀: score = recall * min(1, matched_count/2)  (单离子匹配至多 0.5,
     同 search.py:427-432 单碎片虚高惩罚的思路)
"""
from __future__ import annotations

import json
from functools import lru_cache
from typing import Any, Optional


def parse_pathway(pathway_json: Optional[str]) -> Optional[dict]:
    """惰性解析通路 JSON, lru_cache 缓存 (同化合物跨 query 复用, 91 万行库解析开销摊薄)。"""
    if not pathway_json:
        return None
    return _parse_cached(pathway_json)


@lru_cache(maxsize=8192)
def _parse_cached(pathway_json: str) -> Optional[dict]:
    try:
        obj = json.loads(pathway_json)
        return obj if isinstance(obj, dict) else None
    except (ValueError, TypeError):
        return None


def _trusted(ion: dict, require_recurring: bool) -> bool:
    """离子是否可信: recurring==True 或带 curated key (diag_key/nl_key 非空)。"""
    if require_recurring:
        if ion.get("recurring") is True:
            return True
        key = ion.get("diag_key") or ion.get("nl_key")
        return bool(key)
    return True


def pathway_align_score(
    query_peaks: list,
    query_pmz: float,
    pathway: Optional[dict],
    mz_tol: float = 0.05,
    min_support: float = 0.5,
    require_recurring: bool = True,
) -> tuple:
    """计算候选自身通路 vs query 观测峰的对齐分。

    Returns: (score ∈ [0,1], detail dict)
      detail = {score, recall, n_trusted, matched_count, matched_ions, matched_losses}
    """
    empty_detail = {
        "score": 0.0, "recall": 0.0, "n_trusted": 0, "matched_count": 0,
        "matched_ions": [], "matched_losses": [],
        "candidate_conditioned": True,
        "independent_evidence": False,
        "eligible_for_fusion": False,
        "runtime_influence": False,
        "score_contribution": 0.0,
        "status": "shadow_only_query_exclusion_not_proven",
    }
    if not pathway or not query_peaks:
        return 0.0, empty_detail

    diag = pathway.get("diagnostic_ions") or []
    losses = pathway.get("neutral_losses") or []

    # ── 建候选自身可信特征 m/z (绝对产物离子 m/z, freq 加权) ──
    # diagnostic_ions.mz 与 neutral_losses.fragment_mz 都是绝对 m/z, 统一处理。
    # 按 0.001 Da 桶去重 (诊断离子与其互补丢失产物落在同 m/z 不重复计), freq 求和。
    bucket: dict = {}  # round(mz,3) -> [freq_sum, mz, n]
    for ion in diag:
        try:
            mz = float(ion.get("mz"))
        except (TypeError, ValueError):
            continue
        if not (0 < mz < 5000):
            continue
        if not _trusted(ion, require_recurring):
            continue
        freq = float(ion.get("freq", 0.0)) or 0.0
        if freq < min_support:
            continue
        b = round(mz, 3)
        if b in bucket:
            bucket[b][0] += freq
            bucket[b][2] += 1
        else:
            bucket[b] = [freq, mz, 1]

    for nl in losses:
        # 中性丢失用其绝对产物离子 fragment_mz (挖矿时 = 前体 - 丢失); None 则跳过
        fmz = nl.get("fragment_mz")
        try:
            mz = float(fmz)
        except (TypeError, ValueError):
            continue
        if not (0 < mz < 5000):
            continue
        if not _trusted(nl, require_recurring):
            continue
        freq = float(nl.get("freq", 0.0)) or 0.0
        if freq < min_support:
            continue
        b = round(mz, 3)
        if b in bucket:
            bucket[b][0] += freq
            bucket[b][2] += 1
        else:
            bucket[b] = [freq, mz, 1]

    if not bucket:
        return 0.0, empty_detail

    q_mz = [float(p[0]) for p in query_peaks if p and len(p) >= 2]
    if not q_mz:
        return 0.0, empty_detail

    # ── 匹配 + freq 加权召回 ──
    matched_weight = 0.0
    total_weight = 0.0
    matched_count = 0
    matched_mzs = []
    for b, (freq, mz, _n) in bucket.items():
        total_weight += freq
        if any(abs(q - mz) <= mz_tol for q in q_mz):
            matched_weight += freq
            matched_count += 1
            matched_mzs.append(round(mz, 3))

    recall = (matched_weight / total_weight) if total_weight > 0 else 0.0
    # 防稀疏通胀: 单特征匹配至多 0.5, 两独立特征才达 1.0
    score = recall * min(1.0, matched_count / 2.0)
    score = min(1.0, score)

    detail = {
        "score": round(score, 4),
        "recall": round(recall, 4),
        "n_trusted": len(bucket),
        "matched_count": matched_count,
        "matched_ions": matched_mzs[:8],
        "matched_losses": [],
        "candidate_conditioned": True,
        "independent_evidence": False,
        "eligible_for_fusion": False,
        "runtime_influence": False,
        "score_contribution": 0.0,
        "status": "shadow_only_query_exclusion_not_proven",
    }
    return score, detail


__all__ = ["parse_pathway", "pathway_align_score"]
