"""P-B: 诊断碎片抢救 (diagnostic-fragment rescue)
================================================
utils/spectrum.clean_peaks 在评分前硬过滤 rel>=0.05 + top15, 会砍掉低丰度但具诊断
价值的碎片。本模块在候选循环前, 对 *原始* query 谱抢救"匹配全局诊断离子集"的低强度
峰, 以强度地板注入 query_clean, 让下游 cosine/entropy/fragment/intensity/expert 可见。

纯增量 (USE_DIAGNOSTIC_RESCUE=False 时 search.py 完全不调用 = 当前行为)。

全局诊断 m/z 集: 来自 kb_v3 (12 类, ~134 个去重 m/z), 进程级建一次 (仿
matching/expert_v3.py:_ensure_idf)。匹配容差 0.05 Da (与 expert_v3 MZ_TOL 一致)。

rescue_diagnostic_peaks 返回 *待注入* 的峰列表 (非完整清洗谱); 调用方 (search.py)
负责去重 (丢弃 query_clean 已有 m/z) + 重排。强度地板 = base peak 的 intensity_floor
(默认 5% ≈ REL_THRESHOLD), 保证抢救峰可见但不主导。
"""
from __future__ import annotations

from typing import Optional

# 进程级全局诊断 m/z 集, 按极性分。仿 expert_v3._ION_DF 模式。
_GLOBAL_DIAG_MZ: dict = {}      # polarity -> set(round(mz,3))
_BUILT = False


def _ensure_diag_set() -> None:
    """遍历 kb_v3 全类 × 极性, 建 pos/neg/both 诊断 m/z 集。进程级只建一次。"""
    global _BUILT
    if _BUILT:
        return
    try:
        from kb_v3 import loader as kbv3
    except Exception:
        _BUILT = True
        return
    s = {"pos": set(), "neg": set()}
    try:
        for cls in kbv3.classes():
            for ion in kbv3.get_diagnostic_ions(cls, "pos"):
                try:
                    s["pos"].add(round(float(ion["mz"]), 3))
                except (TypeError, ValueError, KeyError):
                    pass
            for ion in kbv3.get_diagnostic_ions(cls, "neg"):
                try:
                    s["neg"].add(round(float(ion["mz"]), 3))
                except (TypeError, ValueError, KeyError):
                    pass
        s["both"] = s["pos"] | s["neg"]
    except Exception:
        s = {"pos": set(), "neg": set(), "both": set()}
    _GLOBAL_DIAG_MZ.update(s)
    _BUILT = True


def rescue_diagnostic_peaks(
    query_peaks: list,
    query_pmz: float,
    polarity: str = "both",
    mz_tol: float = 0.05,
    intensity_floor: float = 0.05,
) -> list:
    """抢救原始 query 谱中匹配全局诊断离子的低强度峰。

    Args:
        query_peaks: 原始 (clean_peaks 之前) 峰列表 [[mz,int],...]
        polarity: 'pos'/'neg'/'both'
        intensity_floor: 抢救峰相对强度地板 (base peak 的比例)

    Returns:
        待注入峰列表 [[mz, inten], ...] (强度已抬到地板)。空列表 = 无可抢救。
    """
    if not query_peaks:
        return []
    _ensure_diag_set()
    diag = _GLOBAL_DIAG_MZ.get(polarity) or _GLOBAL_DIAG_MZ.get("both")
    if not diag:
        return []

    # base peak (绝对强度最大) 用于相对地板
    try:
        max_int = max(float(p[1]) for p in query_peaks if p and len(p) >= 2 and p[1] is not None)
    except (ValueError, TypeError):
        return []
    if max_int <= 0:
        return []

    floor_int = intensity_floor * max_int
    rescued = []
    for p in query_peaks:
        if not p or len(p) < 2:
            continue
        try:
            mz = float(p[0])
            inten = float(p[1])
        except (TypeError, ValueError):
            continue
        if any(abs(mz - d) <= mz_tol for d in diag):
            # 命中全局诊断离子 → 即便低于 clean_peaks 的 rel/topN 阈值也抢救, 但不低于地板
            rescued.append([mz, max(inten, floor_int)])
    return rescued


def diag_set_sizes() -> dict:
    """调试/验证用: 返回各极性诊断 m/z 集大小。"""
    _ensure_diag_set()
    return {k: len(v) for k, v in _GLOBAL_DIAG_MZ.items()}


__all__ = ["rescue_diagnostic_peaks", "diag_set_sizes"]
