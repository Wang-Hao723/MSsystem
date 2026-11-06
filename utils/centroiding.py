"""智能 peak picking / centroiding 模块

对 mzML 中的 MS2 谱图做防御性预处理, 把 profile 模式连续扫描点压成真正的质心峰。

⚠️ 2026-07-13 重大修正(Fix1+Fix2):
  诊断发现 33 标品对照里很多已知物(Magnoflorine 等)余弦反常偏低, 根因不在匹配算法,
  而在本模块污染了谱图:
    Fix1: 真实 mzML 全部已标记 MS:1000127(centroid spectrum), 但旧代码无视该标记、
          只靠 _is_profile_like 的 gap 占比启发式判断 → 富含碎片的真实质心谱图(>15% 的
          相邻 gap < 0.05)被误判为 profile 残余而错误合并。现由 mzml_io 按 mzML cvParam
          直接跳过本模块(见 _file_is_centroid)。
    Fix2: 旧 _centroid_merge 用单链接邻接聚类(相邻 gap >= merge_da 才断簇) → 100.00 /
          100.04 / 100.08 / 100.12 这种 0.04 间隔的"链"会被串成 1 个簇(span 0.12 Da),
          即把多个不同碎片串成了一个峰。改为"簇总跨度 < merge_da"断簇(complete-linkage),
          杜绝链式合并。

为什么需要它(2026-06-30 实测坐实):
  真实 mzML(orbitrap profile 残余)里单条谱图常有 300-1500 个"峰", 其中大量是同一个
  化学峰的连续 profile 采样点。下游 clean_peaks 只做 top-N 强度截断 + 5% 相对强度过滤,
  不合并近邻点 → top-N 槽被同一峰的多个 profile 点占满, 真正的诊断碎片被挤出匹配窗口,
  cosine 召回受损。

正确做法(仅对真 profile 数据): 把"总跨度 < merge_da"的连续点聚成一簇, 取**强度加权质心**
(mz*int/Σint)作为该化学峰的真实 m/z, 峰高取簇内最大强度(与 clean_peaks/cosine 强度口径一致)。

幂等性: 已质心化的谱图(或被 mzml_io 按 cvParam 判定 centroid 而跳过)不经合并; 即便误入,
      簇内已无 < merge_da 跨度的连续点 → 原样返回。
"""
import logging
import numpy as np

logger = logging.getLogger(__name__)

# ─── 算法版本号 (Fix3): 解析缓存据此失效, 改动 centroiding 即自动重解析, 避免新旧缓存混用 ───
# 历史: v1 = 单链接邻接合并(有链式串峰 bug); centv2_20260713 = cvParam 跳过 + 跨度断簇
CENTROID_ALGO_VERSION = "centv2_20260713"

# profile 残余判定与合并阈值(模块常量, 便于统一调参)
# merge_da: 一簇的总跨度(mz_max - mz_min)小于此值才视为同一化学峰的 profile 采样点
MERGE_DA = 0.05
# profile_ratio: 连续 gap < merge_da 的占比超过此值才判为 profile 数据(否则按已质心化跳过)
PROFILE_RATIO = 0.15


def _gap_profile_ratio(peaks: list, merge_da: float = MERGE_DA) -> float:
    """相邻 m/z gap < merge_da 的占比 — 越高越像 profile 残余。"""
    if len(peaks) < 2:
        return 0.0
    mz = np.array([p[0] for p in peaks], dtype=float)
    mz.sort()
    diffs = np.diff(mz)
    return float(np.mean(diffs < merge_da)) if diffs.size else 0.0


def _is_profile_like(peaks: list, merge_da: float = MERGE_DA,
                     profile_ratio: float = PROFILE_RATIO) -> bool:
    """检测谱图是否含 profile 模式残余。

    特征: 相邻 m/z 间距 < merge_da 的占比 > profile_ratio (默认 15%)。
    旧启发式(">30% 的 gap < 0.01 Da")对真实数据几乎不触发(实测 0.04% 命中),
    因真实 profile 残余 gap 多在 0.01-0.05 之间、且占比 ~50% 而非 <0.01。
    新阈值在 92862 条真实谱图上对 profile 密集谱可靠触发。
    """
    if len(peaks) < 10:
        return False
    return _gap_profile_ratio(peaks, merge_da) > profile_ratio


def _centroid_merge(peaks: list, merge_da: float = MERGE_DA) -> list:
    """把"总跨度 < merge_da"的连续点聚成一簇, 输出每簇的强度加权质心。

    ⚠️ Fix2 (2026-07-13): 改 complete-linkage / 跨度断簇, 杜绝链式串峰。
      旧实现按"相邻 gap >= merge_da 断簇"(single-linkage) → m/z 间隔 0.04 的长链
      (100.00/100.04/100.08/100.12) 会被串成 1 个 span 0.12 的簇, 即多个不同碎片
      被并为一个峰, 污染 cosine。
      新实现: 每簇起点 start_mz 固定, 持续吸收点直到 (mz_j - start_mz) >= merge_da 才断。
      → 任一簇的总跨度 < merge_da (0.05 Da), 单个化学峰的 profile 采样(典型 <0.03 Da)
      仍被正确合并, 而不同碎片绝不被串接。

    每簇质心 m/z = Σ(mz*int)/Σint, 峰高 = 簇内最大强度。
    幂等: 已质心化输入(无 < merge_da 跨度的连续点) → 每点自成一簇 → 原样返回。
    """
    if len(peaks) < 2:
        return [list(p) for p in peaks]
    pts = sorted(peaks, key=lambda p: float(p[0]))
    n = len(pts)
    out = []
    i = 0
    while i < n:
        start_mz = float(pts[i][0])
        mz_acc = 0.0
        int_acc = 0.0
        cint = 0.0
        j = i
        # 持续吸收: 只要当前点相对簇首的跨度 < merge_da 就并入同一簇 (complete-linkage)
        while j < n and (float(pts[j][0]) - start_mz) < merge_da:
            mj = float(pts[j][0])
            ij = float(pts[j][1])
            mz_acc += mj * ij
            int_acc += ij
            if ij > cint:
                cint = ij
            j += 1
        if int_acc > 0:
            out.append([mz_acc / int_acc, cint])
        else:
            # 全零强度簇: 退化取首点
            out.append([start_mz, 0.0])
        i = j  # 跳到下一簇起点(必前进, 防 j==i 死循环)
    return out


def smart_peak_pick(
    peaks: list,
    precursor_mz: float,
    merge_da: float = MERGE_DA,
    profile_ratio: float = PROFILE_RATIO,
) -> list:
    """对可能含 profile 残余的谱图做真正的质心化。

    Args:
        peaks: [[mz, intensity], ...] 原始峰列表
        precursor_mz: 前体离子 m/z (保留入参兼容, 当前质心化与前体无关)
        merge_da: 同一化学峰的 profile 采样点合并阈值(默认 0.05 Da)
        profile_ratio: 判定 profile 的 gap 占比阈值(默认 0.15)

    Returns:
        [[mz, intensity], ...] 质心化后的峰列表; 已质心化输入为 no-op 原样返回。
    """
    if not peaks or len(peaks) < 3:
        return peaks

    # 已质心化的谱图直接跳过(幂等保障)
    if not _is_profile_like(peaks, merge_da, profile_ratio):
        return peaks

    merged = _centroid_merge(peaks, merge_da)
    if len(merged) < len(peaks):
        logger.debug(
            "Profile-like spectrum centroided: %d raw points → %d peaks (pmz=%.2f, "
            "merge_da=%.3f, ratio=%.2f)",
            len(peaks), len(merged), precursor_mz, merge_da,
            _gap_profile_ratio(peaks, merge_da),
        )
    return merged
