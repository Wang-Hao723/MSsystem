"""评分函数：余弦相似度、熵相似度、质量偏差、碎片重叠、强度分布 + 化学专家评分

v3.2 改进:
  - 修复 cosine_score 中库谱 precursor_mz 错误 Bug
  - 新增 entropy_score (基于 ms-entropy 库)
  - fragment_score 改用双向 Dice 系数
  - intensity_score 改用匹配碎片的 Pearson 相关
"""
import math
import numpy as np
from numba import njit
from matchms.similarity.spectrum_similarity_functions import (
    collect_peak_pairs,
    score_best_matches,
)
from config.settings import (
    WEIGHT_COSINE, WEIGHT_ENTROPY, WEIGHT_MASS, WEIGHT_FRAGMENT, WEIGHT_INTENSITY,
    WEIGHT_EXPERT, EXPERT_SYSTEM_ENABLED,
)


@njit(cache=True)
def _batch_modified_cosine_kernel(query_mz, query_intensity, query_pmz,
                                  library_mz, library_intensity, library_lengths,
                                  library_pmz, tolerance):
    """Batch equivalent of matchms ModifiedCosineGreedy for mz_power=0/intensity_power=1."""
    n_library = library_mz.shape[0]
    q_len = query_mz.shape[0]
    output = np.zeros(n_library, dtype=np.float64)
    q_norm_sq = 0.0
    for qi in range(q_len):
        q_norm_sq += query_intensity[qi] * query_intensity[qi]
    if q_norm_sq <= 0.0:
        return output
    q_norm = math.sqrt(q_norm_sq)

    for row in range(n_library):
        l_len = library_lengths[row]
        l_norm_sq = 0.0
        for li in range(l_len):
            l_norm_sq += library_intensity[row, li] * library_intensity[row, li]
        if l_norm_sq <= 0.0:
            continue

        max_pairs = max(1, 2 * q_len * l_len)
        pair_q = np.empty(max_pairs, dtype=np.int64)
        pair_l = np.empty(max_pairs, dtype=np.int64)
        products = np.empty(max_pairs, dtype=np.float64)
        count = 0
        mass_shift = query_pmz - library_pmz[row]

        # matchms concatenates zero-shift pairs first, then precursor-shifted pairs.
        for qi in range(q_len):
            for li in range(l_len):
                if abs(query_mz[qi] - library_mz[row, li]) <= tolerance:
                    pair_q[count] = qi
                    pair_l[count] = li
                    products[count] = query_intensity[qi] * library_intensity[row, li]
                    count += 1
        if abs(mass_shift) > tolerance:
            for qi in range(q_len):
                for li in range(l_len):
                    if abs(query_mz[qi] - (library_mz[row, li] + mass_shift)) <= tolerance:
                        pair_q[count] = qi
                        pair_l[count] = li
                        products[count] = query_intensity[qi] * library_intensity[row, li]
                        count += 1
        if count == 0:
            continue

        order = np.argsort(products[:count], kind="mergesort")
        used_q = np.zeros(q_len, dtype=np.uint8)
        used_l = np.zeros(l_len, dtype=np.uint8)
        score = 0.0
        for pos in range(count - 1, -1, -1):
            pair = order[pos]
            qi = pair_q[pair]
            li = pair_l[pair]
            if used_q[qi] == 0 and used_l[li] == 0:
                score += products[pair]
                used_q[qi] = 1
                used_l[li] = 1
        output[row] = score / (q_norm * math.sqrt(l_norm_sq))
    return output


def batch_modified_cosine_scores(query_peaks, library_peak_lists, query_pmz,
                                 library_pmzs, mz_tol=0.1):
    """Calculate exact ModifiedCosineGreedy scores in one numba boundary crossing."""
    if not query_peaks or not library_peak_lists:
        return np.zeros(len(library_peak_lists), dtype=float)
    query = np.asarray(query_peaks, dtype=float)
    query = query[np.argsort(query[:, 0], kind="mergesort")]
    lengths = np.asarray([len(peaks) for peaks in library_peak_lists], dtype=np.int64)
    max_length = int(lengths.max()) if lengths.size else 0
    if max_length == 0:
        return np.zeros(len(library_peak_lists), dtype=float)
    library_mz = np.zeros((len(library_peak_lists), max_length), dtype=np.float64)
    library_intensity = np.zeros_like(library_mz)
    for row, peaks in enumerate(library_peak_lists):
        if not peaks:
            continue
        array = np.asarray(peaks, dtype=float)
        array = array[np.argsort(array[:, 0], kind="mergesort")]
        library_mz[row, :len(array)] = array[:, 0]
        library_intensity[row, :len(array)] = array[:, 1]
    return _batch_modified_cosine_kernel(
        query[:, 0], query[:, 1], float(query_pmz),
        library_mz, library_intensity, lengths,
        np.asarray(library_pmzs, dtype=np.float64), float(mz_tol),
    )


def _neutral_loss_peaks(peaks, precursor_mz):
    """Project product ions to positive precursor-neutral-loss coordinates."""
    try:
        precursor = float(precursor_mz)
    except (TypeError, ValueError):
        return []
    if not math.isfinite(precursor) or precursor <= 0:
        return []
    output = []
    for peak in peaks or []:
        try:
            mz, intensity = float(peak[0]), float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        loss = precursor - mz
        if (
            math.isfinite(loss)
            and math.isfinite(intensity)
            and loss > 0
            and intensity > 0
        ):
            output.append([loss, intensity])
    output.sort(key=lambda item: item[0])
    return output


def batch_neutral_loss_cosine_scores(
    query_peaks,
    library_peak_lists,
    query_pmz,
    library_pmzs,
    mz_tol=0.1,
):
    """Batch equivalent of matchms ``NeutralLossesCosine``.

    Query and library peaks must already have passed the same preprocessing
    preset.  Transforming both sides to precursor-minus-product-ion coordinates
    makes this an ordinary greedy cosine without a further precursor shift.
    """
    query_losses = _neutral_loss_peaks(query_peaks, query_pmz)
    library_losses = [
        _neutral_loss_peaks(peaks, precursor)
        for peaks, precursor in zip(library_peak_lists, library_pmzs)
    ]
    return batch_modified_cosine_scores(
        query_losses,
        library_losses,
        0.0,
        [0.0] * len(library_losses),
        mz_tol=mz_tol,
    )


def cosine_score(query_peaks, lib_peaks, mz_tol=0.1, query_pmz=None, lib_pmz=None):
    """Modified Cosine 相似度 (0-1)

    修复: 使用实际的 lib_pmz 替代旧的 max(mz)+1.0 哑值。
    ModifiedCosine 使用 precursor_mz 来对齐碎片峰的 m/z shift，
    使用错误值会导致评分不准确。
    """
    if not query_peaks or not lib_peaks:
        return 0.0
    try:
        q_arr = np.asarray(query_peaks, dtype=float).copy()
        l_arr = np.asarray(lib_peaks, dtype=float).copy()
        q_arr = q_arr[np.argsort(q_arr[:, 0], kind="mergesort")]
        l_arr = l_arr[np.argsort(l_arr[:, 0], kind="mergesort")]

        q_ints = q_arr[:, 1]
        q_max = q_ints.max() if q_ints.size > 0 else 1.0
        q_arr[:, 1] = q_ints / q_max if q_max > 0 else q_ints

        l_ints = l_arr[:, 1]
        l_max = l_ints.max() if l_ints.size > 0 else 1.0
        l_arr[:, 1] = l_ints / l_max if l_max > 0 else l_ints

        # 使用实际 precursor_mz；缺省时回退到查询谱最大 m/z
        q_pmz = query_pmz if query_pmz else max(p[0] for p in query_peaks)
        l_pmz_val = lib_pmz if lib_pmz else max(p[0] for p in lib_peaks)

        # Exact ModifiedCosineGreedy calculation on arrays.  The public matchms
        # wrapper constructs Spectrum/Fragments objects for every library candidate;
        # its implementation then delegates to these same numba kernels.
        mass_shift = q_pmz - l_pmz_val
        zero_pairs = collect_peak_pairs(q_arr, l_arr, mz_tol, shift=0.0)
        pair_sets = [] if zero_pairs is None else [zero_pairs]
        if abs(mass_shift) > mz_tol:
            shifted_pairs = collect_peak_pairs(q_arr, l_arr, mz_tol, shift=mass_shift)
            if shifted_pairs is not None:
                pair_sets.append(shifted_pairs)
        if not pair_sets:
            return 0.0
        matching_pairs = pair_sets[0] if len(pair_sets) == 1 else np.concatenate(pair_sets, axis=0)
        matching_pairs = matching_pairs[
            np.argsort(matching_pairs[:, 2], kind="mergesort")[::-1], :
        ]
        val = float(score_best_matches(matching_pairs, q_arr, l_arr, 0.0, 1.0)[0])
        return float(val) if not math.isnan(val) else 0.0
    except Exception:
        return 0.0


def entropy_score(query_peaks, lib_peaks, mz_tol=0.1):
    """熵相似度 (0-1)，基于 ms-entropy 库

    Li et al. (2021) Nature Methods 证明熵相似度系统性优于余弦相似度。
    使用 Shannon 信息熵对 MS/MS 谱加权，自动增强低强度高信息量离子峰。
    """
    if not query_peaks or not lib_peaks:
        return 0.0
    try:
        import ms_entropy
        q = np.array(query_peaks, dtype=float)
        l = np.array(lib_peaks, dtype=float)
        sim = ms_entropy.calculate_entropy_similarity(
            q, l, ms2_tolerance_in_da=mz_tol)
        return float(sim) if not math.isnan(sim) else 0.0
    except ImportError:
        # ms-entropy 未安装时静默降级
        return 0.0
    except Exception:
        return 0.0


def mass_score(pmz, lib_pmz, ppm_tol=None, da_tol=0.5):
    """质量精度得分 (0-1) — 连续、单调的 ppm/Da 双容差。

    旧实现先让 ppm 分量在边界降到 0，再在边界外切换到约 0.45 的
    Da 分量，导致质量误差更大的候选反而得分更高。两个分量本身都
    随误差单调下降；取其最大值既保留低 m/z 的绝对容差兜底，也消除
    ppm/Da 切换处的向上跳变。
    """
    if pmz <= 0 or lib_pmz <= 0:
        return 0.0
    diff_da = abs(pmz - lib_pmz)
    diff_ppm = diff_da / pmz * 1e6

    if ppm_tol is None:
        ppm_tol = 20  # 放宽默认ppm容差

    ppm_component = 0.0
    if ppm_tol and ppm_tol > 0 and diff_ppm <= ppm_tol:
        ppm_component = max(0.0, 1.0 - diff_ppm / ppm_tol)

    da_component = 0.0
    if da_tol and da_tol > 0 and diff_da <= da_tol:
        da_component = max(0.1, 1.0 - diff_da / da_tol) * 0.6

    return max(ppm_component, da_component)


def mass_ppm(pmz, lib_pmz):
    """质量偏差 (ppm)"""
    if pmz <= 0 or lib_pmz <= 0:
        return None
    return abs(pmz - lib_pmz) / pmz * 1e6


def mass_error_penalty(ppm, free_tol=50.0, mid_tol=200.0, floor=0.35):
    """质量误差惩罚因子 ∈ [floor, 1.0]。1.0 = 无惩罚。

    解决"ppm 误差大但 final_score 仍高"的缺陷: final_score 由 cos(0.40)+entropy(0.30)主导,
    mass 仅 0.10 权重 → ppm=300 时 mass_score=0 也只扣 0.10, 分数仍虚高。
    本因子在 final 组装后直接乘入, 让质量明显不符的候选分数诚实下降 (cos 相似≠可信命中)。

    化学依据: 质谱鉴定前提是前体 m/z 匹配; ppm 越大, 命中越可能是错误化合物/错误加合物。
      |ppm| ≤ free_tol(50)    → 1.00  正常仪器质量精度(HRMS<5ppm, 低分辨<100ppm),
                                      自检索真匹配 ppm≈0 → 红线 11/14 完全不受影响。
      free_tol < |ppm| ≤ mid  → 1.00→0.75 线性衰减 (质量轻微不符, 可能校准漂移)
      |ppm| > mid(200)        → 0.75→floor(0.35) 更陡衰减 (质量明显不符, 即便谱图相似不可信)
    floor>0 保留: 仍可能是有效类似物/未知加合物, 且 honest_calibration 已将 confidence 降为 Low。
    """
    if ppm is None:
        return 1.0
    try:
        ap = abs(float(ppm))
    except (TypeError, ValueError):
        return 1.0
    if ap <= free_tol:
        return 1.0
    if ap <= mid_tol:
        return 1.0 - 0.25 * (ap - free_tol) / (mid_tol - free_tol)
    return max(floor, 0.75 - 0.40 * min(ap - mid_tol, 300.0) / 300.0)


def reverse_match_ratio(query_peaks, lib_peaks, mz_tol=0.1,
                        query_pmz=None, lib_pmz=None):
    """反向匹配率 (NIST reverse dot product 思路): lib 谱图的峰(强度加权)在 query 中
    找到匹配的比例。低 = lib 有很多 query 没有的特征峰 = 假阳性信号(如甾体骨架碎片
    通用导致 cos 虚高, 但 lib 的特征峰 query 没有)。

    仅用于 confidence 降级(不改 final_score), 零风险控制假阳性标 High 误导用户。
    与 ModifiedCosine 保持一致：跨加合物/类似物检索时，如果两个前体
    相差超过碎片容差，同时接受 precursor-shifted fragment pair。此前
    cosine 使用 shift、reverse 却只看零位移，会系统性压低跨加合物候选。

    参考: Spectral entropy / NIST MSMatch reverse dot product (假阳性控制标准方法)。
    """
    if not query_peaks or not lib_peaks:
        return 0.0
    try:
        qmz = np.array([p[0] for p in query_peaks], dtype=float)
        if qmz.size == 0:
            return 0.0
        mass_shift = None
        try:
            if query_pmz is not None and lib_pmz is not None:
                candidate_shift = float(query_pmz) - float(lib_pmz)
                if abs(candidate_shift) > mz_tol:
                    mass_shift = candidate_shift
        except (TypeError, ValueError):
            mass_shift = None
        matched_int = 0.0
        total_int = 0.0
        for lmz, lint in lib_peaks:
            la = abs(float(lint))
            total_int += la
            diff = np.abs(qmz - float(lmz))
            direct_match = bool(diff.size and diff.min() <= mz_tol)
            shifted_match = False
            if mass_shift is not None:
                shifted = np.abs(qmz - (float(lmz) + mass_shift))
                shifted_match = bool(shifted.size and shifted.min() <= mz_tol)
            if direct_match or shifted_match:
                matched_int += la
        return float(matched_int / total_int) if total_int > 0 else 0.0
    except Exception:
        return 0.0


def fragment_score(query_clean, lib_clean, mz_tol=0.1):
    """碎片重叠得分 (0-1) — 双向 Dice 系数

    改进: 使用双向匹配 + Dice 系数替代单向匹配率。
    Dice = 2 * |A ∩ B| / (|A| + |B|)，对两边碎片数量差异更鲁棒。
    mz_tol=0.1 与 cosine_score 对齐(旧 0.02 严 5 倍, 致匹配峰过少、分数系统性偏低)。
    """
    if not query_clean or not lib_clean:
        return 0.0
    q_mz = np.array([p[0] for p in query_clean], dtype=float)
    l_mz_arr = np.array([p[0] for p in lib_clean], dtype=float)
    if len(q_mz) == 0 or len(l_mz_arr) == 0:
        return 0.0

    # 双向匹配: 查询→库 和 库→查询
    matched_fwd = _count_matched_peaks(q_mz, l_mz_arr, mz_tol)
    matched_rev = _count_matched_peaks(l_mz_arr, q_mz, mz_tol)
    matched = max(matched_fwd, matched_rev)

    # Dice 系数
    total = len(q_mz) + len(l_mz_arr)
    if total == 0:
        return 0.0
    return min(1.0, 2.0 * matched / total)


def _count_matched_peaks(source_mz, target_mz, tol):
    """统计 source 中有多少峰能在 target 中找到匹配"""
    matched = 0
    used = set()
    for smz in source_mz:
        diffs = np.abs(target_mz - smz)
        sorted_idx = np.argsort(diffs)
        for idx in sorted_idx:
            if diffs[idx] <= tol and idx not in used:
                matched += 1
                used.add(idx)
                break
    return matched


def intensity_score(query_clean, lib_clean, mz_tol=0.1):
    """强度分布相似度 (0-1) — 基于 Pearson 相关

    改进: 对所有匹配碎片的强度做 Pearson 相关，替代旧的 Top-3 峰比较。
    匹配碎片 < 3 个时返回 0.0（样本不足，相关性不可靠）。
    mz_tol=0.1 与 cosine_score 对齐(旧 0.02 严 5 倍, 致半数候选匹配峰<3得0分)。
    """
    if not query_clean or not lib_clean:
        return 0.0

    # 收集匹配碎片的强度
    q_mz = np.array([p[0] for p in query_clean], dtype=float)
    l_mz = np.array([p[0] for p in lib_clean], dtype=float)
    q_ints_all = np.array([p[1] for p in query_clean], dtype=float)
    l_ints_all = np.array([p[1] for p in lib_clean], dtype=float)

    q_matched = []
    l_matched = []
    used_l = set()

    for i, qmz in enumerate(q_mz):
        diffs = np.abs(l_mz - qmz)
        sorted_idx = np.argsort(diffs)
        for idx in sorted_idx:
            if diffs[idx] <= mz_tol and idx not in used_l:
                q_matched.append(q_ints_all[i])
                l_matched.append(l_ints_all[idx])
                used_l.add(idx)
                break

    if len(q_matched) < 3:
        return 0.0

    q_arr = np.array(q_matched, dtype=float)
    l_arr = np.array(l_matched, dtype=float)

    # Pearson 相关
    q_mean = q_arr.mean()
    l_mean = l_arr.mean()
    q_c = q_arr - q_mean
    l_c = l_arr - l_mean
    denom = np.sqrt(float(np.sum(q_c ** 2)) * float(np.sum(l_c ** 2)))
    if denom < 1e-12:
        return 0.0

    pearson = float(np.sum(q_c * l_c)) / denom
    # 映射 [-1, 1] → [0, 1]
    return max(0.0, (pearson + 1.0) / 2.0)


def composite_score(cosine, mass, fragment, intensity):
    """综合得分(基础版,不含专家评分)"""
    score = (cosine * WEIGHT_COSINE + mass * WEIGHT_MASS +
             fragment * WEIGHT_FRAGMENT + intensity * WEIGHT_INTENSITY)
    return max(0.0, min(1.0, score))


def enhanced_composite_score(cosine, mass, fragment, intensity, expert=0.0, entropy=0.0):
    """增强综合得分 — 谱图相似度主导(cos+ent=0.70)。

    v3.5 重校准: 谱图整体形貌相似度(cosine+entropy)是质谱鉴定金标准, 主导综合分;
    mass/frag/inten 受实验谱vs标准品仪器差异影响系统性偏低, 降权; expert 未校准降权。
    统一单一权重公式, 不因 expert 是否>0 切换(消除"权重悬崖": 旧版 expert 略>0 时
    final 反而突降 0.07)。expert=0 时该项自然贡献 0, 无需回退分支。
    """
    score = (cosine * WEIGHT_COSINE +
             entropy * WEIGHT_ENTROPY +
             mass * WEIGHT_MASS +
             fragment * WEIGHT_FRAGMENT +
             intensity * WEIGHT_INTENSITY +
             expert * WEIGHT_EXPERT)
    # 裁剪到 [0, 1] 范围(防御性编程)
    return max(0.0, min(1.0, score))


def confidence(cos, final, expert=0.0):
    """
    置信度分级(增强版,支持专家评分)

    分级逻辑:
    - High: Cosine≥0.9 AND Final≥0.5, 或 Expert≥0.70
    - Medium: Cosine≥0.8 AND Final≥0.4, 或 Expert≥0.45
    - Low: Cosine≥0.7
    - Unknown: 其余
    """
    if EXPERT_SYSTEM_ENABLED and expert > 0:
        if (cos >= 0.9 and final >= 0.5) or expert >= 0.70:
            return "High"
        if (cos >= 0.8 and final >= 0.4) or expert >= 0.45:
            return "Medium"
        if cos >= 0.7 or expert >= 0.25:
            return "Low"
        return "Unknown"
    else:
        if cos >= 0.9 and final >= 0.5:
            return "High"
        if cos >= 0.8 and final >= 0.4:
            return "Medium"
        if cos >= 0.7:
            return "Low"
        return "Unknown"
