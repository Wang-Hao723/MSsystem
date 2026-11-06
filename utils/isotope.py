"""同位素模式识别模块 — 检测 Cl, Br, S 等特征同位素峰簇

基于已知同位素丰度比:
- Cl-35/Cl-37: 3:1 (M : M+2 ≈ 3:1)
- Br-79/Br-81: 1:1 (M : M+2 ≈ 1:1)
- S-32/S-34:   ~22:1 (M : M+2 ≈ 4.4%)
- 含多个卤素: Cl2 → 9:6:1, Br2 → 1:2:1

参考: NIST 同位素丰度表, McLafferty "Interpretation of MS" 第4章
"""
import logging
import math
import re
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

# ── 同位素间距 (Da) ─────────────────────────────────────────────────
# M+1 间距 ≈ 1.003 (C-13 vs C-12)
# M+2 间距 ≈ 1.997 (O-18 vs O-16) 或 1.998 (Cl-37 vs Cl-35)
ISOTOPE_SPACING_M1 = 1.003  # M+1 峰间距
ISOTOPE_SPACING_M2 = 1.997  # M+2 峰间距

# ── 同位素模式定义 ────────────────────────────────────────────────────
# 格式: (名称, 相对质量偏移列表, 理论丰度比列表, 检测用 m/z 容差)
ISOTOPE_PATTERNS = {
    # 单卤素
    "Cl": {
        "offsets": [0.0, 1.997],           # M, M+2
        "ratios":  [3.0, 1.0],             # 3:1
        "tol":     0.05,                   # Da
        "desc":    "1×Cl (M:M+2 ≈ 3:1)",
        "min_ratio": 2.0,                  # M/M+2 比值下限
        "max_ratio": 5.0,                  # M/M+2 比值上限
    },
    "Cl2": {
        "offsets": [0.0, 1.997, 3.994],
        "ratios":  [9.0, 6.0, 1.0],        # 9:6:1
        "tol":     0.05,
        "desc":    "2×Cl (M:M+2:M+4 ≈ 9:6:1)",
        "min_ratio": 1.0,
        "max_ratio": 2.0,
    },
    "Br": {
        "offsets": [0.0, 1.997],
        "ratios":  [1.0, 1.0],              # 1:1
        "tol":     0.05,
        "desc":    "1×Br (M:M+2 ≈ 1:1)",
        "min_ratio": 0.5,
        "max_ratio": 2.0,
    },
    "Br2": {
        "offsets": [0.0, 1.997, 3.994],
        "ratios":  [1.0, 2.0, 1.0],         # 1:2:1
        "tol":     0.05,
        "desc":    "2×Br (M:M+2:M+4 ≈ 1:2:1)",
        "min_ratio": 0.3,
        "max_ratio": 0.7,
    },
    "ClBr": {
        "offsets": [0.0, 1.997, 3.994],
        "ratios":  [3.0, 4.0, 1.0],         # 3:4:1
        "tol":     0.05,
        "desc":    "Cl+Br (M:M+2:M+4 ≈ 3:4:1)",
        "min_ratio": 0.5,
        "max_ratio": 1.5,
    },
    "S": {
        "offsets": [0.0, 1.997],
        "ratios":  [22.0, 1.0],             # ~4.4% M+2
        "tol":     0.05,
        "desc":    "1×S (M:M+2 ≈ 22:1, M+2 ≈ 4.4%)",
        "min_ratio": 10.0,
        "max_ratio": 50.0,
    },
}


def _find_peak_near(peaks_mz: np.ndarray, peaks_int: np.ndarray,
                    target_mz: float, tol: float) -> tuple:
    """
    在 peaks 中找到距离 target_mz ± tol 内最强的峰。

    Returns:
        (mz, intensity) or (None, 0.0)
    """
    if peaks_mz.size == 0:
        return None, 0.0

    mask = np.abs(peaks_mz - target_mz) <= tol
    if not np.any(mask):
        return None, 0.0

    matched_mz = peaks_mz[mask]
    matched_int = peaks_int[mask]
    best_idx = np.argmax(matched_int)
    return float(matched_mz[best_idx]), float(matched_int[best_idx])


def detect_isotope_pattern(peaks: list,
                           precursor_mz: float = None,
                           tol_da: float = 0.05) -> Optional[dict]:
    """
    检测单个谱图中的同位素模式。

    策略:
    1. 如果提供 precursor_mz，在 precursor 附近搜索 M, M+2, M+4
    2. 否则，对所有强度 >10% 的峰尝试 M+2 匹配
    3. 检验丰度比是否符合已知同位素模式

    Args:
        peaks: [[mz, intensity], ...]
        precursor_mz: 前体离子 m/z (可选，优先在此搜索)
        tol_da: m/z 容差

    Returns:
        dict or None:
            {
                "pattern": "Cl",
                "description": "1×Cl (M:M+2 ≈ 3:1)",
                "mz_mono": 300.0,          # 单同位素峰 m/z
                "mz_m2": 301.997,           # M+2 峰 m/z
                "ratio_m_m2": 3.12,         # M/M+2 实际比值
                "ratio_theory": 3.0,        # 理论比值
                "confidence": "high"
            }
    """
    if not peaks or len(peaks) < 2:
        return None

    mz_arr = np.array([p[0] for p in peaks], dtype=float)
    int_arr = np.array([p[1] for p in peaks], dtype=float)

    if mz_arr.size < 2:
        return None

    # 确定搜索起点
    if precursor_mz and precursor_mz > 0:
        # 在 precursor 附近搜索
        search_targets = [precursor_mz]
    else:
        # 对所有强度 > 10% 的峰搜索
        max_int = np.max(int_arr) if int_arr.size > 0 else 0
        if max_int <= 0:
            return None
        rel_int = int_arr / max_int
        search_targets = mz_arr[rel_int >= 0.10].tolist()

    if not search_targets:
        return None

    best_result = None
    best_confidence_score = 0

    for base_mz in search_targets:
        # 对每种模式检查
        for pattern_name, pattern in ISOTOPE_PATTERNS.items():
            offsets = pattern["offsets"]
            theory_ratios = pattern["ratios"]
            tol = pattern["tol"]

            # 找到所有同位素峰
            found_peaks = []
            for offset in offsets:
                target = base_mz + offset
                found_mz, found_int = _find_peak_near(mz_arr, int_arr, target, tol)
                if offset == 0.0:
                    # M 峰必须存在
                    if found_mz is None:
                        break
                    found_peaks.append((found_mz, found_int))
                else:
                    found_peaks.append((found_mz, found_int))

            # 检查是否所有必须的峰都找到了
            if len(found_peaks) < len(offsets):
                continue

            # M 峰缺失则跳过
            if found_peaks[0][0] is None:
                continue

            # 计算丰度比
            m_int = found_peaks[0][1]
            if m_int <= 0:
                continue

            m2_int = found_peaks[1][1] if len(found_peaks) > 1 else 0
            ratio_m_m2 = m_int / m2_int if m2_int > 0 else float("inf")

            # 只检查双峰模式 (M, M+2) 的丰度比
            if len(offsets) == 2:
                min_r = pattern["min_ratio"]
                max_r = pattern["max_ratio"]

                if not (min_r <= ratio_m_m2 <= max_r):
                    continue

                # 计算置信度得分 (越接近理论比值得分越高)
                theory_ratio = theory_ratios[0] / theory_ratios[1]
                deviation = abs(ratio_m_m2 - theory_ratio) / theory_ratio
                conf_score = max(0, 1.0 - deviation)

                # 强度加权: M 峰越强越可靠
                max_int = np.max(int_arr)
                intensity_weight = m_int / max_int if max_int > 0 else 0
                conf_score *= (0.5 + 0.5 * intensity_weight)

                if conf_score > best_confidence_score:
                    best_confidence_score = conf_score
                    confidence = "high" if conf_score > 0.7 else "medium" if conf_score > 0.4 else "low"

                    best_result = {
                        "pattern": pattern_name,
                        "description": pattern["desc"],
                        "mz_mono": round(found_peaks[0][0], 4),
                        "mz_m2": round(found_peaks[1][0], 4) if found_peaks[1][0] else None,
                        "ratio_m_m2": round(ratio_m_m2, 3),
                        "ratio_theory": round(theory_ratio, 3),
                        "confidence": confidence,
                    }

                    # M+4 峰 (如果有)
                    if len(found_peaks) > 2 and found_peaks[2][0] is not None:
                        m4_int = found_peaks[2][1]
                        best_result["mz_m4"] = round(found_peaks[2][0], 4)
                        best_result["ratio_m_m4"] = round(m_int / m4_int, 3) if m4_int > 0 else None

            elif len(offsets) == 3:
                # 三峰模式 (M, M+2, M+4)
                if found_peaks[1][0] is None or found_peaks[2][0] is None:
                    continue

                m2_int = found_peaks[1][1]
                m4_int = found_peaks[2][1]

                ratio_m_m2 = m_int / m2_int if m2_int > 0 else float("inf")
                ratio_m_m4 = m_int / m4_int if m4_int > 0 else float("inf")

                theory_m_m2 = theory_ratios[0] / theory_ratios[1]
                theory_m_m4 = theory_ratios[0] / theory_ratios[2]

                dev1 = abs(ratio_m_m2 - theory_m_m2) / theory_m_m2
                dev2 = abs(ratio_m_m4 - theory_m_m4) / theory_m_m4
                conf_score = max(0, 1.0 - (dev1 + dev2) / 2)

                if conf_score > best_confidence_score:
                    best_confidence_score = conf_score
                    confidence = "high" if conf_score > 0.7 else "medium" if conf_score > 0.4 else "low"

                    best_result = {
                        "pattern": pattern_name,
                        "description": pattern["desc"],
                        "mz_mono": round(found_peaks[0][0], 4),
                        "mz_m2": round(found_peaks[1][0], 4),
                        "mz_m4": round(found_peaks[2][0], 4),
                        "ratio_m_m2": round(ratio_m_m2, 3),
                        "ratio_m_m4": round(ratio_m_m4, 3),
                        "ratio_theory_m2": round(theory_m_m2, 3),
                        "ratio_theory_m4": round(theory_m_m4, 3),
                        "confidence": confidence,
                    }

    return best_result


# Natural-abundance probabilities truncated to nominal M, M+1 and M+2.
# These are used only to compare an observed MS1 precursor envelope with a
# candidate formula.  The score is supporting evidence, never a hard formula
# exclusion, because co-elution and isolation interference can distort ratios.
_NATURAL_ISOTOPES = {
    "C": (0.9893, 0.0107, 0.0),
    "H": (0.999885, 0.000115, 0.0),
    "N": (0.99636, 0.00364, 0.0),
    "O": (0.99757, 0.00038, 0.00205),
    "S": (0.9499, 0.0075, 0.0425),
    "Cl": (0.7578, 0.0, 0.2422),
    "Br": (0.5069, 0.0, 0.4931),
    "Si": (0.92223, 0.04685, 0.03092),
    "F": (1.0, 0.0, 0.0),
    "P": (1.0, 0.0, 0.0),
    "I": (1.0, 0.0, 0.0),
    "Na": (1.0, 0.0, 0.0),
    "K": (0.932581, 0.0, 0.067302),
}
_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def detect_precursor_isotope_envelope(
    ms1_peaks: list,
    precursor_mz: float,
    charge: int | None = 1,
    tol_da: float | None = None,
    tol_ppm: float = 5.0,
    max_isotopes: int = 3,
) -> Optional[dict]:
    """Extract an isotope envelope around a selected precursor from an MS1 scan.

    This function must receive MS1 peaks.  Applying it to product-ion (MS2)
    peaks is scientifically invalid because those peaks are fragments, not the
    precursor isotope envelope.
    """
    if not ms1_peaks or not precursor_mz or precursor_mz <= 0:
        return None
    try:
        z = abs(int(charge or 1))
    except (TypeError, ValueError):
        z = 1
    if z < 1 or z > 8:
        z = 1

    mz_arr = np.asarray([p[0] for p in ms1_peaks], dtype=float)
    int_arr = np.asarray([p[1] for p in ms1_peaks], dtype=float)
    if mz_arr.size == 0 or int_arr.size != mz_arr.size:
        return None
    if tol_da is None:
        # QE-Orbitrap expected-compound extraction is configured at 5 ppm.
        # A fixed 0.01 Da window would be ~43 ppm at m/z 230 and can capture
        # unrelated co-eluting ions as isotope peaks.  The small absolute
        # floor only accommodates centroid values rounded to four decimals.
        tolerance = max(
            0.0005,
            float(precursor_mz) * max(0.0, float(tol_ppm)) * 1e-6,
        )
        tolerance_source = "ppm_with_centroid_rounding_floor"
    else:
        tolerance = float(tol_da)
        if not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("tol_da must be a positive finite number")
        tolerance_source = "explicit_da"
    mono_mz, mono_int = _find_peak_near(mz_arr, int_arr, float(precursor_mz), tolerance)
    if mono_mz is None or mono_int <= 0:
        return None

    spacing = 1.00335483507 / z
    peaks = [{"order": 0, "mz": round(mono_mz, 6), "intensity": round(mono_int, 3),
              "relative_to_m": 1.0}]
    for order in range(1, max(1, int(max_isotopes)) + 1):
        found_mz, found_int = _find_peak_near(
            mz_arr, int_arr, mono_mz + spacing * order, tolerance
        )
        peaks.append({
            "order": order,
            "mz": round(found_mz, 6) if found_mz is not None else None,
            "intensity": round(found_int, 3) if found_mz is not None else 0.0,
            "relative_to_m": round(found_int / mono_int, 5) if found_mz is not None else None,
        })

    observed = sum(1 for peak in peaks[1:] if peak["mz"] is not None)
    return {
        "source": "adjacent_ms1_scan",
        "status": "observed" if observed else "monoisotopic_peak_only",
        "charge": z,
        "spacing_da": round(spacing, 6),
        "mz_tolerance_da": round(tolerance, 6),
        "mz_tolerance_ppm": (
            round(float(tol_ppm), 6) if tol_da is None else None
        ),
        "mz_tolerance_source": tolerance_source,
        "n_isotopes_observed": observed,
        "peaks": peaks,
    }


def theoretical_formula_isotope_ratios(formula: str) -> Optional[dict]:
    """Return approximate natural-abundance M+1/M and M+2/M ratios."""
    text = re.sub(r"\s+", "", str(formula or ""))
    if not text:
        return None
    tokens = _FORMULA_TOKEN_RE.findall(text)
    if not tokens or "".join(f"{e}{n}" for e, n in tokens) != text:
        return None

    distribution = np.array([1.0, 0.0, 0.0], dtype=float)
    counts = {}
    for element, raw_count in tokens:
        if element not in _NATURAL_ISOTOPES:
            return None
        count = int(raw_count or 1)
        if count < 0 or count > 1000:
            return None
        counts[element] = counts.get(element, 0) + count
        atom = np.asarray(_NATURAL_ISOTOPES[element], dtype=float)
        for _ in range(count):
            distribution = np.convolve(distribution, atom)[:3]
    if distribution[0] <= 0:
        return None
    return {
        "m1_over_m": round(float(distribution[1] / distribution[0]), 6),
        "m2_over_m": round(float(distribution[2] / distribution[0]), 6),
        "element_counts": counts,
    }


def evaluate_formula_isotope_consistency(formula: str, envelope: dict | None) -> Optional[dict]:
    """Compare an MS1 envelope with a candidate formula conservatively."""
    if not envelope or envelope.get("source") != "adjacent_ms1_scan":
        return None
    theory = theoretical_formula_isotope_ratios(formula)
    peaks = envelope.get("peaks") or []
    if theory is None or len(peaks) < 2:
        return None

    comparisons = []
    for order, key in ((1, "m1_over_m"), (2, "m2_over_m")):
        if len(peaks) <= order or peaks[order].get("relative_to_m") is None:
            continue
        observed = max(0.0, float(peaks[order]["relative_to_m"]))
        expected = max(0.0, float(theory[key]))
        # A small floor prevents trace-isotope ratios from dominating the score.
        log_error = abs(math.log((observed + 0.01) / (expected + 0.01)))
        comparisons.append((order, observed, expected, log_error))
    if not comparisons:
        return None

    mean_error = sum(item[3] for item in comparisons) / len(comparisons)
    score = float(math.exp(-mean_error / 0.75))
    status = "supporting" if score >= 0.75 else "consistent" if score >= 0.45 else "conflicting"
    return {
        "status": status,
        "score": round(score, 4),
        "formula": str(formula),
        "observed": {
            f"m{order}_over_m": round(observed, 5)
            for order, observed, _, _ in comparisons
        },
        "theoretical": {
            f"m{order}_over_m": round(expected, 5)
            for order, _, expected, _ in comparisons
        },
        "element_counts": theory["element_counts"],
        "interpretation": (
            "MS1 precursor isotope envelope supports the candidate formula"
            if status == "supporting" else
            "MS1 isotope ratios are broadly compatible with the candidate formula"
            if status == "consistent" else
            "MS1 isotope ratios do not closely match the candidate formula; co-elution may interfere"
        ),
    }


def detect_isotopes_batch(spectra_list: list[dict],
                          tol_da: float = 0.05) -> list[dict]:
    """
    批量检测同位素模式。

    Args:
        spectra_list: [{"precursor_mz": float, "peaks": [[mz, int], ...]}, ...]
        tol_da: m/z 容差

    Returns:
        list[dict]: 每个谱图的同位素检测结果
            {
                "index": 0,
                "precursor_mz": 300.0,
                "isotope": {...} or None,
            }
    """
    results = []
    for i, spec in enumerate(spectra_list):
        pmz = spec.get("precursor_mz", 0)
        peaks = spec.get("peaks", [])

        isotope = detect_isotope_pattern(peaks, precursor_mz=pmz, tol_da=tol_da)

        results.append({
            "index": i,
            "precursor_mz": pmz,
            "isotope": isotope,
        })

    # 统计
    detected = sum(1 for r in results if r["isotope"] is not None)
    if detected > 0:
        patterns = {}
        for r in results:
            if r["isotope"]:
                p = r["isotope"]["pattern"]
                patterns[p] = patterns.get(p, 0) + 1
        logger.info("Isotope detection: %d/%d spectra with patterns: %s",
                     detected, len(results), patterns)

    return results
