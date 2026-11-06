#!/usr/bin/env python3
"""Blank-feature assessment for LC-MS/MS analysis.

The production web pipeline does not subtract a global list of blank fragment
ions from every sample spectrum. It instead:

1. aligns each sample MS2 precursor to same-polarity blank MS1 XIC evidence
   using both the configured Da and ppm gates plus a symmetric RT window;
2. compares the same XIC apex statistic in sample and blank;
3. combines replicate blanks within analytical role by detection prevalence
   and median abundance, never pooling carryover and procedural blanks;
4. defaults to searching every feature and annotating sample-enriched,
   blank-associated, blank-dominant and carryover evidence separately;
5. leaves missing or incompatible evidence unevaluated rather than deleting
   the sample feature.

MGF lacks continuous MS1 scans, so it can only use an explicitly labelled,
conservative precursor-intensity fallback. ``subtract_blank`` remains below
for backwards compatibility with direct callers; ``web.server`` uses
``assess_blank_role_groups`` and never performs fragment-global subtraction.
"""
from __future__ import annotations

import logging
import math
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Sequence

import numpy as np

logger = logging.getLogger(__name__)

BLANK_SUBTRACTION_AVAILABLE = True


# ─────────────────────────────────────────────────────────────────
# 配置常量
# ─────────────────────────────────────────────────────────────────
DEFAULT_MZ_TOLERANCE_DA = 0.01  # m/z 匹配容差（Da）
DEFAULT_MZ_TOLERANCE_PPM = 10  # m/z 匹配容差（ppm）
DEFAULT_PRECURSOR_TOLERANCE_DA = 0.01
DEFAULT_PRECURSOR_TOLERANCE_PPM = 10
DEFAULT_RT_TOLERANCE_MINUTES = 0.20
DEFAULT_FEATURE_RT_TOLERANCE_MINUTES = 0.30
DEFAULT_SAMPLE_TO_BLANK_RATIO = 3.0
DEFAULT_BLANK_MIN_DETECTION_FRACTION = 0.50
DEFAULT_BLANK_FILTER_POLICY = "soft_evidence"
DEFAULT_BLANK_CLEAN_MAX_RATIO = 1.0
DEFAULT_BLANK_CLEAN_MIN_REPLICATES = 2
DEFAULT_INTENSITY_THRESHOLD_FACTOR = 0.1  # 空白强度阈值因子（10%）
DEFAULT_SUBTRACTION_MODE: Literal["absolute", "relative", "filter_only"] = "absolute"

BLANK_FILTER_POLICIES = frozenset({"soft_evidence", "hard_legacy"})
BLANK_ROLES = frozenset({
    "unspecified",
    "background_blank",
    "solvent_blank",
    "method_blank",
    "extraction_blank",
    "carryover_blank",
})
BACKGROUND_BLANK_ROLES = frozenset({
    "background_blank",
    "solvent_blank",
    "method_blank",
    "extraction_blank",
})


@dataclass(frozen=True)
class BlankSpectrum:
    """One blank MS/MS spectrum with the acquisition context kept intact.

    Fragment ions are meaningful blank evidence only for a sample spectrum at
    the same precursor and, when recorded, the same chromatographic region.
    Flattening every blank fragment into a global m/z list can otherwise erase
    genuine product ions from unrelated compounds.
    """

    precursor_mz: float
    retention_time_minutes: float | None
    mz_array: np.ndarray
    intensity_array: np.ndarray
    polarity: str | None = None
    precursor_intensity: float | None = None


def _positive_float(value) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if np.isfinite(parsed) and parsed > 0 else None


def _rt_to_minutes(value, unit: str | None) -> float | None:
    """Normalize a retention time to minutes without guessing invalid data."""
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(parsed) or parsed < 0:
        return None
    normalized = str(unit or "minutes").strip().lower()
    if normalized in {"s", "sec", "second", "seconds"}:
        return parsed / 60.0
    if normalized in {"m", "min", "minute", "minutes"}:
        return parsed
    # Unknown units are unsafe for chromatographic alignment.
    return None


def _mass_within_tolerance(
    observed: float,
    reference: float,
    tolerance_da: float,
    tolerance_ppm: float,
    *,
    require_both: bool,
) -> bool:
    diff = abs(float(observed) - float(reference))
    ppm = diff / abs(float(observed)) * 1e6 if observed else float("inf")
    da_ok = tolerance_da <= 0 or diff <= tolerance_da
    ppm_ok = tolerance_ppm <= 0 or ppm <= tolerance_ppm
    return da_ok and ppm_ok if require_both else da_ok or ppm_ok


def _make_blank_spectrum(
    peaks: Iterable[Sequence[float]],
    precursor_mz,
    retention_time=None,
    retention_time_unit: str | None = "minutes",
    polarity: str | None = None,
    precursor_intensity=None,
) -> BlankSpectrum | None:
    precursor = _positive_float(precursor_mz)
    if precursor is None:
        return None
    valid: list[tuple[float, float]] = []
    for peak in peaks or []:
        if len(peak) < 2:
            continue
        mz = _positive_float(peak[0])
        intensity = _positive_float(peak[1])
        if mz is not None and intensity is not None:
            valid.append((mz, intensity))
    if not valid:
        return None
    valid.sort(key=lambda item: item[0])
    return BlankSpectrum(
        precursor_mz=precursor,
        retention_time_minutes=_rt_to_minutes(retention_time, retention_time_unit),
        mz_array=np.asarray([item[0] for item in valid], dtype=float),
        intensity_array=np.asarray([item[1] for item in valid], dtype=float),
        polarity=(str(polarity).strip().upper()
                  if str(polarity or "").strip().upper() in {"POS", "NEG"}
                  else None),
        precursor_intensity=_positive_float(precursor_intensity),
    )


class BlankProfile:
    """
    空白溶剂峰轮廓

    存储从空白样品中提取的峰信息，用于后续扣除。
    """

    def __init__(
        self,
        mz_array: np.ndarray,
        intensity_max: np.ndarray,
        intensity_mean: np.ndarray,
        intensity_std: np.ndarray,
        count: np.ndarray,
        mz_tolerance_da: float = DEFAULT_MZ_TOLERANCE_DA,
        mz_tolerance_ppm: float = DEFAULT_MZ_TOLERANCE_PPM,
        context_spectra: Sequence[BlankSpectrum] | None = None,
        precursor_tolerance_da: float = DEFAULT_PRECURSOR_TOLERANCE_DA,
        precursor_tolerance_ppm: float = DEFAULT_PRECURSOR_TOLERANCE_PPM,
        rt_tolerance_minutes: float = DEFAULT_RT_TOLERANCE_MINUTES,
        query_rt_unit: str = "minutes",
        source_path: str | Path | None = None,
        ms1_scan_count: int = 0,
        available_polarities: Sequence[str] | None = None,
    ):
        """
        Args:
            mz_array: 排序后的 m/z 数组
            intensity_max: 每个 m/z 处的最大强度
            intensity_mean: 每个 m/z 处的平均强度
            intensity_std: 每个 m/z 处的强度标准差
            count: 每个 m/z 处的观测次数（多少个谱图包含此峰）
            mz_tolerance_da: m/z 匹配绝对容差（Da）
            mz_tolerance_ppm: m/z 匹配相对容差（ppm）
        """
        self.mz_array = mz_array
        self.intensity_max = intensity_max
        self.intensity_mean = intensity_mean
        self.intensity_std = intensity_std
        self.count = count
        self.mz_tolerance_da = mz_tolerance_da
        self.mz_tolerance_ppm = mz_tolerance_ppm
        self.context_spectra = tuple(sorted(
            context_spectra or (), key=lambda spectrum: spectrum.precursor_mz,
        ))
        self._context_precursors = np.asarray(
            [spectrum.precursor_mz for spectrum in self.context_spectra],
            dtype=float,
        )
        self.precursor_tolerance_da = float(precursor_tolerance_da)
        self.precursor_tolerance_ppm = float(precursor_tolerance_ppm)
        self.rt_tolerance_minutes = float(rt_tolerance_minutes)
        self.query_rt_unit = str(query_rt_unit or "minutes")
        self.source_path = Path(source_path) if source_path else None
        self.ms1_scan_count = int(ms1_scan_count or 0)
        inferred_polarities = {
            spectrum.polarity for spectrum in self.context_spectra
            if spectrum.polarity in {"POS", "NEG"}
        }
        supplied_polarities = {
            str(value).strip().upper() for value in (available_polarities or ())
            if str(value).strip().upper() in {"POS", "NEG"}
        }
        self.available_polarities = frozenset(
            supplied_polarities or inferred_polarities
        )

        # 预计算阈值数组
        self.threshold = intensity_max * DEFAULT_INTENSITY_THRESHOLD_FACTOR

    def _matching_contexts(
        self,
        precursor_mz,
        retention_time=None,
        retention_time_unit: str | None = None,
        polarity: str | None = None,
    ) -> list[BlankSpectrum]:
        """Return precursor/RT-aligned blank spectra, or no evidence.

        A precursor is mandatory.  If aligned blank spectra carry RT, a query
        RT with a known unit is mandatory too.  Context-free blank profiles are
        intentionally inert rather than globally deleting common fragments.
        """
        precursor = _positive_float(precursor_mz)
        if precursor is None or not self.context_spectra:
            return []
        windows = []
        if self.precursor_tolerance_da > 0:
            windows.append(self.precursor_tolerance_da)
        if self.precursor_tolerance_ppm > 0:
            windows.append(precursor * self.precursor_tolerance_ppm * 1e-6)
        window = min(windows) if windows else 0.0
        left = int(np.searchsorted(
            self._context_precursors, precursor - window, side="left",
        ))
        right = int(np.searchsorted(
            self._context_precursors, precursor + window, side="right",
        ))
        precursor_matches = [
            spectrum for spectrum in self.context_spectra[left:right]
            if _mass_within_tolerance(
                precursor,
                spectrum.precursor_mz,
                self.precursor_tolerance_da,
                self.precursor_tolerance_ppm,
                require_both=True,
            )
        ]
        normalized_polarity = str(polarity or "").strip().upper()
        if normalized_polarity in {"POS", "NEG"}:
            precursor_matches = [
                spectrum for spectrum in precursor_matches
                if spectrum.polarity == normalized_polarity
            ]
        if not precursor_matches:
            return []

        timed = [
            spectrum for spectrum in precursor_matches
            if spectrum.retention_time_minutes is not None
        ]
        if not timed:
            return precursor_matches
        query_minutes = _rt_to_minutes(
            retention_time,
            retention_time_unit or self.query_rt_unit,
        )
        if query_minutes is None:
            return []
        return [
            spectrum for spectrum in timed
            if abs(query_minutes - float(spectrum.retention_time_minutes))
            <= self.rt_tolerance_minutes
        ]

    def has_matching_context(
        self,
        precursor_mz,
        retention_time=None,
        retention_time_unit: str | None = None,
        polarity: str | None = None,
    ) -> bool:
        """Whether this profile can safely be applied to a sample spectrum."""
        return bool(self._matching_contexts(
            precursor_mz,
            retention_time,
            retention_time_unit,
            polarity,
        ))

    @staticmethod
    def _find_peak_in_spectrum(
        spectrum: BlankSpectrum,
        mz: float,
        mz_tolerance_da: float,
        mz_tolerance_ppm: float,
    ) -> tuple[int, float] | None:
        idx = int(np.searchsorted(spectrum.mz_array, mz))
        start = max(0, idx - 2)
        end = min(len(spectrum.mz_array), idx + 3)
        best_match = None
        best_diff = float("inf")
        for i in range(start, end):
            blank_mz = float(spectrum.mz_array[i])
            diff = abs(float(mz) - blank_mz)
            if _mass_within_tolerance(
                float(mz), blank_mz, mz_tolerance_da, mz_tolerance_ppm,
                require_both=False,
            ) and diff < best_diff:
                best_match = (i, blank_mz)
                best_diff = diff
        return best_match

    def find_match(
        self,
        mz: float,
        *,
        precursor_mz=None,
        retention_time=None,
        retention_time_unit: str | None = None,
    ) -> tuple[int, float] | None:
        """
        在空白轮廓中查找匹配的 m/z

        Args:
            mz: 查询的 m/z 值

        Returns:
            (索引, 匹配的空白 m/z) 或 None
        """
        contexts = self._matching_contexts(
            precursor_mz, retention_time, retention_time_unit,
        )
        best_match = None
        best_diff = float("inf")
        for spectrum in contexts:
            match = self._find_peak_in_spectrum(
                spectrum, mz, self.mz_tolerance_da, self.mz_tolerance_ppm,
            )
            if match is not None and abs(float(mz) - match[1]) < best_diff:
                best_match = match
                best_diff = abs(float(mz) - match[1])
        return best_match

    def get_blank_intensity(
        self,
        mz: float,
        mode: Literal["max", "mean", "std"] = "max",
        *,
        precursor_mz=None,
        retention_time=None,
        retention_time_unit: str | None = None,
    ) -> tuple[float, int] | None:
        """
        获取指定 m/z 处的空白强度

        Args:
            mz: 查询的 m/z 值
            mode: 使用哪种强度（max/mean/std）

        Returns:
            (空白强度, 观测次数) 或 None（无匹配）
        """
        contexts = self._matching_contexts(
            precursor_mz, retention_time, retention_time_unit,
        )
        return self._blank_intensity_from_contexts(mz, contexts, mode)

    def _blank_intensity_from_contexts(
        self,
        mz: float,
        contexts: Sequence[BlankSpectrum],
        mode: Literal["max", "mean", "std"] = "max",
    ) -> tuple[float, int] | None:
        """Aggregate one fragment only within pre-aligned blank spectra."""
        intensities: list[float] = []
        for spectrum in contexts:
            match = self._find_peak_in_spectrum(
                spectrum, mz, self.mz_tolerance_da, self.mz_tolerance_ppm,
            )
            if match is not None:
                intensities.append(float(spectrum.intensity_array[match[0]]))
        if not intensities:
            return None
        if mode == "max":
            intensity = max(intensities)
        elif mode == "mean":
            intensity = float(np.mean(intensities))
        elif mode == "std":
            intensity = float(np.std(intensities))
        else:
            intensity = max(intensities)
        return (float(intensity), len(intensities))


def _profile_from_spectra(
    spectra: Sequence[dict],
    mz_tolerance_da: float = DEFAULT_MZ_TOLERANCE_DA,
    mz_tolerance_ppm: float = DEFAULT_MZ_TOLERANCE_PPM,
    *,
    retention_time_unit: str = "minutes",
) -> BlankProfile:
    """Build display aggregates and safe spectrum-level blank contexts."""
    peak_dict: dict[float, dict[str, list[float]]] = {}
    context_spectra: list[BlankSpectrum] = []

    for spec in spectra:
        peaks = spec.get("peaks", [])
        context = _make_blank_spectrum(
            peaks,
            spec.get("precursor_mz"),
            spec.get("retention_time"),
            spec.get("retention_time_unit") or retention_time_unit,
            spec.get("polarity") or spec.get("ion_mode"),
            spec.get("precursor_ms1_intensity") or spec.get("precursor_intensity"),
        )
        if context is not None:
            context_spectra.append(context)
        for mz, intensity in peaks:
            if mz <= 0 or intensity <= 0:
                continue

            # 使用字符串化的 m/z 作为 key（考虑容差）
            # 简单起见，使用 0.01 Da 的 bin
            bin_key = round(mz / 0.01) * 0.01

            if bin_key not in peak_dict:
                peak_dict[bin_key] = {
                    "intensities": [],
                    "mz_values": [],
                }

            peak_dict[bin_key]["intensities"].append(intensity)
            peak_dict[bin_key]["mz_values"].append(mz)

    mz_list: list[float] = []
    intensity_max_list: list[float] = []
    intensity_mean_list: list[float] = []
    intensity_std_list: list[float] = []
    count_list: list[int] = []

    for bin_key, data in sorted(peak_dict.items()):
        intensities = np.array(data["intensities"])
        mz_values = np.array(data["mz_values"])

        mz_avg = np.mean(mz_values)

        mz_list.append(mz_avg)
        intensity_max_list.append(np.max(intensities))
        intensity_mean_list.append(np.mean(intensities))
        intensity_std_list.append(np.std(intensities))
        count_list.append(len(intensities))

    return BlankProfile(
        mz_array=np.asarray(mz_list, dtype=float),
        intensity_max=np.asarray(intensity_max_list, dtype=float),
        intensity_mean=np.asarray(intensity_mean_list, dtype=float),
        intensity_std=np.asarray(intensity_std_list, dtype=float),
        count=np.asarray(count_list, dtype=int),
        mz_tolerance_da=mz_tolerance_da,
        mz_tolerance_ppm=mz_tolerance_ppm,
        context_spectra=context_spectra,
        query_rt_unit=retention_time_unit,
    )


def load_blank_mgf(
    file_path: str | Path,
    mz_tolerance_da: float = DEFAULT_MZ_TOLERANCE_DA,
    mz_tolerance_ppm: float = DEFAULT_MZ_TOLERANCE_PPM,
) -> BlankProfile:
    """Load an MGF blank while preserving precursor and RT per spectrum."""
    file_path = Path(file_path)
    logger.info("加载空白溶剂文件: %s", file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"空白文件不存在: {file_path}")
    spectra = _parse_mgf(file_path)
    logger.info("  空白谱图数: %d", len(spectra))
    if not spectra:
        raise ValueError(f"空白文件中没有找到有效谱图: {file_path}")
    profile = _profile_from_spectra(
        spectra,
        mz_tolerance_da,
        mz_tolerance_ppm,
        retention_time_unit="seconds",
    )
    profile.source_path = file_path
    logger.info(
        "  空白峰数量: %d（有前体上下文的谱图: %d）",
        len(profile.mz_array),
        len(profile.context_spectra),
    )
    return profile


def _parse_mgf(file_path: Path) -> list[dict]:
    """解析 MGF 文件，返回谱图列表"""
    spectra = []
    raw_text = file_path.read_text(encoding="utf-8", errors="replace")

    for block in raw_text.split("BEGIN IONS"):
        if not block.strip():
            continue

        peaks = []
        precursor_mz = None
        retention_time = None
        retention_time_unit = "seconds"
        polarity = None
        precursor_intensity = None

        for line in block.split("\n"):
            line = line.strip()
            if not line:
                continue
            upper = line.upper()
            if upper.startswith("PEPMASS="):
                try:
                    precursor_mz = float(line.split("=", 1)[1].split()[0])
                except (ValueError, IndexError):
                    pass
            elif upper.startswith("RTINSECONDS="):
                try:
                    retention_time = float(line.split("=", 1)[1].split()[0])
                    retention_time_unit = "seconds"
                except (ValueError, IndexError):
                    pass
            elif upper.startswith("RTINMINUTES="):
                try:
                    retention_time = float(line.split("=", 1)[1].split()[0])
                    retention_time_unit = "minutes"
                except (ValueError, IndexError):
                    pass
            elif upper.startswith(("IONMODE=", "ION_MODE=")):
                value = line.split("=", 1)[1].strip().upper()
                if value in {"POS", "POSITIVE", "+"}:
                    polarity = "POS"
                elif value in {"NEG", "NEGATIVE", "-"}:
                    polarity = "NEG"
            elif upper.startswith("CHARGE=") and polarity is None:
                value = line.split("=", 1)[1]
                polarity = "NEG" if "-" in value else "POS" if "+" in value else None
            elif upper.startswith(("PRECURSORINTENSITY=", "PRECURSOR_INTENSITY=")):
                try:
                    precursor_intensity = float(line.split("=", 1)[1].split()[0])
                except (ValueError, IndexError):
                    pass
            elif "END IONS" in line:
                break
            elif "\t" in line or (" " in line and "=" not in line):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        if mz > 0 and intensity > 0:
                            peaks.append([mz, intensity])
                    except ValueError:
                        pass

        if peaks:
            spectra.append({
                "precursor_mz": precursor_mz,
                "retention_time": retention_time,
                "retention_time_unit": retention_time_unit,
                "polarity": polarity,
                "precursor_intensity": precursor_intensity,
                "peaks": peaks,
            })

    return spectra


def subtract_blank(
    peaks: list[list[float]],
    blank_profile: BlankProfile,
    mode: Literal["absolute", "relative", "filter_only"] = DEFAULT_SUBTRACTION_MODE,
    intensity_threshold_factor: float = DEFAULT_INTENSITY_THRESHOLD_FACTOR,
    blank_intensity_mode: Literal["max", "mean"] = "max",
    *,
    precursor_mz=None,
    retention_time=None,
    retention_time_unit: str | None = None,
    polarity: str | None = None,
) -> list[list[float]]:
    """
    从样品峰中扣除空白贡献

    Args:
        peaks: 样品峰列表 [[mz, intensity], ...]
        blank_profile: 空白轮廓对象
        mode: 扣除模式
            - "absolute": 绝对扣除（样品强度 - 空白强度）
            - "relative": 相对扣除（扣除空白强度的比例）
            - "filter_only": 仅过滤，不进行扣除
        intensity_threshold_factor: 强度阈值因子
            如果样品峰强度 < 空白最大强度 * 此因子，则过滤掉
        blank_intensity_mode: 使用空白的哪种强度（max 或 mean）
        precursor_mz: 当前样品 MS/MS 的前体 m/z；缺失时为安全起见不扣除
        retention_time: 当前样品谱图保留时间；空白有 RT 时必须提供
        retention_time_unit: seconds 或 minutes；未知单位时不做 RT 扣除

    Returns:
        扣除后的峰列表
    """
    if not peaks:
        return []
    contexts = blank_profile._matching_contexts(
        precursor_mz,
        retention_time,
        retention_time_unit,
        polarity,
    )
    if not contexts:
        # 旧实现把整个空白文件的所有碎片 m/z 压成一张全局表，导致不同
        # 前体/不同 RT 的真实产物离子被误删。上下文不完整时宁可不扣除。
        return [list(peak) for peak in peaks]

    result = []

    for mz, intensity in peaks:
        blank_info = blank_profile._blank_intensity_from_contexts(
            mz,
            contexts,
            mode=blank_intensity_mode,
        )

        if blank_info is None:
            # 空白中没有此峰，保留
            result.append([mz, intensity])
            continue

        blank_intensity, blank_count = blank_info

        # 计算阈值
        threshold = blank_intensity * intensity_threshold_factor

        if intensity <= threshold:
            # 低于阈值，过滤掉
            continue

        # 高于阈值，进行扣除
        if mode == "absolute":
            new_intensity = intensity - blank_intensity
        elif mode == "relative":
            # 扣除空白强度占总强度的比例
            if intensity > blank_intensity:
                ratio = blank_intensity / intensity
                new_intensity = intensity * (1 - ratio * intensity_threshold_factor)
            else:
                new_intensity = intensity * 0.1  # 保留 10%
        else:  # filter_only
            new_intensity = intensity

        if new_intensity > 0:
            result.append([mz, new_intensity])

    return result


def create_blank_profile_from_list(
    blank_peaks: list[tuple[float, float]],
    mz_tolerance_da: float = DEFAULT_MZ_TOLERANCE_DA,
    mz_tolerance_ppm: float = DEFAULT_MZ_TOLERANCE_PPM,
    *,
    precursor_mz=None,
    retention_time=None,
    retention_time_unit: str = "minutes",
) -> BlankProfile:
    """
    从给定的峰列表创建空白轮廓

    用于手动指定空白峰的场景。

    Args:
        blank_peaks: 空白峰列表 [(mz, intensity), ...]
        mz_tolerance_da: m/z 匹配绝对容差
        mz_tolerance_ppm: m/z 匹配相对容差
        precursor_mz: 可选的空白谱前体；缺失时轮廓仅供展示，不执行扣除
        retention_time: 可选的空白谱 RT
        retention_time_unit: seconds 或 minutes

    Returns:
        BlankProfile 对象
    """
    if not blank_peaks:
        raise ValueError("空白峰列表不能为空")
    return _profile_from_spectra(
        [{
            "precursor_mz": precursor_mz,
            "retention_time": retention_time,
            "retention_time_unit": retention_time_unit,
            "polarity": None,
            "peaks": blank_peaks,
        }],
        mz_tolerance_da,
        mz_tolerance_ppm,
        retention_time_unit=retention_time_unit,
    )


def load_blank_file(
    file_path: str | Path,
    mz_tolerance_da: float = DEFAULT_MZ_TOLERANCE_DA,
    mz_tolerance_ppm: float = DEFAULT_MZ_TOLERANCE_PPM,
) -> BlankProfile:
    """格式无关的空白轮廓加载 (.mgf 或 .mzML/.mzXML)。

    两种格式都保留每张 MS2 的前体和 RT 上下文。全局聚合数组仅用于上传页
    的峰数摘要，不再用于跨谱删峰。

    新增(v3 Phase 1): 原仅支持 MGF, 现兼容用户上传的 mzML 空白(如 BK_POS1.mzML)。
    纯增量, 不改 load_blank_mgf。
    """
    file_path = Path(file_path)
    ext = file_path.suffix.lower()
    if ext == ".mgf":
        return load_blank_mgf(file_path, mz_tolerance_da, mz_tolerance_ppm)
    if ext in (".mzml", ".mzxml", ".xml"):
        from utils.mzml_io import parse_mzml_file
        spectra = parse_mzml_file(str(file_path))
        ms1_scan_count = _count_mzml_ms1_scans(file_path)
        if (not spectra or not any(spec.get("peaks") for spec in spectra)) and not ms1_scan_count:
            raise ValueError(f"空白 mzML 中无有效 MS1/MS2 峰: {file_path}")
        profile = _profile_from_spectra(
            spectra or [],
            mz_tolerance_da,
            mz_tolerance_ppm,
            retention_time_unit="minutes",
        )
        profile.source_path = file_path
        profile.ms1_scan_count = ms1_scan_count
        try:
            from utils.mzml_io import _file_polarity
            file_polarity = _file_polarity(file_path)
            if file_polarity in {"POS", "NEG"}:
                profile.available_polarities = frozenset({file_polarity})
        except Exception:
            pass
        logger.info(
            "空白 mzML %s: %d MS2 → %d 聚合峰（%d 张有前体上下文，%d 张 MS1）",
            file_path.name,
            len(spectra),
            len(profile.mz_array),
            len(profile.context_spectra),
            profile.ms1_scan_count,
        )
        return profile
    raise ValueError(f"不支持的空白文件格式: {ext} (允许 .mgf/.mzML/.mzXML)")


def _count_mzml_ms1_scans(file_path: str | Path) -> int:
    try:
        with _open_ms_reader(Path(file_path)) as reader:
            return sum(1 for spectrum in reader if spectrum.get("ms level") == 1)
    except Exception as exc:
        logger.warning("空白 MS1 计数失败 %s: %s", file_path, exc)
        return 0


def _open_ms_reader(file_path: Path):
    """Open mzML or mzXML with the matching Pyteomics reader."""
    if file_path.suffix.lower() == ".mzxml":
        from pyteomics import mzxml
        return mzxml.read(str(file_path))
    from pyteomics import mzml
    return mzml.read(str(file_path))


def _mzml_rt_minutes(spectrum: dict) -> float | None:
    scans = ((spectrum.get("scanList") or {}).get("scan") or [])
    scan = scans[0] if scans else {}
    raw = (
        scan.get("scan start time")
        if scan.get("scan start time") is not None
        else spectrum.get("retentionTime")
    )
    if raw is None:
        raw = spectrum.get("scan start time")
    if raw is None:
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    unit = str(getattr(raw, "unit_info", "") or "").lower()
    if "second" in unit:
        return value / 60.0
    return value


def _mzml_scan_polarity(spectrum: dict, fallback: str | None = None) -> str | None:
    if "positive scan" in spectrum:
        return "POS"
    if "negative scan" in spectrum:
        return "NEG"
    mzxml_value = str(spectrum.get("polarity") or "").strip()
    if mzxml_value == "+":
        return "POS"
    if mzxml_value == "-":
        return "NEG"
    value = str(fallback or "").strip().upper()
    return value if value in {"POS", "NEG"} else None


def _target_intensity(mz_array, intensity_array, target_mz: float,
                      tolerance_da: float, tolerance_ppm: float):
    if mz_array is None or intensity_array is None or len(mz_array) == 0:
        return None
    window = min(float(tolerance_da), abs(float(target_mz)) * float(tolerance_ppm) * 1e-6)
    left = int(np.searchsorted(mz_array, target_mz - window, side="left"))
    right = int(np.searchsorted(mz_array, target_mz + window, side="right"))
    if right <= left:
        return None
    local = np.asarray(intensity_array[left:right], dtype=float)
    if local.size == 0 or not np.any(np.isfinite(local)):
        return None
    rel_idx = int(np.nanargmax(local))
    intensity = float(local[rel_idx])
    if not np.isfinite(intensity) or intensity <= 0:
        return None
    return float(mz_array[left + rel_idx]), intensity


def _base_assessment(status: str = "no_blank_feature") -> dict:
    return {
        "blank_status": status,
        "blank_feature_matched": False,
        "blank_excluded": False,
        "blank_flagged": False,
        "blank_clean_view_excluded": False,
        "blank_legacy_below_ratio": False,
        "blank_carryover_suspected": False,
        "blank_filter_policy": DEFAULT_BLANK_FILTER_POLICY,
        "blank_evidence_role": "unspecified",
        "blank_evidence_roles": [],
        "blank_evidence_category": "",
        "blank_role_evidence": {},
        "blank_evidence_source": "",
        "blank_sample_intensity": None,
        "blank_intensity": None,
        "blank_sample_ratio": None,
        "blank_mz_error_ppm": None,
        "blank_rt_delta_minutes": None,
        "blank_quantitation_mode": "",
        "blank_sample_xic_area": None,
        "blank_xic_area": None,
        "blank_sample_xic_scans": 0,
        "blank_xic_scans": 0,
        "blank_replicates_total": 1,
        "blank_replicates_applicable": 0,
        "blank_replicates_detected": 0,
        "blank_detection_fraction": None,
        "blank_required_detection_count": 1,
        "blank_replicate_applicable": True,
    }


def normalize_blank_role(value) -> str:
    role = str(value or "unspecified").strip().lower()
    return role if role in BLANK_ROLES else "unspecified"


def _normalize_blank_filter_policy(value) -> str:
    policy = str(value or DEFAULT_BLANK_FILTER_POLICY).strip().lower()
    return policy if policy in BLANK_FILTER_POLICIES else DEFAULT_BLANK_FILTER_POLICY


def _apply_blank_evidence_policy(
    assessment: dict,
    *,
    filter_policy: str,
    blank_role: str,
    ratio_threshold: float,
    clean_max_ratio: float,
    clean_min_replicates: int,
) -> dict:
    """Separate quantitative blank evidence from result visibility.

    ``soft_evidence`` never suppresses library search.  A feature is omitted
    only from the optional clean view when an explicitly typed background blank
    reproducibly contains at least as much signal as the sample.  Carryover and
    untyped blanks remain audit evidence and can never hide an identification.
    ``hard_legacy`` preserves the former sample/blank fold-change gate, except
    that a blank explicitly labelled as carryover is still audit-only.
    """
    policy = _normalize_blank_filter_policy(filter_policy)
    role = normalize_blank_role(blank_role)
    assessment["blank_filter_policy"] = policy
    assessment["blank_evidence_role"] = role
    assessment["blank_evidence_roles"] = (
        [role] if assessment.get("blank_feature_matched") else []
    )
    assessment["blank_clean_view_excluded"] = False
    assessment["blank_flagged"] = False
    assessment["blank_carryover_suspected"] = False

    if not assessment.get("blank_feature_matched"):
        assessment["blank_excluded"] = False
        status = str(assessment.get("blank_status") or "")
        assessment["blank_evidence_category"] = (
            "not_detected" if status == "no_blank_feature" else "unevaluable"
        )
        return assessment

    ratio = _positive_float(assessment.get("blank_sample_ratio"))
    if ratio is None:
        status = str(assessment.get("blank_status") or "")
        assessment["blank_evidence_category"] = (
            "low_blank_prevalence"
            if status == "retained_low_blank_prevalence"
            else "unevaluable"
        )
        assessment["blank_excluded"] = False
        if role == "carryover_blank":
            assessment["blank_flagged"] = True
            assessment["blank_carryover_suspected"] = True
        return assessment

    threshold = max(float(ratio_threshold), 1.0)
    clean_ratio = max(float(clean_max_ratio), 0.0)
    if ratio >= threshold:
        category = "sample_enriched"
    elif ratio > clean_ratio:
        category = "blank_associated"
    else:
        category = "blank_dominant"
    assessment["blank_evidence_category"] = category
    assessment["blank_legacy_below_ratio"] = ratio < threshold

    if role == "carryover_blank":
        assessment["blank_excluded"] = False
        assessment["blank_clean_view_excluded"] = False
        assessment["blank_flagged"] = True
        assessment["blank_carryover_suspected"] = True
        assessment["blank_status"] = "retained_carryover_evidence"
        return assessment

    if policy == "hard_legacy":
        assessment["blank_excluded"] = ratio < threshold
        assessment["blank_flagged"] = assessment["blank_excluded"]
        assessment["blank_status"] = (
            "excluded_below_ratio"
            if assessment["blank_excluded"]
            else "retained_above_ratio"
        )
        return assessment

    # Default V7 path: search everything, then expose evidence-aware views.
    assessment["blank_excluded"] = False
    assessment["blank_flagged"] = category != "sample_enriched"
    if category == "sample_enriched":
        assessment["blank_status"] = "retained_sample_enriched"
    elif category == "blank_associated":
        assessment["blank_status"] = "retained_blank_associated"
    else:
        assessment["blank_status"] = "retained_blank_dominant"
        detected = int(assessment.get("blank_replicates_detected") or 0)
        required = int(assessment.get("blank_required_detection_count") or 0)
        assessment["blank_clean_view_excluded"] = bool(
            role in BACKGROUND_BLANK_ROLES
            and detected >= max(int(clean_min_replicates), 1)
            and detected >= max(required, 1)
        )
    return assessment


def _finalize_blank_assessments(
    assessments: Sequence[dict],
    *,
    filter_policy: str,
    blank_role: str,
    ratio_threshold: float,
    clean_max_ratio: float,
    clean_min_replicates: int,
) -> list[dict]:
    return [
        _apply_blank_evidence_policy(
            assessment,
            filter_policy=filter_policy,
            blank_role=blank_role,
            ratio_threshold=ratio_threshold,
            clean_max_ratio=clean_max_ratio,
            clean_min_replicates=clean_min_replicates,
        )
        for assessment in assessments
    ]


def _sample_intensity(spec: dict) -> float | None:
    return (
        _positive_float(spec.get("precursor_ms1_intensity"))
        or _positive_float(spec.get("precursor_intensity"))
    )


def assess_blank_features(
    sample_spectra: Sequence[dict],
    blank_profile: BlankProfile,
    *,
    sample_path: str | Path | None = None,
    sample_xic: dict[int, dict] | None = None,
    sample_rt_unit: str = "minutes",
    sample_to_blank_ratio: float = DEFAULT_SAMPLE_TO_BLANK_RATIO,
    rt_tolerance_minutes: float = DEFAULT_FEATURE_RT_TOLERANCE_MINUTES,
    precursor_tolerance_da: float = DEFAULT_PRECURSOR_TOLERANCE_DA,
    precursor_tolerance_ppm: float = DEFAULT_PRECURSOR_TOLERANCE_PPM,
    filter_policy: str = "hard_legacy",
    blank_role: str = "unspecified",
    clean_max_ratio: float = DEFAULT_BLANK_CLEAN_MAX_RATIO,
    clean_min_replicates: int = DEFAULT_BLANK_CLEAN_MIN_REPLICATES,
) -> list[dict]:
    """Assess each sample MS2 against a same-polarity blank LC-MS feature.

    mzML uses the same MS1 extracted-ion chromatogram (XIC) statistic on sample
    and blank when the sample raw file is available.  This avoids the former
    asymmetric comparison between one preceding sample MS1 scan and the maximum
    blank signal across an RT window.  MGF has no MS1 and therefore receives a
    clearly labelled, conservative MS2-context fallback.
    """
    assessments = [_base_assessment() for _ in sample_spectra]
    targets = []
    for index, spec in enumerate(sample_spectra):
        precursor = _positive_float(spec.get("precursor_mz"))
        rt = _rt_to_minutes(
            spec.get("retention_time"),
            spec.get("retention_time_unit") or sample_rt_unit,
        )
        polarity = str(spec.get("polarity") or spec.get("ion_mode") or "").strip().upper()
        if precursor is None:
            assessments[index]["blank_status"] = "sample_precursor_missing"
            continue
        if rt is None:
            assessments[index]["blank_status"] = "sample_rt_missing"
            continue
        targets.append((rt, index, precursor, polarity if polarity in {"POS", "NEG"} else None))

    source = blank_profile.source_path
    if source and source.suffix.lower() in {".mzml", ".mzxml", ".xml"} and blank_profile.ms1_scan_count:
        if sample_xic is None and sample_path is not None:
            sample_source = Path(sample_path)
            if sample_source.suffix.lower() in {".mzml", ".mzxml", ".xml"}:
                sample_xic = _extract_mzml_xic_targets(
                    sample_source,
                    targets,
                    rt_tolerance_minutes=rt_tolerance_minutes,
                    precursor_tolerance_da=precursor_tolerance_da,
                    precursor_tolerance_ppm=precursor_tolerance_ppm,
                )
        _assess_mzml_ms1_targets(
            source, targets, sample_spectra, assessments,
            sample_xic=sample_xic,
            sample_to_blank_ratio=sample_to_blank_ratio,
            rt_tolerance_minutes=rt_tolerance_minutes,
            precursor_tolerance_da=precursor_tolerance_da,
            precursor_tolerance_ppm=precursor_tolerance_ppm,
        )
        return _finalize_blank_assessments(
            assessments,
            filter_policy=filter_policy,
            blank_role=blank_role,
            ratio_threshold=sample_to_blank_ratio,
            clean_max_ratio=clean_max_ratio,
            clean_min_replicates=clean_min_replicates,
        )

    blank_profile.rt_tolerance_minutes = float(rt_tolerance_minutes)
    for rt, index, precursor, polarity in targets:
        if (polarity and blank_profile.available_polarities
                and polarity not in blank_profile.available_polarities):
            assessments[index]["blank_status"] = "blank_polarity_not_applicable"
            assessments[index]["blank_replicate_applicable"] = False
            continue
        assessments[index]["blank_replicates_applicable"] = 1
        contexts = blank_profile._matching_contexts(
            precursor, rt, "minutes", polarity,
        )
        if not contexts:
            assessments[index]["blank_status"] = "no_blank_ms2_context"
            continue
        sample_value = _sample_intensity(sample_spectra[index])
        blank_values = [item.precursor_intensity for item in contexts if item.precursor_intensity]
        blank_value = max(blank_values) if blank_values else None
        _finish_assessment(
            assessments[index], sample_value, blank_value,
            sample_to_blank_ratio, "ms2_precursor_intensity",
            mz_error_ppm=min(abs(precursor - item.precursor_mz) / precursor * 1e6 for item in contexts),
            rt_delta=min(abs(rt - item.retention_time_minutes)
                         for item in contexts if item.retention_time_minutes is not None)
                         if any(item.retention_time_minutes is not None for item in contexts) else None,
        )
        assessments[index].update({
            "blank_quantitation_mode": "mgf_ms2_precursor_intensity_fallback",
            "blank_replicates_detected": 1,
            "blank_detection_fraction": 1.0,
        })
    return _finalize_blank_assessments(
        assessments,
        filter_policy=filter_policy,
        blank_role=blank_role,
        ratio_threshold=sample_to_blank_ratio,
        clean_max_ratio=clean_max_ratio,
        clean_min_replicates=clean_min_replicates,
    )


def _trapezoid_area(points: Sequence[tuple[float, float]]) -> float:
    """Return a deterministic XIC area without relying on deprecated np.trapz."""
    if len(points) < 2:
        return 0.0
    ordered = sorted((float(rt), max(float(value), 0.0)) for rt, value in points)
    return float(sum(
        (right_rt - left_rt) * (left_value + right_value) * 0.5
        for (left_rt, left_value), (right_rt, right_value)
        in zip(ordered, ordered[1:])
        if right_rt >= left_rt
    ))


def _extract_mzml_xic_targets(
    mzml_path: Path,
    targets,
    *,
    rt_tolerance_minutes,
    precursor_tolerance_da,
    precursor_tolerance_ppm,
) -> dict[int, dict]:
    """Extract identical apex/area XIC evidence for a collection of MS2 targets."""
    try:
        from utils.mzml_io import _file_polarity
    except Exception as exc:
        raise RuntimeError(f"mzML XIC reader unavailable: {exc}") from exc

    ordered = sorted(targets)
    target_rts = [item[0] for item in ordered]
    file_polarity = _file_polarity(mzml_path)
    traces: dict[int, list[tuple[float, float]]] = {}
    detected_scans: dict[int, int] = {}
    best: dict[int, tuple[float, float, float, float, float]] = {}
    applicable: set[int] = set()

    with _open_ms_reader(mzml_path) as reader:
        for spectrum in reader:
            if spectrum.get("ms level") != 1:
                continue
            rt = _mzml_rt_minutes(spectrum)
            if rt is None:
                continue
            polarity = _mzml_scan_polarity(spectrum, file_polarity)
            left = bisect_left(target_rts, rt - rt_tolerance_minutes)
            right = bisect_right(target_rts, rt + rt_tolerance_minutes)
            if right <= left:
                continue
            mz_array = spectrum.get("m/z array")
            intensity_array = spectrum.get("intensity array")
            for target_rt, index, precursor, target_polarity in ordered[left:right]:
                if target_polarity and polarity and polarity != target_polarity:
                    continue
                if (target_polarity and file_polarity in {"POS", "NEG"}
                        and file_polarity != target_polarity):
                    continue
                applicable.add(index)
                found = _target_intensity(
                    mz_array, intensity_array, precursor,
                    precursor_tolerance_da, precursor_tolerance_ppm,
                )
                intensity = float(found[1]) if found else 0.0
                traces.setdefault(index, []).append((rt, intensity))
                if found is None:
                    continue
                observed_mz = float(found[0])
                detected_scans[index] = detected_scans.get(index, 0) + 1
                previous = best.get(index)
                if previous is None or intensity > previous[0]:
                    best[index] = (
                        intensity, observed_mz, rt, target_rt, precursor,
                    )

    evidence: dict[int, dict] = {}
    for _, index, _, _ in targets:
        points = traces.get(index, [])
        found = best.get(index)
        evidence[index] = {
            "applicable": index in applicable,
            "apex": found[0] if found else None,
            "observed_mz": found[1] if found else None,
            "apex_rt": found[2] if found else None,
            "target_rt": found[3] if found else None,
            "target_mz": found[4] if found else None,
            "area": _trapezoid_area(points),
            "scans": len(points),
            "detected_scans": detected_scans.get(index, 0),
        }
    return evidence


def _assess_mzml_ms1_targets(
    blank_path: Path,
    targets,
    sample_spectra,
    assessments,
    *,
    sample_xic,
    sample_to_blank_ratio,
    rt_tolerance_minutes,
    precursor_tolerance_da,
    precursor_tolerance_ppm,
):
    blank_xic = _extract_mzml_xic_targets(
        blank_path,
        targets,
        rt_tolerance_minutes=rt_tolerance_minutes,
        precursor_tolerance_da=precursor_tolerance_da,
        precursor_tolerance_ppm=precursor_tolerance_ppm,
    )
    for _, index, _, _ in targets:
        blank_found = blank_xic.get(index) or {}
        assessments[index]["blank_replicate_applicable"] = bool(
            blank_found.get("applicable")
        )
        assessments[index]["blank_replicates_applicable"] = int(
            bool(blank_found.get("applicable"))
        )
        if not blank_found.get("applicable"):
            assessments[index]["blank_status"] = "blank_polarity_not_applicable"
            continue
        if not blank_found.get("apex"):
            assessments[index]["blank_status"] = "no_blank_feature"
            continue
        sample_found = (sample_xic or {}).get(index) or {}
        use_paired_xic = bool(sample_xic is not None)
        sample_value = (
            sample_found.get("apex") if use_paired_xic
            else _sample_intensity(sample_spectra[index])
        )
        blank_value = blank_found.get("apex")
        observed_mz = blank_found.get("observed_mz")
        blank_rt = blank_found.get("apex_rt")
        sample_rt = blank_found.get("target_rt")
        precursor = blank_found.get("target_mz")
        _finish_assessment(
            assessments[index], sample_value, blank_value,
            sample_to_blank_ratio,
            "paired_ms1_xic_apex" if use_paired_xic else "blank_ms1",
            mz_error_ppm=abs(precursor - observed_mz) / precursor * 1e6,
            rt_delta=abs(sample_rt - blank_rt),
        )
        assessments[index].update({
            "blank_quantitation_mode": (
                "symmetric_xic_apex" if use_paired_xic
                else "legacy_sample_scan_vs_blank_xic_apex"
            ),
            "blank_sample_xic_area": (
                round(float(sample_found.get("area") or 0.0), 6)
                if use_paired_xic else None
            ),
            "blank_xic_area": round(float(blank_found.get("area") or 0.0), 6),
            "blank_sample_xic_scans": int(sample_found.get("scans") or 0),
            "blank_xic_scans": int(blank_found.get("scans") or 0),
            "blank_replicates_detected": 1,
            "blank_detection_fraction": 1.0,
        })


def _finish_assessment(assessment: dict, sample_value, blank_value,
                       ratio_threshold, source, *, mz_error_ppm, rt_delta):
    assessment.update({
        "blank_feature_matched": True,
        "blank_evidence_source": source,
        "blank_sample_intensity": float(sample_value) if sample_value else None,
        "blank_intensity": float(blank_value) if blank_value else None,
        "blank_mz_error_ppm": round(float(mz_error_ppm), 4) if mz_error_ppm is not None else None,
        "blank_rt_delta_minutes": round(float(rt_delta), 5) if rt_delta is not None else None,
    })
    if not sample_value:
        assessment["blank_status"] = "sample_intensity_missing"
        return
    if not blank_value:
        assessment["blank_status"] = "blank_intensity_missing"
        return
    ratio = float(sample_value) / float(blank_value)
    assessment["blank_sample_ratio"] = round(ratio, 6)
    if ratio < float(ratio_threshold):
        assessment["blank_status"] = "excluded_below_ratio"
        assessment["blank_excluded"] = True
    else:
        assessment["blank_status"] = "retained_above_ratio"


def blank_polarity_coverage(
    sample_spectra: Sequence[dict],
    blank_profiles: Sequence[BlankProfile],
) -> dict:
    """Summarize whether supplied blanks cover every known sample polarity.

    A known POS sample cannot be filtered with known NEG blanks (or vice
    versa).  Profiles without polarity metadata remain conservatively usable
    because rejecting them would change the historical MGF fallback; their
    uncertainty is exposed in ``unknown_blank_profiles``.
    """
    sample_polarities = {
        str(spec.get("polarity") or spec.get("ion_mode") or "").strip().upper()
        for spec in sample_spectra
    } & {"POS", "NEG"}
    blank_polarities: set[str] = set()
    unknown_blank_profiles = 0
    for profile in blank_profiles:
        declared = {
            str(value).strip().upper()
            for value in (getattr(profile, "available_polarities", ()) or ())
        } & {"POS", "NEG"}
        if declared:
            blank_polarities.update(declared)
        else:
            unknown_blank_profiles += 1

    missing = (
        sample_polarities - blank_polarities
        if sample_polarities and unknown_blank_profiles == 0
        else set()
    )
    return {
        "sample_polarities": sorted(sample_polarities),
        "blank_polarities": sorted(blank_polarities),
        "unknown_blank_profiles": unknown_blank_profiles,
        "missing_sample_polarities": sorted(missing),
        "applicable": not missing,
    }


def assign_blank_profiles_by_polarity(
    sample_spectra: Sequence[dict],
    blank_profiles: Sequence[BlankProfile],
) -> dict:
    """Partition uploaded blanks for one sample using resolved ion polarity.

    Known opposite-polarity blanks are excluded from the sample's replicate
    set before any blank statistics are calculated.  A profile without
    polarity metadata is retained for backwards-compatible MGF evidence, but
    is reported separately so callers can surface that uncertainty.

    ``web.server`` calls this only after resolving every sample spectrum to
    POS or NEG.  The unresolved-sample fallback below deliberately assigns all
    profiles rather than silently discarding user-supplied evidence.
    """
    sample_polarities = {
        str(spec.get("polarity") or spec.get("ion_mode") or "").strip().upper()
        for spec in sample_spectra
    } & {"POS", "NEG"}
    assigned_indices: list[int] = []
    skipped_indices: list[int] = []
    unknown_assigned_indices: list[int] = []
    profile_polarities: list[list[str]] = []

    for index, profile in enumerate(blank_profiles):
        declared = {
            str(value).strip().upper()
            for value in (getattr(profile, "available_polarities", ()) or ())
        } & {"POS", "NEG"}
        profile_polarities.append(sorted(declared))
        if declared and sample_polarities and declared.isdisjoint(sample_polarities):
            skipped_indices.append(index)
            continue
        assigned_indices.append(index)
        if not declared:
            unknown_assigned_indices.append(index)

    assigned_polarities = {
        polarity
        for index in assigned_indices
        for polarity in profile_polarities[index]
    }
    if sample_polarities:
        assigned_polarities &= sample_polarities
    skipped_polarities = sorted({
        polarity
        for index in skipped_indices
        for polarity in profile_polarities[index]
    })
    if not sample_polarities:
        assignment_mode = "unresolved_sample_polarity_fallback"
    elif unknown_assigned_indices:
        assignment_mode = (
            "automatic_by_sample_polarity_with_unknown_fallback"
        )
    else:
        assignment_mode = "automatic_by_sample_polarity"
    return {
        "assignment_mode": assignment_mode,
        "sample_polarities": sorted(sample_polarities),
        "profile_polarities": profile_polarities,
        "assigned_indices": assigned_indices,
        "skipped_indices": skipped_indices,
        "unknown_assigned_indices": unknown_assigned_indices,
        "assigned_blank_polarities": sorted(assigned_polarities),
        "skipped_blank_polarities": skipped_polarities,
    }


def assess_blank_replicates(
    sample_spectra: Sequence[dict],
    blank_profiles: Sequence[BlankProfile],
    *,
    sample_path: str | Path | None = None,
    sample_rt_unit: str = "minutes",
    sample_to_blank_ratio: float = DEFAULT_SAMPLE_TO_BLANK_RATIO,
    rt_tolerance_minutes: float = DEFAULT_FEATURE_RT_TOLERANCE_MINUTES,
    precursor_tolerance_da: float = DEFAULT_PRECURSOR_TOLERANCE_DA,
    precursor_tolerance_ppm: float = DEFAULT_PRECURSOR_TOLERANCE_PPM,
    min_detection_fraction: float = DEFAULT_BLANK_MIN_DETECTION_FRACTION,
    filter_policy: str = "hard_legacy",
    blank_role: str = "unspecified",
    clean_max_ratio: float = DEFAULT_BLANK_CLEAN_MAX_RATIO,
    clean_min_replicates: int = DEFAULT_BLANK_CLEAN_MIN_REPLICATES,
) -> list[dict]:
    """Aggregate replicate blanks using prevalence and median blank abundance.

    A one-off contaminant spike in one blank cannot exclude a feature when
    multiple blanks were supplied.  For two replicates both must contain the
    feature; for three or more, at least ``min_detection_fraction`` (default
    50%) and at least two detections are required.
    """
    profiles = list(blank_profiles)
    if not profiles:
        return [_base_assessment("not_requested") for _ in sample_spectra]

    sample_xic = None
    sample_source = Path(sample_path) if sample_path is not None else None
    if sample_source and sample_source.suffix.lower() in {".mzml", ".mzxml", ".xml"}:
        targets = []
        for index, spec in enumerate(sample_spectra):
            precursor = _positive_float(spec.get("precursor_mz"))
            rt = _rt_to_minutes(
                spec.get("retention_time"),
                spec.get("retention_time_unit") or sample_rt_unit,
            )
            polarity = str(
                spec.get("polarity") or spec.get("ion_mode") or ""
            ).strip().upper()
            if precursor is not None and rt is not None:
                targets.append((
                    rt, index, precursor,
                    polarity if polarity in {"POS", "NEG"} else None,
                ))
        if targets:
            sample_xic = _extract_mzml_xic_targets(
                sample_source,
                targets,
                rt_tolerance_minutes=rt_tolerance_minutes,
                precursor_tolerance_da=precursor_tolerance_da,
                precursor_tolerance_ppm=precursor_tolerance_ppm,
            )

    by_profile = [
        assess_blank_features(
            sample_spectra,
            profile,
            sample_path=sample_source,
            sample_xic=sample_xic,
            sample_rt_unit=sample_rt_unit,
            sample_to_blank_ratio=sample_to_blank_ratio,
            rt_tolerance_minutes=rt_tolerance_minutes,
            precursor_tolerance_da=precursor_tolerance_da,
            precursor_tolerance_ppm=precursor_tolerance_ppm,
            filter_policy=filter_policy,
            blank_role=blank_role,
            clean_max_ratio=clean_max_ratio,
            clean_min_replicates=clean_min_replicates,
        )
        for profile in profiles
    ]
    if len(profiles) == 1:
        return by_profile[0]

    combined: list[dict] = []
    min_fraction = min(max(float(min_detection_fraction), 0.0), 1.0)
    for index in range(len(sample_spectra)):
        rows = [assessments[index] for assessments in by_profile]
        applicable = [
            row for row in rows if row.get("blank_replicate_applicable", True)
        ]
        detected = [
            row for row in applicable if row.get("blank_feature_matched")
        ]
        n_applicable = len(applicable)
        n_detected = len(detected)
        required = (
            max(2, int(math.ceil(n_applicable * min_fraction)))
            if n_applicable > 1 else 1
        )
        result = _base_assessment()
        result.update({
            "blank_replicates_total": len(profiles),
            "blank_replicates_applicable": n_applicable,
            "blank_replicates_detected": n_detected,
            "blank_detection_fraction": (
                round(n_detected / n_applicable, 6) if n_applicable else None
            ),
            "blank_required_detection_count": required,
            "blank_replicate_applicable": bool(n_applicable),
        })
        if not n_applicable:
            result["blank_status"] = "blank_polarity_not_applicable"
            combined.append(result)
            continue
        if not detected:
            result["blank_status"] = "no_blank_feature"
            combined.append(result)
            continue

        result["blank_feature_matched"] = True
        sources = sorted({
            str(row.get("blank_evidence_source") or "") for row in detected
            if row.get("blank_evidence_source")
        })
        result["blank_evidence_source"] = "replicate_consensus:" + "+".join(sources)
        result["blank_quantitation_mode"] = "replicate_prevalence_median_apex"
        result["blank_mz_error_ppm"] = min(
            (row["blank_mz_error_ppm"] for row in detected
             if row.get("blank_mz_error_ppm") is not None),
            default=None,
        )
        result["blank_rt_delta_minutes"] = min(
            (row["blank_rt_delta_minutes"] for row in detected
             if row.get("blank_rt_delta_minutes") is not None),
            default=None,
        )
        if n_detected < required:
            result["blank_status"] = "retained_low_blank_prevalence"
            combined.append(result)
            continue

        sample_values = [
            float(row["blank_sample_intensity"]) for row in detected
            if _positive_float(row.get("blank_sample_intensity"))
        ]
        blank_values = [
            float(row["blank_intensity"]) for row in detected
            if _positive_float(row.get("blank_intensity"))
        ]
        sample_areas = [
            float(row["blank_sample_xic_area"]) for row in detected
            if row.get("blank_sample_xic_area") is not None
        ]
        blank_areas = [
            float(row["blank_xic_area"]) for row in detected
            if row.get("blank_xic_area") is not None
        ]
        sample_value = float(np.median(sample_values)) if sample_values else None
        blank_value = float(np.median(blank_values)) if blank_values else None
        _finish_assessment(
            result,
            sample_value,
            blank_value,
            sample_to_blank_ratio,
            result["blank_evidence_source"],
            mz_error_ppm=result["blank_mz_error_ppm"],
            rt_delta=result["blank_rt_delta_minutes"],
        )
        result.update({
            "blank_quantitation_mode": "replicate_prevalence_median_apex",
            "blank_sample_xic_area": (
                round(float(np.median(sample_areas)), 6) if sample_areas else None
            ),
            "blank_xic_area": (
                round(float(np.median(blank_areas)), 6) if blank_areas else None
            ),
            "blank_sample_xic_scans": max(
                (int(row.get("blank_sample_xic_scans") or 0) for row in detected),
                default=0,
            ),
            "blank_xic_scans": max(
                (int(row.get("blank_xic_scans") or 0) for row in detected),
                default=0,
            ),
            "blank_replicates_total": len(profiles),
            "blank_replicates_applicable": n_applicable,
            "blank_replicates_detected": n_detected,
            "blank_detection_fraction": round(n_detected / n_applicable, 6),
            "blank_required_detection_count": required,
        })
        combined.append(result)
    return _finalize_blank_assessments(
        combined,
        filter_policy=filter_policy,
        blank_role=blank_role,
        ratio_threshold=sample_to_blank_ratio,
        clean_max_ratio=clean_max_ratio,
        clean_min_replicates=clean_min_replicates,
    )


def _compact_role_evidence(row: dict) -> dict:
    return {
        key: row.get(key)
        for key in (
            "blank_status",
            "blank_feature_matched",
            "blank_flagged",
            "blank_clean_view_excluded",
            "blank_carryover_suspected",
            "blank_evidence_category",
            "blank_sample_intensity",
            "blank_intensity",
            "blank_sample_ratio",
            "blank_detection_fraction",
            "blank_replicates_total",
            "blank_replicates_applicable",
            "blank_replicates_detected",
            "blank_required_detection_count",
        )
    }


def assess_blank_role_groups(
    sample_spectra: Sequence[dict],
    blank_profiles: Sequence[BlankProfile],
    blank_roles: Sequence[str] | None = None,
    *,
    sample_path: str | Path | None = None,
    sample_rt_unit: str = "minutes",
    sample_to_blank_ratio: float = DEFAULT_SAMPLE_TO_BLANK_RATIO,
    rt_tolerance_minutes: float = DEFAULT_FEATURE_RT_TOLERANCE_MINUTES,
    precursor_tolerance_da: float = DEFAULT_PRECURSOR_TOLERANCE_DA,
    precursor_tolerance_ppm: float = DEFAULT_PRECURSOR_TOLERANCE_PPM,
    min_detection_fraction: float = DEFAULT_BLANK_MIN_DETECTION_FRACTION,
    filter_policy: str = DEFAULT_BLANK_FILTER_POLICY,
    clean_max_ratio: float = DEFAULT_BLANK_CLEAN_MAX_RATIO,
    clean_min_replicates: int = DEFAULT_BLANK_CLEAN_MIN_REPLICATES,
) -> list[dict]:
    """Assess blank replicates within role, then merge role-aware evidence.

    Method, extraction, solvent, carryover and unspecified blanks are never
    pooled as interchangeable replicates.  This prevents a post-injection
    carryover file from acting as a procedural-background replicate.
    """
    profiles = list(blank_profiles)
    if not profiles:
        return [_base_assessment("not_requested") for _ in sample_spectra]
    roles = list(blank_roles or [])
    if len(roles) < len(profiles):
        roles.extend(["unspecified"] * (len(profiles) - len(roles)))
    roles = [normalize_blank_role(role) for role in roles[:len(profiles)]]

    grouped: dict[str, list[BlankProfile]] = {}
    for profile, role in zip(profiles, roles):
        grouped.setdefault(role, []).append(profile)

    by_role = {
        role: assess_blank_replicates(
            sample_spectra,
            role_profiles,
            sample_path=sample_path,
            sample_rt_unit=sample_rt_unit,
            sample_to_blank_ratio=sample_to_blank_ratio,
            rt_tolerance_minutes=rt_tolerance_minutes,
            precursor_tolerance_da=precursor_tolerance_da,
            precursor_tolerance_ppm=precursor_tolerance_ppm,
            min_detection_fraction=min_detection_fraction,
            filter_policy=filter_policy,
            blank_role=role,
            clean_max_ratio=clean_max_ratio,
            clean_min_replicates=clean_min_replicates,
        )
        for role, role_profiles in grouped.items()
    }

    def primary_priority(item: tuple[str, dict]) -> tuple:
        role, row = item
        matched = bool(row.get("blank_feature_matched"))
        clean_excluded = bool(row.get("blank_clean_view_excluded"))
        role_rank = (
            0 if role in BACKGROUND_BLANK_ROLES
            else 1 if role == "carryover_blank"
            else 2
        )
        ratio = _positive_float(row.get("blank_sample_ratio"))
        return (
            0 if clean_excluded else 1,
            0 if matched else 1,
            role_rank,
            ratio if ratio is not None else float("inf"),
            role,
        )

    combined: list[dict] = []
    for index in range(len(sample_spectra)):
        role_rows = {
            role: assessments[index] for role, assessments in by_role.items()
        }
        primary_role, primary = min(
            role_rows.items(), key=primary_priority
        )
        result = dict(primary)
        matched_roles = sorted(
            role for role, row in role_rows.items()
            if row.get("blank_feature_matched")
        )
        result.update({
            "blank_filter_policy": _normalize_blank_filter_policy(filter_policy),
            "blank_evidence_role": primary_role,
            "blank_evidence_roles": matched_roles,
            "blank_role_evidence": {
                role: _compact_role_evidence(row)
                for role, row in sorted(role_rows.items())
            },
            "blank_replicates_total": len(profiles),
            "blank_replicates_applicable": sum(
                int(row.get("blank_replicates_applicable") or 0)
                for row in role_rows.values()
            ),
            "blank_replicates_detected": sum(
                int(row.get("blank_replicates_detected") or 0)
                for row in role_rows.values()
            ),
            "blank_flagged": any(
                bool(row.get("blank_flagged")) for row in role_rows.values()
            ),
            "blank_clean_view_excluded": any(
                bool(row.get("blank_clean_view_excluded"))
                for row in role_rows.values()
            ),
            "blank_carryover_suspected": bool(
                (role_rows.get("carryover_blank") or {}).get(
                    "blank_feature_matched"
                )
            ),
        })
        if result["blank_carryover_suspected"]:
            result["blank_flagged"] = True
        if _normalize_blank_filter_policy(filter_policy) == "hard_legacy":
            result["blank_excluded"] = any(
                bool(row.get("blank_excluded")) for row in role_rows.values()
            )
            if result["blank_excluded"]:
                result["blank_status"] = "excluded_below_ratio"
        else:
            result["blank_excluded"] = False
            if (
                primary_role == "carryover_blank"
                and primary.get("blank_feature_matched")
            ):
                result["blank_status"] = "retained_carryover_evidence"
        combined.append(result)
    return combined
