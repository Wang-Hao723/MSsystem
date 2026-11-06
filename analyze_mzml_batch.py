#!/usr/bin/env python3
"""
分析 mzML / MGF 文件,匹配 local_library.db 质谱库,余弦相似度 >= 0.7。
支持空白溶剂峰扣除(基于 BlankProfile)。

输出格式参考 /Users/zizi/Desktop/lc-ms分析结果/:
- 总体摘要.csv
- <file>_结果.csv   (1 行/匹配,字段含 SELFIES、Confidence)
- <file>_报告.csv   (1 行/匹配,字段含 Peaks_Count、Mass_PPM)
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import re
import sqlite3
import time
import warnings
from pathlib import Path

# Suppress RDKit warnings early
try:
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    pass
from typing import Literal

import numpy as np
from matchms import Spectrum
from matchms.importing import load_from_mgf, load_from_mzml
from matchms.similarity import ModifiedCosineGreedy as ModifiedCosine

from config.settings import DB_PATH
from matching.search import (
    clean_peaks as _canonical_clean_peaks,
    load_library as _canonical_load_library,
    search as _canonical_search,
)

# Blank subtraction
try:
    from matching.blank_subtraction import (
        BlankProfile,
        load_blank_mgf,
        subtract_blank,
    )
    BLANK_SUBTRACTION_AVAILABLE = True
except ImportError:
    BLANK_SUBTRACTION_AVAILABLE = False

# 专家系统(可选)
try:
    from matching.expert_scorer import expert_evaluate
    from matching.expert_logger import ExpertLogger
    EXPERT_AVAILABLE = True
except ImportError:
    EXPERT_AVAILABLE = False

# MS2Classifier DL分类器(可选)
try:
    import os as _os
    _clf_path = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "models", "ms2classifier_v3_best.pt")
    from models.ms2_clf_wrapper import MS2ClfWrapper
    _clf_dev = "cuda" if __import__("torch").cuda.is_available() else "cpu"
    _clf_wrapper = MS2ClfWrapper(_clf_path, device=_clf_dev)
    MS2CLF_AVAILABLE = True
except Exception as _e:
    MS2CLF_AVAILABLE = False
    _clf_wrapper = None
warnings.filterwarnings("ignore")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("analyze_mzml_batch")


# ─────────────────────────────────────────────────────────────────
# 默认配置(CLI 参数可覆盖)
# ─────────────────────────────────────────────────────────────────
DEFAULT_DB_PATH = Path(DB_PATH)
DEFAULT_SAMPLES_DIR = Path("/Users/zizi/Desktop/彭晓扬复方粗提物")
DEFAULT_OUTPUT_DIR = Path("/Users/zizi/Desktop/彭晓扬复方粗提物_results")
DEFAULT_BLANK_DIR = Path("/Users/zizi/Desktop/kb2")
DEFAULT_BLANK_A = DEFAULT_BLANK_DIR / "Blank-A.mgf"
DEFAULT_BLANK_P = DEFAULT_BLANK_DIR / "Blank-P.mgf"

# Matching
MZ_TOLERANCE_DA = 0.02
MZ_TOLERANCE_PPM = 10.0
COSINE_THRESHOLD = 0.6
FINAL_SCORE_THRESHOLD = 0.4
MZ_BUCKET_DALTON = 0.1
TOP_N_PEAKS = 50
REL_THRESHOLD = 0.01
MAX_CANDIDATES_PER_QUERY = 500

# Blank
BLANK_MZ_TOL_DA = 0.01
BLANK_MZ_TOL_PPM = 10
DEFAULT_BLANK_MODE: Literal["absolute", "relative", "filter_only"] = "absolute"
DEFAULT_BLANK_INTENSITY_MODE: Literal["max", "mean"] = "max"
DEFAULT_BLANK_THRESHOLD = 0.1  # 强度阈值因子


# ─────────────────────────────────────────────────────────────────
# 全局状态
# ─────────────────────────────────────────────────────────────────
_lib_records = None
_mz_index = None
# 空白 profile 字典:键为 "A" / "P"(或自定义),值为 BlankProfile
_blank_profiles: dict[str, BlankProfile] = {}


# ─────────────────────────────────────────────────────────────────
# 库加载
# ─────────────────────────────────────────────────────────────────
def _safe_float(v):
    if v is None:
        return None
    try:
        f = float(v)
        return f if not (f != f) else None
    except (TypeError, ValueError):
        return None


def load_library(db_path: Path):
    """Load through the canonical V6 library gate.

    The historical CLI used a second SQL/search implementation that ignored
    ``source_tier`` and therefore could re-enable quarantined or non-
    experimental spectra.  Keeping this wrapper preserves the public CLI API
    while making Web and command-line analyses use exactly the same records.
    """
    global _lib_records, _mz_index
    records, mz_index = _canonical_load_library(Path(db_path))
    _lib_records, _mz_index = records, mz_index
    return records, mz_index


def _query_candidates(pmz, records, mz_index):
    bucket_size = MZ_BUCKET_DALTON
    tol_da = MZ_TOLERANCE_DA
    tol_ppm = MZ_TOLERANCE_PPM
    bucket = int(pmz / bucket_size)
    candidates = []
    for b in range(bucket - 50, bucket + 51):
        if b not in mz_index:
            continue
        for idx in mz_index[b]:
            ref_mz = records[idx]["precursor_mz"]
            diff = abs(pmz - ref_mz)
            ppm = diff / pmz * 1e6 if pmz > 0 else 999999
            if diff <= tol_da or ppm <= tol_ppm:
                candidates.append(idx)
    return candidates


def _clean_peaks(peaks, precursor_mz, top_n=TOP_N_PEAKS):
    if not peaks:
        return []
    mz_arr = np.array([p[0] for p in peaks], dtype=float)
    int_arr = np.array([p[1] for p in peaks], dtype=float)
    if mz_arr.size == 0 or int_arr.size == 0:
        return []
    max_int = float(np.max(int_arr))
    if max_int <= 0:
        return []
    rel = int_arr / max_int
    keep = (rel >= REL_THRESHOLD) & (mz_arr < (precursor_mz - 1.5))
    mz_f = mz_arr[keep]
    int_f = int_arr[keep]
    if mz_f.size == 0:
        return []
    top_idx = np.argsort(int_f)[::-1][:top_n]
    mz_t = mz_f[top_idx]
    int_t = int_f[top_idx]
    order = np.argsort(mz_t)
    return [[float(mz_t[i]), float(int_t[i])] for i in order]


def _score_cosine(query_peaks, lib_peaks, mz_tol=0.1, query_pmz=None):
    """Modified Cosine — 使用实际 precursor_mz"""
    if not query_peaks or not lib_peaks:
        return 0.0
    try:
        q_ints = np.array([p[1] for p in query_peaks], dtype=float)
        q_max = q_ints.max() if q_ints.size > 0 else 1.0
        if q_max > 0:
            q_ints = q_ints / q_max

        l_ints = np.array([p[1] for p in lib_peaks], dtype=float)
        l_max = l_ints.max() if l_ints.size > 0 else 1.0
        if l_max > 0:
            l_ints = l_ints / l_max

        q_pmz = query_pmz if query_pmz else max(p[0] for p in query_peaks)
        l_pmz = max(p[0] for p in lib_peaks) + 1.0

        q_spectrum = Spectrum(
            mz=np.array([p[0] for p in query_peaks], dtype=float),
            intensities=q_ints,
            metadata={"precursor_mz": q_pmz}
        )
        l_spectrum = Spectrum(
            mz=np.array([p[0] for p in lib_peaks], dtype=float),
            intensities=l_ints,
            metadata={"precursor_mz": l_pmz}
        )
        score = ModifiedCosine(tolerance=mz_tol).pair(q_spectrum, l_spectrum)
        arr = score.tolist()
        s_val = arr[0] if isinstance(arr, (list, tuple)) else float(arr)
        return float(s_val) if not math.isnan(s_val) else 0.0
    except Exception:
        return 0.0


def _mass_score(pmz, lib_pmz, ppm_tol=20, da_tol=0.02):
    """质量精度得分 — 自适应双重容差(ppm+Da)"""
    if pmz <= 0 or lib_pmz <= 0:
        return 0.0
    diff_da = abs(pmz - lib_pmz)
    diff_ppm = diff_da / pmz * 1e6
    # 第一层：ppm 精确匹配
    if diff_ppm <= ppm_tol:
        return max(0.0, 1.0 - diff_ppm / ppm_tol)
    # 第二层：Da 容差匹配（衰减得分）
    if diff_da <= da_tol:
        return max(0.1, 1.0 - diff_da / da_tol) * 0.6
    return 0.0


def _mass_ppm(pmz, lib_pmz):
    if pmz <= 0 or lib_pmz <= 0:
        return None
    return abs(pmz - lib_pmz) / pmz * 1e6


def _fragment_score(query_clean, lib_clean, mz_tol=0.02):
    """碎片重叠得分 — 使用容差匹配替代严格精确匹配"""
    if not query_clean or not lib_clean:
        return 0.0
    q_mz_arr = np.array([p[0] for p in query_clean], dtype=float)
    l_mz_arr = np.array([p[0] for p in lib_clean], dtype=float)
    if len(l_mz_arr) == 0:
        return 0.0
    matched = 0
    q_used = set()
    for lmz in l_mz_arr:
        diffs = np.abs(q_mz_arr - lmz)
        sorted_indices = np.argsort(diffs)
        for idx in sorted_indices:
            if diffs[idx] <= mz_tol and idx not in q_used:
                matched += 1
                q_used.add(idx)
                break
    return min(1.0, matched / max(len(l_mz_arr), 1))


def _intensity_score(query_clean, lib_clean):
    if not query_clean or not lib_clean:
        return 0.0
    q_top3 = [p[1] for p in sorted(query_clean, key=lambda x: x[1], reverse=True)[:3]]
    l_top3 = [p[1] for p in sorted(lib_clean, key=lambda x: x[1], reverse=True)[:3]]
    if not l_top3:
        return 0.0
    q_max = max(q_top3) if q_top3 else 1.0
    l_max = max(l_top3) if l_top3 else 1.0
    q_norm = np.array(q_top3) / q_max if q_max > 0 else np.array(q_top3)
    l_norm = np.array(l_top3) / l_max if l_max > 0 else np.array(l_top3)
    n = min(len(q_norm), len(l_norm))
    if n == 0:
        return 0.0
    return 1.0 - min(1.0, float(np.mean(np.abs(q_norm[:n] - l_norm[:n]))))


def _confidence(cos, final, mass, expert=0.0):
    """置信度分级(增强版,含专家评分)"""
    if (cos >= 0.9 and final >= 0.5) or expert >= 0.70:
        return "High"
    if (cos >= 0.8 and final >= 0.4) or expert >= 0.45:
        return "Medium"
    if cos >= 0.7 or expert >= 0.25:
        return "Low"
    return "Unknown"


def _smiles_to_selfies(smiles):
    if not smiles or smiles.upper() in ("N/A", "NONE", "NULL", ""):
        return ""
    try:
        import selfies as sf
        return sf.encoder(smiles)
    except Exception:
        return ""


def _legacy_search_do_not_use(query_peaks, query_pmz, records, mz_index,
                              blank_profile=None,
                              blank_mode=DEFAULT_BLANK_MODE,
                              blank_threshold=DEFAULT_BLANK_THRESHOLD,
                              blank_int_mode=DEFAULT_BLANK_INTENSITY_MODE,
                              query_rt=None, query_rt_unit=None):
    """Return single best match (or None) above cosine threshold.
    If blank_profile is not None, subtract blank from query_peaks first.
    """
    if blank_profile is not None and BLANK_SUBTRACTION_AVAILABLE:
        query_peaks = subtract_blank(
            query_peaks,
            blank_profile,
            mode=blank_mode,
            intensity_threshold_factor=blank_threshold,
            blank_intensity_mode=blank_int_mode,
            precursor_mz=query_pmz,
            retention_time=query_rt,
            retention_time_unit=query_rt_unit,
        )

    candidates = _query_candidates(query_pmz, records, mz_index)
    if not candidates:
        return None
    if len(candidates) > MAX_CANDIDATES_PER_QUERY:
        candidates = sorted(
            candidates,
            key=lambda i: abs(records[i]["precursor_mz"] - query_pmz)
        )[:MAX_CANDIDATES_PER_QUERY]

    query_clean = _clean_peaks(query_peaks, query_pmz)
    best = None
    for idx in candidates:
        lib = records[idx]
        lib_peaks = lib["peaks"]
        if not lib_peaks:
            continue
        cosine = _score_cosine(query_clean, lib_peaks, query_pmz=query_pmz)
        if cosine < COSINE_THRESHOLD:
            continue
        mass = _mass_score(query_pmz, lib["precursor_mz"])
        ppm = _mass_ppm(query_pmz, lib["precursor_mz"])
        frag = _fragment_score(query_clean, lib_peaks)
        inten = _intensity_score(query_clean, lib_peaks)
        # ── 化学专家评分 ──
        expert_score = 0.0
        expert_details = None
        if EXPERT_AVAILABLE:
            try:
                expert_result = expert_evaluate(
                    query_peaks=query_clean,
                    query_pmz=query_pmz,
                    candidate_smiles=lib["smiles"] or "",
                    candidate_formula=lib.get("formula", ""),
                    candidate_classification=lib.get("classification", "Unknown"),
                    candidate_inchikey=lib.get("inchikey", ""),
                )
                expert_score = expert_result.get("expert_total_score", 0.0)
                expert_details = expert_result
            except Exception:
                expert_score = 0.0
                expert_details = None

        final = cosine * 0.30 + mass * 0.15 + frag * 0.15 + inten * 0.10 + expert_score * 0.30
        if final < FINAL_SCORE_THRESHOLD:
            continue
        result = {
            "cosine": cosine,
            "mass_score": mass,
            "mass_ppm": ppm,
            "fragment_score": frag,
            "intensity_score": inten,
            "final_score": final,
            "expert_score": expert_score,
            "expert_details": expert_details,
            # v3.0: expert sub-scores for detailed report
            "nl_score": expert_details.get("neutral_loss_score", 0) if expert_details else 0,
            "di_score": expert_details.get("diagnostic_ion_score", 0) if expert_details else 0,
            "rc_score": expert_details.get("rule_coverage_score", 0) if expert_details else 0,
            "fg_score": expert_details.get("fg_signature_score", 0) if expert_details else 0,
            "sf_score": expert_details.get("structure_fragment_score", 0) if expert_details else 0,
            "be_score": expert_details.get("bond_energy_score", 0) if expert_details else 0,
            "fl_score": expert_details.get("feature_library_score", 0) if expert_details else 0,
            "id_score": expert_details.get("isomer_differentiation_score", 0) if expert_details else 0,
            "chemical_reasoning": expert_details.get("chemical_reasoning", "") if expert_details else "",
            "spectrum_id": lib["spectrum_id"],
            "compound_name": lib["compound_name"] or lib["spectrum_id"],
            "smiles": lib["smiles"],
            "selfies": "",
            "inchikey": lib["inchikey"],
            "formula": lib["formula"],
            "adduct": lib["adduct"],
            "precursor_mz": lib["precursor_mz"],
            "classification": lib["classification"] or "Unknown",
            "lib_peaks": lib_peaks,
            "query_peaks": query_clean,  # for detailed report
        }
        if best is None or final > best["final_score"]:
            best = result
    if best is not None and not best["selfies"]:
        best["selfies"] = _smiles_to_selfies(best["smiles"])
    # MS2Classifier DL分类预测
    if MS2CLF_AVAILABLE and _clf_wrapper is not None:
        try:
            _clf_result = _clf_wrapper.predict(query_peaks, query_pmz)
            best["dl_classification"] = _clf_result["class"]
            best["dl_confidence"] = _clf_result["confidence"]
            best["dl_top3"] = _clf_result["top3"]
        except Exception:
            pass
    return best


def search(query_peaks, query_pmz, records, mz_index, blank_profile=None,
           blank_mode=DEFAULT_BLANK_MODE, blank_threshold=DEFAULT_BLANK_THRESHOLD,
           blank_int_mode=DEFAULT_BLANK_INTENSITY_MODE,
           query_rt=None, query_rt_unit=None, query_polarity=None,
           query_adduct=None, query_adduct_reliable=None):
    """Compatibility wrapper around :mod:`matching.search`.

    Blank subtraction happens once here so the report can retain the exact
    searched peak list.  All candidate eligibility, source-tier, polarity and
    confidence rules are delegated to the canonical V6 search implementation.
    """
    searched_peaks = query_peaks
    if blank_profile is not None and BLANK_SUBTRACTION_AVAILABLE:
        searched_peaks = subtract_blank(
            query_peaks,
            blank_profile,
            mode=blank_mode,
            intensity_threshold_factor=blank_threshold,
            blank_intensity_mode=blank_int_mode,
            precursor_mz=query_pmz,
            retention_time=query_rt,
            retention_time_unit=query_rt_unit,
        )

    match = _canonical_search(
        searched_peaks,
        query_pmz,
        records,
        mz_index,
        blank_profile=None,
        cos_thresh=COSINE_THRESHOLD,
        final_thresh=FINAL_SCORE_THRESHOLD,
        mz_da=MZ_TOLERANCE_DA,
        mz_ppm=MZ_TOLERANCE_PPM,
        query_adduct=query_adduct,
        query_adduct_reliable=query_adduct_reliable,
        query_polarity=query_polarity,
        query_rt=query_rt,
        query_rt_unit=query_rt_unit,
        enforce_polarity=True,
    )
    if match is None:
        return None

    # Adapt canonical result keys to the long-established CSV writer contract.
    result = dict(match)
    result.setdefault("selfies", _smiles_to_selfies(result.get("smiles")))
    result.setdefault("spectrum_id", result.get("lib_spectrum_id", ""))
    result.setdefault("precursor_mz", result.get("lib_precursor_mz"))
    result["query_peaks"] = _canonical_clean_peaks(
        searched_peaks, query_pmz, top_n=TOP_N_PEAKS,
        rel_threshold=REL_THRESHOLD,
    )
    return result


# ─────────────────────────────────────────────────────────────────
# 文件解析
# ─────────────────────────────────────────────────────────────────
_SCAN_RE = re.compile(r"scanId=(\d+)")
_TITLE_POLARITY_RE = re.compile(r"NativeID:\"[^\"]*scanId=\d+[^\"]*\"", re.IGNORECASE)


def _extract_scan_id(title):
    if not title:
        return None
    m = _SCAN_RE.search(str(title))
    return int(m.group(1)) if m else None


def _infer_polarity_from_name(name: str) -> str | None:
    """根据文件名启发式推断极性: 'A'(阴离子) / 'P'(阳离子) / None(未知)。"""
    n = name.upper()
    # Strong negative markers
    if any(tok in n for tok in ("-A", "_A.", "-A.", "NEGATIVE", "NEG-", "NEG_", "M-H")):
        return "A"
    # Strong positive markers
    if any(tok in n for tok in ("-P", "_P.", "-P.", "POSITIVE", "POS-", "POS_", "M+H")):
        return "P"
    return None


def _metadata_polarity(metadata: dict) -> str | None:
    """Normalize common matchms/vendor polarity fields to POS/NEG."""
    for key in ("ionmode", "ion_mode", "polarity", "scan_polarity"):
        value = metadata.get(key)
        if value is None:
            continue
        text = str(value).strip().upper()
        if text.startswith(("NEG", "-")):
            return "NEG"
        if text.startswith(("POS", "+")):
            return "POS"
    charge = metadata.get("charge")
    if charge is not None:
        text = str(charge).strip()
        if "-" in text:
            return "NEG"
        if "+" in text:
            return "POS"
        try:
            numeric = float(text)
            if numeric < 0:
                return "NEG"
            if numeric > 0:
                return "POS"
        except (TypeError, ValueError):
            pass
    return None


def _metadata_adduct(metadata: dict) -> str | None:
    for key in ("adduct", "precursor_type", "precursortype"):
        value = metadata.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def parse_mzml_file(file_path: Path):
    """Parse mzML → list of spectrum dicts."""
    records = []
    spectra = load_from_mzml(str(file_path))
    local_idx = 0
    for s in spectra:
        pmz = s.get("precursor_mz")
        if pmz is None:
            continue
        try:
            pmz = float(pmz)
        except (TypeError, ValueError):
            continue
        if pmz <= 0:
            continue
        peaks_arr = list(zip(s.peaks.mz.tolist(), s.peaks.intensities.tolist()))
        if not peaks_arr:
            continue
        md = dict(s.metadata)
        rt = md.get("retention_time")
        if rt is not None:
            try:
                rt = float(rt)
            except (TypeError, ValueError):
                rt = None
        scan_start = md.get("scan_start_time")
        if scan_start is not None and rt is None:
            try:
                rt = float(scan_start)
            except (TypeError, ValueError):
                rt = None
        title = md.get("title") or ""
        scan_id = _extract_scan_id(title) or local_idx
        records.append({
            "spectrum_id": scan_id,
            "precursor_mz": pmz,
            "peaks": peaks_arr,
            "title": str(title),
            "retention_time": rt,
            "scan_number": md.get("scan_number"),
            "ion_mode": _metadata_polarity(md),
            "adduct": _metadata_adduct(md),
        })
        local_idx += 1
    return records


def parse_mgf_file(file_path: Path):
    """Parse MGF via matchms → list of spectrum dicts (same shape as parse_mzml_file)."""
    records = []
    spectra = load_from_mgf(str(file_path))
    local_idx = 0
    for s in spectra:
        pmz = s.get("precursor_mz")
        if pmz is None:
            continue
        try:
            pmz = float(pmz)
        except (TypeError, ValueError):
            continue
        if pmz <= 0:
            continue
        peaks_arr = list(zip(s.peaks.mz.tolist(), s.peaks.intensities.tolist()))
        if not peaks_arr:
            continue
        md = dict(s.metadata)
        rt = md.get("retention_time")
        if rt is None:
            for key, value in md.items():
                if key and str(key).lower() == "rtinseconds":
                    rt = value
                    break
        if rt is not None:
            try:
                rt = float(rt)
            except (TypeError, ValueError):
                rt = None
        title = md.get("title") or ""
        scan_id = _extract_scan_id(title) or local_idx
        records.append({
            "spectrum_id": scan_id,
            "precursor_mz": pmz,
            "peaks": peaks_arr,
            "title": str(title),
            "retention_time": rt,
            "scan_number": None,
            "ion_mode": _metadata_polarity(md),
            "adduct": _metadata_adduct(md),
        })
        local_idx += 1
    return records


def parse_spectrum_file(file_path: Path):
    """根据后缀分发到 mzML 或 MGF 解析器。"""
    suffix = file_path.suffix.lower()
    if suffix in (".mzml", ".xml"):
        return parse_mzml_file(file_path)
    if suffix in (".mgf",):
        return parse_mgf_file(file_path)
    raise ValueError(f"不支持的文件格式: {file_path}")


# ─────────────────────────────────────────────────────────────────
# 单文件分析
# ─────────────────────────────────────────────────────────────────
def analyze_file(spec_path, records, mz_index, output_dir, blank_profile,
                 blank_mode, blank_threshold, blank_int_mode):
    """Analyze one spectrum file. Writes _结果.csv and _报告.csv.
    Returns dict with total/matched/mz_range/avg_peaks/elapsed/paths.
    """
    t_file_start = time.time()
    logger.info("  分析: %s", spec_path.name)
    spectra = parse_spectrum_file(spec_path)
    t_load = time.time() - t_file_start
    logger.info("    加载谱图: %d (%.1f 秒)", len(spectra), t_load)
    if blank_profile is not None:
        logger.info("    空白扣除: 已启用 (profile=%d 峰, mode=%s, threshold=%.2f)",
                    len(blank_profile.mz_array), blank_mode, blank_threshold)

    matched_results = []
    t_search_start = time.time()
    n_matched = 0
    n_subtracted_peaks = 0
    n_raw_peaks = 0
    rt_unit = "seconds" if spec_path.suffix.lower() == ".mgf" else "minutes"
    file_polarity = {"A": "NEG", "P": "POS"}.get(
        _infer_polarity_from_name(spec_path.name)
    )
    for idx, spec in enumerate(spectra):
        pmz = spec["precursor_mz"]
        peaks = spec["peaks"]
        n_raw_peaks += len(peaks)
        m = search(peaks, pmz, records, mz_index,
                   blank_profile=blank_profile,
                   blank_mode=blank_mode,
                   blank_threshold=blank_threshold,
                   blank_int_mode=blank_int_mode,
                   query_rt=spec.get("retention_time"),
                   query_rt_unit=rt_unit,
                   query_polarity=spec.get("ion_mode") or file_polarity,
                   query_adduct=spec.get("adduct"))
        if m is None:
            continue
        # 估算扣除后剩余峰数(用于 log)
        if blank_profile is not None and BLANK_SUBTRACTION_AVAILABLE:
            sub = subtract_blank(peaks, blank_profile,
                                 mode=blank_mode,
                                 intensity_threshold_factor=blank_threshold,
                                 blank_intensity_mode=blank_int_mode,
                                 precursor_mz=pmz,
                                 retention_time=spec.get("retention_time"),
                                 retention_time_unit=rt_unit)
            n_subtracted_peaks += len(sub)
        n_matched += 1
        matched_results.append({
            "Spectrum": spec["spectrum_id"],
            "Precursor_mz": pmz,
            "Peaks_Count": len(peaks),
            "RT": spec.get("retention_time"),
            "Match_Name": m["compound_name"],
            "Match_SMILES": m["smiles"] or "",
            "Match_SELFIES": m["selfies"] or "",
            "Match_InChIKey": m["inchikey"] or "",
            "Match_Formula": m["formula"] or "",
            "Classification": m["classification"],
            "Adduct": m["adduct"] or "",
            "Cosine_Score": m["cosine"],
            "Final_Score": m["final_score"],
            "Expert_Score": m.get("expert_score", 0.0),
            "Mass_Score": m["mass_score"],
            "Mass_PPM": m["mass_ppm"] if m["mass_ppm"] is not None else "",
            "Fragment_Score": m["fragment_score"],
            "Intensity_Score": m["intensity_score"],
            "Lib_Spectrum_ID": m["spectrum_id"],
            "Lib_Precursor_MZ": m["precursor_mz"],
            # v3.0 expert sub-scores
            "NL_Score": m.get("nl_score", 0.0),
            "DI_Score": m.get("di_score", 0.0),
            "RC_Score": m.get("rc_score", 0.0),
            "FG_Score": m.get("fg_score", 0.0),
            "SF_Score": m.get("sf_score", 0.0),
            "BE_Score": m.get("be_score", 0.0),
            "FL_Score": m.get("fl_score", 0.0),
            "ID_Score": m.get("id_score", 0.0),
            "Chemical_Reasoning": m.get("chemical_reasoning", ""),
            # MS2Classifier DL分类
            "DL_Classification": m.get("dl_classification", "N/A"),
            "DL_Confidence": m.get("dl_confidence", 0.0),
            "DL_Top3": m.get("dl_top3", []),
            # Peaks for HTML report generation
            "Query_Peaks": m.get("query_peaks", []),
            "Lib_Peaks": m.get("lib_peaks", []),
        })
        if (idx + 1) % 200 == 0:
            elapsed = time.time() - t_search_start
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            logger.info("    进度: %d / %d (匹配 %d, %.1f 谱/秒)",
                        idx + 1, len(spectra), n_matched, rate)
    t_total = time.time() - t_file_start
    logger.info("    匹配: %d / %d (耗时 %.1f 秒)", n_matched, len(spectra), t_total)
    if blank_profile is not None and n_raw_peaks > 0:
        sub_pct = 100.0 * (n_raw_peaks - n_subtracted_peaks) / n_raw_peaks
        logger.info("    空白扣除: 原始 %d 峰 → 扣除后 %d 峰 (-%.1f%%)",
                    n_raw_peaks, n_subtracted_peaks, sub_pct)

    if spectra:
        mz_vals = [s["precursor_mz"] for s in spectra]
        mz_range = f"{min(mz_vals):.1f} - {max(mz_vals):.1f}"
        avg_peaks = sum(len(s["peaks"]) for s in spectra) / len(spectra)
    else:
        mz_range = "—"
        avg_peaks = 0.0

    matched_results.sort(key=lambda r: r["Spectrum"])

    # 报告.csv (Detailed report with expert sub-scores and peaks)
    baogao_path = output_dir / f"{spec_path.stem}_报告.csv"
    with baogao_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "Spectrum", "Precursor_mz", "Peaks_Count", "RT",
            "Compound_Name", "Classification", "Adduct",
            "Cosine_Score", "Final_Score", "Expert_Score",
            "Mass_Score", "Mass_PPM", "Fragment_Score", "Intensity_Score",
            "NL_Score", "DI_Score", "RC_Score", "FG_Score", "SF_Score", "BE_Score",
            "FL_Score", "ID_Score",
            "InChIKey", "SMILES",
            "Query_Peaks_JSON", "Lib_Peaks_JSON",
            "Chemical_Reasoning", "DL_Classification", "DL_Confidence", "DL_Top3",
        ])
        for r in matched_results:
            query_peaks_json = json.dumps(r.get("Query_Peaks", []), ensure_ascii=False)
            lib_peaks_json = json.dumps(r.get("Lib_Peaks", []), ensure_ascii=False)
            w.writerow([
                r["Spectrum"],
                f"{r['Precursor_mz']:.4f}",
                r["Peaks_Count"],
                r.get("RT", ""),
                r["Match_Name"],
                r["Classification"],
                r.get("Adduct", ""),
                f"{r['Cosine_Score']:.4f}",
                f"{r['Final_Score']:.4f}",
                f"{r.get('Expert_Score', 0.0):.4f}",
                f"{r['Mass_Score']:.4f}",
                f"{r['Mass_PPM']:.2f}" if r["Mass_PPM"] != "" else "",
                f"{r['Fragment_Score']:.4f}",
                f"{r['Intensity_Score']:.4f}",
                f"{r.get('NL_Score', 0.0):.4f}",
                f"{r.get('DI_Score', 0.0):.4f}",
                f"{r.get('RC_Score', 0.0):.4f}",
                f"{r.get('FG_Score', 0.0):.4f}",
                f"{r.get('SF_Score', 0.0):.4f}",
                f"{r.get('BE_Score', 0.0):.4f}",
                f"{r.get('FL_Score', 0.0):.4f}",
                f"{r.get('ID_Score', 0.0):.4f}",
                r["Match_InChIKey"],
                r["Match_SMILES"],
                query_peaks_json,
                lib_peaks_json,
                r.get("Chemical_Reasoning", ""),
                r.get("DL_Classification", "N/A"),
                f"{r.get("DL_Confidence", 0.0):.4f}",
                str(r.get("DL_Top3", [])),
            ])

    # 结果.csv
    jieguo_path = output_dir / f"{spec_path.stem}_结果.csv"
    with jieguo_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "Spectrum", "Precursor_mz", "Match_Name", "Match_SMILES",
            "Match_SELFIES", "Match_InChIKey", "Match_Formula",
            "Final_Score", "Cosine_Score", "Expert_Score", "Confidence", "Class", "DL_Classification", "DL_Confidence",
        ])
        for r in matched_results:
            confidence = _confidence(r["Cosine_Score"], r["Final_Score"], r["Mass_Score"],
                                    r.get("Expert_Score", 0.0))
            w.writerow([
                r["Spectrum"],
                f"{r['Precursor_mz']:.4f}",
                r["Match_Name"],
                r["Match_SMILES"],
                r["Match_SELFIES"],
                r["Match_InChIKey"],
                r["Match_Formula"],
                f"{r['Final_Score']:.4f}",
                f"{r['Cosine_Score']:.4f}",
                f"{r.get('Expert_Score', 0.0):.4f}",      # 新增
                confidence,
                r["Classification"],
                r.get("DL_Classification", "N/A"),
                f"{r.get("DL_Confidence", 0.0):.4f}",
            ])

    # ── 保存完整匹配数据为JSON (保留所有peak数据，用于HTML报告生成) ──
    match_data_path = output_dir / f"{spec_path.stem}_match_data.json"
    with match_data_path.open("w", encoding="utf-8") as f:
        json.dump(matched_results, f, ensure_ascii=False, default=str)

    return {
        "file_name": spec_path.name,
        "total_spectra": len(spectra),
        "matched": n_matched,
        "match_rate": n_matched / len(spectra) * 100 if spectra else 0.0,
        "mz_range": mz_range,
        "avg_peaks": avg_peaks,
        "elapsed": t_total,
        "blank_used": blank_profile is not None,
        "result_path": str(jieguo_path),
        "report_path": str(baogao_path),
    }


# ─────────────────────────────────────────────────────────────────
# 空白 profile 加载
# ─────────────────────────────────────────────────────────────────
def load_blank_profiles(args) -> dict[str, BlankProfile]:
    """Load blank profiles from CLI args. Returns dict keyed by polarity ("A"/"P")."""
    profiles: dict[str, BlankProfile] = {}
    if not BLANK_SUBTRACTION_AVAILABLE:
        logger.warning("blank_subtraction 模块不可用,空白扣除已禁用")
        return profiles
    if args.no_blank:
        logger.info("空白扣除: --no-blank 已禁用")
        return profiles

    candidates = [("A", args.blank_a), ("P", args.blank_p)]
    for polarity, path in candidates:
        if path is None:
            continue
        p = Path(path)
        if not p.exists():
            logger.warning("空白文件不存在,跳过 [%s]: %s", polarity, p)
            continue
        try:
            logger.info("加载空白 profile [%s]: %s", polarity, p)
            prof = load_blank_mgf(p, mz_tolerance_da=BLANK_MZ_TOL_DA,
                                  mz_tolerance_ppm=BLANK_MZ_TOL_PPM)
            profiles[polarity] = prof
            logger.info("  空白 [%s] 峰数量: %d", polarity, len(prof.mz_array))
        except Exception as e:
            logger.warning("加载空白 [%s] 失败,跳过: %s", polarity, e)
    return profiles


def select_blank_profile(profiles: dict[str, BlankProfile], spec_path: Path):
    """根据文件名启发式选择匹配的 blank profile,未匹配返回 None。"""
    if not profiles:
        return None
    pol = _infer_polarity_from_name(spec_path.name)
    if pol is None:
        logger.warning("    无法从文件名推断极性,跳过空白扣除: %s", spec_path.name)
        return None
    if pol not in profiles:
        logger.warning("    推断极性 %s 但无对应空白 profile,跳过", pol)
        return None
    return profiles[pol]


# ─────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────
def build_arg_parser():
    p = argparse.ArgumentParser(
        description="mzML/MGF 批处理质谱库匹配,支持空白扣除。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--samples", type=Path, default=DEFAULT_SAMPLES_DIR,
                   help="样本文件夹路径")
    p.add_argument("--output", type=Path, default=None,
                   help="输出目录(默认 <samples>_results)")
    p.add_argument("--db", type=Path, default=DEFAULT_DB_PATH,
                   help="质谱库 SQLite 路径")
    p.add_argument("--blank-a", type=Path, default=DEFAULT_BLANK_A,
                   help="阴离子空白 MGF 路径(存在则启用)")
    p.add_argument("--blank-p", type=Path, default=DEFAULT_BLANK_P,
                   help="阳离子空白 MGF 路径(存在则启用)")
    p.add_argument("--no-blank", action="store_true",
                   help="禁用空白扣除")
    p.add_argument("--blank-mode",
                   choices=("absolute", "relative", "filter_only"),
                   default=DEFAULT_BLANK_MODE,
                   help="空白扣除模式")
    p.add_argument("--blank-threshold", type=float, default=DEFAULT_BLANK_THRESHOLD,
                   help="空白强度阈值因子(样品峰 < 空白*threshold 时被过滤)")
    p.add_argument("--blank-int-mode", choices=("max", "mean"),
                   default=DEFAULT_BLANK_INTENSITY_MODE,
                   help="空白强度使用 max 或 mean")
    p.add_argument("--include", nargs="*", default=None,
                   help="只处理文件名包含这些子串的文件(空格分隔)")
    p.add_argument("--exclude", nargs="*", default=None,
                   help="跳过文件名包含这些子串的文件(空格分隔)")
    p.add_argument("--dry-run", action="store_true",
                   help="只加载空白 profile,不做匹配(用于快速验证)")
    return p


def _load_results_for_html(result_csv_path: Path, report_csv_path: Path, sample_path: Path) -> list:
    """Load matched results from JSON match_data file (preserves full peak data for structure images).
    Falls back to CSV-based loading if JSON not found.
    """
    import csv as _csv
    import json as _json

    match_json_path = result_csv_path.parent / f"{sample_path.stem}_match_data.json"

    # ── Primary: load from match_data JSON (full data preservation) ──
    if match_json_path.exists():
        with match_json_path.open("r", encoding="utf-8") as f:
            raw_rows = _json.load(f)
        results = []
        for r in raw_rows:
            pmz = float(r.get("Precursor_mz", 0))
            query_peaks = r.get("Query_Peaks", [])
            lib_peaks = r.get("Lib_Peaks", [])
            results.append({
                "precursor_mz": pmz,
                "spectrum_id": r.get("Spectrum", ""),
                "peaks": query_peaks,
                "query_peaks": query_peaks,  # key used by _gen_rich
                "peaks_count": len(query_peaks),
                "retention_time": r.get("RT", ""),
                "matches": [{
                    "compound_name": r.get("Match_Name", ""),
                    "smiles": r.get("Match_SMILES", ""),
                    "inchikey": r.get("Match_InChIKey", ""),
                    "formula": r.get("Match_Formula", ""),
                    "cosine": float(r.get("Cosine_Score", 0) or 0),
                    "final_score": float(r.get("Final_Score", 0) or 0),
                    "mass_score": float(r.get("Mass_Score", 0) or 0),
                    "mass_ppm": float(r.get("Mass_PPM", 0) or 0) if r.get("Mass_PPM", "") else None,
                    "fragment_score": float(r.get("Fragment_Score", 0) or 0),
                    "intensity_score": float(r.get("Intensity_Score", 0) or 0),
                    "expert_score": float(r.get("Expert_Score", 0) or 0),
                    "classification": r.get("Classification", ""),
                    "adduct": r.get("Adduct", ""),
                    "lib_peaks": lib_peaks,
                    # ── v3.0 expert sub-scores ──
                    "neutral_loss_score": float(r.get("NL_Score", 0) or 0),
                    "diagnostic_ion_score": float(r.get("DI_Score", 0) or 0),
                    "rule_coverage_score": float(r.get("RC_Score", 0) or 0),
                    "fg_signature_score": float(r.get("FG_Score", 0) or 0),
                    "structure_fragment_score": float(r.get("SF_Score", 0) or 0),
                    "bond_energy_score": float(r.get("BE_Score", 0) or 0),
                    "feature_library_score": float(r.get("FL_Score", 0) or 0),
                    "isomer_differentiation_score": float(r.get("ID_Score", 0) or 0),
                    "chemical_reasoning": r.get("Chemical_Reasoning", ""),
                }],
            })
        return results

    # ── Fallback: load from CSV ──
    report_by_spec = {}
    with report_csv_path.open("r", encoding="utf-8") as f:
        reader = _csv.DictReader(f)
        for row in reader:
            sid = row.get("Spectrum", "")
            if sid:
                report_by_spec[sid] = row

    results = []
    with result_csv_path.open("r", encoding="utf-8") as f:
        reader = _csv.DictReader(f)
        for row in reader:
            try:
                pmz = float(row.get("Precursor_mz", 0))
            except (ValueError, TypeError):
                continue
            sid = row.get("Spectrum", "")
            report_row = report_by_spec.get(sid, {})
            results.append({
                "precursor_mz": pmz,
                "spectrum_id": sid,
                "peaks": [],
                "query_peaks": [],
                "peaks_count": 0,
                "retention_time": report_row.get("RT", ""),
                "matches": [{
                    "compound_name": row.get("Match_Name", ""),
                    "smiles": row.get("Match_SMILES", ""),
                    "inchikey": row.get("Match_InChIKey", ""),
                    "formula": row.get("Match_Formula", ""),
                    "cosine": float(row.get("Cosine_Score", 0) or 0),
                    "final_score": float(row.get("Final_Score", 0) or 0),
                    "classification": row.get("Class", ""),
                    "expert_score": float(row.get("Expert_Score", 0) or 0),
                    "lib_peaks": [],
                    "neutral_loss_score": float(report_row.get("NL_Score", 0) or 0),
                    "diagnostic_ion_score": float(report_row.get("DI_Score", 0) or 0),
                    "rule_coverage_score": float(report_row.get("RC_Score", 0) or 0),
                    "fg_signature_score": float(report_row.get("FG_Score", 0) or 0),
                    "structure_fragment_score": float(report_row.get("SF_Score", 0) or 0),
                    "bond_energy_score": float(report_row.get("BE_Score", 0) or 0),
                    "feature_library_score": float(report_row.get("FL_Score", 0) or 0),
                    "isomer_differentiation_score": float(report_row.get("ID_Score", 0) or 0),
                    "chemical_reasoning": report_row.get("Chemical_Reasoning", ""),
                }],
            })
    return results


def _feed_expert_logger(elog, result_csv_path: Path, report_csv_path: Path, sample_name: str):
    """Feed matched results from CSVs into ExpertLogger for proper logging."""
    import csv as _csv
    result_rows = []
    with result_csv_path.open("r", encoding="utf-8") as f:
        reader = _csv.DictReader(f)
        result_rows = list(reader)

    report_rows = {}
    with report_csv_path.open("r", encoding="utf-8") as f:
        reader = _csv.DictReader(f)
        for row in reader:
            report_rows[row.get("Spectrum", "")] = row

    for row in result_rows:
        try:
            pmz = float(row.get("Precursor_mz", 0))
            expert_score = float(row.get("Expert_Score", 0) or 0)
            cosine = float(row.get("Cosine_Score", 0) or 0)
            final = float(row.get("Final_Score", 0) or 0)
        except (ValueError, TypeError):
            continue

        spec_id = row.get("Spectrum", "")
        report_row = report_rows.get(spec_id, {})

        query_info = {
            "spectrum_id": spec_id,
            "precursor_mz": pmz,
            "sample": sample_name,
        }
        match_result = {
            "compound_name": row.get("Match_Name", ""),
            "smiles": row.get("Match_SMILES", ""),
            "inchikey": row.get("Match_InChIKey", ""),
            "formula": row.get("Match_Formula", ""),
            "classification": row.get("Class", ""),
            "cosine": cosine,
            "final_score": final,
            "expert_score": expert_score,
        }
        expert_scores = {
            "expert_total_score": expert_score,
            "neutral_loss_score": float(report_row.get("NeutralLoss_Score", 0) or 0),
            "diagnostic_ion_score": float(report_row.get("DiagIon_Score", 0) or 0),
            "rule_coverage_score": float(report_row.get("RuleCoverage_Score", 0) or 0),
            "confidence_level": (
                "High" if expert_score >= 0.70 else
                "Medium" if expert_score >= 0.45 else
                "Low" if expert_score >= 0.25 else
                "Insufficient"
            ),
        }
        try:
            elog.log_evaluation(query_info, match_result, expert_scores)
        except Exception:
            pass


def main():
    args = build_arg_parser().parse_args()

    samples_dir: Path = args.samples
    output_dir: Path = args.output or samples_dir.parent / f"{samples_dir.name}_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("样本目录: %s", samples_dir)
    logger.info("输出目录: %s", output_dir)

    # Load blank profiles
    global _blank_profiles
    _blank_profiles = load_blank_profiles(args)
    if _blank_profiles:
        logger.info("=" * 60)
        logger.info("空白溶剂峰扣除: 已启用")
        logger.info("  模式: %s", args.blank_mode)
        logger.info("  强度阈值: %.1f%% (样品峰 < 空白*threshold 时过滤)",
                     args.blank_threshold * 100)
        logger.info("  空白强度: %s", args.blank_int_mode)
        logger.info("  已加载 profiles: %s", list(_blank_profiles.keys()))
        logger.info("=" * 60)

    if args.dry_run:
        logger.info("--dry-run: 已加载 %d 个空白 profile,退出", len(_blank_profiles))
        return

    # Load library
    logger.info("加载质谱库...")
    records, mz_index = load_library(args.db)

    # Discover sample files (mzML + MGF)
    all_files: list[Path] = []
    for ext in ("*.mzML", "*.mzml", "*.MGF", "*.mgf"):
        all_files.extend(samples_dir.glob(ext))
    # Deduplicate
    seen = set()
    files = []
    for f in sorted(all_files):
        if f in seen:
            continue
        seen.add(f)
        files.append(f)

    if args.include:
        files = [f for f in files if any(s in f.name for s in args.include)]
    if args.exclude:
        files = [f for f in files if not any(s in f.name for s in args.exclude)]

    if not files:
        logger.error("在 %s 中未找到 mzML/MGF 文件", samples_dir)
        return
    logger.info("找到 %d 个样本文件", len(files))
    for f in files:
        pol = _infer_polarity_from_name(f.name)
        pol_str = f" [极性={pol}]" if pol else ""
        logger.info("  - %s (%.1f MB)%s", f.name, f.stat().st_size / 1e6, pol_str)

    summary_rows = []
    t0 = time.time()

    for fp in files:
        blank_profile = select_blank_profile(_blank_profiles, fp)
        info = analyze_file(
            fp, records, mz_index, output_dir,
            blank_profile=blank_profile,
            blank_mode=args.blank_mode,
            blank_threshold=args.blank_threshold,
            blank_int_mode=args.blank_int_mode,
        )
        summary_rows.append(info)
        logger.info("    → %d / %d 匹配 (%.1f%%)",
                     info["matched"], info["total_spectra"], info["match_rate"])

    elapsed = time.time() - t0
    logger.info("")
    logger.info("=" * 60)
    logger.info("分析完成！")
    total_spec = sum(r["total_spectra"] for r in summary_rows)
    total_matched = sum(r["matched"] for r in summary_rows)
    logger.info("  处理文件: %d", len(files))
    logger.info("  谱图总数: %d", total_spec)
    logger.info("  成功匹配 (cos>=%.2f): %d (%.1f%%)",
                COSINE_THRESHOLD, total_matched,
                total_matched / total_spec * 100 if total_spec else 0)
    logger.info("  输出目录: %s", output_dir)
    logger.info("  总耗时: %.1f 秒", elapsed)
    logger.info("=" * 60)

    # 总体摘要.csv
    sum_path = output_dir / "总体摘要.csv"
    with sum_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "文件名", "总谱图数", "匹配数", "匹配率(%)", "m/z范围",
            "平均峰数", "处理耗时(秒)", "空白扣除", "结果文件", "报告文件",
        ])
        for r in summary_rows:
            w.writerow([
                r["file_name"],
                r["total_spectra"],
                r["matched"],
                f"{r['match_rate']:.1f}",
                r["mz_range"],
                f"{r['avg_peaks']:.1f}",
                f"{r['elapsed']:.1f}",
                "是" if r["blank_used"] else "否",
                r["result_path"],
                r["report_path"],
            ])
    logger.info("  汇总: %s", sum_path)

    # ── 生成 Standard HTML 报告 (匹配网页版 Standard Format) ──
    try:
        from web.report import generate_html
        for fp in files:
            result_path = output_dir / f"{fp.stem}_结果.csv"
            report_csv_path = output_dir / f"{fp.stem}_报告.csv"
            html_path = output_dir / f"{fp.stem}_standard_report.html"
            if result_path.exists() and report_csv_path.exists():
                results_for_html = _load_results_for_html(result_path, report_csv_path, fp)
                generate_html(results_for_html, str(html_path), fp.stem)
                logger.info("  Standard报告: %s", html_path)
    except Exception as e:
        logger.warning("  Standard HTML报告生成失败: %s", e)

    # ── 生成 Detailed HTML 报告 (匹配网页版 Detailed Format) ──
    try:
        from web.app import _gen_rich
        for fp in files:
            result_path = output_dir / f"{fp.stem}_结果.csv"
            report_csv_path = output_dir / f"{fp.stem}_报告.csv"
            rich_html_path = output_dir / f"{fp.stem}_detailed_report.html"
            if result_path.exists() and report_csv_path.exists():
                results_for_rich = _load_results_for_html(result_path, report_csv_path, fp)
                rich_html = _gen_rich(results_for_rich, source_name=fp.name)
                rich_html_path.write_text(rich_html, encoding="utf-8")
                logger.info("  Detailed报告: %s", rich_html_path)
    except Exception as e:
        logger.warning("  Detailed HTML报告生成失败: %s", e)

    # ── 生成 Expert Analysis HTML 报告 (v3.0 新增,不替代现有报告) ──
    try:
        from web.expert_report import generate_expert_report
        for fp in files:
            result_path = output_dir / f"{fp.stem}_结果.csv"
            report_csv_path = output_dir / f"{fp.stem}_报告.csv"
            expert_html_path = output_dir / f"{fp.stem}_expert_analysis_report.html"
            if result_path.exists() and report_csv_path.exists():
                results_for_expert = _load_results_for_html(result_path, report_csv_path, fp)
                generate_expert_report(results_for_expert, str(expert_html_path), fp.stem)
                logger.info("  Expert分析报告: %s", expert_html_path)
    except Exception as e:
        logger.warning("  Expert分析报告生成失败: %s", e)

    # ── 生成专家评估报告 ──
    if EXPERT_AVAILABLE:
        try:
            elog = ExpertLogger()
            # Feed matched results into expert logger
            for fp in files:
                result_path = output_dir / f"{fp.stem}_结果.csv"
                report_path = output_dir / f"{fp.stem}_报告.csv"
                if result_path.exists() and report_path.exists():
                    _feed_expert_logger(elog, result_path, report_path, fp.stem)
            summary = elog.write_summary()
            elog.generate_report()
            logger.info("  专家评估报告: logs/expert_report.md")
            logger.info("  专家评估汇总: logs/expert_summary.json")
            if summary.get("count", 0) > 0:
                logger.info("  评估数量: %d, 平均专家评分: %.3f",
                           summary.get("count", 0),
                           summary.get("expert_scores", {}).get("mean", 0))
        except Exception as e:
            logger.warning("  专家报告生成失败: %s", e)


if __name__ == "__main__":
    main()
