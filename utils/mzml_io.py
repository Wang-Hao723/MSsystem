"""mzML 文件解析"""
import logging
import re
from bisect import bisect_left
from pathlib import Path
from matchms.importing import load_from_mzml

from utils.adduct import NEUTRAL_SHIFTS, adduct_polarity, neutral_mass
from utils.instrument_metadata import extract_mzml_instrument_metadata
from utils.spectrum_contract import is_auxiliary_detector_spectrum

logger = logging.getLogger(__name__)

_SCAN_RE = re.compile(r"scanId=(\d+)")
_SPECIFIC_ADDUCTS = {re.sub(r"\s+", "", key).upper() for key in NEUTRAL_SHIFTS}


def _metadata_value(metadata, *aliases):
    wanted = {re.sub(r"[^a-z0-9]", "", alias.lower()) for alias in aliases}
    for key, value in metadata.items():
        normalized = re.sub(r"[^a-z0-9]", "", str(key or "").lower())
        if normalized in wanted and value is not None and str(value).strip():
            return value
    return None


def _file_polarity(file_path):
    """mzML 文件级极性兜底 (按 cvParam 计数)。

    matchms.load_from_mzml 实测不把 mzML 的极性 cvParam 映射进 metadata.ionmode
    (MS:1000129=negative scan / MS:1000130=positive scan 被丢弃) → spec.polarity 全 None,
    使下游 ICEBERG 兜底 neutral_mass 与检索加合物多假设/诊断峰抢救失效。
    按两 accession 计数取多数: 单极性文件全覆盖, 混合取多数, 无信息返 None 不臆测。
    """
    try:
        neg = pos = 0
        overlap = 16  # 防 accession 串跨 chunk 边界漏计
        prev = b""
        with open(file_path, "rb") as f:
            while True:
                chunk = f.read(1024 * 1024)
                if not chunk:
                    break
                data = prev + chunk
                neg += data.count(b"MS:1000129")
                pos += data.count(b"MS:1000130")
                prev = chunk[-overlap:] if len(chunk) >= overlap else chunk
    except Exception:
        return None
    if neg and pos:
        return "MIXED"
    if neg == pos:
        return None
    return "NEG" if neg > pos else "POS"


def _file_is_centroid(file_path) -> bool | None:
    """mzML 文件级"谱图表示"判定 (Fix1, 按 cvParam 计数)。

    真实仪器导出的 mzML 在每条 spectrum 上标注谱图表示 cvParam:
      MS:1000127 = centroid spectrum (已质心化, 只有峰顶点)
      MS:1000126 = profile spectrum (连续轮廓采样)
    matchms 不暴露该标记, 旧代码靠 centroiding._is_profile_like 的 gap 占比启发式判断 →
    富含碎片的真实质心谱图(>15% 相邻 gap < 0.05)被误判 profile 而错误合并, 污染 cosine
    (33 标品对照里 Magnoflorine 等已知物余弦反常偏低的根因)。

    本函数按两 accession 计数: centroid 多 → True(应跳过 smart_peak_pick); profile 多 →
    False(需质心化); 两者都无/相等 → None(未知, 退回启发式)。
    """
    try:
        cen = prof = 0
        overlap = 16
        prev = b""
        with open(file_path, "rb") as f:
            while True:
                chunk = f.read(1024 * 1024)
                if not chunk:
                    break
                data = prev + chunk
                cen += data.count(b"MS:1000127")
                prof += data.count(b"MS:1000126")
                prev = chunk[-overlap:] if len(chunk) >= overlap else chunk
    except Exception:
        return None
    if cen == prof == 0:
        return None
    if cen == prof:
        return None
    return cen > prof


def _safe_float(value):
    try:
        number = float(value)
        return None if number != number else number
    except (TypeError, ValueError):
        return None


def _retention_time_minutes(raw, *, default_unit="minute"):
    """Normalize a Pyteomics/matchms retention-time value to minutes.

    Pyteomics preserves the mzML CV unit on ``unitfloat.unit_info``.  OpenMS
    serializes scan times in seconds, while many vendor mzML files use minutes;
    coercing either value directly to ``float`` silently loses that distinction.
    Plain numeric values retain the historical mzML default of minutes.
    """
    if raw is None:
        return None
    if isinstance(raw, (list, tuple)):
        if not raw:
            return None
        raw = raw[0]
    value = _safe_float(raw)
    if value is None:
        return None
    unit_text = " ".join((
        str(getattr(raw, "unit_info", "") or ""),
        str(getattr(raw, "unit_name", "") or ""),
    )).strip().lower()
    if not unit_text:
        unit_text = str(default_unit or "minute").strip().lower()
    if "millisecond" in unit_text:
        return value / 60000.0
    if "second" in unit_text and "minute" not in unit_text:
        return value / 60.0
    if "hour" in unit_text:
        return value * 60.0
    return value


def _scan_start_time_minutes(spectrum):
    scans = ((spectrum.get("scanList") or {}).get("scan") or [])
    if not scans:
        return None
    return _retention_time_minutes(scans[0].get("scan start time"))


def _nearest_ms1_peak_intensity(peaks, target_mz, tolerance_da=0.01,
                                tolerance_ppm=10.0):
    """Return the strongest genuine MS1 peak around a selected precursor.

    Selected-ion ``peak intensity`` is optional in mzML and was absent in
    several vendor exports used by this project.  The preceding referenced MS1
    scan is the correct place to recover the precursor abundance used for
    feature-level blank filtering; product-ion TIC is not a substitute.
    """
    if not peaks or target_mz is None:
        return None
    try:
        target = float(target_mz)
    except (TypeError, ValueError):
        return None
    window = min(
        float(tolerance_da),
        abs(target) * float(tolerance_ppm) * 1e-6,
    )
    mzs = [float(item[0]) for item in peaks]
    pos = bisect_left(mzs, target - window)
    best = None
    for idx in range(pos, len(peaks)):
        mz, intensity = peaks[idx]
        if float(mz) > target + window:
            break
        value = _safe_float(intensity)
        if value is not None and value > 0 and (best is None or value > best):
            best = value
    return best


def evaluate_precursor_isolation_purity(
    peaks,
    target_mz,
    lower_offset,
    upper_offset,
    *,
    charge=1,
    tolerance_da=0.01,
    tolerance_ppm=10.0,
):
    """Quantify co-isolation from the genuine precursor MS1 scan.

    The target signal includes M, M+1 and M+2 of the selected charge state so
    that a natural isotope envelope is not mislabeled as contamination.  This
    is audit evidence only: it does not change search ranking or confidence.
    """
    target = _safe_float(target_mz)
    lower = _safe_float(lower_offset)
    upper = _safe_float(upper_offset)
    if not peaks or target is None or lower is None or upper is None:
        return {
            "status": "isolation_window_unavailable",
            "purity": None,
            "target_cluster_intensity": None,
            "isolation_window_intensity": None,
            "coisolated_peak_count": None,
        }
    if lower <= 0 or upper <= 0:
        return {
            "status": "isolation_window_invalid",
            "purity": None,
            "target_cluster_intensity": None,
            "isolation_window_intensity": None,
            "coisolated_peak_count": None,
        }
    try:
        z = max(abs(int(charge or 1)), 1)
    except (TypeError, ValueError):
        z = 1
    start = target - lower
    end = target + upper
    window_peaks = [
        (float(mz), float(intensity))
        for mz, intensity in peaks
        if start <= float(mz) <= end
        and (_safe_float(intensity) or 0.0) > 0
    ]
    total = sum(item[1] for item in window_peaks)
    if total <= 0:
        return {
            "status": "no_ms1_signal_in_isolation_window",
            "purity": None,
            "target_cluster_intensity": 0.0,
            "isolation_window_intensity": 0.0,
            "coisolated_peak_count": 0,
            "window_start_mz": round(start, 6),
            "window_end_mz": round(end, 6),
        }
    isotope_targets = [target + isotope * 1.00335483507 / z for isotope in range(3)]

    def belongs_to_target_cluster(mz):
        return any(
            abs(mz - isotope_mz) <= min(
                float(tolerance_da),
                abs(isotope_mz) * float(tolerance_ppm) * 1e-6,
            )
            for isotope_mz in isotope_targets
            if start <= isotope_mz <= end
        )

    target_signal = sum(
        intensity for mz, intensity in window_peaks
        if belongs_to_target_cluster(mz)
    )
    coisolated = [
        (mz, intensity) for mz, intensity in window_peaks
        if not belongs_to_target_cluster(mz)
    ]
    purity = min(max(target_signal / total, 0.0), 1.0)
    status = "clean" if purity >= 0.80 else "mixed" if purity >= 0.50 else "severe"
    strongest = sorted(coisolated, key=lambda item: item[1], reverse=True)[:5]
    return {
        "status": status,
        "purity": round(purity, 6),
        "target_cluster_intensity": round(target_signal, 6),
        "isolation_window_intensity": round(total, 6),
        "coisolated_peak_count": len(coisolated),
        "strongest_coisolated_peaks": [
            [round(mz, 6), round(intensity, 6)] for mz, intensity in strongest
        ],
        "window_start_mz": round(start, 6),
        "window_end_mz": round(end, 6),
        "charge_used": z,
        "interpretation": (
            "Precursor isolation is dominated by the target isotope cluster."
            if status == "clean" else
            "Mixed precursor isolation can distort fragment and isotope evidence."
            if status == "mixed" else
            "Severe co-isolation: the MS/MS spectrum is likely chimeric."
        ),
    }


def _collision_energy_unit(activation) -> str:
    """Recover an explicit mzML CE unit without inventing one."""
    text = str(activation or "").lower()
    if "normalized collision energy" in text or "nce" in text or "percent" in text:
        return "NCE"
    if "electronvolt" in text or "electron volt" in text or " ev" in text:
        return "eV"
    try:
        raw = activation.get("collision energy")
        unit_info = str(getattr(raw, "unit_info", "") or "").lower()
        if "electronvolt" in unit_info:
            return "eV"
        if "percent" in unit_info:
            return "NCE"
    except Exception:
        pass
    return ""


def _activation_method(activation) -> str:
    """Preserve the explicit mzML activation label for downstream auditing."""
    if not activation:
        return ""
    text = " ".join(
        str(key)
        for key in getattr(activation, "keys", lambda: [])()
        if str(key).strip().lower() != "collision energy"
        and (
            "dissociation" in str(key).lower()
            or str(key).strip().lower() in {"cid", "hcd"}
        )
    ).strip()
    if text:
        return text
    raw = str(activation or "").strip()
    return raw if raw not in {"{}", "None"} else ""


def _collision_energy_mode(activation, collision_energy) -> str:
    text = str(activation or "").lower()
    if "stepped" in text or "step collision" in text:
        return "stepped"
    if "ramp" in text:
        return "ramp"
    return "fixed" if collision_energy is not None else "missing"


def _extract_mzml_precursor_evidence(
    file_path,
    isolation_window_width_da=None,
    isolation_window_width_source=None,
    include_ordinal_index=False,
) -> dict:
    """Read precursor metadata and genuine MS1 isotope evidence by scan ID.

    matchms intentionally exposes a compact metadata subset and, for the
    current Bruker exports, drops selected-ion charge/intensity and collision
    energy.  Pyteomics preserves those fields.  If an MS1 scan is present in
    the mzML, its precursor envelope is attached; MS2-only exports are marked
    unavailable instead of incorrectly treating product ions as isotopes.
    """
    try:
        from pyteomics import mzml as _pmzml
        from utils.isotope import detect_precursor_isotope_envelope
    except Exception:
        return {}

    override_width = _safe_float(isolation_window_width_da)
    if override_width is not None and override_width <= 0:
        raise ValueError("isolation_window_width_da must be positive")
    out = {}
    last_ms1 = None
    last_ms1_id = None
    envelope_cache = {}
    ms2_ordinal = 0
    try:
        with _pmzml.read(str(file_path)) as reader:
            for index, spectrum in enumerate(reader):
                level = spectrum.get("ms level")
                native_id = str(spectrum.get("id") or "")
                scan_id = _extract_scan_id(native_id)
                if scan_id is None:
                    match = re.search(r"scan=(\d+)", native_id)
                    scan_id = int(match.group(1)) if match else index

                if level == 1:
                    mzs = spectrum.get("m/z array")
                    intensities = spectrum.get("intensity array")
                    if mzs is not None and intensities is not None:
                        last_ms1 = list(zip(mzs.tolist(), intensities.tolist()))
                        last_ms1_id = native_id
                        envelope_cache = {}
                    continue
                if level != 2:
                    continue

                if "centroid spectrum" in spectrum:
                    spectrum_representation = "centroid"
                    spectrum_representation_source = "mzml_cv"
                elif "profile spectrum" in spectrum:
                    spectrum_representation = "profile"
                    spectrum_representation_source = "mzml_cv"
                else:
                    spectrum_representation = "unknown"
                    spectrum_representation_source = "missing_cv"

                precursors = ((spectrum.get("precursorList") or {}).get("precursor") or [])
                precursor = precursors[0] if precursors else {}
                selected = (((precursor.get("selectedIonList") or {}).get("selectedIon") or [{}])[0])
                pmz = _safe_float(selected.get("selected ion m/z"))
                charge_raw = selected.get("charge state")
                try:
                    charge = abs(int(charge_raw)) if charge_raw is not None else None
                except (TypeError, ValueError):
                    charge = None
                activation = precursor.get("activation") or {}
                isolation = precursor.get("isolationWindow") or {}
                collision_energy = _safe_float(activation.get("collision energy"))
                collision_energy_unit = _collision_energy_unit(activation)
                activation_method = _activation_method(activation)
                collision_energy_mode = _collision_energy_mode(
                    activation, collision_energy
                )
                ref = str(precursor.get("spectrumRef") or "")

                envelope = None
                isolation_evidence = {
                    "status": "referenced_ms1_not_in_file",
                    "purity": None,
                }
                ms1_precursor_intensity = None
                isotope_status = "referenced_ms1_not_in_file"
                ref_matches_last = bool(
                    last_ms1 is not None
                    and (not ref or ref == last_ms1_id)
                )
                if ref_matches_last and pmz is not None:
                    ms1_precursor_intensity = _nearest_ms1_peak_intensity(
                        last_ms1, pmz,
                    )
                    cache_key = (round(pmz, 6), charge or 1)
                    if cache_key not in envelope_cache:
                        envelope_cache[cache_key] = detect_precursor_isotope_envelope(
                            last_ms1, pmz, charge=charge or 1
                        )
                    envelope = envelope_cache[cache_key]
                    isotope_status = (
                        "observed_from_ms1" if envelope else "precursor_not_observed_in_ms1"
                    )
                    isolation_target = (
                        _safe_float(isolation.get("isolation window target m/z"))
                        or pmz
                    )
                    lower_offset = _safe_float(
                        isolation.get("isolation window lower offset")
                    )
                    upper_offset = _safe_float(
                        isolation.get("isolation window upper offset")
                    )
                    isolation_source = "mzml_metadata"
                    if override_width is not None and (
                        lower_offset is None or upper_offset is None
                    ):
                        lower_offset = override_width / 2.0
                        upper_offset = override_width / 2.0
                        isolation_source = (
                            str(isolation_window_width_source or "").strip()
                            or "user_full_width_override"
                        )
                    isolation_evidence = evaluate_precursor_isolation_purity(
                        last_ms1,
                        isolation_target,
                        lower_offset,
                        upper_offset,
                        charge=charge or 1,
                    )
                    isolation_evidence["isolation_window_source"] = isolation_source

                record = {
                    "native_scan_id": scan_id,
                    "native_id": native_id,
                    "retention_time": _scan_start_time_minutes(spectrum),
                    "retention_time_unit": "minutes",
                    "precursor_charge": charge,
                    "precursor_intensity": _safe_float(selected.get("peak intensity")),
                    "precursor_ms1_intensity": ms1_precursor_intensity,
                    "collision_energy": collision_energy,
                    "collision_energy_unit": collision_energy_unit,
                    "collision_energy_mode": collision_energy_mode,
                    "activation_method": activation_method,
                    "activation_raw": str(activation or ""),
                    "precursor_scan_ref": ref,
                    "isolation_window_target_mz": _safe_float(
                        isolation.get("isolation window target m/z")
                    ),
                    "isolation_window_lower_offset": (
                        lower_offset if ref_matches_last and pmz is not None
                        else _safe_float(isolation.get("isolation window lower offset"))
                    ),
                    "isolation_window_upper_offset": (
                        upper_offset if ref_matches_last and pmz is not None
                        else _safe_float(isolation.get("isolation window upper offset"))
                    ),
                    "precursor_isotope_evidence": envelope,
                    "isotope_evidence_status": isotope_status,
                    "precursor_isolation_evidence": isolation_evidence,
                    "spectrum_representation": spectrum_representation,
                    "spectrum_representation_source": (
                        spectrum_representation_source
                    ),
                }
                out[scan_id] = record
                if include_ordinal_index:
                    # Some ProteoWizard releases omit the optional spectrum-title
                    # cvParam even though the mzML native ID is present. matchms
                    # then exposes title=None, so retain an ordered join key to
                    # keep the MS2 record linked to its genuine preceding MS1.
                    out[("ms2_ordinal", ms2_ordinal)] = record
                ms2_ordinal += 1
    except Exception as exc:
        logger.warning("mzML precursor-evidence extraction failed for %s: %s", file_path, exc)
        return {}
    return out


def parse_mzml_file(
    file_path,
    isolation_window_width_da=None,
    isolation_window_width_source=None,
    allow_unknown_representation=False,
):
    """
    解析 mzML 文件，返回谱图列表。
    每个谱图为 dict: {spectrum_id, precursor_mz, peaks, title, retention_time, scan_number}
    """
    records = []
    file_pol = _file_polarity(file_path)  # mzML 极性 cvParam 文件级兜底 (matchms 丢失 ionmode)
    file_instrument = extract_mzml_instrument_metadata(file_path)
    precursor_evidence = _extract_mzml_precursor_evidence(
        file_path,
        isolation_window_width_da=isolation_window_width_da,
        isolation_window_width_source=isolation_window_width_source,
        include_ordinal_index=True,
    )
    spectra = load_from_mzml(str(file_path))
    local_idx = 0
    for stream_idx, s in enumerate(spectra):
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
        ordinal_meta = precursor_evidence.get(("ms2_ordinal", stream_idx), {})
        # matchms 把 mzML 的 'scan start time' 存成 list/tuple/标量都有可能,
        # 这里统一抽第一个标量并保留 mzML CV 单位语义, 避免
        # OpenMS 输出的 second 被误当成 minute，也避免 float([0.30]) 抛 TypeError。
        # (旧版直接 float(list) 失败 → 整个 mzML 全部 RT 丢失 → 前端 TIC 概览渲染不出来)
        def _coerce_rt(val):
            return _retention_time_minutes(val)
        rt = _safe_float(ordinal_meta.get("retention_time"))
        if rt is None:
            rt = _coerce_rt(md.get("scan_start_time"))
        if rt is None:
            rt = _coerce_rt(md.get("retention_time"))
        title = md.get("title") or ordinal_meta.get("native_id") or ""
        scan_id = _extract_scan_id(title)
        if scan_id is None:
            scan_id = ordinal_meta.get("native_scan_id")
        if scan_id is None:
            scan_id = local_idx
        precursor_meta = precursor_evidence.get(scan_id) or ordinal_meta
        spectrum_representation = str(
            precursor_meta.get("spectrum_representation") or "unknown"
        )
        if spectrum_representation == "profile":
            raise ValueError(
                "mzML contains an unpicked profile MS2 spectrum. Run the "
                "per-spectrum vendor/OpenMS peak-picking stage before matching."
            )
        if (
            spectrum_representation == "unknown"
            and not allow_unknown_representation
        ):
            raise ValueError(
                "mzML omits centroid/profile CV terms for an MS2 spectrum. "
                "Use OpenMS per-spectrum representation detection or explicitly "
                "choose the audited assume-centroid policy."
            )
        # 极性 (准确度修复 #1/#4: 激活加合物中性质量多假设匹配)
        # matchms 从 mzML 极性 cvParam (MS:1000130 positive / MS:1000129 negative) 派生 ionmode
        _ion = str(md.get("ionmode") or "").lower()
        polarity = "POS" if "positive" in _ion or _ion == "pos" else \
                   "NEG" if "negative" in _ion or _ion == "neg" else None
        if polarity is None and file_pol in {"POS", "NEG"}:
            polarity = file_pol  # matchms 丢失 mzML 极性 cvParam 时按文件级多数兜底
        adduct_value = _metadata_value(
            md, "adduct", "precursor_type", "precursortype", "ion_type", "iontype"
        )
        adduct = str(adduct_value).strip() if adduct_value is not None else ""
        adduct_mode = adduct_polarity(adduct)
        polarity_conflict = bool(polarity and adduct_mode and polarity != adduct_mode)
        normalized_adduct = re.sub(r"\s+", "", adduct).upper()
        adduct_reliable = bool(
            not polarity_conflict
            and normalized_adduct in _SPECIFIC_ADDUCTS
            and neutral_mass(pmz, adduct) is not None
        )
        records.append({
            "spectrum_id": scan_id,
            "precursor_mz": pmz,
            "peaks": peaks_arr,
            "title": str(title),
            "retention_time": rt,
            "retention_time_unit": "minutes",
            "scan_number": md.get("scan_number") or precursor_meta.get(
                "native_scan_id"
            ),
            "polarity": polarity,
            "adduct": adduct,
            "adduct_reliable": adduct_reliable,
            "polarity_conflict": polarity_conflict,
            "polarity_source": "scan_cv" if _ion else "file_cv" if polarity else None,
            "file_polarity": file_pol,
            "precursor_charge": precursor_meta.get("precursor_charge"),
            "precursor_intensity": precursor_meta.get("precursor_intensity"),
            "precursor_ms1_intensity": precursor_meta.get("precursor_ms1_intensity"),
            "collision_energy": precursor_meta.get("collision_energy"),
            "collision_energy_unit": precursor_meta.get("collision_energy_unit", ""),
            "collision_energy_mode": precursor_meta.get(
                "collision_energy_mode", "missing"
            ),
            "activation_method": precursor_meta.get("activation_method", ""),
            "activation_raw": precursor_meta.get("activation_raw", ""),
            "collision_cell_region": precursor_meta.get(
                "collision_cell_region", ""
            ),
            "instrument_type": (
                _metadata_value(
                    md, "instrument_type", "instrumenttype",
                    "mass_analyzer", "massanalyzer",
                )
                or file_instrument.get("instrument_type", "")
            ),
            "instrument": (
                _metadata_value(
                    md, "instrument", "instrument_model", "instrumentmodel",
                )
                or file_instrument.get("instrument", "")
            ),
            "precursor_scan_ref": precursor_meta.get("precursor_scan_ref", ""),
            "isolation_window_target_mz": precursor_meta.get("isolation_window_target_mz"),
            "isolation_window_lower_offset": precursor_meta.get("isolation_window_lower_offset"),
            "isolation_window_upper_offset": precursor_meta.get("isolation_window_upper_offset"),
            "precursor_isotope_evidence": precursor_meta.get("precursor_isotope_evidence"),
            "isotope_evidence_status": precursor_meta.get(
                "isotope_evidence_status", "metadata_unavailable"
            ),
            "precursor_isolation_evidence": precursor_meta.get(
                "precursor_isolation_evidence",
                {"status": "metadata_unavailable", "purity": None},
            ),
            "spectrum_representation": spectrum_representation,
            "spectrum_representation_source": precursor_meta.get(
                "spectrum_representation_source", "missing_cv"
            ),
            "raw_peaks_complete": True,
        })
        local_idx += 1
    return records


def _extract_scan_id(title):
    if not title:
        return None
    m = _SCAN_RE.search(str(title))
    return int(m.group(1)) if m else None


def _mzml_chromatogram_time_unit(file_path):
    """Return the explicit mzML time-array unit without decoding binary data."""
    try:
        overlap = 4096
        previous = b""
        with open(file_path, "rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                data = previous + chunk
                position = data.find(b"MS:1000595")  # time array
                if position >= 0:
                    start = data.rfind(b"<cvParam", 0, position)
                    end = data.find(b">", position)
                    if start >= 0 and end >= 0:
                        tag = data[start:end + 1].lower()
                        if b"uo:0000028" in tag or b"millisecond" in tag:
                            return "millisecond"
                        if b"uo:0000010" in tag or b'unitname="second"' in tag:
                            return "second"
                        if b"uo:0000032" in tag or b'unitname="hour"' in tag:
                            return "hour"
                        if b"uo:0000031" in tag or b'unitname="minute"' in tag:
                            return "minute"
                previous = data[-overlap:]
    except Exception:
        return "minute"
    return "minute"


def _chromatogram_kind(chromatogram):
    normalized_keys = {
        re.sub(r"[^a-z0-9]", "", str(key or "").lower())
        for key in chromatogram
    }
    identifier = re.sub(
        r"[^a-z0-9]", "", str(chromatogram.get("id") or "").lower()
    )
    if "totalioncurrentchromatogram" in normalized_keys or identifier == "tic":
        return "tic"
    if "basepeakchromatogram" in normalized_keys or identifier in {"bpc", "bpi"}:
        return "bpc"
    return None


def _numeric_chromatogram_arrays(chromatogram, *, default_time_unit):
    times = chromatogram.get("time array")
    intensities = chromatogram.get("intensity array")
    if times is None or intensities is None:
        return None
    try:
        count = min(len(times), len(intensities))
    except TypeError:
        return None
    if count <= 0:
        return None

    levels = chromatogram.get("ms level")
    try:
        levels_available = levels is not None and len(levels) == count
    except TypeError:
        levels_available = False

    normalized_times = []
    normalized_intensities = []
    normalized_levels = [] if levels_available else None
    for index in range(count):
        rt = _retention_time_minutes(
            times[index], default_unit=default_time_unit,
        )
        intensity = _safe_float(intensities[index])
        if rt is None or intensity is None:
            normalized_times.append(None)
            normalized_intensities.append(None)
        else:
            normalized_times.append(rt)
            normalized_intensities.append(max(0.0, intensity))
        if levels_available:
            level = _safe_float(levels[index])
            normalized_levels.append(int(level) if level is not None else None)
    return {
        "times": normalized_times,
        "intensities": normalized_intensities,
        "levels": normalized_levels,
    }


def _downsample_ordered_indices(indices, limit):
    limit = max(1, int(limit or 1))
    if len(indices) <= limit:
        return indices
    stride = max(1, (len(indices) + limit - 1) // limit)
    sampled = indices[::stride]
    if sampled[-1] != indices[-1]:
        sampled.append(indices[-1])
    return sampled


def _extract_ms1_tic_chromatogram(file_path, max_scans, mzml_module):
    """Recover MS1 TIC points stored only in mzML ``chromatogramList``.

    Some ProteoWizard vendor exports omit every MS1 ``spectrum`` while retaining
    a source TIC chromatogram plus a non-standard ``ms level`` array.  Those
    points are valid for drawing TIC, but they do not contain the peak arrays
    needed to reconstruct BPC.  A source BPC chromatogram is used only when it
    is explicitly present and can be aligned to the MS1-labelled TIC points.
    """
    time_unit = _mzml_chromatogram_time_unit(file_path)
    decoded = {}
    try:
        with mzml_module.MzML(str(file_path), decode_binary=True) as reader:
            for chromatogram in reader.iterfind("chromatogram"):
                kind = _chromatogram_kind(chromatogram)
                if kind and kind not in decoded:
                    arrays = _numeric_chromatogram_arrays(
                        chromatogram, default_time_unit=time_unit,
                    )
                    if arrays:
                        decoded[kind] = arrays
    except Exception as exc:
        logger.debug(
            "mzML source-chromatogram extraction failed for %s: %s",
            file_path,
            exc,
        )
        return None

    tic = decoded.get("tic")
    if not tic or tic.get("levels") is None:
        return None
    ms1_indices = [
        index for index, level in enumerate(tic["levels"])
        if level == 1
        and tic["times"][index] is not None
        and tic["intensities"][index] is not None
    ]
    if not ms1_indices:
        return None

    bpc = decoded.get("bpc")
    bpc_aligned = bool(
        bpc
        and len(bpc["times"]) == len(tic["times"])
        and all(
            tic["times"][index] is not None
            and bpc["times"][index] is not None
            and abs(tic["times"][index] - bpc["times"][index]) <= 1e-7
            for index in (ms1_indices[0], ms1_indices[-1])
        )
    )
    selected_indices = _downsample_ordered_indices(ms1_indices, max_scans)
    points = []
    for index in selected_indices:
        bpc_value = None
        if bpc_aligned:
            bpc_level = bpc["levels"][index] if bpc.get("levels") is not None else 1
            if bpc_level == 1:
                bpc_value = bpc["intensities"][index]
        points.append({
            "rt": tic["times"][index],
            "tic": tic["intensities"][index],
            "bpc": bpc_value,
            "base_peak_mz": None,
            "scan_id": index,
        })
    points.sort(key=lambda point: point["rt"])
    bpc_available = bool(points) and all(
        point.get("bpc") is not None for point in points
    )
    return {
        "points": points,
        "total_ms1_points": len(ms1_indices),
        "bpc_available": bpc_available,
        "bpc_source": (
            "mzml_basepeak_chromatogram" if bpc_available else "unavailable"
        ),
    }


def extract_tic_overview(file_path, max_scans=20000):
    """
    用 pyteomics 直读 mzML 头信息(不解码 peaks 数组)抽取每条 scan 的
    (rt, tic, bpc, base_peak_mz, ms_level)，用于前端 TIC/BPC 全扫描概览渲染。
    MS1 给连续色谱图，MS2 仅做采集位置标记参考。

    matchms.load_from_mzml 只返回带 precursor 的 MS2+, 拿不到 MS1, 所以另起一个轻量通道.
    速度: 2250 MS2 + 892 MS1 的 100MB 文件约 1-2s, 跟分析主流程比可忽略.

    返回: dict {
        "ms1": [{"rt": float (min), "tic": float, "bpc": float|None,
                 "base_peak_mz": float|None, "scan_id": int}, ...],
        "ms2": [{"rt": float, "tic": float, "bpc": float|None,
                 "base_peak_mz": float|None, "scan_id": int,
                 "precursor_mz": float|None}, ...],
        "rt_unit": "minute",
        "bpc_available": bool,
    }
    若文件不是 mzML 或解析失败则返回 None — 调用方据此走 fallback(从 results 兜底).
    """
    ext = Path(file_path).suffix.lower()
    if ext not in (".mzml", ".mzxml", ".xml"):
        return None
    try:
        from pyteomics import mzml as _pmzml
    except Exception:
        return None

    ms1, ms2 = [], []
    try:
        # TIC/BPI are spectrum-level CV values in standards-compliant mzML.
        # Keeping binary arrays undecoded makes lazy backfill practical even for
        # multi-gigabyte profile files and avoids allocating hundreds of millions
        # of intensity values merely to draw chromatograms.
        with _pmzml.read(str(file_path), decode_binary=False) as reader:
            for i, s in enumerate(reader):
                if i >= max_scans:
                    break
                if is_auxiliary_detector_spectrum(s):
                    continue
                lvl = s.get("ms level")
                if lvl not in (1, 2):
                    continue
                # Normalize the explicit mzML CV unit. OpenMS commonly writes
                # seconds, while vendor exports commonly write minutes.
                rt = _scan_start_time_minutes(s)
                if rt is None:
                    continue
                tic = s.get("total ion current")
                try:
                    tic = float(tic) if tic is not None else 0.0
                except (TypeError, ValueError):
                    tic = 0.0
                bpc = s.get("base peak intensity")
                try:
                    bpc = float(bpc) if bpc is not None else None
                except (TypeError, ValueError):
                    bpc = None
                base_peak_mz = s.get("base peak m/z")
                try:
                    base_peak_mz = float(base_peak_mz) if base_peak_mz is not None else None
                except (TypeError, ValueError):
                    base_peak_mz = None
                # Tests and uncommon readers may still provide decoded numeric
                # arrays. Use them only as a safe fallback; pyteomics
                # binary_array_record objects are deliberately not decoded here.
                if tic <= 0:
                    try:
                        intens = s.get("intensity array")
                        if (
                            intens is not None
                            and not isinstance(intens, dict)
                            and not callable(getattr(intens, "decode", None))
                        ):
                            tic = float(sum(float(x) for x in intens))
                    except Exception:
                        tic = 0.0
                if bpc is None or bpc < 0:
                    try:
                        intens = s.get("intensity array")
                        if (
                            intens is not None
                            and not isinstance(intens, dict)
                            and not callable(getattr(intens, "decode", None))
                            and len(intens)
                        ):
                            bpc = float(max(float(x) for x in intens))
                    except Exception:
                        bpc = None
                # scan id 从 native id 提取(用于和 results 里的 spectrum_id 关联)
                native_id = s.get("id") or ""
                m = _SCAN_RE.search(str(native_id))
                if not m:
                    m = re.search(r"scan=(\d+)", str(native_id))
                scan_id = int(m.group(1)) if m else i
                point = {
                    "rt": rt,
                    "tic": tic,
                    "bpc": bpc,
                    "base_peak_mz": base_peak_mz,
                    "scan_id": scan_id,
                }
                if lvl == 1:
                    ms1.append(point)
                else:
                    pmz = None
                    try:
                        pl = s.get("precursorList") or {}
                        precs = pl.get("precursor") or []
                        if precs:
                            sil = precs[0].get("selectedIonList") or {}
                            si = sil.get("selectedIon") or []
                            if si:
                                pmz = float(si[0].get("selected ion m/z"))
                    except (TypeError, ValueError, KeyError):
                        pmz = None
                    ms2.append({**point, "precursor_mz": pmz})
    except Exception:
        return None

    ms1_spectrum_count = len(ms1)
    chromatogram_fallback = None
    if not ms1 and ext == ".mzml":
        chromatogram_fallback = _extract_ms1_tic_chromatogram(
            file_path, max_scans, _pmzml,
        )
        if chromatogram_fallback:
            ms1 = chromatogram_fallback["points"]

    if not ms1 and not ms2:
        return None
    bpc_available = bool(ms1) and all(point.get("bpc") is not None for point in ms1)
    tic_source = (
        "mzml_ms1_spectrum_metadata"
        if ms1_spectrum_count
        else (
            "mzml_tic_chromatogram_ms_level_array"
            if chromatogram_fallback
            else "ms2_spectrum_metadata_approximation"
        )
    )
    bpc_source = (
        "mzml_ms1_spectrum_metadata"
        if bpc_available and ms1_spectrum_count
        else (
            chromatogram_fallback.get("bpc_source", "unavailable")
            if chromatogram_fallback
            else "unavailable"
        )
    )
    return {
        "chromatogram_schema_version": 3,
        "ms1": ms1,
        "ms2": ms2,
        "rt_unit": "minute",
        "tic_source": tic_source,
        "ms1_spectra_available": bool(ms1_spectrum_count),
        "ms1_spectrum_count": ms1_spectrum_count,
        "ms1_tic_point_count": (
            chromatogram_fallback.get("total_ms1_points", len(ms1))
            if chromatogram_fallback else len(ms1)
        ),
        "bpc_available": bpc_available,
        "bpc_source": bpc_source,
        "bpc_scope": (
            "MS1 base peak intensity per survey scan"
            if bpc_available
            else (
                "unavailable_tic_only_without_ms1_peak_arrays_or_bpc_chromatogram"
                if chromatogram_fallback
                else "unavailable_without_continuous_ms1"
            )
        ),
    }
