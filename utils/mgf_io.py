"""MGF 文件解析"""
import re

from matchms.importing import load_from_mgf

from utils.adduct import NEUTRAL_SHIFTS, adduct_polarity, neutral_mass
from utils.collision_energy import parse_collision_energy

_SCAN_RE = re.compile(r"scanId=(\d+)")
_POLARITY_ONLY_ADDUCTS = {"POSITIVE", "POS", "+", "NEGATIVE", "NEG", "-"}
_SUPPORTED_SPECIFIC_ADDUCTS = {
    re.sub(r"\s+", "", key).upper(): value for key, value in NEUTRAL_SHIFTS.items()
}


def _metadata_value(metadata, *aliases):
    """Return the first non-empty metadata value, case/underscore insensitive."""
    wanted = {re.sub(r"[^a-z0-9]", "", alias.lower()) for alias in aliases}
    for key, value in metadata.items():
        normalized = re.sub(r"[^a-z0-9]", "", str(key or "").lower())
        if normalized in wanted and value is not None and str(value).strip():
            return value
    return None


def _charge_magnitude(charge):
    if charge is None or str(charge).strip() == "":
        return None
    try:
        return abs(int(charge))
    except (TypeError, ValueError):
        match = re.search(r"\d+", str(charge))
        return int(match.group()) if match else None


def _float_metadata(metadata, *aliases):
    value = _metadata_value(metadata, *aliases)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        match = re.search(r"[-+]?\d+(?:\.\d+)?", str(value))
        return float(match.group()) if match else None


def _collision_energy_metadata(metadata):
    aliases = (
        "collision_energy", "collisionenergy", "collision_energy_ev", "ce",
        "normalized_collision_energy", "nce",
    )
    raw = _metadata_value(metadata, *aliases)
    value = _float_metadata(metadata, *aliases)
    if value is None:
        return None, "", "missing"
    raw_text = str(raw or "").lower()
    unit = ""
    if re.search(r"\bnce\b|normalized\s*collision\s*energy|%", raw_text):
        unit = "NCE"
    elif re.search(r"\bev\b|electron\s*volt", raw_text):
        unit = "eV"
    if not unit:
        for key in metadata:
            normalized = re.sub(r"[^a-z0-9]", "", str(key or "").lower())
            if normalized in {"nce", "normalizedcollisionenergy"}:
                unit = "NCE"
                break
            if normalized == "collisionenergyev":
                unit = "eV"
                break
    parsed = parse_collision_energy(raw, unit)
    parsed_mode = str(parsed.get("mode") or "")
    mode = "fixed" if parsed_mode == "single" else parsed_mode
    return value, unit, mode


def _activation_metadata(metadata):
    value = _metadata_value(
        metadata,
        "activation",
        "activation_method",
        "activation_type",
        "fragmentation_method",
        "dissociation_method",
    )
    return str(value or "").strip()


def _specific_adduct_is_reliable(adduct, precursor_mz, charge=None):
    """A polarity label is not an adduct and must not create a neutral mass."""
    text = str(adduct or "").strip()
    normalized = re.sub(r"\s+", "", text).upper()
    if normalized in _POLARITY_ONLY_ADDUCTS or normalized not in _SUPPORTED_SPECIFIC_ADDUCTS:
        return False
    expected_charge = abs(int(_SUPPORTED_SPECIFIC_ADDUCTS[normalized][1]))
    observed_charge = _charge_magnitude(charge)
    if observed_charge is not None and observed_charge != expected_charge:
        return False
    return neutral_mass(precursor_mz, text) is not None


def parse_mgf_file(file_path):
    """
    解析 MGF 文件，返回谱图列表。
    每个谱图为 dict: {spectrum_id, precursor_mz, peaks, title, retention_time, scan_number}
    """
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
        # matchms 导入 MGF 时会把 RTINSECONDS= 规范化成 "retention_time" 键(单位秒)。
        # 旧代码找 "rtinsECONDS"/"rtinseconds" 两个不存在的键 → RT 永远 None, TIC 画不出。
        rt = md.get("retention_time")
        rt_unit = "seconds"
        if rt is None:
            rt = _metadata_value(md, "rtinseconds")
        if rt is None:
            rt = _metadata_value(md, "rtinminutes")
            if rt is not None:
                rt_unit = "minutes"
        if rt is not None:
            try:
                rt = float(rt)
            except (TypeError, ValueError):
                rt = None
        title = md.get("title") or ""
        scan_id = _extract_scan_id(title) or local_idx
        # Preserve an explicitly supplied MGF adduct.  matchms normally maps
        # PRECURSOR_TYPE to ``adduct`` but vendor exports also use ADDUCT,
        # IONTYPE, ION_TYPE or PRECURSORTYPE; inspect aliases defensively.
        adduct_value = _metadata_value(
            md, "adduct", "precursor_type", "precursortype", "ion_type", "iontype"
        )
        adduct = str(adduct_value).strip() if adduct_value is not None else ""

        # Infer polarity from IONMODE, CHARGE and ADDUCT independently.  A
        # disagreement is evidence corruption: retain the raw fields but do
        # not silently choose one polarity or enable cross-adduct matching.
        # matchms charge is usually int(-1/1), occasionally "1-"/"2+".
        charge = md.get("charge")
        charge_magnitude = _charge_magnitude(charge)
        (
            collision_energy,
            collision_energy_unit,
            collision_energy_mode,
        ) = _collision_energy_metadata(md)
        activation_method = _activation_metadata(md)
        ionmode = str(md.get("ionmode") or md.get("ion_mode") or "").strip().lower()
        declared_polarity = (
            "POS" if ionmode in {"pos", "positive", "+"} or "positive" in ionmode
            else "NEG" if ionmode in {"neg", "negative", "-"} or "negative" in ionmode
            else None
        )
        charge_polarity = None
        try:
            ci = int(charge)
            charge_polarity = "NEG" if ci < 0 else ("POS" if ci > 0 else None)
        except (TypeError, ValueError):
            cs = str(charge)
            if "-" in cs:
                charge_polarity = "NEG"
            elif "+" in cs:
                charge_polarity = "POS"
        adduct_mode = adduct_polarity(adduct)
        supplied_modes = {
            mode for mode in (declared_polarity, charge_polarity, adduct_mode) if mode
        }
        polarity_conflict = len(supplied_modes) > 1
        polarity = None if polarity_conflict else next(iter(supplied_modes), None)
        adduct_reliable = bool(
            not polarity_conflict
            and _specific_adduct_is_reliable(adduct, pmz, charge)
        )
        records.append({
            "spectrum_id": scan_id,
            "precursor_mz": pmz,
            "peaks": peaks_arr,
            "title": str(title),
            "retention_time": rt,
            "retention_time_unit": rt_unit,
            "charge": charge,
            "precursor_charge": charge_magnitude,
            "collision_energy": collision_energy,
            "collision_energy_unit": collision_energy_unit,
            "collision_energy_mode": collision_energy_mode,
            "activation_method": activation_method,
            "activation_raw": activation_method,
            "collision_cell_region": str(_metadata_value(
                md,
                "collision_cell_region",
                "collision_cell",
                "collision_region",
            ) or ""),
            "instrument_type": str(_metadata_value(
                md, "instrument_type", "instrumenttype",
                "mass_analyzer", "massanalyzer",
            ) or ""),
            "instrument": str(_metadata_value(
                md, "instrument", "instrument_model", "instrumentmodel",
            ) or ""),
            "precursor_intensity": _float_metadata(
                md, "precursor_intensity", "precursorintensity"
            ),
            # An MGF contains isolated MS/MS peak lists, never the associated
            # MS1 envelope.  State this explicitly instead of fabricating an
            # isotope pattern from fragment ions.
            "isotope_evidence_status": "ms1_not_available_in_mgf",
            "precursor_isotope_evidence": None,
            "spectrum_representation": "centroid",
            "spectrum_representation_source": "mgf_peak_list",
            "raw_peaks_complete": True,
            "adduct": adduct,
            "adduct_reliable": adduct_reliable,
            "polarity": polarity,
            "polarity_conflict": polarity_conflict,
            "polarity_source": (
                "IONMODE" if declared_polarity
                else "CHARGE" if charge_polarity
                else "ADDUCT" if adduct_mode
                else None
            ),
            "scan_number": None,
        })
        local_idx += 1
    return records


def _extract_scan_id(title):
    if not title:
        return None
    m = _SCAN_RE.search(str(title))
    return int(m.group(1)) if m else None
