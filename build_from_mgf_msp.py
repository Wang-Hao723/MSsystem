#!/usr/bin/env python3
"""
从本地 MGF/MSP 文件直接构建质谱数据库。

数据来源：
  - MGF_classify_data/    : COMMENT=NAME: 格式（含化合物名、SMILES、分类）
  - MGF_other_data/       : NAME= / SMILES= 顶级字段格式
  - MSP_data/             : 标准 MSP 格式

完全跳过 HDF5，确保所有元数据直接从原始文件提取。
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import re
import sqlite3
import sys
import time
import uuid
import warnings
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from config.settings import (
    DB_PATH as SHARED_DB_PATH,
    DB_PATH_SMALL as SHARED_DB_PATH_SMALL,
    LIBRARY_WRITE_ENABLED,
    USER_CONTRIB_DIR,
)
from utils.storage_guard import (
    estimate_sqlite_rebuild_bytes,
    exclusive_file_lock,
    prune_generated_backups,
    require_storage_capacity,
)
from utils.collision_energy import parse_collision_energy

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("build_from_mgf_msp")

warnings.filterwarnings("ignore")
try:
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    RDLogger.DisableLog("rdChem.*")
except ImportError:
    pass

# ─────────────────────────────────────────────────────────────────
# 数据库 Schema
# ─────────────────────────────────────────────────────────────────
OUTPUT_DIR = Path("output")
DB_PATH = OUTPUT_DIR / "local_library.db"
MAPPING_PATH = OUTPUT_DIR / "library_mapping.csv"
DEFAULT_LIBRARY_ROOT = Path(
    "/Users/zizi/Desktop/MsQCModel_v2/_local_only/mgf_msp_database"
)

BATCH_SIZE = 5000

# These tiers may be retained in the database for traceability or dedicated
# workflows, but are deliberately excluded by matching.search from ordinary
# experimental Top-K search.  Presence-only rows are not spectra at all and
# therefore never enter the Spectra table during a full rebuild.
_NON_SPECTRAL_SOURCE_TIERS = {"literature_presence_only", "quarantined"}
_KNOWN_SOURCE_TIERS = {
    "experimental", "merged", "curated", "predicted", "in_silico",
    "unknown", "user_unverified", "metadata_conflict", "quarantined",
    "literature_presence_only", "massbank_supplement",
}

# ─────────────────────────────────────────────────────────────────
# MGF 解析器（支持两种格式）
# ─────────────────────────────────────────────────────────────────
_RE_COMMENT_SMILES = re.compile(r"SMILES[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_INCHI = re.compile(r"INCHI[:=\s]+(InChI[^\s;]+)", re.IGNORECASE)
_RE_COMMENT_INCHIKEY = re.compile(
    r"\bINCHIKEY\s*[:=]\s*([A-Z]{14}-[A-Z]{10}-[A-Z])\b",
    re.IGNORECASE,
)
_RE_COMMENT_FORMULA = re.compile(r"FORMULA[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_NAME = re.compile(
    r"(?:^|\s)NAME\s*[:=]\s*(.+?)(?="
    r"\s+(?:SMILES|INCHI|INCHIKEY|FORMULA|EXACT_MASS|ADDUCT|ION_MODE|USER_CLASS|COMPOUND_CLASS|INSTRUMENT(?:_TYPE)?)\s*[:=]|;|$)",
    re.IGNORECASE,
)
_RE_COMMENT_EXACT_MASS = re.compile(r"EXACT_MASS[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_ADDUCT = re.compile(r"\[([M][^\]]+)\]", re.IGNORECASE)
_RE_COMMENT_ION_MODE = re.compile(r"ION_MODE[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_USER_CLASS = re.compile(r"USER_CLASS[:=\s]*([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_COMPOUND_CLASS = re.compile(r"COMPOUND_CLASS[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_INSTRUMENT_TYPE = re.compile(
    r"INSTRUMENT(?:_TYPE)?[:=\s]+(.+?)(?=\s+(?:SMILES|INCHI|INCHIKEY|FORMULA|NAME|ADDUCT|ION_MODE|USER_CLASS|COMPOUND_CLASS)[:=]|$)",
    re.IGNORECASE,
)
_RE_COMMENT_SOURCE_TIER = re.compile(r"SOURCE_TIER[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_SOURCE_ACCESSION = re.compile(
    r"(?:DB#|SOURCE_ACCESSION|SPECTRUMID|ACCESSION)\s*[:=]\s*([^\s;\"']+)",
    re.IGNORECASE,
)
_RE_COMMENT_COLLISION_ENERGY = re.compile(
    r"(?:COLLISIONENERGY|COLLISION_ENERGY)\s*[:=][ \t]*?(.*?)"
    r"(?=\s+(?:CONFIDENCE|SOURCE|COVERS_TARGETS|SOURCE_TIER|"
    r"SOURCE_ACCESSION)\s*=|;|$)",
    re.IGNORECASE,
)


def _parse_mgf_comment(comment: str) -> dict[str, Any]:
    if not comment:
        return {}
    block = " ".join(comment.split())

    def _first(pat, text):
        m = pat.search(text)
        if not m:
            return None
        val = m.group(1).strip().strip('"').strip("'")
        return val if val else None

    result: dict[str, Any] = {}
    for pat, key in [
        (_RE_COMMENT_SMILES, "smiles"),
        (_RE_COMMENT_INCHI, "inchi"),
        (_RE_COMMENT_INCHIKEY, "inchikey"),
        (_RE_COMMENT_FORMULA, "formula"),
        (_RE_COMMENT_NAME, "name"),
        (_RE_COMMENT_EXACT_MASS, "exact_mass"),
        (_RE_COMMENT_ADDUCT, "adduct"),
        (_RE_COMMENT_ION_MODE, "ion_mode"),
        (_RE_COMMENT_COLLISION_ENERGY, "collision_energy_raw"),
    ]:
        val = _first(pat, block)
        if val:
            result[key] = val

    uc = _first(_RE_COMMENT_USER_CLASS, block)
    cc = _first(_RE_COMMENT_COMPOUND_CLASS, block)
    result["compound_class"] = cc or uc
    instrument_type = _first(_RE_COMMENT_INSTRUMENT_TYPE, block)
    source_tier = _first(_RE_COMMENT_SOURCE_TIER, block)
    source_accession = _first(_RE_COMMENT_SOURCE_ACCESSION, block)
    if instrument_type:
        result["instrument_type"] = instrument_type
    if source_tier:
        result["source_tier"] = source_tier
    if source_accession:
        result["source_accession"] = source_accession
    return result


def _source_is_in_curated_root(source_file: Path | None) -> bool:
    if source_file is None:
        return False
    try:
        source_file.resolve().relative_to(DEFAULT_LIBRARY_ROOT.resolve())
        return True
    except (OSError, ValueError):
        return False


def _infer_source_tier(fields: dict[str, Any], source_file: Path | None = None) -> str:
    """Infer provenance conservatively; absence is never experimental proof."""
    explicit = str(fields.get("source_tier") or "").strip().lower()
    if explicit:
        return explicit
    source_text = str(source_file or "").lower()
    if any(token in source_text for token in (
        "upload_template_literature", "literature_target", "presence_only",
        "待新增数据库数据",
    )):
        return "literature_presence_only"
    evidence = " ".join(str(fields.get(key) or "") for key in (
        "instrument_type", "instrument", "spectrum_type", "description", "origin",
    )).lower()
    predicted_tokens = (
        "in-silico", "in silico", "predicted", "computed", "theoretical",
        "cfm-id", "cfmid", "sirius", "metfrag",
    )
    if any(token in evidence for token in predicted_tokens):
        return "predicted"
    # This local directory is the reviewed source manifest for a formal full
    # rebuild.  Files elsewhere need an explicit SOURCE_TIER instead of being
    # silently promoted merely because they contain numeric intensities.
    if _source_is_in_curated_root(source_file):
        return "experimental"
    return "user_unverified"


def _stable_spectrum_id(
    file_path: Path,
    local_idx: int,
    source_root: Path | None = None,
) -> str:
    """Include a source-relative path hash so duplicate stems cannot collide."""
    resolved = file_path.resolve()
    root = source_root.resolve() if source_root is not None else None
    if root is None and _source_is_in_curated_root(file_path):
        root = DEFAULT_LIBRARY_ROOT.resolve()
    try:
        identity = resolved.relative_to(root).as_posix() if root else file_path.as_posix()
    except ValueError:
        identity = resolved.as_posix()
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:12]
    return f"{file_path.stem}_{digest}_{local_idx}"


def _derive_inchikey(smiles: str):
    if not smiles:
        return None
    try:
        from rdkit import Chem, rdBase
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
        if mol:
            # MolToInchiKey expects a Mol object, NOT an InChI string
            return Chem.MolToInchiKey(mol)
    except Exception:
        pass
    return None


def _resolve_classification_safely(
    *,
    smiles: str | None,
    raw: Any,
    source_file: str,
) -> str:
    from mssystem.classification import resolve_compound_class

    try:
        from rdkit import rdBase
        with rdBase.BlockLogs():
            return resolve_compound_class(
                smiles=smiles, raw=raw, source_file=source_file
            )
    except ImportError:
        return resolve_compound_class(
            smiles=smiles, raw=raw, source_file=source_file
        )


def _infer_intensity_origin(fields: dict[str, Any], intensities: list[float]) -> str:
    explicit = str(fields.get("intensity_origin") or "").strip().lower()
    if explicit:
        return explicit
    finite_positive = [float(value) for value in intensities if float(value) > 0]
    if len(finite_positive) >= 3:
        rounded = {round(value, 12) for value in finite_positive}
        if len(rounded) == 1:
            # A flat list (commonly every ion = 1000) contains presence
            # information but no measured relative-intensity evidence.
            return "default"
    return "measured"


def _clean_peaks(mz_arr, int_arr, precursor_mz, top_n=50, rel_thresh=0.01, gap=1.5):
    import numpy as np
    if mz_arr.size == 0 or int_arr.size == 0:
        return []
    max_int = float(np.max(int_arr))
    if max_int <= 0:
        return []
    rel = int_arr / max_int
    keep = (rel >= rel_thresh) & (mz_arr < (precursor_mz - gap))
    mz_f = mz_arr[keep]
    int_f = int_arr[keep]
    if mz_f.size == 0:
        return []
    top_idx = np.argsort(int_f)[::-1][:top_n]
    mz_t = mz_f[top_idx]
    int_t = int_f[top_idx]
    order = np.argsort(mz_t)
    return [[float(mz_t[i]), float(int_t[i])] for i in order]


def parse_mgf_file(
    file_path: Path,
    global_idx: int = 0,
    source_root: Path | None = None,
) -> tuple[list[dict], int]:
    """解析单个 MGF 文件，支持两种格式：
    1. Pluskal Lab 格式：NAME=, SMILES= 等为顶级字段
    2. GNPS/MoNA 格式：字段在 COMMENT= 中
    """
    import numpy as np

    records = []
    raw_text = file_path.read_text(encoding="utf-8", errors="replace")
    local_idx = 0

    for block in raw_text.split("BEGIN IONS"):
        if not block.strip() or block.strip() == "END IONS":
            continue

        fields: dict[str, Any] = {}
        peaks_mz = []
        peaks_int = []

        for line in block.split("\n"):
            line = line.rstrip()
            if not line:
                continue

            # ── 顶级字段（Pluskal Lab 格式）──────────────────────
            if line.startswith("NAME="):
                fields["name"] = line[5:].strip()
            elif line.startswith("SMILES="):
                fields["smiles"] = line[7:].strip()
            elif line.startswith("INCHI="):
                fields["inchi"] = line[6:].strip()
            elif line.startswith("INCHIKEY="):
                fields["inchikey"] = line[9:].strip()
            elif line.startswith("INCHIAUX="):
                fields["inchikey"] = line[9:].strip()
            elif line.startswith("EXACTMASS="):
                # ``EXACTMASS=`` is ten characters long.  The historical
                # fixed offset of 11 silently dropped the first digit
                # (408.2 -> 8.2), so split on the delimiter instead.
                fields["exact_mass"] = line.partition("=")[2].strip()
            elif line.startswith("FORMULA="):
                fields["formula"] = line[8:].strip()
            elif line.startswith("ADDUCT="):
                fields["adduct"] = line[7:].strip()
            elif line.startswith("IONMODE="):
                fields["ion_mode"] = line[8:].strip()
            elif line.startswith("COLLISIONENERGY="):
                fields["collision_energy_raw"] = line[16:].strip()
            elif line.startswith("COLLISION_ENERGY="):
                fields["collision_energy_raw"] = line[17:].strip()
            elif line.startswith("DESCRIPTION="):
                fields["description"] = line[12:].strip()
            elif line.startswith("INSTRUMENT_TYPE="):
                fields["instrument_type"] = line[16:].strip()
            elif line.startswith("INSTRUMENT="):
                fields["instrument"] = line[11:].strip()
            elif line.startswith("SPECTRUM_TYPE="):
                fields["spectrum_type"] = line[14:].strip()
            elif line.startswith("SPECTRUMID="):
                fields["source_accession"] = line[11:].strip()
            elif line.startswith("SOURCE_ACCESSION="):
                fields["source_accession"] = line[17:].strip()
            elif line.startswith("ACCESSION="):
                fields["source_accession"] = line[10:].strip()
            elif line.startswith("USI="):
                usi = line[4:].strip()
                fields["usi"] = usi
                if ":accession:" in usi.casefold():
                    fields.setdefault(
                        "source_accession",
                        re.split(r":accession:", usi, flags=re.IGNORECASE, maxsplit=1)[1],
                    )
            elif line.startswith("COMPOUND_NAME="):
                fields["name"] = line[14:].strip()
            elif line.startswith("SOURCE_TIER="):
                fields["source_tier"] = line[12:].strip()
            elif line.startswith("INTENSITY_ORIGIN="):
                fields["intensity_origin"] = line[17:].strip()
            # ── 标准 MGF 字段 ────────────────────────────────────
            elif line.startswith("TITLE="):
                fields["title"] = line[6:].strip()
            elif line.startswith("PEPMASS="):
                try:
                    fields["precursor_mz"] = float(line[8:].split()[0])
                except (ValueError, IndexError):
                    pass
            elif line.startswith("CHARGE="):
                fields["charge_raw"] = line[7:].strip()
            elif line.startswith("RTINSECONDS="):
                try:
                    fields["retention_time"] = float(line[13:].strip())
                except ValueError:
                    pass
            elif line.startswith("COMMENT="):
                parsed = _parse_mgf_comment(line[8:].strip())
                for k, v in parsed.items():
                    if v and k not in fields:
                        fields[k] = v
            elif line.startswith("ION_MODE="):
                fields["ion_mode"] = line[9:].strip()
            # ── 峰数据行 ────────────────────────────────────────
            elif "\t" in line or " " in line.strip():
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        peaks_mz.append(float(parts[0]))
                        peaks_int.append(float(parts[1]))
                    except ValueError:
                        pass

        precursor_mz = fields.get("precursor_mz")
        if precursor_mz is None or not peaks_mz:
            continue
        if re.search(r"\bms\s*1\b", str(fields.get("spectrum_type") or ""), re.IGNORECASE):
            continue

        peaks = _clean_peaks(
            np.array(peaks_mz, dtype=float),
            np.array(peaks_int, dtype=float),
            float(precursor_mz),
        )
        if not peaks:
            continue

        spectrum_id = _stable_spectrum_id(file_path, local_idx, source_root)
        raw_name = fields.get("name") or fields.get("title") or ""
        name = raw_name if raw_name.strip() else spectrum_id

        smiles = fields.get("smiles")
        inchikey = fields.get("inchikey")
        if not inchikey and smiles:
            inchikey = _derive_inchikey(smiles)

        compound_class = fields.get("compound_class") or ""
        # classification 统一用 resolve_compound_class (SMARTS 结构规则 → 文件名归一化 → other),
        # 不再用 file_path.stem 兜底(那是数据集名污染的根源)。
        _raw_cls = compound_class or (file_path.stem if file_path.parent.name == "MGF_classify_data" else "")
        classification = _resolve_classification_safely(
            smiles=smiles, raw=_raw_cls, source_file=str(file_path)
        )
        compound_class = classification
        collision_energy = parse_collision_energy(
            fields.get("collision_energy_raw")
        )

        rec = {
            "spectrum_id": spectrum_id,
            "source_file": str(file_path),
            "precursor_mz": float(precursor_mz),
            "peaks": peaks,
            "name": name,
            "retention_time": _coerce_float(fields.get("retention_time")),
            "charge": _coerce_int(fields.get("charge_raw")),
            "adduct": fields.get("adduct"),
            "ion_mode": fields.get("ion_mode"),
            "collision_energy": collision_energy["value"],
            "collision_energy_raw": collision_energy["raw"],
            "collision_energy_unit": collision_energy["unit"],
            "collision_energy_mode": collision_energy["mode"],
            "collision_energy_min": collision_energy["minimum"],
            "collision_energy_max": collision_energy["maximum"],
            "instrument_type": fields.get("instrument_type"),
            "instrument": fields.get("instrument"),
            "smiles": smiles,
            "inchi": fields.get("inchi"),
            "inchikey": inchikey,
            "formula": fields.get("formula"),
            "exact_mass": _coerce_float(fields.get("exact_mass")),
            "compound_class": compound_class,
            "classification": classification,
            "source_tier": _infer_source_tier(fields, file_path),
            "intensity_origin": _infer_intensity_origin(fields, peaks_int),
            "source_accession": fields.get("source_accession"),
        }
        records.append(rec)
        local_idx += 1
        global_idx += 1

    return records, global_idx


# ─────────────────────────────────────────────────────────────────
# MSP 解析器
# ─────────────────────────────────────────────────────────────────
def parse_msp_file(
    file_path: Path,
    global_idx: int = 0,
    source_root: Path | None = None,
) -> tuple[list[dict], int]:
    """解析标准 MSP 文件。"""
    import numpy as np

    records = []
    raw_text = file_path.read_text(encoding="utf-8", errors="replace")
    local_idx = 0

    for entry in re.split(r"\n\s*\n", raw_text):
        entry = entry.strip()
        if not entry:
            continue

        fields = {}
        peaks_mz = []
        peaks_int = []

        for line in entry.split("\n"):
            line = line.strip()
            if not line:
                continue
            if ":" in line and not line[0].isdigit():
                key, _, val = line.partition(":")
                key = key.strip().upper()
                key_normalized = re.sub(r"[\s_-]+", "", key)
                val = val.strip()
                if key_normalized == "NAME":
                    fields["name"] = val
                elif key_normalized in ("PRECURSORMZ", "PRECURSORM/Z", "PEPMASS"):
                    try:
                        fields["precursor_mz"] = float(val.split()[0])
                    except (ValueError, IndexError):
                        pass
                elif key_normalized == "MW" and "precursor_mz" not in fields:
                    # Fallback: some MoNA HMDB entries use MW instead of PRECURSORMZ
                    try:
                        fields["precursor_mz"] = float(val.split()[0])
                    except (ValueError, IndexError):
                        pass
                elif key_normalized in ("PRECURSORTYPE", "ADDUCT"):
                    fields["adduct"] = val
                elif key_normalized == "IONMODE":
                    fields["ion_mode"] = val
                elif key_normalized == "CHARGE":
                    fields["charge_raw"] = val
                elif key_normalized == "FORMULA":
                    fields["formula"] = val
                elif key_normalized == "SMILES":
                    fields["smiles"] = val
                elif key_normalized == "INCHIKEY":
                    fields["inchikey"] = val
                elif key_normalized == "INCHI":
                    fields["inchi"] = val
                elif key_normalized == "RETENTIONTIME":
                    try:
                        fields["retention_time"] = float(val)
                    except ValueError:
                        pass
                elif key_normalized == "COLLISIONENERGY":
                    fields["collision_energy_raw"] = val
                elif key_normalized == "COLLISIONENERGYUNIT":
                    fields["collision_energy_unit"] = val
                elif key_normalized == "EXACTMASS":
                    try:
                        fields["exact_mass"] = float(val)
                    except ValueError:
                        pass
                elif key_normalized in ("COMMENT", "COMMENTS"):
                    # MoNA MSP: Comments: "SMILES=..." "InChI=..." etc.
                    parsed = _parse_mgf_comment(val)
                    for ck, cv in parsed.items():
                        if cv and ck not in fields:
                            fields[ck] = cv
                elif key_normalized in (
                    "DB#",
                    "ACCESSION",
                    "SOURCEACCESSION",
                    "SPECTRUMID",
                ):
                    fields["source_accession"] = val
                elif key_normalized in ("COMPOUNDCLASS", "CLASS", "CATEGORY"):
                    fields["compound_class"] = val
                elif key_normalized == "IONSOURCE":
                    fields["ion_source"] = val
                elif key_normalized == "INSTRUMENTTYPE":
                    fields["instrument_type"] = val
                elif key_normalized == "INSTRUMENT":
                    fields["instrument"] = val
                elif key_normalized == "SPECTRUMTYPE":
                    fields["spectrum_type"] = val
                elif key_normalized == "SOURCETIER":
                    fields["source_tier"] = val
                elif key_normalized == "INTENSITYORIGIN":
                    fields["intensity_origin"] = val
                elif key_normalized == "NUMPEAKS":
                    pass
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        if mz > 0 and intensity > 0:
                            peaks_mz.append(mz)
                            peaks_int.append(intensity)
                    except ValueError:
                        pass

        precursor_mz = fields.get("precursor_mz")
        if precursor_mz is None or not peaks_mz:
            continue
        if re.search(r"\bms\s*1\b", str(fields.get("spectrum_type") or ""), re.IGNORECASE):
            continue

        peaks = _clean_peaks(
            np.array(peaks_mz, dtype=float),
            np.array(peaks_int, dtype=float),
            float(precursor_mz),
        )
        if not peaks:
            continue

        spectrum_id = _stable_spectrum_id(file_path, local_idx, source_root)
        raw_name = fields.get("name") or ""
        name = raw_name.strip() if raw_name.strip() else spectrum_id

        smiles = fields.get("smiles")
        inchikey = fields.get("inchikey")
        if not inchikey and smiles:
            inchikey = _derive_inchikey(smiles)

        # classification 统一走 resolve_compound_class, 避免 file_path.stem 当类别的污染
        compound_class = _resolve_classification_safely(
            smiles=smiles, raw=fields.get("compound_class"), source_file=str(file_path)
        )

        collision_energy = parse_collision_energy(
            fields.get("collision_energy_raw"), fields.get("collision_energy_unit")
        )

        rec = {
            "spectrum_id": spectrum_id,
            "source_file": str(file_path),
            "precursor_mz": float(precursor_mz),
            "peaks": peaks,
            "name": name,
            "retention_time": _coerce_float(fields.get("retention_time")),
            "charge": _coerce_int(fields.get("charge_raw")),
            "adduct": fields.get("adduct") or fields.get("ion_mode"),
            "ion_mode": fields.get("ion_mode"),
            "smiles": smiles,
            "inchi": fields.get("inchi"),
            "inchikey": inchikey,
            "formula": fields.get("formula"),
            "exact_mass": _coerce_float(fields.get("exact_mass")),
            "collision_energy": collision_energy["value"],
            "collision_energy_raw": collision_energy["raw"],
            "collision_energy_unit": collision_energy["unit"],
            "collision_energy_mode": collision_energy["mode"],
            "collision_energy_min": collision_energy["minimum"],
            "collision_energy_max": collision_energy["maximum"],
            "instrument_type": fields.get("instrument_type"),
            "instrument": fields.get("instrument"),
            "compound_class": compound_class,
            "classification": compound_class,
            "source_tier": _infer_source_tier(fields, file_path),
            "intensity_origin": _infer_intensity_origin(fields, peaks_int),
            "source_accession": fields.get("source_accession"),
        }
        records.append(rec)
        local_idx += 1
        global_idx += 1

    return records, global_idx


# ─────────────────────────────────────────────────────────────────
# 去重
# ─────────────────────────────────────────────────────────────────
def _dedup_key(rec: dict[str, Any]) -> str:
    """Spectrum-level key that never merges different named compounds.

    The historical fallback used precursor m/z alone and collapsed isobars;
    the InChIKey branch also discarded spectra measured at other collision
    energies/adducts.  Identity and spectral acquisition are separate axes, so
    both must participate in a deduplication key.
    """
    inchikey = str(rec.get("inchikey") or "").strip().upper()
    name = re.sub(r"\s+", " ", str(rec.get("name") or "").strip()).casefold()
    formula = re.sub(r"\s+", "", str(rec.get("formula") or "").strip())
    # A literature presence-only row is a target assertion, not a validated
    # structure record.  Keep distinct reported names even when a spreadsheet
    # has accidentally pasted the same SMILES/InChIKey into several rows.
    presence_only = str(rec.get("source_tier") or "").strip().lower() == "literature_presence_only"
    if inchikey and not presence_only:
        identity = f"ik::{inchikey}"
    elif name and formula:
        identity = f"name_formula::{name}|{formula}"
    else:
        identity = f"sid::{rec.get('spectrum_id', 'unknown')}"

    precursor = _coerce_float(rec.get("precursor_mz"))
    precursor_key = f"{precursor:.4f}" if precursor is not None else "missing"
    adduct = re.sub(r"\s+", "", str(rec.get("adduct") or "")).upper()
    peaks = rec.get("peaks") or []
    peak_key = ",".join(
        f"{float(peak[0]):.4f}:{float(peak[1]):.3g}"
        for peak in peaks
        if isinstance(peak, (list, tuple)) and len(peak) >= 2
    )
    return f"{identity}|pmz::{precursor_key}|adduct::{adduct}|peaks::{peak_key}"


# 这些是通用/无信息分类，在做 dedup 时应该被具体分类覆盖
_GENERIC_CLASSIFICATIONS = {
    "MoNA-export-All_Spectra", "MoNA-export-All_LC-MS-MS_QTOF",
    "MoNA_export_LC_MS_Spectra", "MoNA_export_GNPS",
    "MoNA-export-MassBank", "MSMS-Public_experimentspectra-pos-VS19",
    "MSMS-Public_experimentspectra-neg-VS19", "Other", "N/A", "NA",
    "Unknown", "None", "", "CLASS1",
}


def _is_generic_classification(cls: str | None) -> bool:
    if not cls:
        return True
    return cls.strip() in _GENERIC_CLASSIFICATIONS or cls.startswith("MoNA-export")


def _quality_score(rec: dict[str, Any]) -> float:
    """Score record quality for dedup. Higher = better.
    Prioritizes: specific classification > SMILES > real name > formula > peaks.

    中药成分分类（Flavonoid, Alkaloid, Terpenoid 等）获得大幅加分，
    确保 dedup 时不会被 MoNA-export 等通用标签覆盖。
    """
    score = 0.0
    if rec.get("smiles"):
        score += 10.0
    if rec.get("inchikey"):
        score += 3.0
    if rec.get("formula"):
        score += 2.0
    name = rec.get("name", "")
    if name and not name.startswith(("hdf5_", "NCGC")):
        score += 1.0
    # 具体分类加分（Flavonoid, Alkaloid 等 > 通用 MoNA 标签）
    cls = rec.get("compound_class", "")
    if cls and not _is_generic_classification(cls):
        score += 15.0  # 比 SMILES 更高的权重，确保具体分类不被覆盖
    peaks = rec.get("peaks", [])
    if peaks:
        score += min(len(peaks), 50) * 0.01
    return score


def deduplicate(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen = {}
    for rec in records:
        key = _dedup_key(rec)
        score = _quality_score(rec)
        if key not in seen or score > seen[key][1]:
            seen[key] = (rec, score)
    return [v[0] for v in seen.values()]


# ─────────────────────────────────────────────────────────────────
# 构建安全与质量门
# ─────────────────────────────────────────────────────────────────
def _canonical_source_tier(value: Any) -> str:
    """Normalize provenance without upgrading an unknown tier to experimental."""
    raw = str(value or "").strip().lower()
    compact = re.sub(r"[\s_-]+", "", raw)
    if not raw:
        return "user_unverified"
    if raw in _KNOWN_SOURCE_TIERS:
        return raw
    if compact in {"insilico", "predicted", "computed", "theoretical", "cfmid"}:
        return "predicted"
    if compact in {"literaturepresenceonly", "presenceonly", "literaturetarget"}:
        return "literature_presence_only"
    if compact in {"quarantine", "quarantined", "rejected"}:
        return "quarantined"
    if compact in {"unverified", "userunverified", "legacyunresolved"}:
        return "user_unverified"
    # An explicit but unknown provenance label is not proof of an experimental
    # measurement.  Keep it searchable only after a curator maps it explicitly.
    return "user_unverified"


def _quality_gate_records(records: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Apply the shared per-record quality gate without retaining huge audits.

    ``assess_records`` intentionally retains browser-facing audit samples.  A
    full rebuild can contain millions of records, so the builder calls the same
    authoritative single-record gate and only accumulates exact counters.
    """
    from matching.library_quality import assess_record

    accepted: list[dict[str, Any]] = []
    source_tiers: Counter[str] = Counter()
    errors: Counter[str] = Counter()
    warnings_found: Counter[str] = Counter()
    excluded_non_spectra = 0
    rejected = 0

    for index, record in enumerate(records, start=1):
        prepared = dict(record)
        declared_tier = _canonical_source_tier(prepared.get("source_tier"))
        prepared["source_tier"] = declared_tier
        if declared_tier == "literature_presence_only":
            # A presence assertion remains non-spectral even if an upstream
            # converter filled every diagnostic ion with an arbitrary value.
            prepared["intensity_origin"] = "default"

        clean, audit = assess_record(prepared)
        tier = _canonical_source_tier(clean.get("source_tier"))
        clean["source_tier"] = tier
        source_tiers[tier] += 1
        errors.update(audit.get("errors") or [])
        warnings_found.update(audit.get("warnings") or [])

        if not audit.get("accepted") or tier == "quarantined":
            rejected += 1
        elif tier in _NON_SPECTRAL_SOURCE_TIERS:
            excluded_non_spectra += 1
        else:
            accepted.append(clean)

        if index % 100_000 == 0:
            logger.info("  质量门 [%d / %d]", index, len(records))

    summary = {
        "total": len(records),
        "accepted_spectra": len(accepted),
        "rejected": rejected,
        "excluded_non_spectra": excluded_non_spectra,
        "source_tier_counts": dict(source_tiers),
        "error_counts": dict(errors),
        "warning_counts": dict(warnings_found),
    }
    return accepted, summary


def _compound_merge_key(rec: dict[str, Any]) -> tuple[str, ...]:
    """Return a chemically safe Compound identity for a fresh rebuild.

    A complete InChIKey is the strongest key.  Without one, merging requires
    both normalized name and molecular formula.  Name-only rows deliberately
    get a spectrum-specific identity so homonyms cannot share structures or
    classifications through the Compounds table.
    """
    inchikey = str(rec.get("inchikey") or "").strip().upper()
    if inchikey:
        return ("inchikey", inchikey)
    name = re.sub(r"\s+", " ", str(rec.get("name") or "").strip()).casefold()
    formula = re.sub(r"\s+", "", str(rec.get("formula") or "").strip())
    if name and formula:
        return ("name_formula", name, formula)
    return ("spectrum", str(rec.get("spectrum_id") or uuid.uuid4().hex))


def _temporary_sibling(path: Path, purpose: str) -> Path:
    """Create a unique *name* in the target directory (the file is not opened)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.parent / f".{path.name}.{purpose}-{uuid.uuid4().hex}.tmp"


def _cleanup_sqlite_artifacts(path: Path) -> None:
    for candidate in (
        path,
        Path(f"{path}-wal"),
        Path(f"{path}-shm"),
        Path(f"{path}-journal"),
    ):
        try:
            candidate.unlink()
        except FileNotFoundError:
            pass


def _fsync_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


_CORE_TABLES = {"Compounds", "Spectra", "Mapping"}
_CORE_BASE_COLUMNS = {
    "Compounds": {
        "id", "compound_name", "molecular_formula", "smiles", "inchi",
        "inchikey", "exact_mass", "classification", "compound_class",
    },
    "Spectra": {
        "id", "spectrum_id", "source_file", "dataset_source", "source_tier",
        "intensity_origin", "quality_flags_json", "precursor_mz",
        "retention_time", "collision_energy", "collision_energy_raw",
        "collision_energy_unit", "collision_energy_mode",
        "collision_energy_min", "collision_energy_max", "instrument_type",
        "instrument", "ion_mode", "source_accession", "charge", "adduct",
        "peaks_json",
    },
    "Mapping": {"id", "compound_id", "spectrum_id", "score"},
}


def _quote_identifier(value: str) -> str:
    return '"' + str(value).replace('"', '""') + '"'


def _pragma_table_info(conn: sqlite3.Connection, schema: str, table: str) -> list[sqlite3.Row]:
    schema_name = _quote_identifier(schema)
    table_name = _quote_identifier(table)
    return list(conn.execute(f"PRAGMA {schema_name}.table_info({table_name})"))


def _add_extension_column(
    conn: sqlite3.Connection,
    table: str,
    column_info: sqlite3.Row,
) -> None:
    _cid, name, column_type, not_null, default_value, primary_key = column_info[:6]
    if primary_key:
        raise RuntimeError(
            f"cannot preserve unexpected primary-key extension {table}.{name}"
        )
    definition = f"{_quote_identifier(name)} {column_type or ''}".rstrip()
    if not_null:
        if default_value is None:
            raise RuntimeError(
                f"cannot add NOT NULL extension without default: {table}.{name}"
            )
        definition += " NOT NULL"
    if default_value is not None:
        definition += f" DEFAULT {default_value}"
    conn.execute(
        f"ALTER TABLE {_quote_identifier(table)} ADD COLUMN {definition}"
    )


def _preserve_existing_schema_and_metadata(
    existing_db: Path,
    staging_db: Path,
) -> dict[str, Any]:
    """Copy non-core data and safe core extensions into the staged rebuild.

    Core spectra/compound rows are rebuilt from source files.  Every other
    table is independent application state (notably literature targets and
    migration audit history) and is copied byte-for-byte at the SQL row level.
    Core extension values are transferred only through stable identities:
    complete exact InChIKey for Compounds and exact spectrum_id for Spectra.
    """
    stats: dict[str, Any] = {
        "non_core_tables": {},
        "core_extension_columns": {},
        "indexes": [],
        "transferred_non_null": {},
    }
    previous_uri = f"file:{existing_db}?mode=ro"
    with sqlite3.connect(staging_db, uri=True) as conn:
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys=OFF")
        conn.execute("ATTACH DATABASE ? AS previous", (previous_uri,))

        old_tables = list(conn.execute("""
            SELECT name, sql FROM previous.sqlite_master
            WHERE type = 'table' AND name NOT LIKE 'sqlite_%' AND sql IS NOT NULL
            ORDER BY name
        """))
        for row in old_tables:
            table = str(row["name"])
            if table in _CORE_TABLES:
                continue
            conn.execute(str(row["sql"]))
            quoted = _quote_identifier(table)
            conn.execute(
                f"INSERT INTO {quoted} SELECT * FROM previous.{quoted}"
            )
            stats["non_core_tables"][table] = int(
                conn.execute(f"SELECT COUNT(*) FROM {quoted}").fetchone()[0]
            )

        # Add every old extension column that the current builder does not yet
        # know.  This prevents a future migration from disappearing merely
        # because the raw spectral library was rebuilt.
        for table in sorted(_CORE_TABLES):
            old_info = _pragma_table_info(conn, "previous", table)
            if not old_info:
                continue
            new_info = _pragma_table_info(conn, "main", table)
            new_names = {str(row[1]) for row in new_info}
            for column in old_info:
                if str(column[1]) not in new_names:
                    _add_extension_column(conn, table, column)
                    new_names.add(str(column[1]))

            extension_names = [
                str(column[1])
                for column in old_info
                if str(column[1]) not in _CORE_BASE_COLUMNS[table]
            ]
            stats["core_extension_columns"][table] = extension_names

            old_names = {str(row[1]) for row in old_info}
            if table == "Compounds" and "inchikey" in old_names:
                for column in extension_names:
                    quoted_column = _quote_identifier(column)
                    conn.execute(f"""
                        UPDATE Compounds AS n
                        SET {quoted_column} = (
                            SELECT o.{quoted_column}
                            FROM previous.Compounds AS o
                            WHERE o.inchikey = n.inchikey
                              AND o.{quoted_column} IS NOT NULL
                            ORDER BY o.id LIMIT 1
                        )
                        WHERE n.inchikey IS NOT NULL
                          AND LENGTH(n.inchikey) = 27
                          AND SUBSTR(n.inchikey, 15, 1) = '-'
                          AND SUBSTR(n.inchikey, 26, 1) = '-'
                          AND EXISTS (
                            SELECT 1 FROM previous.Compounds AS o
                            WHERE o.inchikey = n.inchikey
                              AND o.{quoted_column} IS NOT NULL
                          )
                    """)
                    stats["transferred_non_null"][f"{table}.{column}"] = int(
                        conn.execute(
                            f"SELECT COUNT(*) FROM Compounds WHERE {quoted_column} IS NOT NULL"
                        ).fetchone()[0]
                    )
            elif table == "Spectra" and "spectrum_id" in old_names:
                for column in extension_names:
                    quoted_column = _quote_identifier(column)
                    conn.execute(f"""
                        UPDATE Spectra AS n
                        SET {quoted_column} = (
                            SELECT o.{quoted_column}
                            FROM previous.Spectra AS o
                            WHERE o.spectrum_id = n.spectrum_id
                              AND o.{quoted_column} IS NOT NULL
                            ORDER BY o.id LIMIT 1
                        )
                        WHERE EXISTS (
                            SELECT 1 FROM previous.Spectra AS o
                            WHERE o.spectrum_id = n.spectrum_id
                              AND o.{quoted_column} IS NOT NULL
                        )
                    """)
                    stats["transferred_non_null"][f"{table}.{column}"] = int(
                        conn.execute(
                            f"SELECT COUNT(*) FROM Spectra WHERE {quoted_column} IS NOT NULL"
                        ).fetchone()[0]
                    )

        # Preserve all named indexes, including those on core extension columns
        # and independent target tables.  SQLite autoindexes have sql=NULL and
        # are recreated by the original CREATE TABLE statements.
        existing_indexes = {
            str(row[0])
            for row in conn.execute(
                "SELECT name FROM main.sqlite_master WHERE type='index'"
            )
        }
        for row in conn.execute("""
            SELECT name, sql FROM previous.sqlite_master
            WHERE type = 'index' AND sql IS NOT NULL
            ORDER BY name
        """):
            name = str(row["name"])
            if name in existing_indexes:
                continue
            conn.execute(str(row["sql"]))
            existing_indexes.add(name)
            stats["indexes"].append(name)

        conn.commit()
        foreign_key_errors = list(conn.execute("PRAGMA foreign_key_check"))
        if foreign_key_errors:
            raise RuntimeError(
                f"preserved schema has foreign-key errors: {foreign_key_errors[:5]}"
            )
        conn.execute("DETACH DATABASE previous")

    return stats


def _audit_staging_database(
    db_path: Path,
    mapping_path: Path,
    *,
    expected_compounds: int,
    expected_spectra: int,
    expected_tiers: Counter[str],
    preservation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run integrity, cardinality, provenance, and mapping audits before swap."""
    uri = f"file:{db_path}?mode=ro"
    with sqlite3.connect(uri, uri=True) as conn:
        quick_check = [str(row[0]) for row in conn.execute("PRAGMA quick_check")]
        if quick_check != ["ok"]:
            raise RuntimeError(f"staging SQLite quick_check failed: {quick_check[:5]}")

        counts = {
            table: int(conn.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0])
            for table in ("Compounds", "Spectra", "Mapping")
        }
        if counts["Compounds"] != expected_compounds:
            raise RuntimeError(
                f"staging Compound cardinality mismatch: {counts['Compounds']} != {expected_compounds}"
            )
        if counts["Spectra"] != expected_spectra:
            raise RuntimeError(
                f"staging Spectrum cardinality mismatch: {counts['Spectra']} != {expected_spectra}"
            )
        if counts["Mapping"] != expected_spectra:
            raise RuntimeError(
                f"staging Mapping cardinality mismatch: {counts['Mapping']} != {expected_spectra}"
            )

        distinct_spectra = int(conn.execute(
            "SELECT COUNT(DISTINCT spectrum_id) FROM Spectra"
        ).fetchone()[0])
        if distinct_spectra != expected_spectra:
            raise RuntimeError("staging contains duplicate spectrum_id values")

        orphan_mappings = int(conn.execute("""
            SELECT COUNT(*) FROM Mapping m
            LEFT JOIN Compounds c ON c.id = m.compound_id
            LEFT JOIN Spectra s ON s.id = m.spectrum_id
            WHERE c.id IS NULL OR s.id IS NULL
        """).fetchone()[0])
        invalid_mapping_counts = int(conn.execute("""
            SELECT COUNT(*) FROM (
                SELECT s.id FROM Spectra s
                LEFT JOIN Mapping m ON m.spectrum_id = s.id
                GROUP BY s.id HAVING COUNT(m.id) != 1
            )
        """).fetchone()[0])
        unmapped_compounds = int(conn.execute("""
            SELECT COUNT(*) FROM Compounds c
            LEFT JOIN Mapping m ON m.compound_id = c.id
            WHERE m.id IS NULL
        """).fetchone()[0])
        if orphan_mappings or invalid_mapping_counts or unmapped_compounds:
            raise RuntimeError(
                "staging mapping audit failed: "
                f"orphans={orphan_mappings}, bad_spectrum_mappings={invalid_mapping_counts}, "
                f"unmapped_compounds={unmapped_compounds}"
            )

        tier_counts = Counter({
            str(tier): int(count)
            for tier, count in conn.execute(
                "SELECT source_tier, COUNT(*) FROM Spectra GROUP BY source_tier"
            )
        })
        if tier_counts != expected_tiers:
            raise RuntimeError(
                f"staging source-tier mismatch: {dict(tier_counts)} != {dict(expected_tiers)}"
            )
        forbidden = int(conn.execute(
            "SELECT COUNT(*) FROM Spectra WHERE source_tier IN (?, ?)",
            tuple(sorted(_NON_SPECTRAL_SOURCE_TIERS)),
        ).fetchone()[0])
        if forbidden:
            raise RuntimeError(f"staging contains {forbidden} non-spectral/quarantined rows")
        malformed_quality = int(conn.execute("""
            SELECT COUNT(*) FROM Spectra
            WHERE source_tier IS NULL OR TRIM(source_tier) = ''
               OR quality_flags_json IS NULL OR NOT json_valid(quality_flags_json)
        """).fetchone()[0])
        if malformed_quality:
            raise RuntimeError(f"staging contains {malformed_quality} unaudited spectra")

        preservation = preservation or {}
        for table, expected_count in preservation.get("non_core_tables", {}).items():
            quoted = _quote_identifier(table)
            actual_count = int(
                conn.execute(f"SELECT COUNT(*) FROM {quoted}").fetchone()[0]
            )
            if actual_count != int(expected_count):
                raise RuntimeError(
                    f"preserved table cardinality mismatch for {table}: "
                    f"{actual_count} != {expected_count}"
                )
        for table, expected_columns in preservation.get(
            "core_extension_columns", {}
        ).items():
            actual_columns = {
                str(row[1]) for row in conn.execute(
                    f"PRAGMA table_info({_quote_identifier(table)})"
                )
            }
            missing = set(expected_columns) - actual_columns
            if missing:
                raise RuntimeError(
                    f"preserved extension columns missing from {table}: {sorted(missing)}"
                )
        actual_indexes = {
            str(row[0])
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='index'"
            )
        }
        missing_indexes = set(preservation.get("indexes", [])) - actual_indexes
        if missing_indexes:
            raise RuntimeError(
                f"preserved indexes missing: {sorted(missing_indexes)}"
            )
        foreign_key_errors = list(conn.execute("PRAGMA foreign_key_check"))
        if foreign_key_errors:
            raise RuntimeError(
                f"staging foreign-key audit failed: {foreign_key_errors[:5]}"
            )

    with mapping_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected_header = {
            "compound_id", "spectrum_row_id", "spectrum_id",
            "compound_name", "inchikey", "precursor_mz",
        }
        if set(reader.fieldnames or []) != expected_header:
            raise RuntimeError("staging mapping CSV header is invalid")
        mapping_rows = sum(1 for _ in reader)
    if mapping_rows != expected_spectra:
        raise RuntimeError(
            f"staging mapping CSV cardinality mismatch: {mapping_rows} != {expected_spectra}"
        )

    return {
        "quick_check": "ok",
        **counts,
        "source_tier_counts": dict(tier_counts),
        "mapping_csv_rows": mapping_rows,
    }


def _backup_existing_database(db_path: Path) -> Path:
    """Create and validate a consistent SQLite snapshot of the old database."""
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    backup_path = db_path.with_name(f"{db_path.name}.prebuild_{timestamp}.bak")
    suffix = 1
    while backup_path.exists():
        backup_path = db_path.with_name(
            f"{db_path.name}.prebuild_{timestamp}_{suffix}.bak"
        )
        suffix += 1
    staging_backup = _temporary_sibling(backup_path, "backup")
    try:
        source_uri = f"file:{db_path}?mode=ro"
        with sqlite3.connect(source_uri, uri=True) as source, sqlite3.connect(staging_backup) as target:
            source.backup(target, pages=16_384, sleep=0.01)
        with sqlite3.connect(f"file:{staging_backup}?mode=ro", uri=True) as conn:
            result = [str(row[0]) for row in conn.execute("PRAGMA quick_check")]
        if result != ["ok"]:
            raise RuntimeError(f"old-database backup quick_check failed: {result[:5]}")
        _fsync_file(staging_backup)
        os.replace(staging_backup, backup_path)
        _fsync_directory(backup_path.parent)
        return backup_path
    except Exception:
        _cleanup_sqlite_artifacts(staging_backup)
        raise


# ─────────────────────────────────────────────────────────────────
# 数据库写入
# ─────────────────────────────────────────────────────────────────
def _write_staging_database(
    records: list[dict[str, Any]],
    db_path: Path,
    mapping_path: Path,
) -> tuple[int, int, Counter[str]]:
    """Write a complete new database at a unique staging path."""
    import sqlalchemy
    from sqlalchemy.orm import Session, declarative_base, relationship

    Base = declarative_base()

    class Compound(Base):
        __tablename__ = "Compounds"
        id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
        compound_name = sqlalchemy.Column(sqlalchemy.String, index=True)
        molecular_formula = sqlalchemy.Column(sqlalchemy.String)
        smiles = sqlalchemy.Column(sqlalchemy.String)
        inchi = sqlalchemy.Column(sqlalchemy.String)
        inchikey = sqlalchemy.Column(sqlalchemy.String, index=True)
        exact_mass = sqlalchemy.Column(sqlalchemy.Float)
        classification = sqlalchemy.Column(sqlalchemy.String, nullable=True)
        compound_class = sqlalchemy.Column(sqlalchemy.String, nullable=True)
        # Extension columns used by the classifier/pathway enrichment pipeline.
        # Unknown future columns from an existing formal DB are added and
        # preserved generically before the atomic swap (see below).
        classification_v3 = sqlalchemy.Column(sqlalchemy.Text, nullable=True)
        classification_source = sqlalchemy.Column(sqlalchemy.Text, nullable=True)
        classification_confidence = sqlalchemy.Column(sqlalchemy.Float, nullable=True)
        fragmentation_pathway = sqlalchemy.Column(sqlalchemy.Text, nullable=True)

    class SpectrumRecord(Base):
        __tablename__ = "Spectra"
        id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
        spectrum_id = sqlalchemy.Column(sqlalchemy.String, unique=True, index=True)
        source_file = sqlalchemy.Column(sqlalchemy.String)
        dataset_source = sqlalchemy.Column(sqlalchemy.String)  # 数据集名[-极性], 由 normalize_source_file 派生
        source_tier = sqlalchemy.Column(sqlalchemy.String)
        intensity_origin = sqlalchemy.Column(sqlalchemy.String)
        quality_flags_json = sqlalchemy.Column(sqlalchemy.String)
        precursor_mz = sqlalchemy.Column(sqlalchemy.Float, index=True)
        retention_time = sqlalchemy.Column(sqlalchemy.Float, nullable=True)
        collision_energy = sqlalchemy.Column(sqlalchemy.Float, nullable=True)
        collision_energy_raw = sqlalchemy.Column(sqlalchemy.Text, nullable=True)
        collision_energy_unit = sqlalchemy.Column(sqlalchemy.String, nullable=True)
        collision_energy_mode = sqlalchemy.Column(sqlalchemy.String, nullable=True)
        collision_energy_min = sqlalchemy.Column(sqlalchemy.Float, nullable=True)
        collision_energy_max = sqlalchemy.Column(sqlalchemy.Float, nullable=True)
        instrument_type = sqlalchemy.Column(sqlalchemy.Text, nullable=True)
        instrument = sqlalchemy.Column(sqlalchemy.Text, nullable=True)
        ion_mode = sqlalchemy.Column(sqlalchemy.String, nullable=True)
        source_accession = sqlalchemy.Column(sqlalchemy.Text, nullable=True)
        charge = sqlalchemy.Column(sqlalchemy.Integer, nullable=True)
        adduct = sqlalchemy.Column(sqlalchemy.String, nullable=True)
        peaks_json = sqlalchemy.Column(sqlalchemy.String)

    class Mapping(Base):
        __tablename__ = "Mapping"
        id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
        compound_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey("Compounds.id"), index=True)
        spectrum_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey("Spectra.id"), index=True)
        score = sqlalchemy.Column(sqlalchemy.Float, nullable=True)

    db_path.parent.mkdir(parents=True, exist_ok=True)
    engine = sqlalchemy.create_engine(f"sqlite:///{db_path}", echo=False)
    Base.metadata.create_all(engine)

    with engine.connect() as conn:
        from sqlalchemy import text
        # Staging is deliberately built in rollback-journal mode.  A detached
        # WAL file must never be paired with the new main database at swap time.
        conn.execute(text("PRAGMA journal_mode=DELETE"))
        conn.execute(text("PRAGMA synchronous=FULL"))
        conn.execute(text("CREATE INDEX IF NOT EXISTS ix_Compounds_inchikey ON Compounds(inchikey)"))
        conn.execute(text("CREATE INDEX IF NOT EXISTS ix_Compounds_compound_name ON Compounds(compound_name)"))
        conn.execute(text("CREATE INDEX IF NOT EXISTS ix_Spectra_spectrum_id ON Spectra(spectrum_id)"))
        conn.execute(text("CREATE INDEX IF NOT EXISTS ix_Spectra_precursor_mz ON Spectra(precursor_mz)"))
        conn.execute(text("CREATE INDEX IF NOT EXISTS ix_Spectra_source_tier ON Spectra(source_tier)"))
        conn.execute(text("CREATE INDEX IF NOT EXISTS ix_Mapping_compound_id ON Mapping(compound_id)"))
        conn.commit()

    all_mapping = []
    compound_cache = {}
    total_compounds = 0
    total_spectra = 0
    t0 = time.time()

    try:
        for batch_start in range(0, len(records), BATCH_SIZE):
            batch = records[batch_start:batch_start + BATCH_SIZE]
            with Session(engine) as session:
                for rec in batch:
                    # ── 插入 Compound ──────────────────────────
                    comp_key = _compound_merge_key(rec)
                    compound_id = compound_cache.get(comp_key)

                    if compound_id is None and rec.get("inchikey"):
                        row = session.execute(
                            sqlalchemy.text(
                                "SELECT id, smiles, molecular_formula FROM Compounds WHERE inchikey = :ik"
                            ),
                            {"ik": rec["inchikey"]},
                        ).fetchone()
                        if row:
                            compound_id = row[0]
                            compound_cache[comp_key] = compound_id
                            # Backfill missing fields only under the exact full
                            # InChIKey identity already chosen above.
                            if rec.get("smiles") and not row[1]:
                                session.execute(
                                    sqlalchemy.text(
                                        "UPDATE Compounds SET smiles = :s, inchi = :i, "
                                        "molecular_formula = :f WHERE id = :id"
                                    ),
                                    {
                                        "s": rec.get("smiles"),
                                        "i": rec.get("inchi"),
                                        "f": rec.get("formula"),
                                        "id": compound_id,
                                    },
                                )
                            if rec.get("formula") and not row[2]:
                                session.execute(
                                    sqlalchemy.text(
                                        "UPDATE Compounds SET molecular_formula = :f WHERE id = :id"
                                    ),
                                    {"f": rec.get("formula"), "id": compound_id},
                                )

                    if compound_id is None:
                        compound = Compound(
                            compound_name=rec.get("name"),
                            molecular_formula=rec.get("formula"),
                            smiles=rec.get("smiles"),
                            inchi=rec.get("inchi"),
                            inchikey=rec.get("inchikey"),
                            exact_mass=rec.get("exact_mass"),
                            classification=rec.get("classification") or rec.get("compound_class"),
                            compound_class=rec.get("compound_class"),
                        )
                        session.add(compound)
                        session.flush()
                        compound_id = compound.id
                        compound_cache[comp_key] = compound_id
                        total_compounds += 1

                    # ── 插入 Spectrum ──────────────────────────
                    from mssystem.classification import normalize_source_file as _norm_src
                    spectrum = SpectrumRecord(
                        spectrum_id=str(rec.get("spectrum_id")),
                        source_file=rec.get("source_file"),
                        dataset_source=_norm_src(rec.get("source_file")),
                        source_tier=rec.get("source_tier"),
                        intensity_origin=rec.get("intensity_origin") or "measured",
                        quality_flags_json=json.dumps(rec.get("quality_flags") or [], ensure_ascii=False),
                        precursor_mz=rec.get("precursor_mz"),
                        retention_time=rec.get("retention_time"),
                        collision_energy=rec.get("collision_energy"),
                        collision_energy_raw=rec.get("collision_energy_raw"),
                        collision_energy_unit=rec.get("collision_energy_unit"),
                        collision_energy_mode=rec.get("collision_energy_mode"),
                        collision_energy_min=rec.get("collision_energy_min"),
                        collision_energy_max=rec.get("collision_energy_max"),
                        instrument_type=rec.get("instrument_type"),
                        instrument=rec.get("instrument"),
                        ion_mode=rec.get("ion_mode"),
                        source_accession=rec.get("source_accession"),
                        charge=rec.get("charge"),
                        adduct=rec.get("adduct"),
                        peaks_json=json.dumps(rec.get("peaks", []), ensure_ascii=False),
                    )
                    session.add(spectrum)
                    session.flush()
                    total_spectra += 1

                    # ── 插入 Mapping ───────────────────────────
                    session.add(Mapping(compound_id=compound_id, spectrum_id=spectrum.id))
                    all_mapping.append({
                        "compound_id": compound_id,
                        "spectrum_row_id": spectrum.id,
                        "spectrum_id": spectrum.spectrum_id,
                        "compound_name": rec.get("name"),
                        "inchikey": rec.get("inchikey"),
                        "precursor_mz": rec.get("precursor_mz"),
                    })

                session.commit()

            elapsed = time.time() - t0
            rate = (batch_start + len(batch)) / elapsed if elapsed > 0 else 0
            logger.info(
                "  [%d / %d] (%.0f rec/sec)  compounds=%d  spectra=%d",
                min(batch_start + BATCH_SIZE, len(records)), len(records),
                rate, total_compounds, total_spectra,
            )
    finally:
        engine.dispose()

    # 导出映射 CSV
    mapping_path.parent.mkdir(parents=True, exist_ok=True)
    with mapping_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "compound_id", "spectrum_row_id", "spectrum_id",
            "compound_name", "inchikey", "precursor_mz",
        ])
        writer.writeheader()
        writer.writerows(all_mapping)

    source_tiers = Counter(str(rec.get("source_tier")) for rec in records)
    return total_compounds, total_spectra, source_tiers


def build_database(records: list[dict], db_path: Path, mapping_path: Path) -> tuple[int, int]:
    """Serialize a build against V6 incremental imports, then rebuild atomically."""
    target = Path(db_path).resolve(strict=False)
    shared_targets = {
        SHARED_DB_PATH.resolve(strict=False),
        SHARED_DB_PATH_SMALL.resolve(strict=False),
    }
    lock_path = (
        USER_CONTRIB_DIR / ".library_import.lock"
        if target in shared_targets
        else Path(db_path).parent / ".library_import.lock"
    )
    with exclusive_file_lock(lock_path):
        return _build_database_locked(records, db_path, mapping_path)


def _build_database_locked(records: list[dict], db_path: Path, mapping_path: Path) -> tuple[int, int]:
    """Build, validate, back up, and atomically replace a formal SQLite library."""
    db_path = Path(db_path)
    mapping_path = Path(mapping_path)
    shared_targets = {
        SHARED_DB_PATH.resolve(strict=False),
        SHARED_DB_PATH_SMALL.resolve(strict=False),
    }
    if db_path.resolve(strict=False) in shared_targets and not LIBRARY_WRITE_ENABLED:
        raise PermissionError("shared library writes are disabled for this instance")
    require_storage_capacity(
        db_path.parent,
        estimate_sqlite_rebuild_bytes(db_path, mapping_path),
        f"atomic rebuild of {db_path.name}",
    )
    db_path.parent.mkdir(parents=True, exist_ok=True)
    mapping_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("应用谱库质量门（%d 条记录）...", len(records))
    quality_records, quality_summary = _quality_gate_records(records)
    logger.info(
        "质量门结果: spectral=%d rejected=%d presence-only=%d tiers=%s",
        quality_summary["accepted_spectra"],
        quality_summary["rejected"],
        quality_summary["excluded_non_spectra"],
        quality_summary["source_tier_counts"],
    )
    if not quality_records:
        raise ValueError("质量门后没有可构建的谱图记录；正式数据库保持不变")

    logger.info("去重中（%d 条质量合格记录）...", len(quality_records))
    deduped = deduplicate(quality_records)
    logger.info("去重后：%d 条记录", len(deduped))
    if not deduped:
        raise ValueError("去重后没有谱图记录；正式数据库保持不变")

    staging_db = _temporary_sibling(db_path, "staging")
    staging_mapping = _temporary_sibling(mapping_path, "staging")
    mapping_rollback: Path | None = None
    backup_path: Path | None = None
    database_replaced = False

    try:
        compounds, spectra, source_tiers = _write_staging_database(
            deduped, staging_db, staging_mapping
        )
        preservation: dict[str, Any] = {}
        if db_path.exists():
            preservation = _preserve_existing_schema_and_metadata(
                db_path, staging_db
            )
            logger.info("已保留现库扩展数据: %s", preservation)
        audit = _audit_staging_database(
            staging_db,
            staging_mapping,
            expected_compounds=compounds,
            expected_spectra=spectra,
            expected_tiers=source_tiers,
            preservation=preservation,
        )
        logger.info("staging 审计通过: %s", audit)
        _fsync_file(staging_db)
        _fsync_file(staging_mapping)

        if db_path.exists():
            backup_path = _backup_existing_database(db_path)
            logger.info("旧数据库一致性备份: %s", backup_path)

        # The CSV is supplementary, but keep its visible version paired with
        # the database.  If the database swap fails, restore the old CSV before
        # propagating the exception; the old database itself is never touched.
        if mapping_path.exists():
            mapping_rollback = _temporary_sibling(mapping_path, "rollback")
            os.replace(mapping_path, mapping_rollback)
        os.replace(staging_mapping, mapping_path)
        _fsync_directory(mapping_path.parent)
        try:
            # Same-directory os.replace is the only operation that changes the
            # formal database path, so readers see either the complete old DB
            # or the complete audited new DB—never a partially populated file.
            os.replace(staging_db, db_path)
            database_replaced = True
        except Exception:
            try:
                mapping_path.unlink()
            except FileNotFoundError:
                pass
            if mapping_rollback and mapping_rollback.exists():
                os.replace(mapping_rollback, mapping_path)
                mapping_rollback = None
            raise

        # The database replacement has succeeded.  Durability housekeeping and
        # removal of the now-obsolete CSV rollback file must not turn this into
        # a reported failure (which would falsely imply the old DB is active).
        try:
            _fsync_directory(db_path.parent)
        except OSError as exc:
            logger.warning("数据库已替换，但目录 fsync 失败: %s", exc)
        if mapping_rollback and mapping_rollback.exists():
            try:
                mapping_rollback.unlink()
            except OSError as exc:
                logger.warning("数据库已替换；旧 mapping 回滚文件保留: %s", exc)
            mapping_rollback = None

        logger.info(
            "数据库原子替换完成: %s (%d compounds, %d spectra)%s",
            db_path,
            compounds,
            spectra,
            f"; backup={backup_path}" if backup_path else "",
        )
        removed_backups = prune_generated_backups(db_path, keep=3)
        if removed_backups:
            logger.info(
                "备份保留策略已清理 %d 个旧 prebuild 快照: %s",
                len(removed_backups),
                ", ".join(path.name for path in removed_backups),
            )
        return compounds, spectra
    except Exception:
        if database_replaced:
            # There is no rollback-by-rewrite after the atomic boundary: the
            # audited new database is already the formal database and the old
            # one is safely retained at backup_path.
            logger.exception(
                "数据库已完成原子替换；替换后的非关键收尾步骤失败，保留新库与旧库备份"
            )
            return compounds, spectra
        # Every failure before the final os.replace leaves db_path untouched.
        # Unique staging siblings and SQLite sidecars are always removed.
        _cleanup_sqlite_artifacts(staging_db)
        try:
            staging_mapping.unlink()
        except FileNotFoundError:
            pass
        if mapping_rollback and mapping_rollback.exists():
            try:
                mapping_path.unlink()
            except FileNotFoundError:
                pass
            os.replace(mapping_rollback, mapping_path)
        raise


# ─────────────────────────────────────────────────────────────────
# 辅助函数
# ─────────────────────────────────────────────────────────────────
def _coerce_float(val):
    if val is None:
        return None
    if isinstance(val, (int, float)):
        return float(val) if not (val != val) else None  # check NaN
    try:
        return float(str(val).strip())
    except (TypeError, ValueError):
        return None


def _coerce_int(val):
    if val is None:
        return None
    text = str(val).strip()
    # Common MGF forms are ``1+``, ``1-``, ``+1`` and ``-1``.  Removing both
    # signs (the historical behaviour) silently converted every negative
    # charge to positive and defeated polarity fallback when adduct is absent.
    match = re.match(r"^([+-])?\s*(\d+(?:\.0*)?)\s*([+-])?", text)
    if not match:
        return None
    try:
        magnitude = int(float(match.group(2)))
    except ValueError:
        return None
    leading, trailing = match.group(1), match.group(3)
    if leading == "-" or trailing == "-":
        return -magnitude
    return magnitude


# ─────────────────────────────────────────────────────────────────
# 主函数
# ─────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="从 MGF/MSP 构建质谱数据库（跳过 HDF5）")
    parser.add_argument("--mgf-dir", type=Path,
                        default=DEFAULT_LIBRARY_ROOT,
                        help="MGF/MSP 根目录")
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR,
                        help="输出目录")
    parser.add_argument("--db-name", default="local_library.db")
    parser.add_argument("--mapping-name", default="library_mapping.csv")
    args = parser.parse_args()

    db_path = args.output / args.db_name
    mapping_path = args.output / args.mapping_name

    all_records = []
    total_files = 0
    total_records = 0

    # ── 解析所有 MGF 文件 ──────────────────────────────────────
    mgf_dir = args.mgf_dir.expanduser().resolve()
    if not mgf_dir.is_dir():
        parser.error(f"MGF/MSP 源目录不存在或不是目录: {mgf_dir}")
    mgf_files = sorted(mgf_dir.rglob("*.mgf"))
    logger.info("扫描 MGF 文件: %d 个", len(mgf_files))

    # Prioritize files with rich metadata (Pluskal Lab, classify data)
    # so dedup keeps their SMILES-rich records over bare MoNA records
    def _file_priority(fp):
        name = fp.name.lower()
        # Highest: Pluskal Lab files (100% SMILES coverage)
        if any(p in name for p in ["nihnp", "mcescaf", "otavapep", "mcebio"]):
            return 0
        # High: classify data (usually has SMILES in COMMENT)
        if "classify_data" in str(fp):
            return 1
        # High: VF-NPL (rich metadata)
        if name.startswith("vf-"):
            return 2
        # Medium: other small files
        return 5

    mgf_files.sort(key=_file_priority)

    global_idx = 0
    for fp in mgf_files:
        records, global_idx = parse_mgf_file(fp, global_idx, source_root=mgf_dir)
        if records:
            all_records.extend(records)
            logger.info("  MGF %s: %d 条", fp.name, len(records))
        total_files += 1

    # ── 解析所有 MSP 文件 ─────────────────────────────────────
    msp_files = sorted(mgf_dir.rglob("*.msp"))
    logger.info("扫描 MSP 文件: %d 个", len(msp_files))

    def _msp_priority(fp):
        name = fp.name.lower()
        # Vaniya-Fiehn and bmdms have good SMILES coverage
        if any(p in name for p in ["vaniya", "bmdms"]):
            return 0
        # FiehnHILIC, GNPS, MassBank have moderate coverage
        if any(p in name for p in ["fiehn", "gnps", "massbank"]):
            return 1
        # MoNA-export-All_Spectra is huge but low SMILES -> process last
        if "all_spectra" in name:
            return 9
        return 3

    msp_files.sort(key=_msp_priority)

    for fp in msp_files:
        records, global_idx = parse_msp_file(fp, global_idx, source_root=mgf_dir)
        if records:
            all_records.extend(records)
            logger.info("  MSP %s: %d 条", fp.name, len(records))
        total_files += 1

    if total_files == 0:
        logger.error("源目录中没有任何 .mgf/.msp 文件: %s", mgf_dir)
        sys.exit(2)
    if not all_records:
        logger.error("找到 %d 个源文件，但没有可解析谱图: %s", total_files, mgf_dir)
        sys.exit(1)

    # ── 构建数据库 ────────────────────────────────────────────
    n_compounds, n_spectra = build_database(all_records, db_path, mapping_path)

    # ── 统计 ─────────────────────────────────────────────────
    with_name = sum(1 for r in all_records if r.get("name") and not r["name"].startswith(
        tuple(["hdf5_", "20231", "20240", "NCGC", "mcebio", "mcescaf", "nihnp", "otavapep"])))
    with_smiles = sum(1 for r in all_records if r.get("smiles"))
    with_inchikey = sum(1 for r in all_records if r.get("inchikey"))

    logger.info("")
    logger.info("=" * 60)
    logger.info("构建完成！")
    logger.info("  数据库: %s", db_path)
    logger.info("  映射:   %s", mapping_path)
    logger.info("  解析文件: %d", total_files)
    logger.info("  总记录数: %d", len(all_records))
    logger.info("  去重后:   %d", n_spectra)
    logger.info("  化合物数: %d", n_compounds)
    logger.info("  有名称:   %d (%.1f%%)", with_name, with_name/len(all_records)*100)
    logger.info("  有SMILES: %d (%.1f%%)", with_smiles, with_smiles/len(all_records)*100)
    logger.info("  有InChIKey: %d (%.1f%%)", with_inchikey, with_inchikey/len(all_records)*100)
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
