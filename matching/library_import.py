"""
matching/library_import.py — 谱库扩增模块
─────────────────────────────────────────
将用户从文献搜集的 CSV/Excel 数据导入质谱库：
    CSV/Excel ──convert_single──► MGF ──parse_mgf_file──► records
    ──归一化/去重──► append_records ──► 追加到 Full / Small 双 SQLite 库

适配文献数据的现实（实测）：
- InChIKey 平均覆盖 <15%、SMILES 0–53%，但化合物名/formula/precursor_mz ≈100%
- formula 格式混乱：ASCII、Unicode 下标、元素乱序并存

核心设计：
- normalize_formula():    Hill 顺序归一化（C₅H₁₁O₂NS → C5H11NO2S）
- compound_match_key():   多级身份键（InChIKey14 → 名称+formula → 名称 → mz）
- append_records():       增量追加到现有库（不删库，复用 inchikey 去重复用 + 名称借用补全）
"""
from __future__ import annotations

import json
import hashlib
import csv
import fcntl
import logging
import os
import re
import sqlite3
import time
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Optional, Sequence

from config import settings as app_settings

logger = logging.getLogger("library_import")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
_MODULE_PROJECT_ROOT = PROJECT_ROOT
SMALL_MASS_THRESHOLD = 800.0  # <800 Da 进 Small 库（与 settings/记忆一致）


def _assert_library_write_enabled() -> None:
    if not app_settings.LIBRARY_WRITE_ENABLED:
        raise PermissionError(
            "This MSsystem instance is configured as a read-only shared-library reader"
        )


def _configured_data_path(filename: str, configured: Path) -> Path:
    """Resolve a shared-data path while preserving isolated test roots."""
    if PROJECT_ROOT.resolve(strict=False) != _MODULE_PROJECT_ROOT.resolve(strict=False):
        return PROJECT_ROOT / "data" / filename
    return Path(configured)


def _user_contrib_dir() -> Path:
    return _configured_data_path("user_contrib", app_settings.USER_CONTRIB_DIR)

# 常见元素精确质量（用于无 SMILES 时从 formula 估算 exact_mass，主要供 Small 库过滤）
_ATOMIC_MASS = {
    "C": 12.0, "H": 1.00782503223, "N": 14.00307400443,
    "O": 15.99491461957, "S": 31.9720711744, "P": 30.97376199842,
    "F": 18.99840316273, "Cl": 34.968852682, "Br": 78.9183376,
    "I": 126.904468, "Na": 22.9897692820, "K": 38.9637064864,
    "Ca": 39.962590863, "Mg": 23.985041697, "Fe": 55.93493633,
    "Si": 27.97692653465, "Se": 79.9165218, "B": 11.00930536,
}


# ════════════════════════════════════════════════════════════════════
# 归一化工具
# ════════════════════════════════════════════════════════════════════
def normalize_formula(f: Any) -> Optional[str]:
    """归一化分子式到 Hill 顺序标准式。

    处理：Unicode 下标（C₆H₁₂→C6H12）、元素乱序（C5H11O2NS→C5H11NO2S）、
    水合物（C6H12O6·H2O→C6H12O6）、空格。
    """
    if f is None:
        return None
    s = str(f).strip()
    if not s:
        return None
    # 1. Unicode 下标/上标归一化（复用 csv_excel_to_mgf）
    try:
        from csv_excel_to_mgf import normalize_unicode
        got = normalize_unicode(s)
        if got:
            s = got
    except Exception:
        pass
    # 2. 去水合物/盐前缀（取第一段）
    s = re.split(r"[·.]", s)[0]
    s = s.replace(" ", "").replace("\t", "")
    # 3. 解析元素 + 原子数
    counts: dict[str, int] = {}
    for elem, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", s):
        if not elem or not re.match(r"^[A-Z][a-z]?$", elem):
            continue
        counts[elem] = counts.get(elem, 0) + (int(cnt) if cnt else 1)
    if not counts:
        return None
    # 4. Hill 顺序：有 C 时 C 在前 H 第二；无 C 时全字母序
    parts = []
    if "C" in counts:
        parts.append(("C", counts.pop("C")))
        if "H" in counts:
            parts.append(("H", counts.pop("H")))
    for elem in sorted(counts):
        parts.append((elem, counts[elem]))
    return "".join(f"{e}{c if c > 1 else ''}" for e, c in parts)


def normalize_compound_name(name: Any) -> Optional[str]:
    """化合物名规范化：去尾部标记符号、折叠空白、转小写（大小写不敏感匹配）。"""
    if name is None:
        return None
    s = str(name).strip()
    if not s:
        return None
    s = re.sub(r"[*./_~]+$", "", s).strip()
    s = re.sub(r"\s+", " ", s)
    return s.lower() if s else None


def _is_real_inchikey(ik: Any) -> bool:
    """校验 InChIKey 合法性（避免误存的 InChI 串/脏数据充当去重键）。"""
    if not ik:
        return False
    s = str(ik).strip()
    if s.startswith("InChI="):
        return False
    # 标准 InChIKey: AAAAAAAAAAAAAA-AAAAAAAAAA-A  （或仅前14位连接层）
    return bool(re.match(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$", s)) or bool(re.match(r"^[A-Z]{14}$", s))


def _full_inchikey(value: Any) -> Optional[str]:
    text = str(value or "").strip().upper()
    return text if re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", text) else None


def _reported_full_inchikey(rec: dict) -> Optional[str]:
    """Return the submitted full key, whether verified or not.

    A reported key is useful as a *conflict guard* and for provenance.  It is
    intentionally not sufficient to establish structural identity.
    """
    return (_full_inchikey(rec.get("reported_inchikey"))
            or _full_inchikey(rec.get("inchikey")))


def _verified_inchikey(rec: dict) -> Optional[str]:
    """Return a full key only when it is backed by a parsed structure."""
    key = _full_inchikey(rec.get("inchikey"))
    if not key:
        return None
    if rec.get("structure_identity_verified") is True:
        return key
    # ``append_records`` is public and older callers may not have passed the
    # new audit marker.  Validate those records directly rather than trusting
    # a user-supplied boolean-like/string value.
    smiles = str(rec.get("smiles") or "").strip()
    inchi = str(rec.get("inchi") or "").strip()
    if not smiles and not inchi:
        return None
    try:
        from matching.library_quality import inchi_metadata, smiles_metadata
        metadata = smiles_metadata(smiles) if smiles else inchi_metadata(inchi)
        derived = _full_inchikey(metadata.get("inchikey"))
    except Exception:
        derived = None
    return key if derived == key else None


def compound_match_key(rec: dict) -> tuple[str, str]:
    """多级化合物身份键：(级别, key)。级别越小优先级越高。

    Level 1: 实验谱使用完整 InChIKey（27 位，保留立体化学层）
    Level 2: 规范化名称 + 归一化 formula   ← 文献数据的主力身份
    Level 3: 规范化名称
    Level 4: precursor_mz（4 位小数，仅兜底）
    """
    ik = _verified_inchikey(rec)
    reported_ik = _reported_full_inchikey(rec)
    source_tier = str(rec.get("source_tier") or "").strip().lower()
    # 文献 presence-only 行的结构尚未由实验谱确认；不能让误粘的同一 SMILES
    # 把不同报道名称合并到一个 compound。
    if ik and _is_real_inchikey(ik) and source_tier != "literature_presence_only":
        return ("ik", str(ik).strip().upper())
    name = normalize_compound_name(rec.get("name"))
    formula = normalize_formula(rec.get("formula"))
    if reported_ik and source_tier != "literature_presence_only":
        # Keep an unverified submitted key as a conflict barrier only.  The
        # name/formula suffix prevents the key from becoming a surrogate
        # structure identity during in-file deduplication.
        return (
            "reported-ik-guard",
            f"{reported_ik}::{name or ''}::{formula or ''}",
        )
    if name and formula:
        return ("name+formula", f"{name}::{formula}")
    if name:
        return ("name", name)
    pmz = rec.get("precursor_mz")
    if pmz:
        return ("mz", f"{float(pmz):.4f}")
    return ("sid", str(rec.get("spectrum_id", "unknown")))


def _exact_mass_from_smiles(smiles: Optional[str]) -> Optional[float]:
    if not smiles:
        return None
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import Descriptors
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
        if mol:
            return float(Descriptors.ExactMolWt(mol))
    except Exception:
        pass
    return None


def _exact_mass_from_formula(formula: Optional[str]) -> Optional[float]:
    """从 formula 用原子量估算精确质量（无 SMILES 时的 fallback）。"""
    nf = normalize_formula(formula)
    if not nf:
        return None
    counts: dict[str, int] = {}
    pos = 0
    for match in re.finditer(r"([A-Z][a-z]?)(\d*)", nf):
        if match.start() != pos:
            return None
        elem, cnt = match.groups()
        if elem not in _ATOMIC_MASS:
            return None
        counts[elem] = counts.get(elem, 0) + (int(cnt) if cnt else 1)
        pos = match.end()
    if pos != len(nf):
        return None
    total = sum(_ATOMIC_MASS[e] * n for e, n in counts.items())
    return total if total > 0 else None


# ════════════════════════════════════════════════════════════════════
# 现有库读取（只读，用于预览比对）
# ════════════════════════════════════════════════════════════════════
def _load_existing_keys(
    db_path: Path, records: Optional[list[dict]] = None
) -> tuple[set[str], set[str], dict[str, bool], dict[str, bool]]:
    """Return exact verified keys and normalized ``name+formula`` identities.

    Incoming records with any full reported InChIKey are only queried by a
    structure-verified exact key.  An exact miss never falls back to a
    name/formula query, because doing so can attach the spectrum to a different
    stereoisomer or constitutional isomer.

    同一化合物多条记录时，SMILES 状态取"有结构"（任一条有就算有）。
    """
    existing_ik: set[str] = set()
    existing_name: set[str] = set()
    ik_has_smi: dict[str, bool] = {}
    name_has_smi: dict[str, bool] = {}
    if not Path(db_path).exists():
        return existing_ik, existing_name, ik_has_smi, name_has_smi
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        if records is None:
            rows = conn.execute(
                "SELECT inchikey, compound_name, molecular_formula, smiles FROM Compounds"
            )
        else:
            # Preview files contain hundreds of rows while the production DB has
            # >1M compounds.  Query only incoming identities instead of loading
            # the entire Compounds table into Python memory.
            rows_list: list[tuple[Any, Any, Any]] = []
            seen_queries: set[tuple[str, str]] = set()
            for rec in records:
                ik = _verified_inchikey(rec)
                reported_ik = _reported_full_inchikey(rec)
                name = normalize_compound_name(rec.get("name"))
                formula = normalize_formula(rec.get("formula"))
                if ik and ("ik", ik) not in seen_queries:
                    seen_queries.add(("ik", ik))
                    rows_list.extend(conn.execute(
                        "SELECT inchikey, compound_name, molecular_formula, smiles FROM Compounds "
                        "WHERE inchikey=?", (ik,)
                    ).fetchall())
                # A full reported key, including an unverified one, is a hard
                # conflict guard.  Do not silently reinterpret it as the
                # different key attached to a same-name/formula DB row.
                if (not ik and not reported_ik and name and formula
                        and ("name_formula", f"{name}|{formula}") not in seen_queries):
                    seen_queries.add(("name_formula", f"{name}|{formula}"))
                    rows_list.extend(conn.execute(
                        "SELECT inchikey, compound_name, molecular_formula, smiles FROM Compounds "
                        "WHERE lower(trim(compound_name))=? AND molecular_formula=?",
                        (name, formula),
                    ).fetchall())
            rows = rows_list
        for ik, nm, formula, smi in rows:
            has = bool(smi and str(smi).strip())
            full_key = _full_inchikey(ik)
            if full_key:
                k = full_key
                existing_ik.add(k)
                ik_has_smi[k] = ik_has_smi.get(k, False) or has
            nn = normalize_compound_name(nm)
            nf = normalize_formula(formula)
            if nn and nf:
                name_formula = f"{nn}::{nf}"
                existing_name.add(name_formula)
                name_has_smi[name_formula] = name_has_smi.get(name_formula, False) or has
        conn.close()
    except Exception as e:
        logger.warning("读取现有库失败 %s: %s", db_path, e)
    return existing_ik, existing_name, ik_has_smi, name_has_smi


# ════════════════════════════════════════════════════════════════════
# 增量追加（核心）
# ════════════════════════════════════════════════════════════════════
_SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS Compounds (
    id INTEGER PRIMARY KEY, compound_name TEXT, molecular_formula TEXT,
    smiles TEXT, inchi TEXT, inchikey TEXT, exact_mass REAL,
    classification TEXT, compound_class TEXT,
    reported_inchikey TEXT, structure_identity_verified INTEGER DEFAULT 0,
    structure_identity_source TEXT
);
CREATE TABLE IF NOT EXISTS Spectra (
    id INTEGER PRIMARY KEY, spectrum_id TEXT UNIQUE, source_file TEXT,
    precursor_mz REAL, retention_time REAL, collision_energy REAL,
    collision_energy_raw TEXT, collision_energy_unit TEXT,
    collision_energy_mode TEXT, collision_energy_min REAL,
    collision_energy_max REAL,
    charge INTEGER, adduct TEXT, peaks_json TEXT,
    instrument_type TEXT, instrument TEXT, ion_mode TEXT,
    source_accession TEXT, activation_method TEXT,
    instrument_domain_id TEXT, collision_cell_region TEXT,
    acquisition_method TEXT, source_asset_sha256 TEXT,
    provenance_json TEXT, reference_material_type TEXT,
    standard_lot TEXT, standard_purity REAL,
    digitization_method TEXT, source_image_file TEXT,
    dataset_source TEXT, source_tier TEXT,
    intensity_origin TEXT, quality_flags_json TEXT
);
CREATE TABLE IF NOT EXISTS Mapping (
    id INTEGER PRIMARY KEY, compound_id INTEGER, spectrum_id INTEGER, score REAL
);
CREATE INDEX IF NOT EXISTS ix_c_name_formula_norm
    ON Compounds(lower(trim(compound_name)), molecular_formula);
"""


def _has_simple_index(
    conn: sqlite3.Connection, table: str, columns: tuple[str, ...]
) -> bool:
    """Return True when any existing index already covers exactly these columns."""
    for row in conn.execute(f'PRAGMA index_list("{table}")'):
        index_name = str(row[1]).replace('"', '""')
        indexed = tuple(
            str(info[2]) for info in conn.execute(f'PRAGMA index_info("{index_name}")')
        )
        if indexed == columns:
            return True
    return False


def _create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(_SCHEMA_SQL)
    # Historical builders used several index names for the same columns.  A
    # name-only ``IF NOT EXISTS`` check created a second 100–400 MB index on
    # production databases during the first incremental import.  Reuse any
    # semantically equivalent index and create one only for a genuinely new DB.
    for table, index_name, columns in (
        ("Compounds", "ix_c_ik", ("inchikey",)),
        ("Compounds", "ix_c_name", ("compound_name",)),
        ("Spectra", "ix_s_sid", ("spectrum_id",)),
        ("Spectra", "ix_s_pmz", ("precursor_mz",)),
        ("Mapping", "ix_m_c", ("compound_id",)),
    ):
        if not _has_simple_index(conn, table, columns):
            quoted_columns = ", ".join(f'"{column}"' for column in columns)
            conn.execute(
                f'CREATE INDEX "{index_name}" ON "{table}" ({quoted_columns})'
            )


def _ensure_quality_columns(conn: sqlite3.Connection) -> None:
    columns = {row[1] for row in conn.execute("PRAGMA table_info(Spectra)")}
    if "dataset_source" not in columns:
        conn.execute("ALTER TABLE Spectra ADD COLUMN dataset_source TEXT")
    if "source_tier" not in columns:
        conn.execute("ALTER TABLE Spectra ADD COLUMN source_tier TEXT")
    if "intensity_origin" not in columns:
        conn.execute("ALTER TABLE Spectra ADD COLUMN intensity_origin TEXT")
    if "quality_flags_json" not in columns:
        conn.execute("ALTER TABLE Spectra ADD COLUMN quality_flags_json TEXT")
    for column, column_type in (
        ("collision_energy_raw", "TEXT"),
        ("collision_energy_unit", "TEXT"),
        ("collision_energy_mode", "TEXT"),
        ("collision_energy_min", "REAL"),
        ("collision_energy_max", "REAL"),
        ("instrument_type", "TEXT"),
        ("instrument", "TEXT"),
        ("ion_mode", "TEXT"),
        ("source_accession", "TEXT"),
        ("activation_method", "TEXT"),
        ("instrument_domain_id", "TEXT"),
        ("collision_cell_region", "TEXT"),
        ("acquisition_method", "TEXT"),
        ("source_asset_sha256", "TEXT"),
        ("provenance_json", "TEXT"),
        ("reference_material_type", "TEXT"),
        ("standard_lot", "TEXT"),
        ("standard_purity", "REAL"),
        ("digitization_method", "TEXT"),
        ("source_image_file", "TEXT"),
    ):
        if column not in columns:
            conn.execute(
                f'ALTER TABLE Spectra ADD COLUMN "{column}" {column_type}'
            )
    compound_columns = {row[1] for row in conn.execute("PRAGMA table_info(Compounds)")}
    if "reported_inchikey" not in compound_columns:
        conn.execute("ALTER TABLE Compounds ADD COLUMN reported_inchikey TEXT")
    if "structure_identity_verified" not in compound_columns:
        conn.execute(
            "ALTER TABLE Compounds ADD COLUMN structure_identity_verified INTEGER DEFAULT 0"
        )
    if "structure_identity_source" not in compound_columns:
        conn.execute("ALTER TABLE Compounds ADD COLUMN structure_identity_source TEXT")
    # Preview/import performs this normalized identity lookup for every incoming
    # compound.  On the production library, omitting the expression index turns
    # each row into a full Compounds-table scan.
    conn.execute(
        "CREATE INDEX IF NOT EXISTS ix_c_name_formula_norm "
        "ON Compounds(lower(trim(compound_name)), molecular_formula)"
    )


def _find_existing_compound(conn: sqlite3.Connection, rec: dict) -> Optional[int]:
    """Find one chemically consistent compound without loading the full DB.

    Full InChIKey is safe for measured spectra.  Presence-only literature rows
    deliberately use name+formula so a copied SMILES cannot merge distinct
    reported compounds.  Name-only reuse is forbidden.
    """
    source_tier = str(rec.get("source_tier") or "").strip().lower()
    ik = _verified_inchikey(rec)
    reported_ik = _reported_full_inchikey(rec)
    if ik and source_tier != "literature_presence_only":
        row = conn.execute(
            "SELECT id FROM Compounds WHERE inchikey=? LIMIT 1", (ik,)
        ).fetchone()
        if row:
            return int(row[0])
        # A complete verified key that is absent from the DB identifies a new
        # compound.  Falling through to name+formula here used to merge an
        # isomer into a different existing InChIKey.
        return None
    if reported_ik and source_tier != "literature_presence_only":
        # Unverified submitted keys are not structural identity.  They are,
        # however, evidence that a name+formula fallback may contradict the
        # submitter's intended compound, so create an unlinked compound row.
        return None

    name = normalize_compound_name(rec.get("name"))
    formula = normalize_formula(rec.get("formula"))
    if name and formula:
        row = conn.execute(
            "SELECT id FROM Compounds WHERE lower(trim(compound_name))=? "
            "AND molecular_formula=? LIMIT 1",
            (name, formula),
        ).fetchone()
        if row:
            return int(row[0])
    return None


def _dataset_source(rec: dict) -> str:
    existing = str(rec.get("dataset_source") or "").strip()
    if existing:
        return existing
    source = Path(str(rec.get("source_file") or "user_contrib")).stem
    return f"user_contrib:{source}"[:240]


def _record_fingerprint(rec: dict) -> str:
    """Stable content fingerprint used for retry-safe spectrum IDs."""
    peaks: list[list[float]] = []
    for peak in rec.get("peaks") or []:
        try:
            peaks.append([round(float(peak[0]), 8), round(float(peak[1]), 8)])
        except (TypeError, ValueError, IndexError):
            continue
    payload = {
        "source_file": str(rec.get("source_file") or ""),
        "source_content_sha256": str(rec.get("source_content_sha256") or ""),
        "source_row_index": rec.get("source_row_index"),
        "dataset_source": _dataset_source(rec),
        "spectrum_id": str(rec.get("spectrum_id") or ""),
        "name": normalize_compound_name(rec.get("name")),
        "formula": normalize_formula(rec.get("formula")),
        "verified_inchikey": _verified_inchikey(rec),
        "reported_inchikey": _reported_full_inchikey(rec),
        "precursor_mz": round(float(rec.get("precursor_mz") or 0), 8),
        "retention_time": rec.get("retention_time"),
        "collision_energy": rec.get("collision_energy"),
        "collision_energy_raw": str(rec.get("collision_energy_raw") or ""),
        "collision_energy_unit": str(rec.get("collision_energy_unit") or ""),
        "collision_energy_mode": str(rec.get("collision_energy_mode") or ""),
        "collision_energy_min": rec.get("collision_energy_min"),
        "collision_energy_max": rec.get("collision_energy_max"),
        "charge": rec.get("charge"),
        "adduct": str(rec.get("adduct") or ""),
        "instrument_type": str(rec.get("instrument_type") or ""),
        "instrument": str(rec.get("instrument") or ""),
        "ion_mode": str(rec.get("ion_mode") or ""),
        "source_accession": str(rec.get("source_accession") or ""),
        "activation_method": str(rec.get("activation_method") or ""),
        "instrument_domain_id": str(rec.get("instrument_domain_id") or ""),
        "collision_cell_region": str(rec.get("collision_cell_region") or ""),
        "acquisition_method": str(rec.get("acquisition_method") or ""),
        "source_asset_sha256": str(rec.get("source_asset_sha256") or ""),
        "reference_material_type": str(rec.get("reference_material_type") or ""),
        "standard_lot": str(rec.get("standard_lot") or ""),
        "standard_purity": rec.get("standard_purity"),
        "digitization_method": str(rec.get("digitization_method") or ""),
        "source_image_file": str(rec.get("source_image_file") or ""),
        "plant_source": str(rec.get("plant_source") or ""),
        "source_tier": str(rec.get("source_tier") or ""),
        "intensity_origin": str(rec.get("intensity_origin") or ""),
        "peaks": peaks,
    }
    encoded = json.dumps(
        payload, sort_keys=True, ensure_ascii=False, separators=(",", ":"), default=str
    )
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _import_batch_id(records: list[dict]) -> str:
    """Order-independent deterministic ID for one sanitized upload payload."""
    fingerprints = sorted(_record_fingerprint(record) for record in records)
    encoded = json.dumps(fingerprints, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:32]


@contextmanager
def _library_import_lock(lock_root: Path | None = None):
    """Serialize imports beside the database set being mutated.

    V6 always placed the lock in the configured shared ``user_contrib``
    directory, even when a caller explicitly supplied an isolated/test
    database. Besides making portable deployments harder, that made a
    read-only shared library block otherwise safe local imports. V7 scopes the
    lock to the target database directory; multi-database production imports
    still share one lock because Full/Small live in the same canonical root.
    """
    root = Path(lock_root) if lock_root is not None else _user_contrib_dir()
    lock_path = root / ".library_import.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+b") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _remove_new_database_files(db_path: Path) -> None:
    for candidate in (
        Path(db_path),
        Path(f"{db_path}-wal"),
        Path(f"{db_path}-shm"),
    ):
        try:
            candidate.unlink(missing_ok=True)
        except OSError:
            logger.exception("无法清理失败导入新建的数据库文件: %s", candidate)


def _retry_spectrum_matches(
    row: sqlite3.Row | tuple[Any, ...], rec: dict
) -> bool:
    dataset_source, precursor_mz, adduct, peaks_json = row
    try:
        existing_peaks = json.loads(peaks_json or "[]")
    except (TypeError, ValueError, json.JSONDecodeError):
        return False
    expected_peaks = rec.get("peaks") or []
    return (
        str(dataset_source or "") == _dataset_source(rec)
        and abs(float(precursor_mz or 0) - float(rec.get("precursor_mz") or 0)) <= 1e-8
        and str(adduct or "") == str(rec.get("adduct") or "")
        and existing_peaks == expected_peaks
    )


def _append_records_with_receipt(
    records: list[dict],
    db_path: Path,
    *,
    ts_prefix: str,
    mass_threshold: Optional[float] = None,
    idempotent: bool = False,
) -> dict[str, Any]:
    """Write one spectral database and return an exact compensation receipt."""
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    existed_before = db_path.exists()
    receipt: dict[str, Any] = {
        "kind": "spectra",
        "db_path": str(db_path),
        "db_existed_before": existed_before,
        "new_compound_ids": [],
        "inserted_spectrum_ids": [],
        "mapping_rows": [],
        "new_compounds": 0,
        "new_spectra": 0,
    }
    conn: Optional[sqlite3.Connection] = None
    committed = False
    try:
        conn = sqlite3.connect(str(db_path))
        conn.execute("PRAGMA busy_timeout=60000")
        conn.execute("PRAGMA journal_mode=WAL")
        # ``CREATE ... IF NOT EXISTS`` also repairs a pre-created but empty DB
        # path (common in tests/admin tooling); relying only on path existence
        # would otherwise fail at the first ALTER TABLE.
        _create_schema(conn)
        _ensure_quality_columns(conn)

        conn.execute("BEGIN IMMEDIATE")
        for i, rec in enumerate(records):
            em = rec.get("exact_mass")
            if mass_threshold is not None:
                # Small is a neutral-mass-bounded library.  Unknown exact mass
                # is not evidence that a record is small; allowing it through
                # used to admit arbitrarily large precursor ions.
                if em is None or float(em) >= mass_threshold:
                    continue

            base = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(rec.get("spectrum_id") or i))
            if idempotent:
                # Per-record content identity stays stable even if the same
                # file is later submitted alongside additional files (which
                # legitimately changes the surrounding batch ID).
                sid = f"user_{_record_fingerprint(rec)[:48]}"
                existing_spectrum = conn.execute(
                    "SELECT dataset_source, precursor_mz, adduct, peaks_json "
                    "FROM Spectra WHERE spectrum_id=?", (sid,)
                ).fetchone()
                if existing_spectrum:
                    if not _retry_spectrum_matches(existing_spectrum, rec):
                        raise RuntimeError(
                            f"确定性 spectrum_id 冲突且内容不同，已拒绝覆盖: {sid}"
                        )
                    continue
            else:
                sid = f"{ts_prefix}_{i}_{base}"[:180]
                while conn.execute(
                    "SELECT 1 FROM Spectra WHERE spectrum_id=?", (sid,)
                ).fetchone():
                    sid = f"{ts_prefix}_{i}_{uuid.uuid4().hex[:10]}_{base}"[:190]

            ik = _verified_inchikey(rec)
            reported_ik = _reported_full_inchikey(rec)
            name = rec.get("name")
            smiles = rec.get("smiles")
            compound_id = _find_existing_compound(conn, rec)

            if compound_id is None:
                conn.execute(
                    "INSERT INTO Compounds (compound_name, molecular_formula, smiles, inchi, "
                    "inchikey, exact_mass, classification, compound_class, reported_inchikey, "
                    "structure_identity_verified, structure_identity_source) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    (
                        name, rec.get("formula"), smiles, rec.get("inchi"), ik, em,
                        rec.get("classification") or rec.get("compound_class") or "Unknown",
                        rec.get("compound_class"), reported_ik,
                        1 if ik else 0,
                        rec.get("structure_identity_source") or (
                            "smiles_derived" if ik else (
                                "reported_unverified" if reported_ik else "none"
                            )
                        ),
                    ),
                )
                compound_id = int(conn.execute("SELECT last_insert_rowid()").fetchone()[0])
                receipt["new_compound_ids"].append(compound_id)
                receipt["new_compounds"] += 1

            conn.execute(
                "INSERT INTO Spectra (spectrum_id, source_file, precursor_mz, retention_time, "
                "collision_energy, collision_energy_raw, collision_energy_unit, "
                "collision_energy_mode, collision_energy_min, collision_energy_max, "
                "charge, adduct, peaks_json, instrument_type, instrument, ion_mode, "
                "source_accession, dataset_source, source_tier, intensity_origin, "
                "quality_flags_json) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    sid, rec.get("source_file"), rec.get("precursor_mz"),
                    rec.get("retention_time"), rec.get("collision_energy"),
                    rec.get("collision_energy_raw"),
                    rec.get("collision_energy_unit"),
                    rec.get("collision_energy_mode"),
                    rec.get("collision_energy_min"),
                    rec.get("collision_energy_max"),
                    rec.get("charge"), rec.get("adduct"),
                    json.dumps(rec.get("peaks", [])),
                    rec.get("instrument_type"), rec.get("instrument"),
                    rec.get("ion_mode"), rec.get("source_accession"),
                    _dataset_source(rec), rec.get("source_tier") or "user_unverified",
                    rec.get("intensity_origin") or "unknown",
                    json.dumps(rec.get("quality_flags") or [], ensure_ascii=False),
                ),
            )
            spec_pk = int(conn.execute("SELECT last_insert_rowid()").fetchone()[0])
            conn.execute(
                "UPDATE Spectra SET activation_method=?, instrument_domain_id=?, "
                "collision_cell_region=?, acquisition_method=?, source_asset_sha256=?, "
                "provenance_json=?, reference_material_type=?, standard_lot=?, "
                "standard_purity=?, digitization_method=?, source_image_file=? "
                "WHERE id=?",
                (
                    rec.get("activation_method"),
                    rec.get("instrument_domain_id"),
                    rec.get("collision_cell_region"),
                    rec.get("acquisition_method"),
                    rec.get("source_asset_sha256"),
                    json.dumps(
                        rec.get("source_provenance") or {}, ensure_ascii=False
                    ),
                    rec.get("reference_material_type"),
                    rec.get("standard_lot"),
                    rec.get("standard_purity"),
                    rec.get("digitization_method"),
                    rec.get("source_image_file"),
                    spec_pk,
                ),
            )
            conn.execute(
                "INSERT INTO Mapping (compound_id, spectrum_id) VALUES (?,?)",
                (compound_id, spec_pk),
            )
            receipt["inserted_spectrum_ids"].append(sid)
            receipt["mapping_rows"].append(
                (compound_id, spec_pk, sid, name, ik, rec.get("precursor_mz"))
            )
            receipt["new_spectra"] += 1

        conn.execute("COMMIT")
        committed = True
        try:
            conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        except sqlite3.Error as checkpoint_error:
            # A concurrent reader can defer a checkpoint.  The transaction is
            # already durably committed; this is not an import failure.
            logger.warning("WAL checkpoint deferred for %s: %s", db_path, checkpoint_error)
    except Exception:
        if conn is not None and not committed:
            try:
                conn.execute("ROLLBACK")
            except sqlite3.Error:
                pass
        if conn is not None:
            conn.close()
            conn = None
        if not existed_before:
            _remove_new_database_files(db_path)
        raise
    finally:
        if conn is not None:
            conn.close()

    return receipt


def append_records(
    records: list[dict],
    db_path: Path,
    mapping_path: Optional[Path] = None,
    ts_prefix: Optional[str] = None,
    mass_threshold: Optional[float] = None,
) -> tuple[int, int]:
    """增量追加 records 到现有库（不删库）。

    - 复用已有化合物（完整 InChIKey / 名称+formula），命中则只加新谱图
    - 禁止按名称借用结构，避免同名/异构体或错误模板造成结构污染
    - spectrum_id 加 ts_prefix 前缀防撞 unique
    - mass_threshold: 若设置，仅追加 exact_mass < threshold 的（Small 库过滤）
    - 原子安全：单库写入使用 SQLite ``BEGIN IMMEDIATE``，失败完整回滚

    返回 (new_compounds, new_spectra)。
    """
    _assert_library_write_enabled()
    if ts_prefix is None:
        ts_prefix = f"user_{int(time.time())}"
    with _library_import_lock(Path(db_path).parent):
        receipt = _append_records_with_receipt(
            records,
            Path(db_path),
            ts_prefix=ts_prefix,
            mass_threshold=mass_threshold,
            idempotent=False,
        )
        try:
            if mapping_path and receipt["mapping_rows"]:
                _append_mapping_csv_atomic(Path(mapping_path), receipt["mapping_rows"])
        except Exception:
            _rollback_database_receipt(receipt)
            raise

    logger.info(
        "追加完成 %s: +%d 化合物, +%d 谱图",
        Path(db_path).name, receipt["new_compounds"], receipt["new_spectra"],
    )
    return int(receipt["new_compounds"]), int(receipt["new_spectra"])


_TARGET_SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS TargetCompounds (
    id INTEGER PRIMARY KEY,
    target_id TEXT UNIQUE NOT NULL,
    source_file TEXT,
    source_content_sha256 TEXT,
    dataset_version TEXT,
    source_row_index INTEGER,
    dataset_source TEXT,
    compound_name TEXT,
    molecular_formula TEXT,
    smiles TEXT,
    inchi TEXT,
    inchikey TEXT,
    reported_inchikey TEXT,
    structure_identity_verified INTEGER DEFAULT 0,
    structure_identity_source TEXT,
    exact_mass REAL,
    precursor_mz REAL,
    adduct TEXT,
    charge INTEGER,
    reported_precursor_mz REAL,
    reported_adduct TEXT,
    reported_retention_time REAL,
    reported_retention_time_unit TEXT,
    reported_collision_energy_raw TEXT,
    instrument TEXT,
    target_searchability TEXT,
    plant_source TEXT,
    provenance TEXT,
    doi TEXT,
    dataset_accession TEXT,
    spectrum_accession TEXT,
    source_asset_sha256 TEXT,
    source_image_file TEXT,
    source_image_sha256 TEXT,
    digitization_method TEXT,
    digitization_reviewer_confirmed INTEGER DEFAULT 0,
    intensity_origin TEXT,
    instrument_type TEXT,
    activation_method TEXT,
    collision_energy REAL,
    collision_energy_unit TEXT,
    collision_energy_mode TEXT,
    acquisition_method TEXT,
    reference_material_type TEXT,
    standard_lot TEXT,
    standard_purity REAL,
    provenance_json TEXT,
    source_tier TEXT NOT NULL DEFAULT 'literature_presence_only',
    quality_flags_json TEXT,
    created_at TEXT
);
CREATE TABLE IF NOT EXISTS DiagnosticIons (
    id INTEGER PRIMARY KEY,
    target_id INTEGER NOT NULL,
    fragment_mz REAL NOT NULL,
    ordinal INTEGER,
    FOREIGN KEY(target_id) REFERENCES TargetCompounds(id)
);
CREATE INDEX IF NOT EXISTS ix_target_dataset ON TargetCompounds(dataset_source);
CREATE INDEX IF NOT EXISTS ix_target_precursor ON TargetCompounds(precursor_mz);
CREATE INDEX IF NOT EXISTS ix_target_name_formula ON TargetCompounds(compound_name, molecular_formula);
CREATE INDEX IF NOT EXISTS ix_diagnostic_target ON DiagnosticIons(target_id);
CREATE INDEX IF NOT EXISTS ix_diagnostic_mz ON DiagnosticIons(fragment_mz);
"""


def _target_identifier(rec: dict) -> str:
    signature = json.dumps({
        "source": _dataset_source(rec),
        "name": normalize_compound_name(rec.get("name")),
        "formula": normalize_formula(rec.get("formula")),
        "precursor_mz": round(float(rec.get("precursor_mz") or 0), 5),
        "adduct": str(rec.get("adduct") or "").strip(),
        # Preserve which formula herb/source reported the target.  The same
        # chemical can legitimately occur in several component herbs.
        "plant_source": str(rec.get("plant_source") or "").strip().casefold(),
        "fragments": [round(mz, 5) for mz in _diagnostic_ions(rec)],
    }, sort_keys=True, ensure_ascii=False)
    return "target_" + hashlib.sha256(signature.encode("utf-8")).hexdigest()[:24]


def _diagnostic_ions(rec: dict) -> list[float]:
    """Return target ions as a stable, sorted, de-duplicated mass list."""
    ions: set[float] = set()
    for peak in rec.get("peaks") or []:
        try:
            mz = round(float(peak[0]), 8)
        except (TypeError, ValueError, IndexError):
            continue
        if mz > 0:
            ions.add(mz)
    return sorted(ions)


def _ensure_target_columns(conn: sqlite3.Connection) -> None:
    columns = {row[1] for row in conn.execute("PRAGMA table_info(TargetCompounds)")}
    additions = {
        "reported_inchikey": "TEXT",
        "structure_identity_verified": "INTEGER DEFAULT 0",
        "structure_identity_source": "TEXT",
        "source_content_sha256": "TEXT",
        "dataset_version": "TEXT",
        "source_row_index": "INTEGER",
        "provenance": "TEXT",
        "doi": "TEXT",
        "dataset_accession": "TEXT",
        "spectrum_accession": "TEXT",
        "reported_precursor_mz": "REAL",
        "reported_adduct": "TEXT",
        "reported_retention_time": "REAL",
        "reported_retention_time_unit": "TEXT",
        "reported_collision_energy_raw": "TEXT",
        "instrument": "TEXT",
        "target_searchability": "TEXT",
        "source_asset_sha256": "TEXT",
        "source_image_file": "TEXT",
        "source_image_sha256": "TEXT",
        "digitization_method": "TEXT",
        "digitization_reviewer_confirmed": "INTEGER DEFAULT 0",
        "intensity_origin": "TEXT",
        "instrument_type": "TEXT",
        "activation_method": "TEXT",
        "collision_energy": "REAL",
        "collision_energy_unit": "TEXT",
        "collision_energy_mode": "TEXT",
        "acquisition_method": "TEXT",
        "reference_material_type": "TEXT",
        "standard_lot": "TEXT",
        "standard_purity": "REAL",
        "provenance_json": "TEXT",
    }
    for column, declaration in additions.items():
        if column not in columns:
            conn.execute(
                f"ALTER TABLE TargetCompounds ADD COLUMN {column} {declaration}"
            )


def _append_target_records_with_receipt(
    records: list[dict],
    db_path: Path,
    *,
    replace_datasets: bool = False,
    protected_tables: Sequence[str] = (),
) -> dict[str, Any]:
    """Store targets transactionally and return their compensation receipt.

    A contributed workbook is a versioned dataset snapshot, not an endless
    append log.  ``replace_datasets=True`` atomically replaces every dataset
    present in ``records`` unless its exact target-ID set is already current.
    This removes superseded, chemically-invalid rows after a corrected upload
    while keeping retries idempotent.  Migration callers can provide
    ``protected_tables``; SQLite's authorizer will then reject any DML or DDL
    against those tables on this write connection.
    """
    db_path = Path(db_path)
    existed_before = db_path.exists()
    receipt: dict[str, Any] = {
        "kind": "targets",
        "db_path": str(db_path),
        "db_existed_before": existed_before,
        "inserted_target_ids": [],
        "new_targets": 0,
        "new_diagnostic_ions": 0,
        "removed_targets": 0,
        "replaced_datasets": [],
        "protected_tables": sorted({
            str(table).strip() for table in protected_tables
            if str(table).strip()
        }),
    }
    if not records:
        return receipt
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn: Optional[sqlite3.Connection] = None
    committed = False
    try:
        conn = sqlite3.connect(str(db_path))
        conn.execute("PRAGMA busy_timeout=60000")
        conn.execute("PRAGMA foreign_keys=ON")
        protected = {
            str(table).strip().casefold() for table in protected_tables
            if str(table).strip()
        }
        if protected:
            protected_actions = {
                getattr(sqlite3, name) for name in (
                    "SQLITE_INSERT",
                    "SQLITE_UPDATE",
                    "SQLITE_DELETE",
                    "SQLITE_CREATE_INDEX",
                    "SQLITE_CREATE_TABLE",
                    "SQLITE_CREATE_TEMP_INDEX",
                    "SQLITE_CREATE_TEMP_TABLE",
                    "SQLITE_CREATE_TEMP_TRIGGER",
                    "SQLITE_CREATE_TEMP_VIEW",
                    "SQLITE_CREATE_TRIGGER",
                    "SQLITE_CREATE_VIEW",
                    "SQLITE_DROP_INDEX",
                    "SQLITE_DROP_TABLE",
                    "SQLITE_DROP_TEMP_INDEX",
                    "SQLITE_DROP_TEMP_TABLE",
                    "SQLITE_DROP_TEMP_TRIGGER",
                    "SQLITE_DROP_TEMP_VIEW",
                    "SQLITE_DROP_TRIGGER",
                    "SQLITE_DROP_VIEW",
                    "SQLITE_ALTER_TABLE",
                    "SQLITE_REINDEX",
                ) if hasattr(sqlite3, name)
            }

            def _authorise(
                action: int,
                argument1: Optional[str],
                argument2: Optional[str],
                _database: Optional[str],
                _trigger: Optional[str],
            ) -> int:
                referenced = {
                    str(argument1 or "").strip().casefold(),
                    str(argument2 or "").strip().casefold(),
                }
                if action in protected_actions and protected.intersection(referenced):
                    return sqlite3.SQLITE_DENY
                return sqlite3.SQLITE_OK

            conn.set_authorizer(_authorise)
        conn.executescript(_TARGET_SCHEMA_SQL)
        _ensure_target_columns(conn)
        conn.execute("BEGIN IMMEDIATE")
        records_to_insert = list(records)
        if replace_datasets:
            grouped: dict[str, list[dict]] = {}
            for rec in records:
                grouped.setdefault(_dataset_source(rec), []).append(rec)
            records_to_insert = []
            for dataset_source, dataset_records in grouped.items():
                desired_ids = {_target_identifier(rec) for rec in dataset_records}
                desired_versions = {
                    str(rec.get("source_content_sha256") or "").strip()
                    for rec in dataset_records
                }
                existing_rows = conn.execute(
                    "SELECT id, target_id, dataset_version FROM TargetCompounds "
                    "WHERE dataset_source=?",
                    (dataset_source,),
                ).fetchall()
                existing_ids = {str(row[1]) for row in existing_rows}
                existing_versions = {str(row[2] or "").strip() for row in existing_rows}
                if existing_ids == desired_ids and existing_versions == desired_versions:
                    continue
                row_ids = [int(row[0]) for row in existing_rows]
                for row_chunk in _chunks(row_ids):
                    if row_chunk:
                        marks = ",".join("?" for _ in row_chunk)
                        conn.execute(
                            f"DELETE FROM DiagnosticIons WHERE target_id IN ({marks})",
                            row_chunk,
                        )
                conn.execute(
                    "DELETE FROM TargetCompounds WHERE dataset_source=?",
                    (dataset_source,),
                )
                receipt["removed_targets"] += len(existing_rows)
                receipt["replaced_datasets"].append(dataset_source)
                records_to_insert.extend(dataset_records)

        for rec in records_to_insert:
            target_id = _target_identifier(rec)
            existing = conn.execute(
                "SELECT id FROM TargetCompounds WHERE target_id=?", (target_id,)
            ).fetchone()
            if existing:
                continue
            verified_ik = _verified_inchikey(rec)
            reported_ik = _reported_full_inchikey(rec)
            conn.execute(
                "INSERT INTO TargetCompounds (target_id, source_file, source_content_sha256, "
                "dataset_version, source_row_index, dataset_source, "
                "compound_name, molecular_formula, smiles, inchi, inchikey, "
                "reported_inchikey, structure_identity_verified, structure_identity_source, "
                "exact_mass, precursor_mz, adduct, charge, reported_precursor_mz, "
                "reported_adduct, reported_retention_time, reported_retention_time_unit, "
                "reported_collision_energy_raw, instrument, target_searchability, "
                "plant_source, provenance, doi, "
                "dataset_accession, spectrum_accession, source_tier, "
                "quality_flags_json, created_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    target_id, rec.get("source_file"), rec.get("source_content_sha256"),
                    rec.get("source_content_sha256"), rec.get("source_row_index"),
                    _dataset_source(rec), rec.get("name"),
                    rec.get("formula"), rec.get("smiles"), rec.get("inchi"), verified_ik,
                    reported_ik, 1 if verified_ik else 0,
                    rec.get("structure_identity_source") or (
                        "smiles_derived" if verified_ik else (
                            "reported_unverified" if reported_ik else "none"
                        )
                    ),
                    rec.get("exact_mass"), rec.get("precursor_mz"), rec.get("adduct"),
                    rec.get("charge"), rec.get("reported_precursor_mz"),
                    rec.get("reported_adduct"), rec.get("reported_retention_time"),
                    rec.get("reported_retention_time_unit"),
                    rec.get("reported_collision_energy_raw"), rec.get("instrument"),
                    rec.get("target_searchability"), rec.get("plant_source"),
                    (
                        json.dumps(rec.get("provenance"), ensure_ascii=False, default=str)
                        if isinstance(rec.get("provenance"), (dict, list, tuple))
                        else rec.get("provenance")
                    ),
                    rec.get("doi"), rec.get("dataset_accession"),
                    rec.get("spectrum_accession"), "literature_presence_only",
                    json.dumps(rec.get("quality_flags") or [], ensure_ascii=False),
                    time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                ),
            )
            row_id = int(conn.execute("SELECT last_insert_rowid()").fetchone()[0])
            conn.execute(
                "UPDATE TargetCompounds SET source_asset_sha256=?, source_image_file=?, "
                "source_image_sha256=?, digitization_method=?, "
                "digitization_reviewer_confirmed=?, intensity_origin=?, "
                "instrument_type=?, activation_method=?, collision_energy=?, "
                "collision_energy_unit=?, collision_energy_mode=?, acquisition_method=?, "
                "reference_material_type=?, standard_lot=?, standard_purity=?, "
                "provenance_json=? WHERE id=?",
                (
                    rec.get("source_asset_sha256"),
                    rec.get("source_image_file"),
                    rec.get("source_image_sha256"),
                    rec.get("digitization_method"),
                    1 if rec.get("digitization_reviewer_confirmed") else 0,
                    rec.get("intensity_origin"),
                    rec.get("instrument_type"),
                    rec.get("activation_method"),
                    rec.get("collision_energy"),
                    rec.get("collision_energy_unit"),
                    rec.get("collision_energy_mode"),
                    rec.get("acquisition_method"),
                    rec.get("reference_material_type"),
                    rec.get("standard_lot"),
                    rec.get("standard_purity"),
                    json.dumps(
                        rec.get("source_provenance") or rec.get("provenance") or {},
                        ensure_ascii=False,
                        default=str,
                    ),
                    row_id,
                ),
            )
            ions = _diagnostic_ions(rec)
            conn.executemany(
                "INSERT INTO DiagnosticIons (target_id, fragment_mz, ordinal) VALUES (?,?,?)",
                [(row_id, mz, ordinal) for ordinal, mz in enumerate(ions, 1)],
            )
            receipt["inserted_target_ids"].append(target_id)
            receipt["new_targets"] += 1
            receipt["new_diagnostic_ions"] += len(ions)
        conn.execute("COMMIT")
        committed = True
    except Exception:
        if conn is not None and not committed:
            try:
                conn.execute("ROLLBACK")
            except sqlite3.Error:
                pass
        if conn is not None:
            conn.close()
            conn = None
        if not existed_before:
            _remove_new_database_files(db_path)
        raise
    finally:
        if conn is not None:
            conn.close()
    return receipt


def append_target_records(records: list[dict], db_path: Path) -> tuple[int, int]:
    """Store presence-only literature targets outside the spectral Top-K tables."""
    if not records:
        return 0, 0
    with _library_import_lock(Path(db_path).parent):
        receipt = _append_target_records_with_receipt(records, db_path)
    logger.info(
        "文献靶标追加完成: +%d targets, +%d diagnostic ions",
        receipt["new_targets"], receipt["new_diagnostic_ions"],
    )
    return int(receipt["new_targets"]), int(receipt["new_diagnostic_ions"])


def _mapping_csv_bytes(rows: list[tuple], include_header: bool) -> bytes:
    import io

    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer)
    if include_header:
        writer.writerow([
            "compound_id", "spectrum_row_id", "spectrum_id",
            "compound_name", "inchikey", "precursor_mz",
        ])
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def _append_mapping_csv_atomic(mapping_path: Path, rows: list[tuple]) -> dict[str, Any]:
    """Append one complete mapping block, truncating it on any write failure."""
    mapping_path = Path(mapping_path)
    mapping_path.parent.mkdir(parents=True, exist_ok=True)
    with mapping_path.open("a+b") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0, os.SEEK_END)
        previous_size = handle.tell()
        payload = _mapping_csv_bytes(rows, include_header=previous_size == 0)
        try:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        except Exception:
            handle.truncate(previous_size)
            handle.flush()
            os.fsync(handle.fileno())
            raise
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    return {
        "kind": "mapping",
        "path": str(mapping_path),
        "previous_size": previous_size,
        "expected_size": previous_size + len(payload),
    }


def _append_mapping_csv(mapping_path: Path, rows: list[tuple]) -> None:
    _append_mapping_csv_atomic(mapping_path, rows)


def _rollback_mapping_receipt(receipt: dict[str, Any]) -> None:
    path = Path(receipt["path"])
    if not path.exists() and int(receipt["previous_size"]) == 0:
        return
    with path.open("r+b") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0, os.SEEK_END)
        current_size = handle.tell()
        expected_size = int(receipt["expected_size"])
        if current_size != expected_size:
            raise RuntimeError(
                f"mapping CSV 在补偿窗口内被其他写入修改，拒绝截断: {path} "
                f"({current_size} != {expected_size})"
            )
        handle.truncate(int(receipt["previous_size"]))
        handle.flush()
        os.fsync(handle.fileno())
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    if int(receipt["previous_size"]) == 0:
        path.unlink(missing_ok=True)


def _chunks(values: list[Any], size: int = 400) -> list[list[Any]]:
    return [values[index:index + size] for index in range(0, len(values), size)]


def _rollback_database_receipt(receipt: dict[str, Any]) -> None:
    db_path = Path(receipt["db_path"])
    if not db_path.exists():
        return
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA busy_timeout=60000")
    try:
        conn.execute("BEGIN IMMEDIATE")
        if receipt["kind"] == "targets":
            target_ids = list(receipt.get("inserted_target_ids") or [])
            for target_chunk in _chunks(target_ids):
                placeholders = ",".join("?" for _ in target_chunk)
                row_ids = [
                    int(row[0]) for row in conn.execute(
                        f"SELECT id FROM TargetCompounds WHERE target_id IN ({placeholders})",
                        target_chunk,
                    )
                ]
                for row_chunk in _chunks(row_ids):
                    if not row_chunk:
                        continue
                    id_marks = ",".join("?" for _ in row_chunk)
                    conn.execute(
                        f"DELETE FROM DiagnosticIons WHERE target_id IN ({id_marks})", row_chunk
                    )
                conn.execute(
                    f"DELETE FROM TargetCompounds WHERE target_id IN ({placeholders})", target_chunk
                )
        elif receipt["kind"] == "spectra":
            spectrum_ids = list(receipt.get("inserted_spectrum_ids") or [])
            for spectrum_chunk in _chunks(spectrum_ids):
                placeholders = ",".join("?" for _ in spectrum_chunk)
                spectrum_pks = [
                    int(row[0]) for row in conn.execute(
                        f"SELECT id FROM Spectra WHERE spectrum_id IN ({placeholders})",
                        spectrum_chunk,
                    )
                ]
                for pk_chunk in _chunks(spectrum_pks):
                    if not pk_chunk:
                        continue
                    pk_marks = ",".join("?" for _ in pk_chunk)
                    conn.execute(
                        f"DELETE FROM Mapping WHERE spectrum_id IN ({pk_marks})", pk_chunk
                    )
                conn.execute(
                    f"DELETE FROM Spectra WHERE spectrum_id IN ({placeholders})", spectrum_chunk
                )
            for compound_id in receipt.get("new_compound_ids") or []:
                conn.execute(
                    "DELETE FROM Compounds WHERE id=? AND NOT EXISTS "
                    "(SELECT 1 FROM Mapping WHERE compound_id=?)",
                    (int(compound_id), int(compound_id)),
                )
        conn.execute("COMMIT")
    except Exception:
        try:
            conn.execute("ROLLBACK")
        except sqlite3.Error:
            pass
        raise
    finally:
        conn.close()
    if not receipt.get("db_existed_before"):
        _remove_new_database_files(db_path)


# ════════════════════════════════════════════════════════════════════
# 高层管线
# ════════════════════════════════════════════════════════════════════
def _stable_upload_identity(file_path: Path) -> tuple[str, str]:
    """Return ``(original_basename, sha256)`` independent of web temp UUIDs.

    The upload endpoint stores raw files as ``<8hex>_<original name>``.  That
    prefix is storage identity only and changes on every retry; letting it flow
    into dataset IDs previously defeated idempotency.  Direct CLI files keep
    their basename unchanged.
    """
    file_path = Path(file_path)
    match = re.fullmatch(r"[0-9a-fA-F]{8}_(.+)", file_path.name)
    original_name = match.group(1) if match else file_path.name
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return original_name, digest.hexdigest()


def _files_to_records(
    file_paths: list[Path], *, archive_mgf: bool = True, strict: bool = False
) -> tuple[list[dict], list[dict]]:
    """Parse tables or audited raw-reference batches into import records.

    CSV/Excel inputs retain the existing direct-table path and optional MGF
    archive. Raw MGF/mzML and image+manifest inputs are parsed by the V15
    reference-intake contract. A spectrum image never becomes a searchable
    experimental spectrum merely because numeric peaks were digitized from it.
    """
    from csv_excel_to_mgf import parse_file, write_mgf

    mgf_dir = _user_contrib_dir() / "mgf"
    if archive_mgf:
        mgf_dir.mkdir(parents=True, exist_ok=True)

    all_records: list[dict] = []
    convert_stats: list[dict] = []
    table_extensions = {".csv", ".xlsx", ".xls"}
    table_paths = [
        Path(path) for path in file_paths
        if Path(path).suffix.lower() in table_extensions
    ]
    reference_paths = [
        Path(path) for path in file_paths
        if Path(path).suffix.lower() not in table_extensions
    ]
    for fp in table_paths:
        fp = Path(fp)
        try:
            original_name, content_sha256 = _stable_upload_identity(fp)
            original_stem = Path(original_name).stem
            table_records, fmt = parse_file(fp, format_hint="auto")
            if strict and not table_records:
                raise ValueError("未解析到任何有效记录")
            if archive_mgf:
                # Content-addressed provenance archive: retrying the same raw
                # workbook reuses one path instead of generating duplicates.
                mgf_out = mgf_dir / f"{original_stem}_{content_sha256[:16]}.mgf"
                mgf_tmp = mgf_out.with_name(
                    f".{mgf_out.name}.{uuid.uuid4().hex}.part"
                )
                try:
                    stat = write_mgf(table_records, mgf_tmp)
                    os.replace(mgf_tmp, mgf_out)
                finally:
                    mgf_tmp.unlink(missing_ok=True)
            else:
                stat = {
                    "total": len(table_records),
                    "written": sum(
                        1 for row in table_records
                        if row.get("precursor_mz") is not None
                        and float(row.get("precursor_mz") or 0) > 0
                    ),
                    "skipped_no_mz": sum(
                        1 for row in table_records
                        if row.get("precursor_mz") is None
                        or float(row.get("precursor_mz") or 0) <= 0
                    ),
                    "skipped_no_frag": sum(
                        1 for row in table_records if not (row.get("fragments") or [])
                    ),
                }
            stat["format"] = fmt
        except Exception as e:
            logger.error("转换失败 %s: %s", fp.name, e)
            convert_stats.append({"file": fp.name, "error": str(e), "total": 0, "written": 0})
            if strict:
                raise ValueError(f"批次文件转换失败，数据库未修改: {fp.name}: {e}") from e
            continue
        convert_stats.append({"file": fp.name, **{k: stat.get(k) for k in
                                ("total", "written", "skipped_no_mz", "skipped_no_frag", "format")}})
        for index, row in enumerate(table_records):
            fragments = row.get("fragments") or []
            intensities = row.get("fragment_intensities") or []
            peaks = []
            for peak_index, fragment_mz in enumerate(fragments):
                intensity = intensities[peak_index] if peak_index < len(intensities) else None
                peaks.append([float(fragment_mz), float(intensity or 1000.0)])
            rec = dict(row)
            rec.update({
                "spectrum_id": f"{original_stem}_{index}",
                "source_file": original_name,
                "source_content_sha256": content_sha256,
                "source_row_index": index,
                "dataset_source": f"user_contrib:{original_stem}"[:240],
                "peaks": peaks,
                "classification": row.get("component_class") or "Unknown",
                "compound_class": row.get("component_class") or "Unknown",
            })
            all_records.append(rec)
        logger.info("转换 %s: %d 条 → MGF", original_name, len(table_records))

    if reference_paths:
        from matching.reference_spectrum_intake import parse_reference_intake

        reference_records, reference_stats = parse_reference_intake(
            reference_paths, strict=strict
        )
        all_records.extend(reference_records)
        convert_stats.extend(reference_stats)

    # 归一化 formula + 估算中性 exact_mass。绝不再用 precursor_mz 冒充
    # 中性质量；两者差一个加合物，会污染 neutral-mass 索引。
    for r in all_records:
        nf = normalize_formula(r.get("formula"))
        if nf:
            r["formula"] = nf
        if not r.get("exact_mass"):
            r["exact_mass"] = (_exact_mass_from_smiles(r.get("smiles"))
                               or _exact_mass_from_formula(r.get("formula")))
    return all_records, convert_stats


def preview_import(file_paths: list[Path], db_path: Optional[Path] = None) -> dict:
    """只转换 + 与现有库比对，不写库。返回预览统计。"""
    if db_path is None:
        db_path = app_settings.DB_PATH

    # Preview is read-only: parsing a workbook must not create a provenance
    # archive or touch the production data directory before user confirmation.
    records, convert_stats = _files_to_records(file_paths, archive_mgf=False)
    from matching.library_quality import assess_import_records
    records, quality = assess_import_records(records)

    from matching.library_quality import PRESENCE_ONLY_SOURCE_TIER, TRUSTED_SOURCE_TIERS
    target_records = [
        record for record in records
        if str(record.get("source_tier") or "").lower() == PRESENCE_ONLY_SOURCE_TIER
    ]
    spectrum_records = [
        record for record in records
        if str(record.get("source_tier") or "").lower() in TRUSTED_SOURCE_TIERS
    ]
    provenance_held_records = [
        record for record in records
        if record not in target_records and record not in spectrum_records
    ]

    # Deduplicate only real spectra.  Literature target rows intentionally keep
    # herb/source provenance and must not be collapsed by a shared ion list.
    try:
        from build_from_mgf_msp import deduplicate
        before = len(spectrum_records)
        spectrum_records = deduplicate(spectrum_records)
        dedup_removed = before - len(spectrum_records)
    except Exception:
        dedup_removed = 0
    records = spectrum_records + target_records

    existing_ik, existing_name_formula, ik_has_smi, name_formula_has_smi = _load_existing_keys(
        Path(db_path), records
    )
    new_comp = 0
    existing_comp = 0
    no_smiles = 0
    existing_with_struct = 0   # 已有且库中含结构
    new_no_struct = 0          # 新化合物且表格无SMILES（导入时将尝试借用补全）
    samples: list[dict] = []
    for r in records:
        verified_ik = _verified_inchikey(r)
        reported_ik = _reported_full_inchikey(r)
        nname = normalize_compound_name(r.get("name"))
        nformula = normalize_formula(r.get("formula"))
        name_formula = (
            f"{nname}::{nformula}" if nname and nformula else None
        )
        if verified_ik:
            is_existing = verified_ik in existing_ik
        elif reported_ik:
            # A submitted key that could not be structure-verified is kept as
            # provenance only.  It also prevents a contradictory
            # name+formula fallback in the preview, matching the write path.
            is_existing = False
        else:
            is_existing = bool(
                name_formula and name_formula in existing_name_formula
            )
        tbl_has_smi = bool(r.get("smiles"))
        db_has_smi: Optional[bool] = None
        if is_existing:
            existing_comp += 1
            if verified_ik and verified_ik in ik_has_smi:
                db_has_smi = ik_has_smi[verified_ik]
            elif name_formula and name_formula in name_formula_has_smi:
                db_has_smi = name_formula_has_smi[name_formula]
            if db_has_smi:
                existing_with_struct += 1
        else:
            new_comp += 1
            if not tbl_has_smi:
                new_no_struct += 1
        if not tbl_has_smi:
            no_smiles += 1
        if len(samples) < 100:
            samples.append({
                "name": r.get("name"),
                "formula": r.get("formula"),
                "precursor_mz": r.get("precursor_mz"),
                "is_new": not is_existing,
                "has_smiles": tbl_has_smi,
                "db_has_smiles": db_has_smi,
                "reported_inchikey": r.get("reported_inchikey"),
                "structure_identity_verified": bool(
                    r.get("structure_identity_verified")
                ),
                "n_peaks": len(r.get("peaks", [])),
                "source_file": r.get("source_file"),
            })

    return {
        "convert": convert_stats,
        "total": len(records),
        "spectral_records": len(spectrum_records),
        "literature_target_records": len(target_records),
        "dedup_removed": dedup_removed,
        "new_compounds": new_comp,
        "existing_compounds": existing_comp,
        "existing_with_structure": existing_with_struct,
        "no_smiles": no_smiles,
        "new_no_structure": new_no_struct,
        "samples": samples,
        "quality": quality,
        "importable_records": len(records),
        "provenance_held_records": len(provenance_held_records),
        "rejected_records": quality["rejected"],
    }


def import_files(
    file_paths: list[Path],
    db_mode: str = "both",
) -> dict:
    """CSV/Excel → 转换 → 归一化去重 → 增量追加到 Full / Small 双库。

    db_mode: "both" | "full" | "small"
    返回统计 dict（含每个库的新增化合物/谱图数）。
    """
    _assert_library_write_enabled()
    full_db_path = Path(app_settings.DB_PATH)
    small_db_path = Path(app_settings.DB_PATH_SMALL)
    mapping_full_path = _configured_data_path(
        "library_mapping.csv", app_settings.LIBRARY_MAPPING_PATH
    )
    mapping_small_path = _configured_data_path(
        "library_mapping_small.csv", app_settings.LIBRARY_MAPPING_SMALL_PATH
    )
    if db_mode not in {"both", "full", "small"}:
        raise ValueError(f"不支持的数据库模式: {db_mode}")

    records, convert_stats = _files_to_records(
        file_paths, archive_mgf=True, strict=True
    )
    from matching.library_quality import assess_import_records
    records, quality = assess_import_records(records)

    if not records:
        return {
            "convert": convert_stats,
            "total_records": 0,
            "dedup_removed": 0,
            "quality": quality,
            "error": "没有通过质量门控的可导入记录；数据库未修改。",
        }

    from matching.library_quality import PRESENCE_ONLY_SOURCE_TIER, TRUSTED_SOURCE_TIERS
    target_records = [
        record for record in records
        if str(record.get("source_tier") or "").lower() == PRESENCE_ONLY_SOURCE_TIER
    ]
    spectrum_records = [
        record for record in records
        if str(record.get("source_tier") or "").lower() in TRUSTED_SOURCE_TIERS
    ]
    provenance_held_records = [
        record for record in records
        if record not in target_records and record not in spectrum_records
    ]

    try:
        from build_from_mgf_msp import deduplicate
        before = len(spectrum_records)
        spectrum_records = deduplicate(spectrum_records)
        dedup_removed = before - len(spectrum_records)
    except Exception:
        dedup_removed = 0

    result: dict[str, Any] = {
        "convert": convert_stats,
        "total_records": len(spectrum_records) + len(target_records),
        "spectral_records": len(spectrum_records),
        "literature_target_records": len(target_records),
        "provenance_held_records": len(provenance_held_records),
        "dedup_removed": dedup_removed,
        "quality": quality,
    }

    if not spectrum_records and not target_records:
        result.update({
            "database_modified": False,
            "error": (
                "通过化学质量门控的记录均缺少可验证实验来源；已保留为预览审计，"
                "但未写入普通检索谱库。请提供 SOURCE_TIER 以及仪器/数据集/DOI 来源。"
            ),
        })
        return result

    import_records = spectrum_records + target_records
    batch_id = _import_batch_id(import_records)
    ts_prefix = f"user_{batch_id[:24]}"
    result["import_batch_id"] = batch_id
    database_receipts: list[dict[str, Any]] = []
    pending_mappings: list[tuple[Path, list[tuple]]] = []
    mapping_receipts: list[dict[str, Any]] = []

    # A process-wide file lock keeps the compensation interval isolated from
    # another web worker.  Each SQLite database still uses its own immediate
    # transaction; if a later database/target/mapping write raises, the exact
    # inserted row IDs are removed in reverse order.  Stable spectrum/target
    # IDs additionally make a retry idempotent after an abrupt process exit.
    with _library_import_lock(full_db_path.parent):
        try:
            if spectrum_records and db_mode in ("both", "full"):
                full_receipt = _append_records_with_receipt(
                    spectrum_records,
                    full_db_path,
                    ts_prefix=ts_prefix,
                    idempotent=True,
                )
                database_receipts.append(full_receipt)
                pending_mappings.append((
                    mapping_full_path,
                    full_receipt["mapping_rows"],
                ))
                result["full"] = {
                    "new_compounds": full_receipt["new_compounds"],
                    "new_spectra": full_receipt["new_spectra"],
                }

            if spectrum_records and db_mode in ("both", "small"):
                small_receipt = _append_records_with_receipt(
                    spectrum_records,
                    small_db_path,
                    ts_prefix=ts_prefix,
                    mass_threshold=SMALL_MASS_THRESHOLD,
                    idempotent=True,
                )
                database_receipts.append(small_receipt)
                pending_mappings.append((
                    mapping_small_path,
                    small_receipt["mapping_rows"],
                ))
                result["small"] = {
                    "new_compounds": small_receipt["new_compounds"],
                    "new_spectra": small_receipt["new_spectra"],
                }

            # Mapping CSV is an audit export, but it must agree with the DB.
            # Append it only after every DB transaction succeeded and retain a
            # byte-offset receipt so later failure cannot leave a partial map.
            for mapping_path, rows in pending_mappings:
                if rows:
                    mapping_receipts.append(
                        _append_mapping_csv_atomic(mapping_path, rows)
                    )

            # Dataset replacement is deliberately last.  Its own transaction
            # rolls back on failure; after it commits there is no later
            # fallible side effect that could require reconstructing deleted
            # superseded target rows.
            if target_records:
                target_receipt = _append_target_records_with_receipt(
                    target_records, full_db_path, replace_datasets=True,
                )
                database_receipts.append(target_receipt)
                result["literature_targets"] = {
                    "new_targets": target_receipt["new_targets"],
                    "new_diagnostic_ions": target_receipt["new_diagnostic_ions"],
                    "removed_superseded_targets": target_receipt["removed_targets"],
                    "replaced_datasets": target_receipt["replaced_datasets"],
                    "database": "full",
                    "ordinary_spectral_search": False,
                }
        except Exception as import_error:
            rollback_errors: list[str] = []
            for receipt in reversed(mapping_receipts):
                try:
                    _rollback_mapping_receipt(receipt)
                except Exception as error:
                    rollback_errors.append(f"mapping: {error}")
            for receipt in reversed(database_receipts):
                try:
                    _rollback_database_receipt(receipt)
                except Exception as error:
                    rollback_errors.append(
                        f"{receipt.get('kind')}@{receipt.get('db_path')}: {error}"
                    )
            if rollback_errors:
                raise RuntimeError(
                    "导入失败，且补偿回滚未完全成功: " + "; ".join(rollback_errors)
                ) from import_error
            raise

    result["database_modified"] = any(
        int(receipt.get("new_spectra", 0)) > 0
        or int(receipt.get("new_targets", 0)) > 0
        for receipt in database_receipts
    )
    return result
