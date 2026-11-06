"""Single-source polarity classification and database coverage auditing."""
from __future__ import annotations

import json
import os
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from utils.adduct import adduct_polarity


POLARITY_RULE_VERSION = (
    "v7-explicit-ion-mode_then-adduct_then-unambiguous-signed-charge-2"
)


def _explicit_mode(value: Any) -> str | None:
    text = str(value or "").strip().upper().replace("_", " ")
    if text in {"POS", "POSITIVE", "P", "+", "ESI+", "ESI POSITIVE"}:
        return "POS"
    if text in {"NEG", "NEGATIVE", "N", "-", "ESI-", "ESI NEGATIVE"}:
        return "NEG"
    return None


def _charge_mode(value: Any) -> str | None:
    if value is None or isinstance(value, bool):
        return None
    text = str(value).strip()
    if not text:
        return None
    if text.endswith("+"):
        return "POS"
    if text.endswith("-"):
        return "NEG"
    if text.startswith("+"):
        return "POS"
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    # Historical library schemas store charge magnitude as a positive integer
    # for both ion modes.  A bare ``1`` therefore cannot prove positive
    # polarity; only an explicit plus sign can.  Negative numbers are
    # unambiguous and remain valid evidence.
    if number < 0:
        return "NEG"
    return None


def polarity_evidence(
    adduct: Any = None,
    charge: Any = None,
    ion_mode: Any = None,
) -> dict[str, Any]:
    """Classify polarity with a deterministic, conflict-aware precedence."""

    explicit = _explicit_mode(ion_mode)
    adduct_mode = adduct_polarity(adduct)
    charge_mode = _charge_mode(charge)
    available = [
        (name, value)
        for name, value in (
            ("explicit_ion_mode", explicit),
            ("adduct", adduct_mode),
            ("signed_charge", charge_mode),
        )
        if value
    ]
    chosen_source, chosen = available[0] if available else ("none", None)
    conflict = len({value for _, value in available}) > 1
    return {
        "polarity": chosen,
        "source": chosen_source,
        "conflict": conflict,
        "explicit_ion_mode": explicit,
        "adduct_polarity": adduct_mode,
        "charge_polarity": charge_mode,
        "rule_version": POLARITY_RULE_VERSION,
    }


def classify_spectrum_polarity(
    adduct: Any = None,
    charge: Any = None,
    ion_mode: Any = None,
) -> str | None:
    return polarity_evidence(adduct, charge, ion_mode)["polarity"]


def _columns(connection: sqlite3.Connection, table: str) -> set[str]:
    return {
        str(row[1])
        for row in connection.execute(f'PRAGMA table_info("{table}")')
    }


def _connect_readonly(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(
        f"file:{path.resolve()}?mode=ro", uri=True, timeout=120
    )
    connection.execute("PRAGMA busy_timeout=120000")
    connection.create_function(
        "mssystem_polarity",
        3,
        lambda adduct, charge, ion_mode: (
            classify_spectrum_polarity(adduct, charge, ion_mode) or "UNKNOWN"
        ),
        deterministic=True,
    )
    connection.create_function(
        "mssystem_polarity_conflict",
        3,
        lambda adduct, charge, ion_mode: int(
            polarity_evidence(adduct, charge, ion_mode)["conflict"]
        ),
        deterministic=True,
    )
    return connection


def recompute_polarity_snapshot(path: str | Path) -> dict[str, Any]:
    """Recompute spectrum and compound coverage from physical SQLite rows."""

    db_path = Path(path).expanduser().resolve()
    stat = db_path.stat()
    connection = _connect_readonly(db_path)
    spectra_columns = _columns(connection, "Spectra")
    ion_mode_expr = "ion_mode" if "ion_mode" in spectra_columns else "NULL"
    classified = (
        "mssystem_polarity(adduct, charge, "
        f"{ion_mode_expr})"
    )
    conflict = (
        "mssystem_polarity_conflict(adduct, charge, "
        f"{ion_mode_expr})"
    )
    try:
        spectrum_row = connection.execute(
            f"""
            SELECT COUNT(*) AS total,
                   SUM({classified}='POS') AS positive,
                   SUM({classified}='NEG') AS negative,
                   SUM({classified}='UNKNOWN') AS unknown,
                   SUM({conflict}) AS conflicts
            FROM Spectra
            """
        ).fetchone()
        compound_row = connection.execute(
            f"""
            WITH compound_modes AS (
                SELECT m.compound_id,
                       MAX({classified}='POS') AS has_positive,
                       MAX({classified}='NEG') AS has_negative,
                       MAX({classified}='UNKNOWN') AS has_unknown
                FROM Mapping m
                JOIN Spectra ON Spectra.id = m.spectrum_id
                GROUP BY m.compound_id
            )
            SELECT COUNT(*) AS mapped,
                   SUM(has_positive) AS with_positive,
                   SUM(has_negative) AS with_negative,
                   SUM(has_positive=1 AND has_negative=1) AS both_modes,
                   SUM(has_positive=1 AND has_negative=0) AS positive_only,
                   SUM(has_positive=0 AND has_negative=1) AS negative_only,
                   SUM(has_positive=0 AND has_negative=0 AND has_unknown=1)
                       AS unknown_only
            FROM compound_modes
            """
        ).fetchone()
        check = connection.execute("PRAGMA quick_check").fetchone()
    finally:
        connection.close()

    total = int(spectrum_row[0] or 0)
    mapped = int(compound_row[0] or 0)
    result = {
        "path": str(db_path),
        "database_size_bytes": int(stat.st_size),
        "database_mtime_ns": int(stat.st_mtime_ns),
        "recomputed_at": datetime.now(timezone.utc).isoformat(),
        "polarity_rule_version": POLARITY_RULE_VERSION,
        "polarity_precedence": [
            "explicit ion_mode when physically present",
            "specific adduct sign",
            "unambiguous signed charge (+1/-1/1+/1-; bare positive magnitude ignored)",
            "UNKNOWN (never guessed from precursor m/z)",
        ],
        "quick_check": check[0] if check else None,
        "total_spectra": total,
        "positive_spectra": int(spectrum_row[1] or 0),
        "negative_spectra": int(spectrum_row[2] or 0),
        "unknown_polarity_spectra": int(spectrum_row[3] or 0),
        "polarity_metadata_conflicts": int(spectrum_row[4] or 0),
        "mapped_compounds": mapped,
        "compounds_with_positive": int(compound_row[1] or 0),
        "compounds_with_negative": int(compound_row[2] or 0),
        "both_polarities": int(compound_row[3] or 0),
        "positive_only_compounds": int(compound_row[4] or 0),
        "negative_only_compounds": int(compound_row[5] or 0),
        "unknown_only_compounds": int(compound_row[6] or 0),
    }
    result["coverage_percent"] = {
        "positive": (
            round(100.0 * result["compounds_with_positive"] / mapped, 3)
            if mapped else None
        ),
        "negative": (
            round(100.0 * result["compounds_with_negative"] / mapped, 3)
            if mapped else None
        ),
        "both": (
            round(100.0 * result["both_polarities"] / mapped, 3)
            if mapped else None
        ),
    }
    result["interpretation"] = (
        "Coverage uses one runtime-identical polarity rule. Unknown and conflict "
        "counts are explicit; no spectrum is assigned from precursor m/z alone."
    )
    return result


def database_signature(path: str | Path) -> tuple[int, int]:
    stat = Path(path).expanduser().resolve().stat()
    return int(stat.st_size), int(stat.st_mtime_ns)


def load_cached_snapshot(
    cache_path: str | Path,
    db_path: str | Path,
) -> dict[str, Any] | None:
    try:
        payload = json.loads(Path(cache_path).read_text(encoding="utf-8"))
    except Exception:
        return None
    size, mtime = database_signature(db_path)
    if (
        int(payload.get("database_size_bytes") or -1) != size
        or int(payload.get("database_mtime_ns") or -1) != mtime
        or payload.get("polarity_rule_version") != POLARITY_RULE_VERSION
    ):
        return None
    return payload


def compound_quality_snapshot(path: str | Path) -> dict[str, Any]:
    """Return database-native structure provenance counts."""

    db_path = Path(path).expanduser().resolve()
    connection = _connect_readonly(db_path)
    columns = _columns(connection, "Compounds")
    reported = (
        "reported_inchikey" if "reported_inchikey" in columns else "NULL"
    )
    verified = (
        "structure_identity_verified"
        if "structure_identity_verified" in columns
        else "0"
    )
    source = (
        "structure_identity_source"
        if "structure_identity_source" in columns
        else "NULL"
    )
    try:
        row = connection.execute(
            f"""
            SELECT COUNT(*) AS total,
                   SUM(smiles IS NULL OR trim(smiles)='') AS missing_smiles,
                   SUM(inchikey IS NULL OR trim(inchikey)='') AS missing_inchikey,
                   SUM(COALESCE({verified},0)=1) AS structure_verified,
                   SUM({source} IS NOT NULL AND trim({source})!='')
                       AS with_structure_source,
                   SUM({reported} IS NOT NULL AND trim({reported})!='')
                       AS with_reported_inchikey
            FROM Compounds
            """
        ).fetchone()
        sources = dict(
            connection.execute(
                f"""
                SELECT COALESCE(NULLIF(trim({source}),''),'(missing)') AS source,
                       COUNT(*)
                FROM Compounds GROUP BY 1 ORDER BY 2 DESC
                """
            ).fetchall()
        )
    finally:
        connection.close()
    return {
        "total_compounds": int(row[0] or 0),
        "missing_smiles": int(row[1] or 0),
        "missing_inchikey": int(row[2] or 0),
        "structure_identity_verified": int(row[3] or 0),
        "with_structure_identity_source": int(row[4] or 0),
        "with_reported_inchikey": int(row[5] or 0),
        "structure_identity_sources": {
            str(key): int(value) for key, value in sources.items()
        },
    }
