#!/usr/bin/env python3
"""Recover per-spectrum MoNA provenance lost by the legacy V6 builder.

The old MSP parser generated ``<stem>_<accepted-index>`` IDs but discarded
``Instrument_type`` and ``Spectrum_type``.  Consequently in-silico and MS1
records could be labelled experimental.  This migration streams the original
MSP (it never loads the multi-GB file into memory), reproduces the old accepted
index exactly, and relabels only IDs that still exist after the old dedup step.
No spectrum or compound row is deleted.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import time
from pathlib import Path
from typing import Iterable


MIGRATION_ID = "v6_20260716_recover_mona_in_silico_ms1_provenance"
PREDICTED_TOKENS = (
    "in-silico", "in silico", "predicted", "theoretical",
    "cfm-id", "cfmid", "sirius", "metfrag",
)


def _entry_classification(lines: list[str]) -> tuple[bool, str | None]:
    """Return ``(legacy_parser_accepted, recovered_tier)`` for one MSP block."""
    precursor: float | None = None
    fallback_mw: float | None = None
    spectrum_type = ""
    provenance_text: list[str] = []
    peaks: list[tuple[float, float]] = []

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if ":" in line and not line[0].isdigit():
            key, _, value = line.partition(":")
            key_upper = key.strip().upper()
            value = value.strip()
            if key_upper in {"PRECURSORMZ", "PRECURSOR M/Z", "PEPMASS"}:
                try:
                    precursor = float(value.split()[0])
                except (ValueError, IndexError):
                    pass
            elif key_upper == "MW" and precursor is None:
                try:
                    fallback_mw = float(value.split()[0])
                except (ValueError, IndexError):
                    pass
            elif key_upper in {"INSTRUMENT_TYPE", "INSTRUMENT TYPE"}:
                provenance_text.append(value.lower())
            elif key_upper in {"SOURCE_TIER", "SOURCE TIER"}:
                provenance_text.append(value.lower())
            elif key_upper in {"SPECTRUM_TYPE", "SPECTRUM TYPE"}:
                spectrum_type = value.lower()
                provenance_text.append(spectrum_type)
            continue

        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            mz, intensity = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        if mz > 0 and intensity > 0:
            peaks.append((mz, intensity))

    precursor = precursor if precursor is not None else fallback_mw
    if precursor is None or not peaks:
        return False, None
    maximum = max(intensity for _, intensity in peaks)
    # Exact acceptance predicate from the legacy _clean_peaks implementation.
    accepted = any(
        intensity / maximum >= 0.01 and mz < precursor - 1.5
        for mz, intensity in peaks
    )
    if not accepted:
        return False, None

    if "ms1" in spectrum_type.replace(" ", ""):
        return True, "quarantined"
    evidence = " ".join(provenance_text)
    if any(token in evidence for token in PREDICTED_TOKENS):
        return True, "predicted"
    return True, None


def iter_recovered_ids(
    msp_path: Path, *, id_offset: int = 0,
) -> Iterable[tuple[str, str]]:
    """Yield legacy spectrum IDs that are provably MS1 or in-silico."""
    stem = msp_path.stem
    accepted_index = 0
    entry: list[str] = []

    def finish(block: list[str]):
        nonlocal accepted_index
        if not block:
            return None
        accepted, tier = _entry_classification(block)
        if not accepted:
            return None
        spectrum_id = f"{stem}_{accepted_index + id_offset}"
        accepted_index += 1
        return (spectrum_id, tier) if tier else None

    with msp_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                entry.append(line)
                continue
            result = finish(entry)
            entry = []
            if result:
                yield result
        result = finish(entry)
        if result:
            yield result


def _backup_database(db_path: Path) -> Path:
    stamp = time.strftime("%Y%m%dT%H%M%S")
    backup = db_path.with_name(f"{db_path.name}.pre_mona_provenance_{stamp}.bak")
    source = sqlite3.connect(f"file:{db_path.resolve()}?mode=ro", uri=True)
    destination = sqlite3.connect(str(backup))
    try:
        source.backup(destination)
        destination.commit()
    except Exception:
        backup.unlink(missing_ok=True)
        raise
    finally:
        destination.close()
        source.close()
    check = sqlite3.connect(f"file:{backup.resolve()}?mode=ro", uri=True)
    try:
        if check.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise RuntimeError("backup quick_check failed")
        if list(check.execute("PRAGMA foreign_key_check")):
            raise RuntimeError("backup foreign_key_check failed")
    except Exception:
        backup.unlink(missing_ok=True)
        raise
    finally:
        check.close()
    return backup


def migrate(
    db_path: Path, msp_path: Path, *, dry_run: bool = False,
    id_offset: int = 0,
) -> dict:
    db_path, msp_path = db_path.resolve(), msp_path.resolve()
    if not db_path.is_file() or not msp_path.is_file():
        raise FileNotFoundError(db_path if not db_path.is_file() else msp_path)

    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA busy_timeout=60000")
    conn.execute(
        "CREATE TEMP TABLE recovered_mona_tier("
        "spectrum_id TEXT PRIMARY KEY, tier TEXT NOT NULL) WITHOUT ROWID"
    )
    scanned = {"predicted": 0, "quarantined": 0}
    batch: list[tuple[str, str]] = []
    for spectrum_id, tier in iter_recovered_ids(msp_path, id_offset=id_offset):
        batch.append((spectrum_id, tier))
        scanned[tier] += 1
        if len(batch) >= 20_000:
            conn.executemany(
                "INSERT OR REPLACE INTO recovered_mona_tier VALUES (?,?)", batch
            )
            batch.clear()
    if batch:
        conn.executemany(
            "INSERT OR REPLACE INTO recovered_mona_tier VALUES (?,?)", batch
        )
    conn.commit()

    matched = dict(conn.execute(
        "SELECT r.tier, COUNT(*) FROM recovered_mona_tier r "
        "JOIN Spectra s ON s.spectrum_id=r.spectrum_id "
        "WHERE s.dataset_source='MoNA-export' GROUP BY r.tier"
    ).fetchall())
    result = {
        "raw_recovered": scanned,
        "database_matched": matched,
        "legacy_id_offset": id_offset,
    }
    if dry_run:
        conn.close()
        result["dry_run"] = True
        return result

    backup = _backup_database(db_path)
    try:
        conn.execute("BEGIN IMMEDIATE")
        conn.execute(
            "CREATE TABLE IF NOT EXISTS LibraryProvenanceMigrations("
            "migration_id TEXT PRIMARY KEY, applied_at TEXT NOT NULL, "
            "backup_path TEXT NOT NULL, details_json TEXT NOT NULL)"
        )
        if conn.execute(
            "SELECT 1 FROM LibraryProvenanceMigrations WHERE migration_id=?",
            (MIGRATION_ID,),
        ).fetchone():
            conn.execute("ROLLBACK")
            backup.unlink(missing_ok=True)
            return {**result, "already_applied": True}

        before_total = int(conn.execute("SELECT COUNT(*) FROM Spectra").fetchone()[0])
        cursor = conn.execute(
            "UPDATE Spectra SET "
            "source_tier=(SELECT tier FROM recovered_mona_tier r "
            "WHERE r.spectrum_id=Spectra.spectrum_id), "
            "quality_flags_json=CASE "
            "WHEN (SELECT tier FROM recovered_mona_tier r WHERE r.spectrum_id=Spectra.spectrum_id)='predicted' "
            "THEN '[\"raw_source_in_silico\"]' "
            "ELSE '[\"raw_source_ms1_not_reference_msms\"]' END "
            "WHERE dataset_source='MoNA-export' AND spectrum_id IN "
            "(SELECT spectrum_id FROM recovered_mona_tier)"
        )
        after_total = int(conn.execute("SELECT COUNT(*) FROM Spectra").fetchone()[0])
        if after_total != before_total:
            raise RuntimeError("spectrum cardinality changed during provenance migration")
        details = {**result, "updated": int(cursor.rowcount)}
        conn.execute(
            "INSERT INTO LibraryProvenanceMigrations VALUES (?,?,?,?)",
            (
                MIGRATION_ID, time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                str(backup), json.dumps(details, sort_keys=True),
            ),
        )
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    finally:
        conn.close()

    verify = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        if verify.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise RuntimeError("post-migration quick_check failed")
        if list(verify.execute("PRAGMA foreign_key_check")):
            raise RuntimeError("post-migration foreign_key_check failed")
        details["after_tiers"] = dict(verify.execute(
            "SELECT source_tier, COUNT(*) FROM Spectra "
            "WHERE dataset_source='MoNA-export' GROUP BY source_tier"
        ).fetchall())
    finally:
        verify.close()
    return {**details, "backup": str(backup)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--msp", type=Path, required=True)
    parser.add_argument(
        "--id-offset", type=int, default=0,
        help="Legacy builder global-index offset for this source file",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    print(json.dumps(migrate(
        args.db, args.msp, dry_run=args.dry_run, id_offset=args.id_offset,
    ), indent=2))


if __name__ == "__main__":
    main()
