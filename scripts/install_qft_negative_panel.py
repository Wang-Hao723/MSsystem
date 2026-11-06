#!/usr/bin/env python3
"""Install the audited Qingfengteng negative-ion candidate panel.

The source Qingfengteng catalogue is a literature-presence list and is almost
entirely recorded as positive ions.  This installer does not manufacture
negative ions for every entry.  It keeps only identities for which the source
catalogue carries a verified full InChIKey and the shared library contains an
identity-matched, experimental [M-H]- reference with a precursor consistent
with the molecular formula.

Only ``TargetCompounds`` rows under the dedicated dataset name are inserted.
The ordinary spectral library and its ranking inputs are not changed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sqlite3
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.expected_compounds import _normalise_adduct  # noqa: E402
from matching.library_quality import (  # noqa: E402
    formula_exact_mass,
    normalize_formula_text,
    normalize_inchikey,
)
from utils.adduct import ion_mz  # noqa: E402


DATASET_SOURCE = "QFT_青风藤_NEG_experimental_reference_20260920"
DATASET_VERSION = "QFT-NEG-panel-v1.0.0"
SOURCE_PLANT = "青风藤"
ALLOWED_LIBRARY_TIERS = {
    "experimental",
    "curated",
    "merged",
    "massbank_supplement",
}
PRECURSOR_TOLERANCE_PPM = 20.0
DEFAULT_DB = Path(
    "/Users/zizi/Desktop/桌面分类/共享与备份/"
    "Mssystem_shared_data/local_library.db"
)
DEFAULT_RECEIPT = (
    PROJECT_ROOT
    / "research/20260920_qft_negative_panel/INSTALL_RECEIPT.json"
)


TARGET_COLUMNS = (
    "target_id",
    "source_file",
    "dataset_source",
    "compound_name",
    "molecular_formula",
    "smiles",
    "inchi",
    "inchikey",
    "exact_mass",
    "precursor_mz",
    "adduct",
    "charge",
    "plant_source",
    "source_tier",
    "quality_flags_json",
    "created_at",
    "reported_inchikey",
    "structure_identity_verified",
    "structure_identity_source",
    "source_content_sha256",
    "dataset_version",
    "source_row_index",
    "provenance",
    "doi",
    "dataset_accession",
    "spectrum_accession",
    "reported_precursor_mz",
    "reported_adduct",
    "reported_retention_time",
    "reported_retention_time_unit",
    "reported_collision_energy_raw",
    "instrument",
    "target_searchability",
)


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _connect(path: Path, *, writable: bool) -> sqlite3.Connection:
    mode = "rw" if writable else "ro"
    connection = sqlite3.connect(f"file:{path.resolve()}?mode={mode}", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys=ON")
    connection.execute("PRAGMA busy_timeout=30000")
    return connection


def _library_reference_priority(row: Mapping[str, Any]) -> tuple[Any, ...]:
    dataset = str(row.get("library_dataset_source") or "").casefold()
    dataset_rank = (
        0 if "ciead" in dataset
        else 1 if "massbank" in dataset
        else 2 if "mona" in dataset
        else 3
    )
    return (
        abs(float(row["precursor_error_ppm"])),
        dataset_rank,
        str(row.get("spectrum_id") or ""),
    )


def build_panel(connection: sqlite3.Connection) -> list[dict[str, Any]]:
    """Return one strict [M-H]- row per verified Qingfengteng identity."""

    target_rows = connection.execute(
        """
        SELECT *
        FROM TargetCompounds
        WHERE source_tier='literature_presence_only'
          AND plant_source LIKE '%青风藤%'
          AND structure_identity_verified=1
          AND trim(coalesce(inchikey,''))<>''
          AND dataset_source<>?
        ORDER BY id
        """,
        (DATASET_SOURCE,),
    ).fetchall()

    # Duplicate source rows occur in the inherited catalogue.  Full InChIKey
    # is the panel identity, so the earliest audited source row is retained.
    source_by_key: dict[str, sqlite3.Row] = {}
    for row in target_rows:
        key = normalize_inchikey(row["inchikey"])
        if key:
            source_by_key.setdefault(key, row)

    reference_rows = connection.execute(
        """
        SELECT
            c.inchikey AS compound_inchikey,
            c.structure_identity_verified AS compound_identity_verified,
            c.molecular_formula AS library_formula,
            s.spectrum_id,
            s.precursor_mz AS reference_precursor_mz,
            s.adduct AS reference_adduct,
            s.ion_mode AS reference_ion_mode,
            s.source_tier AS library_source_tier,
            s.dataset_source AS library_dataset_source,
            s.instrument,
            s.collision_energy_raw
        FROM Compounds c
        JOIN Mapping m ON m.compound_id=c.id
        JOIN Spectra s ON s.id=m.spectrum_id
        WHERE c.structure_identity_verified=1
          AND trim(coalesce(c.inchikey,''))<>''
          AND lower(trim(coalesce(s.source_tier,''))) IN (
              'experimental','curated','merged','massbank_supplement'
          )
        ORDER BY s.id
        """
    ).fetchall()

    references_by_key: dict[str, list[dict[str, Any]]] = {}
    for raw in reference_rows:
        key = normalize_inchikey(raw["compound_inchikey"])
        source = source_by_key.get(key)
        if source is None:
            continue
        adduct = _normalise_adduct(raw["reference_adduct"])
        if adduct != "[M-H]-":
            continue
        tier = str(raw["library_source_tier"] or "").strip().casefold()
        if tier not in ALLOWED_LIBRARY_TIERS:
            continue
        source_formula = normalize_formula_text(source["molecular_formula"])
        library_formula = normalize_formula_text(raw["library_formula"])
        if not source_formula or source_formula != library_formula:
            continue
        neutral_mass = formula_exact_mass(source_formula)
        theoretical_mz = ion_mz(neutral_mass, adduct) if neutral_mass else None
        reference_mz = raw["reference_precursor_mz"]
        if theoretical_mz is None or reference_mz is None:
            continue
        error_ppm = (float(reference_mz) - theoretical_mz) / theoretical_mz * 1e6
        if abs(error_ppm) > PRECURSOR_TOLERANCE_PPM:
            continue
        record = dict(raw)
        record.update(
            theoretical_precursor_mz=theoretical_mz,
            precursor_error_ppm=error_ppm,
        )
        references_by_key.setdefault(key, []).append(record)

    created_at = datetime.now(timezone.utc).isoformat()
    panel: list[dict[str, Any]] = []
    for row_index, key in enumerate(sorted(references_by_key), start=1):
        source = source_by_key[key]
        reference = min(
            references_by_key[key], key=_library_reference_priority
        )
        provenance_payload = {
            "panel_rule": (
                "verified Qingfengteng literature identity intersected with "
                "an identity-matched experimental [M-H]- reference; formula "
                "and precursor required to agree within 20 ppm"
            ),
            "source_target": {
                "row_id": int(source["id"]),
                "target_id": str(source["target_id"] or ""),
                "dataset_source": str(source["dataset_source"] or ""),
                "source_file": str(source["source_file"] or ""),
                "source_content_sha256": str(
                    source["source_content_sha256"] or ""
                ),
            },
            "negative_reference": {
                "spectrum_id": str(reference["spectrum_id"] or ""),
                "dataset_source": str(
                    reference["library_dataset_source"] or ""
                ),
                "source_tier": str(reference["library_source_tier"] or ""),
                "reported_precursor_mz": float(
                    reference["reference_precursor_mz"]
                ),
                "theoretical_precursor_mz": float(
                    reference["theoretical_precursor_mz"]
                ),
                "precursor_error_ppm": float(
                    reference["precursor_error_ppm"]
                ),
                "adduct": "[M-H]-",
            },
            "claim_limit": (
                "candidate panel only; not an authentic-standard marker "
                "confirmation or a complete Qingfengteng quality panel"
            ),
        }
        source_digest = _sha256_json(provenance_payload)
        identity_token = hashlib.sha256(key.encode("ascii")).hexdigest()[:20]
        panel.append(
            {
                "target_id": f"QFTNEG-{identity_token}",
                "source_file": str(source["source_file"] or ""),
                "dataset_source": DATASET_SOURCE,
                "compound_name": str(source["compound_name"] or "").strip(),
                "molecular_formula": normalize_formula_text(
                    source["molecular_formula"]
                ),
                "smiles": str(source["smiles"] or ""),
                "inchi": str(source["inchi"] or ""),
                "inchikey": key,
                "exact_mass": formula_exact_mass(source["molecular_formula"]),
                "precursor_mz": float(reference["theoretical_precursor_mz"]),
                "adduct": "[M-H]-",
                "charge": -1,
                "plant_source": SOURCE_PLANT,
                "source_tier": "literature_presence_only",
                "quality_flags_json": _canonical_json([
                    "candidate_panel_not_authentic_standard_confirmed",
                    "literature_presence_plus_experimental_negative_reference",
                    "not_a_complete_qingfengteng_quality_marker_panel",
                ]),
                "created_at": created_at,
                "reported_inchikey": str(source["reported_inchikey"] or ""),
                "structure_identity_verified": 1,
                "structure_identity_source": (
                    "verified_source_target_and_identity_matched_"
                    "experimental_negative_reference"
                ),
                "source_content_sha256": source_digest,
                "dataset_version": DATASET_VERSION,
                "source_row_index": row_index,
                "provenance": _canonical_json(provenance_payload),
                "doi": str(source["doi"] or ""),
                "dataset_accession": DATASET_SOURCE,
                "spectrum_accession": str(reference["spectrum_id"] or ""),
                "reported_precursor_mz": float(
                    reference["reference_precursor_mz"]
                ),
                "reported_adduct": "[M-H]-",
                "reported_retention_time": None,
                "reported_retention_time_unit": None,
                "reported_collision_energy_raw": str(
                    reference["collision_energy_raw"] or ""
                ),
                "instrument": str(reference["instrument"] or ""),
                "target_searchability": (
                    "MS1_EXACT_MASS_AND_EXPERIMENTAL_LIBRARY_CANDIDATE"
                ),
            }
        )
    return panel


def _panel_projection(rows: list[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [
        {column: row.get(column) for column in TARGET_COLUMNS}
        for row in rows
    ]


def _existing_panel(connection: sqlite3.Connection) -> list[dict[str, Any]]:
    rows = connection.execute(
        f"SELECT {','.join(TARGET_COLUMNS)} FROM TargetCompounds "
        "WHERE dataset_source=? ORDER BY source_row_index,target_id",
        (DATASET_SOURCE,),
    ).fetchall()
    return [dict(row) for row in rows]


def install_panel(
    db_path: Path,
    receipt_path: Path,
    *,
    apply: bool,
) -> dict[str, Any]:
    if apply and os.environ.get("MSSYSTEM_LIBRARY_WRITE_ENABLED", "").lower() not in {
        "1", "true", "yes", "on"
    }:
        raise RuntimeError(
            "apply requires MSSYSTEM_LIBRARY_WRITE_ENABLED=1"
        )
    connection = _connect(db_path, writable=apply)
    try:
        panel = build_panel(connection)
        if not panel:
            raise RuntimeError("strict intersection produced an empty panel")
        names = {str(row["compound_name"]).casefold() for row in panel}
        if "sinomenine" not in names:
            raise RuntimeError("strict panel does not contain sinomenine")
        if len(panel) != len({row["inchikey"] for row in panel}):
            raise RuntimeError("panel contains duplicate full InChIKeys")
        existing = _existing_panel(connection)
        panel_fingerprint = _sha256_json(_panel_projection(panel))
        existing_fingerprint = (
            _sha256_json(_panel_projection(existing)) if existing else None
        )
        status = "preview"
        if apply:
            if existing:
                # created_at is deliberately not part of semantic identity;
                # require every other field to match for an idempotent no-op.
                comparable_columns = tuple(
                    column for column in TARGET_COLUMNS if column != "created_at"
                )
                expected_semantic = [
                    {column: row.get(column) for column in comparable_columns}
                    for row in panel
                ]
                existing_semantic = [
                    {column: row.get(column) for column in comparable_columns}
                    for row in existing
                ]
                if _sha256_json(expected_semantic) != _sha256_json(existing_semantic):
                    raise RuntimeError(
                        "existing QFT negative panel differs; refusing to replace it"
                    )
                status = "already_installed_exact"
            else:
                placeholders = ",".join("?" for _ in TARGET_COLUMNS)
                columns = ",".join(TARGET_COLUMNS)
                connection.execute("BEGIN IMMEDIATE")
                try:
                    connection.executemany(
                        f"INSERT INTO TargetCompounds ({columns}) "
                        f"VALUES ({placeholders})",
                        [tuple(row[column] for column in TARGET_COLUMNS) for row in panel],
                    )
                    inserted = connection.execute(
                        "SELECT COUNT(*) FROM TargetCompounds "
                        "WHERE dataset_source=?",
                        (DATASET_SOURCE,),
                    ).fetchone()[0]
                    if int(inserted) != len(panel):
                        raise RuntimeError(
                            f"insert verification failed: {inserted} != {len(panel)}"
                        )
                    if connection.execute("PRAGMA foreign_key_check").fetchone():
                        raise RuntimeError("foreign_key_check failed")
                    connection.commit()
                except Exception:
                    connection.rollback()
                    raise
                status = "installed"

        payload = {
            "schema_version": "mssystem-qft-negative-panel-install/1.0",
            "status": status,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "database_path": str(db_path.resolve()),
            "dataset_source": DATASET_SOURCE,
            "dataset_version": DATASET_VERSION,
            "panel_row_count": len(panel),
            "unique_full_inchikey_count": len({row["inchikey"] for row in panel}),
            "all_adducts": sorted({row["adduct"] for row in panel}),
            "compound_names": sorted(row["compound_name"] for row in panel),
            "contains_sinomenine": "sinomenine" in names,
            "panel_semantic_sha256": panel_fingerprint,
            "preexisting_panel_semantic_sha256": existing_fingerprint,
            "selection_rule": (
                "verified Qingfengteng literature identities with an exact-"
                "identity experimental [M-H]- reference and <=20 ppm formula-"
                "precursor agreement"
            ),
            "claim_limit": (
                "candidate screening panel; no authentic-standard marker "
                "confirmation and no claim of complete quality evaluation"
            ),
            "records": panel,
        }
        _write_json_atomic(receipt_path, payload)
        return payload
    finally:
        connection.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--apply", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.db.is_file():
        raise SystemExit(f"database does not exist: {args.db}")
    payload = install_panel(
        args.db.expanduser().resolve(),
        args.receipt.expanduser().resolve(),
        apply=args.apply,
    )
    print(json.dumps({
        key: payload[key]
        for key in (
            "status",
            "dataset_source",
            "panel_row_count",
            "unique_full_inchikey_count",
            "compound_names",
            "contains_sinomenine",
            "panel_semantic_sha256",
        )
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
