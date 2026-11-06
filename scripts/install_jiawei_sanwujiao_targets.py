#!/usr/bin/env python3
"""Install a reviewed Jiawei Sanwujiao presence-only target catalogue.

The input is a literature catalogue, not a spectral library.  Every declared
ion becomes one ``TargetCompounds`` row with ``source_tier`` fixed to
``literature_presence_only``.  The command never writes spectra, diagnostic
ions, retention times, or MS/MS evidence.  Dry-run is the default; ``--apply``
also requires ``MSSYSTEM_LIBRARY_WRITE_ENABLED=1``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sqlite3
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.expected_compounds import _normalise_adduct  # noqa: E402
from matching.library_quality import (  # noqa: E402
    formula_exact_mass,
    inchi_metadata,
    normalize_formula_text,
    normalize_inchikey,
    smiles_metadata,
)
from utils.adduct import adduct_polarity, ion_mz, parse_adduct_shift  # noqa: E402


CATALOG_SCHEMA = "mssystem-jiawei-sanwujiao-target-catalog/1.0"
RECEIPT_SCHEMA = "mssystem-presence-target-install/1.0"
EXPECTED_PLANT_SOURCE = "加味三乌胶"
PRESENCE_TIER = "literature_presence_only"
NEUTRAL_MASS_TOLERANCE_PPM = 5.0
PRECURSOR_TOLERANCE_PPM = 20.0
ELECTRON_MASS_DA = 0.000548579909065
DEFAULT_CATALOG = (
    PROJECT_ROOT.parent
    / "research/20261020_jiawei_sanwujiao_targeted/catalog.json"
)
DEFAULT_DB = Path(
    "/Users/zizi/Desktop/桌面分类/共享与备份/"
    "Mssystem_shared_data/local_library.db"
)
DEFAULT_RECEIPT = (
    PROJECT_ROOT.parent
    / "research/20261020_jiawei_sanwujiao_targeted/install_receipt.json"
)


TARGET_COLUMNS = (
    "target_id",
    "source_file",
    "source_content_sha256",
    "dataset_version",
    "source_row_index",
    "dataset_source",
    "compound_name",
    "molecular_formula",
    "smiles",
    "inchi",
    "inchikey",
    "reported_inchikey",
    "structure_identity_verified",
    "structure_identity_source",
    "exact_mass",
    "precursor_mz",
    "adduct",
    "charge",
    "reported_precursor_mz",
    "reported_adduct",
    "reported_retention_time",
    "reported_retention_time_unit",
    "reported_collision_energy_raw",
    "instrument",
    "target_searchability",
    "plant_source",
    "provenance",
    "doi",
    "dataset_accession",
    "spectrum_accession",
    "source_tier",
    "quality_flags_json",
    "created_at",
)

FORBIDDEN_EVIDENCE_KEYS = frozenset({
    "peaks",
    "fragments",
    "diagnostic_ions",
    "ms2",
    "ms2_spectrum",
    "retention_time",
    "reported_retention_time",
    "collision_energy",
    "instrument",
    "spectrum_accession",
})


class CatalogError(ValueError):
    """Raised when the reviewed JSON does not satisfy the target contract."""


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_json(value: Any) -> str:
    return _sha256_bytes(_canonical_json(value).encode("utf-8"))


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _required_text(container: Mapping[str, Any], key: str, *, where: str) -> str:
    value = str(container.get(key) or "").strip()
    if not value:
        raise CatalogError(f"{where}.{key} is required")
    return value


def _finite_positive(value: Any, *, where: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise CatalogError(f"{where} must be a finite positive number") from exc
    if not math.isfinite(number) or number <= 0:
        raise CatalogError(f"{where} must be a finite positive number")
    return number


def _ppm_error(observed: float, expected: float) -> float:
    return (observed - expected) / expected * 1e6


def _validate_no_measured_evidence(
    value: Mapping[str, Any], *, where: str
) -> None:
    present = sorted(FORBIDDEN_EVIDENCE_KEYS.intersection(value))
    if present:
        raise CatalogError(
            f"{where} contains unsupported measured-evidence fields: "
            + ", ".join(present)
        )


def _structure_fields(
    compound: Mapping[str, Any], formula: str, formal_charge: int, *, where: str
) -> dict[str, Any]:
    smiles = str(compound.get("smiles") or "").strip()
    inchi = str(compound.get("inchi") or "").strip()
    reported_key = normalize_inchikey(compound.get("inchikey"))
    verified = compound.get("structure_identity_verified", False)
    if not isinstance(verified, bool):
        raise CatalogError(f"{where}.structure_identity_verified must be boolean")

    metadata: list[tuple[str, dict[str, Any]]] = []
    if smiles:
        metadata.append(("smiles", smiles_metadata(smiles)))
    if inchi:
        metadata.append(("inchi", inchi_metadata(inchi)))
    if (smiles or inchi or reported_key or verified) and not metadata:
        raise CatalogError(
            f"{where} needs SMILES or InChI to validate a supplied structure identity"
        )
    for label, details in metadata:
        if details.get("status") != "valid":
            raise CatalogError(f"{where}.{label} is not a valid RDKit structure")
        raw_formula = str(details.get("formula") or "")
        derived_formula = normalize_formula_text(
            re.sub(r"[+-]\d*$", "", raw_formula)
        )
        if derived_formula != formula:
            raise CatalogError(
                f"{where}.{label} formula {derived_formula!r} differs from {formula!r}"
            )
    if metadata:
        try:
            from rdkit import Chem

            molecule = (
                Chem.MolFromSmiles(smiles)
                if smiles else Chem.MolFromInchi(inchi)
            )
            derived_charge = sum(
                int(atom.GetFormalCharge()) for atom in molecule.GetAtoms()
            ) if molecule is not None else None
        except Exception as exc:
            raise CatalogError(f"{where} formal charge could not be validated") from exc
        if derived_charge != formal_charge:
            raise CatalogError(
                f"{where}.formal_charge {formal_charge} differs from structure "
                f"charge {derived_charge}"
            )

    keys = {
        normalize_inchikey(details.get("inchikey"))
        for _, details in metadata
        if normalize_inchikey(details.get("inchikey"))
    }
    if reported_key:
        keys.add(reported_key)
    if len(keys) > 1:
        raise CatalogError(f"{where} structure representations have different InChIKeys")
    canonical_key = next(iter(keys), "")
    if compound.get("inchikey") and not canonical_key:
        raise CatalogError(f"{where}.inchikey is malformed")
    if verified and not canonical_key:
        raise CatalogError(f"{where} verified structure has no validated InChIKey")

    identity_source = str(compound.get("structure_identity_source") or "").strip()
    if verified and not identity_source:
        raise CatalogError(
            f"{where}.structure_identity_source is required for a verified structure"
        )
    primary = metadata[0][1] if metadata else {}
    return {
        "smiles": str(primary.get("canonical_smiles") or smiles),
        "inchi": str(primary.get("inchi") or inchi),
        "inchikey": canonical_key,
        "reported_inchikey": reported_key,
        "structure_identity_verified": int(verified),
        "structure_identity_source": identity_source,
        "formal_charge": formal_charge,
    }


def _theoretical_precursor_mz(
    exact_mass: float, adduct: str, formal_charge: int
) -> float | None:
    if formal_charge:
        expected = "[M]+" if formal_charge == 1 else "[M]-"
        if adduct != expected:
            return None
        return exact_mass - formal_charge * ELECTRON_MASS_DA
    if adduct in {"[M]+", "[M]-"}:
        return None
    return ion_mz(exact_mass, adduct)


def _source_metadata(payload: Mapping[str, Any]) -> dict[str, str]:
    source = payload.get("source")
    if not isinstance(source, Mapping):
        raise CatalogError("catalog.source must be an object")
    title = _required_text(source, "title", where="catalog.source")
    citation = str(source.get("citation") or "").strip()
    doi = str(source.get("doi") or "").strip()
    url = str(source.get("url") or "").strip()
    if not any((citation, doi, url)):
        raise CatalogError("catalog.source needs citation, doi, or url")
    return {"title": title, "citation": citation, "doi": doi, "url": url}


def _identity_token(
    *,
    inchikey: str,
    name: str,
    formula: str,
    source_row_index: int,
    plant_source: str,
) -> str:
    if inchikey:
        return _canonical_json({
            "inchikey": inchikey,
            "plant_source": plant_source,
        })
    return _canonical_json({
        "name": name.casefold(),
        "formula": formula,
        "source_row_index": source_row_index,
        "plant_source": plant_source,
    })


def load_catalog(catalog_path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Validate the JSON and expand each declared ion into a target row."""

    raw = catalog_path.read_bytes()
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise CatalogError(f"catalog is not valid UTF-8 JSON: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise CatalogError("catalog must be a JSON object")
    if payload.get("schema_version") != CATALOG_SCHEMA:
        raise CatalogError(
            f"catalog.schema_version must equal {CATALOG_SCHEMA!r}"
        )
    for field in (
        "raw_spectra_created",
        "retention_times_created",
        "diagnostic_ions_created",
    ):
        if payload.get(field) not in (None, 0, False):
            raise CatalogError(f"catalog.{field} must be zero/false for this target-only install")
    dataset_source = _required_text(payload, "dataset_source", where="catalog")
    dataset_version = _required_text(payload, "dataset_version", where="catalog")
    plant_source = _required_text(payload, "plant_source", where="catalog")
    if plant_source != EXPECTED_PLANT_SOURCE:
        raise CatalogError(
            f"catalog.plant_source must equal {EXPECTED_PLANT_SOURCE!r}"
        )
    if "source_tier" in payload and payload.get("source_tier") != PRESENCE_TIER:
        raise CatalogError(f"catalog.source_tier must equal {PRESENCE_TIER!r}")
    source = _source_metadata(payload)
    compounds = payload.get("compounds")
    if not isinstance(compounds, list) or not compounds:
        raise CatalogError("catalog.compounds must be a non-empty array")
    source_sha = _sha256_bytes(raw)
    created_at = datetime.now(timezone.utc).isoformat()

    rows: list[dict[str, Any]] = []
    seen_row_indexes: set[int] = set()
    for offset, compound in enumerate(compounds, start=1):
        where = f"catalog.compounds[{offset - 1}]"
        if not isinstance(compound, Mapping):
            raise CatalogError(f"{where} must be an object")
        _validate_no_measured_evidence(compound, where=where)
        name = _required_text(compound, "compound_name", where=where)
        formula = normalize_formula_text(
            _required_text(compound, "molecular_formula", where=where)
        )
        calculated_mass = formula_exact_mass(formula)
        if calculated_mass is None:
            raise CatalogError(f"{where}.molecular_formula is unsupported or invalid")
        stated_mass = _finite_positive(compound.get("exact_mass"), where=f"{where}.exact_mass")
        if abs(_ppm_error(stated_mass, calculated_mass)) > NEUTRAL_MASS_TOLERANCE_PPM:
            raise CatalogError(
                f"{where}.exact_mass differs from formula mass by more than "
                f"{NEUTRAL_MASS_TOLERANCE_PPM:g} ppm"
            )
        row_index = compound.get("source_row_index")
        if isinstance(row_index, bool) or not isinstance(row_index, int) or row_index < 1:
            raise CatalogError(f"{where}.source_row_index must be a positive integer")
        if row_index in seen_row_indexes:
            raise CatalogError(f"duplicate source_row_index: {row_index}")
        seen_row_indexes.add(row_index)
        provenance = compound.get("provenance")
        if not isinstance(provenance, (str, Mapping)) or not provenance:
            raise CatalogError(f"{where}.provenance must be a non-empty string or object")
        compound_plant_source = _required_text(
            compound, "plant_source", where=where
        )
        formal_charge = compound.get("formal_charge", 0)
        if (
            isinstance(formal_charge, bool)
            or not isinstance(formal_charge, int)
            or formal_charge not in {-1, 0, 1}
        ):
            raise CatalogError(f"{where}.formal_charge must be -1, 0, or 1")
        structure = _structure_fields(
            compound, formula, formal_charge, where=where
        )
        ions = compound.get("ions")
        if not isinstance(ions, list) or not ions:
            raise CatalogError(f"{where}.ions must be a non-empty array")
        identity = _identity_token(
            inchikey=structure["inchikey"],
            name=name,
            formula=formula,
            source_row_index=row_index,
            plant_source=compound_plant_source,
        )
        for ion_offset, ion in enumerate(ions, start=1):
            ion_where = f"{where}.ions[{ion_offset - 1}]"
            if not isinstance(ion, Mapping):
                raise CatalogError(f"{ion_where} must be an object")
            _validate_no_measured_evidence(ion, where=ion_where)
            raw_adduct = _required_text(ion, "adduct", where=ion_where)
            adduct = _normalise_adduct(raw_adduct)
            if not adduct or parse_adduct_shift(adduct) is None:
                raise CatalogError(f"{ion_where}.adduct is unsupported: {raw_adduct!r}")
            charge = ion.get("charge")
            if isinstance(charge, bool) or not isinstance(charge, int) or charge == 0:
                raise CatalogError(f"{ion_where}.charge must be a non-zero integer")
            parsed = parse_adduct_shift(adduct)
            assert parsed is not None
            _, magnitude = parsed
            polarity = adduct_polarity(adduct)
            expected_charge = abs(int(magnitude)) * (-1 if polarity == "NEG" else 1)
            if charge != expected_charge:
                raise CatalogError(
                    f"{ion_where}.charge {charge} conflicts with {adduct}"
                )
            precursor = _finite_positive(
                ion.get("precursor_mz"), where=f"{ion_where}.precursor_mz"
            )
            theoretical = _theoretical_precursor_mz(
                calculated_mass, adduct, formal_charge
            )
            if theoretical is None:
                raise CatalogError(
                    f"{ion_where}.adduct is incompatible with formal_charge "
                    f"{formal_charge}"
                )
            precursor_error = _ppm_error(precursor, theoretical)
            if abs(precursor_error) > PRECURSOR_TOLERANCE_PPM:
                raise CatalogError(
                    f"{ion_where}.precursor_mz differs from formula/adduct by more than "
                    f"{PRECURSOR_TOLERANCE_PPM:g} ppm"
                )
            token = _canonical_json({
                "dataset_source": dataset_source,
                "identity": identity,
                "adduct": adduct,
                "charge": charge,
            })
            target_id = "JWSWJ-" + hashlib.sha256(
                token.encode("utf-8")
            ).hexdigest()[:24]
            row_provenance = {
                "catalog_source": source,
                "compound_provenance": provenance,
                "evidence_scope": "literature_presence_and_exact_mass_candidate_only",
                "theoretical_precursor_mz": theoretical,
                "catalog_declared_precursor_mz": precursor,
                "precursor_error_ppm": precursor_error,
                "no_measured_ms2_or_retention_time": True,
            }
            quality_flags = [
                "literature_presence_only",
                "exact_mass_candidate_not_identification",
                "calculated_ion_hypothesis",
                "no_measured_ms2",
                "no_measured_retention_time",
            ]
            if formal_charge:
                quality_flags.append(f"formal_charge:{formal_charge}")
            rows.append({
                "target_id": target_id,
                "source_file": str(catalog_path.resolve()),
                "source_content_sha256": source_sha,
                "dataset_version": dataset_version,
                "source_row_index": row_index,
                "dataset_source": dataset_source,
                "compound_name": name,
                "molecular_formula": formula,
                "smiles": structure["smiles"],
                "inchi": structure["inchi"],
                "inchikey": structure["inchikey"],
                "reported_inchikey": structure["reported_inchikey"],
                "structure_identity_verified": structure["structure_identity_verified"],
                "structure_identity_source": structure["structure_identity_source"],
                "exact_mass": calculated_mass,
                "precursor_mz": theoretical,
                "adduct": adduct,
                "charge": charge,
                "reported_precursor_mz": None,
                "reported_adduct": None,
                "reported_retention_time": None,
                "reported_retention_time_unit": None,
                "reported_collision_energy_raw": None,
                "instrument": None,
                "target_searchability": "MS1_EXACT_MASS_CANDIDATE",
                "plant_source": compound_plant_source,
                "provenance": _canonical_json(row_provenance),
                "doi": str(compound.get("doi") or source["doi"]),
                "dataset_accession": dataset_source,
                "spectrum_accession": None,
                "source_tier": PRESENCE_TIER,
                "quality_flags_json": _canonical_json(quality_flags),
                "created_at": created_at,
            })

    target_ids = [row["target_id"] for row in rows]
    if len(target_ids) != len(set(target_ids)):
        raise CatalogError(
            "catalog contains duplicate compound/adduct identities; merge exact duplicates"
        )
    catalog_info = {
        "dataset_source": dataset_source,
        "dataset_version": dataset_version,
        "plant_source": plant_source,
        "source": source,
        "source_content_sha256": source_sha,
        "compound_count": len(compounds),
        "target_row_count": len(rows),
    }
    return catalog_info, rows


def _connect(path: Path, *, writable: bool) -> sqlite3.Connection:
    mode = "rw" if writable else "ro"
    connection = sqlite3.connect(f"file:{path.resolve()}?mode={mode}", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys=ON")
    connection.execute("PRAGMA busy_timeout=30000")
    if not writable:
        connection.execute("PRAGMA query_only=ON")
    return connection


def _require_schema(connection: sqlite3.Connection) -> None:
    exists = connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='TargetCompounds'"
    ).fetchone()
    if not exists:
        raise RuntimeError("database has no TargetCompounds table")
    columns = {
        str(row[1]) for row in connection.execute("PRAGMA table_info(TargetCompounds)")
    }
    missing = sorted(set(TARGET_COLUMNS).difference(columns))
    if missing:
        raise RuntimeError("TargetCompounds is missing columns: " + ", ".join(missing))


def _database_scope_snapshot(
    connection: sqlite3.Connection, db_path: Path, dataset_source: str
) -> dict[str, Any]:
    """Capture the small proof needed for a later target-only cache rebind."""

    tables = {
        str(row[0])
        for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )
    }
    counts = {
        table: int(connection.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0])
        for table in (
            "Spectra", "Compounds", "Mapping", "TargetCompounds", "DiagnosticIons"
        )
        if table in tables
    }
    schema_rows = [
        list(row)
        for row in connection.execute(
            "SELECT type,name,tbl_name,sql FROM sqlite_master "
            "WHERE name NOT LIKE 'sqlite_%' ORDER BY type,name"
        )
    ]
    digest = hashlib.sha256()
    selected = ",".join(("id", *TARGET_COLUMNS))
    for row in connection.execute(
        f"SELECT {selected} FROM TargetCompounds "
        "WHERE dataset_source IS NULL OR dataset_source<>? ORDER BY id",
        (dataset_source,),
    ):
        digest.update(_canonical_json(list(row)).encode("utf-8"))
        digest.update(b"\n")
    stat = db_path.stat()
    return {
        "file_stat": {
            "path": str(db_path.resolve()),
            "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "ctime_ns": stat.st_ctime_ns,
        },
        "schema_sha256": _sha256_json(schema_rows),
        "schema_version": int(connection.execute("PRAGMA schema_version").fetchone()[0]),
        "table_row_counts": counts,
        "non_dataset_target_rows": counts.get("TargetCompounds", 0)
        - int(connection.execute(
            "SELECT COUNT(*) FROM TargetCompounds WHERE dataset_source=?",
            (dataset_source,),
        ).fetchone()[0]),
        "non_dataset_targets_sha256": digest.hexdigest(),
    }


def _semantic_projection(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    columns = tuple(column for column in TARGET_COLUMNS if column != "created_at")
    return [
        {column: row.get(column) for column in columns}
        for row in sorted(rows, key=lambda item: str(item.get("target_id") or ""))
    ]


def _existing_dataset(
    connection: sqlite3.Connection, dataset_source: str
) -> list[dict[str, Any]]:
    rows = connection.execute(
        f"SELECT {','.join(TARGET_COLUMNS)} FROM TargetCompounds "
        "WHERE dataset_source=? ORDER BY target_id",
        (dataset_source,),
    ).fetchall()
    return [dict(row) for row in rows]


def _target_only_authorizer(
    action: int, arg1: str | None, arg2: str | None, _database: str | None,
    _trigger: str | None,
) -> int:
    write_actions = {
        sqlite3.SQLITE_INSERT,
        sqlite3.SQLITE_UPDATE,
        sqlite3.SQLITE_DELETE,
    }
    ddl_actions = {
        value
        for name in (
            "SQLITE_ALTER_TABLE", "SQLITE_CREATE_INDEX", "SQLITE_CREATE_TABLE",
            "SQLITE_CREATE_TRIGGER", "SQLITE_CREATE_VIEW", "SQLITE_DROP_INDEX",
            "SQLITE_DROP_TABLE", "SQLITE_DROP_TRIGGER", "SQLITE_DROP_VIEW",
        )
        if (value := getattr(sqlite3, name, None)) is not None
    }
    if action in ddl_actions:
        return sqlite3.SQLITE_DENY
    if action in write_actions and arg1 != "TargetCompounds":
        return sqlite3.SQLITE_DENY
    return sqlite3.SQLITE_OK


def install_catalog(
    catalog_path: Path,
    db_path: Path,
    receipt_path: Path,
    *,
    apply: bool,
) -> dict[str, Any]:
    if apply and os.environ.get("MSSYSTEM_LIBRARY_WRITE_ENABLED", "").lower() not in {
        "1", "true", "yes", "on"
    }:
        raise RuntimeError("apply requires MSSYSTEM_LIBRARY_WRITE_ENABLED=1")
    catalog_info, rows = load_catalog(catalog_path)
    connection = _connect(db_path, writable=apply)
    try:
        _require_schema(connection)
        dataset_source = catalog_info["dataset_source"]
        existing = _existing_dataset(connection, dataset_source)
        pre_database_state = _database_scope_snapshot(
            connection, db_path, dataset_source
        )
        expected_semantic = _semantic_projection(rows)
        expected_sha = _sha256_json(expected_semantic)
        existing_sha = _sha256_json(_semantic_projection(existing)) if existing else None
        if existing and existing_sha != expected_sha:
            raise RuntimeError(
                "target dataset already exists with different content; refusing to modify it"
            )
        status = "already_installed_exact" if existing else "preview"
        before_total = int(
            connection.execute("SELECT COUNT(*) FROM TargetCompounds").fetchone()[0]
        )
        inserted = 0
        if apply and not existing:
            conflicting_ids = connection.execute(
                f"SELECT target_id,dataset_source FROM TargetCompounds WHERE target_id IN "
                f"({','.join('?' for _ in rows)})",
                tuple(row["target_id"] for row in rows),
            ).fetchall()
            if conflicting_ids:
                raise RuntimeError("stable target_id collides with an existing dataset")
            prepared_receipt = {
                "schema_version": RECEIPT_SCHEMA,
                "status": "apply_prepared",
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "catalog_path": str(catalog_path.resolve()),
                "database_path": str(db_path.resolve()),
                **catalog_info,
                "target_semantic_sha256": expected_sha,
                "pre_database_state": pre_database_state,
                "target_ids": [row["target_id"] for row in rows],
            }
            _write_json_atomic(receipt_path, prepared_receipt)
            connection.set_authorizer(_target_only_authorizer)
            connection.execute("BEGIN IMMEDIATE")
            try:
                placeholders = ",".join("?" for _ in TARGET_COLUMNS)
                connection.executemany(
                    f"INSERT INTO TargetCompounds ({','.join(TARGET_COLUMNS)}) "
                    f"VALUES ({placeholders})",
                    [tuple(row[column] for column in TARGET_COLUMNS) for row in rows],
                )
                installed = _existing_dataset(connection, dataset_source)
                if _sha256_json(_semantic_projection(installed)) != expected_sha:
                    raise RuntimeError("post-insert semantic verification failed")
                connection.commit()
            except Exception:
                connection.rollback()
                raise
            inserted = len(rows)
            status = "installed"
        after_total = int(
            connection.execute("SELECT COUNT(*) FROM TargetCompounds").fetchone()[0]
        )
        if after_total - before_total != inserted:
            raise RuntimeError("unexpected TargetCompounds row-count change")
        post_database_state = _database_scope_snapshot(
            connection, db_path, dataset_source
        )
        payload = {
            "schema_version": RECEIPT_SCHEMA,
            "status": status,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "catalog_path": str(catalog_path.resolve()),
            "database_path": str(db_path.resolve()),
            **catalog_info,
            "target_semantic_sha256": expected_sha,
            "preexisting_target_semantic_sha256": existing_sha,
            "target_rows_before": before_total,
            "target_rows_after": after_total,
            "inserted_target_rows": inserted,
            "pre_database_state": pre_database_state,
            "post_database_state": post_database_state,
            "only_target_compounds_written": True,
            "diagnostic_ions_inserted": 0,
            "spectra_inserted": 0,
            "claim_scope": "mass-matching candidates only; no measured MS2 or RT",
            "target_ids": [row["target_id"] for row in rows],
        }
        _write_json_atomic(receipt_path, payload)
        return payload
    finally:
        connection.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--apply", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for label, path in (("catalog", args.catalog), ("database", args.db)):
        if not path.is_file():
            raise SystemExit(f"{label} does not exist: {path}")
    payload = install_catalog(
        args.catalog.expanduser().resolve(),
        args.db.expanduser().resolve(),
        args.receipt.expanduser().resolve(),
        apply=args.apply,
    )
    print(json.dumps({
        key: payload[key]
        for key in (
            "status", "dataset_source", "dataset_version", "compound_count",
            "target_row_count", "inserted_target_rows", "target_semantic_sha256",
        )
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
