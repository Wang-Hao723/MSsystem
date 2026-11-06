from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from scripts import install_jiawei_sanwujiao_targets as installer
from scripts import rebind_record_store_after_jiawei_target_only_change as rebind_script
from matching.library_record_store import (
    processing_identity,
    semantic_sha256,
    source_binding,
)
from utils.adduct import ion_mz


def _create_database(path: Path) -> None:
    text_columns = set(installer.TARGET_COLUMNS) - {
        "source_row_index",
        "structure_identity_verified",
        "exact_mass",
        "precursor_mz",
        "charge",
        "reported_precursor_mz",
        "reported_retention_time",
    }
    integer_columns = {"source_row_index", "structure_identity_verified", "charge"}
    definitions = ["id INTEGER PRIMARY KEY"]
    for column in installer.TARGET_COLUMNS:
        if column in text_columns:
            sql_type = "TEXT"
        elif column in integer_columns:
            sql_type = "INTEGER"
        else:
            sql_type = "REAL"
        unique = " UNIQUE" if column == "target_id" else ""
        definitions.append(f"{column} {sql_type}{unique}")
    connection = sqlite3.connect(path)
    try:
        connection.execute(f"CREATE TABLE TargetCompounds ({','.join(definitions)})")
        connection.execute("CREATE TABLE DiagnosticIons (id INTEGER PRIMARY KEY, target_id INTEGER, fragment_mz REAL)")
        connection.execute("CREATE TABLE Spectra (id INTEGER PRIMARY KEY, marker TEXT)")
        connection.execute("CREATE TABLE Compounds (id INTEGER PRIMARY KEY, marker TEXT)")
        connection.execute("CREATE TABLE Mapping (id INTEGER PRIMARY KEY, marker TEXT)")
        connection.execute("INSERT INTO Spectra(marker) VALUES ('unchanged')")
        connection.commit()
    finally:
        connection.close()


def _compound(smiles: str, row: int, plant_source: str) -> dict:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    mass = float(Descriptors.ExactMolWt(molecule))
    return {
        "compound_name": f"isomer-{row}",
        "molecular_formula": formula,
        "exact_mass": mass,
        "smiles": smiles,
        "inchi": Chem.MolToInchi(molecule),
        "inchikey": Chem.MolToInchiKey(molecule),
        "structure_identity_verified": True,
        "structure_identity_source": "reviewed PubChem identity",
        "source_row_index": row,
        "plant_source": plant_source,
        "provenance": {"citation": f"reviewed source row {row}"},
        "ions": [{
            "adduct": "[M+H]+",
            "charge": 1,
            "precursor_mz": ion_mz(mass, "[M+H]+"),
        }],
    }


def _catalog() -> dict:
    # The two structures are formula isomers.  They must remain two targets.
    return {
        "schema_version": installer.CATALOG_SCHEMA,
        "dataset_source": "JWSWJ_加味三乌胶_targeted_20261020",
        "dataset_version": "1.0.0",
        "plant_source": "加味三乌胶",
        "source_tier": "literature_presence_only",
        "source": {
            "title": "Reviewed formula ingredient catalogue",
            "citation": "Local reviewed source",
        },
        "compounds": [
            _compound("CN1C=NC2=C1C(=O)NC(=O)N2C", 1, "川乌"),
            _compound("CN1C2=C(C(=O)N(C1=O)C)NC=N2", 2, "草乌"),
        ],
    }


def _write_catalog(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")


def test_dry_run_is_read_only_and_expands_formula_isomers(tmp_path: Path) -> None:
    database = tmp_path / "library.db"
    catalog = tmp_path / "catalog.json"
    receipt = tmp_path / "receipt.json"
    _create_database(database)
    _write_catalog(catalog, _catalog())

    result = installer.install_catalog(catalog, database, receipt, apply=False)

    assert result["status"] == "preview"
    assert result["compound_count"] == 2
    assert result["target_row_count"] == 2
    assert result["inserted_target_rows"] == 0
    assert len(set(result["target_ids"])) == 2
    with sqlite3.connect(database) as connection:
        assert connection.execute("SELECT COUNT(*) FROM TargetCompounds").fetchone()[0] == 0
        assert connection.execute("SELECT marker FROM Spectra").fetchone()[0] == "unchanged"
    assert json.loads(receipt.read_text(encoding="utf-8"))["status"] == "preview"


def test_apply_is_target_only_and_second_apply_is_exact_noop(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    database = tmp_path / "library.db"
    catalog = tmp_path / "catalog.json"
    receipt = tmp_path / "receipt.json"
    _create_database(database)
    _write_catalog(catalog, _catalog())
    monkeypatch.setenv("MSSYSTEM_LIBRARY_WRITE_ENABLED", "1")

    first = installer.install_catalog(catalog, database, receipt, apply=True)
    second = installer.install_catalog(catalog, database, receipt, apply=True)

    assert first["status"] == "installed"
    assert first["inserted_target_rows"] == 2
    assert second["status"] == "already_installed_exact"
    assert second["inserted_target_rows"] == 0
    with sqlite3.connect(database) as connection:
        rows = connection.execute(
            "SELECT source_tier,reported_precursor_mz,reported_adduct,"
            "reported_retention_time,spectrum_accession,plant_source "
            "FROM TargetCompounds ORDER BY source_row_index"
        ).fetchall()
        assert len(rows) == 2
        assert {row[0] for row in rows} == {"literature_presence_only"}
        assert all(row[1:5] == (None, None, None, None) for row in rows)
        assert [row[5] for row in rows] == ["川乌", "草乌"]
        assert connection.execute("SELECT COUNT(*) FROM DiagnosticIons").fetchone()[0] == 0
        assert connection.execute("SELECT marker FROM Spectra").fetchone()[0] == "unchanged"


def test_existing_dataset_with_different_content_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    database = tmp_path / "library.db"
    catalog = tmp_path / "catalog.json"
    receipt = tmp_path / "receipt.json"
    _create_database(database)
    payload = _catalog()
    _write_catalog(catalog, payload)
    monkeypatch.setenv("MSSYSTEM_LIBRARY_WRITE_ENABLED", "1")
    installer.install_catalog(catalog, database, receipt, apply=True)

    payload["compounds"][0]["compound_name"] = "changed name"
    _write_catalog(catalog, payload)
    with pytest.raises(RuntimeError, match="different content"):
        installer.install_catalog(catalog, database, receipt, apply=True)
    with sqlite3.connect(database) as connection:
        assert connection.execute("SELECT COUNT(*) FROM TargetCompounds").fetchone()[0] == 2


def test_validation_rejects_mass_structure_and_measured_evidence(tmp_path: Path) -> None:
    catalog = tmp_path / "catalog.json"

    wrong_mass = _catalog()
    wrong_mass["compounds"][0]["exact_mass"] += 1.0
    _write_catalog(catalog, wrong_mass)
    with pytest.raises(installer.CatalogError, match="formula mass"):
        installer.load_catalog(catalog)

    wrong_structure = _catalog()
    wrong_structure["compounds"][0]["molecular_formula"] = "C8H10N4O2"
    wrong_structure["compounds"][0]["exact_mass"] = 194.08037556
    _write_catalog(catalog, wrong_structure)
    with pytest.raises(installer.CatalogError, match="formula"):
        installer.load_catalog(catalog)

    fake_measurement = _catalog()
    fake_measurement["compounds"][0]["retention_time"] = 5.2
    _write_catalog(catalog, fake_measurement)
    with pytest.raises(installer.CatalogError, match="measured-evidence"):
        installer.load_catalog(catalog)


def test_permanent_cation_uses_element_mass_and_electron_corrected_m_ion(
    tmp_path: Path,
) -> None:
    molecule = Chem.MolFromSmiles("C[N+](C)(C)C")
    assert molecule is not None
    element_formula = "C4H12N"
    element_mass = installer.formula_exact_mass(element_formula)
    assert element_mass is not None
    payload = _catalog()
    payload["compounds"] = [{
        "compound_name": "tetramethylammonium test identity",
        "molecular_formula": element_formula,
        "formal_charge": 1,
        "exact_mass": element_mass,
        "smiles": Chem.MolToSmiles(molecule, isomericSmiles=True),
        "inchi": Chem.MolToInchi(molecule),
        "inchikey": Chem.MolToInchiKey(molecule),
        "structure_identity_verified": True,
        "structure_identity_source": "reviewed structure",
        "source_row_index": 1,
        "plant_source": "延胡索",
        "provenance": "charged-identity fixture",
        "ions": [{
            "adduct": "[M]+",
            "charge": 1,
            "precursor_mz": element_mass - installer.ELECTRON_MASS_DA,
        }],
    }]
    catalog = tmp_path / "charged.json"
    _write_catalog(catalog, payload)

    _, rows = installer.load_catalog(catalog)

    assert rows[0]["exact_mass"] == pytest.approx(element_mass)
    assert rows[0]["precursor_mz"] == pytest.approx(
        element_mass - installer.ELECTRON_MASS_DA
    )
    assert "formal_charge:1" in json.loads(rows[0]["quality_flags_json"])
    assert rows[0]["reported_precursor_mz"] is None
    assert rows[0]["reported_adduct"] is None

    payload["compounds"][0]["ions"][0]["adduct"] = "[M+H]+"
    payload["compounds"][0]["ions"][0]["precursor_mz"] = ion_mz(
        element_mass, "[M+H]+"
    )
    _write_catalog(catalog, payload)
    with pytest.raises(installer.CatalogError, match="formal_charge"):
        installer.load_catalog(catalog)


def test_trigger_cannot_write_another_table_and_transaction_rolls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    database = tmp_path / "library.db"
    catalog = tmp_path / "catalog.json"
    receipt = tmp_path / "receipt.json"
    _create_database(database)
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE AuditLog (target_id TEXT)")
        connection.execute(
            "CREATE TRIGGER forbidden_side_effect AFTER INSERT ON TargetCompounds "
            "BEGIN INSERT INTO AuditLog(target_id) VALUES (NEW.target_id); END"
        )
    _write_catalog(catalog, _catalog())
    monkeypatch.setenv("MSSYSTEM_LIBRARY_WRITE_ENABLED", "1")

    with pytest.raises(sqlite3.DatabaseError, match="not authorized"):
        installer.install_catalog(catalog, database, receipt, apply=True)
    with sqlite3.connect(database) as connection:
        assert connection.execute("SELECT COUNT(*) FROM TargetCompounds").fetchone()[0] == 0
        assert connection.execute("SELECT COUNT(*) FROM AuditLog").fetchone()[0] == 0


def _semantic_payload(payload: dict) -> dict:
    payload["semantic_sha256"] = semantic_sha256(payload)
    return payload


def test_target_only_rebind_preserves_record_payload_and_updates_bindings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    database = tmp_path / "library.db"
    catalog = tmp_path / "catalog.json"
    install_receipt = tmp_path / "install.json"
    rebind_receipt = tmp_path / "rebind.json"
    _create_database(database)
    _write_catalog(catalog, _catalog())

    cache = tmp_path / "record-store"
    cache.mkdir()
    old_binding = source_binding(database)
    record_files = {
        "records.jsonl": {"sha256": "records", "size_bytes": 7},
        "record_offsets.npy": {"sha256": "offsets", "size_bytes": 8},
        "precursor_masses.npy": {"sha256": "pm", "size_bytes": 9},
        "precursor_record_indices.npy": {"sha256": "pi", "size_bytes": 10},
        "neutral_masses.npy": {"sha256": "nm", "size_bytes": 11},
        "neutral_record_indices.npy": {"sha256": "ni", "size_bytes": 12},
    }
    manifest = _semantic_payload({
        "schema_version": "mssystem.library-record-store/1.0.0",
        "source_db": old_binding,
        "processing_identity": processing_identity(),
        "record_count": 3_153_872,
        "normalized_records_sha256": "normalized-records",
        "record_alignment_sha256": "alignment",
        "mass_counts": {"precursor": 3, "neutral": 2},
        "files": record_files,
    })
    (cache / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    acceptance = _semantic_payload({
        "schema_version": "mssystem.library-record-store-acceptance/1.0.0",
        "manifest_semantic_sha256": manifest["semantic_sha256"],
        "full_record_roundtrip_sha256": manifest["normalized_records_sha256"],
        "verified_record_count": manifest["record_count"],
        "candidate_order_cases": 10_000,
        "candidate_order_mismatches": 0,
        "paired_scoring_queries": 100,
        "score_mismatches": 0,
        "shared_database_unchanged": True,
    })
    (cache / "acceptance.json").write_text(json.dumps(acceptance), encoding="utf-8")

    config_dir = tmp_path / "config"
    config_dir.mkdir()
    config_paths = {
        "record_store": config_dir / "library_record_store.json",
        "metadata_registry": config_dir / "library_metadata_repairs.json",
        "sql_metadata": config_dir / "sql_metadata_repairs.json",
        "mass_registry": config_dir / "library_compound_mass_repairs.json",
        "mass_runtime": config_dir / "source_compound_mass_repairs.json",
    }
    config_paths["record_store"].write_text(json.dumps({
        "cache_dir": str(cache),
        "accepted_cache_manifest_semantic_sha256": manifest["semantic_sha256"],
        "acceptance_receipt_file_sha256": rebind_script.file_sha256(
            cache / "acceptance.json"
        ),
    }), encoding="utf-8")
    for key in ("metadata_registry", "sql_metadata", "mass_registry", "mass_runtime"):
        config_paths[key].write_text(json.dumps(_semantic_payload({
            "source_database_sha256": old_binding["sha256"],
            **({
                "processing_identity": manifest["processing_identity"],
                "cache_manifest_semantic_sha256": manifest["semantic_sha256"],
            } if key == "sql_metadata" else {}),
        })), encoding="utf-8")

    monkeypatch.setenv("MSSYSTEM_LIBRARY_WRITE_ENABLED", "1")
    installer.install_catalog(catalog, database, install_receipt, apply=True)
    preview = rebind_script.rebind(
        database, catalog, install_receipt, rebind_receipt,
        apply=False, project_root=tmp_path, config_paths=config_paths,
    )
    applied = rebind_script.rebind(
        database, catalog, install_receipt, rebind_receipt,
        apply=True, project_root=tmp_path, config_paths=config_paths,
    )

    assert preview["status"] == "preview"
    assert applied["status"] == "applied"
    assert applied["binding_reused_from_preview"] is True
    assert applied["record_store_record_count"] == 3_153_872
    assert applied["normalized_records_sha256"] == "normalized-records"
    installed_manifest = json.loads((cache / "manifest.json").read_text())
    assert installed_manifest["record_count"] == manifest["record_count"]
    assert installed_manifest["normalized_records_sha256"] == manifest["normalized_records_sha256"]
    assert installed_manifest["files"] == record_files
    assert installed_manifest["source_db"]["sha256"] != old_binding["sha256"]
    assert json.loads(config_paths["metadata_registry"].read_text())[
        "source_database_sha256"
    ] == installed_manifest["source_db"]["sha256"]
