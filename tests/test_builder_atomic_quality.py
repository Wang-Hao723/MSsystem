"""Atomic rebuild and provenance-boundary tests for the formal SQLite library."""

from __future__ import annotations

import csv
import json
import os
import sqlite3
import sys
from collections import Counter
from pathlib import Path

import pytest

import build_from_mgf_msp as builder
from matching.library_quality import formula_exact_mass, source_tier_is_searchable


PROTON = 1.007276466621


def _record(
    spectrum_id: str,
    *,
    name: str = "Shared name",
    formula: str | None = "C2H6O",
    source_tier: str = "experimental",
    intensity_origin: str = "measured",
    peak_shift: float = 0.0,
) -> dict:
    neutral = formula_exact_mass(formula) if formula else None
    precursor = (neutral + PROTON) if neutral is not None else 150.0
    return {
        "spectrum_id": spectrum_id,
        "source_file": "synthetic_test.mgf",
        "name": name,
        "formula": formula,
        "precursor_mz": precursor,
        "adduct": "[M+H]+",
        "peaks": [
            [20.0 + peak_shift, 10.0],
            [30.0 + peak_shift, 20.0],
            [40.0 + peak_shift, 30.0],
        ],
        "source_tier": source_tier,
        "intensity_origin": intensity_origin,
        "classification": "other",
        "compound_class": "other",
    }


def _old_database(path: Path, value: str = "old-library") -> bytes:
    with sqlite3.connect(path) as conn:
        conn.execute("CREATE TABLE Sentinel (value TEXT NOT NULL)")
        conn.execute("INSERT INTO Sentinel VALUES (?)", (value,))
        conn.commit()
    return path.read_bytes()


def _assert_no_builder_temporaries(directory: Path) -> None:
    leftovers = [path.name for path in directory.iterdir() if path.name.endswith(".tmp")]
    assert leftovers == []


def test_atomic_build_applies_quality_gate_and_preserves_safe_compound_identity(
    tmp_path: Path,
):
    db_path = tmp_path / "formal.db"
    mapping_path = tmp_path / "mapping.csv"
    _old_database(db_path)
    mapping_path.write_text("old mapping\n", encoding="utf-8")

    records = [
        # Same no-InChIKey name + formula may share one Compound, but both
        # distinct experimental spectra must survive spectrum-level dedup.
        _record("exp_1"),
        _record("exp_2", peak_shift=0.1),
        # Same name with another formula must never inherit the first structure.
        _record("exp_other_formula", formula="C3H8O", peak_shift=0.2),
        # Name-only rows are deliberately spectrum-specific Compounds.
        # Even identical name-only peak lists cannot be deduplicated together.
        _record("name_only_1", formula=None, peak_shift=0.3),
        _record("name_only_2", formula=None, peak_shift=0.3),
        # Predicted provenance is retained, not promoted to experimental.
        _record("predicted_1", formula="C4H10O", source_tier="predicted", peak_shift=0.5),
        _record("predicted_2", formula="C5H12O", source_tier="in-silico", peak_shift=0.6),
        # Unknown explicit provenance is safely retained as unverified.
        _record("unverified", formula="C6H14O", source_tier="mystery_vendor", peak_shift=0.7),
        # A literature target with filled placeholder intensities is still not
        # admitted to the ordinary Spectra table.
        _record(
            "presence_only",
            formula="C7H16O",
            source_tier="literature_presence_only",
            intensity_origin="measured",
            peak_shift=0.8,
        ),
        # Internal formula/adduct/precursor inconsistency is quarantined.
        {
            **_record("bad_mass", formula="C8H18O", peak_shift=0.9),
            "precursor_mz": 500.0,
        },
    ]

    compounds, spectra = builder.build_database(records, db_path, mapping_path)

    assert (compounds, spectra) == (7, 8)
    with sqlite3.connect(db_path) as conn:
        assert conn.execute("PRAGMA quick_check").fetchone() == ("ok",)
        assert conn.execute("SELECT value FROM Sentinel").fetchone() == ("old-library",)
        assert conn.execute("SELECT COUNT(*) FROM Mapping").fetchone()[0] == spectra
        tier_counts = dict(conn.execute(
            "SELECT source_tier, COUNT(*) FROM Spectra GROUP BY source_tier"
        ))
        assert tier_counts == {
            "experimental": 5,
            "predicted": 2,
            "user_unverified": 1,
        }
        assert conn.execute(
            "SELECT COUNT(*) FROM Spectra WHERE spectrum_id IN ('presence_only','bad_mass')"
        ).fetchone()[0] == 0
        shared_ids = {
            row[0]
            for row in conn.execute("""
                SELECT m.compound_id FROM Mapping m
                JOIN Spectra s ON s.id = m.spectrum_id
                WHERE s.spectrum_id IN ('exp_1', 'exp_2')
            """)
        }
        assert len(shared_ids) == 1
        other_formula_id = conn.execute("""
            SELECT m.compound_id FROM Mapping m
            JOIN Spectra s ON s.id = m.spectrum_id
            WHERE s.spectrum_id = 'exp_other_formula'
        """).fetchone()[0]
        assert other_formula_id not in shared_ids
        name_only_ids = {
            row[0]
            for row in conn.execute("""
                SELECT m.compound_id FROM Mapping m
                JOIN Spectra s ON s.id = m.spectrum_id
                WHERE s.spectrum_id IN ('name_only_1', 'name_only_2')
            """)
        }
        assert len(name_only_ids) == 2
        quality_json = conn.execute(
            "SELECT quality_flags_json FROM Spectra WHERE spectrum_id = 'name_only_1'"
        ).fetchone()[0]
        assert "missing_formula_precursor_not_verified" in json.loads(quality_json)

    assert source_tier_is_searchable("experimental") is True
    assert source_tier_is_searchable("predicted") is False
    assert source_tier_is_searchable("user_unverified") is False

    with mapping_path.open(newline="", encoding="utf-8") as handle:
        assert sum(1 for _ in csv.DictReader(handle)) == spectra

    backups = list(tmp_path.glob("formal.db.prebuild_*.bak"))
    assert len(backups) == 1
    with sqlite3.connect(backups[0]) as conn:
        assert conn.execute("PRAGMA quick_check").fetchone() == ("ok",)
        assert conn.execute("SELECT value FROM Sentinel").fetchone() == ("old-library",)
    _assert_no_builder_temporaries(tmp_path)


def test_staging_audit_failure_keeps_old_database_and_mapping(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    db_path = tmp_path / "formal.db"
    mapping_path = tmp_path / "mapping.csv"
    old_bytes = _old_database(db_path, "audit-old")
    mapping_path.write_text("old mapping\n", encoding="utf-8")

    def fail_audit(*_args, **_kwargs):
        raise RuntimeError("injected staging audit failure")

    monkeypatch.setattr(builder, "_audit_staging_database", fail_audit)

    with pytest.raises(RuntimeError, match="injected staging audit failure"):
        builder.build_database([_record("valid")], db_path, mapping_path)

    assert db_path.read_bytes() == old_bytes
    with sqlite3.connect(db_path) as conn:
        assert conn.execute("SELECT value FROM Sentinel").fetchone() == ("audit-old",)
    assert mapping_path.read_text(encoding="utf-8") == "old mapping\n"
    assert list(tmp_path.glob("formal.db.prebuild_*.bak")) == []
    _assert_no_builder_temporaries(tmp_path)


def test_atomic_replace_failure_restores_mapping_and_keeps_old_database(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    db_path = tmp_path / "formal.db"
    mapping_path = tmp_path / "mapping.csv"
    old_bytes = _old_database(db_path, "replace-old")
    mapping_path.write_text("old mapping\n", encoding="utf-8")
    real_replace = os.replace

    def fail_only_formal_database_swap(source, destination):
        source_path = Path(source)
        destination_path = Path(destination)
        if destination_path == db_path and ".staging-" in source_path.name:
            raise OSError("injected atomic replace failure")
        return real_replace(source, destination)

    monkeypatch.setattr(builder.os, "replace", fail_only_formal_database_swap)

    with pytest.raises(OSError, match="injected atomic replace failure"):
        builder.build_database([_record("valid")], db_path, mapping_path)

    assert db_path.read_bytes() == old_bytes
    with sqlite3.connect(db_path) as conn:
        assert conn.execute("SELECT value FROM Sentinel").fetchone() == ("replace-old",)
    assert mapping_path.read_text(encoding="utf-8") == "old mapping\n"
    # The consistent old-library backup may remain, but no partial staging file
    # or rollback mapping is allowed to leak from the failed swap.
    backups = list(tmp_path.glob("formal.db.prebuild_*.bak"))
    assert len(backups) == 1
    with sqlite3.connect(backups[0]) as conn:
        assert conn.execute("SELECT value FROM Sentinel").fetchone() == ("replace-old",)
    _assert_no_builder_temporaries(tmp_path)


def test_mapping_audit_detects_broken_cardinality(tmp_path: Path):
    stage_db = tmp_path / "stage.db"
    stage_mapping = tmp_path / "stage.csv"
    quality_records, _summary = builder._quality_gate_records([_record("valid")])
    records = builder.deduplicate(quality_records)
    compounds, spectra, tiers = builder._write_staging_database(
        records, stage_db, stage_mapping
    )
    with sqlite3.connect(stage_db) as conn:
        conn.execute("DELETE FROM Mapping")
        conn.commit()

    with pytest.raises(RuntimeError, match="Mapping cardinality mismatch"):
        builder._audit_staging_database(
            stage_db,
            stage_mapping,
            expected_compounds=compounds,
            expected_spectra=spectra,
            expected_tiers=Counter(tiers),
        )


def test_atomic_rebuild_preserves_targets_core_extensions_and_indexes(
    tmp_path: Path,
):
    db_path = tmp_path / "formal.db"
    mapping_path = tmp_path / "mapping.csv"
    ethanol_ik = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    with sqlite3.connect(db_path) as conn:
        conn.executescript("""
            CREATE TABLE Compounds (
                id INTEGER PRIMARY KEY,
                compound_name TEXT,
                molecular_formula TEXT,
                smiles TEXT,
                inchi TEXT,
                inchikey TEXT,
                exact_mass REAL,
                classification TEXT,
                compound_class TEXT,
                classification_v3 TEXT,
                classification_source TEXT,
                classification_confidence REAL,
                fragmentation_pathway TEXT,
                custom_review_note TEXT
            );
            CREATE INDEX ix_compound_custom_review
                ON Compounds(custom_review_note);
            CREATE TABLE TargetCompounds (
                id INTEGER PRIMARY KEY,
                target_id TEXT UNIQUE NOT NULL,
                source_file TEXT,
                dataset_source TEXT,
                compound_name TEXT,
                molecular_formula TEXT,
                smiles TEXT,
                inchi TEXT,
                inchikey TEXT,
                exact_mass REAL,
                precursor_mz REAL,
                adduct TEXT,
                charge INTEGER,
                plant_source TEXT,
                source_tier TEXT NOT NULL DEFAULT 'literature_presence_only',
                quality_flags_json TEXT,
                created_at TEXT
            );
            CREATE TABLE DiagnosticIons (
                id INTEGER PRIMARY KEY,
                target_id INTEGER NOT NULL,
                fragment_mz REAL NOT NULL,
                ordinal INTEGER,
                FOREIGN KEY(target_id) REFERENCES TargetCompounds(id)
            );
            CREATE INDEX ix_target_precursor ON TargetCompounds(precursor_mz);
            CREATE INDEX ix_diagnostic_target ON DiagnosticIons(target_id);
            CREATE TABLE LibraryQualityMigrations (
                migration_id TEXT PRIMARY KEY,
                applied_at TEXT NOT NULL,
                backup_path TEXT NOT NULL,
                details_json TEXT NOT NULL
            );
            CREATE TABLE ReviewNotes (id INTEGER PRIMARY KEY, note TEXT);
            CREATE INDEX ix_review_note ON ReviewNotes(note);
        """)
        conn.execute(
            """INSERT INTO Compounds (
                compound_name, molecular_formula, smiles, inchikey,
                classification_v3, classification_source,
                classification_confidence, fragmentation_pathway,
                custom_review_note
            ) VALUES (?,?,?,?,?,?,?,?,?)""",
            (
                "Ethanol", "C2H6O", "CCO", ethanol_ik,
                "alcohol", "curated", 0.99, '{"path":"retained"}',
                "reviewed-by-human",
            ),
        )
        target_id = conn.execute(
            """INSERT INTO TargetCompounds (
                target_id, source_file, dataset_source, compound_name,
                molecular_formula, precursor_mz, adduct, source_tier,
                quality_flags_json
            ) VALUES (?,?,?,?,?,?,?,?,?)""",
            (
                "target-1", "literature.xlsx", "user_contrib:literature",
                "Reported target", "C10H12O3", 181.0859, "[M+H]+",
                "literature_presence_only", "[]",
            ),
        ).lastrowid
        conn.execute(
            "INSERT INTO DiagnosticIons (target_id, fragment_mz, ordinal) VALUES (?,?,?)",
            (target_id, 135.044, 0),
        )
        conn.execute(
            "INSERT INTO LibraryQualityMigrations VALUES (?,?,?,?)",
            ("v6-quality", "2026-07-16", "old.bak", "{}"),
        )
        conn.execute("INSERT INTO ReviewNotes (note) VALUES ('preserve me')")
        conn.commit()

    record = _record("ethanol", name="Ethanol")
    record.update({"smiles": "CCO", "inchikey": ethanol_ik})
    assert builder.build_database([record], db_path, mapping_path) == (1, 1)

    with sqlite3.connect(db_path) as conn:
        assert conn.execute("PRAGMA quick_check").fetchone() == ("ok",)
        assert conn.execute("PRAGMA foreign_key_check").fetchall() == []
        compound_columns = {
            row[1] for row in conn.execute("PRAGMA table_info(Compounds)")
        }
        assert {
            "classification_v3", "classification_source",
            "classification_confidence", "fragmentation_pathway",
            "custom_review_note",
        } <= compound_columns
        enriched = conn.execute(
            """SELECT classification_v3, classification_source,
                      classification_confidence, fragmentation_pathway,
                      custom_review_note
               FROM Compounds WHERE inchikey = ?""",
            (ethanol_ik,),
        ).fetchone()
        assert enriched == (
            "alcohol", "curated", 0.99, '{"path":"retained"}',
            "reviewed-by-human",
        )
        assert conn.execute("SELECT COUNT(*) FROM TargetCompounds").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM DiagnosticIons").fetchone()[0] == 1
        assert conn.execute(
            "SELECT migration_id FROM LibraryQualityMigrations"
        ).fetchone() == ("v6-quality",)
        assert conn.execute("SELECT note FROM ReviewNotes").fetchone() == ("preserve me",)
        indexes = {
            row[0]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='index'"
            )
        }
        assert {
            "ix_compound_custom_review", "ix_target_precursor",
            "ix_diagnostic_target", "ix_review_note",
        } <= indexes


def test_parser_metadata_is_fail_closed_and_preserves_name_charge_and_identity(
    tmp_path: Path,
):
    root = tmp_path / "sources"
    first_dir = root / "one"
    second_dir = root / "two"
    first_dir.mkdir(parents=True)
    second_dir.mkdir(parents=True)
    first = first_dir / "same.mgf"
    second = second_dir / "same.mgf"
    body = """BEGIN IONS
PEPMASS=150
CHARGE=1-
COMMENT=NAME: Multi word compound name
20 10
30 20
40 30
END IONS
"""
    first.write_text(body, encoding="utf-8")
    second.write_text(body, encoding="utf-8")

    first_record = builder.parse_mgf_file(first, source_root=root)[0][0]
    second_record = builder.parse_mgf_file(second, source_root=root)[0][0]

    assert first_record["name"] == "Multi word compound name"
    assert first_record["charge"] == -1
    assert first_record["source_tier"] == "user_unverified"
    assert first_record["spectrum_id"] != second_record["spectrum_id"]
    assert builder._coerce_int("2-") == -2
    assert builder._coerce_int("-3") == -3
    assert builder._coerce_int("2+") == 2


def test_massbank_supplement_mgf_preserves_identity_mode_ce_and_tier(
    tmp_path: Path,
):
    source = tmp_path / "massbank_supplement.mgf"
    source.write_text(
        """BEGIN IONS
TITLE=MassBank exact-key record
NAME=Ethanol
PEPMASS=47.0491
CHARGE=1+
SMILES=CCO
INCHIKEY=LFQSCWFLJHTTHZ-UHFFFAOYSA-N
ION_MODE=POSITIVE
INSTRUMENT_TYPE=LC-ESI-QTOF
INSTRUMENT=Q Exactive
SOURCE_ACCESSION=MSBNK-TEST-0001
SOURCE_TIER=massbank_supplement
COMMENT=INCHIKEY=LFQSCWFLJHTTHZ-UHFFFAOYSA-N COLLISIONENERGY=20-40eV CONFIDENCE=HIGH SOURCE=MassBank_2026.03
20 10
30 20
40 30
END IONS
BEGIN IONS
TITLE=MassBank key derived from structure
NAME=Ethanol
PEPMASS=47.0491
CHARGE=1+
SMILES=CCO
ION_MODE=POSITIVE
SOURCE_ACCESSION=MSBNK-TEST-0002
SOURCE_TIER=massbank_supplement
COMMENT=INCHIKEY= COLLISIONENERGY=20 CONFIDENCE=MED SOURCE=MassBank_2026.03
20 10
30 20
40 30
END IONS
""",
        encoding="utf-8",
    )

    records, _ = builder.parse_mgf_file(source)

    assert len(records) == 2
    first, second = records
    assert first["inchikey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert first["ion_mode"] == "POSITIVE"
    assert first["adduct"] is None
    assert first["source_tier"] == "massbank_supplement"
    assert first["collision_energy_raw"] == "20-40eV"
    assert first["collision_energy_unit"] == "eV"
    assert first["collision_energy_mode"] == "range"
    assert first["collision_energy_min"] == 20.0
    assert first["collision_energy_max"] == 40.0
    assert first["instrument_type"] == "LC-ESI-QTOF"
    assert first["instrument"] == "Q Exactive"
    assert second["inchikey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert second["inchikey"] != "COLLISIONENERGY=20"
    assert builder._canonical_source_tier("massbank_supplement") == (
        "massbank_supplement"
    )
    empty_ce = builder._parse_mgf_comment(
        "INCHIKEY= COLLISIONENERGY= CONFIDENCE=MED "
        "SOURCE=MassBank_2026.03"
    )
    assert not empty_ce.get("collision_energy_raw")

    database = tmp_path / "massbank.db"
    mapping = tmp_path / "massbank.csv"
    builder._write_staging_database([first], database, mapping)
    with sqlite3.connect(database) as connection:
        persisted = connection.execute(
            "SELECT collision_energy_raw, collision_energy_unit, "
            "collision_energy_mode, collision_energy_min, "
            "collision_energy_max, instrument_type, instrument, ion_mode, "
            "source_accession FROM Spectra"
        ).fetchone()
    assert persisted == (
        "20-40eV", "eV", "range", 20.0, 40.0, "LC-ESI-QTOF",
        "Q Exactive", "POSITIVE", "MSBNK-TEST-0001",
    )


def test_flat_placeholder_intensities_and_ms1_are_not_experimental_spectra(
    tmp_path: Path,
):
    flat = tmp_path / "flat.mgf"
    flat.write_text(
        """BEGIN IONS
TITLE=Literature ion list
PEPMASS=150
CHARGE=1+
SOURCE_TIER=experimental
20 1000
30 1000
40 1000
END IONS
""",
        encoding="utf-8",
    )
    flat_record = builder.parse_mgf_file(flat)[0][0]
    assert flat_record["intensity_origin"] == "default"
    accepted, summary = builder._quality_gate_records([flat_record])
    assert accepted == []
    assert summary["excluded_non_spectra"] == 1

    ms1 = tmp_path / "only_ms1.msp"
    ms1.write_text(
        """Name: EI entry
MW: 150
Spectrum_type: MS1
Num Peaks: 3
20 10
30 20
40 30
""",
        encoding="utf-8",
    )
    assert builder.parse_msp_file(ms1)[0] == []

    predicted = tmp_path / "predicted.msp"
    predicted.write_text(
        """Name: Predicted MS2 entry
PrecursorMZ: 101
Precursor_type: [M+H]+
Spectrum_type: MS2
Instrument_type: in-silico QTOF
Ion_mode: P
Num Peaks: 3
20 10
30 20
40 30
""",
        encoding="utf-8",
    )
    predicted_record = builder.parse_msp_file(predicted)[0][0]
    assert predicted_record["adduct"] == "[M+H]+"
    assert predicted_record["source_tier"] == "predicted"


def test_default_source_root_exists_and_empty_cli_fails_before_build(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    assert builder.DEFAULT_LIBRARY_ROOT.is_dir()
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setattr(
        sys,
        "argv",
        ["build_from_mgf_msp.py", "--mgf-dir", str(empty), "--output", str(tmp_path / "out")],
    )
    with pytest.raises(SystemExit) as error:
        builder.main()
    assert error.value.code == 2
    assert not (tmp_path / "out" / "local_library.db").exists()


def test_builder_quality_gate_accepts_only_verified_inchi_identity():
    ethanol_ik = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    record = _record("inchi_only", name="Ethanol")
    record.update({
        "smiles": None,
        "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "inchikey": ethanol_ik,
    })

    accepted, summary = builder._quality_gate_records([record])

    assert summary["accepted_spectra"] == 1
    assert len(accepted) == 1
    assert accepted[0]["inchikey"] == ethanol_ik
    assert accepted[0]["smiles"] == "CCO"
    assert accepted[0]["structure_identity_verified"] is True
    assert accepted[0]["structure_identity_source"] == "inchi_derived"
