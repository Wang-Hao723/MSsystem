"""Regression tests for library provenance, literature targets, and contaminants.

All database tests use temporary SQLite files.  They intentionally avoid the
shared production library so the quality boundary can be tested without any
external state or multi-gigabyte library load.
"""

from __future__ import annotations

import json
import hashlib
import importlib
import sqlite3
from pathlib import Path

import pytest

from build_from_mgf_msp import deduplicate, parse_mgf_file
from matching.compound_annotations import (
    annotate_compound_dict,
    annotation_for,
    compound_annotation_catalog_summary,
    load_compound_annotations,
)
from matching.compound_names import display_compound_name
from matching.honest_calibration import annotate as honest_annotate
from matching.library_import import compound_match_key
from matching.library_quality import (
    PRESENCE_ONLY_SOURCE_TIER,
    QUARANTINED_SOURCE_TIER,
    assess_record,
    assess_import_records,
    assess_records,
    formula_exact_mass,
    neutral_formula_text,
    resolve_structure_formula,
    source_tier_is_searchable,
    validate_structure_metadata,
)
from matching.literature_targets import screen_results


OXIDIZED_IRGAFOS_IK = "AZSKHRTUXHLAHS-UHFFFAOYSA-N"
IRGAFOS_IK = "JKIJEFPNVSHHEI-UHFFFAOYSA-N"


@pytest.fixture(autouse=True)
def _enable_writes_only_for_isolated_test_libraries(monkeypatch):
    """Keep V9 production read-only while allowing this file's temp SQLite tests."""
    from config import settings

    monkeypatch.setattr(settings, "LIBRARY_WRITE_ENABLED", True)


def _glucose_precursor() -> float:
    mass = formula_exact_mass("C6H12O6")
    assert mass is not None
    return mass + 1.007276


def _presence_only_record(name: str = "Reported compound") -> dict:
    return {
        "name": name,
        "source_file": "literature.xlsx",
        "dataset_source": "user_contrib:literature",
        "formula": "C6H12O6",
        "precursor_mz": _glucose_precursor(),
        "adduct": "[M+H]+",
        "charge": 1,
        "peaks": [[89.0, 1000.0], [119.0, 1000.0], [149.0, 1000.0]],
        "intensity_origin": "default",
        "plant_source": "Synthetic test herb",
    }


def _measured_ethanol_record(*, smiles: str = "CCO", name: str = "Ethanol") -> dict:
    return {
        "name": name,
        "source_file": "measured.mgf",
        "dataset_source": "user_contrib:measured",
        "spectrum_id": "scan_1",
        "formula": "C2H6O",
        "smiles": smiles,
        "precursor_mz": 47.049141,
        "adduct": "[M+H]+",
        "charge": 1,
        "peaks": [[19.0, 100.0], [29.0, 75.0], [45.0, 50.0]],
        "intensity_origin": "measured",
        "source_tier": "experimental",
        "instrument_type": "QTOF",
    }


def test_nfkc_numeric_parser_removes_only_whitespace_between_digits():
    from csv_excel_to_mgf import safe_float

    assert safe_float("1 155．276 48") == pytest.approx(1155.27648)
    assert safe_float("1\u00a0155.25") == pytest.approx(1155.25)
    assert safe_float("12 mg") is None


def test_mgf_exactmass_parser_keeps_the_leading_digit(tmp_path: Path):
    source = tmp_path / "exactmass.mgf"
    source.write_text(
        "BEGIN IONS\n"
        "TITLE=Exact-mass regression\n"
        "PEPMASS=407.2072\n"
        "CHARGE=1-\n"
        "ADDUCT=[M-H]-\n"
        "EXACTMASS=408.21480336835\n"
        "59.0127 100\n"
        "135.0802 10\n"
        "149.0962 5\n"
        "END IONS\n",
        encoding="utf-8",
    )

    records, _ = parse_mgf_file(source)

    assert len(records) == 1
    assert records[0]["exact_mass"] == pytest.approx(408.21480336835)


def test_standard_parser_keeps_formula_only_expected_target(tmp_path: Path):
    from csv_excel_to_mgf import parse_standard

    source = tmp_path / "formula_only.csv"
    source.write_text(
        "spectrum_name,formula,precursor_mz,adduct\n"
        "Formula-only target,C6H12O6,,[M+H]+\n",
        encoding="utf-8",
    )
    rows = parse_standard(source)

    assert len(rows) == 1
    assert rows[0]["name"] == "Formula-only target"
    assert rows[0]["formula"] == "C6H12O6"
    assert rows[0]["precursor_mz"] is None
    assert rows[0]["fragments"] == []


def test_expected_target_gate_accepts_formula_row_without_precursor_or_ms2():
    accepted, quality = assess_import_records([
        {
            "name": "Formula-only literature target",
            "source_file": "reported.xlsx",
            "dataset_source": "user_contrib:reported",
            "formula": "C6H12O6",
            "precursor_mz": None,
            "adduct": None,
            "peaks": [],
            "intensity_origin": "default",
        }
    ])

    assert len(accepted) == 1
    target = accepted[0]
    assert target["source_tier"] == PRESENCE_ONLY_SOURCE_TIER
    assert target["target_searchability"] == "theoretical_ion_ready"
    assert "reported_precursor_missing_theoretical_ion_required" in (
        target["quality_flags"]
    )
    assert quality["record_types"] == {
        "expected_target_input": 1,
        "experimental_spectrum_input": 0,
        "expected_targets_accepted": 1,
        "experimental_spectra_accepted": 0,
    }


def test_ao1680_exact_alias_and_external_contaminant_annotation():
    canonical = (
        "Tris(2,4-di-tert-butylphenyl) phosphate (oxidized Irgafos 168)"
    )
    assert display_compound_name("AO1680", OXIDIZED_IRGAFOS_IK) == canonical
    assert display_compound_name("Irgafos 168", IRGAFOS_IK).endswith(
        "phosphite (Irgafos 168)"
    )

    match = {
        "compound_name": "AO1680",
        "inchikey": OXIDIZED_IRGAFOS_IK,
        "classification": "phenolic acid",
        "display_class": "phenolic acid",
    }
    annotate_compound_dict(match)

    assert match["compound_name"] == canonical
    assert match["source_compound_name"] == "AO1680"
    assert match["compound_role"] == "external_contaminant"
    assert match["potential_contaminant"] is True
    assert match["exclude_from_natural_product_count"] is True
    assert match["natural_product_candidate"] is False
    assert match["classification"] == "other"
    assert match["classification_override_source"] == "curated_compound_role"

    # Role curation must never spread by a short connectivity block or name.
    assert annotation_for(OXIDIZED_IRGAFOS_IK[:14]) == {}
    ordinary = {"compound_name": "AO1680", "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C"}
    annotate_compound_dict(ordinary)
    assert ordinary["potential_contaminant"] is False
    assert ordinary["compound_role"] == "unclassified"
    assert ordinary["natural_product_candidate"] is False
    assert ordinary["exclude_from_natural_product_count"] is True


def test_dibutyl_phthalate_is_curated_as_external_plasticizer():
    item = {
        "compound_name": "Spectral Match to Dibutyl phthalate from NIST14",
        "inchikey": "DOIRQSBPFJWKBE-UHFFFAOYSA-N",
        "classification": "phenolic_acid",
    }
    annotate_compound_dict(item)
    assert item["compound_name"] == "Dibutyl phthalate"
    assert item["source_compound_name"].startswith("Spectral Match to")
    assert item["potential_contaminant"] is True
    assert item["contaminant_category"] == "phthalate_plasticizer"
    assert item["classification"] == "other"


def test_v12_alkylbenzene_sulfonates_are_separated_as_external_surfactants():
    exact_keys = {
        "UAZLASMTBCLJKO-UHFFFAOYSA-N",  # decyl
        "FERBTPHUEYEGDN-UHFFFAOYSA-N",  # branched undecyl
        "UCDCOJNNUVYFKJ-UHFFFAOYSA-N",  # undecyl
        "KWXICGTUELOLSQ-UHFFFAOYSA-N",  # dodecyl
    }
    for inchikey in exact_keys:
        annotation = annotation_for(inchikey)
        assert annotation["compound_role"] == "external_contaminant"
        assert annotation["contaminant_category"] == (
            "alkylbenzene_sulfonate_surfactant"
        )
        assert annotation["exclude_from_natural_product_count"] is True

    summary = compound_annotation_catalog_summary()
    assert summary["external_contaminant_count"] == 95
    # Four alkylbenzene sulfonates plus the exact-identity amine-oxide
    # surfactant added for the MS2 586678 regression case.
    assert summary["supplemental_annotation_count"] == 5
    assert summary["supplemental_external_contaminant_count"] == 5


def test_ncgc_formula_identifier_is_replaced_by_readable_exact_structure_name():
    inchikey = "MPYSZNHNOUOLSB-UHFFFAOYSA-N"
    assert display_compound_name(
        "NCGC00385425-01_C21H36O11_", inchikey
    ) == "5-Hydroxybornyl 6-O-apiosyl-beta-D-glucoside"


def test_v9_external_contaminant_baseline_is_preserved_exactly():
    annotations = load_compound_annotations()
    keys = sorted(
        key
        for key, value in annotations.items()
        if value.get("compound_role") == "external_contaminant"
    )
    digest = hashlib.sha256(("\n".join(keys) + "\n").encode("utf-8")).hexdigest()
    summary = compound_annotation_catalog_summary()

    assert len(keys) == 95
    assert digest == "5877216afacc37d4df610777731c101128801fda7e8354defa736c7310148258"
    assert summary["external_contaminant_count"] == 95
    assert summary["class_level_structure_rule_count"] == 66
    assert summary["baseline_contract"]["external_contaminant_count"] == 95

    dehp = annotation_for("BJQHLKABXJIVAM-UHFFFAOYSA-N")
    assert dehp["compound_role"] == "external_contaminant"
    assert dehp["contaminant_category"] == "phthalate_plasticizer"
    assert dehp["role_annotation_scope"] == "exact_identity_curation"


def test_current_catalog_overwrites_stale_cached_role_projection():
    item = {
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "compound_role": "external_contaminant",
        "potential_contaminant": True,
        "contaminant_category": "stale_category",
        "screening_catalog_id": "retired-catalog",
    }

    annotate_compound_dict(item)

    assert item["source_compound_role"] == "external_contaminant"
    assert item["compound_role"] == "unclassified"
    assert item["potential_contaminant"] is False
    assert item["contaminant_category"] == ""
    assert "screening_catalog_id" not in item


def test_diisobutyl_phthalate_is_named_and_curated_as_external_plasticizer():
    item = {
        "compound_name": "Unresolved compound (MGWAVDBGNNKXQV)",
        "inchikey": "MGWAVDBGNNKXQV-UHFFFAOYSA-N",
        "classification": "phenolic_acid",
        "display_class": "phenolic_acid",
    }
    annotate_compound_dict(item)

    assert item["compound_name"] == "Diisobutyl phthalate"
    assert item["compound_name_cn"] == "邻苯二甲酸二异丁酯（DIBP）"
    assert item["compound_role"] == "external_contaminant"
    assert item["potential_contaminant"] is True
    assert item["contaminant_category"] == "phthalate_plasticizer"
    assert item["classification"] == "other"
    assert item["display_class"] == "other"
    assert item["exclude_from_natural_product_count"] is True


@pytest.mark.parametrize(
    ("source_name", "inchikey", "canonical_start"),
    [
        (
            "AO1081",
            "MQWCQFCZUNBTCM-UHFFFAOYSA-N",
            "2,2'-Thiobis",
        ),
        (
            "AO22E46",
            "DXCHWXWXYPEZKM-UHFFFAOYSA-N",
            "2,2'-Ethylidenebis",
        ),
        (
            "AO4426",
            "MDWVSAYEQPLWMX-UHFFFAOYSA-N",
            "4,4'-Methylenebis",
        ),
        (
            "AO4703",
            "VMZVBRIIHDRYGK-UHFFFAOYSA-N",
            "2,6-Di-tert-butyl",
        ),
    ],
)
def test_ao_library_codes_are_named_and_curated_as_synthetic_antioxidants(
    source_name, inchikey, canonical_start
):
    item = {
        "compound_name": source_name,
        "inchikey": inchikey,
        "classification": "lignan",
        "display_class": "lignan",
    }
    annotate_compound_dict(item)

    assert item["compound_name"].startswith(canonical_start)
    assert item["source_compound_name"] == source_name
    assert item["compound_name_cn"]
    assert item["compound_role"] == "external_contaminant"
    assert item["potential_contaminant"] is True
    assert item["contaminant_category"] == "synthetic_phenolic_antioxidant"
    assert item["classification"] == "other"
    assert item["display_class"] == "other"
    assert item["exclude_from_natural_product_count"] is True
    assert item["natural_product_candidate"] is False


def test_missing_intensities_import_only_as_targets(monkeypatch, tmp_path: Path):
    """A diagnostic-ion list must not create a pseudo experimental spectrum."""
    from config import settings
    from matching import library_import

    db_path = tmp_path / "targets_only.db"
    incoming = _presence_only_record()
    incoming.update({
        "retention_time": 4.25,
        "retention_time_unit": "minutes",
        "collision_energy": 20.0,
        "instrument": "LC-ESI-QTOF",
    })
    monkeypatch.setattr(settings, "DB_PATH", db_path)
    monkeypatch.setattr(settings, "DB_PATH_SMALL", tmp_path / "unused_small.db")
    monkeypatch.setattr(
        library_import,
        "_files_to_records",
        lambda _paths, **_kwargs: (
            [incoming], [{"file": "literature.xlsx", "written": 1}]
        ),
    )

    result = library_import.import_files([tmp_path / "literature.xlsx"], db_mode="full")

    assert result["spectral_records"] == 0
    assert result["literature_target_records"] == 1
    assert result["literature_targets"] == {
        "new_targets": 1,
        "new_diagnostic_ions": 3,
        "removed_superseded_targets": 0,
        "replaced_datasets": ["user_contrib:literature"],
        "database": "full",
        "ordinary_spectral_search": False,
    }
    with sqlite3.connect(db_path) as conn:
        tables = {
            row[0]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        assert {"TargetCompounds", "DiagnosticIons"} <= tables
        assert "Spectra" not in tables
        target = conn.execute(
            "SELECT compound_name, source_tier, reported_precursor_mz, "
            "reported_adduct, reported_retention_time, "
            "reported_retention_time_unit, reported_collision_energy_raw, "
            "instrument, target_searchability FROM TargetCompounds"
        ).fetchone()
        assert target == (
            "Reported compound",
            PRESENCE_ONLY_SOURCE_TIER,
            pytest.approx(_glucose_precursor()),
            "[M+H]+",
            pytest.approx(4.25),
            "minutes",
            "20.0",
            "LC-ESI-QTOF",
            "theoretical_ion_ready",
        )
        assert conn.execute("SELECT COUNT(*) FROM DiagnosticIons").fetchone()[0] == 3


def test_incremental_schema_reuses_equivalent_legacy_indexes(tmp_path: Path):
    """Index names may differ across builders; an import must not duplicate them."""
    from matching import library_import

    db_path = tmp_path / "legacy_indexes.db"
    with sqlite3.connect(db_path) as conn:
        conn.executescript(
            """
            CREATE TABLE Compounds (id INTEGER PRIMARY KEY, compound_name TEXT,
                molecular_formula TEXT, inchikey TEXT);
            CREATE TABLE Spectra (id INTEGER PRIMARY KEY, spectrum_id TEXT UNIQUE,
                precursor_mz REAL);
            CREATE TABLE Mapping (id INTEGER PRIMARY KEY, compound_id INTEGER,
                spectrum_id INTEGER, score REAL);
            CREATE INDEX ix_Comp_inchikey ON Compounds(inchikey);
            CREATE INDEX ix_Comp_name ON Compounds(compound_name);
            CREATE INDEX ix_Spec_sid ON Spectra(spectrum_id);
            CREATE INDEX ix_Spec_pmz ON Spectra(precursor_mz);
            CREATE INDEX ix_Map_cid ON Mapping(compound_id);
            """
        )
        library_import._create_schema(conn)
        indexes = {
            row[0]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='index'"
            )
        }

    assert not {"ix_c_ik", "ix_c_name", "ix_s_sid", "ix_s_pmz", "ix_m_c"} & indexes


def test_formula_adduct_precursor_conflict_is_rejected():
    record = _presence_only_record("Impossible glucose row")
    record.update({
        "precursor_mz": 500.0,
        "intensity_origin": "reported_relative",
    })

    clean, audit = assess_record(record, max_precursor_ppm=50.0)

    assert audit["accepted"] is False
    assert audit["source_tier"] == QUARANTINED_SOURCE_TIER
    assert "formula_adduct_precursor_mismatch" in audit["errors"]
    assert abs(audit["precursor_ppm"]) > 50
    assert clean["source_tier"] == QUARANTINED_SOURCE_TIER


def test_formula_precursor_can_repair_a_clear_proton_polarity_typo():
    record = _presence_only_record("Glucose negative ion")
    mass = formula_exact_mass("C6H12O6")
    assert mass is not None
    record.update({
        "precursor_mz": mass - 1.007276,
        "adduct": "[M+H]+",  # deliberately wrong sign
    })

    clean, audit = assess_record(record, max_precursor_ppm=50.0)

    assert audit["accepted"] is True
    assert audit["original_adduct"] == "[M+H]+"
    assert clean["adduct"] == "[M-H]-"
    assert "adduct_repaired_from_formula_precursor" in audit["warnings"]


def test_explicit_deprotonated_formula_is_normalized_to_neutral_parent():
    assert neutral_formula_text("C21H19O12-", "[M-H]-") == "C21H20O12"
    assert neutral_formula_text("C21H20O12", "[M-H]-") == "C21H20O12"
    assert neutral_formula_text("C21H19O12−", "[M-H]-") == "C21H20O12"
    assert neutral_formula_text("C21H19O12-", "[M+H]+") == ""
    assert neutral_formula_text("C21H19O12-", "") == ""


def test_charged_smiles_and_formula_validate_against_neutral_parent():
    charged_smiles = (
        "O=c1c(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)"
        "c(-c2ccc(O)c([O-])c2)oc2cc(O)cc(O)c12"
    )
    resolved = resolve_structure_formula(
        charged_smiles, "C21H19O12-", "[M-H]-"
    )
    valid, status = validate_structure_metadata(
        charged_smiles, "C21H19O12-", "[M-H]-"
    )

    assert resolved == "C21H20O12"
    assert valid is True
    assert status == "validated"


def test_charged_formula_without_supported_adduct_fails_closed():
    charged_smiles = (
        "O=c1c(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)"
        "c(-c2ccc(O)c([O-])c2)oc2cc(O)cc(O)c12"
    )

    assert resolve_structure_formula(charged_smiles, "C21H19O12-", "") == ""
    valid, status = validate_structure_metadata(
        charged_smiles, "C21H19O12-", ""
    )

    assert valid is False
    assert status == "charged_formula_neutralization_unresolved"


def test_charged_library_record_keeps_structure_and_neutral_exact_mass():
    charged_smiles = (
        "O=c1c(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)"
        "c(-c2ccc(O)c([O-])c2)oc2cc(O)cc(O)c12"
    )
    neutral_mass_value = formula_exact_mass("C21H20O12")
    assert neutral_mass_value is not None
    record = {
        "name": "Deprotonated flavonoid glycoside",
        "source_file": "charged.mgf",
        "dataset_source": "user_contrib:charged",
        "spectrum_id": "charged_1",
        "formula": "C21H19O12-",
        "smiles": charged_smiles,
        "precursor_mz": neutral_mass_value - 1.007276,
        "adduct": "[M-H]-",
        "charge": -1,
        "peaks": [[151.0, 100.0], [271.0, 80.0], [300.0, 60.0]],
        "intensity_origin": "measured",
        "source_tier": "experimental",
        "instrument_type": "QTOF",
    }

    clean, audit = assess_record(record)

    assert audit["accepted"] is True
    assert clean["formula"] == "C21H20O12"
    assert clean["smiles"]
    assert clean["exact_mass"] == pytest.approx(neutral_mass_value, abs=1e-8)
    assert "smiles_formula_mismatch_structure_removed" not in audit["warnings"]


def test_shared_precursor_and_fragments_do_not_merge_different_compounds():
    common = {
        "source_tier": PRESENCE_ONLY_SOURCE_TIER,
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",  # copied/incorrect template structure
        "formula": "C20H20O10",
        "precursor_mz": 421.1234,
        "adduct": "[M+H]+",
        "peaks": [[121.0, 1000.0], [151.0, 1000.0], [301.0, 1000.0]],
    }
    first = {**common, "name": "Astragaloside alpha", "spectrum_id": "row_1"}
    second = {**common, "name": "Astragaloside beta", "spectrum_id": "row_2"}

    retained = deduplicate([first, second])

    assert len(retained) == 2
    assert {row["name"] for row in retained} == {
        "Astragaloside alpha",
        "Astragaloside beta",
    }
    assert compound_match_key(first) != compound_match_key(second)


def test_verified_full_inchikey_exact_miss_never_falls_back_to_name_formula(
    tmp_path: Path,
):
    """Same label/formula cannot merge two different structure-derived keys."""
    from matching import library_import
    from matching.library_quality import smiles_metadata

    db_path = tmp_path / "identity_guard.db"
    ethanol_key = smiles_metadata("CCO")["inchikey"]
    ether_key = smiles_metadata("COC")["inchikey"]
    assert ethanol_key != ether_key
    with sqlite3.connect(db_path) as conn:
        library_import._create_schema(conn)
        conn.execute(
            "INSERT INTO Compounds (compound_name, molecular_formula, smiles, inchikey, "
            "exact_mass, classification, structure_identity_verified, "
            "structure_identity_source) VALUES (?,?,?,?,?,?,?,?)",
            (
                "Same submitted label", "C2H6O", "CCO", ethanol_key,
                46.041865, "other", 1, "smiles_derived",
            ),
        )
        conn.commit()

    incoming, audit = assess_record(
        _measured_ethanol_record(smiles="COC", name="Same submitted label")
    )
    assert audit["accepted"] is True
    assert incoming["inchikey"] == ether_key
    assert incoming["structure_identity_verified"] is True

    new_compounds, new_spectra = library_import.append_records(
        [incoming], db_path, ts_prefix="identity_exact_miss"
    )

    assert (new_compounds, new_spectra) == (1, 1)
    with sqlite3.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT inchikey, structure_identity_verified FROM Compounds ORDER BY id"
        ).fetchall()
        mapped_key = conn.execute(
            "SELECT c.inchikey FROM Mapping m "
            "JOIN Compounds c ON c.id=m.compound_id "
            "JOIN Spectra s ON s.id=m.spectrum_id "
            "WHERE s.spectrum_id LIKE 'identity_exact_miss%'"
        ).fetchone()[0]
    assert rows == [(ethanol_key, 1), (ether_key, 1)]
    assert mapped_key == ether_key


def test_unverified_user_inchikey_is_provenance_not_structure_identity(
    tmp_path: Path,
):
    from matching import library_import
    from matching.library_quality import smiles_metadata

    db_path = tmp_path / "unverified_key.db"
    submitted_key = smiles_metadata("CCO")["inchikey"]
    different_key = smiles_metadata("COC")["inchikey"]
    with sqlite3.connect(db_path) as conn:
        library_import._create_schema(conn)
        conn.execute(
            "INSERT INTO Compounds (compound_name, molecular_formula, smiles, inchikey, "
            "exact_mass, classification, structure_identity_verified) "
            "VALUES (?,?,?,?,?,?,?)",
            (
                "Same submitted label", "C2H6O", "COC", different_key,
                46.041865, "other", 1,
            ),
        )
        conn.commit()

    raw = _measured_ethanol_record(smiles="", name="Same submitted label")
    raw["inchikey"] = submitted_key
    incoming, audit = assess_record(raw)

    assert audit["accepted"] is True
    assert incoming["reported_inchikey"] == submitted_key
    assert incoming["inchikey"] is None
    assert incoming["structure_identity_verified"] is False
    assert "unverified_inchikey_removed_from_structure_identity" in audit["warnings"]

    new_compounds, _new_spectra = library_import.append_records(
        [incoming], db_path, ts_prefix="identity_unverified"
    )
    assert new_compounds == 1
    with sqlite3.connect(db_path) as conn:
        inserted = conn.execute(
            "SELECT inchikey, reported_inchikey, structure_identity_verified, "
            "structure_identity_source FROM Compounds ORDER BY id DESC LIMIT 1"
        ).fetchone()
    assert inserted == (None, submitted_key, 0, "reported_unverified")


def test_import_is_atomic_across_full_small_targets_and_retry_idempotent(
    monkeypatch, tmp_path: Path
):
    from config import settings
    from matching import library_import

    full_db = tmp_path / "full.db"
    small_db = tmp_path / "small.db"
    measured = _measured_ethanol_record()
    target = _presence_only_record("Literature target")
    monkeypatch.setattr(settings, "DB_PATH", full_db)
    monkeypatch.setattr(settings, "DB_PATH_SMALL", small_db)
    monkeypatch.setattr(library_import, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(
        library_import,
        "_files_to_records",
        lambda _paths, **_kwargs: (
            [measured, target], [{"file": "mixed.xlsx", "written": 2}]
        ),
    )

    real_mapping_append = library_import._append_mapping_csv_atomic
    calls = 0

    def fail_second_mapping(path, rows):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("injected small mapping failure")
        return real_mapping_append(path, rows)

    monkeypatch.setattr(
        library_import, "_append_mapping_csv_atomic", fail_second_mapping
    )
    with pytest.raises(RuntimeError, match="injected small mapping failure"):
        library_import.import_files([tmp_path / "mixed.xlsx"], db_mode="both")

    # Both databases and the already-appended Full mapping were compensated.
    assert not full_db.exists()
    assert not small_db.exists()
    assert not (tmp_path / "data" / "library_mapping.csv").exists()
    assert not (tmp_path / "data" / "library_mapping_small.csv").exists()

    monkeypatch.setattr(
        library_import, "_append_mapping_csv_atomic", real_mapping_append
    )
    first = library_import.import_files([tmp_path / "mixed.xlsx"], db_mode="both")
    retry = library_import.import_files([tmp_path / "mixed.xlsx"], db_mode="both")

    assert first["full"] == {"new_compounds": 1, "new_spectra": 1}
    assert first["small"] == {"new_compounds": 1, "new_spectra": 1}
    assert first["literature_targets"]["new_targets"] == 1
    assert retry["import_batch_id"] == first["import_batch_id"]
    assert retry["full"] == {"new_compounds": 0, "new_spectra": 0}
    assert retry["small"] == {"new_compounds": 0, "new_spectra": 0}
    assert retry["literature_targets"]["new_targets"] == 0

    with sqlite3.connect(full_db) as conn:
        assert conn.execute("SELECT COUNT(*) FROM Spectra").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM TargetCompounds").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM DiagnosticIons").fetchone()[0] == 3
    with sqlite3.connect(small_db) as conn:
        assert conn.execute("SELECT COUNT(*) FROM Spectra").fetchone()[0] == 1
    for mapping_name in ("library_mapping.csv", "library_mapping_small.csv"):
        mapping_lines = (
            tmp_path / "data" / mapping_name
        ).read_text(encoding="utf-8").splitlines()
        assert len(mapping_lines) == 2


def test_small_import_rejects_unknown_neutral_mass_and_targets_deduplicate_ions(
    tmp_path: Path,
):
    from matching import library_import

    small_db = tmp_path / "small_unknown_mass.db"
    unknown_mass = _measured_ethanol_record(smiles="", name="Unknown huge precursor")
    unknown_mass.update({"formula": None, "exact_mass": None, "precursor_mz": 1500.0})
    nc, ns = library_import.append_records(
        [unknown_mass],
        small_db,
        ts_prefix="small_guard",
        mass_threshold=800.0,
    )
    assert (nc, ns) == (0, 0)
    with sqlite3.connect(small_db) as conn:
        assert conn.execute("SELECT COUNT(*) FROM Spectra").fetchone()[0] == 0

    target_db = tmp_path / "target_ions.db"
    target = _presence_only_record()
    target["peaks"] = [
        [149.0, 1.0], [89.0, 1.0], [149.0, 999.0], [119.0, 1.0]
    ]
    nt, ni = library_import.append_target_records([target], target_db)
    assert (nt, ni) == (1, 3)
    with sqlite3.connect(target_db) as conn:
        ions = [
            row[0] for row in conn.execute(
                "SELECT fragment_mz FROM DiagnosticIons ORDER BY ordinal"
            )
        ]
    assert ions == [89.0, 119.0, 149.0]


def test_corrected_target_dataset_replaces_superseded_snapshot(tmp_path: Path):
    from matching import library_import

    db_path = tmp_path / "targets.db"
    old = _presence_only_record("Old incorrect row")
    old["dataset_source"] = "user_contrib:formula"
    current = _presence_only_record("Correct current row")
    current["dataset_source"] = "user_contrib:formula"

    library_import._append_target_records_with_receipt(
        [old], db_path, replace_datasets=True,
    )
    replaced = library_import._append_target_records_with_receipt(
        [current], db_path, replace_datasets=True,
    )
    retry = library_import._append_target_records_with_receipt(
        [current], db_path, replace_datasets=True,
    )

    assert replaced["removed_targets"] == 1
    assert replaced["new_targets"] == 1
    assert replaced["replaced_datasets"] == ["user_contrib:formula"]
    assert retry["new_targets"] == 0
    assert retry["removed_targets"] == 0
    with sqlite3.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT compound_name FROM TargetCompounds"
        ).fetchall()
        assert rows == [("Correct current row",)]
        assert conn.execute("SELECT COUNT(*) FROM DiagnosticIons").fetchone()[0] == 3


def test_quality_migration_backup_uses_consistent_sqlite_snapshot(tmp_path: Path):
    import importlib.util

    script_path = Path(__file__).resolve().parents[1] / "scripts" / "migrate_v6_library_quality.py"
    spec = importlib.util.spec_from_file_location("migrate_v6_library_quality", script_path)
    assert spec and spec.loader
    migration = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(migration)

    source = tmp_path / "source.db"
    with sqlite3.connect(source) as conn:
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("CREATE TABLE Parent(id INTEGER PRIMARY KEY)")
        conn.execute(
            "CREATE TABLE Child(id INTEGER PRIMARY KEY, parent_id INTEGER "
            "REFERENCES Parent(id))"
        )
        conn.execute("INSERT INTO Parent VALUES (1)")
        conn.execute("INSERT INTO Child VALUES (1, 1)")
    backup = migration._clone_backup(source)

    assert backup.exists()
    assert migration._table_row_counts(backup) == {"Child": 1, "Parent": 1}
    assert migration._foreign_key_errors(backup) == []


def test_legacy_mona_migration_recovers_predicted_and_ms1_by_old_ids(tmp_path: Path):
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    script_path = root / "scripts" / "migrate_legacy_mona_provenance.py"
    spec = importlib.util.spec_from_file_location("migrate_legacy_mona_provenance", script_path)
    assert spec and spec.loader
    migration = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(migration)

    msp = tmp_path / "MoNA-export-All_Spectra.msp"
    msp.write_text(
        "Name: old ms1\nSpectrum_type: MS1\nMW: 100\nNum Peaks: 2\n10 100\n20 50\n\n"
        "Name: predicted\nSpectrum_type: MS2\nPrecursorMZ: 100\n"
        "Instrument_type: in-silico QTOF\nNum Peaks: 2\n10 100\n20 50\n\n"
        # Rejected by the old clean-peaks gap predicate; it must not increment the ID.
        "Name: rejected\nPrecursorMZ: 100\nNum Peaks: 1\n99.5 100\n\n"
        "Name: experimental\nSpectrum_type: MS2\nPrecursorMZ: 100\n"
        "Instrument_type: QTOF\nNum Peaks: 2\n10 100\n20 50\n\n",
        encoding="utf-8",
    )
    assert list(migration.iter_recovered_ids(msp)) == [
        ("MoNA-export-All_Spectra_0", "quarantined"),
        ("MoNA-export-All_Spectra_1", "predicted"),
    ]

    db = tmp_path / "library.db"
    with sqlite3.connect(db) as conn:
        conn.execute(
            "CREATE TABLE Spectra(id INTEGER PRIMARY KEY, spectrum_id TEXT UNIQUE, "
            "dataset_source TEXT, source_tier TEXT, quality_flags_json TEXT)"
        )
        conn.executemany(
            "INSERT INTO Spectra(spectrum_id,dataset_source,source_tier,quality_flags_json) "
            "VALUES (?,'MoNA-export','experimental','[]')",
            [(f"MoNA-export-All_Spectra_{index}",) for index in range(3)],
        )
    result = migration.migrate(db, msp)

    assert result["database_matched"] == {"predicted": 1, "quarantined": 1}
    assert result["updated"] == 2
    with sqlite3.connect(db) as conn:
        tiers = dict(conn.execute("SELECT spectrum_id,source_tier FROM Spectra"))
        assert tiers["MoNA-export-All_Spectra_0"] == "quarantined"
        assert tiers["MoNA-export-All_Spectra_1"] == "predicted"
        assert tiers["MoNA-export-All_Spectra_2"] == "experimental"


def test_web_temp_uuid_does_not_change_upload_identity(tmp_path: Path):
    from matching import library_import

    first = tmp_path / "a1b2c3d4_literature.csv"
    retry = tmp_path / "deadbeef_literature.csv"
    content = b"name,precursor_mz\nexample,101.0\n"
    first.write_bytes(content)
    retry.write_bytes(content)

    first_identity = library_import._stable_upload_identity(first)
    retry_identity = library_import._stable_upload_identity(retry)

    assert first_identity == retry_identity
    assert first_identity[0] == "literature.csv"


def test_numeric_intensities_without_provenance_are_held_out_of_topk_import(
    monkeypatch, tmp_path: Path
):
    from config import settings
    from matching import library_import

    record = _measured_ethanol_record()
    record.pop("source_tier")
    record.pop("instrument_type")
    sanitized, quality = assess_records([record])
    assert len(sanitized) == 1
    assert sanitized[0]["source_tier"] == "user_unverified"
    assert quality["accepted"] == 1
    assert quality["provenance_held"] == 1
    assert "source_provenance_missing_not_experimental" in quality["warning_counts"]

    full_db = tmp_path / "must_not_be_created.db"
    monkeypatch.setattr(settings, "DB_PATH", full_db)
    monkeypatch.setattr(settings, "DB_PATH_SMALL", tmp_path / "small.db")
    monkeypatch.setattr(library_import, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(
        library_import,
        "_files_to_records",
        lambda _paths, **_kwargs: (
            [record], [{"file": "unverified.xlsx", "written": 1}]
        ),
    )

    result = library_import.import_files(
        [tmp_path / "unverified.xlsx"], db_mode="both"
    )
    assert result["spectral_records"] == 0
    assert result["literature_target_records"] == 0
    assert result["provenance_held_records"] == 1
    assert result["database_modified"] is False
    assert "缺少可验证实验来源" in result["error"]
    assert not full_db.exists()


def test_spreadsheet_can_explicitly_supply_trusted_provenance(tmp_path: Path):
    from csv_excel_to_mgf import parse_file

    csv_path = tmp_path / "provenanced.csv"
    csv_path.write_text(
        "spectrum_name,precursor_mz,formula,smiles,adduct,instrument,source_tier,doi,"
        "fragment_mz_1,fragment_intensity_1,fragment_mz_2,fragment_intensity_2,"
        "fragment_mz_3,fragment_intensity_3\n"
        "Ethanol,47.049141,C2H6O,CCO,[M+H]+,QTOF,experimental,10.1000/example,"
        "19,100,29,75,45,50\n",
        encoding="utf-8",
    )
    records, detected = parse_file(csv_path)
    assert detected == "standard"
    assert records[0]["source_tier"] == "experimental"
    assert records[0]["instrument"] == "QTOF"
    assert records[0]["doi"] == "10.1000/example"

    sanitized, quality = assess_records(records)
    assert quality["provenance_held"] == 0
    assert sanitized[0]["source_tier"] == "experimental"


def test_massbank_supplement_remains_trusted_through_quality_gate():
    record = _measured_ethanol_record()
    record["source_tier"] = "massbank_supplement"
    record["source_accession"] = "MSBNK-TEST-0001"

    clean, audit = assess_record(record)
    assert audit["accepted"] is True
    assert clean["source_tier"] == "massbank_supplement"
    assert source_tier_is_searchable(clean["source_tier"]) is True

    sanitized, quality = assess_import_records([record])
    assert sanitized[0]["source_tier"] == "massbank_supplement"
    assert quality["provenance_held"] == 0


def test_nonexperimental_source_tiers_are_excluded_from_library_topk(
    monkeypatch, tmp_path: Path
):
    """Exercise the real SQLite loader, not only the tier predicate helper."""
    from matching import library_import
    search = importlib.import_module("matching.search")

    db_path = tmp_path / "source_tiers.db"
    with sqlite3.connect(db_path) as conn:
        library_import._create_schema(conn)
        for index, tier in enumerate(
            (
                "experimental", "massbank_supplement", "predicted",
                PRESENCE_ONLY_SOURCE_TIER, QUARANTINED_SOURCE_TIER,
            ),
            start=1,
        ):
            compound_id = conn.execute(
                "INSERT INTO Compounds (compound_name, molecular_formula, smiles, inchikey, "
                "exact_mass, classification) VALUES (?,?,?,?,?,?)",
                (
                    f"compound_{tier}", "C2H6O", "CCO",
                    f"AAAAAAAAAAAAA{index}-BBBBBBBBBB-C", 46.041865, "other",
                ),
            ).lastrowid
            spectrum_id = conn.execute(
                "INSERT INTO Spectra (spectrum_id, source_file, precursor_mz, adduct, "
                "peaks_json, dataset_source, source_tier, intensity_origin, "
                "quality_flags_json) VALUES (?,?,?,?,?,?,?,?,?)",
                (
                    f"spec_{tier}", "synthetic.mgf", 47.049141, "[M+H]+",
                    json.dumps([[20.0, 10.0], [30.0, 20.0], [40.0, 30.0]]),
                    "synthetic", tier, "measured", "[]",
                ),
            ).lastrowid
            conn.execute(
                "INSERT INTO Mapping (compound_id, spectrum_id) VALUES (?,?)",
                (compound_id, spectrum_id),
            )
        low_compound_id = conn.execute(
            "INSERT INTO Compounds (compound_name, molecular_formula, smiles, "
            "inchikey, exact_mass, classification) VALUES (?,?,?,?,?,?)",
            (
                "low_information_massbank", "C2H6O", "CCO",
                "CCCCCCCCCCCCCC-DDDDDDDDDD-C", 46.041865, "other",
            ),
        ).lastrowid
        low_spectrum_id = conn.execute(
            "INSERT INTO Spectra (spectrum_id, source_file, precursor_mz, "
            "adduct, peaks_json, dataset_source, source_tier, "
            "intensity_origin, quality_flags_json) VALUES (?,?,?,?,?,?,?,?,?)",
            (
                "spec_low_information_massbank", "synthetic.mgf", 47.049141,
                "[M+H]+", json.dumps([[20.0, 10.0], [30.0, 20.0]]),
                "synthetic", "massbank_supplement", "measured",
                json.dumps(["fewer_than_3_fragment_peaks"]),
            ),
        ).lastrowid
        conn.execute(
            "INSERT INTO Mapping (compound_id, spectrum_id) VALUES (?,?)",
            (low_compound_id, low_spectrum_id),
        )
        coarse_compound_id = conn.execute(
            "INSERT INTO Compounds (compound_name, molecular_formula, smiles, "
            "inchikey, exact_mass, classification) VALUES (?,?,?,?,?,?)",
            (
                "coarse_precursor_massbank", "C2H6O", "CCO",
                "EEEEEEEEEEEEEE-FFFFFFFFFF-E", 46.041865, "other",
            ),
        ).lastrowid
        coarse_spectrum_id = conn.execute(
            "INSERT INTO Spectra (spectrum_id, source_file, precursor_mz, "
            "adduct, peaks_json, dataset_source, source_tier, "
            "intensity_origin, quality_flags_json) VALUES (?,?,?,?,?,?,?,?,?)",
            (
                "spec_coarse_precursor_massbank", "synthetic.mgf", 47.0,
                "[M+H]+",
                json.dumps([[20.0, 10.0], [30.0, 20.0], [40.0, 30.0]]),
                "synthetic", "massbank_supplement", "measured",
                json.dumps([
                    "theoretical_precursor_index_required",
                    "theoretical_precursor_mz=47.049141",
                ]),
            ),
        ).lastrowid
        conn.execute(
            "INSERT INTO Mapping (compound_id, spectrum_id) VALUES (?,?)",
            (coarse_compound_id, coarse_spectrum_id),
        )
        conn.commit()

    monkeypatch.setattr(search, "_lib_records", None)
    monkeypatch.setattr(search, "_mz_index", None)
    monkeypatch.setattr(search, "_nm_index", None)
    monkeypatch.setattr(search, "_lib_db_path", None)
    records, _index = search.load_library(db_path)

    assert [record["source_tier"] for record in records] == [
        "experimental", "massbank_supplement", "massbank_supplement",
    ]
    assert [record["compound_name"] for record in records] == [
        "compound_experimental", "compound_massbank_supplement",
        "coarse_precursor_massbank",
    ]
    coarse = records[-1]
    assert coarse["precursor_mz"] == 47.0
    assert coarse["search_precursor_mz"] == 47.049141
    candidates = search._query_candidates(
        47.049141, records, _index, tol_da=0.0001, tol_ppm=5.0
    )
    assert coarse["spectrum_id"] in {
        records[index]["spectrum_id"] for index in candidates
    }
    assert source_tier_is_searchable("experimental") is True
    assert source_tier_is_searchable("massbank_supplement") is True
    for tier in ("predicted", PRESENCE_ONLY_SOURCE_TIER, QUARANTINED_SOURCE_TIER):
        assert source_tier_is_searchable(tier) is False


def test_massbank_inchi_identity_wins_when_smiles_omits_stereochemistry():
    record = {
        "name": "L-lactic acid",
        "formula": "C3H6O3",
        "smiles": "CC(O)C(=O)O",
        "inchi": (
            "InChI=1S/C3H6O3/c1-2(4)3(5)6/"
            "h2,4H,1H3,(H,5,6)/t2-/m0/s1"
        ),
        "inchikey": "JVTAAEKCZFNVCJ-REOHCLBHSA-N",
        "precursor_mz": 89.024418,
        "adduct": "[M-H]-",
        "ion_mode": "NEGATIVE",
        "peaks": [[43.0, 100.0], [45.0, 80.0], [71.0, 40.0]],
        "source_tier": "massbank_supplement",
        "intensity_origin": "measured",
        "instrument_type": "LC-ESI-QTOF",
        "source_accession": "MSBNK-TEST-STEREO",
    }

    clean, audit = assess_record(
        record,
        require_experimental_provenance=True,
    )

    assert audit["accepted"] is True
    assert clean["source_tier"] == "massbank_supplement"
    assert clean["inchikey"] == "JVTAAEKCZFNVCJ-REOHCLBHSA-N"
    assert clean["structure_identity_source"] == "inchi_derived"
    assert (
        "smiles_stereochemistry_incomplete_inchi_identity_preferred"
        in clean["quality_flags"]
    )


@pytest.mark.parametrize(
    ("source_tier", "compound_cap", "class_cap"),
    [
        ("experimental", "High", "High"),
        ("predicted", "Medium", "Medium"),
        (PRESENCE_ONLY_SOURCE_TIER, "Low", "Unknown"),
        (QUARANTINED_SOURCE_TIER, "Unknown", "Unknown"),
    ],
)
def test_honest_confidence_obeys_source_tier_cap(
    source_tier: str, compound_cap: str, class_cap: str
):
    best = {
        "compound_name": "Synthetic high scorer",
        "smiles": "CCO",
        "has_structure": True,
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "formula": "C2H6O",
        "classification": "flavonoid",
        "display_class": "flavonoid",
        "cosine": 0.99,
        "final_score": 0.90,
        "reverse_match": 0.95,
        "mass_ppm": 1.0,
        "source_tier": source_tier,
    }

    calibrated = honest_annotate(best, [])

    assert calibrated["compound_confidence"] == compound_cap
    assert calibrated["class_confidence"] == class_cap
    if source_tier == PRESENCE_ONLY_SOURCE_TIER:
        assert calibrated["honesty_tier"] == "literature_target_lead"
        assert "无实测强度" in calibrated["honest_headline"]


def test_theoretical_precursor_rescue_is_capped_at_medium():
    best = {
        "compound_name": "Legacy rounded precursor",
        "smiles": "CCO",
        "has_structure": True,
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "formula": "C2H6O",
        "classification": "flavonoid",
        "display_class": "flavonoid",
        "cosine": 0.99,
        "final_score": 0.95,
        "reverse_match": 0.95,
        "mass_ppm": 0.5,
        "source_tier": "massbank_supplement",
        "lib_precursor_mz_basis": "structure_adduct_theoretical_index",
    }

    calibrated = honest_annotate(best, [])

    assert calibrated["compound_confidence"] == "Medium"
    assert calibrated["class_confidence"] == "Medium"

    manual = dict(best)
    manual.pop("lib_precursor_mz_basis")
    manual["quality_flags"] = [
        "precursor_adduct_conflict_manual_review_required"
    ]
    manual_calibrated = honest_annotate(manual, [])
    assert manual_calibrated["compound_confidence"] == "Low"
    assert manual_calibrated["class_confidence"] == "Medium"


def test_explicit_medium_confidence_cap_survives_high_spectral_score():
    best = {
        "compound_name": "Digitized literature figure spectrum",
        "smiles": "CCO",
        "has_structure": True,
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "formula": "C2H6O",
        "classification": "flavonoid",
        "display_class": "flavonoid",
        "cosine": 0.99,
        "final_score": 0.95,
        "reverse_match": 0.95,
        "mass_ppm": 0.5,
        "source_tier": "experimental",
        "quality_flags": [
            "literature_figure_digitized_relative_intensity",
            "confidence_cap_medium",
        ],
    }

    calibrated = honest_annotate(best, [])

    assert calibrated["compound_confidence"] == "Medium"
    assert calibrated["isomer_group_confidence"] == "Medium"
    assert calibrated["class_confidence"] == "Medium"


def test_literature_target_screening_uses_precursor_polarity_and_diagnostic_ions():
    results = [{
        "spectrum_id": "sample_scan_17",
        "precursor_mz": 301.1006,  # 2 ppm from the target
        "retention_time": 4.2,
        "ion_mode": "POS",
        "raw_peaks": [[100.01, 100.0], [150.02, 80.0], [220.0, 20.0]],
    }]
    shared = {
        "formula": "C15H8O7",
        "precursor_mz": 301.1,
        "adduct": "[M+H]+",
        "diagnostic_ions": [100.0, 150.0, 250.0],
        "plant_source": "Synthetic herb",
    }
    targets = [
        {**shared, "target_id": "target_hit", "compound_name": "Target hit", "polarity": "POS"},
        {**shared, "target_id": "wrong_mode", "compound_name": "Wrong mode", "polarity": "NEG"},
        {
            **shared,
            "target_id": "wrong_precursor",
            "compound_name": "Wrong precursor",
            "polarity": "POS",
            "precursor_mz": 302.1,
        },
        {
            **shared,
            "target_id": "no_diagnostic_ions",
            "compound_name": "No diagnostic ions",
            "polarity": "POS",
            "diagnostic_ions": [],
        },
        {
            **shared,
            "target_id": "one_diagnostic_ion",
            "compound_name": "Only one diagnostic ion",
            "polarity": "POS",
            "diagnostic_ions": [100.0],
        },
    ]

    screened = screen_results(
        results,
        targets,
        precursor_ppm=20.0,
        fragment_tolerance_da=0.05,
        min_diagnostic_ions=2,
    )

    assert screened["matched_target_rows"] == 1
    assert screened["matched_unique_compounds"] == 1
    assert screened["matched_distinguishable_groups"] == 1
    assert screened["screening_funnel"]["precursor_matched_unique_compounds"] == 1
    assert screened["interpretation"].startswith("Target-screening evidence only")
    assert [match["target_id"] for match in screened["matches"]] == ["target_hit"]
    hit = screened["matches"][0]
    assert hit["matched_diagnostic_ions"] == [100.0, 150.0]
    assert hit["diagnostic_coverage"] == pytest.approx(0.667, abs=0.001)
    assert hit["evidence_level"] == "literature_target_screening_only"


def test_literature_target_screen_collapses_indistinguishable_isomer_names():
    results = [{
        "spectrum_id": "sample_scan_isomers",
        "precursor_mz": 353.0878,
        "ion_mode": "NEG",
        "query_peaks": [[191.0345, 100.0], [179.0345, 80.0], [135.045, 20.0]],
    }]
    shared = {
        "formula": "C16H18O9",
        "precursor_mz": 353.0878,
        "adduct": "[M-H]-",
        "polarity": "NEG",
        "diagnostic_ions": [191.0345, 179.0345, 135.045],
        "plant_source": "Herb A",
    }
    targets = [
        {**shared, "target_id": "chlorogenic", "compound_name": "Chlorogenic acid"},
        {**shared, "target_id": "cryptochlorogenic", "compound_name": "Cryptochlorogenic acid"},
    ]

    screened = screen_results(results, targets)

    assert screened["matched_target_rows"] == 2
    assert screened["matched_unique_compounds"] == 2
    assert screened["matched_distinguishable_groups"] == 1
    assert len(screened["matches"]) == 1
    group = screened["matches"][0]
    assert group["isomer_ambiguous"] is True
    assert group["nominal_compound_count"] == 2
    assert group["ambiguous_target_names"] == ["Chlorogenic acid", "Cryptochlorogenic acid"]


def test_literature_target_screen_prefers_denoised_query_peaks_over_raw_noise():
    results = [{
        "spectrum_id": "dense_raw_scan",
        "precursor_mz": 301.1,
        "ion_mode": "POS",
        "query_peaks": [[100.0, 100.0]],
        "raw_peaks": [[100.0, 100.0], [150.0, 1.0]],
    }]
    targets = [{
        "target_id": "noise_only_second_ion",
        "compound_name": "Noise-sensitive target",
        "formula": "C15H8O7",
        "precursor_mz": 301.1,
        "adduct": "[M+H]+",
        "polarity": "POS",
        "diagnostic_ions": [100.0, 150.0],
    }]

    screened = screen_results(results, targets, min_diagnostic_ions=2)

    assert screened["matched_target_rows"] == 0
    assert screened["screening_funnel"]["one_diagnostic_ion_unique_compounds"] == 1


def test_one_observed_peak_cannot_satisfy_two_nearby_diagnostic_ions():
    results = [{
        "spectrum_id": "single_centroid",
        "precursor_mz": 301.1,
        "ion_mode": "POS",
        "query_peaks": [[50.04, 100.0]],
    }]
    targets = [{
        "target_id": "two_close_ions",
        "compound_name": "Close-ion target",
        "formula": "C15H8O7",
        "precursor_mz": 301.1,
        "adduct": "[M+H]+",
        "polarity": "POS",
        "diagnostic_ions": [50.00, 50.08],
    }]

    screened = screen_results(
        results,
        targets,
        fragment_tolerance_da=0.05,
        min_diagnostic_ions=2,
    )

    assert screened["matched_target_rows"] == 0
    assert screened["screening_funnel"]["one_diagnostic_ion_unique_compounds"] == 1


def test_legacy_results_backfill_single_task_polarity_but_reject_mixed_mode():
    from web.server import _with_task_polarity_context

    rows = [{"spectrum_id": "legacy", "precursor_mz": 301.1, "ion_mode": None}]
    repaired, source = _with_task_polarity_context(
        rows,
        {"ion_mode": "NEG", "ion_mode_counts": {"POS": 0, "NEG": 10}},
    )
    assert repaired[0]["ion_mode"] == "NEG"
    assert rows[0]["ion_mode"] is None
    assert source == "task_level_single_polarity_backfill"

    with pytest.raises(ValueError, match="mixed/unknown"):
        _with_task_polarity_context(
            rows,
            {"ion_mode": "MIXED", "ion_mode_counts": {"POS": 5, "NEG": 5}},
        )
