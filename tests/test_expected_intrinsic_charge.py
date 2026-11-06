"""Charged Corydalis targets retain ionic chemistry during MS1 screening."""
import sqlite3
from matching.expected_compounds import (
    ELECTRON_MASS_DA, _target_structure_evidence, _formula_only_intrinsic_charge,
    build_expected_ion_hypotheses, load_expected_compound_catalog,
)
from matching.library_quality import formula_exact_mass


def test_quaternary_structure_keeps_elemental_formula_and_charge():
    evidence = _target_structure_evidence("C[N+](C)(C)CCO")
    assert evidence["formula"] == "C5H14NO"
    assert evidence["formal_charge"] == 1


def test_intrinsic_monocation_generates_one_correct_positive_ion_only():
    target = {
        "target_id": "test-choline", "compound_name": "Choline",
        "formula": "C5H14NO", "formal_charge": 1, "adduct": "[M]+",
        "quality_flags": ["calculated_ion_hypothesis"],
        "precursor_mz": formula_exact_mass("C5H14NO") - ELECTRON_MASS_DA,
    }
    positive = build_expected_ion_hypotheses([target], "POS")
    assert len(positive) == 1
    assert positive[0]["adduct"] == "[M]+"
    assert abs(positive[0]["target_mz"] - 104.1069905) < 1e-6
    assert positive[0]["members"][0]["formal_charge"] == 1
    assert build_expected_ion_hypotheses([target], "NEG") == []


def test_intrinsic_anion_adds_electron_and_excludes_positive_mode():
    target = {"compound_name": "Acetate", "formula": "C2H3O2", "formal_charge": -1}
    negative = build_expected_ion_hypotheses([target], "NEG")
    assert len(negative) == 1
    assert abs(negative[0]["target_mz"] - (formula_exact_mass("C2H3O2") + ELECTRON_MASS_DA)) < 1e-8
    assert negative[0]["adduct"] == "[M]-"
    assert build_expected_ion_hypotheses([target], "POS") == []


def test_neutral_targets_keep_existing_adduct_generation():
    target = {"compound_name": "Glucose", "formula": "C6H12O6", "formal_charge": 0}
    assert {x["adduct"] for x in build_expected_ion_hypotheses([target], "POS")} == {
        "[M+H]+", "[M+Na]+", "[M+NH4]+"
    }


def test_formula_only_isomer_group_retains_checked_intrinsic_charge():
    target = {
        "formula": "C21H26NO4", "adduct": "[M]+",
        "source_tier": "literature_presence_only",
        "quality_flags": ["calculated_ion_hypothesis", "formal_charge:1"],
        "precursor_mz": formula_exact_mass("C21H26NO4") - ELECTRON_MASS_DA,
    }
    assert _formula_only_intrinsic_charge(target) == 1
    target["formal_charge"] = _formula_only_intrinsic_charge(target)
    assert build_expected_ion_hypotheses([target], "NEG") == []
    ions = build_expected_ion_hypotheses([target], "POS")
    assert len(ions) == 1 and ions[0]["adduct"] == "[M]+"


def test_formula_only_charge_requires_mass_and_adduct_agreement():
    target = {
        "formula": "C21H26NO4", "adduct": "[M+H]+",
        "source_tier": "literature_presence_only",
        "quality_flags": ["calculated_ion_hypothesis", "formal_charge:1"],
        "precursor_mz": formula_exact_mass("C21H26NO4") - ELECTRON_MASS_DA,
    }
    assert _formula_only_intrinsic_charge(target) == 0
    target["adduct"] = "[M]+"
    target["precursor_mz"] += 1.0
    assert _formula_only_intrinsic_charge(target) == 0


def test_arbitrary_quality_flags_do_not_supply_intrinsic_charge():
    target = {
        "formula": "C21H26NO4", "adduct": "[M]+", "source_tier": "unknown",
        "quality_flags": ["calculated_ion_hypothesis", "formal_charge:1"],
        "precursor_mz": formula_exact_mass("C21H26NO4") - ELECTRON_MASS_DA,
    }
    assert _formula_only_intrinsic_charge(target) == 0


def test_charged_database_structure_is_not_quarantined(tmp_path):
    path = tmp_path / "catalog.db"
    c = sqlite3.connect(path)
    columns = (
        "target_id dataset_source compound_name molecular_formula smiles inchi "
        "inchikey reported_inchikey structure_identity_verified "
        "structure_identity_source exact_mass precursor_mz adduct charge "
        "plant_source source_file source_content_sha256 dataset_version "
        "source_row_index provenance doi dataset_accession spectrum_accession "
        "quality_flags_json source_tier"
    ).split()
    c.execute("CREATE TABLE TargetCompounds (id INTEGER PRIMARY KEY," +
              ",".join(name + " TEXT" for name in columns) + ")")
    c.execute("INSERT INTO TargetCompounds (target_id,dataset_source,compound_name,"
              "molecular_formula,smiles,source_tier) VALUES (?,?,?,?,?,?)",
              ("choline", "charged", "Choline", "C5H14NO", "C[N+](C)(C)CCO",
               "literature_presence_only"))
    c.commit()
    sources = [x[0] for x in c.execute("SELECT DISTINCT dataset_source FROM TargetCompounds")]
    c.close()
    targets = load_expected_compound_catalog(path, sources)
    assert targets
    assert all(x["structure_validation_status"] == "internally_consistent" for x in targets)
    assert all(x["formal_charge"] == 1 for x in targets)
    assert build_expected_ion_hypotheses(targets, "NEG") == []
