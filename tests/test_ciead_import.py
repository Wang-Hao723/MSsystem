import json
import sqlite3

from rdkit import Chem
from rdkit.Chem import Descriptors

from scripts import import_ciead_library as ciead


def _structure(name="cis-Resveratrol"):
    smiles = "OC1=CC(O)=C(C=C1)/C=C/C2=CC(O)=CC=C2"
    molecule = Chem.MolFromSmiles(smiles)
    return ciead.StructureRecord(
        name=name,
        aliases=("Resveratrol",),
        source_formula="C14H12O3",
        formula="C14H12O3",
        smiles=Chem.MolToSmiles(molecule, isomericSmiles=True),
        inchi=Chem.MolToInchi(molecule),
        inchikey=Chem.MolToInchiKey(molecule),
        exact_mass=float(Descriptors.ExactMolWt(molecule)),
        compound_class_raw="Stilbenes",
        superclass_raw="Stilbenoids",
        classification_v3="stilbene",
        classification_source="curated_ciead_npclassifier",
        classification_confidence=0.97,
        structure_source="CIeaD_Table_S2",
    )


def _candidate():
    structure = _structure()
    peaks = [[91.0542, 100.0], [107.0491, 60.0], [135.0440, 40.0]]
    fingerprint = ciead.normalized_peak_fingerprint(peaks)
    return ciead.Candidate(
        source_name=structure.name,
        structure=structure,
        precursor_mz=structure.exact_mass + 1.007276,
        polarity="POS",
        adduct="[M+H]+",
        charge=1,
        collision_energy=20.0,
        collision_energy_raw="20",
        peaks=peaks,
        peak_fingerprint=fingerprint,
        source_accession="CIeaD:test",
        spectrum_id="CIeaD:test",
        activation_method="unresolved",
        source_memberships=("Neg CID MSP", "Pos EAD MSP"),
        quality_flags=["activation_method_unresolved_cid_ead_export_overlap"],
        metadata_match_method="exact_name_formula",
        formula_evidence="C14H12O3",
    )


def test_name_normalization_and_formula_safe_matching():
    hidden = "(-\u200b)\u200b-\u200balpha\u2011Terpineol"
    assert ciead.normalized_name(hidden) == "(-)-alpha-terpineol"
    assert ciead.simple_name("Quercetin 3-O-\u03b2-D-glucoside") == (
        "quercetin3obetadglucoside"
    )
    assert ciead.neutral_formula_for_comparison("C20H20NO4+") == "C20H20NO4"
    assert ciead.neutral_formula_for_comparison("C17H23NO3.1/2H2SO4") == ""


def test_adduct_inference_is_structure_and_polarity_constrained():
    structure = _structure()
    positive, positive_ppm = ciead.infer_adduct(
        structure.exact_mass,
        structure.exact_mass + 1.007276,
        "POS",
    )
    negative, negative_ppm = ciead.infer_adduct(
        structure.exact_mass,
        structure.exact_mass - 1.007276,
        "NEG",
    )
    assert positive == "[M+H]+"
    assert negative == "[M-H]-"
    assert positive_ppm < 0.01
    assert negative_ppm < 0.01
    rejected, _ = ciead.infer_adduct(
        structure.exact_mass,
        structure.exact_mass + 11.0,
        "NEG",
    )
    assert rejected is None


def test_structure_resolver_uses_formula_before_fuzzy_name():
    target = _structure("Quercetin 3-O-beta-D-glucoside")
    other = _structure("Unrelated same-name-like isomer")
    other.source_formula = "C15H10O7"
    other.formula = "C15H10O7"
    resolver = ciead.StructureResolver([target, other])
    resolved, method, details = resolver.resolve(
        "Quercetin-3-O-\u03b2-D-glucoside",
        "C14H12O3",
    )
    assert resolved is target
    assert method in {"normalized_name_formula", "fuzzy_name_formula"}
    assert isinstance(details, dict)


def test_database_novel_selection_rejects_existing_content():
    candidate = _candidate()
    key = (
        candidate.structure.inchikey,
        candidate.polarity,
        candidate.peak_fingerprint,
    )
    index = {key: [candidate.precursor_mz]}
    selected, stats = ciead.select_database_novel([candidate], index)
    assert selected == []
    assert stats["duplicate_existing_database_content"] == 1


def test_apply_import_writes_structure_spectrum_and_mapping(tmp_path, monkeypatch):
    db = tmp_path / "library.db"
    connection = sqlite3.connect(db)
    connection.executescript(
        """
        CREATE TABLE Compounds (
          id INTEGER PRIMARY KEY, compound_name TEXT, molecular_formula TEXT,
          smiles TEXT, inchi TEXT, inchikey TEXT, exact_mass REAL,
          classification TEXT, compound_class TEXT, classification_v3 TEXT,
          classification_source TEXT, classification_confidence REAL,
          fragmentation_pathway TEXT, reported_inchikey TEXT,
          structure_identity_verified INTEGER DEFAULT 0,
          structure_identity_source TEXT
        );
        CREATE TABLE Spectra (
          id INTEGER PRIMARY KEY, spectrum_id TEXT UNIQUE, source_file TEXT,
          precursor_mz REAL, retention_time REAL, collision_energy REAL,
          charge INTEGER, adduct TEXT, peaks_json TEXT, dataset_source TEXT,
          source_tier TEXT, intensity_origin TEXT, quality_flags_json TEXT,
          collision_energy_raw TEXT, collision_energy_unit TEXT,
          collision_energy_mode TEXT, collision_energy_min REAL,
          collision_energy_max REAL, instrument_type TEXT, instrument TEXT,
          ion_mode TEXT, source_accession TEXT
        );
        CREATE TABLE Mapping (
          id INTEGER PRIMARY KEY, compound_id INTEGER, spectrum_id INTEGER,
          score REAL
        );
        """
    )
    connection.commit()
    connection.close()

    monkeypatch.setattr(ciead, "LIBRARY_WRITE_ENABLED", True)
    candidate = _candidate()
    result, ledger = ciead.apply_import(
        db,
        [candidate],
        preferred={},
        rows_by_id={},
        compound_updates=[],
    )
    assert result["inserted_compounds"] == 1
    assert result["inserted_spectra"] == 1
    assert {row["action"] for row in ledger} == {
        "insert_compound",
        "insert_spectrum",
    }

    connection = sqlite3.connect(db)
    spectrum = connection.execute(
        "SELECT source_tier, ion_mode, quality_flags_json FROM Spectra"
    ).fetchone()
    mapping = connection.execute(
        "SELECT COUNT(*) FROM Mapping"
    ).fetchone()[0]
    connection.close()
    assert spectrum[0] == "experimental"
    assert spectrum[1] == "Positive"
    assert "activation_method_unresolved" in json.loads(spectrum[2])[0]
    assert mapping == 1
