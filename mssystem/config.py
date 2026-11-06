from __future__ import annotations

from pathlib import Path

SUPPORTED_INPUT_EXTS = {".mgf", ".msp", ".mzxml"}

METADATA_ALIASES = {
    "precursor_mz": [
        "precursor_mz",
        "PEPMASS",
        "pep_mass",
        "PrecursorMZ",
        "PRECURSOR_MZ",
        "parent_mass",
        "PeptideMass",
        "mz",
    ],
    "retention_time": [
        "retention_time",
        "RTINSECONDS",
        "rt",
        "RETENTIONTIME",
        "scan_start_time",
        "RetTime",
        "rtinseconds",
        "retentionTime",
    ],
    "charge": [
        "charge",
        "CHARGE",
        "z",
        "precursor_charge",
        "ChargeState",
    ],
    "collision_energy": [
        "collision_energy",
        "collisionEnergy",
        "COLLISIONENERGY",
        "Collision_energy",
        "ce",
        "activation_energy",
        "activationenergy",
        "hcd_energy",
        "HCD",
        "collision voltage",
    ],
    "collision_energy_unit": [
        "collision_energy_unit",
        "collisionEnergyUnit",
        "COLLISION_ENERGY_UNIT",
        "ce_unit",
        "energy_unit",
    ],
    "adduct": ["adduct", "Adduct", "IONTYPE", "ion_type", "precursor_type"],
    "name": ["name", "NAME", "compound_name", "title", "TITLE", "Name"],
    "smiles": ["smiles", "SMILES", "canonical_smiles", "structure_smiles"],
    "inchi": ["inchi", "INCHI", "InChI"],
    "inchikey": ["inchikey", "INCHIKEY", "InChIKey", "INCHIAUX"],
    "formula": ["formula", "FORMULA", "molecular_formula", "MolecularFormula"],
    "instrument": ["instrument", "INSTRUMENT", "instrument_model"],
    "instrument_type": [
        "instrument_type",
        "INSTRUMENT_TYPE",
        "mass_analyzer_type",
    ],
    "ion_mode": ["ion_mode", "ionmode", "ION_MODE", "polarity"],
    "source_accession": [
        "source_accession",
        "spectrum_id",
        "accession",
        "db#",
        "DB#",
    ],
    "compound_class": [
        "compound_class",
        "USER_CLASS",
        "user_class",
        "compoundclass",
        "class",
        "Classification",
    ],
}

REQUIRED_COLUMNS = [
    "spectrum_id",
    "source_file",
    "name",
    "precursor_mz",
    "retention_time",
    "charge",
    "collision_energy",
    "collision_energy_raw",
    "collision_energy_unit",
    "collision_energy_mode",
    "adduct",
    "instrument_type",
    "instrument",
    "ion_mode",
    "source_accession",
    "smiles",
    "inchi",
    "inchikey",
    "formula",
    "peaks",
]

DEFAULT_DB_PATH = Path("data/local_library.db")
DEFAULT_MAPPING_CSV = Path("output/library_mapping.csv")
DEFAULT_MODEL_PATH = Path("models/classifier.joblib")
DEFAULT_TARGET_DB = Path("output/known_compound_targets.csv")
