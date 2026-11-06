from pathlib import Path
import sqlite3

from rdkit import Chem

from matching.ce_profiles import (
    lookup_ce_profile_evidence,
    profile_database_summary,
)
from mssystem.data_loader import DataCleaner
from scripts.build_ce_profiles import build_database, resolve_input_files
from utils.collision_energy import parse_collision_energy


def test_collision_energy_parser_preserves_domains_and_modes():
    assert parse_collision_energy("10(NCE)") == {
        "raw": "10(NCE)",
        "value": 10.0,
        "unit": "NCE",
        "mode": "single",
        "minimum": 10.0,
        "maximum": 10.0,
        "values": [10.0],
        "comparable": True,
    }
    assert parse_collision_energy("20 eV")["unit"] == "eV"
    attached_range = parse_collision_energy("20-40eV")
    assert attached_range["unit"] == "eV"
    assert attached_range["mode"] == "range"
    assert attached_range["minimum"] == 20.0
    assert attached_range["maximum"] == 40.0
    assert parse_collision_energy("6V")["unit"] == "V"
    assert not parse_collision_energy("6V")["comparable"]
    nominal = parse_collision_energy("30 % (nominal)")
    assert nominal["unit"] == "percent"
    assert nominal["unit"] != "NCE"
    ramp = parse_collision_energy("Ramp 19.9-29.9 eV")
    assert ramp["mode"] == "ramp"
    assert ramp["minimum"] == 19.9
    assert ramp["maximum"] == 29.9
    assert ramp["value"] is None
    assert parse_collision_energy("40")["unit"] == ""


def _massbank_record(
    name,
    accession,
    inchi,
    inchikey,
    smiles,
    energy,
    peaks,
    *,
    adduct="[M+H]+",
):
    peak_lines = "\n".join(f"{mz} {intensity}" for mz, intensity in peaks)
    return f"""Name: {name}
DB#: {accession}
InChIKey: {inchikey}
InChI: {inchi}
SMILES: {smiles}
Precursor_type: {adduct}
Spectrum_type: MS2
PrecursorMZ: 195.0877
Instrument_type: LC-ESI-ITFT
Instrument: Q Exactive Test
Ion_mode: POSITIVE
Collision_energy: {energy}(NCE)
Formula: C8H10N4O2
Num Peaks: {len(peaks)}
{peak_lines}

"""


def _synthetic_massbank(path: Path):
    smiles = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
    molecule = Chem.MolFromSmiles(smiles)
    inchi = Chem.MolToInchi(molecule)
    inchikey = Chem.MolToInchiKey(molecule)
    records = []
    for energy, peaks in (
        (10, [(55.0, 20), (82.0, 100), (110.0, 40), (138.0, 10)]),
        (20, [(55.0, 30), (82.0, 70), (110.0, 100), (138.0, 30)]),
        (40, [(55.0, 100), (82.0, 20), (110.0, 30), (138.0, 70)]),
    ):
        for replicate in range(2):
            shifted = [
                (mz + replicate * 0.001, intensity)
                for mz, intensity in peaks
            ]
            records.append(
                _massbank_record(
                    "Caffeine",
                    f"MSBNK-TEST-{energy}-{replicate}",
                    inchi,
                    inchikey,
                    smiles,
                    energy,
                    shifted,
                )
            )
    path.write_text("".join(records), encoding="utf-8")
    return inchikey


def test_build_and_lookup_multi_ce_sidecar(tmp_path):
    source = tmp_path / "massbank.msp"
    inchikey = _synthetic_massbank(source)
    output = tmp_path / "massbank_ce_profiles.db"
    summary = build_database([source], output)

    assert summary["counters"]["source_records"] == 6
    assert summary["counters"]["multi_ce_profiles"] == 1
    assert summary["counters"]["multi_ce_identities"] == 1

    connection = sqlite3.connect(output)
    try:
        assert connection.execute("PRAGMA quick_check").fetchone()[0] == "ok"
        assert connection.execute("SELECT COUNT(*) FROM CEProfiles").fetchone()[0] == 1
        assert connection.execute("SELECT COUNT(*) FROM CEBins").fetchone()[0] == 3
        assert connection.execute(
            "SELECT COUNT(*) FROM CEFragmentTrajectories"
        ).fetchone()[0] >= 3
        assert connection.execute(
            "SELECT quality_grade FROM CEProfiles"
        ).fetchone()[0] == "A"
    finally:
        connection.close()

    evidence = lookup_ce_profile_evidence(
        inchikey=inchikey,
        adduct="[M+H]+",
        ion_mode="POS",
        query_collision_energy={"value": 21, "unit": "NCE"},
        query_peaks=[[55, 0.3], [82, 0.7], [110, 1.0], [138, 0.3]],
        query_precursor_mz=195.0877,
        db_path=output,
    )
    assert evidence["status"] == "available_near_energy"
    assert evidence["runtime_influence"] is False
    assert evidence["score_contribution"] == 0.0
    assert evidence["best_profile"]["closest_energy"] == 20.0
    assert evidence["best_profile"]["consensus_cosine"] > 0.95
    assert evidence["best_profile"]["matched_fragment_trajectories"]
    assert evidence["best_profile"]["source_dataset"].endswith("/massbank")
    assert evidence["best_profile"]["identity_verification"] == "structure_verified"

    mismatch = lookup_ce_profile_evidence(
        inchikey=inchikey,
        adduct="[M+H]+",
        ion_mode="POS",
        query_collision_energy={"value": 20, "unit": "eV"},
        db_path=output,
    )
    assert mismatch["status"] == "no_comparable_profile"

    database_summary = profile_database_summary(output)
    assert database_summary["available"] is True
    assert database_summary["profile_count"] == 1
    assert database_summary["identity_count"] == 1
    assert database_summary["energy_bin_count"] == 3
    assert database_summary["source_file_count"] == 1


def _synthetic_mgf(path: Path, *, reported_only: bool = False):
    smiles = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
    molecule = Chem.MolFromSmiles(smiles)
    inchi = Chem.MolToInchi(molecule)
    inchikey = Chem.MolToInchiKey(molecule)
    blocks = []
    for energy in (10, 20, 40):
        identity = (
            ""
            if reported_only
            else f"SMILES={smiles}\nINCHI={inchi}\n"
        )
        blocks.append(
            f"""BEGIN IONS
TITLE=Community caffeine {energy}
PEPMASS=195.0877
INCHIAUX={inchikey}
{identity}FORMULA=C8H10N4O2
ADDUCT=[M+H]+
IONMODE=Positive
MSLEVEL=2
INSTRUMENT_TYPE=Orbitrap
INSTRUMENT=Community Q Exactive
Collision energy={energy}
COLLISION_ENERGY_UNIT=NCE
SPECTRUMID=COMMUNITY-{energy}
55.0 {20 + energy}
82.0 {100 - energy}
110.0 {30 + energy}
138.0 {10 + energy}
END IONS
"""
        )
    path.write_text("\n".join(blocks), encoding="utf-8")
    return inchikey


def test_multi_source_directory_ingestion_and_source_isolation(tmp_path):
    source_dir = tmp_path / "sources"
    source_dir.mkdir()
    massbank = source_dir / "massbank.msp"
    community = source_dir / "community.mgf"
    inchikey = _synthetic_massbank(massbank)
    _synthetic_mgf(community)
    output = tmp_path / "ce_profiles.db"

    summary = build_database([source_dir], output)
    assert summary["counters"]["source_records"] == 9
    assert summary["counters"]["multi_ce_profiles"] == 2
    assert len(resolve_input_files([source_dir, community])) == 2

    database_summary = profile_database_summary(output)
    assert database_summary["schema_version"] == "2"
    assert database_summary["source_file_count"] == 2
    assert sum(database_summary["profiles_by_source"].values()) == 2

    evidence = lookup_ce_profile_evidence(
        inchikey=inchikey,
        adduct="[M+H]+",
        ion_mode="POS",
        query_collision_energy={"value": 20, "unit": "NCE"},
        db_path=output,
    )
    assert evidence["profile_count"] == 2
    assert len({row["source_dataset"] for row in evidence["profiles"]}) == 2

    cleaned = DataCleaner(top_n_peaks=50).clean_files([community])
    assert len(cleaned) == 3
    assert cleaned[0]["collision_energy"] == 10.0
    assert cleaned[0]["collision_energy_raw"] == "10"
    assert cleaned[0]["collision_energy_unit"] == "NCE"
    assert cleaned[0]["instrument"] == "Community Q Exactive"
    assert cleaned[0]["source_accession"] == "COMMUNITY-10"


def test_cross_file_duplicate_spectra_are_not_counted_as_replicates(tmp_path):
    first = tmp_path / "first.msp"
    second = tmp_path / "second.msp"
    _synthetic_massbank(first)
    second.write_bytes(first.read_bytes())
    output = tmp_path / "deduplicated.db"

    summary = build_database([first, second], output)
    assert summary["counters"]["duplicate_spectra_removed"] == 6
    assert summary["counters"]["multi_ce_profiles"] == 1
    connection = sqlite3.connect(output)
    try:
        assert connection.execute(
            "SELECT n_spectra FROM CEProfiles"
        ).fetchone()[0] == 6
        assert connection.execute(
            "SELECT COUNT(*) FROM ProfileSources"
        ).fetchone()[0] == 2
    finally:
        connection.close()


def test_reported_full_inchikey_is_retained_but_provenance_graded(tmp_path):
    source = tmp_path / "reported.mgf"
    _synthetic_mgf(source, reported_only=True)
    output = tmp_path / "reported.db"
    summary = build_database([source], output)
    assert summary["counters"]["multi_ce_profiles"] == 1
    connection = sqlite3.connect(output)
    try:
        row = connection.execute(
            "SELECT identity_verification, structure_verified_spectra "
            "FROM CEProfiles"
        ).fetchone()
        assert row == ("reported_full_inchikey", 0)
    finally:
        connection.close()

    strict = build_database(
        [source],
        tmp_path / "strict.db",
        dry_run=True,
        allow_reported_inchikey=False,
    )
    assert strict["counters"].get("multi_ce_profiles", 0) == 0
    assert strict["counters"]["identity_not_structure_verified"] == 3


def test_data_cleaner_preserves_massbank_ce_metadata(tmp_path):
    source = tmp_path / "massbank.msp"
    _synthetic_massbank(source)
    records = DataCleaner(top_n_peaks=50).clean_files([source])
    assert records
    first = records[0]
    assert first["collision_energy"] == 10.0
    assert first["collision_energy_raw"] == "10(NCE)"
    assert first["collision_energy_unit"] == "NCE"
    assert first["collision_energy_mode"] == "single"
    assert first["instrument_type"] == "LC-ESI-ITFT"
    assert first["instrument"] == "Q Exactive Test"
    assert first["source_accession"] == "MSBNK-TEST-10-0"
