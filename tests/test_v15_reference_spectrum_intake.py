from __future__ import annotations

import json
import importlib
import sqlite3

import pytest

from matching.library_quality import (
    PRESENCE_ONLY_SOURCE_TIER,
    assess_import_records,
)
from matching.reference_spectrum_intake import (
    MANIFEST_SCHEMA_VERSION,
    parse_reference_intake,
)


def _manifest(source_file: str, *, image: bool = False, reviewer: bool = True):
    spectrum = {
        "compound": {
            "name": "Ethanol reference",
            "formula": "C2H6O",
            "smiles": "CCO",
        },
        "precursor_mz": 47.049141,
        "adduct": "[M+H]+",
        "ion_mode": "POS",
        "charge": 1,
        "instrument_metadata": {
            "model": "Agilent 6560 QTOF",
            "analyzer": "qtof",
            "activation": "CID",
            "acquisition_method": "method-20260826",
        },
        "collision_energy": {"value": 20.0, "unit": "eV", "mode": "fixed"},
        "reference_standard": {
            "material_type": "authentic_standard",
            "lot": "LOT-001",
            "purity_percent": 99.5,
        },
        "source_provenance_verified": True,
        "provenance": {"verified": True, "source_accession": "LOCAL-STD-001"},
    }
    if image:
        spectrum.update({
            "image_file": source_file,
            "peaks": [[20.0, 15.0], [30.0, 100.0], [40.0, 35.0]],
            "digitization": {
                "method": "manual_axis_calibrated_peak_pick",
                "reviewer_confirmed": reviewer,
            },
        })
    else:
        spectrum["source_file"] = source_file
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "manifest_id": "unit-test-reference-manifest",
        "spectra": [spectrum],
    }


def test_authentic_standard_mgf_with_manifest_becomes_experimental_record(tmp_path):
    mgf = tmp_path / "a1b2c3d4_standard.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=scan-one\nPEPMASS=47.049141\nCHARGE=1+\n"
        "20 15\n30 100\n40 35\nEND IONS\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "e5f6a7b8_reference.json"
    manifest.write_text(
        json.dumps(_manifest("standard.mgf")), encoding="utf-8"
    )

    records, stats = parse_reference_intake([mgf, manifest], strict=True)

    assert len(records) == 1
    record = records[0]
    assert record["name"] == "Ethanol reference"
    assert record["source_file"] == "standard.mgf"
    assert record["source_tier"] == "experimental"
    assert record["intensity_origin"] == "measured"
    assert record["instrument"] == "Agilent 6560 QTOF"
    assert record["activation_method"] == "CID"
    assert record["collision_energy"] == pytest.approx(20.0)
    assert record["collision_energy_unit"] == "eV"
    assert record["collision_energy_mode"] == "fixed"
    assert record["standard_lot"] == "LOT-001"
    assert record["standard_purity"] == pytest.approx(99.5)
    assert len(record["source_asset_sha256"]) == 64
    assert any(item.get("format") == "raw_mgf_reference" for item in stats)

    accepted, quality = assess_import_records(records)
    assert quality["record_types"]["experimental_spectra_accepted"] == 1
    assert quality["record_types"]["expected_targets_accepted"] == 0
    assert accepted[0]["source_tier"] == "experimental"


def test_authentic_standard_msp_with_manifest_becomes_experimental_record(tmp_path):
    msp = tmp_path / "a1b2c3d4_standard.msp"
    msp.write_text(
        "Name: Ethanol reference\n"
        "PrecursorMZ: 47.049141\n"
        "Precursor_type: [M+H]+\n"
        "Ion_mode: Positive\n"
        "Num Peaks: 3\n"
        "20 15\n30 100\n40 35\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "e5f6a7b8_reference.json"
    manifest.write_text(
        json.dumps(_manifest("standard.msp")), encoding="utf-8"
    )

    records, stats = parse_reference_intake([msp, manifest], strict=True)

    assert len(records) == 1
    assert records[0]["name"] == "Ethanol reference"
    assert records[0]["source_file"] == "standard.msp"
    assert records[0]["source_tier"] == "experimental"
    assert records[0]["collision_energy"] == pytest.approx(20.0)
    assert any(item.get("format") == "raw_msp_reference" for item in stats)


def test_same_condition_mgf_replicates_become_one_auditable_consensus(tmp_path):
    mgf = tmp_path / "a1b2c3d4_replicates.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=replicate-one\nPEPMASS=47.049141\nCHARGE=1+\n"
        "20.000 15\n30.000 100\n40.000 35\nEND IONS\n"
        "BEGIN IONS\nTITLE=replicate-two\nPEPMASS=47.049141\nCHARGE=1+\n"
        "20.005 20\n30.004 80\n40.006 40\n45.000 10\nEND IONS\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "e5f6a7b8_reference.json"
    manifest.write_text(
        json.dumps(_manifest("replicates.mgf")), encoding="utf-8"
    )

    records, stats = parse_reference_intake([mgf, manifest], strict=True)

    assert len(records) == 1
    record = records[0]
    assert record["intensity_origin"] == "measured_consensus"
    assert record["consensus_member_count"] == 2
    assert record["independent_evidence_count"] == 1
    assert record["technical_replicates_count_as_independent_evidence"] is False
    assert len(record["peaks"]) == 3
    assert len(record["source_provenance"]["members"]) == 2
    consensus_stats = next(
        item for item in stats
        if item.get("format") == "technical_replicate_consensus"
    )
    assert consensus_stats["consensus"]["consensus_groups_built"] == 1

    accepted, quality = assess_import_records(records)
    assert quality["record_types"]["experimental_spectra_accepted"] == 1
    assert accepted[0]["intensity_origin"] == "measured_consensus"


def test_digitized_image_is_retained_but_cannot_enter_experimental_search(tmp_path):
    image = tmp_path / "1234abcd_spectrum.png"
    image.write_bytes(b"\x89PNG\r\n\x1a\nunit-test-image")
    manifest = tmp_path / "5678abcd_reference.json"
    manifest.write_text(
        json.dumps(_manifest("spectrum.png", image=True)), encoding="utf-8"
    )

    records, _ = parse_reference_intake([image, manifest], strict=True)

    assert len(records) == 1
    record = records[0]
    assert record["source_file"] == "spectrum.png"
    assert record["source_image_file"] == "spectrum.png"
    assert len(record["source_image_sha256"]) == 64
    assert record["intensity_origin"] == "digitized_image"
    assert record["source_tier"] == "user_unverified"

    accepted, quality = assess_import_records(records)
    assert quality["record_types"]["experimental_spectra_accepted"] == 0
    assert quality["record_types"]["expected_targets_accepted"] == 1
    assert accepted[0]["source_tier"] == PRESENCE_ONLY_SOURCE_TIER


def test_digitized_image_requires_human_review_confirmation(tmp_path):
    image = tmp_path / "spectrum.png"
    image.write_bytes(b"image")
    manifest = tmp_path / "reference.json"
    manifest.write_text(
        json.dumps(_manifest("spectrum.png", image=True, reviewer=False)),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="reviewer_confirmed"):
        parse_reference_intake([image, manifest], strict=True)


def test_mzml_without_compound_manifest_fails_closed(tmp_path):
    path = tmp_path / "reference.mzML"
    path.write_text("<mzML></mzML>", encoding="utf-8")

    with pytest.raises(Exception):
        parse_reference_intake([path], strict=True)


def test_raw_reference_without_auditable_identity_cannot_be_confirmed(tmp_path):
    mgf = tmp_path / "anonymous.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=accession-only\nPEPMASS=47.049141\n"
        "INSTRUMENT=Agilent 6560 QTOF\nACCESSION=PUBLIC-1\n"
        "20 15\n30 100\n40 35\nEND IONS\n",
        encoding="utf-8",
    )

    preview_records, _ = parse_reference_intake([mgf], strict=False)
    assert preview_records[0]["source_tier"] == "user_unverified"
    assert "missing_auditable_compound_identity" in preview_records[0][
        "reference_intake_flags"
    ]
    with pytest.raises(ValueError, match="lacks compound name"):
        parse_reference_intake([mgf], strict=True)


def test_stepped_ce_range_is_preserved_without_midpoint_substitution(tmp_path):
    mgf = tmp_path / "stepped.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=scan-one\nPEPMASS=47.049141\nCHARGE=1+\n"
        "20 15\n30 100\n40 35\nEND IONS\n",
        encoding="utf-8",
    )
    payload = _manifest("stepped.mgf")
    payload["spectra"][0]["collision_energy"] = {
        "value": 20.0,
        "minimum": 10.0,
        "maximum": 30.0,
        "unit": "eV",
        "mode": "stepped",
    }
    manifest = tmp_path / "reference.json"
    manifest.write_text(json.dumps(payload), encoding="utf-8")

    records, _ = parse_reference_intake([mgf, manifest], strict=True)

    assert records[0]["collision_energy"] == pytest.approx(20.0)
    assert records[0]["collision_energy_min"] == pytest.approx(10.0)
    assert records[0]["collision_energy_max"] == pytest.approx(30.0)
    assert records[0]["collision_energy_mode"] == "stepped"


def test_reference_and_image_provenance_survive_sqlite_write(
    tmp_path, monkeypatch
):
    from config import settings
    from matching import library_import

    monkeypatch.setattr(settings, "LIBRARY_WRITE_ENABLED", True)
    mgf = tmp_path / "standard.mgf"
    mgf.write_text(
        "BEGIN IONS\nTITLE=scan-one\nPEPMASS=47.049141\nCHARGE=1+\n"
        "20 15\n30 100\n40 35\nEND IONS\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "reference.json"
    manifest.write_text(json.dumps(_manifest("standard.mgf")), encoding="utf-8")
    raw_records, _ = parse_reference_intake([mgf, manifest], strict=True)
    raw_accepted, _ = assess_import_records(raw_records)

    db_path = tmp_path / "reference_intake.db"
    assert library_import.append_records(
        raw_accepted, db_path, ts_prefix="v15_reference"
    ) == (1, 1)
    with sqlite3.connect(db_path) as conn:
        row = conn.execute(
            "SELECT activation_method, acquisition_method, source_asset_sha256, "
            "reference_material_type, standard_lot, standard_purity, provenance_json "
            "FROM Spectra"
        ).fetchone()
    assert row[:2] == ("CID", "method-20260826")
    assert len(row[2]) == 64
    assert row[3:5] == ("authentic_standard", "LOT-001")
    assert row[5] == pytest.approx(99.5)
    assert json.loads(row[6])["source_accession"] == "LOCAL-STD-001"

    search_module = importlib.import_module("matching.search")

    monkeypatch.setattr(search_module, "_lib_records", None)
    monkeypatch.setattr(search_module, "_mz_index", None)
    monkeypatch.setattr(search_module, "_nm_index", None)
    monkeypatch.setattr(search_module, "_lib_db_path", None)
    loaded, _ = search_module.load_library(db_path)
    assert loaded[0]["activation_method"] == "CID"
    assert loaded[0]["collision_energy_mode"] == "fixed"
    assert loaded[0]["acquisition_method"] == "method-20260826"

    image = tmp_path / "spectrum.png"
    image.write_bytes(b"\x89PNG\r\n\x1a\nunit-test-image")
    image_manifest = tmp_path / "image_reference.json"
    image_manifest.write_text(
        json.dumps(_manifest("spectrum.png", image=True)), encoding="utf-8"
    )
    image_records, _ = parse_reference_intake(
        [image, image_manifest], strict=True
    )
    image_accepted, _ = assess_import_records(image_records)
    assert library_import.append_target_records(image_accepted, db_path) == (1, 3)
    with sqlite3.connect(db_path) as conn:
        target = conn.execute(
            "SELECT source_image_file, source_image_sha256, digitization_method, "
            "digitization_reviewer_confirmed, intensity_origin, provenance_json "
            "FROM TargetCompounds"
        ).fetchone()
        spectral_count = conn.execute("SELECT COUNT(*) FROM Spectra").fetchone()[0]
    assert target[:5] == (
        "spectrum.png",
        image_accepted[0]["source_image_sha256"],
        "manual_axis_calibrated_peak_pick",
        1,
        "digitized_image",
    )
    assert json.loads(target[5])["source_asset_sha256"] == image_accepted[0][
        "source_asset_sha256"
    ]
    assert spectral_count == 1
