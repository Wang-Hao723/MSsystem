from pathlib import Path

import build_from_mgf_msp as builder
import import_neg_spectra as importer


def test_massbank_supplement_is_considered_during_duplicate_detection():
    assert "massbank_supplement" in importer.SEARCHABLE_SOURCE_TIERS


def _candidate(
    *,
    accession: str,
    peaks: list[list[float]],
    precursor_mz: float = 300.0,
) -> dict:
    rec = {
        "spectrum_id": f"test-{accession}-{peaks[0][1]}",
        "name": "Same compound name",
        "formula": "C10H10O2",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "precursor_mz": precursor_mz,
        "adduct": "[M-H]-",
        "ion_mode": "Negative",
        "charge": -1,
        "peaks": peaks,
        "source_accession": accession,
    }
    return {
        "rec": rec,
        "path": Path(f"/tmp/{accession}.msp"),
        "source": accession,
        "identity": importer._identity(rec),
        "inchikey": rec["inchikey"],
        "precursor_mz": precursor_mz,
        "polarity": "NEG",
        "peaks": peaks,
        "fragment_count": len(peaks),
        "peak_fp": importer._peak_fp(peaks),
        "accession": importer._normalize_accession(accession),
    }


def test_peak_fingerprint_keeps_relative_intensity_information():
    low_high = [[100.0, 10.0], [150.0, 100.0], [200.0, 20.0]]
    high_low = [[100.0, 100.0], [150.0, 10.0], [200.0, 20.0]]

    assert importer._peak_fp(low_high) != importer._peak_fp(high_low)
    assert importer._peak_fp(low_high) == importer._peak_fp(
        [[100.0004, 1.0], [150.0004, 10.0], [200.0004, 2.0]]
    )
    assert importer._peak_fp(low_high) == importer._peak_fp(list(reversed(low_high)))


def test_precursor_tolerance_is_numeric_not_rounding_bin():
    assert importer._has_precursor_match([500.0], 500.0149)
    assert not importer._has_precursor_match([500.0], 500.0151)
    assert importer._has_precursor_match([499.994, 500.006], 500.0)


def test_polarity_prefers_terminal_adduct_charge_then_declared_mode():
    assert importer._polarity("[M+FA-H]-", 1, "Positive") == "NEG"
    assert importer._polarity("[M+H]+", -1, "Negative") == "POS"
    assert importer._polarity(None, 1, "Negative") == "NEG"
    assert importer._polarity(None, "1-", None) == "NEG"


def test_batch_dedup_uses_spectrum_accession_not_compound_name():
    first = _candidate(
        accession="VF-1",
        peaks=[[100.0, 10.0], [150.0, 100.0], [200.0, 20.0]],
    )
    different_ce = _candidate(
        accession="VF-2",
        peaks=[[100.0, 100.0], [150.0, 10.0], [200.0, 20.0]],
    )
    same_accession = _candidate(
        accession="VF-1",
        peaks=[[101.0, 100.0], [151.0, 10.0], [201.0, 20.0]],
    )
    same_content = _candidate(
        accession="VS19-copy",
        peaks=[[100.0, 100.0], [150.0, 10.0], [200.0, 20.0]],
    )

    selected, _, summary = importer._select_candidates(
        [first, different_ce, same_accession, same_content],
        {},
    )

    assert [row["accession"] for row in selected] == ["vf-1", "vf-2"]
    assert summary["counts"]["duplicate_batch_accession"] == 1
    assert summary["counts"]["duplicate_batch_content"] == 1


def test_database_duplicate_still_reserves_accession_for_later_copy():
    database_copy = _candidate(
        accession="VF-DB",
        peaks=[[100.0, 10.0], [150.0, 100.0], [200.0, 20.0]],
    )
    redistributed = _candidate(
        accession="VF-DB",
        peaks=[[100.0, 11.0], [150.0, 100.0], [200.0, 20.0]],
    )
    db_index = {
        (
            database_copy["identity"],
            "NEG",
            database_copy["peak_fp"],
        ): [database_copy["precursor_mz"]]
    }

    selected, _, summary = importer._select_candidates(
        [database_copy, redistributed],
        db_index,
    )

    assert selected == []
    assert summary["counts"]["duplicate_database_content"] == 1
    assert summary["counts"]["duplicate_batch_accession"] == 1


def test_msp_parser_preserves_db_number_accession(tmp_path, monkeypatch):
    monkeypatch.setattr(
        builder,
        "_resolve_classification_safely",
        lambda **_: "Other",
    )
    source = tmp_path / "test.msp"
    source.write_text(
        "\n".join(
            [
                "NAME: Test",
                "PRECURSORMZ: 300",
                "PRECURSORTYPE: [M-H]-",
                "IONMODE: Negative",
                "INCHIKEY: AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                "SMILES: CCO",
                "Comment: DB#=VF-NPL-QEHF000001; origin=Vaniya/Fiehn",
                "Num Peaks: 3",
                "100 10",
                "150 100",
                "200 20",
                "",
            ]
        ),
        encoding="utf-8",
    )

    records, _ = builder.parse_msp_file(source)

    assert records[0]["source_accession"] == "VF-NPL-QEHF000001"
    assert records[0]["ion_mode"] == "Negative"


def test_mgf_parser_preserves_spectrum_id_accession(tmp_path, monkeypatch):
    monkeypatch.setattr(
        builder,
        "_resolve_classification_safely",
        lambda **_: "Other",
    )
    source = tmp_path / "test.mgf"
    source.write_text(
        "\n".join(
            [
                "BEGIN IONS",
                "PEPMASS=300",
                "CHARGE=1-",
                "IONMODE=negative",
                "NAME=Test",
                "SMILES=CCO",
                "SPECTRUMID=CCMSLIB00000000001",
                "USI=mzspec:GNPS:TEST:accession:CCMSLIB00000000001",
                "100 10",
                "150 100",
                "200 20",
                "END IONS",
                "",
            ]
        ),
        encoding="utf-8",
    )

    records, _ = builder.parse_mgf_file(source)

    assert records[0]["source_accession"] == "CCMSLIB00000000001"
    assert records[0]["ion_mode"] == "negative"
