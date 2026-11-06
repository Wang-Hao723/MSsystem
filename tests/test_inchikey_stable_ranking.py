import asyncio
import sqlite3
import threading
import time

from matching.compound_annotations import annotate_compound_dict
from matching.library_quality import compound_identity_key, normalize_inchikey
from matching.search import _prescreen_dense_candidates, _select_best_result
from web import server


SINOMENINE_IK = "INYYVPJSBIVGPH-QHRIQVFBSA-N"


def _candidate(**overrides):
    item = {
        "compound_name": "candidate",
        "smiles": "CCO",
        "has_structure": True,
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "formula": "C2H6O",
        "rank_score": 0.8,
        "final_score": 0.8,
        "cosine": 0.9,
        "reverse_match": 0.85,
        "mass_ppm": 2.0,
        "source_tier": "experimental",
        "dataset_source": "MassBank",
        "lib_spectrum_id": "MB000002",
    }
    item.update(overrides)
    return item


def test_full_inchikey_normalization_and_identity_preserve_stereochemistry():
    first = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    second = "AAAAAAAAAAAAAA-DDDDDDDDDD-E"
    assert normalize_inchikey(f" InChIKey={first.lower()} ") == first
    assert normalize_inchikey(first[:14]) == ""
    assert normalize_inchikey("AAAAAAAAAAAAAA-BROKEN") == ""
    assert compound_identity_key({"inchikey": first}) == f"ik::{first}"
    assert compound_identity_key({"inchikey": second}) == f"ik::{second}"
    assert compound_identity_key({"inchikey": first}) != compound_identity_key(
        {"inchikey": second}
    )


def test_top1_is_deterministic_when_scientific_scores_are_identical():
    lower_accession = _candidate(lib_spectrum_id="MB000001")
    higher_accession = _candidate(lib_spectrum_id="MB000002")
    assert _select_best_result([higher_accession, lower_accession])["lib_spectrum_id"] == "MB000001"
    assert _select_best_result([lower_accession, higher_accession])["lib_spectrum_id"] == "MB000001"


def test_dense_prescreen_does_not_merge_distinct_full_inchikeys():
    records = [
        {
            "precursor_mz": 100.0,
            "peaks": [(50.0, 100.0), (70.0, 50.0)],
            "clean_peaks": [(50.0, 100.0), (70.0, 50.0)],
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "source_tier": "experimental",
            "spectrum_id": "A",
        },
        {
            "precursor_mz": 100.0,
            "peaks": [(50.0, 100.0), (70.0, 50.0)],
            "clean_peaks": [(50.0, 100.0), (70.0, 50.0)],
            "inchikey": "AAAAAAAAAAAAAA-DDDDDDDDDD-E",
            "source_tier": "experimental",
            "spectrum_id": "B",
        },
    ]
    selected, _ = _prescreen_dense_candidates(
        [0, 1], records, records[0]["clean_peaks"], 100.0, limit=2
    )
    assert selected == [0, 1]


def test_sinomenine_exact_key_is_curated_and_bilingual():
    item = {"compound_name": "Sinomenine M-H", "inchikey": SINOMENINE_IK}
    annotate_compound_dict(item)
    assert item["compound_name"] == "Sinomenine"
    assert item["compound_name_cn"] == "青藤碱"
    assert item["compound_role"] == "curated_natural"
    assert item["natural_product_candidate"] is True
    assert item["exclude_from_natural_product_count"] is False


def test_server_dedup_uses_full_inchikey_and_stable_representative():
    first = _candidate(
        inchikey="AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        lib_spectrum_id="B",
        final_score=0.7,
        rank_score=0.7,
    )
    better_same = dict(first, lib_spectrum_id="A", final_score=0.8, rank_score=0.8)
    stereo = _candidate(
        inchikey="AAAAAAAAAAAAAA-DDDDDDDDDD-E",
        lib_spectrum_id="C",
        final_score=0.75,
        rank_score=0.75,
    )
    rows = [
        {"spectrum_id": "q1", "match": first},
        {"spectrum_id": "q2", "match": stereo},
        {"spectrum_id": "q3", "match": better_same},
    ]
    deduped = server._dedup_matched(rows)
    assert len(deduped) == 2
    assert deduped[0]["match"]["lib_spectrum_id"] == "A"


def test_library_counts_distinguish_physical_and_loaded_searchable(tmp_path, monkeypatch):
    full = tmp_path / "full.db"
    conn = sqlite3.connect(full)
    conn.executescript(
        """
        CREATE TABLE Spectra (
            id INTEGER PRIMARY KEY, precursor_mz REAL, source_tier TEXT
        );
        CREATE TABLE Compounds (id INTEGER PRIMARY KEY);
        INSERT INTO Compounds VALUES (1), (2);
        INSERT INTO Spectra VALUES
          (1, 100.0, 'experimental'),
          (2, 101.0, 'merged'),
          (3, 102.0, 'predicted'),
          (4, 103.0, 'quarantined');
        """
    )
    conn.commit()
    conn.close()

    monkeypatch.setattr(server, "DB_PATH", full)
    monkeypatch.setattr(server, "DB_PATH_SMALL", tmp_path / "missing.db")
    monkeypatch.setattr(server, "_lib_loaded", True)
    monkeypatch.setattr(server, "_lib_name", "Full (2)")
    monkeypatch.setattr(server, "_lib_spectra_count", 2)

    counts = asyncio.run(server.library_counts())
    assert counts["full"]["physical_spectra"] == 4
    assert counts["full"]["source_searchable_spectra"] == 2
    assert counts["full"]["runtime_searchable_spectra"] == 2
    assert counts["full"]["loaded_searchable_spectra"] == 2
    assert counts["full"]["excluded_by_source"] == 2


def test_library_counts_publish_runtime_structure_gate_before_loading(
    tmp_path, monkeypatch,
):
    full = tmp_path / "full_with_mapping.db"
    conn = sqlite3.connect(full)
    conn.executescript(
        """
        CREATE TABLE Spectra (
            id INTEGER PRIMARY KEY, precursor_mz REAL, source_tier TEXT
        );
        CREATE TABLE Compounds (
            id INTEGER PRIMARY KEY, smiles TEXT
        );
        CREATE TABLE Mapping (
            spectrum_id INTEGER, compound_id INTEGER
        );
        INSERT INTO Compounds VALUES (1, 'CCO'), (2, 'N/A');
        INSERT INTO Spectra VALUES
          (1, 100.0, 'experimental'),
          (2, 101.0, 'experimental');
        INSERT INTO Mapping VALUES (1, 1), (2, 2);
        """
    )
    conn.commit()
    conn.close()

    monkeypatch.setattr(server, "DB_PATH", full)
    monkeypatch.setattr(server, "DB_PATH_SMALL", tmp_path / "missing.db")
    monkeypatch.setattr(server, "_lib_loaded", False)
    monkeypatch.setattr(server, "_lib_name", "")
    monkeypatch.setattr(server, "_lib_spectra_count", 0)

    counts = asyncio.run(server.library_counts())
    assert counts["full"]["source_searchable_spectra"] == 2
    assert counts["full"]["runtime_searchable_spectra"] == 1
    assert counts["full"]["excluded_by_runtime_structure_gate"] == 1
    assert counts["full"]["runtime_projection_exact"] is True
    assert counts["full"]["loaded_searchable_spectra"] is None


def test_library_load_does_not_block_unrelated_async_work(tmp_path, monkeypatch):
    full = tmp_path / "full.db"
    full.touch()
    loader_started = threading.Event()
    loader_release = threading.Event()

    def slow_loader(_path):
        loader_started.set()
        assert loader_release.wait(timeout=2.0)
        return ([{"precursor_mz": 100.0}], {1000: [0]})

    monkeypatch.setattr(server, "DB_PATH", full)
    monkeypatch.setattr(server, "_tasks", {})
    monkeypatch.setattr(server, "_lib_records", [])
    monkeypatch.setattr(server, "_mz_index", {})
    monkeypatch.setattr(server, "_lib_loaded", False)
    monkeypatch.setattr(server, "_lib_name", "")
    monkeypatch.setattr(server, "_lib_spectra_count", 0)
    monkeypatch.setattr(server, "load_library", slow_loader)

    async def scenario():
        started_at = time.monotonic()
        task = asyncio.create_task(server.library_load("full"))
        for _ in range(100):
            if loader_started.is_set():
                break
            await asyncio.sleep(0.005)
        assert loader_started.is_set()
        await asyncio.sleep(0.02)
        # A synchronous loader would freeze the loop until the two-second
        # worker timeout. The ticker must stay responsive during the load.
        assert time.monotonic() - started_at < 0.5
        assert not task.done()
        loader_release.set()
        status = await task
        assert status.loaded is True
        assert status.spectra_count == 1

    try:
        asyncio.run(scenario())
    finally:
        loader_release.set()
