from __future__ import annotations

import sqlite3

from matching.formula_candidate_shadow import (
    build_formula_conditioned_library_shadow,
)


def _database(path):
    connection = sqlite3.connect(path)
    connection.execute("""
        CREATE TABLE Compounds (
            id INTEGER PRIMARY KEY,
            compound_name TEXT,
            molecular_formula TEXT,
            smiles TEXT,
            inchikey TEXT,
            exact_mass REAL,
            classification TEXT,
            reported_inchikey TEXT,
            structure_identity_verified INTEGER
        )
    """)
    connection.executemany(
        "INSERT INTO Compounds VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (
                1,
                "Candidate A",
                "C2H6O",
                "CCO",
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                46.0419,
                "other",
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                1,
            ),
            (
                2,
                "Duplicate A",
                "C2H6O",
                "CCO",
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                46.0419,
                "other",
                "",
                0,
            ),
            (
                3,
                "Candidate B",
                "C2H6O",
                "COC",
                "LCGLNKUTAGEVQW-UHFFFAOYSA-N",
                46.0419,
                "other",
                "",
                0,
            ),
        ],
    )
    connection.commit()
    connection.close()


def test_formula_conditioned_catalog_is_shadow_and_deduplicated(tmp_path):
    database = tmp_path / "library.db"
    _database(database)
    formula_shadow = {
        "task_id": "task",
        "status": "completed_shadow",
        "spectra": [
            {
                "row_index": 7,
                "spectrum_id": "scan=7",
                "matched": False,
                "formula_candidates": [
                    {"formula": "C2H6O", "score": 0.8},
                    {"formula": "C3H8O", "score": 0.1},
                ],
            }
        ],
    }

    payload = build_formula_conditioned_library_shadow(
        formula_shadow,
        db_path=database,
        max_candidates_per_formula=1,
    )

    assert payload["status"] == "completed_shadow"
    assert payload["runtime_influence"] is False
    assert payload["ranking_changed"] is False
    assert payload["summary"]["unmatched_spectra_with_structure_hypotheses"] == 1
    assert payload["summary"]["high_promotions"] == 0
    assert payload["formula_catalog"]["C2H6O"]["unique_structure_candidates"] == 2
    assert payload["formula_catalog"]["C2H6O"]["truncated"] is True
    assert len(payload["formula_catalog"]["C2H6O"]["returned_candidates"]) == 1


def test_formula_conditioned_catalog_fails_closed_without_completed_formula_shadow(
    tmp_path,
):
    payload = build_formula_conditioned_library_shadow(
        {"task_id": "task", "status": "failed", "spectra": []},
        db_path=tmp_path / "missing.db",
    )

    assert payload["status"] == "unavailable"
    assert payload["summary"]["reason"] == "formula_shadow_not_completed"
    assert payload["high_promotions"] == 0
