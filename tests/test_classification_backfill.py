import sqlite3

from scripts.backfill_classification_v3 import (
    apply_updates,
    classify_missing_rows,
    coverage,
    rollback_ledger,
    write_ledger,
)


def _database(tmp_path):
    path = tmp_path / "library.db"
    connection = sqlite3.connect(path)
    connection.execute(
        """
        CREATE TABLE Compounds (
            id INTEGER PRIMARY KEY,
            compound_name TEXT,
            smiles TEXT,
            classification TEXT,
            classification_v3 TEXT,
            classification_source TEXT,
            classification_confidence REAL
        )
        """
    )
    connection.executemany(
        "INSERT INTO Compounds VALUES (?,?,?,?,?,?,?)",
        [
            (1, "Piperidine", "C1CCNCC1", "alkaloid", None, None, None),
            (
                2,
                "Adenosine",
                "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
                "other",
                None,
                None,
                None,
            ),
            (3, "Benzene", "c1ccccc1", "other", None, None, None),
            (4, "Placeholder", "N/A", "other", None, None, None),
            (
                5,
                "Protected",
                "C1CCNCC1",
                "alkaloid",
                "curated_class",
                "curated",
                1.0,
            ),
        ],
    )
    connection.commit()
    connection.row_factory = sqlite3.Row
    return path, connection


def test_backfill_requires_structure_evidence_and_never_overwrites(tmp_path):
    _, connection = _database(tmp_path)

    updates, counters = classify_missing_rows(connection)

    by_id = {update.compound_id: update for update in updates}
    assert set(by_id) == {1, 2}
    assert by_id[1].new_classification_v3 == "alkaloid"
    assert by_id[1].new_classification_source == "smarts_raw_agreement_v6"
    assert by_id[2].new_classification_v3 == "nucleoside"
    assert by_id[2].new_classification_source == "smarts_high_specific_v6"
    assert counters["skipped_no_smarts_class"] == 1
    assert counters["skipped_no_usable_structure"] == 1
    assert counters["eligible_updates"] == 2

    assert apply_updates(connection, updates) == 2
    protected = connection.execute(
        "SELECT classification_v3, classification_source FROM Compounds WHERE id=5"
    ).fetchone()
    assert tuple(protected) == ("curated_class", "curated")
    assert coverage(connection)["with_classification_v3"] == 3
    connection.close()


def test_backfill_ledger_is_reversible_without_clobbering_later_edits(tmp_path):
    db_path, connection = _database(tmp_path)
    updates, _ = classify_missing_rows(connection)
    ledger = tmp_path / "classification.jsonl.gz"
    write_ledger(ledger, db_path, updates)
    apply_updates(connection, updates)

    # A later curator edit must survive rollback because it no longer carries
    # this migration's exact class/source pair.
    connection.execute(
        """
        UPDATE Compounds
        SET classification_v3='curated_nitrogen_heterocycle',
            classification_source='curated'
        WHERE id=1
        """
    )
    connection.commit()

    assert rollback_ledger(connection, ledger) == 1
    rows = {
        row["id"]: (row["classification_v3"], row["classification_source"])
        for row in connection.execute(
            "SELECT id, classification_v3, classification_source FROM Compounds"
        )
    }
    assert rows[1] == ("curated_nitrogen_heterocycle", "curated")
    assert rows[2] == (None, None)
    connection.close()
