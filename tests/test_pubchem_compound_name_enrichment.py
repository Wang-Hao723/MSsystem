import sqlite3

from scripts.enrich_compound_name_aliases_pubchem import (
    _safe_pubchem_name,
    build_alias_payload,
    collect_targets,
)


def _make_compound_db(path):
    with sqlite3.connect(path) as connection:
        connection.execute(
            """
            CREATE TABLE Compounds (
                id INTEGER PRIMARY KEY,
                compound_name TEXT,
                inchikey TEXT,
                structure_identity_verified INTEGER
            )
            """
        )
        connection.executemany(
            """
            INSERT INTO Compounds(compound_name, inchikey, structure_identity_verified)
            VALUES (?, ?, ?)
            """,
            [
                (
                    None,
                    "QECYIHDHCKWWPU-CMDGGOBGSA-N",
                    1,
                ),
                (
                    "Cyclohomonervilasterol",
                    "QECYIHDHCKWWPU-CMDGGOBGSA-N",
                    1,
                ),
                (
                    "NCGC00000001",
                    "BBBBBBBBBBBBBB-CCCCCCCCCC-D",
                    0,
                ),
            ],
        )


def test_collect_targets_requires_verified_structure_and_keeps_local_fallback(
    tmp_path,
):
    database = tmp_path / "names.db"
    _make_compound_db(database)
    targets, local_candidates, status_counts = collect_targets(database, {})

    assert targets == {"QECYIHDHCKWWPU-CMDGGOBGSA-N"}
    assert (
        local_candidates["QECYIHDHCKWWPU-CMDGGOBGSA-N"]
        == "Cyclohomonervilasterol"
    )
    assert status_counts["missing"] == 1


def test_pubchem_name_requires_exact_full_key_and_safe_title():
    key = "QECYIHDHCKWWPU-CMDGGOBGSA-N"
    assert _safe_pubchem_name(
        key,
        {
            "InChIKey": key,
            "Title": "Cyclohomonervilasterol",
            "IUPACName": "long systematic name",
        },
    ) == ("Cyclohomonervilasterol", "PubChem_Title")
    assert _safe_pubchem_name(
        key,
        {
            "InChIKey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "Title": "Wrong structure",
        },
    ) is None


def test_payload_prefers_pubchem_then_uses_exact_key_local_fallback(tmp_path):
    pubchem_key = "QECYIHDHCKWWPU-CMDGGOBGSA-N"
    local_key = "BBBBBBBBBBBBBB-CCCCCCCCCC-D"
    payload = build_alias_payload(
        target_keys={pubchem_key, local_key},
        local_candidates={
            pubchem_key: "Local synonym",
            local_key: "Verified local name",
        },
        properties_by_key={
            pubchem_key: {
                "InChIKey": pubchem_key,
                "Title": "Cyclohomonervilasterol",
            },
            local_key: None,
        },
        database=tmp_path / "library.db",
        status_counts={"missing": 2},
    )

    assert payload["aliases"][pubchem_key] == "Cyclohomonervilasterol"
    assert payload["provenance"][pubchem_key] == "PubChem_Title"
    assert payload["aliases"][local_key] == "Verified local name"
    assert (
        payload["provenance"][local_key]
        == "local_exact_full_inchikey_synonym"
    )
