import json
import sqlite3

from kb_mining.mine_rules_contrastive import MiningConfig, mine_database


def _build_test_db(path):
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        CREATE TABLE Compounds (
            id INTEGER PRIMARY KEY, inchikey TEXT, smiles TEXT,
            classification TEXT, classification_v3 TEXT
        );
        CREATE TABLE Spectra (
            id INTEGER PRIMARY KEY, precursor_mz REAL, peaks_json TEXT,
            adduct TEXT, charge INTEGER, dataset_source TEXT,
            source_file TEXT, source_tier TEXT, collision_energy REAL
        );
        CREATE TABLE Mapping (
            id INTEGER PRIMARY KEY, compound_id INTEGER, spectrum_id INTEGER
        );
        """
    )
    compounds = [
        (1, "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "CCO", "class_a", "class_a"),
        (2, "CCCCCCCCCCCCCC-DDDDDDDDDD-E", "CCN", "class_a", "class_a"),
        (3, "EEEEEEEEEEEEEE-FFFFFFFFFF-G", "CCC", "class_b", "class_b"),
        (4, "GGGGGGGGGGGGGG-HHHHHHHHHH-I", "CCCl", "class_b", "class_b"),
        (5, "IIIIIIIIIIIIII-JJJJJJJJJJ-K", "CCBr", "class_a", "class_a"),
    ]
    connection.executemany("INSERT INTO Compounds VALUES (?,?,?,?,?)", compounds)
    spectrum_id = 1
    mappings = []

    def add(compound_id, source, peaks, tier="experimental"):
        nonlocal spectrum_id
        connection.execute(
            "INSERT INTO Spectra VALUES (?,?,?,?,?,?,?,?,?)",
            (
                spectrum_id, 300.0, json.dumps(peaks), "[M+H]+", 1,
                source, f"{source}.mgf", tier, 20.0,
            ),
        )
        mappings.append((len(mappings) + 1, compound_id, spectrum_id))
        spectrum_id += 1

    # Forty technical/energy replicates still count as one independent compound.
    for _ in range(40):
        add(1, "source_one", [[50.0, 100.0], [100.004, 80.0]])
    add(2, "source_two", [[50.0, 100.0], [100.006, 80.0]])
    add(3, "source_one", [[50.0, 100.0]])
    add(4, "source_two", [[50.0, 100.0]])
    # Predicted spectra must not contribute to empirical support.
    add(5, "source_two", [[50.0, 100.0], [200.0, 90.0]], tier="predicted")
    connection.executemany("INSERT INTO Mapping VALUES (?,?,?)", mappings)
    connection.commit()
    connection.close()


def test_contrastive_miner_counts_structures_not_spectrum_replicates(tmp_path):
    db_path = tmp_path / "library.db"
    _build_test_db(db_path)
    result = mine_database(
        db_path,
        MiningConfig(
            min_independent_compounds=2,
            min_source_families=2,
            min_compounds_per_source=1,
            min_prevalence=0.5,
            min_likelihood_ratio=2.0,
            top_k=10,
        ),
    )

    class_a = result["classes"]["class_a"]["pos"]
    assert class_a["n_independent_compounds"] == 2
    assert class_a["diagnostic_ions"][0]["n_independent_compounds"] == 2
    assert 100.0 <= class_a["diagnostic_ions"][0]["mz"] <= 100.02
    assert not any(abs(rule["mz"] - 50.0) < 0.02 for rule in class_a["diagnostic_ions"])
    assert not any(abs(rule["mz"] - 200.0) < 0.02 for rule in class_a["diagnostic_ions"])
    assert result["meta"]["runtime_enabled"] is False
