from __future__ import annotations

import json

from scripts.build_v11_identity_registry import (
    collect_conflicts,
    validate_resolution,
)


def test_collects_reported_inchikey_vs_structure_conflict(tmp_path):
    artifact = tmp_path / "synthetic_v11_flagship.json"
    artifact.write_text(
        json.dumps(
            {
                "entities": [
                    {
                        "candidates": [
                            {
                                "identity_key": "std2d:AAAAAAAAAAAAAA",
                                "candidate": {
                                    "compound_name": "Caffeine",
                                    "formula": "C8H10N4O2",
                                    "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
                                    "standardized_connectivity_inchikey": "AAAAAAAAAAAAAA",
                                    "reported_connectivity_matches_standardized": False,
                                },
                            }
                        ]
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    conflicts = collect_conflicts([artifact])
    conflict = conflicts[("caffeine", "C8H10N4O2")]

    assert conflict["triggers"] == {
        "reported_inchikey_disagrees_with_stored_structure"
    }
    assert conflict["observed_connectivity_identities"] == {
        "std2d:AAAAAAAAAAAAAA",
        "ik14:RYYVLZVUVIJVGH",
    }
    assert conflict["source_artifacts"] == {artifact.name}


def test_collects_conflict_directly_from_immutable_result_rows(tmp_path):
    artifact = tmp_path / "synthetic_results.json"
    artifact.write_text(
        json.dumps(
            [
                {
                    "spectrum_id": "scan-1",
                    "match": {
                        "compound_name": "Caffeine",
                        "formula": "C8H10N4O2",
                        "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
                        # Ethanol deliberately conflicts with the reported key.
                        "smiles": "CCO",
                        "top_candidates": [],
                    },
                }
            ]
        ),
        encoding="utf-8",
    )

    conflicts = collect_conflicts([artifact])
    conflict = conflicts[("caffeine", "C8H10N4O2")]

    assert conflict["triggers"] == {
        "reported_inchikey_disagrees_with_stored_structure"
    }
    assert "ik14:RYYVLZVUVIJVGH" in conflict["observed_connectivity_identities"]
    assert any(
        token.startswith("std2d:")
        for token in conflict["observed_connectivity_identities"]
    )


def test_pubchem_resolution_requires_formula_and_rdkit_inchikey_agreement():
    conflict = {
        "name_key": "caffeine",
        "query_name": "Caffeine",
        "formula": "C8H10N4O2",
        "observed_connectivity_identities": {
            "std2d:AAAAAAAAAAAAAA",
            "ik14:RYYVLZVUVIJVGH",
        },
        "source_artifacts": {"synthetic_v11_flagship.json"},
        "triggers": {"reported_inchikey_disagrees_with_stored_structure"},
    }
    resolved, reason = validate_resolution(
        conflict,
        {
            "CID": 2519,
            "Title": "Caffeine",
            "MolecularFormula": "C8H10N4O2",
            "SMILES": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "ConnectivitySMILES": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "InChI": "",
            "InChIKey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
            "property_url": "https://pubchem.ncbi.nlm.nih.gov/compound/2519",
        },
    )

    assert reason == "accepted"
    assert resolved is not None
    assert resolved["connectivity_inchikey"] == "RYYVLZVUVIJVGH"
    assert resolved["resolution_triggers"] == [
        "reported_inchikey_disagrees_with_stored_structure"
    ]
