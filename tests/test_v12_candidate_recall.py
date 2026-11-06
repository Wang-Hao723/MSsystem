from __future__ import annotations

import hashlib
import json
import sqlite3
from copy import deepcopy
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import inchi, rdMolDescriptors

from matching.evirank_retrieval import (
    SQLiteCandidateRecallRetriever,
    _reference_content_fingerprint,
    _source_paired_consensus,
    build_candidate_recall_bundle,
    load_evirecall_policy,
)
from pipeline.evirank_retrieval import (
    _sampled_file_sha256,
    read_evirecall_artifact,
    write_evirecall_artifact,
)


QUERY_PEAKS = [[20.0, 100.0], [30.0, 80.0], [40.0, 60.0], [50.0, 40.0]]
DISTRACTOR_PEAKS = [[21.0, 100.0], [31.0, 80.0], [41.0, 60.0], [51.0, 40.0]]


def _structure(smiles: str) -> tuple[str, str]:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return (
        rdMolDescriptors.CalcMolFormula(molecule),
        inchi.MolToInchiKey(molecule),
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _build_library(path: Path) -> dict[str, str]:
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        CREATE TABLE Compounds (
            id INTEGER PRIMARY KEY,
            compound_name TEXT,
            molecular_formula TEXT,
            smiles TEXT,
            inchikey TEXT,
            exact_mass REAL
        );
        CREATE TABLE Spectra (
            id INTEGER PRIMARY KEY,
            spectrum_id TEXT,
            precursor_mz REAL,
            peaks_json TEXT,
            adduct TEXT,
            charge INTEGER,
            ion_mode TEXT,
            source_tier TEXT,
            dataset_source TEXT,
            quality_flags_json TEXT
        );
        CREATE TABLE Mapping (
            id INTEGER PRIMARY KEY,
            compound_id INTEGER,
            spectrum_id INTEGER,
            score REAL
        );
        """
    )
    definitions = [
        # name, structure, precursor, peaks, adduct, ion mode, tier, flags
        ("correct", "CCO", 100.005, QUERY_PEAKS, "[M+H]+", "POS", "experimental", []),
        ("distractor", "CCN", 99.997, DISTRACTOR_PEAKS, "[M+H]+", "POS", "experimental", []),
        ("outside-window", "CC(=O)O", 100.0301, QUERY_PEAKS, "[M+H]+", "POS", "experimental", []),
        ("wrong-polarity", "c1ccccc1", 100.004, QUERY_PEAKS, "[M-H]-", "NEG", "experimental", []),
        ("predicted", "CCCO", 100.003, QUERY_PEAKS, "[M+H]+", "POS", "predicted", []),
        (
            "quarantined-spectrum",
            "CCCC",
            100.002,
            QUERY_PEAKS,
            "[M+H]+",
            "POS",
            "experimental",
            ["excluded_from_ordinary_spectral_search"],
        ),
    ]
    identities: dict[str, str] = {}
    for index, (name, smiles, precursor, peaks, adduct, mode, tier, flags) in enumerate(
        definitions, 1
    ):
        formula, inchikey = _structure(smiles)
        identities[name] = inchikey.split("-", 1)[0]
        connection.execute(
            "INSERT INTO Compounds VALUES (?,?,?,?,?,?)",
            (index, name, formula, smiles, inchikey, 0.0),
        )
        connection.execute(
            "INSERT INTO Spectra VALUES (?,?,?,?,?,?,?,?,?,?)",
            (
                index,
                f"spectrum-{index}",
                precursor,
                json.dumps(peaks),
                adduct,
                1 if mode == "POS" else -1,
                mode,
                tier,
                "synthetic-test",
                json.dumps(flags),
            ),
        )
        connection.execute(
            "INSERT INTO Mapping VALUES (?,?,?,?)",
            (index, index, index, 1.0),
        )

    # A name and formula are not a stable structure identity and must not enter
    # EviRecall even when its spectrum is an exact match.
    index = len(definitions) + 1
    connection.execute(
        "INSERT INTO Compounds VALUES (?,?,?,?,?,?)",
        (index, "name-only", "C5H10O", "", "", 0.0),
    )
    connection.execute(
        "INSERT INTO Spectra VALUES (?,?,?,?,?,?,?,?,?,?)",
        (
            index,
            f"spectrum-{index}",
            100.001,
            json.dumps(QUERY_PEAKS),
            "[M+H]+",
            1,
            "POS",
            "experimental",
            "synthetic-test",
            "[]",
        ),
    )
    connection.execute(
        "INSERT INTO Mapping VALUES (?,?,?,?)",
        (index, index, index, 1.0),
    )
    connection.commit()
    connection.close()
    return identities


def _build_source_paired_library(
    path: Path,
    *,
    sources: tuple[str, ...] = ("Source-A", "Source-B"),
) -> dict[str, str]:
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        CREATE TABLE Compounds (
            id INTEGER PRIMARY KEY,
            compound_name TEXT,
            molecular_formula TEXT,
            smiles TEXT,
            inchikey TEXT,
            exact_mass REAL
        );
        CREATE TABLE Spectra (
            id INTEGER PRIMARY KEY,
            spectrum_id TEXT,
            precursor_mz REAL,
            peaks_json TEXT,
            adduct TEXT,
            charge INTEGER,
            ion_mode TEXT,
            source_tier TEXT,
            dataset_source TEXT,
            quality_flags_json TEXT
        );
        CREATE TABLE Mapping (
            id INTEGER PRIMARY KEY,
            compound_id INTEGER,
            spectrum_id INTEGER,
            score REAL
        );
        """
    )
    formula, baseline_full_key = _structure("CCO")
    challenger_formula, challenger_full_key = _structure("COC")
    assert challenger_formula == formula
    identities = {
        "baseline": baseline_full_key.split("-", 1)[0],
        "baseline_full": baseline_full_key,
        "challenger": challenger_full_key.split("-", 1)[0],
        "challenger_full": challenger_full_key,
        "formula": formula,
    }
    connection.execute(
        "INSERT INTO Compounds VALUES (?,?,?,?,?,?)",
        (1, "baseline-name", formula, "CCO", baseline_full_key, 0.0),
    )
    connection.execute(
        "INSERT INTO Compounds VALUES (?,?,?,?,?,?)",
        (2, "challenger-name", formula, "COC", challenger_full_key, 0.0),
    )
    spectrum_id = 1
    for source in sources:
        repeat_count = 3 if source == "Source-A" else 1
        mirrored_source_a = source in {"Source-A", "Source-A Mirror"}
        baseline_peaks = (
            [[20.2, 100.0], [30.2, 80.0], [40.2, 60.0], [50.2, 40.0]]
            if mirrored_source_a
            else [[20.3, 100.0], [30.3, 75.0], [40.3, 65.0], [50.3, 35.0]]
        )
        challenger_peaks = (
            QUERY_PEAKS
            if mirrored_source_a
            else [[20.0, 100.0], [30.0, 78.0], [40.0, 62.0], [50.0, 42.0]]
        )
        for compound_id, peaks in (
            (1, baseline_peaks),
            (2, challenger_peaks),
        ):
            for repeat in range(repeat_count):
                connection.execute(
                    "INSERT INTO Spectra VALUES (?,?,?,?,?,?,?,?,?,?)",
                    (
                        spectrum_id,
                        f"paired-{source}-{compound_id}-{repeat}",
                        100.0,
                        json.dumps(peaks),
                        "[M+H]+",
                        1,
                        "POS",
                        "experimental",
                        source,
                        "[]",
                    ),
                )
                connection.execute(
                    "INSERT INTO Mapping VALUES (?,?,?,?)",
                    (spectrum_id, compound_id, spectrum_id, 1.0),
                )
                spectrum_id += 1
    connection.commit()
    connection.close()
    return identities


@pytest.fixture()
def synthetic_library(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    path = tmp_path / "evirecall.sqlite"
    return path, _build_library(path)


@pytest.fixture()
def source_paired_library(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    path = tmp_path / "source-paired.sqlite"
    return path, _build_source_paired_library(path)


def _query(query_id: str = "query-1", collision_energy: float = 10.0) -> dict:
    return {
        "spectrum_id": query_id,
        "precursor_mz": 100.0,
        "polarity": "POS",
        "charge": 1,
        "collision_energy": collision_energy,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "absolute",
        "peaks": QUERY_PEAKS,
    }


def _retriever(path: Path, *, policy: dict | None = None) -> SQLiteCandidateRecallRetriever:
    return SQLiteCandidateRecallRetriever(
        path,
        library_snapshot_id="synthetic-library-sha256",
        policy=policy,
    )


def _source_paired_existing(
    identities: dict[str, str],
    *,
    baseline_score: float = 0.55,
    runner_score: float = 0.50,
) -> list[dict]:
    runner_formula, runner_full_key = _structure("CCN")
    return [
        {
            "identity_key": f"ik14:{identities['baseline']}",
            "entity_score": baseline_score,
            "candidate": {
                "formula": identities["formula"],
                "inchikey": identities["baseline_full"],
            },
        },
        {
            "identity_key": f"ik14:{runner_full_key.split('-', 1)[0]}",
            "entity_score": runner_score,
            "candidate": {
                "formula": runner_formula,
                "inchikey": runner_full_key,
            },
        },
    ]


def test_exact_bounded_retrieval_applies_structure_polarity_quality_and_source_gates(
    synthetic_library,
):
    path, identities = synthetic_library
    before = _sha256(path)

    result = _retriever(path).retrieve_scan(_query())

    assert result["status"] == "available"
    assert result["score_is_probability"] is False
    assert result["window_audit"]["effective_tolerance_da"] == pytest.approx(0.02)
    assert result["window_audit"]["blind_window_widening"] is False
    assert result["window_audit"]["blind_adduct_enumeration"] is False
    assert result["window_audit"]["excluded_wrong_polarity"] == 1
    assert result["window_audit"]["excluded_quality"] == 1
    assert result["window_audit"]["excluded_unstable_identity"] == 1
    returned = set(result["candidates"])
    assert returned == {identities["correct"], identities["distractor"]}
    assert identities["outside-window"] not in returned
    assert identities["predicted"] not in returned
    assert result["rankings"]["library_match"][0] == identities["correct"]
    assert result["rankings"]["sensitive"][0] == identities["correct"]
    assert _sha256(path) == before


def test_ground_truth_and_expected_identity_are_rejected_even_when_nested(
    synthetic_library,
):
    path, _ = synthetic_library
    retriever = _retriever(path)
    with pytest.raises(ValueError, match="metadata.ground_truth"):
        retriever.retrieve_scan({**_query(), "metadata": {"ground_truth": "correct"}})
    with pytest.raises(ValueError, match="compound_name"):
        retriever.retrieve_scan({**_query(), "compound_name": "correct"})


def test_custom_policy_cannot_enable_blind_widening_or_production_ranking(
    synthetic_library,
):
    path, _ = synthetic_library
    policy = load_evirecall_policy()
    widened = deepcopy(policy)
    widened["mass_retrieval"]["blind_window_widening_allowed"] = True
    with pytest.raises(ValueError, match="widening"):
        _retriever(path, policy=widened)

    production = deepcopy(policy)
    production["invariants"]["production_ranking_changes"] = True
    with pytest.raises(ValueError, match="production ranking"):
        _retriever(path, policy=production)

    uncollapsed = deepcopy(policy)
    uncollapsed["source_paired_consensus"]["mirrored_provenance_strata"] = (
        "count_labels_as_votes"
    )
    with pytest.raises(ValueError, match="content-collapsed"):
        _retriever(path, policy=uncollapsed)


def test_neutral_loss_view_is_symmetric_and_remains_shadow(synthetic_library):
    path, identities = synthetic_library
    policy = load_evirecall_policy()
    policy["spectral_views"].append(
        {
            "name": "neutral_loss",
            "peak_processing_preset": "library_match",
            "similarity": "matchms_neutral_losses_cosine",
            "fragment_tolerance_da": 0.02,
        }
    )
    result = _retriever(path, policy=policy).retrieve_scan(_query())
    assert result["rankings"]["neutral_loss"][0] == identities["correct"]
    assert "neutral_loss" in result["candidates"][identities["correct"]]["query_peak_audits"]
    entity = _retriever(path, policy=policy).retrieve_entity("e", [_query()])
    assert entity["production_effect"] is False


def test_reference_content_fingerprint_is_name_free_order_and_scale_invariant():
    first = {
        "compound_name": "first-name",
        "dataset_source": "source-a",
        "peaks": [[20.0, 10.0], [30.0, 5.0], [40.0, 1.0]],
    }
    mirrored = {
        "compound_name": "truth-looking-renamed-copy",
        "dataset_source": "source-b",
        "peaks": [[40.0, 10.0], [20.0, 100.0], [30.0, 50.0]],
    }
    distinct = {"peaks": [[20.0, 10.0], [30.1, 5.0], [40.0, 1.0]]}

    assert _reference_content_fingerprint(first) == _reference_content_fingerprint(
        mirrored
    )
    assert _reference_content_fingerprint(first) != _reference_content_fingerprint(
        distinct
    )


def test_dense_precursor_window_exceeds_budget_and_fails_closed(synthetic_library):
    path, _ = synthetic_library
    policy = load_evirecall_policy()
    policy["library_gate"]["maximum_sql_rows_per_precursor_window"] = 1

    result = _retriever(path, policy=policy).retrieve_scan(_query())

    assert result["status"] == "dense_window_budget_exceeded_fail_closed"
    assert result["candidates"] == {}
    assert result["window_audit"]["sql_rows"] > 1


def test_read_only_retrieval_uses_verified_primary_key_mapping_join(
    synthetic_library,
):
    path, identities = synthetic_library

    result = _retriever(path).retrieve_scan(_query())

    assert result["rankings"]["library_match"][0] == identities["correct"]
    strategy = result["window_audit"]["sqlite_query_strategy"]
    assert strategy["mapping_join_strategy"] == (
        "verified_mapping_primary_key_alignment"
    )
    assert strategy["library_read_only"] is True


def test_non_aligned_mapping_without_spectrum_index_fails_closed(
    synthetic_library,
):
    path, _ = synthetic_library
    connection = sqlite3.connect(path)
    connection.execute("UPDATE Mapping SET id=id+100")
    connection.commit()
    connection.close()

    with pytest.raises(RuntimeError, match="Mapping.spectrum_id"):
        _retriever(path)


def test_non_aligned_mapping_uses_explicit_spectrum_index(
    synthetic_library,
):
    path, identities = synthetic_library
    connection = sqlite3.connect(path)
    connection.execute("UPDATE Mapping SET id=id+100")
    connection.execute(
        "CREATE INDEX ix_mapping_spectrum_id ON Mapping(spectrum_id)"
    )
    connection.execute(
        "CREATE INDEX ix_spectra_precursor_mz ON Spectra(precursor_mz)"
    )
    connection.commit()
    connection.close()

    result = _retriever(path).retrieve_scan(_query())

    assert result["rankings"]["library_match"][0] == identities["correct"]
    strategy = result["window_audit"]["sqlite_query_strategy"]
    assert strategy["mapping_join_strategy"] == "indexed_mapping_spectrum_id"
    assert strategy["mapping_spectrum_index"] == "ix_mapping_spectrum_id"
    assert strategy["precursor_index"] == "ix_spectra_precursor_mz"


def test_names_are_display_only_and_do_not_change_identity_ranking(synthetic_library):
    path, _ = synthetic_library
    first = _retriever(path).retrieve_entity(
        "entity-1", [_query()], detection_method_id="lc-ms-method-a"
    )
    first_signature = [
        (
            item["connectivity_identity"],
            item["retrieval_rank"],
            item["retrieval_rrf_utility"],
        )
        for item in first["expanded_candidates"]
    ]

    connection = sqlite3.connect(path)
    connection.execute("UPDATE Compounds SET compound_name='renamed-' || id")
    connection.commit()
    connection.close()
    second = _retriever(path).retrieve_entity(
        "entity-1", [_query()], detection_method_id="lc-ms-method-a"
    )
    second_signature = [
        (
            item["connectivity_identity"],
            item["retrieval_rank"],
            item["retrieval_rrf_utility"],
        )
        for item in second["expanded_candidates"]
    ]

    assert first_signature == second_signature
    assert second["candidate_names_affect_retrieval_or_scoring"] is False


def test_technical_repeats_are_collapsed_but_distinct_collision_energies_sum(
    synthetic_library,
):
    path, identities = synthetic_library
    scans = [
        _query("repeat-a", 10.0),
        _query("repeat-b", 10.0),
        _query("condition-b", 20.0),
    ]

    result = _retriever(path).retrieve_entity(
        "entity-1", scans, detection_method_id="lc-ms-method-a"
    )

    assert result["scan_count"] == 3
    assert result["distinct_condition_count"] == 2
    assert result["technical_repeats_increase_condition_count"] is False
    correct = next(
        item
        for item in result["expanded_candidates"]
        if item["connectivity_identity"] == identities["correct"]
    )
    assert correct["distinct_condition_support"] == 2
    repeats = sorted(
        item["technical_repeat_count"]
        for item in correct["condition_view_contributions"]
    )
    assert repeats == [1, 1, 2, 2]
    assert correct["score_is_probability"] is False
    assert correct["production_effect"] is False


def test_source_paired_consensus_qualifies_only_after_cross_stratum_multice_support(
    source_paired_library,
):
    path, identities = source_paired_library
    result = _retriever(path).retrieve_entity(
        "paired-entity",
        [_query("ce-10", 10.0), _query("ce-20", 20.0)],
        existing_candidates=_source_paired_existing(identities),
        detection_method_id="lc-ms-method-a",
    )

    consensus = result["source_paired_consensus"]
    assert consensus["status"] == "qualified_switch"
    assert consensus["selected_identity_key"] == f"std2d:{identities['challenger']}"
    assert consensus["ground_truth_used_at_inference"] is False
    assert consensus["candidate_names_affect_scoring"] is False
    comparison = next(item for item in consensus["comparisons"] if item["qualified"])
    assert comparison["common_provenance_stratum_count"] == 2
    assert comparison["unique_provenance_evidence_profile_count"] == 2
    assert comparison["unique_provenance_evidence_profile_sizes"] == [1, 1]
    assert comparison["distinct_collision_energy_condition_count"] == 2
    assert comparison["spectral_view_count"] == 2
    assert comparison["paired_observation_count"] == 8
    assert comparison["observation_win_fraction"] == 1.0
    assert comparison["provenance_stratum_win_fraction"] == 1.0
    assert comparison["unique_provenance_profile_win_fraction"] == 1.0
    challenger = next(
        item
        for item in result["expanded_candidates"]
        if item["connectivity_identity"] == identities["challenger"]
    )
    assert challenger["source_paired_selection_qualified"] is True
    scan = _retriever(path).retrieve_scan(_query())
    assert scan["candidates"][identities["challenger"]][
        "source_reference_record_counts"
    ]["library_match"]["source-a"] == 3
    assert comparison["common_provenance_stratum_count"] == 2


def test_source_paired_consensus_fails_closed_with_one_stratum_or_one_ce(
    source_paired_library,
    tmp_path: Path,
):
    path, identities = source_paired_library
    one_ce = _retriever(path).retrieve_entity(
        "one-ce",
        [_query("repeat-a", 10.0), _query("repeat-b", 10.0)],
        existing_candidates=_source_paired_existing(identities),
    )
    one_ce_comparison = one_ce["source_paired_consensus"]["comparisons"][0]
    assert one_ce["source_paired_consensus"]["status"] == "baseline_retained"
    assert "insufficient_distinct_collision_energy_conditions" in one_ce_comparison[
        "failed_gates"
    ]
    assert one_ce_comparison["paired_observation_count"] == 4

    one_source_path = tmp_path / "one-source.sqlite"
    one_source_identities = _build_source_paired_library(
        one_source_path, sources=("Source-A",)
    )
    one_source = _retriever(one_source_path).retrieve_entity(
        "one-source",
        [_query("ce-10", 10.0), _query("ce-20", 20.0)],
        existing_candidates=_source_paired_existing(one_source_identities),
    )
    one_source_comparison = one_source["source_paired_consensus"]["comparisons"][0]
    assert one_source["source_paired_consensus"]["status"] == "baseline_retained"
    assert "insufficient_common_provenance_strata" in one_source_comparison[
        "failed_gates"
    ]


def test_mirrored_source_labels_collapse_to_one_evidence_profile_and_abstain(
    tmp_path: Path,
):
    path = tmp_path / "mirrored-sources.sqlite"
    identities = _build_source_paired_library(
        path, sources=("Source-A", "Source-A Mirror")
    )

    result = _retriever(path).retrieve_entity(
        "mirrored-sources",
        [_query("ce-10", 10.0), _query("ce-20", 20.0)],
        existing_candidates=_source_paired_existing(identities),
    )

    comparison = result["source_paired_consensus"]["comparisons"][0]
    assert result["source_paired_consensus"]["status"] == "baseline_retained"
    assert comparison["common_provenance_stratum_count"] == 2
    assert comparison["unique_provenance_evidence_profile_count"] == 1
    assert comparison["unique_provenance_evidence_profile_sizes"] == [2]
    assert "insufficient_unique_provenance_evidence_profiles" in comparison[
        "failed_gates"
    ]


def test_source_paired_consensus_respects_v11_margin_gate_and_ignores_names(
    source_paired_library,
):
    path, identities = source_paired_library
    scans = [_query("ce-10", 10.0), _query("ce-20", 20.0)]
    high_margin = _retriever(path).retrieve_entity(
        "high-margin",
        scans,
        existing_candidates=_source_paired_existing(
            identities, baseline_score=0.70, runner_score=0.50
        ),
    )
    assert high_margin["source_paired_consensus"]["status"] == "baseline_retained"
    assert high_margin["source_paired_consensus"]["reason"] == "v11_margin_above_gate"

    before = _retriever(path).retrieve_entity(
        "name-blind",
        scans,
        existing_candidates=_source_paired_existing(identities),
    )["source_paired_consensus"]
    connection = sqlite3.connect(path)
    connection.execute("UPDATE Compounds SET compound_name='truth-looking-' || id")
    connection.commit()
    connection.close()
    after = _retriever(path).retrieve_entity(
        "name-blind",
        scans,
        existing_candidates=_source_paired_existing(identities),
    )["source_paired_consensus"]
    assert before == after


def test_native_same_formula_shortlist_is_not_blocked_by_large_retrieval_universe():
    baseline = "AAAAAAAAAAAAAA"
    challenger = "BBBBBBBBBBBBBB"
    expanded = [f"X{index:013d}" for index in range(20)]
    existing = [
        {
            "identity_key": f"std2d:{baseline}",
            "entity_score": 0.55,
            "candidate": {"formula": "C2H6O"},
        },
        {
            "identity_key": f"std2d:{challenger}",
            "entity_score": 0.50,
            "candidate": {"formula": "C2H6O"},
        },
    ]
    scan_results = []
    for condition in ("ce:10|ev|absolute", "ce:20|ev|absolute"):
        candidates = {
            baseline: {
                "candidate": {"formula": "C2H6O"},
                "source_view_scores": {
                    "library_match": {"source-a": 0.2, "source-b": 0.2},
                    "sensitive": {"source-a": 0.2, "source-b": 0.2},
                },
                "source_reference_content_fingerprints": {
                    "library_match": {
                        "source-a": "baseline-a",
                        "source-b": "baseline-b",
                    },
                    "sensitive": {
                        "source-a": "baseline-a",
                        "source-b": "baseline-b",
                    },
                },
            },
            challenger: {
                "candidate": {"formula": "C2H6O"},
                "source_view_scores": {
                    "library_match": {"source-a": 0.8, "source-b": 0.8},
                    "sensitive": {"source-a": 0.8, "source-b": 0.8},
                },
                "source_reference_content_fingerprints": {
                    "library_match": {
                        "source-a": "challenger-a",
                        "source-b": "challenger-b",
                    },
                    "sensitive": {
                        "source-a": "challenger-a",
                        "source-b": "challenger-b",
                    },
                },
            },
        }
        for identity in expanded:
            candidates[identity] = {
                "candidate": {"formula": "C2H6O"},
                "source_view_scores": {
                    "library_match": {"source-a": 0.9, "source-b": 0.9},
                    "sensitive": {"source-a": 0.9, "source-b": 0.9},
                },
            }
        scan_results.append({"condition_id": condition, "candidates": candidates})

    consensus = _source_paired_consensus(
        existing,
        scan_results,
        load_evirecall_policy(),
        expanded_candidate_connectivities=expanded,
    )

    assert consensus["status"] == "qualified_switch"
    assert consensus["candidate_scope_applied"] == "v11_native_same_formula_only"
    assert consensus["native_same_formula_connectivity_count"] == 2
    assert consensus["expanded_same_formula_challenger_count"] == 0
    assert consensus["same_formula_connectivity_count"] == 2
    assert consensus["selected_identity_key"] == f"std2d:{challenger}"


def test_same_connectivity_in_different_detection_methods_is_not_preexcluded(
    synthetic_library,
):
    path, identities = synthetic_library
    retriever = _retriever(path)

    method_a = retriever.retrieve_entity(
        "entity-method-a", [_query("a")], detection_method_id="method-a"
    )
    method_b = retriever.retrieve_entity(
        "entity-method-b", [_query("b")], detection_method_id="method-b"
    )

    assert method_a["detection_method_id"] == "method-a"
    assert method_b["detection_method_id"] == "method-b"
    for result in (method_a, method_b):
        assert identities["correct"] in {
            item["connectivity_identity"]
            for item in result["expanded_candidates"]
        }


def test_existing_connectivity_is_not_added_twice_within_one_entity_method(
    synthetic_library,
):
    path, identities = synthetic_library
    correct_key = identities["correct"]
    existing = [{"identity_key": f"std2d:{correct_key}", "candidate": {}}]

    result = _retriever(path).retrieve_entity(
        "entity-1",
        [_query()],
        existing_candidates=existing,
        detection_method_id="method-a",
    )

    assert result["existing_connectivity_count"] == 1
    assert correct_key not in {
        item["connectivity_identity"] for item in result["expanded_candidates"]
    }


def _minimal_v11_and_rows() -> tuple[dict, list[dict]]:
    rows = [
        {
            "spectrum_id": "query-1",
            "precursor_mz": 100.0,
            "ion_mode": "POS",
            "precursor_charge": 1,
            "collision_energy": 10.0,
            "collision_energy_unit": "eV",
            "raw_peaks": QUERY_PEAKS,
        }
    ]
    v11 = {
        "metadata": {
            "task_id": "candidate-recall-test",
            "algorithm_id": "mssystem-v11-flagship-chem-entity-ranker-1.2.1",
            "policy_sha256": "b" * 64,
            "source_semantic_sha256": "a" * 64,
            "source_row_count": 1,
        },
        "entities": [
            {
                "entity_id": "entity-1",
                "member_result_indices": [0],
                "member_spectrum_ids": ["query-1"],
                "candidates": [],
            }
        ],
    }
    return v11, rows


def test_bundle_requires_exact_v11_raw_membership_partition(synthetic_library):
    path, _ = synthetic_library
    v11, rows = _minimal_v11_and_rows()
    retriever = _retriever(path)

    bundle = build_candidate_recall_bundle(
        v11,
        rows,
        retriever=retriever,
        detection_method_id="method-a",
    )

    assert bundle["metadata"]["source_membership_partition_validated"] is True
    assert bundle["metadata"]["entity_count"] == 1
    assert bundle["entities"][0]["detection_method_id"] == "method-a"
    assert bundle["entities"][0]["expanded_candidate_count"] == 2

    mismatched = deepcopy(rows)
    mismatched[0]["spectrum_id"] = "another-acquisition"
    with pytest.raises(ValueError, match="Raw spectrum id mismatch"):
        build_candidate_recall_bundle(
            v11,
            mismatched,
            retriever=retriever,
            detection_method_id="method-a",
        )


def test_atomic_candidate_recall_sidecar_preserves_all_three_inputs(
    synthetic_library, tmp_path: Path
):
    database, _ = synthetic_library
    v11, rows = _minimal_v11_and_rows()
    v11_path = tmp_path / "candidate-recall-test_v11_flagship.json"
    source_path = tmp_path / "candidate-recall-test_results.json"
    output = tmp_path / "candidate-recall-test_v12_candidate_recall.json"
    v11_path.write_text(json.dumps(v11), encoding="utf-8")
    source_path.write_text(json.dumps(rows), encoding="utf-8")
    before = {
        "v11": _sha256(v11_path),
        "source": _sha256(source_path),
        "database": _sha256(database),
    }

    summary = write_evirecall_artifact(
        "candidate-recall-test",
        v11_path,
        source_path,
        database,
        output,
        library_snapshot_id="synthetic-library-sha256",
        library_sampled_sha256=_sampled_file_sha256(database),
        detection_method_id="method-a",
    )
    page = read_evirecall_artifact(output, limit=10)

    assert summary["database_sampled_fingerprint_unchanged"] is True
    assert summary["production_effect"] is False
    assert page["entities_total"] == 1
    assert page["entities"][0]["expanded_candidate_count"] == 2
    assert before == {
        "v11": _sha256(v11_path),
        "source": _sha256(source_path),
        "database": _sha256(database),
    }
    assert not list(tmp_path.glob(".*.tmp-*"))
