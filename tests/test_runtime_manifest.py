from __future__ import annotations

import json
import os
import sqlite3
from pathlib import Path

from runtime.golden import (
    semantic_result_digest,
    spectrum_projection,
    task_summary_projection,
    verify_golden_baseline,
)
from runtime.manifest import (
    library_snapshot_freshness,
    refresh_runtime_manifest_source,
    role_catalog_snapshot,
    runtime_manifest_freshness,
    source_tree_snapshot,
    sqlite_snapshot,
)


def test_source_tree_snapshot_is_order_invariant_and_excludes_its_output(tmp_path: Path):
    for directory in ("matching", "config", "web"):
        (tmp_path / directory).mkdir()
    (tmp_path / "matching" / "b.py").write_text("B = 2\n", encoding="utf-8")
    (tmp_path / "matching" / "a.py").write_text("A = 1\n", encoding="utf-8")
    (tmp_path / "config" / "runtime_manifest.json").write_text("{}\n", encoding="utf-8")
    (tmp_path / "config" / "golden_baseline.json").write_text("{}\n", encoding="utf-8")
    (tmp_path / "web" / "index.html").write_text("<p>ok</p>\n", encoding="utf-8")

    first = source_tree_snapshot(tmp_path)
    (tmp_path / "config" / "runtime_manifest.json").write_text('{"changed":true}\n', encoding="utf-8")
    (tmp_path / "config" / "golden_baseline.json").write_text('{"changed":true}\n', encoding="utf-8")
    second = source_tree_snapshot(tmp_path)

    assert first["source_tree_sha256"] == second["source_tree_sha256"]
    assert [row["path"] for row in first["files"]] == [
        "matching/a.py",
        "matching/b.py",
        "web/index.html",
    ]


def test_source_snapshot_includes_versioned_contract_and_governance_packages():
    root = Path(__file__).resolve().parents[1]
    snapshot = source_tree_snapshot(root)
    paths = {item["path"] for item in snapshot["files"]}

    assert "schemas/identification.py" in paths
    assert "pipeline/result_contract.py" in paths
    assert "calibration/project.py" in paths
    assert "governance/claim_compiler.py" in paths


def test_source_only_manifest_refresh_detects_subsequent_source_change(tmp_path: Path):
    for directory in ("matching", "config", "web"):
        (tmp_path / directory).mkdir()
    (tmp_path / "web" / "app.py").write_text("VALUE = 1\n", encoding="utf-8")
    (tmp_path / "matching" / "compound_annotations.json").write_text(
        json.dumps(
            {
                "catalog_id": "roles-test",
                "catalog_version": "1",
                "annotations": {
                    "AAAAAAAAAAAAAA-BBBBBBBBBB-C": {
                        "compound_role": "external_contaminant",
                        "role_annotation_scope": "exact_identity_curation",
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "matching" / "screening_catalog.json").write_text(
        json.dumps(
            {
                "catalog_id": "screening-test",
                "released_at": "2026-08-04",
                "entries": [
                    {
                        "full_inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                        "compound_role": "pesticide_screening_target",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    manifest_path = tmp_path / "config" / "runtime_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "library": {
                    "production_full": {"snapshot_id": "sqlite-preserved"},
                    "mapping": {"sha256": "mapping-preserved"},
                }
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "config" / "model_registry.json").write_text(
        json.dumps({"models": []}),
        encoding="utf-8",
    )
    (tmp_path / "config" / "golden_baseline.json").write_text(
        json.dumps({"schema_version": "1.0.0", "tasks": {"frozen": {}}}),
        encoding="utf-8",
    )

    refreshed = refresh_runtime_manifest_source(tmp_path, manifest_path)
    status = runtime_manifest_freshness(tmp_path, manifest_path)

    assert refreshed["library"]["production_full"]["snapshot_id"] == (
        "sqlite-preserved"
    )
    assert refreshed["models"]["registry"] == str(
        (tmp_path / "config" / "model_registry.json").resolve()
    )
    assert refreshed["golden_baseline"]["path"] == str(
        (tmp_path / "config" / "golden_baseline.json").resolve()
    )
    assert refreshed["golden_baseline"]["task_count"] == 1
    assert refreshed["role_catalogs"]["compound_roles"][
        "external_contaminant_count"
    ] == 1
    assert status["fresh"] is True

    (tmp_path / "web" / "app.py").write_text("VALUE = 2\n", encoding="utf-8")
    stale = runtime_manifest_freshness(tmp_path, manifest_path)
    assert stale["fresh"] is False
    assert stale["status"] == "stale"


def test_bundled_role_catalog_snapshot_exposes_release_counts():
    root = Path(__file__).resolve().parents[1]
    catalogs = role_catalog_snapshot(root)

    assert catalogs["compound_roles"]["external_contaminant_count"] == 95
    assert catalogs["compound_roles"]["class_level_structure_rule_count"] == 66
    assert catalogs["screening"]["entry_count"] == 18
    assert catalogs["screening"]["indexed_full_inchikey_count"] == 20
    assert catalogs["screening"]["library_coverage"]["status"] == "current"
    assert catalogs["screening"]["library_coverage"][
        "catalog_full_keys_with_searchable_spectra"
    ] == 20
    assert catalogs["screening"]["library_coverage"][
        "catalog_full_key_coverage_fraction"
    ] == 1.0


def test_online_runtime_snapshot_identity_refreshes_after_atomic_manifest_update(
    tmp_path: Path,
    monkeypatch,
):
    from web import server

    config = tmp_path / "config"
    config.mkdir()
    manifest = config / "runtime_manifest.json"
    manifest.write_text(
        json.dumps({
            "runtime_snapshot_id": "v10-first",
            "source": {"source_tree_sha256": "same-source"},
        }),
        encoding="utf-8",
    )
    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(
        server,
        "RUNTIME_SOURCE_TREE_SHA256",
        "same-source",
    )
    monkeypatch.setattr(
        server,
        "_runtime_snapshot_cache",
        {"signature": None, "value": "v10-bootstrap"},
    )

    assert server._current_runtime_snapshot_id() == "v10-first"

    replacement = config / ".runtime_manifest.json.next"
    replacement.write_text(
        json.dumps({
            "runtime_snapshot_id": "v10-final",
            "source": {"source_tree_sha256": "same-source"},
        }),
        encoding="utf-8",
    )
    os.replace(replacement, manifest)

    assert server._current_runtime_snapshot_id() == "v10-final"


def test_online_runtime_snapshot_identity_requires_restart_after_source_change(
    tmp_path: Path,
    monkeypatch,
):
    from web import server

    config = tmp_path / "config"
    config.mkdir()
    manifest = config / "runtime_manifest.json"
    manifest.write_text(
        json.dumps({
            "runtime_snapshot_id": "v12-loaded",
            "source": {"source_tree_sha256": "loaded-source"},
        }),
        encoding="utf-8",
    )
    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(
        server,
        "RUNTIME_SOURCE_TREE_SHA256",
        "loaded-source",
    )
    monkeypatch.setattr(
        server,
        "_runtime_snapshot_cache",
        {"signature": None, "value": "v12-bootstrap"},
    )

    assert server._current_runtime_snapshot_id() == "v12-loaded"

    replacement = config / ".runtime_manifest.json.next"
    replacement.write_text(
        json.dumps({
            "runtime_snapshot_id": "v12-disk-only",
            "source": {"source_tree_sha256": "edited-source"},
        }),
        encoding="utf-8",
    )
    os.replace(replacement, manifest)

    assert server._current_runtime_snapshot_id() == "v12-loaded"


def test_sqlite_snapshot_changes_when_database_changes(tmp_path: Path):
    path = tmp_path / "library.db"
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE spectra(id INTEGER PRIMARY KEY, name TEXT)")
        connection.execute("INSERT INTO spectra(name) VALUES ('a')")
    first = sqlite_snapshot(path)
    with sqlite3.connect(path) as connection:
        connection.execute("INSERT INTO spectra(name) VALUES ('b')")
    second = sqlite_snapshot(path)

    assert first["row_count"] == 1
    assert second["row_count"] == 2
    assert first["snapshot_id"] != second["snapshot_id"]


def test_library_freshness_detects_tracked_target_content_change(tmp_path: Path):
    path = tmp_path / "library.db"
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE Spectra(id INTEGER PRIMARY KEY, name TEXT)")
        connection.execute("CREATE TABLE Compounds(id INTEGER PRIMARY KEY, name TEXT)")
        connection.execute(
            "CREATE TABLE TargetCompounds(id INTEGER PRIMARY KEY, name TEXT)"
        )
        connection.execute(
            "CREATE TABLE DiagnosticIons(id INTEGER PRIMARY KEY, target_id INTEGER, mz REAL)"
        )
        connection.execute("INSERT INTO Spectra(name) VALUES ('spectrum')")
        connection.execute("INSERT INTO Compounds(name) VALUES ('compound')")
        connection.execute("INSERT INTO TargetCompounds(name) VALUES ('target-a')")
        connection.execute(
            "INSERT INTO DiagnosticIons(target_id, mz) VALUES (1, 123.4)"
        )

    frozen = sqlite_snapshot(path)
    current = library_snapshot_freshness(frozen)
    assert current["matches"] is True
    assert frozen["tracked_table_row_counts"] == {
        "Compounds": 1,
        "DiagnosticIons": 1,
        "Spectra": 1,
        "TargetCompounds": 1,
    }

    with sqlite3.connect(path) as connection:
        connection.execute(
            "UPDATE TargetCompounds SET name='target-b' WHERE id=1"
        )

    stale = library_snapshot_freshness(frozen)
    assert stale["matches"] is False
    assert stale["status"] == "stale"
    assert stale["manifest_table_content_sha256"] != (
        stale["current_table_content_sha256"]
    )


def test_semantic_projection_ignores_peaks_but_detects_identity_change(tmp_path: Path):
    path = tmp_path / "task_results.json"
    record = {
        "spectrum_id": "scan=1",
        "precursor_mz": 123.4,
        "retention_time": 1.2,
        "query_peaks": [[10.0, 999.0]],
        "match": {
            "compound_name": "A",
            "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            "formula": "C1H2",
            "rank_score": 0.8,
            "confidence": "High",
            "isopair_evidence": {"evaluated": True, "runtime_influence": False},
            "evidence_governance": {"claim_status": "isomer_group_only"},
        },
    }
    path.write_text(json.dumps([record]), encoding="utf-8")
    first = semantic_result_digest(path)
    record["query_peaks"] = [[999.0, 1.0]]
    path.write_text(json.dumps([record]), encoding="utf-8")
    second = semantic_result_digest(path)
    record["match"]["inchikey"] = "BBBBBBBBBBBBBB-UHFFFAOYSA-N"
    path.write_text(json.dumps([record]), encoding="utf-8")
    third = semantic_result_digest(path)

    assert first["semantic_sha256"] == second["semantic_sha256"]
    assert second["semantic_sha256"] != third["semantic_sha256"]
    assert first["isopair_evaluated"] == 1


def test_golden_verifier_reports_exact_scientific_path(tmp_path: Path):
    task_id = "abc123"
    task = {
        "filename": "sample.mzML",
        "status": "completed",
        "total_spectra": 1,
        "matched_spectra": 1,
        "unique_compounds": 1,
    }
    baseline = {
        "tasks": {task_id: {"summary": task_summary_projection(task)}}
    }
    baseline_path = tmp_path / "golden.json"
    index_path = tmp_path / "tasks_index.json"
    results_root = tmp_path / "tasks"
    results_root.mkdir()
    baseline_path.write_text(json.dumps(baseline), encoding="utf-8")
    changed = dict(task, matched_spectra=0)
    index_path.write_text(json.dumps({task_id: changed}), encoding="utf-8")

    report = verify_golden_baseline(
        baseline_path,
        index_path,
        results_root,
        verify_semantic_digests=False,
    )

    assert not report["ok"]
    assert any(
        row["path"] == f"tasks.{task_id}.summary.matched_spectra"
        for row in report["differences"]
    )


def test_spectrum_projection_separates_runtime_and_counterfactual_influence():
    projection = spectrum_projection(
        {
            "spectrum_id": "1",
            "match": {
                "inchikey": "A",
                "isopair_evidence": {
                    "evaluated": True,
                    "runtime_influence": False,
                    "ranking_changed": False,
                    "counterfactual_ranking_changed": True,
                },
            },
        }
    )
    assert projection["match"]["isopair"]["runtime_influence"] is False
    assert projection["match"]["isopair"]["counterfactual_ranking_changed"] is True


def test_frozen_runtime_and_golden_files_are_machine_readable():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "config" / "runtime_manifest.json").read_text())
    golden = json.loads((root / "config" / "golden_baseline.json").read_text())
    frozen_root = Path(manifest["project_root"])

    assert manifest["runtime_snapshot_id"].startswith("v16-")
    assert manifest["models"]["registry"] == str(
        (frozen_root / "config" / "model_registry.json").resolve()
    )
    assert manifest["golden_baseline"]["path"] == str(
        (frozen_root / "config" / "golden_baseline.json").resolve()
    )
    assert len(manifest["source"]["source_tree_sha256"]) == 64
    assert manifest["role_catalogs"]["compound_roles"][
        "external_contaminant_count"
    ] == 95
    assert manifest["role_catalogs"]["screening"][
        "indexed_full_inchikey_count"
    ] == 20
    assert manifest["library"]["production_full"]["snapshot_id"].startswith("sqlite-")
    assert len(golden["tasks"]) == 10
    assert all(
        row["result_semantics"]["isopair_runtime_ranking_changed"] == 0
        for task_id, row in golden["tasks"].items()
        if task_id not in {"92becdbf0217", "8a9d1d22d31a"}
    )


def test_frozen_runtime_manifest_matches_current_source_and_catalogs():
    root = Path(__file__).resolve().parents[1]
    status = runtime_manifest_freshness(
        root,
        root / "config" / "runtime_manifest.json",
    )

    assert status["fresh"] is True, status
    assert status["release_ready"] is True, status
    assert status["screening_coverage_current"] is True
    assert status["authentic_reference"]["passed"] is True
    assert status["authentic_reference"]["entry_count"] == 77
    assert status["manifest_source_tree_sha256"] == (
        status["current_source_tree_sha256"]
    )
    assert status["catalog_matches"] == {
        "compound_roles": True,
        "screening": True,
    }
