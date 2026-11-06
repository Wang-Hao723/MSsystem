import copy
import csv
import hashlib
import io
import json
import lzma
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from web import frozen_batch_reader as reader
from web.frozen_batch_tasks import FrozenBatchTasks, install_frozen_task_access


PROTOCOL = "1d09e21e37dad4c856a70c6b09a8ccecf611aa5cadc40567f8acd82767fbe779"


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _write_xz(path, payload):
    compressed = lzma.compress(json.dumps(payload, ensure_ascii=False).encode())
    path.write_bytes(compressed)
    return {"path": str(path), "size": len(compressed), "sha256": _sha(compressed)}


def _summary(sample_key, index):
    return {
        "parsed_searched_ms2_spectra": 3,
        "matched_ms2_spectra": 2,
        "ms1_feature_records": 2,
        "ms1_evidence_level_counts": {"MS1_EXACT_MASS_CANDIDATE": 1, "MS1_ISOTOPE_SUPPORTED": 1},
        "blank_flagged_ms2": 1,
        "blank_clean_view_excluded_ms2": 1,
        "source_sample": {"sample_key": sample_key, "sample_id": f"S{index:03d}", "polarity": "POS", "replicate": 1},
    }


@pytest.fixture
def bridge_fixture(tmp_path):
    native_samples = {}
    conclusion_samples = {}
    for index in range(149):
        key = f"sample-{index:03d}__{index:016x}"
        native_path = f"/data/Mssystem/xbt_real_batch_20261021_ms2/results/{key}.json.xz"
        conclusion_path = f"/data/Mssystem/xbt_real_batch_20261021_ms2/conclusions/{key}.json.xz"
        native_sha = f"{index + 1:064x}"
        native_samples[key] = {
            "status": "completed",
            "source_path": f"/frozen/raw/S{index:03d}.mzML",
            "source_sha256": f"{index + 1000:064x}",
            "started_at": "2026-10-21T00:00:00+00:00",
            "verified_remote_path": native_path,
            "verified_remote_size": 1000 + index,
            "verified_remote_sha256": native_sha,
        }
        conclusion_samples[key] = {
            "status": "completed",
            "finished_at": "2026-10-21T01:00:00+00:00",
            "native_remote_path": native_path,
            "native_remote_sha256": native_sha,
            "verified_remote_path": conclusion_path,
            "verified_remote_size": 2000 + index,
            "verified_remote_sha256": f"{index + 2000:064x}",
            "summary": _summary(key, index),
        }
    native_state = tmp_path / "STATE.json"
    conclusion_state = tmp_path / "CONCLUSIONS_STATE.json"
    terminal_receipt = tmp_path / "CONCLUSIONS_TERMINAL_RECEIPT.json"
    native_state.write_text(json.dumps({"state": "completed", "protocol_sha256": PROTOCOL, "samples": native_samples}))
    conclusion_state.write_text(json.dumps({"state": "completed", "protocol_sha256": "conclusion-protocol", "samples": conclusion_samples}))
    terminal_receipt.write_text(json.dumps({"state": "completed", "sample_count": 149}))
    (tmp_path / "V16_INTEGRATION_STATE.json").write_text(json.dumps({
        "baseline_frozen_artifact_sha256": {
            "STATE.json": _sha(native_state.read_bytes()),
            "CONCLUSIONS_STATE.json": _sha(conclusion_state.read_bytes()),
            "CONCLUSIONS_TERMINAL_RECEIPT.json": _sha(terminal_receipt.read_bytes()),
        }
    }))
    blank_path = "/frozen/raw/BLANK_POS1.mzML"
    manifest = [
        {"path": row["source_path"], "role": "sample", "observed_mode": "POS"}
        for row in native_samples.values()
    ] + [{"path": blank_path, "role": "blank", "observed_mode": "POS"}]
    input_manifest = tmp_path / "INPUT_MANIFEST.json"
    blank_pairing = tmp_path / "BLANK_PAIRING.json"
    input_manifest.write_text(json.dumps(manifest))
    blank_pairing.write_text(json.dumps([{
        "samples": [row["source_path"] for row in native_samples.values()],
        "blanks": [blank_path], "blank_filter_policy": "soft_evidence",
        "method_and_batch_comparability_user_confirmed": True,
        "hard_background_exclusion_confirmed": False,
    }]))

    class Service:
        batch_root = tmp_path
        paths = {"native_state": native_state, "conclusions_state": conclusion_state, "report": tmp_path / "report.html"}
        index = {"xbt_dataset": "XBT_targeted_20260817"}
        remote_host = "unused"
        packaging = {"html_sha256": "unused"}
        input_manifest_path = input_manifest
        blank_pairing_path = blank_pairing

        def detail_payload(self):
            assert set(native_samples) == set(conclusion_samples)
            return {"samples": [{}] * 149}

        def _summary(self):
            return {"report_url": "/report", "download_all_url": "/download-all"}

    bridge = FrozenBatchTasks(Service())
    return SimpleNamespace(bridge=bridge, native=native_samples, conclusions=conclusion_samples)


@pytest.fixture
def reader_fixture(tmp_path, monkeypatch):
    key = "sample-000__0000000000000000"
    results_dir = tmp_path / "results"
    conclusions_dir = tmp_path / "conclusions"
    results_dir.mkdir(); conclusions_dir.mkdir()
    source_sha = "1" * 64
    native = {
        "sample_key": key,
        "source": {"sha256": source_sha, "filename": "S000.mzML"},
        "protocol_sha256": PROTOCOL,
        "parameters": {"db_mode": "full"},
        "native_summary": {"searched": 3},
        "paired_blanks": {"files": 2},
        "runtime": {"version": "16"},
        "native_ms2_results": [
            {"spectrum_id": "scan-10", "match": {"name": "A", "cosine": 0.91}, "blank_flagged": False, "blank_clean_view_excluded": False},
            {"spectrum_id": "scan-11", "match": None, "blank_flagged": True, "blank_clean_view_excluded": True},
            {"spectrum_id": "scan-12", "match": {"name": "B", "cosine": 0.52}, "blank_flagged": True, "blank_clean_view_excluded": False},
        ],
        "ms1_integration": {
            "features": [
                {"feature_id": "ms1-only", "evidence_level": "MS1_EXACT_MASS_CANDIDATE", "associated_ms2_spectrum_index": None},
                {"feature_id": "ms1-ms2", "evidence_level": "MS1_ISOTOPE_SUPPORTED", "associated_ms2_spectrum_index": 0},
            ],
            "groups": [],
        },
    }
    result_record = _write_xz(results_dir / f"{key}.json.xz", native)
    conclusion_counts = {
        "high": 2, "isomer": 22, "candidate": 27, "background": 117,
        "unresolved": 19, "unannotated": 1123,
    }
    native_compounds = []
    for conclusion_key, count in conclusion_counts.items():
        for _ in range(count):
            index = len(native_compounds)
            grouped = conclusion_key in {"high", "isomer", "candidate", "background"}
            unresolved = conclusion_key == "unresolved"
            candidate_members = ([{"compound_name": f"saved member {index}", "cosine": 0.81,
                                   "source_evidence": {"library": "frozen"}}]
                                 if grouped or unresolved else [])
            native_compounds.append({
                "compound_id": f"compound:{index}", "result_id": f"result:{index}",
                "feature_id": f"feature:{index}", "record_unit": "ion_feature",
                "record_type": ("compound_or_isomer_group" if grouped else
                                "unresolved_precursor_event" if unresolved else "unannotated_ion_feature"),
                "identity_key": f"identity:{index}", "identity_scope": "frozen_scope",
                "identity_source": "saved_native_conclusion",
                "display_name": ("Unresolved precursor event" if unresolved else f"saved row {index}"),
                "formula": "C1H2" if grouped else "", "evidence_scope": "saved_native",
                "lead_candidate": ({"compound_name": f"lead {index}", "formula": "C1H2",
                                    "cosine": 0.81, "rank": 1} if candidate_members else {}),
                "candidate_members": candidate_members,
                "conclusion": {"key": conclusion_key, "confirmation_level": "candidate_only",
                               "counted_as_sample_compound": grouped,
                               "counted_as_plant_candidate": conclusion_key in {"high", "isomer", "candidate"}},
                "plant_report": {"excluded_from_plant_candidate_table": conclusion_key == "background",
                                 "review_category": conclusion_key, "reason": "saved_flag"},
                "representative": {"result_index": index, "spectrum_id": f"scan-{index}",
                                   "blank_flagged": conclusion_key == "background",
                                   "blank_clean_view_excluded": False,
                                   "raw_peaks": [[10.0, 20.0]],
                                   "match": {"cosine": 0.81, "matched_peaks": [[10.0, 10.0]]}},
                "observations": [{"result_index": index, "raw_peaks": [[10.0, 20.0]]}],
            })
    native_metadata = {"source": "frozen-native", "compound_count": 1310,
                       "main_result_count": 168, "unresolved_event_count": 19,
                       "unannotated_count": 1123, "preserve_this_exact_field": "yes"}
    conclusion = {
        "sample_key": key,
        "native_result_binding": {
            "compressed_sha256": result_record["sha256"], "compressed_size": result_record["size"],
            "remote_path": result_record["path"], "source_sha256": source_sha,
            "native_protocol_sha256": PROTOCOL,
        },
        "report_views": {
            "xbt_named_candidates": {"metadata": {"selection_rule": "frozen"}, "compounds": [
                copy.deepcopy(native_compounds[168])
            ]},
            "external_or_incompatible_qc": {"metadata": {}, "compounds": []},
        },
        "compound_conclusions": {"metadata": native_metadata, "compounds": native_compounds},
        "summary": {"xbt_named_candidate_rows": 1},
        "primary_identity_artifact": {"entities": [{"entity_id": "e1", "decision": {"status": "candidate"}}]},
    }
    conclusion_record = _write_xz(conclusions_dir / f"{key}.json.xz", conclusion)
    request = {
        "sample_key": key, "source_sha256": source_sha, "native_protocol_sha256": PROTOCOL,
        "results": result_record, "conclusions": conclusion_record,
    }
    monkeypatch.setattr(reader, "ROOT", tmp_path)
    return SimpleNamespace(request=request, native=native, conclusion=conclusion,
                           result_path=results_dir / f"{key}.json.xz",
                           conclusion_path=conclusions_dir / f"{key}.json.xz")


def test_reader_preserves_full_item_indices_and_includes_ms1_only_features(reader_fixture):
    request = reader_fixture.request
    summary = reader.project({**request, "operation": "summary", "row_fields": ["spectrum_id", "blank_flagged"], "match_fields": ["name", "cosine"]})
    assert [row["result_index"] for row in summary["results"]] == [0, 1, 2]
    assert summary["results"][1]["match"] is None
    item = reader.project({**request, "operation": "item", "index": 2})
    assert item == {"result": reader_fixture.native["native_ms2_results"][2], "result_index": 2}
    expected = reader.project({**request, "operation": "expected"})
    assert [row["feature_id"] for row in expected["features"]] == ["ms1-only", "ms1-ms2"]
    assert expected["features"][0]["associated_ms2_spectrum_index"] is None


def test_reader_preserves_saved_compound_conclusion_and_evidence(reader_fixture):
    payload = reader.project({**reader_fixture.request, "operation": "compound_page",
                              "result_view": "xbt_linked", "report_scope": "all",
                              "offset": 0, "limit": 50})
    row = payload["compounds"][0]
    assert row["record_type"] == "unresolved_precursor_event"
    assert row["conclusion"]["key"] == "unresolved"
    assert row["conclusion"]["confirmation_level"] == "candidate_only"
    assert payload["native_metadata"] == reader_fixture.conclusion["compound_conclusions"]["metadata"]
    assert payload["native_total"] == 1310
    assert payload["xbt_linked_total"] == 1
    assert payload["selection"]["partial_view"] is True


def test_native_compound_pages_preserve_all_1310_indices_metadata_and_full_rows(reader_fixture):
    request = reader_fixture.request
    pages = [reader.project({**request, "operation": "compound_page", "result_view": "native",
                             "report_scope": "all", "offset": offset, "limit": 500})
             for offset in (0, 500, 1000)]
    assert [page["native_total"] for page in pages] == [1310, 1310, 1310]
    assert [len(page["compounds"]) for page in pages] == [500, 500, 310]
    assert [row["frozen_compound_index"] for page in pages for row in page["compounds"]] == list(range(1310))
    assert pages[0]["metadata"] == reader_fixture.conclusion["compound_conclusions"]["metadata"]
    assert pages[0]["metadata"]["preserve_this_exact_field"] == "yes"
    assert "observations" not in pages[0]["compounds"][0]
    assert "raw_peaks" not in pages[0]["compounds"][0]["representative"]
    full = reader.project({**request, "operation": "compound_item", "result_view": "native", "index": 0})
    assert full["compound"] == reader_fixture.conclusion["compound_conclusions"]["compounds"][0]
    assert full["compound"]["observations"]
    assert full["compound"]["candidate_members"]


def test_native_compound_saved_flag_filters_and_exact_exports(reader_fixture):
    request = reader_fixture.request
    external = reader.project({**request, "operation": "compound_page", "result_view": "native",
                               "report_scope": "external_review", "offset": 0, "limit": 500})
    assert external["filtered_total"] == 117
    assert all(row["plant_report"]["excluded_from_plant_candidate_table"]
               for row in external["compounds"])
    plant = reader.project({**request, "operation": "compound_page", "result_view": "native",
                            "report_scope": "plant_candidates", "offset": 0, "limit": 500})
    assert plant["filtered_total"] == 51
    csv_bytes = reader.project({**request, "operation": "native_compounds_csv"})
    rows = list(csv.DictReader(io.StringIO(csv_bytes.decode("utf-8-sig"))))
    assert len(rows) == 1310
    assert rows[168]["record_type"] == "unresolved_precursor_event"
    assert json.loads(rows[168]["candidate_members_json"]) == [
        {"compound_name": "saved member 168", "cosine": 0.81,
         "source_evidence": {"library": "frozen"}}
    ]
    native_xz = reader.project({**request, "operation": "native_conclusions_xz"})
    assert native_xz == reader_fixture.conclusion_path.read_bytes()
    assert json.loads(lzma.decompress(native_xz)) == reader_fixture.conclusion


def test_reader_csv_projects_only_identity_and_scan_bound_lead_scores(reader_fixture):
    document = copy.deepcopy(reader_fixture.conclusion)
    accepted = document["compound_conclusions"]["compounds"][0]
    accepted.update({
        "identity_key": "std2d:AAAAAAAAAAAAAA|formula:C6H6",
        "identity_scope": "2d_connectivity", "formula": "C6H6",
    })
    accepted["lead_candidate"] = {
        "compound_name": "A", "formula": "C6H6",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "result_index": 0,
        "lib_spectrum_id": "LIB-0", "cosine": 0.934, "mass_ppm": -1.7,
    }
    accepted["representative"].update({"result_index": 0, "spectrum_id": "scan-0"})
    accepted["representative"]["match"] = {
        "compound_name": "A", "formula": "C6H6",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "lib_spectrum_id": "LIB-0",
        "spectral_evidence_binding": {
            "schema": "same-scan-selected-reference/1", "status": "bound_saved_reference",
            "spectrum_id": "scan-0", "result_index": 0,
            "candidate_inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "lib_spectrum_id": "LIB-0", "saved_confidence": "Low",
            "metrics": {"cosine": 0.934, "mass_ppm": -1.7},
        },
    }
    rejected = document["compound_conclusions"]["compounds"][1]
    rejected.update({
        "identity_key": "std2d:AAAAAAAAAAAAAA|formula:C6H6",
        "identity_scope": "2d_connectivity", "formula": "C6H6",
    })
    rejected["lead_candidate"] = {
        "compound_name": "B", "formula": "C6H6",
        "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E", "result_index": 2,
        "lib_spectrum_id": "LIB-2", "cosine": 0.71, "mass_ppm": 8.2,
    }
    rejected["representative"].update({"result_index": 1, "spectrum_id": "scan-1"})
    rejected["representative"]["match"] = {
        "compound_name": "A", "formula": "C6H6",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "lib_spectrum_id": "LIB-1",
        "cosine": 0.999, "mass_ppm": 0.1,
        "spectral_evidence_binding": {
            "schema": "same-scan-selected-reference/1", "status": "bound_saved_reference",
            "spectrum_id": "scan-1", "result_index": 1,
            "candidate_inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "lib_spectrum_id": "OTHER-LIB", "saved_confidence": "High",
            "metrics": {"cosine": 0.999, "mass_ppm": 0.1},
        },
    }
    record = _write_xz(reader_fixture.conclusion_path, document)
    request = {**reader_fixture.request, "conclusions": record,
               "operation": "native_compounds_csv"}
    rows = list(csv.DictReader(io.StringIO(reader.project(request).decode("utf-8-sig"))))
    assert rows[0]["lead_cosine"] == "0.934"
    assert rows[0]["lead_mass_ppm"] == "-1.7"
    assert rows[0]["lead_score_projection_bound"] == "True"
    assert rows[0]["saved_spectrum_confidence"] == "Low"
    assert rows[1]["lead_cosine"] == ""
    assert rows[1]["lead_mass_ppm"] == ""
    assert rows[1]["lead_cosine_raw"] == "0.71"
    assert rows[1]["lead_mass_ppm_raw"] == "8.2"
    assert rows[1]["lead_score_projection_bound"] == "False"
    assert rows[1]["representative_match_cosine_raw"] == "0.999"


def test_reader_routes_explicit_external_roles_to_review_only():
    external = {
        "conclusion": {"counted_as_plant_candidate": True},
        "plant_report": {}, "lead_candidate": {"compound_role": "pesticide_screening_target"},
        "representative": {"match": {}},
    }
    ordinary = {
        "conclusion": {"counted_as_plant_candidate": True},
        "plant_report": {}, "lead_candidate": {"compound_role": "unclassified"},
        "representative": {"match": {}},
    }
    assert reader._compound_in_scope(external, "external_review") is True
    assert reader._compound_in_scope(external, "plant_candidates") is False
    assert reader._compound_in_scope(ordinary, "external_review") is False
    assert reader._compound_in_scope(ordinary, "plant_candidates") is True


def test_primary_list_is_compact_but_full_saved_entity_and_decisions_remain_exact(reader_fixture):
    document = copy.deepcopy(reader_fixture.conclusion)
    entity = document["primary_identity_artifact"]["entities"][0]
    entity.update({"candidate_recall_support": {"large": ["evidence"] * 100},
                   "candidates": [{"raw": [1] * 100}],
                   "winner": {"entity_score": 0.876, "identity_key": "original",
                              "candidate": {"compound_name": "Saved candidate", "formula": "C6H6"},
                              "fragmentation_explanation": {"original": [1] * 100}}})
    record = _write_xz(reader_fixture.conclusion_path, document)
    request = {**reader_fixture.request, "conclusions": record, "operation": "flagship"}
    compact = reader.project(request)["entities"][0]
    assert compact["decision"] == entity["decision"]
    assert compact["winner"]["entity_score"] == 0.876
    assert compact["_frozen_summary"] is True
    assert "candidates" not in compact
    assert "fragmentation_explanation" not in compact["winner"]
    full = reader.project({**request, "entity_id": "e1"})["entities"][0]
    assert full == entity


def test_reader_rejects_compressed_hash_source_protocol_and_native_binding_tamper(reader_fixture):
    request = reader_fixture.request
    bad_hash = copy.deepcopy(request)
    bad_hash["results"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA/size"):
        reader.project({**bad_hash, "operation": "summary", "row_fields": [], "match_fields": []})

    bad_source = copy.deepcopy(request)
    bad_source["source_sha256"] = "2" * 64
    with pytest.raises(ValueError, match="source SHA"):
        reader.project({**bad_source, "operation": "native_metadata"})

    bad_protocol = copy.deepcopy(request)
    bad_protocol["native_protocol_sha256"] = "3" * 64
    with pytest.raises(ValueError, match="protocol"):
        reader.project({**bad_protocol, "operation": "expected"})

    bad_binding = copy.deepcopy(request)
    bad_binding["results"]["path"] += ".other"
    with pytest.raises(ValueError, match="native binding"):
        reader.project({**bad_binding, "operation": "compounds", "view": "xbt_named_candidates"})


def test_registry_adds_149_without_overwriting_359_and_collision_is_strict(bridge_fixture):
    bridge = bridge_fixture.bridge
    registry = {f"existing-{index}": {"id": f"existing-{index}"} for index in range(359)}
    assert bridge.install(registry) == 149
    assert len(registry) == 508
    snapshot = copy.deepcopy(registry)
    assert bridge.install(registry) == 0
    assert registry == snapshot

    task_id = next(iter(bridge.tasks))
    collision = copy.deepcopy(bridge.tasks[task_id])
    collision["filename"] = "different.mzML"
    registry[task_id] = collision
    with pytest.raises(RuntimeError, match="collid|冲突|identity|metadata changed"):
        bridge.install(registry)


def test_frozen_task_projects_real_blank_pairing_and_explicit_target_selection(bridge_fixture):
    task = next(iter(bridge_fixture.bridge.tasks.values()))
    blank = task["blank_summary"]
    assert blank["assigned_blank_replicate_count"] == 1
    assert blank["assigned_blank_filenames"] == ["BLANK_POS1.mzML"]
    assert blank["blank_operation"] == "comparison_only"
    assert blank["excluded_spectra"] == 0
    assert blank["method_comparability_verified"] is True
    assert blank["background_exclusion_authorized"] is False
    assert task["params"]["expected_compound_dataset_selection_source"] == "user_explicit"
    assert task["expected_compounds_summary"]["target_dataset_selection"] == {
        "selection_source": "user_explicit",
        "selected_datasets": ["XBT_targeted_20260817"],
        "selection_changes_saved_results": False,
    }


def test_legacy_frozen_registration_is_upgraded_without_rejecting_history(bridge_fixture):
    bridge = bridge_fixture.bridge
    task_id = next(iter(bridge.tasks))
    registry = copy.deepcopy(bridge._legacy_tasks)
    for legacy_id, legacy_task in registry.items():
        legacy_task["completed_batch"]["registration_sha256"] = (
            bridge._legacy_registration_sha256[legacy_id]
        )
    assert bridge.install(registry) == 0
    assert registry[task_id]["blank_summary"]["assigned_blank_replicate_count"] == 1
    assert registry[task_id]["completed_batch"]["registration_sha256"] == (
        bridge.tasks[task_id]["completed_batch"]["registration_sha256"]
    )


def test_registry_rejects_changed_frozen_registration_source(bridge_fixture):
    bridge = bridge_fixture.bridge
    integration = bridge.service.batch_root / "V16_INTEGRATION_STATE.json"
    payload = json.loads(integration.read_text())
    payload["baseline_frozen_artifact_sha256"]["STATE.json"] = "0" * 64
    integration.write_text(json.dumps(payload))
    with pytest.raises(RuntimeError, match="登记文件SHA"):
        FrozenBatchTasks(bridge.service)


def test_deleted_imported_task_stays_removed_after_restart_and_can_be_restored(bridge_fixture, monkeypatch, tmp_path):
    monkeypatch.setattr("config.settings.RUNTIME_ROOT", tmp_path)
    bridge = bridge_fixture.bridge
    task_id = next(iter(bridge.tasks))
    audit = tmp_path / "runtime/task_deletion_audit.jsonl"
    audit.parent.mkdir()
    audit.write_text(json.dumps({"event": "batch_task_recycled", "deleted_task_ids": [task_id]}) + "\n")
    registry = {}
    assert bridge.install(registry) == 148
    assert task_id not in registry
    assert bridge.install(registry) == 0
    registry[task_id] = copy.deepcopy(bridge.tasks[task_id])
    assert bridge.install(registry) == 0
    assert task_id in registry


def test_imported_task_removal_uses_native_handlers_and_keeps_rerun_blocked(bridge_fixture, monkeypatch):
    bridge = bridge_fixture.bridge
    app = FastAPI()
    registry = copy.deepcopy(bridge.tasks)
    task_id = next(iter(registry))

    @app.delete("/api/tasks/{key}")
    async def remove(key):
        registry.pop(key)
        return {"deleted": key}

    @app.post("/api/tasks/batch/delete")
    async def batch_remove():
        return {"native_delete_handler_reached": True}

    install_frozen_task_access(app, bridge, lambda: registry, set(), set())
    client = TestClient(app)
    assert client.post("/api/tasks/batch/delete", json={"task_ids": [task_id]}).json() == {"native_delete_handler_reached": True}
    assert client.post(f"/api/tasks/{task_id}/analyze").status_code == 409
    assert client.delete(f"/api/tasks/{task_id}").json() == {"deleted": task_id}
    assert client.get(f"/api/tasks/{task_id}").status_code == 404
    assert task_id in bridge.tasks  # Frozen result bindings remain available to the batch archive.


def _route_client(bridge, monkeypatch):
    app = FastAPI()
    holder = {"tasks": {task_id: copy.deepcopy(task) for task_id, task in bridge.tasks.items()}}
    materialized = {"count": 0}

    @app.get("/api/tasks/{task_id}/dangerous-materializer")
    async def dangerous(task_id):
        materialized["count"] += 1
        return {"materialized": True}

    rows = {
        0: {"spectrum_id": "scan-10", "match": {"name": "A", "cosine": 0.91}, "blank_flagged": False, "blank_clean_view_excluded": False},
        1: {"spectrum_id": "scan-11", "match": None, "blank_flagged": True, "blank_clean_view_excluded": True},
        2: {"spectrum_id": "scan-12", "match": {"name": "B", "cosine": 0.52}, "blank_flagged": True, "blank_clean_view_excluded": False},
    }
    native_compounds = [
        {"compound_id": "c-native", "record_type": "compound_or_isomer_group",
         "conclusion": {"key": "candidate", "counted_as_plant_candidate": True},
         "plant_report": {"excluded_from_plant_candidate_table": False},
         "candidate_members": [{"compound_name": "saved"}],
         "lead_candidate": {"compound_name": "saved"}, "representative": {"result_index": 0}},
        {"compound_id": "c-unresolved", "record_type": "unresolved_precursor_event",
         "conclusion": {"key": "unresolved", "counted_as_plant_candidate": False},
         "plant_report": {"excluded_from_plant_candidate_table": False},
         "candidate_members": [{"compound_name": "ambiguous"}],
         "lead_candidate": {}, "representative": {"result_index": 1}},
    ]

    def fake_read(task_id, operation, **options):
        if operation == "summary":
            return {"results": [{**row, "result_index": index} for index, row in rows.items()], "total": 3}
        if operation == "item":
            return {"result": copy.deepcopy(rows[options["index"]]), "result_index": options["index"]}
        if operation == "items":
            return {"results": [{**copy.deepcopy(rows[index]), "result_index": index} for index in options["indices"]]}
        if operation == "compound_page":
            if options["result_view"] == "xbt_linked" and options["report_scope"] != "all":
                raise ValueError("XBT-linked subset cannot be reclassified by report_scope")
            source = native_compounds if options["result_view"] == "native" else native_compounds[1:]
            selected = list(enumerate(source))
            if options["report_scope"] == "plant_candidates":
                selected = [(index, row) for index, row in selected
                            if row["conclusion"]["counted_as_plant_candidate"]]
            elif options["report_scope"] == "external_review":
                selected = [(index, row) for index, row in selected
                            if row["plant_report"]["excluded_from_plant_candidate_table"]]
            if options.get("conclusion"):
                selected = [(index, row) for index, row in selected
                            if row["conclusion"]["key"] == options["conclusion"]]
            page = selected[options["offset"]:options["offset"] + options["limit"]]
            return {"metadata": {"source": "native-exact"}, "native_metadata": {"source": "native-exact"},
                    "compounds": [{**copy.deepcopy(row), "frozen_compound_index": index}
                                  for index, row in page],
                    "native_total": 2, "filtered_total": len(selected), "xbt_linked_total": 1,
                    "offset": options["offset"], "limit": options["limit"],
                    "has_more": options["offset"] + options["limit"] < len(selected),
                    "selection": {"result_view": options["result_view"],
                                  "report_scope": options["report_scope"],
                                  "partial_view": options["result_view"] == "xbt_linked"},
                    "summary": {"xbt_named_candidate_rows": 1}}
        if operation == "compound_item":
            source = native_compounds if options["result_view"] == "native" else native_compounds[1:]
            index = options.get("index")
            if index is None:
                index = next(i for i, row in enumerate(source)
                             if row["compound_id"] == options["compound_id"])
            return {"compound": copy.deepcopy(source[index]), "frozen_compound_index": index,
                    "result_view": options["result_view"]}
        if operation == "expected":
            return {"features": [
                {"feature_id": "ms1-only", "associated_ms2_spectrum_index": None, "evidence_level": "MS1_EXACT_MASS_CANDIDATE"},
                {"feature_id": "linked", "associated_ms2_spectrum_index": 0, "evidence_level": "MS1_ISOTOPE_SUPPORTED"},
            ], "groups": [], "summary": {"detected_feature_count": 2}, "status": "ok"}
        if operation == "flagship":
            return {"entities": []}
        raise AssertionError(operation)

    monkeypatch.setattr(bridge, "read", fake_read)
    monkeypatch.setattr(bridge, "read_export", lambda task_id, operation: (
        b"sample_key,compound_id\nsample,c-native\n" if operation == "native_compounds_csv"
        else lzma.compress(json.dumps({"sample_key": "sample", "compound_conclusions": {
            "metadata": {"source": "native-exact"}, "compounds": native_compounds}}).encode())
    ))
    install_frozen_task_access(app, bridge, lambda: holder["tasks"], {"spectrum_id", "blank_flagged", "blank_clean_view_excluded"}, {"name", "cosine"})
    return TestClient(app), holder, materialized


def test_middleware_uses_rebound_registry_and_blocks_individual_and_batch_reruns(bridge_fixture, monkeypatch):
    bridge = bridge_fixture.bridge
    client, holder, _ = _route_client(bridge, monkeypatch)
    task_id = next(iter(bridge.tasks))
    rebound = {task_id: copy.deepcopy(bridge.tasks[task_id])}
    rebound[task_id]["message"] = "from-rebound-registry"
    holder["tasks"] = rebound
    status = client.get(f"/api/tasks/{task_id}")
    assert status.status_code == 200
    assert status.json()["message"] == "from-rebound-registry"
    assert client.post(f"/api/tasks/{task_id}/analyze").status_code == 409
    assert client.post(f"/api/tasks/{task_id}/retry").status_code == 409
    batch = client.post("/api/tasks/batch/analyze", json={"task_ids": [task_id]})
    assert batch.status_code == 409


def test_middleware_filters_summary_then_fetches_original_full_indices(bridge_fixture, monkeypatch):
    bridge = bridge_fixture.bridge
    client, _, _ = _route_client(bridge, monkeypatch)
    task_id = next(iter(bridge.tasks))
    response = client.get(
        f"/api/tasks/{task_id}/results",
        params={"blank_view": "clean", "match_state": "matched", "min_cosine": "0.6", "projection": "full"},
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["results"] == [{"result_index": 0, "spectrum_id": "scan-10", "match": {"name": "A", "cosine": 0.91}, "blank_flagged": False, "blank_clean_view_excluded": False}]
    assert payload["total"] == 1
    single = client.get(f"/api/tasks/{task_id}/results/item/2")
    assert single.status_code == 200
    assert single.json()["result"]["spectrum_id"] == "scan-12"
    assert single.json()["result_index"] == 2


def test_middleware_projects_saved_compounds_and_all_ms1_features(bridge_fixture, monkeypatch):
    bridge = bridge_fixture.bridge
    client, _, _ = _route_client(bridge, monkeypatch)
    task_id = next(iter(bridge.tasks))
    compounds = client.get(f"/api/tasks/{task_id}/compounds")
    assert compounds.status_code == 200
    payload = compounds.json()
    assert payload["native_total"] == 2
    assert payload["xbt_linked_total"] == 1
    assert [row["compound_id"] for row in payload["compounds"]] == ["c-native", "c-unresolved"]
    assert payload["metadata"] == {"source": "native-exact"}
    linked = client.get(f"/api/tasks/{task_id}/compounds?result_view=xbt_linked")
    assert linked.status_code == 200
    assert linked.json()["compounds"][0]["record_type"] == "unresolved_precursor_event"
    assert linked.json()["selection"]["partial_view"] is True
    full = client.get(f"/api/tasks/{task_id}/compounds/item/1")
    assert full.status_code == 200
    expected_full = dict(payload["compounds"][1])
    expected_full.pop("frozen_compound_index")
    assert full.json()["compound"] == expected_full
    expected = client.get(f"/api/tasks/{task_id}/expected-compounds")
    assert expected.status_code == 200
    assert [row["feature_id"] for row in expected.json()["features"]] == ["ms1-only", "linked"]
    assert expected.json()["features"][0]["associated_ms2_spectrum_index"] is None


def test_native_compound_csv_and_complete_json_routes_use_saved_exports(bridge_fixture, monkeypatch):
    bridge = bridge_fixture.bridge
    client, _, _ = _route_client(bridge, monkeypatch)
    task_id = next(iter(bridge.tasks))
    csv_response = client.get(f"/api/tasks/{task_id}/frozen-tables/native_compounds")
    assert csv_response.status_code == 200
    assert csv_response.content == b"sample_key,compound_id\nsample,c-native\n"
    assert "native_compounds.csv" in csv_response.headers["content-disposition"]
    json_response = client.get(f"/api/tasks/{task_id}/download/native-conclusions.json.xz")
    assert json_response.status_code == 200
    assert "native_conclusions.json.xz" in json_response.headers["content-disposition"]
    saved = json.loads(lzma.decompress(json_response.content))
    assert len(saved["compound_conclusions"]["compounds"]) == 2
    assert saved["compound_conclusions"]["metadata"] == {"source": "native-exact"}


def test_per_sample_report_routes_call_readonly_renderer_inline_and_attachment(bridge_fixture, monkeypatch):
    bridge = bridge_fixture.bridge
    client, _, _ = _route_client(bridge, monkeypatch)
    task_id = next(iter(bridge.tasks))
    inline = client.get(f"/api/tasks/{task_id}/report")
    assert inline.status_code == 200
    assert "完整原生结论及QC记录 · 2条" in inline.text
    assert "独立MS1靶向特征 · 2条" in inline.text
    assert "完整原生结论与证据JSON.XZ" in inline.text
    assert "content-disposition" not in inline.headers
    download = client.get(f"/api/tasks/{task_id}/download/html")
    assert download.status_code == 200
    assert "attachment" in download.headers["content-disposition"]
    assert download.text == inline.text


def test_unknown_frozen_route_never_falls_through_to_materializer(bridge_fixture, monkeypatch):
    bridge = bridge_fixture.bridge
    client, _, materialized = _route_client(bridge, monkeypatch)
    task_id = next(iter(bridge.tasks))
    response = client.get(f"/api/tasks/{task_id}/dangerous-materializer")
    assert response.status_code == 404
    assert materialized["count"] == 0


def test_invalid_scientific_view_values_are_rejected_not_silently_remapped(bridge_fixture, monkeypatch):
    bridge = bridge_fixture.bridge
    client, _, _ = _route_client(bridge, monkeypatch)
    task_id = next(iter(bridge.tasks))
    cases = [
        (f"/api/tasks/{task_id}/results?blank_view=typo&projection=summary", 422),
        (f"/api/tasks/{task_id}/results?match_state=typo&projection=summary", 422),
        (f"/api/tasks/{task_id}/compounds?report_scope=typo", 422),
        (f"/api/tasks/{task_id}/compounds?result_view=typo", 422),
        (f"/api/tasks/{task_id}/compounds?result_view=xbt_linked&report_scope=plant_candidates", 422),
        (f"/api/tasks/{task_id}/download/csv?result_layer=typo", 422),
    ]
    for path, expected in cases:
        assert client.get(path).status_code == expected
