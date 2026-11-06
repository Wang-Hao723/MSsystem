import copy
import hashlib
import io
import json
import lzma
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from tools import build_frozen_native_delivery as builder
from web import frozen_native_delivery as delivery


PROTOCOL = "1d09e21e37dad4c856a70c6b09a8ccecf611aa5cadc40567f8acd82767fbe779"


def _sha(value):
    return hashlib.sha256(value).hexdigest()


@pytest.fixture
def frozen_inputs(tmp_path):
    remote = tmp_path / "remote"
    (remote / "results").mkdir(parents=True)
    (remote / "conclusions").mkdir()
    samples = []
    entries = {}
    tasks = {}
    for index in range(149):
        sample_key = f"sample-{index:03d}__{index:016x}"
        source_sha = f"{index + 1000:064x}"
        result_path = remote / "results" / f"{sample_key}.json.xz"
        result_document = {"sample_key": sample_key, "source": {"sha256": source_sha},
                           "protocol_sha256": PROTOCOL, "native_ms2_results": []}
        result_raw = lzma.compress(json.dumps(result_document).encode())
        result_path.write_bytes(result_raw)
        conclusion_path = remote / "conclusions" / f"{sample_key}.json.xz"
        compound = {
            "compound_id": f"compound:{index}", "result_id": f"result:{index}",
            "feature_id": f"feature:{index}", "record_unit": "ion_feature",
            "record_type": "unannotated_ion_feature", "identity_key": f"identity:{index}",
            "identity_scope": "unannotated_feature", "identity_source": "saved",
            "display_name": "Unannotated ion feature", "formula": "",
            "lead_candidate": {}, "candidate_members": [],
            "conclusion": {"key": "unannotated", "counted_as_sample_compound": False,
                           "counted_as_plant_candidate": False},
            "plant_report": {"excluded_from_plant_candidate_table": False},
            "representative": {"result_index": index, "spectrum_id": f"scan-{index}"},
        }
        conclusion_document = {
            "sample_key": sample_key,
            "native_result_binding": {"compressed_sha256": _sha(result_raw),
                "compressed_size": len(result_raw), "remote_path": str(result_path),
                "source_sha256": source_sha, "native_protocol_sha256": PROTOCOL},
            "compound_conclusions": {"metadata": {"source": "saved"}, "compounds": [compound]},
            "report_views": {"xbt_named_candidates": {"metadata": {}, "compounds": []}},
        }
        conclusion_raw = lzma.compress(json.dumps(conclusion_document).encode())
        conclusion_path.write_bytes(conclusion_raw)
        task_id = hashlib.sha256(("xbt-20261021-ms1-msms:" + sample_key).encode()).hexdigest()[:12]
        summary = {"compound_result_records": 1, "parsed_searched_ms2_spectra": 2,
                   "matched_ms2_spectra": 1, "ms1_feature_records": 1,
                   "source_sample": {"sample_id": f"S{index:03d}", "polarity": "POS"}}
        entry = {"sample_key": sample_key, "source_sha256": source_sha,
                 "native_protocol_sha256": PROTOCOL,
                 "results": {"path": str(result_path), "size": len(result_raw), "sha256": _sha(result_raw)},
                 "conclusions": {"path": str(conclusion_path), "size": len(conclusion_raw),
                                 "sha256": _sha(conclusion_raw)}}
        entries[task_id] = entry
        tasks[task_id] = {"completed_batch": {"summary": summary}}
        samples.append({**entry, "task_id": task_id, "summary": summary})
    original_index = tmp_path / "FINAL_DELIVERY_INDEX.json"
    original_index.write_text('{"status":"complete"}')
    payload = {"original_delivery_index_path": str(original_index),
               "original_delivery_index_sha256": _sha(original_index.read_bytes()), "samples": samples}
    reader_source = Path(delivery.__file__).with_name("frozen_batch_reader.py").read_text().replace(
        'ROOT = Path("/data/Mssystem/xbt_real_batch_20261021_ms2")', f"ROOT = Path({str(remote)!r})")
    report = tmp_path / "legacy.html"; report.write_text("legacy report")
    packaging = tmp_path / "packaging.json"; packaging.write_text("{}")
    verification = tmp_path / "verification.json"; verification.write_text("{}")
    legacy = tmp_path / "legacy.csv.xz"; legacy.write_bytes(lzma.compress(b"legacy\n"))
    service = SimpleNamespace(
        cache_root=tmp_path / "cache",
        paths={"index": original_index, "report": report, "packaging": packaging,
               "verification": verification}, remote_host="unused",
        exports={"legacy": {"filename": "legacy.csv.xz"}},
        ensure_exports=lambda names: {"legacy": legacy},
    )
    bridge = SimpleNamespace(entries=entries, tasks=tasks, service=service)
    return SimpleNamespace(payload=payload, reader_source=reader_source, bridge=bridge)


def test_builder_exports_all_149_rows_and_binds_every_frozen_input(frozen_inputs, tmp_path, monkeypatch):
    monkeypatch.setattr(builder.shutil, "disk_usage", lambda path: SimpleNamespace(free=100 * 1024**3))
    output = tmp_path / "delivery"
    index = builder._build_delivery(frozen_inputs.payload, output, frozen_inputs.reader_source, workers=2)
    assert index["sample_count"] == 149
    assert index["native_conclusion_rows"] == 149
    assert index["record_type_counts"] == {"unannotated_ion_feature": 149}
    assert index["conclusion_status_counts"] == {"unannotated": 149}
    assert len(index["samples"]) == 149
    assert index["standard_confirmed_count"] == 0
    with lzma.open(output / "all_native_conclusions.csv.xz", "rt", encoding="utf-8-sig") as handle:
        lines = handle.readlines()
    assert len(lines) == 150
    assert "original_row_json" not in lines[0]
    assert "representative_json" not in lines[0]
    assert "XBT数量是附加筛选记录数，不是命名化合物数" in (
        output / "FULL_NATIVE_RESULTS_149.html").read_text()
    saved = json.loads((output / "FULL_NATIVE_DELIVERY_INDEX.json").read_text())
    assert saved["semantic_sha256"] == delivery._semantic_sha256(
        {key: value for key, value in saved.items() if key != "semantic_sha256"})


def _publishable_index(built, bridge, tmp_path):
    value = copy.deepcopy(built)
    csv_path = tmp_path / "all_native_conclusions.csv.xz"
    csv_rows = ["sample_key,compound_id"] + [
        f"{entry['sample_key']},compound:{index}"
        for index, entry in enumerate(bridge.entries.values())
    ]
    csv_path.write_bytes(lzma.compress(("\n".join(csv_rows) + "\n").encode()))
    report_path = tmp_path / "FULL_NATIVE_RESULTS_149.html"
    report_path.write_text("<html>saved native report</html>")
    value["artifacts"] = {
        "native_compounds": {"path": str(delivery.REMOTE_FILES["native_compounds"]),
                             "size": csv_path.stat().st_size, "sha256": _sha(csv_path.read_bytes()), "rows": 149},
        "report": {"path": str(delivery.REMOTE_FILES["report"]),
                   "size": report_path.stat().st_size, "sha256": _sha(report_path.read_bytes())},
    }
    value["semantic_sha256"] = delivery._semantic_sha256(
        {key: item for key, item in value.items() if key != "semantic_sha256"})
    return value, csv_path, report_path


def test_index_validation_and_routes_reject_binding_tamper(frozen_inputs, tmp_path, monkeypatch):
    samples = []
    for task_id, entry in frozen_inputs.bridge.entries.items():
        samples.append({"sample_key": entry["sample_key"], "task_id": task_id,
                        "source_sha256": entry["source_sha256"],
                        "native_result": entry["results"], "conclusion_result": entry["conclusions"],
                        "native_conclusion_rows": 1})
    built = {"schema_version": delivery.SCHEMA, "status": "completed", "sample_count": 149,
             "native_conclusion_rows": 149, "record_type_counts": {"unannotated_ion_feature": 149},
             "conclusion_status_counts": {"unannotated": 149}, "counting_unit": "record_occurrences",
             "standard_confirmed_count": 0,
             "original_delivery_index_sha256": _sha(frozen_inputs.bridge.service.paths["index"].read_bytes()),
             "samples": samples, "artifacts": {}}
    index, csv_path, report_path = _publishable_index(built, frozen_inputs.bridge, tmp_path)
    raw_index = json.dumps(index).encode()
    monkeypatch.setattr(delivery.FrozenNativeDelivery, "_ssh_bytes", lambda self, path, timeout=180: raw_index)
    monkeypatch.setattr(delivery.FrozenNativeDelivery, "_cache_file",
                        lambda self, key: csv_path if key == "native_compounds" else report_path)
    native_index_path = tmp_path / "FULL_NATIVE_DELIVERY_INDEX.json"
    native_index_path.write_bytes(raw_index)
    monkeypatch.setattr(delivery.FrozenNativeDelivery, "index_artifact", lambda self: native_index_path)
    app = FastAPI()
    delivery.register_frozen_native_delivery_routes(app, frozen_inputs.bridge)
    client = TestClient(app)
    base = "/api/completed-batches/xbt-20261021-ms1-msms"
    assert client.get(base + "/native-delivery").json()["native_conclusion_rows"] == 149
    assert client.get(base + "/native-report").text == "<html>saved native report</html>"
    assert lzma.decompress(client.get(base + "/native-exports/native_compounds").content).startswith(b"sample_key")
    assert client.get(base + "/native-exports/native_compounds?decompress=true").content.startswith(b"sample_key")
    assert client.get(base + "/native-index").json()["sample_count"] == 149
    table = client.get(base + "/native-table?offset=1&limit=2").json()
    assert table["total"] == 149 and len(table["rows"]) == 2 and table["has_more"] is True
    sample_key = samples[0]["sample_key"]
    sample_table = client.get(base + "/native-table", params={"sample_key": sample_key}).json()
    assert sample_table["total"] == 1
    assert sample_table["rows"][0]["sample_key"] == sample_key
    assert client.get(base + "/native-table?limit=201").status_code == 422
    archive = client.get(base + "/native-download-all")
    assert archive.status_code == 200
    with zipfile.ZipFile(io.BytesIO(archive.content)) as bundle:
        assert set(bundle.namelist()) == {
            "legacy/legacy.csv.xz", "legacy/既有XBT附加筛选摘要.html",
            "legacy/FINAL_DELIVERY_INDEX.json", "legacy/FINAL_PACKAGING_RECEIPT.json",
            "legacy/FINAL_DELIVERY_VERIFICATION.json",
            "all_native_conclusions.csv.xz", "FULL_NATIVE_RESULTS_149.html",
            "FULL_NATIVE_DELIVERY_INDEX.json",
        }
        assert lzma.decompress(bundle.read("all_native_conclusions.csv.xz")).startswith(b"sample_key")

    damaged = copy.deepcopy(index)
    damaged["samples"][0]["source_sha256"] = "0" * 64
    damaged["semantic_sha256"] = delivery._semantic_sha256(
        {key: value for key, value in damaged.items() if key != "semantic_sha256"})
    checker = delivery.FrozenNativeDelivery(frozen_inputs.bridge)
    with pytest.raises(RuntimeError, match="输入绑定"):
        checker._validate_index(damaged)
