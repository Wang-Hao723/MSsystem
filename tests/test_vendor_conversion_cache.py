from __future__ import annotations

import hashlib
import json
import os
import threading
from copy import deepcopy
from pathlib import Path

import pytest

from utils import vendor_conversion_cache as cache_module
from utils.vendor_conversion_cache import (
    CACHE_MANIFEST_SCHEMA_VERSION,
    CacheEntryInvalidError,
    ConversionLedger,
    LedgerStateError,
    SourceMaterialChangedError,
    VendorConversionCache,
    build_conversion_cache_key,
)
from utils.vendor_raw_conversion import inspect_agilent_d, inspect_thermo_raw


_CONVERTER_IDENTITY = {
    "backend": "docker",
    "image": "proteowizard:3.0.25142",
    "image_identity": "sha256:image-one",
    "wine_runner": "wine64_anyuser",
}
_VALIDATION_RULES_VERSION = "validate-mzml-v1"


def _signature(path: Path) -> dict[str, int | str]:
    content = path.read_bytes()
    metadata = path.stat()
    return {
        "size_bytes": metadata.st_size,
        "mtime_ns": metadata.st_mtime_ns,
        "sha256": hashlib.sha256(content).hexdigest(),
        "sha1_for_mzml_provenance": hashlib.sha1(
            content, usedforsecurity=False
        ).hexdigest(),
    }


def _make_raw_report(root: Path) -> dict:
    raw = root / "sample.d"
    acquisition = raw / "AcqData"
    acquisition.mkdir(parents=True)
    required = {
        "AcqData/MSScan.bin": b"scan-data",
        "AcqData/MSPeak.bin": b"peak-data",
        "AcqData/AcqMethod.xml": b"<method StorageMode='Both'/>",
    }
    for relative, content in required.items():
        (raw / relative).write_bytes(content)
    (raw / "AcqData" / "noncritical.bin").write_bytes(b"auxiliary")
    all_files = [path for path in raw.rglob("*") if path.is_file()]
    return {
        "status": "ok",
        "source_kind": "agilent_d",
        "raw_name": raw.name,
        "source_path": str(raw),
        "file_count": len(all_files),
        "total_bytes": sum(path.stat().st_size for path in all_files),
        "required_file_signatures": {
            relative: _signature(raw / relative) for relative in required
        },
    }


def _key(
    report: dict,
    *,
    converter_identity: dict | None = None,
    peak_picking: str = "vendor",
    validation_rules_version: str = "validate-mzml-v1",
):
    return build_conversion_cache_key(
        report,
        converter_identity=converter_identity or _CONVERTER_IDENTITY,
        resolved_peak_picking=peak_picking,
        validation_rules_version=validation_rules_version,
    )


def _conversion_manifest(path: Path, report: dict) -> dict:
    content = path.read_bytes()
    return {
        "status": "completed",
        "output": {
            "status": "valid",
            "size_bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
            "ms1_count": 3,
            "ms2_count": 7,
        },
        "source": deepcopy(report),
        "converter": deepcopy(_CONVERTER_IDENTITY),
        "conversion_policy": {
            "peak_picking_strategy_resolved": "vendor",
            "peak_picking_filter": "peakPicking true 1-",
        },
    }


def test_cache_key_binds_full_inventory_signatures_converter_policy_and_rules(
    tmp_path,
):
    report = _make_raw_report(tmp_path)
    original = _key(report)
    repeated = _key(deepcopy(report))

    assert repeated.digest == original.digest
    source_material = original.material["source"]
    assert source_material["directory_material"]["file_count"] == 4
    assert source_material["directory_material"]["total_bytes"] > 0
    assert len(source_material["directory_material"]["inventory"]) == 4
    assert report["source_path"] not in json.dumps(original.material)

    auxiliary = Path(report["source_path"]) / "AcqData" / "noncritical.bin"
    auxiliary_stat = auxiliary.stat()
    os.utime(
        auxiliary,
        ns=(auxiliary_stat.st_atime_ns, auxiliary_stat.st_mtime_ns + 1_000_000),
    )


    assert _key(report).digest != original.digest

    auxiliary.write_bytes(b"auxiliary-changed")
    assert _key(report).digest != original.digest

    converter_changed = _key(
        report,
        converter_identity={
            "backend": "docker",
            "image_identity": "sha256:image-two",
        },
    )
    assert converter_changed.digest != original.digest
    assert _key(report, peak_picking="preserve").digest != original.digest
    assert (
        _key(report, validation_rules_version="validate-mzml-v2").digest
        != original.digest
    )


def test_ledger_atomic_claim_allows_only_one_concurrent_owner(tmp_path):
    ledger_path = tmp_path / "ledger.json"
    metadata = {"source_path": "/validated/sample.d", "contract": "fixed"}
    barrier = threading.Barrier(2)
    outcomes: list[tuple[str, str]] = []
    outcome_lock = threading.Lock()

    def attempt(worker_id: str) -> None:
        ledger = ConversionLedger(
            ledger_path,
            job_id="atomic-claim-test",
            worker_id=worker_id,
        )
        barrier.wait()
        try:
            item = ledger.claim(
                "selection-token",
                metadata=metadata,
                details={"worker": worker_id},
            )
        except LedgerStateError as exc:
            result = ("rejected", str(exc))
        else:
            result = ("claimed", str(item["claim_token"]))
        with outcome_lock:
            outcomes.append(result)

    threads = [
        threading.Thread(target=attempt, args=(f"worker-{index}",))
        for index in range(2)
    ]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=5)

    assert all(not thread.is_alive() for thread in threads)
    assert [status for status, _ in outcomes].count("claimed") == 1
    assert [status for status, _ in outcomes].count("rejected") == 1
    current = ConversionLedger(
        ledger_path,
        job_id="atomic-claim-test",
    ).get("selection-token")
    assert current is not None
    assert current["state"] == "converting"
    assert current["attempt"] == 1


def test_cache_key_uses_resolved_not_requested_peak_picking_or_task_path(tmp_path):
    report = _make_raw_report(tmp_path)
    common = {
        "converter_identity": {
            "backend": "docker",
            "image_identity": "sha256:image-one",
        },
        "resolved_peak_picking": "vendor",
        "validation_rules_version": "validate-mzml-v1",
    }
    automatic = build_conversion_cache_key(
        report,
        **common,
        conversion_policy={
            "peak_picking_strategy_requested": "auto",
            "peak_picking_strategy_resolved": "vendor",
            "command_audit_without_source_path": ["--outfile", "task-a.mzML"],
            "peak_picking_filter": "peakPicking true 1-",
        },
    )
    explicit = build_conversion_cache_key(
        report,
        **common,
        conversion_policy={
            "peak_picking_strategy_requested": "vendor",
            "peak_picking_strategy_resolved": "vendor",
            "command_audit_without_source_path": ["--outfile", "task-b.mzML"],
            "peak_picking_filter": "peakPicking true 1-",
        },
    )
    assert automatic.digest == explicit.digest

    with pytest.raises(ValueError, match="immutable image ID"):
        build_conversion_cache_key(
            report,
            converter_identity={"backend": "docker", "image": "mutable:latest"},
            resolved_peak_picking="vendor",
            validation_rules_version="validate-mzml-v1",
        )


def test_cache_key_fails_closed_when_a_signed_source_file_changed(tmp_path):
    report = _make_raw_report(tmp_path)
    scan = Path(report["source_path"]) / "AcqData" / "MSScan.bin"
    scan.write_bytes(b"replacement scan data")

    with pytest.raises(SourceMaterialChangedError, match="size changed"):
        _key(report)

    report["required_file_signatures"]["AcqData/MSScan.bin"] = _signature(scan)
    rebuilt = _key(report)
    assert (
        rebuilt.material["source"]["critical_signatures"]["required_files"][
            "AcqData/MSScan.bin"
        ]["sha256"]
        == hashlib.sha256(scan.read_bytes()).hexdigest()
    )


def test_cache_key_accepts_current_agilent_and_thermo_inspection_reports(tmp_path):
    synthetic_agilent = _make_raw_report(tmp_path)
    agilent_report = inspect_agilent_d(synthetic_agilent["source_path"])
    agilent_key = _key(agilent_report)
    assert agilent_key.material["source"]["source_kind"] == "agilent_d"

    thermo_path = tmp_path / "sample.raw"
    thermo_path.write_bytes(b"synthetic thermo payload")
    thermo_report = inspect_thermo_raw(thermo_path, expected_polarities=["NEG"])
    thermo_key = _key(thermo_report)
    assert thermo_key.material["source"]["source_kind"] == "thermo_raw"
    thermo_pos_report = inspect_thermo_raw(thermo_path, expected_polarities=["POS"])
    assert _key(thermo_pos_report).digest != thermo_key.digest


def test_cache_key_rejects_changed_inspected_source_even_if_mtime_is_restored(
    tmp_path,
):
    thermo_path = tmp_path / "stable.raw"
    thermo_path.write_bytes(b"AAAA")
    report = inspect_thermo_raw(thermo_path)
    original = _key(report)
    metadata = thermo_path.stat()
    thermo_path.write_bytes(b"BBBB")
    os.utime(
        thermo_path,
        ns=(metadata.st_atime_ns, metadata.st_mtime_ns),
    )

    with pytest.raises(SourceMaterialChangedError, match="identity changed"):
        _key(report)

    refreshed = inspect_thermo_raw(thermo_path)
    assert _key(refreshed).digest != original.digest


def test_fresh_inspection_avoids_a_second_full_source_hash(monkeypatch, tmp_path):
    raw = tmp_path / "one-pass.raw"
    raw.write_bytes(b"one-pass vendor payload")
    report = inspect_thermo_raw(raw, expected_polarities=["NEG"])

    def unexpected_second_hash(_path):
        raise AssertionError("fresh inspected source must not be re-read")

    monkeypatch.setattr(cache_module, "_sha256_file", unexpected_second_hash)
    cache_key = _key(report)

    assert cache_key.material["source"]["critical_signatures"][
        "source_file"
    ]["sha256"] == hashlib.sha256(raw.read_bytes()).hexdigest()


def test_cache_store_lookup_and_materialize_independent_task_file(
    monkeypatch,
    tmp_path,
):
    report = _make_raw_report(tmp_path)
    cache_key = _key(report)
    converted = tmp_path / "source.mzML"
    converted.write_bytes(b"<mzML>validated vendor output</mzML>")
    conversion_manifest = _conversion_manifest(converted, report)
    cache = VendorConversionCache(tmp_path / "cache")

    stored = cache.store(
        cache_key,
        converted,
        conversion_manifest,
        validation_rules_version=_VALIDATION_RULES_VERSION,
    )
    assert stored.manifest["schema_version"] == CACHE_MANIFEST_SCHEMA_VERSION
    assert cache.lookup(cache_key) is not None

    # Exercise the hardlink pin fallback even on APFS.  The link must never be
    # published as the task output and the final inode must be independent.
    monkeypatch.setattr(cache_module, "_try_apfs_clone", lambda *_: False)
    task_output = tmp_path / "tasks" / "sample.mzML"
    result = cache.materialize(cache_key, task_output)

    assert result["method"] == "hardlink_then_copy"
    assert task_output.read_bytes() == converted.read_bytes()
    assert task_output.stat().st_ino != stored.payload_path.stat().st_ino
    task_output.write_bytes(b"task-owned mutation")
    assert stored.payload_path.read_bytes() == converted.read_bytes()
    assert cache.lookup(cache_key) is not None
    assert not list(cache.manifest_path(cache_key.digest).parent.glob("*.tmp"))
    assert not list(task_output.parent.glob("*.cachelink"))


@pytest.mark.parametrize("corruption", ["status", "size", "sha256"])
def test_cache_hit_rejects_manifest_status_size_or_hash_corruption(
    tmp_path,
    corruption,
):
    report = _make_raw_report(tmp_path)
    cache_key = _key(report)
    converted = tmp_path / "source.mzML"
    converted.write_bytes(b"stable mzML bytes")
    cache = VendorConversionCache(tmp_path / "cache")
    hit = cache.store(
        cache_key,
        converted,
        _conversion_manifest(converted, report),
        validation_rules_version=_VALIDATION_RULES_VERSION,
    )

    if corruption == "status":
        manifest = json.loads(hit.manifest_path.read_text(encoding="utf-8"))
        manifest["status"] = "writing"
        hit.manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    elif corruption == "size":
        os.chmod(hit.payload_path, hit.payload_path.stat().st_mode | 0o200)
        hit.payload_path.write_bytes(b"different length")
    else:
        original_mode = hit.payload_path.stat().st_mode
        os.chmod(hit.payload_path, original_mode | 0o200)
        data = bytearray(hit.payload_path.read_bytes())
        data[0] ^= 1
        hit.payload_path.write_bytes(bytes(data))

    assert cache.lookup(cache_key) is None


@pytest.mark.parametrize(
    "corruption",
    [
        "audit_status",
        "audit_output_status",
        "audit_output_validation_status",
        "audit_output_size",
        "audit_output_sha256",
    ],
)
def test_cache_hit_rejects_conversion_audit_payload_binding_corruption(
    tmp_path,
    corruption,
):
    report = _make_raw_report(tmp_path)
    cache_key = _key(report)
    converted = tmp_path / "source.mzML"
    converted.write_bytes(b"stable mzML bytes")
    cache = VendorConversionCache(tmp_path / "cache")
    hit = cache.store(
        cache_key,
        converted,
        _conversion_manifest(converted, report),
        validation_rules_version=_VALIDATION_RULES_VERSION,
    )
    manifest = json.loads(hit.manifest_path.read_text(encoding="utf-8"))
    audit = manifest["conversion_audit"]

    if corruption == "audit_status":
        audit["status"] = "failed"
    elif corruption == "audit_output_status":
        audit["output_status"] = "invalid"
    elif corruption == "audit_output_validation_status":
        audit["output"]["status"] = "invalid"
    elif corruption == "audit_output_size":
        audit["output"]["size_bytes"] += 1
    else:
        audit["output"]["sha256"] = "0" * 64
    hit.manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    assert cache.lookup(cache_key) is None
    with pytest.raises(CacheEntryInvalidError):
        cache.materialize(cache_key, tmp_path / "task.mzML")


def test_cache_store_rejects_unvalidated_or_mismatched_output(tmp_path):
    report = _make_raw_report(tmp_path)
    cache_key = _key(report)
    converted = tmp_path / "source.mzML"
    converted.write_bytes(b"mzML")
    cache = VendorConversionCache(tmp_path / "cache")
    manifest = _conversion_manifest(converted, report)
    manifest["output"]["status"] = "invalid"
    with pytest.raises(CacheEntryInvalidError, match="must be valid"):
        cache.store(
            cache_key,
            converted,
            manifest,
            validation_rules_version=_VALIDATION_RULES_VERSION,
        )

    manifest = _conversion_manifest(converted, report)
    manifest["output"]["size_bytes"] += 1
    with pytest.raises(CacheEntryInvalidError, match="size disagrees"):
        cache.store(
            cache_key,
            converted,
            manifest,
            validation_rules_version=_VALIDATION_RULES_VERSION,
        )

    manifest = _conversion_manifest(converted, report)
    manifest["output"]["sha256"] = "0" * 64
    with pytest.raises(CacheEntryInvalidError, match="SHA-256 disagrees"):
        cache.store(
            cache_key,
            converted,
            manifest,
            validation_rules_version=_VALIDATION_RULES_VERSION,
        )


@pytest.mark.parametrize("mismatch", ["converter", "policy", "rules"])
def test_cache_store_rejects_manifest_bound_to_a_different_key(tmp_path, mismatch):
    report = _make_raw_report(tmp_path)
    cache_key = _key(report)
    converted = tmp_path / "source.mzML"
    converted.write_bytes(b"mzML")
    manifest = _conversion_manifest(converted, report)
    rules_version = _VALIDATION_RULES_VERSION
    if mismatch == "converter":
        manifest["converter"]["image_identity"] = "sha256:image-two"
    elif mismatch == "policy":
        manifest["conversion_policy"]["peak_picking_strategy_resolved"] = "preserve"
        manifest["conversion_policy"]["peak_picking_filter"] = None
    else:
        rules_version = "validate-mzml-v2"

    with pytest.raises(CacheEntryInvalidError, match="does not match"):
        VendorConversionCache(tmp_path / "cache").store(
            cache_key,
            converted,
            manifest,
            validation_rules_version=rules_version,
        )


def test_materialize_refuses_overwrite(tmp_path):
    report = _make_raw_report(tmp_path)
    cache_key = _key(report)
    converted = tmp_path / "source.mzML"
    converted.write_bytes(b"mzML")
    cache = VendorConversionCache(tmp_path / "cache")
    cache.store(
        cache_key,
        converted,
        _conversion_manifest(converted, report),
        validation_rules_version=_VALIDATION_RULES_VERSION,
    )
    destination = tmp_path / "already.mzML"
    destination.write_bytes(b"keep")

    with pytest.raises(FileExistsError):
        cache.materialize(cache_key, destination)
    assert destination.read_bytes() == b"keep"


def test_materialize_rejects_dangling_destination_symlink(tmp_path):
    report = _make_raw_report(tmp_path)
    cache_key = _key(report)
    converted = tmp_path / "source.mzML"
    converted.write_bytes(b"mzML")
    cache = VendorConversionCache(tmp_path / "cache")
    cache.store(
        cache_key,
        converted,
        _conversion_manifest(converted, report),
        validation_rules_version=_VALIDATION_RULES_VERSION,
    )
    task_dir = tmp_path / "tasks"
    outside = tmp_path / "outside"
    task_dir.mkdir()
    outside.mkdir()
    destination = task_dir / "escape.mzML"
    outside_target = outside / "outside.mzML"
    destination.symlink_to(outside_target)

    with pytest.raises(CacheEntryInvalidError, match="symbolic link"):
        cache.materialize(cache_key, destination)
    assert destination.is_symlink()
    assert not outside_target.exists()


def test_ledger_persists_valid_state_machine_and_explicit_retry(tmp_path):
    path = tmp_path / "jobs" / "job-1.json"
    ledger = ConversionLedger(path, job_id="job-1")
    ledger.queue("sample-a", metadata={"raw_name": "a.d"})
    converting_a = ledger.mark_converting("sample-a")
    assert converting_a["attempt"] == 1
    validated_a = ledger.mark_validated(
        "sample-a",
        converting_a["claim_token"],
        details={"cache_key": "a" * 64, "sha256": "b" * 64},
    )
    ledger.mark_registered(
        "sample-a",
        converting_a["claim_token"],
        validated_a["registration_token"],
        details={"task_id": "task-a"},
    )

    ledger.queue("sample-b")
    converting_b = ledger.mark_converting("sample-b")
    ledger.mark_failed(
        "sample-b", "converter exited 1", claim_token=converting_b["claim_token"]
    )
    retried = ledger.queue("sample-b")
    assert retried["state"] == "queued"
    assert "error" not in retried
    assert ledger.mark_converting("sample-b")["attempt"] == 2

    reopened = ConversionLedger(path, job_id="job-1")
    assert reopened.get("sample-a")["state"] == "registered"
    assert reopened.get("sample-b")["state"] == "converting"
    with pytest.raises(LedgerStateError, match="invalid conversion transition"):
        reopened.mark_converting("sample-a")
    assert not list(path.parent.glob("*.tmp"))


def test_ledger_restart_requeues_only_interrupted_conversions(tmp_path):
    path = tmp_path / "job-2.json"
    ledger = ConversionLedger(path, job_id="job-2")
    for item_id in ("queued", "converting", "validated", "registered", "failed"):
        ledger.queue(item_id)
    converting_claim = ledger.mark_converting("converting")
    validated_claim = ledger.mark_converting("validated")
    ledger.mark_validated(
        "validated",
        validated_claim["claim_token"],
        details={"output": "ready.mzML"},
    )
    registered_claim = ledger.mark_converting("registered")
    registered_validation = ledger.mark_validated(
        "registered", registered_claim["claim_token"]
    )
    ledger.mark_registered(
        "registered",
        registered_claim["claim_token"],
        registered_validation["registration_token"],
    )
    ledger.mark_queued_failed("failed", "selection rejected")

    restarted = ConversionLedger(path, job_id="job-2")
    assert restarted.recover_interrupted() == []
    assert restarted.get("converting")["claim_token"] == converting_claim["claim_token"]
    assert restarted.recover_interrupted(force=True) == ["converting"]
    plan = restarted.recovery_plan()
    assert plan == {
        "converting": [],
        "failed": ["failed"],
        "queued": ["converting", "queued"],
        "registered": ["registered"],
        "validated": ["validated"],
    }
    recovered = restarted.get("converting")
    assert recovered["attempt"] == 1
    assert recovered["recovered_from"] == "converting"
    assert recovered["history"][-1]["reason"] == "forced_process_recovery"


def test_ledger_claim_token_fences_old_attempt_after_recovery(tmp_path):
    path = tmp_path / "job-fencing.json"
    ledger = ConversionLedger(path, job_id="job-fencing", worker_id="worker-a")
    ledger.queue("sample")
    first = ledger.mark_converting("sample")
    assert ledger.recover_interrupted() == []
    assert ledger.recover_interrupted(force=True) == ["sample"]
    with pytest.raises(LedgerStateError, match="stale claim token"):
        ledger.mark_failed(
            "sample", "late failure from attempt one", claim_token=first["claim_token"]
        )

    second = ledger.mark_converting("sample", owner_id="worker-b")
    assert second["attempt"] == 2
    assert second["claim_token"] != first["claim_token"]
    with pytest.raises(LedgerStateError, match="stale or missing claim token"):
        ledger.mark_validated("sample", first["claim_token"])

    renewed = ledger.renew_lease("sample", second["claim_token"])
    assert renewed["claim_token"] == second["claim_token"]
    validated = ledger.mark_validated("sample", second["claim_token"])
    with pytest.raises(LedgerStateError, match="registration token"):
        ledger.mark_registered(
            "sample", second["claim_token"], "stale-registration-token"
        )
    registered = ledger.mark_registered(
        "sample",
        second["claim_token"],
        validated["registration_token"],
    )
    assert registered["state"] == "registered"


def test_ledger_rejects_corrupt_or_wrong_job_document(tmp_path):
    path = tmp_path / "ledger.json"
    path.write_text("{not-json", encoding="utf-8")
    with pytest.raises(LedgerStateError, match="unreadable"):
        ConversionLedger(path, job_id="job").snapshot()

    path.write_text(
        json.dumps(
            {
                "schema_version": cache_module.LEDGER_SCHEMA_VERSION,
                "job_id": "another-job",
                "items": {},
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(LedgerStateError, match="job_id"):
        ConversionLedger(path, job_id="job").snapshot()
