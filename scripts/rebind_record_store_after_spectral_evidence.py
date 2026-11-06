#!/usr/bin/env python3
"""Audit and rebind accepted record-store metadata after observer-only evidence.

Preview is the default.  ``--apply`` is a maintenance-only transaction: it
requires port 8516 to be closed, takes an exclusive lock, creates a new backup,
stages all four metadata files, and publishes them with guarded rollback.
Cached records, payload arrays, and the source SQLite database are never opened
for writing or copied by this command.
"""

from __future__ import annotations

import argparse
import ast
import copy
import fcntl
import hashlib
import importlib
import json
import os
import shutil
import socket
import sys
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.library_record_store import (  # noqa: E402
    file_sha256,
    processing_identity,
    semantic_sha256,
    validate_acceptance,
)
from runtime.manifest import source_tree_snapshot  # noqa: E402


RECORD_STORE = (
    PROJECT_ROOT / "output/library_record_store/full-20260915-formula-composition-v1"
)
MANIFEST_PATH = RECORD_STORE / "manifest.json"
ACCEPTANCE_PATH = RECORD_STORE / "acceptance.json"
RECORD_CONFIG_PATH = PROJECT_ROOT / "config/library_record_store.json"
SQL_CONFIG_PATH = PROJECT_ROOT / "config/sql_metadata_repairs.json"
METADATA_PATHS = (
    MANIFEST_PATH,
    ACCEPTANCE_PATH,
    RECORD_CONFIG_PATH,
    SQL_CONFIG_PATH,
)
SEARCH_PATH = "matching/search.py"
HELPER_PATH = "matching/spectral_evidence.py"
AUDIT_FIELD = "post_acceptance_spectral_evidence_source_rebind"
HISTORICAL_AUDIT_FIELDS = (
    "post_acceptance_non_normalization_source_rebind",
    "post_acceptance_reliability_source_rebind",
    "post_acceptance_representative_summary_rebind",
    "post_acceptance_target_only_source_rebind",
)
PAYLOAD_NAMES = (
    "records.jsonl",
    "record_offsets.npy",
    "precursor_masses.npy",
    "precursor_record_indices.npy",
    "neutral_masses.npy",
    "neutral_record_indices.npy",
)
PROTOCOL_SCHEMA = "mssystem.spectral-evidence-cache-rebind-protocol/1.0"
PROTOCOL_KEYS = {
    "schema_version",
    "baseline_ref",
    "accepted_source_files_sha256",
    "baseline_search",
    "reviewed_sources",
    "test_receipt",
    "review_receipt",
    "audit_scope",
}
SOURCE_SUFFIXES = {
    ".css", ".html", ".js", ".json", ".jsonl", ".py", ".sh", ".txt",
    ".yaml", ".yml",
}


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"expected JSON object: {path}")
    return payload


def _json_bytes(payload: dict[str, Any]) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def _write_bytes_atomic(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    _write_bytes_atomic(path, _json_bytes(payload))


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _refresh_semantic(payload: dict[str, Any]) -> None:
    payload["semantic_sha256"] = semantic_sha256(
        {key: value for key, value in payload.items() if key != "semantic_sha256"}
    )


def _require_exact_keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    if set(value) != expected:
        raise RuntimeError(
            f"{label} keys differ: missing={sorted(expected - set(value))}, "
            f"extra={sorted(set(value) - expected)}"
        )


def _source_pin(
    value: Any, *, path: str, label: str, include_path: bool = True
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} must be an object")
    expected = {"sha256", "size_bytes"} | ({"path"} if include_path else set())
    _require_exact_keys(value, expected, label)
    if include_path and value["path"] != path:
        raise RuntimeError(f"{label} path differs")
    if not isinstance(value["sha256"], str) or len(value["sha256"]) != 64:
        raise RuntimeError(f"{label} sha256 is invalid")
    if not isinstance(value["size_bytes"], int) or value["size_bytes"] < 0:
        raise RuntimeError(f"{label} size is invalid")
    return value


def _receipt_pin(value: Any, *, label: str) -> dict[str, str]:
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} must be an object")
    _require_exact_keys(value, {"path", "sha256"}, label)
    if not isinstance(value["path"], str) or not value["path"]:
        raise RuntimeError(f"{label} path is invalid")
    if not isinstance(value["sha256"], str) or len(value["sha256"]) != 64:
        raise RuntimeError(f"{label} sha256 is invalid")
    return value


def _resolve_receipt(path_text: str) -> Path:
    path = Path(path_text)
    return path.resolve() if path.is_absolute() else (PROJECT_ROOT / path).resolve()


def _load_protocol(path: Path, expected_sha256: str) -> tuple[dict[str, Any], dict[str, Any]]:
    if len(expected_sha256) != 64 or any(c not in "0123456789abcdef" for c in expected_sha256):
        raise RuntimeError("--protocol-sha256 must be 64 lowercase hexadecimal characters")
    actual = file_sha256(path)
    if actual != expected_sha256:
        raise RuntimeError(f"protocol hash differs: {actual}")
    protocol = _read_json(path)
    _require_exact_keys(protocol, PROTOCOL_KEYS, "protocol")
    if protocol["schema_version"] != PROTOCOL_SCHEMA:
        raise RuntimeError("protocol schema differs")
    if not isinstance(protocol["baseline_ref"], str) or not protocol["baseline_ref"]:
        raise RuntimeError("protocol baseline_ref is invalid")
    if (
        not isinstance(protocol["accepted_source_files_sha256"], str)
        or len(protocol["accepted_source_files_sha256"]) != 64
    ):
        raise RuntimeError("protocol accepted source hash is invalid")
    baseline = _source_pin(
        protocol["baseline_search"], path=SEARCH_PATH, label="baseline_search"
    )
    reviewed = protocol["reviewed_sources"]
    if not isinstance(reviewed, dict):
        raise RuntimeError("reviewed_sources must be an object")
    _require_exact_keys(reviewed, {SEARCH_PATH, HELPER_PATH}, "reviewed_sources")
    _source_pin(
        reviewed[SEARCH_PATH], path=SEARCH_PATH, label="reviewed search",
        include_path=False,
    )
    _source_pin(
        reviewed[HELPER_PATH], path=HELPER_PATH, label="reviewed helper",
        include_path=False,
    )
    test_pin = _receipt_pin(protocol["test_receipt"], label="test_receipt")
    review_pin = _receipt_pin(protocol["review_receipt"], label="review_receipt")
    scope = protocol["audit_scope"]
    if not isinstance(scope, dict):
        raise RuntimeError("audit_scope must be an object")
    _require_exact_keys(scope, {"successful_result_fields", "fresh_fault_isolation"}, "audit_scope")
    if scope != {"successful_result_fields": ["spectral_evidence"], "fresh_fault_isolation": True}:
        raise RuntimeError("audit_scope exceeds or differs from the reviewed observer scope")
    receipts: dict[str, Any] = {}
    for name, pin in (("test", test_pin), ("review", review_pin)):
        receipt_path = _resolve_receipt(pin["path"])
        if not receipt_path.is_file() or file_sha256(receipt_path) != pin["sha256"]:
            raise RuntimeError(f"{name} receipt is missing or altered")
        receipts[name] = {
            "path": str(receipt_path), "sha256": pin["sha256"],
            "payload": _read_json(receipt_path),
        }
    test = receipts["test"]["payload"]
    if (
        test.get("state") != "completed" or test.get("exit_code") != 0
        or test.get("passed") != 61 or test.get("skipped") != 3
        or test.get("source_pins") != reviewed
    ):
        raise RuntimeError("test receipt does not prove the pinned 61/3 successful run")
    review = receipts["review"]["payload"]
    if (
        review.get("state") != "completed"
        or review.get("verdict") != "cache_only_reuse_permitted_with_strict_normalization_equivalence_proof"
        or review.get("source_pins") != reviewed
        or review.get("accepted_source_files_sha256") != protocol["accepted_source_files_sha256"]
    ):
        raise RuntimeError("review receipt does not authorize the pinned cache-only proof")
    del baseline
    return protocol, receipts


def _processing_files() -> list[dict[str, Any]]:
    return [
        dict(item)
        for item in source_tree_snapshot(PROJECT_ROOT)["files"]
        if item["path"].startswith(("matching/", "utils/", "mssystem/"))
        or item["path"] in {"config/settings.py", "config/retrieval_policy.json"}
    ]


def _offset(lines: list[str], line: int, column: int) -> int:
    return sum(len(value) for value in lines[: line - 1]) + column


def _search_reverse_patch(source: bytes, baseline: dict[str, Any]) -> tuple[bytes, dict[str, Any]]:
    text = source.decode("utf-8")
    tree = ast.parse(text)
    searches = [
        node for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == "search"
    ]
    capture_blocks = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.If) and isinstance(node.test, ast.Name)
        and node.test.id == "capture_peak_evidence"
    ]
    if len(searches) != 1 or len(capture_blocks) != 1:
        raise RuntimeError("expected exactly one search definition and capture block")
    search, capture = searches[0], capture_blocks[0]
    if not search.args.args or search.args.args[-1].arg != "capture_peak_evidence":
        raise RuntimeError("capture_peak_evidence is not the trailing search argument")
    if not search.args.defaults or not (
        isinstance(search.args.defaults[-1], ast.Constant)
        and search.args.defaults[-1].value is True
    ):
        raise RuntimeError("capture_peak_evidence default is not True")
    parent: dict[ast.AST, ast.AST] = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parent[child] = node
    owner = parent.get(capture)
    if not (
        isinstance(owner, ast.If)
        and isinstance(owner.test, ast.Compare)
        and ast.unparse(owner.test) == "best is not None"
        and capture in owner.body
    ):
        raise RuntimeError("capture block is not post-selection under best is not None")
    calls = [
        node for node in ast.walk(capture)
        if isinstance(node, ast.Call)
        and ((isinstance(node.func, ast.Name) and node.func.id == "attach_spectral_evidence"))
    ]
    tries = [node for node in capture.body if isinstance(node, ast.Try)]
    handlers = tries[0].handlers if len(tries) == 1 else []
    if len(calls) != 1 or len(handlers) != 1 or ast.unparse(handlers[0].type) != "Exception":
        raise RuntimeError("capture block lacks the reviewed fresh fault isolation")

    lines = text.splitlines(keepends=True)
    argument = search.args.args[-1]
    default = search.args.defaults[-1]
    argument_start = _offset(lines, argument.lineno, argument.col_offset)
    comma = text.rfind(",", 0, argument_start)
    argument_end = _offset(lines, default.end_lineno, default.end_col_offset)
    if comma < 0 or "capture_peak_evidence" not in text[comma:argument_end]:
        raise RuntimeError("could not isolate trailing capture argument")
    block_start = _offset(lines, capture.lineno, 0)
    block_end = _offset(lines, capture.end_lineno + 1, 0)
    without_argument = text[:comma] + text[argument_end:]
    removed_width = argument_end - comma
    reconstructed = (
        without_argument[: block_start - removed_width]
        + without_argument[block_end - removed_width :]
    ).encode("utf-8")
    if (
        len(reconstructed) != baseline["size_bytes"]
        or _sha256_bytes(reconstructed) != baseline["sha256"]
    ):
        raise RuntimeError("AST-scoped reverse patch did not reproduce baseline search bytes")

    transformed = copy.deepcopy(tree)
    transformed_search = next(
        node for node in ast.walk(transformed)
        if isinstance(node, ast.FunctionDef) and node.name == "search"
    )
    if transformed_search.args.args[-1].arg != "capture_peak_evidence":
        raise RuntimeError("transformed search argument changed concurrently")
    transformed_search.args.args.pop()
    transformed_search.args.defaults.pop()
    transformed_capture = next(
        node for node in ast.walk(transformed)
        if isinstance(node, ast.If) and isinstance(node.test, ast.Name)
        and node.test.id == "capture_peak_evidence"
    )
    transformed_owner = next(
        node for node in ast.walk(transformed)
        if isinstance(node, ast.If) and transformed_capture in node.body
    )
    transformed_owner.body.remove(transformed_capture)
    ast.fix_missing_locations(transformed)
    if ast.dump(transformed, include_attributes=False) != ast.dump(
        ast.parse(reconstructed), include_attributes=False
    ):
        raise RuntimeError("whole-module transformed AST differs from baseline AST")

    current_normalizer = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "iter_normalized_library_records"
    )
    baseline_normalizer = next(
        node for node in ast.walk(ast.parse(reconstructed))
        if isinstance(node, ast.FunctionDef) and node.name == "iter_normalized_library_records"
    )
    normalizer_hash = _sha256_bytes(
        ast.dump(current_normalizer, include_attributes=False).encode("utf-8")
    )
    if ast.dump(current_normalizer, include_attributes=False) != ast.dump(
        baseline_normalizer, include_attributes=False
    ):
        raise RuntimeError("iter_normalized_library_records differs from baseline")
    return reconstructed, {
        "removed_trailing_argument": "capture_peak_evidence=True",
        "removed_post_selection_capture_blocks": 1,
        "whole_module_ast_matches_baseline": True,
        "iter_normalized_library_records_ast_sha256": normalizer_hash,
        "iter_normalized_library_records_unchanged": True,
        "capture_after_best_selection": True,
        "fresh_exception_isolation": True,
    }


def _helper_proof(source: bytes) -> dict[str, Any]:
    tree = ast.parse(source.decode("utf-8"))
    allowed = (ast.Expr, ast.Import, ast.ImportFrom, ast.FunctionDef, ast.AsyncFunctionDef)
    if any(not isinstance(node, allowed) for node in tree.body):
        raise RuntimeError("spectral helper has executable top-level statements")
    expressions = [node for node in tree.body if isinstance(node, ast.Expr)]
    if any(not isinstance(node.value, ast.Constant) or not isinstance(node.value.value, str) for node in expressions):
        raise RuntimeError("spectral helper has a non-docstring top-level expression")
    attach = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "attach_spectral_evidence"
    ]
    if len(attach) != 1:
        raise RuntimeError("spectral helper attach function is missing or duplicated")
    stored_subscripts = [
        node for node in ast.walk(attach[0])
        if isinstance(node, ast.Subscript) and isinstance(node.ctx, ast.Store)
    ]
    def allowed_store(node: ast.Subscript) -> bool:
        if isinstance(node.value, ast.Name) and node.value.id in {"by_reference", "captured"}:
            return True
        return (
            isinstance(node.value, ast.Name) and node.value.id == "candidate"
            and isinstance(node.slice, ast.Constant)
            and node.slice.value == "spectral_evidence"
        )

    if not stored_subscripts or any(not allowed_store(node) for node in stored_subscripts):
        raise RuntimeError("spectral helper mutates a reviewed result field other than spectral_evidence")
    return {
        "top_level_only_docstring_imports_and_functions": True,
        "successful_result_fields": ["spectral_evidence"],
        "ranking_influence_declared": False,
    }


def _source_proof(protocol: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    reviewed = protocol["reviewed_sources"]
    source_bytes: dict[str, bytes] = {}
    for relative in (SEARCH_PATH, HELPER_PATH):
        content = (PROJECT_ROOT / relative).read_bytes()
        pin = reviewed[relative]
        if len(content) != pin["size_bytes"] or _sha256_bytes(content) != pin["sha256"]:
            raise RuntimeError(f"reviewed source changed: {relative}")
        source_bytes[relative] = content
    baseline_search, ast_proof = _search_reverse_patch(
        source_bytes[SEARCH_PATH], protocol["baseline_search"]
    )
    helper_proof = _helper_proof(source_bytes[HELPER_PATH])
    files = _processing_files()
    paths = [item["path"] for item in files]
    if paths.count(SEARCH_PATH) != 1 or paths.count(HELPER_PATH) != 1:
        raise RuntimeError("processing source inventory lacks the reviewed files")
    reconstructed = []
    for item in files:
        if item["path"] == SEARCH_PATH:
            reconstructed.append({
                "path": SEARCH_PATH,
                "size_bytes": len(baseline_search),
                "sha256": _sha256_bytes(baseline_search),
            })
        elif item["path"] != HELPER_PATH:
            reconstructed.append(item)
    reconstructed_hash = semantic_sha256(reconstructed)
    if reconstructed_hash != protocol["accepted_source_files_sha256"]:
        raise RuntimeError("covered source reconstruction differs from accepted identity")
    current_hash = semantic_sha256(files)
    current_identity = processing_identity()
    if current_identity.get("source_files_sha256") != current_hash:
        raise RuntimeError("manual proof inventory differs from normal processing_identity")
    return current_identity, {
        "baseline_ref": protocol["baseline_ref"],
        "accepted_source_files_sha256": protocol["accepted_source_files_sha256"],
        "current_source_files_sha256": current_hash,
        "reviewed_sources": reviewed,
        "reconstructed_baseline_search": protocol["baseline_search"],
        "removed_new_helper_from_baseline_snapshot": True,
        "other_covered_source_files_changed": 0,
        **ast_proof,
        **helper_proof,
    }


def _verify_semantic(payload: dict[str, Any], label: str) -> None:
    expected = semantic_sha256({k: v for k, v in payload.items() if k != "semantic_sha256"})
    if payload.get("semantic_sha256") != expected:
        raise RuntimeError(f"{label} semantic hash differs")


def _verify_historical(payloads: list[dict[str, Any]]) -> dict[str, Any]:
    inherited: dict[str, Any] = {}
    for field in HISTORICAL_AUDIT_FIELDS:
        values = [payload.get(field) for payload in payloads]
        if any(value is None for value in values) or any(value != values[0] for value in values[1:]):
            raise RuntimeError(f"historical audit field missing or inconsistent: {field}")
        inherited[field] = values[0]
    return inherited


def _verify_preconditions(
    manifest: dict[str, Any], acceptance: dict[str, Any], record_config: dict[str, Any],
    sql_config: dict[str, Any], current_identity: dict[str, Any], protocol: dict[str, Any],
) -> dict[str, Any]:
    _verify_semantic(manifest, "manifest")
    _verify_semantic(acceptance, "acceptance")
    _verify_semantic(sql_config, "SQL metadata")
    validate_acceptance(acceptance, manifest)
    accepted_identity = manifest.get("processing_identity")
    if not isinstance(accepted_identity, dict):
        raise RuntimeError("manifest processing identity is missing")
    if accepted_identity.get("source_files_sha256") != protocol["accepted_source_files_sha256"]:
        raise RuntimeError("manifest is not bound to the accepted source snapshot")
    if {k: v for k, v in current_identity.items() if k != "source_files_sha256"} != {
        k: v for k, v in accepted_identity.items() if k != "source_files_sha256"
    }:
        raise RuntimeError("numeric runtime or normalization options changed")
    if sql_config.get("processing_identity") != accepted_identity:
        raise RuntimeError("SQL metadata identity differs from accepted manifest")
    if sql_config.get("cache_manifest_semantic_sha256") != manifest["semantic_sha256"]:
        raise RuntimeError("SQL metadata manifest binding differs")
    if record_config.get("accepted_cache_manifest_semantic_sha256") != manifest["semantic_sha256"]:
        raise RuntimeError("record-store manifest binding differs")
    if record_config.get("acceptance_receipt_file_sha256") != file_sha256(ACCEPTANCE_PATH):
        raise RuntimeError("record-store acceptance file binding differs")
    if manifest.get("record_count") != 3_153_872 or manifest.get("complete_library") is not True:
        raise RuntimeError("accepted cache count or completeness differs")
    files = manifest.get("files")
    if not isinstance(files, dict) or set(files) != set(PAYLOAD_NAMES):
        raise RuntimeError("manifest does not retain exactly six accepted payload bindings")
    inherited = _verify_historical([manifest, acceptance, record_config, sql_config])
    return {
        "accepted_identity": accepted_identity,
        "historical_audits": inherited,
        "acceptance_results": {
            key: acceptance[key] for key in (
                "full_record_roundtrip_sha256", "verified_record_count",
                "candidate_order_cases", "candidate_order_mismatches",
                "paired_scoring_queries", "score_mismatches",
                "shared_database_unchanged",
            )
        },
    }


def _only_changes(original: dict[str, Any], proposed: dict[str, Any], allowed: set[str], label: str) -> None:
    keys = set(original) | set(proposed)
    changed = {key for key in keys if original.get(key) != proposed.get(key)}
    if not changed <= allowed:
        raise RuntimeError(f"{label} changed forbidden fields: {sorted(changed - allowed)}")


def _build_proposed(
    manifest: dict[str, Any], acceptance: dict[str, Any], record_config: dict[str, Any],
    sql_config: dict[str, Any], current_identity: dict[str, Any], audit: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    new_manifest = copy.deepcopy(manifest)
    new_manifest["processing_identity"] = current_identity
    new_manifest[AUDIT_FIELD] = audit
    _refresh_semantic(new_manifest)
    new_acceptance = copy.deepcopy(acceptance)
    new_acceptance["manifest_semantic_sha256"] = new_manifest["semantic_sha256"]
    new_acceptance[AUDIT_FIELD] = audit
    _refresh_semantic(new_acceptance)
    acceptance_bytes = _json_bytes(new_acceptance)
    new_record = copy.deepcopy(record_config)
    new_record["accepted_cache_manifest_semantic_sha256"] = new_manifest["semantic_sha256"]
    new_record["acceptance_receipt_file_sha256"] = _sha256_bytes(acceptance_bytes)
    new_record[AUDIT_FIELD] = audit
    new_sql = copy.deepcopy(sql_config)
    new_sql["processing_identity"] = current_identity
    new_sql["cache_manifest_semantic_sha256"] = new_manifest["semantic_sha256"]
    new_sql[AUDIT_FIELD] = audit
    _refresh_semantic(new_sql)
    _only_changes(manifest, new_manifest, {"processing_identity", AUDIT_FIELD, "semantic_sha256"}, "manifest")
    _only_changes(acceptance, new_acceptance, {"manifest_semantic_sha256", AUDIT_FIELD, "semantic_sha256"}, "acceptance")
    _only_changes(record_config, new_record, {"accepted_cache_manifest_semantic_sha256", "acceptance_receipt_file_sha256", AUDIT_FIELD}, "record config")
    _only_changes(sql_config, new_sql, {"processing_identity", "cache_manifest_semantic_sha256", AUDIT_FIELD, "semantic_sha256"}, "SQL metadata")
    validate_acceptance(new_acceptance, new_manifest)
    if new_record["acceptance_receipt_file_sha256"] != _sha256_bytes(acceptance_bytes):
        raise RuntimeError("staged acceptance binding differs")
    if new_sql["cache_manifest_semantic_sha256"] != new_manifest["semantic_sha256"]:
        raise RuntimeError("staged SQL manifest binding differs")
    return new_manifest, new_acceptance, new_record, new_sql


def _port_8516_closed() -> bool:
    for host in ("127.0.0.1", "::1"):
        family = socket.AF_INET6 if ":" in host else socket.AF_INET
        sock = socket.socket(family, socket.SOCK_STREAM)
        sock.settimeout(0.25)
        try:
            if sock.connect_ex((host, 8516)) == 0:
                return False
        finally:
            sock.close()
    return True


@contextmanager
def _exclusive_lock(directory: Path) -> Iterator[None]:
    directory.mkdir(parents=True, exist_ok=True)
    lock_path = directory / ".spectral_evidence_rebind.lock"
    with lock_path.open("a+b") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"another metadata rebind owns {lock_path}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _backup_exclusive(paths: tuple[Path, ...], destination: Path) -> dict[str, Any]:
    destination.mkdir(parents=False, exist_ok=False)
    result: dict[str, Any] = {}
    for path in paths:
        target = destination / path.name
        with path.open("rb") as source, target.open("xb") as output:
            shutil.copyfileobj(source, output)
            output.flush()
            os.fsync(output.fileno())
        result[str(path)] = {
            "backup_path": str(target),
            "sha256": file_sha256(target),
            "size_bytes": target.stat().st_size,
        }
    return result


def _validate_staged(stage: Path) -> None:
    manifest = _read_json(stage / MANIFEST_PATH.name)
    acceptance = _read_json(stage / ACCEPTANCE_PATH.name)
    record = _read_json(stage / RECORD_CONFIG_PATH.name)
    sql = _read_json(stage / SQL_CONFIG_PATH.name)
    validate_acceptance(acceptance, manifest)
    _verify_semantic(manifest, "staged manifest")
    _verify_semantic(acceptance, "staged acceptance")
    _verify_semantic(sql, "staged SQL metadata")
    if record["accepted_cache_manifest_semantic_sha256"] != manifest["semantic_sha256"]:
        raise RuntimeError("staged record config manifest binding differs")
    if record["acceptance_receipt_file_sha256"] != file_sha256(stage / ACCEPTANCE_PATH.name):
        raise RuntimeError("staged record config acceptance binding differs")
    if sql["cache_manifest_semantic_sha256"] != manifest["semantic_sha256"]:
        raise RuntimeError("staged SQL manifest binding differs")
    if sql["processing_identity"] != manifest["processing_identity"]:
        raise RuntimeError("staged SQL processing identity differs")


def rebind(protocol_path: Path, protocol_sha256: str, output_path: Path, *, apply: bool) -> dict[str, Any]:
    protocol_path = protocol_path.resolve()
    output_path = output_path.resolve()
    if output_path in {path.resolve() for path in METADATA_PATHS}:
        raise RuntimeError("output receipt cannot replace a metadata input")
    protocol, receipts = _load_protocol(protocol_path, protocol_sha256)
    current_identity, source_proof = _source_proof(protocol)
    originals = {path: path.read_bytes() for path in METADATA_PATHS}
    original_hashes = {str(path): _sha256_bytes(content) for path, content in originals.items()}
    manifest, acceptance, record_config, sql_config = (
        json.loads(originals[path].decode("utf-8")) for path in METADATA_PATHS
    )
    retained = _verify_preconditions(
        manifest, acceptance, record_config, sql_config, current_identity, protocol
    )
    timestamp = datetime.now(timezone.utc).isoformat()
    audit = {
        "schema_version": "mssystem-spectral-evidence-source-rebind/1.0",
        "rebound_at": timestamp,
        "protocol_path": str(protocol_path),
        "protocol_sha256": protocol_sha256,
        "test_receipt": {k: v for k, v in receipts["test"].items() if k != "payload"},
        "review_receipt": {k: v for k, v in receipts["review"].items() if k != "payload"},
        "source_proof": source_proof,
        "previous_processing_identity": retained["accepted_identity"],
        "current_processing_identity": current_identity,
        "original_metadata_sha256": original_hashes,
        "inherited_historical_acceptance_status": list(HISTORICAL_AUDIT_FIELDS),
        "acceptance_results_retained": retained["acceptance_results"],
        "record_count": manifest["record_count"],
        "complete_library": manifest["complete_library"],
        "normalized_records_sha256": manifest["normalized_records_sha256"],
        "record_alignment_sha256": manifest["record_alignment_sha256"],
        "source_db": manifest["source_db"],
        "payload_bindings": manifest["files"],
        "normalized_record_content_changed": False,
        "payload_files_written": False,
        "sqlite_library_written": False,
        "candidate_scores_changed": False,
        "successful_result_fields": ["spectral_evidence"],
        "fresh_fault_isolation": True,
        "normal_runtime_loader_validation": "pending",
        "scope_excludes_unconditional_all_resource_equivalence": True,
    }
    proposed = _build_proposed(
        manifest, acceptance, record_config, sql_config, current_identity, audit
    )
    proposed_bytes = {path: _json_bytes(payload) for path, payload in zip(METADATA_PATHS, proposed)}
    result: dict[str, Any] = {
        "schema_version": "mssystem-spectral-evidence-cache-rebind/1.0",
        "status": "preview",
        "generated_at": timestamp,
        "protocol_path": str(protocol_path),
        "protocol_sha256": protocol_sha256,
        "record_store": str(RECORD_STORE),
        "source_proof": source_proof,
        "original_metadata_sha256": original_hashes,
        "proposed_metadata_sha256": {
            str(path): _sha256_bytes(content) for path, content in proposed_bytes.items()
        },
        "new_manifest_semantic_sha256": proposed[0]["semantic_sha256"],
        "normal_runtime_loader_validation": "pending",
        "backup_files": {},
    }
    if not apply:
        _write_json_atomic(output_path, result)
        return result
    if not _port_8516_closed():
        raise RuntimeError("maintenance boundary refused: port 8516 is listening")
    with _exclusive_lock(output_path.parent):
        if not _port_8516_closed():
            raise RuntimeError("maintenance boundary changed: port 8516 is listening")
        if output_path.parent.joinpath("pre_spectral_evidence_rebind").exists():
            raise RuntimeError("exclusive backup directory already exists")
        if any(path.read_bytes() != originals[path] for path in METADATA_PATHS):
            raise RuntimeError("metadata changed before backup")
        backup_dir = output_path.parent / "pre_spectral_evidence_rebind"
        result["backup_files"] = _backup_exclusive(METADATA_PATHS, backup_dir)
        stage = Path(tempfile.mkdtemp(prefix=".spectral-evidence-rebind-stage-", dir=output_path.parent))
        published: dict[Path, str] = {}
        try:
            for path in METADATA_PATHS:
                (stage / path.name).write_bytes(proposed_bytes[path])
            _validate_staged(stage)
            _load_protocol(protocol_path, protocol_sha256)
            repeated_identity, repeated_proof = _source_proof(protocol)
            if repeated_identity != current_identity or repeated_proof != source_proof:
                raise RuntimeError("reviewed source changed immediately before publication")
            if any(_sha256_bytes(path.read_bytes()) != original_hashes[str(path)] for path in METADATA_PATHS):
                raise RuntimeError("metadata changed immediately before publication")
            for path in METADATA_PATHS:
                content = (stage / path.name).read_bytes()
                _write_bytes_atomic(path, content)
                published[path] = _sha256_bytes(content)
            _validate_staged_from_installed()
            result["status"] = "metadata_validated_loader_pending"
            result["applied_at"] = datetime.now(timezone.utc).isoformat()
            _write_json_atomic(output_path, result)
        except Exception:
            for path in reversed(METADATA_PATHS):
                expected = published.get(path)
                if expected is not None and path.is_file() and file_sha256(path) == expected:
                    _write_bytes_atomic(path, originals[path])
            raise
        finally:
            shutil.rmtree(stage)
    return result


def _validate_staged_from_installed() -> None:
    manifest = _read_json(MANIFEST_PATH)
    acceptance = _read_json(ACCEPTANCE_PATH)
    record = _read_json(RECORD_CONFIG_PATH)
    sql = _read_json(SQL_CONFIG_PATH)
    validate_acceptance(acceptance, manifest)
    _verify_semantic(manifest, "installed manifest")
    _verify_semantic(acceptance, "installed acceptance")
    _verify_semantic(sql, "installed SQL metadata")
    if record["accepted_cache_manifest_semantic_sha256"] != manifest["semantic_sha256"]:
        raise RuntimeError("installed record config manifest binding differs")
    if record["acceptance_receipt_file_sha256"] != file_sha256(ACCEPTANCE_PATH):
        raise RuntimeError("installed record config acceptance binding differs")
    if sql["cache_manifest_semantic_sha256"] != manifest["semantic_sha256"]:
        raise RuntimeError("installed SQL manifest binding differs")
    if sql["processing_identity"] != manifest["processing_identity"]:
        raise RuntimeError("installed SQL processing identity differs")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = rebind(
        args.protocol, args.protocol_sha256, args.output, apply=args.apply
    )
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
