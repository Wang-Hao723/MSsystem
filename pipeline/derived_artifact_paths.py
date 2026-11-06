"""Explicit, content-bound derived artifact names; legacy paths stay unchanged.

This module performs no writes or policy/source discovery. Callers must provide
the digest verified by the selected implementation's own runtime guard.
"""
from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from pathlib import Path

_NAME = re.compile(r"[a-z0-9][a-z0-9_-]{0,39}\Z")
_TASK = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,99}\Z")
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
_SUFFIXES = {"_v15_flagship.json", "_v16_flagship.json", "_v11_compounds.json"}


def _canonical(value):
    """Strict JSON avoids repr-based identities and lossy numeric-key coercion."""
    def check(item):
        if isinstance(item, Mapping):
            if any(not isinstance(key, str) for key in item):
                raise ValueError("Artifact policy keys must be strings")
            for child in item.values():
                check(child)
        elif isinstance(item, list):
            for child in item:
                check(child)
        elif item is not None and type(item) not in (str, bool, int, float):
            raise ValueError("Artifact binding accepts JSON values only")
    check(value)
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).encode("utf-8")


def implementation_bundle_sha256(component_digests: Mapping[str, str]) -> str:
    if not component_digests or any(
        not isinstance(name, str) or not name or not isinstance(value, str)
        or not _DIGEST.fullmatch(value)
        for name, value in component_digests.items()
    ):
        raise ValueError("Implementation components require named SHA-256 values")
    return hashlib.sha256(_canonical(dict(component_digests))).hexdigest()


def artifact_binding(policy: Mapping, *, actual_implementation_sha256=None):
    """None means the exact legacy path; malformed explicit opt-in is an error."""
    if "artifact_namespace" not in policy:
        return None
    namespace = policy["artifact_namespace"]
    if not isinstance(namespace, Mapping) or set(namespace) != {
        "name", "implementation_sha256"
    }:
        raise ValueError("Explicit artifact_namespace requires name and implementation_sha256")
    name, declared = namespace["name"], namespace["implementation_sha256"]
    if not isinstance(name, str) or not _NAME.fullmatch(name):
        raise ValueError("Invalid artifact namespace name")
    if not isinstance(declared, str) or not _DIGEST.fullmatch(declared):
        raise ValueError("Invalid declared artifact implementation SHA-256")
    if not isinstance(actual_implementation_sha256, str) or not _DIGEST.fullmatch(actual_implementation_sha256):
        raise ValueError("Versioned artifacts require a verified actual implementation digest")
    if declared != actual_implementation_sha256:
        raise ValueError("Declared artifact implementation differs from loaded implementation")
    algorithm = policy.get("algorithm_id")
    if not isinstance(algorithm, str) or not algorithm:
        raise ValueError("Versioned artifacts require an explicit algorithm_id")
    return {
        "schema_version": "derived-artifact-binding/1",
        "namespace": name,
        "policy_sha256": hashlib.sha256(_canonical(policy)).hexdigest(),
        "algorithm_id": algorithm,
        "declared_implementation_sha256": declared,
        "verified_implementation_sha256": actual_implementation_sha256,
        "derived_reason": "explicit_policy_namespace_materialization_from_saved_scan_evidence",
        "legacy_artifact_overwritten": False,
    }


def derived_artifact_path(task_dir, task_id, legacy_suffix, *, policy,
                          actual_implementation_sha256=None):
    if legacy_suffix not in _SUFFIXES:
        raise ValueError("Unsupported derived artifact type")
    binding = artifact_binding(policy, actual_implementation_sha256=actual_implementation_sha256)
    if binding is None:
        return Path(task_dir) / f"{task_id}{legacy_suffix}"
    if not isinstance(task_id, str) or not _TASK.fullmatch(task_id):
        raise ValueError("Invalid task identifier for a versioned artifact")
    return Path(task_dir) / (
        f"{task_id}_derived_{binding['namespace']}_{binding['policy_sha256']}{legacy_suffix}"
    )


def binding_metadata(policy, *, actual_implementation_sha256=None):
    binding = artifact_binding(policy, actual_implementation_sha256=actual_implementation_sha256)
    return {} if binding is None else {"derived_artifact_binding": binding}


def assert_derived_destination(output_path, task_id, legacy_suffix, *, policy,
                               actual_implementation_sha256=None):
    expected = derived_artifact_path(Path(output_path).parent, task_id, legacy_suffix,
        policy=policy, actual_implementation_sha256=actual_implementation_sha256)
    if "artifact_namespace" in policy and Path(output_path) != expected:
        raise ValueError("Versioned output must use its content-bound destination")


def metadata_binding_is_current(metadata, policy, *, actual_implementation_sha256=None):
    expected = artifact_binding(policy, actual_implementation_sha256=actual_implementation_sha256)
    return expected is None or metadata.get("derived_artifact_binding") == expected


def verified_primary_binding(primary, policy, *, actual_implementation_sha256=None):
    """A new compound view cannot silently consume a legacy/different primary."""
    if "artifact_namespace" not in policy:
        return {}
    if isinstance(primary, Mapping):
        payload = primary
        raw = _canonical(primary)
        basis = "canonical_mapping"
    elif isinstance(primary, (str, Path)):
        # Primary artifacts can contain gigabytes of reference evidence. A
        # freshness check needs their bytes and metadata, not every entity.
        with Path(primary).open("rb") as handle:
            primary_sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
            handle.seek(0)
            try:
                import ijson
            except ImportError:
                payload = json.load(handle)
            else:
                try:
                    payload = {"metadata": next(ijson.items(handle, "metadata", use_float=True), {})}
                except ijson.JSONError as exc:
                    raise ValueError("Primary artifact contains invalid JSON") from exc
        basis = "file_bytes"
    else:
        raise ValueError("Versioned compounds require a bound primary artifact")
    if not isinstance(payload, Mapping) or not metadata_binding_is_current(
        payload.get("metadata") or {}, policy,
        actual_implementation_sha256=actual_implementation_sha256,
    ):
        raise ValueError("Compound and primary artifact bindings differ")
    return {"derived_primary_sha256": (primary_sha256 if basis == "file_bytes" else hashlib.sha256(raw).hexdigest()),
            "derived_primary_hash_basis": basis}
