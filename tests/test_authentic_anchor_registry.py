from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from governance.authentic_anchor_registry import (
    audit_active_authentic_anchor_registry,
    install_active_authentic_anchor_library,
)
from matching.authentic_anchor import AuthenticAnchorPolicy, finalize_anchor_library


def _library() -> dict:
    return finalize_anchor_library(
        {
            "library_id": "registry-test",
            "runtime_mode": "shadow_only",
            "construction_policy": {
                "explicit_task_id_only": True,
                "candidate_generation_allowed": False,
                "authentic_library_shadow_search_allowed": True,
            },
            "release_policy": asdict(AuthenticAnchorPolicy().validated()),
            "entries": [
                {
                    "reference_id": "A:NEG",
                    "compound_code": "A",
                    "compound_name": "A",
                    "formula": "C10H10O4",
                    "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                    "neutral_mass": 194.05791,
                    "polarity": "NEG",
                    "instrument_family": "agilent-qtof-6546",
                    "authentic_standard": True,
                    "source_task_id": "standard-a",
                    "spectra": [
                        {
                            "collision_energy_ev": 10.0,
                            "peaks": [[80.0, 100.0], [100.0, 60.0], [120.0, 40.0]],
                        },
                        {
                            "collision_energy_ev": 20.0,
                            "peaks": [[80.0, 80.0], [100.0, 100.0], [140.0, 30.0]],
                        },
                    ],
                }
            ],
        }
    )


def _registry(path: Path, library: dict, *, raw_hash: str = "") -> Path:
    payload = {
        "schema_version": "authentic-anchor-registry/1.0",
        "active_library": "runtime/authentic.json",
        "active_library_content_sha256": library["content_sha256"],
        "state": "active_authentic_reference_display_ranker",
        "instrument_family": "agilent-qtof-6546",
        "formal_level1_release_enabled": False,
    }
    if raw_hash:
        payload["active_library_file_sha256"] = raw_hash
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_active_registry_fails_closed_when_declared_library_is_missing(tmp_path):
    library = _library()
    registry = _registry(tmp_path / "registry.json", library)

    audit = audit_active_authentic_anchor_registry(registry, project_root=tmp_path)

    assert audit["passed"] is False
    assert audit["violations"] == ["active_library_missing"]


def test_explicit_install_is_hash_bound_and_does_not_mutate_source(tmp_path):
    library = _library()
    source = tmp_path / "source.json"
    source.write_text(json.dumps(library), encoding="utf-8")
    import hashlib

    raw_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    registry = _registry(tmp_path / "registry.json", library, raw_hash=raw_hash)
    before = source.read_bytes()

    receipt = install_active_authentic_anchor_library(
        source,
        registry,
        project_root=tmp_path,
    )

    assert receipt["installed"] is True
    assert receipt["source_mutated"] is False
    assert source.read_bytes() == before
    assert receipt["deployment_audit"]["passed"] is True
    assert receipt["deployment_audit"]["entry_count"] == 1


def test_registry_rejects_semantically_valid_but_raw_hash_changed_file(tmp_path):
    library = _library()
    destination = tmp_path / "runtime" / "authentic.json"
    destination.parent.mkdir()
    destination.write_text(json.dumps(library, indent=2), encoding="utf-8")
    registry = _registry(
        tmp_path / "registry.json",
        library,
        raw_hash="0" * 64,
    )

    audit = audit_active_authentic_anchor_registry(registry, project_root=tmp_path)

    assert audit["passed"] is False
    assert audit["actual_content_sha256"] == library["content_sha256"]
    assert audit["violations"] == ["active_library_file_hash_mismatch"]
