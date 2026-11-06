import ast
import hashlib
import json
from pathlib import Path

import pytest

import matching.library_record_store as store
import runtime.library_cache_identity as resolver
import runtime.manifest as runtime_manifest
import runtime.source_compound_mass_repairs as mass_repairs


ROOT = Path(__file__).resolve().parents[1]


def _record(root, relative):
    path = root / relative
    data = path.read_bytes()
    return {"path": relative, "size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def _identity(source, *, top_n=15):
    return {
        "source_files_sha256": source,
        "python": [3, 13],
        "numpy": "2.1.3",
        "load_name_only_leads": False,
        "library_top_n_peaks": top_n,
        "library_relative_threshold": 0.05,
    }


@pytest.fixture
def isolated_snapshot(tmp_path, monkeypatch):
    (tmp_path / "matching").mkdir()
    (tmp_path / "matching/base.py").write_bytes(b"VALUE = 1\n")
    (tmp_path / "matching/isolated.py").write_bytes(b"ISOLATED = 1\n")
    files = [_record(tmp_path, "matching/base.py"), _record(tmp_path, "matching/isolated.py")]
    binding = files[1]
    monkeypatch.setattr(resolver, "_ISOLATED_ADDITIONS", {
        "matching/isolated.py": {
            "size_bytes": binding["size_bytes"],
            "sha256": binding["sha256"],
            "import_token": b"matching.isolated",
        }
    })
    monkeypatch.setattr(resolver, "_processing_files", lambda root: list(files))
    accepted = _identity(resolver._semantic_sha256(files[:1]))
    current = _identity(resolver._semantic_sha256(files))
    return tmp_path, files, accepted, current


def test_reviewed_isolated_addition_must_reconstruct_exact_accepted_fingerprint(isolated_snapshot):
    root, _, accepted, current = isolated_snapshot
    result = resolver.processing_identity_compatibility(accepted, current, project_root=root)
    assert result["compatible"] is True
    assert result["mode"] == "verified_isolated_non_normalizing_additions"
    assert result["reconstructed_source_files_sha256"] == accepted["source_files_sha256"]


def test_changed_reviewed_file_unknown_file_import_or_option_fails_closed(isolated_snapshot, monkeypatch):
    root, files, accepted, current = isolated_snapshot

    changed = [files[0], {**files[1], "sha256": "0" * 64}]
    monkeypatch.setattr(resolver, "_processing_files", lambda project_root: list(changed))
    changed_current = _identity(resolver._semantic_sha256(changed))
    result = resolver.processing_identity_compatibility(accepted, changed_current, project_root=root)
    assert result == {"compatible": False, "reason": "reviewed_isolated_addition_changed", "path": "matching/isolated.py"}

    unknown_path = root / "matching/unknown.py"
    unknown_path.write_bytes(b"UNKNOWN = 1\n")
    unknown = files + [_record(root, "matching/unknown.py")]
    monkeypatch.setattr(resolver, "_processing_files", lambda project_root: list(unknown))
    result = resolver.processing_identity_compatibility(
        accepted, _identity(resolver._semantic_sha256(unknown)), project_root=root,
    )
    assert result["compatible"] is False
    assert result["reason"] == "unreviewed_processing_source_change"

    (root / "matching/base.py").write_bytes(b"import matching.isolated\n")
    imported = [_record(root, "matching/base.py"), files[1]]
    monkeypatch.setattr(resolver, "_processing_files", lambda project_root: list(imported))
    result = resolver.processing_identity_compatibility(
        accepted, _identity(resolver._semantic_sha256(imported)), project_root=root,
    )
    assert result["compatible"] is False
    assert result["reason"] == "isolated_addition_imported_by_processing_source"

    result = resolver.processing_identity_compatibility(
        accepted, {**current, "library_top_n_peaks": 40}, project_root=root,
    )
    assert result == {"compatible": False, "reason": "processing_options_changed"}


def _acceptance(manifest):
    payload = {
        "schema_version": "mssystem.library-record-store-acceptance/1.0.0",
        "manifest_semantic_sha256": manifest["semantic_sha256"],
        "full_record_roundtrip_sha256": manifest["normalized_records_sha256"],
        "verified_record_count": manifest["record_count"],
        "candidate_order_cases": 10000,
        "candidate_order_mismatches": 0,
        "paired_scoring_queries": 100,
        "score_mismatches": 0,
        "shared_database_unchanged": True,
    }
    payload["semantic_sha256"] = store.semantic_sha256(payload)
    return payload


def test_real_loader_uses_same_resolver_without_rebinding_manifest(tmp_path, monkeypatch):
    source = tmp_path / "source.db"
    source.write_bytes(b"source")
    cache = tmp_path / "cache"
    accepted = _identity("accepted")
    current = _identity("current")
    manifest = store.build_record_store(
        [{"spectrum_id": "one", "precursor_mz": 101.0, "neutral_mass": 100.0, "dataset_source": "test"}],
        cache,
        source_db=source,
        identity=accepted,
    )
    (cache / "acceptance.json").write_text(json.dumps(_acceptance(manifest)), encoding="utf-8")
    manifest_before = (cache / "manifest.json").read_bytes()
    seen = []

    def compatible(left, right, *, project_root):
        seen.append((left, right, Path(project_root)))
        return {"compatible": True, "mode": "test_verified"}

    monkeypatch.setattr(resolver, "processing_identity_compatibility", compatible)
    records, precursor, neutral = resolver.load_compatible_record_store(
        cache,
        source_db=source,
        current_identity=current,
        project_root=tmp_path,
        cache_bytes=store.DEFAULT_CACHE_BYTES,
    )
    assert records[0]["spectrum_id"] == "one"
    assert seen == [(accepted, current, tmp_path)]
    assert (cache / "manifest.json").read_bytes() == manifest_before
    records.close()
    del precursor, neutral


def test_active_cache_health_resolves_only_the_two_exact_masscube_additions():
    manifest = json.loads((ROOT / "output/library_record_store/full-20260915-formula-composition-v1/manifest.json").read_text())
    accepted = manifest["processing_identity"]
    current = store.processing_identity()
    result = resolver.processing_identity_compatibility(accepted, current, project_root=ROOT)
    assert result["compatible"] is True
    assert result["mode"] == "verified_isolated_non_normalizing_additions"
    assert [row["path"] for row in result["excluded_additions"]] == [
        "matching/masscube_peak_evidence.py", "matching/masscube_peak_workflow.py"
    ]
    health = runtime_manifest.library_record_store_compatibility(ROOT)
    assert health["compatible"] is True
    assert health["processing_identity_resolution"] == result


def test_source_mass_direct_loader_passes_its_project_root_to_shared_loader():
    tree = ast.parse(Path(mass_repairs.__file__).read_text(encoding="utf-8"))
    calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Name) and node.func.id == "load_compatible_record_store"]
    assert len(calls) == 1
    keywords = {keyword.arg: ast.unparse(keyword.value) for keyword in calls[0].keywords}
    assert keywords["project_root"] == "root"
