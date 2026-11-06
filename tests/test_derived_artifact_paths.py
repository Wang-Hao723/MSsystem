from __future__ import annotations

import ast
import copy
import hashlib
import importlib.util
import json
import os
import sys
import types
import uuid
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("candidate_paths", HERE / "pipeline/derived_artifact_paths.py")
D = importlib.util.module_from_spec(spec)
spec.loader.exec_module(D)
DIGEST = D.implementation_bundle_sha256({"ranker": "a" * 64, "scope": "b" * 64})


def policy():
    return {"algorithm_id": "mssystem-v16-example-1", "decision": {"minimum": 0.5},
            "artifact_namespace": {"name": "exact-tie-v1", "implementation_sha256": DIGEST}}


def path(tmp, suffix="_v16_flagship.json", p=None):
    return D.derived_artifact_path(tmp, "task1", suffix, policy=p or policy(),
                                   actual_implementation_sha256=DIGEST)


def actual_functions(filename, names, p):
    """Execute exact candidate function bodies with unrelated scoring stubbed."""
    tree = ast.parse((HERE / "pipeline" / filename).read_text())
    selected = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert len(selected) == len(names)
    ns = {"Path": Path, "Mapping": Mapping, "json": json, "hashlib": hashlib,
          "datetime": datetime, "timezone": timezone, "os": os, "uuid": uuid,
          "load_v15_flagship_policy": lambda: copy.deepcopy(p),
          "active_implementation_sha256": lambda _: DIGEST,
          "artifact_implementation_is_current": lambda *args: True,
          "V15_FLAGSHIP_ARTIFACT_SUFFIX": "_v15_flagship.json",
          "V16_FLAGSHIP_ARTIFACT_SUFFIX": "_v16_flagship.json",
          "COMPOUND_CONCLUSION_ARTIFACT_SUFFIX": "_v11_compounds.json"}
    for name in ("derived_artifact_path", "binding_metadata", "assert_derived_destination",
                 "metadata_binding_is_current", "verified_primary_binding"):
        ns[name] = getattr(D, name)
    module = ast.Module(body=[ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0), *selected], type_ignores=[])
    exec(compile(ast.fix_missing_locations(module), filename, "exec"), ns)
    return ns


@pytest.mark.parametrize("algorithm,suffix", [("mssystem-v15-example", "_v15_flagship.json"), ("mssystem-v16-example", "_v16_flagship.json")])
def test_legacy_path_helpers_exact(tmp_path, algorithm, suffix):
    p = {"algorithm_id": algorithm}
    f = actual_functions("v15_flagship.py", {"v15_flagship_artifact_path"}, p)
    c = actual_functions("compound_conclusions.py", {"compound_conclusion_artifact_path"}, p)
    assert f["v15_flagship_artifact_path"](tmp_path, "task1") == tmp_path / ("task1" + suffix)
    assert c["compound_conclusion_artifact_path"](tmp_path, "task1") == tmp_path / "task1_v11_compounds.json"


def test_both_actual_path_helpers_share_policy_binding(tmp_path):
    f = actual_functions("v15_flagship.py", {"v15_flagship_artifact_path"}, policy())
    c = actual_functions("compound_conclusions.py", {"compound_conclusion_artifact_path"}, policy())
    primary = f["v15_flagship_artifact_path"](tmp_path, "task1")
    compound = c["compound_conclusion_artifact_path"](tmp_path, "task1")
    assert primary == path(tmp_path)
    assert compound == path(tmp_path, "_v11_compounds.json")
    assert primary.name.removesuffix("_v16_flagship.json") == compound.name.removesuffix("_v11_compounds.json")
    assert D.artifact_binding(policy(), actual_implementation_sha256=DIGEST)["policy_sha256"] in primary.name
    assert not list(tmp_path.iterdir())  # name computation never materializes files


def test_same_algorithm_changed_policy_or_implementation_gets_new_path(tmp_path):
    original = policy(); changed = copy.deepcopy(original); changed["decision"]["minimum"] = 0.6
    assert path(tmp_path, p=original) != path(tmp_path, p=changed)
    changed = copy.deepcopy(original); changed["artifact_namespace"]["implementation_sha256"] = "c" * 64
    newer = D.derived_artifact_path(tmp_path, "task1", "_v16_flagship.json", policy=changed, actual_implementation_sha256="c" * 64)
    assert newer != path(tmp_path)
    assert path(tmp_path, p=dict(reversed(list(original.items())))) == path(tmp_path)


@pytest.mark.parametrize("actual", [None, "", "c" * 64])
def test_missing_or_stale_implementation_cannot_reuse_path(tmp_path, actual):
    with pytest.raises(ValueError):
        D.derived_artifact_path(tmp_path, "task1", "_v16_flagship.json", policy=policy(), actual_implementation_sha256=actual)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("namespace", [None, "", {}, {"name": "x"}, {"name": "../x", "implementation_sha256": DIGEST}, {"name": "x", "implementation_sha256": "old"}, {"name": "x", "implementation_sha256": DIGEST, "ignored": True}])
def test_invalid_explicit_namespace_never_falls_back(tmp_path, namespace):
    p = policy(); p["artifact_namespace"] = namespace
    with pytest.raises(ValueError):
        path(tmp_path, p=p)


@pytest.mark.parametrize("bad", [float("nan"), {1: "numeric key"}, (1, 2)])
def test_policy_digest_cannot_hide_lossy_json_values(bad):
    p = policy(); p["bad"] = bad
    with pytest.raises(ValueError):
        D.artifact_binding(p, actual_implementation_sha256=DIGEST)


def test_full_metadata_and_primary_binding_are_required(tmp_path):
    p = policy(); meta = D.binding_metadata(p, actual_implementation_sha256=DIGEST)
    assert D.metadata_binding_is_current(meta, p, actual_implementation_sha256=DIGEST)
    assert not D.metadata_binding_is_current({}, p, actual_implementation_sha256=DIGEST)
    with pytest.raises(ValueError):
        D.verified_primary_binding({"metadata": {}}, p, actual_implementation_sha256=DIGEST)
    primary = path(tmp_path); primary.write_text(json.dumps({"metadata": meta, "entities": []}))
    bound = D.verified_primary_binding(primary, p, actual_implementation_sha256=DIGEST)
    assert bound["derived_primary_sha256"] == hashlib.sha256(primary.read_bytes()).hexdigest()
    assert bound["derived_primary_hash_basis"] == "file_bytes"
    assert meta["derived_artifact_binding"]["legacy_artifact_overwritten"] is False


def test_actual_writers_preserve_old_files_and_reject_old_destination(tmp_path):
    p = policy(); old = {}
    for suffix in ("_v15_flagship.json", "_v16_flagship.json", "_v11_compounds.json"):
        f = tmp_path / ("task1" + suffix); f.write_bytes(("original " + suffix).encode()); old[f] = f.read_bytes()
    ns = actual_functions("v15_flagship.py", {"write_v15_flagship_artifact", "_write_payload_atomic"}, p)
    ns.update(build_v15_flagship_result=lambda *a, **kw: {"metadata": {"algorithm_id": p["algorithm_id"]}, "policy": p, "entities": []},
              attach_replicate_missingness_audit=lambda x, _: x, apply_v15_identity_resolution=lambda x: x,
              apply_spectral_tie_scope_guard=lambda x: x,
              attach_v15_decision_assurance=lambda x: x, _entity_semantic_sha256=lambda _: "d" * 64)
    with pytest.raises(ValueError): ns["write_v15_flagship_artifact"]("task1", [], tmp_path / "task1_v16_flagship.json")
    primary = path(tmp_path); ns["write_v15_flagship_artifact"]("task1", [], primary)
    nc = actual_functions("compound_conclusions.py", {"write_compound_conclusion_artifact"}, p)
    nc["build_compound_conclusions"] = lambda *a, **kw: {"metadata": {}, "compounds": []}
    with pytest.raises(ValueError): nc["write_compound_conclusion_artifact"]("task1", [], tmp_path / "task1_v11_compounds.json", primary_identity_artifact=primary)
    with pytest.raises(ValueError): nc["write_compound_conclusion_artifact"]("task1", [], path(tmp_path, "_v11_compounds.json"), primary_identity_artifact=tmp_path / "task1_v16_flagship.json")
    compound = path(tmp_path, "_v11_compounds.json")
    nc["write_compound_conclusion_artifact"]("task1", [], compound, primary_identity_artifact=primary)
    assert json.loads(primary.read_text())["metadata"]["derived_artifact_binding"] == json.loads(compound.read_text())["metadata"]["derived_artifact_binding"]
    assert all(f.read_bytes() == content for f, content in old.items())
    assert set(tmp_path.glob("task1_*")) == {*old, primary, compound}  # existing recycle ownership rule


def test_actual_replicate_reader_uses_shared_path(tmp_path, monkeypatch):
    c = types.ModuleType("pipeline.compound_conclusions")
    c.compound_conclusion_artifact_path = lambda directory, task_id: path(directory, "_v11_compounds.json")
    package = types.ModuleType("pipeline"); package.__path__ = []
    monkeypatch.setitem(sys.modules, "pipeline", package)
    monkeypatch.setitem(sys.modules, "pipeline.compound_conclusions", c)
    (tmp_path / "task1_v11_compounds.json").write_text(json.dumps({"compounds": ["legacy"]}))
    path(tmp_path, "_v11_compounds.json").write_text(json.dumps({"compounds": ["versioned"]}))
    ns = actual_functions("replicate_consistency.py", {"_read_compounds"}, policy())
    assert ns["_read_compounds"](tmp_path, "task1")["compounds"] == ["versioned"]


def test_provider_checks_selected_implementation_and_never_guesses(monkeypatch):
    import pipeline
    from pipeline import active_artifact_binding as module

    ranker = types.SimpleNamespace(ALGORITHM_ID="exact-tie", implementation_sha256=lambda: "a" * 64)
    scope = types.SimpleNamespace(implementation_sha256=lambda: "b" * 64)
    # Build the real source registry before substituting only its consumers.
    # A synthetic package with an empty path prevents sibling discovery.
    monkeypatch.setattr(pipeline, "spectral_tie_ranker", ranker)
    monkeypatch.setattr(pipeline, "spectral_tie_scope", scope)
    assert module.active_implementation_sha256({}) is None
    p = policy()
    with pytest.raises(ValueError): module.active_implementation_sha256(p)
    p["spectral_tie_handling"] = {"enabled": True, "algorithm_id": "exact-tie"}
    p["algorithm_id"] = "mssystem-v16-shared-spectral-ties-multilibrary-1.0.0"
    expected = D.implementation_bundle_sha256({**module._LOADED_SOURCE_DIGESTS,
                                               **module._LOADED_DATA_DIGESTS,
                                               "ranker": "a" * 64, "scope": "b" * 64})
    assert module.active_implementation_sha256(p) == expected
    def changed(): raise ValueError("loaded source changed")
    ranker.implementation_sha256 = changed
    with pytest.raises(ValueError, match="loaded source changed"): module.active_implementation_sha256(p)


@pytest.mark.parametrize("missing_required", [False, True])
def test_provider_discovery_in_clean_process_preserves_required_dependencies(missing_required):
    import subprocess

    # Each discovery scenario starts without the parent process's module cache.
    script = '''
import importlib.util
import json
from pathlib import Path
import sys
missing = sys.argv[1] == 'true'
original = importlib.util.find_spec
target = 'pipeline.v16_candidate_missingness'
if missing:
    importlib.util.find_spec = lambda name, *a, **k: None if name == target else original(name, *a, **k)
try:
    from pipeline import active_artifact_binding as provider
except AttributeError as exc:
    if not missing or 'origin' not in str(exc):
        raise
else:
    assert not missing, 'Required missing dependency was accepted'
    assert set(provider._SOURCE_PATHS) == set(provider._DEPENDENCIES)
    assert all(path.is_file() for path in provider._SOURCE_PATHS.values())
    policy = json.loads(Path('config/v15_flagship_policy.json').read_text())
    assert provider.active_implementation_sha256(policy) == policy['artifact_namespace']['implementation_sha256']
'''
    subprocess.run([sys.executable, "-c", script, str(missing_required).lower()],
                   cwd=HERE, check=True, capture_output=True, text=True)


def test_compound_freshness_uses_primary_hash_even_when_size_mtime_match(tmp_path):
    p = policy(); primary = path(tmp_path); source = tmp_path / "task1_results.json"
    source.write_text("[]")
    primary.write_text(json.dumps({"metadata": D.binding_metadata(p, actual_implementation_sha256=DIGEST), "entities": [1]}))
    meta = {"algorithm_id": "compound-test", "report_role_catalog": {},
            "source_artifact_size_bytes": source.stat().st_size,
            "source_artifact_mtime_ns": source.stat().st_mtime_ns,
            "primary_identity_artifact_size_bytes": primary.stat().st_size,
            "primary_identity_artifact_mtime_ns": primary.stat().st_mtime_ns,
            **D.binding_metadata(p, actual_implementation_sha256=DIGEST),
            **D.verified_primary_binding(primary, p, actual_implementation_sha256=DIGEST)}
    compound = path(tmp_path, "_v11_compounds.json"); compound.write_text(json.dumps({"metadata": meta}))
    ns = actual_functions("compound_conclusions.py", {"artifact_is_current"}, p)
    ns.update(COMPOUND_CONCLUSION_ALGORITHM_ID="compound-test", report_role_catalog_summary=lambda: {})
    assert ns["artifact_is_current"](compound, source, primary)
    st = primary.stat(); original = primary.read_text()
    primary.write_text(original.replace('"entities": [1]', '"entities": [2]'))
    os.utime(primary, ns=(st.st_atime_ns, st.st_mtime_ns))
    assert primary.stat().st_size == st.st_size and primary.stat().st_mtime_ns == st.st_mtime_ns
    assert not ns["artifact_is_current"](compound, source, primary)


def test_primary_freshness_rejects_missing_binding_and_runtime_digest_change(tmp_path):
    p = policy(); file = path(tmp_path)
    meta = {"algorithm_id": p["algorithm_id"], "policy_sha256": "policy-test", "identity_resolution_status": "available"}
    file.write_text(json.dumps({"metadata": meta}))
    ns = actual_functions("v15_flagship.py", {"artifact_is_current"}, p)
    ns.update(flagship_policy_sha256=lambda _: "policy-test", _read_artifact_metadata_cached=lambda f, *args: json.loads(Path(f).read_text())["metadata"])
    assert not ns["artifact_is_current"](file)
    meta.update(D.binding_metadata(p, actual_implementation_sha256=DIGEST)); file.write_text(json.dumps({"metadata": meta}))
    assert ns["artifact_is_current"](file)
    ns["active_implementation_sha256"] = lambda _: "f" * 64
    assert not ns["artifact_is_current"](file)
