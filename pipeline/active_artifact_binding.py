"""One explicit implementation provider; no automatic algorithm fallback."""
from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

from pipeline.derived_artifact_paths import implementation_bundle_sha256


_DEPENDENCIES = (
    "pipeline.active_artifact_binding", "pipeline.derived_artifact_paths",
    "pipeline.v15_flagship", "pipeline.v16_candidate_missingness",
    "pipeline.compound_conclusions", "pipeline.entity_confidence",
    "pipeline.reliability_annotations", "pipeline.reliability_contract",
    "pipeline.compound_role_projection",
    "pipeline.replicate_consistency",
    "pipeline.entity_builder", "pipeline.query_input_guards",
    "pipeline.query_search", "pipeline.query_charge_guard", "pipeline.query_record_memo", "pipeline.query_expert_memo",
    "pipeline.mzml_input", "pipeline.mzml_metadata",
    "pipeline.structure_prediction",
    "web.server", "web.app", "web.governance_blocked_projection",
    "matching.flagship_ranker", "matching.v15_flagship_ranker",
    "matching.v15_identity_resolution", "matching.v15_ce_consensus",
    "matching.v15_decision_assurance",
    "matching.conditioned_multilibrary", "matching.library_quality",
    "matching.compound_names", "matching.biosource_context",
    "matching.fragmentation_rules",
    "matching.search", "matching.reference_mass",
    "matching.untargeted_workflow",
    "matching.source_ms_stage", "matching.evirank_retrieval",
    "utils.adduct", "utils.collision_energy", "utils.instrument_metadata",
    "runtime.source_compound_mass_core", "runtime.source_compound_mass_repairs",
    "runtime.library_metadata_repairs", "runtime.verified_sql_retrieval", "runtime.analysis_completion",
    "runtime.legacy_library_session", "runtime.source_file_verification", "runtime.json_encoding",
)
_SOURCE_PATHS = {
    name: Path(importlib.util.find_spec(name).origin)
    for name in _DEPENDENCIES
}
_LOADED_SOURCE_DIGESTS = {
    name: hashlib.sha256(path.read_bytes()).hexdigest()
    for name, path in _SOURCE_PATHS.items()
}

# These small configurations bind the accepted normalized library and the
# verified SQL/metadata overlay. Their content can change while all ranking
# code remains identical; the resulting reports need a different destination.
# Runtime validation separately verifies their referenced large files.
_DATA_PATHS = {
    "data:" + name: Path(__file__).resolve().parents[1] / name
    for name in (
        "config/library_record_store.json",
        "config/sql_metadata_repairs.json",
        "config/source_compound_mass_repairs.json",
    )
}
_LOADED_DATA_DIGESTS = {
    name: hashlib.sha256(path.read_bytes()).hexdigest()
    for name, path in _DATA_PATHS.items()
}


def active_implementation_sha256(policy):
    if "artifact_namespace" not in policy:
        return None
    # Only the existing reviewed exact-tie implementation is connected here.
    # A continuous-score adapter needs its own explicit provider after review.
    from pipeline import spectral_tie_ranker, spectral_tie_scope
    config = policy.get("spectral_tie_handling") or {}
    if (config.get("enabled") is not True
            or config.get("algorithm_id") != spectral_tie_ranker.ALGORITHM_ID
            or policy.get("algorithm_id") != "mssystem-v16-shared-spectral-ties-multilibrary-1.0.0"):
        raise ValueError("No reviewed artifact binding provider for this scoring policy")
    for name, path in _SOURCE_PATHS.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != _LOADED_SOURCE_DIGESTS[name]:
            raise ValueError("Derived artifact implementation changed; restart and rebind: " + name)
    for name, path in _DATA_PATHS.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != _LOADED_DATA_DIGESTS[name]:
            raise ValueError("Derived artifact library binding changed; restart and rebind: " + name)
    return implementation_bundle_sha256({
        **_LOADED_SOURCE_DIGESTS,
        **_LOADED_DATA_DIGESTS,
        "ranker": spectral_tie_ranker.implementation_sha256(),
        "scope": spectral_tie_scope.implementation_sha256(),
    })
