"""
MSsystem Web Server — FastAPI Backend
─────────────────────────────────────
REST API for the MSsystem LC-MS/MS spectral matching system.

Start: bash run.sh
Or:    MSSYSTEM_PORT=8516 python web/server.py
"""
from __future__ import annotations

import csv
import fcntl
import hashlib
import io
import importlib
import json
from runtime.json_encoding import write_json_encoded
import logging
import math
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import time
import uuid
import zipfile
from copy import deepcopy
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from functools import wraps
from pathlib import Path, PurePosixPath
from typing import Any, Literal, Mapping, Optional
from urllib.parse import urlsplit

# ── Project root ────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from fastapi import FastAPI, File, Form, HTTPException, Query, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field, field_validator, model_validator

from runtime.task_state import (
    ACTIVE_TASK_STATES,
    InvalidTaskTransition,
    TaskCancellationRequested,
    TaskState,
    initialize_task_state as _initialize_task_state,
    is_active_task_state as _is_active_task_state,
    transition_task as _transition_task,
)
from runtime.reorderable_executor import ReorderableSerialExecutor

from config.settings import (
    COSINE_THRESHOLD,
    DB_PATH,
    DB_PATH_SMALL,
    LIBRARY_WRITE_ENABLED,
    RUNTIME_ROOT,
    SHARED_DATA_ROOT,
    USER_CONTRIB_DIR,
    VENDOR_CONVERSION_CACHE_ROOT,
    DEFAULT_BLANK_INTENSITY_MODE,
    DEFAULT_BLANK_MODE,
    DEFAULT_BLANK_RT_TOLERANCE_MINUTES,
    DEFAULT_BLANK_MIN_DETECTION_FRACTION,
    DEFAULT_BLANK_SAMPLE_RATIO,
    DEFAULT_BLANK_FILTER_POLICY,
    DEFAULT_BLANK_CLEAN_MAX_RATIO,
    DEFAULT_BLANK_CLEAN_MIN_REPLICATES,
    DEFAULT_BLANK_THRESHOLD,
    FINAL_SCORE_THRESHOLD,
    MAX_CANDIDATES_PER_QUERY,
    MZ_BUCKET_DALTON,
    MZ_TOLERANCE_DA,
    MZ_TOLERANCE_PPM,
    REL_THRESHOLD,
    TOP_N_PEAKS,
    WEIGHT_COSINE,
    WEIGHT_FRAGMENT,
    WEIGHT_INTENSITY,
    WEIGHT_MASS,
    USE_SAMPLE_CONSENSUS_TIER,
    SAMPLE_CONSENSUS_MIN_SPECTRA,
    SAMPLE_CONSENSUS_RT_WINDOW_MIN,
    SAMPLE_CONSENSUS_MIN_COS,
    SAMPLE_CONSENSUS_MIN_FINAL,
    SAMPLE_CONSENSUS_MAX_PPM,
    SAMPLE_CONSENSUS_MIN_CLASS_AGREE,
    USE_STRICT_EVIDENCE_GATE,
    LIMITED_FRAGMENT_STRICT_MAX_MASS_PPM,
    LOAD_NAME_ONLY_LEADS,
    USE_CHIMERIC_DECONVOLUTION,
    USE_ION_IDENTITY_NETWORKING,
    USE_ISOPAIR_MS,
    ION_IDENTITY_RT_TOLERANCE_MIN,
    ION_IDENTITY_MIN_SHAPE_CORRELATION,
    ION_IDENTITY_MIN_OVERLAP,
)
from matching.honest_calibration import annotate as _honest_annotate
from matching.isomer_literature import (
    enrich_isomer_groups_with_pubmed,
    search_pubmed_isomer_evidence,
    write_json_atomic as _write_json_atomic,
)
from matching.compound_names import (
    normalize_match_compound_names as _normalize_match_compound_names,
    normalize_result_compound_names as _normalize_result_compound_names,
)
from pipeline.compound_role_projection import (
    annotate_match_compound_roles as _annotate_match_compound_roles,
    annotate_result_compound_roles as _annotate_result_compound_roles,
    compound_annotation_catalog_summary as _compound_role_catalog_summary,
)
from matching.screening_catalog import catalog_summary as _screening_catalog_summary
from matching.library_quality import (
    compound_identity_key as _compound_identity_key,
    normalize_inchikey as _normalize_inchikey,
    resolve_structure_formula as _resolve_structure_formula,
    validate_structure_metadata as _validate_structure_metadata,
)
from matching.expert_scorer import (
    enrich_neutral_loss_analysis as _enrich_neutral_loss_analysis,
    finalize_expert_evidence as _finalize_expert_evidence,
)
from pipeline.query_input_guards import (
    annotate_precursor_features as _annotate_precursor_features,
    precursor_feature_summary as _precursor_feature_summary,
)
from matching.evidence_governance import (
    COHORT_MZ_PPM as EVIGATE_COHORT_MZ_PPM,
    COHORT_RT_WINDOW_MIN as EVIGATE_COHORT_RT_WINDOW_MIN,
    MIN_IDENTITY_AGREEMENT as EVIGATE_MIN_IDENTITY_AGREEMENT,
    MIN_INDEPENDENT_INJECTIONS as EVIGATE_MIN_INDEPENDENT_INJECTIONS,
    POLICY_ID as EVIGATE_POLICY_ID,
    annotate_evidence_passports as _annotate_evidence_passports,
    build_evidence_passport as _build_evidence_passport,
    build_replicate_consensus as _build_replicate_consensus,
    summarize_evidence_governance as _summarize_evidence_governance,
)
from pipeline.ion_identity_evidence import annotate_ion_identity_networks
from pipeline.mzml_metadata import query_acquisition_projection
from pipeline.reliability_annotations import (
    background_evidence as _background_evidence,
    ce_comparability as _ce_comparability,
    fragmentation_annotation as _fragmentation_annotation,
)
from matching.expected_compounds import (
    list_expected_datasets,
    resolve_expected_dataset_selection,
    run_expected_compounds_workflow,
    unavailable_payload as _expected_compounds_unavailable_payload,
    write_expected_compounds_artifact,
)
from matching.reference_standard_qc import (
    REFERENCE_STANDARD_QC_SCHEMA_VERSION,
    build_reference_standard_qc as _build_reference_standard_qc,
    list_reference_standard_sets as _list_reference_standard_sets,
    read_reference_standard_qc_artifact as _read_reference_standard_qc_artifact,
    reference_standard_qc_artifact_path as _reference_standard_qc_artifact_path,
    validate_reference_standard_set_id as _validate_reference_standard_set_id,
    write_reference_standard_qc_artifact as _write_reference_standard_qc_artifact,
)
from web.expected_compounds_export import (
    ANALYSIS_WARNING_COLUMNS,
    ANALYSIS_WARNING_LABELS,
    EXPECTED_FEATURE_COLUMNS,
    EXPECTED_RAW_FEATURE_COLUMNS,
    analysis_warning_rows,
    append_expected_csv_section,
    expected_feature_column_labels,
    expected_feature_groups,
    expected_group_rows,
    expected_feature_rows,
    expected_raw_feature_column_labels,
    expected_raw_feature_rows,
    expected_pdf_story,
    expected_summary_row,
    inject_expected_html,
    reportlab_font_names,
    render_expected_html_report,
    style_analysis_workbook,
)
from web.reference_standard_qc_export import (
    append_reference_standard_qc_csv_section,
    inject_reference_standard_qc_html,
    reference_standard_qc_csv_text,
    reference_standard_qc_detail_rows,
    reference_standard_qc_pdf_story,
    reference_standard_qc_summary_rows,
)
from web.confidence_guide import (
    CONFIDENCE_GUIDE_POLICY_ID,
    confidence_guide_payload,
    confidence_guide_table_rows,
    confidence_guide_pdf_story,
    format_confidence_guide_worksheet,
    inject_confidence_guide_html,
)
from matching.untargeted_workflow import (
    attach_ms2_to_lcms_features,
    build_lcms_consensus_spectra,
    resolve_search_precursor_mz,
)
from pipeline.reliability_contract import build_reliability_contract
from governance.authentic_anchor_registry import (
    load_active_authentic_anchor_library as _load_active_authentic_anchor_library,
    resolve_active_library_path as _resolve_active_authentic_library_path,
)
from matching.evitwin_ms import (
    build_evitwin_result as _build_evitwin_result,
    read_evitwin_result as _read_evitwin_result,
    write_evitwin_result as _write_evitwin_result,
)
from matching.evitwin_runtime import (
    resolve_evitwin_query_instrument_family as _resolve_evitwin_query_instrument_family,
    resolve_evitwin_query_method_id as _resolve_evitwin_query_method_id,
)
from matching.v16_dual_polarity import (
    build_v16_dual_polarity_shadow as _build_v16_dual_polarity_shadow,
    load_v16_dual_polarity_policy as _load_v16_dual_polarity_policy,
    rt_alignment_semantic_sha256 as _v16_rt_alignment_semantic_sha256,
    semantic_sha256 as _v16_semantic_sha256,
    sidecar_semantic_sha256 as _v16_sidecar_semantic_sha256,
)
from matching.v16_openms_anchor import (
    compact_openms_validation_audit as _v16_compact_openms_validation_audit,
    validate_openms_adduct_group_manifest as _v16_validate_openms_adduct_groups,
)
from matching.flagship_ranker import (
    flagship_policy_sha256 as _flagship_policy_sha256,
)
from matching.v15_flagship_ranker import (
    load_v15_flagship_policy as _load_v15_flagship_policy,
)
from matching.search import (
    _top_candidates_with_best,
    build_deconvolution_run_pool,
    find_matching_fragments,
)
from pipeline.query_search import search
from runtime.source_compound_mass_repairs import source_compound_mass_status
from runtime.library_metadata_repairs import (
    load_library_with_verified_repairs as load_library,
    library_metadata_repair_status,
)
from runtime.verified_sql_retrieval import (
    create_verified_sql_retriever,
    sql_metadata_repair_status,
)
from runtime.verified_sql_artifacts import sql_artifact_is_current
from runtime.analysis_completion import formal_output_check, assert_formal_outputs
from utils.spectrum import (
    aggregate_peak_processing_audits,
    process_peaks,
    resolve_peak_processing_config,
)
from utils.collision_energy import collision_energy_evidence
from utils.storage_guard import (
    InsufficientStorageError,
    estimate_incremental_import_bytes,
    require_storage_capacity,
    storage_status,
)
from pipeline.result_contract import (
    identification_contract_path as _identification_contract_path,
    iter_legacy_result_artifact as _iter_legacy_result_artifact,
    read_identification_contract as _read_identification_contract,
    runtime_snapshot_id_from_manifest as _runtime_snapshot_id_from_manifest,
    write_identification_contract as _write_identification_contract,
)
from pipeline.entity_builder import (
    entity_bundle_path as _entity_bundle_path,
    read_entity_bundle as _read_entity_bundle,
    write_entity_bundle as _write_entity_bundle,
)
from governance.claim_compiler import default_task_claim_context as _default_task_claim_context
from pipeline.claim_contract import (
    claim_contract_path as _claim_contract_path,
    iter_identification_records as _iter_identification_records,
    read_claim_contract as _read_claim_contract,
    write_claim_contract as _write_claim_contract,
)
from pipeline.flagship import (
    flagship_artifact_path as _v11_flagship_artifact_path,
    read_flagship_artifact as _read_v11_flagship_artifact,
    write_flagship_artifact as _write_v11_flagship_artifact,
)
from pipeline.v15_flagship import (
    artifact_is_current as _base_flagship_artifact_is_current,
    read_flagship_artifact as _read_flagship_artifact,
    upgrade_v15_flagship_identity_scope_artifact as _upgrade_flagship_artifact,
    v15_flagship_artifact_path as _flagship_artifact_path,
    write_v15_flagship_artifact as _write_flagship_artifact,
)
from web.v15_multilibrary_export import (
    MULTILIBRARY_COLUMNS,
    multilibrary_column_labels,
    read_multilibrary_evidence,
)
from matching.evirank_retrieval import SQLiteCandidateRecallRetriever
from matching.v15_candidate_recall import load_v15_candidate_recall_policy
from pipeline.evirank_mv import (
    evirank_artifact_path as _evirank_artifact_path,
    read_evirank_artifact as _read_evirank_artifact,
    write_evirank_artifact as _write_evirank_artifact,
)
from pipeline.compound_conclusions import (
    _auxiliary_interpretation as _compound_auxiliary_interpretation,
    artifact_is_current as _compound_conclusion_artifact_is_current,
    compound_conclusion_artifact_path as _compound_conclusion_artifact_path,
    read_compound_conclusion_artifact as _read_compound_conclusion_artifact,
    write_compound_conclusion_artifact as _write_compound_conclusion_artifact,
)
from pipeline.replicate_consistency import (
    build_replicate_consistency as _build_replicate_consistency_view,
    parameter_signature as _replicate_parameter_signature,
    parse_replicate_filename as _parse_replicate_filename,
)
from pipeline.consumer_projection import attach_claim_decisions as _attach_claim_decisions
from pipeline.evograph_shadow import (
    evograph_sidecar_path as _evograph_sidecar_path,
    materialize_evograph_sidecar as _materialize_evograph_sidecar,
    read_evograph_sidecar as _read_evograph_sidecar,
)
from runtime.manifest import (
    runtime_manifest_dependency_signature as _runtime_manifest_dependency_signature,
    runtime_manifest_freshness as _runtime_manifest_freshness,
)
from schemas.acquisition import AcquisitionPlanRequest, FollowupLinkRequest
from schemas.evidiscrim import EviDiscrimArtifact, EviDiscrimRunRequest

# ── Logging ─────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("mssystem.server")

# ── RDKit availability check (for structure images) ──
try:
    from rdkit import Chem as _Chem
    RDKIT_AVAILABLE = True
    logger.info("RDKit available — structure images enabled")
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit NOT available — structure images will be missing in reports!")

# ── FastAPI App ─────────────────────────────────────────────────────
app = FastAPI(
    title="MSsystem · LC-MS/MS Spectral Matching",
    description="REST API for mass spectrometry data analysis",
    version="16.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Paths ───────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
SERVICE_GENERATION = "V16"
DEFAULT_PORT = 8516
OUTPUT_DIR = RUNTIME_ROOT
UPLOAD_DIR = OUTPUT_DIR / "tmp_uploads"
TASK_DIR = OUTPUT_DIR / "tasks"
TASKS_FILE = OUTPUT_DIR / "tasks_index.json"
BLANK_DIR = OUTPUT_DIR / "blank_profiles"   # v3 Phase 1: 上传的空白样本存放处
BLANK_UPLOADS_FILE = OUTPUT_DIR / "blank_uploads_index.json"
VENDOR_RAW_DIR = OUTPUT_DIR / "vendor_raw_staging"
VENDOR_CONVERSION_LEDGER_FILE = (
    OUTPUT_DIR / "runtime" / "vendor_conversion_ledger.json"
)
VENDOR_CONVERSION_VALIDATION_RULES_VERSION = "mssystem.validate-mzml/16.1.0"
RUNTIME_SNAPSHOT_ID = _runtime_snapshot_id_from_manifest(PROJECT_ROOT)
EVITWIN_REGISTRY_FILE = PROJECT_ROOT / "config" / "authentic_anchor_registry.json"


def _evitwin_artifact_path(task_id: str) -> Path:
    return TASK_DIR / f"{task_id}_evitwin_ms.json"


def _evitwin_library_path(registry: Mapping[str, Any]) -> Path:
    return _resolve_active_authentic_library_path(
        registry,
        project_root=PROJECT_ROOT,
    )


def _active_evitwin_library() -> tuple[dict[str, Any], dict[str, Any]]:
    return _load_active_authentic_anchor_library(
        EVITWIN_REGISTRY_FILE,
        project_root=PROJECT_ROOT,
    )


def _task_evitwin_polarity(task: Mapping[str, Any]) -> str:
    filename = str(task.get("filename") or "")
    match = re.search(r"(?:^|_)(POS|NEG)(?:_|\.|$)", filename, re.I)
    if match:
        return match.group(1).upper()
    declared = str(task.get("ion_mode") or "").upper()
    if declared in {"POS", "NEG"}:
        return declared
    counts = task.get("ion_mode_counts") or {}
    observed = [
        mode for mode in ("POS", "NEG") if int(counts.get(mode) or 0) > 0
    ]
    return observed[0] if len(observed) == 1 else ""


def _runtime_manifest_source_tree_sha256(project_root: Path) -> str | None:
    """Read the executable-source identity recorded in a runtime manifest."""

    path = project_root / "config" / "runtime_manifest.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return None
    source = payload.get("source")
    source = source if isinstance(source, dict) else {}
    value = str(source.get("source_tree_sha256") or "").strip()
    return value or None


# A process may adopt an atomically completed manifest after database
# verification only when the executable-source hash is unchanged.  Editing
# Python source on disk must never make an already-loaded process advertise
# the new snapshot before restart.
RUNTIME_SOURCE_TREE_SHA256 = _runtime_manifest_source_tree_sha256(PROJECT_ROOT)
_runtime_snapshot_lock = threading.RLock()
_runtime_manifest_status_check_lock = threading.Lock()
_runtime_snapshot_cache: dict[str, Any] = {
    "signature": None,
    "value": RUNTIME_SNAPSHOT_ID,
}
_runtime_manifest_status_cache: dict[str, Any] = {
    "signature": None,
    "value": None,
}
_VENDOR_LOCAL_SELECTION_TTL_SECONDS = 30 * 60
_vendor_local_selection_lock = threading.RLock()
_vendor_local_selections: dict[str, dict[str, Any]] = {}
_VENDOR_CONVERSION_MAX_WORKERS = 2
_VENDOR_CONVERSION_LEASE_SECONDS = 15 * 60
_VENDOR_CONVERSION_HEARTBEAT_SECONDS = 60.0
_vendor_conversion_slots = threading.BoundedSemaphore(
    _VENDOR_CONVERSION_MAX_WORKERS
)
_VENDOR_CONVERTER_PREFLIGHT_TTL_SECONDS = 30 * 60
_vendor_converter_preflight_lock = threading.Lock()
_vendor_converter_preflight: dict[str, Any] = {
    "prepared_at_monotonic": None,
    "status": None,
}
_vendor_conversion_key_locks_guard = threading.Lock()
_vendor_conversion_key_locks: dict[str, threading.Lock] = {}
_vendor_conversion_service_lock_guard = threading.Lock()
_vendor_conversion_service_lock_handle: Any = None
_vendor_conversion_service_lock_path: Path | None = None
_vendor_conversion_service_lock_refcount = 0
_vendor_conversion_active_claims_lock = threading.Lock()
_vendor_conversion_active_claims: set[str] = set()
_vendor_conversion_accepting_claims = False

for d in [OUTPUT_DIR, UPLOAD_DIR, TASK_DIR, BLANK_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def _current_runtime_snapshot_id() -> str | None:
    """Return the manifest identity attested for this loaded process.

    The final runtime manifest is written atomically after the exact online
    database state has been verified.  The already-loaded Full-library process
    may adopt that identity only when its startup source-tree hash still
    matches; a source-code change requires a restart.
    """
    manifest_path = PROJECT_ROOT / "config" / "runtime_manifest.json"
    try:
        stat = manifest_path.stat()
        signature = (stat.st_mtime_ns, stat.st_size, stat.st_ino)
    except OSError:
        return RUNTIME_SNAPSHOT_ID

    with _runtime_snapshot_lock:
        if _runtime_snapshot_cache.get("signature") != signature:
            value = _runtime_snapshot_id_from_manifest(PROJECT_ROOT)
            manifest_source = _runtime_manifest_source_tree_sha256(PROJECT_ROOT)
            if (
                value
                and RUNTIME_SOURCE_TREE_SHA256
                and manifest_source == RUNTIME_SOURCE_TREE_SHA256
            ):
                _runtime_snapshot_cache["value"] = value
            _runtime_snapshot_cache["signature"] = signature
        return _runtime_snapshot_cache.get("value") or RUNTIME_SNAPSHOT_ID


def _flagship_artifact_is_current(path):
    if not _base_flagship_artifact_is_current(path):
        return False
    try:
        return sql_artifact_is_current(path, project_root=PROJECT_ROOT)
    except (OSError, ValueError, KeyError, TypeError):
        return False


def _v15_ce_evidence_retriever() -> SQLiteCandidateRecallRetriever:
    """Bind V15 CE and Top40 recall to the read-only Full library snapshot."""

    manifest_path = PROJECT_ROOT / "config" / "runtime_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    binding = ((manifest.get("library") or {}).get("production_full") or {})
    snapshot_id = str(binding.get("snapshot_id") or "").strip()
    declared_path = Path(str(binding.get("path") or "")).resolve()
    if not snapshot_id:
        raise ValueError("Runtime manifest has no Full-library snapshot id")
    if declared_path != DB_PATH.resolve():
        raise ValueError("V15 CE evidence library differs from the runtime library")
    return create_verified_sql_retriever(
        DB_PATH,
        library_snapshot_id=snapshot_id,
        policy=load_v15_candidate_recall_policy(),
        project_root=PROJECT_ROOT,
    )


def _materialize_current_v15_flagship(
    task_id: str,
    task: Mapping[str, Any],
    source_artifact: Path,
    output_artifact: Path,
) -> dict[str, Any]:
    """Materialize the active V15 artifact through a verified safe path."""

    if not source_artifact.is_file():
        raise FileNotFoundError("Saved scan-level candidate evidence is unavailable")
    if output_artifact.is_file():
        try:
            if not sql_artifact_is_current(output_artifact, project_root=PROJECT_ROOT):
                raise ValueError('Saved CE/recall evidence uses an older SQL metadata revision')
            return _upgrade_flagship_artifact(
                output_artifact,
                task_id=task_id,
                runtime_snapshot_id=_current_runtime_snapshot_id(),
            )
        except (OSError, TypeError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
            # Only the 1.2 -> 1.3 policy extension is eligible for a no-search
            # upgrade.  Every incompatible or unverifiable artifact falls
            # through to the ordinary full-library rebuild from immutable
            # scan evidence.
            logger.warning(
                "V15 artifact %s is not eligible for verified identity-scope "
                "upgrade; rebuilding from scan evidence: %s",
                output_artifact.name,
                exc,
            )
    v11_artifact = _v11_flagship_artifact_path(TASK_DIR, task_id)
    if not v11_artifact.is_file():
        _write_v11_flagship_artifact(
            task_id,
            _iter_legacy_result_artifact(source_artifact),
            v11_artifact,
            runtime_snapshot_id=_current_runtime_snapshot_id(),
        )
    retriever = _v15_ce_evidence_retriever()
    return _write_flagship_artifact(
        task_id,
        _iter_legacy_result_artifact(source_artifact),
        output_artifact,
        runtime_snapshot_id=_current_runtime_snapshot_id(),
        ce_evidence_retriever=retriever,
        candidate_recall_retriever=retriever,
        detection_method_id=str(
            ((task.get("params") or {}).get("rt_method_id")) or "unspecified"
        ),
    )


def _current_runtime_manifest_status() -> dict[str, Any]:
    """Return a cached source/catalog freshness check for health responses."""

    manifest_path = PROJECT_ROOT / "config" / "runtime_manifest.json"
    try:
        stat = manifest_path.stat()
        manifest_signature = (stat.st_mtime_ns, stat.st_size, stat.st_ino)
    except OSError:
        manifest_signature = None
    try:
        source_signature = _runtime_manifest_dependency_signature(PROJECT_ROOT)
    except OSError:
        source_signature = "unavailable"
    signature = (manifest_signature, source_signature)
    with _runtime_manifest_status_check_lock:
        if _runtime_manifest_status_cache.get("signature") != signature:
            _runtime_manifest_status_cache["value"] = _runtime_manifest_freshness(
                PROJECT_ROOT,
                manifest_path,
            )
            _runtime_manifest_status_cache["signature"] = signature
        value = _runtime_manifest_status_cache.get("value")
        status = dict(value) if isinstance(value, dict) else {
            "status": "unavailable",
            "fresh": False,
        }
        manifest_source = str(
            status.get("manifest_source_tree_sha256") or ""
        ).strip()
        process_matches = bool(
            RUNTIME_SOURCE_TREE_SHA256
            and manifest_source
            and manifest_source == RUNTIME_SOURCE_TREE_SHA256
        )
        current_disk_source = str(
            status.get("current_source_tree_sha256") or ""
        ).strip()
        process_matches_current_disk = bool(
            RUNTIME_SOURCE_TREE_SHA256
            and current_disk_source
            and current_disk_source == RUNTIME_SOURCE_TREE_SHA256
        )
        status.update({
            "loaded_process_source_tree_sha256": RUNTIME_SOURCE_TREE_SHA256,
            "process_source_matches_manifest": process_matches,
            "process_source_matches_current_disk": process_matches_current_disk,
            "restart_required": not process_matches_current_disk,
            "loaded_runtime_snapshot_id": _current_runtime_snapshot_id(),
        })
        if not process_matches_current_disk:
            status["release_ready"] = False
            status["status"] = "restart_required"
        return status


def _resolve_blank_file_token(value: str) -> Path:
    """Resolve a blank upload token, never an arbitrary filesystem path."""
    raw = str(value or "").strip()
    if not raw or "\x00" in raw:
        raise ValueError("missing blank upload token")
    root = BLANK_DIR.resolve()
    supplied = Path(raw)
    # Backward compatibility for already-open browser pages that received the
    # former absolute response, but only if the path is still inside BLANK_DIR.
    candidate = supplied.resolve() if supplied.is_absolute() else (root / supplied.name).resolve()
    if candidate.parent != root:
        raise ValueError("blank upload token is outside the managed blank directory")
    if not candidate.is_file():
        raise ValueError("uploaded blank file no longer exists")
    return candidate


_blank_uploads_file_lock = threading.RLock()
_BLANK_UPLOAD_EXTENSIONS = frozenset({".mgf", ".mzml", ".mzxml", ".xml"})
_BLANK_ROLES = frozenset({
    "unspecified",
    "background_blank",
    "solvent_blank",
    "method_blank",
    "extraction_blank",
    "carryover_blank",
})


def _normalize_blank_role(value: Any) -> str:
    role = str(value or "unspecified").strip().lower()
    return role if role in _BLANK_ROLES else "unspecified"


def _blank_original_filename(token: str) -> str:
    """Recover the user-facing basename from a managed UUID-prefixed token."""
    prefix, separator, remainder = str(token or "").partition("_")
    if (
        separator
        and len(prefix) == 8
        and all(char in "0123456789abcdefABCDEF" for char in prefix)
        and remainder
    ):
        return remainder
    return Path(str(token or "")).name


def _file_content_sha256(path: Path) -> str:
    """Return a byte-level identity for one managed analytical input."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _valid_content_sha256(value: Any) -> str:
    digest = str(value or "").strip().lower()
    return digest if re.fullmatch(r"[0-9a-f]{64}", digest) else ""


def _read_blank_upload_index_unlocked() -> dict[str, dict[str, Any]]:
    if not BLANK_UPLOADS_FILE.is_file():
        return {}
    try:
        payload = json.loads(BLANK_UPLOADS_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        logger.warning("Blank upload index is unreadable; rebuilding from files")
        return {}
    records = payload.get("uploads", payload) if isinstance(payload, dict) else {}
    if not isinstance(records, dict):
        return {}
    return {
        str(token): dict(record)
        for token, record in records.items()
        if isinstance(record, dict)
    }


def _write_blank_upload_index_unlocked(
    records: dict[str, dict[str, Any]],
) -> None:
    BLANK_UPLOADS_FILE.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = BLANK_UPLOADS_FILE.with_name(
        f"{BLANK_UPLOADS_FILE.name}.{os.getpid()}.{threading.get_ident()}.tmp"
    )
    try:
        tmp_path.write_text(
            json.dumps(
                {"schema_version": 2, "uploads": records},
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        tmp_path.replace(BLANK_UPLOADS_FILE)
    finally:
        tmp_path.unlink(missing_ok=True)


def _save_blank_upload_metadata(record: dict[str, Any]) -> None:
    token = Path(str(record.get("blank_path") or "")).name
    if not token:
        return
    with _blank_uploads_file_lock:
        records = _read_blank_upload_index_unlocked()
        records[token] = dict(record)
        _write_blank_upload_index_unlocked(records)


def _register_managed_blank_file(
    path: Path,
    display_name: str,
    *,
    blank_role: str = "background_blank",
    extra_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Validate a managed blank and persist the same metadata for every intake."""

    from matching.blank_subtraction import load_blank_file

    profile = load_blank_file(str(path))
    if (
        not getattr(profile, "context_spectra", None)
        and not getattr(profile, "ms1_scan_count", 0)
    ):
        raise ValueError(
            "blank contains no usable MS1 features or precursor-context MS/MS spectra"
        )
    stat = path.stat()
    record: dict[str, Any] = {
        "blank_path": path.name,
        "blank_filename": display_name,
        "content_sha256": _file_content_sha256(path),
        "size_mb": round(stat.st_size / (1024 * 1024), 2),
        "n_peaks": int(len(profile.mz_array)),
        "n_ms1_scans": int(getattr(profile, "ms1_scan_count", 0)),
        "n_ms2_contexts": len(
            getattr(profile, "context_spectra", ()) or ()
        ),
        "polarities": sorted(
            getattr(profile, "available_polarities", ()) or ()
        ),
        "blank_role": _normalize_blank_role(blank_role),
        "uploaded_at": datetime.now(timezone.utc).isoformat(),
        "metadata_status": "parsed",
        "active_task_ids": [],
        "in_use": False,
    }
    if extra_metadata:
        record.update(extra_metadata)
    _save_blank_upload_metadata(record)
    return record


def _remove_blank_upload_metadata(token: str) -> None:
    token = Path(str(token or "")).name
    with _blank_uploads_file_lock:
        records = _read_blank_upload_index_unlocked()
        if records.pop(token, None) is not None:
            _write_blank_upload_index_unlocked(records)


def _active_blank_task_references() -> dict[str, list[str]]:
    """Return queued/running task references keyed by managed blank token."""
    references: dict[str, list[str]] = defaultdict(list)
    with _tasks_lock:
        for task_id, task in _tasks.items():
            if not _is_active_task_state(task.get("status")):
                continue
            params = task.get("params") or {}
            tokens = list(params.get("blank_files") or [])
            legacy_token = params.get("blank_file")
            if legacy_token:
                tokens.append(legacy_token)
            for value in tokens:
                token = Path(str(value or "")).name
                if token and task_id not in references[token]:
                    references[token].append(task_id)
    return references


def _recover_legacy_blank_metadata(path: Path) -> dict[str, Any]:
    """Inspect a pre-index blank once so polarity switching remains reliable."""
    stat = path.stat()
    record: dict[str, Any] = {
        "blank_path": path.name,
        "blank_filename": _blank_original_filename(path.name),
        "content_sha256": _file_content_sha256(path),
        "size_mb": round(stat.st_size / (1024 * 1024), 2),
        "n_peaks": None,
        "n_ms1_scans": None,
        "n_ms2_contexts": None,
        "polarities": [],
        "blank_role": "unspecified",
        "uploaded_at": datetime.fromtimestamp(
            stat.st_mtime, tz=timezone.utc
        ).isoformat(),
        "metadata_status": "legacy_discovered",
    }
    try:
        from matching.blank_subtraction import load_blank_file

        profile = load_blank_file(path)
        record.update({
            "n_peaks": int(len(profile.mz_array)),
            "n_ms1_scans": int(getattr(profile, "ms1_scan_count", 0)),
            "n_ms2_contexts": len(
                getattr(profile, "context_spectra", ()) or ()
            ),
            "polarities": sorted(
                getattr(profile, "available_polarities", ()) or ()
            ),
            "metadata_status": "legacy_recovered",
        })
    except Exception as exc:
        record["metadata_warning"] = (
            f"Legacy blank metadata could not be recovered: "
            f"{str(exc)[:240]}"
        )
    return record


def _list_managed_blank_uploads() -> list[dict[str, Any]]:
    """List editable blank uploads, including legacy files missing an index."""
    with _blank_uploads_file_lock:
        records = _read_blank_upload_index_unlocked()
        files = {
            path.name: path
            for path in BLANK_DIR.iterdir()
            if path.is_file()
            and not path.name.startswith(".")
            and path.suffix.lower() in _BLANK_UPLOAD_EXTENSIONS
        }
        changed = False
        for token in list(records):
            if token not in files:
                records.pop(token, None)
                changed = True
            elif records[token].get("blank_role") not in _BLANK_ROLES:
                records[token]["blank_role"] = "unspecified"
                changed = True
            if token in files and not _valid_content_sha256(
                records[token].get("content_sha256")
            ):
                records[token]["content_sha256"] = _file_content_sha256(
                    files[token]
                )
                changed = True
        for token, path in files.items():
            if token in records:
                continue
            records[token] = _recover_legacy_blank_metadata(path)
            changed = True
        if changed:
            _write_blank_upload_index_unlocked(records)

    references = _active_blank_task_references()
    output = []
    for token, record in records.items():
        if token not in files:
            continue
        item = dict(record)
        item["blank_path"] = token
        item.setdefault("blank_filename", _blank_original_filename(token))
        item.setdefault(
            "size_mb", round(files[token].stat().st_size / (1024 * 1024), 2)
        )
        item.setdefault("polarities", [])
        item["blank_role"] = _normalize_blank_role(item.get("blank_role"))
        item["active_task_ids"] = references.get(token, [])
        item["in_use"] = bool(item["active_task_ids"])
        output.append(item)
    output.sort(
        key=lambda item: (
            str(item.get("uploaded_at") or ""),
            str(item.get("blank_filename") or ""),
        ),
        reverse=True,
    )
    return output


def _deduplicate_managed_blank_tokens(
    values: list[str],
) -> tuple[list[str], list[dict[str, str]]]:
    """Collapse byte-identical blanks before they can become replicates."""
    tokens = list(dict.fromkeys(
        Path(str(value or "")).name for value in values
        if str(value or "").strip()
    ))
    with _blank_uploads_file_lock:
        records = _read_blank_upload_index_unlocked()
        changed = False
        identities: dict[str, str] = {}
        unique: list[str] = []
        duplicates: list[dict[str, str]] = []
        for token in tokens:
            path = _resolve_blank_file_token(token)
            record = dict(records.get(token) or {})
            digest = _valid_content_sha256(record.get("content_sha256"))
            if not digest:
                digest = _file_content_sha256(path)
                record.update({
                    "blank_path": token,
                    "blank_filename": record.get("blank_filename")
                    or _blank_original_filename(token),
                    "content_sha256": digest,
                })
                records[token] = record
                changed = True
            retained = identities.get(digest)
            if retained:
                duplicates.append({
                    "blank_path": token,
                    "duplicate_of": retained,
                    "content_sha256": digest,
                })
                continue
            identities[digest] = token
            unique.append(token)
        if changed:
            _write_blank_upload_index_unlocked(records)
    return unique, duplicates

# ── Upload Limits ───────────────────────────────────────────────────
MAX_UPLOAD_SIZE_MB = 500          # 单文件最大 500 MB
MAX_UPLOAD_SIZE_BYTES = MAX_UPLOAD_SIZE_MB * 1024 * 1024
MAX_BATCH_UPLOAD_FILES = 100
MAX_BATCH_UPLOAD_SIZE_MB = 5_000
MAX_BATCH_UPLOAD_SIZE_BYTES = MAX_BATCH_UPLOAD_SIZE_MB * 1024 * 1024
MAX_VENDOR_RAW_FILES = 2_000
MAX_VENDOR_RAW_SIZE_MB = 20_000
MAX_VENDOR_RAW_SIZE_BYTES = MAX_VENDOR_RAW_SIZE_MB * 1024 * 1024
SAMPLE_UPLOAD_EXTENSIONS = frozenset({
    ".mgf", ".mzml", ".mzxml", ".xml", ".html", ".htm",
})


def _require_disk_budget(path: Path, planned_bytes: int, operation: str) -> None:
    try:
        require_storage_capacity(path, planned_bytes, operation)
    except InsufficientStorageError as exc:
        raise HTTPException(status_code=507, detail=str(exc)) from exc

# ── Global State ────────────────────────────────────────────────────
_lib_records: list[dict] = []
_mz_index: dict[int, list[int]] = {}
_lib_loaded: bool = False
_lib_name: str = ""
_lib_spectra_count: int = 0
_lib_load_time: float = 0.0
_library_state_lock = threading.RLock()
_library_count_cache_lock = threading.Lock()
_library_count_compute_lock = threading.Lock()
_library_count_cache: dict[str, tuple[tuple[int, int, bool], dict[str, Any]]] = {}

_tasks: dict[str, dict] = {}
_tasks_lock = threading.Lock()
_tasks_file_lock = threading.Lock()
# Presentation migrations may need to read a large persisted result artifact.
# Keep those migrations serial and, critically, never retain the artifact in
# the task registry.  The registry is a metadata index; result rows live on
# disk and are loaded only for the lifetime of one request/worker.
_presentation_refresh_lock = threading.RLock()
_task_cancel_events: dict[str, threading.Event] = {}
_restart_recovery_task_ids: set[str] = set()
_task_delete_preview_lock = threading.Lock()
_task_delete_previews: dict[str, dict[str, Any]] = {}
_task_deletion_audit_lock = threading.Lock()
_task_recycle_lock = threading.RLock()
_task_recycle_purge_preview_lock = threading.Lock()
_task_recycle_purge_previews: dict[str, dict[str, Any]] = {}
# FULL 库单任务会持有大量候选/结构推理对象；并发两份实测可把进程推到
# ~16 GB RSS 并耗尽交换区。FULL 仍严格串行，但等待中的任务可由操作者
# 置顶/上移/下移；正在运行的任务永远不会被重排或抢占。
_executor_full = ReorderableSerialExecutor(thread_name_prefix="ms-full")
_executor_small = ThreadPoolExecutor(max_workers=4, thread_name_prefix="ms-small")
# PubMed enrichment is optional I/O-only post-processing.  It has its own
# serial worker and NCBI rate limiter, so a slow or unavailable literature
# service can never occupy the Full/Small analysis executors.
_executor_literature = ThreadPoolExecutor(
    max_workers=1,
    thread_name_prefix="ms-literature",
)
# EvoGraph-F reads only the committed legacy result and writes a separate
# shadow sidecar.  Its serial worker cannot occupy or delay a search executor.
_executor_evograph = ThreadPoolExecutor(
    max_workers=1,
    thread_name_prefix="ms-evograph",
)
# Report generation can allocate large pandas/openpyxl/reportlab objects.  All
# download path operations are synchronous so FastAPI runs them outside the
# event loop, while this re-entrant gate prevents multiple large exports from
# multiplying memory pressure.  Health/status/cancel APIs remain responsive.
_report_export_lock = threading.RLock()


def _full_queue_snapshot() -> dict[str, Any]:
    """Read the Full executor queue without depending on private internals."""
    snapshot = getattr(_executor_full, "snapshot", None)
    if not callable(snapshot):
        return {"running": None, "queued": [], "reorder_supported": False}
    try:
        payload = dict(snapshot() or {})
    except Exception as exc:
        logger.warning("Could not inspect Full analysis queue: %s", exc)
        return {"running": None, "queued": [], "reorder_supported": False}
    return {
        "running": payload.get("running"),
        "queued": [str(item) for item in payload.get("queued") or []],
        "reorder_supported": True,
    }


def _sync_full_queue_order_metadata(
    snapshot: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Persist restart-stable ordering for active Full-library work."""
    payload = dict(snapshot or _full_queue_snapshot())
    ordered = []
    running = str(payload.get("running") or "").strip()
    if running:
        ordered.append(running)
    ordered.extend(
        task_id for task_id in payload.get("queued") or []
        if task_id not in ordered
    )
    with _tasks_lock:
        for queue_order, task_id in enumerate(ordered):
            task = _tasks.get(task_id)
            if task is not None:
                task["queue_order"] = queue_order
    return payload


def _full_queue_public_payload(
    snapshot: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Project the executor snapshot into stable API rows."""
    payload = dict(snapshot or _full_queue_snapshot())
    waiting_ids = list(payload.get("queued") or [])
    running_id = str(payload.get("running") or "").strip() or None
    with _tasks_lock:
        def row(task_id: str, *, position: Optional[int] = None):
            task = _tasks.get(task_id) or {}
            return {
                "id": task_id,
                "filename": task.get("filename"),
                "status": task.get("status"),
                "position": position,
                "queue_order": task.get("queue_order"),
            }

        running_row = row(running_id) if running_id else None
        queued_rows = [
            row(task_id, position=index)
            for index, task_id in enumerate(waiting_ids, start=1)
        ]
    return {
        "reorder_supported": bool(payload.get("reorder_supported")),
        "running": running_row,
        "queued": queued_rows,
        "queued_count": len(queued_rows),
        "order": waiting_ids,
    }


def _remove_waiting_full_queue_item(task_id: str) -> bool:
    """Discard a deleted Full task if its callable has not started yet."""
    cancel = getattr(_executor_full, "cancel", None)
    if not callable(cancel) or not cancel(task_id):
        return False
    snapshot = _sync_full_queue_order_metadata()
    if _save_tasks() is False:
        logger.error(
            "Full queue item %s was removed but remaining order was not persisted",
            task_id,
        )
    public = _full_queue_public_payload(snapshot)
    _sse_broadcast(
        "queue_reordered",
        {
            "task_id": task_id,
            "action": "remove",
            "order": public["order"],
            "queued_count": public["queued_count"],
        },
    )
    return True


def _serialized_report_export(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        from web.export_jobs import has_progress_callback, report_progress
        # A legacy HTTP request must not silently queue for minutes behind
        # another export. Background jobs can wait with an explicit status.
        acquired = _report_export_lock.acquire(blocking=False)
        if not acquired and not has_progress_callback():
            raise HTTPException(status_code=409, detail="已有导出正在生成，请等待完成或打开导出进度页；没有启动重复导出。")
        try:
            while not acquired:
                report_progress(status="queued", phase="waiting_export_slot")
                acquired = _report_export_lock.acquire(timeout=1)
            report_progress(status="running", phase="preparing")
            return function(*args, **kwargs)
        finally:
            if acquired:
                _report_export_lock.release()

    return wrapped


def _build_excel_export_subprocess(task_id: str) -> bytes:
    """Render a large workbook outside the API process and its Python GIL."""
    if not re.fullmatch(r"[A-Za-z0-9_-]{1,128}", str(task_id or "")):
        raise ValueError("Invalid task identifier for export worker")
    # Test fixtures and a just-created in-memory task may intentionally have no
    # persisted metadata/result pair yet.  Keep that compatibility path local;
    # every normal completed production task uses the isolated worker below.
    try:
        persisted_index = json.loads(TASKS_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        persisted_index = {}
    persisted_result = TASK_DIR / f"{task_id}_results.json"
    if task_id not in persisted_index or not persisted_result.is_file():
        with _tasks_lock:
            task = _tasks.get(task_id)
        if not isinstance(task, dict):
            raise KeyError(f"Task not found: {task_id}")
        return _build_analysis_excel_bytes(
            task_id,
            task,
            _get_task_results(task_id),
            language="zh",
        )
    export_tmp_root = OUTPUT_DIR / "export_tmp"
    export_tmp_root.mkdir(parents=True, exist_ok=True)
    worker_script = PROJECT_ROOT / "scripts" / "render_task_excel.py"
    with tempfile.TemporaryDirectory(
        prefix=f"excel-{task_id[:16]}-",
        dir=export_tmp_root,
    ) as temporary_dir:
        output_path = Path(temporary_dir) / "analysis.xlsx"
        environment = dict(os.environ)
        environment["MSSYSTEM_RUNTIME_ROOT"] = str(OUTPUT_DIR)
        completed = subprocess.run(
            [
                sys.executable,
                str(worker_script),
                "--task-id",
                task_id,
                "--output",
                str(output_path),
            ],
            cwd=PROJECT_ROOT,
            env=environment,
            capture_output=True,
            text=True,
            timeout=30 * 60,
        )
        if completed.returncode != 0 or not output_path.is_file():
            diagnostic = (completed.stderr or completed.stdout or "")[-2000:]
            raise RuntimeError(
                "Excel export worker failed"
                + (f": {diagnostic}" if diagnostic else "")
            )
        return output_path.read_bytes()


def _invalidate_library_state() -> None:
    """Atomically invalidate both Web and canonical-search library caches."""
    global _lib_records, _mz_index, _lib_loaded, _lib_name
    global _lib_spectra_count, _lib_load_time

    search_module = importlib.import_module("matching.search")

    with _library_state_lock:
        _lib_records = []
        _mz_index = {}
        _lib_loaded = False
        _lib_name = ""
        _lib_spectra_count = 0
        _lib_load_time = 0.0
        # No analysis may be queued/running while this helper is called, so it
        # is safe to release the canonical in-memory arrays as one operation.
        with search_module._lib_lock:
            search_module._lib_records = None
            search_module._mz_index = None
            search_module._nm_index = None
            search_module._lib_db_path = None

# ── SSE (Server-Sent Events) Support ─────────────────────────────────
import asyncio
from typing import Callable

_sse_subscribers: list[asyncio.Queue] = []
_sse_lock = threading.Lock()


def _sse_broadcast(event_type: str, data: dict):
    """Broadcast an event to all connected SSE clients."""
    import json as _json
    payload = _json.dumps({"type": event_type, **data}, ensure_ascii=False, default=str)
    with _sse_lock:
        dead = []
        for i, q in enumerate(_sse_subscribers):
            try:
                q.put_nowait(payload)
            except Exception:
                dead.append(i)
        for i in reversed(dead):
            _sse_subscribers.pop(i)

# ── Task Persistence ────────────────────────────────────────────────
def _repair_migrated_upload_path(file_path: str) -> str:
    """Reconnect a cloned task path to this checkout's managed upload copy.

    A version fork may leave both the legacy and the cloned upload files on
    disk.  Testing only ``old_path.exists()`` therefore keeps V12 coupled to
    V11.  Managed ``output/tmp_uploads`` paths are always rebound by basename
    when the local clone exists; arbitrary external input paths are preserved.
    """
    old_path = Path(file_path or "")
    migrated_path = UPLOAD_DIR / old_path.name
    if old_path.name and migrated_path.is_file():
        try:
            old_resolved = old_path.resolve()
            upload_resolved = UPLOAD_DIR.resolve()
            if old_resolved.is_relative_to(upload_resolved):
                return str(old_resolved)
        except (OSError, RuntimeError):
            pass
        if (
            not old_path.exists()
            or (old_path.parent.name == "tmp_uploads" and old_path.parent.parent.name == "output")
        ):
            return str(migrated_path)
    return str(file_path or "")


def _task_index_backup_path() -> Path:
    return TASKS_FILE.with_name(f"{TASKS_FILE.stem}.backup.json")


def _read_task_index(path: Path) -> dict[str, dict]:
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("task index root must be an object")
    invalid = [task_id for task_id, task in payload.items() if not isinstance(task, dict)]
    if invalid:
        raise ValueError(f"task index contains invalid task rows: {invalid[:3]}")
    return payload


def _load_tasks():
    """Load durable tasks and preserve recoverable work across restarts."""
    global _tasks
    _restart_recovery_task_ids.clear()
    primary_error: Exception | None = None
    if TASKS_FILE.exists():
        try:
            _tasks = _read_task_index(TASKS_FILE)
        except Exception as exc:
            primary_error = exc
            logger.error("Primary task index is unreadable: %s", exc)
    else:
        _tasks = {}

    if primary_error is not None:
        backup_path = _task_index_backup_path()
        try:
            _tasks = _read_task_index(backup_path)
        except Exception as backup_exc:
            # Starting with an empty registry would make the next harmless
            # save erase the only evidence that tasks existed.  Refuse startup
            # instead and leave both files untouched for recovery.
            raise RuntimeError(
                "Task registry is unreadable and no valid backup is available; "
                "startup aborted without modifying task data"
            ) from backup_exc
        logger.warning(
            "Recovered %d tasks from backup after primary-index failure",
            len(_tasks),
        )

    for tid, task in _tasks.items():
        _initialize_task_state(task)
        task.setdefault("results", [])
        old_file_path = task.get("file_path") or ""
        repaired_file_path = _repair_migrated_upload_path(old_file_path)
        if repaired_file_path != old_file_path:
            task["file_path"] = repaired_file_path
            logger.info(
                "Repaired migrated upload path for task %s: %s",
                tid, repaired_file_path,
            )

        previous_state = str(task.get("status") or "")
        if previous_state == TaskState.CANCELLING.value:
            _transition_task(
                task,
                TaskState.CANCELLED,
                reason="server_restart_completed_requested_cancellation",
            )
            task["error"] = None
            task["message"] = "停止请求已在服务重启时完成；原始上传仍保留"
            continue
        if previous_state not in {
            TaskState.QUEUED.value,
            TaskState.RUNNING.value,
            TaskState.FINALIZING.value,
        }:
            continue

        upload_path = Path(str(task.get("file_path") or ""))
        if not upload_path.is_file():
            _transition_task(
                task,
                TaskState.FAILED,
                reason="server_restart_recovery_input_missing",
            )
            task["error"] = "原始上传文件缺失，无法在重启后恢复分析"
            task["message"] = task["error"]
            continue

        _transition_task(
            task,
            TaskState.PENDING,
            reason="server_restart_preserved_for_automatic_resume",
        )
        requested_mode = str((task.get("params") or {}).get("db_mode") or "full")
        task["restart_recovery"] = {
            "status": "awaiting_library",
            "previous_state": previous_state,
            "requested_db_mode": requested_mode,
            "detected_at": _now_iso(),
            "automatic_resume": True,
        }
        task["error"] = None
        task["message"] = (
            f"服务重启后任务与原始文件已保留；加载 {requested_mode.upper()} 库后自动恢复排队"
        )
        _restart_recovery_task_ids.add(tid)

    if _frozen_task_bridge is not None:
        added = _frozen_task_bridge.install(_tasks)
        logger.info("Registered %d completed frozen XBT tasks without analysis", added)

    logger.info(
        "Loaded %d tasks from disk (%d awaiting automatic restart recovery)",
        len(_tasks), len(_restart_recovery_task_ids),
    )
    if _save_tasks() is False:
        raise RuntimeError("Could not persist recovered task registry")


def _tasks_to_metadata() -> dict:
    """Extract task metadata (without large results) for persistence."""
    meta = {}
    with _tasks_lock:
        for tid, task in _tasks.items():
            t = {k: v for k, v in task.items() if k != "results"}
            # Only keep a flag indicating results exist on disk
            t["has_results_file"] = bool(task.get("completed_batch")) or bool(task.get("results")) or (
                TASK_DIR / f"{tid}_results.json"
            ).exists()
            meta[tid] = t
    return meta


def _save_tasks(
    *, allow_empty: bool = False, backup_previous: bool = True
) -> bool:
    """Atomically save task metadata while retaining a durable recovery copy."""
    tmp_path: Path | None = None
    backup_tmp_path: Path | None = None
    try:
        with _tasks_file_lock:
            meta = _tasks_to_metadata()
            if not meta and TASKS_FILE.is_file() and not allow_empty:
                try:
                    existing = _read_task_index(TASKS_FILE)
                except Exception:
                    existing = {}
                if existing:
                    raise RuntimeError(
                        "refusing to replace a non-empty task registry with an "
                        "empty registry outside an explicit deletion transaction"
                    )
            tmp_path = TASKS_FILE.with_name(
                f"{TASKS_FILE.name}.{os.getpid()}.{threading.get_ident()}.tmp"
            )
            with open(tmp_path, "w", encoding="utf-8") as f:
                json.dump(meta, f, ensure_ascii=False, indent=2, default=str)
                f.flush()
                os.fsync(f.fileno())

            # Preserve the previous valid registry before publishing the new
            # generation.  Corrupt primaries are never copied over a good backup.
            if backup_previous and TASKS_FILE.is_file():
                try:
                    previous = _read_task_index(TASKS_FILE)
                except Exception:
                    previous = None
                if previous is not None:
                    backup_path = _task_index_backup_path()
                    backup_tmp_path = backup_path.with_name(
                        f"{backup_path.name}.{os.getpid()}.{threading.get_ident()}.tmp"
                    )
                    with open(backup_tmp_path, "w", encoding="utf-8") as backup:
                        json.dump(
                            previous, backup, ensure_ascii=False, indent=2,
                            default=str,
                        )
                        backup.flush()
                        os.fsync(backup.fileno())
                    backup_tmp_path.replace(backup_path)
                    backup_tmp_path = None
            tmp_path.replace(TASKS_FILE)
            tmp_path = None
            return True
    except Exception as e:
        logger.error("Failed to save tasks: %s", e)
        for partial in (tmp_path, backup_tmp_path):
            if partial:
                try:
                    partial.unlink(missing_ok=True)
                except Exception:
                    pass
        return False


_TASK_PROGRESS_DURABLE_INTERVAL_SECONDS = 10.0
_task_progress_persist_lock = threading.Lock()
_task_progress_last_persisted: dict[str, float] = {}


def _save_task_progress_throttled(
    task_id: str,
    *,
    force: bool = False,
    now_monotonic: float | None = None,
) -> bool:
    """Persist volatile progress at most once per interval.

    SSE remains the real-time channel.  Progress-only saves skip rotating the
    recovery backup, while every state transition and final result continues
    to call :func:`_save_tasks` directly.
    """

    now = time.monotonic() if now_monotonic is None else float(now_monotonic)
    with _task_progress_persist_lock:
        previous = _task_progress_last_persisted.get(task_id)
        if (
            not force
            and previous is not None
            and now - previous < _TASK_PROGRESS_DURABLE_INTERVAL_SECONDS
        ):
            return True
        saved = _save_tasks(backup_previous=False)
        if saved:
            _task_progress_last_persisted[task_id] = now
        return saved


# ── Pydantic Models ─────────────────────────────────────────────────
class LibraryStatus(BaseModel):
    loaded: bool
    name: str
    spectra_count: int
    load_time_sec: float


class TaskSummary(BaseModel):
    id: str
    filename: str
    status: str
    progress: int
    message: str
    upload_time: str
    total_spectra: int
    matched_spectra: int
    match_rate: float
    unique_compounds: Optional[int] = None
    contaminant_unique_compounds: Optional[int] = None
    ion_identity_networks: Optional[int] = None
    deconvolved_spectra: Optional[int] = None


class AnalysisParams(BaseModel):
    db_mode: Literal["small", "full"] = "full"
    cosine_threshold: float = Field(COSINE_THRESHOLD, ge=0.0, le=1.0)
    final_score_threshold: float = Field(FINAL_SCORE_THRESHOLD, ge=0.0, le=1.0)
    mz_tolerance_da: float = Field(MZ_TOLERANCE_DA, gt=0.0, le=1.0)
    mz_tolerance_ppm: float = Field(MZ_TOLERANCE_PPM, gt=0.0, le=1000.0)
    top_n_peaks: int = Field(TOP_N_PEAKS, ge=3, le=500)
    peak_processing_preset: Literal[
        "library_match", "sensitive", "network_rich"
    ] = "library_match"
    peak_relative_intensity_threshold: Optional[float] = Field(
        None, ge=0.0, le=1.0
    )
    peak_precursor_exclusion_da: float = Field(0.5, ge=0.0, le=10.0)
    peak_top_n: Optional[int] = Field(None, ge=3, le=500)
    peak_upper_mass_bound_enabled: bool = True
    # mzML spectrum representation is evaluated per MS1/MS2 spectrum.  Profile
    # spectra use vendor/OpenMS peak picking; the old spacing heuristic is not a
    # production backend.
    profile_peak_picking: Literal[
        "auto", "openms", "reject", "assume_centroid"
    ] = "auto"
    profile_peak_picker: Literal["hires", "iterative"] = "hires"
    profile_peak_picker_signal_to_noise: float = Field(1.0, ge=0.0, le=100.0)
    enable_ms1_feature_detection: bool = True
    ms1_mass_error_ppm: float = Field(10.0, gt=0.0, le=100.0)
    # The mzML selected-ion m/z is an isolation setpoint, not always a measured
    # centroid.  This separate, RT/charge-constrained window links it to a
    # continuous-MS1 feature; library search itself still obeys mz_tolerance_ppm.
    ms2_feature_link_tolerance_ppm: float = Field(
        20.0, gt=0.0, le=100.0
    )
    ms1_noise_threshold_intensity: float = Field(1000.0, ge=0.0)
    ms1_chrom_peak_snr: float = Field(3.0, ge=0.0, le=100.0)
    ms1_min_trace_length_seconds: float = Field(5.0, gt=0.0, le=300.0)
    ms1_min_fwhm_seconds: float = Field(1.0, gt=0.0, le=300.0)
    ms1_max_fwhm_seconds: float = Field(60.0, gt=0.0, le=1800.0)
    ms1_remove_single_traces: bool = False
    ms1_charge_low: int = Field(1, ge=1, le=10)
    ms1_charge_high: int = Field(1, ge=1, le=10)
    enable_adduct_grouping: bool = True
    ms1_max_adduct_grouping_features: int = Field(
        10_000, ge=0, le=100_000
    )
    ms1_feature_timeout_seconds: int = Field(300, ge=30, le=3600)
    strict_evidence_gate: bool = USE_STRICT_EVIDENCE_GATE
    enable_chimeric_deconvolution: bool = USE_CHIMERIC_DECONVOLUTION
    enable_ion_identity_networking: bool = USE_ION_IDENTITY_NETWORKING
    # Compound Discoverer-style Expected Compounds is a separate continuous-
    # MS1 target layer.  It never increments ordinary spectral-library matches.
    enable_expected_compounds: bool = False
    expected_compound_datasets: list[str] = Field(
        default_factory=list, max_length=20
    )
    expected_mass_tolerance_ppm: float = Field(5.0, gt=0.0, le=50.0)
    expected_min_peak_intensity: float = Field(1000.0, ge=0.0)
    expected_min_snr: float = Field(1.5, ge=0.0, le=1000.0)
    expected_min_scan_points: int = Field(5, ge=3, le=100)
    expected_max_gap_ratio: float = Field(0.35, ge=0.0, le=1.0)
    expected_max_fwhm_seconds: float = Field(
        60.0, gt=0.0, le=1800.0
    )
    expected_min_relative_valley_depth: float = Field(
        0.10, ge=0.0, le=1.0
    )
    expected_openms_rt_tolerance_minutes: float = Field(
        0.05, gt=0.0, le=1.0
    )
    expected_background_sample_blank_ratio: float = Field(
        5.0, ge=1.0, le=100.0
    )
    expected_generate_common_adducts: bool = True
    expected_compound_ms2_precursor_ppm: float = Field(
        10.0, gt=0.0, le=100.0
    )
    expected_compound_diagnostic_tolerance_da: float = Field(
        0.02, gt=0.0, le=0.5
    )
    # Local RT standards are comparable only inside an explicitly named LC
    # method. Blank means RT remains audit-only and cannot resolve an isomer.
    rt_method_id: Optional[str] = Field(None, max_length=128)
    # Full quadrupole isolation width.  Used only when mzML omits lower/upper
    # offsets; blank means evidence-unavailable rather than silently guessing.
    isolation_window_width_da: Optional[float] = Field(None, gt=0.0, le=20.0)
    # auto=从 mzML CV / MGF IONMODE、CHARGE 自动识别；POS/NEG=用户明确覆盖。
    # 自动模式无法可靠判断时拒绝运行，避免静默放开跨极性匹配。
    ion_mode: Literal["auto", "POS", "NEG"] = "auto"
    # ── V7 blank evidence: identify first, then annotate/filter views ──
    blank_file: Optional[str] = Field(None, max_length=1024)  # BLANK_DIR 内受控 token
    blank_files: list[str] = Field(default_factory=list, max_length=10)
    blank_roles: dict[str, Literal[
        "unspecified",
        "background_blank",
        "solvent_blank",
        "method_blank",
        "extraction_blank",
        "carryover_blank",
    ]] = Field(default_factory=dict)
    blank_filter_policy: Literal[
        "soft_evidence", "hard_legacy"
    ] = DEFAULT_BLANK_FILTER_POLICY
    # Hard background exclusion is intentionally double-opt-in. A role label
    # and matching polarity alone do not establish batch/method comparability.
    blank_method_comparability_confirmed: bool = False
    blank_background_exclusion_confirmed: bool = False
    blank_clean_max_ratio: float = Field(
        DEFAULT_BLANK_CLEAN_MAX_RATIO, ge=0.0, le=1.0,
    )
    blank_clean_min_replicates: int = Field(
        DEFAULT_BLANK_CLEAN_MIN_REPLICATES, ge=2, le=10,
    )
    blank_mode: Literal["absolute", "relative", "filter_only"] = DEFAULT_BLANK_MODE
    blank_threshold: float = Field(DEFAULT_BLANK_THRESHOLD, ge=0.0, le=100.0)
    blank_intensity_mode: Literal["max", "mean"] = DEFAULT_BLANK_INTENSITY_MODE
    # 新任务使用 LC-MS 特征级空白证据；上面三个旧碎片扣除参数仅为历史任务兼容。
    blank_sample_ratio: float = Field(DEFAULT_BLANK_SAMPLE_RATIO, ge=1.0, le=100.0)
    blank_rt_tolerance_minutes: float = Field(
        DEFAULT_BLANK_RT_TOLERANCE_MINUTES, ge=0.01, le=5.0,
    )
    blank_min_detection_fraction: float = Field(
        DEFAULT_BLANK_MIN_DETECTION_FRACTION, ge=0.0, le=1.0,
    )

    @field_validator("blank_files", mode="before")
    @classmethod
    def _normalize_blank_files(cls, value):
        if value is None:
            return []
        if isinstance(value, str):
            return [value]
        return value

    @field_validator("expected_compound_datasets", mode="before")
    @classmethod
    def _normalize_expected_compound_datasets(cls, value):
        if value is None:
            return []
        if isinstance(value, str):
            value = [value]
        return list(dict.fromkeys(
            str(item).strip() for item in value if str(item).strip()
        ))

    @field_validator("blank_roles", mode="before")
    @classmethod
    def _normalize_blank_roles(cls, value):
        if not isinstance(value, dict):
            return {}
        normalized = {}
        for raw_token, raw_role in list(value.items())[:10]:
            token = Path(str(raw_token or "")).name
            if token:
                normalized[token] = _normalize_blank_role(raw_role)
        return normalized

    @field_validator("db_mode", mode="before")
    @classmethod
    def _normalize_db_mode(cls, value):
        return str(value or "full").strip().lower()

    @field_validator("ion_mode", mode="before")
    @classmethod
    def _normalize_ion_mode(cls, value):
        text = str(value or "auto").strip().upper()
        return "auto" if text == "AUTO" else text

    @model_validator(mode="after")
    def _validate_ms1_ranges(self):
        if self.ms1_min_fwhm_seconds > self.ms1_max_fwhm_seconds:
            raise ValueError("ms1_min_fwhm_seconds cannot exceed ms1_max_fwhm_seconds")
        if self.ms1_charge_low > self.ms1_charge_high:
            raise ValueError("ms1_charge_low cannot exceed ms1_charge_high")
        selected_blanks = list(self.blank_files)
        if self.blank_file and self.blank_file not in selected_blanks:
            selected_blanks.append(self.blank_file)
        if self.blank_filter_policy == "hard_legacy" and selected_blanks:
            if not (
                self.blank_method_comparability_confirmed
                and self.blank_background_exclusion_confirmed
            ):
                raise ValueError(
                    "hard_legacy blank exclusion requires explicit confirmation "
                    "of LC-method comparability and background exclusion"
                )
            eligible_roles = {
                "background_blank", "solvent_blank", "method_blank",
                "extraction_blank",
            }
            invalid = [
                Path(token).name for token in selected_blanks
                if self.blank_roles.get(Path(token).name, "unspecified")
                not in eligible_roles
            ]
            if invalid:
                raise ValueError(
                    "hard_legacy blank exclusion accepts only explicitly typed "
                    "solvent/method/extraction blanks; invalid: "
                    + ", ".join(invalid)
                )
        return self


def _peak_config_from_params(params: dict[str, Any]):
    preset = str(params.get("peak_processing_preset") or "library_match")
    explicit_top_n = params.get("peak_top_n")
    if explicit_top_n is None and preset == "library_match":
        explicit_top_n = params.get("top_n_peaks", TOP_N_PEAKS)
    return resolve_peak_processing_config(
        preset,
        relative_intensity_threshold=params.get(
            "peak_relative_intensity_threshold"
        ),
        precursor_exclusion_da=params.get(
            "peak_precursor_exclusion_da", 0.5
        ),
        top_n=explicit_top_n,
        upper_mass_bound_enabled=params.get(
            "peak_upper_mass_bound_enabled", True
        ),
    )


class BatchAnalysisRequest(BaseModel):
    task_ids: list[str] = Field(default_factory=list)
    params: AnalysisParams = Field(default_factory=AnalysisParams)


class TaskQueueReorderRequest(BaseModel):
    task_id: str = Field(..., min_length=1, max_length=128)
    action: Literal["front", "up", "down"]

    @field_validator("task_id")
    @classmethod
    def _clean_task_id(cls, value):
        task_id = str(value or "").strip()
        if not task_id:
            raise ValueError("task_id cannot be empty")
        return task_id


class BatchTaskDeletePreviewRequest(BaseModel):
    task_ids: list[str] = Field(..., min_length=1, max_length=500)

    @field_validator("task_ids")
    @classmethod
    def _unique_task_ids(cls, values):
        unique = list(dict.fromkeys(
            str(value or "").strip()
            for value in values
            if str(value or "").strip()
        ))
        if not unique:
            raise ValueError("Select at least one task")
        return unique


class BatchTaskDeleteRequest(BatchTaskDeletePreviewRequest):
    confirmation_token: str = Field(..., min_length=16, max_length=128)
    confirmation_phrase: str = Field(..., min_length=1, max_length=128)


class TaskRecyclePurgeRequest(BaseModel):
    confirmation_token: str = Field(..., min_length=16, max_length=128)


class BlankRoleUpdate(BaseModel):
    blank_role: Literal[
        "unspecified",
        "background_blank",
        "solvent_blank",
        "method_blank",
        "extraction_blank",
        "carryover_blank",
    ]


class BlankBatchDeleteRequest(BaseModel):
    blank_tokens: list[str] = Field(
        ..., min_length=1, max_length=500
    )

    @field_validator("blank_tokens")
    @classmethod
    def _unique_blank_tokens(cls, values):
        unique = list(dict.fromkeys(
            str(value or "").strip()
            for value in values
            if str(value or "").strip()
        ))
        if not unique:
            raise ValueError("Select at least one blank file")
        return unique


class RTAlignmentRequest(BaseModel):
    task_ids: list[str] = Field(min_length=2, max_length=100)
    rt_tolerance_seconds: float = Field(20.0, gt=0.0, le=300.0)
    mz_tolerance: float = Field(10.0, gt=0.0, le=1000.0)
    mz_unit: Literal["ppm", "Da"] = "ppm"

    @field_validator("task_ids")
    @classmethod
    def _unique_task_ids(cls, value):
        unique = list(dict.fromkeys(str(item).strip() for item in value if str(item).strip()))
        if len(unique) < 2:
            raise ValueError("RT alignment requires at least two distinct tasks")
        return unique


class DualPolarityShadowRequest(BaseModel):
    positive_task_id: str = Field(..., min_length=1, max_length=128)
    negative_task_id: str = Field(..., min_length=1, max_length=128)
    same_sample_confirmed: bool
    same_lc_method_confirmed: bool
    pair_label: str = Field(default="", max_length=160)
    rt_alignment_id: str | None = Field(default=None, max_length=64)

    @model_validator(mode="after")
    def _validate_explicit_pair(self):
        self.positive_task_id = self.positive_task_id.strip()
        self.negative_task_id = self.negative_task_id.strip()
        self.pair_label = self.pair_label.strip()
        if self.positive_task_id == self.negative_task_id:
            raise ValueError("POS and NEG must be two distinct explicit task IDs")
        if not self.same_sample_confirmed:
            raise ValueError("Confirm that POS and NEG were acquired from the same sample")
        if not self.same_lc_method_confirmed:
            raise ValueError("Confirm that POS and NEG used the same LC method")
        if self.rt_alignment_id is not None:
            value = self.rt_alignment_id.strip()
            if not re.fullmatch(r"RTA-[0-9a-f]{12}", value):
                raise ValueError("RT alignment ID is invalid")
            self.rt_alignment_id = value
        return self


class EvidenceGovernanceCohortRequest(BaseModel):
    """Explicit injection cohort for EviGate-MS stability governance.

    Task IDs are never inferred from filenames.  The caller declares blank
    roles and independence groups so technical scans cannot silently masquerade
    as biological replication.
    """

    task_ids: list[str] = Field(min_length=2, max_length=12)
    task_roles: dict[str, Literal[
        "sample", "replicate", "qc", "blank"
    ]] = Field(default_factory=dict)
    independence_groups: dict[str, str] = Field(default_factory=dict)
    limit: int = Field(default=500, ge=1, le=2000)
    offset: int = Field(default=0, ge=0)

    @field_validator("task_ids")
    @classmethod
    def _unique_governance_task_ids(cls, values):
        unique = list(dict.fromkeys(
            str(value or "").strip()
            for value in values
            if str(value or "").strip()
        ))
        if len(unique) < 2:
            raise ValueError(
                "Evidence governance requires at least two distinct task IDs"
            )
        return unique

    @model_validator(mode="after")
    def _validate_governance_labels(self):
        task_ids = set(self.task_ids)
        unknown_roles = set(self.task_roles) - task_ids
        unknown_groups = set(self.independence_groups) - task_ids
        if unknown_roles or unknown_groups:
            unknown = sorted(unknown_roles | unknown_groups)
            raise ValueError(
                f"Cohort labels reference task IDs outside task_ids: {unknown}"
            )
        cleaned_groups = {}
        for task_id, value in self.independence_groups.items():
            label = " ".join(str(value or "").split())[:128]
            if not label:
                raise ValueError(
                    f"independence_groups[{task_id!r}] cannot be blank"
                )
            cleaned_groups[task_id] = label
        self.independence_groups = cleaned_groups
        return self


class IsomerLiteratureSearchRequest(BaseModel):
    candidate_names: list[str] = Field(
        ..., min_length=1, max_length=8
    )
    formula: str = Field(default="", max_length=64)
    max_results: int = Field(default=5, ge=1, le=10)

    @field_validator("candidate_names")
    @classmethod
    def _clean_candidate_names(cls, values):
        cleaned = [
            " ".join(str(value or "").split())[:160]
            for value in values
            if str(value or "").strip()
        ]
        if not cleaned:
            raise ValueError(
                "At least one non-empty candidate name is required"
            )
        return cleaned


class ReferenceStandardQCRequest(BaseModel):
    set_id: str = Field(..., min_length=1, max_length=16)
    mass_tolerance_ppm: float = Field(10.0, gt=0.0, le=100.0)


def _analysis_params_with_task_provenance(
    task: dict[str, Any],
    params: AnalysisParams,
    expected_dataset_selection: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Resolve an isolation-width fallback without disguising its source."""
    payload = params.model_dump()
    if expected_dataset_selection is not None:
        selected = list(
            expected_dataset_selection.get("selected_dataset_sources") or []
        )
        payload["expected_compound_datasets"] = selected
        payload["expected_compound_dataset_selection_source"] = str(
            expected_dataset_selection.get("selection_source") or "none"
        )
        payload["expected_compound_dataset_selection_status"] = str(
            expected_dataset_selection.get("selection_status") or ""
        )
        payload["expected_compound_inferred_formula_code"] = str(
            expected_dataset_selection.get("inferred_formula_code") or ""
        )
    blank_tokens = list(payload.get("blank_files") or [])
    legacy_blank = payload.get("blank_file")
    if legacy_blank and legacy_blank not in blank_tokens:
        blank_tokens.append(legacy_blank)
    if blank_tokens:
        blank_tokens, duplicate_blanks = _deduplicate_managed_blank_tokens(
            blank_tokens
        )
        payload["blank_files"] = blank_tokens
        if legacy_blank and Path(str(legacy_blank)).name not in blank_tokens:
            payload["blank_file"] = None
        payload["duplicate_blank_tokens_ignored"] = duplicate_blanks
        supplied_roles = dict(payload.get("blank_roles") or {})
        with _blank_uploads_file_lock:
            stored_blanks = _read_blank_upload_index_unlocked()
        payload["blank_roles"] = {
            Path(str(token)).name: _normalize_blank_role(
                supplied_roles.get(Path(str(token)).name)
                or (stored_blanks.get(Path(str(token)).name) or {}).get(
                    "blank_role"
                )
            )
            for token in blank_tokens
        }
    requested_width = payload.get("isolation_window_width_da")
    if requested_width is not None:
        payload["isolation_window_width_source"] = "user_full_width_override"
        return payload

    method_width = _safe_float(task.get("recommended_isolation_window_width_da"))
    method_source = str(
        task.get("recommended_isolation_window_width_source") or ""
    ).strip()
    if method_width is not None and 0 < method_width <= 20:
        payload["isolation_window_width_da"] = method_width
        payload["isolation_window_width_source"] = (
            f"raw_method:{method_source}"
            if method_source else "raw_method:Agilent AcqMethod.xml"
        )
    return payload


# ── Helper: load task results ────────────────────────────────────────
def _prepare_results_for_display(results: list[dict]) -> list[dict]:
    """Upgrade cached results to the current naming, role and honesty rules."""
    results = _normalize_result_compound_names(results)
    # Exact-InChIKey role curation is also applied to cached legacy tasks, so
    # AO1680-like contaminants are corrected without rerunning the raw sample.
    results = _annotate_result_compound_roles(results)
    for row in results:
        match = row.get("match") if isinstance(row, dict) else None
        if not isinstance(match, dict):
            continue
        # Public libraries mix neutral molecular formulae with explicitly
        # charged ion compositions (for example C21H19O12- for [M-H]-).
        # Repair the cached presentation view before identity deduplication and
        # isomer calibration; new searches already emit this neutral formula.
        formula_candidates = [match, *(match.get("top_candidates") or [])]
        for candidate in formula_candidates:
            if not isinstance(candidate, dict):
                continue
            resolved_formula = _resolve_structure_formula(
                str(candidate.get("smiles") or ""),
                str(candidate.get("formula") or ""),
                str(candidate.get("adduct") or ""),
            )
            if resolved_formula:
                candidate["formula"] = resolved_formula
        smiles = str(match.get("smiles") or "").strip()
        if match.get("has_structure") is not False and smiles.lower() not in {
            "", "nan", "none", "null", "n/a"
        }:
            has_structure, structure_quality = _validate_structure_metadata(
                smiles,
                str(match.get("formula") or ""),
                str(match.get("adduct") or ""),
            )
            match["has_structure"] = has_structure
            match["structure_quality"] = structure_quality
        elif match.get("has_structure") is False or not smiles:
            match["has_structure"] = False
            match.setdefault("structure_quality", "missing")

        # Confidence calibration needs the candidate-independent precursor
        # evidence stored on the scan row.  Historical artifacts kept these
        # fields outside ``match``; copy them into the transient presentation
        # view so the same High gate applies after restart without rewriting
        # the immutable result artifact.
        match["query_search_precursor_mz_source"] = str(
            row.get("search_precursor_mz_source")
            or match.get("query_search_precursor_mz_source")
            or ""
        )
        match["query_adduct"] = (
            row.get("query_adduct") or match.get("query_adduct")
        )
        match["query_adduct_reliable"] = bool(
            row.get("query_adduct_reliable")
            or match.get("query_adduct_reliable")
        )
        match["isotope_evidence_status"] = str(
            row.get("isotope_evidence_status")
            or match.get("isotope_evidence_status")
            or "not_available"
        )

        # Historical caches could expose a near-tie isomer at Top-K[0] even
        # though ``match`` was the production-selected candidate.  Rebuild the
        # public list before honesty calibration so UI, exports and ambiguity
        # guards all share one Top-1.
        raw_top_candidates = match.get("top_candidates")
        if isinstance(raw_top_candidates, list) and raw_top_candidates:
            selected_candidate = dict(match)
            selected_candidate.pop("top_candidates", None)
            # EviGate passport is result-context metadata, not a candidate
            # property.  Do not copy it into Top-K on repeated presentation
            # upgrades (which would duplicate a full evidence DAG per row).
            selected_candidate.pop("evidence_governance", None)
            selected_name = " ".join(
                str(selected_candidate.get("compound_name") or "")
                .casefold().split()
            )
            selected_formula = str(
                selected_candidate.get("formula") or ""
            ).strip().casefold()
            selected_inchikey = _normalize_inchikey(
                selected_candidate.get("inchikey")
            )

            def same_saved_candidate(candidate: Any) -> bool:
                if not isinstance(candidate, dict):
                    return False
                candidate_inchikey = _normalize_inchikey(
                    candidate.get("inchikey")
                )
                if selected_inchikey and candidate_inchikey:
                    return selected_inchikey == candidate_inchikey
                candidate_name = " ".join(
                    str(candidate.get("compound_name") or "")
                    .casefold().split()
                )
                candidate_formula = str(
                    candidate.get("formula") or ""
                ).strip().casefold()
                return bool(
                    selected_name
                    and selected_formula
                    and candidate_name == selected_name
                    and candidate_formula == selected_formula
                )

            # Legacy Top-K occasionally stored the selected entry again as a
            # name/formula-only row.  If that duplicate is kept, it consumes a
            # Top-K slot and the last real isomer competitor disappears during
            # every export.  Remove only the selected duplicate, then retain
            # every genuinely distinct saved candidate even when the selected
            # identity was absent from the old Top-K artifact.
            retained_candidates = [
                candidate for candidate in raw_top_candidates
                if not same_saved_candidate(candidate)
            ]
            match["top_candidates"] = _top_candidates_with_best(
                selected_candidate,
                [selected_candidate, *retained_candidates],
                limit=max(1, 1 + len(retained_candidates)),
            )

        # Legacy JSON stored a confidence label before source provenance and
        # contaminant roles existed.  Recompute all presentation-only honesty
        # fields so an old AO1680 result cannot remain a natural-product High.
        honest = _honest_annotate(match, match.get("top_candidates"))
        match["confidence"] = honest["compound_confidence"]
        match["confidence_class"] = honest["class_confidence"]
        for key in (
            "exact_structure_confidence", "isomer_group_confidence",
            "identification_scope", "isomer_group_formula",
            "isomer_group_size", "isomer_group_evidence",
            "honesty_tier", "isomer_ambiguous", "isomer_n", "isomer_set",
            "same_formula_alternatives", "same_formula_alternatives_n",
            "noncompetitive_same_formula_n",
            "isomer_ambiguity_type", "isomer_relationship_counts",
            "isomer_variant_subtype_counts",
            "isomer_resolution_status",
            "isomer_literature_evidence",
            "class_consensus", "class_agreement", "honest_headline",
            "confidence_semantics", "confidence_is_probability",
            "confidence_is_project_qvalue",
            "legacy_unvalidated_expert_influence_detected",
            "margin_upgrade_counterfactual_eligible",
            "margin_upgrade_runtime_influence",
            "evidence_tier", "annotation_level", "rt_confirmed", "rt_warning",
            "candidate_comparison_source", "candidate_comparison_cosine",
            "candidate_margin", "deconvolution_score_influence",
            "deconvolution_evidence_state", "chimeric_isolation_warning",
            "two_dimensional_connectivity_confidence",
            "isopair_reranked_2d_connectivity",
            "isopair_high_precision_gate_required",
            "isopair_high_precision_support",
            "isopair_high_precision_gate_result", "isopair_high_eligible",
            "isopair_evidence", "precursor_high_evidence_gate_required",
            "precursor_high_evidence_gate_result",
            "precursor_high_evidence_gate_reason",
            "precursor_high_evidence_gate_reason_zh",
            "direct_matched_product_ions",
        ):
            match[key] = honest.get(key)
        details = match.get("expert_details")
        if not isinstance(details, dict):
            continue
        details = _finalize_expert_evidence(
            details,
            fused_score=match.get("expert_score"),
            classification_context={
                "classification": match.get("classification"),
                "classification_source": match.get("classification_source"),
                "classification_confidence": match.get(
                    "classification_confidence"
                ),
            },
        )
        match["expert_details"] = details
        match["expert_score"] = float(
            details.get("fused_expert_score", 0.0) or 0.0
        )
        match["expert_score_runtime_recomputed"] = True
        fusion = details.get("score_fusion") or {}
        match["expert_score_source"] = fusion.get("selected_source", "none")
        match["expert_score_components"] = fusion.get("components", {})
        analysis = details.get("neutral_loss_analysis")
        if isinstance(analysis, dict):
            classification = (
                match.get("display_class") or match.get("classification")
                or details.get("analysis_class") or "Unknown"
            )
            _enrich_neutral_loss_analysis(analysis, classification)
    # Honesty calibration rebuilds compact isomer/same-formula summaries from
    # Top-K candidates.  Run the now-idempotent exact-name and role passes once
    # more so those newly created nested rows inherit curated names and
    # contaminant roles instead of exposing stale unresolved placeholders.
    results = _normalize_result_compound_names(results)
    results = _annotate_result_compound_roles(results)
    results = _apply_sample_consensus(results)
    _annotate_precursor_features(results)
    # EviGate-MS is a read-only claim-governance layer.  It consumes the final
    # honesty calibration and precursor-event annotations, then records which
    # evidence was applicable/independent.  It never mutates Top-1 or scores.
    _annotate_evidence_passports(results)
    return results


PRESENTATION_CALIBRATION_VERSION = (
    "2026-08-22-v12-complete-isomer-candidate-export-v6"
)


def _project_claims_for_task(task_id: str, results: list[dict]) -> list[dict]:
    """Prepare legacy rows and attach the read-only EviGate 2 decision view."""
    prepared = _prepare_results_for_display(results)
    try:
        _attach_claim_decisions(
            prepared,
            _claim_contract_path(TASK_DIR, task_id),
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        logger.warning("Could not project EviGate 2 claims for %s: %s", task_id, exc)
    return prepared


def _task_feature_inventory(
    task: Mapping[str, Any] | None,
) -> tuple[list[dict[str, Any]] | None, Path | None]:
    """Load the immutable OpenMS feature table for result completeness checks.

    Absence means an independent continuous-MS1 inventory was never produced;
    callers must report partition status as unknown rather than deriving a
    tautological pass from the rows that survived matching.
    """

    if not isinstance(task, Mapping):
        return None, None
    workflow = task.get("untargeted_workflow")
    if not isinstance(workflow, Mapping):
        return None, None
    manifest_value = str(workflow.get("manifest") or "").strip()
    if not manifest_value:
        return None, None
    manifest = Path(manifest_value).expanduser().resolve()
    task_root = TASK_DIR.resolve()
    if not manifest.is_relative_to(task_root):
        raise ValueError("Feature inventory manifest is outside the task store")
    if not manifest.is_file():
        raise FileNotFoundError(f"Feature inventory manifest is missing: {manifest}")
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    features = payload.get("features")
    if not isinstance(features, list):
        raise ValueError("Feature inventory manifest has no features array")
    declared_count = payload.get("feature_count")
    if declared_count is not None and int(declared_count) != len(features):
        raise ValueError(
            "Feature inventory count differs from the manifest features array"
        )
    ion_mode = str(task.get("ion_mode") or "").strip().upper()
    normalized = []
    for feature in features:
        if not isinstance(feature, Mapping):
            raise ValueError("Feature inventory contains a non-object row")
        row = dict(feature)
        if ion_mode and not row.get("ion_mode"):
            row["ion_mode"] = ion_mode
        normalized.append(row)
    return normalized, manifest


def _get_task_results(task_id: str) -> list[dict]:
    """Load task results from memory or disk."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    results = task.get("results", [])
    if not results:
        result_path = TASK_DIR / f"{task_id}_results.json"
        if result_path.exists():
            try:
                with open(result_path, "r", encoding="utf-8") as f:
                    results = json.load(f)
            except Exception as exc:
                raise HTTPException(
                    status_code=500,
                    detail=(
                        "Saved task results are unreadable; the task must be rerun "
                        f"({result_path.name}: {exc})"
                    ),
                ) from exc
    return _project_claims_for_task(task_id, results)


_RESULT_SUMMARY_ROW_FIELDS = {
    "spectrum_id", "precursor_mz", "retention_time", "tic", "ion_mode",
    "search_precursor_mz", "search_precursor_mz_source",
    "search_precursor_delta_ppm", "lcms_feature_used_for_search",
    "isotope_evidence_status",
    "polarity", "precursor_charge", "collision_energy",
    "collision_energy_unit", "collision_energy_mode", "activation_method",
    "instrument_domain_id", "collision_cell_region", "spectrum_representation",
    "raw_peaks_complete", "precursor_feature_id", "lcms_feature_id",
    "limited_fragment_list", "peaks_count", "processed_peaks_count",
    "display_class",
    "dl_confidence", "dl_top3", "input_evidence_note",
    "blank_excluded", "blank_flagged", "blank_clean_view_excluded",
    "blank_feature_matched", "blank_sample_ratio", "blank_status",
    "blank_evidence_category", "blank_evidence_role",
    "ion_identity_network_id", "ion_identity_linked_annotation",
    "query_adduct", "query_adduct_reliable", "lcms_adduct",
    "lcms_neutral_mass",
}

_RESULT_SUMMARY_MATCH_FIELDS = {
    "compound_name", "compound_name_cn", "source_compound_name", "formula",
    "smiles", "has_structure", "structure_quality", "inchikey", "adduct",
    "dataset_source", "source_tier", "lib_spectrum_id", "cosine",
    "raw_cosine", "final_score", "raw_final_score", "rank_score",
    "reverse_match", "mass_ppm", "fragment_score", "mass_score",
    "intensity_score", "expert_score", "expert_score_source",
    "classification", "display_class", "classification_source",
    "display_classification_source", "classification_reliable",
    "classification_status", "classification_reliability_reason",
    "classification_confidence", "display_classification_confidence",
    "confidence", "confidence_class",
    "exact_structure_confidence", "isomer_group_confidence",
    "identification_scope", "isomer_group_formula", "isomer_group_size",
    "honesty_tier", "isomer_ambiguous", "isomer_n",
    "same_formula_alternatives_n", "isomer_ambiguity_type",
    "isomer_variant_subtype_counts", "isomer_resolution_status",
    "honest_headline", "evidence_tier", "annotation_level",
    "candidate_margin", "name_only_lead", "potential_contaminant",
    "contaminant_note", "contaminant_category", "compound_role",
    "compound_role_note", "exclude_from_natural_product_count",
    "natural_product_candidate", "screening_domain", "pesticide_use",
    "chemical_family", "regulatory_status", "identity_source",
    "screening_catalog_id", "screening_catalog_released_at",
    "screening_identity_scope", "screening_identity_note",
    "role_catalog_id", "role_catalog_version", "role_annotation_scope",
    "role_annotation_source", "chimeric_isolation_warning",
    "feature_fragment_evidence",
    "ion_identity_adduct_conflict", "ion_identity_conflict",
    "isopair_high_precision_gate_result",
    "precursor_high_evidence_gate_required",
    "precursor_high_evidence_gate_result",
    "precursor_high_evidence_gate_reason",
    "precursor_high_evidence_gate_reason_zh",
    "direct_matched_product_ions",
    "precursor_metadata_audit", "mass_metadata_consistent",
    "candidate_adduct_hypothesis", "candidate_adduct_hypothesis_source",
    "query_selected_ion_mz", "query_search_precursor_mz",
    "query_search_precursor_mz_source",
    "query_search_precursor_delta_ppm", "query_precursor_mass_audit",
}


def _summary_evidence_flags(
    match: dict[str, Any], row: Optional[dict[str, Any]] = None
) -> dict[str, Any]:
    """Collapse expensive expert/pathway objects to a compact list-view summary.

    The result page must be able to show that classification, biosynthetic
    reasoning and fragmentation analysis actually ran without shipping the
    full expert payload for every table row.  Keep the strings deliberately
    short and capped; the complete reasoning remains available from the item
    endpoint used by the candidate detail dialog.
    """

    from web.report_common import fragmentation_evidence_status

    details = match.get("expert_details") or {}
    row = row or {}
    fragmentation_status = fragmentation_evidence_status(match)
    neutral_loss = details.get("neutral_loss_analysis") or {}
    v3 = details.get("v3") or {}
    pathway_status = str(match.get("pathway_polarity_status") or "")
    has_library_pathway = bool(match.get("pathway_json"))
    reasoning = v3.get("reasoning") or []
    if isinstance(reasoning, str):
        reasoning = [reasoning]
    biosynthetic_pathways: list[str] = []
    for value in reasoning:
        text = " ".join(str(value or "").split())
        if not text or ("生源" not in text and "biosynth" not in text.lower()):
            continue
        text = text[:180]
        if text not in biosynthetic_pathways:
            biosynthetic_pathways.append(text)
        if len(biosynthetic_pathways) >= 3:
            break

    best_pathway = neutral_loss.get("best_pathway") or []
    sequential_losses = neutral_loss.get("sequential_losses") or []
    neutral_losses = neutral_loss.get("neutral_losses") or []
    product_structures = (
        (details.get("structure_fragment_analysis") or {}).get(
            "plausible_fragments"
        )
        or []
    )
    reconstructed_step_count = len(best_pathway) + len(sequential_losses)
    product_structure_count = len(product_structures)
    diagnostic_ions = (details.get("diagnostic_ion_analysis") or {}).get(
        "matched_ions"
    ) or []
    def _specific_ion_form(value: Any) -> str:
        text = " ".join(str(value or "").split())
        return "" if text.casefold() in {
            "", "none", "nan", "null", "unknown", "pos", "positive",
            "+", "neg", "negative", "-",
        } else text

    lcms_adduct = _specific_ion_form(row.get("lcms_adduct"))
    query_adduct = _specific_ion_form(
        row.get("query_adduct") or match.get("query_adduct")
    )
    if lcms_adduct:
        sample_ion_form = lcms_adduct
        sample_ion_form_source = "lcms_ion_identity"
    elif query_adduct and bool(
        row.get("query_adduct_reliable")
        or match.get("query_adduct_reliable")
    ):
        sample_ion_form = query_adduct
        sample_ion_form_source = "explicit_query_metadata"
    else:
        sample_ion_form = ""
        sample_ion_form_source = "unresolved"
    candidate_ion_form_hypothesis = _specific_ion_form(
        match.get("candidate_adduct_hypothesis")
    )
    return {
        "has_neutral_loss": bool(
            neutral_losses
            or sequential_losses
        ),
        "has_fragmentation": bool(
            reconstructed_step_count
            or product_structure_count
        ),
        "has_fragmentation_evidence": bool(
            reconstructed_step_count
            or product_structure_count
            or neutral_losses
            or diagnostic_ions
        ),
        "has_library_pathway_reference": has_library_pathway,
        "has_reasoning": bool(reasoning or details.get("chemical_reasoning")),
        "biosynthetic_pathways": biosynthetic_pathways,
        "reasoning_count": len(reasoning),
        "neutral_loss_count": len(neutral_losses),
        "fragmentation_step_count": reconstructed_step_count,
        "product_structure_count": product_structure_count,
        "diagnostic_ion_count": len(diagnostic_ions),
        "sample_ion_form": sample_ion_form,
        "sample_ion_form_source": sample_ion_form_source,
        "library_adduct": _specific_ion_form(match.get("adduct")),
        # This is inferred from this candidate's structure/formula and the
        # library precursor.  It is useful for audit/display, but it is not an
        # observed sample adduct and must never be presented as independent
        # ion-identity evidence.
        "candidate_ion_form_hypothesis": candidate_ion_form_hypothesis,
        "candidate_ion_form_hypothesis_source": str(
            match.get("candidate_adduct_hypothesis_source") or ""
        ),
        "pathway_status": pathway_status or (
            "available" if has_library_pathway else "missing"
        ),
        "fragmentation_evidence_type": fragmentation_status["evidence_type"],
        "candidate_pathway_eligible_for_fusion": fragmentation_status[
            "candidate_pathway_eligible_for_fusion"
        ],
        "candidate_pathway_runtime_influence": fragmentation_status[
            "candidate_pathway_runtime_influence"
        ],
        "fragmentation_mechanism_confirmed": False,
    }


def _result_summary_projection(
    row: dict[str, Any], result_index: int
) -> dict[str, Any]:
    """Return the browser list projection without peak/evidence bulk data."""

    projected = {
        key: value for key, value in row.items()
        if key in _RESULT_SUMMARY_ROW_FIELDS
    }
    projected["result_index"] = int(result_index)
    if not projected.get("peaks_count"):
        projected["peaks_count"] = len(row.get("query_peaks") or [])
    match = row.get("match")
    if isinstance(match, dict):
        match_projection = {
            key: value for key, value in match.items()
            if key in _RESULT_SUMMARY_MATCH_FIELDS
        }
        passport = match.get("evidence_governance")
        if isinstance(passport, dict):
            match_projection["evidence_governance"] = {
                key: passport.get(key) for key in (
                    "policy_id", "claim_status", "claim_scope",
                    "release_eligible_before_replication", "release_blockers",
                    "applicable_evidence_nodes", "independent_support_groups",
                ) if key in passport
            }
        deconvolution = match.get("chimeric_deconvolution")
        if isinstance(deconvolution, dict):
            match_projection["chimeric_deconvolution"] = {
                key: deconvolution.get(key) for key in (
                    "status", "isolation_status", "target_fraction",
                    "score_contribution", "ranking_changed",
                ) if key in deconvolution
            }
        ion_identity = match.get("ion_identity_evidence")
        if isinstance(ion_identity, dict):
            match_projection["ion_identity_evidence"] = {
                key: ion_identity.get(key) for key in (
                    "network_id", "identity_conflict", "ion_identity",
                    "neutral_mass", "adduct_consistency",
                ) if key in ion_identity
            }
        decision = match.get("final_claim_decision")
        if isinstance(decision, dict):
            calibration = decision.get("calibration") or {}
            match_projection["final_claim_decision"] = {
                key: decision.get(key) for key in (
                    "policy_id", "released", "status", "scope",
                    "independent_injection_count", "reason_codes",
                ) if key in decision
            }
            if isinstance(calibration, dict):
                match_projection["final_claim_decision"]["calibration"] = {
                    key: calibration.get(key) for key in (
                        "status", "q_connectivity", "q_isomer_group",
                        "q_formula", "q_class", "reason",
                    ) if key in calibration
                }
        match_projection["display_evidence_flags"] = _summary_evidence_flags(
            match, row
        )
        projected["match"] = match_projection
    else:
        projected["match"] = None
    return projected


def _with_task_polarity_context(
    results: list[dict], task: dict[str, Any]
) -> tuple[list[dict], str]:
    """Backfill polarity for legacy result JSON only when task mode is singular."""
    missing = [row for row in results if not str(row.get("ion_mode") or "").strip()]
    if not missing:
        return results, "per_spectrum_result"

    mode = str(task.get("ion_mode") or "").strip().upper()
    counts = task.get("ion_mode_counts") or {}
    if mode not in {"POS", "NEG"}:
        pos = int(counts.get("POS") or 0)
        neg = int(counts.get("NEG") or 0)
        if pos and not neg:
            mode = "POS"
        elif neg and not pos:
            mode = "NEG"
    if mode not in {"POS", "NEG"}:
        raise ValueError(
            "Legacy results lack per-spectrum polarity and the task is mixed/unknown; "
            "rerun this task with V6 before literature-target screening"
        )

    output: list[dict] = []
    for row in results:
        copied = dict(row)
        if not str(copied.get("ion_mode") or "").strip():
            copied["ion_mode"] = mode
        output.append(copied)
    return output, "task_level_single_polarity_backfill"


# ── Helper Functions ────────────────────────────────────────────────
def _safe_float(v) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _build_mgf_tic_overview(
    rows: list[dict],
    rt_in_seconds: bool = True,
    source: str = "mgf_ms2_approximation",
) -> Optional[dict]:
    """Build an MS2-only TIC overview from parsed spectra or saved results.

    单位换算由 rt_in_seconds 控制: MGF 的 RTINSECONDS 是秒 (viewer 用分钟, 需 ÷60);
    raw_blank_analysis HTML 报告的 retention_time 已是分钟, 传 False 不再换算。无论哪种
    来源这都是 MS2 近似 (无连续 MS1 色谱), 非真实 MS1 TIC。source 标注来源供前端区分。
    """
    ms2 = []
    for i, row in enumerate(rows):
        rt_val = _safe_float(row.get("retention_time"))
        if rt_val is None:
            continue
        rt_out = (rt_val / 60.0) if rt_in_seconds else rt_val
        peaks = row.get("peaks") or row.get("query_peaks") or []
        tic = _safe_float(row.get("tic"))
        if tic is None:
            try:
                tic = sum(float(p[1]) for p in peaks if len(p) >= 2)
            except (TypeError, ValueError):
                tic = 0.0
        if tic <= 0:
            continue

        scan_id = row.get("spectrum_id", i)
        try:
            scan_id = int(scan_id)
        except (TypeError, ValueError):
            scan_id = str(scan_id)
        bpc = _safe_float(row.get("bpc"))
        base_peak_mz = _safe_float(row.get("base_peak_mz"))
        if bpc is None and peaks:
            try:
                base_peak = max(
                    ((float(p[1]), float(p[0])) for p in peaks if len(p) >= 2),
                    default=None,
                )
                if base_peak is not None:
                    bpc, base_peak_mz = base_peak
            except (TypeError, ValueError):
                bpc = None
                base_peak_mz = None
        ms2.append({
            "rt": rt_out,
            "tic": tic,
            "bpc": bpc,
            "base_peak_mz": base_peak_mz,
            "scan_id": scan_id,
            "precursor_mz": _safe_float(row.get("precursor_mz")),
        })

    if not ms2:
        return None
    ms2.sort(key=lambda p: p["rt"])
    return {
        "chromatogram_schema_version": 3,
        "ms1": [],
        "ms2": ms2,
        "rt_unit": "minute",
        "source": source,
        "tic_source": "ms2_spectrum_metadata_approximation",
        "ms1_spectra_available": False,
        "ms1_spectrum_count": 0,
        "ms1_tic_point_count": 0,
        "bpc_available": False,
        "bpc_source": "unavailable",
        "bpc_scope": "unavailable_without_continuous_ms1",
        "bpc_backfill_attempted": True,
    }


def _tic_overview_has_ms1_bpc(payload: dict[str, Any]) -> bool:
    ms1 = payload.get("ms1") or []
    return bool(ms1) and all(
        point.get("bpc") is not None for point in ms1
        if isinstance(point, dict)
    ) and len(ms1) == sum(isinstance(point, dict) for point in ms1)


def _tic_overview_has_ms1_tic(payload: dict[str, Any]) -> bool:
    points = [
        point for point in (payload.get("ms1") or [])
        if isinstance(point, dict)
    ]
    if len(points) < 2:
        return False
    return all(
        _safe_float(point.get("rt")) is not None
        and _safe_float(point.get("tic")) is not None
        for point in points
    )


def _load_tic_overview_with_bpc_backfill(
    task_id: str,
    task: dict[str, Any],
    tic_path: Path,
) -> dict[str, Any]:
    """Upgrade old chromatogram artifacts from the original mzML on first read.

    Schema v3 also recovers MS1-labelled TIC points that some vendor exports keep
    only in ``chromatogramList``.  Those points improve TIC without pretending
    that an absent MS1 peak array or absent BPC chromatogram can produce BPC.
    """
    with open(tic_path, "r", encoding="utf-8") as handle:
        overview = json.load(handle)
    if _tic_overview_has_ms1_bpc(overview):
        overview["bpc_available"] = True
        return overview
    schema_version = int(overview.get("chromatogram_schema_version") or 1)
    if (
        overview.get("bpc_backfill_attempted") is True
        and schema_version >= 3
    ):
        return overview

    candidate_paths = [task.get("file_path")]
    representation = task.get("spectrum_representation_processing") or {}
    candidate_paths.append(representation.get("output"))
    attempted = False
    seen_paths: set[Path] = set()
    for raw_path in candidate_paths:
        if not raw_path:
            continue
        source_path = Path(str(raw_path))
        if not source_path.exists() or source_path.suffix.lower() not in {
            ".mzml", ".mzxml", ".xml",
        }:
            continue
        resolved_path = source_path.resolve()
        if resolved_path in seen_paths:
            continue
        seen_paths.add(resolved_path)
        from utils.mzml_io import extract_tic_overview

        enriched = extract_tic_overview(source_path)
        attempted = True
        if enriched:
            has_bpc = _tic_overview_has_ms1_bpc(enriched)
            has_tic = _tic_overview_has_ms1_tic(enriched)
            if has_bpc or has_tic:
                enriched["bpc_available"] = has_bpc
                enriched["bpc_backfilled"] = has_bpc
                enriched["bpc_backfill_attempted"] = True
                enriched["bpc_backfill_task_id"] = task_id
                if has_tic and not enriched.get("ms1_spectra_available", True):
                    enriched["tic_chromatogram_backfilled"] = True
                _write_json_atomic(tic_path, enriched)
                logger.info(
                    "[task %s] chromatogram overview upgraded: %d MS1 TIC "
                    "points, BPC available=%s",
                    task_id,
                    len(enriched.get("ms1") or []),
                    has_bpc,
                )
                return enriched

    if attempted:
        overview["chromatogram_schema_version"] = max(
            int(overview.get("chromatogram_schema_version") or 1), 3
        )
        overview["bpc_available"] = False
        overview["bpc_backfill_attempted"] = True
        overview["bpc_scope"] = "unavailable_in_source_mzml_metadata"
        _write_json_atomic(tic_path, overview)
    return overview


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _parse_spectrum_file(
    file_path: str,
    isolation_window_width_da: float | None = None,
    isolation_window_width_source: str | None = None,
    allow_unknown_representation: bool = False,
) -> list[dict]:
    """Parse MGF, mzML/mzXML, or an embedded-MS2 HTML report.

    mzML profile data must already have passed the per-spectrum vendor/OpenMS
    peak-picking stage. MGF/HTML are discrete peak-list formats and are never
    re-centroided by a spacing heuristic.
    """
    from utils.mgf_io import parse_mgf_file as _parse_mgf
    from pipeline.mzml_input import parse_mzml_file as _parse_mzml
    from utils.report_html_io import (
        parse_raw_report_html as _parse_html,
        is_raw_report_html as _is_report_html,
    )
    ext = Path(file_path).suffix.lower()
    if ext in (".mgf",):
        spectra = _parse_mgf(file_path)
        for spec in spectra:
            spec["spectrum_representation"] = "centroid"
            spec["spectrum_representation_source"] = "peak_list_format_contract"
    elif ext in (".mzml", ".mzxml", ".xml"):
        # parse_mzml_file 已按 mzML cvParam (MS:1000127 centroid) 决定是否质心化 (Fix1),
        # 不再二次 smart_peak_pick — 否则会对已质心谱图误判 profile 并错误合并, 污染 cosine。
        spectra = _parse_mzml(
            file_path,
            isolation_window_width_da=isolation_window_width_da,
            isolation_window_width_source=isolation_window_width_source,
            allow_unknown_representation=allow_unknown_representation,
        )
    elif ext in (".html", ".htm"):
        # 外部前处理工具 (raw_blank_analysis_v2) 的 HTML 峰级报告: 每个 UV 峰一张 MS2 谱图
        # (前体+碎片+%+离子模式) 直接嵌在 HTML 里。原始 .d 数据常留在采集仪器上拿不到,
        # 报告本身就是可用输入。先探测签名, 避免把 mssystem 自身导出的结果 HTML 当样本解析。
        if not _is_report_html(file_path):
            raise ValueError(
                "该 .html 不是可解析的分析报告 (缺 MS2 碎片/前体签名)。"
                "仅支持 raw_blank_analysis 风格的峰值报告。"
            )
        # 碎片已是离散质心 top-N 列表 (非 profile), 不做 smart_peak_pick ——
        # 否则 _is_profile_like 启发式会误判并错误合并, 污染 cosine (同 mzML centroid 理由)。
        spectra = _parse_html(file_path)
        for spec in spectra:
            spec["spectrum_representation"] = "centroid"
            spec["spectrum_representation_source"] = "peak_list_format_contract"
    else:
        raise ValueError(f"Unsupported file format: {ext}")

    return spectra


def _resolve_spectrum_polarities(
    spectra: list[dict], requested_mode: str = "auto"
) -> tuple[str, dict[str, int]]:
    """Resolve and audit ion polarity before matching.

    Explicit POS/NEG applies to every spectrum.  Auto mode accepts per-scan
    metadata, fills isolated unknown scans only when the whole file has one
    unambiguous known polarity, and refuses ambiguous/unknown input.  Refusing
    is intentional: query_polarity=None disables the library polarity filter.
    """
    mode = str(requested_mode or "auto").strip().upper()
    if mode == "AUTO":
        mode = "auto"
    if mode not in {"auto", "POS", "NEG"}:
        raise ValueError("ion_mode must be auto, POS, or NEG")

    if mode in {"POS", "NEG"}:
        for spec in spectra:
            spec["polarity"] = mode
            spec["polarity_source"] = "user_override"
        return mode, {mode: len(spectra), "unknown": 0}

    conflicts = sum(bool(spec.get("polarity_conflict")) for spec in spectra)
    if conflicts:
        raise ValueError(
            f"Detected {conflicts} spectra with conflicting IONMODE/CHARGE metadata; "
            "choose positive or negative ion mode explicitly."
        )

    known = [spec.get("polarity") for spec in spectra
             if spec.get("polarity") in {"POS", "NEG"}]
    known_modes = set(known)
    unknown = [spec for spec in spectra if spec.get("polarity") not in {"POS", "NEG"}]
    if not known_modes:
        raise ValueError(
            "Ion mode could not be identified from this file. Choose positive or "
            "negative ion mode explicitly before running."
        )
    if len(known_modes) > 1 and unknown:
        raise ValueError(
            "This is a mixed-polarity file, but some scans have no ion-mode metadata. "
            "Split the file by polarity or choose one ion mode explicitly."
        )

    if len(known_modes) == 1:
        consensus = next(iter(known_modes))
        for spec in unknown:
            spec["polarity"] = consensus
            spec["polarity_source"] = "file_consensus"

    pos = sum(spec.get("polarity") == "POS" for spec in spectra)
    neg = sum(spec.get("polarity") == "NEG" for spec in spectra)
    unresolved = len(spectra) - pos - neg
    detected = "MIXED" if pos and neg else ("POS" if pos else "NEG")
    return detected, {"POS": pos, "NEG": neg, "unknown": unresolved}


def _start_task_stage_heartbeat(
    task_id: str,
    *,
    stage: str,
    message: str,
    timeout_seconds: int,
    interval_seconds: float = 15.0,
) -> tuple[threading.Event, threading.Thread]:
    """Keep a long indivisible stage visibly alive without inventing progress."""
    stop_event = threading.Event()
    started = time.monotonic()

    def _heartbeat() -> None:
        while not stop_event.wait(max(float(interval_seconds), 0.01)):
            elapsed = max(0, int(time.monotonic() - started))
            with _tasks_lock:
                task = _tasks.get(task_id)
                if not task or task.get("status") != "running":
                    return
                task["message"] = (
                    f"{message} Still active ({elapsed}s elapsed; "
                    f"bounded at {int(timeout_seconds)}s)."
                )
                task["stage_heartbeat"] = {
                    "stage": stage,
                    "status": "active",
                    "elapsed_seconds": elapsed,
                    "timeout_seconds": int(timeout_seconds),
                    "updated_at": _now_iso(),
                }
                heartbeat_message = task["message"]
                heartbeat_payload = dict(task["stage_heartbeat"])
            _save_task_progress_throttled(task_id)
            _sse_broadcast(
                "progress",
                {
                    "task_id": task_id,
                    "status": "running",
                    "progress": 3,
                    "message": heartbeat_message,
                    "stage_heartbeat": heartbeat_payload,
                },
            )

    thread = threading.Thread(
        target=_heartbeat,
        name=f"task-heartbeat-{task_id[:8]}",
        daemon=True,
    )
    thread.start()
    return stop_event, thread


def _run_task_isomer_literature(
    task_id: str,
    results: Optional[list[dict]] = None,
) -> None:
    """Build a separate PubMed artifact; analysis success is never coupled to it."""
    artifact_path = TASK_DIR / f"{task_id}_isomer_literature.json"
    try:
        with _tasks_lock:
            task = _tasks.get(task_id)
            if not task:
                return
            task["isomer_literature_summary"] = {
                "status": "running",
                "groups_processed": 0,
                "used_for_scoring": False,
                "updated_at": _now_iso(),
            }
        _save_tasks()
        result_path = TASK_DIR / f"{task_id}_results.json"
        if results is None:
            with result_path.open("r", encoding="utf-8") as handle:
                results = json.load(handle)
        if not isinstance(results, list):
            raise ValueError("Saved task results are not a JSON array")

        def progress(index, total, group, lookup):
            should_save = False
            with _tasks_lock:
                current = _tasks.get(task_id)
                if not current:
                    raise RuntimeError(
                        "task was removed while literature enrichment was running"
                    )
            if index == 1 or index == total or index % 5 == 0:
                with _tasks_lock:
                    current = _tasks.get(task_id)
                    if not current:
                        raise RuntimeError(
                            "task was removed while literature enrichment was running"
                        )
                    current["isomer_literature_summary"] = {
                        "status": "running",
                        "groups_processed": index,
                        "groups_total": total,
                        "current_group": group.get("group_key", ""),
                        "last_lookup_status": lookup.get(
                            "status", "unknown"
                        ),
                        "used_for_scoring": False,
                        "updated_at": _now_iso(),
                    }
                    should_save = True
            if should_save:
                _save_task_progress_throttled(task_id)

        artifact = enrich_isomer_groups_with_pubmed(
            results,
            cache_dir=OUTPUT_DIR / "isomer_literature_cache",
            max_groups=0,
            max_results=5,
            timeout=12.0,
            progress_callback=progress,
        )
        artifact["task_id"] = task_id
        if result_path.is_file():
            digest = hashlib.sha256()
            with result_path.open("rb") as handle:
                for chunk in iter(
                    lambda: handle.read(1024 * 1024), b""
                ):
                    digest.update(chunk)
            artifact["source_results_file"] = str(result_path.resolve())
            artifact["source_results_size_bytes"] = (
                result_path.stat().st_size
            )
            artifact["source_results_sha256"] = digest.hexdigest()
        summary = {
            "status": "completed",
            "groups_processed": artifact["groups_processed"],
            "groups_total": artifact["total_unique_ambiguous_groups"],
            "status_counts": artifact["status_counts"],
            "artifact": artifact_path.name,
            "used_for_scoring": False,
            "updated_at": _now_iso(),
        }
        with _tasks_lock:
            current = _tasks.get(task_id)
            if not current:
                return
            # Hold the task-registry lock through the final atomic write.  A
            # deletion therefore either archives this artifact afterwards or
            # wins first and prevents the artifact from being recreated.
            _write_json_atomic(artifact_path, artifact)
            current["isomer_literature_summary"] = summary
        _save_tasks()
    except Exception as exc:
        logger.warning(
            "Isomer literature enrichment failed for task %s: %s",
            task_id,
            exc,
        )
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current:
                current["isomer_literature_summary"] = {
                    "status": "failed",
                    "error": f"{type(exc).__name__}: {exc}"[:300],
                    "used_for_scoring": False,
                    "updated_at": _now_iso(),
                }
        _save_tasks()


def _maybe_launch_isomer_literature_enrichment(
    task_id: str,
    results: Optional[list[dict]] = None,
) -> bool:
    artifact_path = TASK_DIR / f"{task_id}_isomer_literature.json"
    with _tasks_lock:
        task = _tasks.get(task_id)
        if not task:
            return False
        status = str(
            (task.get("isomer_literature_summary") or {}).get("status")
            or ""
        )
        if _is_active_task_state(status):
            return False
        if status == "completed" and artifact_path.exists():
            return False
        task["isomer_literature_summary"] = {
            "status": "queued",
            "groups_processed": 0,
            "used_for_scoring": False,
            "updated_at": _now_iso(),
        }
    _save_tasks()
    _executor_literature.submit(
        _run_task_isomer_literature,
        task_id,
        results,
    )
    return True


def _evograph_task_summary(metadata: dict[str, Any], *, status: str) -> dict[str, Any]:
    graph_totals = metadata.get("graph_totals") or {}
    edge_counts = metadata.get("edge_status_counts") or {}
    formula_counts = metadata.get("formula_context_counts") or {}
    return {
        "status": status,
        "module_state": "shadow",
        "policy_id": metadata.get("policy_id"),
        "spectra_processed": int(metadata.get("record_count") or 0),
        "truncated_by_task_limit": bool(metadata.get("truncated_by_task_limit")),
        "maximum_spectra_per_task": int(
            metadata.get("maximum_spectra_per_task") or 0
        ),
        "graphs_completed": int(
            (metadata.get("status_counts") or {}).get("completed") or 0
        ),
        "formula_context_available": int(formula_counts.get("available") or 0),
        "nodes": int(graph_totals.get("node_count") or 0),
        "selected_edges": int(graph_totals.get("selected_edge_count") or 0),
        "heuristically_supported_edges": int(
            edge_counts.get("heuristically_supported") or 0
        ),
        "ambiguous_edges": int(edge_counts.get("ambiguous") or 0),
        "unassigned_edges": int(edge_counts.get("unassigned") or 0),
        "cycle_count": int(graph_totals.get("cycle_count") or 0),
        "safety_invariants_passed": bool(metadata.get("safety_invariants_passed")),
        "external_calibration_status": str(
            metadata.get("external_calibration_status") or "unavailable"
        ),
        "heuristic_marginals_are_probabilities": False,
        "log_bayes_factor_available": False,
        "production_ranking_changed_records": 0,
        "legacy_confidence_changed_records": 0,
        "high_promotions": 0,
        "used_for_scoring": False,
        "artifact": metadata.get("artifact"),
        "updated_at": metadata.get("generated_at") or _now_iso(),
    }


def _run_task_evograph_fragmentation(
    task_id: str,
    expected_run_id: str | None,
    force: bool = False,
) -> None:
    """Materialize EvoGraph-F after commit, never inside the search executor."""

    result_path = TASK_DIR / f"{task_id}_results.json"
    artifact_path = _evograph_sidecar_path(TASK_DIR, task_id)
    with _tasks_lock:
        current = _tasks.get(task_id)
        if (
            not current
            or current.get("status") != TaskState.COMPLETED.value
            or current.get("result_run_id") != expected_run_id
        ):
            return
        current["evograph_fragmentation_summary"] = {
            "status": "running",
            "module_state": "shadow",
            "spectra_processed": 0,
            "used_for_scoring": False,
            "production_ranking_changed_records": 0,
            "updated_at": _now_iso(),
        }
    _save_tasks()

    try:
        metadata = _materialize_evograph_sidecar(
            task_id,
            result_path,
            artifact_path,
            policy_path=PROJECT_ROOT / "config" / "evograph_policy.json",
            runtime_snapshot_id=_current_runtime_snapshot_id(),
            force=force,
        )
        with _tasks_lock:
            current = _tasks.get(task_id)
            still_current = bool(
                current
                and current.get("status") == TaskState.COMPLETED.value
                and current.get("result_run_id") == expected_run_id
            )
            if still_current:
                current["evograph_fragmentation_summary"] = _evograph_task_summary(
                    metadata,
                    status="completed",
                )
        if not still_current:
            artifact_path.unlink(missing_ok=True)
        _save_tasks()
    except Exception as exc:
        logger.warning(
            "EvoGraph-F shadow materialization failed for task %s: %s",
            task_id,
            exc,
        )
        with _tasks_lock:
            current = _tasks.get(task_id)
            if (
                current
                and current.get("status") == TaskState.COMPLETED.value
                and current.get("result_run_id") == expected_run_id
            ):
                current["evograph_fragmentation_summary"] = {
                    "status": "failed",
                    "module_state": "shadow",
                    "error": f"{type(exc).__name__}: {exc}"[:500],
                    "used_for_scoring": False,
                    "production_ranking_changed_records": 0,
                    "legacy_confidence_changed_records": 0,
                    "high_promotions": 0,
                    "updated_at": _now_iso(),
                }
        _save_tasks()


def _maybe_launch_evograph_fragmentation(task_id: str, *, force: bool = False) -> bool:
    artifact_path = _evograph_sidecar_path(TASK_DIR, task_id)
    with _tasks_lock:
        task = _tasks.get(task_id)
        if not task or task.get("status") != TaskState.COMPLETED.value:
            return False
        summary_status = str(
            (task.get("evograph_fragmentation_summary") or {}).get("status") or ""
        )
        if summary_status in {"queued", "running"}:
            return False
        if not force and summary_status == "completed" and artifact_path.is_file():
            return False
        expected_run_id = task.get("result_run_id")
        task["evograph_fragmentation_summary"] = {
            "status": "queued",
            "module_state": "shadow",
            "spectra_processed": 0,
            "used_for_scoring": False,
            "production_ranking_changed_records": 0,
            "legacy_confidence_changed_records": 0,
            "high_promotions": 0,
            "updated_at": _now_iso(),
        }
    _save_tasks()
    _executor_evograph.submit(
        _run_task_evograph_fragmentation,
        task_id,
        expected_run_id,
        force,
    )
    return True


def _run_analysis(task_id: str):
    """Execute analysis in background thread with robust progress tracking."""
    global _tasks, _lib_records, _mz_index
    cancelled_before_start = False
    with _tasks_lock:
        task = _tasks.get(task_id)
        if not task:
            return
        cancel_event = _task_cancel_events.setdefault(task_id, threading.Event())
        if cancel_event.is_set() or task.get("status") == TaskState.CANCELLING.value:
            _transition_task(
                task,
                TaskState.CANCELLED,
                reason="cancelled_before_worker_started",
            )
            task["message"] = "Analysis cancelled before execution"
            task["error"] = None
            cancelled_before_start = True
        else:
            _transition_task(
                task,
                TaskState.RUNNING,
                reason="executor_worker_started",
            )
            task["progress"] = 0
            task["message"] = "Starting analysis..."
        task_filename = task.get("filename", "")
    _save_tasks()
    if cancelled_before_start:
        _sse_broadcast(
            "status_change",
            {
                "task_id": task_id,
                "status": "cancelled",
                "progress": 0,
                "message": "Analysis cancelled before execution",
            },
        )
        with _tasks_lock:
            _task_cancel_events.pop(task_id, None)
        return
    _sse_broadcast("status_change", {"task_id": task_id, "status": "running", "progress": 0})
    logger.info(
        "Analysis execution started for task %s (%s)",
        task_id,
        task_filename,
    )

    t_start = time.time()
    try:
        raw_file_path = str(task["file_path"])
        file_path = raw_file_path
        params = task.get("params", {})
        # A missing formal-ranking asset must fail before parsing/searching a
        # complete sample. Full content/binding verification remains owned by
        # the existing verified retriever; this is only its existence preflight.
        from runtime.verified_sql_retrieval import sql_metadata_dependency_paths
        missing_primary_dependencies = [str(path) for path in
            sql_metadata_dependency_paths(PROJECT_ROOT)
            if not str(path).endswith(("-wal", "-journal")) and not path.is_file()]
        with _tasks_lock:
            _tasks[task_id]["primary_dependency_preflight"] = {
                "status": "failed" if missing_primary_dependencies else "passed",
                "missing_paths": missing_primary_dependencies,
                "content_verification": "performed_by_verified_sql_retriever",
            }
        if missing_primary_dependencies:
            raise FileNotFoundError("Required formal-ranking dependencies are missing: "
                + ", ".join(missing_primary_dependencies))
        artifact_dir = TASK_DIR / f"{task_id}_artifacts"
        artifact_dir.mkdir(parents=True, exist_ok=True)
        representation_processing = {
            "status": "not_applicable_peak_list_format",
            "input": raw_file_path,
            "output": raw_file_path,
        }
        input_extension = Path(raw_file_path).suffix.lower()
        if input_extension in {".mzml", ".mzxml", ".xml"}:
            from utils.openms_backend import prepare_mzml_for_analysis

            prepared_path, representation_processing = prepare_mzml_for_analysis(
                raw_file_path,
                artifact_dir,
                policy=params.get("profile_peak_picking", "auto"),
                peak_picker=params.get("profile_peak_picker", "hires"),
                signal_to_noise=params.get(
                    "profile_peak_picker_signal_to_noise", 1.0
                ),
            )
            file_path = str(prepared_path)
        with _tasks_lock:
            _tasks[task_id]["raw_input_file_path"] = raw_file_path
            _tasks[task_id]["analysis_file_path"] = file_path
            _tasks[task_id][
                "spectrum_representation_processing"
            ] = representation_processing

        # ── Phase 1: Parse spectra ──
        logger.info("[task %s] Phase 1/3: Parsing file %s ...", task_id, Path(file_path).name)
        with _tasks_lock:
            _tasks[task_id]["progress"] = 2
            _tasks[task_id]["message"] = "Parsing spectrum file..."
        _save_tasks()

        spectra = _parse_spectrum_file(
            file_path,
            isolation_window_width_da=params.get("isolation_window_width_da"),
            isolation_window_width_source=params.get(
                "isolation_window_width_source"
            ),
            allow_unknown_representation=(
                params.get("profile_peak_picking") == "assume_centroid"
            ),
        )
        total = len(spectra)
        t_parse = time.time() - t_start
        logger.info("[task %s] Phase 1 done: %d spectra parsed in %.1fs", task_id, total, t_parse)

        if total == 0:
            raise ValueError("No valid spectra found in file. Check file format and content.")

        detected_ion_mode, ion_mode_counts = _resolve_spectrum_polarities(
            spectra, params.get("ion_mode", "auto")
        )
        with _tasks_lock:
            _tasks[task_id]["ion_mode"] = detected_ion_mode
            _tasks[task_id]["ion_mode_counts"] = ion_mode_counts
        logger.info(
            "[task %s] Ion mode: %s (%s)", task_id, detected_ion_mode, ion_mode_counts
        )

        # ── Phase 1.25: continuous-MS1 feature workflow ────────────
        # This is deliberately separate from MS2 library peak filtering.
        # Feature IDs come from chromatographic MS1 mass traces; spectra that
        # cannot be assigned remain explicitly marked and later use the legacy
        # precursor-event fallback in the molecular network.
        feature_rows: list[dict[str, Any]] = []
        untargeted_workflow: dict[str, Any] = {
            "requested": bool(params.get("enable_ms1_feature_detection", True)),
            "status": "not_applicable_non_mzml",
            "input": file_path,
            "feature_count": 0,
        }
        if not params.get("enable_ms1_feature_detection", True):
            untargeted_workflow["status"] = "disabled"
        elif input_extension in {".mzml", ".mzxml", ".xml"}:
            if detected_ion_mode not in {"POS", "NEG"}:
                untargeted_workflow.update(
                    {
                        "status": "unavailable_mixed_polarity",
                        "warning": (
                            "OpenMS adduct grouping requires one ion mode per "
                            "feature map; split mixed-polarity acquisitions."
                        ),
                    }
                )
            else:
                try:
                    from utils.openms_backend import detect_untargeted_features

                    feature_timeout = int(
                        params.get("ms1_feature_timeout_seconds", 300)
                    )
                    with _tasks_lock:
                        _tasks[task_id]["progress"] = 3
                        _tasks[task_id]["message"] = (
                            "Detecting MS1 chromatographic features "
                            f"(active; bounded at {feature_timeout}s)..."
                        )
                    heartbeat_stop, heartbeat_thread = (
                        _start_task_stage_heartbeat(
                            task_id,
                            stage="openms_feature_detection",
                            message=(
                                "Detecting MS1 chromatographic features."
                            ),
                            timeout_seconds=feature_timeout,
                        )
                    )
                    try:
                        feature_payload = detect_untargeted_features(
                            file_path,
                            artifact_dir,
                            ion_mode=detected_ion_mode,
                            peak_picker=params.get(
                                "profile_peak_picker", "hires"
                            ),
                            signal_to_noise=params.get(
                                "profile_peak_picker_signal_to_noise", 1.0
                            ),
                            mass_error_ppm=params.get(
                                "ms1_mass_error_ppm", 10.0
                            ),
                            noise_threshold=params.get(
                                "ms1_noise_threshold_intensity", 1000.0
                            ),
                            chrom_peak_snr=params.get(
                                "ms1_chrom_peak_snr", 3.0
                            ),
                            min_trace_length=params.get(
                                "ms1_min_trace_length_seconds", 5.0
                            ),
                            min_fwhm=params.get(
                                "ms1_min_fwhm_seconds", 1.0
                            ),
                            max_fwhm=params.get(
                                "ms1_max_fwhm_seconds", 60.0
                            ),
                            remove_single_traces=params.get(
                                "ms1_remove_single_traces", False
                            ),
                            charge_low=params.get("ms1_charge_low", 1),
                            charge_high=params.get("ms1_charge_high", 1),
                            group_adducts=params.get(
                                "enable_adduct_grouping", True
                            ),
                            max_adduct_grouping_features=params.get(
                                "ms1_max_adduct_grouping_features", 10_000
                            ),
                            timeout_seconds=feature_timeout,
                        )
                    finally:
                        heartbeat_stop.set()
                        heartbeat_thread.join(timeout=1.0)
                    feature_rows = list(feature_payload.pop("features", []) or [])
                    attachment = attach_ms2_to_lcms_features(
                        spectra,
                        feature_rows,
                        mz_tolerance_ppm=params.get(
                            "ms2_feature_link_tolerance_ppm", 20.0
                        ),
                    )
                    untargeted_workflow = {
                        **feature_payload,
                        "requested": True,
                        "status": feature_payload.get("status", "ok"),
                        "feature_count": len(feature_rows),
                        "feature_table_storage": (
                            "OpenMS JSON manifest; omitted from task index"
                        ),
                        "ms2_attachment": attachment,
                    }
                    logger.info(
                        "[task %s] OpenMS features=%d; linked MS2=%d/%d",
                        task_id,
                        len(feature_rows),
                        attachment.get("linked_ms2_spectra", 0),
                        total,
                    )
                except Exception as exc:
                    logger.exception(
                        "[task %s] OpenMS feature workflow unavailable: %s",
                        task_id,
                        exc,
                    )
                    untargeted_workflow.update(
                        {
                            "status": "failed_with_explicit_ms2_fallback",
                            "error": str(exc),
                            "warning": (
                                "Continuous-MS1 feature detection failed. "
                                "Library matching continues, but unlinked MS2 "
                                "network nodes remain precursor-event fallbacks."
                            ),
                        }
                    )
        with _tasks_lock:
            _tasks[task_id]["untargeted_workflow"] = untargeted_workflow

        # ── Phase 1.5: TIC 全扫描概览(MS1+MS2 直读, ~1-2s) ──
        # 单独提取并落盘到 {task_id}_tic_overview.json, 前端按需懒加载, 不塞进 results.json
        try:
            from utils.mzml_io import extract_tic_overview as _extract_tic_overview
            tic_ov = _extract_tic_overview(file_path)
            _ext = Path(file_path).suffix.lower()
            if tic_ov is None and _ext == ".mgf":
                tic_ov = _build_mgf_tic_overview(spectra)  # MGF: RT 秒 → 分
            elif tic_ov is None and _ext in (".html", ".htm"):
                # raw_blank_analysis HTML 报告: 无连续 MS1, RT 已是分钟; 从 10 个 MS2 点
                # 重建近似 TIC (rt_in_seconds=False 不再 ÷60), 前端按 MS2 近似渲染。
                tic_ov = _build_mgf_tic_overview(
                    spectra, rt_in_seconds=False, source="html_report_ms2"
                )
            if tic_ov:
                tic_path = TASK_DIR / f"{task_id}_tic_overview.json"
                with open(tic_path, "w", encoding="utf-8") as _tf:
                    json.dump(tic_ov, _tf, separators=(",", ":"))
                logger.info("[task %s] TIC overview saved: %d MS1 + %d MS2 scans",
                            task_id, len(tic_ov.get("ms1", [])), len(tic_ov.get("ms2", [])))
        except Exception as _e:
            logger.warning("[task %s] TIC overview extraction failed: %s", task_id, _e)

        with _tasks_lock:
            _tasks[task_id]["total_spectra"] = total
            _tasks[task_id]["progress"] = 5
            _tasks[task_id]["message"] = f"Parsed {total} spectra (%.1fs). Starting matching..." % t_parse
        _save_tasks()

        # ── v3 Phase 1: 加载空白轮廓。用户明确请求后，任何加载错误都
        # 使任务失败；静默跳过会把“已扣空白”和“未扣空白”结果混在一起。
        blank_profile = None
        blank_profile_loaded = False
        blank_profiles: list[Any] = []
        blank_roles: list[str] = []
        blank_assessments = []
        blank_summary = {
            "requested": False,
            "feature_matches": 0,
            "excluded_spectra": 0,
            "retained_spectra": 0,
            "unevaluable_spectra": 0,
            "evidence_sources": {},
            "status_counts": {},
        }
        _blank_tokens = list(params.get("blank_files") or [])
        _legacy_blank = params.get("blank_file")
        if _legacy_blank and _legacy_blank not in _blank_tokens:
            _blank_tokens.append(_legacy_blank)
        _blank_tokens = list(dict.fromkeys(_blank_tokens))
        if _blank_tokens:
            from matching.blank_subtraction import (
                assign_blank_profiles_by_polarity,
                assess_blank_role_groups,
                blank_polarity_coverage,
                load_blank_file,
            )

            blank_paths = [_resolve_blank_file_token(token) for token in _blank_tokens]
            blank_role_map = dict(params.get("blank_roles") or {})
            blank_roles = [
                _normalize_blank_role(
                    blank_role_map.get(Path(str(token)).name)
                )
                for token in _blank_tokens
            ]
            blank_profiles = [load_blank_file(path) for path in blank_paths]
            for blank_profile in blank_profiles:
                if (not getattr(blank_profile, "context_spectra", None)
                        and not getattr(blank_profile, "ms1_scan_count", 0)):
                    raise ValueError(
                        "blank contains no usable MS1 features or precursor-context MS/MS spectra"
                    )
            polarity_coverage = blank_polarity_coverage(spectra, blank_profiles)
            blank_assignment = assign_blank_profiles_by_polarity(
                spectra, blank_profiles
            )
            requested_blank_paths = list(blank_paths)
            requested_blank_roles = list(blank_roles)
            requested_blank_profiles = list(blank_profiles)
            assigned_indices = blank_assignment["assigned_indices"]
            skipped_indices = blank_assignment["skipped_indices"]
            assigned_blank_paths = [
                requested_blank_paths[index] for index in assigned_indices
            ]
            skipped_blank_paths = [
                requested_blank_paths[index] for index in skipped_indices
            ]
            assigned_blank_roles = [
                requested_blank_roles[index] for index in assigned_indices
            ]
            unspecified_assigned_count = sum(
                role == "unspecified" for role in assigned_blank_roles
            )
            method_comparability_confirmed = bool(
                params.get("blank_method_comparability_confirmed", False)
            )
            exclusion_confirmed = bool(
                params.get("blank_background_exclusion_confirmed", False)
            )
            eligible_background_roles = {
                "background_blank", "solvent_blank", "method_blank",
                "extraction_blank",
            }
            roles_exclusion_eligible = bool(assigned_blank_roles) and all(
                role in eligible_background_roles
                for role in assigned_blank_roles
            )
            background_exclusion_eligible = bool(
                roles_exclusion_eligible
                and method_comparability_confirmed
                and exclusion_confirmed
            )
            if unspecified_assigned_count:
                blank_pairing_status = (
                    "polarity_only_unverified_method_pairing"
                )
                pairing_warning = (
                    "同极性空白包含未定角色文件，批次与色谱方法配对未经核验；"
                    "这些空白只产生背景警示，不据此剔除样本候选。"
                )
            elif not method_comparability_confirmed:
                blank_pairing_status = (
                    "role_and_polarity_labeled_method_unverified"
                )
                pairing_warning = (
                    "空白角色与极性已标注，但尚未确认与样品使用同批次／同色谱方法；"
                    "因此只作背景比较，不执行扣除。"
                )
            elif not exclusion_confirmed:
                blank_pairing_status = (
                    "method_comparable_background_exclusion_not_authorized"
                )
                pairing_warning = (
                    "方法可比性已确认，但未授权背景硬排除；空白仅作比较。"
                )
            else:
                blank_pairing_status = "explicitly_confirmed_background_pairing"
                pairing_warning = ""
            skipped_blank_details = [
                {
                    "filename": requested_blank_paths[index].name.split(
                        "_", 1
                    )[-1],
                    "polarities": blank_assignment["profile_polarities"][index],
                    "role": requested_blank_roles[index],
                }
                for index in skipped_indices
            ]
            blank_assignment_audit = {
                "blank_assignment_mode": blank_assignment["assignment_mode"],
                "requested_blank_filenames": [
                    path.name.split("_", 1)[-1]
                    for path in requested_blank_paths
                ],
                "requested_blank_replicate_count": len(
                    requested_blank_paths
                ),
                "requested_blank_roles": dict(
                    Counter(requested_blank_roles)
                ),
                "assigned_blank_filenames": [
                    path.name.split("_", 1)[-1]
                    for path in assigned_blank_paths
                ],
                "assigned_blank_replicate_count": len(
                    assigned_blank_paths
                ),
                "assigned_blank_polarities": blank_assignment[
                    "assigned_blank_polarities"
                ],
                "skipped_blank_filenames": [
                    path.name.split("_", 1)[-1]
                    for path in skipped_blank_paths
                ],
                "skipped_blank_replicate_count": len(skipped_blank_paths),
                "skipped_blank_polarities": blank_assignment[
                    "skipped_blank_polarities"
                ],
                "skipped_blank_details": skipped_blank_details,
                "unknown_assigned_blank_profiles": len(
                    blank_assignment["unknown_assigned_indices"]
                ),
                "blank_pairing_status": blank_pairing_status,
                "method_comparability_verified": (
                    method_comparability_confirmed
                ),
                "background_exclusion_authorized": exclusion_confirmed,
                "background_exclusion_eligible": background_exclusion_eligible,
                "unverified_assigned_blank_count": (
                    unspecified_assigned_count
                ),
                "blank_pairing_warning": pairing_warning,
            }
            if not polarity_coverage["applicable"]:
                sample_label = "/".join(
                    polarity_coverage["sample_polarities"]
                ) or "未知"
                blank_label = "/".join(
                    polarity_coverage["blank_polarities"]
                ) or "未知"
                mismatch_warning = (
                    f"空白过滤未执行：样本为 {sample_label}，上传的空白仅为 "
                    f"{blank_label}；缺少同极性空白。请上传 "
                    f"{'/'.join(polarity_coverage['missing_sample_polarities'])} "
                    "空白后重新分析。"
                )
                blank_summary = {
                    "requested": True,
                    "applicable": False,
                    "blank_filename": blank_paths[0].name.split("_", 1)[-1],
                    "blank_filenames": [
                        path.name.split("_", 1)[-1] for path in blank_paths
                    ],
                    "blank_replicate_count": len(blank_paths),
                    "blank_roles": dict(
                        Counter(blank_roles)
                    ),
                    **blank_assignment_audit,
                    "filter_policy": params.get(
                        "blank_filter_policy", DEFAULT_BLANK_FILTER_POLICY,
                    ),
                    "blank_operation": "comparison_not_applicable",
                    "cross_injection_msms_subtraction_performed": False,
                    "selected_blank_files_analyzed_as_samples": False,
                    "sample_polarities": polarity_coverage["sample_polarities"],
                    "blank_polarities": polarity_coverage["blank_polarities"],
                    "missing_sample_polarities": (
                        polarity_coverage["missing_sample_polarities"]
                    ),
                    "unknown_blank_profiles": (
                        polarity_coverage["unknown_blank_profiles"]
                    ),
                    "feature_matches": None,
                    "excluded_spectra": None,
                    "retained_spectra": None,
                    "unevaluable_spectra": total,
                    "status_counts": {
                        "blank_polarity_not_applicable": total,
                    },
                    "warning": mismatch_warning,
                }
                with _tasks_lock:
                    _tasks[task_id]["blank_summary"] = blank_summary
                _save_tasks()
                raise ValueError(mismatch_warning)

            # The browser deliberately sends all blanks selected for the
            # batch.  Once this sample's spectra have resolved POS/NEG
            # metadata, reduce that shared set to the matching task-specific
            # replicates before calculating prevalence or abundance.
            blank_paths = assigned_blank_paths
            blank_roles = assigned_blank_roles
            blank_profiles = [
                requested_blank_profiles[index] for index in assigned_indices
            ]
            blank_profile_loaded = bool(blank_profiles)
            logger.info(
                "[task %s] Automatic blank assignment: sample=%s, assigned=%s, "
                "skipped_opposite=%s, unknown_assigned=%d",
                task_id,
                "/".join(blank_assignment["sample_polarities"]),
                ",".join(path.name for path in blank_paths) or "(none)",
                ",".join(path.name for path in skipped_blank_paths) or "(none)",
                len(blank_assignment["unknown_assigned_indices"]),
            )

            requested_blank_policy = str(params.get(
                "blank_filter_policy", DEFAULT_BLANK_FILTER_POLICY,
            ))
            if (
                requested_blank_policy == "hard_legacy"
                and not background_exclusion_eligible
            ):
                raise ValueError(
                    "空白硬排除被拒绝：必须明确标注合格空白角色，并同时确认"
                    "同批次／同色谱方法可比性与背景排除授权。"
                )

            blank_assessments = assess_blank_role_groups(
                spectra,
                blank_profiles,
                blank_roles,
                sample_path=(
                    file_path if Path(file_path).suffix.lower()
                    in {".mzml", ".mzxml", ".xml"} else None
                ),
                sample_rt_unit=("seconds" if Path(file_path).suffix.lower() == ".mgf" else "minutes"),
                sample_to_blank_ratio=params.get(
                    "blank_sample_ratio", DEFAULT_BLANK_SAMPLE_RATIO,
                ),
                rt_tolerance_minutes=params.get(
                    "blank_rt_tolerance_minutes", DEFAULT_BLANK_RT_TOLERANCE_MINUTES,
                ),
                min_detection_fraction=params.get(
                    "blank_min_detection_fraction",
                    DEFAULT_BLANK_MIN_DETECTION_FRACTION,
                ),
                filter_policy=params.get(
                    "blank_filter_policy", DEFAULT_BLANK_FILTER_POLICY,
                ),
                clean_max_ratio=params.get(
                    "blank_clean_max_ratio", DEFAULT_BLANK_CLEAN_MAX_RATIO,
                ),
                clean_min_replicates=params.get(
                    "blank_clean_min_replicates",
                    DEFAULT_BLANK_CLEAN_MIN_REPLICATES,
                ),
            )
            status_counts = Counter(item.get("blank_status") for item in blank_assessments)
            category_counts = Counter(
                item.get("blank_evidence_category") for item in blank_assessments
                if item.get("blank_evidence_category")
            )
            evidence_sources = Counter(
                item.get("blank_evidence_source") for item in blank_assessments
                if item.get("blank_evidence_source")
            )
            feature_matches = sum(bool(item.get("blank_feature_matched")) for item in blank_assessments)
            excluded = sum(bool(item.get("blank_excluded")) for item in blank_assessments)
            clean_view_excluded = sum(
                bool(item.get("blank_clean_view_excluded"))
                for item in blank_assessments
            )
            flagged = sum(
                bool(item.get("blank_flagged")) for item in blank_assessments
            )
            carryover_suspected = sum(
                bool(item.get("blank_carryover_suspected"))
                for item in blank_assessments
            )
            retained = total - excluded
            unevaluable_statuses = {
                "sample_precursor_missing",
                "sample_rt_missing",
                "blank_polarity_not_applicable",
                "sample_intensity_missing",
                "blank_intensity_missing",
                "no_blank_ms2_context",
            }
            unevaluable = sum(
                int(count) for status, count in status_counts.items()
                if status in unevaluable_statuses
            )
            blank_filter_policy = str(params.get(
                "blank_filter_policy", DEFAULT_BLANK_FILTER_POLICY,
            ))
            blank_summary = {
                "requested": True,
                "applicable": True,
                "blank_filename": blank_paths[0].name.split("_", 1)[-1],
                "blank_filenames": [
                    path.name.split("_", 1)[-1] for path in blank_paths
                ],
                "blank_replicate_count": len(blank_paths),
                "blank_roles": dict(Counter(blank_roles)),
                **blank_assignment_audit,
                "filter_policy": blank_filter_policy,
                "blank_operation": (
                    "comparison_only"
                    if blank_filter_policy == "soft_evidence"
                    else "confirmed_presearch_background_exclusion"
                ),
                "cross_injection_msms_subtraction_performed": False,
                "selected_blank_files_analyzed_as_samples": False,
                "search_behavior": (
                    "identify_all_then_annotate"
                    if blank_filter_policy == "soft_evidence"
                    else "legacy_presearch_hard_filter"
                ),
                "sample_polarities": polarity_coverage["sample_polarities"],
                "blank_polarities": polarity_coverage["blank_polarities"],
                "unknown_blank_profiles": polarity_coverage[
                    "unknown_blank_profiles"
                ],
                "feature_matches": feature_matches,
                "excluded_spectra": excluded,
                "retained_spectra": retained,
                "searched_spectra": total - excluded,
                "flagged_spectra": flagged,
                "clean_view_excluded_spectra": clean_view_excluded,
                "carryover_suspected_spectra": carryover_suspected,
                "legacy_below_ratio_spectra": sum(
                    bool(item.get("blank_legacy_below_ratio"))
                    for item in blank_assessments
                ),
                "retained_without_blank_feature": status_counts.get(
                    "no_blank_feature", 0
                ),
                "unevaluable_spectra": unevaluable,
                "sample_to_blank_ratio_threshold": float(params.get(
                    "blank_sample_ratio", DEFAULT_BLANK_SAMPLE_RATIO,
                )),
                "clean_view_max_ratio": float(params.get(
                    "blank_clean_max_ratio", DEFAULT_BLANK_CLEAN_MAX_RATIO,
                )),
                "clean_view_min_replicates": int(params.get(
                    "blank_clean_min_replicates",
                    DEFAULT_BLANK_CLEAN_MIN_REPLICATES,
                )),
                "rt_tolerance_minutes": float(params.get(
                    "blank_rt_tolerance_minutes", DEFAULT_BLANK_RT_TOLERANCE_MINUTES,
                )),
                "min_detection_fraction": float(params.get(
                    "blank_min_detection_fraction",
                    DEFAULT_BLANK_MIN_DETECTION_FRACTION,
                )),
                "evidence_sources": dict(evidence_sources),
                "status_counts": dict(status_counts),
                "category_counts": dict(category_counts),
                "warning": (
                    "空白文件已读取，但在当前 m/z/RT/极性窗口内没有匹配到任何样品特征。"
                    if feature_matches == 0 else ""
                ),
            }
            with _tasks_lock:
                _tasks[task_id]["blank_summary"] = blank_summary
            logger.info(
                "[task %s] 空白证据评估: %s, policy=%s, 匹配=%d, 检索前排除=%d, "
                "清洁视图隐藏=%d, 标记=%d, source=%s",
                task_id, ",".join(path.name for path in blank_paths),
                blank_filter_policy, feature_matches, excluded,
                clean_view_excluded, flagged,
                dict(evidence_sources),
            )

        # ── Phase 2: Spectral matching ──
        logger.info("[task %s] Phase 2/3: Matching %d spectra...", task_id, total)
        deconvolution_run_pool = ([], {})
        if params.get(
            "enable_chimeric_deconvolution", USE_CHIMERIC_DECONVOLUTION
        ):
            deconvolution_run_pool = build_deconvolution_run_pool(spectra)
            logger.info(
                "[task %s] DecoID-like same-run clean component pool: %d spectra",
                task_id, len(deconvolution_run_pool[0]),
            )
        all_results = []
        peak_config = _peak_config_from_params(params)
        peak_audits: list[dict[str, Any]] = []
        n_matched = 0
        # utils.mgf_io 保留 RTINSECONDS 原单位；mzML/HTML 解析结果为分钟。
        _sample_rt_unit = (
            "seconds" if Path(file_path).suffix.lower() == ".mgf" else "minutes"
        )
        last_progress_update = time.time()
        last_log_update = time.time()
        last_progress_val = 5
        matching_start = time.time()

        for i, spec in enumerate(spectra):
            if cancel_event.is_set():
                raise TaskCancellationRequested(
                    f"Cancellation acknowledged before spectrum {i + 1}/{total}"
                )
            # Perform search
            t_search_start = time.time()
            # V7 默认先检索全部 LC-MS 特征，再把空白证据作为结果分层；
            # 仅显式 hard_legacy 兼容模式仍允许检索前排除。
            # 两种策略都不跨针直接相减未归一化的 MS/MS 碎片强度。
            peaks_for_search = spec["peaks"]
            blank_assessment = (
                blank_assessments[i] if i < len(blank_assessments)
                else {"blank_status": "not_requested", "blank_excluded": False}
            )
            spec_rt_unit = spec.get("retention_time_unit") or _sample_rt_unit
            precursor_mass_audit = resolve_search_precursor_mz(spec)
            search_precursor_mz = precursor_mass_audit.get(
                "search_precursor_mz"
            )
            if search_precursor_mz is None:
                search_precursor_mz = spec["precursor_mz"]
            query_peaks, peak_audit = process_peaks(
                peaks_for_search,
                spec["precursor_mz"],
                config=peak_config,
                precursor_charge=spec.get("precursor_charge"),
            )
            peak_audits.append(peak_audit)
            match = None if blank_assessment.get("blank_excluded") else search(
                peaks_for_search,
                search_precursor_mz,
                records=_lib_records,
                mz_index=_mz_index,
                cos_thresh=params.get("cosine_threshold", COSINE_THRESHOLD),
                final_thresh=params.get("final_score_threshold", FINAL_SCORE_THRESHOLD),
                mz_da=params.get("mz_tolerance_da", MZ_TOLERANCE_DA),
                mz_ppm=params.get("mz_tolerance_ppm", MZ_TOLERANCE_PPM),
                top_n=params.get("top_n_peaks", TOP_N_PEAKS),
                # HTML 摘要每谱仅有 <=5 个碎片。此类输入若仍使用宽松阈值，
                # 会把 cos~0.33/final~0.26/88ppm 的偶然近邻列成候选（实测
                # Arginine-Arginine 假线索）。有限碎片输入强制严格证据门槛；
                # 完整 mzML/MGF 仍尊重用户开关。
                strict_evidence_gate=(
                    True if spec.get("limited_fragment_list") else
                    params.get("strict_evidence_gate", USE_STRICT_EVIDENCE_GATE)
                ),
                strict_max_mass_ppm=(
                    LIMITED_FRAGMENT_STRICT_MAX_MASS_PPM
                    if spec.get("limited_fragment_list") else None
                ),
                query_adduct=spec.get("adduct"),
                query_adduct_reliable=spec.get("adduct_reliable"),
                query_polarity=spec.get("polarity"),
                query_rt=spec.get("retention_time"),   # v5: local RT standard confidence evidence
                query_rt_unit=spec_rt_unit,
                query_charge=spec.get("precursor_charge"),
                query_collision_energy=collision_energy_evidence(spec),
                query_instrument_type=spec.get("instrument_type"),
                query_instrument=spec.get("instrument"),
                query_activation=(
                    spec.get("activation_method") or spec.get("activation")
                ),
                query_collision_energy_mode=spec.get(
                    "collision_energy_mode"
                ),
                query_instrument_domain_id=spec.get(
                    "instrument_domain_id"
                ),
                query_collision_cell_region=spec.get(
                    "collision_cell_region"
                ),
                query_spectrum_representation=spec.get(
                    "spectrum_representation"
                ),
                query_raw_peaks_complete=(
                    spec.get("raw_peaks_complete") is True
                    and not bool(spec.get("limited_fragment_list"))
                ),
                query_rt_method_id=params.get("rt_method_id"),
                query_isotope_evidence=spec.get("precursor_isotope_evidence"),
                isotope_evidence_status=spec.get("isotope_evidence_status"),
                query_isolation_evidence=spec.get("precursor_isolation_evidence"),
                enable_chimeric_deconvolution=params.get(
                    "enable_chimeric_deconvolution", USE_CHIMERIC_DECONVOLUTION,
                ),
                deconvolution_run_pool=deconvolution_run_pool,
                peak_processing_config=peak_config,
            )
            if isinstance(match, dict):
                match["query_precursor_mass_audit"] = dict(
                    precursor_mass_audit
                )
                match["query_selected_ion_mz"] = precursor_mass_audit.get(
                    "selected_ion_mz"
                )
                match["query_search_precursor_mz"] = search_precursor_mz
                match["query_search_precursor_mz_source"] = (
                    precursor_mass_audit.get("search_precursor_mz_source")
                )
                match["query_search_precursor_delta_ppm"] = (
                    precursor_mass_audit.get("search_precursor_delta_ppm")
                )
            t_search = time.time() - t_search_start

            # Peak_Height=base peak(最高碎片强度), Peak_Area=TIC(总离子流,纯MS2近似色谱面积)
            _pks = spec.get("peaks") or []
            tic = sum(p[1] for p in _pks)
            base_peak = max((p[1] for p in _pks), default=0.0)
            all_results.append({
                "spectrum_id": spec.get("spectrum_id", i),
                "precursor_mz": spec["precursor_mz"],
                "search_precursor_mz": search_precursor_mz,
                "search_precursor_mz_source": precursor_mass_audit.get(
                    "search_precursor_mz_source"
                ),
                "search_precursor_delta_ppm": precursor_mass_audit.get(
                    "search_precursor_delta_ppm"
                ),
                "lcms_feature_used_for_search": bool(
                    precursor_mass_audit.get("lcms_feature_used_for_search")
                ),
                "precursor_mass_audit": precursor_mass_audit,
                "peaks_count": len(spec["peaks"]),
                "retention_time": spec.get("retention_time"),
                "retention_time_unit": spec_rt_unit,
                "ion_mode": spec.get("polarity"),
                "precursor_charge": spec.get("precursor_charge"),
                "precursor_intensity": spec.get("precursor_intensity"),
                "precursor_ms1_intensity": spec.get("precursor_ms1_intensity"),
                "collision_energy": spec.get("collision_energy"),
                "collision_energy_unit": spec.get("collision_energy_unit", ""),
                **query_acquisition_projection(spec),
                "instrument_type": spec.get("instrument_type", ""),
                "instrument": spec.get("instrument", ""),
                "rt_method_id": params.get("rt_method_id"),
                "precursor_scan_ref": spec.get("precursor_scan_ref", ""),
                "isotope_evidence_status": spec.get("isotope_evidence_status", "not_available"),
                "precursor_isotope_evidence": spec.get("precursor_isotope_evidence"),
                "precursor_isolation_evidence": spec.get(
                    "precursor_isolation_evidence"
                ),
                "input_format": spec.get("input_format"),
                "spectrum_representation": spec.get(
                    "spectrum_representation"
                ),
                "spectrum_representation_source": spec.get(
                    "spectrum_representation_source"
                ),
                "limited_fragment_list": bool(spec.get("limited_fragment_list", False)),
                "input_evidence_note": spec.get("input_evidence_note", ""),
                "tic": tic,
                "base_peak": base_peak,
                "match": match,
                "blank_used": bool(blank_profile_loaded),
                "blank_removed": 0,
                **blank_assessment,
                "peak_processing_audit": peak_audit,
                "query_peaks": query_peaks,
                # Keep the measured centroid list intact until serialization.
                # Search filtering is represented by query_peaks + its audit;
                # consensus/network processing must not inherit a hidden Top-15
                # truncation.
                "raw_peaks": [list(p) for p in _pks],
                **{
                    key: spec.get(key)
                    for key in (
                        "lcms_feature_id",
                        "lcms_feature_mz",
                        "lcms_feature_rt_minutes",
                        "lcms_feature_rt_start_minutes",
                        "lcms_feature_rt_end_minutes",
                        "lcms_feature_intensity",
                        "lcms_feature_charge",
                        "lcms_isotope_trace_count",
                        "lcms_adduct",
                        "lcms_neutral_mass",
                        "lcms_adduct_group_id",
                        "lcms_feature_link_status",
                        "lcms_feature_candidate_count",
                        "lcms_feature_mz_error_ppm",
                        "lcms_feature_rt_error_seconds",
                    )
                },
            })
            if match is not None:
                n_matched += 1

            # ── Progress tracking (robust) ──
            # Compute progress: 5% → 95% during matching phase
            new_progress = 5 + int((i + 1) / total * 90)
            now = time.time()

            # Update progress at least once per second, or when percentage changes
            should_update = (
                new_progress != last_progress_val or
                (now - last_progress_update) > 1.0
            )

            if should_update:
                matching_elapsed = now - matching_start
                rate = (
                    (i + 1) / matching_elapsed
                    if matching_elapsed > 0 else 0
                )
                eta = (total - i - 1) / rate if rate > 0 else 0

                with _tasks_lock:
                    _tasks[task_id]["progress"] = new_progress
                    _tasks[task_id]["message"] = (
                        f"Matching... {i+1}/{total} spectra "
                        f"({n_matched} matched, {rate:.1f} spec/s, ETA {eta:.0f}s)"
                    )
                last_progress_val = new_progress
                last_progress_update = now
                _save_task_progress_throttled(task_id)
                # SSE broadcast progress
                _sse_broadcast("progress", {
                    "task_id": task_id, "status": "running",
                    "progress": new_progress,
                    "message": _tasks[task_id]["message"],
                    "matched_spectra": n_matched, "total_spectra": total,
                })

            # Log every 10 seconds or every 100 spectra
            if (now - last_log_update) > 10.0 or (i + 1) % 100 == 0:
                elapsed = now - t_start
                matching_elapsed = now - matching_start
                rate = (
                    (i + 1) / matching_elapsed
                    if matching_elapsed > 0 else 0
                )
                logger.info(
                    "[task %s] Progress: %d/%d (%.0f%%) matched=%d rate=%.1f/s elapsed=%.0fs ETA=%.0fs t_search=%.2fs",
                    task_id, i + 1, total, (i + 1) / total * 100,
                    n_matched, rate, elapsed,
                    (total - i - 1) / rate if rate > 0 else 0,
                    t_search,
                )
                last_log_update = now

        # ── Phase 3: Serialize & Finalize ──
        logger.info("[task %s] Phase 3/3: Serializing results...", task_id)
        with _tasks_lock:
            if cancel_event.is_set():
                raise TaskCancellationRequested(
                    "Cancellation acknowledged before finalization"
                )
            _transition_task(
                _tasks[task_id],
                TaskState.FINALIZING,
                reason="spectral_matching_completed",
            )
            _tasks[task_id]["progress"] = 96
            _tasks[task_id]["message"] = f"Serializing {n_matched} results..."
            # The scan-level match count is already final at this point.  Keep
            # it visible while additive entity/claim artifacts are built so
            # the UI does not misleadingly display 0/total during finalization.
            _tasks[task_id]["matched_spectra"] = n_matched
            _tasks[task_id]["matched_count"] = n_matched
        _save_tasks()

        matched = [r for r in all_results if r.get("match")]
        # Select a candidate-independent representative from the complete
        # measured MS/MS arrays. True consensus remains disabled until the
        # required precursor-purity and EIC-quality weights are available.
        # Public result rows are intentionally capped for browser payload size,
        # so consensus construction must happen before serialization.
        consensus_artifacts, consensus_spectrum_summary = (
            build_lcms_consensus_spectra(all_results)
        )
        serialized = _serialize_results(all_results)
        reliability_annotation_summary = {
            "policy": "mssystem-v16-real-sample-reliability-1.0.0",
            "row_count": len(serialized),
            "fragmentation_status_counts": dict(Counter(
                str(row.get("fragmentation_status") or "unknown")
                for row in serialized
            )),
            "ce_comparability_counts": dict(Counter(
                str(
                    ((row.get("match") or {}).get("ce_comparability") or {})
                    .get("status")
                    or "no_candidate"
                )
                for row in serialized
            )),
            "background_B3_counts": dict(Counter(
                str((row.get("background_evidence") or {}).get(
                    "B3_blank_confirmed"
                ) or "unknown")
                for row in serialized
            )),
            "new_background_excluded_peak_count": sum(
                int((row.get("background_evidence") or {}).get(
                    "new_background_excluded_peak_count"
                ) or 0)
                for row in serialized
            ),
            "new_background_soft_downweighted_peak_count": sum(
                int((row.get("background_evidence") or {}).get(
                    "new_background_soft_downweighted_peak_count"
                ) or 0)
                for row in serialized
            ),
            "ranking_influence": False,
            "governance_influence": False,
        }
        peak_processing_summary = aggregate_peak_processing_audits(peak_audits)
        consensus_path = artifact_dir / "lcms_feature_consensus_spectra.json"
        consensus_path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "summary": consensus_spectrum_summary,
                    "spectra": consensus_artifacts,
                },
                ensure_ascii=False,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        untargeted_workflow["consensus_spectra"] = {
            **consensus_spectrum_summary,
            "artifact": str(consensus_path),
        }
        if params.get(
            "enable_ion_identity_networking", USE_ION_IDENTITY_NETWORKING
        ):
            ion_identity_summary = annotate_ion_identity_networks(
                serialized,
                file_path=file_path,
                max_rt_delta_min=ION_IDENTITY_RT_TOLERANCE_MIN,
                min_shape_correlation=ION_IDENTITY_MIN_SHAPE_CORRELATION,
                min_overlap_fraction=ION_IDENTITY_MIN_OVERLAP,
            )
            # IIMN conflicts are confidence evidence.  Re-run the single
            # presentation calibration source after graph annotation; raw
            # scores/ranking remain unchanged.
            serialized = _prepare_results_for_display(serialized)
        else:
            ion_identity_summary = {
                "method": "Ion Identity Molecular Networking",
                "status": "disabled",
                "network_count": 0,
                "networked_features": 0,
                "edge_count": 0,
                "propagated_annotations": 0,
                "runtime_influence": False,
            }

        # ── Independent Expected Compounds layer ───────────────────
        # This happens after ordinary MS/MS matching so associated measured
        # MS2 rows can add diagnostic/library evidence.  Its MS1 candidates
        # remain in a separate artifact and never enter `serialized`,
        # `matched`, match_rate, or unique_compounds.
        with _tasks_lock:
            _tasks[task_id]["progress"] = 97
            _tasks[task_id]["message"] = (
                "Screening Expected Compounds across continuous MS1..."
            )
        _save_tasks()
        try:
            _, expected_compounds_summary = _run_expected_compounds_layer(
                task_id,
                file_path=file_path,
                detected_ion_mode=detected_ion_mode,
                params=params,
                serialized_ms2_rows=serialized,
                openms_features=feature_rows,
                blank_profiles=blank_profiles,
                blank_roles=blank_roles,
            )
        except Exception as exc:
            logger.exception(
                "[task %s] Expected Compounds layer failed: %s",
                task_id,
                exc,
            )
            failure_payload = _expected_compounds_unavailable_payload(
                status="failed",
                input_path=file_path,
                warning=(
                    "Expected Compounds failed independently; ordinary "
                    "spectral-library results remain valid and unchanged."
                ),
            )
            failure_payload.update(
                {
                    "task_id": task_id,
                    "error": f"{type(exc).__name__}: {exc}"[:1000],
                }
            )
            try:
                write_expected_compounds_artifact(
                    _expected_compounds_artifact_path(task_id),
                    failure_payload,
                )
            except Exception:
                logger.exception(
                    "[task %s] Could not persist Expected Compounds failure artifact",
                    task_id,
                )
            expected_compounds_summary = {
                **dict(failure_payload.get("summary") or {}),
                "status": "failed",
                "requested": bool(
                    params.get("enable_expected_compounds", False)
                ),
                "available": False,
                "artifact": str(
                    _expected_compounds_artifact_path(task_id)
                ),
                "error": failure_payload["error"],
                "counted_in_standard_spectral_library_matches": False,
                "standard_spectral_library_match_increment": 0,
                "expected_provider_call_count": int(bool(
                    params.get("enable_expected_compounds", False)
                )),
                "evaluation_panel_access_by_pipeline": 0,
            }
        runtime_reliability_contract = build_reliability_contract(
            target_defaults={
                "target_screening_enabled": False,
                "target_compounds_ref": None,
            },
            expected_defaults={
                "enable_expected_compounds": bool(
                    params.get("enable_expected_compounds", False)
                ),
                "expected_compound_datasets": list(
                    params.get("expected_compound_datasets") or []
                ),
            },
            flagship_policy=_load_v15_flagship_policy(),
        )
        runtime_reliability_contract["runtime_observations"] = {
            "expected_provider_call_count": int(
                expected_compounds_summary.get("expected_provider_call_count", 0)
            ),
            "evaluation_panel_access_by_pipeline": int(
                expected_compounds_summary.get(
                    "evaluation_panel_access_by_pipeline", 0
                )
            ),
            "target_screening_provider_call_count": 0,
            "ordinary_search_completed_before_expected_layer": True,
            "ordinary_candidate_scores_changed_by_expected_layer": False,
        }
        reliability_contract_path = artifact_dir / "reliability_contract.json"
        reliability_contract_path.write_text(
            json.dumps(
                runtime_reliability_contract,
                ensure_ascii=False,
                indent=2,
            ) + "\n",
            encoding="utf-8",
        )
        reliability_annotation_summary["runtime_contract"] = (
            runtime_reliability_contract
        )
        reliability_annotation_summary["runtime_contract_artifact"] = str(
            reliability_contract_path
        )
        _matched_serialized = [r for r in serialized if r.get("match")]
        if blank_summary.get("requested"):
            clean_view_rows = [
                row for row in serialized
                if not bool(row.get("blank_clean_view_excluded"))
            ]
            clean_view_matched = [
                row for row in clean_view_rows if row.get("match")
            ]
            blank_summary.update({
                "all_view_matched_spectra": len(_matched_serialized),
                "clean_view_matched_spectra": len(clean_view_matched),
                "all_view_unique_candidates": len(
                    _dedup_matched(_matched_serialized)
                ),
                "clean_view_unique_candidates": len(
                    _dedup_matched(clean_view_matched)
                ),
                "flagged_matched_spectra": sum(
                    bool(row.get("blank_flagged")) and bool(row.get("match"))
                    for row in serialized
                ),
            })
        unique_compounds = len(_dedup_matched(_matched_serialized))
        feature_summary = _precursor_feature_summary(serialized)
        chimeric_spectra = sum(
            str(((row.get("match") or {}).get("chimeric_deconvolution") or {})
                .get("isolation_status") or "") in {"mixed", "severe"}
            for row in serialized if isinstance(row, dict)
        )
        deconvolved_spectra = sum(
            bool(((row.get("match") or {}).get("chimeric_deconvolution") or {})
                 .get("runtime_influence"))
            for row in serialized if isinstance(row, dict)
        )
        ion_identity_summary = _effective_ion_identity_summary(
            ion_identity_summary
        )
        compound_confidence_summary = _compound_confidence_summary(
            _matched_serialized
        )
        identification_funnel = _identification_funnel_summary(serialized)
        deconvolution_impact = _deconvolution_impact_summary(serialized)
        evidence_governance_summary = _summarize_evidence_governance(
            serialized
        )
        module_availability = _analysis_module_availability(
            params, ion_identity_summary, serialized
        )
        contaminant_unique_compounds = len(_dedup_matched([
            r for r in _matched_serialized
            if r["match"].get("potential_contaminant") is True
        ]))

        with _tasks_lock:
            _tasks[task_id]["progress"] = 99
            _tasks[task_id]["message"] = "Saving results..."

        _save_task_results_given(task_id, serialized)

        # Retain the frozen V11/V14-compatible decision as an audit comparator.
        # V15 is materialized separately and never overwrites either this audit
        # artifact or the scan-level evidence.
        v11_flagship_summary: dict[str, Any]
        try:
            with _tasks_lock:
                _tasks[task_id]["progress"] = 98
                _tasks[task_id]["message"] = (
                    "Building frozen audit entity ranking..."
                )
            _save_tasks()
            v11_flagship_summary = _write_v11_flagship_artifact(
                task_id,
                serialized,
                _v11_flagship_artifact_path(TASK_DIR, task_id),
                runtime_snapshot_id=_current_runtime_snapshot_id(),
            )
        except Exception as v11_flagship_exc:
            logger.exception(
                "[task %s] V11 audit flagship materialization failed: %s",
                task_id,
                v11_flagship_exc,
            )
            v11_flagship_summary = {
                "status": "failed",
                "error": (
                    f"{type(v11_flagship_exc).__name__}: {v11_flagship_exc}"
                )[:1000],
                "production_identity_selector": False,
                "legacy_result_mutated": False,
                "entity_count": 0,
            }

        # V15 primary result: every same-identity library reference is retained,
        # mirrored peak content is collapsed, and CE/instrument-compatible
        # provenance strata enter the robust scan score.
        flagship_summary: dict[str, Any]
        try:
            with _tasks_lock:
                _tasks[task_id]["progress"] = 98
                _tasks[task_id]["message"] = (
                    "Building primary entity ranking with bounded candidate recall..."
                )
            _save_tasks()
            ce_evidence_retriever = _v15_ce_evidence_retriever()
            flagship_summary = _write_flagship_artifact(
                task_id,
                serialized,
                _flagship_artifact_path(TASK_DIR, task_id),
                runtime_snapshot_id=_current_runtime_snapshot_id(),
                ce_evidence_retriever=ce_evidence_retriever,
                candidate_recall_retriever=ce_evidence_retriever,
                detection_method_id=str(
                    params.get("rt_method_id") or "unspecified"
                ),
            )
        except Exception as flagship_exc:
            logger.exception(
                "[task %s] V15 flagship entity ranking failed: %s",
                task_id,
                flagship_exc,
            )
            flagship_summary = {
                "status": "failed",
                "error": f"{type(flagship_exc).__name__}: {flagship_exc}"[:1000],
                "production_identity_selector": True,
                "legacy_result_mutated": False,
                "entity_count": 0,
            }

        # P1 compatibility artifact: one versioned result object for future
        # API/report consumers.  It is additive and never mutates the legacy
        # scan evidence; V11 identity selection is the flagship artifact above.
        identification_contract_summary: dict[str, Any]
        try:
            with _tasks_lock:
                _tasks[task_id]["progress"] = 99
                _tasks[task_id]["message"] = (
                    "Saving identification and evidence sidecars..."
                )
            _save_tasks()
            identification_contract_summary = _write_identification_contract(
                task_id,
                serialized,
                _identification_contract_path(TASK_DIR, task_id),
                runtime_snapshot_id=_current_runtime_snapshot_id(),
            )
            identification_contract_summary["status"] = "available"
        except Exception as contract_exc:
            logger.exception(
                "[task %s] Identification contract sidecar failed: %s",
                task_id,
                contract_exc,
            )
            identification_contract_summary = {
                "status": "failed",
                "error": f"{type(contract_exc).__name__}: {contract_exc}"[:1000],
                "legacy_result_mutated": False,
                "production_ranking_changed_records": 0,
            }

        # P2 candidate-independent entity graph.  Entity construction is
        # completed before candidate aggregation; the latter remains a shadow
        # view and cannot mutate the production scan-level result or Top-1.
        entity_bundle_summary: dict[str, Any]
        try:
            entity_bundle_summary = _write_entity_bundle(
                task_id,
                serialized,
                _entity_bundle_path(TASK_DIR, task_id),
            )
        except Exception as entity_exc:
            logger.exception(
                "[task %s] Entity bundle sidecar failed: %s",
                task_id,
                entity_exc,
            )
            entity_bundle_summary = {
                "status": "failed",
                "error": f"{type(entity_exc).__name__}: {entity_exc}"[:1000],
                "candidate_independent": True,
                "production_ranking_changed": False,
                "feature_consensus_runtime_influence": False,
            }

        # P6 final claim outlet.  A single task supplies only one independent
        # injection group, and no project q-value is fabricated.  Therefore
        # this compiler is expected to withhold ordinary releases until a
        # separately declared, calibrated cohort is available.
        claim_contract_summary: dict[str, Any]
        try:
            claim_contract_summary = _write_claim_contract(
                task_id,
                _iter_identification_records(
                    _identification_contract_path(TASK_DIR, task_id)
                ),
                _claim_contract_path(TASK_DIR, task_id),
                context=_default_task_claim_context(task_id),
            )
            claim_contract_summary["status"] = "available"
        except Exception as claim_exc:
            logger.exception(
                "[task %s] EviGate 2 claim contract failed: %s",
                task_id,
                claim_exc,
            )
            claim_contract_summary = {
                "status": "failed",
                "error": f"{type(claim_exc).__name__}: {claim_exc}"[:1000],
                "released_count": 0,
                "production_ranking_changed_records": 0,
                "legacy_confidence_overwritten_records": 0,
            }

        # Ordinary-user result: one row per deduplicated 2D identity, or one
        # unresolved molecular-formula isomer group.  It is generated only
        # after the formal claim contract exists, so "confirmed" can never be
        # inferred from a library score alone.
        compound_conclusion_summary: dict[str, Any]
        try:
            result_artifact = TASK_DIR / f"{task_id}_results.json"
            primary_artifact = _flagship_artifact_path(TASK_DIR, task_id)
            inventory_manifest = (
                Path(str(untargeted_workflow.get("manifest"))).resolve()
                if untargeted_workflow.get("manifest")
                else None
            )
            projection_feature_rows = [
                {
                    **row,
                    "ion_mode": row.get("ion_mode") or detected_ion_mode,
                }
                for row in feature_rows
            ] if inventory_manifest is not None else None
            if not _flagship_artifact_is_current(primary_artifact):
                raise RuntimeError(
                    "Current V15 primary identity artifact is unavailable"
                )
            compound_conclusion_summary = _write_compound_conclusion_artifact(
                task_id,
                serialized,
                _compound_conclusion_artifact_path(TASK_DIR, task_id),
                claim_decisions_path=_claim_contract_path(TASK_DIR, task_id),
                source_artifact=result_artifact,
                primary_identity_artifact=primary_artifact,
                feature_inventory=projection_feature_rows,
                feature_inventory_source=inventory_manifest,
            )
        except Exception as conclusion_exc:
            logger.exception(
                "[task %s] Compound-level user conclusion failed: %s",
                task_id,
                conclusion_exc,
            )
            compound_conclusion_summary = {
                "status": "failed",
                "error": (
                    f"{type(conclusion_exc).__name__}: {conclusion_exc}"
                )[:1000],
                "confirmed_present_count": 0,
                "main_result_count": 0,
            }

        # Mixed-standard truth is intentionally absent from the ordinary
        # analysis contract.  A completed blind-search task can be audited
        # later only through the dedicated reference-standard QC endpoint.
        _reference_standard_qc_artifact_path(
            TASK_DIR, task_id
        ).unlink(missing_ok=True)
        reference_standard_qc_summary: dict[str, Any] = {
            "status": "not_requested",
            "set_id": None,
            "explicit_posthoc_endpoint_required": True,
            "filename_inference_used": False,
            "used_for_candidate_generation": False,
            "used_for_scoring": False,
            "used_for_ranking": False,
            "production_results_mutated": False,
        }

        primary_artifact_current = False
        primary_artifact_check_error = None
        if flagship_summary.get("status") == "available":
            try:
                primary_artifact_current = _flagship_artifact_is_current(
                    _flagship_artifact_path(TASK_DIR, task_id))
            except (OSError, ValueError, RuntimeError, TypeError, KeyError) as exc:
                primary_artifact_check_error = f"{type(exc).__name__}: {exc}"
        completion_check = formal_output_check({
            "primary_identity": flagship_summary,
            "identification_records": identification_contract_summary,
            "claim_contract": claim_contract_summary,
            "compound_conclusions": compound_conclusion_summary,
        }, primary_artifact_current=primary_artifact_current)
        if primary_artifact_check_error:
            completion_check["primary_artifact_check_error"] = primary_artifact_check_error
        # Preserve the failed modules and committed raw evidence before the
        # outer exception handler records FAILED. An available V11 audit is
        # deliberately not an alternative to the formal identity artifact.
        with _tasks_lock:
            _tasks[task_id].update({
                "flagship_summary": flagship_summary,
                "v11_flagship_audit_summary": v11_flagship_summary,
                "identification_contract_summary": identification_contract_summary,
                "claim_contract_summary": claim_contract_summary,
                "compound_conclusion_summary": compound_conclusion_summary,
                "analysis_completion_check": completion_check,
            })
        _save_tasks()
        assert_formal_outputs(completion_check)

        t_total = time.time() - t_start
        match_rate = len(matched) / total * 100 if total else 0
        logger.info(
            "[task %s] COMPLETED: %d/%d matched (%.1f%%) in %.1fs",
            task_id, len(matched), total, match_rate, t_total,
        )

        with _tasks_lock:
            _tasks[task_id]["progress"] = 100
            _tasks[task_id]["result_run_id"] = _tasks[task_id].get(
                "active_run_id"
            )
            _tasks[task_id]["matched_spectra"] = len(matched)
            _tasks[task_id]["matched_count"] = len(matched)  # 别名, 确保前端/统计读取一致
            # Preserve the scan-level candidate count as an audit diagnostic.
            # The public compound count is replaced below by the scope-aware
            # final table, which excludes unresolved precursor events.
            _tasks[task_id]["legacy_unique_library_candidates"] = (
                unique_compounds
            )
            _tasks[task_id]["unique_compounds"] = unique_compounds
            _tasks[task_id][
                "compound_confidence_summary"
            ] = compound_confidence_summary
            _tasks[task_id]["identification_funnel"] = identification_funnel
            _tasks[task_id]["deconvolution_impact"] = deconvolution_impact
            _tasks[task_id][
                "evidence_governance_summary"
            ] = evidence_governance_summary
            _tasks[task_id][
                "identification_contract_summary"
            ] = identification_contract_summary
            _tasks[task_id]["flagship_summary"] = flagship_summary
            _tasks[task_id]["v11_flagship_audit_summary"] = (
                v11_flagship_summary
            )
            _tasks[task_id]["entity_bundle_summary"] = entity_bundle_summary
            _tasks[task_id]["claim_contract_summary"] = claim_contract_summary
            _tasks[task_id]["compound_conclusion_summary"] = (
                compound_conclusion_summary
            )
            _apply_final_compound_count_summary(
                _tasks[task_id], compound_conclusion_summary
            )
            _tasks[task_id]["reference_standard_qc_summary"] = (
                reference_standard_qc_summary
            )
            _tasks[task_id]["presentation_calibration_version"] = (
                PRESENTATION_CALIBRATION_VERSION
            )
            _tasks[task_id]["precursor_features"] = feature_summary["precursor_features"]
            _tasks[task_id]["matched_precursor_features"] = feature_summary[
                "matched_precursor_features"
            ]
            _tasks[task_id]["high_precursor_features"] = feature_summary[
                "high_precursor_features"
            ]
            _tasks[task_id]["feature_summary"] = feature_summary
            _tasks[task_id]["peak_processing_summary"] = (
                peak_processing_summary
            )
            _tasks[task_id]["untargeted_workflow"] = untargeted_workflow
            _tasks[task_id]["consensus_spectrum_summary"] = (
                consensus_spectrum_summary
            )
            _tasks[task_id]["expected_compounds_summary"] = (
                expected_compounds_summary
            )
            _tasks[task_id]["reliability_annotation_summary"] = (
                reliability_annotation_summary
            )
            boundary_audit = dict(
                _tasks[task_id].get("pipeline_boundary_audit") or {}
            )
            boundary_audit.update({
                "expected_provider_call_count": int(
                    expected_compounds_summary.get(
                        "expected_provider_call_count", 0
                    ) or 0
                ),
                "evaluation_panel_access_by_pipeline": 0,
            })
            _tasks[task_id]["pipeline_boundary_audit"] = boundary_audit
            _tasks[task_id]["ion_identity_summary"] = ion_identity_summary
            _tasks[task_id]["ion_identity_networks"] = int(
                ion_identity_summary.get("network_count") or 0
            )
            _tasks[task_id]["ion_identity_networked_features"] = int(
                ion_identity_summary.get("networked_features") or 0
            )
            _tasks[task_id]["ion_identity_linked_annotations"] = int(
                ion_identity_summary.get("propagated_annotations") or 0
            )
            _tasks[task_id]["module_availability"] = module_availability
            _tasks[task_id]["ms1_scans"] = int(
                module_availability.get("ms1_scan_count") or 0
            )
            _tasks[task_id]["analysis_warnings"] = module_availability["warnings"]
            _tasks[task_id]["chimeric_spectra"] = chimeric_spectra
            _tasks[task_id]["deconvolved_spectra"] = deconvolved_spectra
            _tasks[task_id]["deconvolution_run_component_spectra"] = len(
                deconvolution_run_pool[0]
            )
            _tasks[task_id]["contaminant_unique_compounds"] = contaminant_unique_compounds
            _tasks[task_id]["match_rate"] = round(match_rate, 1)
            # The complete result artifact has already been committed by
            # ``_save_task_results_given``.  Keeping ``serialized`` here made
            # every completed sample permanently resident in the API process;
            # a few hundred LC-MS tasks could therefore grow the server to
            # tens of gigabytes and stall every request in garbage collection.
            _tasks[task_id]["results"] = []
            _tasks[task_id]["blank_summary"] = blank_summary
            _transition_task(
                _tasks[task_id],
                TaskState.COMPLETED,
                reason="formal_identity_and_required_outputs_verified",
            )
            _tasks[task_id]["message"] = (
                f"{_task_message_for_display(_tasks[task_id])} in {t_total:.0f}s"
            )

        # 将最终 completed 状态持久化到 tasks_index.json。
        # results 已在上方 _save_task_results_given 落盘;
        # 但若不在此再 save 一次元数据, 一旦重启, tasks_index.json 里状态会
        # 退回上次保存的 running (即使 results 已安全落盘) —— 见 abe5e53e3273。
        _save_tasks()

        # Literature lookup is a separate, rate-limited I/O worker.  It starts
        # only after the analytical result is safely persisted and cannot
        # delay the Full-library serial queue.
        _maybe_launch_isomer_literature_enrichment(task_id)

        # EvoGraph-F is a read-only research sidecar.  It is deliberately
        # launched after completion and on an executor isolated from ranking.
        _maybe_launch_evograph_fragmentation(task_id)

        # CFM 后台预计算: Docker 可用时非阻塞预灌碎裂缓存, 报告/详情随后可直显 CFM 产物级联
        _maybe_launch_cfm_precompute(task_id)

        # SSE broadcast completion
        _sse_broadcast("status_change", {
            "task_id": task_id, "status": "completed", "progress": 100,
            "matched_spectra": len(matched), "total_spectra": total,
            "unique_compounds": _tasks[task_id].get("unique_compounds"),
            "legacy_unique_library_candidates": _tasks[task_id].get(
                "legacy_unique_library_candidates"
            ),
            "main_result_count": _tasks[task_id].get("main_result_count"),
            "result_record_count": _tasks[task_id].get("result_record_count"),
            "unresolved_precursor_event_count": _tasks[task_id].get(
                "unresolved_precursor_event_count"
            ),
            "compound_count_semantics": _tasks[task_id].get(
                "compound_count_semantics"
            ),
            "precursor_features": feature_summary["precursor_features"],
            "matched_precursor_features": feature_summary["matched_precursor_features"],
            "compound_confidence_summary": compound_confidence_summary,
            "identification_funnel": identification_funnel,
            "deconvolution_impact": deconvolution_impact,
            "evidence_governance_summary": evidence_governance_summary,
            "identification_contract_summary": identification_contract_summary,
            "flagship_summary": flagship_summary,
            "analysis_completion_check": completion_check,
            "v11_flagship_audit_summary": v11_flagship_summary,
            "entity_bundle_summary": entity_bundle_summary,
            "claim_contract_summary": claim_contract_summary,
            "reference_standard_qc_summary": reference_standard_qc_summary,
            "module_availability": module_availability,
            "peak_processing_summary": peak_processing_summary,
            "untargeted_workflow": untargeted_workflow,
            "consensus_spectrum_summary": consensus_spectrum_summary,
            "expected_compounds_summary": expected_compounds_summary,
            "reliability_annotation_summary": (
                reliability_annotation_summary
            ),
            "analysis_warnings": module_availability["warnings"],
            "ion_identity_networks": int(ion_identity_summary.get("network_count") or 0),
            "ion_identity_networked_features": int(
                ion_identity_summary.get("networked_features") or 0
            ),
            "chimeric_spectra": chimeric_spectra,
            "deconvolved_spectra": deconvolved_spectra,
            "contaminant_unique_compounds": contaminant_unique_compounds,
            "filename": task["filename"],
            "message": _tasks[task_id]["message"],
        })

        with _tasks_lock:
            _task_cancel_events.pop(task_id, None)

    except TaskCancellationRequested as exc:
        logger.info("[task %s] CANCELLED: %s", task_id, exc)
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current is not None:
                if current.get("status") == TaskState.RUNNING.value:
                    _transition_task(
                        current,
                        TaskState.CANCELLING,
                        reason="worker_observed_cancel_event",
                    )
                _transition_task(
                    current,
                    TaskState.CANCELLED,
                    reason="worker_stopped_before_artifact_commit",
                )
                current["message"] = "Analysis cancelled safely before result commit"
                current["error"] = None
                current["results"] = []
            _task_cancel_events.pop(task_id, None)
        _save_tasks()
        _sse_broadcast(
            "status_change",
            {
                "task_id": task_id,
                "status": "cancelled",
                "message": "Analysis cancelled safely before result commit",
            },
        )
    except Exception as e:
        logger.exception("[task %s] FAILED: %s", task_id, e)
        with _tasks_lock:
            current = _tasks[task_id]
            _transition_task(
                current,
                TaskState.FAILED,
                reason=f"analysis_exception:{type(e).__name__}",
            )
            current["error"] = str(e)
            current["message"] = f"Error: {e}"
            _task_cancel_events.pop(task_id, None)
        _save_tasks()
        # SSE broadcast failure
        _sse_broadcast("status_change", {
            "task_id": task_id, "status": "failed",
            "error": str(e), "message": f"Error: {e}",
        })


def _msms_fragments(query_peaks, lib_peaks, tol=0.1, top_n=20,
                    query_pmz=None, lib_pmz=None):
    """Build a report string from the same direct/shifted pairs as scoring."""
    if not query_peaks or not lib_peaks:
        return ""
    pairs = find_matching_fragments(
        query_peaks,
        lib_peaks,
        tol=tol,
        query_pmz=query_pmz,
        lib_pmz=lib_pmz,
    )
    matched = []
    seen = set()
    for pair in pairs:
        value = round(float(pair["query_mz"]), 4)
        if value in seen:
            continue
        seen.add(value)
        matched.append(f"{value:.4f}")
    return "; ".join(matched[:top_n]) if matched else ""


def _resolve_formula(smiles, raw_formula):
    """分子式解析: 优先用库原始 formula(合理时), 否则用 RDKit 从 SMILES 计算。
    修复库数据问题: 部分化合物 formula 缺失(列表显示'无分子式'但详情能画结构式)
    或残缺(如 'C5' 缺 H)——这些是源库本身的数据缺陷, 用 SMILES 结构重算补全。"""
    rf = str(raw_formula or "").strip()
    if rf and 'H' in rf and len(rf) >= 4:   # 合理分子式含H且够长(如 C6H12O6)
        return rf
    if smiles:
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is not None:
                calc = rdMolDescriptors.CalcMolFormula(mol)
                if calc:
                    return str(calc)
        except Exception:
            pass
    return rf


def _is_nan(v):
    """float('nan') 检测。nan 是 truthy, 不能用 `not v` 或 `v or x` 判断;
    必须用 `v != v`(nan 唯一性质)。用于清理库/结果中的 NaN 字段。"""
    try:
        return v != v
    except Exception:
        return False


def _result_group_key(m: dict) -> str:
    """Stable compound key for sample-level consensus support."""
    return _compound_identity_key(m)


def _float_or_none(v) -> Optional[float]:
    try:
        f = float(v)
        return None if f != f else f
    except (TypeError, ValueError):
        return None


def _median(vals: list[float]) -> Optional[float]:
    if not vals:
        return None
    vals = sorted(vals)
    n = len(vals)
    mid = n // 2
    if n % 2:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2


def _apply_sample_consensus(results: list[dict]) -> list[dict]:
    """Promote repeated same-sample hits to Strong Candidate when evidence is coherent.

    This is a presentation/evidence-tier upgrade only. It does not change the search
    ranking or promote compound-level Medium to High.
    """
    if not USE_SAMPLE_CONSENSUS_TIER:
        return results

    groups: dict[str, list[dict]] = defaultdict(list)
    for r in results:
        m = r.get("match")
        if not m or not m.get("has_structure"):
            continue
        groups[_result_group_key(m)].append(r)

    # One library candidate can occur at several unrelated RT features.  The
    # old all-or-nothing span check discarded every valid local cluster as soon
    # as that name also appeared elsewhere in the chromatogram.  Split by a
    # bounded RT span and evaluate each chromatographic cluster independently.
    rt_clusters: list[list[dict]] = []
    for rows in groups.values():
        timed = [(_float_or_none(r.get("retention_time")), r) for r in rows]
        timed = [(rt, r) for rt, r in timed if rt is not None]
        timed.sort(key=lambda item: item[0])
        current: list[tuple[float, dict]] = []
        for rt, row in timed:
            if current and rt - current[0][0] > SAMPLE_CONSENSUS_RT_WINDOW_MIN:
                rt_clusters.append([item[1] for item in current])
                current = []
            current.append((rt, row))
        if current:
            rt_clusters.append([item[1] for item in current])

    for cluster_index, rows in enumerate(rt_clusters, start=1):
        if len(rows) < SAMPLE_CONSENSUS_MIN_SPECTRA:
            continue

        rts = [_float_or_none(r.get("retention_time")) for r in rows]
        rts = [x for x in rts if x is not None]
        if len(rts) < 2:
            continue
        # Evidence dimensions must co-occur in each supporting spectrum.  The
        # old max(cosine)/max(final)/min(ppm)/max(class) aggregate could take
        # four different scans, none of which individually passed the gate,
        # then promote every Medium row in the cluster.
        qualified_rows = []
        for row in rows:
            m = row.get("match") or {}
            cos = _float_or_none(m.get("cosine"))
            final = _float_or_none(m.get("final_score"))
            ppm = _float_or_none(m.get("mass_ppm"))
            class_agreement = _float_or_none(m.get("class_agreement"))
            if (
                not m.get("isomer_ambiguous")
                and m.get("mass_metadata_consistent") is not False
                and cos is not None and cos >= SAMPLE_CONSENSUS_MIN_COS
                and final is not None and final >= SAMPLE_CONSENSUS_MIN_FINAL
                and (ppm is None or abs(ppm) <= SAMPLE_CONSENSUS_MAX_PPM)
                and (class_agreement or 0.0) >= SAMPLE_CONSENSUS_MIN_CLASS_AGREE
            ):
                qualified_rows.append(row)
        if len(qualified_rows) < SAMPLE_CONSENSUS_MIN_SPECTRA:
            continue

        matches = [r["match"] for r in qualified_rows]
        qualified_rts = [_float_or_none(r.get("retention_time")) for r in qualified_rows]
        qualified_rts = [x for x in qualified_rts if x is not None]
        cos_vals = [_float_or_none(m.get("cosine")) for m in matches]
        fin_vals = [_float_or_none(m.get("final_score")) for m in matches]
        cos_vals = [x for x in cos_vals if x is not None]
        fin_vals = [x for x in fin_vals if x is not None]
        if not cos_vals or not fin_vals:
            continue

        class_agree = min((_float_or_none(m.get("class_agreement")) or 0.0) for m in matches)

        collision_energies = sorted({
            round(value, 2)
            for value in (_float_or_none(r.get("collision_energy")) for r in qualified_rows)
            if value is not None
        })
        support = {
            "cluster_id": cluster_index,
            "n_spectra": len(qualified_rows),
            "cluster_spectra": len(rows),
            "rt_min": round(min(qualified_rts), 4),
            "rt_max": round(max(qualified_rts), 4),
            "rt_span_min": round(max(qualified_rts) - min(qualified_rts), 4),
            "best_cosine": round(max(cos_vals), 4),
            "median_cosine": round(_median(cos_vals) or 0.0, 4),
            "best_final": round(max(fin_vals), 4),
            "class_agreement": round(class_agree, 3),
            "collision_energies": collision_energies,
            "n_collision_energies": len(collision_energies),
            "multi_energy": len(collision_energies) >= 2,
        }

        for m in matches:
            m["sample_consensus_support"] = support
            m["multi_energy_consensus"] = support["multi_energy"]
            if (m.get("confidence") == "Medium"
                    and m.get("evidence_tier") not in ("High", "RT-confirmed High", "Strong Candidate")):
                m["evidence_tier"] = "Strong Candidate"
                m["annotation_level"] = "Level 2 candidate (sample-consensus supported)"
                m["sample_consensus_promoted"] = True
                energy_text = (
                    f"、覆盖{support['n_collision_energies']}个碰撞能量"
                    if support["multi_energy"] else ""
                )
                suffix = (
                    f" · 样本内{support['n_spectra']}张MS2在"
                    f"{support['rt_span_min']:.2f} min内一致支持{energy_text}"
                )
                base = str(m.get("honest_headline") or m.get("compound_name") or "")
                if "样本内" not in base:
                    m["honest_headline"] = (base + suffix).strip()

    return results


def _serialize_results(results: list[dict]) -> list[dict]:
    """Convert results to JSON-serializable format, including matched fragment data."""
    out = []
    for r in results:
        m = r.get("match")
        fragmentation = _fragmentation_annotation(r)
        background = _background_evidence(r)
        # V7 invariant: molecular networking must retain the same rich MS/MS
        # evidence for matched and unmatched rows. V6 stored raw_peaks only in
        # the matched branch, so the unknown spectra most in need of networking
        # were silently reduced to the Top-N search view.
        _raw = r.get("raw_peaks") or []
        _raw_public = []
        for peak in _raw:
            try:
                mz_value = float(peak[0])
                intensity_value = float(peak[1])
            except (TypeError, ValueError, IndexError, KeyError):
                continue
            if not (
                math.isfinite(mz_value) and math.isfinite(intensity_value)
            ):
                continue
            _raw_public.append([mz_value, intensity_value])
        serialized = {
            "spectrum_id": str(r.get("spectrum_id", "")),
            "precursor_mz": r["precursor_mz"],
            "search_precursor_mz": r.get(
                "search_precursor_mz", r.get("precursor_mz")
            ),
            "search_precursor_mz_source": str(
                r.get("search_precursor_mz_source") or ""
            ),
            "search_precursor_delta_ppm": r.get(
                "search_precursor_delta_ppm"
            ),
            "lcms_feature_used_for_search": bool(
                r.get("lcms_feature_used_for_search", False)
            ),
            "precursor_mass_audit": r.get("precursor_mass_audit") or {},
            "peaks_count": r.get("peaks_count", 0),
            "processed_peaks_count": len(r.get("query_peaks") or []),
            "effective_fragment_count": fragmentation[
                "effective_fragment_count"
            ],
            "fragmentation_status": fragmentation["status"],
            "fragmentation_annotation": fragmentation,
            "background_evidence": background,
            "retention_time": r.get("retention_time"),
            "retention_time_unit": r.get("retention_time_unit"),
            "ion_mode": r.get("ion_mode"),
            "precursor_charge": r.get("precursor_charge"),
            # Preserve only declared query-adduct provenance.  EvoGraph-F uses
            # this field when reliable and never guesses an adduct from polarity.
            "query_adduct": (
                r.get("query_adduct")
                or (m.get("query_adduct") if isinstance(m, dict) else None)
            ),
            "query_adduct_reliable": bool(
                r.get("query_adduct_reliable")
                or (
                    m.get("query_adduct_reliable")
                    if isinstance(m, dict)
                    else False
                )
            ),
            "precursor_intensity": r.get("precursor_intensity"),
            "precursor_ms1_intensity": r.get("precursor_ms1_intensity"),
            "collision_energy": r.get("collision_energy"),
            "collision_energy_unit": str(r.get("collision_energy_unit") or ""),
            **query_acquisition_projection(r),
            "collision_energy_mode": str(
                r.get("collision_energy_mode") or ""
            ),
            "instrument_type": str(r.get("instrument_type") or ""),
            "instrument": str(r.get("instrument") or ""),
            "activation_method": str(
                r.get("activation_method") or r.get("activation") or ""
            ),
            "instrument_domain_id": str(
                r.get("instrument_domain_id") or ""
            ),
            "collision_cell_region": str(
                r.get("collision_cell_region") or ""
            ),
            "rt_method_id": str(r.get("rt_method_id") or ""),
            "precursor_scan_ref": str(r.get("precursor_scan_ref") or ""),
            "isotope_evidence_status": str(r.get("isotope_evidence_status") or "not_available"),
            "precursor_isotope_evidence": r.get("precursor_isotope_evidence"),
            "precursor_isolation_evidence": r.get("precursor_isolation_evidence"),
            # TIC 概览 fallback 用 (mgf 无 MS1 连续色谱时, 前端用 RT+TIC 画 stem 图):
            # 此前漏序列化 → results.json 无 tic/base_peak → fallback 画不出 TIC。
            "tic": float(r.get("tic") or 0.0),
            "base_peak": float(r.get("base_peak") or 0.0),
            # v3 Phase 1: 空白扣除信息 (纯增量; 不扣除时为 False/0)
            "blank_used": bool(r.get("blank_used", False)),
            "blank_removed": int(r.get("blank_removed", 0)),
            "blank_status": str(r.get("blank_status") or "not_requested"),
            "blank_feature_matched": bool(r.get("blank_feature_matched", False)),
            "blank_excluded": bool(r.get("blank_excluded", False)),
            "blank_flagged": bool(r.get("blank_flagged", False)),
            "blank_clean_view_excluded": bool(
                r.get("blank_clean_view_excluded", False)
            ),
            "blank_legacy_below_ratio": bool(
                r.get("blank_legacy_below_ratio", False)
            ),
            "blank_carryover_suspected": bool(
                r.get("blank_carryover_suspected", False)
            ),
            "blank_filter_policy": str(
                r.get("blank_filter_policy") or DEFAULT_BLANK_FILTER_POLICY
            ),
            "blank_evidence_role": str(
                r.get("blank_evidence_role") or "unspecified"
            ),
            "blank_evidence_roles": list(r.get("blank_evidence_roles") or []),
            "blank_evidence_category": str(
                r.get("blank_evidence_category") or ""
            ),
            "blank_role_evidence": r.get("blank_role_evidence") or {},
            "blank_evidence_source": str(r.get("blank_evidence_source") or ""),
            "blank_sample_intensity": r.get("blank_sample_intensity"),
            "blank_intensity": r.get("blank_intensity"),
            "blank_sample_ratio": r.get("blank_sample_ratio"),
            "blank_mz_error_ppm": r.get("blank_mz_error_ppm"),
            "blank_rt_delta_minutes": r.get("blank_rt_delta_minutes"),
            "blank_quantitation_mode": str(r.get("blank_quantitation_mode") or ""),
            "blank_sample_xic_area": r.get("blank_sample_xic_area"),
            "blank_xic_area": r.get("blank_xic_area"),
            "blank_sample_xic_scans": int(r.get("blank_sample_xic_scans") or 0),
            "blank_xic_scans": int(r.get("blank_xic_scans") or 0),
            "blank_replicates_total": int(r.get("blank_replicates_total") or 0),
            "blank_replicates_applicable": int(
                r.get("blank_replicates_applicable") or 0
            ),
            "blank_replicates_detected": int(
                r.get("blank_replicates_detected") or 0
            ),
            "blank_detection_fraction": r.get("blank_detection_fraction"),
            "blank_required_detection_count": int(
                r.get("blank_required_detection_count") or 0
            ),
            "input_format": r.get("input_format"),
            "spectrum_representation": r.get("spectrum_representation"),
            "spectrum_representation_source": r.get(
                "spectrum_representation_source"
            ),
            "limited_fragment_list": bool(r.get("limited_fragment_list", False)),
            "input_evidence_note": str(r.get("input_evidence_note") or ""),
            "peak_processing_audit": r.get("peak_processing_audit"),
            "raw_peaks": _raw_public,
            "raw_peaks_input_count": len(_raw),
            "raw_peaks_saved_count": len(_raw_public),
            "raw_peaks_complete": len(_raw_public) == len(_raw),
            "raw_peaks_storage_policy": "all_valid_numeric_pairs",
            # Continuous-MS1 feature evidence.  These fields are the primary
            # grouping key for V7 networking; precursor_feature_id below is an
            # explicit MS2-event fallback for formats without survey MS1.
            "lcms_feature_id": r.get("lcms_feature_id"),
            "lcms_feature_mz": r.get("lcms_feature_mz"),
            "lcms_feature_rt_minutes": r.get("lcms_feature_rt_minutes"),
            "lcms_feature_rt_start_minutes": r.get(
                "lcms_feature_rt_start_minutes"
            ),
            "lcms_feature_rt_end_minutes": r.get(
                "lcms_feature_rt_end_minutes"
            ),
            "lcms_feature_intensity": r.get("lcms_feature_intensity"),
            "lcms_consensus_spectrum_id": r.get(
                "lcms_consensus_spectrum_id"
            ),
            "lcms_consensus_member_count": int(
                r.get("lcms_consensus_member_count") or 0
            ),
            "lcms_consensus_peak_count": int(
                r.get("lcms_consensus_peak_count") or 0
            ),
            "lcms_consensus_is_representative": bool(
                r.get("lcms_consensus_is_representative", False)
            ),
            "lcms_consensus_peaks": r.get("lcms_consensus_peaks") or [],
            "lcms_feature_charge": r.get("lcms_feature_charge"),
            "lcms_isotope_trace_count": r.get("lcms_isotope_trace_count"),
            "lcms_adduct": r.get("lcms_adduct"),
            "lcms_neutral_mass": r.get("lcms_neutral_mass"),
            "lcms_adduct_group_id": r.get("lcms_adduct_group_id"),
            "lcms_feature_link_status": r.get("lcms_feature_link_status"),
            "lcms_feature_candidate_count": r.get(
                "lcms_feature_candidate_count"
            ),
            "lcms_feature_mz_error_ppm": r.get(
                "lcms_feature_mz_error_ppm"
            ),
            "lcms_feature_rt_error_seconds": r.get(
                "lcms_feature_rt_error_seconds"
            ),
            # Preserve feature/IIMN evidence when cached rows are reserialized.
            "precursor_feature_id": r.get("precursor_feature_id"),
            "precursor_feature_spectra": r.get("precursor_feature_spectra"),
            "ion_identity_status": r.get("ion_identity_status"),
            "ion_identity_network_id": r.get("ion_identity_network_id"),
            "ion_identity": r.get("ion_identity"),
            "ion_identity_neutral_mass": r.get("ion_identity_neutral_mass"),
            "ion_identity_representative_feature_id": r.get(
                "ion_identity_representative_feature_id"
            ),
            "ion_identity_is_representative": r.get(
                "ion_identity_is_representative"
            ),
            "ion_identity_evidence": r.get("ion_identity_evidence"),
            "ion_identity_linked_annotation": r.get(
                "ion_identity_linked_annotation"
            ),
        }
        # 历史无结构命中仍降级为未匹配；只有检索层明确通过严格门控的
        # name_only_lead 才允许作为“谱图线索”进入显示层。
        if m:
            ce_annotation = _ce_comparability(r, m)
            _smi = m.get("smiles")
            if (not _smi or _is_nan(_smi)
                    or str(_smi).strip().lower() in ('', 'nan', 'none', 'null')):
                if not m.get("name_only_lead"):
                    m = None
            else:
                _cn = m.get("compound_name")
                if (_cn is None or _is_nan(_cn)
                        or str(_cn).strip().lower() in ('', 'nan', 'none', 'null')):
                    m["compound_name"] = (str(m.get("inchikey") or "").split("-")[0]
                                          or f"Compound_{r.get('spectrum_id', '')}")
        if m:
            # Compute matched fragments
            qp = r.get("query_peaks", [])
            lp = m.get("lib_peaks", [])
            msms_str = _msms_fragments(
                qp,
                lp,
                query_pmz=(
                    r.get("search_precursor_mz") or r.get("precursor_mz")
                ),
                lib_pmz=m.get("lib_precursor_mz"),
            )
            matched_count = len(msms_str.split("; ")) if msms_str else 0
            # Top query peaks summary
            top_qp = sorted(qp, key=lambda p: -p[1])[:5] if qp else []
            top_qp_str = "; ".join(f"{p[0]:.4f}" for p in top_qp) if top_qp else ""

            # 无结构(匿名谱图/无SMILES)命中降级: 库内15351条无SMILES记录(如匿名Spectraverse谱图)
            # 用spectrum_id冒充化合物名混入命中且cos偏高, 但无法鉴定(不知是何物) → 强制Unknown+标记。
            # 依据: NIST/MoNA 库质量过滤 best practice — 无化合物身份的记录不应作为高置信鉴定。
            _has_struct = (
                m.get("has_structure") is not False
                and bool(m.get("smiles"))
                and str(m.get("smiles")).strip().lower() not in ("", "nan", "none", "null")
            )
            # P4 诚实性校准 (单一真源, 取代原 scores.confidence + FP 降级 + report.py local label):
            # 收编 无结构→Unknown + 假阳性(ppm/reverse_match)降级, 新增 异构不可分/弱匹配诚实降级 + 类级置信。
            # 零排序风险: 只算呈现字段, 不改 final_score → 红线 11/14 不受影响。
            m["query_search_precursor_mz_source"] = str(
                r.get("search_precursor_mz_source")
                or m.get("query_search_precursor_mz_source")
                or ""
            )
            m["query_adduct"] = r.get("query_adduct") or m.get("query_adduct")
            m["query_adduct_reliable"] = bool(
                r.get("query_adduct_reliable") or m.get("query_adduct_reliable")
            )
            m["isotope_evidence_status"] = str(
                r.get("isotope_evidence_status")
                or m.get("isotope_evidence_status")
                or "not_available"
            )
            _ann = _honest_annotate(m, m.get("top_candidates"))
            from web.report_common import classification_evidence_status

            _class_evidence = classification_evidence_status(m)
            _conf = _ann["compound_confidence"]
            # top-5 HTML 摘要即使与某张库谱高度相似，也缺少确认位置异构体所需的
            # 低丰度诊断峰。不得把单一库名显示为化合物级 High；最多 Strong Candidate /
            # Medium，并明确位置异构体未确认。类级置信可保留原校准结果。
            if r.get("limited_fragment_list") and _conf == "High":
                _conf = "Medium"
                _ann.update({
                    "compound_confidence": "Medium",
                    "exact_structure_confidence": "Medium",
                    "honesty_tier": "limited_ms2",
                    "evidence_tier": "Strong Candidate",
                    "annotation_level": (
                        "Strong spectral candidate — limited top-5 MS2; "
                        "positional isomer unconfirmed"
                    ),
                    "honest_headline": (
                        f"{m.get('display_class') or m.get('classification') or 'Class-level'} "
                        "candidate（HTML 仅 top-5 MS2，位置异构体未确认）"
                    ),
                    "rt_confirmed": False,
                })
            if (
                r.get("limited_fragment_list")
                and _ann.get("isomer_group_confidence") == "High"
            ):
                _group_evidence = dict(
                    _ann.get("isomer_group_evidence") or {}
                )
                _group_evidence["confidence"] = "Medium"
                _ann.update({
                    "isomer_group_confidence": "Medium",
                    "isomer_group_evidence": _group_evidence,
                    "evidence_tier": "Strong Candidate",
                    "annotation_level": (
                        "Limited top-5 MS2; isomer-group High not assigned"
                    ),
                })
            _name_only = bool(m.get("name_only_lead") and not _has_struct)
            if _name_only:
                _conf = "Unknown"
                _ann.update({
                    "compound_confidence": "Unknown",
                    "class_confidence": "Unknown",
                    "exact_structure_confidence": "Unknown",
                    "isomer_group_confidence": "Unknown",
                    "identification_scope": "spectral_lead",
                    "isomer_group_formula": "",
                    "isomer_group_size": 0,
                    "isomer_group_evidence": {},
                    "isomer_literature_evidence": {},
                    "honesty_tier": "name_only_lead",
                    "evidence_tier": "Spectral Lead",
                    "annotation_level": "Library spectrum lead — no structure",
                    "honest_headline": f"{m.get('compound_name', '')}（仅谱图线索，无结构）",
                    "rt_confirmed": False,
                })
            serialized["match"] = {
                **{key: m[key] for key in ("query_adduct_charge_evidence",) if key in m},
                "cosine": round(m.get("cosine", 0), 4),
                "entropy": round(m.get("entropy", 0), 4),
                "mass_score": round(m.get("mass_score", 0), 4),
                # 0 ppm is valid evidence, not a missing value.
                "mass_ppm": round(m["mass_ppm"], 2) if m.get("mass_ppm") is not None else None,
                "mass_match_basis": str(m.get("mass_match_basis", "")),
                "mass_metadata_consistent": m.get("mass_metadata_consistent"),
                "mass_query_value": m.get("mass_query_value"),
                "mass_library_value": m.get("mass_library_value"),
                "query_selected_ion_mz": m.get("query_selected_ion_mz"),
                "query_search_precursor_mz": m.get(
                    "query_search_precursor_mz"
                ),
                "query_search_precursor_mz_source": str(
                    m.get("query_search_precursor_mz_source") or ""
                ),
                "query_search_precursor_delta_ppm": m.get(
                    "query_search_precursor_delta_ppm"
                ),
                "query_adduct": m.get("query_adduct"),
                "query_adduct_reliable": bool(
                    m.get("query_adduct_reliable")
                ),
                "isotope_evidence_status": str(
                    m.get("isotope_evidence_status") or "not_available"
                ),
                "precursor_isotope_evidence": m.get(
                    "precursor_isotope_evidence"
                ),
                "query_precursor_mass_audit": m.get(
                    "query_precursor_mass_audit"
                ) or {},
                "lib_neutral_mass": m.get("lib_neutral_mass"),
                "lib_neutral_mass_source": str(m.get("lib_neutral_mass_source", "")),
                "fragment_score": round(m.get("fragment_score", 0), 4),
                "intensity_score": round(m.get("intensity_score", 0), 4),
                "isomer_fragment_score": round(m.get("isomer_fragment_score", 0) or 0, 4),
                "isomer_fragment_bonus": round(m.get("isomer_fragment_bonus", 0) or 0, 4),
                "isomer_fragment_detail": m.get("isomer_fragment_detail"),
                "multimetric_shadow": m.get("multimetric_shadow"),
                "isopair_probability": m.get("isopair_probability"),
                "isopair_fused_score": m.get("isopair_fused_score"),
                "isopair_rank_applied": bool(m.get("isopair_rank_applied")),
                "isopair_selected_connectivity": bool(
                    m.get("isopair_selected_connectivity")
                ),
                "isopair_high_precision_support": _ann.get(
                    "isopair_high_precision_support", False
                ),
                "isopair_high_precision_gate_required": _ann.get(
                    "isopair_high_precision_gate_required", False
                ),
                "isopair_high_precision_gate_result": _ann.get(
                    "isopair_high_precision_gate_result", "not_applicable"
                ),
                "isopair_reranked_2d_connectivity": _ann.get(
                    "isopair_reranked_2d_connectivity", False
                ),
                "isopair_high_eligible": _ann.get(
                    "isopair_high_eligible", False
                ),
                "isopair_evidence": m.get("isopair_evidence") or {},
                "precursor_high_evidence_gate_required": _ann.get(
                    "precursor_high_evidence_gate_required", False
                ),
                "precursor_high_evidence_gate_result": _ann.get(
                    "precursor_high_evidence_gate_result", "not_applicable"
                ),
                "precursor_high_evidence_gate_reason": _ann.get(
                    "precursor_high_evidence_gate_reason"
                ),
                "precursor_high_evidence_gate_reason_zh": _ann.get(
                    "precursor_high_evidence_gate_reason_zh"
                ),
                "direct_matched_product_ions": _ann.get(
                    "direct_matched_product_ions"
                ),
                "final_score": round(m.get("final_score", 0), 4),
                "rank_score": round(m.get("rank_score", m.get("final_score", 0)) or 0, 4),
                "raw_cosine": (
                    round(m["raw_cosine"], 4) if m.get("raw_cosine") is not None else None
                ),
                "raw_entropy": (
                    round(m["raw_entropy"], 4) if m.get("raw_entropy") is not None else None
                ),
                "raw_fragment_score": (
                    round(m["raw_fragment_score"], 4)
                    if m.get("raw_fragment_score") is not None else None
                ),
                "raw_intensity_score": (
                    round(m["raw_intensity_score"], 4)
                    if m.get("raw_intensity_score") is not None else None
                ),
                "raw_reverse_match": (
                    round(m["raw_reverse_match"], 4)
                    if m.get("raw_reverse_match") is not None else None
                ),
                "raw_final_score": (
                    round(m["raw_final_score"], 4)
                    if m.get("raw_final_score") is not None else None
                ),
                "raw_rank_score": (
                    round(m["raw_rank_score"], 4)
                    if m.get("raw_rank_score") is not None else None
                ),
                "chimeric_deconvolution": m.get("chimeric_deconvolution"),
                "search_pass": m.get("search_pass", "primary_top15"),
                "primary_search": m.get("primary_search"),
                "compound_name": str(m.get("compound_name", ""))[:200],
                "compound_name_cn": str(m.get("compound_name_cn", ""))[:100],
                "source_compound_name": str(m.get("source_compound_name", ""))[:500],
                "name_resolution_status": str(m.get("name_resolution_status", "")),
                "smiles": str(m.get("smiles", "")),
                "inchikey": str(m.get("inchikey", "")),
                "reported_inchikey": str(m.get("reported_inchikey", "")),
                "formula": _resolve_formula(m.get("smiles", ""), m.get("formula", "")),
                "adduct": str(m.get("adduct", "")),
                "precursor_metadata_audit": m.get("precursor_metadata_audit") or {},
                "mass_metadata_consistent": m.get("mass_metadata_consistent"),
                "candidate_adduct_hypothesis": str(
                    m.get("candidate_adduct_hypothesis") or ""
                ),
                "candidate_adduct_hypothesis_source": str(
                    m.get("candidate_adduct_hypothesis_source") or ""
                ),
                "lib_retention_time": m.get("lib_retention_time"),
                "rt_evidence": m.get("rt_evidence"),
                "lib_charge": m.get("lib_charge"),
                "lib_collision_energy": m.get("lib_collision_energy"),
                "lib_collision_energy_raw": m.get("lib_collision_energy_raw"),
                "lib_collision_energy_unit": m.get("lib_collision_energy_unit"),
                "lib_collision_energy_mode": m.get("lib_collision_energy_mode"),
                "lib_instrument_type": m.get("lib_instrument_type"),
                "lib_instrument": m.get("lib_instrument"),
                "lib_activation_method": m.get("lib_activation_method"),
                "lib_instrument_domain_id": m.get(
                    "lib_instrument_domain_id"
                ),
                "lib_collision_cell_region": m.get(
                    "lib_collision_cell_region"
                ),
                "query_instrument_type": m.get(
                    "query_instrument_type", ""
                ),
                "query_instrument": m.get("query_instrument", ""),
                "query_rt_method_id": m.get(
                    "query_rt_method_id", ""
                ),
                "lib_source_accession": m.get("lib_source_accession"),
                "ce_comparability": ce_annotation,
                "ce_profile_evidence": m.get("ce_profile_evidence"),
                "feature_fragment_evidence": m.get(
                    "feature_fragment_evidence"
                ),
                "formula_isotope_evidence": m.get("formula_isotope_evidence"),
                # 修复A: 此前序列化丢失 → 564/564 命中 lib_mz=None, 用户无法核对 PPM 来源
                "lib_precursor_mz": m.get("lib_precursor_mz"),
                "lib_spectrum_id": m.get("lib_spectrum_id"),
                # 修复B: 无结构(匿名谱图)标记, 供前端/报告降级展示(否则匿名谱图冒充化合物)
                "has_structure": _has_struct,
                "structure_quality": str(m.get("structure_quality", "")),
                "name_only_lead": _name_only,
                "reverse_match": round(m.get("reverse_match") or 0, 3),
                "classification": str(_class_evidence["label"]),
                "display_class": str(m.get("display_class") or ""),
                "raw_classification": str(m.get("classification", "")),
                "classification_source": str(_class_evidence["source"]),
                "display_classification_source": str(
                    m.get("display_classification_source") or ""
                ),
                "classification_reliable": bool(_class_evidence["reliable"]),
                "classification_status": str(_class_evidence["status"]),
                "classification_reliability_reason": str(
                    _class_evidence["reason"]
                ),
                "classification_confidence": m.get(
                    "display_classification_confidence",
                    m.get("classification_confidence"),
                ),
                "raw_classification_source": str(
                    m.get("classification_source") or "unknown"
                ),
                "raw_classification_confidence": m.get(
                    "classification_confidence"
                ),
                "dataset_source": str(m.get("dataset_source", "")),
                "source_tier": str(m.get("source_tier", "unknown")),
                "intensity_origin": str(m.get("intensity_origin", "unknown")),
                "quality_flags": m.get("quality_flags") or [],
                "compound_role": str(m.get("compound_role", "")),
                "compound_role_note": str(m.get("compound_role_note", "")),
                "potential_contaminant": bool(m.get("potential_contaminant")),
                "contaminant_category": str(m.get("contaminant_category", "")),
                "contaminant_note": str(m.get("contaminant_note", "")),
                "exclude_from_natural_product_count": bool(
                    m.get("exclude_from_natural_product_count")
                ),
                "natural_product_candidate": bool(
                    m.get("natural_product_candidate", False)
                ),
                "dl_classification": str(m.get("dl_classification", "") or ""),
                "dl_confidence": round(m.get("dl_confidence", 0) or 0, 3),
                "confidence": _conf,
                # P4 诚实性校准字段 (供前端/报告诚实呈现: 类级置信/异构集/诚实标题)
                "confidence_class": _ann["class_confidence"],
                "exact_structure_confidence": _ann.get(
                    "exact_structure_confidence", _conf
                ),
                "two_dimensional_connectivity_confidence": _ann.get(
                    "two_dimensional_connectivity_confidence", _conf
                ),
                "isomer_group_confidence": _ann.get(
                    "isomer_group_confidence", "Unknown"
                ),
                "identification_scope": _ann.get(
                    "identification_scope", "putative_candidate"
                ),
                "isomer_group_formula": _ann.get("isomer_group_formula", ""),
                "isomer_group_size": _ann.get("isomer_group_size", 0),
                "isomer_group_evidence": _ann.get(
                    "isomer_group_evidence", {}
                ),
                "honesty_tier": _ann["honesty_tier"],
                "isomer_ambiguous": _ann["isomer_ambiguous"],
                "isomer_n": _ann["isomer_n"],
                "isomer_set": _ann["isomer_set"],
                "same_formula_alternatives": _ann.get("same_formula_alternatives") or [],
                "same_formula_alternatives_n": _ann.get("same_formula_alternatives_n", 0),
                "noncompetitive_same_formula_n": _ann.get("noncompetitive_same_formula_n", 0),
                "isomer_ambiguity_type": _ann.get(
                    "isomer_ambiguity_type", "none"
                ),
                "isomer_relationship_counts": _ann.get(
                    "isomer_relationship_counts", {}
                ),
                "isomer_variant_subtype_counts": _ann.get(
                    "isomer_variant_subtype_counts", {}
                ),
                "isomer_resolution_status": _ann.get(
                    "isomer_resolution_status", "not_applicable"
                ),
                "isomer_literature_evidence": _ann.get(
                    "isomer_literature_evidence", {}
                ),
                "class_consensus": _ann["class_consensus"],
                "class_agreement": _ann["class_agreement"],
                "honest_headline": _ann["honest_headline"],
                "evidence_tier": _ann.get("evidence_tier"),
                "annotation_level": _ann.get("annotation_level"),
                "candidate_comparison_source": _ann.get(
                    "candidate_comparison_source"
                ),
                "candidate_comparison_cosine": _ann.get(
                    "candidate_comparison_cosine"
                ),
                "candidate_margin": _ann.get("candidate_margin"),
                "deconvolution_score_influence": _ann.get(
                    "deconvolution_score_influence", False
                ),
                "deconvolution_evidence_state": _ann.get(
                    "deconvolution_evidence_state", "not_applicable"
                ),
                "chimeric_isolation_warning": _ann.get(
                    "chimeric_isolation_warning"
                ),
                "rt_confirmed": _ann.get("rt_confirmed"),
                "rt_warning": _ann.get("rt_warning"),
                # Fragment ion data
                "ms_ms_fragments": msms_str,
                "fragment_count": matched_count,
                "top_query_peaks": top_qp_str,
                # ── Expert system data (v3.0) ──
                "expert_score": round(m.get("expert_score", 0), 4),
                "expert_details": m.get("expert_details"),
                "expert_score_source": m.get("expert_score_source", "none"),
                "expert_score_components": m.get(
                    "expert_score_components", {}
                ),
                "contrastive_expert_shadow": m.get("contrastive_expert_shadow"),
                # P2-show: 预存馆藏裂解途径 JSON (直接读自 DB, 报告展示用; 评分已在 search 用过)
                "pathway_json": m.get("pathway_json"),
                "pathway_polarity_status": m.get("pathway_polarity_status"),
                "pathway_polarities": m.get("pathway_polarities") or [],
                "pathway_adduct_status": m.get("pathway_adduct_status"),
                "pathway_adducts": m.get("pathway_adducts") or [],
                # v3 Phase 5: Top-K 候选 (近邻/异构体, best 在首位)
                "top_candidates": m.get("top_candidates") or [],
                "spectral_evidence": m.get("spectral_evidence"),
            }
            # Include peaks for spectrum visualization
            serialized["query_peaks"] = [
                [float(p[0]), float(p[1])] for p in qp
            ]
            serialized["lib_peaks"] = [
                [float(p[0]), float(p[1])] for p in lp
            ]
            serialized["deconvolved_query_peaks"] = [
                [float(p[0]), float(p[1])]
                for p in (m.get("deconvolved_query_peaks") or [])
            ]
            # 供详情/导出展示的原始 MS/MS 峰已在公共行字段中完整保存；
            # 报告层不得再做 Top-N 截断。
            serialized["raw_peaks"] = _raw_public
        else:
            serialized["match"] = None
            serialized["query_peaks"] = [
                [float(p[0]), float(p[1])] for p in r.get("query_peaks", [])
            ]
            serialized["lib_peaks"] = []
            # 被空白特征过滤的谱图不是“未命中”，不要再用类别模型把它包装成样品线索。
            if serialized["blank_excluded"]:
                out.append(serialized)
                continue
            # 未命中谱图用 MSClassify 从谱图判类 (补无化合物命中的判类空白)
            try:
                from classifier.msclass_wrapper import get_msclass_wrapper
                _msc_w = get_msclass_wrapper()
                if _msc_w is not None:
                    _msc_r = _msc_w.predict(r.get("query_peaks", []))
                    # 0.15 仅是“保留模型输出”的下限，不是可靠性阈值。softmax 最大概率
                    # 未经当前样本域校准；前端按 >=0.70 强 / 0.45~0.70 中 / <0.45 低展示。
                    if _msc_r and _msc_r["confidence"] >= 0.15:
                        serialized["display_class"] = _msc_r["class"]
                        serialized["msclass_source"] = "unmatched"
                        serialized["dl_confidence"] = round(_msc_r["confidence"], 3)
                        serialized["dl_top3"] = _msc_r["top3"]
            except Exception as _e:
                try:
                    logger.debug("MSClassify unmatched predict failed: %s", _e)
                except Exception:
                    pass
        out.append(serialized)
    return _prepare_results_for_display(out)


def _save_task_results(task_id: str):
    """Save task results as JSON file (reads from _tasks dict)."""
    task = _tasks.get(task_id)
    if not task:
        return
    result_path = TASK_DIR / f"{task_id}_results.json"
    try:
        with open(result_path, "w", encoding="utf-8") as f:
            write_json_encoded(
                task["results"],
                f,
                ensure_ascii=False,
                separators=(",", ":"),
                default=str,
            )
    except Exception as e:
        logger.error("Failed to save results for %s: %s", task_id, e)


def _save_task_results_given(task_id: str, results: list[dict]):
    """Save provided results list to disk for a task."""
    result_path = TASK_DIR / f"{task_id}_results.json"
    try:
        with open(result_path, "w", encoding="utf-8") as f:
            write_json_encoded(
                results,
                f,
                ensure_ascii=False,
                separators=(",", ":"),
                default=str,
            )
    except Exception as e:
        logger.error("Failed to save results for %s: %s", task_id, e)


def _expected_compounds_artifact_path(task_id: str) -> Path:
    return TASK_DIR / f"{task_id}_expected_compounds.json"


def _load_expected_compounds_payload(
    task_id: str,
    task: Optional[dict[str, Any]] = None,
    *,
    tolerate_errors: bool = False,
) -> dict[str, Any]:
    """Load the independent Expected Compounds artifact for UI or export."""

    artifact_path = _expected_compounds_artifact_path(task_id)
    if artifact_path.exists():
        try:
            payload = json.loads(artifact_path.read_text(encoding="utf-8"))
            if isinstance(payload, dict):
                return payload
            raise ValueError("artifact root is not an object")
        except (OSError, ValueError) as exc:
            if not tolerate_errors:
                raise
            return _expected_compounds_unavailable_payload(
                status="failed",
                input_path=str((task or {}).get("file_path") or ""),
                requested=bool(
                    ((task or {}).get("expected_compounds_summary") or {}).get(
                        "requested", True
                    )
                ),
                warning=f"Expected Compounds artifact is unreadable: {exc}",
            )

    task = task or {}
    summary = dict(task.get("expected_compounds_summary") or {})
    return {
        "schema_version": 1,
        "task_id": task_id,
        "method": "Expected Compounds / continuous-MS1 targeted XIC",
        "status": str(summary.get("status") or "not_run"),
        "requested": bool(summary.get("requested", False)),
        "summary": {
            **summary,
            "counted_in_standard_spectral_library_matches": False,
            "standard_spectral_library_match_increment": 0,
        },
        "features": [],
        "scientific_boundaries": {
            "presence_only_rows_are_not_reference_spectra": True,
            "exact_mass_is_not_isomer_specific": True,
            "ordinary_library_matched_count_unchanged": True,
        },
    }


def _load_reference_standard_qc_payload(
    task_id: str, *, tolerate_errors: bool = False
) -> dict[str, Any]:
    artifact = _reference_standard_qc_artifact_path(TASK_DIR, task_id)
    if not artifact.is_file():
        return {}
    try:
        return _read_reference_standard_qc_artifact(artifact)
    except (OSError, ValueError, json.JSONDecodeError):
        if not tolerate_errors:
            raise
        logger.warning(
            "Reference-standard QC artifact is unreadable for task %s", task_id
        )
        return {}


def _annotate_expected_compound_background(
    payload: dict[str, Any],
    *,
    file_path: str,
    params: dict[str, Any],
    blank_profiles: list[Any] | None,
    blank_roles: list[str] | None,
) -> dict[str, Any]:
    """Attach role-aware blank evidence without deleting targeted features."""

    features = [
        feature
        for feature in payload.get("features") or []
        if isinstance(feature, dict)
    ]
    profiles = list(blank_profiles or [])
    if not profiles:
        for feature in features:
            feature.setdefault("blank_status", "not_requested")
            feature.setdefault("blank_feature_matched", False)
            feature.setdefault("blank_flagged", False)
        return {
            "status": "not_requested",
            "targeted_feature_count": len(features),
            "background_flagged_features": 0,
            "identity_confidence_increment": 0,
        }

    from matching.blank_subtraction import assess_blank_role_groups

    proxies = [{
        "spectrum_id": feature.get("feature_id"),
        "precursor_mz": (
            feature.get("observed_apex_mz") or feature.get("target_mz")
        ),
        "retention_time": feature.get("apex_rt_minutes"),
        "retention_time_unit": "minutes",
        "ion_mode": feature.get("polarity"),
        "polarity": feature.get("polarity"),
        "precursor_ms1_intensity": feature.get("apex_intensity"),
    } for feature in features]
    assessments = assess_blank_role_groups(
        proxies,
        profiles,
        list(blank_roles or []),
        sample_path=(
            file_path
            if Path(file_path).suffix.lower()
            in {".mzml", ".mzxml", ".xml"}
            else None
        ),
        sample_rt_unit="minutes",
        sample_to_blank_ratio=params.get(
            "expected_background_sample_blank_ratio", 5.0
        ),
        rt_tolerance_minutes=params.get(
            "expected_openms_rt_tolerance_minutes", 0.05
        ),
        precursor_tolerance_da=1.0,
        precursor_tolerance_ppm=params.get(
            "expected_mass_tolerance_ppm", 5.0
        ),
        min_detection_fraction=params.get(
            "blank_min_detection_fraction",
            DEFAULT_BLANK_MIN_DETECTION_FRACTION,
        ),
        # CD-style background marking is evidence annotation, not identity
        # deletion.  The optional clean view remains a presentation choice.
        filter_policy="soft_evidence",
        clean_max_ratio=params.get(
            "blank_clean_max_ratio", DEFAULT_BLANK_CLEAN_MAX_RATIO
        ),
        clean_min_replicates=params.get(
            "blank_clean_min_replicates",
            DEFAULT_BLANK_CLEAN_MIN_REPLICATES,
        ),
    )
    for feature, assessment in zip(features, assessments):
        feature.update(assessment)
    status_counts = Counter(
        str(item.get("blank_status") or "") for item in assessments
    )
    category_counts = Counter(
        str(item.get("blank_evidence_category") or "")
        for item in assessments
        if item.get("blank_evidence_category")
    )
    return {
        "status": "completed",
        "method": "role-aware paired sample/blank MS1 XIC",
        "targeted_feature_count": len(features),
        "blank_replicate_count": len(profiles),
        "sample_to_blank_ratio_threshold": float(params.get(
            "expected_background_sample_blank_ratio", 5.0
        )),
        "status_counts": dict(status_counts),
        "category_counts": dict(category_counts),
        "background_flagged_features": sum(
            bool(item.get("blank_flagged")) for item in assessments
        ),
        "clean_view_excluded_features": sum(
            bool(item.get("blank_clean_view_excluded"))
            for item in assessments
        ),
        "identity_confidence_increment": 0,
    }


def _run_expected_compounds_layer(
    task_id: str,
    *,
    file_path: str,
    detected_ion_mode: str,
    params: dict[str, Any],
    serialized_ms2_rows: list[dict],
    openms_features: list[dict[str, Any]] | None = None,
    blank_profiles: list[Any] | None = None,
    blank_roles: list[str] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Execute and persist the independent continuous-MS1 target layer."""

    artifact_path = _expected_compounds_artifact_path(task_id)
    requested = bool(params.get("enable_expected_compounds", False))
    provider_call_count = 0
    if not requested:
        payload = _expected_compounds_unavailable_payload(
            status="disabled",
            input_path=file_path,
            requested=False,
            warning="Expected Compounds was disabled for this analysis.",
        )
    else:
        provider_call_count = 1
        payload = run_expected_compounds_workflow(
            file_path,
            # TargetCompounds is a presence-only catalog stored in the
            # production database.  It is independent of whether the user
            # selected the Small or Full spectral-library search arrays.
            db_path=DB_PATH,
            polarity=detected_ion_mode,
            ms2_rows=serialized_ms2_rows,
            dataset_sources=params.get("expected_compound_datasets") or [],
            include_safe_common_adducts=params.get(
                "expected_generate_common_adducts", True
            ),
            mass_ppm=params.get("expected_mass_tolerance_ppm", 5.0),
            minimum_intensity=params.get(
                "expected_min_peak_intensity", 1000.0
            ),
            signal_to_noise=params.get("expected_min_snr", 1.5),
            minimum_scan_points=params.get("expected_min_scan_points", 5),
            maximum_gap_ratio=params.get(
                "expected_max_gap_ratio", 0.35
            ),
            maximum_fwhm_seconds=params.get(
                "expected_max_fwhm_seconds", 60.0
            ),
            minimum_relative_valley_depth=params.get(
                "expected_min_relative_valley_depth", 0.10
            ),
            ms2_precursor_ppm=params.get(
                "expected_compound_ms2_precursor_ppm", 10.0
            ),
            diagnostic_tolerance_da=params.get(
                "expected_compound_diagnostic_tolerance_da", 0.02
            ),
            openms_features=openms_features,
            openms_rt_tolerance_minutes=params.get(
                "expected_openms_rt_tolerance_minutes", 0.05
            ),
        )
        try:
            background_validation = _annotate_expected_compound_background(
                payload,
                file_path=file_path,
                params=params,
                blank_profiles=blank_profiles,
                blank_roles=blank_roles,
            )
        except Exception as exc:
            logger.exception(
                "[task %s] Expected Compounds blank annotation failed: %s",
                task_id,
                exc,
            )
            background_validation = {
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}"[:1000],
                "targeted_feature_count": len(
                    payload.get("features") or []
                ),
                "identity_confidence_increment": 0,
            }
        payload.setdefault("summary", {})[
            "background_validation"
        ] = background_validation
        payload.setdefault("scientific_boundaries", {})[
            "background_evidence_does_not_identify_compounds"
        ] = True
    target_selection = {
        "selected_dataset_sources": list(
            params.get("expected_compound_datasets") or []
        ),
        "selection_source": str(
            params.get("expected_compound_dataset_selection_source") or "none"
        ),
        "selection_status": str(
            params.get("expected_compound_dataset_selection_status") or ""
        ),
        "inferred_formula_code": str(
            params.get("expected_compound_inferred_formula_code") or ""
        ),
        "empty_selection_searches_all_catalogues": False,
    }
    payload["target_dataset_selection"] = target_selection
    payload.setdefault("catalog", {}).update(target_selection)
    payload.setdefault("summary", {})[
        "target_dataset_selection"
    ] = target_selection
    payload["task_id"] = task_id
    write_expected_compounds_artifact(artifact_path, payload)
    summary = dict(payload.get("summary") or {})
    summary.update(
        {
            "status": str(payload.get("status") or summary.get("status") or ""),
            "requested": requested,
            "artifact": str(artifact_path),
            "method": str(payload.get("method") or ""),
            "available": str(payload.get("status") or "") == "ok",
            "counted_in_standard_spectral_library_matches": False,
            "standard_spectral_library_match_increment": 0,
            "expected_provider_call_count": provider_call_count,
            "evaluation_panel_access_by_pipeline": 0,
        }
    )
    return payload, summary


# CFM 后台预计算去重: 同 task 只启动一次 (避免重复分析/重试时多次 Popen)
_cfm_precompute_started: set[str] = set()
_cfm_precompute_lock = threading.Lock()
_auto_cfm_precompute = str(
    os.environ.get("MSSYSTEM_AUTO_PRECOMPUTE_CFM", "0")
).strip().lower() in {"1", "true", "yes", "on"}


def _maybe_launch_cfm_precompute(
    task_id: str,
    *,
    manual: bool = False,
) -> dict[str, Any]:
    """分析完成后, 若 Docker 可用 → 后台非阻塞预灌 CFM 碎裂缓存。

    解决"detailed 报告用 cache_only=True 不跑 docker → 开 Docker 也看不到 CFM"。
    后台跑 precompute_cfm_fragments.py (每化合物 ~4s depth1), 缓存就绪后报告读缓存瞬时显示
    CFM 产物级联; 缓存未就绪时报告退 RDKit 断键推断。不阻塞 analyze 返回与单 worker 事件循环。
    """
    try:
        if not manual and not _auto_cfm_precompute:
            logger.info(
                "[task %s] automatic CFM precompute disabled; use the explicit "
                "task endpoint when the analysis queue is idle",
                task_id,
            )
            return {"status": "disabled_by_default"}
        with _tasks_lock:
            active = [
                key
                for key, value in _tasks.items()
                if _is_active_task_state(value.get("status"))
            ]
        if active:
            logger.info(
                "[task %s] CFM precompute deferred: %d analyses are active",
                task_id,
                len(active),
            )
            return {
                "status": "deferred_analysis_queue_active",
                "active_analysis_count": len(active),
            }
        from web.report_common import docker_available
        if not docker_available():
            logger.info("[task %s] Docker 不可用, 跳过 CFM 后台预计算 (报告将用 RDKit 断键推断)", task_id)
            return {"status": "docker_unavailable"}
        with _cfm_precompute_lock:
            if task_id in _cfm_precompute_started:
                return {"status": "already_started"}
            _cfm_precompute_started.add(task_id)
        ROOT = Path(__file__).resolve().parent.parent
        log_path = ROOT / "logs" / f"cfm_precompute_{task_id}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_f = open(log_path, "a", encoding="utf-8")
        # start_new_session=True → 完全脱离父进程, 服务器退出也不杀; 结果只写缓存无副作用
        subprocess.Popen(
            [sys.executable, str(ROOT / "precompute_cfm_fragments.py"), task_id],
            cwd=str(ROOT), stdout=log_f, stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        logger.info("[task %s] 已后台启动 CFM 预计算 → %s", task_id, log_path)
        return {"status": "started", "log_path": str(log_path)}
    except Exception as e:
        logger.warning("[task %s] 启动 CFM 后台预计算失败(不影响分析): %s", task_id, e)
        return {"status": "error", "error": str(e)}


# ═════════════════════════════════════════════════════════════════════
# API Endpoints
# ═════════════════════════════════════════════════════════════════════

@app.on_event("startup")
async def startup():
    """Load tasks from disk on startup."""
    global _vendor_conversion_accepting_claims
    newly_acquired_ledger_lock = _acquire_vendor_conversion_service_lock()
    try:
        _load_tasks()
        recovered_vendor_items = await asyncio.to_thread(
            _vendor_conversion_ledger().recover_interrupted,
            force=newly_acquired_ledger_lock,
        )
    except Exception as exc:
        _release_vendor_conversion_service_lock()
        logger.error("Vendor conversion ledger recovery failed: %s", exc)
        raise RuntimeError(
            "Vendor conversion ledger recovery failed"
        ) from exc
    else:
        if recovered_vendor_items:
            logger.warning(
                "Recovered %d interrupted vendor conversion item(s)",
                len(recovered_vendor_items),
            )
    with _vendor_conversion_active_claims_lock:
        _vendor_conversion_accepting_claims = True
    logger.info("MSsystem Server started. %d tasks loaded.", len(_tasks))


@app.on_event("shutdown")
async def shutdown():
    """Release the cross-process conversion-ledger ownership lock."""

    global _vendor_conversion_accepting_claims
    with _vendor_conversion_active_claims_lock:
        _vendor_conversion_accepting_claims = False
        active_claims = len(_vendor_conversion_active_claims)
    if active_claims:
        logger.warning(
            "Retaining vendor ledger lock until process exit; %d claim(s) active",
            active_claims,
        )
        return
    _release_vendor_conversion_service_lock()


# ── Health Check ────────────────────────────────────────────────────

@app.get("/api/interpretation/confidence-guide")
def get_confidence_guide(
    language: str = Query("zh", enum=["zh", "en"]),
):
    """Return the canonical, presentation-only confidence interpretation."""

    return confidence_guide_payload(language)


@app.get("/api/health")
async def health_check():
    """Health check endpoint for monitoring and load balancers."""
    disk = storage_status(OUTPUT_DIR)
    runtime_manifest = await asyncio.to_thread(_current_runtime_manifest_status)
    authentic_reference_ready = bool(
        (runtime_manifest.get("authentic_reference") or {}).get("passed")
    )
    screening_catalog = _screening_catalog_summary()
    screening_runtime = (
        (runtime_manifest.get("role_catalogs") or {}).get("screening") or {}
    )
    screening_catalog["library_coverage"] = screening_runtime.get(
        "library_coverage",
        {"status": "not_audited"},
    )
    return {
        "status": "ok" if runtime_manifest.get("release_ready") else "degraded",
        # 前端显式展示当前运行实例，避免连到 V5/V6 或其他端口服务。
        "service_generation": SERVICE_GENERATION,
        "service_version": app.version,
        "runtime_snapshot_id": _current_runtime_snapshot_id(),
        "runtime_manifest": runtime_manifest,
        "role_catalogs": {
            "compound_roles": _compound_role_catalog_summary(),
            "screening": screening_catalog,
        },
        "default_port": DEFAULT_PORT,
        "active_port": int(os.getenv("MSSYSTEM_PORT", str(DEFAULT_PORT))),
        "project_root": str(PROJECT_ROOT),
        "library_loaded": _lib_loaded,
        "library_spectra": _lib_spectra_count,
        "library_metadata_repairs": library_metadata_repair_status(_lib_records, _lib_loaded),
        "sql_metadata_repairs": sql_metadata_repair_status(),
        "source_compound_mass_repairs": source_compound_mass_status(_lib_records, _lib_loaded),
        "tasks_count": len(_tasks),
        "capabilities": {
            "v11_flagship_entity_ranker": True,
            "primary_identity_endpoint": "/api/tasks/{task_id}/flagship",
            "v12_evirank_mv_shadow": True,
            "shadow_identity_endpoint": "/api/tasks/{task_id}/evirank-mv",
            "evitwin_ms_authentic_display_ranker": authentic_reference_ready,
            "evitwin_ms_endpoint": "/api/tasks/{task_id}/evitwin-ms",
            "primary_user_conclusion_endpoint": "/api/tasks/{task_id}/compounds",
            "compound_level_deduplication": True,
            "expert_rank_influence": False,
            "classification_rank_influence": False,
            "candidate_pathway_rank_influence": False,
            "blank_evidence_v2": True,
            "blank_filter_default": DEFAULT_BLANK_FILTER_POLICY,
            "blank_roles": sorted(_BLANK_ROLES),
            "isopair_ms": bool(USE_ISOPAIR_MS),
            "isopair_identity_scope": "2d_connectivity_only",
            "isopair_stereochemistry_claimed": False,
            "versioned_identification_contract": True,
            "candidate_independent_entity_graph": True,
            "project_calibration_fail_closed": True,
            "evigate2_final_claims": True,
            "lightweight_result_projection": True,
            "module_effect_inspector": True,
            "eviacquire_ms_shadow_planner": True,
            "evograph_fragmentation_shadow": True,
            "v16_dual_polarity_shadow_evidence_graph": True,
            "dual_polarity_shadow_endpoint": "/api/workflows/dual-polarity-shadow",
            "dual_polarity_changes_v15_top1": False,
            "v16_cross_occurrence_consensus_ranker": True,
            "cross_occurrence_ranking_endpoint": "/api/workflows/cross-occurrence-ranking",
            "cross_occurrence_activation": "explicit_same_lc_method_cohort_only",
            "cross_occurrence_candidate_scope": "existing_native_same_formula_top5",
            "candidate_independent_followup_linkage": True,
            "auditable_task_state_machine": True,
        },
        "storage": {
            "shared_data_root": str(SHARED_DATA_ROOT),
            "runtime_root": str(RUNTIME_ROOT),
            "vendor_conversion_cache_root": str(VENDOR_CONVERSION_CACHE_ROOT),
            "vendor_conversion_ledger": str(VENDOR_CONVERSION_LEDGER_FILE),
            "vendor_conversion_max_workers": _VENDOR_CONVERSION_MAX_WORKERS,
            "library_write_enabled": LIBRARY_WRITE_ENABLED,
            "free_gib": disk["free_gib"],
            "reserve_gib": disk["reserve_gib"],
            "ok": disk["ok"],
        },
    }


@app.get("/api/desktop/shutdown-status")
async def desktop_shutdown_status():
    """Small live task snapshot; never materialize results or inspect libraries."""
    counts: dict[str, int] = {}
    known = {"pending", "uploaded", "uploading", "converting", "queued", "running", "finalizing",
             "stopping", "cancel_requested", "cancelling", "completed", "failed",
             "cancelled", "canceled", "deleted"}
    with _tasks_lock:
        for task in _tasks.values():
            status = str(task.get("status") or "unknown").casefold()
            counts[status] = counts.get(status, 0) + 1
    return {
        "project_root": str(PROJECT_ROOT), "runtime_root": str(RUNTIME_ROOT),
        "port": int(os.getenv("MSSYSTEM_PORT", str(DEFAULT_PORT))),
        "task_counts": counts, "task_state_unknown": any(key not in known for key in counts),
    }


@app.get("/api/system/info")
async def system_info():
    """Return the single-source runtime identity used by the V16 client."""
    from utils.openms_backend import probe_openms

    host = os.getenv("MSSYSTEM_HOST", "127.0.0.1")
    port = int(os.getenv("MSSYSTEM_PORT", str(DEFAULT_PORT)))
    disk = storage_status(OUTPUT_DIR)
    runtime_status = await asyncio.to_thread(_current_runtime_manifest_status)
    authentic_reference = runtime_status.get("authentic_reference") or {}
    authentic_reference_ready = bool(authentic_reference.get("passed"))
    return {
        "service_generation": SERVICE_GENERATION,
        "service_version": app.version,
        "runtime_snapshot_id": _current_runtime_snapshot_id(),
        "project_root": str(PROJECT_ROOT),
        "host": host,
        "port": port,
        "base_url": f"http://{host}:{port}",
        "api_base": "/api",
        "ui_and_api_same_origin": True,
        "runtime_root": str(RUNTIME_ROOT),
        "shared_data_root": str(SHARED_DATA_ROOT),
        "storage": {
            "free_gib": disk["free_gib"],
            "reserve_gib": disk["reserve_gib"],
            "ok": disk["ok"],
            "vendor_conversion_cache_root": str(VENDOR_CONVERSION_CACHE_ROOT),
            "vendor_conversion_max_workers": _VENDOR_CONVERSION_MAX_WORKERS,
        },
        "database": {
            "full": str(DB_PATH),
            "small": str(DB_PATH_SMALL),
            "loaded": _lib_loaded,
            "loaded_name": _lib_name,
            "loaded_spectra": _lib_spectra_count,
        },
        "network_engine": {
            "name": "V15 condition-aware evidence network",
            "node_mode": (
                "OpenMS LC-MS feature consensus with explicit "
                "MS2 precursor-event fallback"
            ),
            "layers": ["spectral", "analog", "neutral_loss", "ion_identity"],
            "singletons_are_families": False,
        },
        "primary_identity_selector": {
            "name": "V15-Flagship / Conditioned Multi-Library Evidence",
            "decision_unit": "candidate-independent precursor entity",
            "expert_rank_influence": False,
            "classification_rank_influence": False,
            "candidate_pathway_rank_influence": False,
            "stereochemistry_claimed": False,
        },
        "shadow_identity_challenger": {
            "name": "V12-Flagship / EviRank-MV",
            "endpoint": "/api/tasks/{task_id}/evirank-mv",
            "deployment_state": "shadow",
            "production_effect": False,
            "score_is_probability": False,
        },
        "dual_polarity_shadow_graph": {
            "name": "V16 candidate-independent dual-polarity evidence graph",
            "endpoint": "/api/workflows/dual-polarity-shadow",
            "deployment_state": "shadow",
            "pairing": "explicit completed POS task plus completed NEG task",
            "edge_inputs": ["reliable_neutral_mass", "retention_time"],
            "candidate_attachment_after_graph_freeze": True,
            "production_effect": False,
            "compound_count_increment": 0,
        },
        "authentic_reference_ranker": {
            "name": "EviTwin-MS 标品谱联合排序",
            "endpoint": "/api/tasks/{task_id}/evitwin-ms",
            "deployment_state": (
                "active_authentic_reference_display_ranker"
                if authentic_reference_ready
                else "unavailable_fail_closed"
            ),
            "candidate_display_order_active": authentic_reference_ready,
            "deployment_audit": authentic_reference,
            "legacy_v11_v12_artifacts_mutated": False,
            "formal_level1_release_increment": 0,
        },
        "openms_backend": await asyncio.to_thread(probe_openms),
    }


@app.get("/api/modules/governance")
async def get_module_governance():
    """Return machine-validated module states and their runtime assertions."""
    from config import settings as runtime_settings
    from governance.module_registry import audit_module_registry, load_module_registry

    registry_path = PROJECT_ROOT / "config" / "module_governance.json"
    try:
        registry = await asyncio.to_thread(load_module_registry, registry_path)
        audit = await asyncio.to_thread(
            audit_module_registry,
            registry,
            project_root=PROJECT_ROOT,
            runtime_settings=runtime_settings,
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Module governance registry is unavailable: {exc}",
        ) from exc
    return {
        "registry": registry.model_dump(mode="json"),
        "audit": audit,
        "runtime_snapshot_id": _current_runtime_snapshot_id(),
        "claim_compiler": "EviGate 2.0",
    }


@app.get("/api/tasks/{task_id}/module-effects")
async def get_task_module_effects(task_id: str):
    """Explain executed module effects without changing analytical output."""

    from governance.module_effects import build_module_effect_summary

    from acquisition.planner import plan_store_path, plan_store_summary
    from matching.evidiscrim import (
        evidiscrim_artifact_path,
        evidiscrim_artifact_summary,
    )

    task = await asyncio.to_thread(_ensure_task_presentation_current, task_id)
    task_view = dict(task)
    task_view["active_evidence_summary"] = await asyncio.to_thread(
        plan_store_summary, plan_store_path(TASK_DIR, task_id)
    )
    task_view["evidiscrim_summary"] = await asyncio.to_thread(
        evidiscrim_artifact_summary,
        evidiscrim_artifact_path(TASK_DIR, task_id),
    )
    return build_module_effect_summary(task_view)


@app.post("/api/tasks/{task_id}/evidiscrim")
async def create_task_evidiscrim(
    task_id: str,
    request: EviDiscrimRunRequest = EviDiscrimRunRequest(),
):
    """Materialize the FIORA-backed EviDiscrim-MS shadow sidecar."""

    from matching.evidiscrim import (
        evidiscrim_artifact_path,
        load_evidiscrim_policy,
        read_evidiscrim_artifact,
        results_sha256,
        stream_select_evidiscrim_rows,
    )
    from matching.openworld_shadow import (
        FIORA_ALLOWED_ADDUCTS,
        FioraRuntime,
        build_fiora_evidiscrim_shadow,
        write_shadow_artifact,
    )

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="EviDiscrim-MS is available only for completed tasks",
        )
    policy = await asyncio.to_thread(
        load_evidiscrim_policy,
        PROJECT_ROOT / "config" / "evidiscrim_policy.json",
    )
    results_path = TASK_DIR / f"{task_id}_results.json"
    streaming = bool(
        results_path.is_file()
        and results_path.stat().st_size >= 128 * 1024 * 1024
    )
    source_digest_override = None
    selection_audit = None
    source_stat_token = None
    if streaming:
        source_stat = results_path.stat()
        source_stat_token = (source_stat.st_size, source_stat.st_mtime_ns)
        results, source_digest_override, selection_audit = await asyncio.to_thread(
            stream_select_evidiscrim_rows,
            results_path,
            max_spectra=request.max_spectra,
            max_candidates_per_spectrum=request.max_candidates_per_spectrum,
            supported_adducts=FIORA_ALLOWED_ADDUCTS,
            policy=policy,
        )
        selected_stat = results_path.stat()
        if (selected_stat.st_size, selected_stat.st_mtime_ns) != source_stat_token:
            raise HTTPException(
                status_code=409,
                detail="Source result changed during EviDiscrim input selection",
            )
    else:
        results = await asyncio.to_thread(_get_task_results, task_id)
    current_source_digest = source_digest_override or results_sha256(results)
    artifact_path = evidiscrim_artifact_path(TASK_DIR, task_id)
    if artifact_path.is_file() and not request.refresh:
        try:
            existing = await asyncio.to_thread(
                read_evidiscrim_artifact, artifact_path
            )
        except (OSError, ValueError):
            existing = None
        if (
            existing
            and existing.get("source_results_sha256") == current_source_digest
            and existing.get("policy_id") == policy.get("policy_id")
            and existing.get("status") in {"completed_shadow", "not_applicable"}
        ):
            existing["cache_hit"] = True
            existing["artifact"] = str(artifact_path.resolve())
            return existing
    runtime = FioraRuntime(
        executable=Path("/opt/anaconda3/envs/fiora/bin/fiora-predict"),
        model_path=Path(
            "/Users/zizi/Desktop/fiora/fiora/resources/models/"
            "fiora_OS_v1.0.0.pt"
        ),
    )
    try:
        payload = await asyncio.to_thread(
            build_fiora_evidiscrim_shadow,
            results,
            task_id=task_id,
            runtime=runtime,
            policy=policy,
            collision_energy_grid=(request.collision_energies or None),
            max_candidates_per_spectrum=request.max_candidates_per_spectrum,
            max_spectra=request.max_spectra,
            source_results_sha256_override=source_digest_override,
            source_results_immutability_verified=streaming,
        )
        if streaming:
            final_stat = results_path.stat()
            if (final_stat.st_size, final_stat.st_mtime_ns) != source_stat_token:
                raise HTTPException(
                    status_code=409,
                    detail="Source result changed during EviDiscrim model execution",
                )
            payload["runtime"]["input_selection"] = selection_audit
            payload["runtime"]["source_results_path"] = str(
                results_path.resolve()
            )
            payload = EviDiscrimArtifact.model_validate(payload).model_dump(
                mode="json"
            )
        await asyncio.to_thread(write_shadow_artifact, payload, artifact_path)
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"EviDiscrim-MS generation failed: {exc}",
        ) from exc
    payload["cache_hit"] = False
    payload["artifact"] = str(artifact_path.resolve())
    return payload


@app.get("/api/tasks/{task_id}/evidiscrim")
async def get_task_evidiscrim(
    task_id: str,
    limit: int = Query(50, ge=1, le=200),
    offset: int = Query(0, ge=0),
):
    from matching.evidiscrim import (
        evidiscrim_artifact_path,
        evidiscrim_artifact_summary,
        read_evidiscrim_artifact,
    )

    with _tasks_lock:
        if task_id not in _tasks:
            raise HTTPException(status_code=404, detail="Task not found")
    path = evidiscrim_artifact_path(TASK_DIR, task_id)
    if not path.is_file():
        return {
            "status": "not_generated",
            "task_id": task_id,
            "module_id": "evidiscrim_candidate_contrast",
            "deployment_state": "shadow",
            "runtime_influence": False,
            "summary": evidiscrim_artifact_summary(path),
            "spectra": [],
            "total": 0,
            "offset": offset,
            "limit": limit,
            "has_more": False,
        }
    try:
        payload = await asyncio.to_thread(read_evidiscrim_artifact, path)
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"EviDiscrim-MS artifact is invalid: {exc}",
        ) from exc
    spectra = payload.get("spectra") or []
    payload["spectra"] = spectra[offset:offset + limit]
    payload["total"] = len(spectra)
    payload["offset"] = offset
    payload["limit"] = limit
    payload["has_more"] = offset + limit < len(spectra)
    return payload


@app.post("/api/tasks/{task_id}/acquisition-plans")
async def create_task_acquisition_plans(
    task_id: str, request: AcquisitionPlanRequest = AcquisitionPlanRequest()
):
    """Materialize additive shadow plans for unresolved candidate pairs."""

    from pipeline.active_evidence import materialize_active_evidence_plans

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Acquisition planning is available only for completed tasks",
        )
    results = await asyncio.to_thread(_get_task_results, task_id)
    try:
        store = await asyncio.to_thread(
            materialize_active_evidence_plans,
            dict(task),
            results,
            task_dir=TASK_DIR,
            project_root=PROJECT_ROOT,
            runtime_snapshot_id=_current_runtime_snapshot_id(),
            max_plans=request.max_plans,
            max_actions_per_plan=request.max_actions_per_plan,
            refresh=request.refresh,
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500, detail=f"EviAcquire-MS planning failed: {exc}"
        ) from exc
    return store


@app.get("/api/tasks/{task_id}/acquisition-plans")
async def get_task_acquisition_plans(
    task_id: str,
    limit: int = Query(50, ge=1, le=200),
    offset: int = Query(0, ge=0),
):
    from acquisition.planner import plan_store_path, read_plan_store

    with _tasks_lock:
        if task_id not in _tasks:
            raise HTTPException(status_code=404, detail="Task not found")
    store = await asyncio.to_thread(
        read_plan_store, plan_store_path(TASK_DIR, task_id)
    )
    plans = store.get("plans") or []
    return {
        "schema_version": store.get("schema_version"),
        "metadata": store.get("metadata") or {},
        "plans": plans[offset:offset + limit],
        "followups": store.get("followups") or [],
        "total": len(plans),
        "offset": offset,
        "limit": limit,
        "has_more": offset + limit < len(plans),
        "module_state": "shadow",
    }


@app.get("/api/tasks/{task_id}/acquisition-plans/{plan_id}/export")
async def export_task_acquisition_plan(
    task_id: str,
    plan_id: str,
    format: Literal["vendor-neutral", "thermo-qe", "agilent-qtof"] = Query(
        "vendor-neutral"
    ),
):
    from acquisition.exporters.agilent_qtof import export_agilent_qtof_advisory
    from acquisition.exporters.thermo_qe import export_thermo_qe_advisory
    from acquisition.exporters.vendor_neutral import export_plan_csv
    from acquisition.planner import get_plan, plan_store_path

    try:
        plan = await asyncio.to_thread(
            get_plan, plan_store_path(TASK_DIR, task_id), plan_id
        )
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Acquisition plan not found") from exc
    if format == "thermo-qe":
        content = export_thermo_qe_advisory(plan)
    elif format == "agilent-qtof":
        content = export_agilent_qtof_advisory(plan)
    else:
        content = export_plan_csv(plan)
    filename = f"{task_id}_{plan_id}_{format}_advisory.csv"
    return Response(
        content=content,
        media_type="text/csv; charset=utf-8",
        headers={"Content-Disposition": attachment_disposition(filename)},
    )


@app.post("/api/tasks/{task_id}/acquisition-plans/{plan_id}/link-followup")
async def link_task_acquisition_followup(
    task_id: str, plan_id: str, request: FollowupLinkRequest
):
    from acquisition.planner import append_followup, get_plan, plan_store_path
    from acquisition.reconciliation import link_followup_observations

    path = plan_store_path(TASK_DIR, task_id)
    try:
        plan = await asyncio.to_thread(get_plan, path, plan_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Acquisition plan not found") from exc
    with _tasks_lock:
        followup_task = _tasks.get(request.followup_task_id)
    if not followup_task:
        raise HTTPException(status_code=404, detail="Follow-up task not found")
    if followup_task.get("status") != "completed":
        raise HTTPException(status_code=409, detail="Follow-up task is not completed")
    followup_results = await asyncio.to_thread(
        _get_task_results, request.followup_task_id
    )
    delta = await asyncio.to_thread(
        link_followup_observations,
        plan,
        request.followup_task_id,
        followup_results,
        mz_tolerance_ppm=request.mz_tolerance_ppm,
        rt_tolerance_min=request.rt_tolerance_min,
    )
    return await asyncio.to_thread(append_followup, path, delta)


@app.get("/api/tasks/{task_id}/acquisition-plans/{plan_id}")
async def get_task_acquisition_plan(task_id: str, plan_id: str):
    from acquisition.planner import get_plan, plan_store_path

    try:
        return await asyncio.to_thread(
            get_plan, plan_store_path(TASK_DIR, task_id), plan_id
        )
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Acquisition plan not found") from exc


# ── Library ─────────────────────────────────────────────────────────

@app.get("/api/library/status", response_model=LibraryStatus)
async def library_status():
    """Get current library status."""
    return LibraryStatus(
        loaded=_lib_loaded,
        name=_lib_name,
        spectra_count=_lib_spectra_count,
        load_time_sec=_lib_load_time,
    )


@app.get("/api/library/source-registry")
@app.get("/api/library/public-sources", include_in_schema=False)
async def public_source_registry():
    """Return runtime-identical polarity, naming and structure-quality audits."""
    from matching.library_polarity import (
        compound_quality_snapshot,
        load_cached_snapshot,
        recompute_polarity_snapshot,
    )

    registry_path = PROJECT_ROOT / "config" / "public_source_registry.json"
    try:
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Public-source registry is unreadable: {exc}",
        ) from exc
    cache_path = PROJECT_ROOT / "config" / "library_polarity_snapshot.json"
    snapshot = load_cached_snapshot(cache_path, DB_PATH)
    if snapshot is None:
        snapshot = await asyncio.to_thread(
            recompute_polarity_snapshot, DB_PATH
        )
        snapshot["cache_status"] = (
            "recomputed_in_memory; run scripts/recompute_library_polarity.py "
            "--apply to persist"
        )
    else:
        snapshot["cache_status"] = "current"
    registry["database_snapshot"] = snapshot

    structure = registry.get("structure_quality_snapshot")
    if not isinstance(structure, dict) or not structure:
        structure = await asyncio.to_thread(compound_quality_snapshot, DB_PATH)
    registry["structure_quality_snapshot"] = structure

    quarantine_path = (
        PROJECT_ROOT / "matching" / "compound_name_aliases_quarantine.json"
    )
    try:
        alias_quarantine = json.loads(
            quarantine_path.read_text(encoding="utf-8")
        )
        registry["compound_name_alias_audit"] = {
            "quarantined_aliases": len(
                alias_quarantine.get("aliases") or []
            ),
            "reason_counts": alias_quarantine.get("reason_counts") or {},
            "policy": alias_quarantine.get("policy") or "",
        }
    except Exception:
        registry["compound_name_alias_audit"] = {
            "quarantined_aliases": 0,
            "status": "no_persisted_quarantine_audit",
        }

    mapped = int(snapshot.get("mapped_compounds") or 0)
    negative = int(snapshot.get("compounds_with_negative") or 0)
    positive = int(snapshot.get("compounds_with_positive") or 0)
    dual = int(snapshot.get("both_polarities") or 0)
    positive_only = int(snapshot.get("positive_only_compounds") or 0)
    registry["coverage_summary"] = {
        "positive_compound_coverage_pct": (
            round(100.0 * positive / mapped, 2) if mapped else None
        ),
        "negative_compound_coverage_pct": (
            round(100.0 * negative / mapped, 2) if mapped else None
        ),
        "dual_polarity_coverage_pct": (
            round(100.0 * dual / mapped, 2) if mapped else None
        ),
        "positive_only_gap_compounds": positive_only,
        "headline": (
            "V11 is polarity-neutral. Positive and negative coverage are "
            "database quality dimensions, not a system theme; all counts use "
            "the same explicit-ion-mode → adduct → signed-charge rule."
        ),
    }
    return registry


_LIBRARY_COUNTS_SNAPSHOT_PATH = (
    PROJECT_ROOT / "config" / "library_counts_snapshot.json"
)


def _library_counts_snapshot_entry(
    mode: str, path: Path
) -> dict[str, Any] | None:
    """Read a stat-bound count snapshot without opening the multi-GB SQLite DB."""

    try:
        payload = json.loads(
            _LIBRARY_COUNTS_SNAPSHOT_PATH.read_text(encoding="utf-8")
        )
        if payload.get("schema_version") != (
            "mssystem.library-counts-snapshot/1.0.0"
        ):
            return None
        if payload.get("load_name_only_leads") is not bool(LOAD_NAME_ONLY_LEADS):
            return None
        item = (payload.get("libraries") or {}).get(mode) or {}
        stat = path.stat()
        if Path(str(item.get("path") or "")).expanduser().resolve() != path.resolve():
            return None
        if int(item.get("size_bytes") or -1) != int(stat.st_size):
            return None
        if int(item.get("mtime_ns") or -1) != int(stat.st_mtime_ns):
            return None
        declared_snapshot_id = str(item.get("snapshot_id") or "").strip()
        if declared_snapshot_id:
            runtime_manifest = json.loads(
                (PROJECT_ROOT / "config" / "runtime_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            runtime_snapshot_id = str(
                (
                    ((runtime_manifest.get("library") or {}).get("production_full"))
                    or {}
                ).get("snapshot_id")
                or ""
            ).strip()
            if runtime_snapshot_id != declared_snapshot_id:
                return None
        counts = item.get("counts") or {}
        required = {
            "spectra",
            "physical_spectra",
            "source_searchable_spectra",
            "runtime_searchable_spectra",
            "excluded_by_source",
            "excluded_by_runtime_structure_gate",
            "runtime_projection_exact",
            "compounds",
            "source_tiers",
        }
        if not isinstance(counts, Mapping) or not required.issubset(counts):
            return None
        return {
            **dict(counts),
            "count_snapshot_status": "current",
            "count_snapshot_generated_at": payload.get("generated_at"),
        }
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def _library_counts_sync_unlocked():
    """Return physical, provenance-eligible and actually loaded counts.

    A physical SQLite row is not automatically searchable: predicted/MS1/
    quarantined tiers are excluded first, then runtime structure/precursor/name
    gates are applied.  Keeping these counts separate prevents the Full option
    from advertising 4.6M as if every record participated in Top-K search.
    """
    result = {}
    loaded_mode = None
    if _lib_loaded:
        loaded_mode = "small" if _lib_name.lower().startswith("small") else "full"
    for mode, path in (("small", DB_PATH_SMALL), ("full", DB_PATH)):
        entry = {"available": bool(path.exists())}
        if not entry["available"]:
            entry.update(
                spectra=0, physical_spectra=0, source_searchable_spectra=0,
                loaded_searchable_spectra=None, compounds=0, source_tiers={},
            )
            result[mode] = entry
            continue
        try:
            import sqlite3 as _sqlite3

            stat = path.stat()
            signature = (
                int(stat.st_size),
                int(stat.st_mtime_ns),
                bool(LOAD_NAME_ONLY_LEADS),
            )
            cache_key = str(path.resolve())
            with _library_count_cache_lock:
                cached = _library_count_cache.get(cache_key)
            if cached and cached[0] == signature:
                base_entry = dict(cached[1])
            else:
                base_entry = _library_counts_snapshot_entry(mode, path)
            if not cached or cached[0] != signature:
                if base_entry is not None:
                    with _library_count_cache_lock:
                        _library_count_cache[cache_key] = (
                            signature, dict(base_entry)
                        )
                else:
                    conn = _sqlite3.connect(f"file:{path}?mode=ro", uri=True)
                    physical = conn.execute(
                    "SELECT COUNT(*) FROM Spectra "
                    "WHERE precursor_mz IS NOT NULL AND precursor_mz > 0"
                ).fetchone()[0]
                    cp = conn.execute("SELECT COUNT(*) FROM Compounds").fetchone()[0]
                    tables = {
                    row[0]
                    for row in conn.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"
                    )
                }
                    columns = {
                    row[1]
                    for row in conn.execute("PRAGMA table_info(Spectra)").fetchall()
                }
                    compound_columns = {
                    row[1]
                    for row in conn.execute("PRAGMA table_info(Compounds)").fetchall()
                }
                    tier_counts = {}
                    if "source_tier" in columns:
                        tier_counts = {
                        str(tier or "unknown"): int(count)
                        for tier, count in conn.execute(
                            "SELECT source_tier, COUNT(*) FROM Spectra "
                            "WHERE precursor_mz IS NOT NULL AND precursor_mz > 0 "
                            "GROUP BY source_tier"
                        )
                        }
                    searchable_source = sum(
                    tier_counts.get(tier, 0)
                    for tier in (
                        "curated", "experimental", "massbank_supplement", "merged",
                        "legacy_unresolved",
                    )
                    ) if tier_counts else physical

                # The production default rejects placeholder/no-structure
                # compounds before building the in-memory index. Publish the
                # same projected count even before the user loads Full, so the
                # UI cannot appear to "lose" records when 3.81M source-eligible
                # rows become 3.36M runtime-searchable rows.
                    runtime_searchable = searchable_source
                    excluded_no_structure = 0
                    can_project_runtime = (
                    not LOAD_NAME_ONLY_LEADS
                    and {"Mapping", "Compounds"}.issubset(tables)
                    and "smiles" in compound_columns
                    )
                    if can_project_runtime:
                        tier_clause = (
                        "lower(COALESCE(s.source_tier,'')) IN "
                        "('curated','experimental','massbank_supplement',"
                        "'merged','legacy_unresolved')"
                        if "source_tier" in columns else "1=1"
                        )
                        excluded_no_structure = int(conn.execute(
                        f"""
                        WITH missing_compounds AS MATERIALIZED (
                            SELECT id
                            FROM Compounds
                            WHERE lower(trim(COALESCE(smiles,''))) IN
                                  ('','nan','none','null','na','n/a','n.a.')
                        )
                        SELECT COUNT(*)
                        FROM missing_compounds c
                        JOIN Mapping m ON m.compound_id = c.id
                        JOIN Spectra s ON s.id = m.spectrum_id
                        WHERE s.precursor_mz IS NOT NULL
                          AND s.precursor_mz > 0
                          AND {tier_clause}
                        """
                        ).fetchone()[0])
                        runtime_searchable = max(
                        0, searchable_source - excluded_no_structure
                        )
                    conn.close()
                    base_entry = {
                    # ``spectra`` remains a compatibility alias for old clients.
                    "spectra": physical,
                    "physical_spectra": physical,
                    "source_searchable_spectra": searchable_source,
                    "runtime_searchable_spectra": runtime_searchable,
                    "excluded_by_source": max(
                        0, physical - searchable_source
                    ),
                    "excluded_by_runtime_structure_gate": (
                        excluded_no_structure
                    ),
                    "runtime_projection_exact": can_project_runtime,
                    "compounds": cp,
                    "source_tiers": tier_counts,
                    }
                    with _library_count_cache_lock:
                        _library_count_cache[cache_key] = (
                            signature, dict(base_entry)
                        )
            loaded_count = _lib_spectra_count if loaded_mode == mode else None
            entry.update(base_entry)
            entry.update(
                loaded_searchable_spectra=loaded_count,
                loaded=(loaded_mode == mode),
            )
        except Exception as e:
            entry.update(
                spectra=0, physical_spectra=0, source_searchable_spectra=0,
                loaded_searchable_spectra=None, compounds=0,
                source_tiers={}, error=str(e),
            )
        result[mode] = entry
    return result


@app.get("/api/library/counts")
async def library_counts():
    """Return cached library count layers without blocking the API event loop."""
    def compute():
        # Serialize only the first/stale scan. Concurrent dashboard requests
        # then reuse the in-memory signature cache instead of scanning twice.
        with _library_count_compute_lock:
            return _library_counts_sync_unlocked()

    return await asyncio.to_thread(compute)


@app.get("/api/library/ce-profiles/summary")
async def library_ce_profile_summary():
    """Return the read-only multi-source CE sidecar coverage summary."""
    try:
        from matching.ce_profiles import profile_database_summary
        return profile_database_summary()
    except Exception as exc:
        return {
            "available": False,
            "error": str(exc),
        }


async def _resume_recovered_tasks(db_mode: str) -> dict[str, Any]:
    """Re-queue restart-preserved tasks after their requested library loads."""
    def recovery_order(task_id: str) -> tuple[int, int, str]:
        task = _tasks.get(task_id) or {}
        raw_order = task.get("queue_order")
        try:
            queue_order = int(raw_order)
            has_saved_order = 0
        except (TypeError, ValueError):
            queue_order = 0
            has_saved_order = 1
        return (
            has_saved_order,
            queue_order,
            str(task.get("upload_time") or ""),
        )

    with _tasks_lock:
        candidates = sorted(
            (
                task_id for task_id in _restart_recovery_task_ids
                if task_id in _tasks
            ),
            key=recovery_order,
        )

    resumed: list[str] = []
    deferred: list[str] = []
    failed: list[dict[str, str]] = []
    for task_id in candidates:
        with _tasks_lock:
            task = _tasks.get(task_id)
            if not task or task.get("status") != TaskState.PENDING.value:
                _restart_recovery_task_ids.discard(task_id)
                continue
            raw_params = dict(task.get("params") or {})
            requested_mode = str(raw_params.get("db_mode") or "full").lower()
        if requested_mode != db_mode:
            deferred.append(task_id)
            continue
        try:
            recovered_params = AnalysisParams.model_validate(raw_params)
        except Exception as exc:
            with _tasks_lock:
                current = _tasks.get(task_id)
                if current:
                    current["restart_recovery"] = {
                        **dict(current.get("restart_recovery") or {}),
                        "status": "manual_action_required",
                        "error": f"历史参数无法验证：{exc}"[:1000],
                    }
                    current["message"] = (
                        "任务与原始文件已保留，但历史分析参数无法自动恢复；"
                        "请核对参数后点击运行"
                    )
            _restart_recovery_task_ids.discard(task_id)
            failed.append({"id": task_id, "error": str(exc)})
            continue
        try:
            await start_analysis(task_id, recovered_params)
        except Exception as exc:
            detail = str(exc.detail) if isinstance(exc, HTTPException) else str(exc)
            with _tasks_lock:
                current = _tasks.get(task_id)
                if current and current.get("status") == TaskState.PENDING.value:
                    current["restart_recovery"] = {
                        **dict(current.get("restart_recovery") or {}),
                        "status": "manual_action_required",
                        "error": detail[:1000],
                    }
                    current["message"] = (
                        "任务与原始文件已保留，但自动恢复未能排队；"
                        "请核对任务参数后重试"
                    )
            _restart_recovery_task_ids.discard(task_id)
            failed.append({"id": task_id, "error": detail})
            continue
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current:
                current["restart_recovery"] = {
                    **dict(current.get("restart_recovery") or {}),
                    "status": "resumed",
                    "resumed_at": _now_iso(),
                    "automatic_resume": True,
                }
        _restart_recovery_task_ids.discard(task_id)
        resumed.append(task_id)

    if resumed or failed:
        _save_tasks()
    return {
        "resumed": resumed,
        "deferred": deferred,
        "failed": failed,
    }


@app.post("/api/isomer-literature/search")
def isomer_literature_search(
    request: IsomerLiteratureSearchRequest,
):
    """Return cached PubMed leads without changing identification evidence."""
    return search_pubmed_isomer_evidence(
        request.candidate_names,
        formula=request.formula,
        max_results=request.max_results,
        cache_dir=OUTPUT_DIR / "isomer_literature_cache",
    )


@app.post("/api/library/load")
async def library_load(db_mode: str = Query("full", enum=["small", "full"])):
    """Load or reload the reference library without freezing unrelated APIs."""

    db_path = DB_PATH_SMALL if db_mode == "small" else DB_PATH

    if not db_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Database not found: {db_path}. Build it first with build_from_mgf_msp.py",
        )

    def load_in_worker() -> LibraryStatus:
        global _lib_records, _mz_index, _lib_loaded, _lib_name
        global _lib_spectra_count, _lib_load_time

        with _library_state_lock:
            with _tasks_lock:
                active = [
                    task_id for task_id, task in _tasks.items()
                    if _is_active_task_state(task.get("status"))
                ]
            if active:
                raise HTTPException(
                    status_code=409,
                    detail=(
                        "Cannot replace the loaded library while analyses are "
                        "queued or running"
                    ),
                )
            # A failed source/cache verification cannot leave stale records
            # advertised as ready for a new analysis.
            _lib_records = None
            _mz_index = None
            _lib_loaded = False
            _lib_name = ""
            _lib_spectra_count = 0
            _lib_load_time = 0.0
            t0 = time.time()
            try:
                records, mz_index = load_library(db_path)
            except Exception as exc:
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to load library: {exc}",
                ) from exc

            _lib_records = records
            _mz_index = mz_index
            _lib_loaded = True
            _lib_name = (
                f"Small (<800 Da, {len(records):,})"
                if db_mode == "small" else f"Full ({len(records):,})"
            )
            _lib_spectra_count = len(records)
            _lib_load_time = round(time.time() - t0, 1)
            return LibraryStatus(
                loaded=True,
                name=_lib_name,
                spectra_count=_lib_spectra_count,
                load_time_sec=_lib_load_time,
            )

    status = await asyncio.to_thread(load_in_worker)
    restart_recovery = await _resume_recovered_tasks(db_mode)
    logger.info(
        "Library loaded: %s (%d spectra in %.1fs); restart recovery resumed=%d, deferred=%d, failed=%d",
        status.name, status.spectra_count, status.load_time_sec,
        len(restart_recovery["resumed"]),
        len(restart_recovery["deferred"]),
        len(restart_recovery["failed"]),
    )
    return status


# ── Library Augment (tables / raw MS2 / audited images → 追加库) ────
USER_CONTRIB_RAW_DIR = USER_CONTRIB_DIR / "raw"
_IMPORT_ALLOWED_EXT = {
    ".csv", ".xlsx", ".xls",
    ".mgf", ".msp", ".mzml", ".mzxml", ".xml",
    ".json", ".png", ".jpg", ".jpeg", ".tif", ".tiff",
}


async def _save_contrib_files(files: list[UploadFile]) -> list[str]:
    """Atomically save a bounded contribution batch and return final paths."""
    if len(files) > MAX_BATCH_UPLOAD_FILES:
        raise HTTPException(
            status_code=400,
            detail=f"单批最多 {MAX_BATCH_UPLOAD_FILES} 个文件",
        )
    USER_CONTRIB_RAW_DIR.mkdir(parents=True, exist_ok=True)
    _require_disk_budget(
        USER_CONTRIB_RAW_DIR,
        MAX_BATCH_UPLOAD_SIZE_BYTES,
        "library contribution upload",
    )
    saved: list[str] = []
    batch_bytes = 0
    partial: Optional[Path] = None
    try:
        for uf in files:
            display_name = _safe_sample_filename(uf.filename or "")
            if not display_name:
                raise HTTPException(status_code=400, detail="投稿文件缺少文件名")
            ext = Path(display_name).suffix.lower()
            if ext not in _IMPORT_ALLOWED_EXT:
                raise HTTPException(
                    status_code=400,
                    detail=(
                        f"不支持的格式: {ext}（支持 CSV/Excel、MGF、MSP、质心 mzML，"
                        "以及带 V15 JSON 清单的 PNG/JPEG/TIFF 谱图）"
                    ),
                )
            dest = USER_CONTRIB_RAW_DIR / f"{uuid.uuid4().hex[:8]}_{display_name}"
            partial = USER_CONTRIB_RAW_DIR / f".{uuid.uuid4().hex}.part"
            total_bytes = 0
            chunk_size = 4 * 1024 * 1024
            with open(partial, "wb") as handle:
                while True:
                    chunk = await uf.read(chunk_size)
                    if not chunk:
                        break
                    total_bytes += len(chunk)
                    batch_bytes += len(chunk)
                    if total_bytes > MAX_UPLOAD_SIZE_BYTES:
                        raise HTTPException(
                            status_code=413,
                            detail=f"文件过大: {display_name}（上限 {MAX_UPLOAD_SIZE_MB} MB）",
                        )
                    if batch_bytes > MAX_BATCH_UPLOAD_SIZE_BYTES:
                        raise HTTPException(
                            status_code=413,
                            detail=f"投稿批次超过 {MAX_BATCH_UPLOAD_SIZE_MB} MB 总上限",
                        )
                    _require_disk_budget(
                        USER_CONTRIB_RAW_DIR,
                        len(chunk),
                        "library contribution upload",
                    )
                    handle.write(chunk)
            if total_bytes == 0:
                raise HTTPException(status_code=400, detail=f"空文件: {display_name}")
            partial.replace(dest)
            partial = None
            saved.append(str(dest))
            logger.info("库扩增接收: %s (%.2f KB)", display_name, total_bytes / 1024)
    except Exception:
        if partial is not None:
            partial.unlink(missing_ok=True)
        for path in saved:
            Path(path).unlink(missing_ok=True)
        raise
    finally:
        for uf in files:
            try:
                await uf.close()
            except Exception:
                pass
    return saved


@app.post("/api/library/import/preview")
async def library_import_preview(files: list[UploadFile] = File(...)):
    """Preview tables or audited raw-reference batches without a DB write."""
    if not files:
        raise HTTPException(status_code=400, detail="未提供文件")
    saved = await _save_contrib_files(files)
    try:
        from matching.library_import import preview_import
        result = preview_import([Path(p) for p in saved])
        result["saved_files"] = [Path(p).name for p in saved]
        return result
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("库扩增预览失败")
        raise HTTPException(status_code=500, detail=f"预览失败: {e}")
    finally:
        # Preview is read-only; do not accumulate UUID-prefixed temporary raw
        # copies.  The confirmed import request uploads and archives its own
        # content-hashed originals.
        for path in saved:
            Path(path).unlink(missing_ok=True)


@app.post("/api/library/import")
async def library_import(
    db_mode: str = Query("both", enum=["both", "full", "small"]),
    files: list[UploadFile] = File(...),
):
    """Confirm an audited reference batch and invalidate the search cache."""
    if not LIBRARY_WRITE_ENABLED:
        raise HTTPException(status_code=403, detail="This MSsystem instance is read-only")
    if not files:
        raise HTTPException(status_code=400, detail="未提供文件")
    saved = await _save_contrib_files(files)
    try:
        from matching.library_import import import_files
        selected = []
        if db_mode in {"both", "full"}:
            selected.append(DB_PATH)
        if db_mode in {"both", "small"}:
            selected.append(DB_PATH_SMALL)
        source_bytes = sum(Path(path).stat().st_size for path in saved)
        _require_disk_budget(
            SHARED_DATA_ROOT,
            estimate_incremental_import_bytes(selected, source_bytes),
            "shared library import",
        )
        with _library_state_lock:
            with _tasks_lock:
                active = [
                    task_id for task_id, task in _tasks.items()
                    if _is_active_task_state(task.get("status"))
                ]
            if active:
                raise HTTPException(
                    status_code=409,
                    detail=(
                        "Cannot modify the library while analyses are queued or running; "
                        f"active tasks: {', '.join(active[:5])}"
                    ),
                )
            result = import_files([Path(p) for p in saved], db_mode=db_mode)
            modified = bool(result.get("database_modified", True))
            if modified:
                _invalidate_library_state()
            result["library_needs_reload"] = modified
        logger.info("库扩增完成: %s", result)
        return result
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("库扩增导入失败")
        raise HTTPException(status_code=500, detail=f"导入失败: {e}")


# ── File Upload ─────────────────────────────────────────────────────

def _safe_sample_filename(filename: str) -> str:
    """Return a display/storage basename without trusting multipart paths."""
    normalized = str(filename or "").replace("\x00", "").replace("\\", "/")
    return Path(normalized).name.strip()


def _register_pending_task(
    display_name: str,
    destination: Path,
    total_bytes: int,
    *,
    extra_metadata: dict[str, Any] | None = None,
    task_id: str | None = None,
) -> str:
    """Register one managed spectrum file as a pending analysis task."""
    task_id = str(task_id or uuid.uuid4().hex[:12])
    if not re.fullmatch(r"[0-9a-f]{12}", task_id):
        raise ValueError("Managed task id must contain 12 lowercase hex characters")
    with _tasks_lock:
        existing = _tasks.get(task_id)
    if existing is not None:
        if (
            str(existing.get("filename") or "") == display_name
            and Path(str(existing.get("file_path") or "")).resolve()
            == destination.resolve()
        ):
            return task_id
        raise ValueError("Managed task id is already bound to a different upload")
    file_size_mb = total_bytes / (1024 * 1024)
    task = {
        "id": task_id,
        "filename": display_name,
        "file_path": str(destination),
        "file_size_mb": round(file_size_mb, 2),
        "status": "pending",
        "progress": 0,
        "message": "File uploaded, awaiting analysis",
        "upload_time": _now_iso(),
        "start_time": None,
        "end_time": None,
        "total_spectra": 0,
        "matched_spectra": 0,
        "unique_compounds": None,
        "match_rate": 0.0,
        "results": [],
        "error": None,
        "params": {},
    }
    if extra_metadata:
        task.update(extra_metadata)
    _initialize_task_state(task)
    with _tasks_lock:
        _tasks[task_id] = task
    return task_id


def _safe_vendor_relative_path(value: str) -> PurePosixPath:
    """Validate one browser-provided member path for a raw-directory upload."""
    normalized = str(value or "").replace("\x00", "").replace("\\", "/").strip()
    candidate = PurePosixPath(normalized)
    if (
        not normalized
        or candidate.is_absolute()
        or any(part in {"", ".", ".."} for part in candidate.parts)
    ):
        raise ValueError(f"Unsafe raw-directory member path: {value!r}")
    if not candidate.parts[0].lower().endswith(".d"):
        raise ValueError(
            "Every Agilent raw member must be inside one top-level .d directory."
        )
    return candidate


@app.get("/api/vendor/converter/status")
async def vendor_converter_status():
    """Report whether local ProteoWizard or its Docker backend can run."""
    from utils.vendor_raw_conversion import converter_status

    return await asyncio.to_thread(converter_status)


def _vendor_conversion_ledger():
    from utils.vendor_conversion_cache import ConversionLedger

    return ConversionLedger(
        VENDOR_CONVERSION_LEDGER_FILE,
        job_id="v16-local-vendor-conversion",
        lease_seconds=4 * 60 * 60,
    )


def _acquire_vendor_conversion_service_lock() -> bool:
    """Prove single-process ownership before forced ledger recovery."""

    global _vendor_conversion_service_lock_handle
    global _vendor_conversion_service_lock_path
    global _vendor_conversion_service_lock_refcount
    lock_path = VENDOR_CONVERSION_LEDGER_FILE.with_name(
        ".vendor_conversion_service.lock"
    ).resolve(strict=False)
    with _vendor_conversion_service_lock_guard:
        if _vendor_conversion_service_lock_handle is not None:
            if _vendor_conversion_service_lock_path == lock_path:
                _vendor_conversion_service_lock_refcount += 1
                return False
            raise RuntimeError(
                "Vendor conversion service already owns a different ledger lock"
            )
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        handle = lock_path.open("a+b")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            handle.close()
            raise RuntimeError(
                "Another V16 process owns the vendor conversion ledger"
            ) from exc
        _vendor_conversion_service_lock_handle = handle
        _vendor_conversion_service_lock_path = lock_path
        _vendor_conversion_service_lock_refcount = 1
        return True


def _release_vendor_conversion_service_lock() -> None:
    global _vendor_conversion_service_lock_handle
    global _vendor_conversion_service_lock_path
    global _vendor_conversion_service_lock_refcount
    with _vendor_conversion_service_lock_guard:
        if _vendor_conversion_service_lock_refcount > 1:
            _vendor_conversion_service_lock_refcount -= 1
            return
        handle = _vendor_conversion_service_lock_handle
        _vendor_conversion_service_lock_handle = None
        _vendor_conversion_service_lock_path = None
        _vendor_conversion_service_lock_refcount = 0
        if handle is None:
            return
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        finally:
            handle.close()


async def _await_thread_completion(
    operation: Any,
    *args: Any,
    **kwargs: Any,
) -> tuple[Any, asyncio.CancelledError | None]:
    """Finish an already-authorized disk operation before propagating cancel."""

    worker = asyncio.create_task(asyncio.to_thread(operation, *args, **kwargs))
    cancellation: asyncio.CancelledError | None = None
    while True:
        try:
            result = await asyncio.shield(worker)
            return result, cancellation
        except asyncio.CancelledError as exc:
            if cancellation is None:
                cancellation = exc


def _begin_vendor_registration_operation(marker: str) -> None:
    with _vendor_conversion_active_claims_lock:
        if not _vendor_conversion_accepting_claims:
            raise RuntimeError("Vendor conversion service is shutting down")
        if marker in _vendor_conversion_active_claims:
            raise RuntimeError("Vendor conversion registration is already running")
        _vendor_conversion_active_claims.add(marker)


def _end_vendor_registration_operation(marker: str) -> None:
    release_ledger_lock = False
    with _vendor_conversion_active_claims_lock:
        _vendor_conversion_active_claims.discard(marker)
        release_ledger_lock = bool(
            not _vendor_conversion_accepting_claims
            and not _vendor_conversion_active_claims
        )
    if release_ledger_lock:
        _release_vendor_conversion_service_lock()


def _vendor_conversion_key_lock(cache_key: str) -> threading.Lock:
    with _vendor_conversion_key_locks_guard:
        return _vendor_conversion_key_locks.setdefault(
            str(cache_key),
            threading.Lock(),
        )


def _fail_vendor_claim_if_converting(
    ledger: Any,
    item_id: str,
    claim_token: str,
    error: Exception,
) -> None:
    """Keep validated outputs resumable; fail only an unfinished conversion."""

    if not str(claim_token or ""):
        return
    try:
        item = ledger.get(item_id)
        if item and item.get("state") == "converting":
            ledger.mark_failed(item_id, str(error), claim_token)
    except Exception:
        logger.exception("Could not persist vendor conversion failure state")


def _cached_or_prepare_vendor_converter_status(
    *,
    accept_vendor_licenses: bool,
) -> dict[str, Any]:
    """Return one verified converter preflight shared by the two workers."""

    from utils.vendor_raw_conversion import prepare_converter

    with _vendor_converter_preflight_lock:
        prepared_at = _vendor_converter_preflight.get("prepared_at_monotonic")
        prepared_status = _vendor_converter_preflight.get("status")
        if (
            isinstance(prepared_at, (int, float))
            and time.monotonic() - float(prepared_at)
            <= _VENDOR_CONVERTER_PREFLIGHT_TTL_SECONDS
            and isinstance(prepared_status, Mapping)
        ):
            return deepcopy(dict(prepared_status))
        prepared = prepare_converter(
            accept_vendor_licenses=accept_vendor_licenses,
            allow_container_pull=False,
        )
        status = prepared.get("status")
        if prepared.get("ready") is not True or not isinstance(status, Mapping):
            raise RuntimeError("Vendor converter preflight did not return a ready status")
        _vendor_converter_preflight.update(
            {
                "prepared_at_monotonic": time.monotonic(),
                "status": deepcopy(dict(status)),
            }
        )
        return deepcopy(dict(status))


def _vendor_converter_identity(
    status: Mapping[str, Any],
    *,
    requested_backend: str,
    accept_vendor_licenses: bool,
) -> tuple[str, dict[str, Any]]:
    """Build the exact converter identity used by the cache and manifest."""

    from utils.vendor_raw_conversion import (
        ConverterUnavailableError,
        VendorLicenseNotAcceptedError,
    )

    local = status.get("local_msconvert") or {}
    docker = status.get("docker") or {}
    selected = str(requested_backend or "auto").strip().lower()
    if selected == "auto":
        selected = "local" if local.get("path") else "docker"
    if selected == "local":
        executable = str(local.get("path") or "")
        if not executable:
            raise ConverterUnavailableError("Local msconvert was not found.")
        probe = subprocess.run(
            [executable, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=15,
            check=False,
        )
        version = str(probe.stdout or probe.stderr or "").strip()
        if not version:
            raise ConverterUnavailableError(
                "Local msconvert did not report a version; conversion cache is fail-closed."
            )
        return selected, {
            "backend": "local",
            "executable": executable,
            "version": version,
        }
    if selected != "docker":
        raise ValueError("backend must be auto, local, or docker")
    if not accept_vendor_licenses:
        raise VendorLicenseNotAcceptedError(
            "ProteoWizard's Docker image requires explicit acceptance of "
            "the bundled vendor licences."
        )
    if not docker.get("daemon_available") or not docker.get("image_available"):
        raise ConverterUnavailableError(
            "The verified ProteoWizard Docker image is not available."
        )
    image_identity = str(docker.get("image_identity") or "")
    if not image_identity:
        raise ConverterUnavailableError(
            "The ProteoWizard image has no immutable image identity."
        )
    return selected, {
        "backend": "docker",
        "image": str(docker.get("image") or ""),
        "image_identity": image_identity,
        "platform": "linux/amd64",
        "wine_runner": str(docker.get("wine_runner") or "wine"),
        "apple_silicon_compatibility_pinned": bool(
            docker.get("apple_silicon_pinned_image")
        ),
    }


def _vendor_conversion_from_cache(
    *,
    raw_report: Mapping[str, Any],
    cache_hit: Any,
    materialization: Mapping[str, Any],
    destination: Path,
    peak_picking_strategy_requested: str,
    elapsed_seconds: float,
) -> dict[str, Any]:
    audit = cache_hit.manifest.get("conversion_audit") or {}
    output = deepcopy(audit.get("output") or {})
    if output.get("status") != "valid":
        raise ValueError("Cached conversion has no validated output metadata")
    output["path"] = str(destination.resolve())
    policy = deepcopy(audit.get("conversion_policy") or {})
    policy["peak_picking_strategy_requested"] = str(
        peak_picking_strategy_requested
    )
    policy["converter_preflight_reused"] = True
    policy["command_audit_without_source_path"] = [
        "<CONTENT_ADDRESSED_VENDOR_CONVERSION_CACHE>"
    ]
    return {
        "status": "completed",
        "format": f"{raw_report.get('format') or 'Vendor raw'} to mzML",
        "source_kind": raw_report.get("source_kind"),
        "source": deepcopy(dict(raw_report)),
        "output": output,
        "converter": deepcopy(audit.get("converter") or {}),
        "conversion_policy": policy,
        "duration_seconds": round(elapsed_seconds, 3),
        "timing": {
            "converter_seconds": 0.0,
            "cache_lookup_and_materialization_seconds": round(
                elapsed_seconds, 3
            ),
            "source_fingerprinting_passes_per_required_file": 1,
            "mzml_validation_passes": 0,
        },
        "cache": {
            "status": "hit",
            "hit": True,
            "cache_key": cache_hit.key,
            "cache_root": str(VENDOR_CONVERSION_CACHE_ROOT),
            "materialization_method": materialization.get("method"),
            "materialized_sha256": materialization.get("sha256"),
            "task_output_has_independent_inode": True,
        },
    }


class VendorConverterPrepareRequest(BaseModel):
    accept_vendor_licenses: bool = False
    allow_container_pull: bool = False


class VendorLocalConvertRequest(VendorConverterPrepareRequest):
    selection_token: str
    backend: Literal["auto", "local", "docker"] = "auto"
    intake_role: Literal["sample", "blank"] = "sample"
    blank_role: Literal[
        "unspecified",
        "background_blank",
        "solvent_blank",
        "method_blank",
        "extraction_blank",
        "carryover_blank",
    ] = "background_blank"
    peak_picking_strategy: Literal[
        "auto", "preserve", "vendor", "openms"
    ] = "auto"


def _require_loopback_request(request: Request) -> None:
    host = str(request.client.host if request.client else "").strip().lower()
    if host not in {"127.0.0.1", "::1", "localhost", "testclient"}:
        raise HTTPException(
            status_code=403,
            detail="Local raw-directory access is restricted to this Mac.",
        )
    # CORS is intentionally broad for the ordinary read APIs, so additionally
    # reject cross-site browser requests for operations that can open a native
    # picker, read a local directory, or start a multi-GB licensed download.
    origin = str(request.headers.get("origin") or "").strip()
    if origin:
        parsed = urlsplit(origin)
        origin_host = str(parsed.hostname or "").lower()
        origin_port = parsed.port or (443 if parsed.scheme == "https" else 80)
        request_port = request.url.port or (
            443 if request.url.scheme == "https" else 80
        )
        if (
            parsed.scheme != request.url.scheme
            or origin_host not in {"127.0.0.1", "::1", "localhost"}
            or origin_port != request_port
        ):
            raise HTTPException(
                status_code=403,
                detail="Cross-site local raw-data operations are not allowed.",
            )


def _choose_macos_vendor_directory() -> Optional[Path]:
    """Open the real macOS folder chooser, independent of browser support.

    The selected path may be either one complete ``.d`` directory or a parent
    directory whose immediate children are complete ``.d`` directories.
    """

    if sys.platform != "darwin":
        raise RuntimeError(
            "当前系统不支持本机目录选择器。"
        )
    script = (
        'set selectedFolder to choose folder with prompt '
        '"选择一个完整 .d，或选择含多个 .d 的父文件夹" '
        "with invisibles\n"
        "return POSIX path of selectedFolder"
    )
    result = subprocess.run(
        ["osascript", "-e", script],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=30 * 60,
        check=False,
    )
    if result.returncode != 0:
        error = (result.stderr or result.stdout or "").strip()
        if "-128" in error or "User canceled" in error:
            return None
        raise RuntimeError(error or "本机目录选择失败")
    selected = str(result.stdout or "").strip()
    if not selected:
        return None
    return Path(selected).expanduser().resolve()


def _choose_macos_thermo_directory() -> Optional[Path]:
    """Choose a folder containing Thermo RAW files, optionally under POS/NEG."""

    if sys.platform != "darwin":
        raise RuntimeError("当前系统不支持本机目录选择器。")
    script = (
        'set selectedFolder to choose folder with prompt '
        '"选择含 Thermo .raw 的文件夹（可含 POS/NEG 一层子目录）" '
        "with invisibles\n"
        "return POSIX path of selectedFolder"
    )
    result = subprocess.run(
        ["osascript", "-e", script],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=30 * 60,
        check=False,
    )
    if result.returncode != 0:
        error = (result.stderr or result.stdout or "").strip()
        if "-128" in error or "User canceled" in error:
            return None
        raise RuntimeError(error or "Thermo RAW 本机目录选择失败")
    selected = str(result.stdout or "").strip()
    return Path(selected).expanduser().resolve() if selected else None


def _thermo_polarity_hint(path: Path) -> list[str]:
    aliases = {
        "pos": "POS",
        "positive": "POS",
        "positive_mode": "POS",
        "neg": "NEG",
        "negative": "NEG",
        "negative_mode": "NEG",
    }
    for part in reversed(path.parts[:-1]):
        normalized = str(part).strip().lower().replace("-", "_").replace(" ", "_")
        if normalized in aliases:
            return [aliases[normalized]]
    return []


def _validate_local_thermo_raw(raw_file: Path) -> None:
    if (
        not raw_file.is_file()
        or raw_file.is_symlink()
        or raw_file.suffix.lower() != ".raw"
    ):
        raise ValueError("请选择完整的 Thermo Fisher .raw 文件")
    if raw_file.stat().st_size <= 0:
        raise ValueError(f"Thermo RAW 文件为空：{raw_file.name}")


def _discover_local_thermo_raw_files(
    selected_path: Path,
) -> tuple[list[tuple[Path, list[str]]], list[dict[str, str]]]:
    """Discover RAW files at depth 0/1 only; never recursively scan a tree."""

    selected = Path(selected_path).expanduser().resolve()
    if selected.is_file():
        _validate_local_thermo_raw(selected)
        return [(selected, _thermo_polarity_hint(selected))], []
    if not selected.is_dir():
        raise ValueError("请选择一个 Thermo .raw 或包含 .raw 的文件夹")

    candidates: list[Path] = []
    for child in selected.iterdir():
        if child.is_file() and not child.is_symlink() and child.suffix.lower() == ".raw":
            candidates.append(child.resolve())
        elif child.is_dir() and not child.is_symlink():
            candidates.extend(
                grandchild.resolve()
                for grandchild in child.iterdir()
                if grandchild.is_file()
                and not grandchild.is_symlink()
                and grandchild.suffix.lower() == ".raw"
            )
    candidates = sorted(
        set(candidates),
        key=lambda path: str(path.relative_to(selected)).casefold(),
    )
    if not candidates:
        raise ValueError(
            "所选目录及其直接子目录中没有 Thermo .raw 文件"
        )
    if len(candidates) > MAX_BATCH_UPLOAD_FILES:
        raise ValueError(
            f"一次最多选择 {MAX_BATCH_UPLOAD_FILES} 个 Thermo RAW；"
            f"当前发现 {len(candidates)} 个"
        )

    valid: list[tuple[Path, list[str]]] = []
    invalid: list[dict[str, str]] = []
    for candidate in candidates:
        try:
            _validate_local_thermo_raw(candidate)
            valid.append((candidate, _thermo_polarity_hint(candidate)))
        except (OSError, ValueError) as exc:
            invalid.append({"raw_name": candidate.name, "error": str(exc)})
    if not valid:
        raise ValueError("发现的 Thermo .raw 文件均为空或不可读")
    return valid, invalid


def _discover_local_agilent_directories(
    selected_path: Path,
) -> tuple[list[Path], list[dict[str, str]]]:
    """Resolve one selected ``.d`` or the valid ``.d`` children of a parent.

    Discovery is intentionally limited to immediate children.  It avoids an
    unbounded recursive scan while covering the normal batch layout where raw
    directories from one project/run are siblings in the same folder.
    """

    selected = Path(selected_path).expanduser().resolve()
    if not selected.is_dir():
        raise ValueError("请选择一个完整 .d，或选择含多个 .d 的父文件夹")

    if selected.suffix.lower() == ".d":
        _validate_local_agilent_directory(selected)
        return [selected], []

    candidates = sorted(
        (
            child.resolve()
            for child in selected.iterdir()
            if child.is_dir()
            and not child.is_symlink()
            and child.suffix.lower() == ".d"
        ),
        key=lambda path: path.name.casefold(),
    )
    if not candidates:
        raise ValueError(
            "所选父文件夹中没有找到直接子级 .d 目录；"
            "请选择一个 .d，或选择这些 .d 共同所在的父文件夹"
        )
    if len(candidates) > MAX_BATCH_UPLOAD_FILES:
        raise ValueError(
            f"一次最多选择 {MAX_BATCH_UPLOAD_FILES} 个 .d；"
            f"当前父文件夹包含 {len(candidates)} 个"
        )

    valid: list[Path] = []
    invalid: list[dict[str, str]] = []
    for candidate in candidates:
        try:
            _validate_local_agilent_directory(candidate)
            valid.append(candidate)
        except ValueError as exc:
            invalid.append({
                "raw_name": candidate.name,
                "error": str(exc),
            })
    if not valid:
        raise ValueError(
            "所选父文件夹中的 .d 均不完整；每个目录必须包含 "
            "AcqData/MSScan.bin、MSPeak.bin 和 AcqMethod.xml"
        )
    return valid, invalid


def _validate_local_agilent_directory(raw_dir: Path) -> None:
    from utils.vendor_raw_conversion import AGILENT_REQUIRED_FILES

    if not raw_dir.is_dir() or raw_dir.suffix.lower() != ".d":
        raise ValueError("请选择名称以 .d 结尾的完整 Agilent 原始目录")
    missing = [
        relative for relative in AGILENT_REQUIRED_FILES
        if not (raw_dir / relative).is_file()
    ]
    if missing:
        raise ValueError(
            "所选 .d 不完整，缺少：" + "、".join(missing)
        )


def _local_vendor_directory_size(raw_dir: Path) -> tuple[int, int]:
    if raw_dir.is_file():
        return raw_dir.stat().st_size, 1
    total_bytes = 0
    file_count = 0
    for root, directories, filenames in os.walk(raw_dir, followlinks=False):
        directories[:] = [
            name for name in directories
            if not (Path(root) / name).is_symlink()
        ]
        for filename in filenames:
            path = Path(root) / filename
            if path.is_symlink():
                continue
            try:
                total_bytes += path.stat().st_size
                file_count += 1
            except OSError:
                continue
    return total_bytes, file_count


def _remember_local_vendor_source(
    source_path: Path,
    *,
    source_kind: str,
    expected_polarities: list[str] | None = None,
    selection_mode: str,
) -> str:
    now = time.time()
    token = uuid.uuid4().hex
    with _vendor_local_selection_lock:
        expired = [
            key for key, value in _vendor_local_selections.items()
            if now - float(value.get("created_at") or 0)
            > _VENDOR_LOCAL_SELECTION_TTL_SECONDS
        ]
        for key in expired:
            _vendor_local_selections.pop(key, None)
        _vendor_local_selections[token] = {
            "source_path": str(source_path),
            "raw_dir": str(source_path),  # compatibility with unexpired V7 tokens
            "raw_name": source_path.name,
            "source_kind": source_kind,
            "expected_polarities": list(expected_polarities or []),
            "selection_mode": selection_mode,
            "created_at": now,
        }
    return token


def _remember_local_vendor_directory(raw_dir: Path) -> str:
    return _remember_local_vendor_source(
        raw_dir,
        source_kind="agilent_d",
        selection_mode="macos_native_folder_picker",
    )


def _resolve_local_vendor_selection(
    token: str,
) -> tuple[str, Path, str, list[str], str]:
    normalized = str(token or "").strip()
    restored_from_ledger = False
    with _vendor_local_selection_lock:
        selection = _vendor_local_selections.get(normalized)
    if normalized:
        # A browser may retry the same token after the V16 process restarted.
        # Only a picker-created path that was already persisted in the fenced
        # conversion ledger can be restored; arbitrary request paths are never
        # accepted here.  Persisted state takes precedence over an expired
        # same-process picker record after lease recovery.
        try:
            ledger_item = _vendor_conversion_ledger().get(normalized)
        except Exception:
            ledger_item = None
        metadata = (
            ledger_item.get("metadata")
            if isinstance(ledger_item, Mapping)
            else None
        )
        if isinstance(metadata, Mapping):
            restored_from_ledger = True
            selection = {
                "source_path": metadata.get("source_path"),
                "raw_dir": metadata.get("source_path"),
                "raw_name": metadata.get("raw_name"),
                "source_kind": metadata.get("source_kind"),
                "expected_polarities": metadata.get("expected_polarities") or [],
                "selection_mode": metadata.get("selection_mode"),
                "created_at": metadata.get("selection_created_at_epoch"),
            }
    if not selection:
        raise HTTPException(
            status_code=404,
            detail=(
                "本机厂商原始数据选择已失效，"
                "请重新选择 Agilent .d 或 Thermo .raw 批次"
            ),
        )
    if (
        not restored_from_ledger
        and time.time() - float(selection.get("created_at") or 0)
        > _VENDOR_LOCAL_SELECTION_TTL_SECONDS
    ):
        with _vendor_local_selection_lock:
            _vendor_local_selections.pop(normalized, None)
        raise HTTPException(
            status_code=410,
            detail=(
                "本机厂商原始数据选择已超过 30 分钟，"
                "请重新选择"
            ),
        )
    source_path = Path(
        str(selection.get("source_path") or selection.get("raw_dir") or "")
    ).expanduser().resolve()
    source_kind = str(selection.get("source_kind") or "agilent_d")
    try:
        if source_kind == "thermo_raw":
            _validate_local_thermo_raw(source_path)
        elif source_kind == "agilent_d":
            _validate_local_agilent_directory(source_path)
        else:
            raise ValueError(f"未知厂商原始数据类型：{source_kind}")
    except (OSError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    expected_polarities = [
        str(mode).upper()
        for mode in (selection.get("expected_polarities") or [])
        if str(mode).upper() in {"POS", "NEG"}
    ]
    selection_mode = str(
        selection.get("selection_mode") or "macos_native_folder_picker"
    )
    return (
        normalized,
        source_path,
        source_kind,
        expected_polarities,
        selection_mode,
    )


@app.post("/api/vendor/converter/prepare")
async def prepare_vendor_converter(
    payload: VendorConverterPrepareRequest,
    request: Request,
):
    """Explicitly download/verify the one-time ProteoWizard image."""

    _require_loopback_request(request)
    from utils.vendor_raw_conversion import (
        ConverterUnavailableError,
        VendorLicenseNotAcceptedError,
        prepare_converter,
    )

    try:
        prepared = await asyncio.to_thread(
            prepare_converter,
            accept_vendor_licenses=payload.accept_vendor_licenses,
            allow_container_pull=payload.allow_container_pull,
        )
        if prepared.get("ready") and isinstance(prepared.get("status"), Mapping):
            with _vendor_converter_preflight_lock:
                _vendor_converter_preflight.update(
                    {
                        "prepared_at_monotonic": time.monotonic(),
                        "status": deepcopy(prepared["status"]),
                    }
                )
        return prepared
    except VendorLicenseNotAcceptedError as exc:
        raise HTTPException(status_code=428, detail=str(exc)) from exc
    except ConverterUnavailableError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc


@app.post("/api/vendor/select-local-folder")
async def select_local_vendor_folder(request: Request):
    """Select one ``.d`` or a parent folder and return short-lived tokens."""

    _require_loopback_request(request)
    try:
        selected_path = await asyncio.to_thread(
            _choose_macos_vendor_directory
        )
        if selected_path is None:
            return {"cancelled": True}
        raw_directories, invalid_directories = await asyncio.to_thread(
            _discover_local_agilent_directories,
            selected_path,
        )
        from utils.vendor_raw_conversion import parse_agilent_acquisition_method

        acquisitions = await asyncio.gather(*[
            asyncio.to_thread(
                parse_agilent_acquisition_method,
                raw_dir,
            )
            for raw_dir in raw_directories
        ])
        selection_mode = (
            "macos_native_parent_folder_picker"
            if selected_path.suffix.lower() != ".d"
            else "macos_native_folder_picker"
        )
        selections: list[dict[str, Any]] = []
        for raw_dir, acquisition in zip(raw_directories, acquisitions):
            centroid_present = bool(
                acquisition.get("centroid_acquisition_present")
            )
            profile_present = bool(
                acquisition.get("profile_acquisition_present")
            )
            if centroid_present and profile_present:
                spectrum_representation = "mixed"
            elif centroid_present:
                spectrum_representation = "centroid"
            elif profile_present:
                spectrum_representation = "profile"
            else:
                spectrum_representation = "unknown"
            token = _remember_local_vendor_source(
                raw_dir,
                source_kind="agilent_d",
                selection_mode=selection_mode,
            )
            selections.append({
                "selection_token": token,
                "raw_name": raw_dir.name,
                "source_kind": "agilent_d",
                "source_format": "Agilent MassHunter .d",
                "acquisition": {
                    "status": acquisition.get("status"),
                    "instrument": acquisition.get("instrument"),
                    "acquisition_mode": acquisition.get("acquisition_mode"),
                    "polarities": acquisition.get("polarities", []),
                    "storage_modes": acquisition.get("storage_modes", []),
                    "centroid_acquisition_present": centroid_present,
                    "profile_acquisition_present": profile_present,
                },
                "spectrum_representation": spectrum_representation,
                "recommended_peak_picking_strategy": (
                    "preserve"
                    if spectrum_representation == "centroid"
                    else "auto"
                ),
            })

        response: dict[str, Any] = {
            "cancelled": False,
            "selection_mode": selection_mode,
            "selected_root_name": selected_path.name,
            "selection_count": len(selections),
            "selections": selections,
            "invalid_directories": invalid_directories,
            "expires_in_seconds": _VENDOR_LOCAL_SELECTION_TTL_SECONDS,
            "browser_upload_required": False,
        }
        # Preserve the original single-directory response shape for older
        # browser bundles while exposing the new ``selections`` list to V7.
        if len(selections) == 1:
            response.update(selections[0])
        return response
    except (OSError, RuntimeError, ValueError, subprocess.TimeoutExpired) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/vendor/select-local-thermo-folder")
async def select_local_thermo_folder(request: Request):
    """Select a Thermo RAW batch without browser-uploading the binary files."""

    _require_loopback_request(request)
    try:
        selected_path = await asyncio.to_thread(_choose_macos_thermo_directory)
        if selected_path is None:
            return {"cancelled": True}
        discovered, invalid_files = await asyncio.to_thread(
            _discover_local_thermo_raw_files,
            selected_path,
        )
        selection_mode = "macos_native_thermo_parent_folder_picker"
        selections: list[dict[str, Any]] = []
        for raw_file, expected_polarities in discovered:
            token = _remember_local_vendor_source(
                raw_file,
                source_kind="thermo_raw",
                expected_polarities=expected_polarities,
                selection_mode=selection_mode,
            )
            selections.append({
                "selection_token": token,
                "raw_name": raw_file.name,
                "source_kind": "thermo_raw",
                "source_format": "Thermo Fisher RAW",
                "relative_source": (
                    raw_file.name
                    if selected_path.is_file()
                    else str(raw_file.relative_to(selected_path))
                ),
                "expected_polarities": expected_polarities,
                "acquisition": {
                    "status": "verified_after_vendor_conversion",
                    "instrument": "待转换后读取精确型号",
                    "polarities": expected_polarities,
                    "storage_modes": [],
                    "centroid_acquisition_present": None,
                    "profile_acquisition_present": None,
                },
                "spectrum_representation": "unknown",
                "recommended_peak_picking_strategy": "auto",
            })
        response: dict[str, Any] = {
            "cancelled": False,
            "selection_mode": selection_mode,
            "selected_root_name": selected_path.name,
            "selection_count": len(selections),
            "selections": selections,
            "invalid_files": invalid_files,
            "expires_in_seconds": _VENDOR_LOCAL_SELECTION_TTL_SECONDS,
            "browser_upload_required": False,
            "discovery_scope": "selected folder plus one direct subfolder level",
        }
        if len(selections) == 1:
            response.update(selections[0])
        return response
    except (OSError, RuntimeError, ValueError, subprocess.TimeoutExpired) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


def _managed_vendor_destination_path(value: Any, *, intake_role: str) -> Path:
    """Resolve a direct managed output path without following its leaf."""

    requested_destination = Path(str(value or "")).expanduser()
    if not requested_destination.is_absolute():
        raise ValueError("Validated vendor conversion output must be absolute")
    if requested_destination.is_symlink():
        raise ValueError("Validated vendor conversion output must not be a symlink")
    managed_root = (
        BLANK_DIR if intake_role == "blank" else UPLOAD_DIR
    ).expanduser().resolve()
    destination = (
        requested_destination.parent.resolve() / requested_destination.name
    )
    if destination.parent != managed_root:
        raise ValueError(
            "Validated vendor conversion output is outside its managed directory"
        )
    if destination.is_symlink():
        raise ValueError("Validated vendor conversion output must not be a symlink")
    if destination.suffix.lower() != ".mzml":
        raise ValueError("Validated vendor conversion output must end in .mzML")
    return destination


def _register_validated_vendor_conversion(
    details: Mapping[str, Any],
    *,
    registration_token: str,
) -> dict[str, Any]:
    """Idempotently register one ledger-validated converted file."""

    conversion = deepcopy(details.get("conversion") or {})
    selection_audit = deepcopy(details.get("selection_audit") or {})
    display_name = str(details.get("display_name") or "")
    raw_name = str(details.get("raw_name") or "")
    source_kind = str(details.get("source_kind") or "")
    intake_role = str(details.get("intake_role") or "sample")
    blank_role = str(details.get("blank_role") or "background_blank")
    selection_mode = str(details.get("selection_mode") or "")
    destination = _managed_vendor_destination_path(
        details.get("destination"),
        intake_role=intake_role,
    )
    output = conversion.get("output") or {}
    try:
        before = destination.lstat()
    except OSError as exc:
        raise ValueError("Validated vendor conversion output is missing") from exc
    if not stat.S_ISREG(before.st_mode):
        raise ValueError("Validated vendor conversion output is missing")
    expected_size = int(output.get("size_bytes") or -1)
    if expected_size < 0 or before.st_size != expected_size:
        raise ValueError("Validated vendor conversion output size changed")
    expected_sha256 = _valid_content_sha256(output.get("sha256"))
    if not expected_sha256:
        raise ValueError("Validated vendor conversion output SHA-256 is missing")
    digest = hashlib.sha256()
    try:
        descriptor = os.open(destination, os.O_RDONLY | os.O_NOFOLLOW)
    except OSError as exc:
        raise ValueError(
            "Validated vendor conversion output cannot be opened safely"
        ) from exc
    try:
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or (opened.st_dev, opened.st_ino) != (before.st_dev, before.st_ino)
        ):
            raise ValueError("Validated vendor conversion output changed before hashing")
        with os.fdopen(descriptor, "rb", closefd=False) as handle:
            for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                digest.update(chunk)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if (
        (opened.st_dev, opened.st_ino, opened.st_size, opened.st_mtime_ns)
        != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
        or digest.hexdigest() != expected_sha256
    ):
        raise ValueError("Validated vendor conversion output SHA-256 changed")
    source = conversion.get("source") or {}
    acquisition = source.get("acquisition") or {}

    if intake_role == "blank":
        blank_record = _register_managed_blank_file(
            destination,
            display_name,
            blank_role=blank_role,
            extra_metadata={
                "source_raw_filename": raw_name,
                "vendor_conversion": conversion,
                "vendor_local_selection": selection_audit,
            },
        )
        logger.info(
            "Registered local %s blank %s → %s (%d MS1 + %d MS2)",
            source_kind,
            raw_name,
            display_name,
            output.get("ms1_count", 0),
            output.get("ms2_count", 0),
        )
        return {
            "uploaded": [],
            "blank_uploads": [blank_record],
            "count": 0,
            "blank_count": 1,
            "failed": [],
            "failed_count": 0,
            "conversion": conversion,
            "selection_mode": selection_mode,
        }

    stable_task_id = hashlib.sha256(
        f"v16-vendor-registration:{registration_token}".encode("utf-8")
    ).hexdigest()[:12]
    task_id = _register_pending_task(
        display_name,
        destination,
        expected_size,
        task_id=stable_task_id,
        extra_metadata={
            "source_raw_filename": raw_name,
            "vendor_conversion": conversion,
            "vendor_local_selection": selection_audit,
            "recommended_isolation_window_width_da": acquisition.get(
                "isolation_width_da"
            ),
            "recommended_isolation_window_width_source": acquisition.get(
                "isolation_width_source"
            ),
        },
    )
    if _save_tasks() is False:
        raise OSError("Could not persist the registered vendor conversion task")
    logger.info(
        "Registered local %s %s → %s (%d MS1 + %d MS2) → task %s",
        source_kind,
        raw_name,
        display_name,
        output.get("ms1_count", 0),
        output.get("ms2_count", 0),
        task_id,
    )
    return {
        "uploaded": [{
            "task_id": task_id,
            "filename": display_name,
            "source_raw_filename": raw_name,
            "size_mb": round(float(expected_size) / (1024 * 1024), 2),
            "ms1_count": output.get("ms1_count", 0),
            "ms2_count": output.get("ms2_count", 0),
            "isolation_width_da": acquisition.get("isolation_width_da"),
            "conversion_warnings": output.get("warnings", []),
            "peak_picking_strategy": (
                conversion.get("conversion_policy") or {}
            ).get("peak_picking_strategy_resolved"),
        }],
        "blank_uploads": [],
        "count": 1,
        "failed": [],
        "failed_count": 0,
        "conversion": conversion,
        "selection_mode": selection_mode,
    }


@app.post("/api/vendor/convert-local")
async def convert_local_vendor_folder(
    payload: VendorLocalConvertRequest,
    request: Request,
):
    """Convert a native-picker-selected Agilent .d or Thermo .raw locally."""

    _require_loopback_request(request)
    from utils.vendor_conversion_cache import (
        CacheEntryInvalidError,
        LedgerStateError,
        SourceMaterialChangedError,
        VendorConversionCache,
        build_conversion_cache_key,
    )
    from utils.vendor_raw_conversion import (
        ConverterUnavailableError,
        VendorLicenseNotAcceptedError,
        VendorRawConversionError,
        VendorRawValidationError,
        convert_agilent_d,
        convert_thermo_raw,
        inspect_agilent_d,
        inspect_thermo_raw,
        resolve_peak_picking_strategy,
    )

    normalized_token = str(payload.selection_token or "").strip()
    ledger = _vendor_conversion_ledger()
    existing_item = ledger.get(normalized_token)
    if existing_item and existing_item.get("state") == "registered":
        response = (existing_item.get("details") or {}).get("response")
        if isinstance(response, Mapping):
            return deepcopy(dict(response))
        raise HTTPException(
            status_code=409,
            detail="Vendor conversion is registered but its response is unavailable",
        )
    if existing_item and existing_item.get("state") == "validated":
        validated_details = existing_item.get("details") or {}
        registration_token = str(existing_item.get("registration_token") or "")
        claim_token = str(existing_item.get("claim_token") or "")
        registration_marker = (
            f"registration:{normalized_token}:{registration_token}"
        )
        try:
            _begin_vendor_registration_operation(registration_marker)
        except RuntimeError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        try:
            try:
                response, cancellation = await _await_thread_completion(
                    _register_validated_vendor_conversion,
                    validated_details,
                    registration_token=registration_token,
                )
                ledger.mark_registered(
                    normalized_token,
                    claim_token,
                    registration_token,
                    details={"response": response},
                )
                with _vendor_local_selection_lock:
                    _vendor_local_selections.pop(normalized_token, None)
                if cancellation is not None:
                    raise cancellation
                return response
            except (OSError, ValueError, LedgerStateError) as exc:
                raise HTTPException(status_code=409, detail=str(exc)) from exc
        finally:
            _end_vendor_registration_operation(registration_marker)
    if existing_item and existing_item.get("state") == "converting":
        ledger.recover_interrupted()
        existing_item = ledger.get(normalized_token)
        if existing_item and existing_item.get("state") == "converting":
            raise HTTPException(
                status_code=409,
                detail="This vendor conversion is already running",
            )

    (
        token,
        raw_source,
        source_kind,
        expected_polarities,
        selection_mode,
    ) = _resolve_local_vendor_selection(normalized_token)
    raw_name = raw_source.name
    polarity_suffix = (
        f"_{expected_polarities[0]}"
        if len(expected_polarities) == 1
        and not raw_source.stem.upper().endswith(
            (f"_{expected_polarities[0]}", f"-{expected_polarities[0]}")
        )
        else ""
    )
    display_name = f"{raw_source.stem}{polarity_suffix}.mzML"
    storage_name = f"{token}_{display_name}"
    destination_root = BLANK_DIR if payload.intake_role == "blank" else UPLOAD_DIR
    destination = destination_root / storage_name
    request_contract = {
        "backend": payload.backend,
        "intake_role": payload.intake_role,
        "blank_role": payload.blank_role,
        "peak_picking_strategy": payload.peak_picking_strategy,
    }
    old_metadata = (
        existing_item.get("metadata")
        if isinstance(existing_item, Mapping)
        and isinstance(existing_item.get("metadata"), Mapping)
        else {}
    )
    if existing_item and existing_item.get("state") == "queued":
        recovered_destination = old_metadata.get("destination")
        if recovered_destination:
            recovered_name = Path(str(recovered_destination)).name
            if not recovered_name.startswith(f"{token}_"):
                raise HTTPException(
                    status_code=409,
                    detail=(
                        "Recovered vendor output uses a legacy non-unique path; "
                        "select the source again before retrying"
                    ),
                )
            destination = _managed_vendor_destination_path(
                recovered_destination,
                intake_role=payload.intake_role,
            )
    ledger_metadata = {
        "source_path": str(raw_source),
        "raw_name": raw_name,
        "source_kind": source_kind,
        "expected_polarities": expected_polarities,
        "selection_mode": selection_mode,
        "selection_created_at_epoch": old_metadata.get(
            "selection_created_at_epoch",
            time.time(),
        ),
        "display_name": display_name,
        "destination": str(destination),
        "request_contract": request_contract,
    }
    if existing_item and existing_item.get("state") == "queued":
        if (
            str(old_metadata.get("source_path") or "") != str(raw_source)
            or old_metadata.get("request_contract") != request_contract
            or str(old_metadata.get("destination") or "") != str(destination)
        ):
            raise HTTPException(
                status_code=409,
                detail="Recovered vendor conversion request does not match its ledger",
            )
    destination = _managed_vendor_destination_path(
        destination,
        intake_role=payload.intake_role,
    )
    claim_token = ""
    owns_destination = False
    heartbeat_stop = threading.Event()
    heartbeat_errors: list[Exception] = []
    heartbeat_thread: threading.Thread | None = None
    staging_destination: Path | None = None
    try:
        def run_bounded_conversion() -> tuple[dict[str, Any], int, int]:
            nonlocal claim_token, owns_destination
            nonlocal heartbeat_thread, staging_destination
            with _vendor_conversion_slots:
                with _vendor_conversion_active_claims_lock:
                    if not _vendor_conversion_accepting_claims:
                        raise LedgerStateError(
                            "Vendor conversion service is shutting down"
                        )
                    claim = ledger.claim(
                        token,
                        metadata=ledger_metadata,
                        details={"destination": str(destination)},
                        lease_seconds=_VENDOR_CONVERSION_LEASE_SECONDS,
                    )
                    claim_token = str(claim.get("claim_token") or "")
                    _vendor_conversion_active_claims.add(claim_token)
                owns_destination = True
                staging_destination = _managed_vendor_destination_path(
                    destination.parent
                    / f".{token}.{claim_token}.part.mzML",
                    intake_role=payload.intake_role,
                )
                if staging_destination.exists() or staging_destination.is_symlink():
                    if staging_destination.is_symlink():
                        raise ValueError(
                            "Managed vendor conversion staging path is a symlink"
                        )
                    if not staging_destination.is_file():
                        raise ValueError(
                            "Managed vendor conversion staging path is not a file"
                        )
                    staging_destination.unlink()

                def renew_claim() -> None:
                    while not heartbeat_stop.wait(
                        _VENDOR_CONVERSION_HEARTBEAT_SECONDS
                    ):
                        try:
                            ledger.renew_lease(
                                token,
                                claim_token,
                                lease_seconds=_VENDOR_CONVERSION_LEASE_SECONDS,
                            )
                        except Exception as exc:
                            heartbeat_errors.append(exc)
                            heartbeat_stop.set()
                            return

                heartbeat_thread = threading.Thread(
                    target=renew_claim,
                    name=f"vendor-lease-{token[:12]}",
                    daemon=True,
                )
                heartbeat_thread.start()

                def publish_conversion(
                    conversion: dict[str, Any],
                    total_bytes: int,
                    file_count: int,
                ) -> tuple[dict[str, Any], int, int]:
                    if heartbeat_errors:
                        raise LedgerStateError(
                            "Vendor conversion lease renewal failed: "
                            f"{heartbeat_errors[0]}"
                        )
                    ledger.renew_lease(
                        token,
                        claim_token,
                        lease_seconds=_VENDOR_CONVERSION_LEASE_SECONDS,
                    )
                    if (
                        staging_destination is None
                        or staging_destination.is_symlink()
                        or not staging_destination.is_file()
                    ):
                        raise ValueError(
                            "Validated vendor conversion staging output is missing"
                        )
                    if destination.is_symlink():
                        raise ValueError(
                            "Managed vendor conversion destination is a symlink"
                        )
                    if destination.exists() and not destination.is_file():
                        raise ValueError(
                            "Managed vendor conversion destination is not a file"
                        )
                    os.replace(staging_destination, destination)
                    directory_descriptor = os.open(destination.parent, os.O_RDONLY)
                    try:
                        os.fsync(directory_descriptor)
                    finally:
                        os.close(directory_descriptor)
                    output = conversion.get("output")
                    if isinstance(output, dict):
                        output["path"] = str(destination)
                    return conversion, total_bytes, file_count

                raw_report = (
                    inspect_thermo_raw(
                        raw_source,
                        expected_polarities=expected_polarities,
                    )
                    if source_kind == "thermo_raw"
                    else inspect_agilent_d(raw_source)
                )
                total_bytes = int(raw_report.get("total_bytes") or 0)
                file_count = int(raw_report.get("file_count") or 0)
                if total_bytes <= 0:
                    raise ValueError("所选厂商原始数据为空")
                if total_bytes > MAX_VENDOR_RAW_SIZE_BYTES:
                    raise ValueError(
                        f"厂商原始数据超过 {MAX_VENDOR_RAW_SIZE_MB} MB 本机转换上限"
                    )
                status = _cached_or_prepare_vendor_converter_status(
                    accept_vendor_licenses=payload.accept_vendor_licenses,
                )
                selected_backend, converter_identity = _vendor_converter_identity(
                    status,
                    requested_backend=payload.backend,
                    accept_vendor_licenses=payload.accept_vendor_licenses,
                )
                resolved_peak_picking = resolve_peak_picking_strategy(
                    raw_report,
                    payload.peak_picking_strategy,
                )
                cache_key = build_conversion_cache_key(
                    raw_report,
                    converter_identity=converter_identity,
                    resolved_peak_picking=resolved_peak_picking,
                    validation_rules_version=(
                        VENDOR_CONVERSION_VALIDATION_RULES_VERSION
                    ),
                )
                cache = VendorConversionCache(VENDOR_CONVERSION_CACHE_ROOT)
                with _vendor_conversion_key_lock(cache_key.digest):
                    cache_started = time.perf_counter()
                    cache_hit = cache.lookup(cache_key)
                    if cache_hit is not None:
                        cached_size = int(
                            (cache_hit.manifest.get("payload") or {}).get(
                                "size_bytes"
                            )
                            or 0
                        )
                        _require_disk_budget(
                            destination_root,
                            max(cached_size, 64 * 1024 * 1024),
                            "cached vendor mzML materialization",
                        )
                        materialization = cache.materialize(
                            cache_key,
                            staging_destination,
                        )
                        conversion = _vendor_conversion_from_cache(
                            raw_report=raw_report,
                            cache_hit=cache_hit,
                            materialization=materialization,
                            destination=staging_destination,
                            peak_picking_strategy_requested=(
                                payload.peak_picking_strategy
                            ),
                            elapsed_seconds=time.perf_counter() - cache_started,
                        )
                        return publish_conversion(
                            conversion,
                            total_bytes,
                            file_count,
                        )

                    _require_disk_budget(
                        destination_root,
                        max(int(total_bytes * 4.0), 256 * 1024 * 1024),
                        "vendor raw to mzML conversion",
                    )
                    converter = (
                        convert_thermo_raw
                        if source_kind == "thermo_raw"
                        else convert_agilent_d
                    )
                    converter_kwargs: dict[str, Any] = {
                        "backend": selected_backend,
                        "accept_vendor_licenses": payload.accept_vendor_licenses,
                        "allow_container_pull": False,
                        "peak_picking_strategy": payload.peak_picking_strategy,
                        "converter_status_report": status,
                        "raw_report": raw_report,
                    }
                    if source_kind == "thermo_raw":
                        converter_kwargs["expected_polarities"] = expected_polarities
                    conversion = converter(
                        raw_source,
                        staging_destination,
                        **converter_kwargs,
                    )
                    try:
                        stored = cache.store(
                            cache_key,
                            staging_destination,
                            conversion,
                            validation_rules_version=(
                                VENDOR_CONVERSION_VALIDATION_RULES_VERSION
                            ),
                        )
                    except OSError as exc:
                        logger.warning(
                            "Vendor conversion cache store failed; using validated output: %s",
                            exc,
                        )
                        conversion["cache"] = {
                            "status": "store_failed",
                            "hit": False,
                            "cache_key": cache_key.digest,
                            "cache_root": str(VENDOR_CONVERSION_CACHE_ROOT),
                            "error": str(exc),
                        }
                    else:
                        conversion["cache"] = {
                            "status": "stored",
                            "hit": False,
                            "cache_key": cache_key.digest,
                            "cache_root": str(VENDOR_CONVERSION_CACHE_ROOT),
                            "payload_sha256": stored.manifest["payload"]["sha256"],
                        }
                    return publish_conversion(
                        conversion,
                        total_bytes,
                        file_count,
                    )

        worker_result, cancellation = await _await_thread_completion(
            run_bounded_conversion
        )
        conversion, total_bytes, file_count = worker_result
        heartbeat_stop.set()
        if heartbeat_thread is not None:
            heartbeat_thread.join(timeout=5.0)
        if heartbeat_errors:
            raise LedgerStateError(
                "Vendor conversion lease renewal failed: "
                f"{heartbeat_errors[0]}"
            )
        source = conversion.get("source") or {}
        source["source_path"] = raw_name
        acquisition = source.get("acquisition") or {}
        output = conversion.get("output") or {}
        instrument = output.get("instrument") or {}
        models = instrument.get("models") or []
        model_names = [
            str(model.get("name") or "")
            for model in models
            if str(model.get("name") or "")
        ]
        if model_names:
            acquisition["instrument"] = model_names[0]
            acquisition["instrument_models"] = model_names
        observed_polarities = [
            mode
            for mode, count in (output.get("polarities") or {}).items()
            if int(count or 0) > 0
        ]
        if observed_polarities:
            acquisition["polarities_observed"] = observed_polarities
        acquisition["ms1_observed"] = int(output.get("ms1_count") or 0) > 0
        acquisition["ms2_observed"] = int(output.get("ms2_count") or 0) > 0
        converted_acquisition = output.get("acquisition") or {}
        isolation_widths = converted_acquisition.get("isolation_widths_da") or []
        if len(isolation_widths) == 1 and acquisition.get("isolation_width_da") is None:
            acquisition["isolation_width_da"] = isolation_widths[0]
            acquisition["isolation_width_source"] = (
                "converted mzML isolation lower + upper offsets"
            )
        source["acquisition"] = acquisition
        selection_audit = {
            "mode": selection_mode,
            "source_kind": source_kind,
            "expected_polarities": expected_polarities,
            "intake_role": payload.intake_role,
            "file_count": file_count,
            "raw_size_bytes": total_bytes,
        }
        validated_details = {
            "destination": str(destination),
            "display_name": display_name,
            "raw_name": raw_name,
            "source_kind": source_kind,
            "intake_role": payload.intake_role,
            "blank_role": payload.blank_role,
            "selection_mode": selection_mode,
            "selection_audit": selection_audit,
            "conversion": conversion,
        }
        validated_item = ledger.mark_validated(
            token,
            claim_token,
            details=validated_details,
        )
        registration_token = str(validated_item.get("registration_token") or "")
        response, registration_cancellation = await _await_thread_completion(
            _register_validated_vendor_conversion,
            validated_details,
            registration_token=registration_token,
        )
        if cancellation is None:
            cancellation = registration_cancellation
        ledger.mark_registered(
            token,
            claim_token,
            registration_token,
            details={"response": response},
        )
        with _vendor_local_selection_lock:
            _vendor_local_selections.pop(token, None)
        if cancellation is not None:
            raise cancellation
        return response
    except VendorLicenseNotAcceptedError as exc:
        _fail_vendor_claim_if_converting(ledger, token, claim_token, exc)
        raise HTTPException(status_code=428, detail=str(exc)) from exc
    except ConverterUnavailableError as exc:
        _fail_vendor_claim_if_converting(ledger, token, claim_token, exc)
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except (
        CacheEntryInvalidError,
        SourceMaterialChangedError,
        VendorRawValidationError,
        VendorRawConversionError,
        ValueError,
    ) as exc:
        _fail_vendor_claim_if_converting(ledger, token, claim_token, exc)
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except InsufficientStorageError as exc:
        _fail_vendor_claim_if_converting(ledger, token, claim_token, exc)
        raise HTTPException(status_code=507, detail=str(exc)) from exc
    except LedgerStateError as exc:
        _fail_vendor_claim_if_converting(ledger, token, claim_token, exc)
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except Exception as exc:
        _fail_vendor_claim_if_converting(ledger, token, claim_token, exc)
        logger.exception("Local vendor raw conversion failed")
        raise HTTPException(
            status_code=500,
            detail=f"Local vendor raw conversion failed: {exc}",
        ) from exc
    finally:
        heartbeat_stop.set()
        if heartbeat_thread is not None and heartbeat_thread.is_alive():
            heartbeat_thread.join(timeout=5.0)
        release_ledger_lock = False
        if claim_token:
            with _vendor_conversion_active_claims_lock:
                _vendor_conversion_active_claims.discard(claim_token)
                release_ledger_lock = bool(
                    not _vendor_conversion_accepting_claims
                    and not _vendor_conversion_active_claims
                )
        if release_ledger_lock:
            _release_vendor_conversion_service_lock()
        if owns_destination and staging_destination is not None:
            try:
                if not staging_destination.is_symlink():
                    staging_destination.unlink(missing_ok=True)
            except OSError:
                pass


@app.get("/api/tasks/{task_id}/input-audit")
async def task_input_audit(task_id: str):
    """Expose raw-input provenance and MS1/MS2 preservation in one verdict."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    conversion = task.get("vendor_conversion")
    if not isinstance(conversion, dict):
        return {
            "task_id": task_id,
            "filename": task.get("filename"),
            "input_kind": "open_format",
            "vendor_conversion": False,
            "verdict": "not_applicable",
            "note": "The task was not created by managed vendor-raw conversion.",
        }

    source = conversion.get("source") or {}
    acquisition = source.get("acquisition") or {}
    output = conversion.get("output") or {}
    expected_ms1 = acquisition.get("ms1_expected")
    expected_ms2 = acquisition.get("ms2_expected")
    ms1_count = int(output.get("ms1_count") or 0)
    ms2_count = int(output.get("ms2_count") or 0)
    checksum_matches = output.get("source_checksum_matches") or {}
    critical_hashes_match = bool(checksum_matches) and all(
        value is True for value in checksum_matches.values()
    )
    preserves_expected_levels = (
        (expected_ms1 is not True or ms1_count > 0)
        and (expected_ms2 is not True or ms2_count > 0)
    )
    valid = (
        output.get("status") == "valid"
        and preserves_expected_levels
        and critical_hashes_match
        and not output.get("errors")
    )
    return {
        "task_id": task_id,
        "filename": task.get("filename"),
        "source_raw_filename": task.get("source_raw_filename"),
        "input_kind": source.get("format") or conversion.get("source_kind") or "vendor_raw",
        "vendor_conversion": True,
        "verdict": "verified" if valid else "failed_closed",
        "preserves_expected_ms_levels": preserves_expected_levels,
        "paired_raw_provenance": critical_hashes_match,
        "conversion_policy": conversion.get("conversion_policy") or {},
        "analysis_representation_processing": task.get(
            "spectrum_representation_processing"
        ) or {},
        "acquisition": {
            "instrument": (
                acquisition.get("instrument")
                or ((output.get("instrument") or {}).get("models") or [{}])[0].get(
                    "name"
                )
            ),
            "instrument_metadata": output.get("instrument") or {},
            "mode": acquisition.get("acquisition_mode"),
            "polarities_expected": acquisition.get("polarities") or [],
            "polarities_observed": acquisition.get("polarities_observed") or [],
            "ms1_expected": expected_ms1,
            "ms2_expected": expected_ms2,
            "ms1_observed": acquisition.get("ms1_observed"),
            "ms2_observed": acquisition.get("ms2_observed"),
            "isolation_width_da": acquisition.get("isolation_width_da"),
            "collision_energies": acquisition.get("collision_energies") or [],
            "converted_cv_audit": output.get("acquisition") or {},
        },
        "converted_mzml": {
            "status": output.get("status"),
            "sha256": output.get("sha256"),
            "actual_spectrum_count": output.get("actual_spectrum_count"),
            "mass_spectrum_count": output.get("mass_spectrum_count"),
            "auxiliary_spectrum_count": output.get("auxiliary_spectrum_count"),
            "ms1_count": ms1_count,
            "nonempty_ms1_count": output.get("nonempty_ms1_count"),
            "ms2_count": ms2_count,
            "nonempty_ms2_count": output.get("nonempty_ms2_count"),
            "precursor_count": output.get("precursor_count"),
            "polarities": output.get("polarities") or {},
            "representation_by_ms_level": output.get(
                "representation_by_ms_level"
            ) or {},
            "source_checksum_matches": checksum_matches,
            "warnings": output.get("warnings") or [],
            "errors": output.get("errors") or [],
        },
        "converter": conversion.get("converter") or {},
    }


@app.post("/api/vendor/convert")
async def convert_vendor_raw_upload(
    files: list[UploadFile] = File(...),
    manifest: str = Form(...),
    accept_vendor_licenses: bool = Form(False),
    allow_container_pull: bool = Form(False),
    backend: str = Form("auto"),
    peak_picking_strategy: Literal["auto", "preserve", "vendor", "openms"] = Form(
        "auto"
    ),
):
    """Upload one Agilent ``.d`` directory, convert it, and create a task.

    Browsers upload directories as individual multipart files.  ``manifest``
    preserves each ``webkitRelativePath`` so the proprietary directory layout
    is reconstructed only inside the managed staging root.
    """
    from utils.vendor_raw_conversion import (
        ConverterUnavailableError,
        VendorLicenseNotAcceptedError,
        VendorRawConversionError,
        VendorRawValidationError,
        convert_agilent_d,
    )

    if not files:
        raise HTTPException(status_code=400, detail="No .d directory files provided")
    if len(files) > MAX_VENDOR_RAW_FILES:
        raise HTTPException(
            status_code=400,
            detail=f"Agilent .d upload exceeds {MAX_VENDOR_RAW_FILES} files",
        )
    try:
        raw_paths = json.loads(manifest)
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=400, detail="Invalid .d upload manifest") from exc
    if not isinstance(raw_paths, list) or len(raw_paths) != len(files):
        raise HTTPException(
            status_code=400,
            detail="The .d upload manifest does not match the uploaded files",
        )

    try:
        relative_paths = [_safe_vendor_relative_path(value) for value in raw_paths]
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    if len(set(relative_paths)) != len(relative_paths):
        raise HTTPException(status_code=400, detail="Duplicate .d member path")
    raw_roots = {path.parts[0] for path in relative_paths}
    if len(raw_roots) != 1:
        raise HTTPException(
            status_code=400,
            detail="Select exactly one top-level Agilent .d directory",
        )

    job_id = uuid.uuid4().hex[:12]
    VENDOR_RAW_DIR.mkdir(parents=True, exist_ok=True)
    staging_root = VENDOR_RAW_DIR / job_id
    staging_root.mkdir(parents=True, exist_ok=False)
    destination: Path | None = None
    task_id: str | None = None
    total_bytes = 0
    try:
        _require_disk_budget(
            VENDOR_RAW_DIR,
            min(MAX_VENDOR_RAW_SIZE_BYTES, 256 * 1024 * 1024),
            "Agilent .d staging upload",
        )
        for upload, relative in zip(files, relative_paths):
            target = staging_root.joinpath(*relative.parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("wb") as handle:
                while True:
                    chunk = await upload.read(4 * 1024 * 1024)
                    if not chunk:
                        break
                    total_bytes += len(chunk)
                    if total_bytes > MAX_VENDOR_RAW_SIZE_BYTES:
                        raise ValueError(
                            f"Agilent .d upload exceeds {MAX_VENDOR_RAW_SIZE_MB} MB"
                        )
                    if (
                        total_bytes == len(chunk)
                        or total_bytes % (64 * 1024 * 1024) < len(chunk)
                    ):
                        require_storage_capacity(
                            VENDOR_RAW_DIR,
                            len(chunk),
                            "Agilent .d staging upload",
                        )
                    handle.write(chunk)
        if total_bytes == 0:
            raise ValueError("Agilent .d upload is empty")

        raw_name = next(iter(raw_roots))
        raw_dir = staging_root / raw_name
        display_name = f"{Path(raw_name).stem}.mzML"
        storage_name = f"{uuid.uuid4().hex[:8]}_{display_name}"
        destination = UPLOAD_DIR / storage_name
        _require_disk_budget(
            UPLOAD_DIR,
            max(int(total_bytes * 2.5), 256 * 1024 * 1024),
            "Agilent .d to mzML conversion",
        )
        conversion = await asyncio.to_thread(
            convert_agilent_d,
            raw_dir,
            destination,
            backend=backend,
            accept_vendor_licenses=accept_vendor_licenses,
            # Never turn an analysis/upload request into a network operation.
            allow_container_pull=False,
            peak_picking_strategy=peak_picking_strategy,
        )
        # The managed staging path is not useful provenance after cleanup.
        # Raw name + critical-file hashes are sufficient to reproduce identity.
        source = conversion.get("source") or {}
        source["source_path"] = raw_name
        acquisition = source.get("acquisition") or {}
        output = conversion.get("output") or {}
        task_id = _register_pending_task(
            display_name,
            destination,
            int(output.get("size_bytes") or destination.stat().st_size),
            extra_metadata={
                "source_raw_filename": raw_name,
                "vendor_conversion": conversion,
                "recommended_isolation_window_width_da": acquisition.get(
                    "isolation_width_da"
                ),
                "recommended_isolation_window_width_source": acquisition.get(
                    "isolation_width_source"
                ),
            },
        )
        _save_tasks()
        logger.info(
            "Converted Agilent raw %s → %s (%d MS1 + %d MS2) → task %s",
            raw_name,
            display_name,
            output.get("ms1_count", 0),
            output.get("ms2_count", 0),
            task_id,
        )
        return {
            "uploaded": [
                {
                    "task_id": task_id,
                    "filename": display_name,
                    "source_raw_filename": raw_name,
                    "size_mb": round(
                        float(output.get("size_bytes") or 0) / (1024 * 1024), 2
                    ),
                    "ms1_count": output.get("ms1_count", 0),
                    "ms2_count": output.get("ms2_count", 0),
                    "isolation_width_da": acquisition.get("isolation_width_da"),
                    "conversion_warnings": output.get("warnings", []),
                    "peak_picking_strategy": (
                        conversion.get("conversion_policy") or {}
                    ).get("peak_picking_strategy_resolved"),
                }
            ],
            "count": 1,
            "failed": [],
            "failed_count": 0,
            "conversion": conversion,
        }
    except HTTPException:
        raise
    except VendorLicenseNotAcceptedError as exc:
        raise HTTPException(status_code=428, detail=str(exc)) from exc
    except ConverterUnavailableError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except (VendorRawValidationError, VendorRawConversionError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except InsufficientStorageError as exc:
        raise HTTPException(status_code=507, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Agilent .d conversion failed")
        raise HTTPException(
            status_code=500,
            detail=f"Agilent .d conversion failed: {exc}",
        ) from exc
    finally:
        for upload in files:
            try:
                await upload.close()
            except Exception:
                pass
        shutil.rmtree(staging_root, ignore_errors=True)
        if task_id is None and destination is not None:
            try:
                destination.unlink(missing_ok=True)
            except OSError:
                pass


@app.post("/api/upload")
async def upload_files(files: list[UploadFile] = File(...)):
    """Upload samples while isolating failures to the affected file."""
    if not files:
        raise HTTPException(status_code=400, detail="No files provided")
    if len(files) > MAX_BATCH_UPLOAD_FILES:
        raise HTTPException(
            status_code=400,
            detail=f"At most {MAX_BATCH_UPLOAD_FILES} sample files per upload",
        )

    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    uploaded, failed = [], []
    batch_total_bytes = 0
    _require_disk_budget(UPLOAD_DIR, MAX_BATCH_UPLOAD_SIZE_BYTES, "sample upload")
    for uf in files:
        display_name = _safe_sample_filename(uf.filename or "")
        dest: Optional[Path] = None
        partial: Optional[Path] = None
        task_id: Optional[str] = None
        try:
            if not display_name:
                raise ValueError("Missing filename")

            ext = Path(display_name).suffix.lower()
            if ext not in SAMPLE_UPLOAD_EXTENSIONS:
                raise ValueError(
                    f"Unsupported format: {ext or '(none)'}. "
                    "Allowed: .mgf, .mzML, .mzXML, .html (报告)"
                )

            # 先写 .part，验证通过后再原子改名。失败文件不会伪装成可分析样本。
            storage_name = f"{uuid.uuid4().hex[:8]}_{display_name}"
            dest = UPLOAD_DIR / storage_name
            # 保留最终扩展名，HTML 报告验证器会先校验 .html/.htm 后缀。
            partial = UPLOAD_DIR / f".{uuid.uuid4().hex[:8]}.part_{storage_name}"
            total_bytes = 0
            chunk_size = 4 * 1024 * 1024
            with open(partial, "wb") as handle:
                while True:
                    chunk = await uf.read(chunk_size)
                    if not chunk:
                        break
                    total_bytes += len(chunk)
                    batch_total_bytes += len(chunk)
                    if total_bytes > MAX_UPLOAD_SIZE_BYTES:
                        raise ValueError(f"File exceeds {MAX_UPLOAD_SIZE_MB} MB limit")
                    if batch_total_bytes > MAX_BATCH_UPLOAD_SIZE_BYTES:
                        raise ValueError(
                            f"Batch exceeds {MAX_BATCH_UPLOAD_SIZE_MB} MB total limit"
                        )
                    if (
                        total_bytes == len(chunk)
                        or total_bytes % (64 * 1024 * 1024) < len(chunk)
                    ):
                        try:
                            require_storage_capacity(
                                UPLOAD_DIR, len(chunk), "sample upload"
                            )
                        except InsufficientStorageError as exc:
                            raise ValueError(str(exc)) from exc
                    handle.write(chunk)

            if total_bytes == 0:
                raise ValueError("Empty file")

            # HTML 是可执行文档内容，不把任意网页当成样本。
            if ext in (".html", ".htm"):
                from utils.report_html_io import is_raw_report_html
                if not is_raw_report_html(partial):
                    raise ValueError(
                        "Not a parseable MS2 report; HTML must contain "
                        "precursor m/z and an MS2 fragment table"
                    )

            partial.replace(dest)
            partial = None
            file_size_mb = total_bytes / (1024 * 1024)

            task_id = _register_pending_task(
                display_name, dest, total_bytes
            )

            uploaded.append({
                "task_id": task_id,
                "filename": display_name,
                "size_mb": round(file_size_mb, 2),
            })
            logger.info(
                "Uploaded: %s (%.2f MB) → task %s",
                display_name, file_size_mb, task_id,
            )
        except Exception as exc:
            if task_id:
                with _tasks_lock:
                    _tasks.pop(task_id, None)
            for path in (partial, dest):
                if path is not None:
                    try:
                        path.unlink(missing_ok=True)
                    except OSError:
                        pass
            error = str(exc) or exc.__class__.__name__
            failed.append({"filename": display_name or "(unnamed)", "error": error})
            logger.warning(
                "Upload rejected for %s: %s", display_name or "(unnamed)", error,
            )
        finally:
            try:
                await uf.close()
            except Exception:
                pass

    # Save tasks index once after all files are processed (not per-file).
    if uploaded:
        _save_tasks()
    if not uploaded:
        summary = "; ".join(
            f"{item['filename']}: {item['error']}" for item in failed[:5]
        )
        raise HTTPException(status_code=400, detail=summary or "No valid files uploaded")
    return {
        "uploaded": uploaded,
        "count": len(uploaded),
        "failed": failed,
        "failed_count": len(failed),
    }


@app.post("/api/blank/upload")
async def upload_blank(file: UploadFile = File(...)):
    """v3 Phase 1: 上传空白样本 (.mgf/.mzML/.mzXML), 返回服务器路径供 AnalysisParams.blank_file 使用。

    纯增量端点 —— 不创建 task、不影响 /api/upload。立即解析验证可加载并报峰数。
    """
    if not file or not file.filename:
        raise HTTPException(status_code=400, detail="No blank file provided")
    display_name = _safe_sample_filename(file.filename)
    ext = Path(display_name).suffix.lower()
    if ext not in (".mgf", ".mzml", ".mzxml", ".xml"):
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported blank format: {ext}. Allowed: .mgf, .mzML, .mzXML",
        )
    safe_name = f"{uuid.uuid4().hex[:8]}_{display_name}"
    dest = BLANK_DIR / safe_name
    _require_disk_budget(
        BLANK_DIR,
        MAX_UPLOAD_SIZE_BYTES,
        f"blank upload {display_name}",
    )
    total_bytes = 0
    CHUNK = 4 * 1024 * 1024
    with open(dest, "wb") as f:
        while True:
            chunk = await file.read(CHUNK)
            if not chunk:
                break
            total_bytes += len(chunk)
            if total_bytes > MAX_UPLOAD_SIZE_BYTES:
                f.close()
                dest.unlink(missing_ok=True)
                raise HTTPException(
                    status_code=413,
                    detail=f"Blank file exceeds {MAX_UPLOAD_SIZE_MB} MB limit",
                )
            f.write(chunk)
    if total_bytes == 0:
        dest.unlink(missing_ok=True)
        raise HTTPException(status_code=400, detail="Empty blank file")

    # 立即验证可加载 + 报峰数 (前端可展示)
    try:
        record = await asyncio.to_thread(
            _register_managed_blank_file,
            dest,
            display_name,
        )
    except Exception as e:
        dest.unlink(missing_ok=True)
        logger.warning("空白文件 %s 加载验证失败并删除: %s", dest.name, e)
        raise HTTPException(status_code=400, detail=f"Blank file could not be parsed safely: {e}") from e

    logger.info("空白上传: %s (%.2f MB, %s 峰) → %s",
                display_name, total_bytes / (1024 * 1024), record["n_peaks"], dest.name)
    return record


@app.get("/api/blanks")
async def list_uploaded_blanks():
    """List managed blank uploads so the page can edit them after refresh."""
    blanks = await asyncio.to_thread(_list_managed_blank_uploads)
    return {"blanks": blanks, "count": len(blanks)}


@app.post("/api/blanks/{blank_token}/role")
async def update_uploaded_blank_role(
    blank_token: str,
    update: BlankRoleUpdate,
):
    """Persist a blank's analytical role without changing active task snapshots."""
    try:
        blank_path = _resolve_blank_file_token(blank_token)
    except ValueError as exc:
        status_code = 404 if "no longer exists" in str(exc) else 400
        raise HTTPException(status_code=status_code, detail=str(exc)) from exc

    token = blank_path.name
    with _blank_uploads_file_lock:
        records = _read_blank_upload_index_unlocked()
        record = dict(
            records.get(token) or _recover_legacy_blank_metadata(blank_path)
        )
        record["blank_path"] = token
        record["blank_role"] = _normalize_blank_role(update.blank_role)
        records[token] = record
        _write_blank_upload_index_unlocked(records)

    item = dict(record)
    item["active_task_ids"] = _active_blank_task_references().get(token, [])
    item["in_use"] = bool(item["active_task_ids"])
    return item


@app.post("/api/blanks/batch/delete")
async def batch_delete_uploaded_blanks(
    request: BlankBatchDeleteRequest,
):
    """Delete many managed blanks while preserving active-task inputs."""
    active_references = _active_blank_task_references()
    deleted: list[dict[str, Any]] = []
    skipped_in_use: list[dict[str, Any]] = []
    failed: list[dict[str, Any]] = []
    seen_tokens: set[str] = set()

    for supplied_token in request.blank_tokens:
        try:
            blank_path = _resolve_blank_file_token(supplied_token)
        except ValueError as exc:
            token = Path(str(supplied_token or "")).name
            failed.append({
                "blank_path": token or str(supplied_token),
                "blank_filename": _blank_original_filename(token),
                "error": str(exc),
            })
            continue

        token = blank_path.name
        if token in seen_tokens:
            continue
        seen_tokens.add(token)
        task_ids = active_references.get(token, [])
        if task_ids:
            skipped_in_use.append({
                "blank_path": token,
                "blank_filename": _blank_original_filename(token),
                "active_task_ids": task_ids,
            })
            continue

        try:
            blank_path.unlink()
        except FileNotFoundError:
            failed.append({
                "blank_path": token,
                "blank_filename": _blank_original_filename(token),
                "error": "uploaded blank file no longer exists",
            })
            continue
        except OSError as exc:
            failed.append({
                "blank_path": token,
                "blank_filename": _blank_original_filename(token),
                "error": f"Could not delete blank file: {exc}",
            })
            continue
        deleted.append({
            "blank_path": token,
            "blank_filename": _blank_original_filename(token),
        })

    metadata_warning = ""
    if deleted:
        deleted_tokens = {item["blank_path"] for item in deleted}
        try:
            with _blank_uploads_file_lock:
                records = _read_blank_upload_index_unlocked()
                changed = False
                for token in deleted_tokens:
                    if records.pop(token, None) is not None:
                        changed = True
                if changed:
                    _write_blank_upload_index_unlocked(records)
        except OSError as exc:
            metadata_warning = (
                "Blank files were deleted, but the metadata index could not "
                f"be updated: {exc}"
            )
            logger.warning(metadata_warning)

    logger.info(
        "Batch blank deletion: requested=%d deleted=%d in_use=%d failed=%d",
        len(request.blank_tokens),
        len(deleted),
        len(skipped_in_use),
        len(failed),
    )
    return {
        "requested_count": len(request.blank_tokens),
        "deleted": deleted,
        "deleted_count": len(deleted),
        "skipped_in_use": skipped_in_use,
        "skipped_in_use_count": len(skipped_in_use),
        "failed": failed,
        "failed_count": len(failed),
        "metadata_warning": metadata_warning,
    }


@app.delete("/api/blanks/{blank_token}")
async def delete_uploaded_blank(blank_token: str):
    """Delete one managed blank upload when no active task still needs it."""
    try:
        blank_path = _resolve_blank_file_token(blank_token)
    except ValueError as exc:
        status_code = 404 if "no longer exists" in str(exc) else 400
        raise HTTPException(status_code=status_code, detail=str(exc)) from exc

    token = blank_path.name
    active_references = _active_blank_task_references().get(token, [])

    if active_references:
        raise HTTPException(
            status_code=409,
            detail=(
                "Blank file is in use by queued/running task(s): "
                + ", ".join(active_references[:5])
            ),
        )

    try:
        blank_path.unlink()
    except FileNotFoundError as exc:
        raise HTTPException(
            status_code=404, detail="uploaded blank file no longer exists"
        ) from exc
    except OSError as exc:
        raise HTTPException(
            status_code=500, detail=f"Could not delete blank file: {exc}"
        ) from exc
    try:
        _remove_blank_upload_metadata(token)
    except OSError as exc:
        logger.warning(
            "Deleted blank %s but could not update its metadata index: %s",
            token,
            exc,
        )
    logger.info("Deleted managed blank upload %s", token)
    return {"deleted": token}


# ── Multi-sample untargeted workflow ────────────────────────────────

@app.get("/api/workflows/openms/status")
async def openms_workflow_status():
    """Report the exact OpenMS backend available to this V16 runtime."""
    from utils.openms_backend import probe_openms

    return await asyncio.to_thread(probe_openms)


def _v16_exact_task_polarity(task: Mapping[str, Any]) -> str:
    """Read polarity only from parsed task metadata, never from a filename."""

    polarity = str(task.get("ion_mode") or "").strip().upper()
    aliases = {
        "POSITIVE": "POS",
        "+": "POS",
        "NEGATIVE": "NEG",
        "-": "NEG",
    }
    polarity = aliases.get(polarity, polarity)
    if polarity not in {"POS", "NEG"}:
        raise ValueError("Task ion_mode must be exactly POS or NEG")
    counts = task.get("ion_mode_counts") or {}
    if isinstance(counts, Mapping):
        observed = set()
        for key, value in counts.items():
            try:
                populated = int(value or 0) > 0
            except (TypeError, ValueError):
                populated = False
            if not populated:
                continue
            normalised = aliases.get(str(key).strip().upper(), str(key).strip().upper())
            if normalised in {"POS", "NEG"}:
                observed.add(normalised)
        if observed and observed != {polarity}:
            raise ValueError("Task contains mixed or conflicting polarity records")
    return polarity


def _v16_task_binding_projection(
    task_id: str, task: Mapping[str, Any]
) -> dict[str, Any]:
    workflow = task.get("untargeted_workflow") or {}
    if not isinstance(workflow, Mapping):
        workflow = {}
    return {
        "task_id": task_id,
        "status": str(task.get("status") or ""),
        "polarity": _v16_exact_task_polarity(task),
        "result_run_id": str(task.get("result_run_id") or ""),
        "feature_xml": str(
            workflow.get("feature_xml") or workflow.get("output_features") or ""
        ),
    }


def _v16_feature_id_by_openms_unique_id(
    task: Mapping[str, Any],
) -> tuple[dict[str, str], str | None, str | None]:
    workflow = task.get("untargeted_workflow") or {}
    if not isinstance(workflow, Mapping):
        return {}, None, None
    manifest_text = str(workflow.get("manifest") or "").strip()
    if not manifest_text:
        return {}, None, None
    manifest_path = Path(manifest_text).expanduser().resolve()
    if not manifest_path.is_file():
        return {}, str(manifest_path), None
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    mapping = {}
    for feature in payload.get("features") or []:
        if not isinstance(feature, Mapping):
            continue
        unique_id = str(feature.get("openms_unique_id") or "").strip()
        feature_id = str(feature.get("feature_id") or "").strip()
        if unique_id and feature_id:
            mapping[unique_id] = feature_id
    return mapping, str(manifest_path), _file_content_sha256(manifest_path)


def _v16_aligned_rt_map(
    consensus_features: list[Mapping[str, Any]],
    samples: list[Mapping[str, Any]],
    feature_id_maps: list[Mapping[str, str]],
) -> dict[str, dict[str, float]]:
    output = {str(sample["task_id"]): {} for sample in samples}
    for consensus in consensus_features:
        for member in consensus.get("members") or []:
            if not isinstance(member, Mapping):
                continue
            try:
                sample_index = int(member.get("sample_index"))
                rt_minutes = float(member.get("rt_seconds")) / 60.0
            except (TypeError, ValueError, IndexError):
                continue
            if not math.isfinite(rt_minutes) or not 0 <= sample_index < len(samples):
                continue
            task_id = str(samples[sample_index]["task_id"])
            unique_id = str(member.get("feature_unique_id") or "").strip()
            if not unique_id:
                continue
            output[task_id][unique_id] = round(rt_minutes, 10)
            feature_id = str(feature_id_maps[sample_index].get(unique_id) or "")
            if feature_id:
                output[task_id][feature_id] = round(rt_minutes, 10)
    return output


@app.post("/api/workflows/rt-align")
async def align_completed_feature_maps(request: RTAlignmentRequest):
    """Align completed OpenMS feature maps and group cross-sample features."""
    from utils.openms_backend import align_feature_maps

    feature_paths: list[Path] = []
    samples: list[dict[str, Any]] = []
    feature_id_maps: list[dict[str, str]] = []
    with _tasks_lock:
        selected = [_tasks.get(task_id) for task_id in request.task_ids]
    for task_id, task in zip(request.task_ids, selected):
        if not task:
            raise HTTPException(
                status_code=404, detail=f"Task not found: {task_id}"
            )
        if task.get("status") != "completed":
            raise HTTPException(
                status_code=409,
                detail=f"Task {task_id} is not completed",
            )
        workflow = task.get("untargeted_workflow") or {}
        feature_path_text = (
            workflow.get("feature_xml")
            or workflow.get("output_features")
        )
        if not feature_path_text:
            raise HTTPException(
                status_code=409,
                detail=(
                    f"Task {task_id} has no OpenMS featureXML; rerun it with "
                    "MS1 feature detection enabled."
                ),
            )
        feature_path = Path(str(feature_path_text)).expanduser().resolve()
        if not feature_path.is_file():
            raise HTTPException(
                status_code=409,
                detail=f"Task {task_id} featureXML is missing from runtime storage",
            )
        feature_paths.append(feature_path)
        feature_id_map, feature_manifest, feature_manifest_sha256 = (
            _v16_feature_id_by_openms_unique_id(task)
        )
        feature_id_maps.append(feature_id_map)
        result_path = TASK_DIR / f"{task_id}_results.json"
        samples.append(
            {
                "task_id": task_id,
                "filename": str(task.get("filename") or ""),
                "feature_xml": str(feature_path),
                "feature_xml_sha256": _file_content_sha256(feature_path),
                "feature_manifest": feature_manifest,
                "feature_manifest_sha256": feature_manifest_sha256,
                "source_results_file_sha256": (
                    _file_content_sha256(result_path)
                    if result_path.is_file()
                    else None
                ),
                "source_task_metadata_sha256": _v16_semantic_sha256(
                    _v16_task_binding_projection(task_id, task)
                ),
            }
        )

    alignment_id = f"RTA-{uuid.uuid4().hex[:12]}"
    output_dir = TASK_DIR / "rt_alignments" / alignment_id
    try:
        payload = await asyncio.to_thread(
            align_feature_maps,
            feature_paths,
            output_dir,
            rt_tolerance_seconds=request.rt_tolerance_seconds,
            mz_tolerance=request.mz_tolerance,
            mz_unit=request.mz_unit,
        )
    except Exception as exc:
        logger.exception("RT alignment %s failed", alignment_id)
        raise HTTPException(
            status_code=422,
            detail=f"OpenMS RT alignment failed: {exc}",
        ) from exc
    consensus_features = list(payload.get("consensus_features") or [])
    consensus_feature_count = len(consensus_features)
    created_at = _now_iso()
    manifest_path = Path(str(payload.get("manifest") or "")).resolve()
    manifest_payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_payload.update(
        {
            "schema_version": "mssystem.v16.rt-alignment-binding/1.0.0",
            "alignment_id": alignment_id,
            "task_ids": list(request.task_ids),
            "samples": samples,
            "source_feature_xml_sha256": {
                str(sample["task_id"]): str(sample["feature_xml_sha256"])
                for sample in samples
            },
            "source_results_file_sha256": {
                str(sample["task_id"]): sample["source_results_file_sha256"]
                for sample in samples
            },
            "source_task_metadata_sha256": {
                str(sample["task_id"]): str(sample["source_task_metadata_sha256"])
                for sample in samples
            },
            "aligned_rt_minutes": _v16_aligned_rt_map(
                consensus_features, samples, feature_id_maps
            ),
            "created_at": created_at,
            "candidate_identity_used_for_alignment": False,
            "feature_group_mz_used_as_cross_polarity_identity": False,
        }
    )
    manifest_payload["semantic_sha256"] = _v16_rt_alignment_semantic_sha256(
        manifest_payload
    )
    _write_runtime_json_atomic(manifest_path, manifest_payload)
    payload.pop("consensus_features", None)
    payload["consensus_feature_table"] = {
        "row_count": consensus_feature_count,
        "manifest": payload.get("manifest"),
        "consensus_xml": payload.get("output_consensus"),
        "inline_rows_omitted": True,
    }
    payload.update(
        {
            "alignment_id": alignment_id,
            "samples": samples,
            "created_at": created_at,
            "task_ids": list(request.task_ids),
            "source_feature_xml_sha256": manifest_payload[
                "source_feature_xml_sha256"
            ],
            "source_task_metadata_sha256": manifest_payload[
                "source_task_metadata_sha256"
            ],
            "aligned_rt_mapping_counts": {
                task_id: len(values)
                for task_id, values in manifest_payload[
                    "aligned_rt_minutes"
                ].items()
            },
            "semantic_sha256": manifest_payload["semantic_sha256"],
        }
    )
    return payload


_V16_DUAL_PAIR_ID_PATTERN = re.compile(r"^DPS-[0-9a-f]{16}$")
_V16_DUAL_CANDIDATE_FIELDS = (
    "identity_key",
    "connectivity_identity",
    "candidate_connectivity",
    "inchikey",
    "formula",
    "molecular_formula",
    "compound_name",
    "name",
    "dataset_source",
    "candidate_source",
    "source",
    "lib_source_accession",
    "source_accession",
    "lib_spectrum_id",
    "potential_contaminant",
    "compound_role",
    "candidate_role",
    "role",
)
_V16_DUAL_RESULT_FIELDS = (
    "ion_mode",
    "polarity",
    "spectrum_id",
    "lcms_feature_id",
    "precursor_feature_id",
    "lcms_feature_rt_minutes",
    "retention_time",
    "retention_time_unit",
    "lcms_feature_mz",
    "precursor_mz",
    "ion_identity_status",
    "ion_identity_network_id",
    "ion_identity",
    "ion_identity_neutral_mass",
    "lcms_adduct_group_id",
    "lcms_adduct",
    "lcms_neutral_mass",
)


def _v16_dual_sidecar_path(pair_id: str) -> Path:
    value = str(pair_id or "").strip()
    if not _V16_DUAL_PAIR_ID_PATTERN.fullmatch(value):
        raise ValueError("Dual-polarity pair ID is invalid")
    directory = (TASK_DIR / "dual_polarity_shadow").resolve()
    path = (directory / f"{value}.json").resolve()
    if path.parent != directory:
        raise ValueError("Dual-polarity sidecar path escaped runtime storage")
    return path


def _v16_stream_flagship_metadata(path: Path) -> dict[str, Any]:
    """Read only the small V15 metadata object from a potentially huge file."""

    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        return dict(payload.get("metadata") or {})
    with path.open("rb") as handle:
        return dict(next(ijson.items(handle, "metadata"), {}) or {})


def _v16_current_flagship_metadata(
    path: Path, *, task_id: str
) -> dict[str, Any]:
    metadata = _v16_stream_flagship_metadata(path)
    policy = _load_v15_flagship_policy()
    if not (
        str(metadata.get("task_id") or "") == task_id
        and metadata.get("algorithm_id") == policy.get("algorithm_id")
        and metadata.get("policy_sha256") == _flagship_policy_sha256(policy)
        and metadata.get("identity_resolution_status") == "available"
    ):
        raise RuntimeError(
            f"Task {task_id} current V15 flagship artifact is absent or stale"
        )
    return metadata


def _v16_compact_candidate(candidate: Any) -> dict[str, Any] | None:
    if not isinstance(candidate, Mapping):
        return None
    nested = candidate.get("candidate")
    sources = [nested, candidate] if isinstance(nested, Mapping) else [candidate]
    compact: dict[str, Any] = {}
    for source in sources:
        for field in _V16_DUAL_CANDIDATE_FIELDS:
            value = source.get(field)
            if value is not None and value != "":
                compact[field] = value
    return compact or None


def _v16_candidate_independent_result_projection(
    row: Mapping[str, Any]
) -> dict[str, Any]:
    """Retain only fields allowed to construct or audit the frozen graph."""

    compact = {
        field: row.get(field)
        for field in _V16_DUAL_RESULT_FIELDS
        if row.get(field) is not None
    }
    for field in ("ion_identity_evidence", "lcms_adduct_group_evidence"):
        evidence = row.get(field)
        if not isinstance(evidence, Mapping):
            continue
        compact[field] = {
            key: evidence.get(key)
            for key in (
                "status",
                "network_id",
                "neutral_mass",
                "ion_identity",
                "runtime_influence",
                "identity_conflict",
                "neutral_mass_consistent",
                "member_feature_ids",
                "member_feature_count",
            )
            if evidence.get(key) is not None
        }
    return compact


def _v16_stream_flagship_candidate_projection(
    path: Path, *, metadata: Mapping[str, Any]
) -> dict[str, Any]:
    """Stream a minimal member-index/candidate view for shadow annotation."""

    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        entities_source = payload.get("entities") or []
    else:
        handle = path.open("rb")
        entities_source = ijson.items(handle, "entities.item")

    entities = []
    try:
        for entity in entities_source:
            if not isinstance(entity, Mapping):
                continue
            member_indices = sorted(
                {
                    int(value)
                    for value in entity.get("member_result_indices") or []
                    if str(value).strip()
                }
            )
            if not member_indices:
                continue
            winner = _v16_compact_candidate(entity.get("winner"))
            leading = _v16_compact_candidate(
                entity.get("leading_candidate_when_abstained")
            )
            candidates = [
                compact
                for candidate in entity.get("candidates") or []
                if (compact := _v16_compact_candidate(candidate)) is not None
            ]
            recall = entity.get("candidate_recall_support") or {}
            expanded = [
                compact
                for candidate in (
                    recall.get("expanded_candidates") or []
                    if isinstance(recall, Mapping)
                    else []
                )
                if (compact := _v16_compact_candidate(candidate)) is not None
            ]
            entities.append(
                {
                    "entity_id": str(entity.get("entity_id") or ""),
                    "member_result_indices": member_indices,
                    "winner": winner,
                    "leading_candidate_when_abstained": leading,
                    "candidates": candidates,
                    "candidate_recall_support": {
                        "expanded_candidates": expanded
                    },
                }
            )
    finally:
        if "handle" in locals():
            handle.close()
    return {
        "metadata": {
            "task_id": str(metadata.get("task_id") or ""),
            "algorithm_id": str(metadata.get("algorithm_id") or ""),
            "policy_sha256": str(metadata.get("policy_sha256") or ""),
            "entity_semantic_sha256": str(
                metadata.get("entity_semantic_sha256") or ""
            ),
            "identity_resolution_status": str(
                metadata.get("identity_resolution_status") or ""
            ),
            "projection_role": "v16_dual_polarity_candidate_support_only",
        },
        "entities": entities,
    }


def _v16_dual_task_artifacts(
    task_id: str,
    task: Mapping[str, Any],
    *,
    expected_polarity: str,
) -> dict[str, Any]:
    if str(task.get("status") or "") != "completed":
        raise ValueError(f"Task {task_id} is not completed")
    actual_polarity = _v16_exact_task_polarity(task)
    if actual_polarity != expected_polarity:
        raise ValueError(
            f"Task {task_id} ion_mode is {actual_polarity}, not {expected_polarity}"
        )
    result_path = TASK_DIR / f"{task_id}_results.json"
    if not result_path.is_file():
        raise FileNotFoundError(
            f"Task {task_id} saved scan-level results are unavailable"
        )
    flagship_path = _flagship_artifact_path(TASK_DIR, task_id)
    if not flagship_path.is_file():
        raise RuntimeError(
            f"Task {task_id} current V15 flagship artifact is absent or stale"
        )
    flagship_metadata = _v16_current_flagship_metadata(
        flagship_path, task_id=task_id
    )
    workflow = task.get("untargeted_workflow") or {}
    if not isinstance(workflow, Mapping):
        workflow = {}
    feature_text = str(
        workflow.get("feature_xml") or workflow.get("output_features") or ""
    ).strip()
    if not feature_text:
        raise RuntimeError(f"Task {task_id} has no OpenMS featureXML binding")
    feature_path = Path(feature_text).expanduser().resolve()
    if not feature_path.is_file():
        raise RuntimeError(f"Task {task_id} featureXML is unavailable")
    feature_manifest_text = str(workflow.get("manifest") or "").strip()
    feature_manifest_path = (
        Path(feature_manifest_text).expanduser().resolve()
        if feature_manifest_text
        else None
    )
    return {
        "task_id": task_id,
        "polarity": actual_polarity,
        "result_path": result_path,
        "flagship_path": flagship_path,
        "flagship_metadata": flagship_metadata,
        "feature_xml_path": feature_path,
        "feature_manifest_path": feature_manifest_path,
    }


def _v16_dual_file_bindings(
    artifacts: Mapping[str, Mapping[str, Any]],
    flagship_payloads: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    bindings = {}
    for task_id in sorted(artifacts):
        item = artifacts[task_id]
        flagship = flagship_payloads[task_id]
        metadata = flagship.get("metadata") or {}
        bindings[task_id] = {
            "results_path": str(Path(item["result_path"]).resolve()),
            "results_file_sha256": _file_content_sha256(
                Path(item["result_path"])
            ),
            "v15_flagship_path": str(Path(item["flagship_path"]).resolve()),
            "v15_flagship_file_sha256": _file_content_sha256(
                Path(item["flagship_path"])
            ),
            "v15_flagship_entity_semantic_sha256": str(
                metadata.get("entity_semantic_sha256") or ""
            ),
            "v15_flagship_algorithm_id": str(
                metadata.get("algorithm_id") or ""
            ),
            "feature_xml_path": str(Path(item["feature_xml_path"]).resolve()),
            "feature_xml_sha256": _file_content_sha256(
                Path(item["feature_xml_path"])
            ),
            "feature_manifest_path": (
                str(Path(item["feature_manifest_path"]).resolve())
                if item.get("feature_manifest_path") is not None
                else None
            ),
            "feature_manifest_sha256": (
                _file_content_sha256(Path(item["feature_manifest_path"]))
                if item.get("feature_manifest_path") is not None
                and Path(item["feature_manifest_path"]).is_file()
                else None
            ),
        }
    return bindings


def _v16_openms_group_validation(
    artifact: Mapping[str, Any],
    *,
    task_polarity: str,
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate one feature manifest without changing saved result rows."""

    manifest_path = artifact.get("feature_manifest_path")
    payload: Mapping[str, Any] = {}
    manifest_sha256 = None
    if manifest_path is not None and Path(manifest_path).is_file():
        path = Path(manifest_path)
        manifest_sha256 = _file_content_sha256(path)
        loaded = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(loaded, Mapping):
            payload = loaded
    requirements = (
        (policy.get("candidate_independent_anchor") or {}).get(
            "openms_group_required"
        )
        or {}
    )
    return _v16_validate_openms_adduct_groups(
        payload,
        task_polarity=task_polarity,
        manifest_sha256=manifest_sha256,
        feature_xml_sha256=_file_content_sha256(
            Path(artifact["feature_xml_path"])
        ),
        feature_xml_path=artifact["feature_xml_path"],
        requirements=requirements,
    )


def _v16_bind_openms_group_evidence(
    rows: list[dict[str, Any]], validation: Mapping[str, Any]
) -> dict[str, Any]:
    """Bind validated feature evidence to an in-memory shadow projection."""

    evidence_by_feature = validation.get("evidence_by_feature_id") or {}
    bound_rows = 0
    binding_mismatches = 0
    for row in rows:
        feature_id = str(
            row.get("lcms_feature_id")
            or row.get("precursor_feature_id")
            or ""
        ).strip()
        evidence = evidence_by_feature.get(feature_id)
        if not isinstance(evidence, Mapping):
            continue
        row_group = str(row.get("lcms_adduct_group_id") or "").strip()
        row_adduct = str(row.get("lcms_adduct") or "").strip()
        if (
            row_group != str(evidence.get("network_id") or "").strip()
            or row_adduct != str(evidence.get("ion_identity") or "").strip()
        ):
            binding_mismatches += 1
            continue
        row["lcms_adduct_group_evidence"] = {
            key: deepcopy(value)
            for key, value in evidence.items()
            if key
            in {
                "status",
                "network_id",
                "neutral_mass",
                "ion_identity",
                "runtime_influence",
                "production_runtime_influence",
                "identity_conflict",
                "neutral_mass_consistent",
                "member_feature_ids",
                "member_feature_count",
                "manifest_sha256",
                "feature_xml_sha256",
                "candidate_identity_used",
            }
        }
        bound_rows += 1
    audit = _v16_compact_openms_validation_audit(validation)
    audit.update(
        {
            "projected_result_row_count": len(rows),
            "bound_result_row_count": bound_rows,
            "result_binding_mismatch_count": binding_mismatches,
            "saved_result_rows_modified": False,
        }
    )
    return audit


def _v16_prepare_rt_alignment(
    alignment_id: str | None,
    *,
    task_ids: set[str],
    tasks: Mapping[str, Mapping[str, Any]],
    result_semantic_sha256: Mapping[str, str],
    file_bindings: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    if alignment_id is None:
        return None, None
    manifest_path = (
        TASK_DIR / "rt_alignments" / alignment_id / "rt-alignment.json"
    ).resolve()
    alignment_root = (TASK_DIR / "rt_alignments").resolve()
    if alignment_root not in manifest_path.parents or not manifest_path.is_file():
        raise RuntimeError("Requested RT alignment artifact is unavailable")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "mssystem.v16.rt-alignment-binding/1.0.0":
        raise RuntimeError("RT alignment lacks the V16 source-hash binding")
    if set(str(value) for value in manifest.get("task_ids") or []) != task_ids:
        raise RuntimeError("RT alignment is bound to different task IDs")
    declared_semantic = str(manifest.get("semantic_sha256") or "")
    observed_semantic = _v16_rt_alignment_semantic_sha256(manifest)
    if not declared_semantic or declared_semantic != observed_semantic:
        raise RuntimeError("RT alignment semantic hash does not verify")

    expected_feature_hashes = {
        task_id: str(file_bindings[task_id]["feature_xml_sha256"])
        for task_id in sorted(task_ids)
    }
    if dict(manifest.get("source_feature_xml_sha256") or {}) != (
        expected_feature_hashes
    ):
        raise RuntimeError("RT alignment featureXML hashes are stale")
    expected_result_file_hashes = {
        task_id: str(file_bindings[task_id]["results_file_sha256"])
        for task_id in sorted(task_ids)
    }
    if dict(manifest.get("source_results_file_sha256") or {}) != (
        expected_result_file_hashes
    ):
        raise RuntimeError("RT alignment result-file hashes are stale")
    expected_task_hashes = {
        task_id: _v16_semantic_sha256(
            _v16_task_binding_projection(task_id, tasks[task_id])
        )
        for task_id in sorted(task_ids)
    }
    if dict(manifest.get("source_task_metadata_sha256") or {}) != (
        expected_task_hashes
    ):
        raise RuntimeError("RT alignment task-metadata hashes are stale")

    normalised = deepcopy(manifest)
    normalised.pop("semantic_sha256", None)
    normalised.pop("artifact_sha256", None)
    normalised["source_results_sha256"] = dict(result_semantic_sha256)
    normalised["source_task_metadata_sha256"] = expected_task_hashes
    return normalised, {
        "alignment_id": alignment_id,
        "path": str(manifest_path),
        "file_sha256": _file_content_sha256(manifest_path),
        "semantic_sha256": declared_semantic,
    }


def _v16_task_count_snapshot(
    tasks: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    fields = (
        "unique_compounds",
        "main_result_count",
        "result_record_count",
        "isomer_group_result_count",
        "unresolved_precursor_event_count",
    )
    return {
        task_id: {field: task.get(field) for field in fields}
        for task_id, task in sorted(tasks.items())
    }


@app.post("/api/workflows/dual-polarity-shadow")
async def build_dual_polarity_shadow(request: DualPolarityShadowRequest):
    """Build one explicit, read-only POS/NEG evidence sidecar."""

    task_ids = {
        request.positive_task_id,
        request.negative_task_id,
    }
    with _tasks_lock:
        selected_tasks = {
            task_id: deepcopy(_tasks.get(task_id))
            for task_id in task_ids
            if _tasks.get(task_id) is not None
        }
    missing = sorted(task_ids - set(selected_tasks))
    if missing:
        raise HTTPException(
            status_code=404, detail=f"Task not found: {missing[0]}"
        )
    try:
        artifacts = {
            request.positive_task_id: _v16_dual_task_artifacts(
                request.positive_task_id,
                selected_tasks[request.positive_task_id],
                expected_polarity="POS",
            ),
            request.negative_task_id: _v16_dual_task_artifacts(
                request.negative_task_id,
                selected_tasks[request.negative_task_id],
                expected_polarity="NEG",
            ),
        }
    except FileNotFoundError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except (ValueError, RuntimeError) as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc

    try:
        results = {
            task_id: [
                _v16_candidate_independent_result_projection(row)
                for row in _iter_legacy_result_artifact(item["result_path"])
                if isinstance(row, Mapping)
            ]
            for task_id, item in artifacts.items()
        }
        flagship_payloads = {
            task_id: _v16_stream_flagship_candidate_projection(
                Path(item["flagship_path"]),
                metadata=item["flagship_metadata"],
            )
            for task_id, item in artifacts.items()
        }
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise HTTPException(
            status_code=409,
            detail=f"Saved V15 evidence artifact is unreadable: {exc}",
        ) from exc

    for task_id, flagship in flagship_payloads.items():
        if str((flagship.get("metadata") or {}).get("task_id") or "") != task_id:
            raise HTTPException(
                status_code=409,
                detail=f"V15 flagship artifact belongs to another task: {task_id}",
            )

    try:
        policy = _load_v16_dual_polarity_policy()
        openms_validations = {
            task_id: _v16_openms_group_validation(
                item,
                task_polarity=str(item["polarity"]),
                policy=policy,
            )
            for task_id, item in artifacts.items()
        }
        openms_validation_audits = {
            task_id: _v16_bind_openms_group_evidence(
                results[task_id], openms_validations[task_id]
            )
            for task_id in sorted(artifacts)
        }
        before_bindings = _v16_dual_file_bindings(artifacts, flagship_payloads)
        result_semantic_sha256 = {
            task_id: _v16_semantic_sha256(rows)
            for task_id, rows in sorted(results.items())
        }
        rt_alignment, rt_alignment_binding = _v16_prepare_rt_alignment(
            request.rt_alignment_id,
            task_ids=task_ids,
            tasks=selected_tasks,
            result_semantic_sha256=result_semantic_sha256,
            file_bindings=before_bindings,
        )
    except (OSError, ValueError, TypeError, RuntimeError, json.JSONDecodeError) as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc

    request_material = {
        "positive_task_id": request.positive_task_id,
        "negative_task_id": request.negative_task_id,
        "same_sample_confirmed": True,
        "same_lc_method_confirmed": True,
        "pair_label": request.pair_label,
        "rt_alignment": rt_alignment_binding,
        "source_files": before_bindings,
        "policy_sha256": _v16_semantic_sha256(policy),
    }
    request_sha256 = _v16_semantic_sha256(request_material)
    pair_id = f"DPS-{request_sha256[:16]}"
    sidecar_path = _v16_dual_sidecar_path(pair_id)
    if sidecar_path.is_file():
        try:
            existing = json.loads(sidecar_path.read_text(encoding="utf-8"))
            existing_hash = str(
                (existing.get("metadata") or {}).get("request_sha256") or ""
            )
            declared_sidecar_hash = str(
                (existing.get("metadata") or {}).get(
                    "sidecar_semantic_sha256"
                )
                or ""
            )
            if (
                existing_hash == request_sha256
                and declared_sidecar_hash
                == _v16_sidecar_semantic_sha256(existing)
            ):
                return {
                    "status": "available",
                    "pair_id": pair_id,
                    "reused": True,
                    "artifact": str(sidecar_path),
                    "sidecar": existing,
                }
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            pass
        raise HTTPException(
            status_code=409,
            detail="Existing dual-polarity sidecar does not match its source binding",
        )

    task_counts_before = _v16_task_count_snapshot(selected_tasks)
    try:
        sidecar = await asyncio.to_thread(
            _build_v16_dual_polarity_shadow,
            selected_tasks[request.positive_task_id],
            results[request.positive_task_id],
            selected_tasks[request.negative_task_id],
            results[request.negative_task_id],
            rt_alignment=rt_alignment,
            candidate_support_by_task=flagship_payloads,
            policy=policy,
        )
    except (TypeError, ValueError, RuntimeError) as exc:
        raise HTTPException(
            status_code=422,
            detail=f"Dual-polarity graph was not evaluable: {exc}",
        ) from exc

    after_bindings = _v16_dual_file_bindings(artifacts, flagship_payloads)
    if after_bindings != before_bindings:
        raise HTTPException(
            status_code=409,
            detail="Source V15 evidence changed during shadow graph construction",
        )
    with _tasks_lock:
        current_tasks = {
            task_id: deepcopy(_tasks.get(task_id) or {}) for task_id in task_ids
        }
    task_counts_after = _v16_task_count_snapshot(current_tasks)
    if task_counts_after != task_counts_before:
        raise HTTPException(
            status_code=409,
            detail="Task compound counts changed during shadow graph construction",
        )

    metadata = sidecar.setdefault("metadata", {})
    metadata.update(
        {
            "pair_id": pair_id,
            "pair_label": request.pair_label,
            "request_sha256": request_sha256,
            "service_generation": SERVICE_GENERATION,
            "service_version": app.version,
            "runtime_snapshot_id": _current_runtime_snapshot_id(),
            "source_file_sha256": before_bindings,
            "source_results_semantic_projection": (
                "candidate_independent_anchor_fields_only; full immutable files "
                "bound separately by results_file_sha256"
            ),
            "rt_alignment_artifact": rt_alignment_binding,
            "openms_adduct_group_validation": openms_validation_audits,
            "source_file_sha256_unchanged_after_build": True,
            "task_compound_counts_before": task_counts_before,
            "task_compound_counts_after": task_counts_after,
            "task_compound_counts_unchanged": True,
        }
    )
    invariants = sidecar.setdefault("invariants", {})
    invariants.update(
        {
            "production_candidate_scores_changed": False,
            "production_decisions_changed": False,
            "dual_polarity_links_counted_as_compounds": False,
            "candidate_identity_used_for_entity_linking": False,
            "filenames_used_for_pairing": False,
            "cross_polarity_ms2_similarity_affects_top1": False,
        }
    )
    sidecar["production_top1_changed"] = False
    sidecar["production_candidate_scores_changed"] = False
    sidecar["production_decisions_changed"] = False
    sidecar["compound_count_increment"] = 0
    sidecar["dual_polarity_links_counted_as_compounds"] = False
    sidecar["candidate_identity_used_for_entity_linking"] = False
    sidecar["filenames_used_for_pairing"] = False
    sidecar["cross_polarity_ms2_similarity_affects_top1"] = False
    metadata["sidecar_semantic_sha256"] = _v16_sidecar_semantic_sha256(
        sidecar
    )
    try:
        _write_runtime_json_atomic(sidecar_path, sidecar)
    except OSError as exc:
        raise HTTPException(
            status_code=500, detail=f"Could not write dual-polarity sidecar: {exc}"
        ) from exc
    return {
        "status": "available",
        "pair_id": pair_id,
        "reused": False,
        "artifact": str(sidecar_path),
        "sidecar": sidecar,
    }


@app.get("/api/workflows/dual-polarity-shadow/{pair_id}")
async def get_dual_polarity_shadow(pair_id: str):
    """Read one V16 dual-polarity sidecar without materializing anything."""

    try:
        sidecar_path = _v16_dual_sidecar_path(pair_id)
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    if not sidecar_path.is_file():
        raise HTTPException(status_code=404, detail="Dual-polarity sidecar not found")
    try:
        sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise HTTPException(
            status_code=500, detail=f"Dual-polarity sidecar is unreadable: {exc}"
        ) from exc
    declared = str(
        (sidecar.get("metadata") or {}).get("sidecar_semantic_sha256") or ""
    )
    if not declared or declared != _v16_sidecar_semantic_sha256(sidecar):
        raise HTTPException(
            status_code=409, detail="Dual-polarity sidecar semantic hash is invalid"
        )
    return {
        "status": "available",
        "pair_id": pair_id,
        "artifact": str(sidecar_path),
        "sidecar": sidecar,
    }


# ── Task Management ─────────────────────────────────────────────────

@app.get("/api/tasks")
async def list_tasks(
    status: Optional[str] = Query(None, description="Filter by status"),
    sort: str = Query("upload_time_desc", description="Sort order"),
    projection: Literal["full", "minimal"] = Query(
        "full", description="Response field projection"
    ),
):
    """List all analysis tasks."""
    queue_snapshot = _full_queue_snapshot()
    waiting_ids = list(queue_snapshot.get("queued") or [])
    queue_positions = {
        task_id: position
        for position, task_id in enumerate(waiting_ids, start=1)
    }
    with _tasks_lock:
        tasks_list = list(_tasks.values())

    if status:
        tasks_list = [t for t in tasks_list if t["status"] == status]

    # Sort
    if sort == "upload_time_desc":
        tasks_list.sort(key=lambda t: t.get("upload_time", ""), reverse=True)
    elif sort == "upload_time_asc":
        tasks_list.sort(key=lambda t: t.get("upload_time", ""))

    if projection == "minimal":
        return {
            "tasks": [
                {
                    "id": task["id"],
                    "filename": task["filename"],
                    "status": task["status"],
                    "ion_mode": task.get("ion_mode"),
                    "upload_time": task.get("upload_time"),
                }
                for task in tasks_list
            ],
            "count": len(tasks_list),
            "projection": "minimal",
        }

    summaries = []
    for t in tasks_list:
        # 提取匹配化合物名称列表（供前端搜索用，最多20个不重复）
        compound_names = []
        seen_names = set()
        for r in t.get("results", []):
            # Task cards are also user-visible.  Apply the same accession/name
            # cleanup as the result table so a cached Z/AKOS/numeric label
            # cannot leak back through search suggestions after an upgrade.
            _normalize_result_compound_names([r])
            m = r.get("match")
            if m and m.get("compound_name"):
                name = m["compound_name"]
                if name not in seen_names:
                    seen_names.add(name)
                    compound_names.append(name)
                    if len(compound_names) >= 20:
                        break
        summaries.append({
            "id": t["id"],
            "filename": t["filename"],
            "status": t["status"],
            "progress": t["progress"],
            "state_version": t.get("state_version", 1),
            "run_attempt": t.get("run_attempt", 0),
            "active_run_id": t.get("active_run_id"),
            "result_run_id": t.get("result_run_id"),
            "message": _task_message_for_display(t),
            "queue_position": queue_positions.get(t["id"]),
            "queue_length": len(waiting_ids),
            "queue_order": t.get("queue_order"),
            "queue_reorderable": bool(
                queue_snapshot.get("reorder_supported")
                and t.get("status") == TaskState.QUEUED.value
                and t["id"] in queue_positions
            ),
            "upload_time": t["upload_time"],
            "start_time": t.get("start_time"),
            "end_time": t.get("end_time"),
            "total_spectra": t.get("total_spectra", 0),
            "ms1_scans": t.get("ms1_scans"),
            "matched_spectra": t.get("matched_spectra", 0),
            "precursor_features": t.get("precursor_features"),
            "matched_precursor_features": t.get("matched_precursor_features"),
            "high_precursor_features": t.get("high_precursor_features"),
            "unique_compounds": t.get("unique_compounds"),
            "legacy_unique_library_candidates": t.get(
                "legacy_unique_library_candidates"
            ),
            "main_result_count": t.get("main_result_count"),
            "result_record_count": t.get("result_record_count"),
            "unresolved_precursor_event_count": t.get(
                "unresolved_precursor_event_count"
            ),
            "isomer_group_result_count": t.get("isomer_group_result_count"),
            "compound_count_semantics": t.get("compound_count_semantics"),
            "compound_confidence_summary": t.get(
                "compound_confidence_summary", {}
            ),
            "identification_funnel": t.get("identification_funnel", {}),
            "deconvolution_impact": t.get("deconvolution_impact", {}),
            "evidence_governance_summary": t.get(
                "evidence_governance_summary", {}
            ),
            "identification_contract_summary": t.get(
                "identification_contract_summary", {}
            ),
            "flagship_summary": t.get("flagship_summary", {}),
            "analysis_completion_check": t.get("analysis_completion_check", {}),
            "primary_dependency_preflight": t.get("primary_dependency_preflight", {}),
            "evirank_mv_summary": t.get("evirank_mv_summary", {}),
            "entity_bundle_summary": t.get("entity_bundle_summary", {}),
            "claim_contract_summary": t.get("claim_contract_summary", {}),
            "compound_conclusion_summary": t.get(
                "compound_conclusion_summary", {}
            ),
            "reference_standard_qc_summary": t.get(
                "reference_standard_qc_summary", {}
            ),
            "contaminant_unique_compounds": t.get("contaminant_unique_compounds"),
            "module_availability": t.get("module_availability", {}),
            "analysis_warnings": t.get("analysis_warnings", []),
            "peak_processing_summary": t.get("peak_processing_summary", {}),
            "untargeted_workflow": t.get("untargeted_workflow", {}),
            "consensus_spectrum_summary": t.get(
                "consensus_spectrum_summary", {}
            ),
            "expected_compounds_summary": t.get(
                "expected_compounds_summary", {}
            ),
            "reliability_annotation_summary": t.get(
                "reliability_annotation_summary", {}
            ),
            "ion_identity_networks": t.get("ion_identity_networks"),
            "ion_identity_networked_features": t.get("ion_identity_networked_features"),
            "ion_identity_linked_annotations": t.get("ion_identity_linked_annotations"),
            "chimeric_spectra": t.get("chimeric_spectra"),
            "deconvolved_spectra": t.get("deconvolved_spectra"),
            "match_rate": round(t.get("match_rate", 0), 1),
            "file_size_mb": t.get("file_size_mb", 0),
            "params": t.get("params", {}),
            "blank_summary": t.get("blank_summary", {}),
            "isomer_literature_summary": t.get(
                "isomer_literature_summary", {}
            ),
            "evograph_fragmentation_summary": t.get(
                "evograph_fragmentation_summary", {}
            ),
            "ion_mode": t.get("ion_mode"),
            "ion_mode_counts": t.get("ion_mode_counts", {}),
            "compound_names": compound_names,
            "completed_batch": t.get("completed_batch"),
            "ms1_feature_count": t.get("ms1_feature_count"),
        })

    return {"tasks": summaries, "count": len(summaries)}


# ── Task Search (MUST be before {task_id} to avoid route capture) ──

def _task_message_for_display(task: dict) -> str:
    """Hide obsolete natural/literature split wording from completed tasks."""
    if task.get("completed_batch"):
        return str(task.get("message") or "已完成V16批处理结果，只读接入")
    if task.get("status") != "completed":
        return str(task.get("message") or "")
    matched = int(task.get("matched_spectra") or task.get("matched_count") or 0)
    total = int(task.get("total_spectra") or 0)
    ms1_scans = int(
        task.get("ms1_scans")
        or (task.get("module_availability") or {}).get("ms1_scan_count")
        or 0
    )
    unique = task.get("unique_compounds")
    precursor_features = task.get("precursor_features")
    matched_features = task.get("matched_precursor_features")
    rate = float(task.get("match_rate") or (matched / total * 100 if total else 0.0))
    if task.get("compound_count_semantics") == (
        "sample_supported_main_results_excluding_unresolved_events"
    ):
        unique_text = f", {int(unique or 0)} main compound results"
        unresolved = int(task.get("unresolved_precursor_event_count") or 0)
        if unresolved:
            unique_text += f", {unresolved} unresolved precursor events excluded"
    else:
        unique_text = (
            f", {int(unique)} unique library candidates"
            if unique is not None
            else ""
        )
    feature_text = (
        f", {int(matched_features or 0)}/{int(precursor_features)} MS2 precursor features matched"
        if precursor_features is not None else ""
    )
    blank = task.get("blank_summary") or {}
    if not blank.get("requested"):
        blank_text = ""
    elif blank.get("filter_policy") == "soft_evidence":
        blank_text = (
            f", blank-flagged {int(blank.get('flagged_spectra') or 0)}"
            f" (clean-view {int(blank.get('clean_view_excluded_spectra') or 0)})"
        )
    else:
        blank_text = (
            f", blank-filtered {int(blank.get('excluded_spectra') or 0)}"
        )
    ion_networks = int(task.get("ion_identity_networks") or 0)
    ion_text = f", {ion_networks} ion-identity networks" if ion_networks else ""
    chimeric = int(task.get("chimeric_spectra") or 0)
    deconvolved = int(task.get("deconvolved_spectra") or 0)
    deconvolution_text = (
        f", chimeric deconvolution {deconvolved}/{chimeric}" if chimeric else ""
    )
    expected = task.get("expected_compounds_summary") or {}
    expected_status = str(expected.get("status") or "")
    if expected_status == "ok":
        expected_text = (
            f", Expected-MS1 {int(expected.get('detected_feature_count') or 0)} "
            "candidate features (separate layer)"
        )
    elif expected.get("requested") and expected_status:
        expected_text = f", Expected Compounds {expected_status}"
    else:
        expected_text = ""
    availability = task.get("module_availability") or {}
    warnings = task.get("analysis_warnings") or availability.get("warnings") or []
    warning_codes = {
        str(item.get("code") or "")
        for item in warnings if isinstance(item, dict)
    }
    if "continuous_ms1_missing" in warning_codes:
        module_text = ", requested MS1-dependent modules unavailable (0 survey MS1)"
    elif warnings:
        module_text = ", module applicability warning"
    else:
        module_text = ""
    scan_text = (
        f"{ms1_scans} MS1 + {matched}/{total} MS2"
        if ms1_scans > 0 else f"{matched}/{total} spectra"
    )
    return (
        f"Completed: {scan_text} matched{feature_text}"
        f"{unique_text}{blank_text}{deconvolution_text}{ion_text}"
        f"{expected_text}{module_text} ({rate:.1f}%)"
    )


def _apply_final_compound_count_summary(
    task: dict[str, Any], summary: Mapping[str, Any] | None
) -> bool:
    """Make the task-level compound count match the final V15 result table.

    ``unique_compounds`` historically counted distinct scan-level library
    leads.  It included unresolved cross-formula events and therefore did not
    represent the number shown by the final compound table.  Preserve that
    legacy diagnostic explicitly, but make the public task count follow
    ``main_result_count`` from the scope-aware production artifact.
    """

    if not isinstance(summary, Mapping):
        return False
    main = summary.get("main_result_count")
    records = summary.get("unique_report_rows")
    if main is None or records is None:
        return False
    before = {
        key: task.get(key)
        for key in (
            "legacy_unique_library_candidates",
            "unique_compounds",
            "main_result_count",
            "result_record_count",
            "unresolved_precursor_event_count",
            "isomer_group_result_count",
            "compound_count_semantics",
        )
    }
    if task.get("legacy_unique_library_candidates") is None:
        legacy = task.get("unique_compounds")
        if legacy is not None:
            task["legacy_unique_library_candidates"] = int(legacy)
    task.update(
        {
            "unique_compounds": int(main),
            "main_result_count": int(main),
            "result_record_count": int(records),
            "unresolved_precursor_event_count": int(
                summary.get("unresolved_precursor_event_count") or 0
            ),
            "isomer_group_result_count": int(
                summary.get("isomer_group_detected_count") or 0
            ),
            "compound_count_semantics": (
                "sample_supported_main_results_excluding_unresolved_events"
            ),
        }
    )
    after = {key: task.get(key) for key in before}
    return before != after


def _refresh_reference_standard_qc_projection(
    task: dict, prepared_results: list[dict]
) -> None:
    """Rebuild stale post-hoc QC against the current reportable-scope semantics."""
    params = task.get("params") or {}
    set_id = str(params.get("reference_standard_set_id") or "").strip()
    task_id = str(task.get("id") or "").strip()
    if not set_id or not task_id:
        return
    artifact = _reference_standard_qc_artifact_path(TASK_DIR, task_id)
    existing: dict[str, Any] = {}
    if artifact.is_file():
        try:
            existing = _read_reference_standard_qc_artifact(artifact)
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            existing = {}
    if (
        int(existing.get("schema_version") or 0)
        == REFERENCE_STANDARD_QC_SCHEMA_VERSION
        and existing.get("presentation_calibration_version")
        == PRESENTATION_CALIBRATION_VERSION
        and str((existing.get("reference_standard_set") or {}).get("id") or "")
        == set_id
    ):
        return
    polarity = str(task.get("ion_mode") or params.get("ion_mode") or "").upper()
    if polarity not in {"POS", "NEG"}:
        modes = {
            str(row.get("ion_mode") or row.get("polarity") or "").upper()
            for row in prepared_results
            if str(row.get("ion_mode") or row.get("polarity") or "").upper()
            in {"POS", "NEG"}
        }
        if len(modes) != 1:
            return
        polarity = next(iter(modes))
    tolerance = float(
        existing.get("mass_tolerance_ppm_for_observable_event_audit") or 10.0
    )
    payload = _build_reference_standard_qc(
        task_id,
        prepared_results,
        set_id=set_id,
        polarity=polarity,
        filename=str(task.get("filename") or ""),
        mass_tolerance_ppm=tolerance,
        selection_source=str(
            params.get("reference_standard_set_selection_source")
            or "saved_explicit_selection_reprojection"
        ),
    )
    payload["presentation_calibration_version"] = (
        PRESENTATION_CALIBRATION_VERSION
    )
    task["reference_standard_qc_summary"] = (
        _write_reference_standard_qc_artifact(artifact, payload)
    )


def _refresh_blank_pairing_semantics(task: dict) -> None:
    """Upgrade legacy blank summaries without pretending polarity is pairing."""
    summary = task.get("blank_summary")
    if not isinstance(summary, dict) or not summary.get("requested"):
        return
    summary = dict(summary)
    roles = dict(summary.get("blank_roles") or {})
    unspecified = int(roles.get("unspecified") or 0)
    assigned_count = int(summary.get("assigned_blank_replicate_count") or 0)
    params = task.get("params") if isinstance(task.get("params"), dict) else {}
    method_confirmed = bool(
        params.get("blank_method_comparability_confirmed", False)
    )
    exclusion_confirmed = bool(
        params.get("blank_background_exclusion_confirmed", False)
    )
    eligible_role_names = {
        "background_blank", "solvent_blank", "method_blank",
        "extraction_blank",
    }
    nonzero_roles = {
        str(role) for role, count in roles.items() if int(count or 0) > 0
    }
    roles_eligible = bool(nonzero_roles) and nonzero_roles <= eligible_role_names
    exclusion_eligible = bool(
        assigned_count
        and roles_eligible
        and method_confirmed
        and exclusion_confirmed
    )
    policy = str(
        summary.get("filter_policy")
        or params.get("blank_filter_policy")
        or DEFAULT_BLANK_FILTER_POLICY
    )
    if unspecified:
        summary.update({
            "blank_pairing_status": "polarity_only_unverified_method_pairing",
            "unverified_assigned_blank_count": unspecified,
            "blank_pairing_warning": (
                "同极性空白包含未定角色文件，批次与色谱方法配对未经核验；"
                "这些空白只产生背景警示，不据此剔除样本候选。"
            ),
        })
    elif assigned_count and not method_confirmed:
        summary.update({
            "blank_pairing_status": "role_and_polarity_labeled_method_unverified",
            "unverified_assigned_blank_count": assigned_count,
            "blank_pairing_warning": (
                "空白角色与极性已标注，但尚未确认与样品使用同批次／同色谱方法；"
                "因此只作背景比较，不执行扣除。"
            ),
        })
    elif assigned_count and not exclusion_confirmed:
        summary.update({
            "blank_pairing_status": (
                "method_comparable_background_exclusion_not_authorized"
            ),
            "unverified_assigned_blank_count": 0,
            "blank_pairing_warning": (
                "方法可比性已确认，但未授权背景硬排除；空白仅作比较。"
            ),
        })
    elif assigned_count:
        summary.update({
            "blank_pairing_status": "explicitly_confirmed_background_pairing",
            "unverified_assigned_blank_count": 0,
            "blank_pairing_warning": "",
        })
    summary.update({
        "method_comparability_verified": method_confirmed,
        "background_exclusion_authorized": exclusion_confirmed,
        "background_exclusion_eligible": exclusion_eligible,
        "blank_operation": (
            "comparison_only"
            if policy == "soft_evidence"
            else (
                "confirmed_presearch_background_exclusion"
                if exclusion_eligible
                else "hard_exclusion_denied"
            )
        ),
        "cross_injection_msms_subtraction_performed": False,
        "selected_blank_files_analyzed_as_samples": False,
        "background_exclusion_performed": bool(
            policy == "hard_legacy"
            and exclusion_eligible
            and int(summary.get("excluded_spectra") or 0) > 0
        ),
    })
    task["blank_summary"] = summary


def _refresh_completed_task_presentation(
    task: dict, prepared_results: list[dict]
) -> None:
    """Synchronize cached task cards with the current honesty calibration."""
    if task.get("status") != "completed":
        return
    _refresh_blank_pairing_semantics(task)
    matched = [
        row for row in prepared_results
        if isinstance(row, dict) and isinstance(row.get("match"), dict)
    ]
    deduped = _dedup_matched(matched)
    task.setdefault("legacy_unique_library_candidates", len(deduped))
    task["unique_compounds"] = len(deduped)
    task["compound_confidence_summary"] = _compound_confidence_summary(matched)
    # Compound-level High counters must be rebuilt from the same current
    # honesty projection as the row table.  Otherwise an old 1.3 sidecar can
    # report zero/one exact High while the current rows correctly support
    # several unresolved High isomer groups.  This is a presentation-sidecar
    # migration only; the immutable retrieval/ranking artifact is untouched.
    task_id = str(task.get("id") or "").strip()
    source_artifact = TASK_DIR / f"{task_id}_results.json"
    conclusion_artifact = _compound_conclusion_artifact_path(
        TASK_DIR, task_id
    )
    primary_artifact = _flagship_artifact_path(TASK_DIR, task_id)
    if task_id and source_artifact.is_file():
        try:
            feature_inventory, feature_inventory_source = (
                _task_feature_inventory(task)
            )
            if not _compound_conclusion_artifact_is_current(
                conclusion_artifact,
                source_artifact,
                primary_artifact,
                feature_inventory_source,
            ):
                if not _flagship_artifact_is_current(primary_artifact):
                    raise ValueError(
                        "Current V15 primary identity artifact requires materialization"
                    )
                task["compound_conclusion_summary"] = (
                    _write_compound_conclusion_artifact(
                        task_id,
                        prepared_results,
                        conclusion_artifact,
                        claim_decisions_path=_claim_contract_path(
                            TASK_DIR, task_id
                        ),
                        source_artifact=source_artifact,
                        primary_identity_artifact=primary_artifact,
                        feature_inventory=feature_inventory,
                        feature_inventory_source=feature_inventory_source,
                    )
                )
            else:
                current_summary = dict(
                    _read_compound_conclusion_artifact(
                        conclusion_artifact, offset=0, limit=1
                    ).get("metadata") or {}
                )
                current_summary.update(
                    {
                        "status": "available",
                        "artifact": str(conclusion_artifact.resolve()),
                        "size_bytes": conclusion_artifact.stat().st_size,
                    }
                )
                task["compound_conclusion_summary"] = current_summary
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
            logger.warning(
                "Could not refresh compound conclusion projection for %s: %s",
                task_id,
                exc,
            )
    _apply_final_compound_count_summary(
        task, task.get("compound_conclusion_summary")
    )
    try:
        _refresh_reference_standard_qc_projection(task, prepared_results)
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        logger.warning(
            "Could not refresh reference-standard QC projection for %s: %s",
            task_id,
            exc,
        )
    task["evidence_governance_summary"] = (
        _summarize_evidence_governance(prepared_results)
    )
    audit_refresh_needed = (
        task.get("presentation_calibration_version")
        != PRESENTATION_CALIBRATION_VERSION
        or not isinstance(task.get("identification_funnel"), dict)
        or not task.get("identification_funnel")
        or not isinstance(task.get("deconvolution_impact"), dict)
        or not task.get("deconvolution_impact")
    )
    if audit_refresh_needed:
        task["identification_funnel"] = _identification_funnel_summary(
            prepared_results
        )
        task["deconvolution_impact"] = _deconvolution_impact_summary(
            prepared_results
        )
    task["contaminant_unique_compounds"] = sum(
        1 for row in deduped
        if row.get("match", {}).get("potential_contaminant") is True
    )
    feature_summary = _precursor_feature_summary(prepared_results)
    task["feature_summary"] = feature_summary
    task["precursor_features"] = feature_summary["precursor_features"]
    task["matched_precursor_features"] = feature_summary[
        "matched_precursor_features"
    ]
    task["high_precursor_features"] = feature_summary[
        "high_precursor_features"
    ]
    task["ion_identity_summary"] = _effective_ion_identity_summary(
        task.get("ion_identity_summary")
    )
    task["module_availability"] = _analysis_module_availability(
        task.get("params"),
        task.get("ion_identity_summary"),
        prepared_results,
    )
    task["ms1_scans"] = int(
        task["module_availability"].get("ms1_scan_count") or 0
    )
    task["analysis_warnings"] = task["module_availability"]["warnings"]
    task["presentation_calibration_version"] = (
        PRESENTATION_CALIBRATION_VERSION
    )
    task["message"] = _task_message_for_display(task)


def _refresh_task_presentation_from_prepared(
    task_id: str,
    prepared_results: list[dict],
) -> dict:
    """Refresh one task's metadata without caching its heavyweight rows.

    ``prepared_results`` is request-local.  Only fields actually changed by
    the presentation migration are merged back into the shared task registry;
    unrelated background-worker metadata is left untouched.
    """
    with _presentation_refresh_lock:
        with _tasks_lock:
            current = _tasks.get(task_id)
            if not current:
                raise HTTPException(status_code=404, detail="Task not found")
            if (
                current.get("status") != "completed"
                or current.get("presentation_calibration_version")
                == PRESENTATION_CALIBRATION_VERSION
            ):
                return current
            task_view = dict(current)

        before = dict(task_view)
        _refresh_completed_task_presentation(task_view, prepared_results)
        updates = {
            key: value
            for key, value in task_view.items()
            if key != "results" and before.get(key) != value
        }

        with _tasks_lock:
            current = _tasks.get(task_id)
            if not current:
                raise HTTPException(status_code=404, detail="Task not found")
            current.update(updates)
            # A normal production task always has this artifact.  Ensure an
            # older request cannot reintroduce a disk result into RAM.
            if (TASK_DIR / f"{task_id}_results.json").is_file():
                current["results"] = []
            refreshed = current
        _save_tasks()
        return refreshed


def _ensure_task_presentation_current(task_id: str) -> dict:
    """Upgrade cached/in-flight rows, but keep status polling disk-light.

    Completed production rows are deliberately not loaded from disk here.
    Status, module-effect and task-card polling must remain O(metadata).  The
    first result/report request already has the rows locally and performs the
    migration through ``_refresh_task_presentation_from_prepared``.
    """
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if (
        task.get("status") != "completed"
        or task.get("presentation_calibration_version")
        == PRESENTATION_CALIBRATION_VERSION
    ):
        return task

    cached_results = task.get("results") or []
    if not cached_results:
        return task
    prepared_results = _project_claims_for_task(task_id, cached_results)
    return _refresh_task_presentation_from_prepared(task_id, prepared_results)


@app.get("/api/tasks/search")
async def search_tasks(q: str = Query(..., min_length=1, description="Search query (task ID or filename)")):
    """Search tasks by ID or filename (fuzzy match)."""
    with _tasks_lock:
        tasks_list = list(_tasks.values())
    q_lower = q.lower()
    results = [
        t for t in tasks_list
        if q_lower in t["id"].lower() or q_lower in t.get("filename", "").lower()
    ]
    results.sort(key=lambda t: t.get("upload_time", ""), reverse=True)
    summaries = []
    for t in results:
        summaries.append({
            "id": t["id"], "filename": t["filename"], "status": t["status"],
            "progress": t["progress"], "message": _task_message_for_display(t),
            "state_version": t.get("state_version", 1),
            "run_attempt": t.get("run_attempt", 0),
            "active_run_id": t.get("active_run_id"),
            "result_run_id": t.get("result_run_id"),
            "upload_time": t["upload_time"],
            "total_spectra": t.get("total_spectra", 0),
            "ms1_scans": t.get("ms1_scans"),
            "matched_spectra": t.get("matched_spectra", 0),
            "precursor_features": t.get("precursor_features"),
            "matched_precursor_features": t.get("matched_precursor_features"),
            "unique_compounds": t.get("unique_compounds"),
            "legacy_unique_library_candidates": t.get(
                "legacy_unique_library_candidates"
            ),
            "main_result_count": t.get("main_result_count"),
            "result_record_count": t.get("result_record_count"),
            "unresolved_precursor_event_count": t.get(
                "unresolved_precursor_event_count"
            ),
            "isomer_group_result_count": t.get("isomer_group_result_count"),
            "compound_count_semantics": t.get("compound_count_semantics"),
            "compound_confidence_summary": t.get(
                "compound_confidence_summary", {}
            ),
            "identification_funnel": t.get("identification_funnel", {}),
            "deconvolution_impact": t.get("deconvolution_impact", {}),
            "evidence_governance_summary": t.get(
                "evidence_governance_summary", {}
            ),
            "identification_contract_summary": t.get(
                "identification_contract_summary", {}
            ),
            "flagship_summary": t.get("flagship_summary", {}),
            "evirank_mv_summary": t.get("evirank_mv_summary", {}),
            "entity_bundle_summary": t.get("entity_bundle_summary", {}),
            "claim_contract_summary": t.get("claim_contract_summary", {}),
            "compound_conclusion_summary": t.get(
                "compound_conclusion_summary", {}
            ),
            "reference_standard_qc_summary": t.get(
                "reference_standard_qc_summary", {}
            ),
            "contaminant_unique_compounds": t.get("contaminant_unique_compounds"),
            "module_availability": t.get("module_availability", {}),
            "analysis_warnings": t.get("analysis_warnings", []),
            "expected_compounds_summary": t.get(
                "expected_compounds_summary", {}
            ),
            "ion_identity_networks": t.get("ion_identity_networks"),
            "chimeric_spectra": t.get("chimeric_spectra"),
            "deconvolved_spectra": t.get("deconvolved_spectra"),
            "match_rate": round(t.get("match_rate", 0), 1),
            "file_size_mb": t.get("file_size_mb", 0),
        })
    return {"tasks": summaries, "count": len(summaries)}


@app.get("/api/tasks/queue")
async def get_task_queue():
    """Return the live Full-library running item and waiting order."""
    return _full_queue_public_payload()


@app.post("/api/tasks/queue/reorder")
async def reorder_task_queue(request: TaskQueueReorderRequest):
    """Move one not-yet-started Full-library task within the serial queue."""
    with _tasks_lock:
        task = _tasks.get(request.task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")
        if task.get("status") != TaskState.QUEUED.value:
            raise HTTPException(
                status_code=409,
                detail="Only a waiting task can be reordered",
            )
        if str((task.get("params") or {}).get("db_mode") or "full") != "full":
            raise HTTPException(
                status_code=409,
                detail="Only the serial Full-library queue supports reordering",
            )

    reorder = getattr(_executor_full, "reorder", None)
    if not callable(reorder):
        raise HTTPException(
            status_code=503,
            detail="Queue reordering is unavailable in this runtime",
        )
    try:
        snapshot = dict(reorder(request.task_id, request.action) or {})
    except KeyError as exc:
        raise HTTPException(
            status_code=409,
            detail="Task has already started or left the waiting queue",
        ) from exc
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    snapshot["reorder_supported"] = True
    _sync_full_queue_order_metadata(snapshot)
    if _save_tasks() is False:
        raise HTTPException(
            status_code=500,
            detail="Queue changed in memory but its restart order could not be saved",
        )
    response = _full_queue_public_payload(snapshot)
    event_payload = {
        "task_id": request.task_id,
        "action": request.action,
        "order": response["order"],
        "queued_count": response["queued_count"],
    }
    _sse_broadcast("queue_reordered", event_payload)
    logger.info(
        "Full analysis queue reordered: task=%s action=%s order=%s",
        request.task_id,
        request.action,
        response["order"],
    )
    return response


# ── SSE (MUST be before {task_id} to avoid route capture) ──────────

@app.get("/api/tasks/events")
async def task_events():
    """SSE endpoint for real-time task status updates."""
    from fastapi.responses import StreamingResponse

    queue: asyncio.Queue = asyncio.Queue()

    with _sse_lock:
        _sse_subscribers.append(queue)

    async def event_generator():
        try:
            while True:
                try:
                    # Wait for events with timeout (keepalive every 30s)
                    data = await asyncio.wait_for(queue.get(), timeout=30.0)
                    yield f"event: task_update\ndata: {data}\n\n"
                except asyncio.TimeoutError:
                    # Send keepalive
                    yield f": keepalive\n\n"
        except asyncio.CancelledError:
            pass
        finally:
            with _sse_lock:
                try:
                    _sse_subscribers.remove(queue)
                except ValueError:
                    pass

    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no",
        },
    )


@app.get("/api/tasks/{task_id}/isomer-literature")
def get_task_isomer_literature(task_id: str):
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    artifact_path = TASK_DIR / f"{task_id}_isomer_literature.json"
    if artifact_path.exists():
        try:
            return json.loads(artifact_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Literature artifact is unreadable: {exc}",
            ) from exc
    return {
        "task_id": task_id,
        "status": (
            (task.get("isomer_literature_summary") or {}).get("status")
            or "not_started"
        ),
        "summary": task.get("isomer_literature_summary") or {},
        "used_for_scoring": False,
    }


@app.get("/api/expected-compounds/datasets")
async def get_expected_compound_datasets():
    """List target catalogues available to the independent MS1 workflow."""

    try:
        datasets = await asyncio.to_thread(list_expected_datasets, DB_PATH)
    except Exception as exc:
        logger.exception("Could not read Expected Compounds catalogues")
        raise HTTPException(
            status_code=503,
            detail=f"Expected Compounds catalogues are unavailable: {exc}",
        ) from exc
    return {
        "datasets": datasets,
        "dataset_count": len(datasets),
        "target_rows": sum(
            int(item.get("target_rows") or 0) for item in datasets
        ),
        "selection_semantics": (
            "explicit_selection_only; empty_selection_runs_no_target_layer"
        ),
        "filename_formula_dataset_map": {},
        "empty_selection_searches_all_catalogues": False,
        "ordinary_library_match_counted": False,
    }


@app.get("/api/reference-standard-qc/sets")
def get_reference_standard_qc_sets():
    """List explicitly selectable, post-hoc mixed-standard QC truth sets."""

    try:
        sets = _list_reference_standard_sets()
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        logger.exception("Reference-standard QC registry is unavailable")
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    return {
        "sets": sets,
        "set_count": len(sets),
        "selection_semantics": "explicit_only_no_filename_inference",
        "used_for_candidate_generation": False,
        "used_for_scoring": False,
        "used_for_ranking": False,
        "can_project_completed_task_without_reanalysis": True,
    }


@app.get("/api/tasks/{task_id}/reference-standard-qc")
def get_task_reference_standard_qc(task_id: str):
    task = _ensure_task_presentation_current(task_id)
    artifact = _reference_standard_qc_artifact_path(TASK_DIR, task_id)
    if artifact.is_file() and task.get("status") == "completed":
        try:
            return _read_reference_standard_qc_artifact(artifact)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Reference-standard QC artifact is unreadable: {exc}",
            ) from exc
    return {
        "task_id": task_id,
        "status": "not_requested",
        "explicit_selection_required": True,
        "filename_inference_used": False,
        "used_for_candidate_generation": False,
        "used_for_scoring": False,
        "used_for_ranking": False,
        "production_results_mutated": False,
    }


@app.post("/api/tasks/{task_id}/reference-standard-qc")
def project_task_reference_standard_qc(
    task_id: str, request: ReferenceStandardQCRequest
):
    """Attach an explicit QC truth projection without rerunning identification."""

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Reference-standard QC can be projected only after analysis completes",
        )
    try:
        set_id = _validate_reference_standard_set_id(request.set_id)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    results = _get_task_results(task_id)
    polarity = str(
        task.get("ion_mode") or (task.get("params") or {}).get("ion_mode") or ""
    ).strip().upper()
    if polarity not in {"POS", "NEG"}:
        observed_modes = {
            str(row.get("ion_mode") or row.get("polarity") or "").strip().upper()
            for row in results
            if str(row.get("ion_mode") or row.get("polarity") or "").strip().upper()
            in {"POS", "NEG"}
        }
        if len(observed_modes) == 1:
            polarity = next(iter(observed_modes))
        else:
            raise HTTPException(
                status_code=422,
                detail="Task polarity is unresolved or mixed; select a single-polarity task",
            )
    payload = _build_reference_standard_qc(
        task_id,
        results,
        set_id=set_id,
        polarity=polarity,
        filename=str(task.get("filename") or ""),
        mass_tolerance_ppm=request.mass_tolerance_ppm,
        selection_source="explicit_posthoc_request",
    )
    summary = _write_reference_standard_qc_artifact(
        _reference_standard_qc_artifact_path(TASK_DIR, task_id), payload
    )
    with _tasks_lock:
        current = _tasks.get(task_id)
        if current is not None:
            current["reference_standard_set_id"] = set_id
            current["reference_standard_qc_summary"] = summary
            current_params = dict(current.get("params") or {})
            current_params["reference_standard_set_id"] = set_id
            current_params["reference_standard_set_selection_source"] = (
                "explicit_posthoc_request"
            )
            current["params"] = current_params
    _save_tasks()
    return payload


@app.get("/api/tasks/{task_id}/expected-compounds")
def get_task_expected_compounds(
    task_id: str,
    supported_only: bool = Query(False),
    evidence_level: Optional[str] = Query(None),
    limit: int = Query(500, ge=1, le=5000),
    offset: int = Query(0, ge=0),
):
    """Return paginated targeted-MS1 results without mixing them into MS2 hits."""

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    try:
        payload = _load_expected_compounds_payload(task_id, task)
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Expected Compounds artifact is unreadable: {exc}",
        ) from exc

    features = [
        feature for feature in payload.get("features", [])
        if isinstance(feature, dict)
    ]
    if supported_only:
        features = [
            feature for feature in features
            if bool(feature.get("cd_pattern_accepted"))
        ]
    if evidence_level:
        wanted = evidence_level.strip().upper()
        features = [
            feature for feature in features
            if str(feature.get("evidence_level") or "").upper() == wanted
        ]
    total = len(features)
    response = {
        key: value for key, value in payload.items() if key != "features"
    }
    summary = dict(response.get("summary") or {})
    response.update({
        "features": features[offset:offset + limit],
        "features_found": int(
            summary.get("detected_feature_count")
            or summary.get("features_found")
            or len(payload.get("features") or [])
        ),
        "supported_features": int(
            summary.get("supported_features")
            or sum(
                bool(feature.get("cd_pattern_accepted"))
                for feature in payload.get("features") or []
                if isinstance(feature, dict)
            )
        ),
        "filtered_total": total,
        "offset": offset,
        "limit": limit,
        "has_more": offset + limit < total,
        "ordinary_library_match_counted": False,
    })
    return response


@app.get("/api/tasks/{task_id}/expected-compounds/groups")
def get_task_expected_compound_groups(
    task_id: str,
    evidence_view: Literal["all", "supported", "mass_only"] = Query("all"),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
):
    """Return compact, server-deduplicated Expected Compounds groups.

    The full per-XIC artifact remains available from the legacy endpoint and
    all export paths.  This endpoint is a presentation projection designed for
    a responsive result viewer.
    """

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    try:
        payload = _load_expected_compounds_payload(task_id, task)
        groups = expected_group_rows(payload, language="zh")
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Expected Compounds artifact is unreadable: {exc}",
        ) from exc

    supported_total = sum(
        int(group.get("supported_feature_count") or 0) > 0
        for group in groups
    )
    mass_only_total = len(groups) - supported_total
    if evidence_view == "supported":
        filtered = [
            group for group in groups
            if int(group.get("supported_feature_count") or 0) > 0
        ]
    elif evidence_view == "mass_only":
        filtered = [
            group for group in groups
            if int(group.get("supported_feature_count") or 0) == 0
        ]
    else:
        filtered = groups

    response = {
        key: value for key, value in payload.items()
        if key not in {"features", "catalog"}
    }
    response.update({
        "task_id": task_id,
        "groups": filtered[offset:offset + limit],
        "group_total": len(groups),
        "filtered_group_total": len(filtered),
        "supported_group_total": supported_total,
        "mass_only_group_total": mass_only_total,
        "offset": offset,
        "limit": limit,
        "has_more": offset + limit < len(filtered),
        "projection": "deduplicated_group_summary",
        "ordinary_library_match_counted": False,
    })
    return response


@app.post("/api/tasks/{task_id}/isomer-literature/enrich")
def launch_task_isomer_literature(task_id: str):
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Literature enrichment starts only after analysis completes",
        )
    launched = _maybe_launch_isomer_literature_enrichment(task_id)
    with _tasks_lock:
        summary = dict(
            (_tasks.get(task_id) or {}).get(
                "isomer_literature_summary"
            ) or {}
        )
    return {
        "task_id": task_id,
        "scheduled": launched,
        "summary": summary,
        "used_for_scoring": False,
    }


@app.get("/api/tasks/{task_id}")
async def get_task(task_id: str):
    """Get task metadata without embedding the heavyweight result artifact."""
    task = await asyncio.to_thread(_ensure_task_presentation_current, task_id)

    # Return a copy without embedding full results in task list API
    # Frontend uses /tasks/{id}/results for paginated access
    public_task = dict(task)
    public_task.pop("results", None)
    public_task.pop("natural_unique_compounds", None)
    public_task.pop("unclassified_unique_compounds", None)
    return public_task


@app.get("/api/tasks/{task_id}/status")
async def get_task_status(task_id: str):
    """Get task status only (lightweight, for polling)."""
    task = await asyncio.to_thread(_ensure_task_presentation_current, task_id)

    return {
        "id": task["id"],
        "filename": task.get("filename", task["id"]),
        "upload_time": task.get("upload_time"),
        "start_time": task.get("start_time"),
        "end_time": task.get("end_time"),
        "status": task["status"],
        "progress": task["progress"],
        "state_version": task.get("state_version", 1),
        "run_attempt": task.get("run_attempt", 0),
        "active_run_id": task.get("active_run_id"),
        "result_run_id": task.get("result_run_id"),
        "can_cancel": task.get("status") in {
            TaskState.PENDING.value,
            TaskState.QUEUED.value,
            TaskState.RUNNING.value,
        },
        "message": _task_message_for_display(task),
        "params": task.get("params", {}),
        "ion_mode": task.get("ion_mode"),
        "ion_mode_counts": task.get("ion_mode_counts", {}),
        "blank_summary": task.get("blank_summary", {}),
        "spectrum_representation_processing": task.get(
            "spectrum_representation_processing", {}
        ),
        "matched_spectra": task.get("matched_spectra", 0),
        "precursor_features": task.get("precursor_features"),
        "matched_precursor_features": task.get("matched_precursor_features"),
        "high_precursor_features": task.get("high_precursor_features"),
        "unique_compounds": task.get("unique_compounds"),
        "legacy_unique_library_candidates": task.get(
            "legacy_unique_library_candidates"
        ),
        "main_result_count": task.get("main_result_count"),
        "result_record_count": task.get("result_record_count"),
        "unresolved_precursor_event_count": task.get(
            "unresolved_precursor_event_count"
        ),
        "isomer_group_result_count": task.get("isomer_group_result_count"),
        "compound_count_semantics": task.get("compound_count_semantics"),
        "compound_confidence_summary": task.get(
            "compound_confidence_summary", {}
        ),
        "identification_funnel": task.get("identification_funnel", {}),
        "deconvolution_impact": task.get("deconvolution_impact", {}),
        "evidence_governance_summary": task.get(
            "evidence_governance_summary", {}
        ),
        "identification_contract_summary": task.get(
            "identification_contract_summary", {}
        ),
        "flagship_summary": task.get("flagship_summary", {}),
        "analysis_completion_check": task.get("analysis_completion_check", {}),
        "primary_dependency_preflight": task.get("primary_dependency_preflight", {}),
        "evirank_mv_summary": task.get("evirank_mv_summary", {}),
        "entity_bundle_summary": task.get("entity_bundle_summary", {}),
        "claim_contract_summary": task.get("claim_contract_summary", {}),
        "compound_conclusion_summary": task.get(
            "compound_conclusion_summary", {}
        ),
        "reference_standard_qc_summary": task.get(
            "reference_standard_qc_summary", {}
        ),
        "contaminant_unique_compounds": task.get("contaminant_unique_compounds"),
        "module_availability": task.get("module_availability", {}),
        "analysis_warnings": task.get("analysis_warnings", []),
        "peak_processing_summary": task.get("peak_processing_summary", {}),
        "untargeted_workflow": task.get("untargeted_workflow", {}),
        "consensus_spectrum_summary": task.get(
            "consensus_spectrum_summary", {}
        ),
        "expected_compounds_summary": task.get(
            "expected_compounds_summary", {}
        ),
        "reliability_annotation_summary": task.get(
            "reliability_annotation_summary", {}
        ),
        "pipeline_boundary_audit": task.get(
            "pipeline_boundary_audit", {}
        ),
        "ion_identity_networks": task.get("ion_identity_networks"),
        "ion_identity_networked_features": task.get("ion_identity_networked_features"),
        "ion_identity_linked_annotations": task.get("ion_identity_linked_annotations"),
        "chimeric_spectra": task.get("chimeric_spectra"),
        "deconvolved_spectra": task.get("deconvolved_spectra"),
        "isomer_literature_summary": task.get(
            "isomer_literature_summary", {}
        ),
        "evograph_fragmentation_summary": task.get(
            "evograph_fragmentation_summary", {}
        ),
        "ms1_scans": task.get("ms1_scans"),
        "total_spectra": task.get("total_spectra", 0),
        "error": task.get("error"),
    }


@app.post("/api/tasks/batch/analyze")
async def batch_start_analysis(request: BatchAnalysisRequest):
    """Queue multiple tasks with one validated parameter set.

    FULL-library work remains memory-safe because its executor is intentionally
    serial; SMALL-library work uses the existing bounded worker pool.
    """
    task_ids = list(dict.fromkeys(request.task_ids))
    if not task_ids:
        raise HTTPException(status_code=400, detail="Select at least one task")
    if len(task_ids) > 100:
        raise HTTPException(status_code=400, detail="At most 100 tasks per batch")

    queued, skipped, failed = [], [], []
    for task_id in task_ids:
        with _tasks_lock:
            task = _tasks.get(task_id)
            status = task.get("status") if task else None
        if _is_active_task_state(status):
            skipped.append({"id": task_id, "reason": f"already {status}"})
            continue
        try:
            await start_analysis(task_id, request.params)
            queued.append(task_id)
        except HTTPException as exc:
            failed.append({"id": task_id, "error": str(exc.detail)})
        except Exception as exc:
            logger.exception("Batch queue failed for task %s", task_id)
            failed.append({"id": task_id, "error": str(exc) or exc.__class__.__name__})

    return {
        "queued": queued,
        "queued_count": len(queued),
        "skipped": skipped,
        "failed": failed,
    }


@app.post("/api/tasks/{task_id}/analyze")
async def start_analysis(task_id: str, params: AnalysisParams = AnalysisParams()):
    """Start or re-run analysis for a task."""
    global _lib_records, _mz_index

    requested_ion_mode = str(params.ion_mode or "auto").strip().upper()
    if requested_ion_mode not in {"AUTO", "POS", "NEG"}:
        raise HTTPException(status_code=422, detail="ion_mode must be auto, POS, or NEG")

    # Resolve only an explicitly selected presence-only target catalogue before
    # claiming or clearing a task. Sample names and source identities never
    # select a catalogue in the production path.
    with _tasks_lock:
        task_snapshot = _tasks.get(task_id)
        if not task_snapshot:
            raise HTTPException(status_code=404, detail="Task not found")
        target_filename = str(
            task_snapshot.get("filename")
            or Path(str(task_snapshot.get("file_path") or "")).name
        )
    expected_dataset_selection: Optional[dict[str, Any]] = None
    if params.enable_expected_compounds:
        try:
            available_target_datasets = await asyncio.to_thread(
                list_expected_datasets, DB_PATH
            )
        except Exception as exc:
            logger.exception(
                "Could not resolve Expected Compounds dataset for %s", task_id
            )
            raise HTTPException(
                status_code=503,
                detail=f"Expected Compounds catalogues are unavailable: {exc}",
            ) from exc
        expected_dataset_selection = resolve_expected_dataset_selection(
            target_filename,
            params.expected_compound_datasets,
            available_target_datasets,
        )
        unknown_sources = list(
            expected_dataset_selection.get("unknown_requested_sources") or []
        )
        if unknown_sources:
            raise HTTPException(
                status_code=422,
                detail=(
                    "Unknown Expected Compounds dataset(s): "
                    + ", ".join(unknown_sources)
                ),
            )
        if not expected_dataset_selection.get("selected_dataset_sources"):
            raise HTTPException(
                status_code=422,
                detail=(
                    "Expected Compounds requires an explicit target dataset; "
                    "sample filenames never select one automatically"
                ),
            )

    with _library_state_lock:
        if not _lib_loaded:
            raise HTTPException(
                status_code=400,
                detail="Library not loaded. Load it first via /api/library/load",
            )

        actual_db_mode = "small" if str(_lib_name).startswith("Small") else "full"
        if params.db_mode != actual_db_mode:
            raise HTTPException(
                status_code=409,
                detail=(
                    f"Analysis requested {params.db_mode!r}, but the loaded library is "
                    f"{actual_db_mode!r}. Load the selected library before starting analysis."
                ),
            )

        with _tasks_lock:
            task = _tasks.get(task_id)
            if not task:
                raise HTTPException(status_code=404, detail="Task not found")
            if _is_active_task_state(task["status"]):
                raise HTTPException(
                    status_code=409,
                    detail=f"Task is already {task['status']}",
                )

            # Claim the task while the library lock is still held.  A
            # concurrent import will now see "queued" and cannot invalidate
            # the arrays that this analysis is about to use.
            run_id = uuid.uuid4().hex[:16]
            task["active_run_id"] = run_id
            _transition_task(
                task,
                TaskState.QUEUED,
                reason="user_requested_analysis",
                run_id=run_id,
            )
            _task_cancel_events[task_id] = threading.Event()
            task["progress"] = 0
            task.pop("queue_order", None)
            task["message"] = (
                "Queued for Full-library serial analysis..."
                if actual_db_mode == "full"
                else "Queued for analysis..."
            )
            task["start_time"] = None
            task["end_time"] = None
            task["total_spectra"] = 0
            task["ms1_scans"] = None
            task["matched_spectra"] = 0
            task["unique_compounds"] = None
            task["compound_confidence_summary"] = {}
            task["identification_funnel"] = {}
            task["deconvolution_impact"] = {}
            task["evidence_governance_summary"] = {}
            task["identification_contract_summary"] = {}
            task["flagship_summary"] = {}
            task["evirank_mv_summary"] = {}
            task["entity_bundle_summary"] = {}
            task["claim_contract_summary"] = {}
            task["compound_conclusion_summary"] = {}
            task["reference_standard_qc_summary"] = {}
            task["evograph_fragmentation_summary"] = {}
            task["match_rate"] = 0.0
            task["results"] = []
            task["error"] = None
            params_payload = _analysis_params_with_task_provenance(
                task,
                params,
                expected_dataset_selection=expected_dataset_selection,
            )
            task["params"] = params_payload
            task["reference_standard_set_id"] = None
            task["expected_compound_datasets"] = list(
                params_payload.get("expected_compound_datasets") or []
            )
            task["expected_compound_dataset_selection_source"] = str(
                params_payload.get(
                    "expected_compound_dataset_selection_source"
                ) or "none"
            )
            task["expected_compound_dataset_selection_status"] = str(
                params_payload.get(
                    "expected_compound_dataset_selection_status"
                ) or ""
            )
            task["isolation_window_width_source"] = params_payload.get(
                "isolation_window_width_source"
            )
            task["blank_summary"] = {}
            task["ion_identity_summary"] = {}
            task["expected_compounds_summary"] = {}
            task["reliability_annotation_summary"] = {}
            task["pipeline_boundary_audit"] = {
                "expected_enabled": bool(
                    params_payload.get("enable_expected_compounds", False)
                ),
                "target_screening_enabled": bool(
                    params_payload.get("enable_expected_compounds", False)
                ),
                "target_compounds_ref": (
                    list(params_payload.get("expected_compound_datasets") or [])
                    or None
                ),
                "expected_provider_call_count": 0,
                "evaluation_panel_access_by_pipeline": 0,
            }
            task["module_availability"] = {}
            task["analysis_warnings"] = []
            task["ion_identity_networks"] = 0
            task["ion_identity_networked_features"] = 0
            task["ion_identity_linked_annotations"] = 0
            task["chimeric_spectra"] = 0
            task["deconvolved_spectra"] = 0
    _save_tasks()
    try:
        _expected_compounds_artifact_path(task_id).unlink(missing_ok=True)
    except OSError as exc:
        logger.warning(
            "Could not clear stale Expected Compounds artifact for %s: %s",
            task_id,
            exc,
        )
    try:
        _evograph_sidecar_path(TASK_DIR, task_id).unlink(missing_ok=True)
    except OSError as exc:
        logger.warning(
            "Could not clear stale EvoGraph-F artifact for %s: %s",
            task_id,
            exc,
        )
    try:
        _evirank_artifact_path(TASK_DIR, task_id).unlink(missing_ok=True)
    except OSError as exc:
        logger.warning(
            "Could not clear stale V12 EviRank-MV artifact for %s: %s",
            task_id,
            exc,
        )
    try:
        _evitwin_artifact_path(task_id).unlink(missing_ok=True)
    except OSError as exc:
        logger.warning(
            "Could not clear stale EviTwin-MS artifact for %s: %s",
            task_id,
            exc,
        )
    try:
        _compound_conclusion_artifact_path(TASK_DIR, task_id).unlink(
            missing_ok=True
        )
    except OSError as exc:
        logger.warning(
            "Could not clear stale compound conclusion artifact for %s: %s",
            task_id,
            exc,
        )

    # Check file exists
    if not Path(task["file_path"]).exists():
        with _tasks_lock:
            current = _tasks[task_id]
            _transition_task(
                current,
                TaskState.FAILED,
                reason="uploaded_file_missing_before_queue_submit",
            )
            current["error"] = "Uploaded file no longer exists"
            _task_cancel_events.pop(task_id, None)
        _save_tasks()
        raise HTTPException(status_code=404, detail="Uploaded file not found")

    # Submit to executor
    executor = _executor_full if actual_db_mode == "full" else _executor_small
    try:
        executor.submit(_run_analysis, task_id)
    except Exception as exc:
        with _tasks_lock:
            current = _tasks[task_id]
            _transition_task(
                current,
                TaskState.FAILED,
                reason="executor_queue_submission_failed",
            )
            current["error"] = f"Could not queue analysis: {exc}"
            _task_cancel_events.pop(task_id, None)
        _save_tasks()
        raise HTTPException(status_code=503, detail="Could not queue analysis") from exc

    if actual_db_mode == "full":
        queue_snapshot = _sync_full_queue_order_metadata()
        _save_tasks()
        queue_payload = _full_queue_public_payload(queue_snapshot)
        _sse_broadcast(
            "queue_reordered",
            {
                "task_id": task_id,
                "action": "enqueue",
                "order": queue_payload["order"],
                "queued_count": queue_payload["queued_count"],
            },
        )

    logger.info(
        "Analysis queued for task %s (%s; executor=%s)",
        task_id,
        task["filename"],
        "full-serial" if actual_db_mode == "full" else "small-bounded",
    )
    return {"task_id": task_id, "status": "queued"}


def _write_runtime_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(
        f"{path.name}.{os.getpid()}.{threading.get_ident()}.tmp"
    )
    try:
        with open(partial, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, default=str)
            handle.flush()
            os.fsync(handle.fileno())
        partial.replace(path)
    finally:
        partial.unlink(missing_ok=True)


def _task_owned_runtime_paths(task_id: str, task: dict[str, Any]) -> list[Path]:
    """Return managed task paths without ever claiming arbitrary user paths."""
    candidates: list[Path] = []
    if task.get("file_path"):
        candidates.append(Path(str(task["file_path"])))
    if TASK_DIR.is_dir():
        candidates.extend(TASK_DIR.glob(f"{task_id}_*"))
        exact = TASK_DIR / task_id
        if exact.exists() or exact.is_symlink():
            candidates.append(exact)
    unique: list[Path] = []
    seen: set[str] = set()
    for path in candidates:
        token = str(path)
        if token not in seen:
            unique.append(path)
            seen.add(token)
    return unique


def _managed_output_relative(path: Path) -> Path | None:
    try:
        root = OUTPUT_DIR.resolve()
        resolved = path.resolve(strict=False)
        relative = resolved.relative_to(root)
    except (OSError, RuntimeError, ValueError):
        return None
    if relative.parts and relative.parts[0] == "task_recycle_bin":
        return None
    return relative


def _archive_deleted_task_files_unlocked(
    task_id: str,
    task: dict[str, Any],
    *,
    deletion_id: str,
) -> dict[str, Any]:
    """Move task-owned runtime data to a recoverable recycle generation."""
    recycle_root = OUTPUT_DIR / "task_recycle_bin" / deletion_id / task_id
    moved: list[dict[str, str]] = []
    skipped: list[dict[str, str]] = []
    failures: list[dict[str, str]] = []
    for source in _task_owned_runtime_paths(task_id, task):
        if not source.exists() and not source.is_symlink():
            continue
        relative = _managed_output_relative(source)
        if relative is None:
            skipped.append({
                "path": str(source),
                "reason": "outside_managed_runtime_root_not_touched",
            })
            continue
        destination = recycle_root / "files" / relative
        try:
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists() or destination.is_symlink():
                raise FileExistsError(f"recycle destination exists: {destination}")
            shutil.move(str(source), str(destination))
            moved.append({
                "original_path": str(source),
                "recycle_path": str(destination),
            })
        except Exception as exc:
            failures.append({"path": str(source), "error": str(exc)})
            logger.error(
                "Could not archive task-owned path for %s (%s): %s",
                task_id, source, exc,
            )

    _write_runtime_json_atomic(recycle_root / "task.json", task)
    manifest = {
        "schema_version": 1,
        "deletion_id": deletion_id,
        "task_id": task_id,
        "deleted_at": _now_iso(),
        "recoverable": not failures,
        "moved": moved,
        "skipped": skipped,
        "failures": failures,
        "restored_at": None,
    }
    _write_runtime_json_atomic(recycle_root / "manifest.json", manifest)
    return manifest


def _archive_deleted_task_files(
    task_id: str,
    task: dict[str, Any],
    *,
    deletion_id: str,
) -> dict[str, Any]:
    """Archive one task while excluding restore/purge races."""
    with _task_recycle_lock:
        return _archive_deleted_task_files_unlocked(
            task_id,
            task,
            deletion_id=deletion_id,
        )


def _request_audit_context(request: Request) -> dict[str, Any]:
    return {
        "client_host": request.client.host if request.client else None,
        "user_agent": request.headers.get("user-agent", "")[:500],
    }


def _append_task_deletion_audit(payload: dict[str, Any]) -> None:
    path = OUTPUT_DIR / "runtime" / "task_deletion_audit.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(payload, ensure_ascii=False, default=str) + "\n"
    with _task_deletion_audit_lock:
        with open(path, "a", encoding="utf-8") as handle:
            handle.write(line)
            handle.flush()
            os.fsync(handle.fileno())


def _task_recycle_tree_stats(root: Path) -> dict[str, Any]:
    """Measure one recycle tree without following symbolic links."""
    size_bytes = 0
    file_count = 0
    entry_count = 0
    symlinks: list[str] = []
    signature_rows: list[tuple[str, str, int, int]] = []
    stack = [root]
    while stack:
        current = stack.pop()
        try:
            stat = current.lstat()
        except FileNotFoundError:
            continue
        entry_count += 1
        try:
            relative = str(current.relative_to(root)) or "."
        except ValueError:
            relative = str(current)
        if current.is_symlink():
            symlinks.append(str(current))
            signature_rows.append(
                (relative, "link", int(stat.st_size), int(stat.st_mtime_ns))
            )
            continue
        if current.is_file():
            size_bytes += int(stat.st_size)
            file_count += 1
            signature_rows.append(
                (relative, "file", int(stat.st_size), int(stat.st_mtime_ns))
            )
        elif current.is_dir():
            signature_rows.append(
                (relative, "dir", int(stat.st_size), int(stat.st_mtime_ns))
            )
            try:
                with os.scandir(current) as iterator:
                    stack.extend(Path(item.path) for item in iterator)
            except FileNotFoundError:
                continue
    return {
        "size_bytes": size_bytes,
        "file_count": file_count,
        "entry_count": entry_count,
        "symlinks": symlinks,
        "tree_signature": hashlib.sha256(
            json.dumps(
                sorted(signature_rows),
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest(),
    }


def _task_recycle_audit_pairs() -> set[tuple[str, str]]:
    """Read deletion/task pairs from the append-only deletion audit."""
    path = OUTPUT_DIR / "runtime" / "task_deletion_audit.jsonl"
    pairs: set[tuple[str, str]] = set()
    if not path.is_file():
        return pairs
    with _task_deletion_audit_lock:
        try:
            lines = path.read_text(encoding="utf-8").splitlines()
        except OSError:
            return pairs
    for line in lines:
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if row.get("event") not in {
            "single_task_recycled",
            "batch_task_recycled",
        }:
            continue
        deletion_id = str(row.get("deletion_id") or "")
        task_ids: set[str] = set()
        for key in ("task_ids", "deleted_task_ids"):
            task_ids.update(str(value) for value in (row.get(key) or []))
        archives = row.get("archives") or {}
        if isinstance(archives, dict):
            task_ids.update(str(value) for value in archives)
        archive = row.get("archive") or {}
        if isinstance(archive, dict) and archive.get("task_id"):
            task_ids.add(str(archive["task_id"]))
        pairs.update((deletion_id, task_id) for task_id in task_ids)
    return pairs


def _task_recycle_path_is_within(path: Path, parent: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(parent.resolve(strict=False))
        return True
    except (OSError, RuntimeError, ValueError):
        return False


def _task_recycle_inventory_unlocked() -> dict[str, Any]:
    """Return a purge-safe recycle inventory while the recycle lock is held."""
    root = OUTPUT_DIR / "task_recycle_bin"
    entries: list[dict[str, Any]] = []
    structural_errors: list[str] = []
    audit_pairs = _task_recycle_audit_pairs()
    with _tasks_lock:
        current_task_ids = set(_tasks)

    if root.exists() and root.is_symlink():
        structural_errors.append("Recycle root is a symbolic link")
    elif root.is_dir():
        for generation in sorted(root.iterdir(), key=lambda path: path.name):
            if generation.is_symlink() or not generation.is_dir():
                structural_errors.append(
                    f"Unexpected recycle generation entry: {generation.name}"
                )
                continue
            if not re.fullmatch(r"[A-Za-z0-9_-]{1,80}", generation.name):
                structural_errors.append(
                    f"Invalid recycle generation identifier: {generation.name}"
                )
                continue
            for task_root in sorted(
                generation.iterdir(), key=lambda path: path.name
            ):
                if task_root.is_symlink() or not task_root.is_dir():
                    structural_errors.append(
                        f"Unexpected recycled task entry: {task_root.name}"
                    )
                    continue
                deletion_id = generation.name
                task_id = task_root.name
                reasons: list[str] = []
                manifest_path = task_root / "manifest.json"
                task_path = task_root / "task.json"
                manifest: dict[str, Any] = {}
                task: dict[str, Any] = {}

                if not re.fullmatch(r"[A-Za-z0-9_-]{1,80}", task_id):
                    reasons.append("invalid_task_identifier")
                if not manifest_path.is_file() or manifest_path.is_symlink():
                    reasons.append("missing_or_linked_manifest")
                if not task_path.is_file() or task_path.is_symlink():
                    reasons.append("missing_or_linked_task_record")
                if not reasons:
                    try:
                        loaded_manifest = json.loads(
                            manifest_path.read_text(encoding="utf-8")
                        )
                        loaded_task = json.loads(
                            task_path.read_text(encoding="utf-8")
                        )
                        if not isinstance(loaded_manifest, dict):
                            raise ValueError("manifest is not an object")
                        if not isinstance(loaded_task, dict):
                            raise ValueError("task record is not an object")
                        manifest = loaded_manifest
                        task = loaded_task
                    except (OSError, ValueError, json.JSONDecodeError):
                        reasons.append("unreadable_recycle_metadata")

                restored = bool(manifest.get("restored_at"))
                if manifest:
                    if manifest.get("schema_version") != 1:
                        reasons.append("unsupported_manifest_schema")
                    if str(manifest.get("deletion_id") or "") != deletion_id:
                        reasons.append("manifest_deletion_id_mismatch")
                    if str(manifest.get("task_id") or "") != task_id:
                        reasons.append("manifest_task_id_mismatch")
                    if manifest.get("failures"):
                        reasons.append("archive_contains_move_failures")
                if not restored and task_id in current_task_ids:
                    reasons.append("task_id_exists_in_current_registry")
                if (deletion_id, task_id) not in audit_pairs:
                    reasons.append("deletion_pair_missing_from_audit")

                files_root = task_root / "files"
                if not restored:
                    for position, item in enumerate(manifest.get("moved") or []):
                        if not isinstance(item, dict):
                            reasons.append(f"invalid_moved_entry_{position}")
                            continue
                        recycled = Path(str(item.get("recycle_path") or ""))
                        original = Path(str(item.get("original_path") or ""))
                        if not _task_recycle_path_is_within(
                            recycled, files_root
                        ):
                            reasons.append(
                                f"recycled_path_outside_task_root_{position}"
                            )
                        elif not recycled.exists() and not recycled.is_symlink():
                            reasons.append(f"missing_recycled_path_{position}")
                        if _managed_output_relative(original) is None:
                            reasons.append(
                                f"original_path_outside_managed_output_{position}"
                            )

                stats = _task_recycle_tree_stats(task_root)
                if stats["symlinks"]:
                    reasons.append("recycle_tree_contains_symbolic_links")
                manifest_sha256 = None
                if manifest_path.is_file() and not manifest_path.is_symlink():
                    try:
                        manifest_sha256 = hashlib.sha256(
                            manifest_path.read_bytes()
                        ).hexdigest()
                    except OSError:
                        reasons.append("manifest_hash_failed")

                entries.append({
                    "deletion_id": deletion_id,
                    "task_id": task_id,
                    "filename": str(task.get("filename") or task_id),
                    "previous_status": str(task.get("status") or "unknown"),
                    "uploaded_at": task.get("upload_time"),
                    "deleted_at": manifest.get("deleted_at"),
                    "restored_at": manifest.get("restored_at"),
                    "restored": restored,
                    "recoverable": bool(manifest.get("recoverable")),
                    "size_bytes": int(stats["size_bytes"]),
                    "file_count": int(stats["file_count"]),
                    "entry_count": int(stats["entry_count"]),
                    "tree_signature": stats["tree_signature"],
                    "safe_to_purge": not restored and not reasons,
                    "purge_block_reasons": reasons,
                    "manifest_sha256": manifest_sha256,
                })

    entries.sort(
        key=lambda item: str(item.get("deleted_at") or ""), reverse=True
    )
    active = [item for item in entries if not item["restored"]]
    restored = [item for item in entries if item["restored"]]
    unsafe_active = [item for item in active if not item["safe_to_purge"]]
    signature_payload = [
        {
            "deletion_id": item["deletion_id"],
            "task_id": item["task_id"],
            "deleted_at": item["deleted_at"],
            "size_bytes": item["size_bytes"],
            "entry_count": item["entry_count"],
            "tree_signature": item["tree_signature"],
            "manifest_sha256": item["manifest_sha256"],
            "safe_to_purge": item["safe_to_purge"],
        }
        for item in active
    ]
    signature = hashlib.sha256(
        json.dumps(
            {
                "entries": signature_payload,
                "structural_errors": structural_errors,
            },
            ensure_ascii=False,
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()
    return {
        "entries": entries,
        "active_entries": active,
        "restored_entries": restored,
        "unsafe_active_entries": unsafe_active,
        "structural_errors": structural_errors,
        "signature": signature,
        "active_total": len(active),
        "restored_total": len(restored),
        "active_size_bytes": sum(item["size_bytes"] for item in active),
        "active_file_count": sum(item["file_count"] for item in active),
    }


def _task_recycle_inventory() -> dict[str, Any]:
    with _task_recycle_lock:
        return _task_recycle_inventory_unlocked()


def _task_delete_state_signature(task_ids: list[str]) -> str:
    with _tasks_lock:
        rows = [
            {
                "id": task_id,
                "status": (_tasks.get(task_id) or {}).get("status"),
                "state_version": (_tasks.get(task_id) or {}).get("state_version"),
            }
            for task_id in task_ids
        ]
    return hashlib.sha256(
        json.dumps(rows, ensure_ascii=False, sort_keys=True).encode("utf-8")
    ).hexdigest()


@app.post("/api/tasks/batch/delete/preview")
async def preview_batch_delete_tasks(request: BatchTaskDeletePreviewRequest):
    """Issue a short-lived token bound to exact task IDs and state versions."""
    task_ids = request.task_ids
    with _tasks_lock:
        missing = [task_id for task_id in task_ids if task_id not in _tasks]
        total_registered = len(_tasks)
        statuses = Counter(
            str((_tasks.get(task_id) or {}).get("status") or "unknown")
            for task_id in task_ids
        )
    if missing:
        raise HTTPException(
            status_code=404,
            detail=f"Task(s) no longer exist: {', '.join(missing[:5])}",
        )
    token = uuid.uuid4().hex
    expires_at_epoch = time.time() + 120.0
    phrase = f"确认处理 {len(task_ids)} 个任务"
    signature = _task_delete_state_signature(task_ids)
    with _task_delete_preview_lock:
        now = time.time()
        for stale in [
            key for key, value in _task_delete_previews.items()
            if float(value.get("expires_at_epoch") or 0) < now
        ]:
            _task_delete_previews.pop(stale, None)
        _task_delete_previews[token] = {
            "task_ids": task_ids,
            "signature": signature,
            "confirmation_phrase": phrase,
            "expires_at_epoch": expires_at_epoch,
        }
    return {
        "confirmation_token": token,
        "confirmation_phrase": phrase,
        "expires_in_seconds": 120,
        "task_count": len(task_ids),
        "status_counts": dict(statuses),
        "selects_entire_registry": bool(
            total_registered and len(task_ids) == total_registered
        ),
        "typed_confirmation_required": bool(
            len(task_ids) >= 10
            or (total_registered and len(task_ids) == total_registered)
        ),
        "deletion_is_recoverable": True,
    }


def _consume_task_delete_preview(request: BatchTaskDeleteRequest) -> None:
    with _task_delete_preview_lock:
        preview = _task_delete_previews.pop(request.confirmation_token, None)
    if not preview:
        raise HTTPException(
            status_code=409,
            detail="删除确认已失效；请重新预览后再确认",
        )
    if float(preview.get("expires_at_epoch") or 0) < time.time():
        raise HTTPException(status_code=409, detail="删除确认已过期；请重新预览")
    if request.task_ids != preview.get("task_ids"):
        raise HTTPException(status_code=409, detail="任务选择已变化；删除已拒绝")
    if request.confirmation_phrase != preview.get("confirmation_phrase"):
        raise HTTPException(status_code=409, detail="删除确认文字不匹配")
    if _task_delete_state_signature(request.task_ids) != preview.get("signature"):
        raise HTTPException(status_code=409, detail="任务状态已变化；请重新预览")


def _active_delete_conflict_detail(status: str) -> str:
    if status == TaskState.RUNNING.value:
        return "任务正在运行；请先点击“停止”，待状态变为“已停止”后再删除"
    if status == TaskState.CANCELLING.value:
        return "任务正在安全停止；待状态变为“已停止”后再删除"
    if status == TaskState.FINALIZING.value:
        return "任务正在写入结果，当前不能安全删除；请等待写入完成"
    return f"任务当前状态为 {status}，不能安全删除"


@app.delete("/api/tasks/{task_id}")
async def delete_task(task_id: str, request: Request):
    """Remove one safe task from the registry and move its data to recycle."""
    deletion_id = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    with _tasks_lock:
        task = _tasks.get(task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")
        previous_status = str(task.get("status") or "")
        if previous_status in {
            TaskState.RUNNING.value,
            TaskState.CANCELLING.value,
            TaskState.FINALIZING.value,
        }:
            raise HTTPException(
                status_code=409,
                detail=_active_delete_conflict_detail(previous_status),
            )

        # QUEUED is intentionally deletable.  Removal is atomic with the
        # worker's QUEUED -> RUNNING claim: if this lock wins, a later worker
        # sees no task and exits; if the worker wins, status is RUNNING and the
        # request above is refused so no live analysis is torn down.
        _tasks.pop(task_id, None)
        _task_cancel_events.pop(task_id, None)
        _restart_recovery_task_ids.discard(task_id)

    if _save_tasks(allow_empty=True) is False:
        with _tasks_lock:
            _tasks[task_id] = task
        _save_tasks()
        raise HTTPException(
            status_code=500,
            detail="任务注册表保存失败；删除已中止，原始数据未移动",
        )
    if previous_status == TaskState.QUEUED.value:
        _remove_waiting_full_queue_item(task_id)
    with _task_recycle_lock:
        archive = _archive_deleted_task_files(
            task_id, task, deletion_id=deletion_id,
        )
        _append_task_deletion_audit({
            "event": "single_task_recycled",
            "at": _now_iso(),
            "deletion_id": deletion_id,
            "task_ids": [task_id],
            "previous_statuses": {task_id: previous_status},
            "archive": archive,
            **_request_audit_context(request),
        })
    logger.info(
        "Recycled task %s (previous status=%s, deletion_id=%s)",
        task_id, previous_status, deletion_id,
    )
    return {
        "deleted": task_id,
        "previous_status": previous_status,
        "removed_from_queue": previous_status == TaskState.QUEUED.value,
        "deletion_id": deletion_id,
        "recoverable": bool(archive.get("recoverable")),
    }


@app.post("/api/tasks/batch/delete")
async def batch_delete_tasks(
    request_payload: BatchTaskDeleteRequest,
    request: Request,
):
    """Process a state-bound confirmed batch and recycle all safe task data."""
    _consume_task_delete_preview(request_payload)
    task_ids = request_payload.task_ids
    deletion_id = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    deleted: list[str] = []
    stop_requested: list[str] = []
    failed: list[dict[str, str]] = []
    removed_tasks: list[tuple[str, dict[str, Any]]] = []
    status_events: list[dict[str, Any]] = []

    for tid in dict.fromkeys(task_ids):
        try:
            with _tasks_lock:
                task = _tasks.get(tid)
                if not task:
                    failed.append({"id": tid, "error": "Task not found"})
                    continue
                status = str(task.get("status") or "")
                if status == TaskState.RUNNING.value:
                    event = _task_cancel_events.setdefault(tid, threading.Event())
                    event.set()
                    _transition_task(
                        task,
                        TaskState.CANCELLING,
                        reason="batch_delete_requested_safe_stop",
                    )
                    task["message"] = "Cancellation requested; waiting for a safe checkpoint"
                    stop_requested.append(tid)
                    status_events.append({
                        "task_id": tid,
                        "status": task["status"],
                        "progress": task.get("progress", 0),
                        "message": task["message"],
                    })
                    continue
                if status == TaskState.CANCELLING.value:
                    stop_requested.append(tid)
                    continue
                if status == TaskState.FINALIZING.value:
                    failed.append({
                        "id": tid,
                        "error": _active_delete_conflict_detail(status),
                    })
                    continue

                # Includes QUEUED: the lock makes queue removal race-safe in
                # exactly the same way as the single-delete endpoint.
                _tasks.pop(tid, None)
                _task_cancel_events.pop(tid, None)
                _restart_recovery_task_ids.discard(tid)
                removed_tasks.append((tid, task))
            deleted.append(tid)
        except Exception as e:
            failed.append({"id": tid, "error": str(e)})

    if (deleted or stop_requested) and _save_tasks(
        allow_empty=bool(deleted),
    ) is False:
        with _tasks_lock:
            for tid, task in removed_tasks:
                _tasks[tid] = task
        _save_tasks()
        raise HTTPException(
            status_code=500,
            detail="任务注册表保存失败；批量删除已中止，原始数据未移动",
        )
    for tid, task in removed_tasks:
        if task.get("status") == TaskState.QUEUED.value:
            _remove_waiting_full_queue_item(tid)
    archives: dict[str, dict[str, Any]] = {}
    with _task_recycle_lock:
        for tid, task in removed_tasks:
            archives[tid] = _archive_deleted_task_files(
                tid, task, deletion_id=deletion_id,
            )
        _append_task_deletion_audit({
            "event": "batch_task_recycled",
            "at": _now_iso(),
            "deletion_id": deletion_id,
            "requested_task_ids": task_ids,
            "deleted_task_ids": deleted,
            "stop_requested_task_ids": stop_requested,
            "failed": failed,
            "archives": archives,
            **_request_audit_context(request),
        })
    for payload in status_events:
        _sse_broadcast("status_change", payload)

    logger.info(
        "Batch task action: deleted=%d, stop_requested=%d, failed=%d",
        len(deleted),
        len(stop_requested),
        len(failed),
    )
    return {
        "deleted_count": len(deleted),
        "deleted": deleted,
        "stop_requested_count": len(stop_requested),
        "stop_requested": stop_requested,
        "failed": failed,
        "deletion_id": deletion_id,
        "recoverable": all(
            bool(manifest.get("recoverable"))
            for manifest in archives.values()
        ),
    }


@app.get("/api/task-recycle-bin")
async def list_task_recycle_bin(
    limit: int = Query(100, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    inventory = await asyncio.to_thread(_task_recycle_inventory)
    policy = await asyncio.to_thread(app.state.storage_maintenance_manager.load_policy)
    entries = inventory["active_entries"]
    page = entries[offset:offset + limit]
    return {
        "entries": page,
        "total": inventory["active_total"],
        "offset": offset,
        "limit": limit,
        "has_more": offset + len(page) < inventory["active_total"],
        "summary": {
            "active_task_count": inventory["active_total"],
            "restored_metadata_count": inventory["restored_total"],
            "size_bytes": inventory["active_size_bytes"],
            "file_count": inventory["active_file_count"],
            "purge_blocked_count": len(inventory["unsafe_active_entries"]),
        },
        "automatic_permanent_deletion": policy["recycle"]["automatic_cleanup_enabled"],
        "retention_days": policy["recycle"]["retention_days"],
        "retention_enabled_at": policy["recycle"]["last_enabled_at"],
        "permanent_deletion_available": True,
    }


@app.post("/api/task-recycle-bin/purge/preview")
async def preview_task_recycle_bin_purge():
    inventory = await asyncio.to_thread(_task_recycle_inventory)
    if inventory["structural_errors"]:
        raise HTTPException(
            status_code=409,
            detail=(
                "回收站存在异常目录，已拒绝永久清理："
                + "; ".join(inventory["structural_errors"][:3])
            ),
        )
    if inventory["unsafe_active_entries"]:
        blocked = inventory["unsafe_active_entries"][:3]
        detail = "; ".join(
            f"{item['task_id']}: {','.join(item['purge_block_reasons'][:3])}"
            for item in blocked
        )
        raise HTTPException(
            status_code=409,
            detail=f"回收站含有无法安全验证的任务，已拒绝永久清理：{detail}",
        )

    token = uuid.uuid4().hex
    expires_at_epoch = time.time() + 120.0
    with _task_recycle_purge_preview_lock:
        now = time.time()
        for stale in [
            key for key, value in _task_recycle_purge_previews.items()
            if float(value.get("expires_at_epoch") or 0) < now
        ]:
            _task_recycle_purge_previews.pop(stale, None)
        _task_recycle_purge_previews[token] = {
            "signature": inventory["signature"],
            "task_count": inventory["active_total"],
            "size_bytes": inventory["active_size_bytes"],
            "expires_at_epoch": expires_at_epoch,
        }
    return {
        "confirmation_token": token,
        "expires_in_seconds": 120,
        "task_count": inventory["active_total"],
        "size_bytes": inventory["active_size_bytes"],
        "file_count": inventory["active_file_count"],
        "irreversible": True,
        "typed_confirmation_required": False,
    }


def _consume_task_recycle_purge_preview(
    payload: TaskRecyclePurgeRequest,
) -> dict[str, Any]:
    with _task_recycle_purge_preview_lock:
        preview = _task_recycle_purge_previews.pop(
            payload.confirmation_token, None
        )
    if not preview:
        raise HTTPException(
            status_code=409,
            detail="清空回收站确认已失效；请重新预览后确认",
        )
    if float(preview.get("expires_at_epoch") or 0) < time.time():
        raise HTTPException(
            status_code=409,
            detail="清空回收站确认已过期；请重新预览",
        )
    return preview


def _perform_task_recycle_bin_purge(
    preview: dict[str, Any],
    audit_context: dict[str, Any],
) -> dict[str, Any]:
    with _task_recycle_lock:
        inventory = _task_recycle_inventory_unlocked()
        if (
            inventory["signature"] != preview.get("signature")
            or inventory["active_total"] != preview.get("task_count")
            or inventory["active_size_bytes"] != preview.get("size_bytes")
        ):
            raise HTTPException(
                status_code=409,
                detail="回收站内容已变化；本次永久清理已拒绝，请重新预览",
            )
        if inventory["structural_errors"] or inventory["unsafe_active_entries"]:
            raise HTTPException(
                status_code=409,
                detail="回收站安全检查未通过；没有永久删除任何项目",
            )

        purge_id = (
            f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-"
            f"{uuid.uuid4().hex[:8]}"
        )
        receipt_path = (
            OUTPUT_DIR / "runtime" / f"task_recycle_purge_{purge_id}.json"
        )
        receipt: dict[str, Any] = {
            "schema_version": 1,
            "purge_id": purge_id,
            "started_at": _now_iso(),
            "state_signature": inventory["signature"],
            "requested_task_count": inventory["active_total"],
            "requested_size_bytes": inventory["active_size_bytes"],
            "targets": [
                {
                    "deletion_id": item["deletion_id"],
                    "task_id": item["task_id"],
                    "filename": item["filename"],
                    "size_bytes": item["size_bytes"],
                    "file_count": item["file_count"],
                }
                for item in inventory["active_entries"]
            ],
            "removed": [],
            "failures": [],
        }
        _write_runtime_json_atomic(receipt_path, receipt)

        recycle_root = OUTPUT_DIR / "task_recycle_bin"
        for item in inventory["active_entries"]:
            task_root = (
                recycle_root / item["deletion_id"] / item["task_id"]
            )
            try:
                if task_root.is_symlink() or not task_root.is_dir():
                    raise RuntimeError(
                        "validated recycle directory changed before removal"
                    )
                shutil.rmtree(task_root)
                receipt["removed"].append({
                    "deletion_id": item["deletion_id"],
                    "task_id": item["task_id"],
                    "filename": item["filename"],
                    "size_bytes": item["size_bytes"],
                    "file_count": item["file_count"],
                })
                try:
                    task_root.parent.rmdir()
                except OSError:
                    pass
            except Exception as exc:
                receipt["failures"].append({
                    "deletion_id": item["deletion_id"],
                    "task_id": item["task_id"],
                    "error": str(exc),
                })
            _write_runtime_json_atomic(receipt_path, receipt)

        post = _task_recycle_inventory_unlocked()
        receipt.update({
            "finished_at": _now_iso(),
            "removed_task_count": len(receipt["removed"]),
            "removed_size_bytes": sum(
                int(item["size_bytes"]) for item in receipt["removed"]
            ),
            "failure_count": len(receipt["failures"]),
            "remaining_task_count": post["active_total"],
            "remaining_size_bytes": post["active_size_bytes"],
        })
        _write_runtime_json_atomic(receipt_path, receipt)
        _append_task_deletion_audit({
            "event": "task_recycle_bin_permanent_purge",
            "at": receipt["finished_at"],
            "purge_id": purge_id,
            "receipt_path": str(receipt_path),
            "removed_task_count": receipt["removed_task_count"],
            "removed_bytes": receipt["removed_size_bytes"],
            "failure_count": receipt["failure_count"],
            "task_ids": [item["task_id"] for item in receipt["removed"]],
            "deletion_ids": sorted({
                item["deletion_id"] for item in receipt["removed"]
            }),
            **audit_context,
        })
        return {
            "purge_id": purge_id,
            "removed_task_count": receipt["removed_task_count"],
            "removed_size_bytes": receipt["removed_size_bytes"],
            "failure_count": receipt["failure_count"],
            "failures": receipt["failures"],
            "remaining_task_count": receipt["remaining_task_count"],
            "remaining_size_bytes": receipt["remaining_size_bytes"],
            "receipt_path": str(receipt_path),
            "irreversible": True,
        }


@app.delete("/api/task-recycle-bin")
async def purge_task_recycle_bin(
    payload: TaskRecyclePurgeRequest,
    request: Request,
):
    preview = _consume_task_recycle_purge_preview(payload)
    audit_context = _request_audit_context(request)
    return await asyncio.to_thread(
        _perform_task_recycle_bin_purge,
        preview,
        audit_context,
    )


def _restore_recycled_task_unlocked(
    deletion_id: str,
    task_id: str,
    request: Request,
):
    if not re.fullmatch(r"[A-Za-z0-9_-]{1,80}", deletion_id) or not re.fullmatch(
        r"[A-Za-z0-9_-]{1,80}", task_id
    ):
        raise HTTPException(status_code=400, detail="Invalid recycle identifier")
    recycle_root = OUTPUT_DIR / "task_recycle_bin" / deletion_id / task_id
    manifest_path = recycle_root / "manifest.json"
    task_path = recycle_root / "task.json"
    if not manifest_path.is_file() or not task_path.is_file():
        raise HTTPException(status_code=404, detail="Recycled task not found")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    task = json.loads(task_path.read_text(encoding="utf-8"))
    if manifest.get("restored_at"):
        raise HTTPException(status_code=409, detail="Task was already restored")
    with _tasks_lock:
        if task_id in _tasks:
            raise HTTPException(status_code=409, detail="Task ID already exists")

    moves: list[tuple[Path, Path]] = []
    for item in manifest.get("moved") or []:
        recycled = Path(str(item.get("recycle_path") or ""))
        original = Path(str(item.get("original_path") or ""))
        if _managed_output_relative(original) is None:
            raise HTTPException(status_code=409, detail="Unsafe restore destination")
        try:
            recycled.resolve(strict=False).relative_to(OUTPUT_DIR.resolve())
        except (OSError, RuntimeError, ValueError) as exc:
            raise HTTPException(status_code=409, detail="Unsafe recycled path") from exc
        if original.exists() or original.is_symlink():
            raise HTTPException(
                status_code=409,
                detail=f"Restore destination already exists: {original.name}",
            )
        if not recycled.exists() and not recycled.is_symlink():
            raise HTTPException(
                status_code=409,
                detail=f"Recycled data is incomplete: {recycled.name}",
            )
        moves.append((recycled, original))

    completed_moves: list[tuple[Path, Path]] = []
    try:
        for recycled, original in moves:
            original.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(recycled), str(original))
            completed_moves.append((original, recycled))
        _initialize_task_state(task)
        if task.get("status") == TaskState.QUEUED.value:
            _transition_task(
                task,
                TaskState.PENDING,
                reason="restored_from_task_recycle_bin",
            )
        task["message"] = "任务及其原始数据已从回收站恢复"
        task["error"] = None
        task["results"] = []
        with _tasks_lock:
            _tasks[task_id] = task
        if _save_tasks() is False:
            raise RuntimeError("Could not persist restored task registry")
    except Exception as exc:
        with _tasks_lock:
            _tasks.pop(task_id, None)
        for original, recycled in reversed(completed_moves):
            recycled.parent.mkdir(parents=True, exist_ok=True)
            if original.exists() or original.is_symlink():
                shutil.move(str(original), str(recycled))
        raise HTTPException(status_code=500, detail=f"Restore failed: {exc}") from exc

    manifest["restored_at"] = _now_iso()
    _write_runtime_json_atomic(manifest_path, manifest)
    _append_task_deletion_audit({
        "event": "task_restored_from_recycle_bin",
        "at": manifest["restored_at"],
        "deletion_id": deletion_id,
        "task_ids": [task_id],
        **_request_audit_context(request),
    })
    return {
        "restored": task_id,
        "deletion_id": deletion_id,
        "status": task.get("status"),
    }


@app.post("/api/task-recycle-bin/{deletion_id}/{task_id}/restore")
async def restore_recycled_task(
    deletion_id: str,
    task_id: str,
    request: Request,
):
    with _task_recycle_lock:
        return _restore_recycled_task_unlocked(
            deletion_id,
            task_id,
            request,
        )


# ── Results / Export ────────────────────────────────────────────────

@app.get("/api/evidence-governance/policy")
async def get_evidence_governance_policy():
    """Expose the immutable release policy used by UI and paper ablations."""
    return {
        "module": "EviGate-MS",
        "policy_id": EVIGATE_POLICY_ID,
        "purpose": (
            "applicability, provenance, independence and claim-scope governance"
        ),
        "locked_parameters": {
            "cohort_mz_ppm": EVIGATE_COHORT_MZ_PPM,
            "cohort_rt_window_min": EVIGATE_COHORT_RT_WINDOW_MIN,
            "min_independent_injections": EVIGATE_MIN_INDEPENDENT_INJECTIONS,
            "min_identity_agreement": EVIGATE_MIN_IDENTITY_AGREEMENT,
        },
        "invariants": {
            "ranking_changed": False,
            "confidence_field_overwritten": False,
            "stepped_energy_scans_are_independent": False,
            "candidate_set_has_conformal_coverage": False,
        },
    }


@app.get("/api/calibration/policy")
async def get_calibration_policy():
    """Expose P3 q-value policy without implying that a project is calibrated."""
    method_path = PROJECT_ROOT / "config" / "calibration_method.json"
    registry_path = PROJECT_ROOT / "config" / "calibration_registry.json"
    try:
        method = json.loads(method_path.read_text(encoding="utf-8"))
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Calibration policy is unreadable: {exc}",
        ) from exc
    return {
        "method": method,
        "production_calibration_id": registry.get("production_calibration_id"),
        "registered_evidence": registry.get("entries") or [],
        "policy": registry.get("policy") or {},
        "production_calibration_available": bool(
            registry.get("production_calibration_id")
        ),
        "current_legacy_high_is_qvalue_tier": False,
    }


@app.get("/api/calibration/readiness")
async def get_calibration_readiness():
    """Return the machine-audited readiness state for production q-values."""
    path = PROJECT_ROOT / "docs" / "audits" / "V9_P3_CALIBRATION_READINESS.json"
    if not path.is_file():
        raise HTTPException(
            status_code=503,
            detail="Calibration readiness audit has not been generated",
        )
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Calibration readiness audit is unreadable: {exc}",
        ) from exc


@app.post("/api/evidence-governance/cohort")
async def run_evidence_governance_cohort(
    request: EvidenceGovernanceCohortRequest,
):
    """Audit separately submitted injections with feature-collapsed voting."""
    tasks = {}
    with _tasks_lock:
        for task_id in request.task_ids:
            task = _tasks.get(task_id)
            if not task:
                raise HTTPException(
                    status_code=404, detail=f"Task not found: {task_id}"
                )
            if task.get("status") != "completed":
                raise HTTPException(
                    status_code=409,
                    detail=(
                        f"Task {task_id} is {task.get('status')}; cohort governance "
                        "requires completed tasks"
                    ),
                )
            tasks[task_id] = {
                "filename": str(task.get("filename") or ""),
                "ion_mode": str(task.get("ion_mode") or ""),
            }

    loaded = await asyncio.gather(*[
        asyncio.to_thread(_get_task_results, task_id)
        for task_id in request.task_ids
    ])
    task_results = dict(zip(request.task_ids, loaded))
    payload = await asyncio.to_thread(
        _build_replicate_consensus,
        task_results,
        task_roles=request.task_roles,
        independence_groups=request.independence_groups,
    )
    clusters = list(payload.pop("clusters", []) or [])
    # Put released/diagnostic clusters first while keeping deterministic m/z/RT
    # ordering inside each action.  Pagination changes only transport, not audit.
    action_priority = {
        "release_replicate_stable_candidate": 0,
        "release_isomer_group_candidate": 1,
        "abstain_blank_recurrent": 2,
        "abstain_identity_conflict": 3,
        "abstain_scope_or_quality_gate": 4,
        "needs_more_independent_injections": 5,
        "no_candidate": 6,
    }
    clusters.sort(key=lambda cluster: (
        action_priority.get(str(cluster.get("action") or ""), 99),
        str(cluster.get("polarity") or ""),
        float(cluster.get("retention_time_min") or math.inf),
        float(cluster.get("precursor_mz") or math.inf),
        str(cluster.get("cluster_id") or ""),
    ))
    total = len(clusters)
    payload.update({
        "tasks": tasks,
        "clusters": clusters[request.offset:request.offset + request.limit],
        "clusters_total": total,
        "offset": request.offset,
        "limit": request.limit,
        "has_more": request.offset + request.limit < total,
    })
    return payload


@app.get("/api/tasks/{task_id}/tic_overview")
async def get_task_tic_overview(task_id: str):
    """TIC/BPC 概览；旧 mzML 任务首次读取时升级源 TIC/BPC 信息。"""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    tic_path = TASK_DIR / f"{task_id}_tic_overview.json"
    if not tic_path.exists() and Path(task.get("filename", "")).suffix.lower() == ".mgf":
        tic_ov = _build_mgf_tic_overview(_get_task_results(task_id))
        if tic_ov:
            try:
                with open(tic_path, "w", encoding="utf-8") as f:
                    json.dump(tic_ov, f, separators=(",", ":"))
            except Exception as e:
                logger.warning("Failed to cache MGF TIC overview for %s: %s", task_id, e)
    if not tic_path.exists():
        raise HTTPException(status_code=404, detail="TIC overview not available (legacy task — please re-analyze)")
    try:
        return await asyncio.to_thread(
            _load_tic_overview_with_bpc_backfill,
            task_id,
            dict(task),
            tic_path,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to read TIC/BPC overview: {e}")


@app.get("/api/tasks/{task_id}/results")
def get_task_results(
    task_id: str,
    matched_only: bool = Query(False),
    match_state: Literal["all", "matched", "unmatched"] = Query("all"),
    projection: Literal["full", "summary"] = Query("full"),
    min_cosine: Optional[float] = Query(None),
    blank_view: Literal["all", "clean", "flagged"] = Query("all"),
    governance_view: Literal[
        "all", "preeligible", "blocked", "isomer_group"
    ] = Query("all"),
    limit: int = Query(500, ge=1, le=5000),
    offset: int = Query(0, ge=0),
):
    """Get paginated results with an optional lightweight list projection."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    results = task.get("results", [])
    if not results:
        result_path = TASK_DIR / f"{task_id}_results.json"
        if result_path.exists():
            try:
                with open(result_path, "r", encoding="utf-8") as f:
                    results = json.load(f)
            except Exception:
                pass

    results = _prepare_results_for_display(results)
    try:
        claim_projection = _attach_claim_decisions(
            results,
            _claim_contract_path(TASK_DIR, task_id),
        )
    except (OSError, ValueError) as exc:
        claim_projection = {
            "status": "failed",
            "error": str(exc),
            "legacy_artifact_mutated": False,
        }

    if (
        task.get("status") == "completed"
        and task.get("presentation_calibration_version")
        != PRESENTATION_CALIBRATION_VERSION
    ):
        task = _refresh_task_presentation_from_prepared(task_id, results)

    from web.governance_blocked_projection import project_governance_blocked_candidates
    governance_blocked_by_index = {
        row["source_result_index"]: row
        for row in project_governance_blocked_candidates(results)
    }
    indexed_results = list(enumerate(results))

    # Filter while preserving the immutable result-artifact index.  The index
    # is used by the detail endpoint so browser sorting cannot open the wrong
    # spectrum or structure.
    if blank_view == "clean":
        indexed_results = [
            item for item in indexed_results
            if not bool(item[1].get("blank_clean_view_excluded"))
        ]
    elif blank_view == "flagged":
        indexed_results = [
            item for item in indexed_results
            if bool(item[1].get("blank_flagged"))
        ]

    view_counts = {
        "total": len(indexed_results),
        "matched": sum(bool(row.get("match")) for _, row in indexed_results),
        "unmatched": sum(not bool(row.get("match")) for _, row in indexed_results),
        "blank_flagged": sum(
            bool(row.get("blank_flagged")) for _, row in indexed_results
        ),
        "blank_excluded": sum(
            bool(row.get("blank_excluded")) for _, row in indexed_results
        ),
        "blank_clean_view_excluded": sum(
            bool(row.get("blank_clean_view_excluded"))
            for _, row in indexed_results
        ),
    }

    effective_match_state = "matched" if matched_only else match_state
    if effective_match_state == "matched":
        indexed_results = [
            item for item in indexed_results if item[1].get("match")
        ]
    elif effective_match_state == "unmatched":
        indexed_results = [
            item for item in indexed_results if not item[1].get("match")
        ]
    if min_cosine is not None:
        indexed_results = [
            item for item in indexed_results
            if item[1].get("match")
            and item[1]["match"].get("cosine", 0) >= min_cosine
        ]
    if governance_view == "preeligible":
        indexed_results = [
            (result_index, row) for result_index, row in indexed_results
            if bool(((row.get("match") or {}).get(
                "evidence_governance"
            ) or {}).get("release_eligible_before_replication"))
        ]
    elif governance_view == "blocked":
        indexed_results = [
            (result_index, row) for result_index, row in indexed_results
            if isinstance(row.get("match"), dict)
            and not bool((row["match"].get("evidence_governance") or {}).get(
                "release_eligible_before_replication"
            ))
        ]
    elif governance_view == "isomer_group":
        indexed_results = [
            (result_index, row) for result_index, row in indexed_results
            if str(((row.get("match") or {}).get(
                "evidence_governance"
            ) or {}).get("claim_status") or "") == "isomer_group_only"
        ]

    total = len(indexed_results)
    page_items = indexed_results[offset:offset + limit]
    page_result_indices = {result_index for result_index, _ in page_items}
    governance_blocked_candidates = [
        row for result_index, row in governance_blocked_by_index.items()
        if result_index in page_result_indices
    ]
    governance_blocked_count = sum(
        result_index in governance_blocked_by_index
        for result_index, _ in indexed_results
    )
    if projection == "summary":
        paged = [
            _result_summary_projection(row, result_index)
            for result_index, row in page_items
        ]
    else:
        paged = [row for _, row in page_items]

    return {
        "results": paged,
        "total": total,
        "offset": offset,
        "limit": limit,
        "blank_view": blank_view,
        "match_state": effective_match_state,
        "governance_view": governance_view,
        "projection": projection,
        "view_counts": view_counts,
        "has_more": (offset + limit) < total,
        "claim_projection": claim_projection,
        "governance_blocked_candidates": governance_blocked_candidates,
        "governance_blocked_count": governance_blocked_count,
    }


@app.get("/api/tasks/{task_id}/results/item/{result_index}")
def get_task_result_item(task_id: str, result_index: int):
    """Return one complete result row selected by immutable artifact index."""

    source_results = _get_task_results(task_id)
    # Preserve the immutable source-artifact index through presentation
    # deduplication.  Sidecar evidence (for example EvoGraph-F) is keyed by
    # this index, not by the post-dedup table position.
    results = [
        {**row, "_result_index": result_index}
        for result_index, row in enumerate(source_results)
    ]
    if result_index < 0 or result_index >= len(results):
        raise HTTPException(status_code=404, detail="Result index out of range")
    row = deepcopy(results[result_index])
    match = row.get("match")
    if isinstance(match, dict):
        _normalize_match_compound_names(match)
        _annotate_match_compound_roles(match)
        match["auxiliary_interpretation"] = _compound_auxiliary_interpretation(
            match, row
        )
    graph_artifact = _evograph_sidecar_path(TASK_DIR, task_id)
    if graph_artifact.is_file():
        try:
            graph_payload = _read_evograph_sidecar(
                graph_artifact,
                offset=0,
                limit=1,
                result_index=result_index,
            )
            graph_records = graph_payload.get("records") or []
            if graph_records:
                row["evograph_fragmentation"] = graph_records[0]
                if isinstance(match, dict):
                    match["auxiliary_interpretation"]["fragmentation"][
                        "formula_graph"
                    ] = graph_records[0]
        except (OSError, ValueError, TypeError):
            logger.warning(
                "Could not attach EvoGraph record for %s/%s",
                task_id,
                result_index,
            )
    return {
        "task_id": task_id,
        "result_index": result_index,
        "result": row,
        "projection": "full_single_item",
        "presentation_projection_current": True,
        "saved_legacy_result_unchanged": True,
    }


@app.post("/api/tasks/{task_id}/evograph-fragmentation")
async def launch_task_evograph_fragmentation(
    task_id: str,
    refresh: bool = Query(
        False,
        description="Rebuild the shadow artifact even when the source hash is unchanged",
    ),
):
    """Queue an isolated EvoGraph-F shadow materialization."""

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != TaskState.COMPLETED.value:
        raise HTTPException(
            status_code=409,
            detail="EvoGraph-F starts only after immutable results are committed",
        )
    scheduled = _maybe_launch_evograph_fragmentation(task_id, force=refresh)
    with _tasks_lock:
        summary = dict(
            (_tasks.get(task_id) or {}).get("evograph_fragmentation_summary") or {}
        )
    return {
        "task_id": task_id,
        "scheduled": scheduled,
        "summary": summary,
        "module_state": "shadow",
        "used_for_scoring": False,
        "production_ranking_changed": False,
    }


@app.get("/api/tasks/{task_id}/evograph-fragmentation")
async def get_task_evograph_fragmentation(
    task_id: str,
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    result_index: int | None = Query(None, ge=0),
    materialize: bool = Query(
        True,
        description="Build the additive shadow sidecar when it is absent",
    ),
):
    """Read paginated EvoGraph-F hypotheses without touching legacy results."""

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != TaskState.COMPLETED.value:
        raise HTTPException(
            status_code=409,
            detail="EvoGraph-F is available only for completed tasks",
        )

    artifact = _evograph_sidecar_path(TASK_DIR, task_id)
    if not artifact.is_file():
        summary = dict(task.get("evograph_fragmentation_summary") or {})
        if str(summary.get("status") or "") in {"queued", "running"}:
            return {
                "task_id": task_id,
                "metadata": summary,
                "records": [],
                "offset": offset,
                "limit": limit,
                "has_more": False,
                "materialization_in_progress": True,
                "module_state": "shadow",
                "used_for_scoring": False,
            }
        if not materialize:
            raise HTTPException(
                status_code=404,
                detail="EvoGraph-F shadow artifact has not been materialized",
            )
        legacy_path = TASK_DIR / f"{task_id}_results.json"
        if not legacy_path.is_file():
            raise HTTPException(
                status_code=404,
                detail="Saved legacy task results are unavailable",
            )
        try:
            metadata = await asyncio.to_thread(
                _materialize_evograph_sidecar,
                task_id,
                legacy_path,
                artifact,
                policy_path=PROJECT_ROOT / "config" / "evograph_policy.json",
                runtime_snapshot_id=_current_runtime_snapshot_id(),
            )
        except Exception as exc:
            logger.exception("Could not materialize EvoGraph-F for %s", task_id)
            raise HTTPException(
                status_code=500,
                detail=f"EvoGraph-F materialization failed: {exc}",
            ) from exc
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current is not None:
                current["evograph_fragmentation_summary"] = _evograph_task_summary(
                    metadata,
                    status="completed",
                )
        _save_tasks()

    try:
        response = await asyncio.to_thread(
            _read_evograph_sidecar,
            artifact,
            offset=offset,
            limit=limit,
            result_index=result_index,
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"EvoGraph-F shadow artifact is unreadable: {exc}",
        ) from exc
    if result_index is not None and not response.get("records"):
        raise HTTPException(status_code=404, detail="Result index is not in the EvoGraph artifact")
    response.update({
        "task_id": task_id,
        "module_state": "shadow",
        "legacy_results_unchanged": True,
        "used_for_scoring": False,
        "production_ranking_changed": False,
        "legacy_confidence_changed": False,
        "high_promotions": 0,
        "heuristic_marginals_are_probabilities": False,
        "log_bayes_factor_available": False,
    })
    return response


@app.get("/api/tasks/{task_id}/identification-records")
async def get_task_identification_records(
    task_id: str,
    limit: int = Query(500, ge=1, le=5000),
    offset: int = Query(0, ge=0),
    materialize: bool = Query(
        True,
        description=(
            "Build the additive P1 sidecar from the immutable legacy result if absent"
        ),
    ),
):
    """Return the versioned IdentificationRecord DTO without raw peak arrays."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Identification records are available only for completed tasks",
        )

    artifact = _identification_contract_path(TASK_DIR, task_id)
    if not artifact.is_file():
        if not materialize:
            raise HTTPException(
                status_code=404,
                detail="Identification contract has not been materialized",
            )
        legacy_path = TASK_DIR / f"{task_id}_results.json"
        if not legacy_path.is_file():
            raise HTTPException(
                status_code=404,
                detail="Saved legacy task results are unavailable",
            )
        try:
            metadata = await asyncio.to_thread(
                _write_identification_contract,
                task_id,
                _iter_legacy_result_artifact(legacy_path),
                artifact,
                runtime_snapshot_id=_current_runtime_snapshot_id(),
            )
        except Exception as exc:
            logger.exception(
                "Could not materialize identification contract for %s", task_id
            )
            raise HTTPException(
                status_code=500,
                detail=f"Identification contract migration failed: {exc}",
            ) from exc
        metadata["status"] = "available"
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current is not None:
                current["identification_contract_summary"] = metadata
        _save_tasks()

    try:
        response = await asyncio.to_thread(
            _read_identification_contract,
            artifact,
            offset=offset,
            limit=limit,
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Identification contract is unreadable: {exc}",
        ) from exc
    response["task_id"] = task_id
    response["legacy_results_unchanged"] = True
    response["ordinary_ms2_identifications_only"] = True
    response["expected_compounds_counted"] = False
    return response


@app.get("/api/tasks/{task_id}/compounds")
async def get_task_compound_conclusions(
    task_id: str,
    limit: int = Query(2000, ge=1, le=5000),
    offset: int = Query(0, ge=0),
    governance_limit: int = Query(5000, ge=1, le=5000),
    governance_offset: int = Query(0, ge=0),
    report_scope: Literal["all", "plant_candidates", "external_review"] = Query("all"),
    conclusion: Optional[
        Literal[
            "confirmed",
            "high",
            "isomer",
            "candidate",
            "background",
            "unresolved_event",
            "insufficient",
            "unannotated",
        ]
    ] = Query(None),
    materialize: bool = Query(
        True,
        description=(
            "Build the deduplicated user conclusion from immutable scan results "
            "when absent or stale"
        ),
    ),
):
    """Return plain-language, deduplicated compound conclusions for one sample."""

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Compound conclusions are available only for completed tasks",
        )

    source_artifact = TASK_DIR / f"{task_id}_results.json"
    if not source_artifact.is_file():
        raise HTTPException(
            status_code=404,
            detail="Saved scan-level evidence is unavailable",
        )
    primary_artifact = _flagship_artifact_path(TASK_DIR, task_id)
    if not _flagship_artifact_is_current(primary_artifact):
        if not materialize:
            raise HTTPException(
                status_code=404,
                detail="Current V15 primary identity artifact is absent or stale",
            )
        try:
            flagship_metadata = await asyncio.to_thread(
                _materialize_current_v15_flagship,
                task_id,
                task,
                source_artifact,
                primary_artifact,
            )
        except Exception as exc:
            logger.exception(
                "Could not materialize current V15 flagship for %s", task_id
            )
            raise HTTPException(
                status_code=500,
                detail=f"V15 flagship materialization failed: {exc}",
            ) from exc
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current is not None:
                current["flagship_summary"] = flagship_metadata
        _save_tasks()
    try:
        feature_inventory, feature_inventory_source = _task_feature_inventory(task)
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Independent MS1 feature inventory is unreadable: {exc}",
        ) from exc
    artifact = _compound_conclusion_artifact_path(TASK_DIR, task_id)
    if not _compound_conclusion_artifact_is_current(
        artifact, source_artifact, primary_artifact, feature_inventory_source
    ):
        if not materialize:
            raise HTTPException(
                status_code=404,
                detail="Compound conclusion artifact is absent or stale",
            )
        try:
            prepared_rows = await asyncio.to_thread(
                _get_task_results, task_id
            )
            metadata = await asyncio.to_thread(
                _write_compound_conclusion_artifact,
                task_id,
                prepared_rows,
                artifact,
                claim_decisions_path=_claim_contract_path(TASK_DIR, task_id),
                source_artifact=source_artifact,
                primary_identity_artifact=primary_artifact,
                feature_inventory=feature_inventory,
                feature_inventory_source=feature_inventory_source,
            )
        except Exception as exc:
            logger.exception(
                "Could not materialize compound conclusions for %s", task_id
            )
            raise HTTPException(
                status_code=500,
                detail=f"Compound conclusion materialization failed: {exc}",
            ) from exc
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current is not None:
                current["compound_conclusion_summary"] = metadata
        _save_tasks()

    try:
        response = await asyncio.to_thread(
            _read_compound_conclusion_artifact,
            artifact,
            offset=offset,
            limit=limit,
            conclusion=conclusion,
            report_scope=report_scope,
        )
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Compound conclusion artifact is unreadable: {exc}",
        ) from exc
    response_metadata = dict(response.get("metadata") or {})
    with _tasks_lock:
        current = _tasks.get(task_id)
        current_summary = (
            current.get("compound_conclusion_summary", {}) if current else {}
        )
        summary_changed = bool(current is not None and (
            current_summary.get("algorithm_id")
            != response_metadata.get("algorithm_id")
            or current_summary.get("source_semantic_sha256")
            != response_metadata.get("source_semantic_sha256")
        ))
        if current is not None and summary_changed:
            current["compound_conclusion_summary"] = {
                **response_metadata,
                "status": "available",
                "artifact": str(artifact.resolve()),
                "size_bytes": artifact.stat().st_size,
            }
        count_changed = bool(
            current is not None
            and _apply_final_compound_count_summary(current, response_metadata)
        )
    if summary_changed or count_changed:
        _save_tasks()
    response.update({
        "task_id": task_id,
        "result_unit": (
            "deduplicated_compound_unresolved_isomer_group_or_"
            "unresolved_precursor_event"
        ),
        "confirmed_is_level1_only": True,
        "auxiliary_rank_influence": False,
        "legacy_scan_results_unchanged": True,
    })
    from web.governance_blocked_projection import project_governance_blocked_candidates
    saved_results = await asyncio.to_thread(_get_task_results, task_id)
    blocked_candidates = project_governance_blocked_candidates(saved_results)
    response["governance_blocked_candidates"] = blocked_candidates[
        governance_offset:governance_offset + governance_limit
    ]
    response["governance_blocked_count"] = len(blocked_candidates)
    response["governance_blocked_offset"] = governance_offset
    response["governance_blocked_limit"] = governance_limit
    response["governance_blocked_has_more"] = (
        governance_offset + governance_limit < len(blocked_candidates)
    )
    return response


@app.get("/api/tasks/{task_id}/replicate-consistency")
async def get_task_replicate_consistency(task_id: str):
    """Report technical-repeat robustness without changing single-run claims."""

    with _tasks_lock:
        task_snapshot = {
            key: dict(value) for key, value in _tasks.items()
            if isinstance(value, dict)
        }
    current = task_snapshot.get(task_id)
    if not current:
        raise HTTPException(status_code=404, detail="Task not found")
    if current.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Replicate consistency is available only for completed tasks",
        )

    parsed_current = _parse_replicate_filename(current.get("filename"))
    if parsed_current:
        signature = _replicate_parameter_signature(current)
        sibling_ids = []
        for sibling_id, sibling in task_snapshot.items():
            parsed = _parse_replicate_filename(sibling.get("filename"))
            if not parsed or sibling.get("status") != "completed":
                continue
            if (
                parsed["base_key"] == parsed_current["base_key"]
                and parsed["extension"] == parsed_current["extension"]
                and _replicate_parameter_signature(sibling) == signature
            ):
                sibling_ids.append(sibling_id)
        # The repeatability view must consume the current presentation schema.
        # Rebuild only stale compact sidecars; immutable search results remain
        # untouched and the full spectral library is not searched again.
        for sibling_id in sibling_ids:
            source_artifact = TASK_DIR / f"{sibling_id}_results.json"
            artifact = _compound_conclusion_artifact_path(TASK_DIR, sibling_id)
            primary_artifact = _flagship_artifact_path(TASK_DIR, sibling_id)
            sibling_task = task_snapshot.get(sibling_id) or {}
            feature_inventory, feature_inventory_source = (
                _task_feature_inventory(sibling_task)
            )
            if not source_artifact.is_file() or _compound_conclusion_artifact_is_current(
                artifact,
                source_artifact,
                primary_artifact,
                feature_inventory_source,
            ):
                continue
            if not _flagship_artifact_is_current(primary_artifact):
                await asyncio.to_thread(
                    _materialize_current_v15_flagship,
                    sibling_id,
                    sibling_task,
                    source_artifact,
                    primary_artifact,
                )
            prepared_rows = await asyncio.to_thread(
                _get_task_results, sibling_id
            )
            metadata = await asyncio.to_thread(
                _write_compound_conclusion_artifact,
                sibling_id,
                prepared_rows,
                artifact,
                claim_decisions_path=_claim_contract_path(TASK_DIR, sibling_id),
                source_artifact=source_artifact,
                primary_identity_artifact=primary_artifact,
                feature_inventory=feature_inventory,
                feature_inventory_source=feature_inventory_source,
            )
            with _tasks_lock:
                sibling_task = _tasks.get(sibling_id)
                if sibling_task is not None:
                    sibling_task["compound_conclusion_summary"] = metadata
        if sibling_ids:
            _save_tasks()

    with _tasks_lock:
        refreshed_snapshot = {
            key: dict(value) for key, value in _tasks.items()
            if isinstance(value, dict)
        }
    return await asyncio.to_thread(
        _build_replicate_consistency_view,
        task_id,
        refreshed_snapshot,
        TASK_DIR,
    )


@app.get("/api/tasks/{task_id}/flagship")
async def get_task_flagship(
    task_id: str,
    limit: int = Query(500, ge=1, le=5000),
    offset: int = Query(0, ge=0),
    decision_status: Optional[str] = Query(None),
    entity_id: Optional[str] = Query(None),
    materialize: bool = Query(
        True,
        description=(
            "Build the primary V15 conditioned multi-library entity decision from immutable scan results if absent"
        ),
    ),
):
    """Return the primary V15 entity-level identification decisions."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="V15 flagship decisions are available only for completed tasks",
        )

    artifact = _flagship_artifact_path(TASK_DIR, task_id)
    if not _flagship_artifact_is_current(artifact):
        if not materialize:
            raise HTTPException(
                status_code=404,
                detail="Current V15 flagship artifact is absent or stale",
            )
        legacy_path = TASK_DIR / f"{task_id}_results.json"
        if not legacy_path.is_file():
            raise HTTPException(
                status_code=404,
                detail="Saved scan-level candidate evidence is unavailable",
            )
        try:
            metadata = await asyncio.to_thread(
                _materialize_current_v15_flagship,
                task_id,
                task,
                legacy_path,
                artifact,
            )
        except Exception as exc:
            logger.exception("Could not materialize V15 flagship result for %s", task_id)
            raise HTTPException(
                status_code=500,
                detail=f"V15 flagship materialization failed: {exc}",
            ) from exc
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current is not None:
                current["flagship_summary"] = metadata
        _save_tasks()

    try:
        response = await asyncio.to_thread(
            _read_flagship_artifact,
            artifact,
            offset=offset,
            limit=limit,
            decision_status=decision_status,
            entity_id=entity_id,
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"V15 flagship artifact is unreadable: {exc}",
        ) from exc
    if entity_id and not response.get("entities"):
        raise HTTPException(status_code=404, detail="Entity is not in the V15 flagship artifact")
    response.update(
        {
            "task_id": task_id,
            "primary_identity_selector": True,
            "algorithm_generation": "V15",
            "conditioned_multilibrary_rank_influence": True,
            "candidate_recall_target_total_connectivities": 40,
            "candidate_recall_support_affects_primary_ranking": False,
            "v14_audit_baseline_embedded": True,
            "legacy_scan_results_unchanged": True,
            "expert_rank_influence": False,
            "classification_rank_influence": False,
            "candidate_pathway_rank_influence": False,
            "score_is_probability": False,
        }
    )
    return response


@app.get("/api/tasks/{task_id}/flagship/multilibrary")
async def get_task_flagship_multilibrary(
    task_id: str,
    limit: int = Query(5000, ge=1, le=100000),
    offset: int = Query(0, ge=0),
    download: bool = Query(False),
):
    """Return the V15 entity/candidate/CE/library long table as JSON."""

    artifact = _flagship_artifact_path(TASK_DIR, task_id)
    if not _flagship_artifact_is_current(artifact):
        await get_task_flagship(
            task_id,
            limit=1,
            offset=0,
            decision_status=None,
            entity_id=None,
            materialize=True,
        )
    try:
        payload = await asyncio.to_thread(read_multilibrary_evidence, artifact)
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"V15 multi-library evidence is unreadable: {exc}",
        ) from exc
    rows = list(payload.get("rows") or [])
    response_payload = {
        "schema_version": payload.get("schema_version"),
        "task_id": task_id,
        "columns": MULTILIBRARY_COLUMNS,
        "rows": rows[offset : offset + limit],
        "rows_total": len(rows),
        "offset": offset,
        "limit": limit,
        "has_more": offset + limit < len(rows),
        "score_is_probability": False,
        "unconditioned_consensus_affects_ranking": False,
    }
    if not download:
        return response_payload
    return JSONResponse(
        response_payload,
        headers={
            "Content-Disposition": (
                f'attachment; filename="{task_id}_v15_multilibrary_evidence.json"'
            ),
            "Cache-Control": "no-store, max-age=0",
        },
    )


@app.get("/api/tasks/{task_id}/evirank-mv")
async def get_task_evirank_mv(
    task_id: str,
    limit: int = Query(500, ge=1, le=5000),
    offset: int = Query(0, ge=0),
    entity_id: Optional[str] = Query(None),
    detection_method_id: Optional[str] = Query(
        None,
        min_length=1,
        max_length=256,
        description=(
            "Stable acquisition/detection-method identifier. It is retained as "
            "method context and never used as ground-truth identity."
        ),
    ),
    materialize: bool = Query(
        True,
        description=(
            "Build the V12 shadow sidecar from the immutable V11 artifact if absent"
        ),
    ),
):
    """Return V12-Flagship / EviRank-MV shadow decisions without changing V11."""

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="V12 EviRank-MV shadow decisions are available only for completed tasks",
        )

    artifact = _evirank_artifact_path(TASK_DIR, task_id)
    if not artifact.is_file():
        if not materialize:
            raise HTTPException(
                status_code=404,
                detail="V12 EviRank-MV shadow artifact has not been materialized",
            )
        v11_artifact = _v11_flagship_artifact_path(TASK_DIR, task_id)
        if not v11_artifact.is_file():
            legacy_path = TASK_DIR / f"{task_id}_results.json"
            if not legacy_path.is_file():
                raise HTTPException(
                    status_code=404,
                    detail="Saved scan-level candidate evidence is unavailable",
                )
            try:
                flagship_metadata = await asyncio.to_thread(
                    _write_v11_flagship_artifact,
                    task_id,
                    _iter_legacy_result_artifact(legacy_path),
                    v11_artifact,
                    runtime_snapshot_id=_current_runtime_snapshot_id(),
                )
            except Exception as exc:
                logger.exception(
                    "Could not materialize V11 base artifact for V12 task %s",
                    task_id,
                )
                raise HTTPException(
                    status_code=500,
                    detail=f"V11 base materialization for V12 failed: {exc}",
                ) from exc
            with _tasks_lock:
                current = _tasks.get(task_id)
                if current is not None:
                    current["flagship_summary"] = flagship_metadata
            _save_tasks()
        method_id = str(detection_method_id or "").strip() or "unspecified"
        try:
            metadata = await asyncio.to_thread(
                _write_evirank_artifact,
                task_id,
                v11_artifact,
                artifact,
                detection_method_id=method_id,
                runtime_snapshot_id=_current_runtime_snapshot_id(),
            )
        except Exception as exc:
            logger.exception("Could not materialize V12 EviRank-MV result for %s", task_id)
            raise HTTPException(
                status_code=500,
                detail=f"V12 EviRank-MV materialization failed: {exc}",
            ) from exc
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current is not None:
                current["evirank_mv_summary"] = metadata
        _save_tasks()

    try:
        response = await asyncio.to_thread(
            _read_evirank_artifact,
            artifact,
            offset=offset,
            limit=limit,
            entity_id=entity_id,
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"V12 EviRank-MV artifact is unreadable: {exc}",
        ) from exc
    bound_method_id = str((response.get("metadata") or {}).get("detection_method_id") or "")
    requested_method_id = str(detection_method_id or "").strip()
    if requested_method_id and requested_method_id != bound_method_id:
        raise HTTPException(
            status_code=409,
            detail=(
                "The existing V12 artifact is bound to detection_method_id "
                f"{bound_method_id!r}; delete/rematerialize it explicitly instead of "
                "silently rebinding method context"
            ),
        )
    if entity_id and not response.get("entities"):
        raise HTTPException(
            status_code=404,
            detail="Entity is not in the V12 EviRank-MV artifact",
        )
    response.update(
        {
            "task_id": task_id,
            "shadow_challenger": True,
            "production_identity_selector": False,
            "production_identity_endpoint": f"/api/tasks/{task_id}/flagship",
            "v11_production_result_unchanged": True,
            "legacy_scan_results_unchanged": True,
            "candidate_recall_expanded": False,
            "score_is_probability": False,
        }
    )
    return response


@app.get("/api/tasks/{task_id}/flagship/structure")
async def get_task_flagship_structure(
    task_id: str,
    entity_id: str = Query(..., min_length=1),
    format: str = Query("png", enum=["png", "svg"]),
):
    """Render the structure belonging to the V15 entity decision."""

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    artifact = _flagship_artifact_path(TASK_DIR, task_id)
    if not artifact.is_file():
        raise HTTPException(status_code=404, detail="V15 flagship artifact is unavailable")
    try:
        page = await asyncio.to_thread(
            _read_flagship_artifact,
            artifact,
            offset=0,
            limit=1,
            entity_id=entity_id,
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"V15 flagship artifact is unreadable: {exc}",
        ) from exc
    if not page.get("entities"):
        raise HTTPException(status_code=404, detail="Entity is not in the V15 flagship artifact")
    entity = page["entities"][0]
    selected = entity.get("winner") or entity.get("leading_candidate_when_abstained") or {}
    smiles = str((selected.get("candidate") or {}).get("smiles") or "").strip()
    if not smiles:
        raise HTTPException(status_code=404, detail="No structure is available for this entity decision")

    from fastapi.responses import Response
    from web.structure import smiles_to_image_base64, smiles_to_svg

    if format == "svg":
        svg = smiles_to_svg(smiles)
        if not svg:
            raise HTTPException(status_code=500, detail="Failed to render flagship structure")
        return Response(content=svg, media_type="image/svg+xml")
    image_base64 = smiles_to_image_base64(smiles)
    if image_base64:
        import base64

        return Response(content=base64.b64decode(image_base64), media_type="image/png")
    svg = smiles_to_svg(smiles)
    if not svg:
        raise HTTPException(status_code=500, detail="Failed to render flagship structure")
    return Response(content=svg, media_type="image/svg+xml")


@app.get("/api/tasks/{task_id}/evitwin-ms")
async def get_task_evitwin_ms(
    task_id: str,
    materialize: bool = Query(
        True,
        description=(
            "Build the authentic-reference display-ranking sidecar when absent"
        ),
    ),
):
    """Return EviTwin-MS authentic-standard candidate ranking.

    The endpoint may change the candidate order shown in this independent
    authentic-reference lane.  It is bound to the primary V15 entities, never
    overwrites them, and never performs a formal Level-1 release without a
    separate same-method RT contract.
    """

    with _tasks_lock:
        task = dict(_tasks.get(task_id) or {})
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="EviTwin-MS is available only for completed tasks",
        )
    artifact = _evitwin_artifact_path(task_id)
    result_path = TASK_DIR / f"{task_id}_results.json"
    primary_path = _flagship_artifact_path(TASK_DIR, task_id)
    try:
        registry_snapshot = json.loads(
            EVITWIN_REGISTRY_FILE.read_text(encoding="utf-8")
        )
        library_path = _evitwin_library_path(registry_snapshot)
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise HTTPException(
            status_code=503,
            detail=f"EviTwin-MS authentic-reference registry is unavailable: {exc}",
        ) from exc
    source_paths = [
        result_path,
        primary_path,
        EVITWIN_REGISTRY_FILE,
        library_path,
    ]
    analysis_source_path = Path(str(task.get("analysis_file_path") or ""))
    if analysis_source_path.is_file():
        source_paths.append(analysis_source_path)
    invalid_existing_artifact = False
    if artifact.is_file():
        try:
            await asyncio.to_thread(_read_evitwin_result, artifact)
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            invalid_existing_artifact = True
    stale = bool(
        invalid_existing_artifact
        or (
            artifact.is_file()
            and any(
                path.is_file() and path.stat().st_mtime_ns > artifact.stat().st_mtime_ns
                for path in source_paths
            )
        )
    )
    if (not artifact.is_file() or stale) and not materialize:
        raise HTTPException(
            status_code=404,
            detail="EviTwin-MS artifact is absent or stale",
        )
    if not artifact.is_file() or stale:
        if not result_path.is_file():
            raise HTTPException(
                status_code=404,
                detail="Saved result/entity evidence required by EviTwin-MS is unavailable",
            )
        polarity = _task_evitwin_polarity(task)
        if polarity not in {"POS", "NEG"}:
            raise HTTPException(
                status_code=409,
                detail="EviTwin-MS requires one resolved task polarity",
            )
        try:
            if not _flagship_artifact_is_current(primary_path):
                await asyncio.to_thread(
                    _materialize_current_v15_flagship,
                    task_id,
                    task,
                    result_path,
                    primary_path,
                )
            library, registry = await asyncio.to_thread(_active_evitwin_library)
            rows = await asyncio.to_thread(
                lambda: json.loads(result_path.read_text(encoding="utf-8"))
            )
            primary = await asyncio.to_thread(
                lambda: json.loads(primary_path.read_text(encoding="utf-8"))
            )
            instrument_audit = _resolve_evitwin_query_instrument_family(
                task,
                reference_instrument_family=str(
                    registry.get("instrument_family") or "agilent-qtof-6546"
                ),
            )
            method_audit = _resolve_evitwin_query_method_id(task)
            payload = await asyncio.to_thread(
                _build_evitwin_result,
                primary,
                rows,
                library,
                task_id=task_id,
                polarity=polarity,
                query_instrument_family=str(
                    instrument_audit["scoring_instrument_family"]
                ),
                query_method_id=str(
                    method_audit.get("query_method_id") or ""
                ),
                evirank_payload=None,
            )
            payload["metadata"]["instrument_applicability"] = instrument_audit
            payload["metadata"]["query_method_applicability"] = method_audit
            await asyncio.to_thread(_write_evitwin_result, artifact, payload)
        except Exception as exc:
            logger.exception("Could not materialize EviTwin-MS result for %s", task_id)
            raise HTTPException(
                status_code=500,
                detail=f"EviTwin-MS materialization failed: {exc}",
            ) from exc
    try:
        response = await asyncio.to_thread(_read_evitwin_result, artifact)
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"EviTwin-MS artifact is unreadable: {exc}",
        ) from exc
    response.update(
        {
            "task_id": task_id,
            "result_unit": "authentic_reference_ranked_compound_or_isomer_group",
            "candidate_display_order_active": True,
            "primary_v15_artifact_unchanged": True,
            "legacy_v11_v12_artifacts_unchanged": True,
            "formal_level1_release_increment": 0,
        }
    )
    return response


@app.get("/api/tasks/{task_id}/entities")
async def get_task_entities(
    task_id: str,
    node_limit: int = Query(500, ge=1, le=5000),
    node_offset: int = Query(0, ge=0),
    consensus_limit: int = Query(500, ge=1, le=5000),
    consensus_offset: int = Query(0, ge=0),
    materialize: bool = Query(
        True,
        description=(
            "Build the additive P2 entity bundle from the immutable legacy result "
            "if absent"
        ),
    ),
):
    """Return the five-level acquisition entity graph and shadow consensus."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Entity graphs are available only for completed tasks",
        )

    artifact = _entity_bundle_path(TASK_DIR, task_id)
    if not artifact.is_file():
        if not materialize:
            raise HTTPException(
                status_code=404,
                detail="Entity bundle has not been materialized",
            )
        legacy_path = TASK_DIR / f"{task_id}_results.json"
        if not legacy_path.is_file():
            raise HTTPException(
                status_code=404,
                detail="Saved legacy task results are unavailable",
            )
        try:
            metadata = await asyncio.to_thread(
                _write_entity_bundle,
                task_id,
                _iter_legacy_result_artifact(legacy_path),
                artifact,
            )
        except Exception as exc:
            logger.exception("Could not materialize entity bundle for %s", task_id)
            raise HTTPException(
                status_code=500,
                detail=f"Entity bundle migration failed: {exc}",
            ) from exc
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current is not None:
                current["entity_bundle_summary"] = metadata
        _save_tasks()

    try:
        response = await asyncio.to_thread(
            _read_entity_bundle,
            artifact,
            node_offset=node_offset,
            node_limit=node_limit,
            consensus_offset=consensus_offset,
            consensus_limit=consensus_limit,
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Entity bundle is unreadable: {exc}",
        ) from exc
    response["legacy_results_unchanged"] = True
    response["feature_consensus_runtime_influence"] = False
    response["expected_compounds_counted"] = False
    return response


@app.get("/api/tasks/{task_id}/claims")
async def get_task_claims(
    task_id: str,
    limit: int = Query(500, ge=1, le=5000),
    offset: int = Query(0, ge=0),
    materialize: bool = Query(
        True,
        description=(
            "Compile the additive EviGate 2 decision artifact from IdentificationRecords "
            "if absent"
        ),
    ),
):
    """Return final claim decisions; legacy High remains a separate field."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Final claim decisions are available only for completed tasks",
        )
    artifact = _claim_contract_path(TASK_DIR, task_id)
    if not artifact.is_file():
        if not materialize:
            raise HTTPException(
                status_code=404,
                detail="Claim decision contract has not been materialized",
            )
        identification_path = _identification_contract_path(TASK_DIR, task_id)
        if not identification_path.is_file():
            raise HTTPException(
                status_code=409,
                detail="IdentificationRecord contract must be materialized first",
            )
        try:
            metadata = await asyncio.to_thread(
                _write_claim_contract,
                task_id,
                _iter_identification_records(identification_path),
                artifact,
                context=_default_task_claim_context(task_id),
            )
        except Exception as exc:
            logger.exception("Could not compile EviGate 2 claims for %s", task_id)
            raise HTTPException(
                status_code=500,
                detail=f"Claim compilation failed: {exc}",
            ) from exc
        metadata["status"] = "available"
        with _tasks_lock:
            current = _tasks.get(task_id)
            if current is not None:
                current["claim_contract_summary"] = metadata
        _save_tasks()
    try:
        response = await asyncio.to_thread(
            _read_claim_contract,
            artifact,
            offset=offset,
            limit=limit,
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Claim contract is unreadable: {exc}",
        ) from exc
    response["task_id"] = task_id
    response["legacy_results_unchanged"] = True
    response["legacy_high_is_project_qvalue_tier"] = False
    response["qvalues_are_probabilities"] = False
    return response


@app.get("/api/tasks/{task_id}/resources")
async def get_task_resources(task_id: str):
    """List the canonical read-only resources for one task without leaking paths."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    resources = [
        (
            "compound_conclusions",
            _compound_conclusion_artifact_path(TASK_DIR, task_id),
            f"/api/tasks/{task_id}/compounds",
            "deduplicated_plain_language_sample_conclusions",
        ),
        (
            "v15_flagship_entities",
            _flagship_artifact_path(TASK_DIR, task_id),
            f"/api/tasks/{task_id}/flagship",
            "primary_v15_conditioned_multilibrary_identification",
        ),
        (
            "v11_flagship_audit_entities",
            _v11_flagship_artifact_path(TASK_DIR, task_id),
            "",
            "read_only_v11_v14_compatible_audit_baseline",
        ),
        (
            "v12_evirank_mv_shadow",
            _evirank_artifact_path(TASK_DIR, task_id),
            f"/api/tasks/{task_id}/evirank-mv",
            "shadow_v12_multiview_entity_ranking_no_production_effect",
        ),
        (
            "evitwin_ms_authentic_ranking",
            _evitwin_artifact_path(task_id),
            f"/api/tasks/{task_id}/evitwin-ms",
            "active_authentic_reference_candidate_display_ranking",
        ),
        (
            "legacy_results",
            TASK_DIR / f"{task_id}_results.json",
            f"/api/tasks/{task_id}/results",
            "candidate_generation_and_raw_scan_evidence",
        ),
        (
            "identification_records",
            _identification_contract_path(TASK_DIR, task_id),
            f"/api/tasks/{task_id}/identification-records",
            "versioned_identification_contract",
        ),
        (
            "entity_graph",
            _entity_bundle_path(TASK_DIR, task_id),
            f"/api/tasks/{task_id}/entities",
            "candidate_independent_entity_contract",
        ),
        (
            "final_claims",
            _claim_contract_path(TASK_DIR, task_id),
            f"/api/tasks/{task_id}/claims",
            "evigate2_final_claim_contract",
        ),
        (
            "tic_overview",
            TASK_DIR / f"{task_id}_tic_overview.json",
            f"/api/tasks/{task_id}/tic_overview",
            "chromatogram_overview",
        ),
        (
            "expected_compounds",
            _expected_compounds_artifact_path(task_id),
            f"/api/tasks/{task_id}/expected-compounds",
            "independent_targeted_ms1_evidence",
        ),
        (
            "evograph_fragmentation",
            _evograph_sidecar_path(TASK_DIR, task_id),
            f"/api/tasks/{task_id}/evograph-fragmentation",
            "shadow_neutral_loss_graph_hypotheses",
        ),
    ]
    return {
        "task_id": task_id,
        "task_state": task.get("status"),
        "resources": [
            {
                "resource": name,
                "contract_role": role,
                "status": "available" if path.is_file() else "not_materialized",
                "size_bytes": path.stat().st_size if path.is_file() else 0,
                "api_url": url,
            }
            for name, path, url, role in resources
        ],
        "legacy_results_are_immutable_inputs": True,
        "final_claims_are_separate_from_legacy_confidence": True,
    }


@app.get("/api/tasks/{task_id}/evidence-passport/{spectrum_index}")
async def get_task_evidence_passport(task_id: str, spectrum_index: int):
    """Return the full evidence DAG/passport for one immutable result row."""
    results = await asyncio.to_thread(_get_task_results, task_id)
    if spectrum_index < 0 or spectrum_index >= len(results):
        raise HTTPException(status_code=404, detail="Spectrum index out of range")
    row = results[spectrum_index]
    passport = _build_evidence_passport(row)
    if passport is None:
        raise HTTPException(
            status_code=409,
            detail="This spectrum has no library candidate to audit",
        )
    return {
        "task_id": task_id,
        "spectrum_index": spectrum_index,
        "spectrum_id": str(row.get("spectrum_id") or ""),
        "candidate": {
            "compound_name": str((row.get("match") or {}).get(
                "compound_name"
            ) or ""),
            "formula": str((row.get("match") or {}).get("formula") or ""),
            "inchikey": str((row.get("match") or {}).get("inchikey") or ""),
        },
        "passport": passport,
    }


@app.get("/api/tasks/{task_id}/evidence-governance")
async def get_task_evidence_governance(
    task_id: str,
    claim_status: Optional[str] = Query(None),
    limit: int = Query(500, ge=1, le=5000),
    offset: int = Query(0, ge=0),
):
    """Return a reproducible task-level passport ledger without peak arrays."""
    results = await asyncio.to_thread(_get_task_results, task_id)
    ledger = []
    for index, row in enumerate(results):
        match = row.get("match") if isinstance(row, dict) else None
        if not isinstance(match, dict):
            continue
        passport = match.get("evidence_governance")
        if not isinstance(passport, dict):
            # Old cached result payloads may predate passport annotation.
            # Build the immutable audit view lazily so the ledger remains
            # complete without requiring a task rerun or changing rankings.
            passport = _build_evidence_passport(row)
            if not isinstance(passport, dict):
                continue
            match["evidence_governance"] = passport
        if claim_status and str(passport.get("claim_status") or "") != claim_status:
            continue
        ledger.append({
            "spectrum_index": index,
            "spectrum_id": str(row.get("spectrum_id") or ""),
            "precursor_feature_id": str(
                row.get("lcms_feature_id")
                or row.get("precursor_feature_id")
                or ""
            ),
            "precursor_mz": row.get("precursor_mz"),
            "retention_time_min": row.get("retention_time"),
            "compound_name": str(match.get("compound_name") or ""),
            "formula": str(match.get("formula") or ""),
            "inchikey": str(match.get("inchikey") or ""),
            "claim_status": passport.get("claim_status"),
            "claim_scope": passport.get("claim_scope"),
            "release_eligible_before_replication": passport.get(
                "release_eligible_before_replication"
            ),
            "release_blockers": passport.get("release_blockers") or [],
            "applicable_evidence_nodes": passport.get(
                "applicable_evidence_nodes"
            ),
            "independent_support_groups": passport.get(
                "independent_support_groups"
            ),
            "counterfactuals": passport.get("counterfactuals") or {},
            "passport_url": (
                f"/api/tasks/{task_id}/evidence-passport/{index}"
            ),
        })
    total = len(ledger)
    return {
        "task_id": task_id,
        "policy_id": EVIGATE_POLICY_ID,
        "summary": _summarize_evidence_governance(results),
        "claim_status_filter": claim_status,
        "ledger": ledger[offset:offset + limit],
        "total": total,
        "offset": offset,
        "limit": limit,
        "has_more": offset + limit < total,
        "ranking_changed": False,
    }


def _dedup_matched(matched_results: list[dict]) -> list[dict]:
    """对匹配结果按化合物身份去重，保留 final_score 最高的。

    去重键: 有效完整 InChIKey 优先，否则结构、规范名+分子式、谱图 ID 分层降级。
    完整键不同的立体异构体绝不再因名称相同而合并。同一化合物的多条谱图只保留
    排名最高的一条；完全同分时使用来源和 accession 稳定排序。
    """

    def priority(row: dict) -> tuple:
        match = row.get("match", {}) or {}
        ppm = match.get("mass_ppm")
        try:
            ppm_abs = abs(float(ppm)) if ppm is not None else 999999.0
        except (TypeError, ValueError):
            ppm_abs = 999999.0
        return (
            # Presentation ordering keeps reviewed external/background
            # chemicals in the exported audit, but never places them ahead of
            # ordinary sample constituents solely because their library cosine
            # is high.  Core search scores and Top-1 selection are unchanged.
            1 if match.get("potential_contaminant") is True else 0,
            -float(match.get("rank_score", match.get("final_score", 0.0)) or 0.0),
            -float(match.get("final_score", 0.0) or 0.0),
            -float(match.get("cosine", 0.0) or 0.0),
            ppm_abs,
            _compound_identity_key(match),
            str(match.get("dataset_source") or ""),
            str(match.get("lib_spectrum_id") or ""),
            str(row.get("spectrum_id") or ""),
        )

    best: dict[str, dict] = {}
    for r in matched_results:
        m = r.get("match", {}) or {}
        # 排除历史 ICEBERG 结构推理兜底结果(已下线): cosine/final=0.70 为占位假值非真实相似度,
        # 每条都是假 InChIKey 无法去重 → 旧任务 results 残留的 putative 结果一律剔除不导出。
        if m.get("dataset_source") == "iceberg_fallback":
            continue
        # 过滤无分子结构的无效匹配(SMILES 空/nan): MoNA 虚拟代谢转化条目等无结构无法鉴定,
        # formula 会残缺(如 C6/C8/C9), 是噪声结果 → 剔除不导出。
        _smi = m.get("smiles")
        if not _smi or _is_nan(_smi) or str(_smi).strip().lower() in ('', 'nan', 'none', 'null'):
            continue
        # 清理 compound_name 的 NaN: float('nan') truthy 绕过 `or` 兜底 → 显示 'nan'
        _cn = m.get("compound_name")
        if _cn is None or _is_nan(_cn) or str(_cn).strip().lower() in ('', 'nan', 'none', 'null'):
            m["compound_name"] = (str(m.get("inchikey") or "").split("-")[0]
                                  or f"Compound_{r.get('spectrum_id', '')}")
        normalized_ik = _normalize_inchikey(m.get("inchikey"))
        if normalized_ik:
            m["inchikey"] = normalized_ik
        key = _compound_identity_key(m)
        if not key:
            continue  # 无身份标识的罕见情况，跳过(不入去重，但仍会出现在原列表)
        prev = best.get(key)
        if prev is None or priority(r) < priority(prev):
            best[key] = r
    return sorted(best.values(), key=priority)


def _compound_confidence_summary(matched_results: list[dict]) -> dict:
    """Return explicit, non-interchangeable compound-level High counts.

    A compound can have a lower-ranked scan calibrated High while its
    highest-ranked representative scan remains Medium (for example because the
    representative has a same-formula near-tie).  Both facts are useful, but
    silently choosing one as "compound-level High" makes task summaries appear
    inconsistent.  Ranking is never changed here.
    """
    matched = [
        row for row in matched_results
        if isinstance(row, dict) and isinstance(row.get("match"), dict)
    ]
    representatives = _dedup_matched(matched)
    high_rows = [
        row for row in matched
        if str((row.get("match") or {}).get("confidence") or "") == "High"
    ]
    high_supported = _dedup_matched(high_rows)
    representative_high = [
        row for row in representatives
        if str((row.get("match") or {}).get("confidence") or "") == "High"
    ]
    supported_keys = {
        _compound_identity_key(row.get("match") or {})
        for row in high_supported
    }
    representative_high_keys = {
        _compound_identity_key(row.get("match") or {})
        for row in representative_high
    }
    supported_keys.discard("")
    representative_high_keys.discard("")
    representative_2d_connectivity_high = [
        row for row in representatives
        if (
            str(
                (row.get("match") or {}).get(
                    "two_dimensional_connectivity_confidence"
                )
                or (row.get("match") or {}).get("confidence")
                or ""
            ) == "High"
            and str(
                (row.get("match") or {}).get("identification_scope") or ""
            ) != "isomer_group"
        )
    ]
    representative_isomer_group_high = [
        row for row in representatives
        if (
            (
                str(
                    (row.get("match") or {}).get("identification_scope") or ""
                ) == "isomer_group"
                or (
                    not str(
                        (row.get("match") or {}).get("identification_scope") or ""
                    )
                    and bool((row.get("match") or {}).get("isomer_ambiguous"))
                )
            )
            and str(
                (row.get("match") or {}).get("isomer_group_confidence") or ""
            ) == "High"
        )
    ]
    representative_class_only_high = [
        row for row in representatives
        if (
            str((row.get("match") or {}).get("confidence_class") or "")
            == "High"
            and str((row.get("match") or {}).get("confidence") or "")
            != "High"
            and not (
                str(
                    (row.get("match") or {}).get("identification_scope") or ""
                ) == "isomer_group"
                and str(
                    (row.get("match") or {}).get(
                        "isomer_group_confidence"
                    ) or ""
                ) == "High"
            )
            and not (
                str(
                    (row.get("match") or {}).get(
                        "two_dimensional_connectivity_confidence"
                    ) or ""
                ) == "High"
                and str(
                    (row.get("match") or {}).get("identification_scope") or ""
                ) != "isomer_group"
            )
        )
    ]
    representative_high_reportable_scope = (
        representative_2d_connectivity_high
        + representative_isomer_group_high
    )
    return {
        "unique_library_candidates": len(representatives),
        # Explicit reader-facing count.  Keep the historical keys below for
        # regression compatibility, but never require a client to infer that
        # their unqualified "High" means full structure identity only.
        "high_reportable_scope_unique_compounds": len(
            representative_high_reportable_scope
        ),
        "representative_high_reportable_scope_unique_compounds": len(
            representative_high_reportable_scope
        ),
        "full_identity_high_supported_unique_compounds": len(high_supported),
        "representative_full_identity_high_unique_compounds": len(
            representative_high
        ),
        "full_identity_high_spectra": len(high_rows),
        "high_supported_unique_compounds": len(high_supported),
        "representative_high_unique_compounds": len(representative_high),
        "representative_2d_connectivity_high_unique_compounds": len(
            representative_2d_connectivity_high
        ),
        "representative_isomer_group_high_unique_compounds": len(
            representative_isomer_group_high
        ),
        "representative_class_only_high_unique_compounds": len(
            representative_class_only_high
        ),
        "high_spectra": len(high_rows),
        "supported_but_representative_not_high": len(
            supported_keys - representative_high_keys
        ),
        "high_supported_policy": "at_least_one_high_ms2_per_identity",
        "representative_high_policy": (
            "highest_raw_rank_score_representative_is_high"
        ),
        "high_reportable_scope_policy": (
            "representative_exact_2d_connectivity_high_plus_"
            "representative_unresolved_isomer_group_high"
        ),
        "legacy_unqualified_high_fields_scope": "full_structure_identity",
        "ranking_changed": False,
    }


_COUNTERFACTUAL_SCORE_FIELDS = (
    "cosine",
    "entropy",
    "fragment_score",
    "intensity_score",
    "reverse_match",
    "final_score",
    "rank_score",
)
_HONEST_PRESENTATION_FIELDS = (
    "exact_structure_confidence",
    "two_dimensional_connectivity_confidence",
    "isomer_group_confidence",
    "identification_scope",
    "isomer_group_formula",
    "isomer_group_size",
    "isomer_group_evidence",
    "honesty_tier",
    "isomer_ambiguous",
    "isomer_n",
    "isomer_set",
    "same_formula_alternatives",
    "same_formula_alternatives_n",
    "noncompetitive_same_formula_n",
    "isomer_ambiguity_type",
    "isomer_relationship_counts",
    "isomer_resolution_status",
    "class_consensus",
    "class_agreement",
    "honest_headline",
    "evidence_tier",
    "annotation_level",
    "rt_confirmed",
    "rt_warning",
    "candidate_comparison_source",
    "candidate_comparison_cosine",
    "candidate_margin",
    "deconvolution_score_influence",
    "deconvolution_evidence_state",
    "chimeric_isolation_warning",
)


def _raw_spectrum_counterfactual_match(match: dict) -> dict:
    """Rebuild one read-only no-deconvolution confidence counterfactual.

    Accepted deconvolution stores every replaced score as ``raw_*``.  Restoring
    those fields and rerunning the same honesty calibration quantifies the
    module's net effect without reranking Top-1, widening a margin, or mutating
    the saved result.
    """
    counterfactual = dict(match)
    for field in _COUNTERFACTUAL_SCORE_FIELDS:
        raw_value = counterfactual.get(f"raw_{field}")
        if raw_value is not None:
            counterfactual[field] = raw_value
    candidates = []
    for candidate in counterfactual.get("top_candidates") or []:
        copied = dict(candidate)
        for field in _COUNTERFACTUAL_SCORE_FIELDS:
            raw_value = copied.get(f"raw_{field}")
            if raw_value is not None:
                copied[field] = raw_value
        candidates.append(copied)
    counterfactual["top_candidates"] = candidates
    deconvolution = counterfactual.get("chimeric_deconvolution")
    if isinstance(deconvolution, dict):
        counterfactual["chimeric_deconvolution"] = {
            **deconvolution,
            "accepted": False,
            "runtime_influence": False,
            "status": "counterfactual_disabled",
            "score_contribution": 0.0,
        }
    honest = _honest_annotate(counterfactual, candidates)
    counterfactual["confidence"] = honest["compound_confidence"]
    counterfactual["confidence_class"] = honest["class_confidence"]
    for field in _HONEST_PRESENTATION_FIELDS:
        counterfactual[field] = honest.get(field)
    return counterfactual


def _deconvolution_impact_summary(results: list[dict]) -> dict:
    """Compare current results with an immutable raw-spectrum counterfactual."""
    current_rows = [
        row for row in results
        if (
            isinstance(row, dict)
            and not bool(row.get("blank_excluded"))
            and isinstance(row.get("match"), dict)
        )
    ]
    baseline_rows = []
    for row in current_rows:
        baseline_rows.append({
            "spectrum_id": row.get("spectrum_id"),
            "precursor_feature_id": row.get("precursor_feature_id"),
            "match": _raw_spectrum_counterfactual_match(row["match"]),
        })
    top1_violations = 0
    margin_violations = 0
    rank_score_violations = 0
    for current_row, baseline_row in zip(current_rows, baseline_rows):
        match = current_row["match"]
        baseline_match = baseline_row["match"]
        evidence = match.get("chimeric_deconvolution") or {}
        if not evidence.get("runtime_influence"):
            continue
        current_identity = _compound_identity_key(match)
        selected_before = str(
            evidence.get("selected_identity_before_deconvolution") or ""
        )
        if (
            evidence.get("top1_identity_preserved") is False
            or (selected_before and selected_before != current_identity)
        ):
            top1_violations += 1
        current_margin = match.get("candidate_margin")
        baseline_margin = baseline_match.get("candidate_margin")
        margin_changed = evidence.get("deconvolved_candidate_margin_used") is True
        if current_margin is None or baseline_margin is None:
            margin_changed = margin_changed or current_margin != baseline_margin
        else:
            try:
                margin_changed = margin_changed or not math.isclose(
                    float(current_margin),
                    float(baseline_margin),
                    rel_tol=0.0,
                    abs_tol=1e-6,
                )
            except (TypeError, ValueError):
                margin_changed = True
        if margin_changed:
            margin_violations += 1
        try:
            rank_changed = not math.isclose(
                float(match.get("rank_score") or 0.0),
                float(match.get("raw_rank_score") or 0.0),
                rel_tol=0.0,
                abs_tol=1e-9,
            )
        except (TypeError, ValueError):
            rank_changed = True
        if rank_changed:
            rank_score_violations += 1
    current = _compound_confidence_summary(current_rows)
    baseline = _compound_confidence_summary(baseline_rows)
    current_features = _precursor_feature_summary(current_rows)
    baseline_features = _precursor_feature_summary(baseline_rows)
    accepted = sum(
        bool(((row.get("match") or {}).get("chimeric_deconvolution") or {}).get(
            "runtime_influence"
        ))
        for row in current_rows
    )
    metric_pairs = {
        "high_supported_unique_compounds": (
            baseline["high_supported_unique_compounds"],
            current["high_supported_unique_compounds"],
        ),
        "representative_high_unique_compounds": (
            baseline["representative_high_unique_compounds"],
            current["representative_high_unique_compounds"],
        ),
        "high_precursor_features": (
            baseline_features["high_precursor_features"],
            current_features["high_precursor_features"],
        ),
        "high_spectra": (
            baseline["high_spectra"],
            current["high_spectra"],
        ),
    }
    return {
        "method": "frozen_raw_spectrum_counterfactual",
        "accepted_score_increase_spectra": accepted,
        "top1_changed": top1_violations > 0,
        "candidate_margin_changed": margin_violations > 0,
        "raw_rank_score_changed": rank_score_violations > 0,
        "invariant_violations": {
            "top1_identity": top1_violations,
            "candidate_margin": margin_violations,
            "raw_rank_score": rank_score_violations,
        },
        "baseline": {key: pair[0] for key, pair in metric_pairs.items()},
        "current": {key: pair[1] for key, pair in metric_pairs.items()},
        "delta": {key: pair[1] - pair[0] for key, pair in metric_pairs.items()},
        "interpretation": (
            "Only spectra passing every reconstruction, regularisation, target-"
            "fraction, non-redundancy and spectral-gain gate contribute."
        ),
    }


def _identification_funnel_summary(results: list[dict]) -> dict:
    """Explain scan, event, candidate and confidence counts without conflation."""
    rows = [row for row in results if isinstance(row, dict)]
    search_rows = [row for row in rows if not bool(row.get("blank_excluded"))]
    matched = [
        row for row in search_rows if isinstance(row.get("match"), dict)
    ]
    representatives = _dedup_matched(matched)
    feature_summary = _precursor_feature_summary(search_rows)
    confidence_counts = Counter(
        str((row.get("match") or {}).get("confidence") or "Unknown")
        for row in matched
    )
    representative_confidence = Counter(
        str(row["match"].get("confidence") or "Unknown")
        for row in representatives
    )
    blocker_counts: Counter[str] = Counter()
    primary_counts: Counter[str] = Counter()

    def blockers(match: dict) -> list[str]:
        found: list[str] = []
        smiles = str(match.get("smiles") or "").strip().lower()
        if match.get("has_structure") is False or smiles in {
            "", "nan", "none", "null", "n/a"
        }:
            found.append("no_validated_structure")
        try:
            cosine = float(match.get("cosine") or 0.0)
            final = float(match.get("final_score") or 0.0)
        except (TypeError, ValueError):
            cosine, final = 0.0, 0.0
        if cosine < 0.70:
            found.append("cosine_below_0_70")
        elif cosine < 0.80:
            found.append("cosine_below_0_80")
        elif cosine < 0.90:
            found.append("cosine_below_high_threshold")
        if final < 0.50:
            found.append("final_score_below_high_threshold")
        if match.get("mass_metadata_consistent") is False:
            found.append("precursor_formula_adduct_conflict")
        try:
            if (
                match.get("mass_ppm") is not None
                and abs(float(match["mass_ppm"])) > 50.0
            ):
                found.append("mass_error_above_50_ppm")
        except (TypeError, ValueError):
            pass
        try:
            if (
                match.get("reverse_match") is not None
                and float(match["reverse_match"]) < 0.60
            ):
                found.append("reverse_coverage_below_0_60")
        except (TypeError, ValueError):
            pass
        if match.get("isomer_ambiguous"):
            found.append("competitive_same_formula_isomer")
        try:
            if float(match.get("class_agreement") or 0.0) < 0.50:
                found.append("class_consensus_below_0_50")
        except (TypeError, ValueError):
            found.append("class_consensus_below_0_50")
        if match.get("ion_identity_conflict") or match.get(
            "ion_identity_adduct_conflict"
        ):
            found.append("ion_identity_conflict")
        if str(match.get("source_tier") or "").strip().lower() not in {
            "experimental", "merged", "curated", "massbank_supplement",
        }:
            found.append("library_source_confidence_cap")
        return found

    for row in representatives:
        match = row["match"]
        reportable_high = (
            (
                str(match.get("identification_scope") or "") == "isomer_group"
                and str(match.get("isomer_group_confidence") or "") == "High"
            )
            or (
                str(match.get("identification_scope") or "") != "isomer_group"
                and str(
                    match.get("two_dimensional_connectivity_confidence")
                    or match.get("confidence")
                    or ""
                ) == "High"
            )
        )
        if reportable_high:
            continue
        reasons = blockers(match)
        if not reasons:
            reasons = ["other_calibration_guard"]
        blocker_counts.update(reasons)
        primary_counts[reasons[0]] += 1

    unmatched = [
        row for row in search_rows if not isinstance(row.get("match"), dict)
    ]
    unmatched_peak_bins = Counter()
    for row in unmatched:
        count = len(row.get("query_peaks") or [])
        label = (
            "0" if count == 0
            else "1-2" if count <= 2
            else "3-5" if count <= 5
            else "6-10" if count <= 10
            else "11+"
        )
        unmatched_peak_bins[label] += 1
    return {
        "method": "explicit_identification_unit_funnel",
        "ms2_scans": len(rows),
        "blank_excluded_ms2_scans": len(rows) - len(search_rows),
        "ms2_precursor_events_approximate": feature_summary["precursor_features"],
        "matched_ms2_scans": len(matched),
        "unmatched_ms2_scans": len(unmatched),
        "matched_ms2_precursor_events": feature_summary[
            "matched_precursor_features"
        ],
        "unique_library_candidates": len(representatives),
        "high_supported_unique_identities": len(_dedup_matched([
            row for row in matched
            if str(row["match"].get("confidence") or "") == "High"
        ])),
        "representative_exact_high_identities": representative_confidence["High"],
        "representative_2d_connectivity_high_identities": sum(
            str(
                row["match"].get("two_dimensional_connectivity_confidence")
                or row["match"].get("confidence")
                or ""
            ) == "High"
            and str(row["match"].get("identification_scope") or "")
            != "isomer_group"
            for row in representatives
        ),
        "representative_isomer_group_high_identities": sum(
            str(row["match"].get("identification_scope") or "")
            == "isomer_group"
            and str(
                row["match"].get("isomer_group_confidence") or ""
            ) == "High"
            for row in representatives
        ),
        "representative_class_high_identities": sum(
            str(row["match"].get("confidence_class") or "") == "High"
            for row in representatives
        ),
        "representative_class_only_high_identities": sum(
            str(row["match"].get("confidence_class") or "") == "High"
            and str(row["match"].get("confidence") or "") != "High"
            and not (
                str(row["match"].get("identification_scope") or "")
                == "isomer_group"
                and str(
                    row["match"].get("isomer_group_confidence") or ""
                ) == "High"
            )
            and not (
                str(
                    row["match"].get(
                        "two_dimensional_connectivity_confidence"
                    ) or ""
                ) == "High"
                and str(row["match"].get("identification_scope") or "")
                != "isomer_group"
            )
            for row in representatives
        ),
        "representative_isomer_ambiguous_identities": sum(
            bool(row["match"].get("isomer_ambiguous"))
            for row in representatives
        ),
        "matched_spectrum_confidence_counts": dict(confidence_counts),
        "representative_confidence_counts": dict(representative_confidence),
        "non_high_primary_blockers": dict(primary_counts),
        "non_high_all_blockers": dict(blocker_counts),
        "unmatched_query_peak_bins": dict(unmatched_peak_bins),
        "unmatched_with_nonzero_query_peaks": sum(
            len(row.get("query_peaks") or []) > 0 for row in unmatched
        ),
        "unmatched_reason_scope": (
            "No library candidate passed configured retrieval and score "
            "thresholds; absence of a match is not absence of an MS1 feature."
        ),
        "units_are_not_interchangeable": True,
    }


def _effective_ion_identity_summary(summary: dict | None) -> dict:
    """Repair applicability status in legacy cached IIMN summaries."""
    effective = dict(summary) if isinstance(summary, dict) else {}
    extraction = effective.get("trace_extraction")
    extraction = extraction if isinstance(extraction, dict) else {}
    extraction_status = str(extraction.get("status") or "")
    current_status = str(effective.get("status") or "not_run")
    if (
        extraction_status
        and extraction_status not in {"ok", "provided_traces"}
        and current_status in {"ok", "not_run"}
    ):
        effective["status"] = extraction_status
    elif (
        extraction_status in {"ok", "provided_traces"}
        and current_status == "ok"
        and int(effective.get("usable_trace_count") or 0) == 0
    ):
        effective["status"] = "no_usable_ms1_traces"
    return effective


def _analysis_module_availability(
    params: dict | None,
    ion_identity_summary: dict | None,
    results: list[dict] | None,
) -> dict:
    """Describe whether MS1-dependent modules were actually applicable.

    Enabled is not the same as executed.  DecoID-like deconvolution needs MS1
    isolation context and IIMN needs continuous MS1 XICs.  This summary keeps a
    missing-MS1 input from being presented as a successful run with zero
    findings.
    """
    params = params if isinstance(params, dict) else {}
    rows = [row for row in (results or []) if isinstance(row, dict)]
    ion_summary = _effective_ion_identity_summary(ion_identity_summary)
    extraction = ion_summary.get("trace_extraction")
    extraction = extraction if isinstance(extraction, dict) else {}
    trace_status = str(extraction.get("status") or ion_summary.get("status") or "")
    raw_ms1_count = extraction.get("ms1_scan_count")
    try:
        ms1_scan_count = int(raw_ms1_count) if raw_ms1_count is not None else None
    except (TypeError, ValueError):
        ms1_scan_count = None

    isolation_status_counts = Counter()
    deconvolution_status_counts = Counter()
    for row in rows:
        isolation = row.get("precursor_isolation_evidence")
        if isinstance(isolation, dict):
            status = str(isolation.get("status") or "")
            if status:
                isolation_status_counts[status] += 1
        match = row.get("match")
        deconvolution = (
            match.get("chimeric_deconvolution")
            if isinstance(match, dict) else None
        )
        if isinstance(deconvolution, dict):
            status = str(deconvolution.get("status") or "")
            if status:
                deconvolution_status_counts[status] += 1

    ms1_positive_isolation_statuses = {"clean", "mixed", "severe"}
    no_ms1_isolation_statuses = {"referenced_ms1_not_in_file"}
    if ms1_scan_count is not None:
        continuous_ms1_available: bool | None = ms1_scan_count > 0
    elif any(
        isolation_status_counts.get(status, 0)
        for status in ms1_positive_isolation_statuses
    ):
        continuous_ms1_available = True
    elif any(
        isolation_status_counts.get(status, 0)
        for status in no_ms1_isolation_statuses
    ) or trace_status == "no_ms1_scans":
        continuous_ms1_available = False
        ms1_scan_count = 0
    else:
        continuous_ms1_available = None

    iimn_requested = bool(params.get(
        "enable_ion_identity_networking", USE_ION_IDENTITY_NETWORKING
    ))
    usable_trace_count = int(ion_summary.get("usable_trace_count") or 0)
    if not iimn_requested:
        iimn_status = "disabled"
        iimn_applicable = False
        iimn_executed = False
    elif continuous_ms1_available is False or trace_status == "no_ms1_scans":
        iimn_status = "unavailable_no_ms1_scans"
        iimn_applicable = False
        iimn_executed = False
    elif str(ion_summary.get("status") or "") == "no_usable_ms1_traces":
        iimn_status = "unavailable_no_usable_ms1_traces"
        iimn_applicable = False
        iimn_executed = False
    elif str(ion_summary.get("status") or "") == "ok":
        iimn_status = "evaluated"
        iimn_applicable = True
        iimn_executed = True
    else:
        detail = str(ion_summary.get("status") or "not_run")
        iimn_status = f"unavailable_{detail}"
        iimn_applicable = False
        iimn_executed = False

    deconvolution_requested = bool(params.get(
        "enable_chimeric_deconvolution", USE_CHIMERIC_DECONVOLUTION
    ))
    accepted_deconvolution = int(
        deconvolution_status_counts.get("accepted_score_applied", 0)
    )
    mixed_or_severe = sum(
        isolation_status_counts.get(status, 0) for status in ("mixed", "severe")
    )
    solver_precheck_statuses = {
        "disabled",
        "module_unavailable",
        "isolation_evidence_unavailable",
        "clean_isolation_no_deconvolution",
        "isolation_window_unavailable",
        "no_coisolated_library_component",
        "evaluation_failed",
    }
    solver_runs = sum(
        count for status, count in deconvolution_status_counts.items()
        if status not in solver_precheck_statuses
    )
    if not deconvolution_requested:
        deconvolution_status = "disabled"
        deconvolution_applicable = False
        deconvolution_executed = False
    elif continuous_ms1_available is False:
        deconvolution_status = "unavailable_no_ms1_isolation_context"
        deconvolution_applicable = False
        deconvolution_executed = False
    elif not deconvolution_status_counts:
        deconvolution_status = "not_evaluated_no_library_candidates"
        deconvolution_applicable = continuous_ms1_available is True
        deconvolution_executed = False
    elif mixed_or_severe:
        deconvolution_status = "evaluated"
        deconvolution_applicable = True
        deconvolution_executed = solver_runs > 0
    elif isolation_status_counts.get("clean", 0):
        deconvolution_status = "evaluated_no_chimeric_isolation"
        deconvolution_applicable = True
        deconvolution_executed = False
    elif any(
        isolation_status_counts.get(status, 0)
        for status in (
            "isolation_window_unavailable",
            "isolation_window_invalid",
            "no_ms1_signal_in_isolation_window",
        )
    ) or deconvolution_status_counts.get("isolation_window_unavailable", 0):
        deconvolution_status = "unavailable_isolation_window_evidence"
        deconvolution_applicable = False
        deconvolution_executed = False
    else:
        deconvolution_status = "not_evaluated_unknown_isolation_state"
        deconvolution_applicable = False
        deconvolution_executed = False

    warnings = []
    requested_ms1_modules = []
    feature_detection_requested = bool(
        params.get("enable_ms1_feature_detection", True)
    )
    linked_lcms_rows = sum(bool(row.get("lcms_feature_id")) for row in rows)
    linked_lcms_features = len({
        str(row.get("lcms_feature_id"))
        for row in rows if row.get("lcms_feature_id")
    })
    if feature_detection_requested:
        requested_ms1_modules.append("连续 MS1 feature finding")
    if deconvolution_requested:
        requested_ms1_modules.append("DecoID 类嵌合谱去卷积")
    if iimn_requested:
        requested_ms1_modules.append("IIMN")
    if continuous_ms1_available is False and requested_ms1_modules:
        warnings.append({
            "code": "continuous_ms1_missing",
            "severity": "warning",
            "modules": requested_ms1_modules,
            "message": (
                f"当前输入含 {len(rows)} 张 MS2，但连续 MS1 扫描为 0；"
                f"{'、'.join(requested_ms1_modules)}未实际执行。"
            ),
            "action": (
                "请从原始采集文件重新导出同时包含 MS1/MS2 的 mzML，"
                "不要使用仅保留 MS2 的过滤器；若原始采集本身无 survey MS1，"
                "则需重新采集。仅填写隔离窗宽度不能替代缺失的 MS1。"
            ),
        })
    elif iimn_requested and iimn_status == "unavailable_no_usable_ms1_traces":
        warnings.append({
            "code": "iimn_no_usable_ms1_traces",
            "severity": "warning",
            "modules": ["IIMN"],
            "message": "检测到 MS1，但目标前体没有足够的非零 XIC 点，IIMN 未构网。",
            "action": "请检查色谱峰覆盖、RT 单位、MS1 扫描密度及 XIC 质量。",
        })
    if (
        deconvolution_requested
        and continuous_ms1_available is not False
        and deconvolution_status == "unavailable_isolation_window_evidence"
    ):
        warnings.append({
            "code": "deconvolution_isolation_evidence_missing",
            "severity": "warning",
            "modules": ["DecoID 类嵌合谱去卷积"],
            "message": "缺少可审计的 MS1 隔离纯度/隔离窗证据，去卷积未运行。",
            "action": "请保留 mzML precursor isolation window 元数据，或填写真实的全隔离窗宽度。",
        })

    return {
        "input_ms2_spectra": len(rows),
        "ms1_scan_count": ms1_scan_count,
        "continuous_ms1_available": continuous_ms1_available,
        "ms1_feature_detection": {
            "requested": feature_detection_requested,
            "applicable": (
                feature_detection_requested
                and continuous_ms1_available is not False
            ),
            "executed_with_ms2_links": linked_lcms_rows > 0,
            "status": (
                "disabled"
                if not feature_detection_requested
                else "unavailable_no_ms1_scans"
                if continuous_ms1_available is False
                else "evaluated_with_ms2_links"
                if linked_lcms_rows
                else "evaluated_no_ms2_links_or_backend_unavailable"
            ),
            "linked_ms2_spectra": linked_lcms_rows,
            "lcms_features_with_ms2": linked_lcms_features,
        },
        "iimn": {
            "requested": iimn_requested,
            "applicable": iimn_applicable,
            "executed": iimn_executed,
            "status": iimn_status,
            "detail_status": str(ion_summary.get("status") or "not_run"),
            "usable_trace_count": usable_trace_count,
            "network_count": int(ion_summary.get("network_count") or 0),
        },
        "chimeric_deconvolution": {
            "requested": deconvolution_requested,
            "applicable": deconvolution_applicable,
            "executed": deconvolution_executed,
            "status": deconvolution_status,
            "solver_runs": int(solver_runs),
            "accepted_spectra": accepted_deconvolution,
            "mixed_or_severe_spectra": int(mixed_or_severe),
            "status_counts": dict(deconvolution_status_counts),
        },
        "isolation_status_counts": dict(isolation_status_counts),
        "warnings": warnings,
    }


def _display_class(m: dict, language: str = "zh") -> str:
    """Return a reader-facing class without promoting a model fallback."""
    from web.report_common import classification_evidence_status

    evidence = classification_evidence_status(m)
    label = str(evidence["label"])
    if evidence["status"] == "unresolved":
        return "Classification unresolved" if language == "en" else "类别尚未可靠确定"
    if evidence["status"] == "tentative":
        return f"Tentative: {label}" if language == "en" else f"待核查：{label}"
    return label


def _dataset_source(m: dict) -> str:
    """数据集来源 token (数据集名[-极性])。旧 results 无此字段 → 空串。"""
    return str(m.get("dataset_source") or "")


def _fmt_csv_number(value, digits: int = 4) -> str:
    """Format optional numeric audit values without turning missing data into zero."""
    if value is None:
        return ""
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return ""


def _readable_dataset_text(value: Any) -> str:
    values = value if isinstance(value, (list, tuple, set)) else [value]
    output = []
    for item in values:
        text = str(item or "").strip()
        if not text:
            continue
        for prefix in ("user_contrib:", "literature:", "database:"):
            if text.startswith(prefix):
                text = text[len(prefix):]
                break
        if text not in output:
            output.append(text)
    return " | ".join(output)


def _yes_no(value: Any, language: str = "zh") -> str:
    if language == "en":
        return "Yes" if bool(value) else "No"
    return "是" if bool(value) else "否"


MS2_RESULT_LABELS = {
    "zh": {
        "Compound_Name": "化合物",
        "Formula": "分子式",
        "Record_Type": "记录类型",
        "Identification_Group_ID": "鉴定组 ID",
        "Candidate_Record_Status": "候选记录状态",
        "Candidate_Rank": "组内候选排名",
        "Candidate_Role": "组内候选角色",
        "Is_Ranking_Lead": "是否当前最佳候选",
        "Candidate_Compound_Name": "候选化合物名称",
        "Candidate_Formula": "候选分子式",
        "Candidate_InChIKey": "候选 InChIKey",
        "Candidate_SMILES": "候选 SMILES",
        "Candidate_Library_Adduct": "候选谱库加合物",
        "Candidate_Cosine": "候选 MS/MS 相似度",
        "Candidate_Entropy": "候选谱熵相似度",
        "Candidate_Reverse_Match": "候选反向匹配覆盖",
        "Candidate_Mass_PPM": "候选质量误差 (ppm)",
        "Candidate_Final_Score": "候选综合匹配分",
        "Candidate_Rank_Score": "候选排序分",
        "Candidate_Confidence": "候选自身置信度",
        "Candidate_Classification": "候选化学类别",
        "Candidate_Relationship": "与当前最佳候选关系",
        "Candidate_Dataset_Source": "候选谱库来源",
        "Candidate_Library_Spectrum": "候选馆藏谱图编号",
        "Isomer_Competition_Status": "异构体竞争状态",
        "Isomer_Candidate_Count": "异构体组候选总数",
        "Isomer_Competitive_Candidate_Count": "其余竞争候选数",
        "Isomer_Best_Candidate": "当前最佳候选",
        "Isomer_Competing_Candidates": "其余异构体竞争候选",
        "Isomer_Candidate_Ranking_Details": "异构体候选完整排名与证据",
        "Isomer_Candidate_Details_JSON": "异构体候选机器可读明细 (JSON)",
        "Isomer_Group_Confidence": "异构体组置信度",
        "Isomer_Resolution_Status": "异构体解析状态",
        "Adduct": "加合物",
        "Adduct_Source": "加合物依据",
        "Adduct_Inference_Error_PPM": "加合物推断误差 (ppm)",
        "Ion_Mode": "离子模式",
        "Precursor_mz": "MS2 隔离目标 m/z",
        "Search_Precursor_mz": "检索用实测/回退 m/z",
        "Search_Precursor_Source": "检索质量轴来源",
        "Selected_to_Search_Delta_PPM": "隔离目标至检索质量差 (ppm)",
        "MS1_Feature_Used_For_Search": "使用连续 MS1 特征质心检索",
        "Precursor_Charge": "前体电荷",
        "Retention_Time": "RT (min)",
        "Query_Collision_Energy": "样本 CE",
        "Query_Collision_Energy_Unit": "样本 CE 单位",
        "Query_Instrument_Type": "样本仪器类型",
        "Query_Instrument": "样本仪器",
        "Confidence": "可报告身份置信度",
        "Reportable_Identity_Scope": "可报告身份层级",
        "Reportable_Identity_Label": "可报告身份置信度说明",
        "Final_Claim_Status": "鉴定结论状态码",
        "Final_Claim_Scope": "当前结论层级",
        "Final_Claim_Released": "可作结论",
        "Project_Q_Value": "项目级 q 值",
        "Calibration_Status": "项目校准状态",
        "Independent_Injection_Count": "独立进样数",
        "Final_Claim_Reason_Codes": "鉴定结论条件理由代码",
        "EviGate2_Policy_ID": "技术规则版本",
        "Final_Score": "综合匹配分",
        "Cosine_Score": "MS/MS 相似度",
        "Mass_PPM": "质量误差 (ppm)",
        "Fragment_Count": "匹配碎片数",
        "MS_MS_Fragments": "匹配碎片 m/z",
        "Matched_Fragment_Details": "匹配碎片 (m/z:强度)",
        "MS_MS_Spectrum": "样本 MS2 峰 (m/z:强度)",
        "Peaks_Count": "用于比对的 MS2 峰数",
        "Classification": "化学类别",
        "Display_Classification": "展示类别",
        "Classification_Source": "分类来源",
        "Display_Classification_Source": "展示分类来源",
        "Classification_Reliable": "分类证据可靠",
        "Classification_Status": "分类证据状态",
        "Classification_Reliability_Reason": "分类可靠性说明",
        "Classification_Confidence": "分类置信度",
        "Compound_Role": "来源角色",
        "Potential_Contaminant": "可能外源或实验背景物",
        "Contaminant_Category": "外源物类别",
        "Natural_Product_Candidate": "已策展天然产物",
        "Fragmentation_Evidence_Type": "碎裂证据类型",
        "Feature_Fragment_Status": "特征碎片计分状态",
        "Feature_Fragment_Independent_Count": "独立特征碎裂证据数",
        "Feature_Fragment_Present_Count": "已观测独立特征碎片数",
        "Feature_Fragment_Score_Nats": "特征碎片证据分 (lnLR, nats)",
        "Feature_Fragment_Score_Bits": "特征碎片证据分 (bits)",
        "Feature_Fragment_Score_Decibans": "特征碎片证据分 (decibans)",
        "Feature_Fragment_Production_Score_Contribution": "特征碎片生产排序贡献",
        "Candidate_Pathway_Eligible_For_Fusion": "候选通路可参与融合",
        "Candidate_Pathway_Runtime_Influence": "候选通路运行时影响",
        "Dataset_Source": "谱库来源",
        "Source_Tier": "谱库证据类型",
        "Library_Spectrum": "馆藏谱图编号",
        "Library_Reported_Adduct": "馆藏原始加合物",
        "Library_Collision_Energy": "馆藏 CE",
        "Library_Collision_Energy_Unit": "馆藏 CE 单位",
        "Library_CE_Metadata_Status": "馆藏 CE 元数据状态",
        "Library_Instrument_Type": "馆藏仪器类型",
        "Library_Instrument": "馆藏仪器",
        "InChIKey": "InChIKey",
        "SMILES": "SMILES",
        "Spectrum": "样本谱图编号",
    },
    "en": {},
}

MS2_RESULT_VALUE_LABELS = {
    "zh": {
        "Adduct": {
            "Unresolved": "未能可靠判定",
        },
        "Record_Type": {
            "single_identification": "单一鉴定结果",
            "isomer_candidate": "异构体组逐候选记录",
        },
        "Isomer_Competition_Status": {
            "unresolved_competition": "存在未解析的异构体竞争候选",
            "unresolved_candidate_records_incomplete": "异构体组未解析，但旧任务候选明细不完整",
            "not_applicable": "不适用（无未解析异构体竞争）",
        },
        "Adduct_Source": {
            "lcms_ion_identity": "样品离子身份网络",
            "explicit_query_metadata": "样品显式元数据",
            "library_candidate_hypothesis": "谱库候选离子形式（非样品独立测定）",
            "candidate_structure_formula_precursor_hypothesis": (
                "候选结构／分子式与前体质量共同支持（候选依赖）"
            ),
            "candidate_formula_precursor_unique_10ppm": (
                "候选分子式 + 前体 m/z 唯一支持（10 ppm；候选依赖）"
            ),
            "unsupported_precursor_charge": "当前电荷态不支持安全推断",
            "missing_precursor_mz": "缺少前体 m/z，未猜测",
            "missing_formula_or_polarity": "缺少分子式或极性，未猜测",
            "no_unique_formula_mz_solution": "无法唯一判定（未猜测）",
        },
        "Search_Precursor_Source": {
            "measured_continuous_ms1_feature_centroid": (
                "连续 MS1 特征实测质心"
            ),
            "ms2_selected_ion_mz_fallback": (
                "mzML selected-ion / 隔离目标回退"
            ),
        },
        "Reportable_Identity_Scope": {
            "exact_structure": "完整身份",
            "exact_2d_connectivity": "单一二维连接结构",
            "isomer_group": "未解异构体组",
            "compound_class": "化学类别",
            "putative_candidate": "谱库候选",
            "spectral_lead": "谱图线索",
            "literature_target_lead": "文献线索",
        },
        "Classification_Status": {
            "reliable": "可靠",
            "tentative": "待核查",
            "unresolved": "未可靠确定",
        },
        "Library_CE_Metadata_Status": {
            "reported": "已报告（含单位）",
            "unit_missing": "已报告但单位缺失",
        },
        "Compound_Role": {
            "external_contaminant": "可能外源物",
            "curated_contaminant": "可能外源物",
            "industrial_contaminant": "可能工业来源物",
            "environmental_contaminant": "可能环境来源物",
            "curated_natural": "已策展天然产物",
            "unclassified": "来源未分类",
        },
        "Contaminant_Category": {
            "alkylbenzene_sulfonate_surfactant": "烷基苯磺酸类表面活性剂",
            "amine_oxide_surfactant": "叔胺 N-氧化物表面活性剂",
            "phthalate_plasticizer": "邻苯二甲酸酯类增塑剂",
            "synthetic_phenolic_antioxidant": "合成酚类抗氧化剂",
        },
        "Final_Claim_Released": {"Yes": "是", "No": "否"},
        "Classification_Reliable": {"Yes": "是", "No": "否"},
        "Candidate_Pathway_Eligible_For_Fusion": {"Yes": "是", "No": "否"},
        "Candidate_Pathway_Runtime_Influence": {"Yes": "是", "No": "否"},
        "MS1_Feature_Used_For_Search": {"Yes": "是", "No": "否"},
    },
}


def _localized_ms2_rows(
    rows: list[dict[str, Any]],
    language: str,
) -> list[dict[str, Any]]:
    mappings = MS2_RESULT_VALUE_LABELS.get(language, {})
    localized = []
    for source in rows:
        row = dict(source)
        for column, value in row.items():
            if isinstance(value, bool):
                row[column] = (
                    "是" if value else "否"
                ) if language == "zh" else (
                    "Yes" if value else "No"
                )
        for column, values in mappings.items():
            value = row.get(column)
            if value in values:
                row[column] = values[value]
        localized.append(row)
    return localized


def _ms2_column_labels(language: str = "zh") -> list[str]:
    from web.report_common import RESULTS_COLUMN_NAMES

    if language == "en":
        return list(RESULTS_COLUMN_NAMES)
    labels = MS2_RESULT_LABELS["zh"]
    return [labels[column] for column in RESULTS_COLUMN_NAMES]


def _ms2_unified_column_labels(language: str = "zh") -> list[str]:
    from web.report_common import MS2_UNIFIED_COLUMN_NAMES

    if language == "en":
        return list(MS2_UNIFIED_COLUMN_NAMES)
    labels = MS2_RESULT_LABELS["zh"]
    return [labels[column] for column in MS2_UNIFIED_COLUMN_NAMES]


def _focused_csv_text(
    rows: list[dict[str, Any]],
    columns: list[str],
    labels: list[str],
) -> str:
    if len(columns) != len(labels):
        raise ValueError(
            "CSV schema error: column keys and display labels differ in length "
            f"({len(columns)} != {len(labels)})"
        )

    def safe_cell(value: Any) -> Any:
        if not isinstance(value, str):
            return value
        # Standard result tables do not require physical line breaks inside a
        # cell.  Normalising all line separators also protects spreadsheet
        # programs that inconsistently handle quoted multiline CSV fields.
        text = re.sub(r"[\r\n\u2028\u2029]+", " ", value).strip()
        # Excel's hard cell limit is 32,767 characters.  Keep headroom for an
        # explicit note instead of allowing the tail of a spectrum to appear
        # as a new worksheet row.
        limit = 32000
        if len(text) > limit:
            marker = (
                " …[单元格内容过长，标准表已截断；"
                "完整数据请查阅详细审计或任务原始数据]"
            )
            text = text[: limit - len(marker)] + marker
        return text

    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow([safe_cell(label) for label in labels])
    for index, row in enumerate(rows, 1):
        if not isinstance(row, dict):
            raise TypeError(f"CSV row {index} is not a mapping")
        values = [safe_cell(row.get(column, "")) for column in columns]
        writer.writerow(values)

    text = buffer.getvalue()
    # Parse back with the standard CSV reader before download.  This catches a
    # future quoting or schema regression at the producer instead of leaving a
    # user to discover shifted columns in Excel.
    parsed = list(csv.reader(io.StringIO(text)))
    bad = [
        (index, len(row))
        for index, row in enumerate(parsed, 1)
        if len(row) != len(columns)
    ]
    if bad:
        raise ValueError(f"CSV parse-back validation failed: {bad[:5]}")
    return text


def _targeted_csv_text(
    expected_payload: dict,
    *,
    language: str = "zh",
) -> str:
    rows = expected_feature_rows(expected_payload, language=language)
    return _focused_csv_text(
        rows,
        EXPECTED_FEATURE_COLUMNS,
        expected_feature_column_labels(language),
    )


def _targeted_feature_audit_csv_text(
    expected_payload: dict,
    *,
    language: str = "zh",
) -> str:
    rows = expected_raw_feature_rows(
        expected_payload, language=language
    )
    return _focused_csv_text(
        rows,
        EXPECTED_RAW_FEATURE_COLUMNS,
        expected_raw_feature_column_labels(language),
    )


def _ms2_csv_text(
    matched: list[dict],
    *,
    language: str = "zh",
    compounds: list[dict[str, Any]] | None = None,
    task_id: str = "",
    results: list[dict] | None = None,
    candidate_view: str = "all",
) -> str:
    from web.lossless_peak_export import prepare_reader_peak_cells
    from web.report_common import (
        MS2_READER_COLUMN_NAMES,
        ms2_reader_column_labels,
        ms2_reader_rows,
    )

    rows = ms2_reader_rows(
        matched,
        compounds=compounds,
        language=language,
    )
    if candidate_view != "all":
        from web.plant_candidate_export import plant_reader_views
        plant_rows, external_rows = plant_reader_views(
            matched, compounds=compounds, language=language,
        )
        rows = plant_rows if candidate_view == "plant_candidates" else external_rows
    rows = _localized_ms2_rows(rows, language)
    rows = prepare_reader_peak_cells(
        rows, task_id=task_id, results=results, language=language,
    )
    csv_text = _focused_csv_text(
        rows,
        MS2_READER_COLUMN_NAMES,
        ms2_reader_column_labels(language),
    )
    if results is not None and candidate_view == "all":
        from web.governance_blocked_projection import (
            append_governance_blocked_csv_section, project_governance_blocked_candidates,
        )
        csv_text = append_governance_blocked_csv_section(
            csv_text, project_governance_blocked_candidates(results), language=language,
        )
    return csv_text


def _isomer_candidate_csv_text(
    matched: list[dict],
    *,
    language: str = "zh",
    compounds: list[dict[str, Any]] | None = None,
) -> str:
    """Return the normalised one-row-per-candidate isomer competition table."""
    from web.report_common import (
        ISOMER_CANDIDATE_EXPORT_COLUMNS,
        isomer_candidate_export_labels,
        isomer_candidate_export_rows,
        isomer_candidate_export_rows_from_compounds,
    )

    rows = (
        isomer_candidate_export_rows_from_compounds(
            compounds, language=language
        )
        if compounds is not None
        else isomer_candidate_export_rows(matched, language=language)
    )
    rows = _localized_ms2_rows(rows, language)
    return _focused_csv_text(
        rows,
        ISOMER_CANDIDATE_EXPORT_COLUMNS,
        isomer_candidate_export_labels(language),
    )


def _analysis_summary_rows(
    task: dict,
    results: list[dict],
    matched: list[dict],
    matched_clean: list[dict],
    blank_excluded: list[dict],
    blank_evidence: list[dict],
    expected_summary: dict,
    confidence_summary: dict,
    compound_summary: dict | None = None,
    *,
    language: str = "zh",
) -> list[dict[str, Any]]:
    compound_evidence = dict(
        expected_summary.get("Compound_Evidence_Counts") or {}
    )
    module_warning_count = len(
        analysis_warning_rows(task.get("analysis_warnings") or [])
    )
    compound_summary = dict(compound_summary or {})
    legacy_high = int(
        confidence_summary.get("high_supported_unique_compounds", 0) or 0
    )
    high_scope = int(
        compound_summary.get(
            "high_confidence_scope_group_count", legacy_high
        ) or 0
    )
    high_exact = int(
        compound_summary.get(
            "high_confidence_2d_connectivity_count",
            compound_summary.get(
                "high_confidence_detected_count",
                compound_summary.get(
                    "high_confidence_exact_structure_count", legacy_high
                ),
            ),
        ) or 0
    )
    high_isomer = int(
        compound_summary.get(
            "high_confidence_isomer_group_count",
            max(0, high_scope - high_exact),
        ) or 0
    )
    unverified_blank = int(
        compound_summary.get("unverified_blank_warning_count", 0) or 0
    )
    if language == "en":
        return [
            {"Item": "START HERE", "Value": "00_Confidence_Guide", "Note": "Plain-language rules for Confirmed, High, Medium, Low, identity scope and every main score"},
            {"Item": "READ FIRST", "Value": "02_Target_High", "Note": "Targeted database hits supported by an experimental library or measured diagnostic MS2"},
            {"Item": "READ SECOND", "Value": "03_MS2_Library", "Note": "Reader-facing MS2 results; every unresolved isomer candidate is already expanded in this sheet"},
            {"Item": "OPTIONAL ISOMER AUDIT", "Value": "MS2_Isomer_Candidates", "Note": "Candidate-only audit view; not required to read the complete MS2 result"},
            {"Item": "READ THIRD", "Value": "04_Target_Isotope", "Note": "Targeted candidates with isotope support; not confirmations"},
            {"Item": "READ LAST", "Value": "05_Target_MassOnly", "Note": "Exact-mass candidates only; do not treat them as confirmed compounds"},
            {"Item": "Sample file", "Value": task.get("filename", ""), "Note": ""},
            {"Item": "MS1 survey scans", "Value": expected_summary.get("MS1_Scan_Count", 0), "Note": ""},
            {"Item": "MS2 spectra", "Value": task.get("total_spectra", len(results)), "Note": ""},
            {"Item": "MS2 library matches", "Value": task.get("matched_spectra", 0), "Note": "Spectrum-level count"},
            {"Item": "MS2 match rate (%)", "Value": round(task.get("match_rate", 0), 1), "Note": ""},
            {"Item": "Deduplicated MS2 results", "Value": len(matched), "Note": "Best spectrum retained per candidate"},
            {"Item": "High-evidence identification groups", "Value": high_scope, "Note": "High single 2D connectivities plus unresolved High isomer groups"},
            {"Item": "High single 2D connectivities", "Value": high_exact, "Note": "One 2D connectivity is supportable; full stereochemical identity is reported separately"},
            {"Item": "High unresolved isomer groups", "Value": high_isomer, "Note": "High evidence supports the group, not one unique isomer"},
            {"Item": "Unverified blank-association warnings", "Value": unverified_blank, "Note": "Signals also observed in polarity-matched blanks whose role/method pairing was not verified; warning only, not automatic exclusion"},
            {"Item": "Targeted-search status", "Value": expected_summary.get("Status_Label", ""), "Note": "Independent continuous-MS1 targeted layer"},
            {"Item": "Target database", "Value": _readable_dataset_text(str(expected_summary.get("Selected_Datasets") or "").split(" | ")), "Note": ""},
            {"Item": "Database target rows", "Value": expected_summary.get("Target_Catalog_Rows", 0), "Note": ""},
            {"Item": "Targeted compound/ambiguity summaries", "Value": expected_summary.get("Compound_Summary_Rows", 0), "Note": "One routine row per identity or unresolved ambiguity group"},
            {"Item": "Raw targeted XIC peak segments", "Value": expected_summary.get("Detected_Features", 0), "Note": "See 10_Target_Peak_Audit or the separate targeted peak-audit CSV"},
            {"Item": "Targeted MS2-supported candidate groups", "Value": compound_evidence.get("EXPERIMENTAL_LIBRARY_CONFIRMED", 0) + compound_evidence.get("MS1_MS2_DIAGNOSTIC_PROBABLE", 0), "Note": "Experimental-library candidate or reported diagnostic-ion screening support; not standard-confirmed"},
            {"Item": "Targeted isotope-supported groups", "Value": compound_evidence.get("MS1_ISOTOPE_SUPPORTED", 0), "Note": "Supportive but not confirmatory"},
            {"Item": "Targeted exact-mass-only groups", "Value": compound_evidence.get("MS1_EXACT_MASS_CANDIDATE", 0), "Note": "Requires additional evidence"},
            {"Item": "Targeted experimental-library candidate groups", "Value": compound_evidence.get("EXPERIMENTAL_LIBRARY_CONFIRMED", 0), "Note": "MS/MS library candidate support; not authentic-standard confirmation"},
            {"Item": "Counting rule", "Value": "Independent", "Note": "Targeted MS1 rows are not added to ordinary MS2 library-match counts"},
            {"Item": "MS1 module applicability notes", "Value": module_warning_count, "Note": "Full text is isolated in 11_Module_Warnings when non-zero; it is not repeated across result rows"},
            {"Item": "Blank-filtered spectra", "Value": len(blank_excluded), "Note": "See appendix sheet when non-zero"},
            {"Item": "Clean-view MS2 results", "Value": len(matched_clean), "Note": "Blank-related rows removed"},
            {"Item": "Blank-evidence rows", "Value": len(blank_evidence), "Note": "See appendix sheet when non-zero"},
        ]
    return [
        {"项目": "开始前：先读懂评分", "结果": "00_置信度说明", "说明": "用通俗语言说明已确证、High、Medium、Low、身份范围和主要评分能支持什么、不能说明什么"},
        {"项目": "第一步：先看", "结果": "02_靶向_高证据", "说明": "实验谱库或实测诊断 MS2 支持的数据库靶向结果"},
        {"项目": "第二步：再看", "结果": "03_MS2谱库鉴定", "说明": "面向阅读的 MS2 结果；未解异构体的最佳候选与全部竞争候选已在本表逐行展开"},
        {"项目": "可选：异构体审计", "结果": "MS2_异构体候选明细", "说明": "仅候选审计视图；完整 MS2 结果无需切换工作表即可阅读"},
        {"项目": "第三步：辅助看", "结果": "04_靶向_同位素", "说明": "有同位素支持，但仍不是确证结果"},
        {"项目": "最后：谨慎看", "结果": "05_靶向_仅质量", "说明": "只有精确质量证据，不能直接当作已鉴定化合物"},
        {"项目": "样本文件", "结果": task.get("filename", ""), "说明": ""},
        {"项目": "MS1 扫描数", "结果": expected_summary.get("MS1_Scan_Count", 0), "说明": "用于独立靶向 XIC"},
        {"项目": "MS2 谱图数", "结果": task.get("total_spectra", len(results)), "说明": ""},
        {"项目": "MS2 谱库匹配数", "结果": task.get("matched_spectra", 0), "说明": "谱图层计数"},
        {"项目": "MS2 匹配率（%）", "结果": round(task.get("match_rate", 0), 1), "说明": ""},
        {"项目": "MS2 去重结果数", "结果": len(matched), "说明": "每个候选保留最佳谱图"},
        {"项目": "High 级鉴定组数", "结果": high_scope, "说明": "High 单一二维连接与 High 未解异构体组之和"},
        {"项目": "High 单一二维连接数", "结果": high_exact, "说明": "当前证据可支持到一个二维连接结构；完整立体身份另行报告"},
        {"项目": "High 未解异构体组数", "结果": high_isomer, "说明": "High 证据支持该组，不能指定组内唯一异构体"},
        {"项目": "未核验空白关联警示数", "结果": unverified_blank, "说明": "同极性空白中亦有信号，但空白角色／方法配对未经核验；仅警示，不自动剔除"},
        {"项目": "靶向检索状态", "结果": expected_summary.get("Status_Label", ""), "说明": "独立连续 MS1 靶向层"},
        {"项目": "靶向数据库", "结果": _readable_dataset_text(str(expected_summary.get("Selected_Datasets") or "").split(" | ")), "说明": ""},
        {"项目": "数据库目标行数", "结果": expected_summary.get("Target_Catalog_Rows", 0), "说明": ""},
        {"项目": "靶向检出特征数", "结果": expected_summary.get("Compound_Summary_Rows", 0), "说明": "标准表每个化合物/未解歧义组一行"},
        {"项目": "靶向原始 XIC 峰段数", "结果": expected_summary.get("Detected_Features", 0), "说明": "完整记录见 10_靶向峰审计 或单独的靶向峰审计 CSV"},
        {"项目": "靶向 MS2 支持候选组", "结果": compound_evidence.get("EXPERIMENTAL_LIBRARY_CONFIRMED", 0) + compound_evidence.get("MS1_MS2_DIAGNOSTIC_PROBABLE", 0), "说明": "实验谱库候选或诊断离子筛查支持；不是标准品确证"},
        {"项目": "靶向同位素支持组", "结果": compound_evidence.get("MS1_ISOTOPE_SUPPORTED", 0), "说明": "支持性证据，仍非确证"},
        {"项目": "仅精确质量候选组", "结果": compound_evidence.get("MS1_EXACT_MASS_CANDIDATE", 0), "说明": "仍需进一步证据确认"},
        {"项目": "实验谱库候选支持组", "结果": compound_evidence.get("EXPERIMENTAL_LIBRARY_CONFIRMED", 0), "说明": "MS/MS 谱库候选支持；不是同方法标准品确证"},
        {"项目": "计数规则", "结果": "独立统计", "说明": "靶向 MS1 结果不叠加到普通 MS2 谱库匹配数"},
        {"项目": "MS1 模块适用性提示", "结果": module_warning_count, "说明": "非零时完整文字独立放在 11_模块提示，不在结果行重复铺开"},
        {"项目": "空白过滤排除谱图", "结果": len(blank_excluded), "说明": "非零时见附表"},
        {"项目": "空白清洁视图 MS2 结果", "结果": len(matched_clean), "说明": "已移除空白相关结果"},
        {"项目": "空白证据记录", "结果": len(blank_evidence), "说明": "非零时见附表"},
    ]


def _build_analysis_excel_bytes(
    task_id: str,
    task: dict,
    results: list[dict],
    *,
    language: str = "zh",
    output_path: Path | None = None,
    evidence_prevalidated: bool = False,
) -> bytes | Path:
    """Build the focused workbook used by both single and batch exports."""
    import resource
    profile_started = time.perf_counter()
    profile_cpu_started = time.process_time()

    def profile(phase: str, **fields):
        usage = resource.getrusage(resource.RUSAGE_SELF)
        _export_progress(
            phase=phase,
            export_profile={
                "phase": phase,
                "wall_seconds": time.perf_counter() - profile_started,
                "process_cpu_seconds": time.process_time() - profile_cpu_started,
                "max_rss_kib": int(usage.ru_maxrss // 1024
                                   if sys.platform == "darwin"
                                   else usage.ru_maxrss),
                **fields,
            },
        )

    from web.lossless_peak_export import (
        add_raw_peak_sheets, peak_sheet_names, prepare_reader_peak_cells,
    )

    all_matched = [result for result in results if result.get("match")]
    matched = sorted(
        _dedup_matched(all_matched),
        key=lambda result: float(
            (result.get("match") or {}).get("final_score") or 0
        ),
        reverse=True,
    )
    matched_clean = sorted(_dedup_matched([
        result for result in all_matched
        if not bool(result.get("blank_clean_view_excluded"))
    ]), key=lambda result: float(
        (result.get("match") or {}).get("final_score") or 0
    ), reverse=True)
    blank_excluded = [
        result for result in results if result.get("blank_excluded")
    ]
    blank_evidence = [
        result for result in results
        if result.get("blank_feature_matched")
        and (
            result.get("blank_flagged")
            or result.get("blank_clean_view_excluded")
            or result.get("blank_carryover_suspected")
        )
    ]
    unmatched = [
        result for result in results
        if not result.get("match") and not result.get("blank_excluded")
    ]

    from web.report_common import (
        ISOMER_CANDIDATE_EXPORT_COLUMNS,
        MS2_READER_COLUMN_NAMES,
        isomer_candidate_export_labels,
        isomer_candidate_export_rows,
        isomer_candidate_export_rows_from_compounds,
        ms2_reader_column_labels,
        ms2_reader_rows,
    )

    compound_conclusions = (
        _compound_conclusions_for_export(task_id, task, prevalidated=True)
        if evidence_prevalidated
        else _compound_conclusions_for_export(task_id, task)
    )
    match_rows = ms2_reader_rows(
        matched,
        compounds=compound_conclusions,
        language=language,
    )
    # A clean-view subset cannot safely reuse all-task final groups.  Rebuild
    # its reader rows from the filtered spectra so blank-related identities do
    # not leak back into the clean worksheet.
    clean_match_rows = ms2_reader_rows(
        matched_clean,
        compounds=None,
        language=language,
    )
    match_rows = _localized_ms2_rows(match_rows, language)
    clean_match_rows = _localized_ms2_rows(
        clean_match_rows, language
    )
    match_rows = prepare_reader_peak_cells(
        match_rows, task_id=task_id, results=results, language=language,
    )
    clean_match_rows = prepare_reader_peak_cells(
        clean_match_rows, task_id=task_id, results=results, language=language,
    )
    from web.plant_candidate_export import plant_reader_views
    plant_candidate_rows, external_review_rows = plant_reader_views(
        matched, compounds=compound_conclusions, language=language,
    )
    plant_candidate_rows = prepare_reader_peak_cells(
        _localized_ms2_rows(plant_candidate_rows, language), task_id=task_id,
        results=results, language=language,
    )
    external_review_rows = prepare_reader_peak_cells(
        _localized_ms2_rows(external_review_rows, language), task_id=task_id,
        results=results, language=language,
    )
    isomer_candidate_rows = (
        isomer_candidate_export_rows_from_compounds(
            compound_conclusions, language=language
        )
        if compound_conclusions is not None
        else isomer_candidate_export_rows(matched, language=language)
    )
    expected_payload = _load_expected_compounds_payload(
        task_id, task, tolerate_errors=True
    )
    reference_qc_payload = _load_reference_standard_qc_payload(
        task_id, tolerate_errors=True
    )
    flagship_path = _flagship_artifact_path(TASK_DIR, task_id)
    multilibrary_rows = (
        list(read_multilibrary_evidence(flagship_path).get("rows") or [])
        if flagship_path.is_file()
        else []
    )
    target_groups = expected_feature_groups(
        expected_payload, language=language
    )
    target_raw_rows = expected_raw_feature_rows(
        expected_payload, language=language
    )
    module_warning_rows = analysis_warning_rows(
        task.get("analysis_warnings") or []
    )
    target_summary = expected_summary_row(
        expected_payload, language=language
    )
    confidence_summary = _compound_confidence_summary(all_matched)
    compound_summary = dict(task.get("compound_conclusion_summary") or {})
    summary_rows = _analysis_summary_rows(
        task,
        results,
        matched,
        matched_clean,
        blank_excluded,
        blank_evidence,
        target_summary,
        confidence_summary,
        compound_summary,
        language=language,
    )
    summary_rows.extend(
        reference_standard_qc_summary_rows(
            reference_qc_payload, language=language
        )
    )
    reference_qc_rows = reference_standard_qc_detail_rows(
        reference_qc_payload, language=language
    )

    summary_columns = (
        ["Item", "Value", "Note"]
        if language == "en"
        else ["项目", "结果", "说明"]
    )
    sheet_names = (
        {
            "guide": "00_Confidence_Guide",
            "summary": "01_Read_Me",
            "target_high": "02_Target_High",
            "ms2": "03_MS2_Library",
            "isomer_candidates": "MS2_Isomer_Candidates",
            "target_isotope": "04_Target_Isotope",
            "target_mass": "05_Target_MassOnly",
            "unmatched": "06_Unmatched",
            "ms2_clean": "07_MS2_Clean",
            "blank_filtered": "08_Blank_Filtered",
            "blank_evidence": "09_Blank_Evidence",
            "target_audit": "10_Target_Peak_Audit",
            "module_warnings": "11_Module_Warnings",
            "reference_qc": "12_Reference_QC",
            "multilibrary": "13_MultiLibrary_Evidence",
        }
        if language == "en"
        else {
            "guide": "00_置信度说明",
            "summary": "01_结果导航",
            "target_high": "02_靶向_高证据",
            "ms2": "03_MS2谱库鉴定",
            "isomer_candidates": "MS2_异构体候选明细",
            "target_isotope": "04_靶向_同位素",
            "target_mass": "05_靶向_仅质量",
            "unmatched": "06_未匹配",
            "ms2_clean": "07_MS2清洁视图",
            "blank_filtered": "08_空白过滤",
            "blank_evidence": "09_空白证据",
            "target_audit": "10_靶向峰审计",
            "module_warnings": "11_模块提示",
            "reference_qc": "12_混标回收审计",
            "multilibrary": "13_多谱库证据",
        }
    )
    target_labels = expected_feature_column_labels(language)
    target_raw_labels = expected_raw_feature_column_labels(language)
    ms2_labels = ms2_reader_column_labels(language)
    unmatched_rows = [{
        "Spectrum": result.get("spectrum_id", ""),
        "Precursor_mz": result.get("precursor_mz"),
        "Retention_Time": result.get("retention_time"),
        "Peaks_Count": result.get("peaks_count", 0),
    } for result in unmatched]
    blank_rows = [{
        "Spectrum": result.get("spectrum_id", ""),
        "Precursor_mz": result.get("precursor_mz"),
        "Retention_Time_Min": result.get("retention_time"),
        "Sample_Intensity": result.get("blank_sample_intensity"),
        "Blank_Intensity": result.get("blank_intensity"),
        "Sample_Blank_Ratio": result.get("blank_sample_ratio"),
        "Evidence_Source": result.get("blank_evidence_source", ""),
        "MZ_Error_PPM": result.get("blank_mz_error_ppm"),
        "RT_Delta_Min": result.get("blank_rt_delta_minutes"),
        "Status": result.get("blank_status", ""),
    } for result in blank_excluded]
    blank_evidence_rows = [{
        "Spectrum": result.get("spectrum_id", ""),
        "Precursor_mz": result.get("precursor_mz"),
        "Retention_Time_Min": result.get("retention_time"),
        "Compound": (result.get("match") or {}).get("compound_name", ""),
        "Sample_Blank_Ratio": result.get("blank_sample_ratio"),
        "Evidence_Category": result.get("blank_evidence_category", ""),
        "Blank_Role": result.get("blank_evidence_role", ""),
        "Flagged": _yes_no(result.get("blank_flagged"), language),
        "Carryover_Suspected": _yes_no(
            result.get("blank_carryover_suspected"), language
        ),
        "Clean_View_Excluded": _yes_no(
            result.get("blank_clean_view_excluded"), language
        ),
        "Status": result.get("blank_status", ""),
    } for result in blank_evidence]

    def export_frame(
        rows: list[dict[str, Any]],
        columns: list[str],
        labels: list[str],
    ) -> pd.DataFrame:
        frame = pd.DataFrame(rows, columns=columns)
        frame.columns = labels
        return frame

    peak_index_name, peak_json_name = peak_sheet_names(language)
    if language == "en":
        summary_rows.append({"Item": "Plant candidate reading table",
                             "Value": "16_Plant_Candidates / 17_External_Review",
                             "Note": "03_MS2_Library is the unfiltered technical table. Reviewed external or conflicting candidates are excluded from sheet 16 and retained in 17. Unclassified origin is not proof of plant origin; original ranks and unresolved group identities remain unchanged."})
        summary_rows.append({"Item": "Complete saved raw peaks",
                             "Value": f"{peak_index_name} / {peak_json_name}",
                             "Note": "All saved spectra, including unmatched and blank rows. Join by zero-based result index; concatenate JSON chunks by chunk number and verify SHA256 of ASCII JSON. Missing or incomplete historical arrays remain marked. Also downloadable as Complete raw peaks JSONL.gz."})
    else:
        summary_rows.append({"项目": "植物候选阅读与外源复核",
                             "结果": "16_植物候选 / 17_外源与冲突复核",
                             "说明": "03_MS2谱库鉴定是未按来源筛选的技术全表。已核验外源或身份冲突候选排除于16表并保留在17表；来源未分类不等于植物内源。保留原候选排名及未解异构体限定。"})
        summary_rows.append({"项目": "完整已保存原始峰",
                             "结果": f"{peak_index_name} / {peak_json_name}",
                             "说明": "包含未匹配和空白排除谱。按结果数组索引关联，按片段序号拼接 JSON，核对 ASCII JSON 的 SHA256；历史未保存或不完整数组仍标注缺失。也可下载完整原始峰 JSONL.gz。"})
    peak_source_path = TASK_DIR / f"{task_id}_results.json"
    peak_source_sha256 = None
    if peak_source_path.is_file():
        import hashlib
        with peak_source_path.open("rb") as source:
            peak_source_sha256 = hashlib.file_digest(source, "sha256").hexdigest()

    profile("projecting_results", result_rows=len(results),
            projected_ms2_rows=len(match_rows))
    excel_target = Path(output_path) if output_path is not None else io.BytesIO()
    profile("rendering_worksheets")
    with pd.ExcelWriter(excel_target, engine="openpyxl") as writer:
        pd.DataFrame(confidence_guide_table_rows(language)).to_excel(
            writer, sheet_name=sheet_names["guide"], index=False
        )
        pd.DataFrame(
            summary_rows, columns=summary_columns
        ).to_excel(writer, sheet_name=sheet_names["summary"], index=False)
        export_frame(
            target_groups["high"],
            EXPECTED_FEATURE_COLUMNS,
            target_labels,
        ).to_excel(
            writer, sheet_name=sheet_names["target_high"], index=False
        )
        export_frame(
            match_rows, MS2_READER_COLUMN_NAMES, ms2_labels
        ).to_excel(writer, sheet_name=sheet_names["ms2"], index=False)
        export_frame(
            isomer_candidate_rows,
            ISOMER_CANDIDATE_EXPORT_COLUMNS,
            isomer_candidate_export_labels(language),
        ).to_excel(
            writer,
            sheet_name=sheet_names["isomer_candidates"],
            index=False,
        )
        export_frame(
            target_groups["isotope"],
            EXPECTED_FEATURE_COLUMNS,
            target_labels,
        ).to_excel(
            writer, sheet_name=sheet_names["target_isotope"], index=False
        )
        export_frame(
            target_groups["mass_only"],
            EXPECTED_FEATURE_COLUMNS,
            target_labels,
        ).to_excel(
            writer, sheet_name=sheet_names["target_mass"], index=False
        )
        unmatched_columns = [
            "Spectrum", "Precursor_mz", "Retention_Time", "Peaks_Count"
        ]
        unmatched_labels = (
            unmatched_columns
            if language == "en"
            else ["谱图编号", "前体 m/z", "RT (min)", "峰数"]
        )
        export_frame(
            unmatched_rows, unmatched_columns, unmatched_labels
        ).to_excel(
            writer, sheet_name=sheet_names["unmatched"], index=False
        )
        if any(
            bool(result.get("blank_clean_view_excluded"))
            for result in all_matched
        ):
            export_frame(
                clean_match_rows, MS2_READER_COLUMN_NAMES, ms2_labels
            ).to_excel(
                writer, sheet_name=sheet_names["ms2_clean"], index=False
            )
        if blank_rows:
            pd.DataFrame(blank_rows).to_excel(
                writer, sheet_name=sheet_names["blank_filtered"], index=False
            )
        if blank_evidence_rows:
            pd.DataFrame(blank_evidence_rows).to_excel(
                writer, sheet_name=sheet_names["blank_evidence"], index=False
            )
        export_frame(
            target_raw_rows,
            EXPECTED_RAW_FEATURE_COLUMNS,
            target_raw_labels,
        ).to_excel(
            writer, sheet_name=sheet_names["target_audit"], index=False
        )
        if module_warning_rows:
            export_frame(
                module_warning_rows,
                ANALYSIS_WARNING_COLUMNS,
                ANALYSIS_WARNING_LABELS[language],
            ).to_excel(
                writer,
                sheet_name=sheet_names["module_warnings"],
                index=False,
            )
        if reference_qc_rows:
            pd.DataFrame(reference_qc_rows).to_excel(
                writer,
                sheet_name=sheet_names["reference_qc"],
                index=False,
            )
        export_frame(
            multilibrary_rows,
            MULTILIBRARY_COLUMNS,
            multilibrary_column_labels(language),
        ).to_excel(
            writer,
            sheet_name=sheet_names["multilibrary"],
            index=False,
        )
        from web.governance_blocked_projection import (
            GOVERNANCE_BLOCKED_EXPORT_COLUMNS, governance_blocked_export_rows,
            project_governance_blocked_candidates,
        )
        pd.DataFrame(
            governance_blocked_export_rows(
                project_governance_blocked_candidates(results)
            ),
            columns=GOVERNANCE_BLOCKED_EXPORT_COLUMNS,
        ).to_excel(
            writer, sheet_name=("Governance_Blocked" if language == "en" else "治理未放行候选"),
            index=False,
        )
        profile("styling_workbook", worksheet_count=len(writer.book.worksheets))
        style_analysis_workbook(writer.book)
        format_confidence_guide_worksheet(
            writer.book[sheet_names["guide"]], language
        )
        profile("rendering_raw_peaks")
        add_raw_peak_sheets(
            writer.book, results, task_id, language=language,
            source_results_sha256=peak_source_sha256,
        )
        for view_rows, view_name in (
            (plant_candidate_rows, "16_Plant_Candidates" if language == "en" else "16_植物候选"),
            (external_review_rows, "17_External_Review" if language == "en" else "17_外源与冲突复核"),
        ):
            export_frame(view_rows, MS2_READER_COLUMN_NAMES, ms2_labels).to_excel(
                writer, sheet_name=view_name, index=False,
            )
            worksheet = writer.book[view_name]
            worksheet.freeze_panes = "C2"
            worksheet.auto_filter.ref = worksheet.dimensions
            # Preserve the reader-table formatting, including column widths.
            import copy
            source_sheet = writer.book[sheet_names["ms2"]]
            for key, dimension in source_sheet.column_dimensions.items():
                worksheet.column_dimensions[key].width = dimension.width
            for target, original in zip(worksheet[1], source_sheet[1]):
                target._style = copy.copy(original._style)
            worksheet.row_dimensions[1].height = source_sheet.row_dimensions[1].height

        profile("finalizing_xlsx", worksheet_count=len(writer.book.worksheets))

    if output_path is not None:
        with Path(output_path).open("rb") as handle:
            os.fsync(handle.fileno())
        profile("xlsx_ready", size_bytes=Path(output_path).stat().st_size)
        return Path(output_path)
    excel_target.seek(0)
    payload = excel_target.getvalue()
    profile("xlsx_ready", size_bytes=len(payload))
    return payload


_NATIVE_MSMS_EVENT_COLUMNS = (
    "result_index", "spectrum_id", "precursor_mz", "search_precursor_mz",
    "search_precursor_mz_source", "search_precursor_delta_ppm",
    "retention_time", "retention_time_unit", "ion_mode", "precursor_charge",
    "collision_energy", "collision_energy_unit", "collision_energy_mode",
    "activation_method", "instrument_configuration_id", "lcms_feature_id",
    "precursor_feature_id", "lcms_feature_mz", "lcms_feature_rt_minutes",
    "lcms_feature_used_for_search", "blank_status", "blank_flagged",
    "blank_excluded", "blank_clean_view_excluded", "blank_sample_ratio",
    "blank_sample_intensity", "blank_intensity", "match_present",
    "top1_compound_name", "top1_confidence", "top1_evidence_tier",
    "top1_annotation_level", "top1_final_score", "top1_cosine",
    "top1_entropy", "top1_mass_ppm", "top1_formula", "top1_inchikey",
    "top1_smiles", "top1_adduct", "top1_dataset_source",
    "top1_library_spectrum_id", "native_record_json",
)


def _native_msms_blank_status(row: Mapping[str, Any]) -> str:
    saved = str(row.get("blank_status") or "").strip()
    if saved:
        return saved
    if bool(row.get("blank_excluded")):
        return "excluded"
    if bool(row.get("blank_flagged")):
        return "flagged"
    if bool(row.get("blank_clean_view_excluded")):
        return "clean_view_excluded"
    return "not_flagged"


def _native_msms_event_row(row: Mapping[str, Any]) -> dict[str, Any]:
    """Project one saved MS/MS result without replacing its scan-level evidence."""
    match = row.get("match") if isinstance(row.get("match"), Mapping) else {}
    saved_row = {key: value for key, value in row.items() if key != "_result_index"}
    return {
        "result_index": row.get("_result_index"),
        "spectrum_id": row.get("spectrum_id"),
        "precursor_mz": row.get("precursor_mz"),
        "search_precursor_mz": row.get("search_precursor_mz"),
        "search_precursor_mz_source": row.get("search_precursor_mz_source"),
        "search_precursor_delta_ppm": row.get("search_precursor_delta_ppm"),
        "retention_time": row.get("retention_time"),
        "retention_time_unit": row.get("retention_time_unit"),
        "ion_mode": row.get("ion_mode"),
        "precursor_charge": row.get("precursor_charge"),
        "collision_energy": row.get("collision_energy"),
        "collision_energy_unit": row.get("collision_energy_unit"),
        "collision_energy_mode": row.get("collision_energy_mode"),
        "activation_method": row.get("activation_method"),
        "instrument_configuration_id": row.get("instrument_configuration_id"),
        "lcms_feature_id": row.get("lcms_feature_id"),
        "precursor_feature_id": row.get("precursor_feature_id"),
        "lcms_feature_mz": row.get("lcms_feature_mz"),
        "lcms_feature_rt_minutes": row.get("lcms_feature_rt_minutes"),
        "lcms_feature_used_for_search": row.get("lcms_feature_used_for_search"),
        "blank_status": _native_msms_blank_status(row),
        "blank_flagged": row.get("blank_flagged"),
        "blank_excluded": row.get("blank_excluded"),
        "blank_clean_view_excluded": row.get("blank_clean_view_excluded"),
        "blank_sample_ratio": row.get("blank_sample_ratio"),
        "blank_sample_intensity": row.get("blank_sample_intensity"),
        "blank_intensity": row.get("blank_intensity"),
        "match_present": bool(match),
        "top1_compound_name": match.get("compound_name"),
        "top1_confidence": match.get("confidence"),
        "top1_evidence_tier": match.get("evidence_tier"),
        "top1_annotation_level": match.get("annotation_level"),
        "top1_final_score": match.get("final_score"),
        "top1_cosine": match.get("cosine"),
        "top1_entropy": match.get("entropy"),
        "top1_mass_ppm": match.get("mass_ppm"),
        "top1_formula": match.get("formula"),
        "top1_inchikey": match.get("inchikey"),
        "top1_smiles": match.get("smiles"),
        "top1_adduct": match.get("adduct"),
        "top1_dataset_source": match.get("dataset_source"),
        "top1_library_spectrum_id": match.get("lib_spectrum_id"),
        "native_record_json": json.dumps(
            saved_row, ensure_ascii=False, sort_keys=True,
            separators=(",", ":"), allow_nan=True,
        ),
    }


def _native_msms_event_csv_chunks(rows: list[dict[str, Any]]):
    """Stream one reversible CSV record per saved MS/MS result row."""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=_NATIVE_MSMS_EVENT_COLUMNS)
    writer.writeheader()
    yield "\ufeff" + buffer.getvalue()
    for row in rows:
        buffer.seek(0)
        buffer.truncate(0)
        writer.writerow(_native_msms_event_row(row))
        yield buffer.getvalue()


@app.get("/api/tasks/{task_id}/download/csv")
@_serialized_report_export
def download_csv(
    task_id: str,
    format_type: str = Query("results", enum=["results", "report"]),
    result_layer: Literal[
        "ms2", "isomer_candidates", "targeted", "targeted_features",
        "reference_qc", "multilibrary", "plant_candidates", "external_review"
    ] = Query("ms2"),
    ms2_schema: Literal["reader_v3", "complete_v2"] = Query("reader_v3"),
    blank_view: Literal["all", "clean", "flagged"] = Query("all"),
):
    """Download results as CSV."""
    from web.lossless_peak_export import attachment_disposition
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    # Complete events come from the immutable saved artifact, before the
    # display reader attaches final-claim annotations to match dictionaries.
    native_source = TASK_DIR / f"{task_id}_results.json"
    source_rows = (
        _iter_legacy_result_artifact(native_source)
        if ms2_schema == "complete_v2" and native_source.is_file()
        else _get_task_results(task_id)
    )
    results = [
        {**row, "_result_index": index}
        for index, row in enumerate(source_rows)
    ]
    all_matched = [r for r in results if r.get("match")]
    selected_results = results
    selected_matched = all_matched
    if blank_view == "clean":
        selected_results = [
            r for r in results
            if not bool(r.get("blank_clean_view_excluded"))
        ]
        selected_matched = [
            r for r in all_matched
            if not bool(r.get("blank_clean_view_excluded"))
        ]
    elif blank_view == "flagged":
        selected_results = [r for r in results if bool(r.get("blank_flagged"))]
        selected_matched = [
            r for r in all_matched if bool(r.get("blank_flagged"))
        ]
    from fastapi.responses import StreamingResponse

    if ms2_schema == "complete_v2":
        if format_type != "results" or result_layer != "ms2":
            raise HTTPException(
                status_code=400,
                detail=(
                    "complete_v2 is available only for format_type=results "
                    "and result_layer=ms2"
                ),
            )
        filename = (
            f"{Path(task['filename']).stem}_raw_msms_events_{blank_view}.csv"
        )
        return StreamingResponse(
            _native_msms_event_csv_chunks(selected_results),
            media_type="text/csv; charset=utf-8",
            headers={
                "Content-Disposition": attachment_disposition(filename),
                "Cache-Control": "no-store, max-age=0",
                "X-MSsystem-Result-Layer": "ms2",
                "X-MSsystem-MS2-Schema": "complete_v2",
                "X-MSsystem-Record-Unit": "raw-msms-event",
            },
        )

    matched = _dedup_matched(selected_matched)
    expected_payload = _load_expected_compounds_payload(
        task_id, task, tolerate_errors=True
    )
    reference_qc_payload = _load_reference_standard_qc_payload(
        task_id, tolerate_errors=True
    )

    report_mode = (format_type == "report")
    if result_layer in {"plant_candidates", "external_review"}:
        compound_conclusions = (
            _compound_conclusions_for_export(task_id, task)
            if blank_view == "all" else None
        )
        csv_text = _ms2_csv_text(
            matched, language="zh", compounds=compound_conclusions,
            task_id=task_id, results=results, candidate_view=result_layer,
        )
        filename = f"{Path(task['filename']).stem}_{result_layer}.csv"
        return StreamingResponse(
            iter(["\ufeff" + csv_text]), media_type="text/csv; charset=utf-8",
            headers={"Content-Disposition": attachment_disposition(filename),
                     "Cache-Control": "no-store, max-age=0",
                     "X-MSsystem-Result-Layer": result_layer},
        )
    if result_layer == "multilibrary":
        artifact = _flagship_artifact_path(TASK_DIR, task_id)
        if not artifact.is_file():
            raise HTTPException(
                status_code=404,
                detail=(
                    "V15 flagship artifact is unavailable; open the V15 flagship "
                    "result once to materialize it"
                ),
            )
        payload = read_multilibrary_evidence(artifact)
        csv_text = _focused_csv_text(
            list(payload.get("rows") or []),
            MULTILIBRARY_COLUMNS,
            multilibrary_column_labels("zh"),
        )
        filename = f"{Path(task['filename']).stem}_v15_multilibrary_evidence.csv"
        return StreamingResponse(
            iter(["\ufeff" + csv_text]),
            media_type="text/csv; charset=utf-8",
            headers={
                "Content-Disposition": attachment_disposition(filename),
                "Cache-Control": "no-store, max-age=0",
                "X-MSsystem-Result-Layer": "multilibrary",
            },
        )
    if result_layer == "targeted_features":
        csv_text = _targeted_feature_audit_csv_text(
            expected_payload, language="zh"
        )
        filename = (
            f"{Path(task['filename']).stem}_targeted_peak_audit.csv"
        )
        return StreamingResponse(
            iter(["\ufeff" + csv_text]),
            media_type="text/csv; charset=utf-8",
            headers={
                "Content-Disposition": attachment_disposition(filename),
                "Cache-Control": "no-store, max-age=0",
                "X-MSsystem-Result-Layer": "targeted_features",
            },
        )

    if result_layer == "reference_qc":
        csv_text = reference_standard_qc_csv_text(
            reference_qc_payload, language="zh"
        )
        if not csv_text:
            raise HTTPException(
                status_code=404,
                detail="No explicit reference-standard QC set is attached to this task",
            )
        filename = f"{Path(task['filename']).stem}_reference_standard_qc.csv"
        return StreamingResponse(
            iter(["\ufeff" + csv_text]),
            media_type="text/csv; charset=utf-8",
            headers={
                "Content-Disposition": attachment_disposition(filename),
                "Cache-Control": "no-store, max-age=0",
                "X-MSsystem-Result-Layer": "reference_qc",
            },
        )

    if result_layer == "isomer_candidates":
        compound_conclusions = (
            _compound_conclusions_for_export(task_id, task)
            if blank_view == "all"
            else None
        )
        csv_text = _isomer_candidate_csv_text(
            matched,
            language="zh",
            compounds=compound_conclusions,
        )
        filename = (
            f"{Path(task['filename']).stem}_isomer_competing_candidates.csv"
        )
        return StreamingResponse(
            iter(["\ufeff" + csv_text]),
            media_type="text/csv; charset=utf-8",
            headers={
                "Content-Disposition": attachment_disposition(filename),
                "Cache-Control": "no-store, max-age=0",
                "X-MSsystem-Result-Layer": "isomer_candidates",
            },
        )

    if not report_mode:
        if result_layer == "targeted":
            csv_text = _targeted_csv_text(
                expected_payload, language="zh"
            )
            filename = (
                f"{Path(task['filename']).stem}_targeted_compounds.csv"
            )
        else:
            compound_conclusions = (
                _compound_conclusions_for_export(task_id, task)
                if blank_view == "all"
                else None
            )
            csv_text = _ms2_csv_text(
                matched,
                language="zh",
                compounds=compound_conclusions,
                task_id=task_id, results=selected_results,
            )
            filename = (
                f"{Path(task['filename']).stem}_ms2_{blank_view}.csv"
            )
        return StreamingResponse(
            iter(["\ufeff" + csv_text]),
            media_type="text/csv; charset=utf-8",
            headers={
                "Content-Disposition": attachment_disposition(filename),
                "Cache-Control": "no-store, max-age=0",
                "X-MSsystem-Result-Layer": result_layer,
                "X-MSsystem-MS2-Schema": (
                    "reader_v3" if result_layer == "ms2" else ""
                ),
            },
        )

    # report 态是明确选择的完整审计表；标准结果已在上方返回单一扁平表。
    filename = (
        f"{Path(task['filename']).stem}_{format_type}_{blank_view}.csv"
    )
    from web.report_common import REPORT_COLUMN_NAMES, row_to_dict

    cols = REPORT_COLUMN_NAMES

    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(cols)
    for r in matched:
        d = row_to_dict(r, report_mode=True)
        w.writerow([d.get(c, "") for c in cols])

    # 特征级空白过滤审计。被排除的谱图不会进入库检索，因此必须单列导出，
    # 不能混入“未命中”或从 matched 表里悄然消失。
    _blank_excluded = [r for r in results if r.get("blank_excluded")]
    if _blank_excluded:
        w.writerow([])
        w.writerow([f"=== 空白特征过滤排除（共 {len(_blank_excluded)} 条） ==="])
        w.writerow([
            "Spectrum", "Precursor_mz", "Retention_Time_Min", "Sample_Intensity",
            "Blank_Intensity", "Sample_Blank_Ratio", "Evidence_Source",
            "MZ_Error_PPM", "RT_Delta_Min", "Quantitation_Mode",
            "Sample_XIC_Area", "Blank_XIC_Area", "Blank_Detected",
            "Blank_Applicable", "Blank_Detection_Fraction", "Status",
        ])
        for r in _blank_excluded:
            w.writerow([
                r.get("spectrum_id", ""), _fmt_csv_number(r.get("precursor_mz"), 4),
                _fmt_csv_number(r.get("retention_time"), 4),
                _fmt_csv_number(r.get("blank_sample_intensity"), 4),
                _fmt_csv_number(r.get("blank_intensity"), 4),
                _fmt_csv_number(r.get("blank_sample_ratio"), 4),
                r.get("blank_evidence_source", ""),
                _fmt_csv_number(r.get("blank_mz_error_ppm"), 4),
                _fmt_csv_number(r.get("blank_rt_delta_minutes"), 4),
                r.get("blank_quantitation_mode", ""),
                _fmt_csv_number(r.get("blank_sample_xic_area"), 4),
                _fmt_csv_number(r.get("blank_xic_area"), 4),
                int(r.get("blank_replicates_detected") or 0),
                int(r.get("blank_replicates_applicable") or 0),
                _fmt_csv_number(r.get("blank_detection_fraction"), 4),
                r.get("blank_status", ""),
            ])

    _blank_flagged = [
        r for r in results
        if r.get("blank_flagged") and not r.get("blank_excluded")
    ]
    if _blank_flagged:
        w.writerow([])
        w.writerow([
            f"=== 空白证据标记（仍保留鉴定；共 {len(_blank_flagged)} 条） ==="
        ])
        w.writerow([
            "Spectrum", "Precursor_mz", "Retention_Time_Min", "Compound",
            "Sample_Intensity", "Blank_Intensity", "Sample_Blank_Ratio",
            "Evidence_Category", "Blank_Role", "Carryover_Suspected",
            "Clean_View_Excluded", "Detection_Fraction", "Status",
        ])
        for r in _blank_flagged:
            w.writerow([
                r.get("spectrum_id", ""),
                _fmt_csv_number(r.get("precursor_mz"), 4),
                _fmt_csv_number(r.get("retention_time"), 4),
                (r.get("match") or {}).get("compound_name", ""),
                _fmt_csv_number(r.get("blank_sample_intensity"), 4),
                _fmt_csv_number(r.get("blank_intensity"), 4),
                _fmt_csv_number(r.get("blank_sample_ratio"), 4),
                r.get("blank_evidence_category", ""),
                r.get("blank_evidence_role", ""),
                bool(r.get("blank_carryover_suspected")),
                bool(r.get("blank_clean_view_excluded")),
                _fmt_csv_number(r.get("blank_detection_fraction"), 4),
                r.get("blank_status", ""),
            ])


    # 未命中谱图的 MSClassify 谱图判类 (display_class) — 让 CSV 导出也包含未命中谱图的分类
    # (此前 CSV 只写 matched, 未命中谱图完全不出现 → 用户在 CSV 里看不到这批谱图的任何分类)
    _um_classified = [r for r in results if not r.get("match") and r.get("display_class")]
    if _um_classified:
        w.writerow([])
        w.writerow([f"=== 未命中谱图判类 (MSClassify 从 MS2 谱图推断, 非库匹配; 共 {len(_um_classified)} 条; 阈值>=0.15) ==="])
        w.writerow(["Spectrum", "Precursor_mz", "Retention_Time",
                    "Predicted_Class", "DL_Confidence", "Confidence_Level"])
        for r in _um_classified:
            _conf = r.get("dl_confidence")
            _cf = float(_conf) if _conf is not None else 0.0
            w.writerow([
                r.get("spectrum_id", ""), f"{r['precursor_mz']:.4f}",
                r.get("retention_time") if r.get("retention_time") is not None else "",
                r.get("display_class", ""), f"{_cf:.3f}",
                ("强信号(>=0.70)" if _cf >= 0.70 else
                 "中等信号(0.45~0.70)" if _cf >= 0.45 else
                 "低信号(0.15~0.45)"),
            ])

    append_expected_csv_section(w, expected_payload, language="zh")
    append_reference_standard_qc_csv_section(
        w, reference_qc_payload, language="zh"
    )

    from web.governance_blocked_projection import (
        append_governance_blocked_csv_section, project_governance_blocked_candidates,
    )
    csv_text = append_governance_blocked_csv_section(
        buf.getvalue(), project_governance_blocked_candidates(selected_results),
        language="zh",
    )
    return StreamingResponse(
        iter(["\ufeff" + csv_text]),
        media_type="text/csv; charset=utf-8",
        headers={
            "Content-Disposition": attachment_disposition(filename),
            "Cache-Control": "no-store, max-age=0",
            "X-MSsystem-Result-Layer": "ms2_report",
        },
    )


@app.get("/api/tasks/{task_id}/download/raw-peaks")
@_serialized_report_export
def download_raw_peaks(task_id: str):
    """Download every saved raw array, including unmatched and blank spectra."""
    import hashlib
    from urllib.parse import quote
    from web.lossless_peak_export import iter_raw_peak_jsonl_gzip
    from fastapi.responses import StreamingResponse

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.get("status") != "completed":
        raise HTTPException(status_code=409, detail="Raw peak export requires a completed task")
    source_path = TASK_DIR / f"{task_id}_results.json"
    source_sha256 = None
    if source_path.is_file():
        with source_path.open("rb") as source:
            source_sha256 = hashlib.file_digest(source, "sha256").hexdigest()
        rows = _iter_legacy_result_artifact(source_path)
    elif isinstance(task.get("results"), list):
        rows = task["results"]
    else:
        raise HTTPException(status_code=404, detail="Saved result arrays are unavailable")
    filename = f"{Path(task['filename']).stem}_raw_peaks.jsonl.gz"
    return StreamingResponse(
        iter_raw_peak_jsonl_gzip(rows, task_id, source_sha256),
        media_type="application/gzip",
        headers={
            "Content-Disposition": f"attachment; filename=\"raw_peaks.jsonl.gz\"; filename*=UTF-8''{quote(filename)}",
            "Cache-Control": "no-store, max-age=0",
            "X-MSsystem-Result-Layer": "complete_saved_raw_peaks",
            "X-MSsystem-Source-Results-SHA256": source_sha256 or "in-memory-no-file-hash",
        },
    )


@app.get("/api/tasks/{task_id}/download/excel")
@_serialized_report_export
def download_excel(task_id: str):
    return _download_excel_impl(task_id)


def _download_excel_impl(task_id: str, *, evidence_prevalidated: bool = False):
    """Download results as Excel (.xlsx)."""
    from web.lossless_peak_export import attachment_disposition
    from fastapi.responses import StreamingResponse
    from starlette.background import BackgroundTask

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    # Refresh stale compound projections before the isolated workbook worker
    # reads the persisted task index.
    if not evidence_prevalidated:
        _compound_scope_summary_for_export(task_id, task)
    try:
        from web.export_jobs import has_progress_callback
        if has_progress_callback():
            # Background exports already run off the HTTP event loop. Reuse
            # this process's fully verified library instead of verifying the
            # entire reference cache in a fresh process for every workbook.
            import resource
            read_started = time.perf_counter()
            read_cpu_started = time.process_time()
            results = _get_task_results(task_id)
            read_usage = resource.getrusage(resource.RUSAGE_SELF)
            _export_progress(
                phase="reading_results",
                read_wall_seconds=time.perf_counter() - read_started,
                export_profile={
                    "phase": "reading_results",
                    "wall_seconds": time.perf_counter() - read_started,
                    "process_cpu_seconds": time.process_time() - read_cpu_started,
                    "max_rss_kib": int(read_usage.ru_maxrss // 1024
                                       if sys.platform == "darwin"
                                       else read_usage.ru_maxrss),
                    "result_rows": len(results),
                },
            )
            export_tmp_root = OUTPUT_DIR / "export_tmp"
            export_tmp_root.mkdir(parents=True, exist_ok=True)
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f"excel-{task_id[:16]}-", suffix=".xlsx",
                dir=export_tmp_root,
            )
            os.close(descriptor)
            temporary_path = Path(temporary_name)
            try:
                _build_analysis_excel_bytes(
                    task_id, task, results, language="zh",
                    output_path=temporary_path,
                    evidence_prevalidated=True,
                )
            except Exception:
                temporary_path.unlink(missing_ok=True)
                raise
            filename = f"{Path(task['filename']).stem}_analysis.xlsx"
            return FileResponse(
                temporary_path,
                media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                headers={"Content-Disposition": attachment_disposition(filename)},
                background=BackgroundTask(temporary_path.unlink, missing_ok=True),
            )
        else:
            excel_bytes = _build_excel_export_subprocess(task_id)
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
        logger.exception("Excel export worker failed for task %s", task_id)
        raise HTTPException(
            status_code=500,
            detail=f"Excel generation failed: {exc}",
        ) from exc
    filename = f"{Path(task['filename']).stem}_analysis.xlsx"
    return StreamingResponse(
        iter([excel_bytes]),
        media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={"Content-Disposition": attachment_disposition(filename)},
    )


def _compound_scope_summary_for_export(
    task_id: str, task: dict[str, Any]
) -> dict[str, Any]:
    """Validate primary and compound artifacts before returning report counters.

    Persisted tasks never fall back to stale counters after SQL validation or
    rebuilding fails. Legacy in-memory tasks without a saved source retain
    their existing test/report compatibility.
    """
    compound_summary = dict(task.get("compound_conclusion_summary") or {})
    compound_source = TASK_DIR / f"{task_id}_results.json"
    compound_artifact = _compound_conclusion_artifact_path(TASK_DIR, task_id)
    primary_artifact = _flagship_artifact_path(TASK_DIR, task_id)
    if compound_source.is_file():
        try:
            feature_inventory, feature_inventory_source = (
                _task_feature_inventory(task)
            )
            refreshed_flagship_summary = None
            if not _flagship_artifact_is_current(primary_artifact):
                refreshed_flagship_summary = _materialize_current_v15_flagship(
                    task_id,
                    task,
                    compound_source,
                    primary_artifact,
                )
                if not _flagship_artifact_is_current(primary_artifact):
                    raise RuntimeError("Primary identity artifact remains stale after rebuilding")
            if not _compound_conclusion_artifact_is_current(
                compound_artifact,
                compound_source,
                primary_artifact,
                feature_inventory_source,
            ):
                prepared_rows = _get_task_results(task_id)
                compound_summary = _write_compound_conclusion_artifact(
                    task_id,
                    prepared_rows,
                    compound_artifact,
                    claim_decisions_path=_claim_contract_path(
                        TASK_DIR, task_id
                    ),
                    source_artifact=compound_source,
                    primary_identity_artifact=primary_artifact,
                    feature_inventory=feature_inventory,
                    feature_inventory_source=feature_inventory_source,
                )
                if not _compound_conclusion_artifact_is_current(
                    compound_artifact,
                    compound_source,
                    primary_artifact,
                    feature_inventory_source,
                ):
                    raise RuntimeError("Compound conclusion artifact remains stale after rebuilding")
            else:
                current_payload = _read_compound_conclusion_artifact(
                    compound_artifact, offset=0, limit=1
                )
                if not isinstance(current_payload.get("metadata"), dict):
                    raise RuntimeError("Current compound artifact metadata is missing")
                compound_summary = dict(current_payload["metadata"])
            if not compound_artifact.is_file():
                raise RuntimeError("Current compound conclusion artifact is missing")
            if refreshed_flagship_summary is not None:
                task["flagship_summary"] = refreshed_flagship_summary
            task["compound_conclusion_summary"] = compound_summary
            _apply_final_compound_count_summary(task, compound_summary)
            with _tasks_lock:
                current = _tasks.get(task_id)
                if current is not None:
                    current["compound_conclusion_summary"] = compound_summary
                    _apply_final_compound_count_summary(
                        current, compound_summary
                    )
            _save_tasks()
        except (OSError, ValueError, TypeError, RuntimeError, json.JSONDecodeError) as exc:
            logger.exception(
                "Export refused: current primary/compound evidence unavailable for %s: %s",
                task_id,
                exc,
            )
            raise HTTPException(
                status_code=500,
                detail=f"Export refused: current primary/compound evidence could not be validated ({exc})",
            ) from exc
    return compound_summary


def _compound_conclusions_for_export(
    task_id: str, task: dict[str, Any], *, prevalidated: bool = False
) -> list[dict[str, Any]] | None:
    """Return the current final compound rows, or ``None`` for legacy fallback.

    ``None`` is reserved for tasks without a persisted result source. A saved
    source requires both primary and compound artifacts to pass their current
    dependency checks, including the active SQL metadata repair binding.
    """
    source_artifact = TASK_DIR / f"{task_id}_results.json"
    if not source_artifact.is_file():
        return None
    artifact = _compound_conclusion_artifact_path(TASK_DIR, task_id)
    # Compound freshness alone cannot detect a stale SQL binding inside an
    # otherwise unchanged primary artifact. Validate upstream first on every
    # persisted export, before opening any cached compound rows.
    if not prevalidated:
        _compound_scope_summary_for_export(task_id, task)
    try:
        if not artifact.is_file():
            raise RuntimeError("Current compound conclusion artifact is missing")
        payload = _read_compound_conclusion_artifact(
            artifact, offset=0, limit=5000
        )
        compounds = list(payload.get("compounds") or [])
        while payload.get("has_more"):
            payload = _read_compound_conclusion_artifact(
                artifact,
                offset=len(compounds),
                limit=5000,
            )
            page = list(payload.get("compounds") or [])
            if not page:
                raise RuntimeError("Compound artifact pagination ended before all rows were read")
            compounds.extend(page)
        return compounds
    except (OSError, ValueError, TypeError, RuntimeError, json.JSONDecodeError) as exc:
        logger.exception(
            "Export refused: current compound rows unavailable for %s: %s",
            task_id,
            exc,
        )
        raise HTTPException(
            status_code=500,
            detail=f"Export refused: current compound rows could not be loaded ({exc})",
        ) from exc


@app.get("/api/tasks/{task_id}/download/html")
@_serialized_report_export
def download_html_report(
    task_id: str,
    format_type: str = Query(
        "standard", enum=["standard", "detailed", "targeted"]
    ),
    language: str = Query("zh", enum=["zh", "en"]),
):
    """Download results as HTML report."""
    from fastapi.responses import HTMLResponse

    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    results = _get_task_results(task_id)
    from web.governance_blocked_projection import (
        inject_governance_blocked_html, project_governance_blocked_candidates,
    )
    governance_blocked_candidates = project_governance_blocked_candidates(results)
    # Detailed/English exports must use the same current compound-scope
    # counters as the interactive result page.  Older tasks may carry a stale
    # sidecar whose summary counted only exact-structure High rows.
    compound_summary = _compound_scope_summary_for_export(task_id, task)
    # Sidecars such as EvoGraph-F are keyed to the immutable scan-row index in
    # ``*_results.json``.  Presentation deduplication must retain that key;
    # otherwise the right graph exists on disk but cannot be joined back to
    # its representative compound in a detailed report.
    for result_index, row in enumerate(results):
        row["_result_index"] = result_index
    expected_payload = _load_expected_compounds_payload(
        task_id, task, tolerate_errors=True
    )
    reference_qc_payload = _load_reference_standard_qc_payload(
        task_id, tolerate_errors=True
    )
    if format_type == "targeted":
        html_content = render_expected_html_report(
            expected_payload,
            task.get("filename", "Analysis"),
            language=language,
            analysis_warnings=task.get("analysis_warnings") or [],
        )
        html_content = inject_reference_standard_qc_html(
            html_content, reference_qc_payload, language=language
        )
        html_content = inject_confidence_guide_html(
            html_content, language=language
        )
        suffix = "_en" if language == "en" else ""
        filename = (
            f"{Path(task['filename']).stem}_targeted_library_report"
            f"{suffix}.html"
        )
        return HTMLResponse(
            content=html_content,
            headers={
                "Content-Disposition": (
                    attachment_disposition(filename)
                ),
                "Cache-Control": "no-store, max-age=0",
                "X-MSsystem-Report-Layer": "targeted",
            },
        )

    # 报告同时保留两个不同口径：原始 MS2 扫描统计来自未去重结果，
    # 主表则按化合物身份去重。旧实现先去重再统计，曾把 2494 张原始
    # MS2 误报为 2265，且把 356 张命中谱误报为 127 个候选。
    total_spectra = len(results)
    matched_spectra = sum(1 for row in results if row.get("match"))
    confidence_summary = _compound_confidence_summary(
        [row for row in results if row.get("match")]
    )
    # 去重 matched (同化合物多谱图只留最高分), 保留全部 unmatched
    results = _dedup_matched([r for r in results if r.get("match")]) + \
              [r for r in results if not r.get("match")]

    if format_type == "detailed":
        graph_artifact = _evograph_sidecar_path(TASK_DIR, task_id)
        if graph_artifact.is_file():
            try:
                graph_payload = _read_evograph_sidecar(
                    graph_artifact,
                    offset=0,
                    limit=max(1, total_spectra),
                )
                graph_by_index = {
                    int(record.get("result_index")): record
                    for record in (graph_payload.get("records") or [])
                    if record.get("result_index") is not None
                }
                for row in results:
                    graph = graph_by_index.get(
                        int(row.get("_result_index", -1))
                    )
                    if graph is not None:
                        row["evograph_fragmentation"] = graph
            except (OSError, ValueError, TypeError):
                logger.warning(
                    "Could not attach EvoGraph records to detailed report %s",
                    task_id,
                )

    if language == "en":
        from web.report_en import generate_english_html
        html_content = generate_english_html(
            results,
            task.get("filename", "Analysis"),
            detailed=format_type == "detailed",
            total_spectra=total_spectra,
            matched_spectra=matched_spectra,
            high_supported_compounds=confidence_summary.get(
                "high_supported_unique_compounds", 0
            ),
            high_confidence_scope_groups=compound_summary.get(
                "high_confidence_scope_group_count"
            ),
            high_confidence_exact_structures=compound_summary.get(
                "high_confidence_2d_connectivity_count",
                compound_summary.get(
                    "high_confidence_detected_count",
                    compound_summary.get("high_confidence_exact_structure_count"),
                ),
            ),
            high_confidence_isomer_groups=compound_summary.get(
                "high_confidence_isomer_group_count"
            ),
        )
        html_content = inject_expected_html(
            html_content,
            expected_payload,
            language="en",
            analysis_warnings=task.get("analysis_warnings") or [],
        )
        html_content = inject_reference_standard_qc_html(
            html_content, reference_qc_payload, language="en"
        )
        html_content = inject_confidence_guide_html(
            html_content, language="en"
        )
        html_content = inject_governance_blocked_html(
            html_content, governance_blocked_candidates, language="en"
        )
        filename = f"{Path(task['filename']).stem}_{format_type}_report_en.html"
        return HTMLResponse(
            content=html_content,
            headers={"Content-Disposition": attachment_disposition(filename)},
        )

    # Convert to format expected by report generator
    html_results = []
    for r in results:
        m = r.get("match")
        if m:
            # lib_peaks is stored at top level alongside match, move it in
            m["lib_peaks"] = r.get("lib_peaks", [])
            # ── Extract expert sub-scores from expert_details ──
            ed = m.get("expert_details") or {}
            if ed:
                m["neutral_loss_score"] = ed.get("neutral_loss_score", 0)
                m["diagnostic_ion_score"] = ed.get("diagnostic_ion_score", 0)
                m["rule_coverage_score"] = ed.get("rule_coverage_score", 0)
                m["fg_signature_score"] = ed.get("fg_signature_score", 0)
                m["structure_fragment_score"] = ed.get("structure_fragment_score", 0)
                m["bond_energy_score"] = ed.get("bond_energy_score", 0)
                m["feature_library_score"] = ed.get("feature_library_score", 0)
                m["isomer_differentiation_score"] = ed.get("isomer_differentiation_score", 0)
                m["chemical_reasoning"] = ed.get("chemical_reasoning", "")
        html_results.append({
            "precursor_mz": r["precursor_mz"],
            "search_precursor_mz": r.get(
                "search_precursor_mz", r.get("precursor_mz")
            ),
            "search_precursor_mz_source": r.get(
                "search_precursor_mz_source", ""
            ),
            "search_precursor_delta_ppm": r.get(
                "search_precursor_delta_ppm"
            ),
            "lcms_feature_used_for_search": bool(
                r.get("lcms_feature_used_for_search", False)
            ),
            "result_index": r.get("_result_index"),
            "spectrum_id": r.get("spectrum_id", ""),
            "peaks_count": r.get("peaks_count", 0),
            "retention_time": r.get("retention_time"),
            "retention_time_unit": r.get("retention_time_unit"),
            "precursor_scan_ref": r.get("precursor_scan_ref"),
            "collision_energy": r.get("collision_energy"),
            "collision_energy_unit": r.get("collision_energy_unit"),
            "ion_mode": r.get("ion_mode"),
            "precursor_charge": r.get("precursor_charge"),
            "tic": r.get("tic"),
            "base_peak": r.get("base_peak"),
            "query_peaks": r.get("query_peaks", []),
            "raw_peaks": r.get("raw_peaks", []),
            "raw_peaks_input_count": r.get("raw_peaks_input_count"),
            "raw_peaks_saved_count": r.get("raw_peaks_saved_count"),
            "raw_peaks_complete": r.get("raw_peaks_complete"),
            "raw_peaks_storage_policy": r.get("raw_peaks_storage_policy"),
            "evograph_fragmentation": r.get("evograph_fragmentation"),
            "peaks": r.get("query_peaks", []),  # alias for batch compat
            "matches": [m] if m else [],
            "best_cosine": m["cosine"] if m else 0.0,
            # 未命中谱图的 MSClassify 判类 (供 report.py 未命中区展示分类; 命中谱图为 None)
            "display_class": r.get("display_class"),
            "dl_confidence": r.get("dl_confidence"),
        })

    if format_type == "detailed":
        # Generate rich HTML with dual spectra
        from web.app import _gen_rich
        html_content = _gen_rich(
            html_results,
            governance_blocked_candidates=governance_blocked_candidates,
            source_name=task.get("filename", "Analysis"),
            total_spectra=total_spectra,
            matched_spectra=matched_spectra,
            high_supported_compounds=confidence_summary.get(
                "high_supported_unique_compounds", 0
            ),
            high_confidence_scope_groups=compound_summary.get(
                "high_confidence_scope_group_count"
            ),
            high_confidence_exact_structures=compound_summary.get(
                "high_confidence_2d_connectivity_count",
                compound_summary.get(
                    "high_confidence_detected_count",
                    compound_summary.get("high_confidence_exact_structure_count"),
                ),
            ),
            high_confidence_isomer_groups=compound_summary.get(
                "high_confidence_isomer_group_count"
            ),
        )
    else:
        from web.report import generate_html
        tmp = tempfile.NamedTemporaryFile(suffix=".html", delete=False, mode="w")
        try:
            generate_html(
                html_results,
                tmp.name,
                task.get("filename", "Analysis"),
                high_confidence_scope_groups=compound_summary.get(
                    "high_confidence_scope_group_count"
                ),
                high_confidence_exact_structures=compound_summary.get(
                    "high_confidence_2d_connectivity_count",
                    compound_summary.get(
                        "high_confidence_detected_count",
                        compound_summary.get("high_confidence_exact_structure_count"),
                    ),
                ),
                high_confidence_isomer_groups=compound_summary.get(
                    "high_confidence_isomer_group_count"
                ),
            )
            with open(tmp.name, "r", encoding="utf-8") as f:
                html_content = f.read()
        finally:
            try:
                os.unlink(tmp.name)
            except Exception:
                pass

    html_content = inject_expected_html(
        html_content,
        expected_payload,
        language="zh",
        analysis_warnings=task.get("analysis_warnings") or [],
    )
    html_content = inject_reference_standard_qc_html(
        html_content, reference_qc_payload, language="zh"
    )
    html_content = inject_confidence_guide_html(
        html_content, language="zh"
    )
    if format_type == "standard":
        html_content = inject_governance_blocked_html(
            html_content, governance_blocked_candidates, language="zh"
        )
    filename = f"{Path(task['filename']).stem}_{format_type}_report.html"
    return HTMLResponse(
        content=html_content,
        headers={"Content-Disposition": attachment_disposition(filename)},
    )


@app.get("/api/tasks/{task_id}/download/structure/result/{result_index}")
async def download_structure_image_by_result(
    task_id: str,
    result_index: int,
    format: str = Query("png", enum=["png", "svg"]),
):
    """Render a structure by immutable result-artifact index."""

    results = _get_task_results(task_id)
    if result_index < 0 or result_index >= len(results):
        raise HTTPException(status_code=404, detail="Result index out of range")
    match = results[result_index].get("match") or {}
    smiles = str(match.get("smiles") or "")
    if not smiles:
        raise HTTPException(status_code=404, detail="No SMILES available for this match")

    from web.structure import smiles_to_image_base64, smiles_to_svg
    from fastapi.responses import Response

    if format == "svg":
        svg = smiles_to_svg(smiles)
        if not svg:
            raise HTTPException(status_code=500, detail="Failed to generate SVG")
        return Response(content=svg, media_type="image/svg+xml")

    img_b64 = smiles_to_image_base64(smiles)
    if img_b64:
        import base64
        return Response(
            content=base64.b64decode(img_b64), media_type="image/png"
        )
    svg = smiles_to_svg(smiles)
    if not svg:
        raise HTTPException(status_code=500, detail="Failed to generate structure image")
    return Response(content=svg, media_type="image/svg+xml")


@app.get("/api/tasks/{task_id}/download/structure/{spectrum_index}")
async def download_structure_image(task_id: str, spectrum_index: int, format: str = Query("png", enum=["png", "svg"])):
    """Download structure image for a specific matched spectrum."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    results = _get_task_results(task_id)

    # ⚠️ 必须与前端 results.js 的 _currentMatched 排序方式完全一致:
    # 前端表格行号 = _currentMatched 按 final_score 降序后的 index, 并把同一 index 传给本接口。
    # 此前这里用"原始 results 顺序"取 matched[index], 导致 96% 行取到错的化合物/撞空 SMILES (实测151/157错位)。
    # 现统一按 final_score 降序; Python sorted 与现代 JS Array.sort 都稳定, final_score 相等者保持原始相对顺序, 两端一致。
    matched = sorted(
        [r for r in results if r.get("match")],
        key=lambda r: r["match"].get("final_score", 0.0),
        reverse=True,
    )
    if spectrum_index < 0 or spectrum_index >= len(matched):
        raise HTTPException(status_code=404, detail="Spectrum index out of range")

    smiles = matched[spectrum_index]["match"].get("smiles", "")
    if not smiles:
        raise HTTPException(status_code=404, detail="No SMILES available for this match")

    from web.structure import smiles_to_image_base64, smiles_to_svg

    if format == "svg":
        svg = smiles_to_svg(smiles)
        if not svg:
            raise HTTPException(status_code=500, detail="Failed to generate SVG")
        from fastapi.responses import Response
        return Response(content=svg, media_type="image/svg+xml")

    img_b64 = smiles_to_image_base64(smiles)
    if img_b64:
        import base64
        from fastapi.responses import Response
        img_bytes = base64.b64decode(img_b64)
        return Response(content=img_bytes, media_type="image/png")

    # PNG not available (no Cairo) — fall back to SVG
    svg = smiles_to_svg(smiles)
    if not svg:
        raise HTTPException(status_code=500, detail="Failed to generate structure image")
    from fastapi.responses import Response
    return Response(content=svg, media_type="image/svg+xml")


@app.get("/api/tasks/{task_id}/download/fragmentation/{spectrum_index}")
async def download_fragmentation_cascade(
    task_id: str,
    spectrum_index: int,
    result_index: int | None = Query(
        None,
        ge=0,
        description="Immutable source-result index; preferred over sorted display index",
    ),
):
    """碎裂产物结构级联 (CFM-ID): 前体→[-Δm]→碎片产物 2D 结构, 供网页详情异步加载。
    与 _gen_rich._pathway_reaction_svg 同源逻辑, 给 web UI 用。"""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    results = _get_task_results(task_id)
    if result_index is not None:
        if result_index >= len(results) or not results[result_index].get("match"):
            raise HTTPException(status_code=404, detail="Result index out of range or unmatched")
        selected_row = results[result_index]
    else:
        matched = sorted([r for r in results if r.get("match")],
                         key=lambda r: r["match"].get("final_score", 0.0), reverse=True)
        if spectrum_index < 0 or spectrum_index >= len(matched):
            raise HTTPException(status_code=404, detail="Spectrum index out of range")
        selected_row = matched[spectrum_index]
    mm = selected_row["match"]
    smiles = mm.get("smiles", "")
    if not smiles:
        raise HTTPException(status_code=404, detail="No SMILES available")
    adduct = str(mm.get("adduct", "") or "")
    # 极性派发统一走 adduct_polarity(), 优先用 query 真实 ion_mode (任务解析时持久化, 见 result["ion_mode"]),
    # 库 adduct 仅兜底。修旧启发式 bug: [M+HCOO]-/[M+FA-H]-/[M+Cl]- 含 '+' 被误判成正离子 → 送错正离子 CFM。
    from utils.adduct import adduct_polarity as _adduct_polarity
    _q_pol = (selected_row.get("ion_mode")
              or task.get("ion_mode")
              or _adduct_polarity(adduct))
    ion = "-" if _q_pol == "NEG" else "+"
    ed = mm.get("expert_details") or {}
    nl = ed.get("neutral_loss_analysis") or {}
    # 收集中性丢失 Δm (供 RDKit 兜底断键推断)
    losses = []
    for s in (nl.get("best_pathway") or []) + (nl.get("sequential_losses") or []):
        try:
            d = float(s.get("mass_difference", 0))
            if d > 0:
                losses.append(d)
        except Exception:
            pass
    for it in (nl.get("neutral_losses") or []):
        try:
            d = float(it.get("mass_difference", 0))
            if d > 0:
                losses.append(d)
        except Exception:
            pass

    query_peaks = selected_row.get("query_peaks") or []
    product_candidates = (
        (ed.get("structure_fragment_analysis") or {}).get(
            "plausible_fragments"
        ) or []
    )

    # 解耦 CFM (问题3): 单化合物详情 allow_docker=True → 缓存 miss 即时跑 docker, 不再被
    # 0.5Da 观测门卡死 (CFM 跑出预测产物就展示)。阻塞调用放线程池, 避免卡住事件循环。
    from web.report_common import docker_available, render_fragmentation_viz
    from fastapi.responses import HTMLResponse
    html = await asyncio.to_thread(
        render_fragmentation_viz,
        smiles,
        ion,
        losses,
        query_peaks,
        True,
        product_candidates=product_candidates,
    )
    if html:
        return HTMLResponse(content=html)
    if not docker_available():
        return HTMLResponse(content=(
            '<div style="padding:0.5rem 0.6rem;border:1px solid #fcd34d;background:#fffbeb;'
            'border-radius:5px;color:#92400e;font-size:0.72rem;line-height:1.5;">'
            '<b>CFM-ID 附加结构级联尚未生成：</b>Docker Desktop 当前未运行，无法调用 '
            'CFM 碎裂模型。上方的中性丢失、碎裂步骤、诊断离子和规则结构推理不受影响；'
            '启动 Docker 后重新打开详情即可重试。</div>'
        ))
    return HTMLResponse(content=(
        '<div style="padding:0.5rem 0.6rem;border:1px solid #cbd5e1;background:#f8fafc;'
        'border-radius:5px;color:#475569;font-size:0.72rem;line-height:1.5;">'
        '<b>CFM-ID 已运行，但未生成可用产物结构。</b>基础中性丢失、碎裂步骤、诊断离子和'
        '规则结构推理仍可在本面板下方查看。</div>'
    ))


@app.post("/api/tasks/{task_id}/precompute-cfm")
async def precompute_cfm(task_id: str):
    """手动触发 CFM 碎裂后台预计算 (给已完成/旧任务补 CFM 缓存用)。

    Docker 可用时非阻塞启动 precompute_cfm_fragments.py; 缓存就绪后报告/详情即显示
    CFM 产物级联。重复调用安全 (同 task 去重)。
    """
    with _tasks_lock:
        task = _tasks.get(task_id)
    rp = TASK_DIR / f"{task_id}_results.json"
    if not task and not rp.exists():
        raise HTTPException(status_code=404, detail="Task not found")
    launch = _maybe_launch_cfm_precompute(task_id, manual=True)
    from web.report_common import docker_available
    return JSONResponse({
        "task_id": task_id,
        "docker_available": docker_available(),
        "launch": launch,
        "message": (
            "CFM 后台预计算已启动；全局最多一个任务，完成后刷新报告即可见。"
            if launch.get("status") == "started"
            else (
                "分析队列仍在运行；为避免拖慢鉴定，CFM 预计算未启动。"
                if launch.get("status") == "deferred_analysis_queue_active"
                else (
                    "Docker 不可用 — 请先开启 Docker Desktop (含 cfmid 镜像), 再重试。"
                    if launch.get("status") == "docker_unavailable"
                    else "CFM 预计算未重复启动；请查看 launch 状态。"
                )
            )
        ),
    })


@app.post("/api/tasks/{task_id}/retry")
async def retry_task(task_id: str, params: AnalysisParams = AnalysisParams()):
    """Re-run a previously completed/failed task with optionally different params."""
    return await start_analysis(task_id, params)


@app.post("/api/tasks/{task_id}/cancel")
async def cancel_task(task_id: str):
    """Request cooperative cancellation before final artifact commit."""
    with _tasks_lock:
        task = _tasks.get(task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")
        state = task.get("status")
        if state == TaskState.PENDING.value:
            _transition_task(
                task,
                TaskState.CANCELLED,
                reason="user_cancelled_pending_task",
            )
            task["message"] = "Pending task cancelled"
            task["error"] = None
            accepted = False
        elif state in {TaskState.QUEUED.value, TaskState.RUNNING.value}:
            event = _task_cancel_events.setdefault(task_id, threading.Event())
            event.set()
            _transition_task(
                task,
                TaskState.CANCELLING,
                reason="user_requested_cooperative_cancellation",
            )
            task["message"] = "Cancellation requested; waiting for a safe checkpoint"
            accepted = True
        elif state == TaskState.CANCELLING.value:
            accepted = True
        elif state == TaskState.FINALIZING.value:
            raise HTTPException(
                status_code=409,
                detail="Task is committing result artifacts and can no longer be cancelled safely",
            )
        else:
            raise HTTPException(
                status_code=409,
                detail=f"Task is already {state} and is not cancellable",
            )
        payload = {
            "task_id": task_id,
            "status": task["status"],
            "cancellation_requested": accepted,
            "state_version": task.get("state_version"),
            "message": task.get("message"),
        }
    _save_tasks()
    _sse_broadcast("status_change", payload)
    return payload


@app.get("/api/tasks/{task_id}/state-history")
async def get_task_state_history(task_id: str):
    """Return the persisted state transition ledger for one task."""
    with _tasks_lock:
        task = _tasks.get(task_id)
        if not task:
            raise HTTPException(status_code=404, detail="Task not found")
        _initialize_task_state(task)
        return {
            "task_id": task_id,
            "status": task.get("status"),
            "state_version": task.get("state_version"),
            "run_attempt": task.get("run_attempt"),
            "active_run_id": task.get("active_run_id"),
            "result_run_id": task.get("result_run_id"),
            "history": list(task.get("state_history") or []),
        }


class TaskStateReconciliationRequest(BaseModel):
    expected_state_version: int = Field(ge=1)
    dry_run: bool = True


@app.post("/api/tasks/{task_id}/reconcile-completed-state")
async def reconcile_completed_task_state(task_id: str, request: TaskStateReconciliationRequest):
    """Recover an unrecorded status overwrite without re-running the analysis."""
    from runtime.task_reconciliation import (
        reconciled_completed_task,
        verify_completed_task_artifacts,
    )
    with _tasks_lock:
        current = _tasks.get(task_id)
        if current is None:
            raise HTTPException(404, "Task not found")
        if current.get("state_version") != request.expected_state_version:
            raise HTTPException(409, "Task state version changed")
        snapshot = deepcopy({key: value for key, value in current.items() if key != "results"})
    try:
        evidence = await asyncio.to_thread(verify_completed_task_artifacts, snapshot, TASK_DIR)
        if request.dry_run:
            return {"task_id": task_id, "dry_run": True, "eligible": True, "evidence": evidence}
        with _tasks_lock:
            original = _tasks.get(task_id)
            if original is None:
                raise ValueError("Task was removed during verification")
            repaired = reconciled_completed_task(original, evidence)
            _tasks[task_id] = repaired
        if _save_tasks() is False:
            with _tasks_lock:
                if _tasks.get(task_id) is repaired:
                    _tasks[task_id] = original
            raise HTTPException(500, "Reconciled state could not be persisted")
    except (ValueError, TypeError, KeyError, OSError) as exc:
        raise HTTPException(409, str(exc)) from exc
    payload = {"task_id": task_id, "status": repaired["status"],
               "state_version": repaired["state_version"], "analysis_rerun": False,
               "evidence": evidence}
    _sse_broadcast("status_change", payload)
    return payload


def _load_task_results_data(task_id: str) -> tuple[Optional[dict], list[dict]]:
    """Load task metadata and results (from memory or disk)."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        # 内存里没有该任务(reload/重启后内存丢失, 或前端传了旧/失效 id)时,
        # 也尝试从磁盘 results 文件恢复, 避免"结果已落盘但内存丢失"时
        # 批量导出/下载误报 404。
        rp = TASK_DIR / f"{task_id}_results.json"
        if rp.exists():
            try:
                with open(rp, "r", encoding="utf-8") as f:
                    results = json.load(f)
                return (
                    {"id": task_id, "filename": task_id, "status": "completed"},
                    _project_claims_for_task(task_id, results),
                )
            except Exception:
                pass
        return None, []
    results = task.get("results", [])
    if not results:
        rp = TASK_DIR / f"{task_id}_results.json"
        if rp.exists():
            try:
                with open(rp, "r", encoding="utf-8") as f:
                    results = json.load(f)
            except Exception:
                pass
    return task, _project_claims_for_task(task_id, results)


def _build_html_report(
    task: dict,
    results: list[dict],
    compound_summary: dict[str, Any] | None = None,
) -> str:
    """PDF 失败兜底 HTML 报告 — 直接复用 standard 报告生成器 (report.py:generate_html)。

    消除原先独立的第 4 套不一致模板 (裸 confidence / 无结构 / 无推理 / 配色异于 detailed),
    使 PDF 失败兜底与 standard HTML 报告完全同源 (honest 双置信度 + 同列 + 同配色)。"""
    import tempfile
    html_results = []
    for r in results:
        m = r.get("match")
        if m:
            html_results.append({
                "precursor_mz": r["precursor_mz"],
                "spectrum_id": r.get("spectrum_id", ""),
                "peaks_count": r.get("peaks_count", 0),
                "retention_time": r.get("retention_time"),
                "retention_time_unit": r.get("retention_time_unit"),
                "precursor_scan_ref": r.get("precursor_scan_ref"),
                "collision_energy": r.get("collision_energy"),
                "collision_energy_unit": r.get("collision_energy_unit"),
                "tic": r.get("tic", 0),
                "base_peak": r.get("base_peak", 0),
                "query_peaks": r.get("query_peaks", []),
                "raw_peaks": r.get("raw_peaks", []),
                "peaks": r.get("query_peaks", []),
                "matches": [m],
                "best_cosine": m.get("cosine", 0.0),
                "display_class": r.get("display_class"),
                "dl_confidence": r.get("dl_confidence"),
            })
    from web.report import generate_html
    tmp = tempfile.NamedTemporaryFile(suffix=".html", delete=False, mode="w")
    try:
        compound_summary = dict(
            compound_summary
            or task.get("compound_conclusion_summary")
            or {}
        )
        generate_html(
            html_results,
            tmp.name,
            task.get("filename", "Analysis"),
            high_confidence_scope_groups=compound_summary.get(
                "high_confidence_scope_group_count"
            ),
            high_confidence_exact_structures=compound_summary.get(
                "high_confidence_2d_connectivity_count",
                compound_summary.get(
                    "high_confidence_detected_count",
                    compound_summary.get("high_confidence_exact_structure_count"),
                ),
            ),
            high_confidence_isomer_groups=compound_summary.get(
                "high_confidence_isomer_group_count"
            ),
        )
        with open(tmp.name, "r", encoding="utf-8") as f:
            html_content = f.read()
        task_id = str(task.get("id") or "")
        expected_payload = _load_expected_compounds_payload(
            task_id, task, tolerate_errors=True
        )
        html_content = inject_expected_html(
            html_content,
            expected_payload,
            language="zh",
            analysis_warnings=task.get("analysis_warnings") or [],
        )
        reference_qc_payload = _load_reference_standard_qc_payload(
            task_id, tolerate_errors=True
        )
        html_content = inject_reference_standard_qc_html(
            html_content, reference_qc_payload, language="zh"
        )
        return inject_confidence_guide_html(
            html_content, language="zh"
        )
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def _governance_blocked_pdf_story(
    results: list[dict],
    *,
    language: str,
    regular_font: str,
    bold_font: str,
) -> list[Any]:
    """Render the same non-release projection used by API, HTML and tables."""
    from html import escape as html_escape

    from reportlab.lib import colors
    from reportlab.lib.units import cm
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.platypus import Paragraph, Spacer, Table, TableStyle

    from web.governance_blocked_projection import (
        project_governance_blocked_candidates,
    )

    rows = project_governance_blocked_candidates(results)
    if not rows:
        return []
    english = str(language).lower().startswith("en")
    styles = getSampleStyleSheet()
    cell = ParagraphStyle(
        f"GovernanceBlockedCell-{language}",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=5.8,
        leading=7.0,
        wordWrap="CJK",
    )
    header = ParagraphStyle(
        f"GovernanceBlockedHeader-{language}",
        parent=cell,
        fontName=bold_font,
        textColor=colors.white,
        alignment=1,
    )
    heading_style = ParagraphStyle(
        f"GovernanceBlockedHeading-{language}",
        parent=styles["Heading2"],
        fontName=bold_font,
    )
    heading = (
        f"Governance not released ({len(rows)} saved scan candidate records)"
        if english
        else f"治理未放行候选（{len(rows)} 条保存的扫描候选记录）"
    )
    explanation = (
        "These candidates retain their original policy and blockers. They are "
        "shown for review and are not counted as released identifications."
        if english else
        "这些候选保留原治理策略与拦截原因，仅供复核，不计入已发布鉴定数。"
    )
    labels = (
        ["Index", "Spectrum", "Feature", "Candidate", "Formula", "Status", "Original blockers"]
        if english else
        ["索引", "谱图", "特征", "候选名称", "分子式", "治理状态", "原始拦截原因"]
    )
    table_rows = [[Paragraph(html_escape(label), header) for label in labels]]
    for row in rows[:200]:
        name = str(row.get("compound_name") or "")
        aliases = "·".join(map(str, row.get("compound_name_cn") or []))
        if aliases:
            name = f"{name} ({aliases})"
        values = [
            row.get("source_result_index", ""),
            row.get("spectrum_id", ""),
            row.get("source_feature_id", ""),
            name,
            row.get("formula", ""),
            row.get("claim_status", ""),
            row.get("reason_text", ""),
        ]
        table_rows.append([
            Paragraph(html_escape(str(value or "")), cell) for value in values
        ])
    table = Table(
        table_rows,
        colWidths=[0.7 * cm, 1.7 * cm, 1.5 * cm, 3.2 * cm,
                   1.6 * cm, 2.8 * cm, 5.1 * cm],
        repeatRows=1,
    )
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#7C2D12")),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#D1D5DB")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [
            colors.white, colors.HexColor("#FFF7ED")
        ]),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ("TOPPADDING", (0, 0), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]))
    story: list[Any] = [
        Spacer(1, 0.25 * cm),
        Paragraph(html_escape(heading), heading_style),
        Paragraph(html_escape(explanation), cell),
        Spacer(1, 0.1 * cm),
        table,
    ]
    if len(rows) > 200:
        note = (
            f"PDF shows the first 200 of {len(rows)} records; CSV and Excel retain all rows."
            if english else
            f"PDF 显示前 200/{len(rows)} 条；CSV 与 Excel 保留全部记录。"
        )
        story.append(Paragraph(html_escape(note), cell))
    return story


def _build_pdf_report(
    task: dict,
    results: list[dict],
    language: str = "zh",
    *,
    compound_summary: dict[str, Any] | None = None,
    reference_standard_qc_payload: dict[str, Any] | None = None,
) -> bytes:
    """Generate a PDF report using reportlab (pure Python, no system deps)."""
    if language == "en":
        return _build_pdf_report_en(
            task,
            results,
            compound_summary=compound_summary,
            reference_standard_qc_payload=reference_standard_qc_payload,
        )
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.units import mm, cm
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from web.report_common import (
        _standard_compound_name,
        exact_structure_confidence,
        final_claim_decision,
        isomer_candidate_export_rows,
        isomer_candidate_export_rows_from_compounds,
        reportable_identity_confidence,
        reportable_identity_label,
    )
    from web.raw_spectrum_locator import raw_spectrum_locator

    all_matched = [r for r in results if r.get("match")]
    matched = _dedup_matched(all_matched)
    task_id = str(task.get("id") or "").strip()
    compound_conclusions = (
        _compound_conclusions_for_export(task_id, task)
        if task_id
        else None
    )
    confidence_summary = _compound_confidence_summary(all_matched)
    total = len(results)
    matched_scan_count = len(all_matched)
    rate = matched_scan_count / total * 100 if total else 0
    filename = task.get("filename", "Unknown")

    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, leftMargin=1.5*cm, rightMargin=1.5*cm,
                            topMargin=1.5*cm, bottomMargin=1.5*cm)
    styles = getSampleStyleSheet()
    regular_font, bold_font = reportlab_font_names("zh")
    for style_name in ("Normal", "Title", "Heading1", "Heading2"):
        styles[style_name].fontName = regular_font
    locator_cell_style = ParagraphStyle(
        "RawSpectrumLocatorCellZh",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=6.2,
        leading=7.6,
    )
    result_cell_style = ParagraphStyle(
        "ResultTableCellZh",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=6.0,
        leading=7.1,
        wordWrap="CJK",
        splitLongWords=True,
    )
    result_number_style = ParagraphStyle(
        "ResultTableNumberZh",
        parent=result_cell_style,
        alignment=1,
    )
    result_header_style = ParagraphStyle(
        "ResultTableHeaderZh",
        parent=result_cell_style,
        fontName=bold_font,
        fontSize=5.8,
        leading=6.8,
        alignment=1,
        textColor=colors.white,
    )
    elements = []

    # Title
    elements.append(Paragraph("质谱分析报告", styles["Title"]))
    elements.append(Paragraph(
        f"文件: {filename}  |  获得谱库候选的 MS2: "
        f"{matched_scan_count}/{total} ({rate:.1f}%)  |  去重候选: {len(matched)}",
        styles["Normal"],
    ))
    elements.append(Spacer(1, 0.5*cm))

    # 正常报告只呈现一个化合物置信度；完整内部校准留在审计导出。
    compound_summary = dict(
        compound_summary
        or task.get("compound_conclusion_summary")
        or {}
    )
    legacy_high = int(
        confidence_summary["high_supported_unique_compounds"] or 0
    )
    n_high_scope = int(
        compound_summary.get(
            "high_confidence_scope_group_count", legacy_high
        ) or 0
    )
    n_high_exact = int(
        compound_summary.get(
            "high_confidence_2d_connectivity_count",
            compound_summary.get(
                "high_confidence_detected_count",
                compound_summary.get(
                    "high_confidence_exact_structure_count", legacy_high
                ),
            ),
        ) or 0
    )
    n_high_isomer = int(
        compound_summary.get(
            "high_confidence_isomer_group_count",
            max(0, n_high_scope - n_high_exact),
        ) or 0
    )
    final_claims = [
        final_claim_decision(row.get("match")) for row in matched
        if final_claim_decision(row.get("match"))
    ]
    released_final_claims = sum(bool(item.get("released")) for item in final_claims)
    summary_data = [
        ["总 MS2 谱图", "有谱库候选的 MS2", "候选获得率", "High 级鉴定组"],
        [
            str(total),
            str(matched_scan_count),
            f"{rate:.1f}%",
            str(n_high_scope),
        ],
    ]
    summary_table = Table(summary_data, colWidths=[3.5*cm]*4)
    summary_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#4299e1")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), bold_font),
        ("FONTNAME", (0, 1), (-1, -1), regular_font),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7fafc")]),
    ]))
    elements.append(summary_table)
    elements.append(Paragraph(
        f"High 级鉴定组组成：单一二维连接 {n_high_exact}；"
        f"未解异构体组 {n_high_isomer}。异构体组不等同于唯一结构确证。",
        styles["Normal"],
    ))
    elements.append(Paragraph(
        "说明：化合物名称为当前结果保留的最佳谱库候选，置信度列表示"
        "该名称的鉴定证据强弱。鉴定结论与证据等级分开；"
        "项目级 q 值或独立进样不足时，结果仍作为候选。报告生成不会重新打分或改变候选排序。",
        styles["Normal"],
    ))
    elements.append(Paragraph(
        f"鉴定结论：{released_final_claims}/{len(final_claims)} 条已评估代表候选达到当前结论条件；其余未达到当前结论条件。",
        styles["Normal"],
    ))
    elements.append(Spacer(1, 0.5*cm))
    elements.extend(confidence_guide_pdf_story(
        "zh", regular_font=regular_font, bold_font=bold_font
    ))

    # 靶向库结果置于 MS2 表格之前，避免藏在报告末尾。
    expected_payload = _load_expected_compounds_payload(
        str(task.get("id") or ""), task, tolerate_errors=True
    )
    elements.extend(expected_pdf_story(
        expected_payload,
        language="zh",
        analysis_warnings=task.get("analysis_warnings") or [],
    ))
    elements.extend(reference_standard_qc_pdf_story(
        reference_standard_qc_payload or {}, language="zh"
    ))
    elements.extend(_governance_blocked_pdf_story(
        results,
        language="zh",
        regular_font=regular_font,
        bold_font=bold_font,
    ))

    # Results table
    if matched:
        header = [
            Paragraph(_esc(label), result_header_style)
            for label in [
                "#", "原始 MS2 定位", "化合物名称", "分子式", "Class",
                "来源", "MS/MS 相似度", "综合匹配分", "质量误差",
                "可报告置信度", "鉴定结论",
            ]
        ]
        table_data = [header]
        for i, r in enumerate(matched[:200]):  # cap at 200 rows
            m = r["match"]
            locator = raw_spectrum_locator(r, filename)
            rt_text = (
                f"{locator['retention_time_min']:.3f} min"
                if locator["retention_time_min"] is not None
                else "未记录"
            )
            spectrum_text = str(locator["spectrum_id"] or "未记录")
            reportable_confidence_text = reportable_identity_label(
                m, language="zh"
            )
            full_identity_confidence = exact_structure_confidence(m)
            if reportable_identity_confidence(m) != full_identity_confidence:
                reportable_confidence_text += (
                    f"<br/><font color='#667085'>完整身份 "
                    f"{_esc(full_identity_confidence)}</font>"
                )
            claim_text = (
                "可作结论" if final_claim_decision(m).get("released")
                else "未达到当前结论条件" if final_claim_decision(m)
                else "尚未评估"
            )
            table_data.append([
                Paragraph(str(i + 1), result_number_style),
                Paragraph(
                    f"m/z {float(r['precursor_mz']):.5f}<br/>RT "
                    f"{_esc(rt_text)}<br/>MS2 {_esc(spectrum_text)}",
                    result_cell_style,
                ),
                Paragraph(_esc(_standard_compound_name(
                    m.get("compound_name", ""), m
                )[:42]), result_cell_style),
                Paragraph(_esc(str(m.get("formula", ""))), result_cell_style),
                Paragraph(_esc(str(_display_class(m))[:24]), result_cell_style),
                Paragraph(_esc(str(_dataset_source(m) or "—")[:28]), result_cell_style),
                Paragraph(f"{m.get('cosine', 0):.4f}", result_number_style),
                Paragraph(f"{m.get('final_score', 0):.4f}", result_number_style),
                Paragraph(_esc(str(m.get("mass_ppm", "—"))), result_number_style),
                Paragraph(
                    reportable_confidence_text,
                    result_cell_style,
                ),
                Paragraph(_esc(claim_text), result_cell_style),
            ])

        col_widths = [
            0.45*cm, 2.2*cm, 2.9*cm, 1.35*cm, 1.45*cm, 1.75*cm,
            1.1*cm, 1.1*cm, 0.85*cm, 2.1*cm, 2.1*cm,
        ]
        results_table = Table(table_data, colWidths=col_widths, repeatRows=1)
        results_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2b6cb0")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), bold_font),
            ("FONTNAME", (0, 1), (-1, -1), regular_font),
            ("FONTSIZE", (0, 0), (-1, -1), 7),
            ("ALIGN", (0, 0), (0, -1), "CENTER"),
            ("ALIGN", (1, 0), (1, -1), "RIGHT"),
            ("ALIGN", (4, 0), (6, -1), "RIGHT"),
            ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#cbd5e0")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7fafc")]),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
            ("TOPPADDING", (0, 0), (-1, -1), 2.5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
        ]))
        elements.append(Paragraph("匹配结果", styles["Heading2"]))
        elements.append(results_table)

        isomer_rows = (
            isomer_candidate_export_rows_from_compounds(
                compound_conclusions, language="zh"
            )
            if compound_conclusions is not None
            else isomer_candidate_export_rows(matched, language="zh")
        )
        if isomer_rows:
            elements.append(Spacer(1, 0.35 * cm))
            elements.append(Paragraph(
                "未解析异构体竞争候选（完整逐候选清单）",
                styles["Heading2"],
            ))
            elements.append(Paragraph(
                "每组第 1 名是当前排序最佳候选，随后为全部保留的竞争异构体；"
                "排序第 1 不等于已经唯一解析。该清单独立于结构推理文字是否可用，"
                "因此不会因缺少专家文字而从 PDF 中消失。",
                styles["Normal"],
            ))
            grouped_isomers: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
            for candidate in isomer_rows:
                key = (
                    candidate.get("Result_Spectrum"),
                    candidate.get("Precursor_mz"),
                    candidate.get("Retention_Time"),
                )
                grouped_isomers.setdefault(key, []).append(candidate)
            for group_index, candidates in enumerate(
                grouped_isomers.values(), 1
            ):
                first = candidates[0]
                lead = next(
                    (
                        candidate for candidate in candidates
                        if candidate.get("Is_Ranking_Lead")
                    ),
                    first,
                )
                elements.append(Paragraph(
                    f"<b>第 {group_index} 组：</b>"
                    f"{_esc(lead.get('Compound_Name') or '未命名候选')}（当前最佳候选）；"
                    f"共 {int(first.get('Candidate_Count') or len(candidates))} 个候选；"
                    f"MS2 {_esc(first.get('Result_Spectrum') or '未记录')}；"
                    f"m/z {_esc(first.get('Precursor_mz') or '—')}；"
                    f"RT {_esc(first.get('Retention_Time') or '—')} min",
                    styles["Normal"],
                ))
                candidate_table_data = [[
                    Paragraph(_esc(label), result_header_style)
                    for label in [
                        "排名／角色", "候选化合物", "分子式／InChIKey",
                        "Cosine", "Final／Rank", "加合物",
                        "候选关系", "谱库来源／谱图",
                    ]
                ]]
                for candidate in candidates:
                    candidate_table_data.append([
                        Paragraph(
                            f"{_esc(candidate.get('Candidate_Rank') or '—')}<br/>"
                            f"{_esc(candidate.get('Candidate_Role') or '—')}",
                            result_cell_style,
                        ),
                        Paragraph(
                            f"{_esc(candidate.get('Compound_Name') or '—')}<br/>"
                            f"<font color='#667085'>{_esc(candidate.get('Candidate_Record_Status') or '')}</font>",
                            result_cell_style,
                        ),
                        Paragraph(
                            f"{_esc(candidate.get('Formula') or '—')}<br/>"
                            f"{_esc(candidate.get('InChIKey') or '—')}",
                            result_cell_style,
                        ),
                        Paragraph(_esc(candidate.get("Cosine") or "—"), result_number_style),
                        Paragraph(
                            f"{_esc(candidate.get('Final_Score') or '—')} / "
                            f"{_esc(candidate.get('Rank_Score') or '—')}",
                            result_number_style,
                        ),
                        Paragraph(_esc(candidate.get("Adduct") or "—"), result_cell_style),
                        Paragraph(_esc(candidate.get("Relationship") or "—"), result_cell_style),
                        Paragraph(
                            f"{_esc(candidate.get('Dataset_Source') or '—')}<br/>"
                            f"{_esc(candidate.get('Library_Spectrum') or '—')}",
                            result_cell_style,
                        ),
                    ])
                candidate_table = Table(
                    candidate_table_data,
                    colWidths=[
                        1.7 * cm, 3.0 * cm, 3.1 * cm, 1.1 * cm,
                        1.45 * cm, 1.25 * cm, 2.4 * cm, 2.7 * cm,
                    ],
                    repeatRows=1,
                )
                candidate_table.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#173F5F")),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("FONTNAME", (0, 0), (-1, 0), bold_font),
                    ("FONTNAME", (0, 1), (-1, -1), regular_font),
                    ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#CBD5E0")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7FAFC")]),
                    ("LEFTPADDING", (0, 0), (-1, -1), 2),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 2),
                    ("TOPPADDING", (0, 0), (-1, -1), 2.5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
                ]))
                elements.append(candidate_table)
                elements.append(Spacer(1, 0.2 * cm))

        # v3 结构推理 + Top-K 候选详情 (仅 High/Medium 命中, 限量30避免PDF过长)
        _rs = ParagraphStyle("rs", parent=styles["Normal"], fontSize=7, leading=9, textColor=colors.HexColor("#234e52"))
        _ts = ParagraphStyle("ts", parent=styles["Normal"], fontSize=7, leading=9, textColor=colors.HexColor("#2a4365"))
        _n_det = 0
        for r in matched[:200]:
            if _n_det >= 30:
                break
            m = r["match"]
            if reportable_identity_confidence(m) not in ("High", "Medium"):
                continue
            _v3 = (m.get("expert_details") or {}).get("v3") or {}
            _reason = "; ".join(
                str(x)[:50]
                for x in (_v3.get("reasoning") or [])[:3]
                if "异构" not in str(x)
                and "isomer" not in str(x).lower()
            )
            if not _reason:
                continue
            if _n_det == 0:
                elements.append(Spacer(1, 0.4*cm))
                elements.append(Paragraph(
                    "结构推理与 Top-K 候选（High/Medium 候选）",
                    styles["Heading2"],
                ))
            _n_det += 1
            _nm = _standard_compound_name(
                m.get("compound_name", ""), m
            )[:30].replace("&", "&amp;").replace("<", "&lt;")
            locator = raw_spectrum_locator(r, filename)
            rt_text = (
                f"{locator['retention_time_min']:.3f} min"
                if locator["retention_time_min"] is not None
                else "未记录"
            )
            elements.append(Paragraph(
                f"<b>{_nm}</b>（原始 MS2：RT {rt_text}；"
                f"前体 m/z {r['precursor_mz']:.5f}；"
                f"谱图编号 {str(locator['spectrum_id'] or '未记录')}；"
                f"cos {m.get('cosine',0):.3f}）",
                styles["Normal"],
            ))
            elements.append(Paragraph(
                f"<b>结构推理：</b>{_esc(_reason)}", _rs
            ))
            # 馆藏预存裂解途径摘要 (诊断离子 + 主要中性丢失, 直接读 DB fragmentation_pathway)
            try:
                from web.report_common import _parse_pathway_json
                _pw = _parse_pathway_json(m.get("pathway_json"))
                if _pw and (_pw.get("diagnostic_ions") or _pw.get("neutral_losses")):
                    _diags = ", ".join(f"{float(d.get('mz',0)):.1f}" for d in (_pw.get("diagnostic_ions") or [])[:4] if d.get("mz"))
                    _nls = ", ".join(f"-{float(n.get('loss_da',0)):.0f}" for n in (_pw.get("neutral_losses") or [])[:4] if n.get("loss_da"))
                    elements.append(Paragraph(
                        (
                            f"<b>馆藏碎裂证据：</b>诊断离子 {_diags or '—'}；"
                            f"中性丢失 {_nls or '—'}（馆藏 "
                            f"{_pw.get('n_spectra','?')} 张谱图）"
                        ),
                        _rs,
                    ))
            except Exception:
                pass
            _seen = set(); _tc = []
            for _c in (m.get("top_candidates") or [])[:6]:
                _rn = (_c.get("compound_name") or "").lower().split("(")[0].split("[")[0].strip()
                _nn = "".join(ch for ch in _rn if ch.isalnum())
                _kk = f"{_nn}|{_c.get('formula')}"
                if _kk in _seen:
                    continue
                _seen.add(_kk)
                _candidate_name = _standard_compound_name(
                    _c.get("compound_name", ""), _c
                )[:22]
                _tc.append(
                    f"{_esc(_candidate_name)}"
                    f"({float(_c.get('cosine', 0) or 0):.2f})"
                )
            if _tc:
                elements.append(Paragraph(
                    f"<b>其他谱库候选：</b>{' · '.join(_tc)}",
                    _ts,
                ))
            elements.append(Spacer(1, 0.15*cm))

        if len(matched) > 200:
            elements.append(Spacer(1, 0.3*cm))
            elements.append(Paragraph(f"(仅显示前 200 条，共 {len(matched)} 条匹配)", styles["Normal"]))

    doc.build(elements)
    return buf.getvalue()


def _build_pdf_report_en(
    task: dict,
    results: list[dict],
    *,
    compound_summary: dict[str, Any] | None = None,
    reference_standard_qc_payload: dict[str, Any] | None = None,
) -> bytes:
    """Generate a system-authored fully English PDF from structured evidence."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import cm
    from reportlab.platypus import (
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from web.report_en import class_name_en
    from web.report_common import (
        _standard_compound_name,
        exact_structure_confidence,
        final_claim_decision,
        isomer_candidate_export_rows,
        isomer_candidate_export_rows_from_compounds,
        reportable_identity_confidence,
        reportable_identity_label,
    )
    from web.raw_spectrum_locator import raw_spectrum_locator

    all_matched = [row for row in results if row.get("match")]
    matched = _dedup_matched(all_matched)
    task_id = str(task.get("id") or "").strip()
    compound_conclusions = (
        _compound_conclusions_for_export(task_id, task)
        if task_id
        else None
    )
    confidence_summary = _compound_confidence_summary(all_matched)
    total = len(results)
    matched_scan_count = len(all_matched)
    unmatched_scan_count = max(0, total - matched_scan_count)
    rate = matched_scan_count / total * 100 if total else 0.0
    filename = str(task.get("filename") or "Unknown")
    compound_summary = dict(
        compound_summary
        or task.get("compound_conclusion_summary")
        or {}
    )
    legacy_high = int(
        confidence_summary["high_supported_unique_compounds"] or 0
    )
    n_high_scope = int(
        compound_summary.get(
            "high_confidence_scope_group_count", legacy_high
        ) or 0
    )
    n_high_exact = int(
        compound_summary.get(
            "high_confidence_2d_connectivity_count",
            compound_summary.get(
                "high_confidence_detected_count",
                compound_summary.get(
                    "high_confidence_exact_structure_count", legacy_high
                ),
            ),
        ) or 0
    )
    n_high_isomer = int(
        compound_summary.get(
            "high_confidence_isomer_group_count",
            max(0, n_high_scope - n_high_exact),
        ) or 0
    )
    final_claims = [
        final_claim_decision(row.get("match")) for row in matched
        if final_claim_decision(row.get("match"))
    ]
    released_final_claims = sum(bool(item.get("released")) for item in final_claims)

    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        leftMargin=1.25 * cm,
        rightMargin=1.25 * cm,
        topMargin=1.25 * cm,
        bottomMargin=1.25 * cm,
    )
    styles = getSampleStyleSheet()
    regular_font, bold_font = reportlab_font_names("en")
    for style_name in ("Normal", "Title", "Heading1", "Heading2"):
        styles[style_name].fontName = regular_font
    small = ParagraphStyle(
        "EvidenceSmall",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=7,
        leading=9,
        textColor=colors.HexColor("#344054"),
    )
    result_cell = ParagraphStyle(
        "ResultTableCellEn",
        parent=small,
        fontSize=5.8,
        leading=6.8,
        wordWrap="CJK",
        splitLongWords=True,
    )
    result_number = ParagraphStyle(
        "ResultTableNumberEn",
        parent=result_cell,
        alignment=1,
    )
    result_header = ParagraphStyle(
        "ResultTableHeaderEn",
        parent=result_cell,
        fontName=bold_font,
        fontSize=5.6,
        leading=6.5,
        alignment=1,
        textColor=colors.white,
    )
    note = ParagraphStyle(
        "EvidenceNote",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=7,
        leading=9,
        textColor=colors.HexColor("#475467"),
        backColor=colors.HexColor("#F8FAFC"),
        borderColor=colors.HexColor("#CBD5E1"),
        borderWidth=0.5,
        borderPadding=5,
    )
    elements = [
        Paragraph("LC-MS/MS Identification Report", styles["Title"]),
        Paragraph(
            f"Source: {_esc(filename)} &nbsp;|&nbsp; "
            f"Matched MS/MS spectra: {matched_scan_count}/{total} ({rate:.1f}%) "
            f"&nbsp;|&nbsp; Unique matched compounds: {len(matched)}",
            styles["Normal"],
        ),
        Spacer(1, 0.35 * cm),
    ]
    summary = Table(
        [
            [
                "Total spectra",
                "Matched MS/MS spectra",
                "Unique matched compounds",
                "Unmatched MS/MS spectra",
                "High-evidence identification groups",
            ],
            [
                str(total),
                f"{matched_scan_count} ({rate:.1f}%)",
                str(len(matched)),
                str(unmatched_scan_count),
                str(n_high_scope),
            ],
        ],
        colWidths=[3.2 * cm] * 5,
    )
    summary.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1A5276")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), bold_font),
                ("FONTNAME", (0, 1), (-1, -1), regular_font),
                ("FONTSIZE", (0, 0), (-1, -1), 7),
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#CBD5E1")),
            ]
        )
    )
    elements.extend(
        [
            summary,
            Spacer(1, 0.25 * cm),
            Paragraph(
                f"<b>High-group composition.</b> {n_high_exact} single 2D "
                f"structure(s) and {n_high_isomer} unresolved isomer group(s). "
                "An isomer-group result does not establish one unique structure.",
                note,
            ),
            Spacer(1, 0.18 * cm),
            Paragraph(
                "<b>Interpretation.</b> A library-spectrum match is normally a "
                "putative Level 2 annotation. The displayed compound name is the "
                "best library candidate retained by the saved analysis, and the "
                "confidence column communicates the strength of that assignment. "
                f"The separate identification-conclusion check released {released_final_claims}/"
                f"{len(final_claims)} assessed representative candidates; "
                "candidates lacking project q-values or independent injections are withheld. "
                "Level 1 confirmation requires an authentic standard under the same "
                "analytical method. Stored fragmentation "
                "signatures are empirical recurring ions and neutral-loss sequences, "
                "not atom-mapped bond-cleavage mechanisms.",
                note,
            ),
            Spacer(1, 0.3 * cm),
        ]
    )
    elements.extend(confidence_guide_pdf_story(
        "en", regular_font=regular_font, bold_font=bold_font
    ))

    expected_payload = _load_expected_compounds_payload(
        str(task.get("id") or ""), task, tolerate_errors=True
    )
    elements.extend(expected_pdf_story(
        expected_payload,
        language="en",
        analysis_warnings=task.get("analysis_warnings") or [],
    ))
    elements.extend(reference_standard_qc_pdf_story(
        reference_standard_qc_payload or {}, language="en"
    ))
    elements.extend(_governance_blocked_pdf_story(
        results,
        language="en",
        regular_font=regular_font,
        bold_font=bold_font,
    ))

    if matched:
        table_data = [[
            Paragraph(_esc(label), result_header)
            for label in [
                "#", "Raw MS2 locator", "Compound", "Formula", "Class",
                "Source", "Cosine", "Final", "ppm",
                "Reportable confidence", "Identification conclusion",
            ]
        ]]
        for index, row in enumerate(matched[:200], 1):
            match = row["match"]
            locator = raw_spectrum_locator(row, filename)
            rt_text = (
                f"{locator['retention_time_min']:.3f} min"
                if locator["retention_time_min"] is not None
                else "not recorded"
            )
            reportable_confidence_text = reportable_identity_label(
                match, language="en"
            )
            full_identity_confidence = exact_structure_confidence(match)
            if (
                reportable_identity_confidence(match)
                != full_identity_confidence
            ):
                reportable_confidence_text += (
                    f"<br/><font color='#667085'>Full identity: "
                    f"{_esc(full_identity_confidence)}</font>"
                )
            claim_text = (
                "Released" if final_claim_decision(match).get("released")
                else "Criteria not met" if final_claim_decision(match)
                else "Not assessed"
            )
            table_data.append(
                [
                    Paragraph(str(index), result_number),
                    Paragraph(
                        f"m/z {float(row.get('precursor_mz') or 0):.5f}<br/>"
                        f"RT {rt_text}<br/>MS2 "
                        f"{_esc(locator['spectrum_id'] or 'not recorded')}",
                        result_cell,
                    ),
                    Paragraph(_esc(_standard_compound_name(
                        match.get("compound_name", ""), match
                    )[:42]), result_cell),
                    Paragraph(_esc(str(match.get("formula") or "")[:22]), result_cell),
                    Paragraph(_esc(class_name_en(
                        _display_class(match, language="en")
                    )[:24]), result_cell),
                    Paragraph(_esc(str(
                        _dataset_source(match) or "Not recorded"
                    )[:28]), result_cell),
                    Paragraph(
                        f"{float(match.get('cosine') or 0):.4f}",
                        result_number,
                    ),
                    Paragraph(
                        f"{float(match.get('final_score') or 0):.4f}",
                        result_number,
                    ),
                    Paragraph(_esc(
                        f"{float(match['mass_ppm']):.1f}"
                        if match.get("mass_ppm") is not None
                        else "—"
                    ), result_number),
                    Paragraph(
                        reportable_confidence_text,
                        result_cell,
                    ),
                    Paragraph(_esc(claim_text), result_cell),
                ]
            )
        results_table = Table(
            table_data,
            colWidths=[
                0.6 * cm,
                2.25 * cm,
                2.85 * cm,
                1.75 * cm,
                1.65 * cm,
                1.75 * cm,
                1.15 * cm,
                1.15 * cm,
                0.9 * cm,
                1.75 * cm,
                1.3 * cm,
            ],
            repeatRows=1,
        )
        results_table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1A5276")),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("FONTNAME", (0, 0), (-1, 0), bold_font),
                    ("FONTNAME", (0, 1), (-1, -1), regular_font),
                    ("FONTSIZE", (0, 0), (-1, -1), 6.5),
                    ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#D0D5DD")),
                    (
                        "ROWBACKGROUNDS",
                        (0, 1),
                        (-1, -1),
                        [colors.white, colors.HexColor("#F8FAFC")],
                    ),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 2),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 2),
                    ("TOPPADDING", (0, 0), (-1, -1), 2.5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
                ]
            )
        )
        elements.extend(
            [
                Paragraph("Matched compounds", styles["Heading2"]),
                results_table,
            ]
        )

        isomer_rows = (
            isomer_candidate_export_rows_from_compounds(
                compound_conclusions, language="en"
            )
            if compound_conclusions is not None
            else isomer_candidate_export_rows(matched, language="en")
        )
        if isomer_rows:
            elements.extend([
                Spacer(1, 0.35 * cm),
                Paragraph(
                    "Unresolved isomer competition — complete candidate list",
                    styles["Heading2"],
                ),
                Paragraph(
                    "The first row in each group is the current ranking lead, "
                    "followed by every retained competitive isomer. A ranking "
                    "lead is not a uniquely resolved structure.",
                    note,
                ),
            ])
            grouped_isomers: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
            for candidate in isomer_rows:
                key = (
                    candidate.get("Result_Spectrum"),
                    candidate.get("Precursor_mz"),
                    candidate.get("Retention_Time"),
                )
                grouped_isomers.setdefault(key, []).append(candidate)
            for group_index, candidates in enumerate(
                grouped_isomers.values(), 1
            ):
                first = candidates[0]
                lead = next(
                    (
                        candidate for candidate in candidates
                        if candidate.get("Is_Ranking_Lead")
                    ),
                    first,
                )
                elements.append(Paragraph(
                    f"<b>Group {group_index}:</b> "
                    f"{_esc(lead.get('Compound_Name') or 'Unnamed candidate')} "
                    f"(ranking lead); {int(first.get('Candidate_Count') or len(candidates))} candidates; "
                    f"MS2 {_esc(first.get('Result_Spectrum') or 'not recorded')}; "
                    f"m/z {_esc(first.get('Precursor_mz') or '—')}; "
                    f"RT {_esc(first.get('Retention_Time') or '—')} min",
                    styles["Normal"],
                ))
                candidate_table_data = [[
                    Paragraph(_esc(label), result_header)
                    for label in [
                        "Rank / role", "Candidate", "Formula / InChIKey",
                        "Cosine", "Final / rank", "Adduct",
                        "Relationship", "Library source / spectrum",
                    ]
                ]]
                for candidate in candidates:
                    candidate_table_data.append([
                        Paragraph(
                            f"{_esc(candidate.get('Candidate_Rank') or '—')}<br/>"
                            f"{_esc(candidate.get('Candidate_Role') or '—')}",
                            result_cell,
                        ),
                        Paragraph(
                            f"{_esc(candidate.get('Compound_Name') or '—')}<br/>"
                            f"<font color='#667085'>{_esc(candidate.get('Candidate_Record_Status') or '')}</font>",
                            result_cell,
                        ),
                        Paragraph(
                            f"{_esc(candidate.get('Formula') or '—')}<br/>"
                            f"{_esc(candidate.get('InChIKey') or '—')}",
                            result_cell,
                        ),
                        Paragraph(_esc(candidate.get("Cosine") or "—"), result_number),
                        Paragraph(
                            f"{_esc(candidate.get('Final_Score') or '—')} / "
                            f"{_esc(candidate.get('Rank_Score') or '—')}",
                            result_number,
                        ),
                        Paragraph(_esc(candidate.get("Adduct") or "—"), result_cell),
                        Paragraph(_esc(candidate.get("Relationship") or "—"), result_cell),
                        Paragraph(
                            f"{_esc(candidate.get('Dataset_Source') or '—')}<br/>"
                            f"{_esc(candidate.get('Library_Spectrum') or '—')}",
                            result_cell,
                        ),
                    ])
                candidate_table = Table(
                    candidate_table_data,
                    colWidths=[
                        1.8 * cm, 2.8 * cm, 3.0 * cm, 1.1 * cm,
                        1.4 * cm, 1.2 * cm, 2.4 * cm, 2.7 * cm,
                    ],
                    repeatRows=1,
                )
                candidate_table.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#173F5F")),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("FONTNAME", (0, 0), (-1, 0), bold_font),
                    ("FONTNAME", (0, 1), (-1, -1), regular_font),
                    ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#D0D5DD")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F8FAFC")]),
                    ("LEFTPADDING", (0, 0), (-1, -1), 2),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 2),
                    ("TOPPADDING", (0, 0), (-1, -1), 2.5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
                ]))
                elements.append(candidate_table)
                elements.append(Spacer(1, 0.2 * cm))

        detail_count = 0
        for row in matched:
            if detail_count >= 30:
                break
            match = row["match"]
            if match.get("confidence") not in ("High", "Medium"):
                continue
            detail_count += 1
            if detail_count == 1:
                elements.extend(
                    [
                        Spacer(1, 0.35 * cm),
                        Paragraph("Structured evidence details", styles["Heading2"]),
                    ]
                )
            name = _esc(
                _standard_compound_name(
                    match.get("compound_name", ""), match
                )[:80]
            )
            locator = raw_spectrum_locator(row, filename)
            rt_text = (
                f"{locator['retention_time_min']:.3f} min "
                f"({locator['retention_time_clock']})"
                if locator["retention_time_min"] is not None
                else "not recorded"
            )
            elements.append(
                Paragraph(
                    f"<b>{name}</b> — raw MS2 locator: RT {rt_text}; "
                    f"precursor m/z {float(row.get('precursor_mz') or 0):.5f}; "
                    f"spectrum/native ID {_esc(locator['spectrum_id'] or 'not recorded')}; cosine "
                    f"{float(match.get('cosine') or 0):.3f}; final "
                    f"{float(match.get('final_score') or 0):.3f}",
                    styles["Normal"],
                )
            )
            elements.append(
                Paragraph(
                    "<b>Classification evidence:</b> "
                    f"{_esc(str(match.get('classification_source') or 'unknown').replace('_', ' '))}; "
                    f"confidence "
                    f"{float(match.get('classification_confidence') or 0):.3f}.",
                    small,
                )
            )

            try:
                from web.report_common import _parse_pathway_json

                pathway = _parse_pathway_json(match.get("pathway_json")) or {}
            except Exception:
                pathway = {}
            diagnostic_items = []
            for item in (pathway.get("diagnostic_ions") or [])[:6]:
                key = str(item.get("diag_key") or "").replace("_", " ").title()
                diagnostic_items.append(
                    f"{float(item.get('mz') or 0):.3f}"
                    + (f" ({_esc(key)})" if key else "")
                )
            loss_items = []
            for item in (pathway.get("neutral_losses") or [])[:6]:
                annotation = str(
                    item.get("nl_key") or item.get("annotation") or ""
                ).replace("_", " ").title()
                loss_items.append(
                    f"{float(item.get('loss_da') or 0):.3f} Da"
                    + (f" ({_esc(annotation)})" if annotation else "")
                )
            if diagnostic_items or loss_items:
                elements.append(
                    Paragraph(
                        "<b>Empirical library fragmentation signature:</b> "
                        f"diagnostic products {', '.join(diagnostic_items) or 'none'}; "
                        f"neutral losses {', '.join(loss_items) or 'none'}; "
                        f"aggregated from {_esc(pathway.get('n_spectra', 'unknown'))} "
                        "reference spectrum/s.",
                        small,
                    )
                )

            isotope = match.get("formula_isotope_evidence") or {}
            if isotope:
                elements.append(
                    Paragraph(
                        "<b>Adjacent-MS1 isotope evidence:</b> "
                        f"{_esc(isotope.get('status', 'unknown'))}; consistency "
                        f"{float(isotope.get('score') or 0):.3f}. This is a bounded "
                        "formula-level tie-breaker, not a unique-structure identifier.",
                        small,
                    )
                )
            query_adduct = str(match.get("query_adduct") or "").strip()
            reliable = bool(match.get("query_adduct_reliable"))
            elements.append(
                Paragraph(
                    "<b>Adduct evidence:</b> query "
                    f"{_esc(query_adduct or 'not assigned')}; library "
                    f"{_esc(match.get('adduct') or 'not recorded')}; "
                    f"mass basis {_esc(str(match.get('mass_match_basis') or 'precursor m/z').replace('_', ' '))}. "
                    + (
                        "The query adduct came from explicit metadata."
                        if reliable
                        else "A specific query adduct was not inferred from one precursor m/z."
                    ),
                    small,
                )
            )
            elements.append(Spacer(1, 0.12 * cm))

        if len(matched) > 200:
            elements.append(
                Paragraph(
                    f"The table displays the first 200 of {len(matched)} unique matches.",
                    styles["Normal"],
                )
            )

    doc.build(elements)
    return buffer.getvalue()


def _esc(s: str) -> str:
    """HTML escape."""
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


_BATCH_FORMAT_LABELS = {
    "csv": "CSV 完整数据包（MS2 鉴定结果 + 靶向汇总 + 原始峰审计）",
    "csv_ms2": "CSV MS2 鉴定结果表（含异构体候选）",
    "csv_isomer_candidates": "CSV 异构体竞争候选明细",
    "csv_targeted": "CSV 靶向化合物汇总",
    "csv_targeted_features": "CSV 靶向峰审计",
    "csv_audit": "CSV MS2 详细审计",
    "excel": "Excel 分层工作簿",
    "report": "HTML 详细综合报告（兼容）",
    "html_standard": "HTML 标准综合报告",
    "html_detailed": "HTML 详细综合报告",
    "html_targeted": "HTML 靶向库专用报告",
    "pdf": "PDF 综合报告",
}

_BATCH_FORMAT_TOKENS = {
    "csv": "CSV",
    "csv_ms2": "CSV_MS2",
    "csv_isomer_candidates": "CSV_ISOMER_CANDIDATES",
    "csv_targeted": "CSV_TARGETED",
    "csv_targeted_features": "CSV_TARGETED_AUDIT",
    "csv_audit": "CSV_MS2_AUDIT",
    "excel": "EXCEL_LAYERED",
    "report": "HTML_DETAILED",
    "html_standard": "HTML_STANDARD",
    "html_detailed": "HTML_DETAILED",
    "html_targeted": "HTML_TARGETED",
    "pdf": "PDF",
}


def _safe_export_component(value: Any, fallback: str = "item", max_length: int = 120) -> str:
    """Return a cross-platform-safe filename/path component."""
    text = str(value or "").strip()
    text = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", text)
    text = re.sub(r"\s+", "_", text).strip(" ._-")
    return (text[:max_length].rstrip(" ._-") or fallback)


def _extract_export_project(filename: Any) -> str:
    """Mirror the frontend project grouping rule and return a safe project name."""
    raw = Path(str(filename or "")).name
    name = Path(raw).stem
    name = re.sub(r"^\d{8,14}[-_]", "", name)
    # F编号: 开头的 F+1~6位数字 (F002 / F004021); 短编号只认开头, 避免误吞。
    fcode = re.match(r"^(F\d{1,6})", name, re.IGNORECASE)
    if fcode:
        return fcode.group(1).upper()
    # 兼容旧约定: 6位长F编号出现在名称中间。
    fcode_mid = re.search(r"(F\d{6})", name, re.IGNORECASE)
    if fcode_mid:
        return fcode_mid.group(1).upper()
    prefix = re.match(r"^([A-Za-z]+)", name)
    if prefix:
        return prefix.group(1).upper()
    return _safe_export_component(name, "Other", 60)


def _batch_zip_filename(
    sample_filenames: list[str],
    format_name: str,
    exported: int,
    date_stamp: Optional[str] = None,
    language: str = "zh",
) -> str:
    """Build a deterministic ZIP name from the exported sample filenames."""
    unique_samples = list(dict.fromkeys(
        _safe_export_component(
            Path(str(filename)).stem, "sample", 72
        )
        for filename in sample_filenames if str(filename or "").strip()
    ))
    if not unique_samples:
        sample_part = "MSsystem"
    elif len(unique_samples) <= 4:
        sample_part = "__".join(unique_samples)
    else:
        sample_part = (
            f"{unique_samples[0]}__{unique_samples[1]}"
            f"__等{len(unique_samples)}个样本"
        )
    sample_part = _safe_export_component(sample_part, "MSsystem", 170)
    token = _BATCH_FORMAT_TOKENS.get(
        format_name,
        _safe_export_component(format_name, "DATA").upper(),
    )
    if language == "en" and (
        format_name in {"report", "pdf"}
        or format_name.startswith("html_")
    ):
        token = f"{token}_EN"
    stamp = date_stamp or datetime.now().astimezone().strftime("%Y%m%d")
    count_token = f"_{exported}samples" if exported > 1 else ""
    return f"{sample_part}_{token}{count_token}_{stamp}.zip"


def _unique_zip_member(path: str, task_id: str, used: set[str]) -> str:
    """Prevent duplicate source filenames from silently overwriting ZIP members."""
    candidate = path
    if candidate in used:
        p = PurePosixPath(path)
        candidate = str(p.with_name(f"{p.stem}__{task_id[:8]}{p.suffix}"))
        serial = 2
        while candidate in used:
            candidate = str(p.with_name(f"{p.stem}__{task_id[:8]}_{serial}{p.suffix}"))
            serial += 1
    used.add(candidate)
    return candidate


@app.post("/api/export/batch")
def batch_export_query(format: str = Query(..., description="分层批量导出格式"),
                       task_ids: list[str] = None,
                       language: str = Query("zh", enum=["zh", "en"])):
    """兼容旧版/浏览器缓存版前端以 ?format=xxx 查询参数调用批量导出。

    新版前端用 /api/export/batch/{format} (路径参数); 二者等价, 均委托给 batch_export。
    背景: 旧前端发 /api/export/batch?format=csv 时, URL 路径为 /api/export/batch
    (无 {format} 段), 匹配不到 /api/export/batch/{format} 路由 → 404(路由不存在,
    非"无匹配结果")。此路由补上查询参数入口, 无需用户清浏览器缓存。
    """
    return batch_export(format, task_ids or [], language)


from web.export_jobs import (
    ExportJobs, canonical_sql_signature, rebind_request as _rebind_export_request,
    report_progress as _export_progress,
)
from web.lossless_peak_export import attachment_disposition

_export_jobs = ExportJobs(OUTPUT_DIR / "export_jobs")


def _export_source_identity(task_ids, request_parameters=None):
    from runtime.verified_sql_retrieval import sql_metadata_dependency_signature
    from web.export_identity import build_export_identity
    artifacts = {}
    task_bindings = {}
    def derived_path(builder, fallback):
        try:
            return builder()
        except ValueError:
            # A source-guard mismatch is still refused by the renderer.  The
            # identity endpoint remains inspectable and binds the missing
            # fallback instead of concealing a changed raw result.
            return fallback

    for tid in sorted(task_ids):
        with _tasks_lock:
            task = dict(_tasks.get(tid) or {})
        artifacts[f"{tid}:results"] = TASK_DIR / f"{tid}_results.json"
        artifacts[f"{tid}:compounds"] = derived_path(
            lambda: _compound_conclusion_artifact_path(TASK_DIR, tid),
            TASK_DIR / f"{tid}_v11_compounds.json",
        )
        artifacts[f"{tid}:primary"] = derived_path(
            lambda: _flagship_artifact_path(TASK_DIR, tid),
            TASK_DIR / f"{tid}_v16_flagship.json",
        )
        artifacts[f"{tid}:claims"] = _claim_contract_path(TASK_DIR, tid)
        artifacts[f"{tid}:identification"] = _identification_contract_path(TASK_DIR, tid)
        artifacts[f"{tid}:expected"] = _expected_compounds_artifact_path(tid)
        artifacts[f"{tid}:reference_qc"] = _reference_standard_qc_artifact_path(TASK_DIR, tid)
        workflow = task.get("untargeted_workflow")
        feature_path = str(workflow.get("manifest") or "").strip() if isinstance(workflow, dict) else ""
        if feature_path:
            artifacts[f"{tid}:feature_manifest"] = Path(feature_path).expanduser().resolve()
        task_bindings[tid] = {
            key: task.get(key) for key in (
                "filename", "ion_mode", "total_spectra", "matched_spectra",
                "match_rate", "params", "analysis_warnings",
                "expected_compounds_summary", "compound_conclusion_summary",
                "feature_summary", "untargeted_workflow",
            )
        }
    for relative in (
        "web/server.py", "web/export_jobs.py", "web/export_identity.py",
        "web/lossless_peak_export.py", "web/report_common.py",
        "web/plant_candidate_export.py", "web/governance_blocked_projection.py",
        "scripts/render_task_excel.py", "config/v15_flagship_policy.json",
    ):
        artifacts[f"exporter:{relative}"] = PROJECT_ROOT / relative
    library_identity = canonical_sql_signature(
        sql_metadata_dependency_signature(), PROJECT_ROOT
    )
    return build_export_identity(artifacts, bindings={
        "runtime_snapshot_id": _current_runtime_snapshot_id(),
        "runtime_source_tree_sha256": RUNTIME_SOURCE_TREE_SHA256,
        "library_identity": library_identity,
        "tasks": task_bindings,
        "request": dict(request_parameters or {}),
    })


def _export_source_signature(task_ids, request_parameters=None):
    return _export_source_identity(task_ids, request_parameters)["sha256"]


@app.post("/api/export/jobs/{format}", status_code=202)
def create_export_job(format: str, task_ids: list[str], language: str = Query("zh", enum=["zh", "en"])):
    if format not in _BATCH_FORMAT_LABELS or not task_ids:
        raise HTTPException(status_code=400, detail="请选择有效的导出格式和样本")
    task_ids = list(dict.fromkeys(task_ids))
    with _tasks_lock:
        if any(tid not in _tasks for tid in task_ids):
            raise HTTPException(status_code=404, detail="选中的分析任务不存在")
    parameters = {"format": format, "language": language, "mode": "batch"}
    identity = _export_source_identity(task_ids, parameters)

    def render():
        if format != "excel":
            return batch_export(format, task_ids, language)

        def locked_excel_batch():
            import resource
            preparation_started = time.perf_counter()
            preparation_cpu_started = time.process_time()
            with _tasks_lock:
                selected = {tid: _tasks[tid] for tid in task_ids}
            for tid, task in selected.items():
                _compound_scope_summary_for_export(tid, task)
            usage = resource.getrusage(resource.RUSAGE_SELF)
            preparation_seconds = time.perf_counter() - preparation_started
            _export_progress(phase="validating_sources", export_profile={
                "phase": "validating_sources",
                "wall_seconds": preparation_seconds,
                "process_cpu_seconds": time.process_time() - preparation_cpu_started,
                "max_rss_kib": int(usage.ru_maxrss // 1024
                                   if sys.platform == "darwin" else usage.ru_maxrss),
            })
            final_identity = _export_source_identity(task_ids, parameters)
            final_request = {
                **parameters, "task_ids": task_ids,
                "source_signature": final_identity["sha256"],
                "_identity_wall_seconds": final_identity["identity_wall_seconds"],
                "_preparation_wall_seconds": preparation_seconds,
            }
            _rebind_export_request(final_request)
            response = batch_export(format, task_ids, language)
            after = _export_source_identity(task_ids, parameters)
            if after["sha256"] != final_identity["sha256"]:
                raise RuntimeError("Excel export inputs changed while the workbook was rendering")
            return response

        return _serialized_report_export(locked_excel_batch)()

    return _export_jobs.submit(
        {"mode": "batch", "format": format, "task_ids": task_ids, "language": language,
         "source_signature": identity["sha256"],
         "_identity_wall_seconds": identity["identity_wall_seconds"]},
        render,
    )


@app.post("/api/tasks/{task_id}/export/jobs/{format}", status_code=202)
def create_single_export_job(
    task_id: str, format: Literal["csv", "excel", "html", "pdf", "raw-peaks", "multilibrary", "structure"],
    format_type: Literal["results", "report", "standard", "detailed", "targeted"] = "results",
    result_layer: Literal["ms2", "isomer_candidates", "targeted", "targeted_features", "reference_qc", "multilibrary", "plant_candidates", "external_review"] = "ms2",
    ms2_schema: Literal["reader_v3", "complete_v2"] = "reader_v3",
    blank_view: Literal["all", "clean", "flagged"] = "all",
    language: Literal["zh", "en"] = "zh",
    spectrum_index: int = Query(0, ge=0), image_format: Literal["png", "svg"] = "png",
):
    with _tasks_lock:
        if task_id not in _tasks:
            raise HTTPException(404, "分析任务不存在")
        task = _tasks[task_id]
    if format == "html" and format_type == "results":
        format_type = "standard"
    if (format == "csv" and format_type not in {"results", "report"}) or (format == "html" and format_type not in {"standard", "detailed", "targeted"}):
        raise HTTPException(400, "该文件类型不支持所选报告内容")
    def render():
        _export_progress(current_task=task_id, phase="reading_results")
        if format == "csv":
            return download_csv(task_id, format_type, result_layer, ms2_schema, blank_view)
        if format == "excel":
            def locked_excel():
                import resource
                preparation_started = time.perf_counter()
                preparation_cpu_started = time.process_time()
                _compound_scope_summary_for_export(task_id, task)
                usage = resource.getrusage(resource.RUSAGE_SELF)
                preparation_seconds = time.perf_counter() - preparation_started
                _export_progress(phase="validating_sources", export_profile={
                    "phase": "validating_sources",
                    "wall_seconds": preparation_seconds,
                    "process_cpu_seconds": time.process_time() - preparation_cpu_started,
                    "max_rss_kib": int(usage.ru_maxrss // 1024
                                       if sys.platform == "darwin" else usage.ru_maxrss),
                })
                final_identity = _export_source_identity([task_id], parameters)
                final_request = {
                    **parameters, "task_ids": [task_id],
                    "source_signature": final_identity["sha256"],
                    "_identity_wall_seconds": final_identity["identity_wall_seconds"],
                    "_preparation_wall_seconds": preparation_seconds,
                }
                _rebind_export_request(final_request)
                response = _download_excel_impl(
                    task_id, evidence_prevalidated=True
                )
                after = _export_source_identity([task_id], parameters)
                if after["sha256"] != final_identity["sha256"]:
                    if getattr(response, "path", None):
                        Path(response.path).unlink(missing_ok=True)
                        response.background = None
                    raise RuntimeError(
                        "Excel export inputs changed while the workbook was rendering"
                    )
                return response

            return _serialized_report_export(locked_excel)()
        if format == "html":
            return download_html_report(task_id, format_type, language)
        if format == "pdf":
            return download_pdf(task_id, language)
        if format == "raw-peaks":
            return download_raw_peaks(task_id)
        if format == "multilibrary":
            return get_task_flagship_multilibrary(task_id, 100000, 0, True)
        async def structure():
            response = await download_structure_image(task_id, spectrum_index, image_format)
            actual = "svg" if "svg" in response.headers.get("content-type", "") else "png"
            response.headers["Content-Disposition"] = f'attachment; filename="{task_id}_structure_{spectrum_index + 1}.{actual}"'
            return response
        return structure()
    parameters = {"mode": "single", "format": format,
        "format_type": format_type, "result_layer": result_layer, "ms2_schema": ms2_schema,
        "blank_view": blank_view, "language": language, "spectrum_index": spectrum_index,
        "image_format": image_format}
    identity = _export_source_identity([task_id], parameters)
    return _export_jobs.submit({**parameters, "task_ids": [task_id],
        "source_signature": identity["sha256"],
        "_identity_wall_seconds": identity["identity_wall_seconds"]}, render)


@app.get("/api/export/jobs")
def active_export_jobs():
    with _export_jobs.lock:
        return {"jobs": [state for job_id in list(_export_jobs.owned)
                         if (state := _export_jobs.get(job_id))["status"] in {"queued", "running"}]}


@app.get("/api/export/jobs/{job_id}")
def get_export_job(job_id: str):
    try:
        state = _export_jobs.get(job_id)
        if state["status"] == "completed" and state["request"].get("source_signature"):
            try:
                request_parameters = {
                    key: value for key, value in state["request"].items()
                    if key not in {"task_ids", "source_signature"}
                    and not str(key).startswith("_")
                }
                current_signature = _export_source_signature(
                    state["request"]["task_ids"], request_parameters
                )
            except (OSError, ValueError, TypeError, KeyError):
                current_signature = None
            if state["request"]["source_signature"] != current_signature:
                # An immutable historical file remains downloadable. It must
                # never be reused as a newly generated current-version export.
                return {**state, "stale": True, "can_reuse": False,
                        "notice": "历史导出文件；重新导出会使用当前输入与系统版本。"}

        return state
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/api/export/jobs/{job_id}/download")
def download_export_job(job_id: str):
    try:
        state = get_export_job(job_id)
        if state["status"] != "completed":
            raise HTTPException(status_code=409, detail=state.get("error", "导出尚未完成"))
        archive = _export_jobs.archive_path(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    headers = {key: value for key, value in state["headers"].items()
               if key.lower() in {"content-disposition", "x-mssystem-batch-format",
                                  "x-mssystem-exported-tasks", "x-mssystem-failed-tasks"}}
    headers["X-MSsystem-Export-Stale"] = str(bool(state.get("stale"))).lower()
    return FileResponse(archive, media_type=state.get("media_type", "application/zip"), headers=headers)


@app.post("/api/export/batch/{format}")
@_serialized_report_export
def batch_export(
    format: str,
    task_ids: list[str],
    language: str = Query("zh", enum=["zh", "en"]),
):
    """Batch export: package multiple task results into a ZIP.
    format: explicit CSV layer, layered Excel, HTML report type, or PDF
    task_ids: JSON array of task ID strings in request body.
    """
    if format not in _BATCH_FORMAT_LABELS:
        raise HTTPException(status_code=400, detail=f"不支持的批量导出格式: {format}")
    if not task_ids:
        raise HTTPException(status_code=400, detail="未选择任何任务")

    from runtime.verified_sql_retrieval import sql_metadata_repair_status, sql_metadata_repair_preflight
    sql_status = sql_metadata_repair_status()
    if (sql_status.get("status") == "not_initialized"
            and any((TASK_DIR / f"{tid}_results.json").is_file() for tid in task_ids)):
        _export_progress(phase="verifying_library")
        verification = sql_metadata_repair_preflight()
        if verification.get("compatible") is not True:
            raise HTTPException(status_code=503, detail=(
                "谱库核验未通过，尚未开始逐样本导出。原因：" + str(verification.get("error", "unknown"))))
    if sql_status.get("status") == "invalid_binding":
        raise HTTPException(status_code=503, detail=(
            "谱库证据校验失效，导出已停止；须重启服务并通过完整谱库核验后重试。原因："
            + str(sql_status.get("error", ""))))

    # 项目分组不再只存在于前端。后端按同一规则识别项目，保证下载名和ZIP目录可信。
    with _tasks_lock:
        selected_project_by_id = {
            tid: _extract_export_project((_tasks.get(tid) or {}).get("filename"))
            for tid in task_ids if tid in _tasks
        }
    selected_projects = list(dict.fromkeys(selected_project_by_id.values()))
    multi_project = len(selected_projects) > 1

    # Batch ZIPs can contain dozens of large result tables.  A BytesIO plus
    # ``getvalue()`` held two complete copies in API memory and the frontend
    # then created a third Blob copy.  Spool after 8 MiB so large archives are
    # assembled on disk and streamed in bounded chunks.
    export_tmp_root = OUTPUT_DIR / "export_tmp"
    export_tmp_root.mkdir(parents=True, exist_ok=True)
    buf = tempfile.SpooledTemporaryFile(
        max_size=8 * 1024 * 1024,
        mode="w+b",
        prefix="batch-export-",
        suffix=".zip",
        dir=export_tmp_root,
    )
    exported = 0
    exported_projects: list[str] = []
    exported_sample_filenames: list[str] = []
    manifest_rows: list[list[Any]] = []
    used_members: set[str] = set()
    failed = []  # 记录失败/跳过的任务, 便于诊断
    exported_at = datetime.now().astimezone().isoformat(timespec="seconds")

    failure_details = []
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for position, tid in enumerate(task_ids):
            _export_progress(phase="reading_results", completed=position, current_task=tid)
            # 异常隔离: 单个 task 数据异常或导出失败不影响其他 task, 仅记录跳过。
            try:
                task, results = _load_task_results_data(tid)
                if not task:
                    logger.warning("batch_export: task %s not found, skipped", tid)
                    failed.append(tid)
                    continue
                results = results or []

                total_spectra = len(results)
                matched_spectra = sum(
                    1 for row in results if row.get("match")
                )
                confidence_summary = _compound_confidence_summary(
                    [row for row in results if row.get("match")]
                )
                _export_progress(phase="checking_evidence")
                compound_summary = _compound_scope_summary_for_export(
                    tid, task
                )
                compound_conclusions = _compound_conclusions_for_export(
                    tid, task
                )
                matched = _dedup_matched([r for r in results if r.get("match")])
                expected_payload = _load_expected_compounds_payload(
                    tid, task, tolerate_errors=True
                )
                reference_qc_payload = _load_reference_standard_qc_payload(
                    tid, tolerate_errors=True
                )
                expected_rows = expected_feature_rows(
                    expected_payload, language=language
                )
                if not matched and not expected_rows:
                    continue

                source_filename = str(task.get("filename") or tid or "task")
                project = _extract_export_project(source_filename)
                safe_name = _safe_export_component(Path(source_filename).stem, tid, 140)
                folder = f"{project}/" if multi_project else ""
                archive_name = ""

                _export_progress(phase="writing_csv" if format.startswith("csv") else "writing_report")
                if format == "csv_ms2":
                    ms2_name = _unique_zip_member(
                        f"{folder}{safe_name}_ms2_identifications.csv",
                        f"{tid}-ms2",
                        used_members,
                    )
                    zf.writestr(
                        ms2_name,
                        _ms2_csv_text(
                            matched,
                            language=language,
                            compounds=compound_conclusions,
                            task_id=tid, results=results,
                        ).encode("utf-8-sig"),
                    )
                    archive_name = ms2_name

                elif format == "csv_isomer_candidates":
                    archive_name = _unique_zip_member(
                        f"{folder}{safe_name}_isomer_competing_candidates.csv",
                        f"{tid}-isomer-candidates",
                        used_members,
                    )
                    zf.writestr(
                        archive_name,
                        _isomer_candidate_csv_text(
                            matched,
                            language=language,
                            compounds=compound_conclusions,
                        ).encode("utf-8-sig"),
                    )

                elif format == "csv_targeted":
                    archive_name = _unique_zip_member(
                        f"{folder}{safe_name}_targeted_compounds.csv",
                        f"{tid}-targeted",
                        used_members,
                    )
                    zf.writestr(
                        archive_name,
                        _targeted_csv_text(
                            expected_payload, language=language
                        ).encode("utf-8-sig"),
                    )

                elif format == "csv_targeted_features":
                    archive_name = _unique_zip_member(
                        f"{folder}{safe_name}_targeted_peak_audit.csv",
                        f"{tid}-targeted-audit",
                        used_members,
                    )
                    zf.writestr(
                        archive_name,
                        _targeted_feature_audit_csv_text(
                            expected_payload, language=language
                        ).encode("utf-8-sig"),
                    )

                elif format == "csv_audit":
                    from web.report_common import (
                        REPORT_COLUMN_NAMES,
                        row_to_dict,
                    )

                    audit_buffer = io.StringIO()
                    audit_writer = csv.writer(audit_buffer)
                    audit_writer.writerow(REPORT_COLUMN_NAMES)
                    for result in matched:
                        audit_row = row_to_dict(result, report_mode=True)
                        audit_writer.writerow([
                            audit_row.get(column, "")
                            for column in REPORT_COLUMN_NAMES
                        ])
                    archive_name = _unique_zip_member(
                        f"{folder}{safe_name}_ms2_detailed_audit.csv",
                        f"{tid}-ms2-audit",
                        used_members,
                    )
                    from web.governance_blocked_projection import (
                        append_governance_blocked_csv_section,
                        project_governance_blocked_candidates,
                    )
                    audit_text = append_governance_blocked_csv_section(
                        audit_buffer.getvalue(),
                        project_governance_blocked_candidates(results),
                        language=language,
                    )
                    zf.writestr(
                        archive_name,
                        audit_text.encode("utf-8-sig"),
                    )

                elif format == "csv":
                    ms2_name = _unique_zip_member(
                        f"{folder}{safe_name}_ms2.csv",
                        f"{tid}-ms2",
                        used_members,
                    )
                    targeted_name = _unique_zip_member(
                        f"{folder}{safe_name}_targeted_ms1.csv",
                        f"{tid}-targeted",
                        used_members,
                    )
                    targeted_audit_name = _unique_zip_member(
                        f"{folder}{safe_name}_targeted_peak_audit.csv",
                        f"{tid}-targeted-audit",
                        used_members,
                    )
                    reference_qc_name = ""
                    zf.writestr(
                        ms2_name,
                        _ms2_csv_text(
                            matched,
                            language=language,
                            compounds=compound_conclusions,
                            task_id=tid, results=results,
                        ).encode("utf-8-sig"),
                    )
                    zf.writestr(
                        targeted_name,
                        _targeted_csv_text(
                            expected_payload, language=language
                        ).encode("utf-8-sig"),
                    )
                    zf.writestr(
                        targeted_audit_name,
                        _targeted_feature_audit_csv_text(
                            expected_payload, language=language
                        ).encode("utf-8-sig"),
                    )
                    reference_qc_text = reference_standard_qc_csv_text(
                        reference_qc_payload, language=language
                    )
                    if reference_qc_text:
                        reference_qc_name = _unique_zip_member(
                            f"{folder}{safe_name}_reference_standard_qc.csv",
                            f"{tid}-reference-qc",
                            used_members,
                        )
                        zf.writestr(
                            reference_qc_name,
                            reference_qc_text.encode("utf-8-sig"),
                        )
                    archive_members = [
                        ms2_name,
                        targeted_name,
                        targeted_audit_name,
                    ]
                    if reference_qc_name:
                        archive_members.append(reference_qc_name)
                    archive_name = " | ".join(archive_members)

                elif format == "excel":
                    excel_bytes = _build_analysis_excel_bytes(
                        tid, task, results, language=language
                    )
                    archive_name = _unique_zip_member(
                        f"{folder}{safe_name}_results.xlsx", tid, used_members
                    )
                    zf.writestr(archive_name, excel_bytes)

                elif format == "html_standard":
                    if language == "en":
                        from web.report_en import generate_english_html

                        html_content = generate_english_html(
                            results,
                            task.get("filename", "Analysis"),
                            detailed=False,
                            total_spectra=total_spectra,
                            matched_spectra=matched_spectra,
                            high_supported_compounds=confidence_summary.get(
                                "high_supported_unique_compounds", 0
                            ),
                            high_confidence_scope_groups=compound_summary.get(
                                "high_confidence_scope_group_count"
                            ),
                            high_confidence_exact_structures=compound_summary.get(
                                "high_confidence_2d_connectivity_count",
                                compound_summary.get(
                                    "high_confidence_detected_count",
                                    compound_summary.get("high_confidence_exact_structure_count"),
                                ),
                            ),
                            high_confidence_isomer_groups=compound_summary.get(
                                "high_confidence_isomer_group_count"
                            ),
                        )
                        html_content = inject_expected_html(
                            html_content,
                            expected_payload,
                            language="en",
                            analysis_warnings=(
                                task.get("analysis_warnings") or []
                            ),
                        )
                    else:
                        html_content = _build_html_report(
                            task, results, compound_summary
                        )
                    html_content = inject_reference_standard_qc_html(
                        html_content,
                        reference_qc_payload,
                        language=language,
                    )
                    html_content = inject_confidence_guide_html(
                        html_content, language=language
                    )
                    from web.governance_blocked_projection import (
                        inject_governance_blocked_html,
                        project_governance_blocked_candidates,
                    )
                    html_content = inject_governance_blocked_html(
                        html_content,
                        project_governance_blocked_candidates(results),
                        language=language,
                    )
                    suffix = "_en" if language == "en" else ""
                    archive_name = _unique_zip_member(
                        f"{folder}{safe_name}_standard_report{suffix}.html",
                        tid,
                        used_members,
                    )
                    zf.writestr(archive_name, html_content)

                elif format in {"report", "html_detailed"}:
                    if language == "en":
                        from web.report_en import generate_english_html
                        html_content = generate_english_html(
                            results,
                            task.get("filename", "Analysis"),
                            detailed=True,
                            total_spectra=total_spectra,
                            matched_spectra=matched_spectra,
                            high_supported_compounds=confidence_summary.get(
                                "high_supported_unique_compounds", 0
                            ),
                            high_confidence_scope_groups=compound_summary.get(
                                "high_confidence_scope_group_count"
                            ),
                            high_confidence_exact_structures=compound_summary.get(
                                "high_confidence_2d_connectivity_count",
                                compound_summary.get(
                                    "high_confidence_detected_count",
                                    compound_summary.get("high_confidence_exact_structure_count"),
                                ),
                            ),
                            high_confidence_isomer_groups=compound_summary.get(
                                "high_confidence_isomer_group_count"
                            ),
                        )
                        from web.governance_blocked_projection import (
                            inject_governance_blocked_html,
                            project_governance_blocked_candidates,
                        )
                        html_content = inject_governance_blocked_html(
                            html_content,
                            project_governance_blocked_candidates(results),
                            language="en",
                        )
                        html_content = inject_expected_html(
                            html_content,
                            expected_payload,
                            language="en",
                            analysis_warnings=(
                                task.get("analysis_warnings") or []
                            ),
                        )
                    else:
                        # 详细报告：支持搜索、排序、点击展开匹配详情。
                        # ``matched`` 已按与单报告/PDF/Excel 相同的完整身份键
                        # 去重。旧批量分支虽然算出了它，却误遍历原始
                        # ``results``，导致同一化合物因多张 MS2 重复成多行。
                        report_results = matched + [
                            row for row in results if not row.get("match")
                        ]
                        html_results = []
                        for r in report_results:
                            m = r.get("match")
                            if m:
                                m["lib_peaks"] = r.get("lib_peaks", [])
                                ed = m.get("expert_details") or {}
                                if ed:
                                    for ek in ["neutral_loss_score", "diagnostic_ion_score", "rule_coverage_score",
                                               "fg_signature_score", "structure_fragment_score", "bond_energy_score",
                                               "feature_library_score", "isomer_differentiation_score", "chemical_reasoning"]:
                                        if ek in ed and ek not in m:
                                            m[ek] = ed[ek]
                            html_results.append({
                                "precursor_mz": r["precursor_mz"],
                                "spectrum_id": r.get("spectrum_id", ""),
                                "peaks_count": r.get("peaks_count", 0),
                                "retention_time": r.get("retention_time"),
                                "retention_time_unit": r.get(
                                    "retention_time_unit"
                                ),
                                "precursor_scan_ref": r.get(
                                    "precursor_scan_ref"
                                ),
                                "collision_energy": r.get(
                                    "collision_energy"
                                ),
                                "collision_energy_unit": r.get(
                                    "collision_energy_unit"
                                ),
                                "ion_mode": r.get("ion_mode"),
                                "precursor_charge": r.get(
                                    "precursor_charge"
                                ),
                                "tic": r.get("tic"),
                                "base_peak": r.get("base_peak"),
                                "query_peaks": r.get("query_peaks", []),
                                "raw_peaks": r.get("raw_peaks", []),
                                "raw_peaks_input_count": r.get(
                                    "raw_peaks_input_count"
                                ),
                                "raw_peaks_saved_count": r.get(
                                    "raw_peaks_saved_count"
                                ),
                                "raw_peaks_complete": r.get(
                                    "raw_peaks_complete"
                                ),
                                "raw_peaks_storage_policy": r.get(
                                    "raw_peaks_storage_policy"
                                ),
                                "peaks": r.get("query_peaks", []),
                                "matches": [m] if m else [],
                                "best_cosine": m["cosine"] if m else 0.0,
                                "display_class": r.get("display_class"),
                                "dl_confidence": r.get("dl_confidence"),
                            })
                        from web.app import _gen_rich
                        from web.governance_blocked_projection import project_governance_blocked_candidates
                        html_content = _gen_rich(
                            html_results,
                            governance_blocked_candidates=project_governance_blocked_candidates(results),
                            source_name=task.get("filename", "Analysis"),
                            total_spectra=total_spectra,
                            matched_spectra=matched_spectra,
                            high_supported_compounds=confidence_summary.get(
                                "high_supported_unique_compounds", 0
                            ),
                            high_confidence_scope_groups=compound_summary.get(
                                "high_confidence_scope_group_count"
                            ),
                            high_confidence_exact_structures=compound_summary.get(
                                "high_confidence_2d_connectivity_count",
                                compound_summary.get(
                                    "high_confidence_detected_count",
                                    compound_summary.get("high_confidence_exact_structure_count"),
                                ),
                            ),
                            high_confidence_isomer_groups=compound_summary.get(
                                "high_confidence_isomer_group_count"
                            ),
                        )
                        html_content = inject_expected_html(
                            html_content,
                            expected_payload,
                            language="zh",
                            analysis_warnings=(
                                task.get("analysis_warnings") or []
                            ),
                        )
                    html_content = inject_reference_standard_qc_html(
                        html_content,
                        reference_qc_payload,
                        language=language,
                    )
                    html_content = inject_confidence_guide_html(
                        html_content, language=language
                    )
                    suffix = "_en" if language == "en" else ""
                    archive_name = _unique_zip_member(
                        f"{folder}{safe_name}_detailed_report{suffix}.html", tid, used_members
                    )
                    zf.writestr(archive_name, html_content)

                elif format == "html_targeted":
                    html_content = render_expected_html_report(
                        expected_payload,
                        task.get("filename", "Analysis"),
                        language=language,
                        analysis_warnings=(
                            task.get("analysis_warnings") or []
                        ),
                    )
                    html_content = inject_reference_standard_qc_html(
                        html_content,
                        reference_qc_payload,
                        language=language,
                    )
                    html_content = inject_confidence_guide_html(
                        html_content, language=language
                    )
                    suffix = "_en" if language == "en" else ""
                    archive_name = _unique_zip_member(
                        f"{folder}{safe_name}_targeted_library_report{suffix}.html",
                        tid,
                        used_members,
                    )
                    zf.writestr(archive_name, html_content)

                elif format == "pdf":
                    try:
                        pdf_bytes = _build_pdf_report(
                            task,
                            results,
                            language=language,
                            compound_summary=compound_summary,
                            reference_standard_qc_payload=reference_qc_payload,
                        )
                        suffix = "_en" if language == "en" else ""
                        archive_name = _unique_zip_member(
                            f"{folder}{safe_name}_report{suffix}.pdf", tid, used_members
                        )
                        zf.writestr(archive_name, pdf_bytes)
                    except Exception as e:
                        raise RuntimeError(f"PDF 生成失败：{e}") from e

                if format in {"csv", "csv_ms2", "csv_audit"}:
                    _export_progress(phase="packing_raw_peaks")
                    from web.lossless_peak_export import iter_raw_peak_jsonl_gzip
                    raw_name = _unique_zip_member(
                        f"{folder}{safe_name}_raw_peaks.jsonl.gz",
                        f"{tid}-raw-peaks", used_members,
                    )
                    raw_path = TASK_DIR / f"{tid}_results.json"
                    raw_source_sha = None
                    if raw_path.is_file():
                        import hashlib
                        with raw_path.open("rb") as raw_source:
                            raw_source_sha = hashlib.file_digest(raw_source, "sha256").hexdigest()
                    # gzip is already compressed; recompressing it in ZIP
                    # adds CPU time without changing the exported evidence.
                    raw_info = zipfile.ZipInfo(raw_name)
                    raw_info.compress_type = zipfile.ZIP_STORED
                    with zf.open(raw_info, "w") as raw_member:
                        for block in iter_raw_peak_jsonl_gzip(results, tid, raw_source_sha):
                            raw_member.write(block)
                    archive_name += f" | {raw_name}"

                exported += 1
                if project not in exported_projects:
                    exported_projects.append(project)
                exported_sample_filenames.append(source_filename)
                matched_spectra = task.get("matched_spectra")
                if matched_spectra is None:
                    matched_spectra = sum(1 for r in results if r.get("match"))
                manifest_rows.append([
                    project,
                    source_filename,
                    tid,
                    archive_name,
                    _BATCH_FORMAT_LABELS[format] + (
                        " English"
                        if language == "en" and (
                            format in {"report", "pdf"}
                            or format.startswith("html_")
                        )
                        else ""
                    ),
                    matched_spectra,
                    len(matched),
                    task.get("match_rate", ""),
                    json.dumps(task.get("params") or {}, ensure_ascii=False, sort_keys=True),
                    exported_at,
                ])
            except Exception as e:
                logger.error("batch_export: task %s failed (skipped, others continue): %s", tid, e, exc_info=True)
                failed.append(tid)
                failure_details.append({"task_id": tid, "error": str(getattr(e, "detail", e))})
                continue
            finally:
                _export_progress(completed=position + 1, exported=exported, failed=len(failed))

        if failed:
            zf.writestr("00_export_failures.json", json.dumps(
                {"failed_task_ids": failed, "errors": failure_details}, ensure_ascii=False, indent=2))
        if manifest_rows:
            guide_rows = confidence_guide_table_rows(language)
            if guide_rows:
                guide_buffer = io.StringIO()
                guide_fields = [*guide_rows[0].keys(), "Policy_ID"]
                guide_writer = csv.DictWriter(
                    guide_buffer, fieldnames=guide_fields
                )
                guide_writer.writeheader()
                for guide_row in guide_rows:
                    guide_writer.writerow({
                        **guide_row,
                        "Policy_ID": CONFIDENCE_GUIDE_POLICY_ID,
                    })
                zf.writestr(
                    "00_confidence_and_score_guide.csv",
                    guide_buffer.getvalue().encode("utf-8-sig"),
                )
            manifest_buf = io.StringIO()
            manifest_writer = csv.writer(manifest_buf)
            manifest_writer.writerow([
                "Project", "Source_File", "Task_ID", "Output_File", "Format",
                "Matched_Spectra", "Unique_Compounds", "Match_Rate_%",
                "Analysis_Parameters", "Exported_At",
            ])
            manifest_writer.writerows(manifest_rows)
            zf.writestr("00_export_manifest.csv", manifest_buf.getvalue().encode("utf-8-sig"))

    if exported == 0:
        buf.close()
        if failed:
            reason = failure_details[0]["error"] if failure_details else "分析任务不存在或无法读取"
            raise HTTPException(status_code=500, detail=f"导出失败：{len(failed)} 个任务未能生成文件。{reason}")
        raise HTTPException(status_code=404, detail="所选任务均无可导出的匹配结果")

    buf.seek(0, os.SEEK_END)
    archive_size = buf.tell()
    buf.seek(0)
    zip_name = _batch_zip_filename(
        exported_sample_filenames, format, exported, language=language
    )
    ascii_zip_name = re.sub(r"[^A-Za-z0-9._-]+", "_", zip_name).strip("._")
    if not ascii_zip_name.lower().endswith(".zip"):
        ascii_zip_name = f"MSsystem_{format}_{exported}samples.zip"
    from fastapi.responses import StreamingResponse
    from urllib.parse import quote

    def stream_archive():
        try:
            while True:
                chunk = buf.read(1024 * 1024)
                if not chunk:
                    break
                yield chunk
        finally:
            buf.close()

    return StreamingResponse(
        stream_archive(),
        media_type="application/zip",
        headers={
            "Content-Disposition": (
                f"attachment; filename=\"{ascii_zip_name}\"; "
                f"filename*=UTF-8''{quote(zip_name, safe='')}"
            ),
            "Cache-Control": "no-store, max-age=0",
            "X-MSsystem-Batch-Format": format,
            "X-MSsystem-Exported-Tasks": str(exported),
            "X-MSsystem-Failed-Tasks": str(len(failed)),
            "Content-Length": str(archive_size),
        },
    )


# ── PDF Export ──────────────────────────────────────────────────────

@app.get("/api/tasks/{task_id}/download/pdf")
@_serialized_report_export
def download_pdf(
    task_id: str,
    language: str = Query("zh", enum=["zh", "en"]),
):
    """Generate and download a PDF report (reportlab, no system deps)."""
    task, results = _load_task_results_data(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    matched = [r for r in results if r.get("match")]
    expected_payload = _load_expected_compounds_payload(
        task_id, task, tolerate_errors=True
    )
    if not matched and not expected_feature_rows(expected_payload):
        raise HTTPException(
            status_code=404,
            detail="No MS2 matches or Expected Compounds results to export",
        )

    compound_summary = _compound_scope_summary_for_export(task_id, task)
    reference_qc_payload = _load_reference_standard_qc_payload(
        task_id, tolerate_errors=True
    )
    try:
        pdf_bytes = _build_pdf_report(
            task,
            results,
            language=language,
            compound_summary=compound_summary,
            reference_standard_qc_payload=reference_qc_payload,
        )
    except Exception as e:
        logger.error("PDF generation failed for task %s: %s", task_id, e)
        raise HTTPException(status_code=500, detail=f"PDF generation failed: {e}")

    suffix = "_en" if language == "en" else ""
    filename = f"{Path(task['filename']).stem}_report{suffix}.pdf"
    from fastapi.responses import Response
    return Response(
        content=pdf_bytes,
        media_type="application/pdf",
        headers={"Content-Disposition": attachment_disposition(filename)},
    )


# ═════════════════════════════════════════════════════════════════════
# Molecular Networking (GNPS-style Feature-Based Molecular Networking)
# ═════════════════════════════════════════════════════════════════════

def _build_network_for_task(task_id: str,
                            cos_thresh: float = 0.7, min_peaks: int = 6,
                            topk: int = 10, mz_tol: float = 0.1,
                            precursor_window: float = 200.0,
                            max_component_size: int = 0,
                            top_n: int = 70, rel_threshold: float = 0.01,
                            topk_policy: str = "mutual"):
    """Load one completed task and build its molecular network."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    results = _get_task_results(task_id)
    if not results:
        raise HTTPException(status_code=404, detail="No results to network")

    from matching.molecular_network import build_network
    network = build_network(
        results,
        cos_thresh=float(cos_thresh),
        min_peaks=int(min_peaks),
        topk=int(topk),
        mz_tol=float(mz_tol),
        precursor_window=float(precursor_window),
        max_component_size=int(max_component_size),
        top_n=int(top_n),
        rel_threshold=float(rel_threshold),
        topk_policy=str(topk_policy),
    )
    return task, results, network


@app.get("/api/tasks/{task_id}/network")
def get_task_network(
    task_id: str,
    cos_thresh: float = Query(0.7, ge=0.0, le=1.0, description="Edge cosine threshold"),
    min_peaks: int = Query(6, ge=1, le=50, description="Minimum matched fragment peaks"),
    topk: int = Query(10, ge=0, le=50, description="Top neighbours retained per node"),
    topk_policy: Literal["either", "mutual", "off"] = Query(
        "mutual", description="Top-K edge retention policy"
    ),
    mz_tol: float = Query(0.1, ge=0.01, le=2.0, description="Fragment bin width (Da)"),
    precursor_window: float = Query(200.0, ge=0.0, le=2000.0, description="Precursor m/z window (Da)"),
    max_component_size: int = Query(0, ge=0, le=5000, description="Maximum molecular-family size; 0 disables artificial splitting"),
):
    """Build a GNPS-style molecular network for a completed task."""
    _, _, network = _build_network_for_task(
        task_id,
        cos_thresh=cos_thresh,
        min_peaks=min_peaks,
        topk=topk,
        mz_tol=mz_tol,
        precursor_window=precursor_window,
        max_component_size=max_component_size,
        topk_policy=topk_policy,
    )
    from matching.molecular_network import edge_public, node_public
    return {
        "task_id": task_id,
        "nodes": [node_public(node) for node in network["nodes"]],
        "edges": [
            edge_public(edge)
            for edge in (network.get("edge_details") or network["edges"])
        ],
        "families": network["families"],
        "stats": network["stats"],
        "params": network["params"],
    }


@app.get("/api/tasks/{task_id}/network/node/{node_id}")
def get_task_network_node(
    task_id: str,
    node_id: int,
    cos_thresh: float = Query(0.7, ge=0.0, le=1.0),
    min_peaks: int = Query(6, ge=1, le=50),
    topk: int = Query(10, ge=0, le=50),
    topk_policy: Literal["either", "mutual", "off"] = Query("mutual"),
    mz_tol: float = Query(0.1, ge=0.01, le=2.0),
    max_component_size: int = Query(0, ge=0, le=5000),
):
    """Return the exact measured spectrum selected for one network node."""
    with _tasks_lock:
        task = _tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    results = _get_task_results(task_id)
    from matching.evidence_network import prepare_feature_nodes

    nodes, _ = prepare_feature_nodes(
        results,
        top_n=70,
        rel_threshold=0.01,
        consensus_tolerance_da=max(0.01, min(float(mz_tol), 0.10)),
    )
    node = next(
        (
            candidate
            for candidate in nodes
            if int(candidate.get("idx", -1)) == node_id
        ),
        None,
    )
    if node is None or node_id < 0 or node_id >= len(results):
        raise HTTPException(status_code=404, detail="Network node not found")
    result = results[node_id]
    peak_evidence = node.get("peak_evidence") or {}
    return {
        "idx": node_id,
        "peaks": node.get("peaks") or [],
        "peak_source": node.get("peak_source"),
        "peak_evidence": node.get("peak_evidence"),
        "annotation_result_index": peak_evidence.get(
            "annotation_result_index"
        ),
        "annotation_spectrum_id": node.get("spectrum_id"),
        "feature_source": node.get("feature_source"),
        "feature_id": node.get("feature_id"),
        "spectrum_count": node.get("spectrum_count"),
        "precursor_mz": peak_evidence.get("precursor_mz", node.get("precursor_mz")),
        "retention_time": peak_evidence.get(
            "retention_time", node.get("retention_time")
        ),
        "feature_precursor_mz": node.get("precursor_mz"),
        "feature_retention_time": node.get("retention_time"),
        "spectrum_id": peak_evidence.get("spectrum_id") or node.get("spectrum_id"),
        "lcms_consensus_spectrum_id": (
            peak_evidence.get("artifact_id")
            or node.get("lcms_consensus_spectrum_id")
        ),
        "collision_energy": peak_evidence.get(
            "collision_energy", node.get("collision_energy")
        ),
        "collision_energy_unit": peak_evidence.get("collision_energy_unit"),
        "collision_energy_mode": peak_evidence.get("collision_energy_mode"),
        "activation_method": peak_evidence.get("activation_method"),
        "instrument_configuration_id": peak_evidence.get(
            "instrument_configuration_id"
        ),
        "match": node.get("match") or result.get("match"),
        "in_network": True,
    }


@app.get("/api/tasks/{task_id}/network/graphml")
def download_network_graphml(
    task_id: str,
    cos_thresh: float = Query(0.7, ge=0.0, le=1.0),
    min_peaks: int = Query(6, ge=1, le=50),
    topk: int = Query(10, ge=0, le=50),
    topk_policy: Literal["either", "mutual", "off"] = Query("mutual"),
    mz_tol: float = Query(0.1, ge=0.01, le=2.0),
    precursor_window: float = Query(200.0, ge=0.0, le=2000.0),
    max_component_size: int = Query(0, ge=0, le=5000),
):
    """Download the current network as Cytoscape-compatible GraphML."""
    task, _, network = _build_network_for_task(
        task_id,
        cos_thresh=cos_thresh,
        min_peaks=min_peaks,
        topk=topk,
        mz_tol=mz_tol,
        precursor_window=precursor_window,
        max_component_size=max_component_size,
        topk_policy=topk_policy,
    )
    from matching.molecular_network import to_graphml
    graphml = to_graphml(network)
    if not graphml:
        raise HTTPException(status_code=500, detail="GraphML export unavailable")
    filename = f"{Path(task.get('filename', task_id)).stem}_network.graphml"
    from fastapi.responses import Response
    return Response(
        content=graphml.encode("utf-8"),
        media_type="application/graphml+xml",
        headers={"Content-Disposition": attachment_disposition(filename)},
    )


@app.get("/api/tasks/{task_id}/network/nodes.csv")
def download_network_nodes_csv(
    task_id: str,
    cos_thresh: float = Query(0.7, ge=0.0, le=1.0),
    min_peaks: int = Query(6, ge=1, le=50),
    topk: int = Query(10, ge=0, le=50),
    topk_policy: Literal["either", "mutual", "off"] = Query("mutual"),
    mz_tol: float = Query(0.1, ge=0.01, le=2.0),
    precursor_window: float = Query(200.0, ge=0.0, le=2000.0),
    max_component_size: int = Query(0, ge=0, le=5000),
):
    """Download the exact displayed network's node table."""
    task, _, network = _build_network_for_task(
        task_id,
        cos_thresh=cos_thresh,
        min_peaks=min_peaks,
        topk=topk,
        mz_tol=mz_tol,
        precursor_window=precursor_window,
        max_component_size=max_component_size,
        topk_policy=topk_policy,
    )
    from matching.molecular_network import to_nodes_csv
    from fastapi.responses import Response
    filename = f"{Path(task.get('filename', task_id)).stem}_network_nodes.csv"
    return Response(
        content=to_nodes_csv(network).encode("utf-8-sig"),
        media_type="text/csv; charset=utf-8",
        headers={"Content-Disposition": attachment_disposition(filename)},
    )


@app.get("/api/tasks/{task_id}/network/edges.csv")
def download_network_edges_csv(
    task_id: str,
    cos_thresh: float = Query(0.7, ge=0.0, le=1.0),
    min_peaks: int = Query(6, ge=1, le=50),
    topk: int = Query(10, ge=0, le=50),
    topk_policy: Literal["either", "mutual", "off"] = Query("mutual"),
    mz_tol: float = Query(0.1, ge=0.01, le=2.0),
    precursor_window: float = Query(200.0, ge=0.0, le=2000.0),
    max_component_size: int = Query(0, ge=0, le=5000),
):
    """Download the exact displayed network's edge table."""
    task, _, network = _build_network_for_task(
        task_id,
        cos_thresh=cos_thresh,
        min_peaks=min_peaks,
        topk=topk,
        mz_tol=mz_tol,
        precursor_window=precursor_window,
        max_component_size=max_component_size,
        topk_policy=topk_policy,
    )
    from matching.molecular_network import to_edges_csv
    from fastapi.responses import Response
    filename = f"{Path(task.get('filename', task_id)).stem}_network_edges.csv"
    return Response(
        content=to_edges_csv(network).encode("utf-8-sig"),
        media_type="text/csv; charset=utf-8",
        headers={"Content-Disposition": attachment_disposition(filename)},
    )


# ═════════════════════════════════════════════════════════════════════
# Serve Static Frontend
# ═════════════════════════════════════════════════════════════════════

from web.v16_candidate_comparison import (
    CandidateComparisonService,
    register_candidate_comparison_routes,
)


def _candidate_comparison_task_snapshot(task_id):
    with _tasks_lock:
        task = _tasks.get(task_id)
        return dict(task) if task else None


_candidate_comparison_service = CandidateComparisonService(
    task_directory=lambda: TASK_DIR,
    task_lookup=_candidate_comparison_task_snapshot,
    artifact_path=lambda directory, task_id: _flagship_artifact_path(directory, task_id),
    artifact_current=lambda artifact: _flagship_artifact_is_current(artifact),
    read_entity=lambda artifact, entity_id: _read_flagship_artifact(
        artifact, entity_id=entity_id, offset=0, limit=1),
    retriever_factory=lambda: _v15_ce_evidence_retriever(),
    runtime_snapshot=lambda: _current_runtime_snapshot_id(),
)
register_candidate_comparison_routes(app, _candidate_comparison_service)


from web.v16_cross_occurrence_api import (
    CrossOccurrenceService,
    register_cross_occurrence_routes,
)


def _cross_occurrence_task_snapshot(task_id):
    with _tasks_lock:
        task = _tasks.get(task_id)
        return dict(task) if task else None


def _cross_occurrence_runtime_idle():
    queue = _full_queue_snapshot()
    if queue.get("running") or queue.get("queued"):
        return False
    with _tasks_lock:
        return not any(
            str(task.get("status") or "") in {"queued", "running", "finalizing"}
            for task in _tasks.values()
        )


def _cross_occurrence_ensure_artifact(task_id, task, source, artifact):
    metadata = _materialize_current_v15_flagship(
        task_id,
        task,
        Path(source),
        Path(artifact),
    )
    with _tasks_lock:
        current = _tasks.get(task_id)
        if current is not None:
            current["flagship_summary"] = metadata
    _save_tasks()
    return metadata


_cross_occurrence_service = CrossOccurrenceService(
    task_directory=lambda: TASK_DIR,
    task_lookup=_cross_occurrence_task_snapshot,
    artifact_path=lambda directory, task_id: _flagship_artifact_path(directory, task_id),
    artifact_current=lambda artifact: _flagship_artifact_is_current(artifact),
    retriever_factory=lambda: _v15_ce_evidence_retriever(),
    runtime_snapshot=lambda: _current_runtime_snapshot_id(),
    runtime_idle=_cross_occurrence_runtime_idle,
    ensure_artifact=_cross_occurrence_ensure_artifact,
)
register_cross_occurrence_routes(app, _cross_occurrence_service)


from web.completed_batches import register_completed_batch_routes

_completed_batch_service = register_completed_batch_routes(
    app,
    batch_root=PROJECT_ROOT.parent / "research" / "20261021_xbt_ms2_integrated",
    cache_root=Path("/Volumes/LaCie/Mssystem_V16/completed_batch_exports"),
    remote_host="msqc-server",
)

from web.frozen_batch_tasks import FrozenBatchTasks, install_frozen_task_access

_frozen_task_bridge = (
    FrozenBatchTasks(_completed_batch_service) if _completed_batch_service else None
)
install_frozen_task_access(
    app, _frozen_task_bridge, lambda: _tasks,
    _RESULT_SUMMARY_ROW_FIELDS, _RESULT_SUMMARY_MATCH_FIELDS,
)

from web.frozen_native_delivery import register_frozen_native_delivery_routes

_frozen_native_delivery = register_frozen_native_delivery_routes(app, _frozen_task_bridge)


STATIC_DIR = Path(__file__).resolve().parent / "static"
STATIC_DIR.mkdir(parents=True, exist_ok=True)


@app.get("/")
async def root():
    """Serve the main web client."""
    index_path = STATIC_DIR / "index.html"
    if index_path.exists():
        return HTMLResponse(
            content=index_path.read_text(encoding="utf-8"),
            headers={
                "Cache-Control": "no-store, max-age=0",
                "Pragma": "no-cache",
            },
        )
    # Fallback: return a simple redirect message
    return HTMLResponse(content="""
    <html><body style="font-family:sans-serif;text-align:center;padding:3rem;">
    <h1>MSsystem Server Running</h1>
    <p>Frontend not found. Place index.html in web/static/</p>
    <p>API docs: <a href="/docs">/docs</a></p>
    </body></html>
    """)


def _storage_tasks_snapshot():
    with _tasks_lock:
        return {task_id: dict(task) for task_id, task in _tasks.items()}


from web.storage_management import install_storage_management

storage_maintenance_manager = install_storage_management(
    app,
    output_dir=OUTPUT_DIR,
    tasks_snapshot=_storage_tasks_snapshot,
    recycle_inventory_unlocked=_task_recycle_inventory_unlocked,
    recycle_lock=_task_recycle_lock,
    cache_lock=_report_export_lock,
    append_audit=_append_task_deletion_audit,
)


# Mount static files (for any additional assets)
if STATIC_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


# ═════════════════════════════════════════════════════════════════════
# Main Entry Point
# ═════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import uvicorn
    logger.info("Starting MSsystem Web Server...")
    uvicorn.run(
        "web.server:app",
        host=os.getenv("MSSYSTEM_HOST", "127.0.0.1"),
        port=int(os.getenv("MSSYSTEM_PORT", str(DEFAULT_PORT))),
        # 与 run.sh 一致；热重载会丢失内存中的任务/谱库状态。
        reload=False,
        log_level="info",
    )
