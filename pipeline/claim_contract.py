"""Compact EviGate 2.0 decision sidecar derived from IdentificationRecord DTOs."""

from __future__ import annotations

import hashlib
import json
from runtime.json_encoding import write_json_encoded
import os
from collections import Counter
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from governance.claim_compiler import (
    EVIGATE2_POLICY_ID,
    ClaimCompilationContext,
    compile_claim_decision,
)
from schemas.decision import CalibrationStatus, ClaimScope
from schemas.identification import IdentificationRecord


CLAIM_CONTRACT_SCHEMA_VERSION = "1.0.0"


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def claim_contract_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}_claim_decisions.json"


def iter_identification_records(path: Path) -> Iterator[IdentificationRecord]:
    try:
        import ijson

        with path.open("rb") as handle:
            for payload in ijson.items(handle, "records.item", use_float=True):
                if isinstance(payload, Mapping):
                    yield IdentificationRecord.model_validate(payload)
        return
    except ImportError:
        pass
    payload = json.loads(path.read_text(encoding="utf-8"))
    for item in payload.get("records") or []:
        if isinstance(item, Mapping):
            yield IdentificationRecord.model_validate(item)


def write_claim_contract(
    task_id: str,
    records: Iterable[IdentificationRecord],
    output_path: Path,
    *,
    context: ClaimCompilationContext,
) -> dict[str, Any]:
    decisions = []
    digest = hashlib.sha256()
    status_counts: Counter[str] = Counter()
    scope_counts: Counter[str] = Counter()
    released = 0
    abstained = 0
    calibrated = 0
    ranking_changed = 0
    confidence_overwritten = 0
    for record in records:
        decision = compile_claim_decision(record, context)
        item = {
            "record_id": record.record_id,
            "spectrum_id": record.observation.spectrum_id,
            "candidate_identity_key": record.candidate.identity_key,
            "legacy_confidence": record.scores.legacy_confidence,
            "decision": decision.model_dump(mode="json", exclude_none=True),
        }
        decisions.append(item)
        digest.update(_canonical(item))
        digest.update(b"\n")
        status_counts[decision.status] += 1
        scope_counts[decision.scope.value] += 1
        released += int(decision.released)
        abstained += int(decision.scope == ClaimScope.ABSTAIN)
        calibrated += int(decision.calibration.status == CalibrationStatus.CALIBRATED)
        ranking_changed += int(decision.ranking_changed)
        confidence_overwritten += int(decision.confidence_field_overwritten)
    metadata = {
        "contract_schema_version": CLAIM_CONTRACT_SCHEMA_VERSION,
        "policy_id": EVIGATE2_POLICY_ID,
        "task_id": task_id,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "context": context.model_dump(mode="json"),
        "record_count": len(decisions),
        "released_count": released,
        "abstained_count": abstained,
        "calibrated_count": calibrated,
        "status_counts": dict(sorted(status_counts.items())),
        "scope_counts": dict(sorted(scope_counts.items())),
        "production_ranking_changed_records": ranking_changed,
        "legacy_confidence_overwritten_records": confidence_overwritten,
        "decision_semantic_sha256": digest.hexdigest(),
        "qvalues_treated_as_probabilities": False,
        "coverage_guarantee": False,
        "legacy_results_mutated": False,
        "expected_compounds_counted": False,
    }
    payload = {"metadata": metadata, "decisions": decisions}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        write_json_encoded(payload, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, output_path)
    metadata["artifact"] = str(output_path.resolve())
    metadata["size_bytes"] = output_path.stat().st_size
    return metadata


def read_claim_contract(
    path: Path,
    *,
    offset: int = 0,
    limit: int = 500,
) -> dict[str, Any]:
    if offset < 0 or limit < 1:
        raise ValueError("offset must be >=0 and limit must be >=1")
    try:
        import ijson
    except ImportError:
        payload = json.loads(path.read_text(encoding="utf-8"))
        decisions = payload.get("decisions") or []
        return {
            "metadata": payload.get("metadata") or {},
            "decisions": decisions[offset : offset + limit],
            "offset": offset,
            "limit": limit,
            "has_more": offset + limit < len(decisions),
        }
    with path.open("rb") as handle:
        metadata = next(ijson.items(handle, "metadata", use_float=True), {})
    decisions = []
    with path.open("rb") as handle:
        for index, decision in enumerate(
            ijson.items(handle, "decisions.item", use_float=True)
        ):
            if index < offset:
                continue
            if len(decisions) >= limit:
                break
            decisions.append(decision)
    total = int(metadata.get("record_count") or 0)
    return {
        "metadata": metadata,
        "decisions": decisions,
        "offset": offset,
        "limit": limit,
        "has_more": offset + len(decisions) < total,
    }
