"""Join versioned contracts onto legacy rows for staged API/report migration."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict, deque
from pathlib import Path
from typing import Any


def _identity_key(match: dict[str, Any]) -> str:
    inchikey = str(match.get("inchikey") or "").strip()
    if inchikey:
        return f"inchikey:{inchikey.upper()}"
    formula = str(match.get("formula") or "unknown_formula").strip()
    name = str(match.get("compound_name") or "unknown_compound").strip()
    digest = hashlib.sha256(f"{formula}|{name}".encode("utf-8")).hexdigest()[:20]
    return f"legacy:{digest}"


def attach_claim_decisions(
    results: list[dict[str, Any]],
    claim_contract: Path,
) -> dict[str, Any]:
    """Attach EviGate 2 decisions in memory; never rewrite legacy artifacts."""
    if not claim_contract.is_file():
        return {
            "status": "unavailable",
            "attached": 0,
            "unmatched_contract_decisions": 0,
            "legacy_artifact_mutated": False,
        }
    payload = json.loads(claim_contract.read_text(encoding="utf-8"))
    queues: dict[tuple[str, str], deque] = defaultdict(deque)
    for item in payload.get("decisions") or []:
        if not isinstance(item, dict):
            continue
        key = (
            str(item.get("spectrum_id") or ""),
            str(item.get("candidate_identity_key") or ""),
        )
        queues[key].append(item)
    attached = 0
    already_attached = 0
    for row in results:
        match = row.get("match") if isinstance(row, dict) else None
        if not isinstance(match, dict):
            continue
        if isinstance(match.get("final_claim_decision"), dict):
            already_attached += 1
            continue
        key = (str(row.get("spectrum_id") or ""), _identity_key(match))
        if not queues[key]:
            continue
        item = queues[key].popleft()
        match["final_claim_decision"] = item.get("decision") or {}
        match["final_claim_record_id"] = item.get("record_id")
        match["legacy_confidence_at_claim_compilation"] = item.get(
            "legacy_confidence"
        )
        attached += 1
    remaining = sum(len(queue) for queue in queues.values())
    return {
        "status": "available",
        "policy_id": (payload.get("metadata") or {}).get("policy_id"),
        "attached": attached,
        "already_attached": already_attached,
        "unmatched_contract_decisions": remaining,
        "legacy_artifact_mutated": False,
    }
