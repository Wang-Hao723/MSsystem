#!/usr/bin/env python3
"""Upgrade completed standard runs to the current V15 reporting scope.

The migration reuses the saved V15 ranked shortlist and immutable scan-level
results.  It does not repeat library retrieval or alter candidate order,
candidate score, scientific database names, formulae, structures, or source
accessions.  Each rewritten artifact is written atomically by the production
artifact writers, and a machine-readable before/after receipt is persisted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.compound_conclusions import (  # noqa: E402
    read_compound_conclusion_artifact,
    write_compound_conclusion_artifact,
)
from pipeline.v15_flagship import (  # noqa: E402
    upgrade_v15_flagship_identity_scope_artifact,
)
from web import server  # noqa: E402


DEFAULT_RECEIPT = (
    ROOT
    / "docs/audits/v15_production_identity_scope_migration_20260829"
    / "receipt.json"
)


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    )
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _ranking_projection(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for entity in payload.get("entities") or []:
        if not isinstance(entity, Mapping):
            continue
        selected = entity.get("winner") or entity.get(
            "leading_candidate_when_abstained"
        ) or {}
        rows.append(
            {
                "entity_id": entity.get("entity_id"),
                "ranking_lead": selected.get("identity_key"),
                "candidates": [
                    {
                        "identity_key": candidate.get("identity_key"),
                        "entity_rank": candidate.get("entity_rank"),
                        "entity_score": candidate.get("entity_score"),
                    }
                    for candidate in entity.get("candidates") or []
                    if isinstance(candidate, Mapping)
                ],
            }
        )
    return rows


def _scientific_identity_projection(
    payload: Mapping[str, Any],
) -> list[dict[str, Any]]:
    rows = []
    for entity in payload.get("entities") or []:
        if not isinstance(entity, Mapping):
            continue
        for candidate in entity.get("candidates") or []:
            if not isinstance(candidate, Mapping):
                continue
            display = candidate.get("candidate") or {}
            rows.append(
                {
                    "entity_id": entity.get("entity_id"),
                    "identity_key": candidate.get("identity_key"),
                    "compound_name": display.get("compound_name"),
                    "formula": display.get("formula"),
                    "inchikey": display.get("inchikey"),
                    "smiles": display.get("smiles"),
                    "dataset_source": display.get("dataset_source"),
                    "lib_spectrum_id": display.get("lib_spectrum_id"),
                    "lib_source_accession": display.get(
                        "lib_source_accession"
                    ),
                }
            )
    return rows


def _projection_hash(payload: Mapping[str, Any], field: str) -> str:
    projection = (
        _ranking_projection(payload)
        if field == "ranking"
        else _scientific_identity_projection(payload)
    )
    return hashlib.sha256(_canonical(projection)).hexdigest()


def _compound_counts(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"status": "missing"}
    metadata = read_compound_conclusion_artifact(
        path, offset=0, limit=1
    ).get("metadata") or {}
    return {
        "status": "available",
        "algorithm_id": metadata.get("algorithm_id"),
        "unique_report_rows": int(metadata.get("unique_report_rows") or 0),
        "main_result_count": int(metadata.get("main_result_count") or 0),
        "isomer_group_detected_count": int(
            metadata.get("isomer_group_detected_count") or 0
        ),
        "unresolved_precursor_event_count": int(
            metadata.get("unresolved_precursor_event_count") or 0
        ),
    }


def migrate(receipt_path: Path) -> dict[str, Any]:
    server._load_tasks()
    task_rows = [
        (task_id, dict(task))
        for task_id, task in sorted(server._tasks.items())
        if task.get("status") == "completed"
    ]
    forbidden = [
        {"task_id": task_id, "filename": task.get("filename")}
        for task_id, task in task_rows
        if str(task.get("filename") or "").lower().startswith("automs")
    ]
    if forbidden:
        raise RuntimeError(
            "AutoMS runs are forbidden in this production migration: "
            + json.dumps(forbidden, ensure_ascii=False)
        )

    receipt: dict[str, Any] = {
        "schema_version": "mssystem.v15.production-identity-scope-migration/1.0.0",
        "started_at": datetime.now(timezone.utc).isoformat(),
        "project_root": str(ROOT),
        "completed_task_count": len(task_rows),
        "failed_task_count_not_migrated": sum(
            task.get("status") == "failed" for task in server._tasks.values()
        ),
        "automs_task_count": 0,
        "library_retrieval_repeated": False,
        "candidate_ranking_changed": False,
        "candidate_scores_changed": False,
        "database_scientific_identity_fields_changed": False,
        "tasks": [],
        "failures": [],
    }
    _atomic_json(receipt_path, receipt)

    for position, (task_id, task_snapshot) in enumerate(task_rows, start=1):
        filename = str(task_snapshot.get("filename") or task_id)
        flagship_path = server._flagship_artifact_path(
            server.TASK_DIR, task_id
        )
        result_path = server.TASK_DIR / f"{task_id}_results.json"
        compound_path = server._compound_conclusion_artifact_path(
            server.TASK_DIR, task_id
        )
        row: dict[str, Any] = {
            "task_id": task_id,
            "filename": filename,
            "position": position,
            "status": "started",
        }
        try:
            before_flagship = json.loads(
                flagship_path.read_text(encoding="utf-8")
            )
            before_metadata = dict(before_flagship.get("metadata") or {})
            before_resolution_algorithm = (
                ((before_flagship.get("policy") or {}).get(
                    "identity_resolution"
                ) or {}).get("algorithm_id")
            )
            before_ranking = _projection_hash(before_flagship, "ranking")
            before_identity = _projection_hash(
                before_flagship, "scientific_identity"
            )
            before_compounds = _compound_counts(compound_path)
            before_flagship_sha = _sha256(flagship_path)
            del before_flagship

            flagship_summary = upgrade_v15_flagship_identity_scope_artifact(
                flagship_path,
                task_id=task_id,
                runtime_snapshot_id=server._current_runtime_snapshot_id(),
            )
            after_flagship = json.loads(
                flagship_path.read_text(encoding="utf-8")
            )
            after_ranking = _projection_hash(after_flagship, "ranking")
            after_identity = _projection_hash(
                after_flagship, "scientific_identity"
            )
            if before_ranking != after_ranking:
                raise RuntimeError("candidate order, lead, rank, or score changed")
            if before_identity != after_identity:
                raise RuntimeError(
                    "database name/formula/structure/provenance fields changed"
                )

            prepared_rows = server._get_task_results(task_id)
            compound_summary = write_compound_conclusion_artifact(
                task_id,
                prepared_rows,
                compound_path,
                claim_decisions_path=server._claim_contract_path(
                    server.TASK_DIR, task_id
                ),
                source_artifact=result_path,
                primary_identity_artifact=flagship_path,
            )
            del prepared_rows
            after_compounds = _compound_counts(compound_path)
            invariant = (
                read_compound_conclusion_artifact(
                    compound_path, offset=0, limit=1
                ).get("metadata")
                or {}
            ).get("precursor_event_main_result_invariant") or {}
            if invariant.get("satisfied") is not True:
                raise RuntimeError(
                    "one-primary-result-per-precursor-event invariant failed"
                )

            with server._tasks_lock:
                current = server._tasks.get(task_id)
                if current is None:
                    raise RuntimeError("task disappeared from the registry")
                current["flagship_summary"] = flagship_summary
                current["compound_conclusion_summary"] = compound_summary
                server._apply_final_compound_count_summary(
                    current, compound_summary
                )
            if not server._save_tasks():
                raise RuntimeError("task index persistence failed")

            row.update(
                {
                    "status": "completed",
                    "source_flagship_algorithm_id": before_metadata.get(
                        "algorithm_id"
                    ),
                    "target_flagship_algorithm_id": (
                        (after_flagship.get("metadata") or {}).get(
                            "algorithm_id"
                        )
                    ),
                    "source_identity_resolution_algorithm_id": (
                        before_resolution_algorithm
                    ),
                    "target_identity_resolution_algorithm_id": (
                        ((after_flagship.get("policy") or {}).get(
                            "identity_resolution"
                        ) or {}).get("algorithm_id")
                    ),
                    "before_ranking_sha256": before_ranking,
                    "after_ranking_sha256": after_ranking,
                    "ranking_unchanged": True,
                    "before_scientific_identity_sha256": before_identity,
                    "after_scientific_identity_sha256": after_identity,
                    "database_scientific_identity_fields_unchanged": True,
                    "before_flagship_sha256": before_flagship_sha,
                    "after_flagship_sha256": _sha256(flagship_path),
                    "before_decision_scope_counts": before_metadata.get(
                        "decision_scope_counts"
                    ) or {},
                    "after_decision_scope_counts": (
                        (after_flagship.get("metadata") or {}).get(
                            "decision_scope_counts"
                        )
                        or {}
                    ),
                    "before_compound_counts": before_compounds,
                    "after_compound_counts": after_compounds,
                    "main_result_count_delta": int(
                        after_compounds.get("main_result_count") or 0
                    )
                    - int(before_compounds.get("main_result_count") or 0),
                    "unresolved_event_count_delta": int(
                        after_compounds.get(
                            "unresolved_precursor_event_count"
                        )
                        or 0
                    )
                    - int(
                        before_compounds.get(
                            "unresolved_precursor_event_count"
                        )
                        or 0
                    ),
                    "isomer_group_count_delta": int(
                        after_compounds.get("isomer_group_detected_count") or 0
                    )
                    - int(
                        before_compounds.get("isomer_group_detected_count") or 0
                    ),
                    "precursor_event_main_result_invariant": invariant,
                }
            )
            del after_flagship
        except Exception as exc:
            row.update(
                {
                    "status": "failed",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
            receipt["failures"].append(
                {
                    "task_id": task_id,
                    "filename": filename,
                    "error": row["error"],
                }
            )
        receipt["tasks"].append(row)
        receipt["last_updated_at"] = datetime.now(timezone.utc).isoformat()
        _atomic_json(receipt_path, receipt)
        print(
            f"[{position}/{len(task_rows)}] {task_id} {filename}: "
            f"{row['status']}",
            flush=True,
        )

    completed = [row for row in receipt["tasks"] if row["status"] == "completed"]
    receipt.update(
        {
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "migrated_task_count": len(completed),
            "migration_failure_count": len(receipt["failures"]),
            "tasks_with_changed_decision_scope_counts": sum(
                row.get("before_decision_scope_counts")
                != row.get("after_decision_scope_counts")
                for row in completed
            ),
            "tasks_with_changed_main_result_count": sum(
                int(row.get("main_result_count_delta") or 0) != 0
                for row in completed
            ),
            "tasks_with_changed_unresolved_event_count": sum(
                int(row.get("unresolved_event_count_delta") or 0) != 0
                for row in completed
            ),
            "tasks_with_changed_isomer_group_count": sum(
                int(row.get("isomer_group_count_delta") or 0) != 0
                for row in completed
            ),
            "aggregate_main_result_count_before": sum(
                int((row.get("before_compound_counts") or {}).get(
                    "main_result_count"
                ) or 0)
                for row in completed
            ),
            "aggregate_main_result_count_after": sum(
                int((row.get("after_compound_counts") or {}).get(
                    "main_result_count"
                ) or 0)
                for row in completed
            ),
            "aggregate_unresolved_event_count_before": sum(
                int((row.get("before_compound_counts") or {}).get(
                    "unresolved_precursor_event_count"
                ) or 0)
                for row in completed
            ),
            "aggregate_unresolved_event_count_after": sum(
                int((row.get("after_compound_counts") or {}).get(
                    "unresolved_precursor_event_count"
                ) or 0)
                for row in completed
            ),
            "aggregate_isomer_group_count_before": sum(
                int((row.get("before_compound_counts") or {}).get(
                    "isomer_group_detected_count"
                ) or 0)
                for row in completed
            ),
            "aggregate_isomer_group_count_after": sum(
                int((row.get("after_compound_counts") or {}).get(
                    "isomer_group_detected_count"
                ) or 0)
                for row in completed
            ),
        }
    )
    _atomic_json(receipt_path, receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    receipt = migrate(args.receipt.resolve())
    print(json.dumps({
        key: receipt.get(key)
        for key in (
            "migrated_task_count",
            "migration_failure_count",
            "tasks_with_changed_decision_scope_counts",
            "tasks_with_changed_main_result_count",
            "aggregate_main_result_count_before",
            "aggregate_main_result_count_after",
            "aggregate_unresolved_event_count_before",
            "aggregate_unresolved_event_count_after",
            "aggregate_isomer_group_count_before",
            "aggregate_isomer_group_count_after",
        )
    }, ensure_ascii=False, indent=2), flush=True)
    return 0 if not receipt["failures"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
