#!/usr/bin/env python3
"""Split the audited zytw/BZL payload into two independent datasets.

This is a deterministic lineage-preserving transformation.  It does not
resolve or alter chemistry: records are copied from the already audited
combined payload, assigned a dataset-local row index, and given the stable
dataset keys requested for V10 (``user_contrib:zytw`` and
``user_contrib:bzl``).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any


OLD_DATASET_SOURCE = "user_contrib:literature_template_zytw_bzl"
SOURCE_TIER = "literature_presence_only"
DATASETS = {
    "zytw": {
        "dataset_source": "user_contrib:zytw",
        "rows": 118,
        "structures": 87,
        "formula_rows": 118,
        "plant_source": "Ligularia cymbulifera",
        "payload_name": "zytw_target_payload_v5.json",
        "audit_name": "zytw_structure_audit_v5.json",
        "csv_name": "literature_template_zytw_v5.csv",
    },
    "bzl": {
        "dataset_source": "user_contrib:bzl",
        "rows": 228,
        "structures": 118,
        "formula_rows": 227,
        "plant_source": "Scutellaria barbata",
        "payload_name": "bzl_target_payload_v3.json",
        "audit_name": "bzl_structure_audit_v3.json",
        "csv_name": "literature_template_bzl_v3.csv",
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json_atomic(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    fragment_columns = [f"fragment_mz_{index}" for index in range(1, 16)]
    columns = [
        "spectrum_name", "precursor_mz", "charge", "collision_energy",
        "retention_time", "instrument", "adduct", "component_class",
        "plant_source", "smiles", "inchi", "inchikey", "formula",
        "source_tier", "dataset_source", "structure_identity_verified",
        "structure_identity_source", "confidence", "resolve_reason",
        "source_row_index", "source_dataset", "reported_mass", "cas",
        "provenance",
    ] + fragment_columns
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for record in records:
                row = {
                    "spectrum_name": record["name"],
                    "precursor_mz": record.get("precursor_mz") or "",
                    "charge": record.get("charge") or "",
                    "collision_energy": "",
                    "retention_time": "",
                    "instrument": record.get("instrument") or "Literature-curated",
                    "adduct": record.get("adduct") or "",
                    "component_class": "",
                    "plant_source": record.get("plant_source") or "",
                    "smiles": record.get("smiles") or "",
                    "inchi": record.get("inchi") or "",
                    "inchikey": record.get("inchikey") or "",
                    "formula": record.get("formula") or "",
                    "source_tier": record.get("source_tier") or SOURCE_TIER,
                    "dataset_source": record["dataset_source"],
                    "structure_identity_verified": int(
                        bool(record.get("structure_identity_verified"))
                    ),
                    "structure_identity_source": (
                        record.get("structure_identity_source") or "none"
                    ),
                    "confidence": record.get("confidence") or "",
                    "resolve_reason": record.get("resolve_reason") or "",
                    "source_row_index": record["source_row_index"],
                    "source_dataset": record["dataset"],
                    "reported_mass": record.get("source_mass") or "",
                    "cas": record.get("cas") or "",
                    "provenance": record.get("provenance") or "",
                }
                row.update({column: "" for column in fragment_columns})
                writer.writerow(row)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _summary(records: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "source_rows": len(records),
        "resolved_total": sum(bool(row.get("smiles")) for row in records),
        "unresolved_total": sum(not row.get("smiles") for row in records),
        "formula_ready": sum(bool(row.get("formula")) for row in records),
        "identity_only": sum(not row.get("formula") for row in records),
        "identity_source_counts": dict(Counter(
            str(row.get("structure_identity_source") or "none")
            for row in records
        )),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--combined-payload", type=Path, required=True)
    parser.add_argument("--combined-audit", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()

    combined_bytes = args.combined_payload.read_bytes()
    combined_sha256 = hashlib.sha256(combined_bytes).hexdigest()
    combined = json.loads(combined_bytes)
    audit = json.loads(args.combined_audit.read_text(encoding="utf-8"))
    records = combined.get("records")
    audit_rows = audit.get("rows")
    if combined.get("schema_version") != 1:
        raise ValueError("expected combined payload schema_version=1")
    if combined.get("dataset_source") != OLD_DATASET_SOURCE:
        raise ValueError("combined payload has an unexpected dataset_source")
    if not isinstance(records, list) or len(records) != 346:
        raise ValueError("expected exactly 346 combined payload rows")
    if not isinstance(audit_rows, list) or len(audit_rows) != 346:
        raise ValueError("expected exactly 346 combined audit rows")
    if Counter(str(row.get("dataset") or "") for row in records) != Counter(
        {"zytw": 118, "bzl": 228}
    ):
        raise ValueError("unexpected dataset membership in combined payload")

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "kind": "zytw_bzl_payload_split_manifest",
        "source_payload": str(args.combined_payload.resolve()),
        "source_payload_sha256": combined_sha256,
        "source_audit": str(args.combined_audit.resolve()),
        "source_audit_sha256": _sha256(args.combined_audit),
        "old_dataset_source": OLD_DATASET_SOURCE,
        "datasets": {},
    }

    for dataset, expected in DATASETS.items():
        selected = [dict(row) for row in records if row.get("dataset") == dataset]
        selected_audit = [
            dict(row) for row in audit_rows if row.get("dataset") == dataset
        ]
        if len(selected) != expected["rows"] or len(selected_audit) != expected["rows"]:
            raise ValueError(f"{dataset} row count changed during split")
        for local_index, record in enumerate(selected):
            if record.get("plant_source") != expected["plant_source"]:
                raise ValueError(f"{dataset}[{local_index}] plant_source mismatch")
            if record.get("source_tier") != SOURCE_TIER or record.get("peaks"):
                raise ValueError(f"{dataset}[{local_index}] is not presence-only")
            record["dataset_source"] = expected["dataset_source"]
            record["source_row_index"] = local_index
            record["source_file"] = expected["csv_name"]
        for row in selected_audit:
            row["dataset_source"] = expected["dataset_source"]

        summary = _summary(selected)
        if summary["resolved_total"] != expected["structures"]:
            raise ValueError(f"{dataset} structure coverage changed during split")
        if summary["formula_ready"] != expected["formula_rows"]:
            raise ValueError(f"{dataset} formula coverage changed during split")

        payload_path = args.output_dir / expected["payload_name"]
        audit_path = args.output_dir / expected["audit_name"]
        csv_path = args.output_dir / expected["csv_name"]
        split_payload = {
            "schema_version": 1,
            "dataset_source": expected["dataset_source"],
            "source_tier": SOURCE_TIER,
            "expected_record_count": expected["rows"],
            "split_from": {
                "dataset_source": OLD_DATASET_SOURCE,
                "payload_path": str(args.combined_payload.resolve()),
                "payload_sha256": combined_sha256,
            },
            "input_sha256": combined.get("input_sha256") or {},
            "records": selected,
        }
        split_audit = {
            "schema_version": 1,
            "dataset_source": expected["dataset_source"],
            "split_from": {
                "dataset_source": OLD_DATASET_SOURCE,
                "audit_path": str(args.combined_audit.resolve()),
                "audit_sha256": _sha256(args.combined_audit),
            },
            "summary": summary,
            "rows": selected_audit,
        }
        _write_json_atomic(payload_path, split_payload)
        _write_json_atomic(audit_path, split_audit)
        _write_csv(csv_path, selected)
        manifest["datasets"][dataset] = {
            "dataset_source": expected["dataset_source"],
            "payload_path": str(payload_path.resolve()),
            "payload_sha256": _sha256(payload_path),
            "audit_path": str(audit_path.resolve()),
            "audit_sha256": _sha256(audit_path),
            "csv_path": str(csv_path.resolve()),
            "csv_sha256": _sha256(csv_path),
            **summary,
        }

    manifest_path = args.manifest or (
        args.output_dir / "zytw_bzl_split_manifest_20260806.json"
    )
    _write_json_atomic(manifest_path, manifest)
    print(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True))
    print(f"manifest={manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
