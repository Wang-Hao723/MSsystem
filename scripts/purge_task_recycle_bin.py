#!/usr/bin/env python3
"""Audit and permanently purge task data already moved to the recycle bin.

The default mode is read-only.  ``--execute`` only removes a recycled task when
all of the following are true:

* its manifest agrees with the directory identifiers;
* its task ID is absent from the current task registry;
* the deletion audit contains the same deletion/task pair;
* every recorded recycled path is inside that task's recycle directory; and
* the recycle directory contains no symbolic links.

An execution receipt is written outside the recycle tree before any removal and
updated after the purge.  The append-only task deletion audit is retained.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


IDENTIFIER_RE = re.compile(r"[A-Za-z0-9_-]{1,80}\Z")
RECYCLE_EVENTS = {"single_task_recycled", "batch_task_recycled"}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with partial.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        partial.replace(path)
    finally:
        partial.unlink(missing_ok=True)


def _append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def _lexically_within(path: Path, parent: Path) -> bool:
    try:
        path.absolute().relative_to(parent.absolute())
        return True
    except ValueError:
        return False


def _tree_stats(root: Path) -> tuple[int, int, list[str]]:
    """Return allocated-ish bytes, entry count, and symlinks without following links."""
    total_bytes = 0
    entries = 0
    symlinks: list[str] = []
    stack = [root]
    while stack:
        current = stack.pop()
        try:
            stat = current.lstat()
        except FileNotFoundError:
            continue
        entries += 1
        if current.is_symlink():
            symlinks.append(str(current))
            continue
        total_bytes += int(stat.st_size)
        if current.is_dir():
            with os.scandir(current) as iterator:
                stack.extend(Path(item.path) for item in iterator)
    return total_bytes, entries, symlinks


def _read_audit_pairs(path: Path) -> tuple[set[tuple[str, str]], list[str]]:
    pairs: set[tuple[str, str]] = set()
    errors: list[str] = []
    if not path.is_file():
        return pairs, [f"missing deletion audit: {path}"]
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                errors.append(f"audit line {line_number}: invalid JSON ({exc})")
                continue
            if row.get("event") not in RECYCLE_EVENTS:
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
    return pairs, errors


def _task_ids_from_registry(path: Path) -> set[str]:
    payload = _load_json(path)
    if not isinstance(payload, dict):
        raise ValueError(f"task registry is not a JSON object: {path}")
    return {str(task_id) for task_id in payload}


def _inspect_task(
    task_root: Path,
    *,
    output_dir: Path,
    current_ids: set[str],
    audit_pairs: set[tuple[str, str]],
) -> dict[str, Any]:
    deletion_id = task_root.parent.name
    task_id = task_root.name
    reasons: list[str] = []
    manifest_path = task_root / "manifest.json"
    task_path = task_root / "task.json"
    manifest: dict[str, Any] = {}

    if not IDENTIFIER_RE.fullmatch(deletion_id):
        reasons.append("invalid deletion directory identifier")
    if not IDENTIFIER_RE.fullmatch(task_id):
        reasons.append("invalid task directory identifier")
    if task_root.is_symlink():
        reasons.append("task recycle root is a symbolic link")
    if not manifest_path.is_file() or manifest_path.is_symlink():
        reasons.append("manifest.json is missing or is a symbolic link")
    if not task_path.is_file() or task_path.is_symlink():
        reasons.append("task.json is missing or is a symbolic link")

    if not reasons:
        try:
            loaded = _load_json(manifest_path)
            if not isinstance(loaded, dict):
                raise ValueError("manifest is not an object")
            manifest = loaded
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            reasons.append(f"manifest cannot be read: {exc}")

    if manifest:
        if manifest.get("schema_version") != 1:
            reasons.append("unsupported manifest schema")
        if str(manifest.get("deletion_id") or "") != deletion_id:
            reasons.append("manifest deletion_id does not match directory")
        if str(manifest.get("task_id") or "") != task_id:
            reasons.append("manifest task_id does not match directory")
        if manifest.get("restored_at"):
            reasons.append("task has been restored")

    if task_id in current_ids:
        reasons.append("task ID is present in the current task registry")
    if (deletion_id, task_id) not in audit_pairs:
        reasons.append("deletion/task pair is absent from the deletion audit")

    moved_paths: list[str] = []
    files_root = task_root / "files"
    for position, item in enumerate(manifest.get("moved") or []):
        if not isinstance(item, dict):
            reasons.append(f"moved[{position}] is not an object")
            continue
        recycled_text = str(item.get("recycle_path") or "")
        original_text = str(item.get("original_path") or "")
        recycled = Path(recycled_text)
        original = Path(original_text)
        moved_paths.append(recycled_text)
        if not recycled_text or not _lexically_within(recycled, files_root):
            reasons.append(f"moved[{position}] recycled path is outside task recycle files")
        elif not recycled.exists() and not recycled.is_symlink():
            reasons.append(f"moved[{position}] recycled path is missing")
        if not original_text or not _lexically_within(original, output_dir):
            reasons.append(f"moved[{position}] original path is outside managed output")

    bytes_total, entry_count, symlinks = _tree_stats(task_root)
    if symlinks:
        reasons.append(f"recycle tree contains {len(symlinks)} symbolic link(s)")

    return {
        "deletion_id": deletion_id,
        "task_id": task_id,
        "deleted_at": manifest.get("deleted_at"),
        "recoverable": manifest.get("recoverable"),
        "moved_path_count": len(moved_paths),
        "bytes": bytes_total,
        "entry_count": entry_count,
        "eligible": not reasons,
        "reasons": reasons,
        "task_root": str(task_root),
    }


def _inventory(output_dir: Path) -> dict[str, Any]:
    registry_path = output_dir / "tasks_index.json"
    audit_path = output_dir / "runtime" / "task_deletion_audit.jsonl"
    recycle_root = output_dir / "task_recycle_bin"
    current_ids = _task_ids_from_registry(registry_path)
    audit_pairs, audit_errors = _read_audit_pairs(audit_path)
    rows: list[dict[str, Any]] = []
    orphan_rows: list[dict[str, Any]] = []
    structural_errors = list(audit_errors)

    if recycle_root.exists() and recycle_root.is_symlink():
        structural_errors.append("recycle root is a symbolic link")
    elif recycle_root.is_dir():
        for generation in sorted(recycle_root.iterdir(), key=lambda path: path.name):
            if generation.is_symlink() or not generation.is_dir():
                structural_errors.append(
                    f"unexpected recycle generation entry: {generation}"
                )
                continue
            if not IDENTIFIER_RE.fullmatch(generation.name):
                structural_errors.append(
                    f"invalid recycle generation identifier: {generation.name}"
                )
                continue
            for task_root in sorted(generation.iterdir(), key=lambda path: path.name):
                if task_root.is_symlink() or not task_root.is_dir():
                    structural_errors.append(f"unexpected task entry: {task_root}")
                    continue
                rows.append(
                    _inspect_task(
                        task_root,
                        output_dir=output_dir,
                        current_ids=current_ids,
                        audit_pairs=audit_pairs,
                    )
                )

    # A post-processing worker from an older runtime could finish after task
    # deletion and recreate a task-prefixed sidecar.  Only deletion-audited,
    # non-current task IDs are eligible; unrelated files are never inferred to
    # be disposable merely because they are old.
    task_dir = output_dir / "tasks"
    deleted_ids = {task_id for _deletion_id, task_id in audit_pairs}
    if task_dir.is_dir():
        for task_id in sorted(deleted_ids - current_ids):
            candidates = list(task_dir.glob(f"{task_id}_*"))
            exact = task_dir / task_id
            if exact.exists() or exact.is_symlink():
                candidates.append(exact)
            for path in sorted(set(candidates), key=lambda item: item.name):
                reasons: list[str] = []
                if not _lexically_within(path, task_dir):
                    reasons.append("path is outside the managed task directory")
                bytes_total, entry_count, symlinks = _tree_stats(path)
                if symlinks:
                    reasons.append(
                        f"artifact contains {len(symlinks)} symbolic link(s)"
                    )
                orphan_rows.append({
                    "task_id": task_id,
                    "path": str(path),
                    "bytes": bytes_total,
                    "entry_count": entry_count,
                    "eligible": not reasons,
                    "reasons": reasons,
                })

    eligible = [row for row in rows if row["eligible"]]
    held = [row for row in rows if not row["eligible"]]
    orphan_eligible = [row for row in orphan_rows if row["eligible"]]
    orphan_held = [row for row in orphan_rows if not row["eligible"]]
    return {
        "schema_version": 1,
        "generated_at": _utc_now(),
        "output_dir": str(output_dir),
        "task_registry": {
            "path": str(registry_path),
            "sha256": _sha256(registry_path),
            "current_task_count": len(current_ids),
        },
        "deletion_audit": {
            "path": str(audit_path),
            "sha256": _sha256(audit_path) if audit_path.is_file() else None,
            "confirmed_pair_count": len(audit_pairs),
        },
        "structural_errors": structural_errors,
        "summary": {
            "recycled_task_count": len(rows),
            "eligible_task_count": len(eligible),
            "held_task_count": len(held),
            "eligible_bytes": sum(int(row["bytes"]) for row in eligible),
            "held_bytes": sum(int(row["bytes"]) for row in held),
            "audited_orphan_artifact_count": len(orphan_rows),
            "eligible_orphan_artifact_count": len(orphan_eligible),
            "held_orphan_artifact_count": len(orphan_held),
            "eligible_orphan_bytes": sum(
                int(row["bytes"]) for row in orphan_eligible
            ),
        },
        "eligible": eligible,
        "held": held,
        "orphan_eligible": orphan_eligible,
        "orphan_held": orphan_held,
    }


def _execute(inventory: dict[str, Any], receipt_path: Path) -> dict[str, Any]:
    if inventory["structural_errors"]:
        raise RuntimeError("structural/audit errors make execution unsafe")
    output_dir = Path(inventory["output_dir"])
    registry_path = output_dir / "tasks_index.json"
    if _sha256(registry_path) != inventory["task_registry"]["sha256"]:
        raise RuntimeError("task registry changed after audit; rerun the command")

    receipt = dict(inventory)
    receipt["mode"] = "execute"
    receipt["execution_started_at"] = _utc_now()
    receipt["removed"] = []
    receipt["removed_orphan_artifacts"] = []
    receipt["failures"] = []
    _atomic_json(receipt_path, receipt)

    for row in inventory["eligible"]:
        task_root = Path(row["task_root"])
        try:
            if not task_root.is_dir() or task_root.is_symlink():
                raise RuntimeError("validated task recycle directory changed")
            shutil.rmtree(task_root)
            receipt["removed"].append(
                {
                    "deletion_id": row["deletion_id"],
                    "task_id": row["task_id"],
                    "bytes": row["bytes"],
                    "entry_count": row["entry_count"],
                }
            )
            generation = task_root.parent
            try:
                generation.rmdir()
            except OSError:
                pass
        except Exception as exc:  # preserve a machine-readable partial receipt
            receipt["failures"].append(
                {
                    "deletion_id": row["deletion_id"],
                    "task_id": row["task_id"],
                    "error": str(exc),
                }
            )
        _atomic_json(receipt_path, receipt)

    for row in inventory["orphan_eligible"]:
        path = Path(row["path"])
        try:
            if path.is_symlink() or not path.exists():
                raise RuntimeError("validated orphan artifact changed")
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
            receipt["removed_orphan_artifacts"].append({
                "task_id": row["task_id"],
                "path": row["path"],
                "bytes": row["bytes"],
                "entry_count": row["entry_count"],
            })
        except Exception as exc:
            receipt["failures"].append({
                "task_id": row["task_id"],
                "path": row["path"],
                "error": str(exc),
            })
        _atomic_json(receipt_path, receipt)

    receipt["execution_finished_at"] = _utc_now()
    receipt["removed_task_count"] = len(receipt["removed"])
    receipt["removed_orphan_artifact_count"] = len(
        receipt["removed_orphan_artifacts"]
    )
    receipt["removed_bytes"] = (
        sum(int(row["bytes"]) for row in receipt["removed"])
        + sum(
            int(row["bytes"])
            for row in receipt["removed_orphan_artifacts"]
        )
    )
    receipt["failure_count"] = len(receipt["failures"])
    receipt["post_inventory"] = _inventory(output_dir)["summary"]
    _atomic_json(receipt_path, receipt)

    audit_path = output_dir / "runtime" / "task_deletion_audit.jsonl"
    _append_jsonl(
        audit_path,
        {
            "event": "task_recycle_bin_permanent_purge",
            "at": receipt["execution_finished_at"],
            "receipt_path": str(receipt_path),
            "removed_task_count": receipt["removed_task_count"],
            "removed_orphan_artifact_count": receipt[
                "removed_orphan_artifact_count"
            ],
            "removed_bytes": receipt["removed_bytes"],
            "failure_count": receipt["failure_count"],
            "task_ids": [row["task_id"] for row in receipt["removed"]],
            "orphan_task_ids": sorted({
                row["task_id"]
                for row in receipt["removed_orphan_artifacts"]
            }),
            "deletion_ids": sorted({row["deletion_id"] for row in receipt["removed"]}),
        },
    )
    return receipt


def _human_bytes(value: int) -> str:
    amount = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if amount < 1024.0 or unit == "TiB":
            return f"{amount:.2f} {unit}"
        amount /= 1024.0
    raise AssertionError("unreachable")


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=project_root / "output",
        help="managed V15 output directory",
    )
    parser.add_argument(
        "--receipt",
        type=Path,
        help="audit receipt path (default: output/runtime/timestamped JSON)",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="permanently remove eligible recycled task directories",
    )
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    receipt_path = (
        args.receipt.resolve()
        if args.receipt
        else output_dir / "runtime" / f"task_recycle_purge_{stamp}.json"
    )

    try:
        inventory = _inventory(output_dir)
        if args.execute:
            result = _execute(inventory, receipt_path)
            print(
                f"removed={result['removed_task_count']} "
                f"orphans={result['removed_orphan_artifact_count']} "
                f"bytes={result['removed_bytes']} "
                f"human={_human_bytes(result['removed_bytes'])} "
                f"failures={result['failure_count']}"
            )
            print(f"receipt={receipt_path}")
            return 0 if not result["failure_count"] else 2

        inventory["mode"] = "dry-run"
        _atomic_json(receipt_path, inventory)
        summary = inventory["summary"]
        print(
            f"eligible={summary['eligible_task_count']} "
            f"bytes={summary['eligible_bytes']} "
            f"human={_human_bytes(summary['eligible_bytes'])} "
            f"orphans={summary['eligible_orphan_artifact_count']} "
            f"orphan_bytes={summary['eligible_orphan_bytes']} "
            f"held={summary['held_task_count']} "
            f"structural_errors={len(inventory['structural_errors'])}"
        )
        print(f"receipt={receipt_path}")
        return 0 if not inventory["structural_errors"] else 2
    except Exception as exc:
        print(f"purge refused: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
