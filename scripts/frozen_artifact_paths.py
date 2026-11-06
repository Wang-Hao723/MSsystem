"""Resolve frozen V15 inputs without silently changing their identity."""

import hashlib
import json
from pathlib import Path


def _sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def _candidate_path(recorded_path, project_root):
    original = Path(recorded_path)
    root = Path(project_root).resolve()
    if original.is_file():
        return original

    # Only the explicitly declared sibling baseline can be relocated.
    # Do not search data directories or substitute a similarly named file.
    old_root = root.parent / 'mssystem_V15'
    try:
        relative = original.relative_to(old_root)
    except ValueError:
        raise FileNotFoundError(f'Frozen artifact unavailable: {original}') from None
    candidate = (root / relative).resolve()
    if not candidate.is_relative_to(root) or not candidate.is_file():
        raise FileNotFoundError(f'Frozen artifact unavailable: {original}')
    return candidate


def resolve_frozen_v15_artifact(recorded_path, expected_sha256, project_root):
    if not isinstance(expected_sha256, str) or len(expected_sha256) != 64:
        raise ValueError('A frozen artifact requires its SHA-256 binding')
    candidate = _candidate_path(recorded_path, project_root)
    actual = _sha256(candidate)
    if actual != expected_sha256:
        raise ValueError(f'Frozen artifact SHA-256 mismatch: {candidate}')
    return candidate


def resolve_identity_scope_migrated_v15_artifact(
    recorded_path,
    expected_sha256,
    project_root,
    *,
    task_id,
    migration_receipt_path,
    migration_receipt_sha256,
):
    """Resolve an exact frozen file or its receipt-bound scope migration.

    The fallback is deliberately limited to the V15 identity-scope migration.
    Its receipt must bind the requested pre-migration file hash to the current
    file hash and prove that candidate order, rank, score, and scientific
    identity projections did not change.
    """
    if not isinstance(expected_sha256, str) or len(expected_sha256) != 64:
        raise ValueError('A frozen artifact requires its SHA-256 binding')
    if (
        not isinstance(migration_receipt_sha256, str)
        or len(migration_receipt_sha256) != 64
    ):
        raise ValueError('A migration receipt requires its SHA-256 binding')

    candidate = _candidate_path(recorded_path, project_root)
    actual = _sha256(candidate)
    if actual == expected_sha256:
        return candidate

    receipt_path = Path(migration_receipt_path)
    if not receipt_path.is_file():
        raise FileNotFoundError(
            f'Frozen migration receipt unavailable: {receipt_path}'
        )
    if _sha256(receipt_path) != migration_receipt_sha256:
        raise ValueError(f'Frozen migration receipt SHA-256 mismatch: {receipt_path}')

    receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
    if receipt.get('schema_version') != (
        'mssystem.v15.production-identity-scope-migration/1.0.0'
    ):
        raise ValueError('Unexpected frozen migration receipt schema')
    required_global_invariants = {
        'library_retrieval_repeated': False,
        'candidate_ranking_changed': False,
        'candidate_scores_changed': False,
        'database_scientific_identity_fields_changed': False,
        'migration_failure_count': 0,
    }
    if any(
        receipt.get(key) != value
        for key, value in required_global_invariants.items()
    ):
        raise ValueError('Frozen migration receipt does not preserve ranking inputs')

    matches = [
        row for row in (receipt.get('tasks') or [])
        if str(row.get('task_id') or '') == str(task_id or '')
    ]
    if len(matches) != 1:
        raise ValueError(
            f'Frozen migration receipt has no unique task binding: {task_id}'
        )
    row = matches[0]
    ranking_before = row.get('before_ranking_sha256')
    ranking_after = row.get('after_ranking_sha256')
    identity_before = row.get('before_scientific_identity_sha256')
    identity_after = row.get('after_scientific_identity_sha256')
    preserved = (
        row.get('status') == 'completed'
        and row.get('before_flagship_sha256') == expected_sha256
        and row.get('after_flagship_sha256') == actual
        and row.get('ranking_unchanged') is True
        and isinstance(ranking_before, str)
        and len(ranking_before) == 64
        and ranking_before == ranking_after
        and row.get('database_scientific_identity_fields_unchanged') is True
        and isinstance(identity_before, str)
        and len(identity_before) == 64
        and identity_before == identity_after
    )
    if not preserved:
        raise ValueError(f'Frozen migration binding or ranking invariant failed: {task_id}')
    return candidate
