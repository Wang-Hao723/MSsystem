"""Reject derived identity evidence bound to a previous SQL metadata revision."""
from pathlib import Path
import json

import ijson

from runtime.verified_sql_retrieval import verified_sql_snapshot_id, sql_metadata_repair_preflight


def expected_sql_snapshot(project_root):
    root = Path(project_root)
    if not (root / 'config/sql_metadata_repairs.json').is_file():
        if (root / 'config/library_metadata_repairs.json').is_file():
            raise ValueError('SQL metadata configuration is missing beside active reference repairs')
        return None
    manifest = json.loads((root / 'config/runtime_manifest.json').read_text())
    base = manifest['library']['production_full']['snapshot_id']
    return verified_sql_snapshot_id(base, project_root=root)


def sql_artifact_is_current(path, *, project_root):
    """Both recall and CE must use the same verified metadata as live retrieval.

    An absent configuration leaves non-production fixtures unchanged. A broken
    configured binding is an error, never permission to use legacy evidence.
    """
    expected = expected_sql_snapshot(project_root)
    if expected is None:
        return True
    verification = sql_metadata_repair_preflight(project_root)
    if not verification.get('compatible'):
        raise ValueError('SQL metadata evidence dependencies are invalid: ' + str(verification.get('error')))
    with Path(path).open('rb') as handle:
        metadata = next(ijson.items(handle, 'metadata'), {})
    return all(metadata.get(field) == expected for field in (
        'candidate_recall_support_library_snapshot_id',
        'source_paired_ce_gate_library_snapshot_id',
    ))
