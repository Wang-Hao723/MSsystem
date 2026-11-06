"""Compatibility-first orchestration for the staged V10 architecture split."""

from .result_contract import (
    identification_contract_path,
    iter_legacy_result_artifact,
    read_identification_contract,
    runtime_snapshot_id_from_manifest,
    write_identification_contract,
)

__all__ = [
    "identification_contract_path",
    "iter_legacy_result_artifact",
    "read_identification_contract",
    "runtime_snapshot_id_from_manifest",
    "write_identification_contract",
]
