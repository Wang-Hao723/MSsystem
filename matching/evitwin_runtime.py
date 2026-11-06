"""Runtime applicability helpers for EviTwin-MS.

The reference registry describes the reference instrument.  It must never be
copied into the query metadata: the query family is resolved independently
from the saved acquisition record or, as a fallback, from mzML metadata.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping

from utils.instrument_metadata import (
    extract_mzml_instrument_metadata,
    normalize_instrument_family,
)


_ACQ_METHOD_RE = re.compile(
    r'<sourceFile id="AcqMethod\.xml" name="AcqMethod\.xml"[^>]*>'
    r'.*?<cvParam[^>]+name="SHA-1" value="([0-9a-fA-F]{40})"',
    re.S,
)


def _task_instrument_labels(task: Mapping[str, Any]) -> list[str]:
    conversion = task.get("vendor_conversion") or {}
    source = conversion.get("source") or {}
    acquisition = source.get("acquisition") or conversion.get("acquisition") or {}
    output_instrument = (conversion.get("output") or {}).get("instrument") or {}
    labels: list[str] = []
    for value in (
        task.get("instrument_family"),
        task.get("instrument_type"),
        task.get("instrument"),
        acquisition.get("instrument"),
    ):
        if str(value or "").strip():
            labels.append(str(value))
    for value in acquisition.get("instrument_models") or ():
        if str(value or "").strip():
            labels.append(str(value))
    for model in output_instrument.get("models") or ():
        if isinstance(model, Mapping) and str(model.get("name") or "").strip():
            labels.append(str(model.get("name")))
    for value in acquisition.get("devices") or ():
        if str(value or "").strip():
            labels.append(str(value))
    return labels


def resolve_evitwin_query_instrument_family(
    task: Mapping[str, Any],
    *,
    reference_instrument_family: str,
) -> dict[str, Any]:
    """Resolve the query analyzer independently and return an audit record."""

    labels = _task_instrument_labels(task)
    query_family = normalize_instrument_family(*labels)
    source = "saved_task_acquisition_metadata" if query_family else ""
    inspected_path = ""
    if not query_family:
        for key in ("analysis_file_path", "raw_input_file_path", "file_path"):
            value = str(task.get(key) or "").strip()
            if not value:
                continue
            path = Path(value)
            if not path.is_file():
                continue
            metadata = extract_mzml_instrument_metadata(path)
            query_family = normalize_instrument_family(
                metadata.get("instrument_type"),
                metadata.get("instrument"),
            )
            inspected_path = str(path.resolve())
            if query_family:
                source = "mzml_instrument_configuration"
                break

    reference_family = normalize_instrument_family(reference_instrument_family)
    compatible = bool(
        query_family and reference_family and query_family == reference_family
    )
    # AuthenticAnchor entries currently carry the registry's specific family
    # token.  Reuse that token only after independent family compatibility has
    # been established; otherwise pass the observed family so scoring abstains.
    scoring_family = (
        str(reference_instrument_family)
        if compatible
        else query_family or "unknown"
    )
    return {
        "query_instrument_family": query_family or "unknown",
        "reference_instrument_family": str(reference_instrument_family),
        "instrument_family_compatible": compatible,
        "scoring_instrument_family": scoring_family,
        "resolution_source": source or "unavailable",
        "inspected_mzml_path": inspected_path,
        "saved_instrument_labels": labels,
    }


def resolve_evitwin_query_method_id(task: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve an explicit LC-MS method id or exact embedded method hash."""

    declared = str(((task.get("params") or {}).get("rt_method_id")) or "").strip()
    if declared:
        return {
            "query_method_id": declared,
            "resolution_source": "saved_task_rt_method_id",
            "inspected_mzml_path": None,
            "acq_method_xml_sha1": None,
            "method_inferred_from_filename": False,
        }
    for key in ("analysis_file_path", "raw_input_file_path", "file_path"):
        value = str(task.get(key) or "").strip()
        if not value:
            continue
        path = Path(value)
        if not path.is_file() or path.suffix.casefold() != ".mzml":
            continue
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            header = handle.read(400_000)
        match = _ACQ_METHOD_RE.search(header)
        if match:
            digest = match.group(1).lower()
            return {
                "query_method_id": f"acqmethod-sha1:{digest}",
                "resolution_source": "mzml_source_file_acq_method_sha1",
                "inspected_mzml_path": str(path.resolve()),
                "acq_method_xml_sha1": digest,
                "method_inferred_from_filename": False,
            }
    return {
        "query_method_id": "",
        "resolution_source": "unavailable",
        "inspected_mzml_path": None,
        "acq_method_xml_sha1": None,
        "method_inferred_from_filename": False,
    }


__all__ = [
    "resolve_evitwin_query_instrument_family",
    "resolve_evitwin_query_method_id",
]
