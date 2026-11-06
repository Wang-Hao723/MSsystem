"""Technical-repeat consistency view for completed LC-MS/MS tasks.

This projection does not rescore spectra and does not promote a single-injection
candidate.  It groups explicit numeric filename repeats only when the saved
analysis parameters are identical, then reports how often each already-saved
compound group was detected and reached High evidence.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from itertools import combinations
from pathlib import Path
from typing import Any, Mapping


REPLICATE_CONSISTENCY_SCHEMA_VERSION = "1.0.0"
REPLICATE_CONSISTENCY_ALGORITHM_ID = (
    "mssystem-v12-technical-repeat-consistency-1.0.0"
)

_REPLICATE_FILENAME = re.compile(
    r"^(?P<base>.+?)[_-](?:rep(?:licate)?[_-]?)?(?P<rep>\d+)"
    r"(?P<extension>\.mzml|\.mzxml|\.mgf)$",
    re.IGNORECASE,
)


def _text(value: Any) -> str:
    return " ".join(str(value or "").split())


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def parse_replicate_filename(filename: Any) -> dict[str, Any] | None:
    """Parse only an explicit terminal repeat number such as ``sample_2.mzML``."""

    name = Path(_text(filename)).name
    match = _REPLICATE_FILENAME.fullmatch(name)
    if not match:
        return None
    return {
        "filename": name,
        "base": match.group("base"),
        "base_key": match.group("base").casefold(),
        "replicate_number": int(match.group("rep")),
        "extension": match.group("extension").casefold(),
    }


def _parameter_payload(task: Mapping[str, Any]) -> dict[str, Any]:
    params = task.get("params")
    return dict(params) if isinstance(params, Mapping) else {}


def parameter_signature(task: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        _parameter_payload(task),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _read_compounds(task_dir: Path, task_id: str) -> dict[str, Any] | None:
    from pipeline.compound_conclusions import compound_conclusion_artifact_path
    artifact = compound_conclusion_artifact_path(task_dir, task_id)
    try:
        payload = json.loads(artifact.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    if not isinstance(payload, Mapping) or not isinstance(
        payload.get("compounds"), list
    ):
        return None
    return dict(payload)


def _is_main(compound: Mapping[str, Any]) -> bool:
    conclusion = compound.get("conclusion")
    return bool(
        isinstance(conclusion, Mapping)
        and conclusion.get("counted_as_sample_compound")
    )


def _is_high(compound: Mapping[str, Any]) -> bool:
    conclusion = compound.get("conclusion")
    return bool(
        isinstance(conclusion, Mapping)
        and conclusion.get("key") in {"high", "isomer"}
        and conclusion.get("high_confidence_group_detected")
    )


def _lead(compound: Mapping[str, Any]) -> dict[str, Any]:
    lead = compound.get("lead_candidate")
    lead = lead if isinstance(lead, Mapping) else {}
    representative = compound.get("representative")
    representative = representative if isinstance(representative, Mapping) else {}
    match = representative.get("match")
    match = match if isinstance(match, Mapping) else {}
    return {
        "compound_name": _text(
            lead.get("compound_name")
            or match.get("compound_name")
            or compound.get("display_name")
        ),
        "inchikey": _text(lead.get("inchikey") or match.get("inchikey")),
        "formula": _text(
            compound.get("formula") or lead.get("formula") or match.get("formula")
        ),
        "rank_score": _finite(
            lead.get("rank_score")
            or match.get("rank_score")
            or match.get("final_score")
        ),
        "cosine": _finite(lead.get("cosine") or match.get("cosine")),
        "retention_time": _finite(representative.get("retention_time")),
        "precursor_mz": _finite(representative.get("precursor_mz")),
    }


def _jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return round(len(left & right) / len(union), 4) if union else 1.0


def _unavailable(
    task_id: str,
    reason: str,
    reason_zh: str,
) -> dict[str, Any]:
    return {
        "schema_version": REPLICATE_CONSISTENCY_SCHEMA_VERSION,
        "algorithm_id": REPLICATE_CONSISTENCY_ALGORITHM_ID,
        "task_id": task_id,
        "status": "unavailable",
        "reason": reason,
        "reason_zh": reason_zh,
        "ranking_changed": False,
        "confidence_changed": False,
    }


def build_replicate_consistency(
    task_id: str,
    tasks: Mapping[str, Mapping[str, Any]],
    task_dir: Path,
) -> dict[str, Any]:
    """Build a deterministic repeatability projection from compound sidecars."""

    current = tasks.get(task_id)
    if not isinstance(current, Mapping):
        return _unavailable(task_id, "task_not_found", "任务不存在")
    parsed_current = parse_replicate_filename(current.get("filename"))
    if not parsed_current:
        return _unavailable(
            task_id,
            "filename_has_no_explicit_replicate_suffix",
            "文件名末尾没有可识别的重复针编号",
        )
    signature = parameter_signature(current)
    candidates: dict[int, tuple[str, Mapping[str, Any], dict[str, Any]]] = {}
    for sibling_id, sibling in tasks.items():
        if not isinstance(sibling, Mapping) or sibling.get("status") != "completed":
            continue
        parsed = parse_replicate_filename(sibling.get("filename"))
        if not parsed:
            continue
        if (
            parsed["base_key"] != parsed_current["base_key"]
            or parsed["extension"] != parsed_current["extension"]
            or parameter_signature(sibling) != signature
        ):
            continue
        repeat_number = int(parsed["replicate_number"])
        existing = candidates.get(repeat_number)
        if sibling_id == task_id:
            candidates[repeat_number] = (sibling_id, sibling, parsed)
            continue
        if existing and existing[0] == task_id:
            continue
        if not existing or _text(sibling.get("end_time")) > _text(
            existing[1].get("end_time")
        ):
            candidates[repeat_number] = (sibling_id, sibling, parsed)

    selected = [candidates[key] for key in sorted(candidates)]
    if len(selected) < 2:
        return _unavailable(
            task_id,
            "fewer_than_two_completed_compatible_repeats",
            "没有找到至少两针参数一致的已完成重复分析",
        )

    injections: list[dict[str, Any]] = []
    payloads: dict[str, dict[str, Any]] = {}
    for sibling_id, sibling, parsed in selected:
        payload = _read_compounds(task_dir, sibling_id)
        if payload is None:
            continue
        payloads[sibling_id] = payload
        compounds = [
            item for item in payload.get("compounds") or []
            if isinstance(item, Mapping) and _text(item.get("identity_key"))
        ]
        injections.append({
            "task_id": sibling_id,
            "filename": _text(sibling.get("filename")),
            "replicate_number": parsed["replicate_number"],
            "compound_group_count": len(compounds),
            "main_group_count": sum(_is_main(item) for item in compounds),
            "high_group_count": sum(_is_high(item) for item in compounds),
        })
    if len(injections) < 2:
        return _unavailable(
            task_id,
            "fewer_than_two_readable_compound_artifacts",
            "至少一针的化合物结论文件不可用",
        )

    injection_count = len(injections)
    robustness_threshold = math.ceil(2 * injection_count / 3)
    by_task: dict[str, dict[str, Mapping[str, Any]]] = {}
    for injection in injections:
        sibling_id = injection["task_id"]
        by_task[sibling_id] = {
            _text(item.get("identity_key")): item
            for item in payloads[sibling_id].get("compounds") or []
            if isinstance(item, Mapping) and _text(item.get("identity_key"))
        }

    identity_union = set().union(*(set(items) for items in by_task.values()))
    compound_rows: list[dict[str, Any]] = []
    for identity_key in sorted(identity_union):
        per_injection: list[dict[str, Any]] = []
        observed: list[Mapping[str, Any]] = []
        for injection in injections:
            sibling_id = injection["task_id"]
            compound = by_task[sibling_id].get(identity_key)
            if compound is not None:
                observed.append(compound)
            per_injection.append({
                "task_id": sibling_id,
                "filename": injection["filename"],
                "replicate_number": injection["replicate_number"],
                "detected": compound is not None,
                "main_result": _is_main(compound) if compound else False,
                "high": _is_high(compound) if compound else False,
                "conclusion": _text(
                    ((compound or {}).get("conclusion") or {}).get("key")
                ),
                "lead_candidate": _lead(compound) if compound else None,
            })
        detected_count = sum(item["detected"] for item in per_injection)
        main_count = sum(item["main_result"] for item in per_injection)
        high_count = sum(item["high"] for item in per_injection)
        current_compound = by_task.get(task_id, {}).get(identity_key)
        representative = current_compound or max(
            observed,
            key=lambda item: _lead(item).get("rank_score") or -1.0,
        )
        current_lead = _lead(representative)
        compound_rows.append({
            "identity_key": identity_key,
            "identity_scope": _text(representative.get("identity_scope")),
            "formula": _text(representative.get("formula")),
            "display_name": _text(representative.get("display_name")),
            "representative_lead_candidate": current_lead,
            "detected_injection_count": detected_count,
            "main_injection_count": main_count,
            "high_injection_count": high_count,
            "injection_count": injection_count,
            "robust_main": main_count >= robustness_threshold,
            "robust_high": high_count >= robustness_threshold,
            "all_injection_high": high_count == injection_count,
            "single_injection_high_only": high_count == 1,
            "per_injection": per_injection,
        })

    high_sets = {
        task_key: {
            identity for identity, compound in compounds.items() if _is_high(compound)
        }
        for task_key, compounds in by_task.items()
    }
    main_sets = {
        task_key: {
            identity for identity, compound in compounds.items() if _is_main(compound)
        }
        for task_key, compounds in by_task.items()
    }
    pairwise = []
    for left, right in combinations(injections, 2):
        left_id = left["task_id"]
        right_id = right["task_id"]
        pairwise.append({
            "left_task_id": left_id,
            "right_task_id": right_id,
            "left_filename": left["filename"],
            "right_filename": right["filename"],
            "high_intersection_count": len(high_sets[left_id] & high_sets[right_id]),
            "high_union_count": len(high_sets[left_id] | high_sets[right_id]),
            "high_jaccard": _jaccard(high_sets[left_id], high_sets[right_id]),
            "main_intersection_count": len(main_sets[left_id] & main_sets[right_id]),
            "main_union_count": len(main_sets[left_id] | main_sets[right_id]),
            "main_jaccard": _jaccard(main_sets[left_id], main_sets[right_id]),
        })

    compound_rows.sort(key=lambda item: (
        -int(item["robust_high"]),
        -int(item["high_injection_count"]),
        -int(item["main_injection_count"]),
        _text(item["display_name"]).casefold(),
        item["identity_key"],
    ))
    high_union = set().union(*high_sets.values())
    main_union = set().union(*main_sets.values())
    return {
        "schema_version": REPLICATE_CONSISTENCY_SCHEMA_VERSION,
        "algorithm_id": REPLICATE_CONSISTENCY_ALGORITHM_ID,
        "task_id": task_id,
        "status": "available",
        "replicate_group_key": parsed_current["base_key"],
        "replicate_group_label": parsed_current["base"],
        "grouping_basis": "terminal_filename_repeat_number_and_identical_parameters",
        "parameter_signature_sha256": signature,
        "injection_count": injection_count,
        "robustness_threshold": robustness_threshold,
        "injections": injections,
        "summary": {
            "high_union_count": len(high_union),
            "robust_high_count": sum(item["robust_high"] for item in compound_rows),
            "all_injection_high_count": sum(
                item["all_injection_high"] for item in compound_rows
            ),
            "single_injection_high_count": sum(
                item["single_injection_high_only"] for item in compound_rows
            ),
            "main_union_count": len(main_union),
            "robust_main_count": sum(item["robust_main"] for item in compound_rows),
        },
        "pairwise": pairwise,
        "compounds": compound_rows,
        "interpretation": (
            "Technical-repeat robustness is a separate cohort projection. It does "
            "not rescore spectra, change per-injection ranking, or constitute "
            "independent identity evidence."
        ),
        "interpretation_zh": (
            "重复针稳健性是独立展示层：不重算谱图、不改动单针候选排序或置信等级，"
            "也不等同于独立的结构确证证据。"
        ),
        "ranking_changed": False,
        "confidence_changed": False,
    }


__all__ = [
    "REPLICATE_CONSISTENCY_ALGORITHM_ID",
    "REPLICATE_CONSISTENCY_SCHEMA_VERSION",
    "build_replicate_consistency",
    "parameter_signature",
    "parse_replicate_filename",
]
