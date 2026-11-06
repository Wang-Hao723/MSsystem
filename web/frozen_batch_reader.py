"""Stdlib-only reader executed over SSH; reads saved JSON, never analysis code."""

from __future__ import annotations

import hashlib
import csv
import io
import json
import lzma
import re
import sys
from pathlib import Path

ROOT = Path("/data/Mssystem/xbt_real_batch_20261021_ms2")
MAX_RESPONSE = 64 * 1024**2
MAX_COMPOUND_PAGE = 500
_EXTERNAL_COMPOUND_ROLES = frozenset({
    "external_contaminant", "curated_contaminant", "pesticide_screening_target",
})


def _native_compounds(document):
    payload = document.get("compound_conclusions")
    if not isinstance(payload, dict) or not isinstance(payload.get("compounds"), list):
        raise ValueError("Frozen native compound conclusions are missing")
    return payload


def _explicit_external_review(row):
    representative = row.get("representative") or {}
    match = representative.get("match") or {}
    lead = row.get("lead_candidate") or {}
    plant = row.get("plant_report") or {}
    roles = {
        str(value or "").strip()
        for value in (
            row.get("compound_role"), lead.get("compound_role"),
            match.get("compound_role"), plant.get("role"),
        )
    }
    return bool(
        plant.get("excluded_from_plant_candidate_table") is True
        or match.get("potential_contaminant") is True
        or roles.intersection(_EXTERNAL_COMPOUND_ROLES)
    )


def _compound_in_scope(row, report_scope):
    if report_scope == "all":
        return True
    if report_scope == "plant_candidates":
        return bool((row.get("conclusion") or {}).get("counted_as_plant_candidate")) \
            and not _explicit_external_review(row)
    if report_scope == "external_review":
        return _explicit_external_review(row)
    raise ValueError("Unknown compound report scope")


def _compact_compound(row, index):
    """Drop bulky display arrays while retaining saved identity/evidence fields."""
    compact = dict(row)
    compact["frozen_compound_index"] = index
    omitted = []
    for key in ("observations", "occurrences"):
        if key in compact:
            compact.pop(key)
            omitted.append(key)
    representative = compact.get("representative")
    if isinstance(representative, dict):
        representative = dict(representative)
        for key in ("raw_peaks", "processed_peaks", "query_peaks", "library_peaks", "peaks"):
            if key in representative:
                representative.pop(key)
                omitted.append("representative." + key)
        match = representative.get("match")
        if isinstance(match, dict):
            match = dict(match)
            for key in ("matched_peaks", "peak_pairs", "query_peaks", "library_peaks", "peaks"):
                if key in match:
                    match.pop(key)
                    omitted.append("representative.match." + key)
            representative["match"] = match
        compact["representative"] = representative
    compact["_frozen_projection"] = {
        "kind": "compact_saved_compound_conclusion",
        "omitted_heavy_fields": omitted,
        "full_row_operation": "compound_item",
    }
    return compact


def _candidate_identity_tokens(candidate, identity_scope=""):
    tokens = set()
    supplied = str(candidate.get("identity_key") or "").strip()
    if supplied:
        tokens.add("identity:" + supplied)
        connectivity = re.search(r"(?:std2d:|ik::)?([A-Z]{14})(?:[-|:]|$)", supplied, re.I)
        if connectivity:
            tokens.add("connectivity:" + connectivity.group(1).upper())
    inchikey = str(candidate.get("inchikey") or "").strip().upper().removeprefix("INCHIKEY=")
    if re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", inchikey):
        tokens.add("inchikey:" + inchikey)
        if re.search(r"2d|connectivity|isomer_group|std2d", identity_scope, re.I):
            tokens.add("connectivity:" + inchikey[:14])
    smiles = str(candidate.get("smiles") or "").strip()
    if smiles:
        tokens.add("smiles:" + smiles)
    name = str(candidate.get("compound_name") or "").strip().casefold()
    formula = str(candidate.get("formula") or "").strip()
    if name and formula:
        tokens.add(f"name-formula:{name}|{formula}")
    return tokens


def _lead_scores_are_bound(row):
    representative = row.get("representative") or {}
    match = representative.get("match") or {}
    lead = row.get("lead_candidate") or {}
    rep_index = representative.get("result_index")
    stored_index = row.get("representative_result_index")
    lead_index = lead.get("result_index", lead.get("representative_result_index"))
    if stored_index is not None and rep_index is not None and int(stored_index) != int(rep_index):
        return False
    if lead_index is not None and rep_index is not None and int(lead_index) != int(rep_index):
        return False
    if (lead.get("spectrum_id") and representative.get("spectrum_id")
            and str(lead["spectrum_id"]) != str(representative["spectrum_id"])):
        return False
    scope = str(row.get("identity_scope") or row.get("identity_key") or "")
    representative_identity = {
        **match,
        "identity_key": match.get("identity_key") or row.get("identity_key"),
        "formula": match.get("formula") or row.get("formula"),
    }
    lead_tokens = _candidate_identity_tokens(lead, scope)
    representative_tokens = _candidate_identity_tokens(representative_identity, scope)
    structural = ("identity:", "connectivity:", "inchikey:", "smiles:")
    lead_structural = {token for token in lead_tokens if token.startswith(structural)}
    representative_structural = {
        token for token in representative_tokens if token.startswith(structural)
    }
    if lead_structural and representative_structural:
        return bool(lead_structural & representative_structural)
    return bool(lead_tokens & representative_tokens)


def _full_inchikey(value):
    value = str(value or "").strip().upper().removeprefix("INCHIKEY=")
    return value if re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", value) else ""


def _bound_spectral_evidence(row):
    representative = row.get("representative") or {}
    match = representative.get("match") or {}
    lead = row.get("lead_candidate") or {}
    binding = match.get("spectral_evidence_binding") or {}
    if (binding.get("schema") != "same-scan-selected-reference/1"
            or binding.get("status") != "bound_saved_reference"):
        return None
    if str(binding.get("spectrum_id") or "") != str(representative.get("spectrum_id") or ""):
        return None
    if (binding.get("result_index") is None or representative.get("result_index") is None
            or int(binding["result_index"]) != int(representative["result_index"])):
        return None
    binding_key = _full_inchikey(binding.get("candidate_inchikey"))
    current_keys = [
        value for value in (_full_inchikey(match.get("inchikey")),
                            _full_inchikey(lead.get("inchikey"))) if value
    ]
    if not binding_key or not current_keys or any(value != binding_key for value in current_keys):
        return None
    binding_reference = str(binding.get("lib_spectrum_id") or "").strip()
    current_references = [
        str(value).strip() for value in (match.get("lib_spectrum_id"), lead.get("lib_spectrum_id"))
        if str(value or "").strip()
    ]
    if (not binding_reference or not current_references
            or any(value != binding_reference for value in current_references)):
        return None
    return binding


def _bound_score(row, field):
    representative = row.get("representative") or {}
    match = representative.get("match") or {}
    declared_binding = match.get("spectral_evidence_binding") or {}
    binding = _bound_spectral_evidence(row)
    if declared_binding.get("status") == "bound_saved_reference" and binding is None:
        return None
    if binding is not None:
        value = (binding.get("metrics") or {}).get(field)
        return value if value is not None and str(value).strip() != "" else None
    value = match.get(field)
    if value is not None and str(value).strip() != "":
        return value
    lead = row.get("lead_candidate") or {}
    return lead.get(field) if _lead_scores_are_bound(row) else None


def _compound_csv(document, request):
    payload = _native_compounds(document)
    output = io.StringIO(newline="")
    fields = [
        "sample_key", "source_sha256", "frozen_compound_index", "compound_id",
        "result_id", "feature_id", "record_unit", "record_type", "identity_key",
        "identity_scope", "identity_source", "display_name", "formula",
        "conclusion_key", "confirmation_level", "counted_as_sample_compound",
        "counted_as_plant_candidate", "evidence_scope", "representative_result_index",
        "representative_spectrum_id", "lead_compound_name", "lead_compound_name_cn",
        "lead_formula", "lead_smiles", "lead_inchikey", "lead_adduct",
        "lead_confidence", "lead_cosine", "lead_mass_ppm", "lead_rank",
        "lead_score_projection_bound", "lead_cosine_raw", "lead_mass_ppm_raw",
        "representative_match_cosine_raw", "representative_match_mass_ppm_raw", "blank_flagged",
        "saved_spectrum_confidence", "spectral_evidence_binding_json",
        "blank_clean_view_excluded", "plant_report_excluded",
        "plant_report_review_category", "plant_report_reason", "conclusion_json",
        "plant_report_json", "lead_candidate_json", "candidate_members_json",
    ]
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for index, row in enumerate(payload["compounds"]):
        conclusion = row.get("conclusion") or {}
        plant = row.get("plant_report") or {}
        lead = row.get("lead_candidate") or {}
        representative = row.get("representative") or {}
        match = representative.get("match") or {}
        writer.writerow({
            "sample_key": request["sample_key"], "source_sha256": request["source_sha256"],
            "frozen_compound_index": index, "compound_id": row.get("compound_id"),
            "result_id": row.get("result_id"), "feature_id": row.get("feature_id"),
            "record_unit": row.get("record_unit"), "record_type": row.get("record_type"),
            "identity_key": row.get("identity_key"), "identity_scope": row.get("identity_scope"),
            "identity_source": row.get("identity_source"), "display_name": row.get("display_name"),
            "formula": row.get("formula"), "conclusion_key": conclusion.get("key"),
            "confirmation_level": conclusion.get("confirmation_level"),
            "counted_as_sample_compound": conclusion.get("counted_as_sample_compound"),
            "counted_as_plant_candidate": conclusion.get("counted_as_plant_candidate"),
            "evidence_scope": row.get("evidence_scope") or conclusion.get("evidence_scope"),
            "representative_result_index": representative.get("result_index"),
            "representative_spectrum_id": representative.get("spectrum_id"),
            "lead_compound_name": lead.get("compound_name"),
            "lead_compound_name_cn": lead.get("compound_name_cn"), "lead_formula": lead.get("formula"),
            "lead_smiles": lead.get("smiles"), "lead_inchikey": lead.get("inchikey"),
            "lead_adduct": lead.get("adduct"), "lead_confidence": lead.get("confidence"),
            "lead_cosine": _bound_score(row, "cosine"),
            "lead_mass_ppm": _bound_score(row, "mass_ppm"),
            "lead_rank": lead.get("rank", lead.get("candidate_rank")),
            "lead_score_projection_bound": bool(
                _bound_spectral_evidence(row) is not None or _lead_scores_are_bound(row)
            ),
            "lead_cosine_raw": lead.get("cosine"),
            "lead_mass_ppm_raw": lead.get("mass_ppm"),
            "representative_match_cosine_raw": match.get("cosine"),
            "representative_match_mass_ppm_raw": match.get("mass_ppm"),
            "saved_spectrum_confidence": (
                (_bound_spectral_evidence(row) or {}).get("saved_confidence")
                or match.get("native_spectrum_confidence")
            ),
            "spectral_evidence_binding_json": json.dumps(
                match.get("spectral_evidence_binding") or {},
                ensure_ascii=False, separators=(",", ":"),
            ),
            "blank_flagged": representative.get("blank_flagged"),
            "blank_clean_view_excluded": representative.get("blank_clean_view_excluded"),
            "plant_report_excluded": plant.get("excluded_from_plant_candidate_table"),
            "plant_report_review_category": plant.get("review_category"),
            "plant_report_reason": plant.get("reason"),
            "conclusion_json": json.dumps(conclusion, ensure_ascii=False, separators=(",", ":")),
            "plant_report_json": json.dumps(plant, ensure_ascii=False, separators=(",", ":")),
            "lead_candidate_json": json.dumps(lead, ensure_ascii=False, separators=(",", ":")),
            "candidate_members_json": json.dumps(row.get("candidate_members") or [], ensure_ascii=False, separators=(",", ":")),
        })
    return ("\ufeff" + output.getvalue()).encode("utf-8")


def read_document(request, kind, *, return_raw=False):
    key = request["sample_key"]
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+__[a-f0-9]{16}", key):
        raise ValueError("Invalid frozen sample key")
    record = request[kind]
    path = ROOT / kind / (key + ".json.xz")
    if str(path) != record["path"]:
        raise ValueError("Frozen path binding mismatch")
    raw = path.read_bytes()
    if len(raw) != record["size"] or hashlib.sha256(raw).hexdigest() != record["sha256"]:
        raise ValueError("Frozen compressed result SHA/size mismatch")
    value = json.loads(lzma.decompress(raw))
    if value.get("sample_key") != key:
        raise ValueError("Frozen sample identity mismatch")
    if kind == "results":
        if value["source"]["sha256"] != request["source_sha256"]:
            raise ValueError("Frozen source SHA mismatch")
        if value.get("protocol_sha256") != request["native_protocol_sha256"]:
            raise ValueError("Frozen native protocol mismatch")
    else:
        binding = value["native_result_binding"]
        expected = request["results"]
        if (binding["compressed_sha256"] != expected["sha256"]
                or binding["compressed_size"] != expected["size"]
                or binding["remote_path"] != expected["path"]
                or binding["source_sha256"] != request["source_sha256"]
                or binding["native_protocol_sha256"] != request["native_protocol_sha256"]):
            raise ValueError("Frozen conclusion/native binding mismatch")
    return raw if return_raw else value


def project(request):
    operation = request["operation"]
    if operation in {"summary", "item", "items", "expected", "native_metadata"}:
        document = read_document(request, "results")
        rows = document["native_ms2_results"]
        if operation == "item":
            index = request["index"]
            if not isinstance(index, int) or index < 0 or index >= len(rows):
                raise IndexError("Saved spectrum index out of range")
            return {"result": rows[index], "result_index": index}
        if operation == "items":
            indices = request["indices"]
            if len(indices) > 500 or any(not isinstance(i, int) or i < 0 or i >= len(rows) for i in indices):
                raise ValueError("Full saved result page exceeds 500 records or invalid index")
            return {"results": [{**rows[index], "result_index": index} for index in indices]}
        if operation == "summary":
            fields = set(request["row_fields"])
            match_fields = set(request["match_fields"])
            result = []
            for index, row in enumerate(rows):
                selected = {key: value for key, value in row.items() if key in fields}
                selected["result_index"] = index
                selected["match"] = ({key: value for key, value in row["match"].items()
                                      if key in match_fields} if row.get("match") else None)
                result.append(selected)
            return {"results": result, "total": len(rows)}
        if operation == "expected":
            layer = document["ms1_integration"]
            return {**layer, "status": "ok", "summary": {
                "requested": True, "status": "ok",
                "detected_feature_count": len(layer["features"]),
                "standard_confirmed_count": 0,
            }}
        return {key: document[key] for key in (
            "source", "parameters", "native_summary", "paired_blanks", "runtime")}
    if operation == "native_conclusions_xz":
        return read_document(request, "conclusions", return_raw=True)
    document = read_document(request, "conclusions")
    if operation in {"compound_page", "compound_item", "native_compounds_csv"}:
        native = _native_compounds(document)
        if operation == "native_compounds_csv":
            return _compound_csv(document, request)
        result_view = request.get("result_view", "native")
        if result_view not in {"native", "xbt_linked"}:
            raise ValueError("Unknown saved compound result view")
        if result_view == "xbt_linked":
            selected = document["report_views"].get("xbt_named_candidates")
            if not isinstance(selected, dict) or not isinstance(selected.get("compounds"), list):
                raise ValueError("Frozen XBT-linked conclusion view is missing")
            source_rows = selected["compounds"]
            source_metadata = selected.get("metadata") or {}
        else:
            source_rows = native["compounds"]
            source_metadata = native.get("metadata") or {}
        if operation == "compound_item":
            index = request.get("index")
            compound_id = request.get("compound_id")
            if index is not None:
                if not isinstance(index, int) or index < 0 or index >= len(source_rows):
                    raise IndexError("Saved compound index out of range")
                return {"compound": source_rows[index], "frozen_compound_index": index,
                        "result_view": result_view}
            matches = [(i, row) for i, row in enumerate(source_rows)
                       if row.get("compound_id") == compound_id]
            if len(matches) != 1:
                raise IndexError("Saved compound ID is missing or non-unique")
            index, row = matches[0]
            return {"compound": row, "frozen_compound_index": index, "result_view": result_view}
        offset = request.get("offset", 0)
        limit = request.get("limit", 200)
        if (not isinstance(offset, int) or offset < 0 or not isinstance(limit, int)
                or limit < 1 or limit > MAX_COMPOUND_PAGE):
            raise ValueError("Frozen compound page is invalid or exceeds 500 records")
        scope = request.get("report_scope", "all")
        if result_view == "xbt_linked" and scope != "all":
            raise ValueError("XBT-linked subset cannot be reclassified by report_scope")
        filtered = list(enumerate(source_rows))
        if result_view == "native":
            filtered = [(index, row) for index, row in filtered if _compound_in_scope(row, scope)]
        conclusion = request.get("conclusion")
        if conclusion:
            filtered = [(index, row) for index, row in filtered
                        if (row.get("conclusion") or {}).get("key") == conclusion]
        page = filtered[offset:offset + limit]
        return {
            "metadata": source_metadata,
            "native_metadata": native.get("metadata") or {},
            "compounds": [_compact_compound(row, index) for index, row in page],
            "native_total": len(native["compounds"]),
            "filtered_total": len(filtered),
            "xbt_linked_total": len((document.get("report_views", {}).get("xbt_named_candidates") or {}).get("compounds") or []),
            "offset": offset, "limit": limit, "has_more": offset + limit < len(filtered),
            "selection": {
                "result_view": result_view, "report_scope": scope,
                "source": ("compound_conclusions" if result_view == "native"
                           else "report_views.xbt_named_candidates"),
                "partial_view": result_view == "xbt_linked",
                "selection_changes_spectral_confidence": False,
            },
            "summary": document.get("summary") or {},
        }
    if operation == "flagship":
        artifact = document["primary_identity_artifact"]
        if request.get("entity_id"):
            return {**artifact, "entities": [row for row in artifact["entities"]
                                              if row.get("entity_id") == request["entity_id"]]}
        def candidate_summary(value):
            if value is None:
                return None
            keys = {"identity_key", "entity_score", "confidence", "high_evidence", "score_components"}
            result = {key: item for key, item in value.items() if key in keys}
            candidate_keys = {"compound_name", "compound_name_cn", "formula", "smiles", "inchikey", "adduct",
                              "confidence", "dataset_source", "lib_spectrum_id", "compound_role",
                              "natural_product_candidate", "potential_contaminant"}
            result["candidate"] = {key: item for key, item in (value.get("candidate") or {}).items()
                                   if key in candidate_keys}
            return result
        entities = []
        for row in artifact["entities"]:
            projected = {key: value for key, value in row.items()
                         if key not in {"candidate_recall_support", "candidates", "winner", "leading_candidate_when_abstained"}}
            projected["winner"] = candidate_summary(row.get("winner"))
            projected["leading_candidate_when_abstained"] = candidate_summary(row.get("leading_candidate_when_abstained"))
            projected["_frozen_summary"] = True
            entities.append(projected)
        return {**artifact, "entities": entities, "projection": "summary_full_entity_available_by_id"}
    raise ValueError("Unsupported frozen read operation")


def main():
    result = project(json.load(sys.stdin))
    if isinstance(result, bytes):
        sys.stdout.buffer.write(result)
        return
    encoded = json.dumps(result, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    if len(encoded.encode("utf-8")) > MAX_RESPONSE:
        raise ValueError("Frozen read response exceeds 64 MiB transport limit")
    sys.stdout.write(encoded)


if __name__ == "__main__":
    main()
