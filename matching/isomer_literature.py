"""Provenance-preserving literature support for unresolved isomer groups.

The module has two deliberately separate layers:

* a small, versioned set of manually scoped methodological/pair-specific
  sources; and
* an optional PubMed bibliographic lookup that returns leads, not spectral
  evidence.

Neither layer changes search scores, Top-1, confidence or exact-identity
status.  Literature text is never converted into measured fragment
intensities.  Exact identification still requires comparable experimental
evidence such as an authenticated standard, validated CCS, or a demonstrated
stereoselective method.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable


SOURCE_FILE = (
    Path(__file__).resolve().parent.parent
    / "config"
    / "isomer_literature_sources.json"
)
PUBMED_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_FORMULA_RE = re.compile(r"^(?:[A-Z][a-z]?\d*){2,}$")
_LIBRARY_SUFFIX_RE = re.compile(
    r"\s*\((?:known\s+isomers|isobaric\s+peaks).*?\)\s*$",
    re.IGNORECASE,
)
_CCS_SUFFIX_RE = re.compile(r"\s*\[\s*CCS\s*=.*?\]\s*$", re.IGNORECASE)
_UNRESOLVED_RE = re.compile(
    r"^(?:unresolved\s+compound|unknown|unnamed|compound[_\s-]*\d*)",
    re.IGNORECASE,
)
_NCBI_RATE_LOCK = threading.Lock()
_NCBI_LAST_REQUEST = 0.0


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _one_line(value: Any) -> str:
    return " ".join(str(value or "").split())


def normalize_literature_name(value: Any) -> str:
    text = _one_line(value)
    text = _LIBRARY_SUFFIX_RE.sub("", text)
    text = _CCS_SUFFIX_RE.sub("", text)
    return text.strip(" ,;")


def _usable_names(values: Iterable[Any], *, limit: int = 4) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    for value in values:
        name = normalize_literature_name(value)
        folded = name.casefold()
        if (
            not name
            or len(name) < 3
            or len(name) > 120
            or _UNRESOLVED_RE.match(name)
            or name.count("_") >= 2
            or folded in seen
        ):
            continue
        seen.add(folded)
        names.append(name)
        if len(names) >= limit:
            break
    return names


@lru_cache(maxsize=1)
def load_curated_sources() -> dict[str, Any]:
    try:
        payload = json.loads(SOURCE_FILE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"schema_version": 1, "policy": {}, "sources": []}
    if not isinstance(payload, dict):
        return {"schema_version": 1, "policy": {}, "sources": []}
    payload.setdefault("policy", {})
    payload.setdefault("sources", [])
    return payload


def _source_applies(source: dict[str, Any], ambiguity_type: str) -> bool:
    applies = {
        str(value)
        for value in (source.get("applies_to") or [])
    }
    return (
        "all_isomer_groups" in applies
        or ambiguity_type in applies
    )


def _pair_source_matches(
    source: dict[str, Any],
    candidate_names: list[str],
) -> bool:
    pair = _usable_names(source.get("pair_names") or [], limit=8)
    if not pair:
        return False
    available = {name.casefold() for name in candidate_names}
    return all(name.casefold() in available for name in pair)


def _pubmed_query_terms(
    candidate_names: list[str],
    formula: str = "",
) -> list[str]:
    method = (
        '("mass spectrometry"[Title/Abstract] OR '
        '"LC-MS"[Title/Abstract] OR "MS/MS"[Title/Abstract] OR '
        'chromatography[Title/Abstract] OR '
        '"ion mobility"[Title/Abstract])'
    )
    queries: list[str] = []
    if len(candidate_names) >= 2:
        pair = " AND ".join(
            f'"{name.replace(chr(34), "")}"[Title/Abstract]'
            for name in candidate_names[:2]
        )
        queries.append(f"({pair}) AND {method}")
    if candidate_names:
        selected = candidate_names[0].replace('"', "")
        queries.append(
            f'("{selected}"[Title/Abstract]) AND '
            f"(isomer*[Title/Abstract] OR stereoisomer*[Title/Abstract] "
            f"OR separation[Title/Abstract] OR "
            f"differentiation[Title/Abstract]) AND {method}"
        )
    formula = str(formula or "").strip()
    if _FORMULA_RE.fullmatch(formula):
        queries.append(
            f'("{formula}"[Title/Abstract]) AND '
            f"(isomer*[Title/Abstract] OR separation[Title/Abstract]) "
            f"AND {method}"
        )
    return list(dict.fromkeys(queries))


def pubmed_search_url(
    candidate_names: Iterable[Any],
    formula: str = "",
) -> str:
    names = _usable_names(candidate_names)
    queries = _pubmed_query_terms(names, formula)
    query = queries[0] if queries else (
        '"mass spectrometry"[Title/Abstract] AND isomer*[Title/Abstract]'
    )
    return "https://pubmed.ncbi.nlm.nih.gov/?" + urllib.parse.urlencode({
        "term": query,
    })


def build_isomer_literature_context(
    selected: dict[str, Any],
    alternatives: Iterable[dict[str, Any]],
    *,
    ambiguity_type: str,
    variant_subtype_counts: dict[str, int] | None = None,
    isomer_detail: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return local sources, explicit gaps and a reproducible PubMed query."""
    alternatives = [
        item for item in alternatives if isinstance(item, dict)
    ]
    names = _usable_names([
        selected.get("compound_name"),
        *[item.get("compound_name") for item in alternatives],
    ])
    formula = str(selected.get("formula") or "").strip()
    catalog = load_curated_sources()
    variant_subtype_counts = variant_subtype_counts or {}
    sources = [
        dict(source)
        for source in catalog.get("sources") or []
        if isinstance(source, dict)
        and _source_applies(source, ambiguity_type)
    ]
    has_structural_isomer_competition = bool(
        ambiguity_type in {
            "constitutional_or_positional_isomer_group",
            "mixed_isomer_group",
        }
        or variant_subtype_counts.get(
            "stereochemical_isotopic_or_tautomeric_variant"
        )
        or variant_subtype_counts.get(
            "same_connectivity_variant_unresolved"
        )
    )
    if not has_structural_isomer_competition:
        sources = [
            source for source in sources
            if "ion_mobility_method_lead"
            not in (source.get("usable_for") or [])
        ]
    pair_sources = [
        source for source in sources
        if _pair_source_matches(source, names)
    ]
    methodological_sources = [
        source for source in sources
        if not source.get("pair_names")
    ]
    # Keep the row payload compact while always retaining the strongest
    # pair-specific evidence and at least one reporting/method source.
    selected_sources = [
        *pair_sources,
        *methodological_sources[:3],
    ]
    seen_ids: set[str] = set()
    selected_sources = [
        source for source in selected_sources
        if not (
            str(source.get("source_id") or "") in seen_ids
            or seen_ids.add(str(source.get("source_id") or ""))
        )
    ]

    detail = isomer_detail or {}
    missing: list[str] = []
    replication = detail.get("reference_replication") or {}
    if replication.get("status") != "complete":
        missing.append("replicated_condition-matched_reference_spectra")
    if detail.get("metadata_comparability") != "complete":
        missing.append("comparable_adduct_collision_energy_and_instrument")
    if not int(detail.get("independent_supported_features") or 0):
        missing.append("replicated_candidate-specific_product_ions_or_losses")
    missing.extend([
        "same-method_authenticated_standard_retention_time",
        "validated_orthogonal_CCS_or_stereoselective_separation",
    ])
    if variant_subtype_counts.get(
        "structure_definition_specificity_variant"
    ):
        missing.append("complete_stereochemical_or_isotopic_structure_definition")
    if variant_subtype_counts.get(
        "protonation_or_charge_layer_variant"
    ):
        missing.append("standardized_parent_structure_and_ionization_state")

    pair_available = bool(pair_sources)
    return {
        "schema_version": 1,
        "status": (
            "pair_specific_method_evidence_available"
            if pair_available else "methodological_evidence_only"
        ),
        "candidate_names": names,
        "formula": formula,
        "sources": selected_sources,
        "pair_specific_source_count": len(pair_sources),
        "methodological_source_count": len(methodological_sources),
        "missing_experimental_information": list(dict.fromkeys(missing)),
        "pubmed_lookup_required": not pair_available,
        "pubmed_query_url": pubmed_search_url(names, formula),
        "runtime_network_access": False,
        "used_for_scoring": False,
        "used_for_ranking": False,
        "eligible_for_exact_identity_promotion": False,
        "interpretation": (
            "Curated literature defines method leads and evidence limits only. "
            "It is not a measured library spectrum or an authenticated "
            "same-method standard."
        ),
    }


def _request_bytes(url: str, timeout: float) -> bytes:
    global _NCBI_LAST_REQUEST
    email = str(os.getenv("MSSYSTEM_NCBI_EMAIL") or "").strip()
    api_key = str(os.getenv("NCBI_API_KEY") or "").strip()
    separator = "&" if "?" in url else "?"
    params = {"tool": "MSsystem_V7"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    request = urllib.request.Request(
        url + separator + urllib.parse.urlencode(params),
        headers={
            "User-Agent": (
                "MSsystem-V7/1.0 isomer-literature-audit "
                "(NCBI E-utilities)"
            ),
            "Accept": "application/json, application/xml, text/xml",
        },
    )
    # NCBI permits 3 requests/s without an API key and 10 requests/s with one.
    # Serialize this small optional channel so UI and background enrichment
    # cannot accidentally exceed that policy together.
    minimum_interval = 0.11 if api_key else 0.34
    with _NCBI_RATE_LOCK:
        delay = minimum_interval - (time.monotonic() - _NCBI_LAST_REQUEST)
        if delay > 0:
            time.sleep(delay)
        _NCBI_LAST_REQUEST = time.monotonic()
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def _esearch(query: str, max_results: int, timeout: float) -> list[str]:
    url = PUBMED_EUTILS + "/esearch.fcgi?" + urllib.parse.urlencode({
        "db": "pubmed",
        "retmode": "json",
        "retmax": max_results,
        "sort": "relevance",
        "term": query,
    })
    payload = json.loads(_request_bytes(url, timeout).decode("utf-8"))
    ids = ((payload.get("esearchresult") or {}).get("idlist") or [])
    return [str(value) for value in ids if str(value).isdigit()]


def _article_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return _one_line("".join(element.itertext()))


def _parse_pubmed_articles(xml_payload: bytes) -> list[dict[str, Any]]:
    root = ET.fromstring(xml_payload)
    records: list[dict[str, Any]] = []
    for article in root.findall(".//PubmedArticle"):
        pmid = _article_text(article.find("./MedlineCitation/PMID"))
        article_node = article.find("./MedlineCitation/Article")
        if article_node is None or not pmid:
            continue
        title = _article_text(article_node.find("./ArticleTitle"))
        abstract = " ".join(
            _article_text(node)
            for node in article_node.findall("./Abstract/AbstractText")
        ).strip()
        journal = _article_text(article_node.find("./Journal/Title"))
        year = _article_text(
            article_node.find("./Journal/JournalIssue/PubDate/Year")
        )
        if not year:
            year = _article_text(
                article_node.find("./Journal/JournalIssue/PubDate/MedlineDate")
            )[:4]
        ids: dict[str, str] = {}
        for identifier in article.findall(
            "./PubmedData/ArticleIdList/ArticleId"
        ):
            kind = str(identifier.attrib.get("IdType") or "")
            value = _article_text(identifier)
            if kind and value:
                ids[kind] = value
        records.append({
            "pmid": pmid,
            "pmcid": ids.get("pmc", ""),
            "doi": ids.get("doi", ""),
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "title": title,
            "journal": journal,
            "year": year,
            "_abstract_for_scope_only": abstract,
        })
    return records


def _efetch(pmids: list[str], timeout: float) -> list[dict[str, Any]]:
    if not pmids:
        return []
    url = PUBMED_EUTILS + "/efetch.fcgi?" + urllib.parse.urlencode({
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
    })
    return _parse_pubmed_articles(_request_bytes(url, timeout))


def _scope_pubmed_record(
    record: dict[str, Any],
    candidate_names: list[str],
) -> dict[str, Any]:
    abstract = str(record.pop("_abstract_for_scope_only", "") or "")
    corpus = f"{record.get('title', '')} {abstract}".casefold()
    # Match longer/nested names first and mask their spans.  Otherwise
    # "phellodendrine" would be counted as a second candidate merely because
    # the text contains only "phellodendrine chloride"; the same trap applies
    # to ephedrine/pseudoephedrine.
    remaining = corpus
    name_hits: list[str] = []
    for name in sorted(candidate_names, key=len, reverse=True):
        pattern = re.compile(
            rf"(?<![a-z0-9]){re.escape(name.casefold())}(?![a-z0-9])"
        )
        if pattern.search(remaining):
            name_hits.append(name)
            remaining = pattern.sub(
                lambda match: " " * len(match.group(0)),
                remaining,
            )
    name_hits = [
        name for name in candidate_names if name in set(name_hits)
    ]
    tag_terms = {
        "tandem_ms": ("ms/ms", "tandem mass"),
        "chromatographic_separation": (
            "chromatograph", "stationary phase", "column"
        ),
        "chiral_separation": ("chiral", "enantio", "diastereo"),
        "ion_mobility_or_ccs": (
            "ion mobility", "collision cross section", "ccs"
        ),
        "retention_time": ("retention time",),
        "diagnostic_fragment": (
            "diagnostic ion", "diagnostic fragment", "product ion"
        ),
        "authenticated_standard": (
            "authentic standard", "reference standard"
        ),
    }
    tags = [
        tag for tag, terms in tag_terms.items()
        if any(term in corpus for term in terms)
    ]
    # A one-name group can only yield a compound-specific/general lead. It
    # must never be labelled pair-specific merely because its sole readable
    # name was found in the paper.
    pair_specific = len(candidate_names) >= 2 and len(name_hits) >= 2
    record.update({
        "candidate_name_hits": name_hits,
        "evidence_tags": tags,
        "evidence_scope": (
            "pair_specific_bibliographic_method_lead"
            if pair_specific
            else "single_candidate_or_general_bibliographic_lead"
        ),
        "abstract_screened_for_scope": bool(abstract),
        "diagnostic_fragment_data_extracted": False,
        "used_for_scoring": False,
        "used_for_ranking": False,
        "eligible_for_exact_identity_promotion": False,
        "full_text_and_experimental_conditions_required": True,
    })
    return record


def _cache_file(
    cache_dir: Path | None,
    payload: dict[str, Any],
) -> Path | None:
    if cache_dir is None:
        return None
    digest = hashlib.sha256(
        json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    return Path(cache_dir) / f"{digest}.json"


def _cache_is_fresh(payload: dict[str, Any], days: int = 30) -> bool:
    try:
        searched_at = datetime.fromisoformat(
            str(payload.get("searched_at") or "")
        )
        if searched_at.tzinfo is None:
            searched_at = searched_at.replace(tzinfo=timezone.utc)
    except (TypeError, ValueError):
        return False
    age_seconds = (
        datetime.now(timezone.utc) - searched_at.astimezone(timezone.utc)
    ).total_seconds()
    return 0 <= age_seconds <= days * 86400


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        temporary = Path(handle.name)
    temporary.replace(path)


def search_pubmed_isomer_evidence(
    candidate_names: Iterable[Any],
    *,
    formula: str = "",
    max_results: int = 5,
    timeout: float = 12.0,
    cache_dir: Path | None = None,
) -> dict[str, Any]:
    """Search PubMed without treating bibliographic hits as MS/MS evidence."""
    names = _usable_names(candidate_names)
    formula = str(formula or "").strip()
    max_results = max(1, min(int(max_results), 10))
    timeout = max(2.0, min(float(timeout), 30.0))
    queries = _pubmed_query_terms(names, formula)
    request_key = {
        "schema_version": 3,
        "candidate_names": names,
        "formula": formula,
        "max_results": max_results,
        "queries": queries,
    }
    cache_path = _cache_file(cache_dir, request_key)
    if cache_path and cache_path.exists():
        try:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if isinstance(cached, dict) and _cache_is_fresh(cached):
                cached["cache_hit"] = True
                return cached
        except (OSError, json.JSONDecodeError):
            pass

    base = {
        **request_key,
        "searched_at": _now_iso(),
        "search_provider": "NCBI PubMed E-utilities",
        "search_url": pubmed_search_url(names, formula),
        "attempts": [],
        "records": [],
        "cache_hit": False,
        "used_for_scoring": False,
        "used_for_ranking": False,
        "eligible_for_exact_identity_promotion": False,
        "interpretation": (
            "Results are bibliographic leads only. Full text, exact structures "
            "and experimental conditions must be checked before manual "
            "curation; no reported ion is imported as a measured spectrum."
        ),
    }
    if not queries:
        base["status"] = "search_terms_insufficient"
        base["error"] = ""
        return base

    try:
        pmids: list[str] = []
        selected_query = ""
        for query in queries:
            ids = _esearch(query, max_results, timeout)
            base["attempts"].append({
                "query": query,
                "result_count_returned": len(ids),
            })
            if ids:
                pmids = ids
                selected_query = query
                break
        records = [
            _scope_pubmed_record(record, names)
            for record in _efetch(pmids, timeout)
        ]
        base.update({
            "status": (
                "bibliographic_leads_found"
                if records else "searched_no_bibliographic_leads"
            ),
            "selected_query": selected_query,
            "records": records,
            "error": "",
        })
    except (
        OSError,
        TimeoutError,
        urllib.error.URLError,
        urllib.error.HTTPError,
        json.JSONDecodeError,
        ET.ParseError,
        ValueError,
    ) as exc:
        base.update({
            "status": "lookup_temporarily_unavailable",
            "error": f"{type(exc).__name__}: {_one_line(exc)}"[:240],
        })
    if cache_path and base.get("status") != "lookup_temporarily_unavailable":
        try:
            write_json_atomic(cache_path, base)
        except OSError:
            pass
    return base


def collect_isomer_literature_groups(
    results: Iterable[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Deduplicate current, spectrally competitive ambiguity groups."""
    # Local import avoids a module-import cycle: honest_calibration calls the
    # lightweight context builder above.
    from matching.honest_calibration import annotate

    groups: dict[str, dict[str, Any]] = {}
    for row_index, row in enumerate(results):
        if not isinstance(row, dict):
            continue
        match = row.get("match")
        if not isinstance(match, dict):
            continue
        annotation = annotate(match, match.get("top_candidates"))
        if not annotation.get("isomer_ambiguous"):
            continue
        literature = annotation.get("isomer_literature_evidence") or {}
        names = literature.get("candidate_names") or _usable_names([
            match.get("compound_name"),
            *[
                item.get("compound_name")
                for item in annotation.get("isomer_set") or []
                if isinstance(item, dict)
            ],
        ])
        inchikeys = sorted({
            str(value or "").strip()
            for value in [
                match.get("inchikey"),
                *[
                    item.get("inchikey")
                    for item in annotation.get("isomer_set") or []
                    if isinstance(item, dict)
                ],
            ]
            if str(value or "").strip()
        })
        formula = str(
            annotation.get("isomer_group_formula")
            or match.get("formula")
            or ""
        ).strip()
        identity = {
            "formula": formula,
            "inchikeys": inchikeys,
            "candidate_names": names,
        }
        group_key = hashlib.sha256(
            json.dumps(
                identity,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()[:20]
        group = groups.setdefault(group_key, {
            "group_key": group_key,
            **identity,
            "bibliographic_name_scope": (
                "multiple_distinct_candidate_names_available"
                if len(names) >= 2
                else "single_distinct_candidate_name_only"
            ),
            "bibliographic_scope_warning": (
                ""
                if len(names) >= 2
                else (
                    "Only one distinct readable candidate name is available "
                    "for this ambiguity group. PubMed hits are single-compound "
                    "or general method leads, not candidate-pair evidence."
                )
            ),
            "ambiguity_type": annotation.get(
                "isomer_ambiguity_type", "none"
            ),
            "variant_subtype_counts": annotation.get(
                "isomer_variant_subtype_counts", {}
            ),
            "literature_evidence": literature,
            "result_row_indices": [],
        })
        group["result_row_indices"].append(row_index)
    return list(groups.values())


def enrich_isomer_groups_with_pubmed(
    results: Iterable[dict[str, Any]],
    *,
    cache_dir: Path,
    max_groups: int = 0,
    max_results: int = 5,
    timeout: float = 12.0,
    progress_callback=None,
) -> dict[str, Any]:
    """Build a task-level literature artifact without touching match scores."""
    groups = collect_isomer_literature_groups(results)
    total_groups = len(groups)
    selected_groups = (
        groups[:max_groups]
        if int(max_groups or 0) > 0 else groups
    )
    enriched: list[dict[str, Any]] = []
    for index, group in enumerate(selected_groups, 1):
        local = group.get("literature_evidence") or {}
        if local.get("pubmed_lookup_required") is False:
            lookup = {
                "schema_version": 3,
                "status": "curated_pair_specific_source_available",
                "records": [],
                "search_url": local.get("pubmed_query_url", ""),
                "used_for_scoring": False,
                "used_for_ranking": False,
                "eligible_for_exact_identity_promotion": False,
                "interpretation": (
                    "The curated pair-specific source is retained with its "
                    "method scope; no network text was converted to spectra."
                ),
            }
        else:
            lookup = search_pubmed_isomer_evidence(
                group.get("candidate_names") or [],
                formula=group.get("formula") or "",
                max_results=max_results,
                timeout=timeout,
                cache_dir=cache_dir,
            )
        enriched.append({
            **group,
            "pubmed_lookup": lookup,
        })
        if progress_callback is not None:
            progress_callback(index, len(selected_groups), group, lookup)

    status_counts: dict[str, int] = {}
    for group in enriched:
        status = str(
            (group.get("pubmed_lookup") or {}).get("status")
            or "unknown"
        )
        status_counts[status] = status_counts.get(status, 0) + 1
    return {
        "schema_version": 1,
        "generated_at": _now_iso(),
        "search_provider": "NCBI PubMed E-utilities",
        "total_unique_ambiguous_groups": total_groups,
        "groups_processed": len(enriched),
        "groups_deferred_by_limit": max(0, total_groups - len(enriched)),
        "status_counts": dict(sorted(status_counts.items())),
        "used_for_scoring": False,
        "used_for_ranking": False,
        "eligible_for_exact_identity_promotion": False,
        "groups": enriched,
    }
