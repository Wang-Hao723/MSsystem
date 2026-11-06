"""MSsystem V11 flagship entity ranker.

The ranker has one and only one identity-selection path:

1. construct a candidate-independent precursor entity;
2. score library candidates from measured-spectrum and precursor-mass fields;
3. collapse technical repeats inside one collision-energy bin;
4. aggregate support across distinct collision-energy bins;
5. veto only candidate-independent formula/adduct/precursor contradictions; and
6. emit one bounded-scope entity decision.

Expert rules, predicted classes, biosynthetic labels, legacy confidence tiers and
candidate pathway annotations are intentionally inaccessible to the scoring
functions below.  They may still be displayed by downstream consumers, but they
cannot select an identity in V11.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import re
import statistics
import unicodedata
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from functools import lru_cache
from pathlib import Path
from typing import Any

from matching.conditioned_multilibrary import aggregate_references
from matching.library_quality import formula_exact_mass, normalize_formula_text
from utils.adduct import adduct_polarity, ion_mz


FLAGSHIP_SCHEMA_VERSION = "1.0.0"
MODULE_ID = "v11_flagship_chem_entity_ranker"
_FULL_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_CONNECTIVITY_BLOCK = re.compile(r"^[A-Z]{14}$")

try:  # RDKit is optional for a conservative identity fallback.
    from rdkit import Chem as _Chem
    from rdkit import rdBase as _rdBase
    from rdkit.Chem.MolStandardize import rdMolStandardize as _rdMolStandardize
except Exception:  # pragma: no cover - exercised in minimal deployments
    _Chem = None
    _rdBase = None
    _rdMolStandardize = None


def default_flagship_policy_path() -> Path:
    return Path(__file__).resolve().parent.parent / "config" / "v11_flagship_policy.json"


def default_identity_registry_path() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "config"
        / "canonical_name_identity_registry.json"
    )


def load_flagship_policy(path: Path | None = None) -> dict[str, Any]:
    policy_path = path or default_flagship_policy_path()
    payload = json.loads(policy_path.read_text(encoding="utf-8"))
    if payload.get("deployment_state") != "production":
        raise ValueError("V11 flagship policy must declare production deployment")
    forbidden = set(payload.get("forbidden_ranking_inputs") or [])
    required = {
        "expert_score",
        "classification",
        "biosynthetic_pathway",
        "candidate_pathway_alignment",
        "legacy_confidence_label",
    }
    if not required <= forbidden:
        raise ValueError("V11 flagship policy does not forbid every auxiliary rank input")
    registry_policy = payload.get("identity_registry") or {}
    if registry_policy.get("required"):
        registry_path = default_identity_registry_path()
        if not registry_path.is_file():
            raise ValueError("V11 canonical identity registry is required but missing")
        actual = hashlib.sha256(registry_path.read_bytes()).hexdigest()
        if actual != str(registry_policy.get("sha256") or ""):
            raise ValueError("V11 canonical identity registry hash does not match policy")
    return payload


def flagship_policy_sha256(policy: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical(policy)).hexdigest()


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def _text(value: Any) -> str:
    return str(value or "").strip()


def _strict_name_key(value: Any) -> str:
    return " ".join(
        unicodedata.normalize("NFKC", _text(value)).casefold().split()
    )


@lru_cache(maxsize=8)
def _identity_registry(expected_sha256: str) -> tuple[str, dict[tuple[str, str], dict[str, Any]]]:
    path = default_identity_registry_path()
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected_sha256:
        raise ValueError("canonical identity registry changed after policy freeze")
    payload = json.loads(raw)
    entries = {
        (
            _strict_name_key(entry.get("name_key") or entry.get("query_name")),
            normalize_formula_text(entry.get("formula")),
        ): dict(entry)
        for entry in (payload.get("entries") or [])
        if _strict_name_key(entry.get("name_key") or entry.get("query_name"))
        and normalize_formula_text(entry.get("formula"))
    }
    return str(payload.get("registry_id") or ""), entries


def _registered_identity(
    candidate: Mapping[str, Any], policy: Mapping[str, Any]
) -> tuple[str, dict[str, Any] | None]:
    registry_policy = policy.get("identity_registry") or {}
    if registry_policy.get("enabled") is False:
        return "", None
    expected = str(registry_policy.get("sha256") or "")
    if not expected:
        return "", None
    _, entries = _identity_registry(expected)
    key = (
        _strict_name_key(candidate.get("compound_name")),
        normalize_formula_text(candidate.get("formula")),
    )
    entry = entries.get(key)
    if not entry:
        return "", None
    block = _text(entry.get("connectivity_inchikey")).upper()
    if not _CONNECTIVITY_BLOCK.fullmatch(block):
        return "", None
    return block, entry


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _search_precursor_mz(row: Mapping[str, Any]) -> float | None:
    """Return the mass axis actually used for high-resolution retrieval."""

    search_mz = _finite(row.get("search_precursor_mz"))
    if search_mz is not None and search_mz > 0:
        return search_mz
    return _finite(row.get("precursor_mz"))


def _unit_interval(value: Any, *, missing: float = 0.0) -> float:
    number = _finite(value)
    if number is None:
        return missing
    return max(0.0, min(1.0, number))


def _first_present(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        value = mapping.get(key)
        if value is not None and value != "":
            return value
    return None


@lru_cache(maxsize=32768)
def _standardized_structure(text: str) -> tuple[str, str]:
    """Return a bounded tautomer-standardized 2-D SMILES and InChIKey block."""

    if not text:
        return "", ""
    if _Chem is None:
        return "", ""
    try:
        blocker = _rdBase.BlockLogs() if _rdBase is not None else None
        molecule = _Chem.MolFromSmiles(text)
        if molecule is None:
            return "", ""
        standardized = molecule
        if _rdMolStandardize is not None:
            enumerator = _rdMolStandardize.TautomerEnumerator()
            # The bound prevents pathological structures from dominating task
            # runtime while still resolving common keto-enol representations.
            enumerator.SetMaxTautomers(32)
            enumerator.SetMaxTransforms(100)
            standardized = enumerator.Canonicalize(molecule)
        canonical = _Chem.MolToSmiles(
            standardized, canonical=True, isomericSmiles=False
        )
        inchikey = _Chem.MolToInchiKey(standardized)
        block = inchikey.split("-", 1)[0] if _FULL_INCHIKEY.fullmatch(inchikey) else ""
        del blocker
        return canonical, block
    except Exception:
        return "", ""


def _canonical_nonisomeric_smiles(value: Any) -> str:
    return _standardized_structure(_text(value))[0]


def connectivity_identity(
    candidate: Mapping[str, Any], *, standardize_tautomers: bool = True
) -> tuple[str, str]:
    """Return a 2-D connectivity identity and its evidence source.

    A bounded RDKit canonical-tautomer transform is attempted before the
    reported InChIKey fallback.  Both paths deliberately collapse stereo layers,
    making it impossible for ordinary LC-MS/MS similarity to silently claim
    stereochemistry merely because library structure encodings differ.
    """

    if standardize_tautomers:
        standardized_smiles, standardized_block = _standardized_structure(
            _text(candidate.get("smiles"))
        )
        if standardized_block:
            return f"std2d:{standardized_block}", "tautomer_standardized_connectivity"
        if standardized_smiles:
            digest = hashlib.sha256(standardized_smiles.encode("utf-8")).hexdigest()[:24]
            return f"std2d-smiles:{digest}", "tautomer_standardized_nonisomeric_smiles"
    inchikey = _text(candidate.get("inchikey")).upper().replace("INCHIKEY=", "")
    if _FULL_INCHIKEY.fullmatch(inchikey):
        return f"ik14:{inchikey[:14]}", "inchikey_connectivity_block"
    first_block = inchikey.split("-", 1)[0]
    if _CONNECTIVITY_BLOCK.fullmatch(first_block):
        return f"ik14:{first_block}", "reported_inchikey_connectivity_block"
    formula = normalize_formula_text(candidate.get("formula"))
    name = re.sub(r"[^a-z0-9]+", "", _text(candidate.get("compound_name")).lower())
    fallback = f"{formula}|{name}" if formula or name else "unidentified"
    digest = hashlib.sha256(fallback.encode("utf-8")).hexdigest()[:24]
    return f"fallback:{digest}", "name_formula_fallback"


def _entity_key(task_id: str, row: Mapping[str, Any], row_index: int) -> tuple[str, str, str]:
    lcms = _text(row.get("lcms_feature_id"))
    precursor = _text(row.get("precursor_feature_id"))
    spectrum = _text(row.get("spectrum_id")) or f"row-{row_index}"
    if lcms:
        source, raw = "lcms_feature_id", lcms
    elif precursor:
        source, raw = "precursor_feature_id", precursor
    else:
        # No candidate identity is used to merge fallback scans.  In the
        # absence of an acquisition feature identifier, one row is one entity.
        source, raw = "spectrum_fallback", f"{spectrum}|{row_index}"
    digest = hashlib.sha256(f"{source}|{raw}".encode("utf-8")).hexdigest()[:18]
    return f"entity:{task_id}:{digest}", source, raw


def _collision_energy_bin(row: Mapping[str, Any], policy: Mapping[str, Any]) -> tuple[str, float | None, str]:
    value = _finite(row.get("collision_energy"))
    unit = _text(row.get("collision_energy_unit")).lower() or "unspecified"
    width = max(1e-9, float((policy.get("collision_energy") or {}).get("bin_width") or 0.5))
    if value is None:
        return f"unknown:{unit}", None, unit
    binned = round(value / width) * width
    return f"{binned:.6f}:{unit}", binned, unit


def _mass_accuracy(candidate: Mapping[str, Any]) -> tuple[float, str]:
    explicit = _finite(candidate.get("mass_score"))
    if explicit is not None:
        return max(0.0, min(1.0, explicit)), "mass_score"
    ppm = _finite(candidate.get("mass_ppm"))
    if ppm is None:
        # Retrieval itself already enforced a precursor window.  Missing the
        # displayed ppm value is not evidence against the candidate.
        return 0.5, "retrieval_only_missing_ppm"
    return math.exp(-0.5 * (abs(ppm) / 25.0) ** 2), "ppm_gaussian_25"


def _spectral_score(candidate: Mapping[str, Any], policy: Mapping[str, Any]) -> tuple[float, dict[str, Any]]:
    weights = policy.get("scan_score") or {}
    cosine = _unit_interval(
        _first_present(candidate, "raw_cosine", "candidate_comparison_cosine", "cosine")
    )
    entropy = _unit_interval(_first_present(candidate, "raw_entropy", "entropy"))
    reverse = _unit_interval(_first_present(candidate, "raw_reverse_match", "reverse_match"))
    fragment = _unit_interval(_first_present(candidate, "raw_fragment_score", "fragment_score"))
    intensity = _unit_interval(_first_present(candidate, "raw_intensity_score", "intensity_score"))
    mass, mass_source = _mass_accuracy(candidate)
    components = {
        "cosine": cosine,
        "entropy": entropy,
        "reverse_match": reverse,
        "mass_accuracy": mass,
        "fragment_coverage": fragment,
        "intensity_agreement": intensity,
    }
    score = sum(float(weights.get(name) or 0.0) * value for name, value in components.items())
    return max(0.0, min(1.0, score)), {
        **{key: round(value, 6) for key, value in components.items()},
        "mass_accuracy_source": mass_source,
        "expert_score_used": False,
        "classification_used": False,
        "pathway_used": False,
    }


_DISPLAY_FIELDS = (
    "compound_name",
    "compound_name_cn",
    "source_compound_name",
    "formula",
    "smiles",
    "inchikey",
    "reported_inchikey",
    "adduct",
    "has_structure",
    "structure_quality",
    "name_only_lead",
    "dataset_source",
    "source_tier",
    "lib_spectrum_id",
    "lib_source_accession",
    "classification",
    "classification_source",
    "classification_reliable",
    "classification_reliability_reason",
    "feature_fragment_evidence",
    "potential_contaminant",
    "contaminant_category",
    "compound_role",
)


def _display_candidate(
    candidate: Mapping[str, Any], *, standardize_tautomers: bool = True
) -> dict[str, Any]:
    display = {
        key: candidate.get(key)
        for key in _DISPLAY_FIELDS
        if candidate.get(key) is not None
    }
    standardized_smiles, standardized_block = (
        _standardized_structure(_text(candidate.get("smiles")))
        if standardize_tautomers else ("", "")
    )
    if standardized_smiles:
        display["standardized_nonisomeric_smiles"] = standardized_smiles
    if standardized_block:
        display["standardized_connectivity_inchikey"] = standardized_block
        reported = _text(candidate.get("inchikey")).upper().split("-", 1)[0]
        display["reported_connectivity_matches_standardized"] = (
            reported == standardized_block
        )
    display["structure_standardization"] = (
        "bounded_rdkit_tautomer_nonisomeric"
        if standardized_smiles and _rdMolStandardize is not None
        else "disabled_for_ablation"
        if not standardize_tautomers
        else "unavailable_or_raw_fallback"
    )
    return display


def _reference_source_family(candidate: Mapping[str, Any]) -> str:
    """Return a conservative dataset-level provenance family.

    Spectrum accessions are deliberately excluded: multiple scans from one
    deposited dataset are technical/library replication, not independent
    provenance.  Polarity and version suffixes are collapsed so they cannot
    manufacture extra source votes.
    """

    raw = _text(candidate.get("dataset_source")).lower()
    if not raw:
        return ""
    normalized = re.sub(r"[^a-z0-9]+", "-", raw).strip("-")
    normalized = re.sub(r"-(?:pos|neg|positive|negative)$", "", normalized)
    normalized = re.sub(r"-(?:v)?\d+(?:-\d+)*$", "", normalized)
    families = (
        ("spectraverse", "spectraverse"),
        ("massbank", "massbank"),
        ("bmdms", "bmdms"),
        ("vaniya", "vaniya-fiehn"),
        ("mona", "mona"),
        ("gnps", "gnps"),
        ("nist", "nist"),
    )
    return next(
        (family for marker, family in families if marker in normalized),
        normalized,
    )


def _merge_projected_references(records):
    """Same reference ID replaces its fallback only; other IDs always survive."""
    grouped = {}
    for index, record in enumerate(records):
        sid = record.get('lib_spectrum_id')
        key = ('id', str(sid)) if sid else ('unknown', index)
        grouped.setdefault(key, []).append(record)
    result = []
    for members in grouped.values():
        compact = [r for r in members if r.get('_reference_projection') == 'compact_search_evidence']
        if compact:
            # Keep original compact records; downstream fingerprint handling is unchanged.
            result.extend(compact)
        else:
            # Duplicate exported slots are not additional reference evidence.
            # Prefer an explicitly recorded mass score over a missing-field fallback.
            result.append(max(members, key=lambda r: (
                (r.get('score_components') or {}).get('mass_accuracy_source') == 'mass_score',
                bool(r.get('library_peak_fingerprint') or r.get('lib_peaks')),
                float(r.get('spectral_score') or 0),
            )))
    return result


def _restore_reference_mass_scores(row):
 new=copy.deepcopy(row);m=new.get('match') or {};cs=[m]+(m.get('top_candidates') or []);by=defaultdict(list);changes=[]
 for c in cs:
  for ref in c.get('reference_evidence_records') or []:
   sid=ref.get('lib_spectrum_id') or ref.get('spectrum_id');v=ref.get('mass_score')
   if sid and v is not None and math.isfinite(float(v)):by[sid].append(float(v))
 for i,c in enumerate(cs):
  sid=c.get('lib_spectrum_id');vals=by.get(sid,[])
  if c.get('mass_score') is None and vals and max(vals)-min(vals)<=1e-10:
   c['mass_score']=vals[0];changes.append({'position':i,'reference_id':sid,'restored_mass_score':vals[0],'old_mass_accuracy':_mass_accuracy(([row.get('match') or {}]+((row.get('match') or {}).get('top_candidates') or []))[i])[0]})
 return new,changes


def _candidate_rows(row: Mapping[str, Any], policy: Mapping[str, Any]) -> list[dict[str, Any]]:
    row, _ = _restore_reference_mass_scores(row)
    match = row.get("match") if isinstance(row.get("match"), Mapping) else None
    if not match:
        return []
    candidates: list[Mapping[str, Any]] = [match]
    candidates.extend(
        candidate for candidate in (match.get("top_candidates") or [])
        if isinstance(candidate, Mapping)
    )
    maximum = max(1, int((policy.get("candidate_pool") or {}).get("maximum_candidates_per_scan") or 20))
    standardize_tautomers = bool(
        ((policy.get("candidate_pool") or {}).get("tautomer_standardization") or {}).get(
            "enabled", True
        )
    )
    by_identity: dict[str, dict[str, Any]] = {}
    references_by_identity: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for candidate in candidates:
        registry_block, registry_entry = _registered_identity(candidate, policy)
        resolved_candidate = dict(candidate)
        if registry_entry is not None:
            reported_inchikey = _text(candidate.get("inchikey"))
            reported_smiles = _text(candidate.get("smiles"))
            resolved_candidate["inchikey"] = registry_entry.get("inchikey")
            resolved_candidate["smiles"] = registry_entry.get("isomeric_smiles")
            identity = f"std2d:{registry_block}"
            identity_source = "pubchem_exact_name_formula_registry"
        else:
            reported_inchikey = ""
            reported_smiles = ""
            identity, identity_source = connectivity_identity(
                candidate, standardize_tautomers=standardize_tautomers
            )
        spectral, components = _spectral_score(candidate, policy)
        compact_references = [
            reference
            for reference in (candidate.get("reference_evidence_records") or [])
            if isinstance(reference, Mapping)
        ]
        reference_inputs = compact_references + [candidate]
        for reference in reference_inputs:
            reference_spectral, reference_components = _spectral_score(
                reference, policy
            )
            references_by_identity[identity].append({
                # Names and biological annotations are intentionally not
                # projected into the V15 reference aggregator.
                "spectral_score": reference_spectral,
                "score_components": reference_components,
                "dataset_source": reference.get("dataset_source"),
                "source_tier": reference.get("source_tier"),
                "lib_spectrum_id": (
                    reference.get("lib_spectrum_id")
                    or reference.get("spectrum_id")
                ),
                "lib_source_accession": (
                    reference.get("lib_source_accession")
                    or reference.get("source_accession")
                ),
                "library_peak_fingerprint": reference.get(
                    "library_peak_fingerprint"
                ),
                "lib_peaks": (
                    reference.get("lib_peaks")
                    or reference.get("_score_library_peaks")
                    or (
                        row.get("lib_peaks")
                        if reference is candidate and candidate is match
                        else []
                    )
                ),
                "lib_collision_energy": reference.get(
                    "lib_collision_energy"
                ),
                "lib_collision_energy_unit": reference.get(
                    "lib_collision_energy_unit"
                ),
                "lib_collision_energy_mode": reference.get(
                    "lib_collision_energy_mode"
                ),
                "lib_instrument_type": reference.get("lib_instrument_type"),
                "lib_instrument": reference.get("lib_instrument"),
                "lib_activation_method": reference.get(
                    "lib_activation_method"
                ),
                "lib_instrument_domain_id": reference.get(
                    "lib_instrument_domain_id"
                ),
                "lib_collision_cell_region": reference.get(
                    "lib_collision_cell_region"
                ),
                "_reference_projection": (
                    "compact_search_evidence"
                    if reference is not candidate
                    else "candidate_fallback"
                ),
            })
        display = _display_candidate(
            resolved_candidate,
            standardize_tautomers=standardize_tautomers,
        )
        if registry_entry is not None:
            ranker_standardized_before_resolution = _text(
                display.get("standardized_connectivity_inchikey")
            )
            display["standardized_connectivity_inchikey"] = registry_block
            display["reported_connectivity_matches_standardized"] = bool(
                reported_inchikey
                and reported_inchikey.upper().split("-", 1)[0] == registry_block
            )
            if (
                ranker_standardized_before_resolution
                and ranker_standardized_before_resolution != registry_block
            ):
                display[
                    "ranker_standardized_connectivity_before_identity_resolution"
                ] = ranker_standardized_before_resolution
            display.update({
                "reported_inchikey_before_identity_resolution": reported_inchikey or None,
                "reported_smiles_before_identity_resolution": reported_smiles or None,
                "canonical_identity_resolution": {
                    "status": "resolved",
                    "registry_id": _identity_registry(
                        str((policy.get("identity_registry") or {}).get("sha256") or "")
                    )[0],
                    "source": registry_entry.get("identity_source"),
                    "pubchem_cid": registry_entry.get("pubchem_cid"),
                    "pubchem_title": registry_entry.get("pubchem_title"),
                    "exact_normalized_name_and_formula": True,
                    "ground_truth_inputs_used": False,
                },
            })
        item = {
            "identity_key": identity,
            "identity_source": identity_source,
            "candidate": display,
            "spectral_score": spectral,
            "components": components,
            "reference_source_families": sorted(filter(None, {
                _reference_source_family(candidate)
            })),
        }
        previous = by_identity.get(identity)
        if previous is None:
            by_identity[identity] = item
            continue
        previous_display = dict(previous["candidate"])
        item_display = dict(item["candidate"])
        merged_sources = sorted(set(previous.get("reference_source_families") or []) | set(
            item.get("reference_source_families") or []
        ))
        preferred_display = min(
            (
                (previous_display, float(previous["spectral_score"])),
                (item_display, spectral),
            ),
            key=lambda option: (
                -int(bool(option[0].get("canonical_identity_resolution"))),
                -int(bool(option[0].get("reported_connectivity_matches_standardized"))),
                -int(bool(option[0].get("has_structure"))),
                -float(option[1]),
                _text(option[0].get("inchikey")),
                _text(option[0].get("compound_name")),
            ),
        )[0]
        if (
            -spectral,
            _text(item_display.get("inchikey")),
            _text(item_display.get("compound_name")),
        ) < (
            -float(previous["spectral_score"]),
            _text(previous_display.get("inchikey")),
            _text(previous_display.get("compound_name")),
        ):
            item["candidate"] = preferred_display
            item["reference_source_families"] = merged_sources
            by_identity[identity] = item
        else:
            previous["candidate"] = preferred_display
            previous["reference_source_families"] = merged_sources

    for identity, records in references_by_identity.items():
        references_by_identity[identity] = _merge_projected_references(records)

    multilibrary_policy = policy.get("conditioned_multilibrary") or {}
    multilibrary_enabled = bool(multilibrary_policy.get("enabled"))
    if multilibrary_enabled:
        for identity, item in by_identity.items():
            evidence = aggregate_references(
                row,
                references_by_identity[identity],
                multilibrary_policy,
            )
            legacy_score = float(item["spectral_score"])
            item["legacy_spectral_score"] = round(legacy_score, 8)
            item["spectral_score"] = float(
                evidence["conditioned_multilibrary_score"]
            )
            item["conditioned_multilibrary_evidence"] = evidence
            item["reference_source_families"] = list(
                evidence.get("source_families") or []
            )

        legacy_ordered = sorted(
            by_identity.values(),
            key=lambda item: (
                -float(item.get("legacy_spectral_score") or 0.0),
                item["identity_key"],
            ),
        )[:maximum]
        rank_fraction = float(
            (policy.get("scan_score") or {}).get("reciprocal_rank_fraction")
            or 0.30
        )
        spectral_fraction = float(
            (policy.get("scan_score") or {}).get("spectral_fraction") or 0.70
        )
        for rank, item in enumerate(legacy_ordered, start=1):
            reciprocal = 1.0 / math.log2(rank + 1.0)
            item["legacy_scan_rank"] = rank
            item["legacy_scan_support"] = round(
                max(
                    0.0,
                    min(
                        1.0,
                        spectral_fraction
                        * float(item.get("legacy_spectral_score") or 0.0)
                        + rank_fraction * reciprocal,
                    ),
                ),
                6,
            )
    ordered = sorted(
        by_identity.values(),
        key=lambda item: (-float(item["spectral_score"]), item["identity_key"]),
    )[:maximum]
    rank_fraction = float((policy.get("scan_score") or {}).get("reciprocal_rank_fraction") or 0.30)
    spectral_fraction = float((policy.get("scan_score") or {}).get("spectral_fraction") or 0.70)
    for rank, item in enumerate(ordered, start=1):
        reciprocal = 1.0 / math.log2(rank + 1.0)
        item["scan_rank"] = rank
        item["reciprocal_rank_support"] = round(reciprocal, 6)
        item["scan_support"] = round(
            max(0.0, min(1.0, spectral_fraction * float(item["spectral_score"]) + rank_fraction * reciprocal)),
            6,
        )
        if not multilibrary_enabled:
            item["legacy_spectral_score"] = round(
                float(item["spectral_score"]), 8
            )
            item["legacy_scan_rank"] = rank
            item["legacy_scan_support"] = item["scan_support"]
    return ordered


def _compact_row(row: Mapping[str, Any], row_index: int) -> dict[str, Any]:
    match = row.get("match") if isinstance(row.get("match"), Mapping) else {}
    precursor_mass_audit = (
        row.get("precursor_mass_audit")
        if isinstance(row.get("precursor_mass_audit"), Mapping)
        else {}
    )
    reliability_value = _first_present(
        row,
        "search_precursor_mz_reliable",
    )
    reliability_source = "result_row"
    if reliability_value is None:
        reliability_value = precursor_mass_audit.get(
            "search_precursor_mz_reliable"
        )
        reliability_source = (
            "precursor_mass_audit"
            if reliability_value is not None
            else "legacy_unspecified"
        )
    chimeric = (
        match.get("chimeric_deconvolution")
        if isinstance(match.get("chimeric_deconvolution"), Mapping)
        else {}
    )
    query_peaks = row.get("query_peaks") or []
    raw_peaks = row.get("raw_peaks") or []
    quality_peaks = query_peaks or raw_peaks
    return {
        "result_index": row_index,
        "spectrum_id": _text(row.get("spectrum_id")) or str(row_index),
        "precursor_mz": _finite(row.get("precursor_mz")),
        "search_precursor_mz": _search_precursor_mz(row),
        "search_precursor_mz_source": _text(
            row.get("search_precursor_mz_source")
            or "ms2_selected_ion_mz_fallback"
        ),
        "search_precursor_mz_reliable": (
            bool(reliability_value)
            if reliability_value is not None
            else True
        ),
        "search_precursor_mz_reliability_status": (
            "reliable"
            if reliability_value is True
            else "unreliable"
            if reliability_value is False
            else "legacy_unspecified_not_failed"
        ),
        "search_precursor_mz_reliability_source": reliability_source,
        "search_precursor_delta_ppm": _finite(
            row.get("search_precursor_delta_ppm")
        ),
        "lcms_feature_used_for_search": bool(
            row.get("lcms_feature_used_for_search")
        ),
        "precursor_charge": row.get("precursor_charge"),
        "precursor_intensity": _finite(row.get("precursor_intensity")),
        "tic": _finite(row.get("tic")),
        "retention_time": _finite(row.get("retention_time")),
        "retention_time_unit": _text(row.get("retention_time_unit")),
        "ion_mode": _text(row.get("ion_mode")).upper(),
        "collision_energy": _finite(row.get("collision_energy")),
        "collision_energy_unit": _text(row.get("collision_energy_unit")),
        "collision_energy_mode": _text(row.get("collision_energy_mode")),
        "instrument": _text(row.get("instrument")),
        "instrument_type": _text(row.get("instrument_type")),
        "activation_method": _text(
            row.get("activation_method") or row.get("activation")
        ),
        "instrument_domain_id": _text(row.get("instrument_domain_id")),
        "collision_cell_region": _text(row.get("collision_cell_region")),
        "query_adduct": _text(row.get("query_adduct") or match.get("query_adduct")),
        "query_adduct_reliable": bool(
            row.get("query_adduct_reliable") or match.get("query_adduct_reliable")
        ),
        "lcms_adduct": _text(row.get("lcms_adduct")),
        "lcms_neutral_mass": _finite(row.get("lcms_neutral_mass")),
        "lcms_adduct_group_id": _text(row.get("lcms_adduct_group_id")),
        "lcms_adduct_group_evidence": (
            dict(row.get("lcms_adduct_group_evidence") or {})
            if isinstance(row.get("lcms_adduct_group_evidence"), Mapping)
            else {}
        ),
        "ion_identity": _text(row.get("ion_identity")),
        "ion_identity_status": _text(row.get("ion_identity_status")).lower(),
        "ion_identity_network_id": _text(row.get("ion_identity_network_id")),
        "ion_identity_neutral_mass": _finite(
            row.get("ion_identity_neutral_mass")
        ),
        "ion_identity_evidence": (
            dict(row.get("ion_identity_evidence") or {})
            if isinstance(row.get("ion_identity_evidence"), Mapping)
            else {}
        ),
        "parent_formula": _text(row.get("parent_formula")),
        "parent_formula_independent": bool(row.get("parent_formula_independent")),
        "input_format": _text(row.get("input_format")),
        "spectrum_representation": _text(row.get("spectrum_representation")),
        "limited_fragment_list": bool(row.get("limited_fragment_list")),
        "query_evidence_peak_count": len(quality_peaks),
        "query_evidence_peak_count_source": (
            "processed_query_peaks" if query_peaks else "raw_peaks_fallback"
        ),
        "chimeric_deconvolution": {
            key: chimeric.get(key)
            for key in (
                "status",
                "isolation_status",
                "isolation_purity",
                "target_fraction",
                "score_contribution",
                "ranking_changed",
            )
            if key in chimeric
        },
        "raw_peaks": raw_peaks,
        "query_peaks": query_peaks,
    }


def _source_projection(row: Mapping[str, Any], candidates: list[Mapping[str, Any]], entity_id: str) -> dict[str, Any]:
    return {
        "entity_id": entity_id,
        "spectrum_id": row.get("spectrum_id"),
        "precursor_mz": row.get("precursor_mz"),
        "search_precursor_mz": row.get("search_precursor_mz"),
        "search_precursor_mz_source": row.get("search_precursor_mz_source"),
        "retention_time": row.get("retention_time"),
        "collision_energy": row.get("collision_energy"),
        "collision_energy_unit": row.get("collision_energy_unit"),
        "collision_energy_mode": row.get("collision_energy_mode"),
        "instrument": row.get("instrument"),
        "instrument_type": row.get("instrument_type"),
        "activation_method": row.get("activation_method"),
        "instrument_domain_id": row.get("instrument_domain_id"),
        "collision_cell_region": row.get("collision_cell_region"),
        "ion_mode": row.get("ion_mode"),
        "candidates": [
            {
                "identity_key": candidate.get("identity_key"),
                "candidate": candidate.get("candidate"),
                "spectral_score": candidate.get("spectral_score"),
                "scan_support": candidate.get("scan_support"),
            }
            for candidate in candidates
        ],
    }


def _normalize_adduct(value: Any) -> str:
    text = re.sub(r"\s+", "", _text(value))
    aliases = {
        "M+H": "[M+H]+",
        "M-H": "[M-H]-",
        "M+": "[M]+",
        "M-": "[M]-",
        "M+NA": "[M+Na]+",
        "M+K": "[M+K]+",
        "M+NH4": "[M+NH4]+",
        "M+CL": "[M+Cl]-",
    }
    return aliases.get(text.upper(), text)


def _mass_tolerance(mz: float, policy: Mapping[str, Any]) -> float:
    gate = policy.get("chemistry_gate") or {}
    return max(
        float(gate.get("parent_mass_tolerance_da") or 0.02),
        abs(mz) * float(gate.get("parent_mass_tolerance_ppm") or 25.0) * 1e-6,
    )


def _formula_adduct_check(
    candidate: Mapping[str, Any],
    row: Mapping[str, Any],
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    formula = normalize_formula_text(candidate.get("formula"))
    precursor = _search_precursor_mz(row)
    if not formula:
        return {
            "status": "unavailable",
            "reason": "candidate_formula_unavailable",
            "veto": False,
            "missing_is_negative_evidence": False,
        }
    neutral_mass = formula_exact_mass(formula)
    if neutral_mass is None or precursor is None:
        return {
            "status": "unavailable",
            "reason": "formula_mass_or_precursor_unavailable",
            "formula": formula,
            "veto": False,
            "missing_is_negative_evidence": False,
        }

    mode = _text(row.get("ion_mode")).upper()
    library_adduct = _normalize_adduct(candidate.get("adduct"))
    library_expected = ion_mz(neutral_mass, library_adduct) if library_adduct else None
    library_polarity = adduct_polarity(library_adduct) if library_adduct else None
    library_audit: dict[str, Any] = {
        "reported_library_adduct": library_adduct or None,
        "library_adduct_is_candidate_conditioned": bool(library_adduct),
        "library_adduct_can_veto": False,
    }
    if library_adduct:
        if library_polarity and mode.startswith(("POS", "NEG")) and library_polarity != mode[:3]:
            library_audit["library_adduct_audit"] = "polarity_mismatch_not_vetoed"
        elif library_expected is not None:
            library_error = precursor - library_expected
            library_audit.update({
                "library_adduct_audit": (
                    "mass_consistent"
                    if abs(library_error) <= _mass_tolerance(library_expected, policy)
                    else "mass_mismatch_not_vetoed"
                ),
                "library_adduct_expected_precursor_mz": round(library_expected, 8),
                "library_adduct_mass_error_da": round(library_error, 8),
                "library_adduct_mass_error_ppm": round(
                    library_error / library_expected * 1e6, 3
                ),
            })
        else:
            library_audit["library_adduct_audit"] = "unsupported_or_nonspecific"

    # A hard veto requires an observation that exists independently of the
    # library candidate.  Library adduct labels describe the reference spectrum
    # and are not necessarily the adduct measured in the query injection.
    # OpenMS may attach a polarity-default [M+H]+/[M-H]- label to a singleton
    # feature without observing an adduct edge.  V15 accepts LC-MS adduct and
    # neutral-mass evidence for a hard veto only when a consistency-validated
    # adduct group or a non-conflicting runtime ion-identity network supports
    # it.  A raw OpenMS group identifier alone is an inferred hypothesis.
    gate = policy.get("chemistry_gate") or {}
    require_reliable_lcms_evidence = bool(
        gate.get(
            "lcms_adduct_hard_veto_requires_validated_network_or_group_evidence",
            False,
        )
    )
    ion_identity_evidence = (
        row.get("ion_identity_evidence")
        if isinstance(row.get("ion_identity_evidence"), Mapping)
        else {}
    ) or {}
    ion_identity = _normalize_adduct(row.get("ion_identity"))
    ion_identity_reliable = bool(
        ion_identity
        and _text(row.get("ion_identity_status")).lower() == "networked"
        and _text(row.get("ion_identity_network_id"))
        and not bool(ion_identity_evidence.get("identity_conflict"))
        and _text(ion_identity_evidence.get("status") or "networked").lower()
        == "networked"
        and ion_identity_evidence.get("runtime_influence") is True
    )
    lcms_group_evidence = (
        row.get("lcms_adduct_group_evidence")
        if isinstance(row.get("lcms_adduct_group_evidence"), Mapping)
        else {}
    ) or {}
    lcms_group_reliable = bool(
        _text(row.get("lcms_adduct_group_id"))
        and _text(lcms_group_evidence.get("status")).lower()
        in {"supported", "validated"}
        and lcms_group_evidence.get("neutral_mass_consistent") is True
        and not bool(lcms_group_evidence.get("identity_conflict"))
        and int(lcms_group_evidence.get("member_feature_count") or 0) >= 2
        and lcms_group_evidence.get("runtime_influence") is True
    )

    independent_checks: list[dict[str, Any]] = []
    independent_adducts: list[tuple[str, str]] = []
    lcms_adduct = _normalize_adduct(row.get("lcms_adduct"))
    if require_reliable_lcms_evidence and ion_identity_reliable:
        independent_adducts.append(
            ("ion_identity_network_adduct", ion_identity)
        )
    elif lcms_adduct and (
        not require_reliable_lcms_evidence or lcms_group_reliable
    ):
        independent_adducts.append(("lcms_feature_adduct", lcms_adduct))
    query_adduct = _normalize_adduct(row.get("query_adduct"))
    if query_adduct and bool(row.get("query_adduct_reliable")):
        independent_adducts.append(("reliable_query_adduct", query_adduct))
    for source, adduct in dict.fromkeys(independent_adducts):
        expected = ion_mz(neutral_mass, adduct)
        polarity = adduct_polarity(adduct)
        if polarity and mode.startswith(("POS", "NEG")) and polarity != mode[:3]:
            independent_checks.append({
                "source": source,
                "adduct": adduct,
                "consistent": False,
                "reason": "candidate_independent_adduct_polarity_conflict",
            })
            continue
        if expected is None:
            continue
        error = precursor - expected
        independent_checks.append({
            "source": source,
            "adduct": adduct,
            "consistent": abs(error) <= _mass_tolerance(expected, policy),
            "reason": (
                "candidate_independent_formula_adduct_mass_consistent"
                if abs(error) <= _mass_tolerance(expected, policy)
                else "candidate_independent_formula_adduct_mass_mismatch"
            ),
            "expected_precursor_mz": expected,
            "mass_error_da": error,
            "mass_error_ppm": error / expected * 1e6,
        })

    neutral_source = "lcms_feature_neutral_mass"
    if require_reliable_lcms_evidence and ion_identity_reliable:
        observed_neutral = _finite(row.get("ion_identity_neutral_mass"))
        neutral_source = "ion_identity_network_neutral_mass"
    elif not require_reliable_lcms_evidence or lcms_group_reliable:
        observed_neutral = _finite(row.get("lcms_neutral_mass"))
    else:
        observed_neutral = None
    if observed_neutral is not None and observed_neutral > 0:
        neutral_error = neutral_mass - observed_neutral
        independent_checks.append({
            "source": neutral_source,
            "adduct": None,
            "consistent": abs(neutral_error) <= _mass_tolerance(observed_neutral, policy),
            "reason": (
                "candidate_independent_formula_neutral_mass_consistent"
                if abs(neutral_error) <= _mass_tolerance(observed_neutral, policy)
                else "candidate_independent_formula_neutral_mass_mismatch"
            ),
            "expected_neutral_mass": neutral_mass,
            "observed_neutral_mass": observed_neutral,
            "mass_error_da": neutral_error,
            "mass_error_ppm": neutral_error / observed_neutral * 1e6,
        })

    lcms_evidence_audit = {
        "lcms_adduct_evidence_policy": (
            "validated_network_or_group_required"
            if require_reliable_lcms_evidence
            else "legacy_any_lcms_feature_annotation"
        ),
        "ion_identity_network_evidence_available": ion_identity_reliable,
        "lcms_validated_adduct_group_evidence_available": lcms_group_reliable,
        "unvalidated_lcms_adduct_group_ignored": bool(
            require_reliable_lcms_evidence
            and _text(row.get("lcms_adduct_group_id"))
            and not lcms_group_reliable
            and not ion_identity_reliable
        ),
        "unlinked_lcms_feature_adduct_ignored": bool(
            require_reliable_lcms_evidence
            and lcms_adduct
            and not ion_identity_reliable
            and not lcms_group_reliable
        ),
        "unlinked_lcms_feature_neutral_mass_ignored": bool(
            require_reliable_lcms_evidence
            and _finite(row.get("lcms_neutral_mass")) is not None
            and not ion_identity_reliable
            and not lcms_group_reliable
        ),
    }

    evaluable_independent = [check for check in independent_checks if "consistent" in check]
    consistent_independent = [check for check in evaluable_independent if check["consistent"]]
    if evaluable_independent and not consistent_independent:
        selected = min(
            evaluable_independent,
            key=lambda check: abs(float(check.get("mass_error_ppm") or 1e12)),
        )
        configured_transformations = {
            _normalize_adduct(value)
            for value in (
                gate.get("precursor_transformation_ambiguity_adducts") or []
            )
            if _normalize_adduct(value)
        }
        eligible_independent_sources = {
            _text(value)
            for value in (
                gate.get(
                    "precursor_transformation_ambiguity_independent_sources"
                ) or []
            )
            if _text(value)
        }
        independent_sources = {
            _text(check.get("source")) for check in evaluable_independent
        }
        transformation_ambiguity = bool(
            gate.get(
                "mass_consistent_library_neutral_loss_can_mark_precursor_transformation_ambiguity",
                False,
            )
            and library_adduct in configured_transformations
            and library_audit.get("library_adduct_audit") == "mass_consistent"
            and independent_sources
            and independent_sources.issubset(eligible_independent_sources)
        )
        if transformation_ambiguity:
            return {
                "status": "ambiguous",
                "reason": (
                    "candidate_independent_ion_identity_conflicts_with_"
                    "mass_consistent_library_neutral_loss_adduct"
                ),
                "formula": formula,
                "adduct": library_adduct,
                "adduct_source": "candidate_library_neutral_loss_annotation",
                "expected_precursor_mz": round(float(library_expected), 8),
                "observed_precursor_mz": round(precursor, 8),
                "mass_error_da": library_audit.get(
                    "library_adduct_mass_error_da"
                ),
                "mass_error_ppm": library_audit.get(
                    "library_adduct_mass_error_ppm"
                ),
                "candidate_independent_checks": independent_checks,
                "candidate_independent_contradiction_preserved": True,
                "candidate_conditioned": True,
                "precursor_transformation_ambiguity": True,
                "precursor_transformation_identity_confirmed": False,
                "veto": False,
                **lcms_evidence_audit,
                **library_audit,
            }
        return {
            "status": "contradicted",
            "reason": selected["reason"],
            "formula": formula,
            "adduct": selected.get("adduct"),
            "adduct_source": selected["source"],
            "expected_precursor_mz": (
                round(float(selected["expected_precursor_mz"]), 8)
                if selected.get("expected_precursor_mz") is not None else None
            ),
            "observed_precursor_mz": round(precursor, 8),
            "mass_error_da": round(float(selected.get("mass_error_da") or 0.0), 8),
            "mass_error_ppm": round(float(selected.get("mass_error_ppm") or 0.0), 3),
            "candidate_independent_checks": independent_checks,
            "candidate_conditioned": False,
            "veto": True,
            **lcms_evidence_audit,
            **library_audit,
        }

    independent_adduct_match = next(
        (check for check in consistent_independent if check.get("adduct")),
        None,
    )
    if independent_adduct_match is not None:
        return {
            "status": "consistent",
            "reason": independent_adduct_match["reason"],
            "formula": formula,
            "adduct": independent_adduct_match["adduct"],
            "adduct_source": independent_adduct_match["source"],
            "expected_precursor_mz": round(
                float(independent_adduct_match["expected_precursor_mz"]), 8
            ),
            "observed_precursor_mz": round(precursor, 8),
            "mass_error_da": round(float(independent_adduct_match["mass_error_da"]), 8),
            "mass_error_ppm": round(float(independent_adduct_match["mass_error_ppm"]), 3),
            "candidate_independent_checks": independent_checks,
            "candidate_conditioned": False,
            "veto": False,
            "ion_identity_network_evidence_used": bool(
                independent_adduct_match["source"]
                == "ion_identity_network_adduct"
            ),
            **lcms_evidence_audit,
            **library_audit,
        }

    hypotheses = (
        gate.get("candidate_mass_inference_adducts_positive") or []
        if mode.startswith("POS")
        else gate.get("candidate_mass_inference_adducts_negative") or []
        if mode.startswith("NEG")
        else []
    )
    matches = []
    for adduct in hypotheses:
        expected = ion_mz(neutral_mass, adduct)
        if expected is None:
            continue
        error = precursor - expected
        if abs(error) <= _mass_tolerance(expected, policy):
            matches.append((abs(error), str(adduct), expected, error))
    matches.sort()
    if len(matches) == 1:
        _, adduct, expected, error = matches[0]
        return {
            "status": "consistent",
            "reason": "unique_candidate_formula_mass_adduct_hypothesis",
            "formula": formula,
            "adduct": adduct,
            "adduct_source": "candidate_formula_mass_unique_inference",
            "expected_precursor_mz": round(expected, 8),
            "observed_precursor_mz": round(precursor, 8),
            "mass_error_da": round(error, 8),
            "mass_error_ppm": round(error / expected * 1e6, 3),
            "veto": False,
            "candidate_conditioned": True,
            "candidate_independent_checks": independent_checks,
            **lcms_evidence_audit,
            **library_audit,
        }
    if len(matches) > 1:
        return {
            "status": "ambiguous",
            "reason": "multiple_candidate_formula_mass_adduct_hypotheses",
            "formula": formula,
            "adduct_hypotheses": [item[1] for item in matches],
            "veto": False,
            "candidate_conditioned": True,
            "candidate_independent_checks": independent_checks,
            **lcms_evidence_audit,
            **library_audit,
        }
    return {
        "status": "unavailable",
        "reason": "no_unique_supported_adduct_hypothesis",
        "formula": formula,
        "reported_adduct": library_adduct or None,
        "veto": False,
        "missing_is_negative_evidence": False,
        "candidate_conditioned": True,
        "candidate_independent_checks": independent_checks,
        **lcms_evidence_audit,
        **library_audit,
    }


def _aggregate_chemistry(checks: list[dict[str, Any]]) -> dict[str, Any]:
    statuses = Counter(str(check.get("status") or "unavailable") for check in checks)
    consistent = [check for check in checks if check.get("status") == "consistent"]
    contradicted = [check for check in checks if check.get("status") == "contradicted"]
    transformation_ambiguities = [
        check for check in checks if check.get("precursor_transformation_ambiguity")
    ]
    evaluable = consistent + contradicted
    if transformation_ambiguities:
        selected = min(
            transformation_ambiguities,
            key=lambda item: abs(float(item.get("mass_error_ppm") or 1e12)),
        )
        status, veto = "ambiguous", False
    elif consistent:
        selected = min(
            consistent,
            key=lambda item: (
                abs(float(item.get("mass_error_ppm") or 0.0)),
                _text(item.get("adduct")),
            ),
        )
        status, veto = "consistent", False
    elif evaluable and len(contradicted) == len(evaluable):
        selected = min(
            contradicted,
            key=lambda item: abs(float(item.get("mass_error_ppm") or 1e12)),
        )
        status, veto = "contradicted", True
    else:
        selected = next((check for check in checks if check.get("status") == "ambiguous"), checks[0] if checks else {})
        status, veto = str(selected.get("status") or "unavailable"), False
    return {
        **selected,
        "status": status,
        "veto": veto,
        "status_counts": dict(sorted(statuses.items())),
        "evaluated_collision_energy_bins": len(checks),
        "precursor_transformation_ambiguity": bool(
            transformation_ambiguities
        ),
        "precursor_transformation_ambiguity_bin_count": len(
            transformation_ambiguities
        ),
        "consistent_metadata_overrides_isolated_conflicting_library_metadata": bool(consistent and contradicted),
    }


def _fragmentation_explanation(
    candidate: Mapping[str, Any],
    chemistry: Mapping[str, Any],
    representative_row: Mapping[str, Any] | None,
) -> dict[str, Any]:
    if representative_row is None:
        return {
            "status": "unavailable",
            "reason": "representative_spectrum_unavailable",
            "formula_conserving_only": True,
            "mass_only_mechanisms_emitted": 0,
        }
    # Only an adduct resolved by the chemistry check may seed the graph.  A
    # reference-library adduct may differ from the query injection and is never
    # silently substituted here.
    adduct = _text(chemistry.get("adduct"))
    formula = normalize_formula_text(candidate.get("formula"))
    if not formula or not adduct:
        return {
            "status": "unavailable",
            "reason": "formula_or_unique_adduct_unavailable",
            "formula": formula or None,
            "adduct": adduct or None,
            "formula_conserving_only": True,
            "mass_only_mechanisms_emitted": 0,
            "mechanism_claimed": False,
        }
    try:
        from matching.evograph_fragmentation import build_evograph_graph

        graph_row = dict(representative_row)
        graph_row["match"] = {**dict(candidate), "formula": formula, "adduct": adduct}
        graph = build_evograph_graph(
            graph_row,
            result_index=int(representative_row.get("result_index") or 0),
        )
    except Exception as exc:
        return {
            "status": "failed",
            "reason": f"{type(exc).__name__}: {exc}"[:300],
            "formula": formula,
            "adduct": adduct,
            "formula_conserving_only": True,
            "mass_only_mechanisms_emitted": 0,
            "mechanism_claimed": False,
        }

    nodes = []
    for node in graph.get("nodes") or []:
        formula_candidates = [
            {
                key: item.get(key)
                for key in ("formula", "mass_error_da", "mass_error_ppm", "dbe")
                if item.get(key) is not None
            }
            for item in (node.get("formula_candidates") or [])[:3]
        ]
        nodes.append({
            "node_id": node.get("node_id"),
            "node_type": node.get("node_type"),
            "mz": node.get("mz"),
            "relative_intensity": node.get("relative_intensity"),
            "formula_candidates": formula_candidates,
            "assignment_status": (
                "formula_candidates" if formula_candidates else "unassigned"
            ),
        })
    edges = [
        {
            key: edge.get(key)
            for key in (
                "parent_id",
                "child_id",
                "kind",
                "observed_loss_da",
                "loss_formula",
                "theoretical_loss_da",
                "mass_error_da",
                "mass_error_ppm",
                "parent_formula",
                "child_formula",
                "formula_conservation",
            )
        }
        for edge in (graph.get("edges") or [])
        if edge.get("kind") == "formula_conserving"
        and edge.get("formula_conservation") is True
    ]
    assigned_fragments = sum(
        bool(node.get("formula_candidates"))
        for node in nodes
        if node.get("node_type") == "fragment"
    )
    fragment_nodes = sum(node.get("node_type") == "fragment" for node in nodes)
    status = "completed" if edges else "partial" if assigned_fragments else "unavailable"
    return {
        "status": status,
        "reason": None if status == "completed" else (
            "fragment_formula_assignments_without_conserving_edges"
            if status == "partial"
            else graph.get("reason") or "no_formula_conserving_fragment_assignment"
        ),
        "formula": formula,
        "adduct": adduct,
        "candidate_conditioned": True,
        "formula_conserving_only": True,
        "mass_only_mechanisms_emitted": 0,
        "mechanism_claimed": False,
        "fragment_nodes": fragment_nodes,
        "assigned_fragment_nodes": assigned_fragments,
        "formula_conserving_edge_count": len(edges),
        "nodes": nodes,
        "edges": edges,
    }


def _candidate_variant(candidate: Mapping[str, Any]) -> tuple[str, ...]:
    return (
        _text(candidate.get("inchikey")),
        _text(candidate.get("smiles")),
        normalize_formula_text(candidate.get("formula")),
        _text(candidate.get("adduct")),
        _text(candidate.get("compound_name")),
    )


def _query_spectrum_release_audit(
    scans: list[dict[str, Any]],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit candidate-independent query quality before a 2-D release.

    The gate never changes retrieval, candidate order, or any score.  A
    technical replicate can rescue another replicate within the same CE bin,
    but it cannot create an additional measurement condition.
    """

    minimum_peaks = max(
        1, int(config.get("minimum_valid_fragment_peaks_per_ce_bin") or 3)
    )
    minimum_bins = max(
        1, int(config.get("minimum_quality_qualified_ce_bins") or 2)
    )
    guarded_isolation = {
        _text(value).lower()
        for value in (
            config.get("chimeric_isolation_statuses_requiring_deconvolution")
            or ("mixed", "severe")
        )
    }
    accepted_deconvolution = {
        _text(value).lower()
        for value in (
            config.get("accepted_chimeric_deconvolution_statuses")
            or ("accepted_score_applied",)
        )
    }
    fail_unreliable_precursor = bool(
        config.get("explicit_unreliable_search_precursor_is_failure", True)
    )
    fail_limited_list = bool(
        config.get("limited_fragment_list_is_failure", True)
    )

    bins: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for scan in scans:
        bins[str(scan.get("ce_key") or "unknown")].append(scan)

    bin_audits = []
    quality_bins = 0
    failure_counts: Counter[str] = Counter()
    for ce_key in sorted(bins):
        scan_audits = []
        for scan in bins[ce_key]:
            row = scan.get("row") or {}
            failures = []
            peak_count = int(row.get("query_evidence_peak_count") or 0)
            if peak_count < minimum_peaks:
                failures.append("too_few_valid_fragment_peaks")
            if (
                fail_unreliable_precursor
                and row.get("search_precursor_mz_reliability_status")
                == "unreliable"
            ):
                failures.append("search_precursor_mz_explicitly_unreliable")
            if fail_limited_list and row.get("limited_fragment_list") is True:
                failures.append("limited_fragment_list_input")
            deconvolution = row.get("chimeric_deconvolution") or {}
            isolation_status = _text(
                deconvolution.get("isolation_status")
            ).lower()
            deconvolution_status = _text(
                deconvolution.get("status")
            ).lower()
            if (
                isolation_status in guarded_isolation
                and deconvolution_status not in accepted_deconvolution
            ):
                failures.append(
                    "mixed_or_severe_isolation_without_accepted_deconvolution"
                )
            failure_counts.update(failures)
            scan_audits.append(
                {
                    "result_index": int(row.get("result_index") or 0),
                    "spectrum_id": row.get("spectrum_id"),
                    "query_evidence_peak_count": peak_count,
                    "query_evidence_peak_count_source": row.get(
                        "query_evidence_peak_count_source"
                    ),
                    "search_precursor_mz_reliability_status": row.get(
                        "search_precursor_mz_reliability_status"
                    ),
                    "search_precursor_mz_reliability_source": row.get(
                        "search_precursor_mz_reliability_source"
                    ),
                    "limited_fragment_list": bool(
                        row.get("limited_fragment_list")
                    ),
                    "chimeric_isolation_status": isolation_status or None,
                    "chimeric_deconvolution_status": (
                        deconvolution_status or None
                    ),
                    "qualified": not failures,
                    "failed_checks": failures,
                }
            )
        qualified = any(item["qualified"] for item in scan_audits)
        quality_bins += int(qualified)
        bin_audits.append(
            {
                "collision_energy_bin": ce_key,
                "technical_replicate_count": len(scan_audits),
                "qualified": qualified,
                "technical_replicates_create_extra_conditions": False,
                "scans": scan_audits,
            }
        )

    failed_gates = []
    if quality_bins < minimum_bins:
        failed_gates.append("insufficient_quality_qualified_ce_bins")
    affects_release_scope = bool(config.get("affects_release_scope"))
    return {
        "enabled": bool(config.get("enabled")),
        "algorithm_id": config.get("algorithm_id"),
        "deployment_state": config.get("deployment_state"),
        "affects_release_scope": affects_release_scope,
        "status": (
            "demoted"
            if failed_gates and affects_release_scope
            else "failed_audit_only"
            if failed_gates
            else "passed"
        ),
        "failed_gates": failed_gates,
        "failure_reason_counts": dict(sorted(failure_counts.items())),
        "observed_quality_qualified_ce_bins": quality_bins,
        "observed_total_ce_bins": len(bin_audits),
        "minimum_valid_fragment_peaks_per_ce_bin": minimum_peaks,
        "minimum_quality_qualified_ce_bins": minimum_bins,
        "explicit_unreliable_search_precursor_is_failure": (
            fail_unreliable_precursor
        ),
        "limited_fragment_list_is_failure": fail_limited_list,
        "chimeric_isolation_statuses_requiring_deconvolution": sorted(
            guarded_isolation
        ),
        "accepted_chimeric_deconvolution_statuses": sorted(
            accepted_deconvolution
        ),
        "collision_energy_evidence": bin_audits,
        "winner_rank_changed": False,
        "winner_score_changed": False,
        "score_is_probability": False,
    }


def _replicate_support_median(values: list[float | None], *, omit_unobserved: bool) -> float:
    """Keep measured zeroes, but do not invent them for truncated Top-K absence."""
    observed = [value for value in values if value is not None]
    basis = observed if omit_unobserved else [0.0 if value is None else value for value in values]
    return statistics.median(basis) if basis else 0.0


def _finalize_entity(entity: Mapping[str, Any], policy: Mapping[str, Any]) -> dict[str, Any]:
    pool_policy = policy.get("candidate_pool") or {}
    replicate_rule = pool_policy.get("same_energy_replicate_support", "legacy_zero_fill")
    if replicate_rule not in {"legacy_zero_fill", "observed_candidates_only_v16"}:
        raise ValueError("Unknown same-energy replicate aggregation policy")
    topk_absence_is_negative = bool(pool_policy.get("missing_from_truncated_topk_is_negative_evidence", False))
    omit_unobserved_replicates = replicate_rule == "observed_candidates_only_v16" and not topk_absence_is_negative
    scans = sorted(
        list(entity.get("scans") or []),
        key=lambda scan: (
            -float(scan["row"].get("tic") or 0.0),
            str(scan["row"].get("spectrum_id") or ""),
            int(scan["row"]["result_index"]),
        ),
    )
    bins: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for scan in scans:
        bins[str(scan["ce_key"])].append(scan)
    bin_keys = sorted(bins)
    candidate_ce: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "bins": {},
            "instances": [],
            "variants": set(),
            "reference_source_families": set(),
        }
    )
    ce_winners: dict[str, str] = {}
    legacy_ce_winners: dict[str, str] = {}

    for ce_key in bin_keys:
        replicate_scans = bins[ce_key]
        identities = sorted({
            item["identity_key"]
            for scan in replicate_scans
            for item in scan.get("candidates") or []
        })
        ce_records = []
        legacy_ce_records = []
        for identity in identities:
            supports = []
            legacy_supports = []
            instances = []
            top1_count = 0
            legacy_top1_count = 0
            for scan in replicate_scans:
                item = next(
                    (candidate for candidate in scan.get("candidates") or [] if candidate["identity_key"] == identity),
                    None,
                )
                supports.append(float(item.get("scan_support") or 0.0) if item else None)
                legacy_supports.append(
                    float(item.get("legacy_scan_support") or 0.0)
                    if item else 0.0
                )
                if item:
                    instances.append((float(item.get("scan_support") or 0.0), scan, item))
                    top1_count += int(int(item.get("scan_rank") or 999999) == 1)
                    legacy_top1_count += int(
                        int(item.get("legacy_scan_rank") or 999999) == 1
                    )
                    candidate_ce[identity]["variants"].add(_candidate_variant(item["candidate"]))
                    candidate_ce[identity]["reference_source_families"].update(
                        item.get("reference_source_families") or []
                    )
            ce_support = _replicate_support_median(supports, omit_unobserved=omit_unobserved_replicates)
            legacy_ce_support = (
                statistics.median(legacy_supports) if legacy_supports else 0.0
            )
            representative = max(
                instances,
                key=lambda value: (value[0], -int(value[1]["row"]["result_index"])),
                default=None,
            )
            record = {
                "collision_energy_bin": ce_key,
                "support": round(ce_support, 6),
                "legacy_v14_support": round(legacy_ce_support, 6),
                "technical_replicates": len(replicate_scans),
                "replicates_with_candidate": len(instances),
                "replicates_missing_from_truncated_topk": len(replicate_scans) - len(instances),
                "truncated_topk_absence_imputed_as_zero": not omit_unobserved_replicates,
                "technical_replicates_count_as_independent_votes": False,
                "scan_top1_fraction": round(top1_count / len(replicate_scans), 6) if replicate_scans else 0.0,
                "legacy_v14_scan_top1_fraction": (
                    round(legacy_top1_count / len(replicate_scans), 6)
                    if replicate_scans else 0.0
                ),
                "representative_result_index": (
                    int(representative[1]["row"]["result_index"]) if representative else None
                ),
                "conditioned_multilibrary_evidence": (
                    representative[2].get("conditioned_multilibrary_evidence")
                    if representative else None
                ),
            }
            candidate_ce[identity]["bins"][ce_key] = record
            if representative:
                candidate_ce[identity]["instances"].append(representative)
            ce_records.append((ce_support, identity))
            legacy_ce_records.append((legacy_ce_support, identity))
        if ce_records:
            ce_winners[ce_key] = min(ce_records, key=lambda value: (-value[0], value[1]))[1]
        if legacy_ce_records:
            legacy_ce_winners[ce_key] = min(
                legacy_ce_records, key=lambda value: (-value[0], value[1])
            )[1]

    weights = policy.get("entity_score") or {}
    candidates = []
    total_bins = len(bin_keys)
    for identity, data in sorted(candidate_ce.items()):
        supports = [float((data["bins"].get(key) or {}).get("support") or 0.0) for key in bin_keys]
        legacy_supports = [
            float(
                (data["bins"].get(key) or {}).get("legacy_v14_support")
                or 0.0
            )
            for key in bin_keys
        ]
        # The versioned V16 rule also retains a genuinely measured zero. Only
        # a CE bin with no observed candidate is absent; score > 0 is not an
        # observation test. Positive support still controls coverage below.
        observed_supports = (
            [float(data["bins"][key]["support"]) for key in bin_keys if key in data["bins"]]
            if omit_unobserved_replicates else [support for support in supports if support > 0]
        )
        legacy_observed_supports = [
            support for support in legacy_supports if support > 0
        ]
        observed_bins = sum(support > 0 for support in observed_supports)
        # The scan artifact exposes a truncated Top-K, not an exhaustive score
        # vector.  Absence outside that list is therefore missing evidence, not
        # a measured zero.  Coverage remains an explicit score component and a
        # release-scope constraint, so sparse support cannot masquerade as a
        # replicated multi-energy identification.
        support_basis = supports if topk_absence_is_negative else observed_supports
        legacy_support_basis = (
            legacy_supports
            if topk_absence_is_negative
            else legacy_observed_supports
        )
        mean_support = statistics.fmean(support_basis) if support_basis else 0.0
        median_support = statistics.median(support_basis) if support_basis else 0.0
        legacy_mean_support = (
            statistics.fmean(legacy_support_basis)
            if legacy_support_basis else 0.0
        )
        legacy_median_support = (
            statistics.median(legacy_support_basis)
            if legacy_support_basis else 0.0
        )
        coverage = observed_bins / total_bins if total_bins else 0.0
        # The embedded historical lane retains the old padded-replicate result,
        # including its coverage. A V16 correction must not rewrite that audit.
        legacy_coverage = (
            len(legacy_observed_supports) / total_bins if total_bins else 0.0
        ) if omit_unobserved_replicates else coverage
        top1_fraction = (
            sum(ce_winners.get(key) == identity for key in bin_keys) / total_bins
            if total_bins else 0.0
        )
        legacy_top1_fraction = (
            sum(legacy_ce_winners.get(key) == identity for key in bin_keys)
            / total_bins
            if total_bins else 0.0
        )
        source_families = sorted(data["reference_source_families"])
        source_saturation = max(
            1,
            int((policy.get("reference_provenance") or {}).get(
                "source_diversity_saturation", 3
            )),
        )
        source_diversity = min(1.0, len(source_families) / source_saturation)
        entity_score = (
            float(weights.get("mean_collision_energy_support") or 0.45) * mean_support
            + float(weights.get("median_collision_energy_support") or 0.30) * median_support
            + float(weights.get("collision_energy_coverage") or 0.15) * coverage
            + float(weights.get("collision_energy_top1_fraction") or 0.10) * top1_fraction
        )
        legacy_entity_score = (
            float(weights.get("mean_collision_energy_support") or 0.45)
            * legacy_mean_support
            + float(weights.get("median_collision_energy_support") or 0.30)
            * legacy_median_support
            + float(weights.get("collision_energy_coverage") or 0.15)
            * legacy_coverage
            + float(weights.get("collision_energy_top1_fraction") or 0.10)
            * legacy_top1_fraction
        )
        instances = sorted(
            data["instances"],
            key=lambda value: (-value[0], int(value[1]["row"]["result_index"])),
        )
        checks = [
            _formula_adduct_check(item[2]["candidate"], item[1]["row"], policy)
            for item in instances
        ]
        chemistry = _aggregate_chemistry(checks)
        display = dict(instances[0][2]["candidate"]) if instances else {}
        candidate = {
            "identity_key": identity,
            "identity_scope": "2d_connectivity",
            "identity_source": instances[0][2].get("identity_source") if instances else "unknown",
            "candidate": display,
            "entity_score": round(max(0.0, min(1.0, entity_score)), 6),
            "legacy_v14_entity_score": round(
                max(0.0, min(1.0, legacy_entity_score)), 6
            ),
            "score_is_probability": False,
            "score_components": {
                "mean_collision_energy_support": round(mean_support, 6),
                "median_collision_energy_support": round(median_support, 6),
                "collision_energy_coverage": round(coverage, 6),
                "collision_energy_top1_fraction": round(top1_fraction, 6),
                "topk_absence_imputed_as_zero": topk_absence_is_negative,
                "reference_source_diversity": round(source_diversity, 6),
                "reference_source_diversity_affects_score": False,
                "source_count_bonus_applied": False,
            },
            "reference_source_family_count": len(source_families),
            "reference_source_families": source_families,
            "reference_source_accessions_count_as_independent": False,
            "collision_energy_bins_supported": observed_bins,
            "collision_energy_bins_total": total_bins,
            "collision_energy_evidence": [data["bins"][key] for key in bin_keys if key in data["bins"]],
            "representative_result_index": (
                int(instances[0][1]["row"]["result_index"]) if instances else None
            ),
            "representative_scan_score": round(float(instances[0][0]), 6) if instances else None,
            "chemistry_gate": chemistry,
            "vetoed": bool(chemistry.get("veto")),
            "variant_count_collapsed_to_connectivity": len(data["variants"]),
            "stereochemistry_claimed": False,
            "ranking_inputs": [
                "cosine",
                "entropy",
                "reverse_match",
                "mass_accuracy",
                "fragment_coverage",
                "intensity_agreement",
                "within_scan_rank",
                "collision_energy_consistency",
            ] + (
                [
                    "mirror_collapsed_reference_profiles",
                    "provenance_stratum_median_similarity",
                    "compatible_instrument_domain_near_ce_preference",
                ]
                if (policy.get("conditioned_multilibrary") or {}).get("enabled")
                else []
            ),
            "conditioned_multilibrary_runtime_influence": any(
                bool(
                    (
                        (data["bins"].get(key) or {}).get(
                            "conditioned_multilibrary_evidence"
                        )
                        or {}
                    ).get("runtime_influence")
                )
                for key in bin_keys
            ),
        }
        candidate["_representative_row"] = instances[0][1]["row"] if instances else None
        candidates.append(candidate)

    legacy_order = sorted(
        candidates,
        key=lambda candidate: (
            bool(candidate["vetoed"]),
            -float(candidate["legacy_v14_entity_score"]),
            candidate["identity_key"],
        ),
    )
    legacy_rank_by_identity = {
        candidate["identity_key"]: rank
        for rank, candidate in enumerate(legacy_order, start=1)
    }
    for candidate in candidates:
        candidate["legacy_v14_entity_rank"] = legacy_rank_by_identity[
            candidate["identity_key"]
        ]
    candidates.sort(
        key=lambda candidate: (
            bool(candidate["vetoed"]),
            -float(candidate["entity_score"]),
            candidate["identity_key"],
        )
    )
    valid = [candidate for candidate in candidates if not candidate["vetoed"]]
    decision_policy = policy.get("decision") or {}
    minimum = float(decision_policy.get("minimum_entity_score") or 0.50)
    margin_required = float(decision_policy.get("minimum_connectivity_margin") or 0.04)
    minimum_ce = int(decision_policy.get("minimum_distinct_collision_energies_for_connectivity") or 2)
    minimum_supported_ce = int(decision_policy.get("minimum_winner_collision_energy_support") or 2)
    isomer_margin = float(decision_policy.get("isomer_group_margin") or 0.05)

    leading = valid[0] if valid else None
    runner = valid[1] if len(valid) > 1 else None
    margin = (
        float(leading["entity_score"]) - float(runner["entity_score"])
        if leading and runner else float(leading["entity_score"]) if leading else 0.0
    )
    winner = leading if leading and float(leading["entity_score"]) >= minimum else None
    if winner is None:
        decision_status = "unmatched" if not candidates else "abstained"
        decision_scope = "none"
        blockers = ["no_library_candidate"] if not candidates else ["entity_score_below_minimum"]
        if candidates and not valid:
            blockers = ["all_candidates_failed_candidate_independent_chemistry_gate"]
    else:
        transformation_ambiguity_candidates = [
            candidate
            for candidate in valid
            if candidate.get("chemistry_gate", {}).get(
                "precursor_transformation_ambiguity"
            )
            and float(candidate.get("entity_score") or 0.0) >= minimum
        ]
        same_formula_competitors = [
            candidate
            for candidate in valid[1:]
            if normalize_formula_text(candidate["candidate"].get("formula"))
            and normalize_formula_text(candidate["candidate"].get("formula"))
            == normalize_formula_text(winner["candidate"].get("formula"))
            and float(winner["entity_score"]) - float(candidate["entity_score"]) <= isomer_margin
        ]
        if transformation_ambiguity_candidates:
            decision_status = "selected_candidate_only"
            decision_scope = "precursor_transformation_ambiguity_candidate"
            blockers = [
                "candidate_independent_ion_identity_conflicts_with_"
                "mass_consistent_library_neutral_loss_adduct"
            ]
        elif same_formula_competitors:
            decision_status = "selected_isomer_group"
            decision_scope = "molecular_formula_isomer_group"
            blockers = ["competitive_same_formula_connectivity"]
        elif (
            total_bins >= minimum_ce
            and int(winner["collision_energy_bins_supported"]) >= minimum_supported_ce
            and margin >= margin_required
        ):
            decision_status = "selected"
            decision_scope = "2d_connectivity_candidate"
            blockers = []
        elif total_bins < minimum_ce:
            decision_status = "selected_candidate_only"
            decision_scope = "single_collision_energy_candidate"
            blockers = ["insufficient_distinct_collision_energies"]
        else:
            decision_status = "selected_candidate_only"
            decision_scope = "ambiguous_connectivity_candidate"
            blockers = ["insufficient_entity_margin_or_energy_coverage"]

    base_formal_release_eligible = bool(
        winner is not None
        and decision_status == "selected"
        and decision_scope == "2d_connectivity_candidate"
    )
    query_quality_config = policy.get("query_spectrum_release_gate") or {}
    query_quality_audit: dict[str, Any] = {
        "enabled": bool(query_quality_config.get("enabled")),
        "algorithm_id": query_quality_config.get("algorithm_id"),
        "status": "not_applicable",
        "failed_gates": [],
        "winner_rank_changed": False,
        "winner_score_changed": False,
        "score_is_probability": False,
    }
    query_quality_failures: list[str] = []
    if (
        base_formal_release_eligible
        and query_quality_config.get("enabled") is True
    ):
        query_quality_audit = _query_spectrum_release_audit(
            scans, query_quality_config
        )
        if query_quality_config.get("affects_release_scope") is True:
            query_quality_failures = list(
                query_quality_audit.get("failed_gates") or []
            )

    abstention_config = policy.get("cross_platform_abstention_gate") or {}
    abstention_audit: dict[str, Any] = {
        "enabled": bool(abstention_config.get("enabled")),
        "algorithm_id": abstention_config.get("algorithm_id"),
        "status": "not_applicable",
        "failed_gates": [],
        "winner_rank_changed": False,
        "winner_score_changed": False,
        "score_is_probability": False,
        "external_independent_validation_complete": bool(
            abstention_config.get("external_independent_validation_complete")
        ),
    }
    if (
        base_formal_release_eligible
        and abstention_config.get("enabled") is True
    ):
        minimum_release_score = float(
            abstention_config.get("minimum_entity_score") or 0.0
        )
        minimum_release_margin = float(
            abstention_config.get("minimum_winner_margin") or 0.0
        )
        minimum_release_ce = int(
            abstention_config.get(
                "minimum_collision_energy_bins_supported"
            ) or 0
        )
        release_failures: list[str] = []
        if float(winner.get("entity_score") or 0.0) < minimum_release_score:
            release_failures.append("entity_score_below_cross_platform_release_gate")
        if margin < minimum_release_margin:
            release_failures.append("winner_margin_below_cross_platform_release_gate")
        if int(winner.get("collision_energy_bins_supported") or 0) < minimum_release_ce:
            release_failures.append("winner_ce_support_below_cross_platform_release_gate")
        abstention_audit.update(
            {
                "status": "demoted" if release_failures else "passed",
                "failed_gates": release_failures,
                "observed_entity_score": round(
                    float(winner.get("entity_score") or 0.0), 6
                ),
                "observed_winner_margin": round(margin, 6),
                "observed_collision_energy_bins_supported": int(
                    winner.get("collision_energy_bins_supported") or 0
                ),
                "minimum_entity_score": minimum_release_score,
                "minimum_winner_margin": minimum_release_margin,
                "minimum_collision_energy_bins_supported": minimum_release_ce,
                "calibration_state": abstention_config.get("calibration_state"),
            }
        )

    abstention_failures = list(abstention_audit.get("failed_gates") or [])
    additional_release_failures = list(
        dict.fromkeys([*query_quality_failures, *abstention_failures])
    )
    if additional_release_failures:
        decision_status = "selected_candidate_only"
        if query_quality_failures and abstention_failures:
            decision_scope = "multi_gate_abstention_candidate"
        elif query_quality_failures:
            decision_scope = "query_spectrum_quality_candidate"
        else:
            decision_scope = "cross_platform_abstention_candidate"
        blockers = additional_release_failures

    if winner is not None:
        winner["fragmentation_explanation"] = _fragmentation_explanation(
            winner["candidate"],
            winner["chemistry_gate"],
            winner.get("_representative_row"),
        )

    shortlist_size = max(1, int((policy.get("candidate_pool") or {}).get("shortlist_size") or 5))
    public_candidates = []
    for rank, candidate in enumerate(candidates[:shortlist_size], start=1):
        clean = {key: value for key, value in candidate.items() if not key.startswith("_")}
        clean["entity_rank"] = rank
        public_candidates.append(clean)
    public_winner = None
    if winner is not None:
        public_winner = next(
            (candidate for candidate in public_candidates if candidate["identity_key"] == winner["identity_key"]),
            {key: value for key, value in winner.items() if not key.startswith("_")},
        )
    public_leading = None
    if leading is not None:
        public_leading = next(
            (candidate for candidate in public_candidates if candidate["identity_key"] == leading["identity_key"]),
            {key: value for key, value in leading.items() if not key.startswith("_")},
        )
    legacy_leading = next(
        (candidate for candidate in legacy_order if not candidate["vetoed"]),
        None,
    )

    # The scan-level comparator is an audit baseline, not an input-order
    # accident.  Use the highest-TIC member and deterministic tie-breakers.
    baseline_scan = scans[0] if scans else None
    baseline_identity = (
        baseline_scan["candidates"][0]["identity_key"]
        if baseline_scan and baseline_scan.get("candidates")
        else None
    )

    precursor_values = [
        float(scan["row"]["search_precursor_mz"])
        for scan in scans
        if scan["row"].get("search_precursor_mz") is not None
    ]
    selected_ion_values = [
        float(scan["row"]["precursor_mz"])
        for scan in scans
        if scan["row"].get("precursor_mz") is not None
    ]
    precursor_sources = sorted({
        _text(scan["row"].get("search_precursor_mz_source"))
        for scan in scans
        if _text(scan["row"].get("search_precursor_mz_source"))
    })
    rt_values = [
        float(scan["row"]["retention_time"])
        for scan in scans
        if scan["row"].get("retention_time") is not None
    ]
    tic_values = [
        float(scan["row"]["tic"])
        for scan in scans
        if scan["row"].get("tic") is not None
    ]
    # These two arrays form a positional binding consumed by EviRecall and
    # EviTwin.  Sorting them independently corrupts that binding whenever the
    # lexical spectrum-id order differs from the immutable result-row order.
    member_pairs = sorted(
        (
            int(scan["row"]["result_index"]),
            str(scan["row"]["spectrum_id"]),
        )
        for scan in scans
    )
    return {
        "entity_id": entity["entity_id"],
        "entity_source": entity["entity_source"],
        "entity_source_value": entity["entity_source_value"],
        "member_scan_count": len(scans),
        "member_result_indices": [index for index, _spectrum_id in member_pairs],
        "member_spectrum_ids": [spectrum_id for _index, spectrum_id in member_pairs],
        "representative_precursor_mz": round(statistics.median(precursor_values), 8) if precursor_values else None,
        "representative_search_precursor_mz": round(statistics.median(precursor_values), 8) if precursor_values else None,
        "representative_selected_ion_mz": round(statistics.median(selected_ion_values), 8) if selected_ion_values else None,
        "search_precursor_mz_sources": precursor_sources,
        "representative_search_precursor_mz_source": (
            precursor_sources[0] if len(precursor_sources) == 1 else
            "mixed" if precursor_sources else None
        ),
        "lcms_feature_used_for_search": any(
            bool(scan["row"].get("lcms_feature_used_for_search"))
            for scan in scans
        ),
        "representative_retention_time": round(statistics.median(rt_values), 8) if rt_values else None,
        "representative_tic": max(tic_values) if tic_values else None,
        "ion_modes": sorted({_text(scan["row"].get("ion_mode")) for scan in scans if _text(scan["row"].get("ion_mode"))}),
        "distinct_collision_energy_bins": total_bins,
        "collision_energy_bins": [
            {
                "key": key,
                "value": bins[key][0].get("ce_value"),
                "unit": bins[key][0].get("ce_unit"),
                "technical_replicates": len(bins[key]),
                "condition_evidence_count": 1,
                "independent_sample_vote_count": 0,
            }
            for key in bin_keys
        ],
        "decision": {
            "status": decision_status,
            "scope": decision_scope,
            "winner_margin": round(margin, 6),
            "minimum_score": minimum,
            "minimum_margin": margin_required,
            "blockers": blockers,
            "query_spectrum_release_gate": query_quality_audit,
            "cross_platform_abstention_gate": abstention_audit,
            "stereochemistry_claimed": False,
            "score_is_probability": False,
        },
        "winner": public_winner,
        "leading_candidate_when_abstained": public_leading if winner is None else None,
        "candidates": public_candidates,
        "candidate_count": len(candidates),
        "vetoed_candidate_count": sum(bool(candidate["vetoed"]) for candidate in candidates),
        "ranking_changed_from_scan_top1": bool(
            public_winner
            and baseline_identity
            and public_winner["identity_key"] != baseline_identity
        ),
        "ranking_changed_from_v14": bool(
            public_winner
            and legacy_leading
            and public_winner["identity_key"] != legacy_leading["identity_key"]
        ),
        "v14_audit_baseline": {
            "identity_key": (
                legacy_leading["identity_key"] if legacy_leading else None
            ),
            "entity_score": (
                legacy_leading["legacy_v14_entity_score"]
                if legacy_leading else None
            ),
            "rank": 1 if legacy_leading else None,
        },
        "scan_top1_audit_baseline": {
            "selection": "highest_tic_then_spectrum_id_then_lowest_result_index",
            "result_index": (
                int(baseline_scan["row"]["result_index"]) if baseline_scan else None
            ),
            "spectrum_id": (
                str(baseline_scan["row"].get("spectrum_id") or "")
                if baseline_scan else None
            ),
            "identity_key": baseline_identity,
        },
    }


def build_flagship_result(
    task_id: str,
    rows: Iterable[Mapping[str, Any]],
    *,
    policy: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build one deterministic V11 result without mutating legacy rows."""

    selected_policy = dict(policy or load_flagship_policy())
    groups: dict[str, dict[str, Any]] = {}
    source_digest = hashlib.sha256()
    source_rows = 0
    matched_rows = 0
    for row_index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            continue
        source_rows += 1
        entity_id, source, source_value = _entity_key(task_id, row, row_index)
        group = groups.setdefault(
            entity_id,
            {
                "entity_id": entity_id,
                "entity_source": source,
                "entity_source_value": source_value,
                "scans": [],
            },
        )
        candidates = _candidate_rows(row, selected_policy)
        matched_rows += int(bool(candidates))
        ce_key, ce_value, ce_unit = _collision_energy_bin(row, selected_policy)
        compact = _compact_row(row, row_index)
        group["scans"].append({
            "row": compact,
            "ce_key": ce_key,
            "ce_value": ce_value,
            "ce_unit": ce_unit,
            "candidates": candidates,
        })
        source_digest.update(_canonical(_source_projection(compact, candidates, entity_id)))
        source_digest.update(b"\n")

    entities = [
        _finalize_entity(groups[entity_id], selected_policy)
        for entity_id in sorted(groups)
    ]
    decisions = Counter(entity["decision"]["status"] for entity in entities)
    scopes = Counter(entity["decision"]["scope"] for entity in entities)
    fragmentation = Counter(
        str(((entity.get("winner") or {}).get("fragmentation_explanation") or {}).get("status") or "not_applicable")
        for entity in entities
    )
    abstention_gate_statuses = Counter(
        str(
            (entity.get("decision") or {})
            .get("cross_platform_abstention_gate", {})
            .get("status")
            or "not_applicable"
        )
        for entity in entities
    )
    query_quality_gate_statuses = Counter(
        str(
            (entity.get("decision") or {})
            .get("query_spectrum_release_gate", {})
            .get("status")
            or "not_applicable"
        )
        for entity in entities
    )
    registry_policy = selected_policy.get("identity_registry") or {}
    registry_id = ""
    if registry_policy.get("sha256"):
        registry_id = _identity_registry(str(registry_policy["sha256"]))[0]
    return {
        "metadata": {
            "schema_version": str(
                selected_policy.get("schema_version") or FLAGSHIP_SCHEMA_VERSION
            ),
            "module_id": str(selected_policy.get("module_id") or MODULE_ID),
            "algorithm_id": selected_policy["algorithm_id"],
            "algorithm_name": str(
                selected_policy.get("display_name") or "V11-Flagship"
            ),
            "policy_sha256": flagship_policy_sha256(selected_policy),
            "task_id": task_id,
            "source_row_count": source_rows,
            "source_matched_row_count": matched_rows,
            "source_semantic_sha256": source_digest.hexdigest(),
            "entity_count": len(entities),
            "entities_with_candidate": sum(bool(entity.get("winner")) for entity in entities),
            "entities_abstained": sum(entity["decision"]["status"] == "abstained" for entity in entities),
            "entities_unmatched": sum(entity["decision"]["status"] == "unmatched" for entity in entities),
            "entities_reranked_from_first_scan": sum(bool(entity.get("ranking_changed_from_scan_top1")) for entity in entities),
            "entities_reranked_from_v14": sum(
                bool(entity.get("ranking_changed_from_v14")) for entity in entities
            ),
            "entities_with_conditioned_multilibrary_influence": sum(
                any(
                    bool(candidate.get("conditioned_multilibrary_runtime_influence"))
                    for candidate in entity.get("candidates") or []
                )
                for entity in entities
            ),
            "decision_status_counts": dict(sorted(decisions.items())),
            "decision_scope_counts": dict(sorted(scopes.items())),
            "fragmentation_status_counts": dict(sorted(fragmentation.items())),
            "query_spectrum_release_gate_status_counts": dict(
                sorted(query_quality_gate_statuses.items())
            ),
            "entities_demoted_by_query_spectrum_release_gate": int(
                query_quality_gate_statuses.get("demoted", 0)
            ),
            "entities_failed_query_spectrum_audit_only": int(
                query_quality_gate_statuses.get("failed_audit_only", 0)
            ),
            "cross_platform_abstention_gate_status_counts": dict(
                sorted(abstention_gate_statuses.items())
            ),
            "entities_demoted_by_cross_platform_abstention_gate": int(
                abstention_gate_statuses.get("demoted", 0)
            ),
            "score_is_probability": False,
            "legacy_result_mutated": False,
            "production_identity_selector": True,
            "identity_registry_id": registry_id,
            "identity_registry_sha256": registry_policy.get("sha256"),
            "v14_audit_baseline_retained": bool(
                (selected_policy.get("conditioned_multilibrary") or {}).get(
                    "enabled"
                )
            ),
        },
        "policy": selected_policy,
        "invariants": {
            "one_decision_per_candidate_independent_entity": True,
            "technical_replicates_are_independent_votes": False,
            "distinct_collision_energies_are_independent_samples": False,
            "distinct_collision_energies_are_separate_measurement_conditions": True,
            "expert_score_affects_ranking": False,
            "classification_affects_ranking": False,
            "feature_fragment_evidence_affects_ranking": False,
            "biosynthetic_pathway_affects_ranking": False,
            "candidate_pathway_affects_ranking": False,
            "conditioned_multilibrary_affects_ranking": bool(
                (selected_policy.get("conditioned_multilibrary") or {}).get(
                    "enabled"
                )
            ),
            "source_count_bonus_applied": False,
            "mirrored_reference_sources_are_independent_votes": False,
            "source_accessions_are_independent_votes": False,
            "cross_unit_collision_energy_conversion_applied": False,
            "missing_chemistry_is_negative_evidence": False,
            "missing_from_truncated_topk_is_negative_evidence": bool(
                (selected_policy.get("candidate_pool") or {}).get(
                    "missing_from_truncated_topk_is_negative_evidence", False
                )
            ),
            "candidate_independent_formula_adduct_contradiction_can_veto": True,
            "candidate_library_adduct_contradiction_can_veto": False,
            "query_spectrum_release_gate_changes_top1": False,
            "query_spectrum_release_gate_changes_top5": False,
            "query_spectrum_release_gate_changes_entity_score": False,
            "query_spectrum_release_gate_affects_release_scope": bool(
                (selected_policy.get("query_spectrum_release_gate") or {}).get(
                    "affects_release_scope"
                )
            ),
            "technical_replicates_create_query_quality_conditions": False,
            "cross_platform_abstention_gate_changes_top1": False,
            "cross_platform_abstention_gate_changes_top5": False,
            "cross_platform_abstention_gate_changes_entity_score": False,
            "cross_platform_abstention_gate_independently_validated": False,
            "identity_registry_uses_ground_truth_inputs": False,
            "tautomer_standardization_enabled": bool(
                ((selected_policy.get("candidate_pool") or {}).get(
                    "tautomer_standardization"
                ) or {}).get("enabled", True)
            ),
            "mass_only_neutral_loss_called_mechanism": False,
            "stereochemistry_claimed": False,
        },
        "entities": entities,
    }


__all__ = [
    "FLAGSHIP_SCHEMA_VERSION",
    "MODULE_ID",
    "build_flagship_result",
    "connectivity_identity",
    "default_flagship_policy_path",
    "flagship_policy_sha256",
    "load_flagship_policy",
]
