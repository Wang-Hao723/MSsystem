"""AuthenticAnchor-MS: conservative authentic-standard MS/MS evidence.

The module offers two explicit shadow paths.  Candidate-conditional evaluation
can compare identities already recalled by the ordinary search.  A separate
precursor-mass-filtered authentic-library search can recover identities absent
from public libraries, but it is never injected into the primary candidate
pool.  Both paths remain shadow-only until a frozen, identity-disjoint cohort
has validated a release policy.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Mapping, Sequence


AUTHENTIC_ANCHOR_SCHEMA_VERSION = "authentic-anchor-ms/1.0"

PROTON_MASS = 1.007276466621
SODIUM_ION_MASS = 22.989218
POTASSIUM_ION_MASS = 38.963158
AMMONIUM_ION_MASS = 18.033823
CHLORIDE_ION_MASS = 34.968853
FORMIC_ACID_MASS = 46.005479
ACETIC_ACID_MASS = 60.021129


@dataclass(frozen=True)
class AuthenticAnchorPolicy:
    """Frozen conservative gate for candidate-conditional anchor evidence."""

    fragment_mz_tolerance_da: float = 0.02
    collision_energy_tolerance_ev: float = 1.0
    minimum_collision_energy_votes: int = 2
    minimum_median_similarity: float = 0.50
    minimum_winner_margin: float = 0.05
    minimum_query_product_ions: int = 3
    require_baseline_separation: bool = True
    allow_primary_rank_change: bool = False

    def validated(self) -> "AuthenticAnchorPolicy":
        if not 0 < float(self.fragment_mz_tolerance_da) <= 0.2:
            raise ValueError("fragment_mz_tolerance_da must be in (0, 0.2]")
        if not 0 <= float(self.collision_energy_tolerance_ev) <= 10:
            raise ValueError("collision_energy_tolerance_ev must be in [0, 10]")
        if int(self.minimum_collision_energy_votes) < 2:
            raise ValueError("at least two collision-energy votes are required")
        if not 0 <= float(self.minimum_median_similarity) <= 1:
            raise ValueError("minimum_median_similarity must be in [0, 1]")
        if not 0 <= float(self.minimum_winner_margin) <= 1:
            raise ValueError("minimum_winner_margin must be in [0, 1]")
        if int(self.minimum_query_product_ions) < 3:
            raise ValueError("at least three measured product ions are required")
        if self.allow_primary_rank_change:
            raise ValueError(
                "AuthenticAnchor-MS v1 is shadow-only; primary rank changes are disabled"
            )
        return self


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _full_inchikey(value: Any) -> str:
    token = str(value or "").strip().upper().replace("INCHIKEY=", "")
    parts = token.split("-")
    if (
        len(parts) == 3
        and len(parts[0]) == 14
        and len(parts[1]) == 10
        and len(parts[2]) == 1
        and all(part.isalpha() for part in parts)
    ):
        return token
    return ""


def _connectivity(value: Any) -> str:
    full = _full_inchikey(value)
    if full:
        return full[:14]
    token = str(value or "").strip().upper().split("-", 1)[0]
    return token if len(token) == 14 and token.isalpha() else ""


def _normalized_peaks(
    peaks: Iterable[Sequence[Any]] | None,
) -> tuple[list[float], list[float]]:
    pairs: list[tuple[float, float]] = []
    for peak in peaks or ():
        try:
            mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
            pairs.append((mz, intensity))
    pairs.sort()
    norm = math.sqrt(sum(intensity * intensity for _, intensity in pairs))
    if not pairs or norm <= 0:
        return [], []
    return (
        [mz for mz, _ in pairs],
        [intensity / norm for _, intensity in pairs],
    )


def greedy_cosine_similarity(
    left: Iterable[Sequence[Any]] | None,
    right: Iterable[Sequence[Any]] | None,
    *,
    tolerance_da: float,
) -> float:
    """Return a deterministic one-to-one cosine similarity in ``[0, 1]``."""

    left_mz, left_intensity = _normalized_peaks(left)
    right_mz, right_intensity = _normalized_peaks(right)
    if not left_mz or not right_mz:
        return 0.0
    pairs: list[tuple[float, float, int, int]] = []
    for left_index, mass in enumerate(left_mz):
        for right_index, candidate_mass in enumerate(right_mz):
            delta = abs(mass - candidate_mass)
            if delta <= tolerance_da:
                pairs.append(
                    (
                        left_intensity[left_index] * right_intensity[right_index],
                        delta,
                        left_index,
                        right_index,
                    )
                )
    used_left: set[int] = set()
    used_right: set[int] = set()
    score = 0.0
    for product, _, left_index, right_index in sorted(
        pairs, key=lambda item: (-item[0], item[1], item[2], item[3])
    ):
        if left_index in used_left or right_index in used_right:
            continue
        used_left.add(left_index)
        used_right.add(right_index)
        score += product
    return max(0.0, min(1.0, score))


def _canonical_payload(payload: Mapping[str, Any]) -> bytes:
    clean = dict(payload)
    clean.pop("content_sha256", None)
    return json.dumps(
        clean,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def content_sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_payload(payload)).hexdigest()


def finalize_anchor_library(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Validate and bind a library artifact to a deterministic content hash."""

    clean = dict(payload)
    clean["schema_version"] = AUTHENTIC_ANCHOR_SCHEMA_VERSION
    clean["content_sha256"] = content_sha256(clean)
    validate_anchor_library(clean)
    return clean


def _policy_from_library(library: Mapping[str, Any]) -> AuthenticAnchorPolicy:
    raw = library.get("release_policy") or {}
    fields = set(AuthenticAnchorPolicy.__dataclass_fields__)
    policy = AuthenticAnchorPolicy(**{key: raw[key] for key in fields if key in raw})
    return policy.validated()


def validate_anchor_library(payload: Mapping[str, Any]) -> None:
    if payload.get("schema_version") != AUTHENTIC_ANCHOR_SCHEMA_VERSION:
        raise ValueError("unsupported AuthenticAnchor-MS schema")
    if payload.get("runtime_mode") != "shadow_only":
        raise ValueError("AuthenticAnchor-MS v1 artifacts must remain shadow_only")
    policy = _policy_from_library(payload)
    if asdict(policy) != dict(payload.get("release_policy") or {}):
        raise ValueError("release_policy must be complete and canonical")
    construction = payload.get("construction_policy") or {}
    if construction.get("explicit_task_id_only") is not True:
        raise ValueError("reference spectra require explicit task IDs")
    if construction.get("candidate_generation_allowed") is not False:
        raise ValueError("anchor evidence cannot alter the primary candidate pool")
    if construction.get("authentic_library_shadow_search_allowed") is not True:
        raise ValueError("artifact must declare its shadow-library search policy")
    entries = payload.get("entries") or []
    if not isinstance(entries, list) or not entries:
        raise ValueError("authentic-standard library is empty")
    seen: set[str] = set()
    for entry in entries:
        reference_id = str(entry.get("reference_id") or "").strip()
        if not reference_id or reference_id in seen:
            raise ValueError("reference_id values must be non-empty and unique")
        seen.add(reference_id)
        if entry.get("authentic_standard") is not True:
            raise ValueError(f"{reference_id}: authentic_standard must be true")
        if not str(entry.get("source_task_id") or "").strip():
            raise ValueError(f"{reference_id}: source_task_id is required")
        if not _full_inchikey(entry.get("inchikey")):
            raise ValueError(f"{reference_id}: a valid full InChIKey is required")
        if str(entry.get("polarity") or "").upper() not in {"POS", "NEG"}:
            raise ValueError(f"{reference_id}: polarity must be POS or NEG")
        spectra = entry.get("spectra") or []
        if not spectra:
            raise ValueError(f"{reference_id}: no reference spectra")
        for spectrum in spectra:
            if _finite(spectrum.get("collision_energy_ev")) is None:
                raise ValueError(f"{reference_id}: collision energy is required")
            if len(_normalized_peaks(spectrum.get("peaks"))[0]) < 3:
                raise ValueError(f"{reference_id}: reference spectrum has <3 peaks")
    expected_hash = str(payload.get("content_sha256") or "")
    if expected_hash and expected_hash != content_sha256(payload):
        raise ValueError("authentic-standard library content hash mismatch")


def write_anchor_library(path: str | Path, payload: Mapping[str, Any]) -> Path:
    """Write a validated artifact atomically."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    bound = finalize_anchor_library(payload)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(bound, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    return destination


def read_anchor_library(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_anchor_library(payload)
    return payload


def _query_peaks(spectrum: Mapping[str, Any]) -> Sequence[Sequence[Any]]:
    return (
        spectrum.get("query_peaks")
        or spectrum.get("raw_peaks")
        or spectrum.get("peaks")
        or ()
    )


def _score_entry(
    queries: Sequence[Mapping[str, Any]],
    entry: Mapping[str, Any],
    policy: AuthenticAnchorPolicy,
) -> dict[str, Any]:
    references_by_ce: dict[float, list[Mapping[str, Any]]] = defaultdict(list)
    for reference in entry.get("spectra") or ():
        ce = _finite(reference.get("collision_energy_ev"))
        if ce is not None:
            references_by_ce[round(ce, 3)].append(reference)
    per_ce: dict[float, list[float]] = defaultdict(list)
    for query in queries:
        peaks = _query_peaks(query)
        if len(_normalized_peaks(peaks)[0]) < policy.minimum_query_product_ions:
            continue
        query_ce = _finite(
            query.get("collision_energy_ev")
            if query.get("collision_energy_ev") is not None
            else query.get("collision_energy")
        )
        if query_ce is None or not references_by_ce:
            continue
        nearest = min(references_by_ce, key=lambda ce: abs(ce - query_ce))
        if abs(nearest - query_ce) > policy.collision_energy_tolerance_ev:
            continue
        similarities = [
            greedy_cosine_similarity(
                peaks,
                reference.get("peaks") or (),
                tolerance_da=policy.fragment_mz_tolerance_da,
            )
            for reference in references_by_ce[nearest]
        ]
        if similarities:
            per_ce[round(query_ce, 3)].append(max(similarities))
    ce_scores = {
        str(ce): median(values)
        for ce, values in sorted(per_ce.items())
        if values
    }
    return {
        "reference_id": entry.get("reference_id"),
        "compound_code": entry.get("compound_code"),
        "compound_name": entry.get("compound_name"),
        "inchikey": entry.get("inchikey"),
        "connectivity": _connectivity(entry.get("inchikey")),
        "median_similarity": median(ce_scores.values()) if ce_scores else 0.0,
        "collision_energy_scores": ce_scores,
        "collision_energy_count": len(ce_scores),
    }


def _actual_ce_similarity(
    query_peaks: Sequence[Sequence[Any]],
    references_by_ce: Mapping[float, Sequence[Mapping[str, Any]]],
    query_ce: float,
    policy: AuthenticAnchorPolicy,
) -> tuple[float | None, dict[str, Any]]:
    """Interpolate similarity scalars at one recorded query CE.

    This deliberately does not interpolate spectra and does not extrapolate
    beyond the measured reference-CE range.  It is an audit-only way to use a
    ramp-acquisition CE without pretending that one query spectrum supplied
    multiple independent CE votes.
    """

    points = {
        float(reference_ce): max(
            greedy_cosine_similarity(
                query_peaks,
                reference.get("peaks") or (),
                tolerance_da=policy.fragment_mz_tolerance_da,
            )
            for reference in references
        )
        for reference_ce, references in references_by_ce.items()
        if references
    }
    ordered = sorted(points.items())
    audit = {
        "query_collision_energy_ev": float(query_ce),
        "reference_collision_energy_scores": {
            str(ce): score for ce, score in ordered
        },
        "outside_reference_range_clamped": False,
    }
    if not ordered:
        return None, {**audit, "status": "no_reference_collision_energy"}
    exact = next(
        (
            (ce, score)
            for ce, score in ordered
            if math.isclose(ce, query_ce, rel_tol=0.0, abs_tol=1e-6)
        ),
        None,
    )
    if exact is not None:
        return exact[1], {
            **audit,
            "status": "exact_reference_collision_energy",
            "bracketing_reference_collision_energies_ev": [exact[0]],
        }
    if len(ordered) < 2 or not ordered[0][0] < query_ce < ordered[-1][0]:
        return None, {
            **audit,
            "status": "outside_measured_reference_collision_energy_range",
        }
    for (left_ce, left_score), (right_ce, right_score) in zip(
        ordered, ordered[1:]
    ):
        if left_ce < query_ce < right_ce:
            fraction = (query_ce - left_ce) / (right_ce - left_ce)
            return left_score + fraction * (right_score - left_score), {
                **audit,
                "status": "scalar_similarity_interpolated",
                "bracketing_reference_collision_energies_ev": [
                    left_ce,
                    right_ce,
                ],
                "interpolation_fraction": fraction,
            }
    raise AssertionError("recorded query CE was not bracketed")


def _score_entry_actual_ce_shadow(
    queries: Sequence[Mapping[str, Any]],
    entry: Mapping[str, Any],
    policy: AuthenticAnchorPolicy,
) -> dict[str, Any]:
    references_by_ce: dict[float, list[Mapping[str, Any]]] = defaultdict(list)
    for reference in entry.get("spectra") or ():
        ce = _finite(reference.get("collision_energy_ev"))
        if ce is not None:
            references_by_ce[round(ce, 3)].append(reference)
    per_ce: dict[float, list[float]] = defaultdict(list)
    alignment_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for query in queries:
        peaks = _query_peaks(query)
        if len(_normalized_peaks(peaks)[0]) < policy.minimum_query_product_ions:
            continue
        query_ce = _finite(
            query.get("collision_energy_ev")
            if query.get("collision_energy_ev") is not None
            else query.get("collision_energy")
        )
        if query_ce is None:
            continue
        similarity, alignment = _actual_ce_similarity(
            peaks,
            references_by_ce,
            query_ce,
            policy,
        )
        alignment_rows[str(round(query_ce, 3))].append(alignment)
        if similarity is not None:
            per_ce[round(query_ce, 3)].append(similarity)
    ce_scores = {
        str(ce): median(values)
        for ce, values in sorted(per_ce.items())
        if values
    }
    return {
        "reference_id": entry.get("reference_id"),
        "compound_code": entry.get("compound_code"),
        "compound_name": entry.get("compound_name"),
        "inchikey": entry.get("inchikey"),
        "connectivity": _connectivity(entry.get("inchikey")),
        "median_similarity": median(ce_scores.values()) if ce_scores else 0.0,
        "collision_energy_scores": ce_scores,
        "collision_energy_count": len(ce_scores),
        "collision_energy_alignment": dict(sorted(alignment_rows.items())),
    }


def _restrict_scores_to_shared_collision_energies(
    scores: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Make multi-candidate scores comparable on the same query CE rows.

    A missing candidate/CE comparison is not a zero similarity observation.
    Requiring the intersection prevents a partially covered reference from
    becoming an artificial low-scoring opponent.  The full per-candidate CE
    map remains available for audit, while ranking fields contain only shared
    query conditions.
    """

    if not scores:
        return [], []
    shared = set(str(key) for key in scores[0].get("collision_energy_scores", {}))
    for score in scores[1:]:
        shared.intersection_update(
            str(key) for key in score.get("collision_energy_scores", {})
        )
    shared_values = sorted(shared, key=float)
    aligned = []
    for score in scores:
        available = dict(score.get("collision_energy_scores") or {})
        comparable = {key: available[key] for key in shared_values}
        aligned.append({
            **dict(score),
            "available_collision_energy_scores": available,
            "collision_energy_scores": comparable,
            "collision_energy_count": len(comparable),
            "median_similarity": median(comparable.values()) if comparable else 0.0,
        })
    return aligned, shared_values


def _candidate_requested(entry: Mapping[str, Any], requested: set[str]) -> bool:
    keys = {
        str(entry.get("reference_id") or "").strip().upper(),
        str(entry.get("compound_code") or "").strip().upper(),
        _full_inchikey(entry.get("inchikey")),
        _connectivity(entry.get("inchikey")),
    }
    return bool(keys & requested)


def _ion_hypotheses(neutral_mass: float, polarity: str) -> list[tuple[str, float]]:
    if polarity == "NEG":
        return [
            ("[M-H]-", neutral_mass - PROTON_MASS),
            ("[M+Cl]-", neutral_mass + CHLORIDE_ION_MASS),
            ("[M+HCOO]-", neutral_mass + FORMIC_ACID_MASS - PROTON_MASS),
            ("[M+CH3COO]-", neutral_mass + ACETIC_ACID_MASS - PROTON_MASS),
            ("[M-2H]2-", (neutral_mass - 2 * PROTON_MASS) / 2),
        ]
    return [
        ("[M+H]+", neutral_mass + PROTON_MASS),
        ("[M+Na]+", neutral_mass + SODIUM_ION_MASS),
        ("[M+K]+", neutral_mass + POTASSIUM_ION_MASS),
        ("[M+NH4]+", neutral_mass + AMMONIUM_ION_MASS),
        ("[M+2H]2+", (neutral_mass + 2 * PROTON_MASS) / 2),
    ]


def _mass_filtered_reference_ids(
    library: Mapping[str, Any],
    *,
    precursor_mz: float,
    polarity: str,
    tolerance_ppm: float,
) -> tuple[list[str], list[dict[str, Any]]]:
    selected: list[str] = []
    audit: list[dict[str, Any]] = []
    for entry in library.get("entries") or ():
        if str(entry.get("polarity") or "").upper() != polarity:
            continue
        neutral_mass = _finite(entry.get("neutral_mass"))
        if neutral_mass is None:
            continue
        best = min(
            (
                (
                    abs((precursor_mz - theoretical) / theoretical * 1e6),
                    label,
                    theoretical,
                )
                for label, theoretical in _ion_hypotheses(neutral_mass, polarity)
                if theoretical > 0
            ),
            default=(math.inf, "", 0.0),
        )
        if best[0] <= tolerance_ppm:
            reference_id = str(entry.get("reference_id") or "")
            selected.append(reference_id)
            audit.append({
                "reference_id": reference_id,
                "adduct_hypothesis": best[1],
                "theoretical_mz": round(best[2], 8),
                "mass_error_abs_ppm": round(best[0], 6),
            })
    return selected, audit


def evaluate_authentic_anchor(
    query_spectra: Sequence[Mapping[str, Any]],
    library: Mapping[str, Any],
    *,
    candidate_identity_keys: Sequence[str],
    query_polarity: str,
    query_instrument_family: str,
    baseline_separated: bool | None,
) -> dict[str, Any]:
    """Evaluate authentic-standard evidence for an existing candidate set.

    ``candidate_identity_keys`` is mandatory by design.  This keeps the
    authentic-standard registry out of candidate generation and makes leakage
    at the call site auditable.
    """

    validate_anchor_library(library)
    policy = _policy_from_library(library)
    requested = {
        str(value or "").strip().upper()
        for value in candidate_identity_keys
        if str(value or "").strip()
    }
    if not requested:
        raise ValueError("candidate_identity_keys must come from an existing candidate set")
    polarity = str(query_polarity or "").strip().upper()
    instrument = str(query_instrument_family or "").strip().casefold()
    if polarity not in {"POS", "NEG"}:
        raise ValueError("query_polarity must be POS or NEG")
    if not instrument:
        raise ValueError("query_instrument_family is required")
    compatible_entries = []
    instrument_mismatch_count = 0
    for entry in library.get("entries") or ():
        if not _candidate_requested(entry, requested):
            continue
        if str(entry.get("polarity") or "").upper() != polarity:
            continue
        reference_instrument = str(entry.get("instrument_family") or "").strip().casefold()
        if not reference_instrument or reference_instrument != instrument:
            instrument_mismatch_count += 1
            continue
        compatible_entries.append(entry)
    scores, shared_ce_values = _restrict_scores_to_shared_collision_energies([
        _score_entry(query_spectra, entry, policy)
        for entry in compatible_entries
    ])
    # Structural keys are used only to serialize exact score ties.  A tied
    # score cannot receive a CE vote or pass the margin gate below, so neither
    # a database accession nor library insertion order can select a winner.
    scores.sort(key=lambda item: (
        -item["median_similarity"],
        str(item.get("connectivity") or ""),
        _full_inchikey(item.get("inchikey")),
    ))
    winner = scores[0] if scores else None
    runner_up = scores[1] if len(scores) > 1 else None
    margin = (
        winner["median_similarity"] - runner_up["median_similarity"]
        if winner and runner_up else 0.0
    )
    ce_votes: Counter[str] = Counter()
    ce_values = sorted({
        ce
        for score in scores
        for ce in score["collision_energy_scores"]
    }, key=float)
    collision_energy_tie_count = 0
    for ce in ce_values:
        eligible = [
            (score["reference_id"], score["collision_energy_scores"].get(ce, -1.0))
            for score in scores
        ]
        if eligible:
            best = max(item[1] for item in eligible)
            winners = [
                reference_id
                for reference_id, similarity in eligible
                if math.isclose(similarity, best, rel_tol=0.0, abs_tol=1e-12)
            ]
            if len(winners) == 1:
                ce_votes[winners[0]] += 1
            else:
                collision_energy_tie_count += 1
    winner_votes = ce_votes.get(winner["reference_id"], 0) if winner else 0
    connectivities = {
        score["connectivity"] for score in scores if score["connectivity"]
    }
    full_keys = {
        _full_inchikey(score["inchikey"]) for score in scores if score["inchikey"]
    }
    same_2d_stereoisomer_set = bool(
        len(full_keys) >= 2 and len(connectivities) < len(full_keys)
    )

    reasons: list[str] = []
    decision = "abstained"
    if len(scores) < 2:
        reasons.append("fewer_than_two_compatible_candidate_anchors")
        if instrument_mismatch_count:
            reasons.append("instrument_family_mismatch")
    elif same_2d_stereoisomer_set:
        decision = "isomer_group_only"
        reasons.append("same_2d_stereoisomers_not_resolved_by_achiral_msms")
    elif policy.require_baseline_separation and baseline_separated is not True:
        decision = "isomer_group_only"
        reasons.append("chromatography_not_proven_baseline_separated")
    elif winner["collision_energy_count"] < policy.minimum_collision_energy_votes:
        reasons.append("insufficient_collision_energy_coverage")
    elif winner_votes < policy.minimum_collision_energy_votes:
        reasons.append("insufficient_consistent_collision_energy_votes")
    elif winner["median_similarity"] < policy.minimum_median_similarity:
        reasons.append("winner_similarity_below_frozen_threshold")
    elif margin < policy.minimum_winner_margin:
        reasons.append("winner_margin_below_frozen_threshold")
    else:
        decision = "anchor_connectivity_supported"
        reasons.append("multi_ce_authentic_anchor_and_chromatography_gate_passed")

    return {
        "schema_version": AUTHENTIC_ANCHOR_SCHEMA_VERSION,
        "module": "AuthenticAnchor-MS",
        "runtime_mode": "shadow_only",
        "primary_rank_changed": False,
        "candidate_generation_used": False,
        "identity_scope": "2d_connectivity_only",
        "decision": decision,
        "reason_codes": reasons,
        "winner_reference_id": winner.get("reference_id") if winner else None,
        "winner_connectivity": winner.get("connectivity") if winner else None,
        "winner_margin": round(margin, 6),
        "winner_collision_energy_votes": int(winner_votes),
        "collision_energy_tie_count": collision_energy_tie_count,
        "shared_candidate_collision_energy_conditions": shared_ce_values,
        "baseline_separated": baseline_separated,
        "same_2d_stereoisomer_set": same_2d_stereoisomer_set,
        "instrument_mismatch_count": instrument_mismatch_count,
        "candidate_scores": scores,
        "collision_energy_votes": dict(ce_votes),
        "release_policy": asdict(policy),
        "library_content_sha256": library.get("content_sha256"),
    }


def evaluate_authentic_anchor_actual_ce_shadow(
    query_spectra: Sequence[Mapping[str, Any]],
    library: Mapping[str, Any],
    *,
    candidate_identity_keys: Sequence[str],
    query_polarity: str,
    query_instrument_family: str,
) -> dict[str, Any]:
    """Score recorded ramp-acquisition CE against measured reference CEs.

    The result is an audit/display sidecar only.  Scalar similarity may be
    interpolated between two measured reference CEs, but one query CE remains
    one vote.  The function cannot change either the primary rank or the
    EviTwin display rank.
    """

    validate_anchor_library(library)
    policy = _policy_from_library(library)
    requested = {
        str(value or "").strip().upper()
        for value in candidate_identity_keys
        if str(value or "").strip()
    }
    if not requested:
        raise ValueError("candidate_identity_keys must come from an existing candidate set")
    polarity = str(query_polarity or "").strip().upper()
    instrument = str(query_instrument_family or "").strip().casefold()
    if polarity not in {"POS", "NEG"}:
        raise ValueError("query_polarity must be POS or NEG")
    if not instrument:
        raise ValueError("query_instrument_family is required")
    compatible_entries = []
    instrument_mismatch_count = 0
    for entry in library.get("entries") or ():
        if not _candidate_requested(entry, requested):
            continue
        if str(entry.get("polarity") or "").upper() != polarity:
            continue
        reference_instrument = str(
            entry.get("instrument_family") or ""
        ).strip().casefold()
        if not reference_instrument or reference_instrument != instrument:
            instrument_mismatch_count += 1
            continue
        compatible_entries.append(entry)
    scores, shared_ce_values = _restrict_scores_to_shared_collision_energies([
        _score_entry_actual_ce_shadow(query_spectra, entry, policy)
        for entry in compatible_entries
    ])
    scores = [score for score in scores if score["collision_energy_count"] > 0]
    # Structural keys serialize exact score ties only.  Ties remain unresolved
    # below; names, accessions and library insertion order never select a lead.
    scores.sort(key=lambda item: (
        -item["median_similarity"],
        str(item.get("connectivity") or ""),
        _full_inchikey(item.get("inchikey")),
    ))
    winner = scores[0] if scores else None
    runner_up = scores[1] if len(scores) > 1 else None
    margin = (
        winner["median_similarity"] - runner_up["median_similarity"]
        if winner and runner_up
        else 0.0
    )
    ce_votes: Counter[str] = Counter()
    ce_values = sorted({
        ce
        for score in scores
        for ce in score["collision_energy_scores"]
    }, key=float)
    collision_energy_tie_count = 0
    for ce in ce_values:
        eligible = [
            (score["reference_id"], score["collision_energy_scores"].get(ce, -1.0))
            for score in scores
        ]
        best = max((item[1] for item in eligible), default=-1.0)
        winners = [
            reference_id
            for reference_id, similarity in eligible
            if math.isclose(similarity, best, rel_tol=0.0, abs_tol=1e-12)
        ]
        if len(winners) == 1:
            ce_votes[winners[0]] += 1
        elif winners:
            collision_energy_tie_count += 1
    winner_votes = ce_votes.get(winner["reference_id"], 0) if winner else 0
    full_keys = {
        _full_inchikey(score.get("inchikey"))
        for score in scores
        if score.get("inchikey")
    }
    connectivities = {
        score.get("connectivity") for score in scores if score.get("connectivity")
    }
    same_2d = bool(len(full_keys) >= 2 and len(connectivities) < len(full_keys))
    exact_tie = bool(
        winner
        and runner_up
        and math.isclose(
            float(winner["median_similarity"]),
            float(runner_up["median_similarity"]),
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    )
    if len(scores) < 2:
        decision = "abstained"
        reasons = ["fewer_than_two_actual_ce_scoreable_candidate_anchors"]
    elif same_2d:
        decision = "unresolved_actual_ce_candidate_group"
        reasons = ["same_2d_stereoisomers_not_resolved_by_achiral_msms"]
    elif exact_tie:
        decision = "unresolved_actual_ce_candidate_group"
        reasons = ["actual_ce_similarity_tie"]
    elif int(winner.get("collision_energy_count") or 0) < 2:
        decision = "single_ce_shadow_candidate_group"
        reasons = ["one_recorded_query_ce_cannot_supply_multi_ce_evidence"]
    elif winner_votes < 2:
        decision = "unresolved_actual_ce_candidate_group"
        reasons = ["fewer_than_two_consistent_actual_ce_votes"]
    else:
        decision = "multi_ce_actual_ce_shadow_candidate_group"
        reasons = ["actual_ce_shadow_scored_but_release_not_enabled"]
    return {
        "schema_version": AUTHENTIC_ANCHOR_SCHEMA_VERSION,
        "module": "AuthenticAnchor-MS",
        "runtime_mode": "shadow_only",
        "scoring_mode": "recorded_query_ce_scalar_similarity_interpolation",
        "decision": decision,
        "reason_codes": reasons,
        "winner_reference_id": winner.get("reference_id") if winner else None,
        "winner_connectivity": winner.get("connectivity") if winner else None,
        "winner_margin": round(margin, 6),
        "winner_collision_energy_votes": int(winner_votes),
        "collision_energy_tie_count": collision_energy_tie_count,
        "shared_candidate_collision_energy_conditions": shared_ce_values,
        "same_2d_stereoisomer_set": same_2d,
        "instrument_mismatch_count": instrument_mismatch_count,
        "candidate_scores": scores,
        "collision_energy_votes": dict(ce_votes),
        "query_collision_energy_condition_count": len(shared_ce_values),
        "reference_spectra_interpolated": False,
        "scalar_similarity_interpolation_used": any(
            row.get("status") == "scalar_similarity_interpolated"
            for score in scores
            for rows in (score.get("collision_energy_alignment") or {}).values()
            for row in rows
        ),
        "outside_reference_range_clamped": False,
        "candidate_generation_used": False,
        "primary_rank_changed": False,
        "evitwin_display_rank_changed": False,
        "ranking_use_allowed": False,
        "formal_level1_released": False,
        "compound_names_used_for_scoring": False,
        "database_accessions_used_for_scoring": False,
        "database_source_order_used_for_scoring": False,
        "truth_used_for_scoring": False,
        "library_content_sha256": library.get("content_sha256"),
    }


def search_authentic_anchor_library(
    query_spectra: Sequence[Mapping[str, Any]],
    library: Mapping[str, Any],
    *,
    precursor_mz: float,
    query_polarity: str,
    query_instrument_family: str,
    baseline_separated: bool | None,
    precursor_tolerance_ppm: float = 10.0,
) -> dict[str, Any]:
    """Search all mass-compatible authentic anchors as a separate shadow view.

    This path is useful when a public library failed to recall a compound.  It
    remains visibly separate from the primary candidate pool and therefore
    cannot silently convert a post-reveal laboratory optimization into a
    confirmatory result.
    """

    validate_anchor_library(library)
    observed = _finite(precursor_mz)
    if observed is None or observed <= 0:
        raise ValueError("precursor_mz must be a positive finite number")
    if not 0 < float(precursor_tolerance_ppm) <= 100:
        raise ValueError("precursor_tolerance_ppm must be in (0, 100]")
    polarity = str(query_polarity or "").strip().upper()
    if polarity not in {"POS", "NEG"}:
        raise ValueError("query_polarity must be POS or NEG")
    reference_ids, mass_audit = _mass_filtered_reference_ids(
        library,
        precursor_mz=observed,
        polarity=polarity,
        tolerance_ppm=float(precursor_tolerance_ppm),
    )
    if not reference_ids:
        return {
            "schema_version": AUTHENTIC_ANCHOR_SCHEMA_VERSION,
            "module": "AuthenticAnchor-MS",
            "runtime_mode": "shadow_only",
            "search_mode": "precursor_mass_filtered_authentic_library",
            "primary_rank_changed": False,
            "primary_candidate_pool_changed": False,
            "candidate_generation_used": True,
            "candidate_generation_source": "authentic_anchor_mass_filter_only",
            "identity_scope": "2d_connectivity_only",
            "decision": "abstained",
            "reason_codes": ["no_mass_compatible_authentic_anchor"],
            "mass_compatible_anchor_count": 0,
            "mass_filter_audit": [],
            "library_content_sha256": library.get("content_sha256"),
        }
    result = evaluate_authentic_anchor(
        query_spectra,
        library,
        candidate_identity_keys=reference_ids,
        query_polarity=polarity,
        query_instrument_family=query_instrument_family,
        baseline_separated=baseline_separated,
    )
    result.update({
        "search_mode": "precursor_mass_filtered_authentic_library",
        "primary_candidate_pool_changed": False,
        "candidate_generation_used": True,
        "candidate_generation_source": "authentic_anchor_mass_filter_only",
        "mass_compatible_anchor_count": len(reference_ids),
        "mass_filter_audit": mass_audit,
        "precursor_mz": observed,
        "precursor_tolerance_ppm": float(precursor_tolerance_ppm),
    })
    return result


def search_authentic_anchor_feature(
    query_spectra: Sequence[Mapping[str, Any]],
    library: Mapping[str, Any],
    *,
    query_polarity: str,
    query_instrument_family: str,
    baseline_separated: bool | None,
    precursor_tolerance_ppm: float = 10.0,
    precursor_group_tolerance_da: float = 0.01,
) -> dict[str, Any]:
    """Search each measured precursor/adduct group for one LC feature.

    A feature may yield MS/MS at ``[M-H]-`` and an adduct such as
    ``[M+HCOO]-``.  Averaging those precursor masses creates a fictitious ion;
    this helper clusters the measured precursor groups before any search.
    """

    if not 0 < float(precursor_group_tolerance_da) <= 0.2:
        raise ValueError("precursor_group_tolerance_da must be in (0, 0.2]")
    observed_rows = []
    for spectrum in query_spectra:
        observed = _finite(
            spectrum.get("search_precursor_mz")
            if spectrum.get("search_precursor_mz") is not None
            else spectrum.get("precursor_mz")
        )
        if observed is not None and observed > 0:
            observed_rows.append((observed, spectrum))
    if not observed_rows:
        raise ValueError("query_spectra contain no measured precursor m/z")
    groups: list[list[tuple[float, Mapping[str, Any]]]] = []
    for observed, spectrum in sorted(observed_rows, key=lambda item: item[0]):
        if (
            not groups
            or abs(observed - median(value for value, _ in groups[-1]))
            > precursor_group_tolerance_da
        ):
            groups.append([])
        groups[-1].append((observed, spectrum))
    searches = []
    for group in groups:
        searches.append(search_authentic_anchor_library(
            [spectrum for _, spectrum in group],
            library,
            precursor_mz=median(value for value, _ in group),
            query_polarity=query_polarity,
            query_instrument_family=query_instrument_family,
            baseline_separated=baseline_separated,
            precursor_tolerance_ppm=precursor_tolerance_ppm,
        ))
    decision_priority = {
        "anchor_connectivity_supported": 2,
        "isomer_group_only": 1,
        "abstained": 0,
    }

    def strength(item: Mapping[str, Any]) -> tuple[Any, ...]:
        scores = item.get("candidate_scores") or []
        best_similarity = (
            float(scores[0].get("median_similarity") or 0.0) if scores else 0.0
        )
        return (
            decision_priority.get(str(item.get("decision") or ""), -1),
            best_similarity,
            float(item.get("winner_margin") or 0.0),
            int(item.get("winner_collision_energy_votes") or 0),
            -float(item.get("precursor_mz") or 0.0),
        )

    selected = dict(max(searches, key=strength))
    selected.update({
        "search_mode": "feature_precursor_grouped_authentic_library",
        "precursor_group_count": len(searches),
        "selected_precursor_group_mz": selected.get("precursor_mz"),
        "precursor_group_searches": [
            {
                "precursor_mz": item.get("precursor_mz"),
                "query_spectrum_count": len(group),
                "decision": item.get("decision"),
                "winner_reference_id": item.get("winner_reference_id"),
                "winner_margin": item.get("winner_margin"),
            }
            for item, group in zip(searches, groups)
        ],
    })
    return selected


__all__ = [
    "AUTHENTIC_ANCHOR_SCHEMA_VERSION",
    "AuthenticAnchorPolicy",
    "content_sha256",
    "evaluate_authentic_anchor",
    "evaluate_authentic_anchor_actual_ce_shadow",
    "finalize_anchor_library",
    "greedy_cosine_similarity",
    "read_anchor_library",
    "search_authentic_anchor_feature",
    "search_authentic_anchor_library",
    "validate_anchor_library",
    "write_anchor_library",
]
