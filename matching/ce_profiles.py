"""Read-only multi-source collision-energy profile evidence.

The module is intentionally shadow-only.  It summarizes whether a candidate
has a structure/adduct/instrument/unit-isolated multi-energy profile and how a
query compares with the nearest CE consensus; it never changes production
scores or ranking.
"""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
import sqlite3
from typing import Any, Iterable

from config.settings import CE_PROFILE_DB, CE_PROFILE_MAX_NEAR_DELTA
from matching.library_quality import normalize_inchikey
from matching.rule_quality import (
    collision_energy_value_unit,
    normalize_adduct,
    normalize_polarity,
)
from matching.scores import cosine_score, entropy_score
from utils.instrument_metadata import normalize_instrument_family


def _path_and_stamp(path: Path | None) -> tuple[str, int]:
    resolved = Path(path or CE_PROFILE_DB).resolve()
    try:
        stamp = resolved.stat().st_mtime_ns
    except OSError:
        stamp = 0
    return str(resolved), stamp


def _json(value: Any, fallback: Any) -> Any:
    try:
        parsed = json.loads(value or "")
    except (TypeError, ValueError, json.JSONDecodeError):
        return fallback
    return parsed


@lru_cache(maxsize=20_000)
def _profile_rows(
    path_text: str,
    stamp: int,
    inchikey: str,
    adduct: str,
    ion_mode: str,
    ce_unit: str,
) -> tuple[tuple[Any, ...], ...]:
    del stamp
    path = Path(path_text)
    if not path.is_file():
        return ()
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        columns = {
            str(row[1])
            for row in connection.execute("PRAGMA table_info(CEProfiles)")
        }
        source_column = (
            "source_dataset"
            if "source_dataset" in columns
            else "'MassBank_NIST'"
        )
        verification_column = (
            "identity_verification"
            if "identity_verification" in columns
            else "'structure_verified'"
        )
        verified_count_column = (
            "structure_verified_spectra"
            if "structure_verified_spectra" in columns
            else "n_spectra"
        )
        rows = connection.execute(
            f"""
            SELECT id, instrument_type, instrument_model, ce_unit,
                   n_energy_levels, n_spectra, ce_min, ce_max, quality_grade,
                   {source_column}, {verification_column},
                   {verified_count_column}
            FROM CEProfiles
            WHERE inchikey=? AND adduct=? AND ion_mode=? AND ce_unit=?
            ORDER BY
                CASE quality_grade WHEN 'A' THEN 0 WHEN 'B' THEN 1 ELSE 2 END,
                n_energy_levels DESC, n_spectra DESC, id
            """,
            (inchikey, adduct, ion_mode, ce_unit),
        ).fetchall()
        return tuple(tuple(row) for row in rows)
    finally:
        connection.close()


@lru_cache(maxsize=20_000)
def _profile_sources(
    path_text: str,
    stamp: int,
    profile_id: int,
) -> tuple[tuple[Any, ...], ...]:
    del stamp
    path = Path(path_text)
    if not path.is_file():
        return ()
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        has_table = connection.execute(
            "SELECT 1 FROM sqlite_master "
            "WHERE type='table' AND name='ProfileSources'"
        ).fetchone()
        if not has_table:
            return ()
        rows = connection.execute(
            """
            SELECT source_dataset, source_file, n_spectra
            FROM ProfileSources
            WHERE profile_id=?
            ORDER BY n_spectra DESC, source_dataset, source_file
            """,
            (profile_id,),
        ).fetchall()
        return tuple(tuple(row) for row in rows)
    finally:
        connection.close()


@lru_cache(maxsize=50_000)
def _nearest_bin(
    path_text: str,
    stamp: int,
    profile_id: int,
    query_energy: float,
) -> tuple[Any, ...] | None:
    del stamp
    path = Path(path_text)
    if not path.is_file():
        return None
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        row = connection.execute(
            """
            SELECT ce_value, n_spectra, median_precursor_mz,
                   consensus_peaks_json, accessions_json
            FROM CEBins
            WHERE profile_id=?
            ORDER BY ABS(ce_value - ?), n_spectra DESC, ce_value
            LIMIT 1
            """,
            (profile_id, query_energy),
        ).fetchone()
        return tuple(row) if row else None
    finally:
        connection.close()


@lru_cache(maxsize=20_000)
def _trajectory_rows(
    path_text: str,
    stamp: int,
    profile_id: int,
) -> tuple[tuple[Any, ...], ...]:
    del stamp
    path = Path(path_text)
    if not path.is_file():
        return ()
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        rows = connection.execute(
            """
            SELECT fragment_mz, energy_intensity_json, energy_frequency_json,
                   n_levels, trend, apex_energy
            FROM CEFragmentTrajectories
            WHERE profile_id=?
            ORDER BY n_levels DESC, fragment_mz
            """,
            (profile_id,),
        ).fetchall()
        return tuple(tuple(row) for row in rows)
    finally:
        connection.close()


def _matched_trajectories(
    query_peaks: Iterable[Iterable[Any]],
    trajectories: Iterable[tuple[Any, ...]],
    *,
    tolerance: float = 0.05,
    limit: int = 12,
) -> list[dict[str, Any]]:
    query_mz: list[float] = []
    for peak in query_peaks or []:
        try:
            query_mz.append(float(peak[0]))
        except (TypeError, ValueError, IndexError):
            continue
    matches: list[dict[str, Any]] = []
    for row in trajectories:
        fragment_mz = float(row[0])
        if not any(abs(fragment_mz - observed) <= tolerance for observed in query_mz):
            continue
        matches.append(
            {
                "fragment_mz": fragment_mz,
                "energy_intensity": _json(row[1], {}),
                "energy_frequency": _json(row[2], {}),
                "n_levels": int(row[3]),
                "trend": str(row[4]),
                "apex_energy": row[5],
            }
        )
        if len(matches) >= limit:
            break
    return matches


def lookup_ce_profile_evidence(
    *,
    inchikey: Any,
    adduct: Any,
    ion_mode: Any,
    query_collision_energy: Any,
    query_peaks: Iterable[Iterable[Any]] | None = None,
    query_precursor_mz: float | None = None,
    query_instrument_type: Any = None,
    query_instrument: Any = None,
    db_path: Path | None = None,
) -> dict[str, Any]:
    """Return auditable nearest-CE consensus evidence for one candidate."""
    query_peak_list = [] if query_peaks is None else list(query_peaks)
    evidence: dict[str, Any] = {
        "status": "unavailable",
        "runtime_influence": False,
        "score_contribution": 0.0,
    }
    path_text, stamp = _path_and_stamp(db_path)
    if stamp == 0:
        evidence["status"] = "profile_database_missing"
        return evidence
    exact_key = normalize_inchikey(inchikey)
    if not exact_key:
        evidence["status"] = "candidate_identity_unverified"
        return evidence
    normalized_adduct = normalize_adduct(adduct)
    polarity = {
        "positive": "pos",
        "negative": "neg",
    }.get(normalize_polarity(ion_mode) or "", "")
    if not normalized_adduct or polarity not in {"pos", "neg"}:
        evidence["status"] = "candidate_adduct_or_polarity_unknown"
        return evidence
    query_energy, query_unit = collision_energy_value_unit(
        query_collision_energy
    )
    if query_energy is None:
        evidence["status"] = "query_collision_energy_missing"
        return evidence
    if query_unit not in {"eV", "NCE"}:
        evidence["status"] = "query_collision_energy_unit_unknown"
        evidence["query_collision_energy"] = query_energy
        return evidence

    rows = _profile_rows(
        path_text,
        stamp,
        exact_key,
        normalized_adduct,
        polarity,
        query_unit,
    )
    if not rows:
        evidence.update(
            {
                "status": "no_comparable_profile",
                "query_collision_energy": query_energy,
                "query_collision_energy_unit": query_unit,
            }
        )
        return evidence

    candidates: list[dict[str, Any]] = []
    query_instrument_family = normalize_instrument_family(
        query_instrument_type, query_instrument
    )
    for profile in rows[:8]:
        profile_id = int(profile[0])
        energy_bin = _nearest_bin(
            path_text, stamp, profile_id, round(float(query_energy), 6)
        )
        if energy_bin is None:
            continue
        consensus_payload = _json(energy_bin[3], [])
        consensus_peaks = [
            [float(item["mz"]), float(item["intensity"])]
            for item in consensus_payload
            if isinstance(item, dict)
            and item.get("mz") is not None
            and item.get("intensity") is not None
        ]
        closest_energy = float(energy_bin[0])
        delta = abs(closest_energy - float(query_energy))
        library_precursor = (
            float(energy_bin[2]) if energy_bin[2] is not None else None
        )
        candidate = {
            "profile_id": profile_id,
            "instrument_type": str(profile[1]),
            "instrument_model": str(profile[2]),
            "instrument_family": normalize_instrument_family(
                profile[1], profile[2]
            ) or "unknown",
            "collision_energy_unit": str(profile[3]),
            "n_energy_levels": int(profile[4]),
            "n_profile_spectra": int(profile[5]),
            "energy_min": float(profile[6]),
            "energy_max": float(profile[7]),
            "quality_grade": str(profile[8]),
            "source_dataset": str(profile[9]),
            "identity_verification": str(profile[10]),
            "structure_verified_spectra": int(profile[11]),
            "closest_energy": closest_energy,
            "energy_delta": delta,
            "near_energy": delta <= CE_PROFILE_MAX_NEAR_DELTA,
            "consensus_spectra": int(energy_bin[1]),
            "consensus_cosine": (
                cosine_score(
                    query_peak_list,
                    consensus_peaks,
                    mz_tol=0.05,
                    query_pmz=query_precursor_mz,
                    lib_pmz=library_precursor,
                )
                if consensus_peaks and query_peak_list
                else None
            ),
            "consensus_entropy": (
                entropy_score(
                    query_peak_list, consensus_peaks, mz_tol=0.05
                )
                if consensus_peaks and query_peak_list
                else None
            ),
            "source_accessions": _json(energy_bin[4], [])[:12],
            "source_provenance": [
                {
                    "source_dataset": str(source[0]),
                    "source_file": str(source[1]),
                    "n_spectra": int(source[2]),
                }
                for source in _profile_sources(
                    path_text, stamp, profile_id
                )[:12]
            ],
            "matched_fragment_trajectories": _matched_trajectories(
                query_peak_list,
                _trajectory_rows(path_text, stamp, profile_id),
            ),
        }
        candidate["same_instrument_family_as_query"] = bool(
            query_instrument_family
            and candidate["instrument_family"] == query_instrument_family
        )
        candidates.append(candidate)

    if not candidates:
        evidence["status"] = "profile_bins_missing"
        return evidence
    candidates.sort(
        key=lambda item: (
            (
                0
                if item["same_instrument_family_as_query"] and item["near_energy"]
                else 1
                if item["same_instrument_family_as_query"]
                else 2
                if item["near_energy"]
                else 3
            )
            if query_instrument_family else 0,
            item["energy_delta"],
            {"A": 0, "B": 1, "C": 2}.get(item["quality_grade"], 3),
            -item["n_energy_levels"],
            -item["n_profile_spectra"],
            item["profile_id"],
        )
    )
    best = candidates[0]
    evidence.update(
        {
            "status": (
                "available_near_energy"
                if best["near_energy"]
                else "available_distant_energy"
            ),
            "query_collision_energy": float(query_energy),
            "query_collision_energy_unit": query_unit,
            "query_instrument_family": query_instrument_family or "unknown",
            "profile_count": len(rows),
            "best_profile": best,
            "profiles": candidates[:5],
            "interpretation": (
                "Multi-source CE consensus evidence; shadow-only and not "
                "used to change final score, rank, or confidence."
            ),
        }
    )
    return evidence


def profile_database_summary(db_path: Path | None = None) -> dict[str, Any]:
    path_text, stamp = _path_and_stamp(db_path)
    path = Path(path_text)
    if stamp == 0:
        return {"available": False, "path": path_text}
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        meta = dict(connection.execute("SELECT key, value FROM Meta"))
        profile_count = connection.execute(
            "SELECT COUNT(*) FROM CEProfiles"
        ).fetchone()[0]
        identity_count = connection.execute(
            "SELECT COUNT(DISTINCT inchikey) FROM CEProfiles"
        ).fetchone()[0]
        bin_count = connection.execute("SELECT COUNT(*) FROM CEBins").fetchone()[0]
        trajectory_count = connection.execute(
            "SELECT COUNT(*) FROM CEFragmentTrajectories"
        ).fetchone()[0]
        by_unit = dict(
            connection.execute(
                "SELECT ce_unit, COUNT(*) FROM CEProfiles GROUP BY ce_unit"
            )
        )
        by_quality = dict(
            connection.execute(
                "SELECT quality_grade, COUNT(*) FROM CEProfiles "
                "GROUP BY quality_grade"
            )
        )
        profile_columns = {
            str(row[1])
            for row in connection.execute("PRAGMA table_info(CEProfiles)")
        }
        by_source = (
            dict(
                connection.execute(
                    "SELECT source_dataset, COUNT(*) FROM CEProfiles "
                    "GROUP BY source_dataset"
                )
            )
            if "source_dataset" in profile_columns
            else {"MassBank_NIST": int(profile_count)}
        )
        by_identity_verification = (
            dict(
                connection.execute(
                    "SELECT identity_verification, COUNT(*) FROM CEProfiles "
                    "GROUP BY identity_verification"
                )
            )
            if "identity_verification" in profile_columns
            else {"structure_verified": int(profile_count)}
        )
        source_file_count = connection.execute(
            "SELECT COUNT(*) FROM sqlite_master "
            "WHERE type='table' AND name='SourceFiles'"
        ).fetchone()[0]
        if source_file_count:
            source_file_count = connection.execute(
                "SELECT COUNT(*) FROM SourceFiles"
            ).fetchone()[0]
        return {
            "available": True,
            "path": path_text,
            "schema_version": meta.get("schema_version"),
            "profile_count": int(profile_count),
            "identity_count": int(identity_count),
            "energy_bin_count": int(bin_count),
            "fragment_trajectory_count": int(trajectory_count),
            "profiles_by_unit": by_unit,
            "profiles_by_quality": by_quality,
            "profiles_by_source": by_source,
            "profiles_by_identity_verification": by_identity_verification,
            "source_file_count": int(source_file_count),
            "build_summary": _json(meta.get("summary_json"), {}),
        }
    finally:
        connection.close()
