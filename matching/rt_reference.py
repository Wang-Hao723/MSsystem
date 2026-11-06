"""Local RT standard evidence for v5 confidence calibration.

Retention time is only comparable inside the same LC method.  This module
therefore consumes a user-curated local standards table rather than treating
public-library RT as universal truth.
"""
from __future__ import annotations

import csv
import re
from functools import lru_cache
from pathlib import Path
from typing import Optional

from matching.library_quality import normalize_inchikey

try:
    from config.settings import (
        RT_STANDARDS_PATH, QUERY_RT_UNIT, RT_DEFAULT_TOL_SEC, RT_STRONG_TOL_SEC,
    )
except Exception:
    RT_STANDARDS_PATH = Path("data/rt_standards.csv")
    QUERY_RT_UNIT = "auto"
    RT_DEFAULT_TOL_SEC = 30.0
    RT_STRONG_TOL_SEC = 12.0


_NA = {"", "nan", "none", "null", "na", "n/a"}


def _clean(v) -> str:
    if v is None:
        return ""
    try:
        if v != v:
            return ""
    except Exception:
        pass
    s = str(v).strip()
    return "" if s.lower() in _NA else s


def _ik14(v) -> str:
    key = normalize_inchikey(v)
    return key.split("-", 1)[0] if key else ""


def _full_inchikey(v) -> str:
    return normalize_inchikey(v)


def _name_key(v) -> str:
    s = _clean(v).lower()
    s = re.sub(r"\([^)]*\)|\[[^]]*\]", " ", s)
    return re.sub(r"[^a-z0-9]+", "", s)


def _method_key(v) -> str:
    return re.sub(r"[^a-z0-9]+", "", _clean(v).lower())


def _bool(v) -> bool:
    return _clean(v).lower() in {"1", "true", "yes", "y", "on"}


def _float(v) -> Optional[float]:
    s = _clean(v)
    if not s:
        return None
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _seconds(value, unit: str = "auto") -> Optional[float]:
    x = _float(value)
    if x is None:
        return None
    u = (unit or "auto").strip().lower()
    if u.startswith("sec") or u in {"s", "second", "seconds"}:
        return x
    if u.startswith("min") or u in {"m", "minute", "minutes"}:
        return x * 60.0
    # LC retention times under ~80 are almost always minutes; MGF RTINSECONDS
    # and vendor exports are usually hundreds of seconds.
    return x if x > 80.0 else x * 60.0


def normalize_query_rt(rt) -> Optional[float]:
    return _seconds(rt, QUERY_RT_UNIT)


def _first(row: dict, names: tuple[str, ...]) -> str:
    lower = {str(k).strip().lower(): v for k, v in row.items()}
    for n in names:
        if n in lower:
            return _clean(lower[n])
    return ""


@lru_cache(maxsize=4)
def load_rt_standards(path: str | Path | None = None) -> list[dict]:
    p = Path(path or RT_STANDARDS_PATH)
    if not p.exists():
        return []
    rows: list[dict] = []
    with p.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            unit = _first(raw, ("rt_unit", "unit")) or "auto"
            rt_raw = _first(raw, ("rt_min", "retention_time_min", "rt", "retention_time"))
            rt_sec_raw = _first(raw, ("rt_sec", "retention_time_sec"))
            rt_sec = _seconds(rt_sec_raw, "second") if rt_sec_raw else _seconds(rt_raw, unit)
            if rt_sec is None:
                continue
            tol_raw = _first(raw, ("tolerance_sec", "rt_tolerance_sec", "window_sec"))
            tol_min_raw = _first(raw, ("tolerance_min", "rt_tolerance_min", "window_min"))
            tol_sec = _seconds(tol_raw, "second") if tol_raw else _seconds(tol_min_raw, "minute")
            if tol_sec is None:
                tol_sec = float(RT_DEFAULT_TOL_SEC)
            name = _first(raw, ("compound_name", "name", "standard_name", "compound"))
            full_inchikey = _full_inchikey(
                _first(raw, ("inchikey", "inchi_key"))
            )
            method_id = _first(
                raw,
                (
                    "method_id", "lc_method_id", "chromatography_method_id",
                ),
            )
            separation_mode = _first(
                raw,
                (
                    "separation_mode", "lc_separation_mode",
                    "stationary_phase_mode",
                ),
            )
            rows.append({
                "compound_name": name,
                "name_key": _name_key(name),
                "inchikey": full_inchikey,
                "inchikey14": (
                    full_inchikey.split("-", 1)[0]
                    if full_inchikey else ""
                ),
                "formula": _first(raw, ("formula", "molecular_formula")).upper(),
                "rt_sec": float(rt_sec),
                "tolerance_sec": float(tol_sec),
                "method_id": method_id,
                "method_key": _method_key(method_id),
                "separation_mode": separation_mode,
                "stereochemistry_capable": bool(
                    _bool(_first(raw, (
                        "stereochemistry_capable",
                        "chiral_separation",
                    )))
                    or "chiral" in separation_mode.lower()
                ),
                "source": _first(raw, ("source", "batch", "method")) or "local_rt_standard",
                "note": _first(raw, ("note", "notes", "comment")),
            })
    return rows


def _candidate_keys(match: dict) -> dict:
    full_inchikey = _full_inchikey(match.get("inchikey"))
    return {
        "inchikey": full_inchikey,
        "inchikey14": (
            full_inchikey.split("-", 1)[0]
            if full_inchikey else ""
        ),
        "name_key": _name_key(match.get("compound_name")),
        "formula": _clean(match.get("formula")).upper(),
    }


def match_rt_evidence(
    match: dict,
    query_rt,
    path: str | Path | None = None,
    *,
    query_method_id: str | None = None,
) -> Optional[dict]:
    """Return method-scoped local-standard RT evidence.

    RT is not portable across LC methods.  Exact-structure promotion therefore
    requires both an exact full InChIKey and an explicit matching method ID.
    Connectivity/name-only matches remain useful audit evidence but cannot
    resolve stereochemical identity.
    """
    q_rt_sec = normalize_query_rt(query_rt)
    if q_rt_sec is None or not match:
        return None
    keys = _candidate_keys(match)
    query_method_key = _method_key(query_method_id)
    best = None
    for row in load_rt_standards(path):
        matched_by = None
        if (
            keys["inchikey"]
            and row["inchikey"]
            and keys["inchikey"] == row["inchikey"]
        ):
            matched_by = "full_inchikey"
        elif (
            keys["inchikey14"]
            and row["inchikey14"]
            and keys["inchikey14"] == row["inchikey14"]
        ):
            matched_by = "connectivity_inchikey14"
        elif keys["name_key"] and row["name_key"] and keys["name_key"] == row["name_key"]:
            if (
                not row["formula"]
                or not keys["formula"]
                or row["formula"] == keys["formula"]
            ):
                matched_by = "exact_name_formula"
        if not matched_by:
            continue
        delta = q_rt_sec - row["rt_sec"]
        abs_delta = abs(delta)
        chromatographic_status = (
            "confirmed" if abs_delta <= float(RT_STRONG_TOL_SEC) else (
            "within_window" if abs_delta <= row["tolerance_sec"] else "mismatch"
        )
        )
        if query_method_key and row["method_key"]:
            method_status = (
                "same_method"
                if query_method_key == row["method_key"]
                else "method_mismatch"
            )
        else:
            method_status = "method_metadata_incomplete"
        method_comparable = method_status == "same_method"
        within_window = chromatographic_status in {
            "confirmed", "within_window",
        }
        identity_scope = (
            "exact_structure"
            if matched_by == "full_inchikey"
            else (
                "connectivity"
                if matched_by == "connectivity_inchikey14"
                else "name_formula"
            )
        )
        eligible_for_exact_identity = bool(
            within_window
            and method_comparable
            and identity_scope == "exact_structure"
        )
        if method_status != "same_method":
            status = method_status
        else:
            status = chromatographic_status
        ev = {
            "matched": bool(within_window and method_comparable),
            "status": status,
            "chromatographic_status": chromatographic_status,
            "method_status": method_status,
            "method_comparable": method_comparable,
            "query_method_id": _clean(query_method_id),
            "reference_method_id": row["method_id"],
            "matched_by": matched_by,
            "identity_scope": identity_scope,
            "eligible_for_exact_identity": eligible_for_exact_identity,
            "stereochemistry_capable": bool(
                row["stereochemistry_capable"]
            ),
            "separation_mode": row["separation_mode"],
            "query_rt_sec": round(q_rt_sec, 3),
            "reference_rt_sec": round(row["rt_sec"], 3),
            "rt_delta_sec": round(delta, 3),
            "abs_delta_sec": round(abs_delta, 3),
            "tolerance_sec": round(row["tolerance_sec"], 3),
            "standard_name": row["compound_name"],
            "source": row["source"],
            "note": row["note"],
        }
        priority = (
            int(eligible_for_exact_identity),
            int(method_comparable),
            int(identity_scope == "exact_structure"),
            -ev["abs_delta_sec"],
        )
        if (
            best is None
            or priority > best.get(
                "_priority", (-1, -1, -1, float("-inf"))
            )
        ):
            ev["_priority"] = priority
            best = ev
    if best is not None:
        best.pop("_priority", None)
    return best
