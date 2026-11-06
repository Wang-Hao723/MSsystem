#!/usr/bin/env python3
"""Summarize MSsystem v5 evidence tiers and suggest RT/curation priorities."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def _rows(data):
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        for k in ("results", "data", "spectra", "items"):
            if isinstance(data.get(k), list):
                return data[k]
    return []


def _match(row):
    if not isinstance(row, dict):
        return None
    if isinstance(row.get("match"), dict):
        return row["match"]
    ms = row.get("matches")
    if isinstance(ms, list) and ms:
        return ms[0]
    return None


def _as_float(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _evidence_tier(m: dict) -> str:
    existing = m.get("evidence_tier")
    if existing:
        return existing
    conf = m.get("confidence") or "Unknown"
    cconf = m.get("confidence_class") or "Unknown"
    iso = bool(m.get("isomer_ambiguous") or int(m.get("isomer_n") or 0) > 0)
    cos = _as_float(m.get("cosine"))
    final = _as_float(m.get("final_score"))
    ppm = abs(_as_float(m.get("mass_ppm"), 999.0))
    cls_ag = _as_float(m.get("class_agreement"), 1.0)
    if conf == "High":
        return "High"
    if cconf == "High" and iso:
        return "Class High"
    if conf == "Medium" and (not iso) and cos >= 0.85 and final >= 0.55 and ppm <= 50 and cls_ag >= 0.5:
        return "Strong Candidate"
    if cconf == "High":
        return "Class High"
    return "Putative"


def summarize(path: Path, output: Path | None = None):
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = _rows(data)
    matched = [(r, _match(r)) for r in rows if _match(r)]
    conf = Counter((m.get("confidence") or "Unknown") for _, m in matched)
    cconf = Counter((m.get("confidence_class") or "Unknown") for _, m in matched)
    etier = Counter(_evidence_tier(m) for _, m in matched)
    htier = Counter((m.get("honesty_tier") or "Unknown") for _, m in matched)

    print(f"File: {path}")
    print(f"Spectra: {len(rows)}")
    print(f"Matched spectra: {len(matched)}")
    print(f"Compound confidence: {dict(conf)}")
    print(f"Class confidence: {dict(cconf)}")
    print(f"Evidence tiers: {dict(etier)}")
    print(f"Honesty tiers: {dict(htier)}")

    candidates = []
    for r, m in matched:
        ev = _evidence_tier(m)
        need_rt = ev in {"Strong Candidate", "Class High"} or m.get("isomer_ambiguous")
        if not need_rt:
            continue
        candidates.append({
            "spectrum_id": r.get("spectrum_id", ""),
            "query_rt": r.get("retention_time", ""),
            "compound_name": m.get("compound_name", ""),
            "formula": m.get("formula", ""),
            "inchikey": m.get("inchikey", ""),
            "classification": m.get("classification", ""),
            "cosine": m.get("cosine", ""),
            "final_score": m.get("final_score", ""),
            "mass_ppm": m.get("mass_ppm", ""),
            "compound_confidence": m.get("confidence", ""),
            "class_confidence": m.get("confidence_class", ""),
            "evidence_tier": ev,
            "isomer_n": m.get("isomer_n", 0),
            "priority_score": round(
                _as_float(m.get("cosine")) * 2
                + _as_float(m.get("final_score"))
                + (1 if ev == "Strong Candidate" else 0)
                + min(int(m.get("isomer_n") or 0), 3) * 0.2,
                4,
            ),
        })
    candidates.sort(key=lambda x: -x["priority_score"])
    print(f"RT/curation priority candidates: {len(candidates)}")
    for c in candidates[:20]:
        print(
            f"  {c['evidence_tier']:<16} cos={c['cosine']} final={c['final_score']} "
            f"iso={c['isomer_n']} rt={c['query_rt']} name={c['compound_name'][:60]}"
        )

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        fields = list(candidates[0].keys()) if candidates else [
            "spectrum_id", "query_rt", "compound_name", "formula", "inchikey",
            "classification", "cosine", "final_score", "mass_ppm",
            "compound_confidence", "class_confidence", "evidence_tier",
            "isomer_n", "priority_score",
        ]
        with output.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(candidates)
        print(f"Wrote: {output}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_json", type=Path)
    ap.add_argument("--output", "-o", type=Path)
    args = ap.parse_args()
    summarize(args.results_json, args.output)


if __name__ == "__main__":
    main()
