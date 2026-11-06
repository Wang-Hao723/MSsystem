"""Run the competing-candidate entity example: multi-energy scans -> entity ranking.

For each query scan, the shipped scan-level matcher searches the example
library; the per-scan results then drive the entity-level flagship layer
(matching.flagship_ranker), which groups scans by collision energy, aggregates
candidate support, ranks connectivities and records the reporting decision
with its margins. Evidence records keep the scan -> reference -> matched-peak
links behind every rank.
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path


def ce_number_of(query):
    match = re.match(r"(\d+(?:\.\d+)?)", str(query.get("collisionenergy", "")))
    return float(match.group(1)) if match else 0.0


def parse_mgf(path):
    block, out = {}, []
    for line in Path(path).read_text().splitlines():
        if line == "END IONS":
            if block:
                out.append(block)
            block = {}
        elif line == "BEGIN IONS":
            block = {"peaks": []}
        elif not block or not line.strip():
            continue
        elif "=" in line and not line[:1].isdigit():
            key, value = line.split("=", 1)
            block[key.lower()] = value
        else:
            parts = line.split()
            if len(parts) >= 2:
                block["peaks"].append([float(parts[0]), float(parts[1])])
    return out


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--source", required=True, type=Path)
    cli.add_argument("--db", required=True, type=Path)
    cli.add_argument("--out", required=True, type=Path)
    cli.add_argument("--task-id", default="entity-demo-flavonoid")
    args = cli.parse_args()
    sys.path.insert(0, str(args.source.resolve()))
    from matching.search import load_library, search
    from matching.flagship_ranker import build_flagship_result

    demo_dir = Path(__file__).resolve().parent
    records, index = load_library(args.db)
    if not records:
        raise RuntimeError("The example library contains no usable reference spectra")

    rows, evidence = [], []
    queries = sorted(parse_mgf(demo_dir / "entity_queries.mgf"), key=ce_number_of)
    for position, query in enumerate(queries):
        peaks = query["peaks"]
        precursor = float(query["pepmass"].split()[0])
        ce = ce_number_of(query)
        result = search(
            peaks, precursor, records, index,
            query_polarity="POS", enforce_polarity=True,
            query_adduct="[M+H]+", query_adduct_reliable=True, query_charge=1,
            query_collision_energy=ce, top_n=10,
        )
        if result is None:
            raise RuntimeError(f"No eligible candidate for {query.get('title')}")
        scan_evidence = []
        for candidate in [result, *(result.get("top_candidates") or [])][1:] or [result]:
            view = candidate.get("spectral_evidence") or {}
            accession = candidate.get("lib_source_accession") or (view.get("reference") or {}).get("source_accession")
            scan_evidence.append({
                "candidate": candidate.get("compound_name"),
                "inchikey": candidate.get("inchikey"),
                "reference_accession": accession,
                "final_score": candidate.get("final_score"),
                "matched_peak_pairs": view.get("matched_pairs") or [],
            })
        evidence.append({"scan": str(query.get("title")), "collision_energy_eV": ce,
                         "candidates": scan_evidence})
        rows.append({
            "spectrum_id": str(query.get("title")),
            "precursor_feature_id": "PF-FLAVONOID-001",
            "precursor_mz": precursor,
            "retention_time": position + 1.0,
            "collision_energy": ce,
            "collision_energy_unit": "eV",
            "ion_mode": "POS",
            "tic": float(sum(intensity for _, intensity in peaks)),
            "raw_peaks": peaks,
            "query_peaks": peaks,
            "result_index": position,
            "match": result,
        })

    payload = build_flagship_result(args.task_id, rows)
    args.out.mkdir(parents=True, exist_ok=True)

    ranking_rows = []
    for entity in payload.get("entities") or []:
        for candidate in entity.get("candidates") or []:
            inner = candidate.get("candidate") or {}
            identity = str(candidate.get("identity_key") or "")
            ik14 = identity.split(":", 1)[-1][:14]
            ranking_rows.append({
                "entity_id": entity.get("entity_id"),
                "rank": candidate.get("entity_rank"),
                "candidate": inner.get("compound_name"),
                "identity_scope": candidate.get("identity_scope"),
                "connectivity_inchikey14": ik14,
                "formula": inner.get("formula"),
                "entity_score": candidate.get("entity_score"),
                "ce_bins_supported": candidate.get("collision_energy_bins_supported"),
                "ce_bins_total": candidate.get("collision_energy_bins_total"),
                "decision_status": (entity.get("decision") or {}).get("status"),
                "decision_scope": (entity.get("decision") or {}).get("scope"),
                "winner_margin": (entity.get("decision") or {}).get("winner_margin"),
                "member_scans": ";".join(entity.get("member_spectrum_ids") or []),
                "entity_ce_bins": ";".join(str(b) for b in (entity.get("collision_energy_bins") or [])),
            })
    with (args.out / "entity_ranking.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ranking_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(ranking_rows)

    (args.out / "entity_decision.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    (args.out / "scan_evidence.json").write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")

    decision = next((e.get("decision") for e in payload.get("entities") or []), None)
    print(f"entity ranking rows: {len(ranking_rows)}; decision: {json.dumps(decision, ensure_ascii=False)[:200] if decision else None}")


if __name__ == "__main__":
    main()
