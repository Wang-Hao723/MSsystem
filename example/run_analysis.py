"""Export MSsystem's scan-level ranking and the engine's spectral evidence."""
import argparse
import csv
import json
from pathlib import Path
import sys

from matchms.importing import load_from_mgf


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--source", required=True, type=Path)
    cli.add_argument("--db", required=True, type=Path)
    cli.add_argument("--out", required=True, type=Path)
    args = cli.parse_args()
    sys.path.insert(0, str(args.source.resolve()))
    from matching.search import load_library, search

    query_file = Path(__file__).resolve().parent / "example.mgf"
    records, index = load_library(args.db)
    if not records:
        raise RuntimeError("The example library contains no usable reference spectra")
    rows, evidence_records = [], []
    for spectrum in load_from_mgf(str(query_file)):
        query_id = spectrum.get("title")
        result = search(
            list(zip(spectrum.peaks.mz, spectrum.peaks.intensities)),
            float(spectrum.get("precursor_mz")), records, index,
            query_polarity="NEG", enforce_polarity=True,
            query_adduct="[M-H]-", query_adduct_reliable=True, query_charge=-1,
            query_collision_energy=spectrum.get("collisionenergy"),
        )
        if result is None:
            raise RuntimeError(f"No eligible candidate for {query_id}")
        for rank, candidate in enumerate(result.get("top_candidates") or [result], 1):
            view = candidate.get("spectral_evidence") or {}
            accession = candidate.get("lib_source_accession") or (view.get("reference") or {}).get("source_accession")
            pairs = view.get("matched_pairs") or []
            rows.append({
                "query": query_id, "rank": rank, "candidate": candidate["compound_name"],
                "inchikey": candidate["inchikey"], "reference_accession": accession,
                "score": candidate["final_score"], "matched_peak_pairs": len(pairs),
            })
            evidence_records.append({"query": query_id, "rank": rank, "spectral_evidence": view})
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "ranking.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.out / "spectral_evidence.json").write_text(
        json.dumps(evidence_records, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    print(f"Wrote {len(rows)} ranked entries from {len(records)} reference spectra.")


if __name__ == "__main__":
    main()
