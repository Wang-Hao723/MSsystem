"""Check the identity, reference accession and peak pairs in each result."""
import csv
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main(output):
    actual = read_tsv(output / "ranking.tsv")
    expected = read_tsv(ROOT / "expected" / "ranking.tsv")
    actual_views = json.loads((output / "spectral_evidence.json").read_text())
    expected_views = json.loads((ROOT / "expected" / "spectral_evidence.json").read_text())
    if len(actual) != len(expected) or len(actual_views) != len(expected_views):
        raise SystemExit("The number of ranked entries or evidence records differs.")
    for row, baseline, view, saved in zip(actual, expected, actual_views, expected_views):
        for field in ("query", "rank", "candidate", "inchikey", "reference_accession", "matched_peak_pairs"):
            if row[field] != baseline[field]:
                raise SystemExit(f"{field} differs for {baseline['query']}: {row[field]}")
        if not math.isclose(float(row["score"]), float(baseline["score"]), rel_tol=0, abs_tol=1e-4):
            raise SystemExit(f"The score differs for {baseline['query']}: {row['score']}")
        evidence = view["spectral_evidence"]
        reference = saved["spectral_evidence"]
        if (view["query"] != row["query"] or str(view["rank"]) != row["rank"]
                or evidence.get("status") != "captured"
                or not evidence.get("replay_matches_comparison")
                or evidence["reference"]["source_accession"] != row["reference_accession"]
                or evidence["matched_pairs"] != reference["matched_pairs"]):
            raise SystemExit(f"Peak evidence does not reproduce for {row['query']}")
    print(f"PASS: {len(actual)}/{len(expected)} rankings and all saved peak pairs reproduced.")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
