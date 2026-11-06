"""Verify the competing-candidate entity example outputs.

Checks that the entity pipeline consumed scans at multiple collision energies,
that the true connectivity ranks first, that the structurally distinct
competitor connectivity is retained among the ranked candidates, and that
every scan's evidence carries reference accessions and matched peak pairs.
"""
import csv
import json
import sys
from pathlib import Path

TRUTH_CONNECTIVITY = "OVSQVDMCBVZWGM"      # isoquercitrin / hyperoside (quercetin O-hexoside)
COMPETITOR_CONNECTIVITY = "OIUBYZLTFSLSBY"  # spiraeoside (quercetin 4'-O-glucoside)


def main():
    result_dir = Path(sys.argv[1])
    ranking = list(csv.DictReader((result_dir / "entity_ranking.tsv").open(encoding="utf-8"), delimiter="\t"))
    decision = json.loads((result_dir / "entity_decision.json").read_text(encoding="utf-8"))
    evidence = json.loads((result_dir / "scan_evidence.json").read_text(encoding="utf-8"))

    failures = []

    def check(condition, message):
        print(("PASS " if condition else "FAIL ") + message)
        if not condition:
            failures.append(message)

    entities = decision.get("entities") or []
    check(len(entities) == 1, "one precursor entity aggregates the query scans")
    entity = entities[0] if entities else {}
    check(len(entity.get("member_spectrum_ids") or []) >= 3,
          f"at least 3 member scans recorded (got {len(entity.get('member_spectrum_ids') or [])})")
    energy_bins = [str(b) for b in (entity.get("collision_energy_bins") or [])]
    check(len(set(energy_bins)) >= 2,
          f"scans span at least 2 collision-energy bins (got {sorted(set(energy_bins))})")

    connectivities = {row["connectivity_inchikey14"]: row for row in ranking}
    first = next((row for row in ranking if str(row["rank"]) == "1"), None)
    check(first is not None and first["connectivity_inchikey14"] == TRUTH_CONNECTIVITY,
          f"rank-1 connectivity is the truth {TRUTH_CONNECTIVITY} (got {first and first['connectivity_inchikey14']})")
    check(COMPETITOR_CONNECTIVITY in connectivities,
          f"distinct competitor connectivity {COMPETITOR_CONNECTIVITY} retained among ranked candidates")
    winner_inner = ((entity.get("winner") or {}).get("candidate")) if isinstance(entity.get("winner"), dict) else None
    check(bool(winner_inner or entity.get("winner")),
          "winning candidate recorded with the decision")
    first_row_support = (int(first["ce_bins_supported"]), int(first["ce_bins_total"])) if first and first.get("ce_bins_supported") else None
    check(first_row_support is not None and first_row_support[0] == first_row_support[1] and first_row_support[1] >= 2,
          f"rank-1 connectivity supported at every collision-energy bin (got {first_row_support})")

    status = (entity.get("decision") or {}).get("status")
    check(bool(status), f"reporting decision recorded (status={status})")

    for scan in evidence:
        with_pairs = [c for c in scan["candidates"] if c["matched_peak_pairs"]]
        with_accessions = [c for c in scan["candidates"] if c["reference_accession"]]
        check(bool(with_pairs), f"scan {scan['scan']} carries matched peak pairs")
        check(bool(with_accessions), f"scan {scan['scan']} carries reference accessions")

    if failures:
        print(f"\n{len(failures)} check(s) failed")
        sys.exit(1)
    print("\nall entity-example checks passed")


if __name__ == "__main__":
    main()
