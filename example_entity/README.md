# MSsystem competing-candidate entity example

This example exercises the centre of the MSsystem workflow: several scans of one
precursor at multiple collision energies, experimental references for competing
structures, entity-level candidate ranking, a reporting decision with margins,
and the scan → reference → matched-peak evidence behind every rank. It
complements `../example/`, which covers installation and scan-level retrieval
for a single compound.

## Composition

| Role | Records | Structure | Source, licence |
| --- | --- | --- | --- |
| Query scans (3) | MSBNK-BGC_Munich-RP017101/02/03, positive mode, 10/20/40 eV | Isoquercitrin (truth connectivity `OVSQVDMCBVZWGM`) | BGC Munich, MassBank 2026.03, CC BY |
| References (4) | MSBNK-Fiocruz-FIO00160–163, 10/20/30/40 eV | Isoquercitrin | Fiocruz, CC BY-NC |
| References (3) | MSBNK-RIKEN_ReSpect-PS046501–03, 10/20/30 eV | Hyperoside — same 2D connectivity, different stereochemistry | RIKEN ReSpect, CC BY-NC |
| References (3) | MSBNK-RIKEN_ReSpect-PS086901–03, 10/20/30 eV | Spiraeoside — different 2D connectivity (`OIUBYZLTFSLSBY`) | RIKEN ReSpect, CC BY-NC |

All records come from the pinned MassBank-data release 2026.03 (commit
`705afb7bccc3b2c42410a744eef73674716a60ef`). `prepare_entity_inputs.py`
converts the bundled records into query and reference MGF files and writes
per-record provenance with licence, authors, record URL and SHA-256. Where a
record's `CH$EXACT_MASS` deviates from the formula monoisotopic mass by more
than 50 ppm, the formula mass is used and the substitution is disclosed in
`provenance.json`.

## Run

From the repository root (see `../README_EN.md` for environment setup):

```bash
uv venv --python python3.13 .venv
uv --no-cache pip install --python .venv/bin/python -r requirements.lock
bash example_entity/run_entity_demo.sh "$PWD/.venv/bin/python" "$PWD"
```

The script builds a local SQLite library from the reference MGF, searches each
query scan with the shipped matcher, feeds the per-scan results to the
entity-level flagship layer (`matching.flagship_ranker.build_flagship_result`),
and writes `entity_ranking.tsv`, `entity_decision.json` and `scan_evidence.json`
to a fresh work directory. `verify_entity_outputs.py` then checks the result
(15 checks: entity aggregation, ≥2 collision-energy bins, truth connectivity at
rank 1, retention of the distinct competitor connectivity, per-bin support,
reporting decision, and per-scan evidence).

## Expected result

`expected/` holds the outputs of the validated run (macOS arm64, CPython
3.13.5, `verification/requirements-macos-arm64-py313-tested.txt`; the run is
deterministic — two independent runs produced byte-identical ranking and
evidence files):

- rank 1 — connectivity `OVSQVDMCBVZWGM` (quercetin O-hexoside;
  isoquercitrin and hyperoside collapse to one connectivity identity),
  entity score 0.893473, supported at 3/3 collision-energy bins;
- rank 2 — connectivity `OIUBYZLTFSLSBY` (spiraeoside), entity score 0.709793;
- decision — status `selected`, scope `2d_connectivity_candidate`,
  winner margin 0.18368 (minimum score 0.50, minimum margin 0.04).

The example does not separate the same-connectivity sugar stereoisomers
(isoquercitrin vs hyperoside): at the connectivity level they are one
candidate, and the sugar-identity distinction remains with standards and
same-condition comparison, as discussed in the manuscript. Distinguishing the
two 2D connectivities (quercetin O-hexoside vs quercetin 4′-O-glucoside) is
what the entity ranking demonstrates here.
