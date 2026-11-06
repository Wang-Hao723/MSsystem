# Recomputation chain for the manuscript figures and tables

This document maps every main-text figure and table of the manuscript to the
materials that reproduce it: where the saved records live (Additional file 2 of
the manuscript, attached byte-identically to the v1.0.0 GitHub release), which
command replays them, which inputs are publicly re-obtainable, and where the
current boundary of reproduction lies. It implements the "figure → query
inputs → reference inputs → configuration → command → expected output"
inventory discussed for the submission, and states gaps as gaps.

`data_sources/PUBLIC_INPUTS.json` (frozen URLs and hashes of the public inputs)
and `data_sources/ACCESSION_REGISTER.tsv` (the full data-supplement inventory
with per-dataset access states) are copied verbatim from Additional file 2.

| Manuscript item | Materials | How to reproduce | Public inputs | Current boundary |
| --- | --- | --- | --- | --- |
| Fig. 3 + Table 2 — complete workflow, 64 standard records (50/64 vs 56/64, intervals) | AF2 `formal32/` (per-record paired outputs, resampling) | `python formal32/replay_formal32.py` (needs NumPy) | none (laboratory acquisitions) | Statistical replay from saved per-record outputs. Re-scoring from raw scans needs the four original vendor acquisitions and the production-library snapshot (not publicly archived). |
| Fig. 4 — MSMLS common-reference equal-output curves (155 vs 144 at 208 outputs, −5.29 pp interval) | AF2 `display_data/statistics_ALL_COMMON_OUTPUT_COUNTS.csv` (357-row curve), `statistics_EQUAL_OUTPUT_RESULT.json` | curve and interval recomputation from the saved matched counts | MSMLS spectra are commercial; 322 queries and the 216,128-record reference snapshot are not in the score-only replay | Score-only replay. Native-score regeneration from spectra requires the MSMLS inputs and frozen program configurations. |
| Fig. 5 + Table 4 — member retention (270/288/289 of 494) and four-arm reference processing (535/535/546/548 of 1,066) | AF2 `retention/` (S3 predictions, replay scripts) and `four_branch_replay.zip` (validated portable four-branch package) | `python retention/replay_from_bundle.py --output-dir <empty>`; four-branch per its in-package instructions | none (laboratory mixed samples) | Statistical replay from saved predictions. The original four-branch per-unit archive files and part of the member-level inputs are not recovered; upstream re-scoring is not possible from the released materials. |
| Table 3 (rows 1–2) — MSMLS comparisons | see Fig. 4 | see Fig. 4 | see Fig. 4 | see Fig. 4 |
| Table 3 (row 3) — MassBank common-reference comparison (30/66 vs 31/66) | `isomer_controls/` package and input snapshot (registered in `ACCESSION_REGISTER.tsv`) | `build_massbank202603_query_windows.py` rebuilds query windows from the pinned MassBank release | MassBank-data 2026.03, commit `705afb7b…`, URL and hash in `PUBLIC_INPUTS.json` | Rebuild path exists for inputs; frozen prediction replay per package instructions. |
| Table 3 (row 4) — SIRIUS common-query comparison (77/132 vs 60/132) | AF2 saved predictions (registered) | replay of saved predictions | the 107 anonymised queries are laboratory spectra | Program comparison under fixed queries; candidate resources differ per program. |
| Fig. 6 — blinded rounds recovery (82/88; 74/86/93 of 98) | AF2 `display_data/statistics_blind_mixture_recovery.csv` | statistics recomputation from the saved counts | none (blinded laboratory mixtures) | Statistics provided; raw blind injections and their input manifests are not yet publicly archived. |
| Fig. 7 — replicate repeatability, injection amount, Orbitrap replay | analytical-robustness summary tables (registered in `ACCESSION_REGISTER.tsv`) | saved tables; R sources in the figure delivery | none | Result tables retained; raw replicate/Orbitrap acquisitions pending public deposition. |
| Fig. 8 — phillyrin case | `figures/phillyrin_member_peaks.csv` + `Phillyrin_provenance.json` (registered) | R source from the figure delivery replots from the CSV | none (laboratory scans; measured peak lists retained) | The reference peak array of the matched reference is not recovered; the mirror comparison is presented through stored scores, and no spectrum is invented. |
| Fig. 9 — flavonoid interpretation | knowledge illustration (mass positions from rules) | n/a (no dataset) | n/a | — |
| Fig. 10 + reporting calibration (69→75/76, 26→30/39, 94/21/23→88/5/45) | AF2 `display_data/statistics_candidate_review_reporting.csv`; figure delivery R sources | statistics from the saved counts | none | Same-cohort calibration; thresholds and outcomes stored together. |
| Entity example (repository) — `example_entity/` | this repository, main branch | `bash example_entity/run_entity_demo.sh …` | MassBank 2026.03 records bundled with attribution (CC BY / CC BY-NC) | Fully executable from public inputs; deterministic outputs in `example_entity/expected/`. |
| Scan-level example — `example/` (Additional file 3) | this repository + release asset | `bash example/run_demo.sh …` | MassBank records bundled (CC BY-SA) | Fully executable from public inputs. |

## Status summary

- **Statistically replayable from saved records now**: Fig. 3/Table 2, Fig. 4/Table 3 (rows 1–2), Fig. 5/Table 4, Fig. 6, Fig. 10 (all materials in Additional file 2; scripts included).
- **Fully executable from public inputs now**: the two repository examples (`example/`, `example_entity/`).
- **Input-level regeneration**: available for the MassBank common-reference experiment via the pinned public release; not available for laboratory/vendor and MSMLS inputs, whose raw files remain undeposited.
- **Known unrecovered items** (recorded in `ACCESSION_REGISTER.tsv`, not silently dropped): four-branch original per-unit files, part of the member-level inputs, the phillyrin reference peak array, blind-round and Orbitrap raw acquisitions, and the full operational spectral library.

Adding a DOI to these materials does not by itself close the input-level gaps;
closing them requires either public deposition of the corresponding raw
inputs or replacing the affected claims with replayable evaluations, as stated
in the manuscript's availability section.
