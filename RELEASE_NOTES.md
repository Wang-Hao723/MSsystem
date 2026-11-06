# Public release notes — v1.0.0 (2026-11-02 source snapshot)

This repository publishes the reviewed 2026-11-02 MSsystem source snapshot that is also delivered as Additional File 1 of the manuscript. The source payload is byte-identical to that archive; the following public-release adjustments were made at repository level and are recorded here.

| # | Change | Reason |
| --- | --- | --- |
| 1 | `CITATION.cff`: title changed from "MSsystem V16 source supplement: …" to "MSsystem: …"; `version` changed from `16-source-supplement-2026-11-02` to `1.0.0` | Public citation metadata uses the release name "MSsystem" and semantic version `1.0.0`; internal working-tree generation labels are not part of the public versioning. This is the only modification to a snapshot file, so the per-file SHA-256 for `CITATION.cff` in `SOURCE_MANIFEST.json` refers to the archived version. |
| 2 | Added root `README.md` and this `RELEASE_NOTES.md` | Repository landing page and change record; not part of the snapshot payload. |
| 3 | Added `example/` containing Additional File 3 (executable example) minus the nested copy of the Additional File 1 source archive | The source already lives at the repository root; `example/run_demo.sh` receives the repository root as its source-directory argument. All other example files are unmodified, and `example/SHA256SUMS.txt` entries continue to match the retained files. |

No Python source, configuration, test, documentation or data file of the snapshot payload was modified. `SOURCE_MANIFEST.json`, `SOURCE_PRIVACY_SCAN.json` and `CURRENT_LICENSE_SCOPE.json` describe the archived snapshot exactly as delivered with the manuscript.

`example_entity/` (added 2026-11-07, repository main branch) contains the competing-candidate entity example: three multi-energy public scans of isoquercitrin are aggregated into one precursor entity and ranked against experimental references for hyperoside (same 2D connectivity, different stereochemistry) and spiraeoside (different 2D connectivity). The example exercises multi-scan entity aggregation, connectivity-level competition, the reporting decision with margins and per-scan evidence export, using the shipped matcher and entity layer. It is additive to the reviewed snapshot and is described in RELEASE_NOTES rather than in the manuscript's Additional file 3.

The manuscript additional files are attached byte-identically to the `v1.0.0` GitHub release (Additional File 1: source archive; Additional File 2: evaluation results and recomputation materials; Additional File 3: executable example bundle).
