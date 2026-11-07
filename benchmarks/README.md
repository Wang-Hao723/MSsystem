# Search-scaling benchmark

`search_scaling.py` measures library-build time, per-scan retrieval latency and
peak memory as the reference library grows. Libraries of 10/100/1,000/5,000
positive-mode MS2 records are drawn as a deterministic sample from the pinned
MassBank-data 2026.03 release, with the ten query-compound records guaranteed
present in every library; ten fixed query scans are then searched with the
shipped matcher. Each size runs in a fresh process.

Recorded environment: Apple M4 Pro, 48 GB RAM, macOS arm64, CPython 3.13.5,
repo `.venv` from `requirements.lock` (2026-11-07).

| Library refs (input) | Distinct compounds loaded | Build time (s) | Median scan latency (ms) | 95th pct (ms) | Queries matched | Peak RSS (MB) |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | 3 | 2.1 | 22.5 | 23.0 | 10/10 | 436 |
| 100 | 84 | 2.2 | 22.2 | 22.5 | 10/10 | 443 |
| 1,000 | 866 | 2.9 | 22.7 | 23.2 | 10/10 | 447 |
| 5,000 | 4,355 | 5.6 | 42.4 | 43.4 | 10/10 | 469 |

Notes: repeated measurements of the same structure consolidate into one
compound entry with its spectra retained, so loaded-compound counts are below
input-record counts. The first query in a process adds ≈1 s of initialisation
(reflected in `search_ms_max`); steady-state scans are the median/p95 columns.
Peak RSS is dominated by the interpreter and numerical dependencies and grows
by ~8 % over a 500-fold library increase.

## Reproduce

```bash
# place the extracted MassBank-data 2026.03 tree (commit 705afb7…) here:
ln -s /path/to/MassBank-MassBank-data-705afb7 benchmarks/massbank_source
for n in 10 100 1000 5000; do python benchmarks/search_scaling.py --size $n; done
```

`massbank_source/` (a symlink to the external extracted tree) and the
generated `massbank_sample/` are intentionally not part of the repository; the
records are re-obtainable from the pinned public release.
