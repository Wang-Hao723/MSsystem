# Test-suite execution record — 2026-11-07

Environment: macOS arm64, CPython 3.13.5 (repo `.venv`, `requirements.lock`),
`MSSYSTEM_SHARED_DATA_ROOT`/`MSSYSTEM_RUNTIME_ROOT` pointing at empty scratch
directories, `MSSYSTEM_LIBRARY_WRITE_ENABLED=0`.

Command:

```bash
python -m pytest tests/ -q --tb=no -p no:cacheprovider --continue-on-collection-errors
```

Result: **2328 passed, 96 failed, 8 skipped, 55 errors** (41 s).

Triage of the 96 failures and 55 errors (sampled with `--tb=line`):

- the majority import or load the `eval/` package (for example
  `tests/test_regression_smoke.py:708` → `ModuleNotFoundError: No module named
  'eval'`, `eval/eval_standards.py` not found) — `eval/` is one of the
  evaluation assets deliberately excluded from the reviewed source snapshot
  (see `DEPENDENCY_EXCLUSIONS.json`);
- the remainder fail in the artifact-integrity guard
  (`pipeline/derived_artifact_paths.py:65`, "Declared artifact implementation
  differs from loaded implementation"), which compares declared artifact
  bindings against files that are likewise part of the excluded evaluation
  holdings;
- no failure was traced to the scan-level matcher, the entity layer, or the
  shipped examples.

This is an executed summary of the archived snapshot on one machine, not a
continuous-integration guarantee. The two executable examples carry their own
verification (`example/verify_outputs.py`, `example_entity/verify_entity_outputs.py`),
both passing in the recorded environments.
