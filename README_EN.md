# MSsystem source supplement — 2026-11-02

This archive contains a reviewed snapshot of 855 selected source files from the current V16 working tree. It is a current dirty-working-tree capture, not a clean export of Git commit `e03b36cc5ebd287df9e273c24403b644999d2daa`.

## Authorization and licenses

The user explicitly confirmed on 2026-11-02: “确认有权，按MIT和CC BY 4.0准备投稿补充文件”. This extends the earlier MIT authorization to all self-developed source among the 855 files listed in `CURRENT_LICENSE_SCOPE.json`. The project-owned source is provided under the MIT License in root `LICENSE`.

Third-party source and bundles are not relicensed as project code. Cytoscape 3.28.1, Plotly.js 2.35.2, SmilesDrawer 2.0.3, and the MassCube-derived/fixture material retain their original MIT notices and copyright attributions; see `THIRD_PARTY_NOTICES.md` and embedded headers. `CURRENT_LICENSE_SCOPE.json` marks these paths separately.

The same user reply authorizes CC BY 4.0 preparation for author-owned data, while third-party data retain their original terms. Data licensing is supplied with the separate data supplement, not by this source archive.

No public upload, repository publication, DOI registration, or archive identifier was created by this preparation step. `CITATION.cff` therefore contains no invented DOI.

## Runtime and installation

The recorded runtime was CPython 3.13.5 on macOS arm64. `requirements.lock` records exact direct dependency versions, including SciPy 1.16.3. The older `environment.yml` requests Python 3.10 and is retained as historical project material; it is not the installation recommendation for this snapshot.

Suggested isolated setup, not tested as a fresh installation during packaging:

```bash
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.lock
```

Some packages are platform-sensitive. The lock lists direct dependencies rather than a complete transitive environment export.

## External inputs

This source archive is not a complete scientific-reproduction package. The production SQLite spectral library, mapping tables, large repair registries, model weights, raw LC-MS/MS inputs, and frozen evaluation artifacts are separate. See `DEPENDENCY_EXCLUSIONS.json` and the manuscript data supplement.

At minimum, configure compatible external shared data and a writable runtime directory:

```bash
export MSSYSTEM_SHARED_DATA_ROOT=/absolute/path/to/Mssystem_shared_data
export MSSYSTEM_RUNTIME_ROOT=/absolute/path/to/a/writable/runtime
export MSSYSTEM_LIBRARY_WRITE_ENABLED=0
```

Environment variables in `config/settings.py` override workstation-specific defaults. `config/runtime_manifest.json` contains the same sanitized payload as the reviewed pre-authorization package: absolute local path values were replaced by `${MSSYSTEM_LOCAL_PATH}/<basename>`, with modifications recorded in `SOURCE_PRIVACY_SCAN.json`.

## Launch and API example

After installing dependencies and obtaining compatible external inputs:

```bash
MSSYSTEM_PYTHON="$PWD/.venv/bin/python" bash run.sh --foreground --no-open --port 8516
curl -fsS http://127.0.0.1:8516/api/health
curl -fsS http://127.0.0.1:8516/openapi.json -o openapi.json
curl -fsS -F 'files=@example.mgf' http://127.0.0.1:8516/api/upload
```

The upload endpoint uses multipart field `files` and accepts `.mgf`, `.mzML`, `.mzXML`, or supported raw-report `.html` inputs. Startup, database loading, upload, and scientific analysis were not executed during source licensing preparation.

## Verification

The 855 source payloads were verified byte-for-byte against the pre-confirmation archive. No Python source changed during licensing. ZIP CRC, path safety, archive membership, and every per-file SHA-256 were checked after rebuilding. The earlier AST results remain recorded but were not rerun because the Python payload was unchanged. See `INSTALLATION_VERIFICATION.json` outside the ZIP.

