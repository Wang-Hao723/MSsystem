#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 3 )); then
  printf 'Usage: bash %s /absolute/path/to/venv/bin/python /path/to/MSsystem_source_20261102 [writable-work-parent]\n' "$0" >&2
  exit 2
fi
demo_dir=$(cd -- "$(dirname -- "$0")" && pwd -P)
python_bin=$1
source_dir=$(cd -- "$2" && pwd -P)
work_parent=${3:-${TMPDIR:-/tmp}}
work=$(mktemp -d "$work_parent/mssystem-example.XXXXXX")
printf 'Run directory: %s\n' "$work"
mkdir "$work/library" "$work/runtime" "$work/result"
cd "$source_dir"

env -i PATH="$PATH" MSSYSTEM_SHARED_DATA_ROOT="$work/unused" \
  MSSYSTEM_RUNTIME_ROOT="$work/runtime" MSSYSTEM_LIBRARY_WRITE_ENABLED=0 \
  "$python_bin" build_from_mgf_msp.py --mgf-dir "$demo_dir/reference" \
  --output "$work/library" 2>&1 | tee "$work/build.log"

env -i PATH="$PATH" MSSYSTEM_SHARED_DATA_ROOT="$work/library" \
  MSSYSTEM_RUNTIME_ROOT="$work/runtime" MSSYSTEM_LIBRARY_WRITE_ENABLED=0 \
  "$python_bin" "$demo_dir/run_analysis.py" --source "$source_dir" \
  --db "$work/library/local_library.db" --out "$work/result" 2>&1 | tee "$work/search.log"

"$python_bin" "$demo_dir/verify_outputs.py" "$work/result"
cat "$work/result/ranking.tsv"
