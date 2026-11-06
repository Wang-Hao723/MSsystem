#!/usr/bin/env python3
"""Fast startup storage/layout check (no hashing and no cleanup)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.settings import (  # noqa: E402
    DB_PATH,
    DB_PATH_SMALL,
    LIBRARY_WRITE_ENABLED,
    RUNTIME_ROOT,
    SHARED_DATA_ROOT,
)

_GUARD_SPEC = importlib.util.spec_from_file_location(
    "mssystem_storage_guard", ROOT / "utils" / "storage_guard.py"
)
if _GUARD_SPEC is None or _GUARD_SPEC.loader is None:
    raise RuntimeError("cannot load storage guard")
_GUARD = importlib.util.module_from_spec(_GUARD_SPEC)
_GUARD_SPEC.loader.exec_module(_GUARD)


def main() -> int:
    for database in (DB_PATH, DB_PATH_SMALL):
        if not database.is_file():
            raise FileNotFoundError(f"shared library database not found: {database}")
    status = _GUARD.require_storage_capacity(
        RUNTIME_ROOT, 2 * 1024 ** 3, "MSsystem startup"
    )
    print(
        f"Storage OK: shared={SHARED_DATA_ROOT}; runtime={RUNTIME_ROOT}; "
        f"writer={LIBRARY_WRITE_ENABLED}; free={status['free_gib']:.2f} GiB; "
        f"reserve={status['reserve_gib']:.2f} GiB"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
