#!/usr/bin/env python3
"""Cross-machine installation audit for MSsystem V10."""

from __future__ import annotations

import argparse
import importlib
import json
import platform
import sys
from importlib import metadata
from pathlib import Path
from typing import Callable, Iterable

from packaging.version import InvalidVersion, Version


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REQUIRED_MODULES = (
    "numpy",
    "pandas",
    "scipy",
    "sklearn",
    "matchms",
    "pyteomics",
    "rdkit",
    "fastapi",
    "uvicorn",
    "sqlalchemy",
    "ijson",
)
OPTIONAL_MODULES = (
    "pyopenms",
    "streamlit",
    "openpyxl",
    "reportlab",
    "weasyprint",
)
REQUIRED_PROJECT_FILES = (
    "config/settings.py",
    "config/module_governance.json",
    "web/server.py",
    "matching/search.py",
    "run.sh",
    "requirements.lock",
)
IMPORT_DISTRIBUTIONS = {
    "sklearn": "scikit-learn",
    "pyteomics": "pyteomics",
    "rdkit": "rdkit",
    "pyopenms": "pyopenms",
}


def _module_version(module: object) -> str:
    return str(getattr(module, "__version__", "installed"))


def _read_lock_versions(path: Path) -> dict[str, str]:
    versions: dict[str, str] = {}
    if not path.is_file():
        return versions
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "==" not in line:
            continue
        name, version = line.split("==", 1)
        versions[name.strip().lower().replace("_", "-")] = version.strip()
    return versions


def _same_version(left: str, right: str) -> bool:
    try:
        return Version(left) == Version(right)
    except InvalidVersion:
        return left == right


def check_installation(
    root: Path,
    *,
    required_modules: Iterable[str] = REQUIRED_MODULES,
    optional_modules: Iterable[str] = OPTIONAL_MODULES,
    importer: Callable[[str], object] = importlib.import_module,
    version_info: tuple[int, int, int] | None = None,
    include_local_paths: bool = False,
) -> dict:
    """Return a serializable report and never modify the source tree."""
    root = root.resolve()
    current_version = version_info or (
        sys.version_info.major,
        sys.version_info.minor,
        sys.version_info.micro,
    )
    required: dict[str, dict] = {}
    optional: dict[str, dict] = {}
    lock_versions = _read_lock_versions(root / "requirements.lock")
    version_drift: list[dict[str, str]] = []
    for name, target in (
        *[(name, required) for name in required_modules],
        *[(name, optional) for name in optional_modules],
    ):
        try:
            module = importer(name)
        except Exception as exc:
            target[name] = {
                "available": False,
                "error": f"{type(exc).__name__}: {exc}",
            }
        else:
            runtime_version = _module_version(module)
            distribution_name = IMPORT_DISTRIBUTIONS.get(name, name)
            normalized_distribution = distribution_name.lower().replace("_", "-")
            locked_version = lock_versions.get(normalized_distribution)
            try:
                metadata_version = metadata.version(distribution_name)
            except metadata.PackageNotFoundError:
                metadata_version = "unavailable"
            runtime_matches_lock = (
                None
                if not locked_version or runtime_version == "installed"
                else _same_version(runtime_version, locked_version)
            )
            metadata_matches_lock = (
                None
                if not locked_version or metadata_version == "unavailable"
                else _same_version(metadata_version, locked_version)
            )
            runtime_matches_metadata = (
                None
                if runtime_version == "installed" or metadata_version == "unavailable"
                else _same_version(runtime_version, metadata_version)
            )
            target[name] = {
                "available": True,
                "runtime_version": runtime_version,
                "package_metadata_version": metadata_version,
                "locked_version": locked_version,
                "runtime_matches_lock": runtime_matches_lock,
                "metadata_matches_lock": metadata_matches_lock,
                "runtime_matches_package_metadata": runtime_matches_metadata,
            }
            if runtime_matches_lock is False or metadata_matches_lock is False or runtime_matches_metadata is False:
                version_drift.append(
                    {
                        "module": name,
                        "runtime_version": runtime_version,
                        "package_metadata_version": metadata_version,
                        "locked_version": str(locked_version or "unrecorded"),
                    }
                )
    files = {
        relative: (root / relative).is_file()
        for relative in REQUIRED_PROJECT_FILES
    }
    python_supported = current_version >= (3, 10, 0)
    required_ok = all(item["available"] for item in required.values())
    files_ok = all(files.values())
    core_ready = python_supported and required_ok and files_ok
    lock_consistent = not version_drift
    status = (
        "ready"
        if core_ready and lock_consistent
        else "ready_with_version_drift"
        if core_ready
        else "not_ready"
    )
    return {
        "schema_version": "1.0.0",
        "status": status,
        "runtime_ready": core_ready,
        "strict_reproducibility_ready": core_ready and lock_consistent,
        "version_drift": version_drift,
        "python": {
            "version": ".".join(str(value) for value in current_version),
            "executable": (
                sys.executable if include_local_paths else Path(sys.executable).name
            ),
            "supported": python_supported,
            "minimum": "3.10.0",
        },
        "platform": {
            "system": platform.system(),
            "machine": platform.machine(),
            "release": platform.release(),
        },
        "required_modules": required,
        "optional_modules": optional,
        "required_project_files": files,
        "optional_modules_block_core_startup": False,
        "root": str(root) if include_local_paths else ".",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional JSON evidence path; no file is written when omitted",
    )
    parser.add_argument(
        "--strict-optional",
        action="store_true",
        help="Also fail when an optional integration is unavailable",
    )
    parser.add_argument(
        "--include-local-paths",
        action="store_true",
        help="Include absolute interpreter and project paths in diagnostic output",
    )
    args = parser.parse_args()
    report = check_installation(
        args.root,
        include_local_paths=args.include_local_paths,
    )
    if args.strict_optional and not all(
        item["available"] for item in report["optional_modules"].values()
    ):
        report["status"] = "not_ready_optional"
    rendered = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if report["status"] == "ready" else 2


if __name__ == "__main__":
    raise SystemExit(main())
