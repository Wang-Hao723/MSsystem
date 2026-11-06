from pathlib import Path
from types import SimpleNamespace

from scripts.check_installation import check_installation


ROOT = Path(__file__).resolve().parents[1]


def _importer(name):
    if name == "missing_optional":
        raise ModuleNotFoundError(name)
    return SimpleNamespace(__version__="1.2.3")


def test_optional_integration_does_not_block_core_readiness():
    report = check_installation(
        ROOT,
        required_modules=("available",),
        optional_modules=("missing_optional",),
        importer=_importer,
        version_info=(3, 10, 0),
    )

    assert report["status"] == "ready"
    assert report["required_modules"]["available"]["available"] is True
    assert report["optional_modules"]["missing_optional"]["available"] is False
    assert report["optional_modules_block_core_startup"] is False


def test_python_below_310_fails_readiness():
    report = check_installation(
        ROOT,
        required_modules=("available",),
        optional_modules=(),
        importer=_importer,
        version_info=(3, 9, 19),
    )

    assert report["status"] == "not_ready"
    assert report["python"]["supported"] is False


def test_runtime_module_and_package_metadata_drift_is_reported():
    report = check_installation(
        ROOT,
        required_modules=("rdkit",),
        optional_modules=(),
        importer=_importer,
        version_info=(3, 10, 0),
    )

    assert report["runtime_ready"] is True
    assert report["strict_reproducibility_ready"] is False
    assert report["status"] == "ready_with_version_drift"
    assert report["version_drift"][0]["module"] == "rdkit"
