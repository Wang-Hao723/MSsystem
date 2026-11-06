from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace


def test_verified_preflight_executes_immutable_docker_image(monkeypatch, tmp_path):
    from utils import vendor_raw_conversion

    mutable_tag = "proteowizard/test:mutable-tag"
    immutable_identity = f"sha256:{'a' * 64}"
    raw = tmp_path / "sample.raw"
    raw.write_bytes(b"vendor raw")
    output = tmp_path / "sample.mzML"
    raw_report = {
        "status": "ok",
        "source_kind": "thermo_raw",
        "source_path": str(raw.resolve()),
        "format": "Thermo Fisher RAW",
        "acquisition": {
            "ms1_expected": None,
            "ms2_expected": True,
            "polarities": ["NEG"],
        },
        "source_file_signature": {},
    }
    status = {
        "local_msconvert": {"available": False, "path": None},
        "docker": {
            "path": "/usr/local/bin/docker",
            "cli_available": True,
            "daemon_available": True,
            "image": mutable_tag,
            "image_available": True,
            "image_identity": immutable_identity,
            "wine_runner": "wine64_anyuser",
            "wine_runner_verified": True,
        },
        "recommended_backend": "docker",
    }
    executed_commands: list[list[str]] = []

    def fake_run(command, **kwargs):
        executed_commands.append(list(command))
        output_mount = next(
            token for token in command if str(token).endswith(":/output")
        )
        host_output = Path(str(output_mount)[: -len(":/output")])
        output_name = command[command.index("--outfile") + 1]
        (host_output / output_name).write_bytes(b"converted mzML")
        return SimpleNamespace(returncode=0, stdout="converted", stderr="")

    def fake_validate(path, **kwargs):
        selected = Path(path)
        return {
            "status": "valid",
            "path": str(selected),
            "size_bytes": selected.stat().st_size,
            "sha256": hashlib.sha256(selected.read_bytes()).hexdigest(),
        }

    monkeypatch.setattr(vendor_raw_conversion.subprocess, "run", fake_run)
    monkeypatch.setattr(vendor_raw_conversion, "validate_mzml", fake_validate)

    manifest = vendor_raw_conversion.convert_thermo_raw(
        raw,
        output,
        backend="docker",
        docker_image=mutable_tag,
        accept_vendor_licenses=True,
        converter_status_report=status,
        raw_report=raw_report,
    )

    assert len(executed_commands) == 1
    command = executed_commands[0]
    assert immutable_identity in command
    assert mutable_tag not in command
    assert command[command.index(immutable_identity) + 1] == (
        vendor_raw_conversion._docker_wine_runner()
    )
    assert manifest["converter"]["image"] == mutable_tag
    assert manifest["converter"]["image_identity"] == immutable_identity
    assert manifest["converter"]["execution_image"] == immutable_identity
    audit = manifest["conversion_policy"]["command_audit_without_source_path"]
    assert immutable_identity in audit
    assert mutable_tag not in audit
    assert manifest["conversion_policy"]["converter_preflight_reused"] is True

