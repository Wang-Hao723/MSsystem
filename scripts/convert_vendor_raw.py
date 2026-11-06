#!/usr/bin/env python3
"""Convert Agilent MassHunter .d or Thermo Fisher .raw to validated mzML."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.vendor_raw_conversion import (  # noqa: E402
    VendorRawConversionError,
    convert_vendor_raw,
    converter_status,
    inspect_vendor_raw,
    validate_mzml,
    write_conversion_manifest,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert proprietary Agilent .d or Thermo .raw data with ProteoWizard, "
            "then verify that MS1 and MS2 were retained."
        )
    )
    parser.add_argument(
        "input", type=Path, nargs="?", help="Agilent .d directory or Thermo .raw file"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output .mzML path (default: <input parent>/<input stem>.mzML)",
    )
    parser.add_argument(
        "--backend",
        choices=("auto", "local", "docker"),
        default="auto",
    )
    parser.add_argument("--msconvert-bin", type=Path)
    parser.add_argument("--docker-image")
    parser.add_argument(
        "--accept-vendor-licenses",
        action="store_true",
        help="Confirm acceptance required by the ProteoWizard vendor-reader image.",
    )
    parser.add_argument(
        "--pull",
        action="store_true",
        help="Allow a one-time pull of the ProteoWizard Docker image.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--inspect",
        action="store_true",
        help="Inspect and fingerprint the vendor source without converting.",
    )
    parser.add_argument(
        "--validate-mzml",
        type=Path,
        help="Validate an existing mzML instead of converting.",
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        help=(
            "Matching Agilent .d or Thermo .raw for --validate-mzml; enables "
            "acquisition expectations and paired source checksums."
        ),
    )
    parser.add_argument(
        "--expected-polarity",
        choices=("POS", "NEG"),
        action="append",
        default=[],
        help="Optional expected polarity; repeat only for an intentionally mixed run.",
    )
    parser.add_argument(
        "--peak-picking-strategy",
        choices=("auto", "preserve", "vendor", "openms"),
        default="auto",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="Show converter backend availability.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.status:
            payload = converter_status(
                msconvert_bin=args.msconvert_bin,
                docker_image=args.docker_image,
            )
        elif args.validate_mzml:
            paired_raw = args.raw_dir or args.input
            if args.raw_dir and args.input:
                if args.raw_dir.expanduser().resolve() != args.input.expanduser().resolve():
                    raise ValueError(
                        "Specify the paired raw source once: positional input or --raw-dir."
                    )
            if paired_raw:
                raw_report = inspect_vendor_raw(
                    paired_raw,
                    expected_polarities=args.expected_polarity,
                )
                acquisition = raw_report.get("acquisition") or {}
                expected_source_sha1 = {
                    Path(relative).name: details.get(
                        "sha1_for_mzml_provenance", ""
                    )
                    for relative, details in (
                        raw_report.get("required_file_signatures") or {}
                    ).items()
                }
                source_signature = raw_report.get("source_file_signature") or {}
                if source_signature:
                    expected_source_sha1[Path(paired_raw).name] = (
                        source_signature.get("sha1_for_mzml_provenance", "")
                    )
                payload = validate_mzml(
                    args.validate_mzml,
                    expected_ms1=acquisition.get("ms1_expected"),
                    expected_ms2=bool(acquisition.get("ms2_expected", True)),
                    expected_polarities=acquisition.get("polarities") or [],
                    isolation_width_fallback_da=acquisition.get(
                        "isolation_width_da"
                    ),
                    expected_source_sha1=expected_source_sha1,
                )
                payload["paired_raw"] = raw_report
            else:
                payload = validate_mzml(args.validate_mzml)
        elif args.inspect:
            if args.raw_dir:
                raise ValueError("--raw-dir is only used with --validate-mzml")
            if not args.input:
                raise ValueError("vendor raw input is required with --inspect")
            payload = inspect_vendor_raw(
                args.input,
                expected_polarities=args.expected_polarity,
            )
        else:
            if args.raw_dir:
                raise ValueError("--raw-dir is only used with --validate-mzml")
            if not args.input:
                raise ValueError("Agilent .d or Thermo .raw input is required")
            output = args.output or args.input.with_suffix(".mzML")
            payload = convert_vendor_raw(
                args.input,
                output,
                expected_polarities=args.expected_polarity,
                backend=args.backend,
                msconvert_bin=args.msconvert_bin,
                docker_image=args.docker_image,
                accept_vendor_licenses=args.accept_vendor_licenses,
                allow_container_pull=args.pull,
                overwrite=args.overwrite,
                peak_picking_strategy=args.peak_picking_strategy,
            )
            manifest_path = output.with_suffix(output.suffix + ".conversion.json")
            write_conversion_manifest(payload, manifest_path)
            payload["manifest_path"] = str(manifest_path.resolve())
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    except (VendorRawConversionError, OSError, ValueError) as exc:
        print(
            json.dumps(
                {"status": "error", "error": str(exc)},
                ensure_ascii=False,
                indent=2,
            ),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
