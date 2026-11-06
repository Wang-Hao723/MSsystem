from __future__ import annotations

from matching.evitwin_runtime import (
    resolve_evitwin_query_instrument_family,
    resolve_evitwin_query_method_id,
)


def test_qtof_query_maps_to_specific_registry_token_only_after_family_match():
    result = resolve_evitwin_query_instrument_family(
        {
            "vendor_conversion": {
                "source": {
                    "acquisition": {
                        "instrument_models": ["Agilent instrument model"],
                        "devices": ["Q-TOF"],
                    }
                }
            }
        },
        reference_instrument_family="agilent-qtof-6546",
    )

    assert result["query_instrument_family"] == "qtof"
    assert result["instrument_family_compatible"] is True
    assert result["scoring_instrument_family"] == "agilent-qtof-6546"


def test_orbitrap_query_never_inherits_qtof_reference_identity():
    result = resolve_evitwin_query_instrument_family(
        {
            "vendor_conversion": {
                "source": {
                    "acquisition": {"instrument_models": ["Q Exactive"]}
                }
            }
        },
        reference_instrument_family="agilent-qtof-6546",
    )

    assert result["query_instrument_family"] == "orbitrap"
    assert result["instrument_family_compatible"] is False
    assert result["scoring_instrument_family"] == "orbitrap"


def test_unknown_query_instrument_fails_closed():
    result = resolve_evitwin_query_instrument_family(
        {"filename": "unknown.mzML"},
        reference_instrument_family="agilent-qtof-6546",
    )

    assert result["query_instrument_family"] == "unknown"
    assert result["instrument_family_compatible"] is False
    assert result["scoring_instrument_family"] == "unknown"


def test_query_method_uses_exact_mzml_acq_method_hash(tmp_path):
    path = tmp_path / "blind.mzML"
    path.write_text(
        '<mzML><sourceFile id="AcqMethod.xml" name="AcqMethod.xml" '
        'location="file:///batch/sample.d/AcqData">\n'
        '<cvParam name="SHA-1" '
        'value="d9fb1be0bfc0adb646f329ad5d8bbfd5d1e043a0"/>'
        '</sourceFile></mzML>',
        encoding="utf-8",
    )

    result = resolve_evitwin_query_method_id({"analysis_file_path": str(path)})

    assert result["query_method_id"] == (
        "acqmethod-sha1:d9fb1be0bfc0adb646f329ad5d8bbfd5d1e043a0"
    )
    assert result["resolution_source"] == "mzml_source_file_acq_method_sha1"
    assert result["method_inferred_from_filename"] is False


def test_declared_query_method_takes_precedence_over_mzml_hash(tmp_path):
    path = tmp_path / "blind.mzML"
    path.write_text("<mzML/>", encoding="utf-8")

    result = resolve_evitwin_query_method_id({
        "analysis_file_path": str(path),
        "params": {"rt_method_id": "validated-c18-method"},
    })

    assert result["query_method_id"] == "validated-c18-method"
    assert result["resolution_source"] == "saved_task_rt_method_id"
