from __future__ import annotations

import hashlib

from matching.feature_fragment_benchmark_v3_1.raw_metadata import (
    extract_mzml_spectrum_metadata,
    extract_mzxml_scan_metadata,
    load_psi_ms_mapping,
    validate_psi_ms_mapping,
)


def test_correct_psi_ms_terms_and_mzml_origin(tmp_path) -> None:
    content = """<?xml version="1.0" encoding="UTF-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="IC1">
      <cvParam accession="MS:1000484" name="orbitrap"/>
    </instrumentConfiguration>
  </instrumentConfigurationList>
  <run id="run1">
    <spectrumList count="1">
      <spectrum id="controllerType=0 controllerNumber=1 scan=7" index="0">
        <cvParam accession="MS:1000511" name="ms level" value="2"/>
        <cvParam accession="MS:1000130" name="positive scan"/>
        <scanList count="1"><scan instrumentConfigurationRef="IC1">
          <scanWindowList count="1"><scanWindow>
            <cvParam accession="MS:1000501" name="scan window lower limit" value="50"/>
            <cvParam accession="MS:1000500" name="scan window upper limit" value="1000"/>
          </scanWindow></scanWindowList>
        </scan></scanList>
        <precursorList count="1"><precursor>
          <isolationWindow>
            <cvParam accession="MS:1000827" name="isolation window target m/z" value="400.2"/>
            <cvParam accession="MS:1000828" name="isolation window lower offset" value="0.4"/>
            <cvParam accession="MS:1000829" name="isolation window upper offset" value="0.7"/>
          </isolationWindow>
          <selectedIonList count="1"><selectedIon>
            <cvParam accession="MS:1000744" name="selected ion m/z" value="400.2"/>
            <cvParam accession="MS:1000041" name="charge state" value="2"/>
          </selectedIon></selectedIonList>
          <activation>
            <cvParam accession="MS:1002481" name="higher energy beam-type collision-induced dissociation"/>
            <cvParam accession="MS:1002482" name="statistical threshold" value="0.95"/>
            <cvParam accession="MS:1000045" name="collision energy" value="30" unitAccession="UO:0000266" unitName="electronvolt"/>
          </activation>
        </precursor></precursorList>
      </spectrum>
    </spectrumList>
  </run>
</mzML>
"""
    path = tmp_path / "probe.mzML"
    path.write_text(content, encoding="utf-8")
    row = extract_mzml_spectrum_metadata(path, index=0)
    assert row["value_origin"] == "mzml_cv_param"
    assert row["source_file_sha256"] == hashlib.sha256(content.encode()).hexdigest()
    assert row["instrument_configuration_ref"] == "IC1"
    assert row["activation"]["methods"][0]["accession"] == "MS:1002481"
    assert row["activation"]["non_activation_terms_seen"][0]["accession"] == "MS:1002482"
    assert row["activation"]["collision_energy"]["unit"] == "eV"
    assert row["activation"]["collision_energy"]["values"] == [30.0]
    assert row["precursor_isolation"]["status"] == "verified"
    assert row["precursor_isolation"]["lower_mz"] == 399.8
    assert row["precursor_isolation"]["upper_mz"] == 400.9
    assert row["scan_range"]["minimum_mz"] == 50.0
    assert row["scan_range"]["maximum_mz"] == 1000.0
    assert row["selected_ion"]["charge"] == 2
    assert row["conversion_applied"] is False


def test_psi_mapping_rejects_statistical_threshold_as_activation() -> None:
    mapping = load_psi_ms_mapping()
    validate_psi_ms_mapping(mapping)
    assert "MS:1000422" in mapping["mzml"]["activation"]
    assert "MS:1002481" in mapping["mzml"]["activation"]
    assert "MS:1002472" in mapping["mzml"]["activation"]
    assert "MS:1002482" not in mapping["mzml"]["activation"]
    assert mapping["mzml"]["terms_not_activation"]["MS:1002482"]["canonical_name"] == "statistical threshold"


def test_mzxml_attributes_never_claim_mzml_cv_origin(tmp_path) -> None:
    content = """<?xml version="1.0"?>
<mzXML><msRun><scan num="7" msLevel="2" polarity="+" lowMz="50" highMz="1000">
  <precursorMz precursorCharge="1" activationMethod="CID" collisionEnergy="25" collisionEnergyUnit="eV" windowWideness="1.2">300.2</precursorMz>
</scan></msRun></mzXML>
"""
    path = tmp_path / "probe.mzXML"
    path.write_text(content, encoding="utf-8")
    row = extract_mzxml_scan_metadata(path, scan_number=7)
    assert row["value_origin"] == "mzxml_attribute"
    assert row["activation"]["value_origin"] == "mzxml_attribute"
    assert row["activation"]["collision_energy"]["value_origin"] == "mzxml_attribute"
    assert row["activation"]["collision_energy"]["status"] == "resolved"
    assert row["precursor_isolation"]["lower_offset_mz"] == 0.6
    assert row["precursor_isolation"]["upper_offset_mz"] == 0.6
    assert row["scan_range"]["status"] == "resolved"
    assert row["conversion_applied"] is False


def test_mzxml_numeric_collision_energy_without_unit_stays_unresolved(tmp_path) -> None:
    content = """<mzXML><msRun><scan num="3" msLevel="2" lowMz="10" highMz="500"><precursorMz collisionEnergy="35">250</precursorMz></scan></msRun></mzXML>"""
    path = tmp_path / "unitless.mzXML"
    path.write_text(content, encoding="utf-8")
    row = extract_mzxml_scan_metadata(path, scan_number="3")
    assert row["activation"]["collision_energy"]["status"] == "unit_unresolved"
    assert row["activation"]["collision_energy"]["unit"] is None
    assert "collision_energy_unit" in row["blocking_missing_fields"]
