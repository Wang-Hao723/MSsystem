"""Convert the fixed MassBank records into query and reference MGF files."""
from decimal import Decimal
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
QUERY = ("MSBNK-IPB_Halle-PB002410", "MSBNK-IPB_Halle-PB002412")
REFERENCE = ("MSBNK-IPB_Halle-PB002411", "MSBNK-RIKEN-PR040049")
COMMIT = "705afb7bccc3b2c42410a744eef73674716a60ef"


def convert(accession, role):
    raw = (ROOT / "massbank_records" / f"{accession}.txt").read_text(encoding="utf-8")
    header, peak_block = raw.split("PK$PEAK: m/z int. rel.int.\n", maxsplit=1)
    fields = dict(line.split(": ", maxsplit=1) for line in header.splitlines()
                  if ": " in line and not line.startswith("CH$LINK: "))
    inchi = next(line.split("INCHIKEY ", 1)[1] for line in header.splitlines()
                 if line.startswith("CH$LINK: INCHIKEY "))
    assert fields["ACCESSION"] == accession
    assert fields["LICENSE"] == "CC BY-SA"
    assert "AC$MASS_SPECTROMETRY: ION_MODE NEGATIVE" in header
    assert any(line.endswith(" [M-H]-") for line in header.splitlines()
               if line.startswith("MS$FOCUSED_ION: "))
    exact = Decimal(fields["CH$EXACT_MASS"])
    precursor = exact - Decimal("1.007276466621")
    ce = next(line.split("COLLISION_ENERGY ", 1)[1] for line in header.splitlines()
              if line.startswith("AC$MASS_SPECTROMETRY: COLLISION_ENERGY "))
    peaks = [line.split()[:2] for line in peak_block.split("//", 1)[0].splitlines() if line.strip()]
    assert len(peaks) == int(fields["PK$NUM_PEAK"])
    metadata = [
        "BEGIN IONS", f"TITLE={accession}", f"PEPMASS={precursor:.7f}", "CHARGE=1-",
        "IONMODE=NEGATIVE", "ADDUCT=[M-H]-", f"COLLISIONENERGY={ce}",
    ]
    if role == "reference":
        metadata += [
            "NAME=Quercetin", f"SMILES={fields['CH$SMILES']}", f"INCHIKEY={inchi}",
            f"FORMULA={fields['CH$FORMULA']}", f"EXACTMASS={exact}", f"SOURCE_ACCESSION={accession}",
            "SOURCE_TIER=experimental", "INTENSITY_ORIGIN=measured",
            f"INSTRUMENT_TYPE={fields['AC$INSTRUMENT_TYPE']}", f"INSTRUMENT={fields['AC$INSTRUMENT']}",
            "SPECTRUM_TYPE=MS2",
        ]
    contributor = accession.split("-")[1]
    provenance = {
        "accession": accession, "role": role, "authors": fields["AUTHORS"],
        "licence_as_recorded": fields["LICENSE"], "collision_energy": ce,
        "peaks": len(peaks), "precursor_mz": float(f"{precursor:.7f}"),
        "precursor_mz_basis": "recorded neutral monoisotopic mass minus proton mass",
        "record": f"https://massbank.eu/MassBank/RecordDisplay?id={accession}",
        "source": f"https://raw.githubusercontent.com/MassBank/MassBank-data/{COMMIT}/{contributor}/{accession}.txt",
    }
    return "\n".join([*metadata, *(" ".join(pair) for pair in peaks), "END IONS", ""]), provenance


def main():
    provenance = []
    for accessions, role, relative_file in (
        (QUERY, "query", "example.mgf"),
        (REFERENCE, "reference", "reference/MassBank_quercetin.mgf"),
    ):
        target = ROOT / relative_file
        target.parent.mkdir(exist_ok=True)
        converted = [convert(acc, role) for acc in accessions]
        target.write_text("\n".join(block for block, _ in converted), encoding="utf-8")
        provenance.extend(row for _, row in converted)
    (ROOT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
