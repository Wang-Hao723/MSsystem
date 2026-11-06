"""Convert fixed MassBank records into the competing-candidate entity example inputs.

Query: three positive-mode scans of isoquercitrin at 10/20/40 eV (BGC Munich, CC BY).
References: isoquercitrin (Fiocruz, CC BY-NC), hyperoside (RIKEN ReSpect, CC BY-NC,
same 2D connectivity, different stereochemistry) and spiraeoside (RIKEN ReSpect,
CC BY-NC, different 2D connectivity) within the same precursor window.
"""
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
COMMIT = "705afb7bccc3b2c42410a744eef73674716a60ef"
PROTON = 1.007276466621

QUERIES = {
    "MSBNK-BGC_Munich-RP017101": "Isoquercitrin",
    "MSBNK-BGC_Munich-RP017102": "Isoquercitrin",
    "MSBNK-BGC_Munich-RP017103": "Isoquercitrin",
}
REFERENCES = {
    "MSBNK-Fiocruz-FIO00160": "Isoquercitrin",
    "MSBNK-Fiocruz-FIO00161": "Isoquercitrin",
    "MSBNK-Fiocruz-FIO00162": "Isoquercitrin",
    "MSBNK-Fiocruz-FIO00163": "Isoquercitrin",
    "MSBNK-RIKEN_ReSpect-PS046501": "Hyperoside",
    "MSBNK-RIKEN_ReSpect-PS046502": "Hyperoside",
    "MSBNK-RIKEN_ReSpect-PS046503": "Hyperoside",
    "MSBNK-RIKEN_ReSpect-PS086901": "Spiraeoside",
    "MSBNK-RIKEN_ReSpect-PS086902": "Spiraeoside",
    "MSBNK-RIKEN_ReSpect-PS086903": "Spiraeoside",
}
DATASET_SOURCE = {
    "BGC_Munich": "BGC_Munich",
    "Fiocruz": "Fiocruz",
    "RIKEN_ReSpect": "RIKEN_ReSpect",
}

MONOISOTOPIC = {"C": 12.0, "H": 1.00782503223, "O": 15.99491461957, "N": 14.00307400443,
                "S": 31.9720711744, "P": 30.97376199842, "Na": 22.9897692820}


def formula_mass(formula):
    total = 0.0
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        total += MONOISOTOPIC[element] * (int(count) if count else 1)
    return total


def ce_number(text):
    match = re.match(r"\s*(\d+(?:\.\d+)?)", str(text))
    return float(match.group(1)) if match else None


def convert(accession, display_name, role):
    raw_path = next((ROOT / "massbank_records").glob(f"{accession}.txt"))
    raw = raw_path.read_text(encoding="utf-8")
    header, peak_block = raw.split("PK$PEAK: m/z int. rel.int.\n", maxsplit=1)
    fields = dict(line.split(": ", maxsplit=1) for line in header.splitlines()
                  if ": " in line and not line.startswith(("CH$LINK:", "CH$NAME:", "AC$", "MS$", "COMMENT")))
    names = [line.split("CH$NAME: ", 1)[1] for line in header.splitlines() if line.startswith("CH$NAME: ")]
    inchi = next(line.split("INCHIKEY ", 1)[1] for line in header.splitlines()
                 if line.startswith("CH$LINK: INCHIKEY "))
    exact = formula_mass(fields["CH$FORMULA"])
    recorded = float(fields.get("CH$EXACT_MASS", "0") or 0)
    mass_basis = "formula monoisotopic"
    if recorded and abs(recorded - exact) / exact <= 5e-5:
        exact = recorded
        mass_basis = "recorded CH$EXACT_MASS"
    precursor = exact + PROTON
    ce_raw = next(line.split("COLLISION_ENERGY ", 1)[1].strip() for line in header.splitlines()
                  if line.startswith("AC$MASS_SPECTROMETRY: COLLISION_ENERGY "))
    ce = ce_number(ce_raw)
    assert fields["ACCESSION"] == accession
    assert "ION_MODE POSITIVE" in header
    peaks = [line.split()[:2] for line in peak_block.split("//", 1)[0].splitlines() if line.strip()]
    assert len(peaks) == int(fields["PK$NUM_PEAK"]), accession
    assert ce is not None, accession
    contributor = accession.split("-")[1]
    dataset_source = DATASET_SOURCE[contributor]
    metadata = [
        "BEGIN IONS", f"TITLE={accession}", f"PEPMASS={precursor:.7f}", "CHARGE=1+",
        "IONMODE=POSITIVE", "ADDUCT=[M+H]+", f"COLLISIONENERGY={ce:g} eV", f"RT={fields.get('AC$CHROMATOGRAPHY: RETENTION_TIME', '').split(': ')[-1] or 0}",
    ]
    if role == "reference":
        metadata += [
            f"NAME={display_name}", f"SMILES={fields['CH$SMILES']}", f"INCHIKEY={inchi}",
            f"FORMULA={fields['CH$FORMULA']}", f"EXACTMASS={exact:.4f}", f"SOURCE_ACCESSION={accession}",
            f"DATASET_SOURCE={dataset_source}", "SOURCE_TIER=experimental", "INTENSITY_ORIGIN=measured",
            f"INSTRUMENT_TYPE={fields.get('AC$INSTRUMENT_TYPE', '')}", f"INSTRUMENT={fields.get('AC$INSTRUMENT', '')}",
            "SPECTRUM_TYPE=MS2",
        ]
    provenance = {
        "accession": accession, "role": role, "display_name": display_name,
        "record_synonyms": names, "authors": fields.get("AUTHORS", ""),
        "licence_as_recorded": fields["LICENSE"], "collision_energy_as_recorded": ce_raw,
        "collision_energy_eV": ce, "peaks": len(peaks),
        "precursor_mz": round(precursor, 7), "precursor_adduct": "[M+H]+",
        "exact_mass_used": round(exact, 4), "exact_mass_basis": mass_basis,
        "formula": fields["CH$FORMULA"], "inchikey": inchi,
        "record_sha256": hashlib.sha256(raw.encode("utf-8")).hexdigest(),
        "record": f"https://massbank.eu/MassBank/RecordDisplay?id={accession}",
        "source": f"https://raw.githubusercontent.com/MassBank/MassBank-data/{COMMIT}/{contributor}/{accession}.txt",
        "massbank_data_release": "2026.03",
    }
    if mass_basis == "formula monoisotopic" and recorded:
        provenance["recorded_exact_mass_discarded"] = recorded
        provenance["discard_reason"] = "differs from formula monoisotopic mass by more than 50 ppm"
    return "\n".join([*metadata, *(" ".join(pair) for pair in peaks), "END IONS", ""]), provenance


def main():
    provenance = []
    for records, role, target in (
        (QUERIES, "query", ROOT / "entity_queries.mgf"),
        (REFERENCES, "reference", ROOT / "reference" / "flavonoid_pos_refs.mgf"),
    ):
        blocks = []
        for accession, display_name in records.items():
            block, record_provenance = convert(accession, display_name, role)
            blocks.append(block)
            provenance.append(record_provenance)
        target.parent.mkdir(exist_ok=True)
        target.write_text("\n".join(blocks), encoding="utf-8")
    (ROOT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
