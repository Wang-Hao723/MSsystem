"""Measure library-build time, scan-search latency and peak memory by library size.

Libraries of increasing size are built from a deterministic sample of real
positive-mode MassBank 2026.03 MS2 records; ten fixed query scans are then
searched against each library with the shipped matcher. Each size runs in a
fresh process; results append as one CSV row per size.
"""
import argparse
import csv
import json
import random
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MB_ROOT = Path(__file__).resolve().parent / "massbank_sample"
PROTON = 1.007276466621
MONOISOTOPIC = {"C": 12.0, "H": 1.00782503223, "O": 15.99491461957, "N": 14.00307400443,
                "S": 31.9720711744, "P": 30.97376199842, "F": 18.99840316273, "Cl": 34.968852682,
                "Br": 78.9183376, "I": 126.9044719, "Si": 27.97692653446, "Se": 79.9165218}
SAMPLE_SEED = 20261103
QUERY_ACCESSIONS = [
    "MSBNK-BGC_Munich-RP017101", "MSBNK-BGC_Munich-RP017102", "MSBNK-BGC_Munich-RP017103",
    "MSBNK-Fiocruz-FIO00160", "MSBNK-Fiocruz-FIO00161", "MSBNK-Fiocruz-FIO00162",
    "MSBNK-RIKEN_ReSpect-PS086903", "MSBNK-RIKEN_ReSpect-PS046501",
    "MSBNK-RIKEN_ReSpect-PS046502", "MSBNK-RIKEN_ReSpect-PS046503",
]


def formula_mass(formula):
    total = 0.0
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        total += MONOISOTOPIC[element] * (int(count) if count else 1)
    return total


def parse_record(path, min_peaks=5):
    raw = path.read_text(encoding="utf-8", errors="ignore")
    if "ION_MODE POSITIVE" not in raw or "MS_TYPE MS2" not in raw:
        return None
    if "PK$PEAK: m/z int. rel.int." not in raw:
        return None
    header, peak_block = raw.split("PK$PEAK: m/z int. rel.int.\n", maxsplit=1)
    fields = {}
    inchikey = formula = smiles = None
    for line in header.splitlines():
        if line.startswith("CH$LINK: INCHIKEY ") and inchikey is None:
            inchikey = line.split("INCHIKEY ", 1)[1].strip()
        elif line.startswith("CH$FORMULA: ") and formula is None:
            formula = line.split(": ", 1)[1].strip()
        elif line.startswith("CH$SMILES: ") and smiles is None:
            smiles = line.split(": ", 1)[1].strip()
        elif ": " in line and not line.startswith(("CH$NAME:", "CH$LINK:", "AC$", "MS$", "COMMENT", "CH$IUPAC")):
            key, value = line.split(": ", 1)
            fields.setdefault(key, value)
    if not (inchikey and formula and smiles):
        return None
    try:
        peaks = [line.split()[:2] for line in peak_block.split("//", 1)[0].splitlines() if line.strip()]
    except IndexError:
        return None
    if len(peaks) < min_peaks:
        return None
    try:
        precursor = formula_mass(formula) + PROTON
    except KeyError:
        return None
    ce_raw = next((line.split("COLLISION_ENERGY ", 1)[1].strip() for line in header.splitlines()
                   if line.startswith("AC$MASS_SPECTROMETRY: COLLISION_ENERGY ")), "0")
    ce = re.match(r"\s*(\d+(?:\.\d+)?)", ce_raw)
    accession = fields.get("ACCESSION", path.stem)
    block = ["BEGIN IONS", f"TITLE={accession}", f"PEPMASS={precursor:.7f}", "CHARGE=1+",
             "IONMODE=POSITIVE", "ADDUCT=[M+H]+", f"COLLISIONENERGY={ce.group(1) if ce else 0} eV",
             f"NAME={accession}", f"SMILES={smiles}", f"INCHIKEY={inchikey}", f"FORMULA={formula}",
             f"EXACTMASS={formula_mass(formula):.4f}", f"SOURCE_ACCESSION={accession}",
             "DATASET_SOURCE=MassBank202603", "SOURCE_TIER=experimental", "INTENSITY_ORIGIN=measured",
             f"INSTRUMENT_TYPE={fields.get('AC$INSTRUMENT_TYPE', '')}",
             f"INSTRUMENT={fields.get('AC$INSTRUMENT', '')}", "SPECTRUM_TYPE=MS2",
             *((" ".join(pair)) for pair in peaks), "END IONS", ""]
    return "\n".join(block)


def ensure_sample(max_records):
    index_file = MB_ROOT / "sample_index.txt"
    if not index_file.exists():
        source = Path(sys.argv[0]).resolve().parent / "massbank_source"
        if not source.exists():
            raise SystemExit("Place the extracted MassBank-data tree under benchmarks/massbank_source first")
        candidates = []
        for path in source.rglob("*.txt"):
            if len(candidates) >= 60000:
                break
            candidates.append(path)
        rng = random.Random(SAMPLE_SEED)
        rng.shuffle(candidates)
        MB_ROOT.mkdir(parents=True, exist_ok=True)
        blocks = []
        for accession in QUERY_ACCESSIONS:
            hit = next(source.rglob(f"{accession}.txt"))
            block = parse_record(hit, min_peaks=2)
            assert block, accession
            blocks.append(block)
        for path in candidates:
            if len(blocks) >= max_records + len(QUERY_ACCESSIONS):
                break
            block = parse_record(path)
            if block:
                blocks.append(block)
        (MB_ROOT / "refs_all.mgf").write_text("\n".join(blocks), encoding="utf-8")
        index_file.write_text(f"{len(blocks)}\n", encoding="utf-8")
    return int(index_file.read_text().strip())


def run_size(size, out_csv):
    total = ensure_sample(6000)
    assert size <= total, f"requested {size} > sampled {total}"
    work = Path(tempfile.mkdtemp(prefix=f"bench-{size}-"))
    (work / "library").mkdir()
    (work / "runtime").mkdir()
    (work / "unused").mkdir()
    ref_dir = work / "refs"
    ref_dir.mkdir()
    text = (MB_ROOT / "refs_all.mgf").read_text(encoding="utf-8")
    blocks = [b for b in text.split("\nEND IONS\n") if b.strip()]
    (ref_dir / "refs.mgf").write_text(
        "\n".join(blocks[:size]).replace("BEGIN IONS", "BEGIN IONS") + ("\nEND IONS\n" if blocks else ""),
        encoding="utf-8")

    env_cmd = ["env", f"PATH={__import__('os').environ['PATH']}",
               f"MSSYSTEM_SHARED_DATA_ROOT={work}/unused",
               f"MSSYSTEM_RUNTIME_ROOT={work}/runtime", "MSSYSTEM_LIBRARY_WRITE_ENABLED=0"]
    t0 = time.perf_counter()
    subprocess.run(env_cmd + [str(REPO / ".venv/bin/python"), "build_from_mgf_msp.py",
                              "--mgf-dir", str(ref_dir), "--output", str(work / "library")],
                   cwd=REPO, check=True, capture_output=True)
    build_seconds = time.perf_counter() - t0

    script = f'''
import sys, time, json, resource
sys.path.insert(0, {str(REPO)!r})
from matching.search import load_library, search
records, index = load_library({str(work / "library" / "local_library.db")!r})
import re
def parse(path):
    block, out = {{"peaks": []}}, []
    started = False
    for line in open(path):
        line = line.strip()
        if line == "BEGIN IONS":
            block = {{"peaks": []}}; started = True; continue
        if line == "END IONS":
            out.append(block); continue
        if not started or not line: continue
        if "=" in line and not line[0].isdigit():
            k, v = line.split("=", 1); block[k.lower()] = v
        else:
            parts = line.split()
            if len(parts) >= 2: block["peaks"].append([float(parts[0]), float(parts[1])])
    return out
queries = parse({str(MB_ROOT / "queries.mgf")!r})
times = []
matched = 0
for q in queries:
    peaks = q["peaks"]; pmz = float(q["pepmass"].split()[0])
    ce = float(re.match(r"(\\d+(?:\\.\\d+)?)", q.get("collisionenergy", "0")).group(1))
    t0 = time.perf_counter()
    r = search(peaks, pmz, records, index, query_polarity="POS", enforce_polarity=True,
           query_adduct="[M+H]+", query_adduct_reliable=True, query_charge=1,
           query_collision_energy=ce, top_n=10)
    if r is not None:
        matched += 1
    times.append((time.perf_counter() - t0) * 1000.0)
times.sort()
raw_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
rss = raw_rss / 1048576.0 if sys.platform == "darwin" else raw_rss / 1024.0
print(json.dumps({{"n_records": len(records), "search_ms_median": times[len(times)//2],
                   "search_ms_p95": times[int(len(times)*0.95)-1], "search_ms_min": times[0],
                   "search_ms_max": times[-1], "queries_matched": matched, "peak_rss_mb": rss}}))
'''
    (work / "probe.py").write_text(script, encoding="utf-8")
    probe = subprocess.run(env_cmd + [str(REPO / ".venv/bin/python"), str(work / "probe.py")],
                           cwd=REPO, check=True, capture_output=True, text=True)
    stats = json.loads(probe.stdout.strip().splitlines()[-1])
    row = {"library_refs": size, "records_loaded": stats["n_records"],
           "build_seconds": round(build_seconds, 2),
           "search_ms_median": round(stats["search_ms_median"], 2),
           "search_ms_p95": round(stats["search_ms_p95"], 2),
           "search_ms_min": round(stats["search_ms_min"], 2),
           "search_ms_max": round(stats["search_ms_max"], 2),
           "queries_matched": stats["queries_matched"],
           "peak_rss_mb": round(stats["peak_rss_mb"], 1)}
    new = not out_csv.exists()
    with out_csv.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        if new:
            writer.writeheader()
        writer.writerow(row)
    print(row)
    shutil.rmtree(work, ignore_errors=True)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--size", type=int, required=True)
    cli.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "search_scaling_results.csv")
    args = cli.parse_args()
    if not (MB_ROOT / "queries.mgf").exists():
        MB_ROOT.mkdir(parents=True, exist_ok=True)
        source = Path(__file__).resolve().parent / "massbank_source"
        blocks = []
        for accession in QUERY_ACCESSIONS:
            hit = next(source.rglob(f"{accession}.txt"))
            block = parse_record(hit, min_peaks=2)
            assert block, accession
            blocks.append(block)
        (MB_ROOT / "queries.mgf").write_text("\n".join(blocks), encoding="utf-8")
    run_size(args.size, args.out)


if __name__ == "__main__":
    main()
