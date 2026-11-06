# MSsystem executable example for JCI Software review

This bundle contains the November 2, 2026 MSsystem source snapshot and a small executable example. The example builds a SQLite library from two experimental references, runs the shipped MSsystem matcher on two experimental queries, and checks the ranked results and matched peaks against the bundled outputs. It demonstrates scan-level retrieval and evidence export for one compound, quercetin.

## Install and run

The example was tested in a fresh Python 3.12.3 virtual environment on Linux x86-64. Start from the directory containing this README. Installation requires uv and access to PyPI. Analysis uses the included spectra and a newly built local reference library.

    python3.12 -m zipfile -e Additional_File_1_MSsystem_source_20261102.zip .
    uv venv --python python3.12 .venv
    uv --no-cache pip install --python .venv/bin/python --only-binary :all: -r MSsystem_source_20261102/requirements.lock
    bash run_demo.sh "$PWD/.venv/bin/python" "$PWD/MSsystem_source_20261102"

The full resolved dependency versions are listed in verification/requirements-linux-x86_64-py312-tested.txt. The last command creates a fresh work directory and prints its path. Its optional third argument selects the parent directory. The validation run used /dev/shm as the parent because the database builder enforces a free-space reserve. Use a writable volume satisfying that reserve. The script points MSsystem at its own library and runtime paths. The installed environment occupied approximately 1.4 GiB during validation.

The result/ directory contains two files:

- ranking.tsv records query accession, rank, candidate name, InChIKey, reference accession, score and matched-peak count.
- spectral_evidence.json records the selected references, query and library scoring views, and matched peak pairs. The test verifies that the stored comparison cosine can be replayed from those views.

The expected/ directory contains the validated output. This CPU example covers the installed software's scan-level matcher. The manuscript benchmarks, web interface, entity-level aggregation and optional model assets require separate validation.

## Input data and attribution

The four source records in massbank_records/ come from [MassBank-data release 2026.03](https://github.com/MassBank/MassBank-data/tree/705afb7bccc3b2c42410a744eef73674716a60ef). The original headers state **CC BY-SA**. The converted MGF files retain those terms. The included MSsystem source archive carries its MIT licence and third-party notices.

| Role | MassBank accession | Contributor | Collision energy |
| --- | --- | --- | --- |
| Query | [MSBNK-IPB_Halle-PB002410](https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-IPB_Halle-PB002410) | C. Boettcher, Institute of Plant Biochemistry, Halle | 25 eV |
| Query | [MSBNK-IPB_Halle-PB002412](https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-IPB_Halle-PB002412) | C. Boettcher, Institute of Plant Biochemistry, Halle | 55 eV |
| Reference | [MSBNK-IPB_Halle-PB002411](https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-IPB_Halle-PB002411) | C. Boettcher, Institute of Plant Biochemistry, Halle | 40 eV |
| Reference | [MSBNK-RIKEN-PR040049](https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-RIKEN-PR040049) | F. Matsuda, Plant Science Center, RIKEN | Ramp 5–60 eV |

The script prepare_inputs.py converts the four records to MGF and preserves their peak masses and intensities. The MGF precursor value is calculated from the reported neutral exact mass and the negative-mode adduct [M−H]⁻. Query titles use the original accessions. The reference library excludes both query accessions. Full input provenance is in provenance.json.

## Verification

On Linux x86-64, both queries returned quercetin at rank 1. The 25 eV query matched MSBNK-RIKEN-PR040049 with 5 peak pairs. The 55 eV query matched MSBNK-IPB_Halle-PB002411 with 10 peak pairs. The verification script reproduced both rankings and peak-pair records after rebuilding the library in a second empty work directory.
