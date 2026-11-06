"""Strict, source-listed MassBank corpus construction for V15.

The builder reads only files named in an explicit manifest.  Every accepted
record retains its repository accession, license label, chemical identity
provenance, exact acquisition domain, unconverted collision-energy quantity,
quality limitations, and full normalized peak list.  Rejected records are
written with stable reason codes.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import urllib.request
from collections import Counter, defaultdict
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any, TextIO

from .mirror import exact_peak_hash


_FULL_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)"
_FIXED_CE = re.compile(rf"^\s*({_NUMBER})\s*(eV|V|NCE)\s*$", re.IGNORECASE)
_RAMP_CE = re.compile(
    rf"^\s*Ramp\s+({_NUMBER})\s*-\s*({_NUMBER})\s*(eV|V)\s*$",
    re.IGNORECASE,
)
_NOMINAL_PERCENT = re.compile(
    rf"^\s*({_NUMBER})\s*%\s*\(nominal\)\s*$", re.IGNORECASE
)
_DB_ID_ONLY = re.compile(
    r"^(?:NCGC[-_A-Z0-9]*|CCMSLIB[A-Z0-9_-]*|HMDB\d+|CHEBI:\d+|"
    r"CID:\d+|DTXSID\d+|MSBNK[-_A-Z0-9]+|LM[A-Z0-9]+)$",
    re.IGNORECASE,
)
_DB_ID_PREFIX = re.compile(r"^(?:NCGC|CCMSLIB|HMDB\d|MSBNK[-_])", re.IGNORECASE)
_EMPTY_NAME = re.compile(r"^(?:n/?a|unknown|unnamed|no name|not available)$", re.I)


class CorpusBuildError(ValueError):
    """Raised for an invalid corpus manifest or cache."""


class RecordReject(CorpusBuildError):
    """Expected fail-closed rejection of one source record."""

    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}" if detail else code)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def iter_json_array(path: str | Path, *, chunk_size: int = 1 << 20) -> Iterator[Any]:
    """Stream a top-level JSON array without loading the file into memory."""
    decoder = json.JSONDecoder()
    buffer = ""
    position = 0
    started = False
    ended = False
    with Path(path).open("r", encoding="utf-8") as handle:
        while True:
            if position >= len(buffer) and not ended:
                block = handle.read(chunk_size)
                if block:
                    buffer = block
                    position = 0
                else:
                    ended = True
            while position < len(buffer) and buffer[position].isspace():
                position += 1
            if not started:
                if position >= len(buffer):
                    if ended:
                        raise CorpusBuildError(f"empty JSON file: {path}")
                    buffer += handle.read(chunk_size)
                    continue
                if buffer[position] != "[":
                    raise CorpusBuildError(f"top-level JSON value is not an array: {path}")
                started = True
                position += 1
                continue
            while position < len(buffer) and (buffer[position].isspace() or buffer[position] == ","):
                position += 1
            if position < len(buffer) and buffer[position] == "]":
                position += 1
                break
            if position >= len(buffer):
                if ended:
                    raise CorpusBuildError(f"unterminated JSON array: {path}")
                buffer = ""
                position = 0
                continue
            try:
                value, next_position = decoder.raw_decode(buffer, position)
            except json.JSONDecodeError:
                if ended:
                    raise CorpusBuildError(f"invalid JSON array: {path}")
                buffer = buffer[position:] + handle.read(chunk_size)
                position = 0
                continue
            yield value
            position = next_position
            if position > chunk_size:
                buffer = buffer[position:]
                position = 0
        trailing = buffer[position:] + handle.read()
        if trailing.strip():
            raise CorpusBuildError(f"unexpected content after JSON array: {path}")


def iter_jsonl(path: str | Path) -> Iterator[dict[str, Any]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise CorpusBuildError(f"{path}:{line_number}: JSONL value is not an object")
            yield value


def _tags(items: Sequence[Mapping[str, Any]] | None) -> dict[str, str]:
    return {
        str(item.get("subtag")): str(item.get("value"))
        for item in (items or [])
        if item.get("subtag") and item.get("value") is not None
    }


def _links(items: Sequence[Mapping[str, Any]] | None) -> dict[str, list[str]]:
    links: dict[str, list[str]] = defaultdict(list)
    for item in items or []:
        database = str(item.get("database") or "").strip()
        identifier = str(item.get("identifier") or "").strip()
        if database and identifier:
            links[database].append(identifier)
    return dict(links)


def _first_link(links: Mapping[str, Sequence[str]], database: str) -> str:
    values = links.get(database) or []
    return str(values[0]) if values else ""


def stable_compound_name(names: Sequence[Any] | None) -> tuple[str, int]:
    """Select a non-identifier display name and return its source-list index."""
    for index, value in enumerate(names or []):
        name = re.sub(r"\s+", " ", str(value or "")).strip(" ;|")
        if not name or _EMPTY_NAME.fullmatch(name):
            continue
        if _DB_ID_ONLY.fullmatch(name) or _DB_ID_PREFIX.match(name):
            continue
        if re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", name):
            continue
        return name, index
    raise RecordReject("stable_compound_name_missing")


def normalize_formula(formula: Any) -> str:
    text = re.sub(r"\s+", "", str(formula or ""))
    if not text or any(marker in text for marker in (".", "·", "[", "]", "+", "-")):
        raise RecordReject("formula_missing_or_ambiguous", text)
    counts: dict[str, int] = defaultdict(int)
    consumed = ""
    for element, count_text in _FORMULA_TOKEN.findall(text):
        count = int(count_text or "1")
        if count < 1:
            raise RecordReject("formula_invalid_count", text)
        counts[element] += count
        consumed += element + count_text
    if not counts or consumed != text:
        raise RecordReject("formula_unparseable", text)
    order: list[str] = []
    if "C" in counts:
        order.append("C")
    if "H" in counts:
        order.append("H")
    order.extend(sorted(element for element in counts if element not in {"C", "H"}))
    return "".join(
        element + (str(counts[element]) if counts[element] != 1 else "")
        for element in order
    )


def parse_collision_energy(raw_value: Any) -> dict[str, Any]:
    raw = str(raw_value or "").strip()
    match = _NOMINAL_PERCENT.fullmatch(raw)
    if match:
        value = float(match.group(1))
        return {
            "raw": raw,
            "value": value,
            "minimum": value,
            "maximum": value,
            "unit": "percent_nominal",
            "mode": "fixed",
            "conversion_applied": False,
            "unit_provenance": "MassBank COLLISION_ENERGY literal",
            "cross_vendor_or_cross_model_equivalence_assumed": False,
        }
    match = _RAMP_CE.fullmatch(raw)
    if match:
        minimum, maximum = sorted((float(match.group(1)), float(match.group(2))))
        unit = match.group(3)
        unit = "eV" if unit.lower() == "ev" else "V"
        return {
            "raw": raw,
            "value": None,
            "minimum": minimum,
            "maximum": maximum,
            "unit": unit,
            "mode": "ramp",
            "conversion_applied": False,
            "unit_provenance": "MassBank COLLISION_ENERGY literal",
            "cross_vendor_or_cross_model_equivalence_assumed": False,
        }
    match = _FIXED_CE.fullmatch(raw)
    if match:
        value = float(match.group(1))
        raw_unit = match.group(2).lower()
        unit = {"ev": "eV", "v": "V", "nce": "NCE"}[raw_unit]
        return {
            "raw": raw,
            "value": value,
            "minimum": value,
            "maximum": value,
            "unit": unit,
            "mode": "fixed",
            "conversion_applied": False,
            "unit_provenance": "MassBank COLLISION_ENERGY literal",
            "cross_vendor_or_cross_model_equivalence_assumed": False,
        }
    raise RecordReject("collision_energy_unresolved", raw)


def _charge_from_adduct(adduct: str) -> int:
    match = re.search(r"\](\d*)([+-])$", adduct)
    if not match:
        raise RecordReject("adduct_charge_unresolved", adduct)
    magnitude = int(match.group(1) or "1")
    return magnitude if match.group(2) == "+" else -magnitude


def _first_number(value: Any) -> float:
    match = re.search(_NUMBER, str(value or ""))
    if not match:
        raise RecordReject("numeric_value_missing", str(value or ""))
    number = float(match.group(0))
    if not math.isfinite(number):
        raise RecordReject("numeric_value_nonfinite", str(value or ""))
    return number


def _precursor_mz(
    reported_value: Any,
    *,
    neutral_exact_mass: Any,
    adduct: str,
    source: Mapping[str, Any],
) -> tuple[float, dict[str, Any]]:
    if str(reported_value or "").strip():
        return _first_number(reported_value), {
            "basis": "MassBank focused-ion precursor m/z",
            "computed": False,
        }
    if not source.get("allow_precursor_mz_from_neutral_exact_mass"):
        raise RecordReject("precursor_mz_missing")
    try:
        neutral_mass = float(neutral_exact_mass)
    except (TypeError, ValueError):
        raise RecordReject("neutral_exact_mass_missing_for_precursor_inference")
    proton_mass_da = 1.007276466621
    offsets = {"[M+H]+": proton_mass_da, "[M-H]-": -proton_mass_da}
    if adduct not in offsets:
        raise RecordReject("precursor_inference_adduct_unsupported", adduct)
    return neutral_mass + offsets[adduct], {
        "basis": "computed from MassBank neutral exact mass and reported adduct",
        "computed": True,
        "neutral_exact_mass": neutral_mass,
        "mass_offset_da": offsets[adduct],
        "constant": "proton mass in unified atomic mass units",
        "constant_value_da": proton_mass_da,
    }


def parse_scan_range(raw: Any) -> tuple[float, float] | None:
    """Parse a positive m/z scan interval without treating ``-`` as a sign.

    MassBank commonly serializes scan limits as ``55-1254``.  Applying a
    signed-number regular expression to that literal yields ``55`` and
    ``-1254``.  This parser instead requires an explicit range separator and
    accepts the hyphen-minus, Unicode en/em dashes, or the word ``to``.
    Ambiguous, non-positive, or zero-width intervals fail closed.
    """
    text = str(raw or "").strip()
    if not text:
        return None
    if re.search(r"(?<![\d.])-\s*\d", text):
        return None
    unsigned = r"(?:\d+(?:\.\d*)?|\.\d+)"
    match = re.search(
        rf"(?<![\d.])({unsigned})\s*(?:-|–|—|\bto\b)\s*({unsigned})(?![\d.])",
        text,
        flags=re.IGNORECASE,
    )
    if not match:
        return None
    lower, upper = float(match.group(1)), float(match.group(2))
    if not all(math.isfinite(value) and value > 0 for value in (lower, upper)):
        return None
    minimum, maximum = sorted((lower, upper))
    if minimum == maximum:
        return None
    return minimum, maximum


def _mass_range(raw: Any) -> tuple[float, float] | None:
    """Backward-compatible private alias used by the V2 corpus builder."""
    return parse_scan_range(raw)


def _chemont(raw: str) -> dict[str, Any]:
    parts = [part.strip() for part in str(raw or "").split(";") if part.strip()]
    if len(parts) < 2 or not parts[0].upper().startswith("CHEMONTID:"):
        raise RecordReject("chemont_leaf_class_missing", raw)
    return {
        "ontology": "ChemOnt",
        "ontology_id": parts[0],
        "path": parts[1:],
        "leaf_class": parts[-1],
        "raw": raw,
    }


def _normalize_peaks(record: Mapping[str, Any]) -> list[list[float]]:
    values = record.get("peak", {}).get("peak", {}).get("values", [])
    peaks: list[tuple[float, float]] = []
    for item in values or []:
        try:
            mz = float(item["mz"])
            intensity = float(item["intensity"])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
            peaks.append((mz, intensity))
    if not peaks:
        raise RecordReject("positive_peak_list_missing")
    base = max(intensity for _, intensity in peaks)
    return [[round(mz, 6), round(intensity / base, 8)] for mz, intensity in sorted(peaks)]


def _structure_checks(inchi: str, inchikey: str, formula: str) -> dict[str, Any]:
    result = {
        "inchi_reported": bool(inchi and inchi.startswith("InChI=")),
        "inchikey_matches_reported_inchi": None,
        "formula_matches_reported_inchi": None,
        "validation_library": "RDKit when parseable",
    }
    if not result["inchi_reported"]:
        return result
    try:
        from rdkit.Chem import rdMolDescriptors
        from rdkit.Chem.inchi import InchiToInchiKey, MolFromInchi

        derived_key = InchiToInchiKey(inchi)
        molecule = MolFromInchi(inchi)
        result["inchikey_matches_reported_inchi"] = derived_key == inchikey
        if molecule is not None:
            derived_formula = normalize_formula(rdMolDescriptors.CalcMolFormula(molecule))
            result["formula_matches_reported_inchi"] = derived_formula == formula
            result["derived_formula"] = derived_formula
        if derived_key and derived_key != inchikey:
            raise RecordReject("inchikey_inchi_mismatch", f"{inchikey}!={derived_key}")
    except RecordReject:
        raise
    except Exception as exc:
        result["validation_error"] = type(exc).__name__
    return result


def _domain_for_record(
    instrument: str,
    activation_raw: str,
    ce: Mapping[str, Any],
    source: Mapping[str, Any],
    domains: Mapping[str, Mapping[str, Any]],
) -> Mapping[str, Any]:
    candidates = []
    for domain_id in source.get("domain_ids", []):
        domain = domains.get(str(domain_id))
        if not domain:
            raise CorpusBuildError(f"source references unknown domain: {domain_id}")
        if instrument not in domain.get("models", []):
            continue
        if activation_raw not in domain.get("activation_raw", []):
            continue
        if ce["unit"] not in domain.get("allowed_ce_units", []):
            continue
        if ce["mode"] not in domain.get("allowed_ce_modes", []):
            continue
        candidates.append(domain)
    if not candidates:
        raise RecordReject(
            "exact_instrument_activation_ce_domain_unresolved",
            f"instrument={instrument!r};activation={activation_raw!r};"
            f"ce={ce['unit']}/{ce['mode']}",
        )
    if len(candidates) > 1:
        raise RecordReject("exact_instrument_domain_ambiguous", instrument)
    return candidates[0]


def _convert_record(
    record: Mapping[str, Any],
    *,
    source: Mapping[str, Any],
    domains: Mapping[str, Mapping[str, Any]],
    repository: Mapping[str, Any],
    source_file_sha256: str,
    accepted_licenses: set[str],
) -> dict[str, Any]:
    accession = str(record.get("accession") or "").strip()
    if not accession:
        raise RecordReject("accession_missing")
    if str(record.get("license") or "") not in accepted_licenses:
        raise RecordReject("license_not_allowed", str(record.get("license") or ""))
    if record.get("deprecated") not in ({}, None):
        raise RecordReject("deprecated_record")

    compound = record.get("compound") or {}
    acquisition = record.get("acquisition") or {}
    acquisition_ms = acquisition.get("mass_spectrometry") or {}
    tags = _tags(acquisition_ms.get("subtags"))
    focused = _tags((record.get("mass_spectrometry") or {}).get("focused_ion"))
    links = _links(compound.get("link"))
    name, name_index = stable_compound_name(compound.get("names"))
    inchikey = _first_link(links, "INCHIKEY")
    if not _FULL_INCHIKEY.fullmatch(inchikey):
        raise RecordReject("full_inchikey_missing", inchikey)
    formula = normalize_formula(compound.get("formula"))
    chemont = _chemont(_first_link(links, "ChemOnt"))
    inchi = str(compound.get("inchi") or "")
    structure_validation = _structure_checks(inchi, inchikey, formula)

    activation_raw = str(tags.get("FRAGMENTATION_MODE") or "").strip()
    ce = parse_collision_energy(tags.get("COLLISION_ENERGY"))
    instrument = str(acquisition.get("instrument") or "").strip()
    domain = _domain_for_record(instrument, activation_raw, ce, source, domains)
    polarity_raw = str(acquisition_ms.get("ion_mode") or "").strip().upper()
    if polarity_raw not in {"POSITIVE", "NEGATIVE"}:
        raise RecordReject("polarity_unresolved", polarity_raw)
    polarity = polarity_raw.lower()
    adduct = str(focused.get("PRECURSOR_TYPE") or focused.get("ION_TYPE") or "").strip()
    if not adduct:
        raise RecordReject("precursor_adduct_missing")
    charge = _charge_from_adduct(adduct)
    if (charge > 0) != (polarity == "positive"):
        raise RecordReject("adduct_polarity_inconsistent", f"{adduct}/{polarity}")
    precursor_raw = focused.get("PRECURSOR_M/Z") or focused.get("FULL_SCAN_FRAGMENT_ION_PEAK")
    precursor_mz, precursor_mz_provenance = _precursor_mz(
        precursor_raw,
        neutral_exact_mass=compound.get("mass"),
        adduct=adduct,
        source=source,
    )
    peaks = _normalize_peaks(record)
    peak_hash = exact_peak_hash(peaks)
    reported_range = _mass_range(tags.get("MASS_RANGE_M/Z"))
    peak_min = min(peak[0] for peak in peaks)
    peak_max = max(peak[0] for peak in peaks)
    record_hash = sha256_bytes(canonical_bytes(record))
    source_url = f"{repository['api_root']}/records/{accession}"
    return {
        "schema_version": "2.0.0",
        "spectrum_id": f"MassBank:{accession}",
        "evidence_stage": "raw_msms_spectrum",
        "experimental_status": "public_repository_experimental_spectrum",
        "reference_identity_status": "repository_reported_reference_identity",
        "query_confirmation_status": "not_applicable",
        "source": {
            "repository": repository["name"],
            "repository_release": repository["release"],
            "source_id": source["source_id"],
            "contributor": source["contributor"],
            "laboratory_id": source["laboratory_id"],
            "research_id": source["research_id"],
            "accession": accession,
            "record_url": source_url,
            "query_url": source["query_url"],
            "license_label": str(record.get("license") or ""),
            "source_file_sha256": source_file_sha256,
            "record_canonical_json_sha256": record_hash,
            "record_hash_scope": "canonical parsed API record",
            "record_created": (record.get("date") or {}).get("created"),
            "record_updated": (record.get("date") or {}).get("updated"),
        },
        "compound": {
            "display_name": name,
            "name_provenance": {
                "source": "MassBank compound.names",
                "list_index": name_index,
                "accession": accession,
            },
            "database_identifiers": links,
            "inchi": inchi or None,
            "inchikey": inchikey,
            "connectivity": inchikey[:14],
            "formula": formula,
            "smiles": compound.get("smiles"),
            "structure_validation": structure_validation,
            "chemical_class": chemont,
        },
        "acquisition": {
            "instrument_vendor": domain["vendor"],
            "instrument_family": domain["family"],
            "instrument_model": instrument,
            "instrument_type": acquisition.get("instrument_type"),
            "analyzer": domain["analyzer"],
            "activation_raw": activation_raw,
            "activation": domain["activation"],
            "activation_evidence": domain.get(
                "activation_evidence", "MassBank contributor text"
            ),
            "psi_beam_type_cid_confirmed": domain.get(
                "psi_beam_type_cid_confirmed"
            ),
            "psi_hcd_cv_confirmed": domain.get("psi_hcd_cv_confirmed"),
            "activation_architecture_interpretation": domain.get(
                "architecture_interpretation"
            ),
            "collision_cell_region": domain["collision_cell_region"],
            "instrument_domain_id": domain["domain_id"],
            "collision_energy": ce,
            "adduct": adduct,
            "polarity": polarity,
            "charge": charge,
            "precursor_mz": precursor_mz,
            "precursor_mz_provenance": precursor_mz_provenance,
            "reported_mass_range": (
                {"minimum": reported_range[0], "maximum": reported_range[1]}
                if reported_range
                else None
            ),
        },
        "independence": {
            "exact_peak_hash": peak_hash,
            "repository_splash": (record.get("peak") or {}).get("splash"),
            "mirror_cluster_id": f"exact:{peak_hash}",
            "two_dimensional_identity_id": inchikey[:14],
            "formula_family_id": formula,
            "source_origin_id": source["source_id"],
            "laboratory_id": source["laboratory_id"],
            "research_id": source["research_id"],
            "independent_fragmentation_group_id": sha256_bytes(
                canonical_bytes(
                    {
                        "source": source["source_id"],
                        "identity": inchikey[:14],
                        "domain": domain["domain_id"],
                        "ce": ce,
                        "adduct": adduct,
                        "polarity": polarity,
                        "peak_hash": peak_hash,
                    }
                )
            )[:24],
        },
        "quality": {
            "centroid_peak_list": True,
            "positive_finite_peaks": True,
            "peak_count": len(peaks),
            "at_least_three_peaks": len(peaks) >= 3,
            "collision_energy_verified": True,
            "cross_ce_unit_conversion_applied": False,
            "instrument_domain_verified": True,
            "full_inchikey_present": True,
            "chemont_leaf_present": True,
            "reported_mass_range_present": reported_range is not None,
            "reported_mass_range_contains_precursor": (
                reported_range[0] <= precursor_mz <= reported_range[1]
                if reported_range
                else None
            ),
            "observed_peak_minimum_mz": peak_min,
            "observed_peak_maximum_mz": peak_max,
            "precursor_isolation_window_reported": False,
            "absence_evidence_limit": (
                "reported_mass_range_available_but_isolation_window_unreported"
                if reported_range
                else "mass_range_and_isolation_window_unreported"
            ),
        },
        "peaks": peaks,
    }


def _download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".partial")
    request = urllib.request.Request(
        url, headers={"User-Agent": "Mssystem-V15-feature-fragment-benchmark/2.0"}
    )
    with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    temporary.replace(target)


def _write_jsonl(handle: TextIO, value: Mapping[str, Any]) -> bytes:
    line = canonical_bytes(value) + b"\n"
    handle.write(line.decode("utf-8"))
    return line


def build_massbank_corpus(
    manifest_path: str | Path,
    *,
    cache_dir: str | Path,
    output_jsonl: str | Path,
    rejects_jsonl: str | Path,
    download_missing: bool = False,
    max_records: int | None = None,
) -> dict[str, Any]:
    """Build a normalized corpus and return a structured construction receipt."""
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    repository = manifest.get("repository") or {}
    sources = manifest.get("sources") or []
    domain_list = manifest.get("instrument_domains") or []
    if not sources or not domain_list:
        raise CorpusBuildError("manifest requires sources and instrument_domains")
    domains = {str(domain["domain_id"]): domain for domain in domain_list}
    if len(domains) != len(domain_list):
        raise CorpusBuildError("duplicate instrument domain id")
    accepted_licenses = set(
        (manifest.get("license_policy") or {}).get("accepted_labels_for_local_analysis")
        or []
    )
    if not accepted_licenses:
        raise CorpusBuildError("manifest has no accepted local-analysis license labels")

    cache_dir = Path(cache_dir)
    output_jsonl = Path(output_jsonl)
    rejects_jsonl = Path(rejects_jsonl)
    output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    rejects_jsonl.parent.mkdir(parents=True, exist_ok=True)
    corpus_digest = hashlib.sha256()
    reject_digest = hashlib.sha256()
    accepted = 0
    rejected = 0
    raw_total = 0
    source_counts: dict[str, dict[str, int]] = {}
    rejection_reasons: Counter[str] = Counter()
    domain_counts: Counter[str] = Counter()
    license_counts: Counter[str] = Counter()
    precursor_mz_basis_counts: Counter[str] = Counter()
    structure_check_counts: Counter[str] = Counter()
    reported_mass_range_counts: Counter[str] = Counter()
    identity_sets: dict[str, set[str]] = defaultdict(set)
    formula_sets: dict[str, set[str]] = defaultdict(set)
    source_file_hashes: dict[str, str] = {}

    with output_jsonl.open("w", encoding="utf-8") as corpus_handle, rejects_jsonl.open(
        "w", encoding="utf-8"
    ) as reject_handle:
        for source in sources:
            source_id = str(source["source_id"])
            path = cache_dir / str(source["cache_filename"])
            if not path.exists():
                if not download_missing:
                    raise CorpusBuildError(f"missing explicit cache file: {path}")
                _download(str(source["query_url"]), path)
            source_hash = file_sha256(path)
            source_file_hashes[source_id] = source_hash
            source_raw = 0
            source_accepted = 0
            source_rejected = 0
            for record in iter_json_array(path):
                if max_records is not None and raw_total >= max_records:
                    break
                source_raw += 1
                raw_total += 1
                accession = str(record.get("accession") or "") if isinstance(record, Mapping) else ""
                try:
                    if not isinstance(record, Mapping):
                        raise RecordReject("source_record_not_object")
                    normalized = _convert_record(
                        record,
                        source=source,
                        domains=domains,
                        repository=repository,
                        source_file_sha256=source_hash,
                        accepted_licenses=accepted_licenses,
                    )
                except RecordReject as exc:
                    rejected += 1
                    source_rejected += 1
                    rejection_reasons[exc.code] += 1
                    line = _write_jsonl(
                        reject_handle,
                        {
                            "source_id": source_id,
                            "accession": accession,
                            "reason_code": exc.code,
                            "detail": exc.detail,
                        },
                    )
                    reject_digest.update(line)
                    continue
                line = _write_jsonl(corpus_handle, normalized)
                corpus_digest.update(line)
                accepted += 1
                source_accepted += 1
                domain_id = normalized["acquisition"]["instrument_domain_id"]
                domain_counts[domain_id] += 1
                license_counts[normalized["source"]["license_label"]] += 1
                precursor_mz_basis_counts[
                    "computed" if normalized["acquisition"]["precursor_mz_provenance"]["computed"] else "reported"
                ] += 1
                structure_validation = normalized["compound"]["structure_validation"]
                for field in (
                    "inchi_reported",
                    "inchikey_matches_reported_inchi",
                    "formula_matches_reported_inchi",
                ):
                    value = structure_validation.get(field)
                    label = "unknown" if value is None else str(bool(value)).lower()
                    structure_check_counts[f"{field}:{label}"] += 1
                reported_mass_range_counts[
                    "reported"
                    if normalized["quality"]["reported_mass_range_present"]
                    else "unreported"
                ] += 1
                identity_sets[source_id].add(normalized["compound"]["connectivity"])
                formula_sets[source_id].add(normalized["compound"]["formula"])
            if max_records is None:
                expected = int(source["expected_record_count"])
                if source_raw != expected:
                    raise CorpusBuildError(
                        f"{source_id}: expected {expected} records, observed {source_raw}"
                    )
            source_counts[source_id] = {
                "raw_records": source_raw,
                "accepted_records": source_accepted,
                "rejected_records": source_rejected,
                "independent_2d_identities": len(identity_sets[source_id]),
                "formula_families": len(formula_sets[source_id]),
            }
            if max_records is not None and raw_total >= max_records:
                break

    all_identities: set[str] = set()
    all_formulas: set[str] = set()
    for values in identity_sets.values():
        all_identities.update(values)
    for values in formula_sets.values():
        all_formulas.update(values)
    return {
        "schema_version": "2.0.0",
        "status": "built",
        "evidence_stage": "raw_msms_spectrum",
        "query_confirmation_status": "not_applicable",
        "manifest_path": str(manifest_path),
        "manifest_sha256": file_sha256(manifest_path),
        "repository_release": repository.get("release"),
        "raw_record_count": raw_total,
        "accepted_record_count": accepted,
        "rejected_record_count": rejected,
        "independent_2d_identity_count": len(all_identities),
        "formula_family_count": len(all_formulas),
        "source_counts": source_counts,
        "domain_counts": dict(sorted(domain_counts.items())),
        "license_counts": dict(sorted(license_counts.items())),
        "precursor_mz_basis_counts": dict(sorted(precursor_mz_basis_counts.items())),
        "structure_check_counts": dict(sorted(structure_check_counts.items())),
        "reported_mass_range_counts": dict(sorted(reported_mass_range_counts.items())),
        "rejection_reasons": dict(sorted(rejection_reasons.items())),
        "source_file_sha256": source_file_hashes,
        "corpus_jsonl": str(output_jsonl),
        "corpus_jsonl_sha256": corpus_digest.hexdigest(),
        "rejects_jsonl": str(rejects_jsonl),
        "rejects_jsonl_sha256": reject_digest.hexdigest(),
        "scientific_scope": "public_reference_spectra_only",
        "runtime_effect": "none_shadow_benchmark_only",
    }
