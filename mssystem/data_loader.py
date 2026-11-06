from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
from matchms.importing import load_from_mgf, load_from_msp
from pyteomics import mzxml

from .config import METADATA_ALIASES, SUPPORTED_INPUT_EXTS
from utils.collision_energy import parse_collision_energy

logger = logging.getLogger(__name__)


def _compact_metadata_key(key: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(key or "").lower())


_RAW_METADATA_ALIAS_LOOKUP = {
    _compact_metadata_key(alias): field
    for field, aliases in METADATA_ALIASES.items()
    for alias in aliases
}

# Patterns for extracting structured data from MGF COMMENT fields
_RE_COMMENT_SMILES = re.compile(r"SMILES[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_INCHI = re.compile(r"INCHI[:=\s]+(InChI[^\s;]+)", re.IGNORECASE)
_RE_COMMENT_INCHIKEY = re.compile(r"INCHIKEY[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_FORMULA = re.compile(r"FORMULA[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_NAME = re.compile(
    r"NAME[:=\s]+([^\s;]+)", re.IGNORECASE
)
_RE_COMMENT_EXACT_MASS = re.compile(r"EXACT_MASS[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_ADDUCT = re.compile(r"\[([M][^\]]+)\]", re.IGNORECASE)
_RE_COMMENT_ION_MODE = re.compile(r"ION_MODE[:=\s]+([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_USER_CLASS = re.compile(r"USER_CLASS[:=\s]*([^\s;]+)", re.IGNORECASE)
_RE_COMMENT_COMPOUND_CLASS = re.compile(
    r"COMPOUND_CLASS[:=\s]+([^\s;]+)", re.IGNORECASE
)
_RE_COMMENT_COLLISION_ENERGY = re.compile(
    r"(?:COLLISION[_\s]?ENERGY|ACTIVATION[_\s]?ENERGY|HCD[_\s]?ENERGY|CE)"
    r"\s*[:=]\s*"
    r"((?:\d+(?:\.\d*)?|\.\d+)"
    r"(?:\s*(?:eV|NCE|V|%|percent))?)",
    re.IGNORECASE,
)
_RE_COMMENT_COLLISION_ENERGY_UNIT = re.compile(
    r"(?:COLLISION[_\s]?ENERGY[_\s]?UNIT|CE[_\s]?UNIT|ENERGY[_\s]?UNIT)"
    r"\s*[:=]\s*([A-Za-z%]+)",
    re.IGNORECASE,
)


def _parse_comment_field(comment: str) -> Dict[str, Optional[str]]:
    """Extract SMILES, InChI, InChIKey, Formula, Name, ExactMass, Adduct,
    ION_MODE, USER_CLASS from an MGF COMMENT block."""
    if not comment:
        return {}

    # Join multi-line comments
    block = " ".join(comment.split())

    def _first(pat: re.Pattern, text: str) -> Optional[str]:
        m = pat.search(text)
        if not m:
            return None
        return m.group(1).strip().strip('"').strip("'")

    result: Dict[str, Optional[str]] = {}
    val = _first(_RE_COMMENT_SMILES, block)
    if val:
        result["smiles"] = val
    val = _first(_RE_COMMENT_INCHI, block)
    if val:
        result["inchi"] = val
    val = _first(_RE_COMMENT_INCHIKEY, block)
    if val:
        result["inchikey"] = val
    val = _first(_RE_COMMENT_FORMULA, block)
    if val:
        result["formula"] = val
    val = _first(_RE_COMMENT_NAME, block)
    if val:
        result["name"] = val
    val = _first(_RE_COMMENT_EXACT_MASS, block)
    if val:
        result["exact_mass"] = val
    val = _first(_RE_COMMENT_ADDUCT, block)
    if val:
        result["adduct"] = val
    val = _first(_RE_COMMENT_ION_MODE, block)
    if val:
        result["ion_mode"] = val
    val = _first(_RE_COMMENT_USER_CLASS, block)
    if val:
        result["compound_class"] = val
    val = _first(_RE_COMMENT_COMPOUND_CLASS, block)
    if val:
        result["compound_class"] = val  # prefer compound_class if present
    val = _first(_RE_COMMENT_COLLISION_ENERGY, block)
    if val:
        result["collision_energy"] = val
    val = _first(_RE_COMMENT_COLLISION_ENERGY_UNIT, block)
    if val:
        result["collision_energy_unit"] = val

    return result


class DataCleaner:
    def __init__(
        self,
        top_n_peaks: int = 50,
        relative_intensity_threshold: float = 0.01,
        precursor_fragment_gap: float = 1.5,
        normalize_intensities: bool = True,
    ) -> None:
        self.top_n_peaks = top_n_peaks
        self.relative_intensity_threshold = relative_intensity_threshold
        self.precursor_fragment_gap = precursor_fragment_gap
        self.normalize_intensities = normalize_intensities

    def clean_directory(self, input_dir: str | Path) -> List[Dict[str, Any]]:
        input_dir = Path(input_dir)
        if not input_dir.exists():
            raise FileNotFoundError(f"Input directory not found: {input_dir}")

        files = [
            p
            for p in input_dir.rglob("*")
            if p.is_file() and p.suffix.lower() in SUPPORTED_INPUT_EXTS
        ]
        return self.clean_files(files)

    def clean_files(self, files: Iterable[str | Path]) -> List[Dict[str, Any]]:
        cleaned: List[Dict[str, Any]] = []
        for file_path in files:
            path = Path(file_path)
            if path.suffix.lower() not in SUPPORTED_INPUT_EXTS:
                logger.warning("Skipping unsupported format: %s", path)
                continue

            try:
                if path.suffix.lower() == ".mgf":
                    cleaned.extend(self._parse_matchms_file(path, "mgf"))
                elif path.suffix.lower() == ".msp":
                    cleaned.extend(self._parse_matchms_file(path, "msp"))
                elif path.suffix.lower() == ".mzxml":
                    cleaned.extend(self._parse_mzxml(path))
            except Exception as exc:  # noqa: BLE001
                logger.exception("Failed parsing %s: %s", path, exc)
        return cleaned

    def _parse_matchms_file(
        self, file_path: Path, file_type: str
    ) -> List[Dict[str, Any]]:
        if file_type == "mgf":
            return self._parse_mgf_raw(file_path)
        # MSP still uses matchms
        loader = load_from_msp
        return self._parse_msp_matchms(file_path, loader)

    def _parse_mgf_raw(self, file_path: Path) -> List[Dict[str, Any]]:
        """Parse MGF file by reading raw lines directly, so each COMMENT=
        line is captured individually rather than being collapsed by matchms."""
        records: List[Dict[str, Any]] = []
        raw_text = file_path.read_text(encoding="utf-8", errors="replace")

        i = 0
        for block in raw_text.split("BEGIN IONS"):
            if not block.strip() or block.strip() == "END IONS":
                continue

            fields: Dict[str, Any] = {}
            peaks_mz: List[float] = []
            peaks_int: List[float] = []

            for line in block.split("\n"):
                line = line.rstrip()
                if not line:
                    continue
                if "=" in line:
                    raw_key, raw_value = line.split("=", 1)
                    normalized_key = self._normalize_key(raw_key)
                    alias_field = _RAW_METADATA_ALIAS_LOOKUP.get(normalized_key)
                    if alias_field == "precursor_mz":
                        try:
                            fields["precursor_mz"] = float(
                                raw_value.strip().split()[0]
                            )
                        except (ValueError, IndexError):
                            pass
                        continue
                    if alias_field == "retention_time":
                        try:
                            fields["retention_time"] = float(raw_value.strip())
                        except ValueError:
                            pass
                        continue
                    if alias_field == "charge":
                        fields["charge_raw"] = raw_value.strip()
                        continue
                    if normalized_key == "comment":
                        parsed = _parse_comment_field(raw_value.strip())
                        for key, value in parsed.items():
                            if value is not None and key not in fields:
                                fields[key] = value
                        continue
                    if alias_field:
                        fields[alias_field] = raw_value.strip()
                        continue
                if line.startswith("TITLE="):
                    fields["title"] = line[6:].strip()
                elif line.startswith("PEPMASS="):
                    try:
                        fields["precursor_mz"] = float(line[8:].split()[0])
                    except (ValueError, IndexError):
                        pass
                elif line.startswith("CHARGE="):
                    fields["charge_raw"] = line[7:].strip()
                elif line.startswith("RTINSECONDS="):
                    try:
                        fields["retention_time"] = float(line[13:].strip())
                    except ValueError:
                        pass
                elif line.startswith("COMMENT="):
                    comment_val = line[8:].strip()
                    # Parse key=value or key:value pairs from the comment
                    parsed = _parse_comment_field(comment_val)
                    for k, v in parsed.items():
                        if v is not None and k not in fields:
                            fields[k] = v
                elif line.startswith("ION_MODE="):
                    fields["ion_mode"] = line[10:].strip()
                elif line.upper().startswith("COLLISION_ENERGY="):
                    fields["collision_energy"] = line.split("=", 1)[1].strip()
                elif line.upper().startswith("INSTRUMENT_TYPE="):
                    fields["instrument_type"] = line.split("=", 1)[1].strip()
                elif line.upper().startswith("INSTRUMENT="):
                    fields["instrument"] = line.split("=", 1)[1].strip()
                elif line.startswith("NAME="):
                    fields["name"] = line[5:].strip()
                elif line.startswith("SMILES="):
                    fields["smiles"] = line[7:].strip()
                elif line.startswith("INCHI="):
                    fields["inchi"] = line[6:].strip()
                elif line.startswith("INCHIAUX="):
                    fields["inchikey"] = line[9:].strip()
                elif line.startswith("EXACTMASS="):
                    fields["exact_mass"] = line[10:].strip()
                elif line.startswith("FORMULA="):
                    fields["formula"] = line[8:].strip()
                elif line.startswith("ADDUCT="):
                    fields["adduct"] = line[7:].strip()
                elif line.startswith("IONMODE="):
                    fields["ion_mode"] = line[8:].strip()
                elif line.startswith("DESCRIPTION="):
                    fields["description"] = line[12:].strip()
                elif line.startswith("COMPOUND_NAME="):
                    fields["name"] = line[14:].strip()
                elif "\t" in line or " " in line.strip():
                    # Peak line: m/z[tab/space]intensity
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            mz = float(parts[0])
                            inten = float(parts[1])
                            peaks_mz.append(mz)
                            peaks_int.append(inten)
                        except ValueError:
                            pass

            precursor_mz = fields.get("precursor_mz")
            if precursor_mz is None or not peaks_mz:
                continue

            peaks = self._clean_peaks(
                np.array(peaks_mz, dtype=float),
                np.array(peaks_int, dtype=float),
                float(precursor_mz),
            )
            if not peaks:
                continue

            standardized = self._standardize_metadata(fields)
            standardized.update(
                {
                    "spectrum_id": f"{file_path.stem}_{i}",
                    "source_file": str(file_path),
                    "precursor_mz": float(precursor_mz),
                    "peaks": peaks,
                }
            )
            # Ensure COMMENT-extracted fields override defaults
            if fields.get("smiles"):
                standardized["smiles"] = fields["smiles"]
            if fields.get("inchi"):
                standardized["inchi"] = fields["inchi"]
            if fields.get("inchikey"):
                standardized["inchikey"] = fields["inchikey"]
            if fields.get("formula"):
                standardized["formula"] = fields["formula"]
            if fields.get("name"):
                standardized["name"] = fields["name"]
            if fields.get("exact_mass"):
                em = self._coerce_float(fields["exact_mass"])
                if em:
                    standardized["exact_mass"] = em
            if fields.get("adduct"):
                standardized["adduct"] = fields["adduct"]
            if fields.get("compound_class"):
                standardized["compound_class"] = fields["compound_class"]
                if not standardized.get("classification"):
                    standardized["classification"] = fields["compound_class"]
            # Adduct from [M+H]+ style in TITLE
            if not standardized.get("adduct") and fields.get("ion_mode"):
                standardized["adduct"] = fields["ion_mode"]

            records.append(standardized)
            i += 1

        return records

    def _parse_msp_matchms(
        self, file_path: Path, loader: Any
    ) -> List[Dict[str, Any]]:
        records: List[Dict[str, Any]] = []

        for i, spec in enumerate(loader(str(file_path))):
            metadata = dict(spec.metadata or {})
            precursor_mz = self._extract_float(metadata, "precursor_mz")

            mz_arr = (
                np.array(spec.peaks.mz, dtype=float)
                if spec.peaks is not None
                else np.array([])
            )
            int_arr = (
                np.array(spec.peaks.intensities, dtype=float)
                if spec.peaks is not None
                else np.array([])
            )

            if precursor_mz is None or mz_arr.size == 0:
                continue

            peaks = self._clean_peaks(mz_arr, int_arr, precursor_mz)
            if not peaks:
                continue

            # Parse COMMENT block for MSP
            comment_str = metadata.get("comment") or ""
            parsed_comment = _parse_comment_field(comment_str)

            full_meta: Dict[str, Any] = dict(metadata)
            for k, v in parsed_comment.items():
                if v is not None and k not in full_meta:
                    full_meta[k] = v

            standardized = self._standardize_metadata(full_meta)
            standardized.update(
                {
                    "spectrum_id": f"{file_path.stem}_{i}",
                    "source_file": str(file_path),
                    "precursor_mz": precursor_mz,
                    "peaks": peaks,
                }
            )
            if parsed_comment.get("smiles"):
                standardized["smiles"] = parsed_comment["smiles"]
            if parsed_comment.get("inchi"):
                standardized["inchi"] = parsed_comment["inchi"]
            if parsed_comment.get("inchikey"):
                standardized["inchikey"] = parsed_comment["inchikey"]
            if parsed_comment.get("formula"):
                standardized["formula"] = parsed_comment["formula"]
            if parsed_comment.get("name"):
                standardized["name"] = parsed_comment["name"]
            if parsed_comment.get("exact_mass"):
                em = self._coerce_float(parsed_comment["exact_mass"])
                if em:
                    standardized["exact_mass"] = em
            if parsed_comment.get("compound_class"):
                standardized["compound_class"] = parsed_comment["compound_class"]
                standardized["classification"] = parsed_comment["compound_class"]

            records.append(standardized)

        return records

    def _parse_mzxml(self, file_path: Path) -> List[Dict[str, Any]]:
        records: List[Dict[str, Any]] = []

        with mzxml.MzXML(str(file_path)) as reader:
            for scan in reader:
                try:
                    ms_level = int(scan.get("msLevel", 1))
                    if ms_level != 2:
                        continue

                    precursor_mz = self._extract_precursor_from_scan(scan)
                    if precursor_mz is None:
                        continue

                    mz_arr = np.array(scan.get("m/z array", []), dtype=float)
                    int_arr = np.array(scan.get("intensity array", []), dtype=float)
                    if mz_arr.size == 0:
                        continue

                    peaks = self._clean_peaks(mz_arr, int_arr, precursor_mz)
                    if not peaks:
                        continue

                    metadata = {
                        "retention_time": scan.get("retentionTime")
                        or scan.get("retention_time")
                        or scan.get("scan start time"),
                        "charge": scan.get("precursorCharge"),
                        "collision_energy": scan.get("collisionEnergy"),
                        "name": f"scan_{scan.get('num', 'unknown')}",
                    }

                    standardized = self._standardize_metadata(metadata)
                    standardized.update(
                        {
                            "spectrum_id": f"{file_path.stem}_{scan.get('num', len(records))}",
                            "source_file": str(file_path),
                            "precursor_mz": precursor_mz,
                            "peaks": peaks,
                        }
                    )
                    records.append(standardized)
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        "Skipping problematic scan in %s: %s", file_path, exc
                    )

        return records

    def _extract_precursor_from_scan(self, scan: Dict[str, Any]) -> Optional[float]:
        precursor = scan.get("precursorMz")
        if isinstance(precursor, list) and precursor:
            if isinstance(precursor[0], dict):
                return self._coerce_float(
                    precursor[0].get("precursorMz") or precursor[0].get("value")
                )
            return self._coerce_float(precursor[0])
        if isinstance(precursor, dict):
            return self._coerce_float(
                precursor.get("precursorMz") or precursor.get("value")
            )
        return self._coerce_float(precursor)

    def _clean_peaks(
        self,
        mz_values: np.ndarray,
        intensity_values: np.ndarray,
        precursor_mz: float,
    ) -> List[List[float]]:
        if mz_values.size == 0 or intensity_values.size == 0:
            return []

        if mz_values.size != intensity_values.size:
            min_len = min(mz_values.size, intensity_values.size)
            mz_values = mz_values[:min_len]
            intensity_values = intensity_values[:min_len]

        max_int = float(np.max(intensity_values)) if intensity_values.size else 0.0
        if max_int <= 0:
            return []

        # Step 1: Normalize intensities to 0-1 BEFORE any filtering.
        # This matches the library peak format and ensures scale-invariant matching.
        if self.normalize_intensities:
            intensity_values = intensity_values / max_int
            max_int = 1.0

        rel_int = intensity_values  # already normalized if normalize=True
        if not self.normalize_intensities:
            rel_int = intensity_values / max_int

        # Step 2: Remove the precursor ion peak (and any satellite peaks near it).
        # Using a wider gap (3 Da) catches [M+H-H2O]+, [M+Na]+, etc.
        # This prevents the dominant precursor peak from dominating the cosine score.
        keep = rel_int >= self.relative_intensity_threshold
        keep &= mz_values < (precursor_mz - self.precursor_fragment_gap)

        mz_filtered = mz_values[keep]
        int_filtered = intensity_values[keep]  # use raw/normalized intensities

        if mz_filtered.size == 0:
            return []

        # Step 3: Sort by intensity, take top N peaks
        top_idx = np.argsort(int_filtered)[::-1][: self.top_n_peaks]
        mz_top = mz_filtered[top_idx]
        int_top = int_filtered[top_idx]

        # Step 4: Sort by m/z for consistent output
        order = np.argsort(mz_top)
        return [[float(mz_top[i]), float(int_top[i])] for i in order]

    def _standardize_metadata(self, metadata: Dict[str, Any]) -> Dict[str, Any]:
        raw_collision_energy = self._extract_alias_value(
            metadata, "collision_energy"
        )
        collision_energy = parse_collision_energy(
            raw_collision_energy,
            self._extract_alias_value(metadata, "collision_energy_unit"),
        )
        std = {
            "name": self._extract_str(metadata, "name") or "Unknown",
            "retention_time": self._extract_float(metadata, "retention_time"),
            "charge": self._extract_charge(metadata),
            "collision_energy": collision_energy["value"],
            "collision_energy_raw": collision_energy["raw"],
            "collision_energy_unit": collision_energy["unit"],
            "collision_energy_mode": collision_energy["mode"],
            "collision_energy_min": collision_energy["minimum"],
            "collision_energy_max": collision_energy["maximum"],
            "adduct": self._extract_str(metadata, "adduct"),
            "smiles": self._extract_str(metadata, "smiles"),
            "inchi": self._extract_str(metadata, "inchi"),
            "inchikey": self._extract_str(metadata, "inchikey"),
            "formula": self._extract_str(metadata, "formula"),
            "instrument": self._extract_str(metadata, "instrument"),
            "instrument_type": self._extract_str(metadata, "instrument_type"),
            "ion_mode": self._extract_str(metadata, "ion_mode"),
            "source_accession": self._extract_str(metadata, "source_accession"),
            "compound_class": self._extract_str(metadata, "compound_class"),
            "exact_mass": self._extract_float(metadata, "exact_mass"),
        }

        # Generate InChIKey from SMILES/InChI if missing
        if not std["inchikey"] and (std["smiles"] or std["inchi"]):
            try:
                from rdkit import Chem
                if std["smiles"]:
                    mol = Chem.MolFromSmiles(std["smiles"])
                elif std["inchi"]:
                    mol = Chem.MolFromInchi(std["inchi"])
                else:
                    mol = None

                if mol:
                    std["inchikey"] = Chem.MolToInchiKey(mol)
            except Exception:
                pass  # Keep inchikey as None if generation fails

        # Generate formula from SMILES if missing
        if not std["formula"] and std["smiles"]:
            try:
                from rdkit import Chem
                from rdkit.Chem import rdMolDescriptors
                mol = Chem.MolFromSmiles(std["smiles"])
                if mol:
                    std["formula"] = rdMolDescriptors.CalcMolFormula(mol)
            except Exception:
                pass

        # Calculate exact mass from formula or SMILES if missing
        if not std["exact_mass"]:
            try:
                from rdkit import Chem
                from rdkit.Chem import rdMolDescriptors
                if std["smiles"]:
                    mol = Chem.MolFromSmiles(std["smiles"])
                    if mol:
                        std["exact_mass"] = rdMolDescriptors.CalcExactMolWt(mol)
                elif std["formula"]:
                    # Try to calculate from formula
                    pass  # This would require additional parsing
            except Exception:
                pass

        return std

    def _extract_alias_value(self, metadata: Dict[str, Any], field: str) -> Any:
        aliases = METADATA_ALIASES.get(field, [field])
        normalized = {self._normalize_key(k): v for k, v in metadata.items()}

        for alias in aliases:
            key_norm = self._normalize_key(alias)
            if key_norm in normalized:
                return normalized[key_norm]

        for key, value in normalized.items():
            if re.search(field.replace("_", ".*"), key, re.IGNORECASE):
                return value
        return None

    @staticmethod
    def _normalize_key(key: str) -> str:
        return _compact_metadata_key(key)

    def _extract_float(self, metadata: Dict[str, Any], field: str) -> Optional[float]:
        raw = self._extract_alias_value(metadata, field)
        return self._coerce_float(raw)

    def _extract_str(self, metadata: Dict[str, Any], field: str) -> Optional[str]:
        raw = self._extract_alias_value(metadata, field)
        if raw is None:
            return None
        if isinstance(raw, (list, tuple)):
            raw = raw[0] if raw else None
        if raw is None:
            return None
        return str(raw).strip() or None

    def _extract_charge(self, metadata: Dict[str, Any]) -> Optional[int]:
        raw = self._extract_alias_value(metadata, "charge")
        if raw is None:
            return None
        if isinstance(raw, (list, tuple)) and raw:
            raw = raw[0]
        if isinstance(raw, str):
            raw = raw.strip().replace("+", "").replace("-", "")
        try:
            return int(float(str(raw)))
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _coerce_float(value: Any) -> Optional[float]:
        if value is None:
            return None
        if isinstance(value, (list, tuple)):
            value = value[0] if value else None
        if value is None:
            return None
        if isinstance(value, (int, float)):
            return float(value)

        value_str = str(value).strip()
        if not value_str:
            return None

        if " " in value_str:
            value_str = value_str.split()[0]

        match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", value_str)
        if not match:
            return None
        try:
            return float(match.group(0))
        except ValueError:
            return None

    @staticmethod
    def export_cleaned_json(
        cleaned_data: List[Dict[str, Any]], out_file: str | Path
    ) -> None:
        out_file = Path(out_file)
        out_file.parent.mkdir(parents=True, exist_ok=True)
        with out_file.open("w", encoding="utf-8") as handle:
            json.dump(cleaned_data, handle, ensure_ascii=False, indent=2)


class BlankSubtractor:
    """Identify and subtract blank/solvent peaks from sample spectra."""

    # Common solvent/adduct peaks (m/z values)
    COMMON_SOLVENT_PEAKS = {
        # Water clusters
        18.0344, 19.0184, 36.0687,
        # Methanol
        32.0339, 33.0184, 50.0523,
        # Acetonitrile
        42.0338, 59.0497,
        # Ethanol
        46.0398, 47.0417, 64.0606,
        # Isopropanol
        60.0554, 61.0598,
        # Acetone
        58.0393, 59.0453,
        # Formic acid
        45.0094, 46.0173,
        # Acetic acid
        60.0211, 61.0290,
        # Ammonium
        18.0344,
        # Sodium adducts (Na+ = 22.9898)
        # Common plasticizers
        149.0239, 149.0915, 149.1213, 207.0848,
        # Column bleed
        73.0473, 75.0576, 89.0596, 103.0752, 117.0909,
        # PEG fragments
        89.0596, 133.0853, 177.1110, 221.1367,
    }

    def __init__(
        self,
        mz_tolerance: float = 0.01,
        min_intensity_ratio: float = 0.8,
        min_occurrence_rate: float = 0.5,
    ):
        """
        Args:
            mz_tolerance: m/z tolerance for peak matching (Da)
            min_intensity_ratio: Min intensity ratio (blank/sample) to consider peak as blank
            min_occurrence_rate: Min fraction of blank spectra containing the peak
        """
        self.mz_tolerance = mz_tolerance
        self.min_intensity_ratio = min_intensity_ratio
        self.min_occurrence_rate = min_occurrence_rate
        self._blank_peaks: List[Tuple[float, float]] = []
        self._blank_peaks_set: set = set()

    def fit_blank(self, blank_file: str | Path) -> "BlankSubtractor":
        """Learn blank peaks from a blank/solvent MGF file."""
        cleaner = DataCleaner(top_n_peaks=100, relative_intensity_threshold=0.001,
                             precursor_fragment_gap=1.5, normalize_intensities=True)
        blank_records = cleaner.clean_files([blank_file])

        if not blank_records:
            logger.warning("No spectra found in blank file: %s", blank_file)
            return self

        # Collect all peaks from blank with occurrence counting
        peak_occurrence: Dict[float, List[float]] = {}

        for rec in blank_records:
            peaks = rec.get("peaks", [])
            for mz, inten in peaks:
                # Round to 4 decimal places for grouping
                mz_rounded = round(mz, 4)
                if mz_rounded not in peak_occurrence:
                    peak_occurrence[mz_rounded] = []
                peak_occurrence[mz_rounded].append(inten)

        # Identify blank peaks: those appearing in >= min_occurrence_rate of spectra
        n_spectra = len(blank_records)
        min_spectra = int(n_spectra * self.min_occurrence_rate)

        for mz, intensities in peak_occurrence.items():
            if len(intensities) >= min_spectra:
                # Use median intensity as representative
                avg_inten = np.median(intensities)
                self._blank_peaks.append((mz, avg_inten))
                self._blank_peaks_set.add(mz)

        # Sort by m/z
        self._blank_peaks.sort(key=lambda x: x[0])

        logger.info(
            "Blank subtraction: %d peaks identified from %d spectra",
            len(self._blank_peaks), n_spectra
        )
        return self

    def fit_from_peak_list(self, peak_list: List[Tuple[float, float]]) -> "BlankSubtractor":
        """Set blank peaks from a provided list of (mz, intensity) tuples."""
        self._blank_peaks = list(peak_list)
        self._blank_peaks.sort(key=lambda x: x[0])
        self._blank_peaks_set = {round(mz, 4) for mz, _ in peak_list}
        logger.info("Blank subtraction: %d peaks set manually", len(self._blank_peaks))
        return self

    def subtract_blank(
        self,
        spectra: List[Dict[str, Any]],
        subtract: bool = True
    ) -> List[Dict[str, Any]]:
        """Remove blank peaks from spectra."""
        if not subtract or not self._blank_peaks:
            return spectra

        result = []
        for rec in spectra:
            new_rec = rec.copy()
            peaks = rec.get("peaks", [])

            # Filter out peaks that match blank
            filtered_peaks = []
            removed_count = 0

            for mz, inten in peaks:
                if self._is_blank_peak(mz, inten):
                    removed_count += 1
                    continue
                filtered_peaks.append((mz, inten))

            # Re-apply cleaning if too many peaks were removed
            if removed_count > 0 and len(filtered_peaks) > 0:
                precursor_mz = rec.get("precursor_mz")
                if precursor_mz:
                    mz_arr = np.array([p[0] for p in filtered_peaks])
                    int_arr = np.array([p[1] for p in filtered_peaks])
                    cleaner = DataCleaner()
                    cleaned = cleaner._clean_peaks(mz_arr, int_arr, precursor_mz)
                    if cleaned:
                        filtered_peaks = cleaned

            if filtered_peaks:
                new_rec["peaks"] = filtered_peaks
                new_rec["blank_peaks_removed"] = removed_count
                result.append(new_rec)

        logger.info(
            "Blank subtraction: processed %d spectra, avg %.1f peaks removed per spectrum",
            len(spectra),
            sum(r.get("blank_peaks_removed", 0) for r in result) / len(result) if result else 0
        )
        return result

    def _is_blank_peak(self, mz: float, intensity: float) -> bool:
        """Check if a peak is a blank peak."""
        # Check against learned blank peaks
        for blank_mz, blank_inten in self._blank_peaks:
            if abs(mz - blank_mz) <= self.mz_tolerance:
                # Check intensity ratio
                if blank_inten > 0 and intensity / blank_inten <= self.min_intensity_ratio:
                    return True
                return True

        # Check against common solvent peaks
        for solvent_mz in self.COMMON_SOLVENT_PEAKS:
            if abs(mz - solvent_mz) <= self.mz_tolerance:
                return True

        return False

    def get_blank_peaks(self) -> List[Tuple[float, float]]:
        """Return the list of blank peaks (mz, intensity)."""
        return self._blank_peaks.copy()

    def get_common_solvent_peaks(self) -> List[float]:
        """Return the list of common solvent m/z values."""
        return sorted(self.COMMON_SOLVENT_PEAKS)
