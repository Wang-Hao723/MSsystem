from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, cast

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sqlalchemy import (
    JSON,
    Column,
    Float,
    ForeignKey,
    Integer,
    String,
    create_engine,
)
from sqlalchemy.orm import Session, declarative_base, relationship

from .config import DEFAULT_DB_PATH, DEFAULT_MAPPING_CSV

Base = declarative_base()


class Compound(Base):
    __tablename__ = "Compounds"

    id = Column(Integer, primary_key=True)
    compound_name = Column(String, index=True)
    molecular_formula = Column(String)
    smiles = Column(String)
    inchi = Column(String)
    inchikey = Column(String, index=True)
    exact_mass = Column(Float)
    classification = Column(String, nullable=True)
    compound_class = Column(String, nullable=True)  # original class from MGF

    mappings = relationship("Mapping", back_populates="compound")


class SpectrumRecord(Base):
    __tablename__ = "Spectra"

    id = Column(Integer, primary_key=True)
    spectrum_id = Column(String, unique=True, index=True)
    source_file = Column(String)
    precursor_mz = Column(Float, index=True)
    retention_time = Column(Float, nullable=True)
    collision_energy = Column(Float, nullable=True)
    collision_energy_raw = Column(String, nullable=True)
    collision_energy_unit = Column(String, nullable=True)
    collision_energy_mode = Column(String, nullable=True)
    collision_energy_min = Column(Float, nullable=True)
    collision_energy_max = Column(Float, nullable=True)
    charge = Column(Integer, nullable=True)
    adduct = Column(String, nullable=True)
    instrument_type = Column(String, nullable=True)
    instrument = Column(String, nullable=True)
    ion_mode = Column(String, nullable=True)
    source_accession = Column(String, nullable=True)
    peaks_json = Column(String)

    mappings = relationship("Mapping", back_populates="spectrum")


class Mapping(Base):
    __tablename__ = "Mapping"

    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey("Compounds.id"), index=True)
    spectrum_id = Column(Integer, ForeignKey("Spectra.id"), index=True)
    score = Column(Float, nullable=True)

    compound = relationship("Compound", back_populates="mappings")
    spectrum = relationship("SpectrumRecord", back_populates="mappings")


@dataclass
class DedupKey:
    inchikey: Optional[str]
    precursor_mz: Optional[float]
    retention_time: Optional[float]


class LibraryBuilder:
    def __init__(
        self,
        db_path: str | Path = DEFAULT_DB_PATH,
        mapping_csv_path: str | Path = DEFAULT_MAPPING_CSV,
        rt_tolerance: float = 0.2,
    ) -> None:
        self.db_path = Path(db_path)
        self.mapping_csv_path = Path(mapping_csv_path)
        self.rt_tolerance = rt_tolerance

    def build_library(
        self,
        cleaned_data: List[Dict[str, Any]],
        reference_data: Optional[List[Dict[str, Any]]] = None,
        overwrite: bool = True,
    ) -> Tuple[Path, Path]:
        merged = list(cleaned_data)
        if reference_data:
            merged.extend(reference_data)

        deduped = self._deduplicate_records(merged)
        enriched = [self._enrich_structure_fields(rec) for rec in deduped]

        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        if overwrite and self.db_path.exists():
            self.db_path.unlink()
        engine = create_engine(f"sqlite:///{self.db_path}")
        Base.metadata.create_all(engine)

        with Session(engine) as session:
            mapping_rows = self._insert_records(session, enriched)
            session.commit()

        self._export_mapping_csv(mapping_rows)
        return self.db_path, self.mapping_csv_path

    def load_library_frames(self) -> Dict[str, pd.DataFrame]:
        engine = create_engine(f"sqlite:///{self.db_path}")
        return {
            "Compounds": pd.read_sql_table("Compounds", engine),
            "Spectra": pd.read_sql_table("Spectra", engine),
            "Mapping": pd.read_sql_table("Mapping", engine),
        }

    def _deduplicate_records(
        self, records: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        selected: Dict[str, Dict[str, Any]] = {}

        for rec in records:
            key = self._build_dedup_key(rec)
            key_str = self._serialize_dedup_key(key)
            score = self._quality_score(rec)

            if key_str not in selected or score > self._quality_score(
                selected[key_str]
            ):
                selected[key_str] = rec

        return list(selected.values())

    def _build_dedup_key(self, rec: Dict[str, Any]) -> DedupKey:
        return DedupKey(
            inchikey=(rec.get("inchikey") or rec.get("InChIKey")),
            precursor_mz=rec.get("precursor_mz"),
            retention_time=rec.get("retention_time"),
        )

    def _serialize_dedup_key(self, key: DedupKey) -> str:
        if key.inchikey:
            return f"inchikey::{key.inchikey}"

        if key.precursor_mz is None:
            return "unknown::none"

        rt_bucket = "none"
        if key.retention_time is not None:
            rt_bucket = f"{round(key.retention_time / self.rt_tolerance) * self.rt_tolerance:.2f}"

        return f"mzrt::{key.precursor_mz:.4f}::{rt_bucket}"

    @staticmethod
    def _quality_score(rec: Dict[str, Any]) -> float:
        peaks = rec.get("peaks", [])
        if not peaks:
            return 0.0
        intensities = [float(p[1]) for p in peaks if len(p) == 2]
        return sum(intensities)

    def _enrich_structure_fields(self, rec: Dict[str, Any]) -> Dict[str, Any]:
        out = dict(rec)
        smiles = out.get("smiles")
        inchi = out.get("inchi")
        mol = None

        if smiles:
            mol = Chem.MolFromSmiles(smiles)
        elif inchi:
            try:
                mol = Chem.MolFromInchi(inchi)
            except Exception:
                mol = None

        if mol is None:
            # Still preserve compound_class even without structure
            return out

        # Use explicit conditional — setdefault() won't overwrite existing None keys
        if not out.get("smiles"):
            out["smiles"] = Chem.MolToSmiles(mol, canonical=True)

        # Generate InChI -> InChIKey
        inchi_val = None
        try:
            inchi_val = Chem.MolToInchi(mol)
            if inchi_val:
                if not out.get("inchi"):
                    out["inchi"] = inchi_val
                if not out.get("inchikey"):
                    out["inchikey"] = Chem.InchiToInchiKey(inchi_val)
        except Exception:
            # Fallback: derive InChIKey from canonical SMILES
            try:
                canonical_smiles = out.get("smiles")
                if canonical_smiles:
                    tmp_mol = Chem.MolFromSmiles(canonical_smiles)
                    if tmp_mol:
                        tmp_inchi = Chem.MolToInchi(tmp_mol)
                        if tmp_inchi and not out.get("inchikey"):
                            out["inchikey"] = Chem.InchiToInchiKey(tmp_inchi)
            except Exception:
                pass

        if not out.get("formula"):
            out["formula"] = rdMolDescriptors.CalcMolFormula(mol)
        if not out.get("exact_mass"):
            out["exact_mass"] = float(rdMolDescriptors.CalcExactMolWt(mol))
        return out

    def _insert_records(
        self,
        session: Session,
        records: Iterable[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        compounds_cache: Dict[str, int] = {}
        mapping_rows: List[Dict[str, Any]] = []

        for rec in records:
            compound_id = self._upsert_compound(session, compounds_cache, rec)

            spectrum = SpectrumRecord(
                spectrum_id=str(rec.get("spectrum_id")),
                source_file=rec.get("source_file"),
                precursor_mz=self._to_float(rec.get("precursor_mz")),
                retention_time=self._to_float(rec.get("retention_time")),
                collision_energy=self._to_float(rec.get("collision_energy")),
                collision_energy_raw=rec.get("collision_energy_raw"),
                collision_energy_unit=rec.get("collision_energy_unit"),
                collision_energy_mode=rec.get("collision_energy_mode"),
                collision_energy_min=self._to_float(rec.get("collision_energy_min")),
                collision_energy_max=self._to_float(rec.get("collision_energy_max")),
                charge=self._to_int(rec.get("charge")),
                adduct=rec.get("adduct"),
                instrument_type=rec.get("instrument_type"),
                instrument=rec.get("instrument"),
                ion_mode=rec.get("ion_mode"),
                source_accession=rec.get("source_accession"),
                peaks_json=json.dumps(rec.get("peaks", []), ensure_ascii=False),
            )
            session.add(spectrum)
            session.flush()

            mapping = Mapping(compound_id=compound_id, spectrum_id=spectrum.id)
            session.add(mapping)
            session.flush()

            mapping_rows.append(
                {
                    "compound_id": compound_id,
                    "spectrum_row_id": spectrum.id,
                    "spectrum_id": spectrum.spectrum_id,
                    "compound_name": rec.get("name"),
                    "inchikey": rec.get("inchikey"),
                    "precursor_mz": rec.get("precursor_mz"),
                }
            )

        return mapping_rows

    def _upsert_compound(
        self,
        session: Session,
        cache: Dict[str, int],
        rec: Dict[str, Any],
    ) -> int:
        key = (
            rec.get("inchikey") or rec.get("name") or f"anon::{rec.get('spectrum_id')}"
        )
        if key in cache:
            return cache[key]

        existing = None
        if rec.get("inchikey"):
            existing = (
                session.query(Compound)
                .filter(Compound.inchikey == rec.get("inchikey"))
                .first()
            )

        if existing:
            cache[key] = cast(int, existing.id)
            return cast(int, existing.id)

        compound = Compound(
            compound_name=rec.get("name"),
            molecular_formula=rec.get("formula"),
            smiles=rec.get("smiles"),
            inchi=rec.get("inchi"),
            inchikey=rec.get("inchikey"),
            exact_mass=self._to_float(rec.get("exact_mass")),
            classification=rec.get("classification") or rec.get("compound_class"),
            compound_class=rec.get("compound_class"),
        )
        session.add(compound)
        session.flush()

        cache[key] = cast(int, compound.id)
        return cast(int, compound.id)

    def _export_mapping_csv(self, mapping_rows: List[Dict[str, Any]]) -> None:
        self.mapping_csv_path.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = [
            "compound_id",
            "spectrum_row_id",
            "spectrum_id",
            "compound_name",
            "inchikey",
            "precursor_mz",
        ]
        with self.mapping_csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(mapping_rows)

    @staticmethod
    def _to_float(value: Any) -> Optional[float]:
        if value is None:
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _to_int(value: Any) -> Optional[int]:
        if value is None:
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            return None
