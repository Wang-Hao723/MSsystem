#!/usr/bin/env python3
"""Train and evaluate the formula-conditional IsoPair-MS isomer ranker.

The task is deliberately narrower than generic molecular-similarity learning:

* positive pair: two experimental spectra assigned to the same 2-D molecular
  connectivity (the first 14 InChIKey characters);
* hard negative: spectra from a different connectivity in the same molecular
  formula and adduct stratum;
* split unit: the complete molecular-formula family, so no isomer family can
  occur in more than one of train/validation/test;
* evaluation unit: leave-one-spectrum-out retrieval among same-formula
  candidate identities.

The spectral backbone is the self-supervised DreaMS checkpoint.  It remains
frozen; only a small symmetric pair head is trained.  This avoids pretending
that a few local standards can support training a 96-million-parameter model
from scratch while still allowing the local library's isomer confusions to
shape the decision boundary.

Generated embeddings and weights are artifacts, not source files.  They are
written below ``models/isopair_ms`` (ignored by Git).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import re
import sqlite3
import sys
import time
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from matching.library_quality import normalize_formula_text
DEFAULT_DB = PROJECT_ROOT.parent.parent / "Mssystem_shared_data" / "local_library.db"
DEFAULT_DREAMS_ROOT = Path("/Users/zizi/Desktop/DreaMS")
DEFAULT_OUT = PROJECT_ROOT / "models" / "isopair_ms"


@dataclass(frozen=True)
class SpectrumRow:
    spectrum_id: int
    accession: str
    formula: str
    connectivity: str
    adduct: str
    precursor_mz: float
    collision_energy: float | None
    collision_energy_unit: str
    instrument_type: str
    dataset_source: str
    peaks: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class QueryCase:
    query_index: int
    correct_connectivity: str
    candidate_indices: tuple[tuple[str, tuple[int, ...]], ...]
    formula: str
    adduct: str


def _stable_fraction(text: str, seed: int) -> float:
    digest = hashlib.sha256(f"{seed}|{text}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") / float(2**64)


def split_formula_key(formula: str) -> str:
    """Normalize bracketed/charged formula spellings before split hashing."""
    text = re.sub(r"\s+", "", str(formula or ""))
    bracketed = re.fullmatch(r"\[([^\]]+)\](?:\^?\d*[+-])?", text)
    if bracketed:
        text = bracketed.group(1)
    else:
        text = re.sub(r"(?:\^?\d*[+-])$", "", text)
    normalized = normalize_formula_text(text)
    return normalized or text.upper()


def formula_split(formula: str, seed: int) -> str:
    """Return a deterministic formula-disjoint 70/15/15 split."""
    value = _stable_fraction(split_formula_key(formula), seed)
    if value < 0.15:
        return "test"
    if value < 0.30:
        return "validation"
    return "train"


def enforce_formula_and_connectivity_disjoint_splits(
    spectra: Sequence[SpectrumRow],
    splits: Sequence[str],
    *,
    seed: int,
) -> list[str]:
    """Assign connected formula/connectivity components to exactly one split.

    Formula normalization alone is insufficient when one connectivity is stored
    under ionic and neutral formula spellings.  Unioning both identifiers makes
    transitive leakage impossible: if either normalized formula or 2-D
    connectivity overlaps, the complete component moves together.
    """
    if len(spectra) != len(splits):
        raise ValueError("spectra and split arrays must have equal length")
    parent = list(range(len(spectra)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    formula_owner: dict[str, int] = {}
    connectivity_owner: dict[str, int] = {}
    for index, row in enumerate(spectra):
        formula = split_formula_key(row.formula)
        connectivity = str(row.connectivity or "").strip().upper()
        if formula:
            if formula in formula_owner:
                union(index, formula_owner[formula])
            else:
                formula_owner[formula] = index
        if connectivity:
            if connectivity in connectivity_owner:
                union(index, connectivity_owner[connectivity])
            else:
                connectivity_owner[connectivity] = index

    members: dict[int, list[int]] = defaultdict(list)
    for index in range(len(spectra)):
        members[find(index)].append(index)
    output = list(splits)
    for indices in members.values():
        component_formulas = sorted({
            split_formula_key(spectra[index].formula) for index in indices
        })
        component_connectivities = sorted({
            str(spectra[index].connectivity or "").strip().upper()
            for index in indices
        })
        component_key = "formula=" + ",".join(component_formulas)
        component_key += "|connectivity=" + ",".join(component_connectivities)
        assigned = formula_split(component_key, seed)
        for index in indices:
            output[index] = assigned
    return output


def split_overlap_audit(
    spectra: Sequence[SpectrumRow],
    splits: Sequence[str],
) -> dict[str, Any]:
    formulas: dict[str, set[str]] = defaultdict(set)
    connectivities: dict[str, set[str]] = defaultdict(set)
    for row, split in zip(spectra, splits):
        formulas[str(split)].add(split_formula_key(row.formula))
        connectivities[str(split)].add(str(row.connectivity or "").upper())
    pairwise = {}
    names = ("train", "validation", "test")
    for offset, left in enumerate(names):
        for right in names[offset + 1:]:
            pairwise[f"{left}_vs_{right}"] = {
                "formula_overlap": len(formulas[left] & formulas[right]),
                "connectivity_overlap": len(
                    connectivities[left] & connectivities[right]
                ),
            }
    return {
        "policy": "normalized_formula_and_connectivity_union_components",
        "pairwise": pairwise,
        "zero_overlap": all(
            values["formula_overlap"] == 0
            and values["connectivity_overlap"] == 0
            for values in pairwise.values()
        ),
    }


def _device(name: str) -> torch.device:
    if name != "auto":
        return torch.device(name)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def _valid_peaks(value: Any, precursor_mz: float) -> tuple[tuple[float, float], ...]:
    try:
        parsed = json.loads(value) if isinstance(value, str) else value
    except (TypeError, ValueError, json.JSONDecodeError):
        return ()
    cleaned: dict[float, float] = {}
    for peak in parsed or []:
        try:
            mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if (
            math.isfinite(mz)
            and math.isfinite(intensity)
            and 0.0 < mz < precursor_mz + 2.0
            and intensity > 0.0
        ):
            cleaned[mz] = cleaned.get(mz, 0.0) + intensity
    return tuple(sorted(cleaned.items()))


def _float_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _connect_readonly(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute("PRAGMA cache_size=-200000")
    return connection


def select_family_strata(
    connection: sqlite3.Connection,
    *,
    seed: int,
    train_families: int,
    validation_families: int,
    test_families: int,
    min_identities: int,
    max_identities: int,
) -> dict[str, list[tuple[str, str]]]:
    """Select one well-replicated adduct stratum per molecular formula."""
    query = """
        WITH strata AS (
            SELECT trim(c.molecular_formula) AS formula,
                   lower(trim(s.adduct)) AS adduct,
                   substr(trim(c.inchikey), 1, 14) AS connectivity,
                   COUNT(*) AS n_spectra,
                   COUNT(DISTINCT coalesce(s.dataset_source, '')) AS n_sources
            FROM Compounds c
            JOIN Mapping m ON m.compound_id = c.id
            JOIN Spectra s ON s.id = m.spectrum_id
            WHERE s.source_tier IN (
                    'experimental', 'curated', 'massbank_supplement'
                  )
              AND length(trim(coalesce(c.inchikey, ''))) = 27
              AND trim(coalesce(c.molecular_formula, '')) <> ''
              AND trim(coalesce(s.adduct, '')) <> ''
              AND s.precursor_mz > 0 AND s.precursor_mz <= 1000
              AND length(coalesce(s.peaks_json, '')) > 10
            GROUP BY 1, 2, 3
        ), family_strata AS (
            SELECT formula, adduct,
                   COUNT(*) AS n_identities,
                   SUM(CASE WHEN n_spectra >= 2 THEN 1 ELSE 0 END)
                       AS n_replicated_identities,
                   SUM(CASE WHEN n_sources >= 2 THEN 1 ELSE 0 END)
                       AS n_multisource_identities,
                   SUM(n_spectra) AS n_spectra
            FROM strata
            GROUP BY formula, adduct
        )
        SELECT formula, adduct, n_identities, n_replicated_identities,
               n_multisource_identities, n_spectra
        FROM family_strata
        WHERE n_identities BETWEEN ? AND ?
          AND n_replicated_identities >= 2
        ORDER BY formula, n_multisource_identities DESC,
                 n_replicated_identities DESC, n_spectra DESC, adduct
    """
    best_per_formula: dict[str, sqlite3.Row] = {}
    for row in connection.execute(query, (min_identities, max_identities)):
        best_per_formula.setdefault(str(row["formula"]), row)

    limits = {
        "train": train_families,
        "validation": validation_families,
        "test": test_families,
    }
    selected: dict[str, list[tuple[str, str]]] = {name: [] for name in limits}
    ordered = sorted(
        best_per_formula.values(),
        key=lambda row: (
            _stable_fraction(f"order|{row['formula']}", seed),
            str(row["formula"]),
        ),
    )
    for row in ordered:
        split = formula_split(str(row["formula"]), seed)
        if len(selected[split]) >= limits[split]:
            continue
        selected[split].append((str(row["formula"]), str(row["adduct"])))
        if all(len(selected[name]) >= limits[name] for name in limits):
            break
    return selected


def fetch_spectra(
    connection: sqlite3.Connection,
    selected: dict[str, list[tuple[str, str]]],
    *,
    max_spectra_per_identity: int,
) -> tuple[list[SpectrumRow], list[str]]:
    """Fetch a bounded deterministic subset with a single database pass."""
    split_by_stratum = {
        (formula, adduct): split
        for split, strata in selected.items()
        for formula, adduct in strata
    }
    connection.execute(
        "CREATE TEMP TABLE selected_isopair_strata ("
        "formula TEXT NOT NULL, adduct TEXT NOT NULL, split TEXT NOT NULL, "
        "PRIMARY KEY (formula, adduct))"
    )
    connection.executemany(
        "INSERT INTO selected_isopair_strata VALUES (?, ?, ?)",
        [
            (formula, adduct, split)
            for (formula, adduct), split in split_by_stratum.items()
        ],
    )
    query = """
        WITH source_ranked AS (
            SELECT s.id AS spectrum_row_id,
                   s.spectrum_id AS accession,
                   trim(c.molecular_formula) AS formula,
                   substr(trim(c.inchikey), 1, 14) AS connectivity,
                   lower(trim(s.adduct)) AS adduct,
                   s.precursor_mz, s.collision_energy,
                   coalesce(s.collision_energy_unit, '') AS collision_energy_unit,
                   coalesce(s.instrument_type, '') AS instrument_type,
                   coalesce(s.dataset_source, '') AS dataset_source,
                   s.peaks_json, chosen.split,
                   ROW_NUMBER() OVER (
                       PARTITION BY trim(c.molecular_formula), lower(trim(s.adduct)),
                                    substr(trim(c.inchikey), 1, 14),
                                    coalesce(s.dataset_source, '')
                       ORDER BY s.spectrum_id, s.id
                   ) AS within_source_rank
            FROM selected_isopair_strata chosen
            JOIN Compounds c ON trim(c.molecular_formula) = chosen.formula
            JOIN Mapping m ON m.compound_id = c.id
            JOIN Spectra s ON s.id = m.spectrum_id
                          AND lower(trim(s.adduct)) = chosen.adduct
            WHERE s.source_tier IN (
                    'experimental', 'curated', 'massbank_supplement'
                  )
              AND length(trim(coalesce(c.inchikey, ''))) = 27
              AND s.precursor_mz > 0 AND s.precursor_mz <= 1000
              AND length(coalesce(s.peaks_json, '')) > 10
        ), eligible AS (
            SELECT *, ROW_NUMBER() OVER (
                       PARTITION BY formula, adduct, connectivity
                       ORDER BY within_source_rank, dataset_source,
                                accession, spectrum_row_id
                   ) AS spectrum_rank
            FROM source_ranked
        )
        SELECT * FROM eligible
        WHERE spectrum_rank <= ?
        ORDER BY split, formula, adduct, connectivity,
                 dataset_source, accession, spectrum_row_id
    """
    spectra: list[SpectrumRow] = []
    splits: list[str] = []
    for row in connection.execute(query, (max_spectra_per_identity,)):
        precursor_mz = float(row["precursor_mz"])
        peaks = _valid_peaks(row["peaks_json"], precursor_mz)
        if len(peaks) < 3:
            continue
        spectra.append(SpectrumRow(
            spectrum_id=int(row["spectrum_row_id"]),
            accession=str(row["accession"] or ""),
            formula=str(row["formula"]),
            connectivity=str(row["connectivity"]),
            adduct=str(row["adduct"]),
            precursor_mz=precursor_mz,
            collision_energy=_float_or_none(row["collision_energy"]),
            collision_energy_unit=str(row["collision_energy_unit"] or ""),
            instrument_type=str(row["instrument_type"] or ""),
            dataset_source=str(row["dataset_source"] or ""),
            peaks=peaks,
        ))
        splits.append(str(row["split"]))
    return spectra, splits


def filter_usable_families(
    spectra: Sequence[SpectrumRow],
    splits: Sequence[str],
) -> tuple[list[SpectrumRow], list[str]]:
    """Keep families containing at least two replicated connectivities."""
    by_family: dict[tuple[str, str], dict[str, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for index, row in enumerate(spectra):
        by_family[(row.formula, row.adduct)][row.connectivity].append(index)
    usable = {
        family for family, identities in by_family.items()
        if sum(len(indices) >= 2 for indices in identities.values()) >= 2
    }
    keep = [
        index for index, row in enumerate(spectra)
        if (row.formula, row.adduct) in usable
    ]
    return [spectra[index] for index in keep], [splits[index] for index in keep]


def _import_dreams(dreams_root: Path):
    root = str(dreams_root.resolve())
    if root not in sys.path:
        sys.path.insert(0, root)
    from dreams.api import PreTrainedModel
    from dreams.models.dreams.dreams import DreaMS as DreaMSModel
    from dreams.utils.data import SpectrumPreprocessor
    from dreams.utils.dformats import DataFormatA
    return PreTrainedModel, DreaMSModel, SpectrumPreprocessor, DataFormatA


def compute_ssl_embeddings(
    spectra: Sequence[SpectrumRow],
    *,
    dreams_root: Path,
    checkpoint: Path,
    batch_size: int,
    device: torch.device,
) -> np.ndarray:
    """Extract frozen precursor-token embeddings from self-supervised DreaMS."""
    PreTrainedModel, DreaMSModel, SpectrumPreprocessor, DataFormatA = (
        _import_dreams(dreams_root)
    )
    wrapper = PreTrainedModel.from_ckpt(
        checkpoint,
        DreaMSModel,
        n_highest_peaks=100,
    )
    model = wrapper.model.eval().to(device)
    preprocessor = SpectrumPreprocessor(DataFormatA(), n_highest_peaks=100)
    output = np.zeros((len(spectra), int(model.d_model)), dtype=np.float32)
    started = time.time()
    for start in range(0, len(spectra), batch_size):
        rows = spectra[start:start + batch_size]
        batch = np.stack([
            preprocessor(
                np.asarray(row.peaks, dtype=np.float32),
                prec_mz=row.precursor_mz,
                high_form=True,
            )
            for row in rows
        ])
        tensor = torch.from_numpy(batch).to(device)
        with torch.inference_mode():
            encoded = model(tensor)[:, 0, :]
            encoded = F.normalize(encoded, dim=-1)
        output[start:start + len(rows)] = encoded.float().cpu().numpy()
        if device.type == "mps":
            torch.mps.empty_cache()
        done = start + len(rows)
        if done == len(rows) or done % (batch_size * 20) == 0 or done == len(spectra):
            elapsed = max(time.time() - started, 1e-6)
            print(
                f"DreaMS embeddings: {done}/{len(spectra)} "
                f"({done / elapsed:.1f} spectra/s)",
                flush=True,
            )
    del model, wrapper
    if device.type == "mps":
        torch.mps.empty_cache()
    return output


def pair_metadata(left: SpectrumRow, right: SpectrumRow) -> np.ndarray:
    left_ce = left.collision_energy
    right_ce = right.collision_energy
    same_ce_unit = bool(
        left.collision_energy_unit
        and left.collision_energy_unit == right.collision_energy_unit
    )
    ce_comparable = left_ce is not None and right_ce is not None and same_ce_unit
    ce_delta = (
        min(abs(left_ce - right_ce) / 100.0, 2.0)
        if ce_comparable else 0.0
    )
    return np.asarray([
        float(ce_comparable),
        float(ce_delta),
        float(bool(left.instrument_type) and (
            left.instrument_type.strip().lower()
            == right.instrument_type.strip().lower()
        )),
        float(bool(left.dataset_source) and (
            left.dataset_source == right.dataset_source
        )),
        min(abs(left.precursor_mz - right.precursor_mz) / 0.05, 2.0),
    ], dtype=np.float32)


class IsoPairHead(nn.Module):
    """Symmetric pair classifier over frozen spectral representations."""

    def __init__(self, embedding_dim: int, metadata_dim: int = 5):
        super().__init__()
        self.embedding_dim = int(embedding_dim)
        self.metadata_dim = int(metadata_dim)
        self.network = nn.Sequential(
            nn.Linear(2 * embedding_dim + metadata_dim + 1, 384),
            nn.LayerNorm(384),
            nn.GELU(),
            nn.Dropout(0.12),
            nn.Linear(384, 96),
            nn.GELU(),
            nn.Dropout(0.08),
            nn.Linear(96, 1),
        )

    def forward(
        self,
        left: torch.Tensor,
        right: torch.Tensor,
        metadata: torch.Tensor,
    ) -> torch.Tensor:
        left = F.normalize(left, dim=-1)
        right = F.normalize(right, dim=-1)
        difference = torch.abs(left - right)
        product = left * right
        cosine = torch.sum(product, dim=-1, keepdim=True)
        features = torch.cat([difference, product, cosine, metadata], dim=-1)
        return self.network(features).squeeze(-1)


def _family_index(
    spectra: Sequence[SpectrumRow],
    splits: Sequence[str],
    split: str,
) -> dict[tuple[str, str], dict[str, list[int]]]:
    families: dict[tuple[str, str], dict[str, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for index, (row, row_split) in enumerate(zip(spectra, splits)):
        if row_split == split:
            families[(split_formula_key(row.formula), row.adduct)][
                row.connectivity
            ].append(index)
    return {
        family: dict(identities)
        for family, identities in families.items()
        if sum(len(indices) >= 2 for indices in identities.values()) >= 2
    }


def _spectrum_signature(row: SpectrumRow) -> tuple[tuple[float, float], ...]:
    maximum = max((intensity for _, intensity in row.peaks), default=1.0)
    return tuple(
        (round(mz, 3), round(intensity / maximum, 3))
        for mz, intensity in row.peaks
    )


def _prefer_cross_source_pair(
    indices: Sequence[int],
    spectra: Sequence[SpectrumRow],
    rng: random.Random,
) -> tuple[int, int]:
    pairs = [
        (left, right)
        for offset, left in enumerate(indices)
        for right in indices[offset + 1:]
    ]
    nonduplicate = [
        pair for pair in pairs
        if _spectrum_signature(spectra[pair[0]])
        != _spectrum_signature(spectra[pair[1]])
    ]
    cross_source = [
        pair for pair in pairs
        if spectra[pair[0]].dataset_source
        and spectra[pair[0]].dataset_source != spectra[pair[1]].dataset_source
        and _spectrum_signature(spectra[pair[0]])
        != _spectrum_signature(spectra[pair[1]])
    ]
    return rng.choice(cross_source or nonduplicate or pairs)


def make_training_pairs(
    spectra: Sequence[SpectrumRow],
    embeddings: np.ndarray,
    splits: Sequence[str],
    *,
    seed: int,
    positives_per_identity: int,
    hard_negatives_per_positive: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mine same-formula negatives that are hardest in the frozen embedding."""
    rng = random.Random(seed)
    families = _family_index(spectra, splits, "train")
    left_indices: list[int] = []
    right_indices: list[int] = []
    labels: list[float] = []

    for identities in families.values():
        identity_items = [
            (identity, indices)
            for identity, indices in sorted(identities.items())
            if len(indices) >= 2
        ]
        if len(identity_items) < 2:
            continue
        all_other = {
            identity: [
                index
                for other_identity, other_indices in identity_items
                if other_identity != identity
                for index in other_indices
            ]
            for identity, _ in identity_items
        }
        for identity, indices in identity_items:
            for _ in range(positives_per_identity):
                anchor, positive = _prefer_cross_source_pair(indices, spectra, rng)
                anchor_signature = _spectrum_signature(spectra[anchor])
                candidates = [
                    index for index in all_other[identity]
                    if _spectrum_signature(spectra[index]) != anchor_signature
                ]
                if not candidates:
                    continue
                left_indices.append(anchor)
                right_indices.append(positive)
                labels.append(1.0)
                anchor_embedding = embeddings[anchor]
                similarities = embeddings[candidates] @ anchor_embedding
                hardest_order = np.argsort(-similarities)
                n_hard = min(hard_negatives_per_positive, len(candidates))
                for position in hardest_order[:n_hard]:
                    left_indices.append(anchor)
                    right_indices.append(candidates[int(position)])
                    labels.append(0.0)

    order = list(range(len(labels)))
    rng.shuffle(order)
    return (
        np.asarray([left_indices[index] for index in order], dtype=np.int64),
        np.asarray([right_indices[index] for index in order], dtype=np.int64),
        np.asarray([labels[index] for index in order], dtype=np.float32),
    )


def train_head(
    spectra: Sequence[SpectrumRow],
    embeddings: np.ndarray,
    pairs: tuple[np.ndarray, np.ndarray, np.ndarray],
    *,
    device: torch.device,
    seed: int,
    epochs: int,
    batch_size: int,
) -> IsoPairHead:
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = IsoPairHead(embeddings.shape[1]).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=2e-4, weight_decay=1e-4)
    left_indices, right_indices, labels = pairs
    rng = np.random.default_rng(seed)
    embeddings_tensor = torch.from_numpy(embeddings)
    positive_count = max(int((labels >= 0.5).sum()), 1)
    negative_count = max(int((labels < 0.5).sum()), 1)
    positive_weight = torch.tensor(
        negative_count / positive_count,
        dtype=torch.float32,
        device=device,
    )

    for epoch in range(epochs):
        model.train()
        order = rng.permutation(len(labels))
        running_loss = 0.0
        running_correct = 0
        seen = 0
        for start in range(0, len(order), batch_size):
            batch_order = order[start:start + batch_size]
            left_np = left_indices[batch_order]
            right_np = right_indices[batch_order]
            left = embeddings_tensor[left_np].to(device)
            right = embeddings_tensor[right_np].to(device)
            metadata = torch.from_numpy(np.stack([
                pair_metadata(spectra[int(a)], spectra[int(b)])
                for a, b in zip(left_np, right_np)
            ])).to(device)
            target = torch.from_numpy(labels[batch_order]).to(device)
            optimizer.zero_grad(set_to_none=True)
            logits = model(left, right, metadata)
            loss = F.binary_cross_entropy_with_logits(
                logits, target, pos_weight=positive_weight
            )
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            running_loss += float(loss.detach().cpu()) * len(batch_order)
            running_correct += int(
                ((logits.detach() >= 0) == (target >= 0.5)).sum().cpu()
            )
            seen += len(batch_order)
        print(
            f"epoch {epoch + 1:02d}/{epochs}: "
            f"loss={running_loss / max(seen, 1):.5f} "
            f"pair_acc={running_correct / max(seen, 1):.4f}",
            flush=True,
        )
    return model.eval()


def build_query_cases(
    spectra: Sequence[SpectrumRow],
    splits: Sequence[str],
    split: str,
    *,
    max_queries_per_identity: int = 1,
) -> list[QueryCase]:
    families = _family_index(spectra, splits, split)
    cases: list[QueryCase] = []
    for (formula, adduct), identities in sorted(families.items()):
        eligible = {
            identity: tuple(indices)
            for identity, indices in identities.items()
            if len(indices) >= 2
        }
        if len(eligible) < 2:
            continue
        for correct_identity, indices in sorted(eligible.items()):
            for query_index in indices[:max_queries_per_identity]:
                candidates: list[tuple[str, tuple[int, ...]]] = []
                for identity, references in sorted(eligible.items()):
                    available = tuple(
                        index for index in references if index != query_index
                    )
                    if available:
                        candidates.append((identity, available))
                if len(candidates) >= 2:
                    cases.append(QueryCase(
                        query_index=query_index,
                        correct_connectivity=correct_identity,
                        candidate_indices=tuple(candidates),
                        formula=formula,
                        adduct=adduct,
                    ))
    return cases


def _cosine_score(
    left: SpectrumRow,
    right: SpectrumRow,
) -> float:
    """Standalone equivalent of V7's production Modified Cosine score."""
    from matchms.similarity.spectrum_similarity_functions import (
        collect_peak_pairs,
        score_best_matches,
    )
    def production_view(row: SpectrumRow) -> np.ndarray:
        peaks = np.asarray(row.peaks, dtype=float)
        if peaks.size == 0:
            return peaks.reshape(0, 2)
        relative = peaks[:, 1] >= float(peaks[:, 1].max()) * 0.05
        peaks = peaks[relative]
        peaks = peaks[np.abs(peaks[:, 0] - row.precursor_mz) > 0.5]
        peaks = peaks[peaks[:, 0] <= row.precursor_mz + 1.5]
        if len(peaks) > 15:
            strongest = np.argsort(peaks[:, 1], kind="stable")[::-1][:15]
            peaks = peaks[strongest]
        return peaks[np.argsort(peaks[:, 0], kind="stable")]

    query = production_view(left).copy()
    library = production_view(right).copy()
    if query.size == 0 or library.size == 0:
        return 0.0
    query = query[np.argsort(query[:, 0], kind="mergesort")]
    library = library[np.argsort(library[:, 0], kind="mergesort")]
    query[:, 1] /= max(float(query[:, 1].max()), 1e-12)
    library[:, 1] /= max(float(library[:, 1].max()), 1e-12)
    tolerance = 0.10
    shift = left.precursor_mz - right.precursor_mz
    zero_pairs = collect_peak_pairs(query, library, tolerance, shift=0.0)
    pair_sets = [] if zero_pairs is None else [zero_pairs]
    if abs(shift) > tolerance:
        shifted_pairs = collect_peak_pairs(query, library, tolerance, shift=shift)
        if shifted_pairs is not None:
            pair_sets.append(shifted_pairs)
    if not pair_sets:
        return 0.0
    pairs = pair_sets[0] if len(pair_sets) == 1 else np.concatenate(pair_sets)
    pairs = pairs[np.argsort(pairs[:, 2], kind="mergesort")[::-1], :]
    value = float(score_best_matches(pairs, query, library, 0.0, 1.0)[0])
    return value if math.isfinite(value) else 0.0


def score_cases(
    model: IsoPairHead,
    spectra: Sequence[SpectrumRow],
    embeddings: np.ndarray,
    cases: Sequence[QueryCase],
    *,
    device: torch.device,
) -> list[dict[str, Any]]:
    model.eval()
    output: list[dict[str, Any]] = []
    with torch.inference_mode():
        for case in cases:
            query = spectra[case.query_index]
            candidate_rows: list[dict[str, Any]] = []
            for identity, reference_indices in case.candidate_indices:
                baseline_values = [
                    _cosine_score(query, spectra[index])
                    for index in reference_indices
                ]
                left = torch.from_numpy(np.repeat(
                    embeddings[case.query_index][None, :],
                    len(reference_indices), axis=0,
                )).to(device)
                right = torch.from_numpy(embeddings[list(reference_indices)]).to(device)
                metadata = torch.from_numpy(np.stack([
                    pair_metadata(query, spectra[index])
                    for index in reference_indices
                ])).to(device)
                probability = torch.sigmoid(model(left, right, metadata))
                pair_probabilities = probability.float().cpu().numpy()
                candidate_rows.append({
                    "connectivity": identity,
                    "baseline": float(max(baseline_values)),
                    "isopair": float(max(pair_probabilities)),
                    "frozen_ssl": float(max(
                        embeddings[list(reference_indices)]
                        @ embeddings[case.query_index]
                    )),
                })
            output.append({
                "formula": case.formula,
                "adduct": case.adduct,
                "query_index": case.query_index,
                "correct_connectivity": case.correct_connectivity,
                "candidates": candidate_rows,
            })
    return output


def _rank_metrics(rows: Sequence[dict[str, Any]], score_key: str) -> dict[str, float]:
    correct = 0
    reciprocal_ranks: list[float] = []
    margins: list[float] = []
    winner_scores: list[float] = []
    winner_correct: list[bool] = []
    for row in rows:
        ranked = sorted(
            row["candidates"],
            key=lambda candidate: (
                -float(candidate[score_key]),
                str(candidate["connectivity"]),
            ),
        )
        correct_identity = row["correct_connectivity"]
        rank = next(
            index + 1 for index, candidate in enumerate(ranked)
            if candidate["connectivity"] == correct_identity
        )
        is_correct = rank == 1
        correct += int(is_correct)
        reciprocal_ranks.append(1.0 / rank)
        winner_scores.append(float(ranked[0][score_key]))
        margins.append(
            float(ranked[0][score_key]) - float(ranked[1][score_key])
            if len(ranked) > 1 else float(ranked[0][score_key])
        )
        winner_correct.append(is_correct)
    return {
        "n_queries": len(rows),
        "top1_accuracy": correct / max(len(rows), 1),
        "mean_reciprocal_rank": float(np.mean(reciprocal_ranks)) if rows else 0.0,
        "winner_scores": winner_scores,
        "margins": margins,
        "winner_correct": winner_correct,
    }


def apply_fusion(rows: Sequence[dict[str, Any]], alpha: float) -> None:
    for row in rows:
        for candidate in row["candidates"]:
            candidate["fused"] = (
                alpha * float(candidate["baseline"])
                + (1.0 - alpha) * float(candidate["isopair"])
            )


def tune_fusion(validation_rows: Sequence[dict[str, Any]]) -> float:
    best = (-1.0, -1.0, 1.0)
    best_alpha = 1.0
    for alpha in np.linspace(0.0, 1.0, 21):
        apply_fusion(validation_rows, float(alpha))
        metrics = _rank_metrics(validation_rows, "fused")
        key = (
            float(metrics["top1_accuracy"]),
            float(metrics["mean_reciprocal_rank"]),
            -abs(float(alpha) - 0.5),
        )
        if key > best:
            best = key
            best_alpha = float(alpha)
    return best_alpha


def _wilson_lower(correct: int, total: int, z: float = 1.96) -> float:
    if total <= 0:
        return 0.0
    proportion = correct / total
    denominator = 1.0 + z * z / total
    centre = proportion + z * z / (2.0 * total)
    spread = z * math.sqrt(
        proportion * (1.0 - proportion) / total
        + z * z / (4.0 * total * total)
    )
    return max(0.0, (centre - spread) / denominator)


def tune_selective_high(
    validation_rows: Sequence[dict[str, Any]],
    *,
    target_precision: float,
    min_calls: int,
) -> dict[str, float]:
    metrics = _rank_metrics(validation_rows, "fused")
    scores = np.asarray(metrics["winner_scores"], dtype=float)
    margins = np.asarray(metrics["margins"], dtype=float)
    correct = np.asarray(metrics["winner_correct"], dtype=bool)
    score_grid = sorted(set(np.quantile(scores, np.linspace(0.25, 0.95, 15))))
    margin_grid = sorted(set(np.quantile(margins, np.linspace(0.10, 0.95, 18))))
    best: tuple[int, float, float, float, float] | None = None
    for score_threshold in score_grid:
        for margin_threshold in margin_grid:
            selected = (scores >= score_threshold) & (margins >= margin_threshold)
            total = int(selected.sum())
            if total < min_calls:
                continue
            n_correct = int(correct[selected].sum())
            precision = n_correct / total
            lower = _wilson_lower(n_correct, total)
            # Qualification is based on the uncertainty-aware lower bound,
            # not the optimistic point estimate selected on this same cohort.
            if lower < target_precision:
                continue
            key = (total, lower, precision, -score_threshold, -margin_threshold)
            if best is None or key > best:
                best = key
                chosen = {
                    "score_threshold": float(score_threshold),
                    "margin_threshold": float(margin_threshold),
                    "validation_calls": total,
                    "validation_precision": precision,
                    "validation_wilson_lower_95": lower,
                }
    if best is None:
        return {
            "score_threshold": 1.0,
            "margin_threshold": 1.0,
            "validation_calls": 0,
            "validation_precision": 0.0,
            "validation_wilson_lower_95": 0.0,
        }
    return chosen


def selective_metrics(
    rows: Sequence[dict[str, Any]],
    thresholds: dict[str, float],
) -> dict[str, float]:
    metrics = _rank_metrics(rows, "fused")
    scores = np.asarray(metrics["winner_scores"], dtype=float)
    margins = np.asarray(metrics["margins"], dtype=float)
    correct = np.asarray(metrics["winner_correct"], dtype=bool)
    selected = (
        (scores >= thresholds["score_threshold"])
        & (margins >= thresholds["margin_threshold"])
    )
    total = int(selected.sum())
    n_correct = int(correct[selected].sum())
    return {
        "calls": total,
        "coverage": total / max(len(rows), 1),
        "precision": n_correct / max(total, 1),
        "wilson_lower_95": _wilson_lower(n_correct, total),
    }


def baseline_high_metrics(rows: Sequence[dict[str, Any]]) -> dict[str, float]:
    """Approximate V7's exact-identity High guard inside isomer families."""
    metrics = _rank_metrics(rows, "baseline")
    scores = np.asarray(metrics["winner_scores"], dtype=float)
    margins = np.asarray(metrics["margins"], dtype=float)
    correct = np.asarray(metrics["winner_correct"], dtype=bool)
    selected = (scores >= 0.90) & (margins > 0.05)
    total = int(selected.sum())
    n_correct = int(correct[selected].sum())
    return {
        "calls": total,
        "coverage": total / max(len(rows), 1),
        "precision": n_correct / max(total, 1),
        "wilson_lower_95": _wilson_lower(n_correct, total),
        "definition": "winner modified cosine >=0.90 and same-formula margin >0.05",
    }


def baseline_ambiguous_rows(
    rows: Sequence[dict[str, Any]],
    *,
    maximum_margin: float = 0.05,
) -> list[dict[str, Any]]:
    """Cases for which production cosine leaves two identities competitive."""
    selected = []
    for row in rows:
        scores = sorted(
            (float(candidate["baseline"]) for candidate in row["candidates"]),
            reverse=True,
        )
        if len(scores) >= 2 and scores[0] - scores[1] <= maximum_margin:
            selected.append(row)
    return selected


def high_similarity_isomer_rows(
    rows: Sequence[dict[str, Any]],
    *,
    minimum_wrong_similarity: float = 0.70,
) -> list[dict[str, Any]]:
    """Cases with at least one strongly similar wrong same-formula identity."""
    return [
        row for row in rows
        if max(
            (
                float(candidate["baseline"])
                for candidate in row["candidates"]
                if candidate["connectivity"] != row["correct_connectivity"]
            ),
            default=0.0,
        ) >= minimum_wrong_similarity
    ]


def subset_report(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    return {
        score_key: _jsonable_metrics(_rank_metrics(rows, score_key))
        for score_key in ("baseline", "frozen_ssl", "isopair", "fused")
    }


def _jsonable_metrics(metrics: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value for key, value in metrics.items()
        if key not in {"winner_scores", "margins", "winner_correct"}
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--dreams-root", type=Path, default=DEFAULT_DREAMS_ROOT)
    parser.add_argument("--ssl-checkpoint", type=Path)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--seed", type=int, default=3407)
    parser.add_argument("--train-families", type=int, default=600)
    parser.add_argument("--validation-families", type=int, default=150)
    parser.add_argument("--test-families", type=int, default=150)
    parser.add_argument("--min-identities", type=int, default=2)
    parser.add_argument("--max-identities", type=int, default=12)
    parser.add_argument("--spectra-per-identity", type=int, default=4)
    parser.add_argument("--embedding-batch-size", type=int, default=64)
    parser.add_argument("--training-batch-size", type=int, default=256)
    parser.add_argument("--epochs", type=int, default=12)
    parser.add_argument("--positives-per-identity", type=int, default=2)
    parser.add_argument("--hard-negatives", type=int, default=2)
    parser.add_argument("--target-high-precision", type=float, default=0.95)
    parser.add_argument("--min-high-calls", type=int, default=20)
    parser.add_argument("--reuse-cache", action="store_true")
    args = parser.parse_args()

    checkpoint = args.ssl_checkpoint or (
        args.dreams_root / "dreams" / "models" / "pretrained" / "ssl_model.ckpt"
    )
    if not args.db.exists():
        raise SystemExit(f"Database does not exist: {args.db}")
    if not checkpoint.exists():
        raise SystemExit(f"DreaMS SSL checkpoint does not exist: {checkpoint}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    device = _device(args.device)
    print(f"device={device} db={args.db} checkpoint={checkpoint}")

    cache = args.out_dir / "pilot_embeddings.npz"
    metadata_cache = args.out_dir / "pilot_spectra.json"
    if args.reuse_cache and cache.exists() and metadata_cache.exists():
        payload = np.load(cache)
        embeddings = payload["embeddings"].astype(np.float32)
        split_values = payload["splits"].astype(str).tolist()
        raw_rows = json.loads(metadata_cache.read_text(encoding="utf-8"))
        spectra = [SpectrumRow(
            **{
                **row,
                "peaks": tuple(tuple(peak) for peak in row["peaks"]),
            }
        ) for row in raw_rows]
        split_values = enforce_formula_and_connectivity_disjoint_splits(
            spectra, split_values, seed=args.seed
        )
        print(f"reused cache: {len(spectra)} spectra")
    else:
        connection = _connect_readonly(args.db)
        selected = select_family_strata(
            connection,
            seed=args.seed,
            train_families=args.train_families,
            validation_families=args.validation_families,
            test_families=args.test_families,
            min_identities=args.min_identities,
            max_identities=args.max_identities,
        )
        print("selected families:", {key: len(value) for key, value in selected.items()})
        spectra, split_values = fetch_spectra(
            connection,
            selected,
            max_spectra_per_identity=args.spectra_per_identity,
        )
        connection.close()
        spectra, split_values = filter_usable_families(spectra, split_values)
        split_values = enforce_formula_and_connectivity_disjoint_splits(
            spectra, split_values, seed=args.seed
        )
        print(
            "usable spectra:", len(spectra),
            "split counts:", {
                split: split_values.count(split)
                for split in ("train", "validation", "test")
            },
        )
        embeddings = compute_ssl_embeddings(
            spectra,
            dreams_root=args.dreams_root,
            checkpoint=checkpoint,
            batch_size=args.embedding_batch_size,
            device=device,
        )
        np.savez_compressed(
            cache,
            embeddings=embeddings.astype(np.float16),
            splits=np.asarray(split_values),
        )
        metadata_cache.write_text(
            json.dumps([asdict(row) for row in spectra], ensure_ascii=False),
            encoding="utf-8",
        )

    split_audit = split_overlap_audit(spectra, split_values)
    if not split_audit["zero_overlap"]:
        raise RuntimeError(f"IsoPair split leakage remains: {split_audit}")

    pairs = make_training_pairs(
        spectra,
        embeddings,
        split_values,
        seed=args.seed,
        positives_per_identity=args.positives_per_identity,
        hard_negatives_per_positive=args.hard_negatives,
    )
    positives = int((pairs[2] == 1).sum())
    negatives = int((pairs[2] == 0).sum())
    print(f"training pairs: {len(pairs[2])} positive={positives} negative={negatives}")
    model = train_head(
        spectra,
        embeddings,
        pairs,
        device=device,
        seed=args.seed,
        epochs=args.epochs,
        batch_size=args.training_batch_size,
    )

    validation_cases = build_query_cases(spectra, split_values, "validation")
    test_cases = build_query_cases(spectra, split_values, "test")
    print(f"retrieval cases: validation={len(validation_cases)} test={len(test_cases)}")
    validation_rows = score_cases(
        model, spectra, embeddings, validation_cases, device=device
    )
    test_rows = score_cases(model, spectra, embeddings, test_cases, device=device)

    fusion_alpha = tune_fusion(validation_rows)
    apply_fusion(validation_rows, fusion_alpha)
    apply_fusion(test_rows, fusion_alpha)
    high_thresholds = tune_selective_high(
        validation_rows,
        target_precision=args.target_high_precision,
        min_calls=args.min_high_calls,
    )

    report = {
        "schema_version": 1,
        "method": "IsoPair-MS formula-conditional hard-negative pair ranker",
        "created_at_unix": time.time(),
        "seed": args.seed,
        "split_policy": (
            "SHA-256 normalized-formula/connectivity-union-component "
            "disjoint 70/15/15"
        ),
        "split_audit": split_audit,
        "identity_definition": "first 14 InChIKey characters (2-D connectivity)",
        "stereochemistry_claimed": False,
        "source_tiers": ["experimental", "curated", "massbank_supplement"],
        "n_spectra": len(spectra),
        "n_training_pairs": len(pairs[2]),
        "n_validation_queries": len(validation_rows),
        "n_test_queries": len(test_rows),
        "fusion_alpha_baseline": fusion_alpha,
        "fusion_alpha_isopair": 1.0 - fusion_alpha,
        "validation": {
            "baseline": _jsonable_metrics(_rank_metrics(validation_rows, "baseline")),
            "frozen_ssl": _jsonable_metrics(_rank_metrics(validation_rows, "frozen_ssl")),
            "isopair": _jsonable_metrics(_rank_metrics(validation_rows, "isopair")),
            "fused": _jsonable_metrics(_rank_metrics(validation_rows, "fused")),
            "selective_high": selective_metrics(validation_rows, high_thresholds),
            "baseline_ambiguous_subset": subset_report(
                baseline_ambiguous_rows(validation_rows)
            ),
            "high_similarity_isomer_subset": subset_report(
                high_similarity_isomer_rows(validation_rows)
            ),
        },
        "test": {
            "baseline": _jsonable_metrics(_rank_metrics(test_rows, "baseline")),
            "frozen_ssl": _jsonable_metrics(_rank_metrics(test_rows, "frozen_ssl")),
            "isopair": _jsonable_metrics(_rank_metrics(test_rows, "isopair")),
            "fused": _jsonable_metrics(_rank_metrics(test_rows, "fused")),
            "baseline_high": baseline_high_metrics(test_rows),
            "selective_high": selective_metrics(test_rows, high_thresholds),
            "baseline_ambiguous_subset": subset_report(
                baseline_ambiguous_rows(test_rows)
            ),
            "high_similarity_isomer_subset": subset_report(
                high_similarity_isomer_rows(test_rows)
            ),
        },
        "selective_high_thresholds": high_thresholds,
        "limitations": [
            "Pilot validation uses local public-library records and is not a substitute for an external acquisition-standard isomer panel.",
            "The self-supervised DreaMS backbone may have seen unlabeled spectra from overlapping public repositories; the supervised pair head and all evaluation families are formula-disjoint.",
            "Ordinary achiral CID/HCD cannot support general stereoisomer claims; labels and metrics target 2-D connectivity only.",
        ],
    }

    artifact = {
        "schema_version": 1,
        "model_name": "IsoPair-MS",
        "embedding_dim": embeddings.shape[1],
        "metadata_dim": 5,
        "state_dict": {key: value.detach().cpu() for key, value in model.state_dict().items()},
        "fusion_alpha_baseline": fusion_alpha,
        "fusion_alpha_isopair": 1.0 - fusion_alpha,
        "selective_high_thresholds": high_thresholds,
        "training_report": report,
        "ssl_checkpoint_name": checkpoint.name,
    }
    torch.save(artifact, args.out_dir / "isopair_ms_v1.pt")
    (args.out_dir / "locked_evaluation.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (args.out_dir / "test_predictions.json").write_text(
        json.dumps(test_rows, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
