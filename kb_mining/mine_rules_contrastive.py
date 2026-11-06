#!/usr/bin/env python3
"""Mine conservative class-specific MS/MS rule candidates.

Unlike the legacy miner, this implementation:

* streams rows instead of ``fetchall``;
* counts an InChIKey/structure at most once per polarity and stratum;
* excludes predicted/quarantined spectra;
* uses relative-intensity and precursor/charge guards;
* requires support from multiple independent dataset sources; and
* compares class prevalence with the out-of-class background.

The output is a candidate catalogue and is never consumed automatically by the
runtime scorer.  Candidates must still receive formula/structure/literature
curation and pass a benchmark before promotion.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "local_library.db"
DEFAULT_OUT = ROOT / "kb_v3" / "empirical_rules_candidates_v2.json"

VALID_SOURCE_TIERS = {"experimental", "curated", "merged"}
EXCLUDED_CLASSES = {"", "other", "unknown", "n/a", "na", "none", "blank", "generic"}
KNOWN_ADDUCTS = {
    "[M]+", "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+",
    "[M]-", "[M-H]-", "[M+Cl]-", "[M+FA-H]-", "[M+Hac-H]-", "[M-H-H2O]-",
}


@dataclass(frozen=True)
class MiningConfig:
    bin_width: float = 0.02
    min_relative_intensity: float = 0.01
    min_independent_compounds: int = 20
    min_source_families: int = 2
    min_compounds_per_source: int = 3
    max_source_dominance: float = 0.80
    min_prevalence: float = 0.08
    min_likelihood_ratio: float = 2.0
    max_background_prevalence: float = 0.35
    top_k: int = 30
    min_fragment_mz: float = 30.0
    max_fragment_mz: float = 1000.0
    min_neutral_loss: float = 1.0
    max_neutral_loss: float = 400.0
    precursor_exclusion_da: float = 1.0


def _polarity(adduct: str, charge: Any) -> str:
    text = str(adduct or "").strip().upper().replace(" ", "")
    if text.endswith("+"):
        return "pos"
    if text.endswith("-"):
        return "neg"
    if text in {"POS", "POSITIVE", "+"}:
        return "pos"
    if text in {"NEG", "NEGATIVE", "-"}:
        return "neg"
    try:
        value = int(charge)
    except (TypeError, ValueError):
        return "unknown"
    return "pos" if value > 0 else "neg" if value < 0 else "unknown"


def _normalize_adduct(adduct: str) -> str:
    try:
        from matching.rule_quality import normalize_adduct
        value = normalize_adduct(adduct)
    except Exception:
        value = str(adduct or "").strip()
    return value if value in KNOWN_ADDUCTS else ""


def _identity(row: sqlite3.Row) -> str:
    inchikey = str(row["inchikey"] or "").strip().upper()
    if inchikey:
        return f"inchikey:{inchikey}"
    smiles = str(row["smiles"] or "").strip()
    if smiles:
        return f"smiles:{smiles}"
    return f"compound_id:{row['compound_id']}"


def _source_family(row: sqlite3.Row) -> str:
    source = str(row["dataset_source"] or "").strip()
    if source:
        return source
    source_file = Path(str(row["source_file"] or "unknown")).name
    return source_file or "unknown"


def _bin(value: float, width: float) -> int:
    return int(round(float(value) / width))


def _new_group(row: sqlite3.Row, already_seen: bool) -> dict[str, Any]:
    return {
        "compound_id": row["compound_id"],
        "identity": _identity(row),
        "class": str(row["cls"] or "").strip().lower(),
        "skip": already_seen,
        "polarity": defaultdict(lambda: {
            "diag": defaultdict(lambda: [0.0, 0]),
            "loss": defaultdict(lambda: [0.0, 0]),
            "diag_sources": defaultdict(set),
            "loss_sources": defaultdict(set),
            "sources": set(),
            "adducts": defaultdict(lambda: {
                "diag": defaultdict(lambda: [0.0, 0]),
                "loss": defaultdict(lambda: [0.0, 0]),
                "diag_sources": defaultdict(set),
                "loss_sources": defaultdict(set),
            }),
        }),
    }


def _accumulate_spectrum(group: dict[str, Any], row: sqlite3.Row, cfg: MiningConfig) -> bool:
    if group["skip"]:
        return False
    pol = _polarity(row["adduct"], row["charge"])
    if pol == "unknown":
        return False
    try:
        charge_magnitude = abs(int(row["charge"] or 1))
    except (TypeError, ValueError):
        match = re.search(r"\](\d+)[+-]$", str(row["adduct"] or ""))
        charge_magnitude = int(match.group(1)) if match else 1
    # Broad small-molecule rules must not be contaminated by doubly/triply
    # charged peptide fragments.  Multi-charge mining requires a separate,
    # charge-normalized workflow.
    if charge_magnitude != 1:
        return False
    try:
        precursor = float(row["precursor_mz"])
        peaks = json.loads(row["peaks_json"] or "[]")
    except (TypeError, ValueError, json.JSONDecodeError):
        return False
    parsed = []
    for peak in peaks:
        try:
            mz, intensity = float(peak[0]), float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if mz > 0 and intensity > 0 and math.isfinite(mz) and math.isfinite(intensity):
            parsed.append((mz, intensity))
    if precursor <= 0 or not parsed:
        return False

    maximum = max(intensity for _, intensity in parsed)
    threshold = maximum * cfg.min_relative_intensity
    source = _source_family(row)
    adduct = _normalize_adduct(row["adduct"])
    state = group["polarity"][pol]
    state["sources"].add(source)
    seen_diag: dict[int, float] = {}
    seen_loss: dict[int, float] = {}
    for mz, intensity in parsed:
        if intensity < threshold or mz >= precursor - cfg.precursor_exclusion_da:
            continue
        if cfg.min_fragment_mz <= mz <= cfg.max_fragment_mz:
            seen_diag[_bin(mz, cfg.bin_width)] = mz
        loss = precursor - mz
        if cfg.min_neutral_loss <= loss <= cfg.max_neutral_loss:
            seen_loss[_bin(loss, cfg.bin_width)] = loss

    for bin_id, mz in seen_diag.items():
        state["diag"][bin_id][0] += mz
        state["diag"][bin_id][1] += 1
        state["diag_sources"][bin_id].add(source)
        if adduct:
            state["adducts"][adduct]["diag"][bin_id][0] += mz
            state["adducts"][adduct]["diag"][bin_id][1] += 1
            state["adducts"][adduct]["diag_sources"][bin_id].add(source)
    for bin_id, loss in seen_loss.items():
        state["loss"][bin_id][0] += loss
        state["loss"][bin_id][1] += 1
        state["loss_sources"][bin_id].add(source)
        if adduct:
            state["adducts"][adduct]["loss"][bin_id][0] += loss
            state["adducts"][adduct]["loss"][bin_id][1] += 1
            state["adducts"][adduct]["loss_sources"][bin_id].add(source)
    return True


def _finalize_group(group: dict[str, Any], aggregate: dict[str, Any]) -> None:
    if group["skip"] or group["class"] in EXCLUDED_CLASSES:
        return
    cls = group["class"]
    for pol, state in group["polarity"].items():
        if not state["sources"]:
            continue
        key = (cls, pol)
        aggregate["class_totals"][key] += 1
        aggregate["global_totals"][pol] += 1
        for kind in ("diag", "loss"):
            for bin_id, (value_sum, value_n) in state[kind].items():
                aggregate[f"{kind}_counts"][key][bin_id] += 1
                aggregate[f"global_{kind}_counts"][pol][bin_id] += 1
                aggregate[f"{kind}_centroids"][key][bin_id][0] += value_sum / value_n
                aggregate[f"{kind}_centroids"][key][bin_id][1] += 1
                aggregate[f"{kind}_sources"][key][bin_id].update(
                    state[f"{kind}_sources"][bin_id]
                )
                for source in state[f"{kind}_sources"][bin_id]:
                    aggregate[f"{kind}_source_support"][key][bin_id][source] += 1
        for adduct, adduct_state in state["adducts"].items():
            adduct_key = (cls, pol, adduct)
            aggregate["adduct_totals"][adduct_key] += 1
            for kind in ("diag", "loss"):
                for bin_id, (value_sum, value_n) in adduct_state[kind].items():
                    aggregate[f"adduct_{kind}_counts"][adduct_key][bin_id] += 1
                    aggregate[f"adduct_{kind}_centroids"][adduct_key][bin_id][0] += value_sum / value_n
                    aggregate[f"adduct_{kind}_centroids"][adduct_key][bin_id][1] += 1
                    aggregate[f"adduct_{kind}_sources"][adduct_key][bin_id].update(
                        adduct_state[f"{kind}_sources"][bin_id]
                    )
                    for source in adduct_state[f"{kind}_sources"][bin_id]:
                        aggregate[f"adduct_{kind}_source_support"][adduct_key][bin_id][source] += 1


def _aggregate_template() -> dict[str, Any]:
    nested_counter = lambda: defaultdict(Counter)
    nested_centroid = lambda: defaultdict(lambda: defaultdict(lambda: [0.0, 0]))
    nested_sources = lambda: defaultdict(lambda: defaultdict(set))
    nested_source_support = lambda: defaultdict(lambda: defaultdict(Counter))
    return {
        "class_totals": Counter(), "global_totals": Counter(),
        "diag_counts": nested_counter(), "loss_counts": nested_counter(),
        "global_diag_counts": nested_counter(), "global_loss_counts": nested_counter(),
        "diag_centroids": nested_centroid(), "loss_centroids": nested_centroid(),
        "diag_sources": nested_sources(), "loss_sources": nested_sources(),
        "diag_source_support": nested_source_support(),
        "loss_source_support": nested_source_support(),
        "adduct_totals": Counter(),
        "adduct_diag_counts": nested_counter(), "adduct_loss_counts": nested_counter(),
        "adduct_diag_centroids": nested_centroid(), "adduct_loss_centroids": nested_centroid(),
        "adduct_diag_sources": nested_sources(), "adduct_loss_sources": nested_sources(),
        "adduct_diag_source_support": nested_source_support(),
        "adduct_loss_source_support": nested_source_support(),
    }


def _rule_candidates(
    key: tuple[str, str],
    kind: str,
    aggregate: dict[str, Any],
    cfg: MiningConfig,
) -> list[dict[str, Any]]:
    class_total = aggregate["class_totals"][key]
    pol = key[1]
    other_total = max(aggregate["global_totals"][pol] - class_total, 0)
    rows = []
    for bin_id, support in aggregate[f"{kind}_counts"][key].items():
        sources = aggregate[f"{kind}_sources"][key][bin_id]
        source_support = aggregate[f"{kind}_source_support"][key][bin_id]
        qualified_sources = {
            source: count for source, count in source_support.items()
            if count >= cfg.min_compounds_per_source
        }
        source_dominance = max(source_support.values(), default=0) / max(support, 1)
        if (
            support < cfg.min_independent_compounds
            or len(qualified_sources) < cfg.min_source_families
            or source_dominance > cfg.max_source_dominance
        ):
            continue
        prevalence = support / class_total if class_total else 0.0
        background_support = max(
            aggregate[f"global_{kind}_counts"][pol][bin_id] - support, 0
        )
        background_prevalence = background_support / other_total if other_total else 0.0
        smoothing = 0.5 / max(class_total, other_total, 1)
        likelihood_ratio = (prevalence + smoothing) / (background_prevalence + smoothing)
        if (
            prevalence < cfg.min_prevalence
            or likelihood_ratio < cfg.min_likelihood_ratio
            or background_prevalence > cfg.max_background_prevalence
        ):
            continue
        value_sum, value_n = aggregate[f"{kind}_centroids"][key][bin_id]
        centroid = value_sum / value_n
        score = prevalence * math.log2(max(likelihood_ratio, 1.0)) * math.log1p(support)
        rows.append({
            "mz" if kind == "diag" else "loss_da": round(centroid, 5),
            "prevalence": round(prevalence, 4),
            "background_prevalence": round(background_prevalence, 4),
            "likelihood_ratio": round(likelihood_ratio, 3),
            "n_independent_compounds": int(support),
            "n_background_compounds": int(background_support),
            "n_source_families": len(qualified_sources),
            "source_support": dict(sorted(qualified_sources.items())),
            "max_source_dominance": round(source_dominance, 4),
            "candidate_score": round(score, 4),
            "evidence_level": "empirical_candidate_not_enabled",
        })
    return sorted(rows, key=lambda row: (-row["candidate_score"], row.get("mz", row.get("loss_da"))))[:cfg.top_k]


def _adduct_candidates(
    key: tuple[str, str, str],
    kind: str,
    aggregate: dict[str, Any],
    cfg: MiningConfig,
) -> list[dict[str, Any]]:
    total = aggregate["adduct_totals"][key]
    rows = []
    for bin_id, support in aggregate[f"adduct_{kind}_counts"][key].items():
        sources = aggregate[f"adduct_{kind}_sources"][key][bin_id]
        source_support = aggregate[f"adduct_{kind}_source_support"][key][bin_id]
        qualified_sources = {
            source: count for source, count in source_support.items()
            if count >= cfg.min_compounds_per_source
        }
        source_dominance = max(source_support.values(), default=0) / max(support, 1)
        prevalence = support / total if total else 0.0
        if (
            support < cfg.min_independent_compounds
            or len(qualified_sources) < cfg.min_source_families
            or source_dominance > cfg.max_source_dominance
            or prevalence < cfg.min_prevalence
        ):
            continue
        value_sum, value_n = aggregate[f"adduct_{kind}_centroids"][key][bin_id]
        rows.append({
            "mz" if kind == "diag" else "loss_da": round(value_sum / value_n, 5),
            "prevalence_within_adduct": round(prevalence, 4),
            "n_independent_compounds": int(support),
            "n_source_families": len(qualified_sources),
            "source_support": dict(sorted(qualified_sources.items())),
            "max_source_dominance": round(source_dominance, 4),
            "evidence_level": "empirical_adduct_candidate_not_enabled",
        })
    return sorted(rows, key=lambda row: (-row["prevalence_within_adduct"], row.get("mz", row.get("loss_da"))))[:cfg.top_k]


def mine_database(db_path: Path, cfg: MiningConfig = MiningConfig()) -> dict[str, Any]:
    connection = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    query = """
        SELECT c.id AS compound_id, c.inchikey, c.smiles,
               COALESCE(NULLIF(c.classification_v3, ''), c.classification) AS cls,
               s.precursor_mz, s.peaks_json, s.adduct, s.charge,
               s.dataset_source, s.source_file, s.source_tier, s.collision_energy
        FROM Compounds c
        JOIN Mapping m ON m.compound_id = c.id
        JOIN Spectra s ON s.id = m.spectrum_id
        WHERE s.source_tier IN ('experimental', 'curated', 'merged')
          AND COALESCE(NULLIF(c.classification_v3, ''), c.classification) IS NOT NULL
        ORDER BY c.id
    """
    aggregate = _aggregate_template()
    seen_identities: set[str] = set()
    current_group = None
    rows_seen = rows_used = parse_failures = 0
    ce_audit: dict[str, list[float]] = defaultdict(lambda: [0.0, float("inf"), float("-inf")])

    for row in connection.execute(query):
        rows_seen += 1
        if current_group is None or row["compound_id"] != current_group["compound_id"]:
            if current_group is not None:
                _finalize_group(current_group, aggregate)
            identity = _identity(row)
            current_group = _new_group(row, identity in seen_identities)
            seen_identities.add(identity)
        if _accumulate_spectrum(current_group, row, cfg):
            rows_used += 1
        else:
            parse_failures += 1
        if row["collision_energy"] is not None:
            try:
                ce = float(row["collision_energy"])
                source = _source_family(row)
                ce_audit[source][0] += 1
                ce_audit[source][1] = min(ce_audit[source][1], ce)
                ce_audit[source][2] = max(ce_audit[source][2], ce)
            except (TypeError, ValueError):
                pass
    if current_group is not None:
        _finalize_group(current_group, aggregate)
    connection.close()

    classes: dict[str, Any] = {}
    for key, total in sorted(aggregate["class_totals"].items()):
        if total < cfg.min_independent_compounds:
            continue
        cls, pol = key
        entry = {
            "n_independent_compounds": int(total),
            "diagnostic_ions": _rule_candidates(key, "diag", aggregate, cfg),
            "neutral_losses": _rule_candidates(key, "loss", aggregate, cfg),
            "adduct_strata": {},
        }
        for adduct_key, adduct_total in sorted(aggregate["adduct_totals"].items()):
            if adduct_key[:2] != key or adduct_total < cfg.min_independent_compounds:
                continue
            adduct = adduct_key[2]
            entry["adduct_strata"][adduct] = {
                "n_independent_compounds": int(adduct_total),
                "diagnostic_ions": _adduct_candidates(adduct_key, "diag", aggregate, cfg),
                "neutral_losses": _adduct_candidates(adduct_key, "loss", aggregate, cfg),
            }
        classes.setdefault(cls, {})[pol] = entry

    ce_sources = {
        source: {"n_spectra": int(values[0]), "min": values[1], "max": values[2]}
        for source, values in sorted(ce_audit.items())
    }
    return {
        "meta": {
            "schema_version": "2.0-candidate",
            "database": str(db_path),
            "rows_seen": rows_seen,
            "rows_used": rows_used,
            "rows_skipped_or_unusable": parse_failures,
            "unique_identities_seen": len(seen_identities),
            "counting_unit": "one InChIKey/SMILES per polarity",
            "included_source_tiers": sorted(VALID_SOURCE_TIERS),
            "runtime_enabled": False,
            "collision_energy_policy": (
                "Database CE is unitless; native values are audited but never treated as eV/NCE "
                "or used for portable production rules."
            ),
            "config": cfg.__dict__,
        },
        "collision_energy_source_audit": ce_sources,
        "classes": classes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--min-compounds", type=int, default=20)
    parser.add_argument("--min-sources", type=int, default=2)
    parser.add_argument("--min-compounds-per-source", type=int, default=3)
    parser.add_argument("--min-prevalence", type=float, default=0.08)
    parser.add_argument("--min-lr", type=float, default=2.0)
    args = parser.parse_args()
    cfg = MiningConfig(
        min_independent_compounds=args.min_compounds,
        min_source_families=args.min_sources,
        min_compounds_per_source=args.min_compounds_per_source,
        min_prevalence=args.min_prevalence,
        min_likelihood_ratio=args.min_lr,
    )
    result = mine_database(args.db.resolve(), cfg)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result["meta"], ensure_ascii=False, indent=2))
    print(f"classes={len(result['classes'])}; wrote {args.output}")


if __name__ == "__main__":
    main()
