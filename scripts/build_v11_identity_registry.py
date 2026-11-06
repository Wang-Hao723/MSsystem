#!/usr/bin/env python3
"""Build a frozen PubChem name/formula identity-repair registry.

The input is candidate artifacts, never a ground-truth list.  Only names that
are bound to multiple 2-D connectivities inside one entity, or whose reported
InChIKey disagrees with the ranker's standardized structure, are queried.  A record is
accepted when PubChem resolves the exact submitted name, the returned formula
matches, and RDKit reproduces the returned full InChIKey.

The second trigger is important for real spectral libraries: a candidate can
carry a correct deposited InChIKey while an unrelated SMILES was attached
during aggregation.  Silently preferring either field changes the evaluated
connectivity.  The registry resolves that conflict from an external authority
without using sample truth or expected-compound labels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from matching.library_quality import normalize_formula_text  # noqa: E402
from matching.flagship_ranker import connectivity_identity  # noqa: E402

try:  # noqa: E402
    from rdkit import Chem
except Exception as exc:  # pragma: no cover - deployment dependency check
    raise SystemExit(f"RDKit is required: {exc}")


DEFAULT_OUTPUT = ROOT / "config" / "canonical_name_identity_registry.json"
PUBCHEM_PROPERTIES = (
    "Title,MolecularFormula,IsomericSMILES,ConnectivitySMILES,InChI,InChIKey"
)
FULL_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
CONNECTIVITY_BLOCK = re.compile(r"^[A-Z]{14}$")


def strict_name_key(value: Any) -> str:
    return " ".join(
        unicodedata.normalize("NFKC", str(value or "")).casefold().split()
    )


def _eligible_name(value: str) -> bool:
    lower = value.casefold()
    return bool(
        3 <= len(value) <= 100
        and "_" not in value
        and "known isomer" not in lower
        and "suspect related" not in lower
    )


def _connectivity_block(value: Any) -> str:
    token = str(value or "").strip().upper().replace("INCHIKEY=", "")
    if FULL_INCHIKEY.fullmatch(token):
        return token[:14]
    first = token.split("-", 1)[0]
    return first if CONNECTIVITY_BLOCK.fullmatch(first) else ""


def _candidate_display(candidate: Mapping[str, Any]) -> Mapping[str, Any]:
    nested = candidate.get("candidate")
    return nested if isinstance(nested, Mapping) else candidate


def _smiles_connectivity(value: Any) -> str:
    smiles = str(value or "").strip()
    if not smiles:
        return ""
    try:
        identity, _ = connectivity_identity(
            {"smiles": smiles}, standardize_tautomers=True
        )
        if identity.startswith("std2d:"):
            return _connectivity_block(identity.split(":", 1)[1])
        return ""
    except Exception:
        return ""


def _observed_identity_tokens(candidate: dict[str, Any]) -> set[str]:
    display = _candidate_display(candidate)
    values: set[str] = set()
    identity = str(candidate.get("identity_key") or "").strip()
    if identity:
        values.add(identity)
    standardized = (
        _connectivity_block(display.get("standardized_connectivity_inchikey"))
        or _smiles_connectivity(display.get("smiles"))
    )
    reported = _connectivity_block(
        display.get("reported_inchikey_before_identity_resolution")
        or display.get("reported_inchikey")
        or display.get("inchikey")
    )
    if standardized:
        values.add(f"std2d:{standardized}")
    if reported:
        values.add(f"ik14:{reported}")
    return values


def _candidate_has_unresolved_structure_conflict(candidate: dict[str, Any]) -> bool:
    display = _candidate_display(candidate)
    resolution = display.get("canonical_identity_resolution") or {}
    if resolution.get("status") == "resolved":
        return False
    if display.get("reported_connectivity_matches_standardized") is False:
        return True
    standardized = (
        _connectivity_block(display.get("standardized_connectivity_inchikey"))
        or _smiles_connectivity(display.get("smiles"))
    )
    reported = _connectivity_block(
        display.get("reported_inchikey") or display.get("inchikey")
    )
    return bool(standardized and reported and standardized != reported)


def _candidate_groups(payload: Any) -> Iterable[list[dict[str, Any]]]:
    """Yield candidate pools from flagship entities or immutable result rows."""

    if isinstance(payload, Mapping):
        for entity in payload.get("entities") or []:
            candidates = entity.get("candidates") if isinstance(entity, Mapping) else None
            if isinstance(candidates, list):
                yield [item for item in candidates if isinstance(item, dict)]
        return
    if not isinstance(payload, list):
        return
    for row in payload:
        if not isinstance(row, Mapping):
            continue
        match = row.get("match")
        if not isinstance(match, dict):
            continue
        candidates = [match]
        candidates.extend(
            item for item in (match.get("top_candidates") or [])
            if isinstance(item, dict)
        )
        yield candidates


def _merge_conflict(
    collected: dict[tuple[str, str], dict[str, Any]],
    *,
    key: tuple[str, str],
    query_name: str,
    identities: set[str],
    artifact_name: str,
    trigger: str,
) -> None:
    item = collected.setdefault(key, {
        "query_name": query_name or key[0],
        "name_key": key[0],
        "formula": key[1],
        "observed_connectivity_identities": set(),
        "source_artifacts": set(),
        "triggers": set(),
    })
    item["observed_connectivity_identities"].update(identities)
    item["source_artifacts"].add(artifact_name)
    item["triggers"].add(trigger)


def collect_conflicts(paths: list[Path]) -> dict[tuple[str, str], dict[str, Any]]:
    collected: dict[tuple[str, str], dict[str, Any]] = {}
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        for candidate_group in _candidate_groups(payload):
            groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
            for candidate in candidate_group:
                display = _candidate_display(candidate)
                name = str(display.get("compound_name") or "").strip()
                key = strict_name_key(name)
                formula = normalize_formula_text(display.get("formula"))
                if key and formula and _eligible_name(name):
                    group_key = (key, formula)
                    groups[group_key].append(candidate)
                    if _candidate_has_unresolved_structure_conflict(candidate):
                        _merge_conflict(
                            collected,
                            key=group_key,
                            query_name=name,
                            identities=_observed_identity_tokens(candidate),
                            artifact_name=path.name,
                            trigger="reported_inchikey_disagrees_with_stored_structure",
                        )
            for key, candidates in groups.items():
                identities = sorted({
                    str(candidate.get("identity_key") or "")
                    for candidate in candidates
                    if candidate.get("identity_key")
                })
                if len(identities) < 2:
                    continue
                _merge_conflict(
                    collected,
                    key=key,
                    query_name=str(
                        _candidate_display(candidates[0]).get("compound_name")
                        or key[0]
                    ),
                    identities=set(identities),
                    artifact_name=path.name,
                    trigger="same_name_formula_bound_to_multiple_connectivities",
                )
    return collected


def fetch_pubchem(name: str, timeout: float) -> dict[str, Any]:
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{urllib.parse.quote(name, safe='')}/property/{PUBCHEM_PROPERTIES}/JSON"
    )
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "MSsystem-V14-identity-registry/2.0"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = json.load(response)
    rows = list((payload.get("PropertyTable") or {}).get("Properties") or [])
    if len(rows) != 1:
        raise ValueError(f"expected_one_pubchem_resolution_got_{len(rows)}")
    return {**rows[0], "property_url": url}


def validate_resolution(conflict: dict[str, Any], row: dict[str, Any]) -> tuple[dict[str, Any] | None, str]:
    formula = normalize_formula_text(row.get("MolecularFormula"))
    if not formula or formula != conflict["formula"]:
        return None, "pubchem_formula_mismatch"
    inchikey = str(row.get("InChIKey") or "").upper()
    if not FULL_INCHIKEY.fullmatch(inchikey):
        return None, "pubchem_full_inchikey_unavailable"
    smiles = str(row.get("SMILES") or "")
    molecule = Chem.MolFromSmiles(smiles) if smiles else None
    if molecule is None:
        return None, "pubchem_smiles_unparseable"
    derived = Chem.MolToInchiKey(molecule)
    if derived != inchikey:
        return None, "rdkit_full_inchikey_mismatch"
    return {
        "name_key": conflict["name_key"],
        "query_name": conflict["query_name"],
        "formula": formula,
        "pubchem_cid": int(row["CID"]),
        "pubchem_title": str(row.get("Title") or ""),
        "inchikey": inchikey,
        "connectivity_inchikey": inchikey[:14],
        "isomeric_smiles": smiles,
        "connectivity_smiles": str(row.get("ConnectivitySMILES") or ""),
        "inchi": str(row.get("InChI") or ""),
        "identity_source": "pubchem_exact_name_formula_rdkit_verified",
        "property_url": row["property_url"],
        "observed_connectivity_identities": sorted(
            conflict["observed_connectivity_identities"]
        ),
        "source_artifacts": sorted(conflict["source_artifacts"]),
        "resolution_triggers": sorted(conflict.get("triggers") or []),
    }, "accepted"


def load_base_registry(path: Path | None) -> tuple[dict[str, Any], str | None]:
    if path is None:
        return {}, None
    resolved = path.expanduser().resolve()
    payload = json.loads(resolved.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("entries"), list):
        raise ValueError("base registry is invalid")
    return payload, hashlib.sha256(resolved.read_bytes()).hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(
        payload, ensure_ascii=False, indent=2, sort_keys=True
    ).encode("utf-8") + b"\n"
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_bytes(encoded)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "artifacts", nargs="*", type=Path,
        help=(
            "V11 flagship artifacts and/or immutable *_results.json files; "
            "defaults to output/v11_flagship/mix33-*"
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--base-registry",
        type=Path,
        help="Optional frozen registry whose accepted entries are retained.",
    )
    parser.add_argument("--timeout", type=float, default=20.0)
    args = parser.parse_args()
    paths = [path.resolve() for path in args.artifacts] or sorted(
        (ROOT / "output" / "v11_flagship").glob("mix33-*_v11_flagship.json")
    )
    if not paths or any(not path.is_file() for path in paths):
        parser.error("one or more input artifacts are missing")

    conflicts = collect_conflicts(paths)
    base_registry, base_registry_sha256 = load_base_registry(args.base_registry)
    entries_by_key = {
        (
            strict_name_key(entry.get("name_key") or entry.get("query_name")),
            normalize_formula_text(entry.get("formula")),
        ): dict(entry)
        for entry in (base_registry.get("entries") or [])
        if strict_name_key(entry.get("name_key") or entry.get("query_name"))
        and normalize_formula_text(entry.get("formula"))
    }
    exclusions = []
    for key, conflict in sorted(conflicts.items()):
        existing = entries_by_key.get(key)
        if existing is not None:
            existing["observed_connectivity_identities"] = sorted(
                set(existing.get("observed_connectivity_identities") or [])
                | set(conflict["observed_connectivity_identities"])
            )
            existing["source_artifacts"] = sorted(
                set(existing.get("source_artifacts") or [])
                | set(conflict["source_artifacts"])
            )
            existing["resolution_triggers"] = sorted(
                set(existing.get("resolution_triggers") or [])
                | set(conflict.get("triggers") or [])
            )
            continue
        try:
            row = fetch_pubchem(conflict["query_name"], args.timeout)
            entry, reason = validate_resolution(conflict, row)
        except (urllib.error.URLError, TimeoutError, ValueError, KeyError) as exc:
            entry, reason = None, f"{type(exc).__name__}: {exc}"
        if entry is not None:
            entries_by_key[key] = entry
        else:
            exclusions.append({
                "name_key": conflict["name_key"],
                "query_name": conflict["query_name"],
                "formula": conflict["formula"],
                "reason": reason,
            })
        time.sleep(0.22)

    artifact_hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths
    }
    payload = {
        "schema_version": "1.0.0",
        "registry_id": "mssystem-v14-pubchem-name-formula-identity-repair-2.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": "PubChem PUG REST exact submitted name",
        "selection_contract": {
            "ground_truth_inputs_used": False,
            "triggers": [
                "same exact normalized name and formula bound to multiple candidate connectivities",
                "reported InChIKey disagrees with ranker-standardized stored structure connectivity",
            ],
            "formula_must_match": True,
            "rdkit_must_reproduce_full_inchikey": True,
            "unresolved_or_ambiguous_is_neutral": True,
        },
        "base_registry": {
            "path": str(args.base_registry.expanduser().resolve())
            if args.base_registry is not None else None,
            "sha256": base_registry_sha256,
            "retained_entry_count": len(base_registry.get("entries") or []),
        },
        "source_artifact_sha256": artifact_hashes,
        "entries": sorted(
            entries_by_key.values(),
            key=lambda item: (item["name_key"], item["formula"]),
        ),
        "exclusions": sorted(exclusions, key=lambda item: (item["name_key"], item["formula"])),
    }
    output = args.output.resolve()
    atomic_json(output, payload)
    print(json.dumps({
        "output": str(output),
        "conflicts": len(conflicts),
        "accepted": len(entries_by_key),
        "excluded": len(exclusions),
        "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
