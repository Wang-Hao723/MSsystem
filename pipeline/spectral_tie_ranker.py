"""Query-stage shared spectral support; normalized library code stays immutable.

Only the scan/entity aggregators differ from the accepted base implementation.
Their code is bound to a private namespace; no global matching module is patched.
The borrowed helpers are tied to an exact source version and must be reviewed
again when that base changes. Original scan results and library records are inputs.
"""
from __future__ import annotations

import hashlib
import importlib
from pathlib import Path
from types import CodeType, FunctionType

ALGORITHM_ID = "mssystem-v16-shared-spectral-ties-1.0.2"
BASE_SOURCE_SHA256 = "810715c682a7ce468e6563e152edb6f1197081626eded0d7b6a51fb3b2c22023"
LOADED_IMPLEMENTATION_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _assert_loaded_functions(source, filename, functions):
    compiled = compile(source, filename, "exec", dont_inherit=True)
    expected = {value.co_name: value for value in compiled.co_consts
                if isinstance(value, CodeType)}
    for function in functions:
        # Code equality compares bytecode/constants recursively. Marshal bytes
        # also encode incidental string interning, which is not code identity.
        if function.__code__ != expected[function.__name__]:
            raise ValueError("Loaded spectral ranking code differs from its source binding")


def implementation_sha256():
    source = Path(__file__).read_bytes()
    if hashlib.sha256(source).hexdigest() != LOADED_IMPLEMENTATION_SHA256:
        raise ValueError("Spectral tie adapter source changed; reload the reviewed implementation")
    _assert_loaded_functions(source, __file__, [
        _ranker_namespace, build_flagship_result, build_v15_flagship_result,
        _candidate_rows_shared_ties, _finalize_entity_shared_ties,
        _original_base_entity_audit, _restore_original_legacy_fields,
    ])
    return LOADED_IMPLEMENTATION_SHA256


def _ranker_namespace():
    implementation_sha256()
    base = importlib.import_module("matching.flagship_ranker")
    source = Path(base.__file__).read_bytes()
    actual = hashlib.sha256(source).hexdigest()
    if actual != BASE_SOURCE_SHA256:
        raise ValueError("Spectral tie adapter requires a reviewed base ranker version")
    _assert_loaded_functions(source, base.__file__, [
        base._candidate_rows, base._finalize_entity, base.build_flagship_result,
    ])
    namespace = dict(vars(base))
    namespace["_original_base_candidate_rows"] = base._candidate_rows
    namespace["_original_base_finalize_entity"] = base._finalize_entity
    for name, implementation in [
        ("_original_base_entity_audit", _original_base_entity_audit),
        ("_restore_original_legacy_fields", _restore_original_legacy_fields),
        ("_candidate_rows", _candidate_rows_shared_ties),
        ("_finalize_entity", _finalize_entity_shared_ties),
    ]:
        bound = FunctionType(implementation.__code__, namespace, name,
                             implementation.__defaults__, implementation.__closure__)
        bound.__kwdefaults__ = implementation.__kwdefaults__
        namespace[name] = bound
    build = FunctionType(base.build_flagship_result.__code__, namespace,
                         "build_flagship_result", base.build_flagship_result.__defaults__,
                         base.build_flagship_result.__closure__)
    build.__kwdefaults__ = base.build_flagship_result.__kwdefaults__
    namespace["build_flagship_result"] = build
    return namespace


def build_flagship_result(task_id, rows, *, policy=None):
    base = importlib.import_module("matching.flagship_ranker")
    selected = policy if policy is not None else base.load_flagship_policy()
    config = selected.get("spectral_tie_handling") or {}
    if config.get("enabled") is not True:
        return base.build_flagship_result(task_id, rows, policy=selected)
    expected = {
        "algorithm_id": ALGORITHM_ID,
        "scan_score_decimal_places": 8,
        "ce_support_decimal_places": 6,
        "rank_support": "mean_over_tied_ordinal_positions",
        "top1_support": "equal_share_among_tied_leaders",
        "identity_key_affects_evidence_support": False,
        "library_records_changed": False,
    }
    if any(config.get(key) != value for key, value in expected.items()):
        raise ValueError("Spectral tie policy differs from the validated implementation")
    result = _ranker_namespace()["build_flagship_result"](task_id, rows, policy=selected)
    result["metadata"].update({
        "spectral_tie_handling_algorithm_id": ALGORITHM_ID,
        "spectral_tie_base_source_sha256": BASE_SOURCE_SHA256,
        "spectral_tie_implementation_sha256": implementation_sha256(),
        "spectral_tie_handler_phase": "query_entity_aggregation",
    })
    return result


def build_v15_flagship_result(task_id, rows, *, policy=None):
    from matching.v15_flagship_ranker import load_v15_flagship_policy
    selected = policy if policy is not None else load_v15_flagship_policy()
    return build_flagship_result(task_id, rows, policy=selected)


def artifact_implementation_is_current(metadata, policy):
    if (policy.get("spectral_tie_handling") or {}).get("enabled") is not True:
        return True
    from pipeline.spectral_tie_scope import (
        ALGORITHM_ID as scope_algorithm_id,
        implementation_sha256 as scope_implementation_sha256,
    )
    return bool(
        metadata.get("spectral_tie_handling_algorithm_id") == ALGORITHM_ID
        and metadata.get("spectral_tie_base_source_sha256") == BASE_SOURCE_SHA256
        and metadata.get("spectral_tie_implementation_sha256") == implementation_sha256()
        and metadata.get("spectral_tie_scope_guard_algorithm_id") == scope_algorithm_id
        and metadata.get("spectral_tie_scope_guard_source_sha256") == scope_implementation_sha256()
    )


def _candidate_rows_shared_ties(row: Mapping[str, Any], policy: Mapping[str, Any]) -> list[dict[str, Any]]:
    original_base_candidates = _original_base_candidate_rows(row, policy)
    row, _ = _restore_reference_mass_scores(row)
    match = row.get("match") if isinstance(row.get("match"), Mapping) else None
    if not match:
        return []
    candidates: list[Mapping[str, Any]] = [match]
    candidates.extend(
        candidate for candidate in (match.get("top_candidates") or [])
        if isinstance(candidate, Mapping)
    )
    maximum = max(1, int((policy.get("candidate_pool") or {}).get("maximum_candidates_per_scan") or 20))
    standardize_tautomers = bool(
        ((policy.get("candidate_pool") or {}).get("tautomer_standardization") or {}).get(
            "enabled", True
        )
    )
    by_identity: dict[str, dict[str, Any]] = {}
    references_by_identity: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for candidate in candidates:
        registry_block, registry_entry = _registered_identity(candidate, policy)
        resolved_candidate = dict(candidate)
        if registry_entry is not None:
            reported_inchikey = _text(candidate.get("inchikey"))
            reported_smiles = _text(candidate.get("smiles"))
            resolved_candidate["inchikey"] = registry_entry.get("inchikey")
            resolved_candidate["smiles"] = registry_entry.get("isomeric_smiles")
            identity = f"std2d:{registry_block}"
            identity_source = "pubchem_exact_name_formula_registry"
        else:
            reported_inchikey = ""
            reported_smiles = ""
            identity, identity_source = connectivity_identity(
                candidate, standardize_tautomers=standardize_tautomers
            )
        spectral, components = _spectral_score(candidate, policy)
        compact_references = [
            reference
            for reference in (candidate.get("reference_evidence_records") or [])
            if isinstance(reference, Mapping)
        ]
        reference_inputs = compact_references + [candidate]
        for reference in reference_inputs:
            reference_spectral, reference_components = _spectral_score(
                reference, policy
            )
            references_by_identity[identity].append({
                # Names and biological annotations are intentionally not
                # projected into the V15 reference aggregator.
                "spectral_score": reference_spectral,
                "score_components": reference_components,
                "dataset_source": reference.get("dataset_source"),
                "source_tier": reference.get("source_tier"),
                "lib_spectrum_id": (
                    reference.get("lib_spectrum_id")
                    or reference.get("spectrum_id")
                ),
                "lib_source_accession": (
                    reference.get("lib_source_accession")
                    or reference.get("source_accession")
                ),
                "library_peak_fingerprint": reference.get(
                    "library_peak_fingerprint"
                ),
                "lib_peaks": (
                    reference.get("lib_peaks")
                    or reference.get("_score_library_peaks")
                    or (
                        row.get("lib_peaks")
                        if reference is candidate and candidate is match
                        else []
                    )
                ),
                "lib_collision_energy": reference.get(
                    "lib_collision_energy"
                ),
                "lib_collision_energy_unit": reference.get(
                    "lib_collision_energy_unit"
                ),
                "lib_collision_energy_mode": reference.get(
                    "lib_collision_energy_mode"
                ),
                "lib_instrument_type": reference.get("lib_instrument_type"),
                "lib_instrument": reference.get("lib_instrument"),
                "lib_activation_method": reference.get(
                    "lib_activation_method"
                ),
                "lib_instrument_domain_id": reference.get(
                    "lib_instrument_domain_id"
                ),
                "lib_collision_cell_region": reference.get(
                    "lib_collision_cell_region"
                ),
                "_reference_projection": (
                    "compact_search_evidence"
                    if reference is not candidate
                    else "candidate_fallback"
                ),
            })
        display = _display_candidate(
            resolved_candidate,
            standardize_tautomers=standardize_tautomers,
        )
        if registry_entry is not None:
            ranker_standardized_before_resolution = _text(
                display.get("standardized_connectivity_inchikey")
            )
            display["standardized_connectivity_inchikey"] = registry_block
            display["reported_connectivity_matches_standardized"] = bool(
                reported_inchikey
                and reported_inchikey.upper().split("-", 1)[0] == registry_block
            )
            if (
                ranker_standardized_before_resolution
                and ranker_standardized_before_resolution != registry_block
            ):
                display[
                    "ranker_standardized_connectivity_before_identity_resolution"
                ] = ranker_standardized_before_resolution
            display.update({
                "reported_inchikey_before_identity_resolution": reported_inchikey or None,
                "reported_smiles_before_identity_resolution": reported_smiles or None,
                "canonical_identity_resolution": {
                    "status": "resolved",
                    "registry_id": _identity_registry(
                        str((policy.get("identity_registry") or {}).get("sha256") or "")
                    )[0],
                    "source": registry_entry.get("identity_source"),
                    "pubchem_cid": registry_entry.get("pubchem_cid"),
                    "pubchem_title": registry_entry.get("pubchem_title"),
                    "exact_normalized_name_and_formula": True,
                    "ground_truth_inputs_used": False,
                },
            })
        item = {
            "identity_key": identity,
            "identity_source": identity_source,
            "candidate": display,
            "spectral_score": spectral,
            "components": components,
            "reference_source_families": sorted(filter(None, {
                _reference_source_family(candidate)
            })),
        }
        previous = by_identity.get(identity)
        if previous is None:
            by_identity[identity] = item
            continue
        previous_display = dict(previous["candidate"])
        item_display = dict(item["candidate"])
        merged_sources = sorted(set(previous.get("reference_source_families") or []) | set(
            item.get("reference_source_families") or []
        ))
        preferred_display = min(
            (
                (previous_display, float(previous["spectral_score"])),
                (item_display, spectral),
            ),
            key=lambda option: (
                -int(bool(option[0].get("canonical_identity_resolution"))),
                -int(bool(option[0].get("reported_connectivity_matches_standardized"))),
                -int(bool(option[0].get("has_structure"))),
                -float(option[1]),
                _text(option[0].get("inchikey")),
                _text(option[0].get("compound_name")),
            ),
        )[0]
        if (
            -spectral,
            _text(item_display.get("inchikey")),
            _text(item_display.get("compound_name")),
        ) < (
            -float(previous["spectral_score"]),
            _text(previous_display.get("inchikey")),
            _text(previous_display.get("compound_name")),
        ):
            item["candidate"] = preferred_display
            item["reference_source_families"] = merged_sources
            by_identity[identity] = item
        else:
            previous["candidate"] = preferred_display
            previous["reference_source_families"] = merged_sources

    for identity, records in references_by_identity.items():
        references_by_identity[identity] = _merge_projected_references(records)

    multilibrary_policy = policy.get("conditioned_multilibrary") or {}
    multilibrary_enabled = bool(multilibrary_policy.get("enabled"))
    if multilibrary_enabled:
        for identity, item in by_identity.items():
            evidence = aggregate_references(
                row,
                references_by_identity[identity],
                multilibrary_policy,
            )
            legacy_score = float(item["spectral_score"])
            item["legacy_spectral_score"] = round(legacy_score, 8)
            item["spectral_score"] = float(
                evidence["conditioned_multilibrary_score"]
            )
            item["conditioned_multilibrary_evidence"] = evidence
            item["reference_source_families"] = list(
                evidence.get("source_families") or []
            )

        legacy_ordered = sorted(
            by_identity.values(),
            key=lambda item: (
                -float(item.get("legacy_spectral_score") or 0.0),
                item["identity_key"],
            ),
        )[:maximum]
        rank_fraction = float(
            (policy.get("scan_score") or {}).get("reciprocal_rank_fraction")
            or 0.30
        )
        spectral_fraction = float(
            (policy.get("scan_score") or {}).get("spectral_fraction") or 0.70
        )
        for rank, item in enumerate(legacy_ordered, start=1):
            reciprocal = 1.0 / math.log2(rank + 1.0)
            item["legacy_scan_rank"] = rank
            item["legacy_scan_support"] = round(
                max(
                    0.0,
                    min(
                        1.0,
                        spectral_fraction
                        * float(item.get("legacy_spectral_score") or 0.0)
                        + rank_fraction * reciprocal,
                    ),
                ),
                6,
            )
    ordered = sorted(
        by_identity.values(),
        key=lambda item: (-round(float(item["spectral_score"]), 8), item["identity_key"]),
    )
    rank_fraction = float((policy.get("scan_score") or {}).get("reciprocal_rank_fraction", 0.30))
    spectral_fraction = float((policy.get("scan_score") or {}).get("spectral_fraction", 0.70))
    tie_groups: dict[float, list[dict[str, Any]]] = defaultdict(list)
    for item in ordered:
        # The scorer/reference aggregator already exports eight decimal places.
        # This is equality at that recorded precision, not a tuned near-tie gate.
        tie_groups[round(float(item["spectral_score"]), 8)].append(item)
    first_rank = 1
    for spectral, tied in tie_groups.items():
        last_rank = first_rank + len(tied) - 1
        reciprocal = math.fsum(
            1.0 / math.log2(rank + 1.0)
            for rank in range(first_rank, last_rank + 1)
        ) / len(tied)
        members = [
            {"identity_key": item["identity_key"],
             "formula": normalize_formula_text(item["candidate"].get("formula"))}
            for item in tied
        ]
        retained = max(0, min(last_rank, maximum) - first_rank + 1)
        for item in tied:
            item["spectral_tie_count"] = len(tied)
            item["spectral_tie_rank_start"] = first_rank
            item["spectral_tie_rank_end"] = last_rank
            item["spectral_tie_members"] = members if len(tied) > 1 else []
            item["spectral_tie_members_truncated"] = len(tied) - retained
            item["reciprocal_rank_support"] = round(reciprocal, 6)
            item["scan_top1_share"] = 1.0 / len(tied) if first_rank == 1 else 0.0
            item["scan_support"] = round(
                max(0.0, min(1.0, spectral_fraction * spectral + rank_fraction * reciprocal)),
                6,
            )
        first_rank = last_rank + 1
    for rank, item in enumerate(ordered, start=1):
        # This ordinal is display-only. Evidence support and first-place votes
        # above are shared over the complete tie group, before scan truncation.
        item["scan_rank"] = rank
        if not multilibrary_enabled:
            item["legacy_spectral_score"] = round(
                float(item["spectral_score"]), 8
            )
            item["legacy_scan_rank"] = rank
            item["legacy_scan_support"] = round(
                max(0.0, min(1.0, spectral_fraction * float(item["spectral_score"])
                    + rank_fraction / math.log2(rank + 1.0))), 6,
            )
    selected = ordered[:maximum]
    if selected:
        selected[0]["_original_base_projection_for_legacy_audit"] = original_base_candidates
    return selected


def _finalize_entity_shared_ties(entity: Mapping[str, Any], policy: Mapping[str, Any]) -> dict[str, Any]:
    original_base_audit = _original_base_entity_audit(entity, policy)
    pool_policy = policy.get("candidate_pool") or {}
    replicate_rule = pool_policy.get("same_energy_replicate_support", "legacy_zero_fill")
    if replicate_rule not in {"legacy_zero_fill", "observed_candidates_only_v16"}:
        raise ValueError("Unknown same-energy replicate aggregation policy")
    topk_absence_is_negative = bool(pool_policy.get("missing_from_truncated_topk_is_negative_evidence", False))
    omit_unobserved_replicates = replicate_rule == "observed_candidates_only_v16" and not topk_absence_is_negative
    scans = sorted(
        list(entity.get("scans") or []),
        key=lambda scan: (
            -float(scan["row"].get("tic") or 0.0),
            str(scan["row"].get("spectrum_id") or ""),
            int(scan["row"]["result_index"]),
        ),
    )
    bins: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for scan in scans:
        bins[str(scan["ce_key"])].append(scan)
    bin_keys = sorted(bins)
    candidate_ce: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "bins": {},
            "instances": [],
            "variants": set(),
            "reference_source_families": set(),
        }
    )
    ce_winners: dict[str, dict[str, float]] = {}
    legacy_ce_winners: dict[str, str] = {}

    for ce_key in bin_keys:
        replicate_scans = bins[ce_key]
        identities = sorted({
            item["identity_key"]
            for scan in replicate_scans
            for item in scan.get("candidates") or []
        })
        ce_records = []
        legacy_ce_records = []
        # A scan limit can cut through a tie. Its already computed equal-score
        # members remain observed for vote accounting, without adding display
        # candidates or reconstructing missing reference spectra.
        truncated_tie_supports = []
        for scan in replicate_scans:
            observed = {}
            for item in scan.get("candidates") or []:
                if item.get("spectral_tie_members_truncated", 0):
                    for member in item.get("spectral_tie_members") or []:
                        observed[member["identity_key"]] = (
                            float(item["scan_support"]), float(item["scan_top1_share"])
                        )
            truncated_tie_supports.append(observed)
        for identity in identities:
            supports = []
            legacy_supports = []
            instances = []
            top1_count = 0
            legacy_top1_count = 0
            tied_scores_outside_limit = 0
            for scan, truncated in zip(replicate_scans, truncated_tie_supports):
                item = next(
                    (candidate for candidate in scan.get("candidates") or [] if candidate["identity_key"] == identity),
                    None,
                )
                supports.append(
                    float(item.get("scan_support") or 0.0) if item else
                    truncated[identity][0] if identity in truncated else None
                )
                legacy_supports.append(
                    float(item.get("legacy_scan_support") or 0.0)
                    if item else 0.0
                )
                if item:
                    instances.append((float(item.get("scan_support") or 0.0), scan, item))
                    top1_count += float(item.get(
                        "scan_top1_share",
                        int(int(item.get("scan_rank") or 999999) == 1),
                    ))
                    legacy_top1_count += int(
                        int(item.get("legacy_scan_rank") or 999999) == 1
                    )
                    candidate_ce[identity]["variants"].add(_candidate_variant(item["candidate"]))
                    candidate_ce[identity]["reference_source_families"].update(
                        item.get("reference_source_families") or []
                    )
                elif identity in truncated:
                    top1_count += truncated[identity][1]
                    tied_scores_outside_limit += 1
            ce_support = _replicate_support_median(supports, omit_unobserved=omit_unobserved_replicates)
            legacy_ce_support = (
                statistics.median(legacy_supports) if legacy_supports else 0.0
            )
            representative = max(
                instances,
                key=lambda value: (value[0], -int(value[1]["row"]["result_index"])),
                default=None,
            )
            record = {
                "collision_energy_bin": ce_key,
                "support": round(ce_support, 6),
                "legacy_v14_support": round(legacy_ce_support, 6),
                "technical_replicates": len(replicate_scans),
                "replicates_with_candidate": len(instances),
                "replicates_with_tied_score_outside_scan_limit": tied_scores_outside_limit,
                "replicates_missing_from_truncated_topk": len(replicate_scans) - len(instances) - tied_scores_outside_limit,
                "truncated_topk_absence_imputed_as_zero": not omit_unobserved_replicates,
                "technical_replicates_count_as_independent_votes": False,
                "scan_top1_fraction": round(top1_count / len(replicate_scans), 6) if replicate_scans else 0.0,
                "legacy_v14_scan_top1_fraction": (
                    round(legacy_top1_count / len(replicate_scans), 6)
                    if replicate_scans else 0.0
                ),
                "representative_result_index": (
                    int(representative[1]["row"]["result_index"]) if representative else None
                ),
                "conditioned_multilibrary_evidence": (
                    representative[2].get("conditioned_multilibrary_evidence")
                    if representative else None
                ),
            }
            candidate_ce[identity]["bins"][ce_key] = record
            if representative:
                candidate_ce[identity]["instances"].append(representative)
            ce_records.append((ce_support, identity))
            legacy_ce_records.append((legacy_ce_support, identity))
        if ce_records:
            # CE support is recorded to six decimal places. A lexical display
            # winner cannot receive an extra evidence vote when support ties.
            vote_records = list(ce_records)
            omitted_identities = set().union(*(set(values) for values in truncated_tie_supports)) - set(identities)
            for identity in sorted(omitted_identities):
                supports = [values[identity][0] if identity in values else None
                            for values in truncated_tie_supports]
                vote_records.append((_replicate_support_median(
                    supports, omit_unobserved=omit_unobserved_replicates), identity))
            best_support = max(round(value[0], 6) for value in vote_records)
            tied = [identity for support, identity in vote_records
                    if round(support, 6) == best_support]
            ce_winners[ce_key] = {identity: 1.0 / len(tied) for identity in tied}
            for _support, identity in ce_records:
                candidate_ce[identity]["bins"][ce_key]["bin_top1_share"] = (
                    ce_winners[ce_key].get(identity, 0.0)
                )
        if legacy_ce_records:
            legacy_ce_winners[ce_key] = min(
                legacy_ce_records, key=lambda value: (-value[0], value[1])
            )[1]

    weights = policy.get("entity_score") or {}
    candidates = []
    total_bins = len(bin_keys)
    for identity, data in sorted(candidate_ce.items()):
        supports = [float((data["bins"].get(key) or {}).get("support") or 0.0) for key in bin_keys]
        legacy_supports = [
            float(
                (data["bins"].get(key) or {}).get("legacy_v14_support")
                or 0.0
            )
            for key in bin_keys
        ]
        # The versioned V16 rule also retains a genuinely measured zero. Only
        # a CE bin with no observed candidate is absent; score > 0 is not an
        # observation test. Positive support still controls coverage below.
        observed_supports = (
            [float(data["bins"][key]["support"]) for key in bin_keys if key in data["bins"]]
            if omit_unobserved_replicates else [support for support in supports if support > 0]
        )
        legacy_observed_supports = [
            support for support in legacy_supports if support > 0
        ]
        observed_bins = sum(support > 0 for support in observed_supports)
        # The scan artifact exposes a truncated Top-K, not an exhaustive score
        # vector.  Absence outside that list is therefore missing evidence, not
        # a measured zero.  Coverage remains an explicit score component and a
        # release-scope constraint, so sparse support cannot masquerade as a
        # replicated multi-energy identification.
        support_basis = supports if topk_absence_is_negative else observed_supports
        legacy_support_basis = (
            legacy_supports
            if topk_absence_is_negative
            else legacy_observed_supports
        )
        mean_support = statistics.fmean(support_basis) if support_basis else 0.0
        median_support = statistics.median(support_basis) if support_basis else 0.0
        legacy_mean_support = (
            statistics.fmean(legacy_support_basis)
            if legacy_support_basis else 0.0
        )
        legacy_median_support = (
            statistics.median(legacy_support_basis)
            if legacy_support_basis else 0.0
        )
        coverage = observed_bins / total_bins if total_bins else 0.0
        # The embedded historical lane retains the old padded-replicate result,
        # including its coverage. A V16 correction must not rewrite that audit.
        legacy_coverage = (
            len(legacy_observed_supports) / total_bins if total_bins else 0.0
        ) if omit_unobserved_replicates else coverage
        top1_fraction = (
            math.fsum(ce_winners.get(key, {}).get(identity, 0.0) for key in bin_keys) / total_bins
            if total_bins else 0.0
        )
        legacy_top1_fraction = (
            sum(legacy_ce_winners.get(key) == identity for key in bin_keys)
            / total_bins
            if total_bins else 0.0
        )
        source_families = sorted(data["reference_source_families"])
        source_saturation = max(
            1,
            int((policy.get("reference_provenance") or {}).get(
                "source_diversity_saturation", 3
            )),
        )
        source_diversity = min(1.0, len(source_families) / source_saturation)
        entity_score = (
            float(weights.get("mean_collision_energy_support", 0.45)) * mean_support
            + float(weights.get("median_collision_energy_support", 0.30)) * median_support
            + float(weights.get("collision_energy_coverage", 0.15)) * coverage
            + float(weights.get("collision_energy_top1_fraction", 0.10)) * top1_fraction
        )
        legacy_entity_score = (
            float(weights.get("mean_collision_energy_support") or 0.45)
            * legacy_mean_support
            + float(weights.get("median_collision_energy_support") or 0.30)
            * legacy_median_support
            + float(weights.get("collision_energy_coverage") or 0.15)
            * legacy_coverage
            + float(weights.get("collision_energy_top1_fraction") or 0.10)
            * legacy_top1_fraction
        )
        instances = sorted(
            data["instances"],
            key=lambda value: (-value[0], int(value[1]["row"]["result_index"])),
        )
        checks = [
            _formula_adduct_check(item[2]["candidate"], item[1]["row"], policy)
            for item in instances
        ]
        chemistry = _aggregate_chemistry(checks)
        display = dict(instances[0][2]["candidate"]) if instances else {}
        candidate = {
            "identity_key": identity,
            "identity_scope": "2d_connectivity",
            "identity_source": instances[0][2].get("identity_source") if instances else "unknown",
            "candidate": display,
            "entity_score": round(max(0.0, min(1.0, entity_score)), 6),
            "legacy_v14_entity_score": round(
                max(0.0, min(1.0, legacy_entity_score)), 6
            ),
            "score_is_probability": False,
            "score_components": {
                "mean_collision_energy_support": round(mean_support, 6),
                "median_collision_energy_support": round(median_support, 6),
                "collision_energy_coverage": round(coverage, 6),
                "collision_energy_top1_fraction": round(top1_fraction, 6),
                "topk_absence_imputed_as_zero": topk_absence_is_negative,
                "reference_source_diversity": round(source_diversity, 6),
                "reference_source_diversity_affects_score": False,
                "source_count_bonus_applied": False,
            },
            "reference_source_family_count": len(source_families),
            "reference_source_families": source_families,
            "reference_source_accessions_count_as_independent": False,
            "collision_energy_bins_supported": observed_bins,
            "collision_energy_bins_total": total_bins,
            "collision_energy_evidence": [data["bins"][key] for key in bin_keys if key in data["bins"]],
            "representative_result_index": (
                int(instances[0][1]["row"]["result_index"]) if instances else None
            ),
            "representative_scan_score": round(float(instances[0][0]), 6) if instances else None,
            "chemistry_gate": chemistry,
            "vetoed": bool(chemistry.get("veto")),
            "variant_count_collapsed_to_connectivity": len(data["variants"]),
            "stereochemistry_claimed": False,
            "ranking_inputs": [
                "cosine",
                "entropy",
                "reverse_match",
                "mass_accuracy",
                "fragment_coverage",
                "intensity_agreement",
                "within_scan_rank",
                "collision_energy_consistency",
            ] + (
                [
                    "mirror_collapsed_reference_profiles",
                    "provenance_stratum_median_similarity",
                    "compatible_instrument_domain_near_ce_preference",
                ]
                if (policy.get("conditioned_multilibrary") or {}).get("enabled")
                else []
            ),
            "conditioned_multilibrary_runtime_influence": any(
                bool(
                    (
                        (data["bins"].get(key) or {}).get(
                            "conditioned_multilibrary_evidence"
                        )
                        or {}
                    ).get("runtime_influence")
                )
                for key in bin_keys
            ),
        }
        candidate["_representative_row"] = instances[0][1]["row"] if instances else None
        candidates.append(candidate)

    legacy_order = sorted(
        candidates,
        key=lambda candidate: (
            bool(candidate["vetoed"]),
            -float(candidate["legacy_v14_entity_score"]),
            candidate["identity_key"],
        ),
    )
    legacy_rank_by_identity = {
        candidate["identity_key"]: rank
        for rank, candidate in enumerate(legacy_order, start=1)
    }
    for candidate in candidates:
        candidate["legacy_v14_entity_rank"] = legacy_rank_by_identity[
            candidate["identity_key"]
        ]
    candidates.sort(
        key=lambda candidate: (
            bool(candidate["vetoed"]),
            -float(candidate["entity_score"]),
            candidate["identity_key"],
        )
    )
    valid = [candidate for candidate in candidates if not candidate["vetoed"]]
    decision_policy = policy.get("decision") or {}
    minimum = float(decision_policy.get("minimum_entity_score") or 0.50)
    margin_required = float(decision_policy.get("minimum_connectivity_margin") or 0.04)
    minimum_ce = int(decision_policy.get("minimum_distinct_collision_energies_for_connectivity") or 2)
    minimum_supported_ce = int(decision_policy.get("minimum_winner_collision_energy_support") or 2)
    isomer_margin = float(decision_policy.get("isomer_group_margin") or 0.05)

    leading = valid[0] if valid else None
    runner = valid[1] if len(valid) > 1 else None
    margin = (
        float(leading["entity_score"]) - float(runner["entity_score"])
        if leading and runner else float(leading["entity_score"]) if leading else 0.0
    )
    winner = leading if leading and float(leading["entity_score"]) >= minimum else None
    truncated_top_spectral_ties = []
    for scan in scans:
        item = next((item for item in scan.get("candidates") or []
                     if item.get("spectral_tie_rank_start") == 1
                     and item.get("spectral_tie_members_truncated", 0) > 0), None)
        if item is not None:
            truncated_top_spectral_ties.append({
                "result_index": scan["row"]["result_index"],
                "collision_energy_bin": scan["ce_key"],
                "members": item["spectral_tie_members"],
                "members_truncated": item["spectral_tie_members_truncated"],
            })
    if winner is None:
        decision_status = "unmatched" if not candidates else "abstained"
        decision_scope = "none"
        blockers = ["no_library_candidate"] if not candidates else ["entity_score_below_minimum"]
        if candidates and not valid:
            blockers = ["all_candidates_failed_candidate_independent_chemistry_gate"]
    else:
        transformation_ambiguity_candidates = [
            candidate
            for candidate in valid
            if candidate.get("chemistry_gate", {}).get(
                "precursor_transformation_ambiguity"
            )
            and float(candidate.get("entity_score") or 0.0) >= minimum
        ]
        same_formula_competitors = [
            candidate
            for candidate in valid[1:]
            if normalize_formula_text(candidate["candidate"].get("formula"))
            and normalize_formula_text(candidate["candidate"].get("formula"))
            == normalize_formula_text(winner["candidate"].get("formula"))
            and float(winner["entity_score"]) - float(candidate["entity_score"]) <= isomer_margin
        ]
        if transformation_ambiguity_candidates:
            decision_status = "selected_candidate_only"
            decision_scope = "precursor_transformation_ambiguity_candidate"
            blockers = [
                "candidate_independent_ion_identity_conflicts_with_"
                "mass_consistent_library_neutral_loss_adduct"
            ]
        elif same_formula_competitors:
            decision_status = "selected_isomer_group"
            decision_scope = "molecular_formula_isomer_group"
            blockers = ["competitive_same_formula_connectivity"]
        elif truncated_top_spectral_ties:
            decision_status = "selected_candidate_only"
            decision_scope = "truncated_spectral_tie_candidate"
            blockers = ["equally_supported_scan_candidates_outside_candidate_limit"]
        elif (
            total_bins >= minimum_ce
            and int(winner["collision_energy_bins_supported"]) >= minimum_supported_ce
            and margin >= margin_required
        ):
            decision_status = "selected"
            decision_scope = "2d_connectivity_candidate"
            blockers = []
        elif total_bins < minimum_ce:
            decision_status = "selected_candidate_only"
            decision_scope = "single_collision_energy_candidate"
            blockers = ["insufficient_distinct_collision_energies"]
        else:
            decision_status = "selected_candidate_only"
            decision_scope = "ambiguous_connectivity_candidate"
            blockers = ["insufficient_entity_margin_or_energy_coverage"]

    base_formal_release_eligible = bool(
        winner is not None
        and decision_status == "selected"
        and decision_scope == "2d_connectivity_candidate"
    )
    query_quality_config = policy.get("query_spectrum_release_gate") or {}
    query_quality_audit: dict[str, Any] = {
        "enabled": bool(query_quality_config.get("enabled")),
        "algorithm_id": query_quality_config.get("algorithm_id"),
        "status": "not_applicable",
        "failed_gates": [],
        "winner_rank_changed": False,
        "winner_score_changed": False,
        "score_is_probability": False,
    }
    query_quality_failures: list[str] = []
    if (
        base_formal_release_eligible
        and query_quality_config.get("enabled") is True
    ):
        query_quality_audit = _query_spectrum_release_audit(
            scans, query_quality_config
        )
        if query_quality_config.get("affects_release_scope") is True:
            query_quality_failures = list(
                query_quality_audit.get("failed_gates") or []
            )

    abstention_config = policy.get("cross_platform_abstention_gate") or {}
    abstention_audit: dict[str, Any] = {
        "enabled": bool(abstention_config.get("enabled")),
        "algorithm_id": abstention_config.get("algorithm_id"),
        "status": "not_applicable",
        "failed_gates": [],
        "winner_rank_changed": False,
        "winner_score_changed": False,
        "score_is_probability": False,
        "external_independent_validation_complete": bool(
            abstention_config.get("external_independent_validation_complete")
        ),
    }
    if (
        base_formal_release_eligible
        and abstention_config.get("enabled") is True
    ):
        minimum_release_score = float(
            abstention_config.get("minimum_entity_score") or 0.0
        )
        minimum_release_margin = float(
            abstention_config.get("minimum_winner_margin") or 0.0
        )
        minimum_release_ce = int(
            abstention_config.get(
                "minimum_collision_energy_bins_supported"
            ) or 0
        )
        release_failures: list[str] = []
        if float(winner.get("entity_score") or 0.0) < minimum_release_score:
            release_failures.append("entity_score_below_cross_platform_release_gate")
        if margin < minimum_release_margin:
            release_failures.append("winner_margin_below_cross_platform_release_gate")
        if int(winner.get("collision_energy_bins_supported") or 0) < minimum_release_ce:
            release_failures.append("winner_ce_support_below_cross_platform_release_gate")
        abstention_audit.update(
            {
                "status": "demoted" if release_failures else "passed",
                "failed_gates": release_failures,
                "observed_entity_score": round(
                    float(winner.get("entity_score") or 0.0), 6
                ),
                "observed_winner_margin": round(margin, 6),
                "observed_collision_energy_bins_supported": int(
                    winner.get("collision_energy_bins_supported") or 0
                ),
                "minimum_entity_score": minimum_release_score,
                "minimum_winner_margin": minimum_release_margin,
                "minimum_collision_energy_bins_supported": minimum_release_ce,
                "calibration_state": abstention_config.get("calibration_state"),
            }
        )

    abstention_failures = list(abstention_audit.get("failed_gates") or [])
    additional_release_failures = list(
        dict.fromkeys([*query_quality_failures, *abstention_failures])
    )
    if additional_release_failures:
        decision_status = "selected_candidate_only"
        if query_quality_failures and abstention_failures:
            decision_scope = "multi_gate_abstention_candidate"
        elif query_quality_failures:
            decision_scope = "query_spectrum_quality_candidate"
        else:
            decision_scope = "cross_platform_abstention_candidate"
        blockers = additional_release_failures

    if winner is not None:
        winner["fragmentation_explanation"] = _fragmentation_explanation(
            winner["candidate"],
            winner["chemistry_gate"],
            winner.get("_representative_row"),
        )

    _restore_original_legacy_fields(candidates, original_base_audit)
    shortlist_size = max(1, int((policy.get("candidate_pool") or {}).get("shortlist_size") or 5))
    # Preserve the full leading entity tie in the decision even when the
    # deterministic display shortlist ends inside it.
    leading_score_ties = [candidate for candidate in valid
                          if leading is not None
                          and candidate["entity_score"] == leading["entity_score"]]
    leading_tie_identities = [candidate["identity_key"] for candidate in leading_score_ties]
    public_identity_keys = {candidate["identity_key"] for candidate in candidates[:shortlist_size]}
    public_candidates = []
    for rank, candidate in enumerate(candidates[:shortlist_size], start=1):
        clean = {key: value for key, value in candidate.items() if not key.startswith("_")}
        clean["entity_rank"] = rank
        public_candidates.append(clean)
    public_winner = None
    if winner is not None:
        public_winner = next(
            (candidate for candidate in public_candidates if candidate["identity_key"] == winner["identity_key"]),
            {key: value for key, value in winner.items() if not key.startswith("_")},
        )
    public_leading = None
    if leading is not None:
        public_leading = next(
            (candidate for candidate in public_candidates if candidate["identity_key"] == leading["identity_key"]),
            {key: value for key, value in leading.items() if not key.startswith("_")},
        )
    legacy_leading = next(
        (candidate for candidate in legacy_order if not candidate["vetoed"]),
        None,
    )

    # The scan-level comparator is an audit baseline, not an input-order
    # accident.  Use the highest-TIC member and deterministic tie-breakers.
    baseline_scan = scans[0] if scans else None
    baseline_identity = (
        baseline_scan["candidates"][0]["identity_key"]
        if baseline_scan and baseline_scan.get("candidates")
        else None
    )

    precursor_values = [
        float(scan["row"]["search_precursor_mz"])
        for scan in scans
        if scan["row"].get("search_precursor_mz") is not None
    ]
    selected_ion_values = [
        float(scan["row"]["precursor_mz"])
        for scan in scans
        if scan["row"].get("precursor_mz") is not None
    ]
    precursor_sources = sorted({
        _text(scan["row"].get("search_precursor_mz_source"))
        for scan in scans
        if _text(scan["row"].get("search_precursor_mz_source"))
    })
    rt_values = [
        float(scan["row"]["retention_time"])
        for scan in scans
        if scan["row"].get("retention_time") is not None
    ]
    tic_values = [
        float(scan["row"]["tic"])
        for scan in scans
        if scan["row"].get("tic") is not None
    ]
    # These two arrays form a positional binding consumed by EviRecall and
    # EviTwin.  Sorting them independently corrupts that binding whenever the
    # lexical spectrum-id order differs from the immutable result-row order.
    member_pairs = sorted(
        (
            int(scan["row"]["result_index"]),
            str(scan["row"]["spectrum_id"]),
        )
        for scan in scans
    )
    return {
        "entity_id": entity["entity_id"],
        "entity_source": entity["entity_source"],
        "entity_source_value": entity["entity_source_value"],
        "member_scan_count": len(scans),
        "member_result_indices": [index for index, _spectrum_id in member_pairs],
        "member_spectrum_ids": [spectrum_id for _index, spectrum_id in member_pairs],
        "representative_precursor_mz": round(statistics.median(precursor_values), 8) if precursor_values else None,
        "representative_search_precursor_mz": round(statistics.median(precursor_values), 8) if precursor_values else None,
        "representative_selected_ion_mz": round(statistics.median(selected_ion_values), 8) if selected_ion_values else None,
        "search_precursor_mz_sources": precursor_sources,
        "representative_search_precursor_mz_source": (
            precursor_sources[0] if len(precursor_sources) == 1 else
            "mixed" if precursor_sources else None
        ),
        "lcms_feature_used_for_search": any(
            bool(scan["row"].get("lcms_feature_used_for_search"))
            for scan in scans
        ),
        "representative_retention_time": round(statistics.median(rt_values), 8) if rt_values else None,
        "representative_tic": max(tic_values) if tic_values else None,
        "ion_modes": sorted({_text(scan["row"].get("ion_mode")) for scan in scans if _text(scan["row"].get("ion_mode"))}),
        "distinct_collision_energy_bins": total_bins,
        "collision_energy_bins": [
            {
                "key": key,
                "value": bins[key][0].get("ce_value"),
                "unit": bins[key][0].get("ce_unit"),
                "technical_replicates": len(bins[key]),
                "condition_evidence_count": 1,
                "independent_sample_vote_count": 0,
            }
            for key in bin_keys
        ],
        "decision": {
            "status": decision_status,
            "scope": decision_scope,
            "winner_margin": round(margin, 6),
            "minimum_score": minimum,
            "minimum_margin": margin_required,
            "blockers": blockers,
            "query_spectrum_release_gate": query_quality_audit,
            "cross_platform_abstention_gate": abstention_audit,
            "stereochemistry_claimed": False,
            "score_is_probability": False,
            "score_tie_handling": {
                "scan_spectral_decimal_places": 8,
                "collision_energy_support_decimal_places": 6,
                "scan_rank_support": "mean_over_tied_ordinal_positions",
                "top1_votes": "equal_share_among_tied_leaders",
                "identity_key_affects_evidence_score": False,
                "leading_entity_tie_identities": leading_tie_identities,
                "leading_entity_tie_count": len(leading_tie_identities),
                "leading_entity_tie_outside_display": [identity for identity in leading_tie_identities
                                                        if identity not in public_identity_keys],
                "truncated_top_spectral_ties": truncated_top_spectral_ties,
            },
        },
        "winner": public_winner,
        "leading_candidate_when_abstained": public_leading if winner is None else None,
        "candidates": public_candidates,
        "candidate_count": len(candidates),
        "vetoed_candidate_count": sum(bool(candidate["vetoed"]) for candidate in candidates),
        "ranking_changed_from_scan_top1": bool(
            public_winner
            and baseline_identity
            and public_winner["identity_key"] != baseline_identity
        ),
        "ranking_changed_from_v14": bool(
            public_winner
            and (original_base_audit.get("v14_audit_baseline") or {}).get("identity_key")
            and public_winner["identity_key"] != original_base_audit["v14_audit_baseline"]["identity_key"]
        ),
        "v14_audit_baseline": copy.deepcopy(original_base_audit["v14_audit_baseline"]),
        "legacy_v14_audit_integrity": {
            "status": "restored_from_original_base_projection_and_aggregation",
            "original_scan_support_and_top20_used": True,
            "tie_imputed_support_used_for_legacy_audit": False,
            "original_legacy_audit_is_not_the_new_tie_ranking": True,
        },
        "scan_top1_audit_baseline": {
            "selection": "highest_tic_then_spectrum_id_then_lowest_result_index",
            "result_index": (
                int(baseline_scan["row"]["result_index"]) if baseline_scan else None
            ),
            "spectrum_id": (
                str(baseline_scan["row"].get("spectrum_id") or "")
                if baseline_scan else None
            ),
            "identity_key": baseline_identity,
        },
    }


def _original_base_entity_audit(entity, policy):
    """Replay original support/Top20, with display-only expansion for all audits.

    The original finalizer uses shortlist_size only when constructing public
    copies, after all scores, chemistry, ranks and decision calculations.
    Increasing that display limit does not change original candidate admission.
    """
    original = dict(entity)
    original["scans"] = []
    identities = set()
    for scan in entity.get("scans") or []:
        copied = dict(scan)
        candidates = scan.get("candidates") or []
        if candidates:
            key = "_original_base_projection_for_legacy_audit"
            if key not in candidates[0]:
                raise ValueError("Exact-tie legacy audit requires original base member projection")
            copied["candidates"] = candidates[0][key]
        else:
            copied["candidates"] = []
        identities.update(c["identity_key"] for c in copied["candidates"])
        original["scans"].append(copied)
    audit_policy = copy.deepcopy(policy)
    audit_policy.setdefault("candidate_pool", {})["shortlist_size"] = max(1, len(identities))
    return _original_base_finalize_entity(original, audit_policy)


def _restore_original_legacy_fields(candidates, original_base_audit):
    originals = {c["identity_key"]: c for c in original_base_audit["candidates"]}
    for candidate in candidates:
        prior = originals.get(candidate["identity_key"])
        if prior is None:
            # A numeric rounding boundary can change which identity fits into
            # the existing cap. Never invent an original audit for that identity.
            candidate["legacy_v14_entity_score"] = None
            candidate["legacy_v14_entity_rank"] = None
            candidate["legacy_v14_audit_status"] = "identity_not_in_original_base_entity_pool"
            for evidence in candidate["collision_energy_evidence"]:
                evidence["legacy_v14_support"] = None
                evidence["legacy_v14_scan_top1_fraction"] = None
            continue
        for key in ("legacy_v14_entity_score", "legacy_v14_entity_rank"):
            candidate[key] = prior[key]
        prior_bins = {e["collision_energy_bin"]: e for e in prior["collision_energy_evidence"]}
        for evidence in candidate["collision_energy_evidence"]:
            old = prior_bins.get(evidence["collision_energy_bin"])
            for key in ("legacy_v14_support", "legacy_v14_scan_top1_fraction"):
                evidence[key] = old[key] if old is not None else None
        candidate["legacy_v14_audit_status"] = "exact_original_base_aggregation"
