"""The single scientific result object shared by future API and report layers."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .decision import ClaimDecision
from .evidence import EvidenceNode, ModuleState


IDENTIFICATION_SCHEMA_VERSION = "1.0.0"


class ObservationReference(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    task_id: str
    spectrum_id: str
    precursor_event_id: str | None = None
    chromatographic_feature_id: str | None = None
    ion_identity_group_id: str | None = None
    neutral_analyte_id: str | None = None
    injection_group_id: str | None = None
    polarity: str | None = None
    precursor_mz: float | None = None
    retention_time_minutes: float | None = None


class CandidateIdentity(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    identity_key: str
    compound_name: str
    inchikey: str | None = None
    connectivity_key: str | None = None
    formula: str | None = None
    smiles: str | None = None
    adduct: str | None = None
    source_tier: str | None = None
    library_spectrum_id: str | None = None
    library_accession: str | None = None


class ScoreBundle(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    selected_rank: int = Field(default=1, ge=1)
    cosine: float | None = None
    entropy: float | None = None
    reverse_match: float | None = None
    mass_ppm: float | None = None
    final_score: float | None = None
    rank_score: float | None = None
    raw_rank_score: float | None = None
    candidate_margin: float | None = None
    legacy_confidence: str | None = None
    ranking_source: str = "legacy_production_search"
    ranking_changed: bool = False


class IdentificationRecord(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: str = IDENTIFICATION_SCHEMA_VERSION
    record_id: str
    record_status: str = "matched_library_candidate"
    observation: ObservationReference
    candidate: CandidateIdentity
    scores: ScoreBundle
    evidence: tuple[EvidenceNode, ...]
    decision: ClaimDecision
    module_states: dict[str, ModuleState] = Field(default_factory=dict)
    provenance: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def enforce_record_invariants(self) -> "IdentificationRecord":
        node_ids = [node.node_id for node in self.evidence]
        if len(node_ids) != len(set(node_ids)):
            raise ValueError("Evidence node_id values must be unique within a record")
        if self.decision.ranking_changed != self.scores.ranking_changed:
            raise ValueError("Decision and score ranking_changed fields must agree")
        return self
