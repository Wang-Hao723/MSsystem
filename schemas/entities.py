"""Candidate-independent acquisition entities and their measured relationships."""

from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator


ENTITY_SCHEMA_VERSION = "1.0.0"


class EntityLevel(str, Enum):
    MS2_SCAN = "ms2_scan"
    PRECURSOR_EVENT = "precursor_event"
    CHROMATOGRAPHIC_FEATURE = "chromatographic_feature"
    ION_IDENTITY_GROUP = "ion_identity_group"
    NEUTRAL_ANALYTE = "neutral_analyte"


class EntityNode(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    entity_id: str
    level: EntityLevel
    task_id: str
    polarity: str | None = None
    representative_mz: float | None = None
    retention_time_min: float | None = None
    rt_start_min: float | None = None
    rt_end_min: float | None = None
    member_count: int = Field(default=1, ge=1)
    construction_source: str
    attributes: dict[str, Any] = Field(default_factory=dict)


class EntityEdge(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    source_id: str
    target_id: str
    relation: str
    construction_basis: str
    candidate_independent: bool = True


class EntityGraph(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: str = ENTITY_SCHEMA_VERSION
    task_id: str
    nodes: tuple[EntityNode, ...]
    edges: tuple[EntityEdge, ...]
    summary: dict[str, Any]
    construction_inputs: tuple[str, ...]
    candidate_independent: bool = True
    graph_sha256: str

    @model_validator(mode="after")
    def validate_graph(self) -> "EntityGraph":
        node_ids = [node.entity_id for node in self.nodes]
        if len(node_ids) != len(set(node_ids)):
            raise ValueError("Entity IDs must be unique")
        known = set(node_ids)
        dangling = [
            (edge.source_id, edge.target_id)
            for edge in self.edges
            if edge.source_id not in known or edge.target_id not in known
        ]
        if dangling:
            raise ValueError(f"Entity graph contains dangling edges: {dangling[:3]}")
        if not self.candidate_independent or any(
            not edge.candidate_independent for edge in self.edges
        ):
            raise ValueError("Entity construction must remain candidate-independent")
        return self
