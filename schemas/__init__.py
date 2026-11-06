"""Versioned scientific result contracts used across MSsystem V10."""

from .decision import CalibrationVector, ClaimDecision, ClaimScope
from .evidence import EvidenceNode, EvidenceProvenance, ModuleState
from .identification import (
    CandidateIdentity,
    IdentificationRecord,
    ObservationReference,
    ScoreBundle,
)
from .legacy import legacy_row_to_identification

__all__ = [
    "CalibrationVector",
    "CandidateIdentity",
    "ClaimDecision",
    "ClaimScope",
    "EvidenceNode",
    "EvidenceProvenance",
    "IdentificationRecord",
    "ModuleState",
    "ObservationReference",
    "ScoreBundle",
    "legacy_row_to_identification",
]
