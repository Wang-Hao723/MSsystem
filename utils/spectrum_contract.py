"""Shared detector and spectrum-type contracts for mzML consumers."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


AUXILIARY_DETECTOR_SPECTRUM_TERMS = (
    "electromagnetic radiation spectrum",
    "absorption spectrum",
    "emission spectrum",
    "wavelength array",
)


def is_auxiliary_detector_spectrum(spectrum: Mapping[str, Any]) -> bool:
    """Return true for DAD/UV/optical spectra even if MS level was synthesized.

    Some OpenMS round trips assign MS level 1 to detector traces that originally
    had no mass-spectrometry level.  Optical CV terms therefore take precedence
    over ``ms level`` in every consumer, including representation audits and
    TIC/BPC extraction.
    """

    return any(term in spectrum for term in AUXILIARY_DETECTOR_SPECTRUM_TERMS)
