"""Utility package with lazy compatibility exports.

Importing :mod:`utils` used to import matchms, numba and matplotlib
immediately.  Lightweight tools such as vendor-raw inspection and storage
preflight should not need that scientific runtime (and its writable caches),
so the historical convenience symbols are resolved only when requested.
"""

__all__ = ["clean_peaks", "parse_mgf_file", "parse_mzml_file"]


def __getattr__(name):
    if name == "clean_peaks":
        from .spectrum import clean_peaks

        return clean_peaks
    if name == "parse_mgf_file":
        from .mgf_io import parse_mgf_file

        return parse_mgf_file
    if name == "parse_mzml_file":
        from .mzml_io import parse_mzml_file

        return parse_mzml_file
    raise AttributeError(name)
