MassCube 1.0.22 reference implementation for formula regression tests.

`masscube_1_0_22_feature_evaluation.py` is the unchanged
`masscube/feature_evaluation.py` extracted from the PyPI 1.0.22 wheel.
The original MIT license is in `MASSCUBE_LICENSE.txt` in this directory.
The per-file whitespace attribute retains the unchanged vendor file's five
original trailing-space lines; authored files use the repository defaults.

Paper: https://doi.org/10.1038/s41467-025-60640-5

Wheel SHA256:
`9276d3cff2a73f67b634c8ccedf797e206d86292998e9614eb31121a5d5bf9d8`

These functions are used as independent expected values on eligible traces.
The implementation deliberately returns explicit invalid status for degenerate
traces instead of adopting the author's short/constant-trace ISF fallback.
