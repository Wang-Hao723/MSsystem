import random

import numpy as np
from matchms import Spectrum
from matchms.similarity import ModifiedCosineGreedy, NeutralLossesCosine

from matching.scores import (
    batch_modified_cosine_scores,
    batch_neutral_loss_cosine_scores,
    cosine_score,
)


def _matchms_reference(query, library, query_pmz, library_pmz, tolerance=0.1):
    q_intensity = np.asarray([peak[1] for peak in query], dtype=float)
    l_intensity = np.asarray([peak[1] for peak in library], dtype=float)
    q_intensity /= q_intensity.max()
    l_intensity /= l_intensity.max()
    q = Spectrum(
        mz=np.asarray([peak[0] for peak in query], dtype=float),
        intensities=q_intensity,
        metadata={"precursor_mz": query_pmz},
    )
    lib = Spectrum(
        mz=np.asarray([peak[0] for peak in library], dtype=float),
        intensities=l_intensity,
        metadata={"precursor_mz": library_pmz},
    )
    return float(ModifiedCosineGreedy(tolerance=tolerance).pair(q, lib)["score"])


def test_array_cosine_is_numerically_equivalent_to_matchms_wrapper():
    rng = random.Random(20260715)
    for case in range(100):
        query_pmz = rng.uniform(150, 1000)
        library_pmz = query_pmz + rng.uniform(-0.5, 0.5)
        query = sorted([[rng.uniform(20, query_pmz - 1), rng.random()] for _ in range(15)])
        library = sorted([[rng.uniform(20, library_pmz - 1), rng.random()] for _ in range(15)])
        # Mix exact and precursor-shifted matches into the random background.
        if case % 2 == 0:
            library[0][0] = query[0][0]
        if case % 3 == 0:
            library[1][0] = query[1][0] - (query_pmz - library_pmz)
        library.sort()

        expected = _matchms_reference(query, library, query_pmz, library_pmz)
        observed = cosine_score(
            query, library, query_pmz=query_pmz, lib_pmz=library_pmz,
        )
        assert observed == pytest.approx(expected, abs=1e-12)


def test_batch_cosine_is_numerically_equivalent_to_matchms_wrapper():
    rng = random.Random(20260716)
    query_pmz = 500.0
    query = sorted([[rng.uniform(20, 490), rng.random()] for _ in range(15)])
    libraries, precursor_mz, expected = [], [], []
    for case in range(100):
        pmz = query_pmz + rng.uniform(-0.5, 0.5)
        peaks = sorted([[rng.uniform(20, max(21, pmz - 1)), rng.random()] for _ in range(15)])
        if case % 2 == 0:
            peaks[0][0] = query[0][0]
        if case % 3 == 0:
            peaks[1][0] = query[1][0] - (query_pmz - pmz)
        peaks.sort()
        libraries.append(peaks)
        precursor_mz.append(pmz)
        expected.append(_matchms_reference(query, peaks, query_pmz, pmz))

    observed = batch_modified_cosine_scores(query, libraries, query_pmz, precursor_mz)
    assert observed == pytest.approx(expected, abs=1e-12)


def test_batch_neutral_loss_cosine_is_equivalent_to_matchms_wrapper():
    query_pmz = 300.0
    query = [[100.0, 1.0], [150.0, 0.8], [220.0, 0.4]]
    libraries = [
        [[90.0, 1.0], [140.0, 0.8], [210.0, 0.4]],
        [[90.5, 1.0], [141.0, 0.8], [205.0, 0.4]],
    ]
    library_pmz = [290.0, 290.0]
    query_spectrum = Spectrum(
        mz=np.asarray([peak[0] for peak in query]),
        intensities=np.asarray([peak[1] for peak in query]),
        metadata={"precursor_mz": query_pmz},
    )
    similarity = NeutralLossesCosine(tolerance=0.02)
    expected = []
    for peaks, precursor in zip(libraries, library_pmz):
        reference = Spectrum(
            mz=np.asarray([peak[0] for peak in peaks]),
            intensities=np.asarray([peak[1] for peak in peaks]),
            metadata={"precursor_mz": precursor},
        )
        expected.append(float(similarity.pair(reference, query_spectrum)["score"]))
    observed = batch_neutral_loss_cosine_scores(
        query,
        libraries,
        query_pmz,
        library_pmz,
        mz_tol=0.02,
    )
    assert observed == pytest.approx(expected, abs=1e-12)


# Imported last so the test helper remains visually focused on the numerical path.
import pytest  # noqa: E402
