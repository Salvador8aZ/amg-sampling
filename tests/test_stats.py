import math
import random

import pytest

from tests.stats import chi2_sf, chi_square_gof


@pytest.mark.parametrize(
    "x, df, p",
    [
        (3.841458820694124, 1, 0.05),
        (6.634896601021214, 1, 0.01),
        (18.307038053275146, 10, 0.05),
        (124.34211340400407, 100, 0.05),
        (0.5, 3, 0.9188914116),
    ],
)
def test_chi2_sf_known_quantiles(x, df, p):
    assert chi2_sf(x, df) == pytest.approx(p, rel=1e-8)


@pytest.mark.parametrize("x", [0.1, 1.0, 5.0, 40.0])
def test_chi2_sf_two_degrees_of_freedom_is_exponential(x):
    assert chi2_sf(x, 2) == pytest.approx(math.exp(-x / 2), rel=1e-10)


def test_gof_accepts_a_fair_die_and_rejects_a_loaded_one():
    rng = random.Random(1)
    probs = {k: 1 / 6 for k in range(6)}
    fair = {k: 0 for k in range(6)}
    for _ in range(6000):
        fair[rng.randrange(6)] += 1
    assert chi_square_gof(fair, probs)[2] > 1e-3
    loaded = {k: 0 for k in range(6)}
    for _ in range(6000):
        loaded[min(rng.randrange(6), rng.randrange(6))] += 1
    assert chi_square_gof(loaded, probs)[2] < 1e-10


def test_gof_rejects_categories_with_zero_probability():
    with pytest.raises(AssertionError):
        chi_square_gof({"a": 3, "b": 1}, {"a": 1.0})
