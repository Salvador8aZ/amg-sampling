"""Small, dependency-free statistical tests used to validate samplers.

Pearson's chi-square goodness-of-fit test with the upper tail of the
chi-square distribution computed from the regularised incomplete gamma
function (series / continued fraction, as in Numerical Recipes §6.2).
"""

from __future__ import annotations

import math
from collections.abc import Mapping

# Significance level used by sampler tests. Tests are seeded, so they are
# deterministic; under the null hypothesis each test would fail with
# probability ALPHA for a random seed.
ALPHA = 1e-3


def _gamma_p_series(a: float, x: float) -> float:
    term = total = 1.0 / a
    ap = a
    for _ in range(10_000):
        ap += 1.0
        term *= x / ap
        total += term
        if abs(term) < abs(total) * 1e-15:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gamma_q_continued_fraction(a: float, x: float) -> float:
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 10_000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def chi2_sf(statistic: float, df: int) -> float:
    """P(χ²_df ≥ statistic)."""
    if df <= 0:
        raise ValueError("df must be positive.")
    if statistic <= 0:
        return 1.0
    a, x = df / 2.0, statistic / 2.0
    if x < a + 1.0:
        return 1.0 - _gamma_p_series(a, x)
    return _gamma_q_continued_fraction(a, x)


def chi_square_gof(
    observed: Mapping, probabilities: Mapping, *, min_expected: float = 5.0
) -> tuple[float, int, float]:
    """Pearson goodness-of-fit of ``observed`` counts against ``probabilities``.

    Categories with expected count below ``min_expected`` are pooled into one
    bin (the usual validity condition for the χ² approximation). Every observed
    category must have positive probability. Returns (statistic, df, p-value).
    """
    unknown = set(observed) - set(probabilities)
    if unknown:
        raise AssertionError(
            "observed categories with zero probability: "
            f"{sorted(map(str, unknown))[:5]}"
        )
    total = sum(observed.values())
    bins: list[tuple[float, float]] = []
    pooled_obs = pooled_exp = 0.0
    for key, p in probabilities.items():
        expected = total * float(p)
        count = observed.get(key, 0)
        if expected < min_expected:
            pooled_obs += count
            pooled_exp += expected
        else:
            bins.append((count, expected))
    if pooled_exp > 0:
        bins.append((pooled_obs, pooled_exp))
    if len(bins) < 2:
        raise ValueError("Need at least two bins; increase the sample size.")
    statistic = sum((o - e) ** 2 / e for o, e in bins)
    df = len(bins) - 1
    return statistic, df, chi2_sf(statistic, df)


def binomial_z(successes: int, trials: int, p: float) -> float:
    """Standardised deviation of a binomial count from its mean ``trials * p``."""
    return (successes - trials * p) / math.sqrt(trials * p * (1 - p))
