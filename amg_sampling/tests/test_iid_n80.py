"""Validation of the IID sampler at n = 80 against exact summary distributions."""

import math
import random
from collections import Counter
from fractions import Fraction

import pytest

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import cycle_structure
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.formulas import num_all_states
from amg_sampling.exact.summaries import (
    cycle_count_distribution,
    cycle_length_multiplicity_distribution,
    largest_cycle_distribution,
)
from amg_sampling.samplers.rejection import iter_samples
from amg_sampling.tests.stats import ALPHA, chi_square_gof

N_DSBS = 80
SAMPLES = 8_000
Z_LIMIT = 4.0


def normalise(counts):
    total = sum(counts.values())
    return {k: Fraction(v, total) for k, v in counts.items()}


def sampled_structures(breaks, space, seed):
    theta = InitialConfiguration(breaks)
    return [
        cycle_structure(theta, r)
        for r in iter_samples(theta, space, random.Random(seed), SAMPLES)
    ]


@pytest.fixture(scope="module")
def single_chromosome_proper():
    # Θ(1,(80)): PROPER = DERANGED, so the DERANGED formulas are exact here.
    return sampled_structures((N_DSBS,), StateSpace.PROPER, 801)


@pytest.fixture(scope="module")
def all_states():
    return sampled_structures((N_DSBS,), StateSpace.ALL, 802)


@pytest.mark.parametrize(
    "name, exact, summary",
    [
        (
            "cycles",
            lambda: cycle_count_distribution(N_DSBS, StateSpace.DERANGED),
            lambda c: c.num_cycles,
        ),
        (
            "largest",
            lambda: largest_cycle_distribution(N_DSBS, StateSpace.DERANGED),
            lambda c: c.parts[0],
        ),
        (
            "C2",
            lambda: cycle_length_multiplicity_distribution(
                N_DSBS, 2, StateSpace.DERANGED
            ),
            lambda c: c.count(2),
        ),
    ],
)
def test_single_chromosome_summaries(single_chromosome_proper, name, exact, summary):
    observed = Counter(summary(c) for c in single_chromosome_proper)
    _, _, p = chi_square_gof(observed, normalise(exact()))
    assert p > ALPHA, name


@pytest.mark.parametrize(
    "name, exact, summary",
    [
        (
            "cycles",
            lambda: cycle_count_distribution(N_DSBS, StateSpace.ALL),
            lambda c: c.num_cycles,
        ),
        (
            "largest",
            lambda: largest_cycle_distribution(N_DSBS, StateSpace.ALL),
            lambda c: c.parts[0],
        ),
        (
            "C1",
            lambda: cycle_length_multiplicity_distribution(N_DSBS, 1, StateSpace.ALL),
            lambda c: c.count(1),
        ),
        (
            "C2",
            lambda: cycle_length_multiplicity_distribution(N_DSBS, 2, StateSpace.ALL),
            lambda c: c.count(2),
        ),
    ],
)
def test_all_summaries(all_states, name, exact, summary):
    observed = Counter(summary(c) for c in all_states)
    _, _, p = chi_square_gof(observed, normalise(exact()))
    assert p > ALPHA, name


def ewens_half_moments(n):
    """Derived: number of cycles under uniform ALL is Σ Bernoulli(1/(2i+1)), i < n."""
    probs = [Fraction(1, 2 * i + 1) for i in range(n)]
    mean = sum(probs)
    var = sum(p * (1 - p) for p in probs)
    return mean, var


def test_ewens_moments_match_exact_distribution():
    dist = cycle_count_distribution(N_DSBS, StateSpace.ALL)
    total = num_all_states(N_DSBS)
    mean = Fraction(sum(c * m for c, m in dist.items()), total)
    var = Fraction(sum(c * c * m for c, m in dist.items()), total) - mean**2
    assert (mean, var) == ewens_half_moments(N_DSBS)


def test_sampled_cycle_count_mean_and_variance(all_states):
    mean, var = (float(x) for x in ewens_half_moments(N_DSBS))
    ks = [c.num_cycles for c in all_states]
    m = len(ks)
    sample_mean = sum(ks) / m
    assert abs(sample_mean - mean) / math.sqrt(var / m) < Z_LIMIT
    # Variance: the sample variance has standard error sqrt((mu4 - var^2) / m),
    # with the exact fourth central moment taken from the exact distribution.
    dist = cycle_count_distribution(N_DSBS, StateSpace.ALL)
    total = num_all_states(N_DSBS)
    mu4 = float(
        sum(Fraction(m_ * (c - Fraction(mean)) ** 4, total) for c, m_ in dist.items())
    )
    sample_var = sum((k - sample_mean) ** 2 for k in ks) / (m - 1)
    assert abs(sample_var - var) / math.sqrt((mu4 - var**2) / m) < Z_LIMIT
