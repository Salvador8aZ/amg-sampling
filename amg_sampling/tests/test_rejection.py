import math
import random
from collections import Counter
from fractions import Fraction

import pytest

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import cycle_structure
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.enumerate import (
    count_states,
    exact_cycle_distribution,
    iter_states,
)
from amg_sampling.exact.formulas import (
    cycle_distribution_all,
    cycle_distribution_deranged,
    cycle_distribution_proper,
    num_all_states,
    num_proper_states,
)
from amg_sampling.samplers.rejection import (
    SamplingStats,
    iter_samples,
    sample_state,
    theoretical_acceptance,
)
from amg_sampling.tests.reference import compositions
from amg_sampling.tests.stats import ALPHA, chi_square_gof

ALL, DERANGED, PROPER = StateSpace.ALL, StateSpace.DERANGED, StateSpace.PROPER

# Two-sided normal tail probability for |z| > 4 is 6.3e-5.
Z_LIMIT = 4.0


def normalise(counts):
    total = sum(counts.values())
    return {k: Fraction(v, total) for k, v in counts.items()}


# -- level 1: every state is equally likely ------------------------------------


@pytest.mark.parametrize(
    "breaks, space, seed",
    [
        ((2, 2), ALL, 1),
        ((4,), DERANGED, 2),
        ((3,), DERANGED, 3),
        ((2, 2), PROPER, 4),
        ((2, 1, 1), PROPER, 5),
        ((1, 1, 1, 1), PROPER, 6),
        ((3, 1), PROPER, 7),
        ((2, 1), PROPER, 8),
    ],
)
def test_state_level_uniformity(breaks, space, seed):
    theta = InitialConfiguration(breaks)
    states = set(iter_states(theta, space))
    per_state = 150
    samples = list(
        iter_samples(theta, space, random.Random(seed), per_state * len(states))
    )
    observed = Counter(samples)
    assert set(observed) == states  # every sample is a member; every member appears
    _, _, p = chi_square_gof(observed, {s: Fraction(1, len(states)) for s in states})
    assert p > ALPHA


# -- level 2: cycle-structure distribution matches the exact one ----------------


def exact_distribution(theta, space):
    if space is ALL:
        return cycle_distribution_all(theta.num_dsbs)
    if space is DERANGED:
        return cycle_distribution_deranged(theta.num_dsbs)
    return cycle_distribution_proper(theta)


@pytest.mark.parametrize(
    "breaks, space, seed",
    [
        ((6,), ALL, 21),
        ((3, 2, 2), ALL, 22),
        ((7,), DERANGED, 23),
        ((2, 2, 1, 1), DERANGED, 24),
        ((6,), PROPER, 25),
        ((3, 3), PROPER, 26),
        ((2, 2, 2), PROPER, 27),
        ((2, 1, 1, 1, 1), PROPER, 28),
        ((3, 2, 2), PROPER, 29),
        ((4, 4, 4), PROPER, 30),
        ((6, 3, 2, 1), PROPER, 31),
        ((3, 3, 3, 3), PROPER, 32),
        ((2,) * 6, PROPER, 33),
    ],
)
def test_cycle_structure_distribution(breaks, space, seed):
    theta = InitialConfiguration(breaks)
    exact = normalise(exact_distribution(theta, space))
    observed = Counter(
        cycle_structure(theta, r)
        for r in iter_samples(theta, space, random.Random(seed), 20_000)
    )
    _, _, p = chi_square_gof(observed, exact)
    assert p > ALPHA


def test_exact_references_agree_for_the_small_cases():
    # The formulas used above are themselves checked against enumeration.
    theta = InitialConfiguration((2, 2, 1, 1))
    for space in StateSpace:
        formula = {c: m for c, m in exact_distribution(theta, space).items() if m}
        assert formula == exact_cycle_distribution(theta, space)


# -- level 3: acceptance rate ----------------------------------------------------


@pytest.mark.parametrize(
    "breaks, space, seed",
    [
        ((4,), DERANGED, 41),
        ((10,), DERANGED, 42),
        ((2, 2), PROPER, 43),
        ((1,) * 6, PROPER, 44),
        ((3, 3), PROPER, 45),
        ((4, 4, 4), PROPER, 46),
        ((1,) * 12, PROPER, 47),
        ((20, 20), PROPER, 48),
    ],
)
def test_acceptance_rate(breaks, space, seed):
    # The sampler stops after a fixed number `a` of acceptances, so the number
    # of rejections F is negative binomial: E[F] = a(1-p)/p, Var[F] = a(1-p)/p^2.
    theta = InitialConfiguration(breaks)
    p = float(theoretical_acceptance(theta, space))
    stats = SamplingStats()
    accepted = 5_000
    for _ in iter_samples(theta, space, random.Random(seed), accepted, stats=stats):
        pass
    assert stats.accepted == accepted
    mean = accepted * (1 - p) / p
    sd = math.sqrt(accepted * (1 - p)) / p
    assert abs(stats.rejected - mean) / sd < Z_LIMIT


def test_all_never_rejects():
    stats = SamplingStats()
    list(
        iter_samples(
            InitialConfiguration((3, 2)), ALL, random.Random(0), 200, stats=stats
        )
    )
    assert (stats.proposals, stats.accepted, stats.rejected) == (200, 200, 0)
    assert stats.acceptance_rate == 1.0


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_theoretical_acceptance_matches_enumeration(n):
    for breaks in compositions(n):
        theta = InitialConfiguration(breaks)
        for space in StateSpace:
            assert theoretical_acceptance(theta, space) == Fraction(
                count_states(theta, space), num_all_states(n)
            )


@pytest.mark.parametrize("n", range(2, 11))
def test_all_ones_acceptance_matches_recursion(n):
    theta = InitialConfiguration([1] * n)
    assert theoretical_acceptance(theta, PROPER, max_chromosomes=n) == Fraction(
        num_proper_states(theta), num_all_states(n)
    )


def test_theoretical_acceptance_refuses_infeasible_recursion():
    theta = InitialConfiguration([2] * 20)
    with pytest.raises(ValueError, match="3\\^k"):
        theoretical_acceptance(theta, PROPER)
    # Closed forms are still available for this Θ.
    assert theoretical_acceptance(theta, DERANGED) > 0
    assert theoretical_acceptance(InitialConfiguration([1] * 80), PROPER) > 0


def test_stats_before_any_proposal():
    assert SamplingStats().acceptance_rate is None


# -- level 4: independence (sanity check; the argument is in sampling.md §3) -----


@pytest.mark.parametrize(
    "breaks, space, seed", [((2,), PROPER, 51), ((3,), DERANGED, 52)]
)
def test_non_overlapping_pairs_are_uniform_on_the_product(breaks, space, seed):
    theta = InitialConfiguration(breaks)
    states = list(iter_states(theta, space))
    per_cell = 100
    draws = list(
        iter_samples(theta, space, random.Random(seed), 2 * per_cell * len(states) ** 2)
    )
    pairs = Counter(zip(draws[0::2], draws[1::2]))
    cells = {(a, b): Fraction(1, len(states) ** 2) for a in states for b in states}
    _, _, p = chi_square_gof(pairs, cells)
    assert p > ALPHA


# -- reproducibility, errors -------------------------------------------------------


def test_same_seed_same_samples_and_stats():
    theta = InitialConfiguration((3, 2, 2))
    s1, s2 = SamplingStats(), SamplingStats()
    a = list(iter_samples(theta, PROPER, random.Random(99), 100, stats=s1))
    b = list(iter_samples(theta, PROPER, random.Random(99), 100, stats=s2))
    assert a == b and s1 == s2
    rng1, rng2 = random.Random(7), random.Random(7)
    assert [sample_state(theta, PROPER, rng1) for _ in range(20)] == [
        sample_state(theta, PROPER, rng2) for _ in range(20)
    ]


def test_global_random_state_is_untouched():
    before = random.getstate()
    list(iter_samples(InitialConfiguration((2, 2)), PROPER, random.Random(3), 50))
    assert random.getstate() == before


@pytest.mark.parametrize("space", [DERANGED, PROPER])
def test_empty_state_space_is_refused(space):
    with pytest.raises(ValueError, match="empty"):
        sample_state(InitialConfiguration([1]), space, random.Random(0))
    with pytest.raises(ValueError, match="empty"):
        next(iter_samples(InitialConfiguration([1]), space, random.Random(0)))


def test_max_proposals_is_enforced():
    theta = InitialConfiguration([1] * 10)  # acceptance about 0.28
    raised = 0
    for seed in range(50):
        try:
            sample_state(theta, PROPER, random.Random(seed), max_proposals=1)
        except RuntimeError:
            raised += 1
    assert 0 < raised < 50
