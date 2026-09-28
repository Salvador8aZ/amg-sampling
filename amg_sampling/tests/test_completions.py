"""Fixed (observed) rejoins: partial matchings, exact completions, uniform completion sampling."""

import random
from collections import Counter
from fractions import Fraction
from itertools import combinations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import InvalidMatchingError, PartialMatching, RejoinMatching
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.enumerate import count_states, iter_states
from amg_sampling.exact.formulas import num_completions
from amg_sampling.samplers.rejection import SamplingStats, iter_samples
from amg_sampling.samplers.uniform import sample_completion
from amg_sampling.tests.reference import compositions
from amg_sampling.tests.stats import ALPHA, chi_square_gof


# -- PartialMatching ------------------------------------------------------------


def test_partial_matching_basics():
    f = PartialMatching(8, [(5, 2), (0, 7)])
    assert f.pairs() == ((0, 7), (2, 5))
    assert f.free_ends == (1, 3, 4, 6)
    assert f.partner(2) == 5 and f.partner(1) is None
    assert len(f) == 2
    assert f == PartialMatching(8, [(0, 7), (2, 5)])
    assert f.is_completed_by(RejoinMatching.from_pairs(8, [(0, 7), (2, 5), (1, 3), (4, 6)]))
    assert not f.is_completed_by(RejoinMatching.from_pairs(8, [(0, 1), (2, 5), (7, 3), (4, 6)]))


def test_partial_matching_accepts_a_one_shot_iterator():
    assert PartialMatching(4, zip([0], [1])).pairs() == ((0, 1),)


@pytest.mark.parametrize(
    "num_ends, pairs",
    [(3, []), (4, [(0, 0)]), (4, [(0, 1), (1, 2)]), (4, [(0, 4)]), (-2, [])],
)
def test_partial_matching_rejects_invalid_input(num_ends, pairs):
    with pytest.raises(InvalidMatchingError):
        PartialMatching(num_ends, pairs)


@pytest.mark.parametrize("f", [0, 2, 4, 6, 10])
def test_num_completions(f):
    assert num_completions(f) == [1, 1, 3, 15, 105, 945][f // 2]


# -- exact completions -------------------------------------------------------------


def fixed_edge_sets(num_ends, max_edges=2):
    """Every set of up to ``max_edges`` disjoint fixed edges."""
    edges = list(combinations(range(num_ends), 2))
    yield PartialMatching(num_ends)
    for k in range(1, max_edges + 1):
        for chosen in combinations(edges, k):
            ends = [v for e in chosen for v in e]
            if len(set(ends)) == len(ends):
                yield PartialMatching(num_ends, chosen)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_completion_enumeration_equals_filtered_enumeration(n):
    # Listing completions directly must give exactly the states that contain
    # the fixed edges, for every layout, state space and small fixed set.
    for breaks in compositions(n):
        theta = InitialConfiguration(breaks)
        for space in StateSpace:
            everything = list(iter_states(theta, space))
            for fixed in fixed_edge_sets(theta.num_ends):
                listed = list(iter_states(theta, space, fixed))
                assert len(listed) == len(set(listed))
                assert set(listed) == {r for r in everything if fixed.is_completed_by(r)}


@given(st.data())
def test_all_completions_count_is_double_factorial(data):
    n = data.draw(st.integers(1, 5))
    theta = InitialConfiguration([n])
    ends = data.draw(st.permutations(range(2 * n)))
    k = data.draw(st.integers(0, n))
    fixed = PartialMatching(2 * n, zip(ends[: 2 * k : 2], ends[1 : 2 * k : 2]))
    assert count_states(theta, StateSpace.ALL, fixed) == num_completions(2 * (n - k))


def test_fixed_dsb_edge_leaves_no_deranged_completion():
    theta = InitialConfiguration([3])
    fixed = PartialMatching(6, [(0, 1)])
    assert count_states(theta, StateSpace.ALL, fixed) == 3
    assert count_states(theta, StateSpace.DERANGED, fixed) == 0


def test_mismatched_sizes_raise():
    with pytest.raises(ValueError):
        list(iter_states(InitialConfiguration([2]), StateSpace.ALL, PartialMatching(6)))


# -- uniform completion sampling -----------------------------------------------------


def test_sample_completion_contains_fixed_edges_and_is_reproducible():
    fixed = PartialMatching(12, [(0, 5), (7, 10)])
    rng1, rng2 = random.Random(3), random.Random(3)
    a = [sample_completion(fixed, rng1) for _ in range(30)]
    b = [sample_completion(fixed, rng2) for _ in range(30)]
    assert a == b
    assert all(fixed.is_completed_by(r) for r in a)


def test_sample_completion_with_no_free_ends_returns_the_fixed_matching():
    fixed = PartialMatching(4, [(0, 2), (1, 3)])
    assert sample_completion(fixed, random.Random(0)).pairs() == ((0, 2), (1, 3))


def test_completion_sampling_is_uniform_over_completions():
    theta = InitialConfiguration((3, 1))
    fixed = PartialMatching(8, [(1, 6)])  # 6 free ends: 15 completions
    completions = set(iter_states(theta, StateSpace.ALL, fixed))
    assert len(completions) == 15
    rng = random.Random(21)
    observed = Counter(sample_completion(fixed, rng) for _ in range(150 * 15))
    assert set(observed) == completions
    assert chi_square_gof(observed, {c: Fraction(1, 15) for c in completions})[2] > ALPHA


@pytest.mark.parametrize(
    "breaks, pairs, seed",
    [((3, 2), [(1, 6)], 31), ((2, 2, 1), [(0, 9), (3, 5)], 32), ((4, 1), [(8, 2)], 33)],
)
def test_conditioned_rejection_is_uniform_on_proper_completions(breaks, pairs, seed):
    theta = InitialConfiguration(breaks)
    fixed = PartialMatching(theta.num_ends, pairs)
    target = set(iter_states(theta, StateSpace.PROPER, fixed))
    stats = SamplingStats()
    draws = list(
        iter_samples(theta, StateSpace.PROPER, random.Random(seed), 150 * len(target), stats=stats, fixed=fixed)
    )
    observed = Counter(draws)
    assert set(observed) == target
    assert chi_square_gof(observed, {s: Fraction(1, len(target)) for s in target})[2] > ALPHA
    # Proposals come only from completions, so acceptance is |target| / #completions.
    assert stats.proposals >= stats.accepted


def test_conditioned_rejection_respects_max_proposals():
    theta = InitialConfiguration([3])
    fixed = PartialMatching(6, [(0, 1)])  # no deranged completion exists
    with pytest.raises(RuntimeError):
        next(iter_samples(theta, StateSpace.DERANGED, random.Random(0), fixed=fixed, max_proposals_per_sample=50))
