import random
from collections import Counter
from itertools import permutations
from math import factorial

import pytest

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.enumerate import iter_states
from amg_sampling.samplers.uniform import matching_from_ordering, sample_matching
from amg_sampling.tests.stats import ALPHA, chi_square_gof

# -- the counting argument, checked exhaustively ------------------------------


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_every_matching_arises_from_exactly_2n_nfactorial_orderings(n):
    # Deterministic check of the uniformity proof: pairing consecutive
    # positions maps the (2n)! orderings onto the (2n-1)!! matchings,
    # each matching having exactly 2^n n! preimages.
    counts = Counter(matching_from_ordering(p) for p in permutations(range(2 * n)))
    all_states = set(iter_states(InitialConfiguration([n]), StateSpace.ALL))
    assert set(counts) == all_states
    assert set(counts.values()) == {2**n * factorial(n)}


def test_matching_from_ordering_rejects_bad_input():
    with pytest.raises(ValueError):
        matching_from_ordering([0, 0])
    with pytest.raises(ValueError):
        matching_from_ordering([0, 1, 2])


# -- basic behaviour -----------------------------------------------------------


@pytest.mark.parametrize("n", [0, 1, 5, 80])
def test_sample_is_a_perfect_matching(n):
    r = sample_matching(n, random.Random(0))
    assert r.num_ends == 2 * n
    assert all(r[r[v]] == v and r[v] != v for v in range(2 * n))


def test_invalid_arguments():
    with pytest.raises(TypeError):
        sample_matching(3, 42)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        sample_matching(-1, random.Random(0))


# -- reproducibility -----------------------------------------------------------


def test_same_seed_gives_same_sequence():
    a = [sample_matching(10, random.Random(123)) for _ in range(3)]
    rng1, rng2 = random.Random(123), random.Random(123)
    seq1 = [sample_matching(10, rng1) for _ in range(50)]
    seq2 = [sample_matching(10, rng2) for _ in range(50)]
    assert seq1 == seq2
    assert a[0] == a[1] == a[2] == seq1[0]


def test_different_seeds_give_different_sequences():
    seq1 = [sample_matching(10, random.Random(1)) for _ in range(5)]
    seq2 = [sample_matching(10, random.Random(2)) for _ in range(5)]
    assert seq1 != seq2


def test_global_random_state_is_untouched():
    before = random.getstate()
    sample_matching(20, random.Random(5))
    assert random.getstate() == before


# -- state-level uniformity ----------------------------------------------------


@pytest.mark.parametrize("n, seed", [(2, 11), (3, 12), (4, 13)])
def test_state_frequencies_are_uniform_over_all(n, seed):
    states = list(iter_states(InitialConfiguration([n]), StateSpace.ALL))
    samples_per_state = 200
    rng = random.Random(seed)
    observed = Counter(
        sample_matching(n, rng) for _ in range(samples_per_state * len(states))
    )
    assert set(observed) == set(states)
    _, _, p = chi_square_gof(observed, {s: 1 / len(states) for s in states})
    assert p > ALPHA


def test_the_uniformity_test_detects_a_biased_sampler():
    # Negative control: redraw half of the matchings that do not pair ends 0
    # and 1, raising P(r[0] = 1) from 1/5 to 0.28.
    n, rng = 3, random.Random(14)
    states = list(iter_states(InitialConfiguration([n]), StateSpace.ALL))

    def biased():
        r = sample_matching(n, rng)
        if r[0] != 1 and rng.random() < 0.5:
            r = sample_matching(n, rng)
        return r

    observed = Counter(biased() for _ in range(200 * len(states)))
    assert chi_square_gof(observed, {s: 1 / len(states) for s in states})[2] < ALPHA
