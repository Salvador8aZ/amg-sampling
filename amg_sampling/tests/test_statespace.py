from collections import Counter

import networkx as nx
import pytest
from hypothesis import given

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import cycle_structure
from amg_sampling.core.matching import RejoinMatching
from amg_sampling.core.statespace import (
    StateSpace,
    chromosome_components,
    is_connected,
    is_deranged,
    is_proper,
)
from amg_sampling.tests.reference import (
    all_perfect_matchings,
    amg_multigraph,
    compositions,
    exchange_multigraph,
    reference_is_proper,
)
from amg_sampling.tests.paper_tables import TABLE_1, TABLE_3
from amg_sampling.tests.strategies import states


# -- deterministic -----------------------------------------------------------


def test_example_2_inversion_and_ring_are_both_proper():
    theta = InitialConfiguration([2])
    assert is_proper(theta, RejoinMatching.from_pairs(4, [(0, 2), (1, 3)]))
    assert is_proper(theta, RejoinMatching.from_pairs(4, [(0, 3), (1, 2)]))  # ring


def test_repaired_dsb_is_not_deranged():
    theta = InitialConfiguration([2])
    r = RejoinMatching.from_pairs(4, [(0, 1), (2, 3)])
    assert StateSpace.ALL.contains(theta, r)
    assert not StateSpace.DERANGED.contains(theta, r)
    assert not StateSpace.PROPER.contains(theta, r)


def test_deranged_but_disconnected_is_not_proper():
    # Θ(2,(2,2)) with only intra-chromosomal rejoins.
    theta = InitialConfiguration((2, 2))
    r = RejoinMatching.from_pairs(8, [(0, 2), (1, 3), (4, 6), (5, 7)])
    assert is_deranged(theta, r)
    assert not is_connected(theta, r)
    assert StateSpace.DERANGED.contains(theta, r)
    assert not StateSpace.PROPER.contains(theta, r)
    assert chromosome_components(theta, r) == (frozenset({0}), frozenset({1}))


def test_connected_but_not_deranged_is_not_proper():
    theta = InitialConfiguration((2, 1))
    r = RejoinMatching.from_pairs(6, [(0, 1), (2, 4), (3, 5)])
    assert is_connected(theta, r)
    assert not is_proper(theta, r)


def test_chromosome_components():
    theta = InitialConfiguration((1, 1, 1, 1))
    r = RejoinMatching.from_pairs(8, [(0, 6), (1, 7), (2, 5), (3, 4)])
    assert chromosome_components(theta, r) == (frozenset({0, 3}), frozenset({1, 2}))
    assert not is_connected(theta, r)


def test_incompatible_sizes_raise():
    with pytest.raises(ValueError):
        StateSpace.PROPER.contains(InitialConfiguration([2]), RejoinMatching([1, 0]))


# -- paper counts, via brute force over all matchings ------------------------
# These check the PROPER predicate (and cycle structures) against the paper.
# They use a brute-force test helper, not a library enumerator.

def proper_distribution(breaks):
    theta = InitialConfiguration(breaks)
    dist = Counter()
    for pairs in all_perfect_matchings(theta.num_ends):
        r = RejoinMatching.from_pairs(theta.num_ends, pairs)
        if is_proper(theta, r):
            dist[str(cycle_structure(theta, r))] += 1
    return dict(dist)


@pytest.mark.parametrize("breaks", list(TABLE_1))
def test_table_1_proper_counts(breaks):
    assert proper_distribution(breaks) == TABLE_1[breaks]


@pytest.mark.parametrize("breaks", list(TABLE_3))
def test_table_3_proper_counts(breaks):
    assert proper_distribution(breaks) == TABLE_3[breaks]


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_all_and_deranged_counts(n):
    # |ALL| = (2n-1)!! (paper theorem 1); |DERANGED| is OEIS A053871.
    theta = InitialConfiguration([n])
    counts = Counter()
    for pairs in all_perfect_matchings(2 * n):
        r = RejoinMatching.from_pairs(2 * n, pairs)
        for space in StateSpace:
            counts[space] += space.contains(theta, r)
    double_factorial = [1, 1, 3, 15, 105, 945][n]
    assert counts[StateSpace.ALL] == double_factorial
    assert counts[StateSpace.DERANGED] == [1, 0, 2, 8, 60, 544][n]
    # Θ(1,(n)) has one chromosome, so every deranged state is proper.
    assert counts[StateSpace.PROPER] == counts[StateSpace.DERANGED]


# -- chromosome-connectivity criterion vs explicit graphs ---------------------


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, pytest.param(6, marks=pytest.mark.slow)])
def test_criterion_agrees_with_explicit_graph_exhaustively(n):
    # Every DSB distribution with n <= 6 and every perfect matching.
    for breaks in compositions(n):
        theta = InitialConfiguration(breaks)
        for pairs in all_perfect_matchings(2 * n):
            r = RejoinMatching.from_pairs(2 * n, pairs)
            assert is_connected(theta, r) == nx.is_connected(amg_multigraph(breaks, pairs))
            assert is_proper(theta, r) == reference_is_proper(breaks, pairs)


@given(states(max_chromosomes=8, max_breaks=5))
def test_criterion_agrees_with_explicit_graph(state):
    theta, r = state
    pairs = r.pairs()
    assert is_connected(theta, r) == nx.is_connected(amg_multigraph(theta.breaks, pairs))
    assert is_proper(theta, r) == reference_is_proper(theta.breaks, pairs)


@given(states())
def test_deranged_iff_no_c1(state):
    theta, r = state
    assert is_deranged(theta, r) == (cycle_structure(theta, r).count(1) == 0)


@given(states())
def test_state_spaces_are_nested(state):
    theta, r = state
    proper = StateSpace.PROPER.contains(theta, r)
    deranged = StateSpace.DERANGED.contains(theta, r)
    assert StateSpace.ALL.contains(theta, r)
    assert not proper or deranged


@given(states())
def test_chromosome_components_partition_the_chromosomes(state):
    theta, r = state
    comps = chromosome_components(theta, r)
    assert sorted(i for c in comps for i in c) == list(range(theta.num_chromosomes))
    assert (len(comps) == 1) == is_connected(theta, r)
    for v in range(theta.num_ends):
        same = [c for c in comps if theta.chromosome_of(v) in c]
        assert theta.chromosome_of(r[v]) in same[0]
