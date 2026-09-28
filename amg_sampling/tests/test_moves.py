"""The 2-switch move: neighbour counts, symmetry and connectivity of state spaces.

The connectivity checks are the ones the Markov chain relies on: a chain that
cannot reach every state of its state space would be silently wrong.
"""

from collections import Counter
from itertools import combinations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import (
    InvalidMatchingError,
    PartialMatching,
    RejoinMatching,
)
from amg_sampling.core.moves import (
    movable_ends,
    num_switch_neighbours,
    switch,
    switch_in_place,
    switch_neighbours,
)
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.enumerate import iter_states
from amg_sampling.tests.reference import compositions
from amg_sampling.tests.strategies import matchings, states

SPACES = (StateSpace.ALL, StateSpace.DERANGED, StateSpace.PROPER)


def edges(matching):
    return set(matching.pairs())


def ordered_pairs(matching, movable):
    r = matching.partners
    return [(u, v) for u in movable for v in movable if v not in (u, r[u])]


def num_components(states_, movable_of):
    """Connected components of the switch graph restricted to ``states_``."""
    remaining = set(states_)
    count = 0
    while remaining:
        count += 1
        stack = [remaining.pop()]
        while stack:
            for y in switch_neighbours(stack.pop(), movable_of):
                if y in remaining:
                    remaining.remove(y)
                    stack.append(y)
    return count


def partial_matchings(ends, k):
    """Every set of ``k`` disjoint pairs of ``ends``."""
    if k == 0:
        yield ()
        return
    if len(ends) < 2 * k:
        return
    first, rest = ends[0], ends[1:]
    yield from partial_matchings(rest, k)
    for i, w in enumerate(rest):
        for pm in partial_matchings(rest[:i] + rest[i + 1 :], k - 1):
            yield ((first, w), *pm)


# -- the move ---------------------------------------------------------------------


def test_switch_example():
    r = RejoinMatching.from_pairs(6, [(0, 2), (1, 4), (3, 5)])
    s = switch(r, 0, 1)
    assert edges(s) == {(0, 1), (2, 4), (3, 5)}
    assert r.pairs() == ((0, 2), (1, 4), (3, 5))  # the original is unchanged


def test_switch_rejects_ends_of_the_same_edge():
    r = RejoinMatching.from_pairs(4, [(0, 2), (1, 3)])
    for v in (0, 2):
        with pytest.raises(ValueError, match="same rejoin edge"):
            switch(r, 0, v)
    with pytest.raises(IndexError):
        switch(r, 0, 4)


@given(matchings(12), st.data())
def test_switch_is_reversible(r, data):
    u = data.draw(st.sampled_from(range(12)))
    v = data.draw(st.sampled_from([x for x in range(12) if x not in (u, r[u])]))
    s = switch(r, u, v)
    assert s[u] == v and s[r[u]] == r[v]
    assert len(edges(r) - edges(s)) == 2
    assert switch(s, u, r[u]) == r


@given(matchings(10), st.data())
def test_switch_in_place_matches_switch(r, data):
    u = data.draw(st.sampled_from(range(10)))
    v = data.draw(st.sampled_from([x for x in range(10) if x not in (u, r[u])]))
    partners = list(r.partners)
    switch_in_place(partners, u, v)
    assert RejoinMatching(partners) == switch(r, u, v)


# -- neighbour counts and proposal symmetry (docs/theory/mcmc.md §1) ------------


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_every_matching_has_m_m_minus_1_distinct_neighbours(n):
    for r in iter_states(InitialConfiguration([n]), StateSpace.ALL):
        neighbours = list(switch_neighbours(r))
        assert len(neighbours) == len(set(neighbours)) == num_switch_neighbours(n)
        # Each neighbour differs in exactly two edges, so it is never r itself.
        assert all(len(edges(r) - edges(s)) == 2 for s in neighbours)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_each_neighbour_arises_from_exactly_four_ordered_pairs(n):
    # This is what makes the proposal "u uniform, then v uniform among the
    # ends not on u's edge" symmetric: P(r -> s) = 4 / (2n (2n - 2)).
    for r in iter_states(InitialConfiguration([n]), StateSpace.ALL):
        movable = range(2 * n)
        counts = Counter(switch(r, u, v) for u, v in ordered_pairs(r, movable))
        assert set(counts) == set(switch_neighbours(r))
        assert set(counts.values()) == {4}


@pytest.mark.parametrize("n", [3, 4, 5])
def test_moves_with_fixed_edges_keep_them(n):
    theta = InitialConfiguration([n])
    fixed = PartialMatching(2 * n, [(0, 3)])
    for r in iter_states(theta, StateSpace.ALL):
        if not fixed.is_completed_by(r):
            continue
        movable = movable_ends(r, fixed)
        assert movable == fixed.free_ends
        neighbours = list(switch_neighbours(r, movable))
        assert len(set(neighbours)) == num_switch_neighbours(n - 1)
        assert all(fixed.is_completed_by(s) for s in neighbours)
        counts = Counter(switch(r, u, v) for u, v in ordered_pairs(r, movable))
        assert set(counts.values()) == {4}


def test_movable_ends_requires_a_completion():
    fixed = PartialMatching(4, [(0, 1)])
    with pytest.raises(InvalidMatchingError):
        movable_ends(RejoinMatching.from_pairs(4, [(0, 2), (1, 3)]), fixed)
    r = RejoinMatching.from_pairs(4, [(0, 2), (1, 3)])
    assert movable_ends(r) == (0, 1, 2, 3)
    with pytest.raises(ValueError, match="union of rejoin edges"):
        list(switch_neighbours(r, [0, 1]))


def test_num_switch_neighbours():
    assert [num_switch_neighbours(m) for m in range(5)] == [0, 0, 2, 6, 12]
    with pytest.raises(ValueError):
        num_switch_neighbours(-1)


# -- connectivity of the restricted switch graphs (docs/theory/mcmc.md §2) -------


def layouts(max_n):
    return [b for n in range(1, max_n + 1) for b in compositions(n)]


@pytest.mark.parametrize("breaks", layouts(5))
def test_every_state_space_is_connected_under_switches(breaks):
    theta = InitialConfiguration(breaks)
    for space in SPACES:
        states_ = list(iter_states(theta, space))
        if states_:
            assert num_components(states_, None) == 1, space


@pytest.mark.slow
@pytest.mark.parametrize("breaks", list(compositions(6)))
def test_every_state_space_is_connected_under_switches_n6(breaks):
    theta = InitialConfiguration(breaks)
    for space in SPACES:
        assert num_components(list(iter_states(theta, space)), None) == 1, space


def completion_spaces(max_n):
    """(Θ, fixed, space, states) for every layout, fixed set and space."""
    for breaks in layouts(max_n):
        theta = InitialConfiguration(breaks)
        n = theta.num_dsbs
        all_states = list(iter_states(theta, StateSpace.ALL))
        for k in range(1, n - 1):  # leave at least two free edges
            for pairs in partial_matchings(list(range(2 * n)), k):
                fixed = PartialMatching(2 * n, pairs)
                for space in (StateSpace.DERANGED, StateSpace.PROPER):
                    yield (
                        theta,
                        fixed,
                        space,
                        [
                            r
                            for r in all_states
                            if fixed.is_completed_by(r) and space.contains(theta, r)
                        ],
                    )


def test_completion_spaces_are_connected_under_switches():
    checked = 0
    for theta, fixed, space, states_ in completion_spaces(4):
        if states_:
            checked += 1
            assert num_components(states_, fixed.free_ends) == 1, (
                theta,
                fixed,
                space,
            )
    assert checked > 1000


@pytest.mark.slow
def test_completion_spaces_are_connected_under_switches_n5():
    for theta, fixed, space, states_ in completion_spaces(5):
        if states_:
            assert num_components(states_, fixed.free_ends) == 1


# -- the fast membership test used by the chain --------------------------------


@given(states(max_chromosomes=4, max_breaks=4))
def test_contains_partners_agrees_with_contains(state):
    theta, r = state
    for space in SPACES:
        assert space.contains_partners(theta, list(r.partners)) == space.contains(
            theta, r
        )


def test_pairs_helper_counts():
    # The helper that generates fixed-edge sets, against a direct count of pairs
    # of disjoint pairs of 8 ends.
    got = sum(1 for _ in partial_matchings(list(range(8)), 2))
    expected = sum(
        1 for a, b in combinations(combinations(range(8), 2), 2) if not set(a) & set(b)
    )
    assert got == expected
