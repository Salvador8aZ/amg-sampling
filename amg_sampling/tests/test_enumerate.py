from fractions import Fraction

import pytest

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import CycleStructure
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.enumerate import (
    count_states,
    exact_cycle_distribution,
    iter_states,
    probabilities,
)
from amg_sampling.tests.paper_tables import TABLE_1, TABLE_3
from amg_sampling.tests.reference import compositions, reference_states

SPACES = list(StateSpace)


def as_pair_sets(states):
    return [s.pairs() for s in states]


def by_string(dist):
    return {str(c): m for c, m in dist.items()}


# -- the state set equals the independent brute-force reference ---------------


@pytest.mark.parametrize("space", SPACES, ids=lambda s: s.value)
@pytest.mark.parametrize(
    "n", [1, 2, 3, 4, 5, pytest.param(6, marks=pytest.mark.slow)]
)
def test_state_set_matches_reference(n, space):
    # Every DSB distribution with n DSBs; compare the SET of states, not only counts.
    for breaks in compositions(n):
        theta = InitialConfiguration(breaks)
        produced = as_pair_sets(iter_states(theta, space))
        assert len(produced) == len(set(produced)), "a state was produced twice"
        assert set(produced) == reference_states(breaks, space.value)


@pytest.mark.parametrize("space", SPACES, ids=lambda s: s.value)
@pytest.mark.parametrize("breaks", [(4,), (2, 1, 1), (3, 2), (1, 1, 1, 1, 1)])
def test_every_state_is_a_member(breaks, space):
    theta = InitialConfiguration(breaks)
    for r in iter_states(theta, space):
        assert space.contains(theta, r)


def test_order_is_deterministic():
    theta = InitialConfiguration((2, 2))
    assert list(iter_states(theta, StateSpace.ALL)) == list(iter_states(theta, StateSpace.ALL))


def test_iter_states_is_lazy():
    theta = InitialConfiguration([7])  # 135135 states in ALL
    first = next(iter_states(theta, StateSpace.ALL))
    assert first.pairs() == ((0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11), (12, 13))


# -- counts --------------------------------------------------------------------


@pytest.mark.parametrize(
    "breaks, sizes",
    [
        # (|ALL|, |DERANGED|, |PROPER|)
        ((1,), (1, 0, 0)),
        ((2,), (3, 2, 2)),
        ((1, 1), (3, 2, 2)),
        ((4,), (105, 60, 60)),
        ((2, 2), (105, 60, 56)),
        ((1, 1, 1, 1), (105, 60, 48)),
        ((2, 2, 1), (945, 544, 512)),
        ((6,), (10395, 6040, 6040)),
        ((2, 2, 1, 1), (10395, 6040, 5568)),
    ],
)
def test_state_space_sizes(breaks, sizes):
    theta = InitialConfiguration(breaks)
    assert tuple(count_states(theta, s) for s in SPACES) == sizes


# -- paper tables, via the production enumerator --------------------------------


@pytest.mark.parametrize("breaks", list(TABLE_1) + list(TABLE_3))
def test_paper_tables(breaks):
    theta = InitialConfiguration(breaks)
    expected = {**TABLE_1, **TABLE_3}[breaks]
    assert by_string(exact_cycle_distribution(theta, StateSpace.PROPER)) == expected


# -- exact distributions -------------------------------------------------------


@pytest.mark.parametrize("space", SPACES, ids=lambda s: s.value)
@pytest.mark.parametrize("breaks", [(3,), (2, 2), (3, 1, 1), (2, 2, 1, 1)])
def test_distribution_sums_to_count_and_partitions_n(breaks, space):
    theta = InitialConfiguration(breaks)
    dist = exact_cycle_distribution(theta, space)
    assert sum(dist.values()) == count_states(theta, space)
    for c, m in dist.items():
        assert isinstance(m, int) and m > 0
        assert c.n == theta.num_dsbs


def test_distribution_includes_c1_in_all():
    theta = InitialConfiguration([4])
    dist = by_string(exact_cycle_distribution(theta, StateSpace.ALL))
    assert dist == {"C4": 48, "C3+C1": 32, "2C2": 12, "C2+2C1": 12, "4C1": 1}


def test_probabilities_are_exact_fractions():
    theta = InitialConfiguration((2, 2))
    p = probabilities(exact_cycle_distribution(theta, StateSpace.PROPER))
    assert p == {CycleStructure([4]): Fraction(6, 7), CycleStructure([2, 2]): Fraction(1, 7)}
    assert sum(p.values()) == 1


def test_probabilities_of_empty_space_raise():
    theta = InitialConfiguration([1])
    with pytest.raises(ValueError):
        probabilities(exact_cycle_distribution(theta, StateSpace.PROPER))


# -- does the cycle-structure distribution depend on the chromosome layout? ----

LAYOUTS = {
    4: [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)],
    5: [(5,), (3, 2), (3, 1, 1), (2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)],
}


@pytest.mark.parametrize("space", [StateSpace.ALL, StateSpace.DERANGED], ids=lambda s: s.value)
@pytest.mark.parametrize("n", [4, 5])
def test_all_and_deranged_do_not_depend_on_layout(n, space):
    # PROVED: ALL and DERANGED are the same set of matchings for every Θ with
    # the same n, and the cycle structure depends only on (d, r).
    dists = [exact_cycle_distribution(InitialConfiguration(b), space) for b in LAYOUTS[n]]
    assert all(d == dists[0] for d in dists)
    states = [set(iter_states(InitialConfiguration(b), space)) for b in LAYOUTS[n]]
    assert all(s == states[0] for s in states)


@pytest.mark.parametrize("n", [4, 5])
def test_proper_depends_on_layout(n):
    dists = {
        b: by_string(exact_cycle_distribution(InitialConfiguration(b), StateSpace.PROPER))
        for b in LAYOUTS[n]
    }
    assert len({tuple(sorted(d.items())) for d in dists.values()}) > 1
    # Normalised probabilities differ too, not only the totals.
    p = {b: probabilities(exact_cycle_distribution(InitialConfiguration(b), StateSpace.PROPER))
         for b in LAYOUTS[n]}
    assert len({tuple(sorted((str(c), q) for c, q in v.items())) for v in p.values()}) > 1


@pytest.mark.parametrize("n", [4, 5])
def test_single_chromosome_proper_equals_deranged(n):
    theta = InitialConfiguration([n])
    assert set(iter_states(theta, StateSpace.PROPER)) == set(iter_states(theta, StateSpace.DERANGED))
