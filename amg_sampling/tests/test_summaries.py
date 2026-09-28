from collections import Counter

import pytest

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.enumerate import exact_cycle_distribution
from amg_sampling.exact.formulas import (
    cycle_count_all,
    cycle_distribution_all,
    cycle_distribution_deranged,
    num_all_states,
    num_deranged_states,
)
from amg_sampling.exact.summaries import (
    cycle_count_distribution,
    cycle_length_multiplicity_distribution,
    largest_cycle_distribution,
)

ALL, DERANGED = StateSpace.ALL, StateSpace.DERANGED


def theorem_2(n, space):
    return cycle_distribution_all(n) if space is ALL else cycle_distribution_deranged(n)


def aggregate(dist, key):
    out = Counter()
    for c, m in dist.items():
        out[key(c)] += m
    return {k: v for k, v in out.items() if v}


SUMMARIES = {
    "cycles": (cycle_count_distribution, lambda c: c.num_cycles),
    "largest": (largest_cycle_distribution, lambda c: c.parts[0]),
}


@pytest.mark.parametrize("space", [ALL, DERANGED], ids=lambda s: s.value)
@pytest.mark.parametrize("n", range(1, 19))
def test_against_theorem_2_aggregated_over_partitions(n, space):
    dist = theorem_2(n, space)
    for func, key in SUMMARIES.values():
        assert func(n, space) == aggregate(dist, key)
    for length in range(2 if space is DERANGED else 1, n + 1):
        assert cycle_length_multiplicity_distribution(n, length, space) == aggregate(
            dist, lambda c: c.count(length)
        )


@pytest.mark.parametrize("space", [ALL, DERANGED], ids=lambda s: s.value)
@pytest.mark.parametrize("n", range(2, 7))
def test_against_enumeration(n, space):
    dist = exact_cycle_distribution(InitialConfiguration([n]), space)
    assert cycle_count_distribution(n, space) == aggregate(dist, lambda c: c.num_cycles)
    assert largest_cycle_distribution(n, space) == aggregate(dist, lambda c: c.parts[0])
    assert cycle_length_multiplicity_distribution(n, 2, space) == aggregate(dist, lambda c: c.count(2))


@pytest.mark.parametrize("n", [1, 10, 40, 80])
def test_cycle_count_all_agrees_with_stirling_formula(n):
    assert cycle_count_distribution(n, ALL) == {c: cycle_count_all(n, c) for c in range(1, n + 1)}


@pytest.mark.parametrize("n", [2, 25, 80])
def test_totals(n):
    for space, total in [(ALL, num_all_states(n)), (DERANGED, num_deranged_states(n))]:
        assert sum(cycle_count_distribution(n, space).values()) == total
        assert sum(largest_cycle_distribution(n, space).values()) == total
        assert sum(cycle_length_multiplicity_distribution(n, 2, space).values()) == total


def test_proper_is_rejected():
    with pytest.raises(ValueError):
        cycle_count_distribution(5, StateSpace.PROPER)
    with pytest.raises(ValueError):
        cycle_length_multiplicity_distribution(5, 1, DERANGED)
