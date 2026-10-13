"""Exact weighted cycle-structure distributions (exact/weighted.py)."""

import math
import random
from collections import Counter

import pytest

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import cycle_structure
from amg_sampling.core.statespace import StateSpace
from amg_sampling.data.sheth import convert, load_dataset
from amg_sampling.exact.enumerate import (
    exact_cycle_distribution,
    iter_states,
    probabilities,
)
from amg_sampling.exact.weighted import total_variation, weighted_cycle_distribution
from amg_sampling.samplers.targets import ProximityTarget
from amg_sampling.tests.paper_patient import TABLE_2
from amg_sampling.tests.test_app import REAL
from amg_sampling.tests.test_targets import synthetic_target, weighted_distribution

PROPER = StateSpace.PROPER


@pytest.mark.parametrize("breaks", [(4,), (2, 2), (3, 2), (2, 2, 1)])
def test_uniform_matches_the_exact_counts(breaks):
    theta = InitialConfiguration(list(breaks))
    got = weighted_cycle_distribution(theta, PROPER, None)
    exact = probabilities(exact_cycle_distribution(theta, PROPER))
    assert got.probabilities == pytest.approx({c: float(p) for c, p in exact.items()})
    assert got.num_states == sum(exact_cycle_distribution(theta, PROPER).values())
    assert got.num_map_states == got.num_states  # every state ties
    assert got.log_partition == pytest.approx(math.log(got.num_states))


@pytest.mark.parametrize(
    "breaks, alpha", [((5,), 1.0), ((3, 2), 2.0), ((2, 2, 2), 0.5)]
)
def test_weighted_matches_brute_force(breaks, alpha):
    theta = InitialConfiguration(list(breaks))
    target = synthetic_target(theta, alpha, seed=7)
    got = weighted_cycle_distribution(theta, PROPER, target)
    brute = weighted_distribution(theta, PROPER, target)
    total = sum(brute.values())
    assert got.probabilities == pytest.approx({c: w / total for c, w in brute.items()})
    assert got.log_partition == pytest.approx(math.log(total))
    # the MAP state is a state of maximum weight
    weights = {s: target.log_weight(s.partners) for s in iter_states(theta, PROPER)}
    best = max(weights.values())
    assert weights[got.map_state] == pytest.approx(best)
    assert got.map_probability == pytest.approx(math.exp(best) / total)
    assert got.num_map_states == sum(
        1 for w in weights.values() if abs(w - best) < 1e-9
    )


def test_probabilities_are_sorted_and_sum_to_one():
    theta = InitialConfiguration([3, 3])
    got = weighted_cycle_distribution(theta, PROPER, synthetic_target(theta, 1.5))
    values = list(got.probabilities.values())
    assert values == sorted(values, reverse=True)
    assert sum(values) == pytest.approx(1.0)
    assert got.most_probable_structure == next(iter(got.probabilities))


def test_large_alpha_does_not_overflow():
    theta = InitialConfiguration([4])
    got = weighted_cycle_distribution(theta, PROPER, synthetic_target(theta, 400.0))
    assert sum(got.probabilities.values()) == pytest.approx(1.0)
    assert got.map_probability <= 1.0


def test_most_probable_structure_can_differ_from_the_map_structure():
    # Hand-made weights: Θ(1,(4)) in ALL has 105 states, 48 of them C4 and 12
    # of them 2C2. One 2C2 state weighted 3x is the MAP state, but C4 keeps
    # the most mass.
    theta = InitialConfiguration([4])
    states = list(iter_states(theta, StateSpace.ALL))
    by_structure = Counter(str(cycle_structure(theta, s)) for s in states)
    heavy = next(s for s in states if str(cycle_structure(theta, s)) == "2C2")

    class OneHeavy:
        def log_weight(self, partners):
            return math.log(3.0) if tuple(partners) == heavy.partners else 0.0

    got = weighted_cycle_distribution(theta, StateSpace.ALL, OneHeavy())
    assert str(got.map_structure(theta)) == "2C2"
    assert got.num_map_states == 1
    assert by_structure["2C2"] + 2 < by_structure["C4"]
    assert str(got.most_probable_structure) == "C4"


def test_total_variation():
    assert total_variation({"a": 0.5, "b": 0.5}, {"a": 0.5, "b": 0.5}) == 0
    assert total_variation({"a": 1.0}, {"b": 1.0}) == 1
    assert total_variation({"a": 0.7, "b": 0.3}, {"a": 0.4, "c": 0.6}) == pytest.approx(
        0.6
    )


def test_inter_chromosomal_weight_changes_the_distribution():
    theta = InitialConfiguration([2, 2])
    base = synthetic_target(theta, 1.0, seed=2)
    inter = ProximityTarget(
        base.chromosomes, base.coordinates, 1.0, inter_log_weight=-3.0
    )
    a = weighted_cycle_distribution(theta, PROPER, base)
    b = weighted_cycle_distribution(theta, PROPER, inter)
    assert (
        total_variation(
            {str(c): p for c, p in a.probabilities.items()},
            {str(c): p for c, p in b.probabilities.items()},
        )
        > 0.01
    )


@pytest.mark.skipif(not REAL.is_file(), reason=f"Sheth dataset not available at {REAL}")
def test_p05_1657_uniform_reproduces_table_2():
    pc = convert(load_dataset(REAL).patients["P05-1657"], (8, 12))
    got = weighted_cycle_distribution(pc.theta, PROPER, None, pc.fixed)
    assert got.num_states == 945
    assert {str(c): round(p * 945) for c, p in got.probabilities.items()} == TABLE_2


def test_chain_agrees_with_weighted_exact_with_an_inter_weight():
    from amg_sampling.analysis.null import sample_null
    from amg_sampling.tests.test_mcmc import assert_matches

    theta = InitialConfiguration([3, 2])
    base = synthetic_target(theta, 1.0, seed=5)
    target = ProximityTarget(
        base.chromosomes, base.coordinates, 1.0, inter_log_weight=-2.0
    )
    exact = weighted_distribution(theta, PROPER, target)
    sample = sample_null(
        theta,
        PROPER,
        random.Random(3),
        20_000,
        method="mcmc",
        burn_in=2_000,
        thin=5,
        target=target,
    )
    ess = sample.ess_by_key()
    assert_matches(
        sample.cycle_structures, exact, lambda k: ess.get(k, sample.reference_ess())
    )
