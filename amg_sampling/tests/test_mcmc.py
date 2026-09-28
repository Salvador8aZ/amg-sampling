"""Validation of the Markov chain sampler (docs/theory/mcmc.md).

Level 0: the exact transition kernel, for every state of small state spaces.
Level 1: the fast chain draws its steps from that kernel.
Level 2: chain averages agree with exact distributions, using the batch-means
standard errors that the analyses report.
"""

import random
from collections import Counter
from fractions import Fraction

import pytest

from amg_sampling import directories
from amg_sampling.analysis.estimates import (
    ProportionEstimate,
    batch_means_ess,
    batch_means_standard_error,
)
from amg_sampling.analysis.null import num_cycles, sample_null
from amg_sampling.analysis.runs import Problem, SamplerSettings, run_observed
from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import cycle_structure
from amg_sampling.core.matching import (
    InvalidMatchingError,
    PartialMatching,
    RejoinMatching,
)
from amg_sampling.core.statespace import StateSpace
from amg_sampling.data.sheth import convert, default_data_path, load_dataset
from amg_sampling.exact.enumerate import exact_cycle_distribution, iter_states
from amg_sampling.exact.summaries import cycle_count_distribution
from amg_sampling.samplers.mcmc import iter_chain, transition_probabilities
from amg_sampling.samplers.rejection import SamplingStats
from amg_sampling.tests.paper_patient import TABLE_2
from amg_sampling.tests.reference import compositions
from amg_sampling.tests.stats import ALPHA, chi_square_gof

ALL, DERANGED, PROPER = StateSpace.ALL, StateSpace.DERANGED, StateSpace.PROPER
SPACES = (ALL, DERANGED, PROPER)
Z_LIMIT = 4.0
REAL = default_data_path()


def layouts(max_n):
    return [b for n in range(2, max_n + 1) for b in compositions(n)]


def check_kernel(theta, space, states, fixed=None):
    """Rows sum to 1, the kernel is symmetric, so columns sum to 1 as well."""
    in_space = set(states)
    kernel = {s: transition_probabilities(theta, space, s, fixed) for s in states}
    for s, row in kernel.items():
        assert sum(row.values()) == 1
        assert set(row) <= in_space
    for s, row in kernel.items():
        for t, p in row.items():
            assert kernel[t].get(s, 0) == p  # detailed balance for the uniform target
    for t in states:
        assert sum(kernel[s].get(t, 0) for s in states) == 1  # uniform is stationary


# -- level 0: the exact kernel ------------------------------------------------------


@pytest.mark.parametrize("breaks", layouts(4))
def test_kernel_is_symmetric_and_preserves_uniform(breaks):
    theta = InitialConfiguration(breaks)
    for space in SPACES:
        states = list(iter_states(theta, space))
        if states:
            check_kernel(theta, space, states)


@pytest.mark.parametrize(
    "breaks, pairs",
    [((3, 1), [(1, 6)]), ((2, 2), [(0, 5)]), ((4,), [(1, 2)]), ((2, 1, 1), [(0, 7)])],
)
def test_kernel_with_fixed_rejoins(breaks, pairs):
    theta = InitialConfiguration(breaks)
    fixed = PartialMatching(theta.num_ends, pairs)
    for space in SPACES:
        states = list(iter_states(theta, space, fixed=fixed))
        if states:
            check_kernel(theta, space, states, fixed)
            for s in states:
                assert all(
                    fixed.is_completed_by(t)
                    for t in transition_probabilities(theta, space, s, fixed)
                )


def test_kernel_of_a_single_state_chain():
    theta = InitialConfiguration((2,))
    r = RejoinMatching.from_pairs(4, [(0, 2), (1, 3)])
    fixed = PartialMatching(4, [(0, 2)])
    assert transition_probabilities(theta, ALL, r, fixed) == {r: Fraction(1)}
    with pytest.raises(ValueError, match="not in"):
        transition_probabilities(
            theta, DERANGED, RejoinMatching.from_pairs(4, [(0, 1), (2, 3)])
        )


# -- level 1: the chain's steps follow the kernel -----------------------------------


@pytest.mark.parametrize(
    "breaks, space, seed", [((2, 2), PROPER, 1), ((3, 1), DERANGED, 2), ((4,), ALL, 3)]
)
def test_one_step_distribution_matches_the_kernel(breaks, space, seed):
    theta = InitialConfiguration(breaks)
    start = next(iter(iter_states(theta, space)))
    kernel = transition_probabilities(theta, space, start)
    rng = random.Random(seed)
    counts = Counter(
        next(iter_chain(theta, space, rng, 1, initial=start)) for _ in range(20_000)
    )
    _, _, p_value = chi_square_gof(counts, kernel)
    assert p_value > ALPHA


# -- chain invariants ---------------------------------------------------------------


def test_chain_stays_in_the_space_and_keeps_fixed_edges():
    theta = InitialConfiguration((3, 2, 2))
    fixed = PartialMatching(theta.num_ends, [(0, 5), (7, 12)])
    stats = SamplingStats()
    states = list(
        iter_chain(theta, PROPER, random.Random(0), 500, fixed=fixed, stats=stats)
    )
    assert all(PROPER.contains(theta, r) for r in states)
    assert all(fixed.is_completed_by(r) for r in states)
    assert stats.proposals == 500 and 0 < stats.accepted < 500


def test_burn_in_and_thinning_count_steps():
    theta = InitialConfiguration((3, 3))
    stats = SamplingStats()
    list(
        iter_chain(
            theta, DERANGED, random.Random(1), 7, burn_in=11, thin=3, stats=stats
        )
    )
    assert stats.proposals == 11 + 7 * 3


def test_chain_is_reproducible():
    theta = InitialConfiguration((4, 2))
    a = list(iter_chain(theta, PROPER, random.Random(5), 50, thin=2))
    b = list(iter_chain(theta, PROPER, random.Random(5), 50, thin=2))
    assert a == b
    assert len(set(a)) > 1


def test_chain_argument_checks():
    theta = InitialConfiguration((2, 2))
    rng = random.Random(0)
    with pytest.raises(ValueError, match="burn_in"):
        next(iter_chain(theta, PROPER, rng, 1, burn_in=-1))
    with pytest.raises(ValueError, match="thin"):
        next(iter_chain(theta, PROPER, rng, 1, thin=0))
    with pytest.raises(TypeError):
        next(iter_chain(theta, PROPER, 0, 1))
    not_proper = RejoinMatching.from_pairs(8, [(0, 1), (2, 3), (4, 6), (5, 7)])
    with pytest.raises(ValueError, match="not in"):
        next(iter_chain(theta, PROPER, rng, 1, initial=not_proper))
    fixed = PartialMatching(8, [(0, 2)])
    start = RejoinMatching.from_pairs(8, [(0, 4), (1, 5), (2, 6), (3, 7)])
    with pytest.raises(InvalidMatchingError):
        next(iter_chain(theta, PROPER, rng, 1, initial=start, fixed=fixed))


def test_single_state_chain_repeats_it():
    theta = InitialConfiguration((2, 1))
    fixed = PartialMatching(6, [(0, 2), (1, 4)])
    states = set(iter_chain(theta, ALL, random.Random(0), 20, fixed=fixed))
    assert states == {RejoinMatching.from_pairs(6, [(0, 2), (1, 4), (3, 5)])}


# -- batch means --------------------------------------------------------------------


def batches(series, size):
    full = len(series) // size
    return [sum(series[b * size : (b + 1) * size]) for b in range(full)]


def test_batch_means_ess_of_iid_samples_is_close_to_n():
    rng = random.Random(3)
    n, size = 40_000, 200
    series = [rng.random() < 0.3 for _ in range(n)]
    ess = batch_means_ess(batches(series, size), size, sum(series), n)
    assert 0.7 * n < ess <= n


def test_batch_means_ess_of_a_sticky_chain():
    # Two-state chain that flips with probability q: the integrated
    # autocorrelation time is (1 - q) / q, so ESS ≈ N q / (1 - q).
    rng = random.Random(4)
    n, size, q = 200_000, 1_000, 0.1
    state, series = False, []
    for _ in range(n):
        if rng.random() < q:
            state = not state
        series.append(state)
    ess = batch_means_ess(batches(series, size), size, sum(series), n)
    assert ess == pytest.approx(n * q / (1 - q), rel=0.25)


def test_batch_means_edge_cases():
    assert batch_means_ess([5], 10, 5, 10) is None  # one batch
    assert batch_means_ess([0, 0], 10, 0, 20) is None  # no hits
    assert batch_means_ess([5, 5], 10, 10, 20) is None  # zero batch variance
    assert batch_means_standard_error([1.0]) is None
    assert batch_means_standard_error([1.0, 3.0]) == pytest.approx(1.0)


def test_proportion_estimate_with_effective_trials():
    iid = ProportionEstimate(300, 1000)
    mcmc = ProportionEstimate(300, 1000, effective_trials=250)
    assert mcmc.standard_error == pytest.approx(2 * iid.standard_error)
    assert "ESS 250" in ProportionEstimate(3, 1000, effective_trials=250).describe()
    assert mcmc.to_dict()["effective_samples"] == 250
    with pytest.raises(ValueError, match="effective_trials"):
        ProportionEstimate(3, 10, effective_trials=11)


# -- level 2: chain averages against exact distributions ----------------------------


def assert_matches(counts, exact, ess_of):
    total = sum(counts.values())
    exact_total = sum(exact.values())
    for key, m in exact.items():
        p = m / exact_total
        est = ProportionEstimate(
            counts.get(key, 0), total, effective_trials=ess_of(key)
        )
        se = (p * (1 - p) / est.n_eff) ** 0.5
        assert abs(est.estimate - p) < Z_LIMIT * se, (key, est.estimate, p)
    assert set(counts) <= set(exact)


@pytest.mark.parametrize(
    "breaks, space, seed",
    [
        ((3, 2), PROPER, 11),
        ((2, 2, 1), PROPER, 12),
        ((5,), DERANGED, 13),
        ((4, 1), ALL, 14),
    ],
)
def test_chain_reproduces_exact_cycle_distributions(breaks, space, seed):
    theta = InitialConfiguration(breaks)
    sample = sample_null(
        theta, space, random.Random(seed), 20_000, method="mcmc", burn_in=500, thin=3
    )
    ess = sample.ess_by_key()
    default = sample.reference_ess()
    assert_matches(
        sample.cycle_structures,
        exact_cycle_distribution(theta, space),
        lambda c: ess.get(c, default),
    )


def test_chain_converges_from_an_atypical_start():
    # The single-cycle state C_6 of Θ(1,(6)) has probability 1/10395 under
    # Uniform(ALL); after burn-in the chain has forgotten it.
    theta = InitialConfiguration((6,))
    start = RejoinMatching.from_pairs(
        12, [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 0)]
    )
    assert str(cycle_structure(theta, start)) == "C6"
    counts = Counter(
        num_cycles(cycle_structure(theta, r))
        for r in iter_chain(
            theta, ALL, random.Random(21), 20_000, burn_in=2_000, thin=5, initial=start
        )
    )
    exact = cycle_count_distribution(6, ALL)
    assert_matches(counts, exact, lambda k: None)


def test_observed_analysis_with_mcmc_matches_enumeration():
    ds = load_dataset(directories.test_data("synthetic_sheth.csv"))
    pc = convert(ds.patients["SYN-1"], (1, 2))
    problem = Problem("patient", pc.theta, "SYN-1", pc, ds.source)
    sampler = SamplerSettings("mcmc", 20_000, burn_in=500, thin=3)
    results = run_observed(problem, PROPER, sampler, seed=7)
    assert results["null"]["sampling"]["method"] == "mcmc"
    for row in results["null_cycle_distribution"]:
        p = row["exact_probability"]["value"]
        se = (p * (1 - p) / row["effective_samples"]) ** 0.5
        assert abs(row["sampled_probability"] - p) < Z_LIMIT * se
    for row in results["rejoins"]:
        est = row["null_probability"]
        assert est["effective_samples"] is not None
        assert est["effective_samples"] <= est["samples"]


def test_chain_at_n40_matches_the_exact_number_of_cycles():
    theta = InitialConfiguration((40,))
    sample = sample_null(
        theta,
        DERANGED,
        random.Random(40),
        20_000,
        method="mcmc",
        burn_in=2_000,
        thin=20,
    )
    counts = Counter()
    for c, m in sample.cycle_structures.items():
        counts[c.num_cycles] += m
    exact = cycle_count_distribution(40, DERANGED)
    ess = sample.ess_by_key(num_cycles)
    default = sample.reference_ess()
    # Only values with a non-negligible probability are compared.
    total = sum(exact.values())
    likely = {k: m for k, m in exact.items() if m / total > 1e-3}
    counts = Counter({k: m for k, m in counts.items() if k in likely})
    rest = sample.num_samples - sum(counts.values())
    likely["rest"] = total - sum(likely.values())
    counts["rest"] = rest
    assert_matches(counts, likely, lambda k: ess.get(k, default))


@pytest.mark.skipif(not REAL.is_file(), reason=f"Sheth dataset not available at {REAL}")
def test_mcmc_reproduces_the_p05_1657_completions():
    pc = convert(load_dataset(REAL).patients["P05-1657"], (8, 12))
    sample = sample_null(
        pc.theta,
        PROPER,
        random.Random(1657),
        30_000,
        fixed=pc.fixed,
        method="mcmc",
        burn_in=1_000,
        thin=3,
    )
    counts = {str(c): m for c, m in sample.cycle_structures.items()}
    ess = {str(c): e for c, e in sample.ess_by_key().items()}
    default = sample.reference_ess()
    assert_matches(counts, TABLE_2, lambda k: ess.get(k, default))
