"""Validation of the proximity target for the Markov chain (docs/theory/mcmc.md §5).

Level 0: the exact kernel satisfies detailed balance for the weighted target.
Level 1: alpha = 0 replays the uniform chain step for step.
Level 2: chain frequencies match exact weighted enumeration.
"""

import math
import random
from collections import Counter

import pytest

from amg_sampling.analysis.null import sample_null
from amg_sampling.analysis.report import format_report
from amg_sampling.app import build_problem, execute
from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import cycle_structure
from amg_sampling.core.matching import PartialMatching
from amg_sampling.core.statespace import StateSpace
from amg_sampling.data.sheth import convert, load_dataset
from amg_sampling.exact.enumerate import iter_states
from amg_sampling.samplers.mcmc import iter_chain, transition_probabilities
from amg_sampling.samplers.targets import ProximityTarget
from amg_sampling.tests.test_app import REAL, SYNTHETIC, config
from amg_sampling.tests.test_mcmc import assert_matches

ALL, DERANGED, PROPER = StateSpace.ALL, StateSpace.DERANGED, StateSpace.PROPER


def synthetic_target(theta, alpha, seed=0):
    """Random breakpoint positions (bp) for every DSB of ``theta``."""
    rng = random.Random(seed)
    chromosomes, coordinates = [], []
    for chrom, breaks in enumerate(theta.breaks, start=1):
        positions = sorted(rng.sample(range(10_000, 50_000_000), breaks))
        for p in positions:
            chromosomes += [chrom, chrom]
            coordinates += [p, p + 1]
    return ProximityTarget(tuple(chromosomes), tuple(coordinates), alpha)


def weighted_distribution(theta, space, target, fixed=None):
    dist = Counter()
    for s in iter_states(theta, space, fixed):
        dist[cycle_structure(theta, s)] += math.exp(target.log_weight(s.partners))
    return dist


# -- the weights --------------------------------------------------------------------


def test_edge_weights():
    t = ProximityTarget((1, 1, 1, 1, 2, 2), (0, 1, 2_000_000, 2_000_001, 5, 6), 1.0)
    assert t.log_edge_weight(0, 2) == pytest.approx(-math.log(2.0))  # 2 Mb apart
    assert t.log_edge_weight(0, 4) == 0.0  # different chromosomes
    assert t.log_edge_weight(0, 1) == pytest.approx(-math.log(1_000 / 1_000_000))
    assert t.log_weight([2, 4, 0, 5, 1, 3]) == pytest.approx(
        t.log_edge_weight(0, 2) + t.log_edge_weight(1, 4) + t.log_edge_weight(3, 5)
    )
    zero = ProximityTarget(t.chromosomes, t.coordinates, 0.0)
    assert all(zero.log_edge_weight(u, v) == 0 for u in range(6) for v in range(6))


def test_invalid_parameters():
    with pytest.raises(ValueError, match="alpha"):
        ProximityTarget((1, 1), (0, 1), -1.0)
    with pytest.raises(ValueError, match="equal length"):
        ProximityTarget((1,), (0, 1), 1.0)
    with pytest.raises(ValueError, match="positive"):
        ProximityTarget((1, 1), (0, 1), 1.0, scale_bp=0)


# -- level 0: detailed balance --------------------------------------------------------


@pytest.mark.parametrize("breaks", [(4,), (2, 2), (3, 1), (2, 1, 1), (5,), (3, 2)])
@pytest.mark.parametrize("space", [ALL, DERANGED, PROPER])
def test_kernel_satisfies_detailed_balance(breaks, space):
    theta = InitialConfiguration(list(breaks))
    target = synthetic_target(theta, alpha=1.5)
    states = list(iter_states(theta, space))
    pi = {s: math.exp(target.log_weight(s.partners)) for s in states}
    kernel = {
        s: transition_probabilities(theta, space, s, target=target) for s in states
    }
    for s, row in kernel.items():
        assert sum(row.values()) == pytest.approx(1.0)
        for t, p in row.items():
            assert pi[s] * p == pytest.approx(pi[t] * kernel[t].get(s, 0))


def test_kernel_with_fixed_rejoins():
    theta = InitialConfiguration([3, 2])
    target = synthetic_target(theta, alpha=1.0)
    fixed = PartialMatching(theta.num_ends, [(1, 6)])
    states = list(iter_states(theta, PROPER, fixed))
    pi = {s: math.exp(target.log_weight(s.partners)) for s in states}
    for s in states:
        row = transition_probabilities(theta, PROPER, s, fixed, target)
        assert set(row) <= set(states)
        for t, p in row.items():
            row_t = transition_probabilities(theta, PROPER, t, fixed, target)
            assert pi[s] * p == pytest.approx(pi[t] * row_t.get(s, 0))


# -- level 1: alpha = 0 is the uniform chain ------------------------------------------


@pytest.mark.parametrize("space", [ALL, DERANGED, PROPER])
def test_alpha_zero_replays_the_uniform_chain(space):
    theta = InitialConfiguration([3, 3, 2])
    target = synthetic_target(theta, alpha=0.0)
    uniform = list(iter_chain(theta, space, random.Random(7), 500, burn_in=50))
    weighted = list(
        iter_chain(theta, space, random.Random(7), 500, burn_in=50, target=target)
    )
    assert weighted == uniform


# -- level 2: chain against exact weighted enumeration ---------------------------------


@pytest.mark.parametrize(
    "breaks, space, alpha",
    [((5,), PROPER, 1.0), ((3, 3), PROPER, 1.0), ((2, 2, 2), DERANGED, 2.0)],
)
def test_chain_matches_exact_weighted_distribution(breaks, space, alpha):
    theta = InitialConfiguration(list(breaks))
    target = synthetic_target(theta, alpha, seed=3)
    exact = weighted_distribution(theta, space, target)
    sample = sample_null(
        theta,
        space,
        random.Random(11),
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


def test_sample_null_rejects_iid_with_a_target():
    theta = InitialConfiguration([2, 2])
    with pytest.raises(ValueError, match="mcmc"):
        sample_null(
            theta, PROPER, random.Random(0), 10, target=synthetic_target(theta, 1.0)
        )


@pytest.mark.skipif(not REAL.is_file(), reason=f"Sheth dataset not available at {REAL}")
def test_proximity_chain_on_the_p05_1657_completions():
    pc = convert(load_dataset(REAL).patients["P05-1657"], (8, 12))
    target = ProximityTarget.from_patient(pc, alpha=1.0)
    exact = weighted_distribution(pc.theta, PROPER, target, pc.fixed)
    sample = sample_null(
        pc.theta,
        PROPER,
        random.Random(1657),
        30_000,
        fixed=pc.fixed,
        method="mcmc",
        burn_in=2_000,
        thin=5,
        target=target,
    )
    ess = sample.ess_by_key()
    assert_matches(
        sample.cycle_structures, exact, lambda k: ess.get(k, sample.reference_ess())
    )


# -- the application ------------------------------------------------------------------


PATIENT = (
    "problem=patient",
    "problem.patient_id=SYN-1",
    "problem.chromosomes=[1,2]",
    f"data.sheth_csv={SYNTHETIC}",
)


def test_default_target_is_uniform():
    cfg = config()
    assert cfg.target.name == "uniform"


def test_proximity_run_on_synthetic_data(tmp_path):
    cfg = config(
        *PATIENT,
        "sampler=mcmc",
        "sampler.num_samples=2000",
        "sampler.burn_in=500",
        "target=proximity",
        "target.alpha=1.5",
        "analysis=observed",
        "output.plots=false",
    )
    results, _ = execute(cfg, tmp_path)
    assert results["target"]["name"] == "proximity"
    assert results["target"]["alpha"] == 1.5
    assert results["null_model"].startswith("Proximity(PROPER")
    assert results["null"]["exact"]["unavailable"] == {
        "all": "exact results are computed for the uniform null only"
    }
    report = format_report(results)
    assert "proximity-weighted target" in report
    assert "under the proximity target for Θ" in report
    assert "uniform combinatorial null" not in report


@pytest.mark.parametrize(
    "overrides, message",
    [
        (("sampler=iid", "target=proximity", *PATIENT), "sampler=mcmc"),
        (("sampler=mcmc", "target=proximity"), "problem=patient"),
    ],
)
def test_proximity_target_errors(tmp_path, overrides, message):
    cfg = config(*overrides, "sampler.num_samples=10", "output.plots=false")
    with pytest.raises(ValueError, match=message):
        execute(cfg, tmp_path)


def test_uniform_results_are_unchanged(tmp_path):
    cfg = config(*PATIENT, "sampler=mcmc", "sampler.num_samples=500")
    assert build_problem(cfg).patient is not None
    results, _ = execute(cfg, tmp_path)
    assert results["target"] == {"name": "uniform"}
    assert results["null_model"] == "Uniform(PROPER(Θ))"
    assert results["exact"]["total_states"] is not None  # uniform exact reference kept
