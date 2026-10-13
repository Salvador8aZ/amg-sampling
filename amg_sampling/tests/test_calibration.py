"""Fitting the proximity target to observed junctions (analysis/calibration.py)."""

import math
import random

import pytest

from amg_sampling.analysis.calibration import (
    PARAMETER_LIMIT,
    Choice,
    fit,
    log_likelihood,
    patient_choices,
    split_patients,
)
from amg_sampling.data.sheth import load_dataset
from amg_sampling.tests.test_app import SYNTHETIC

IDS = [f"P{i:02d}" for i in range(57)]


# -- the split ---------------------------------------------------------------------


def test_split_is_deterministic_disjoint_and_sized():
    train, test = split_patients(IDS, 0.2, seed=1)
    assert (train, test) == split_patients(reversed(IDS), 0.2, seed=1)
    assert len(test) == round(0.2 * 57) == 11
    assert set(train).isdisjoint(test)
    assert sorted(train + test) == IDS
    assert split_patients(IDS, 0.2, seed=2)[1] != test


def test_split_forces_patients_into_the_test_set():
    for seed in range(20):
        train, test = split_patients(IDS, 0.2, seed=seed, force_test=["P05"])
        assert "P05" in test and "P05" not in train
        assert len(test) == 11


def test_split_rejects_bad_input():
    with pytest.raises(ValueError, match="test_fraction"):
        split_patients(IDS, 1.0)
    with pytest.raises(ValueError, match="Unknown"):
        split_patients(IDS, 0.2, force_test=["nope"])


# -- choices from a patient ----------------------------------------------------------


def test_patient_choices_on_synthetic_data():
    ds = load_dataset(SYNTHETIC)
    syn1 = patient_choices(ds.patients["SYN-1"])
    # Two junctions, each seen from both ends.
    assert len(syn1.choices) == 4 and syn1.skipped_components == 0
    # chr1 has DSBs at 100, 300, 700 and chr2 one at 500: 8 ends. An end on
    # chr1 has 4 other chr1 ends and 2 chr2 ends as candidates.
    first = syn1.choices[0]
    assert first.chosen_inter and first.chosen_s == 0.0
    assert (len(first.intra), first.num_inter) == (4, 2)
    # Distances below 1 kb are floored: s = -log(1 kb / 1 Mb) = log(1000).
    assert max(first.intra) == pytest.approx(math.log(1000))


def test_refused_components_are_counted():
    ds = load_dataset(SYNTHETIC)
    # SYN-3 joins one breakpoint end to two others: conversion refuses it.
    assert patient_choices(ds.patients["SYN-3"]).skipped_components == 1


# -- the fit -------------------------------------------------------------------------


def simulate(alpha, beta, n, seed=0):
    """Choices drawn from the model itself."""
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        intra = tuple(rng.uniform(-4, 7) for _ in range(rng.randint(3, 12)))
        num_inter = rng.randint(0, 10)
        weights = [math.exp(alpha * s) for s in intra] + [math.exp(beta)] * num_inter
        k = rng.choices(range(len(weights)), weights)[0]
        if k < len(intra):
            out.append(Choice(intra[k], False, intra, num_inter))
        else:
            out.append(Choice(0.0, True, intra, num_inter))
    return out


@pytest.mark.parametrize("alpha, beta", [(0.3, -2.5), (0.0, 0.0), (-0.4, 1.0)])
def test_fit_recovers_known_parameters(alpha, beta):
    choices = simulate(alpha, beta, 4_000, seed=11)
    f = fit(choices, fit_inter=True)
    assert f.converged
    assert abs(f.alpha - alpha) < 4 * f.alpha_se
    assert abs(f.inter_log_weight - beta) < 4 * f.inter_log_weight_se
    assert f.log_likelihood >= f.log_likelihood_uniform


def test_fit_is_the_maximum():
    choices = simulate(0.5, -1.0, 1_000, seed=3)
    f = fit(choices, fit_inter=True)
    best = log_likelihood(choices, f.alpha, f.inter_log_weight)
    for da, db in [(0.01, 0), (-0.01, 0), (0, 0.01), (0, -0.01)]:
        assert log_likelihood(choices, f.alpha + da, f.inter_log_weight + db) < best


def test_alpha_only_fixes_the_inter_weight():
    choices = simulate(0.3, -2.0, 1_000, seed=4)
    f = fit(choices, fit_inter=False)
    assert f.inter_log_weight == 0.0 and f.inter_log_weight_se is None
    best = log_likelihood(choices, f.alpha)
    assert log_likelihood(choices, f.alpha + 0.01) < best
    assert log_likelihood(choices, f.alpha - 0.01) < best


def test_separation_is_reported_as_not_converged():
    # Every end picks its nearest candidate: the likelihood rises without bound.
    choices = [Choice(5.0, False, (5.0, 1.0, 0.0), 0) for _ in range(20)]
    f = fit(choices, fit_inter=False)
    # The fit runs off until the curvature underflows or the limit is hit.
    assert not f.converged and 5 < f.alpha <= PARAMETER_LIMIT


def test_fit_needs_inter_candidates_for_the_inter_weight():
    with pytest.raises(ValueError, match="another chromosome"):
        fit([Choice(1.0, False, (1.0, 2.0), 0)], fit_inter=True)
    with pytest.raises(ValueError, match="No choices"):
        fit([])


def test_target_kwargs_build_the_fitted_target():
    from amg_sampling.samplers.targets import ProximityTarget

    f = fit(simulate(0.3, -2.0, 500, seed=5), fit_inter=True)
    t = ProximityTarget((1, 1, 2, 2), (0, 1, 0, 1), **f.target_kwargs())
    assert t.log_edge_weight(0, 2) == f.inter_log_weight
    assert t.alpha == f.alpha
