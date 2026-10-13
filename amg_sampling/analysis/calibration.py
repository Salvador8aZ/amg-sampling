"""Fit the proximity target to observed junctions, and test it out of sample.

Model
-----
The proximity target (:class:`~amg_sampling.samplers.targets.ProximityTarget`)
gives a rejoin of ends ``x`` and ``y`` the log-weight

    alpha · s(x, y)     with s = −log(max(d, d_min) / scale), if x and y lie
                        d base pairs apart on one chromosome;
    beta                if they lie on different chromosomes

(``beta`` is ``inter_log_weight``). Two models are fitted: ``alpha`` alone with
``beta = 0`` (the original target), and ``alpha`` and ``beta`` together.

The exact likelihood of a patient's observed rejoins under the target needs
its normalising constant, a sum over all rejoin matchings, which cannot be
computed for realistic components. The fit therefore uses a **composite
(pseudo-)likelihood**: each end of an observed junction is treated as choosing
its partner among the admissible ends of its component (every end except
itself and its own DSB's other end) with probability proportional to the
weight. This is a conditional-logit model. Its log-likelihood is concave in
the parameters (the Hessian is minus a sum of covariance matrices), so the
maximum is unique, and Newton's method with step halving finds it. Standard
errors come from the observed information. Choices are not independent
(observed ends are matched in pairs), so standard errors are a guide, not
exact intervals.

Components that conversion refuses (an end used by two junctions) are skipped
and counted.

Out-of-sample test
------------------
Patients are split 80/20 into training and test sets, by patient, never by
junction, with a fixed seed. The parameters are fitted on the training
patients and evaluated on the test patients by the mean log pseudo-likelihood
per choice, compared with the uniform target (all parameters 0). A positive
difference means the fitted model predicts unseen patients' junctions better
than uniform.
"""

from __future__ import annotations

import math
import random
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from amg_sampling.data.sheth import ConversionError, chromosome_components, convert
from amg_sampling.data.sheth.models import PatientRecord
from amg_sampling.samplers.targets import ProximityTarget

# A fitted parameter at this bound means the likelihood was still rising (e.g.
# every junction joins the nearest candidate): no finite maximum exists.
PARAMETER_LIMIT = 50.0


@dataclass(frozen=True)
class Choice:
    """One end's choice of partner among the candidates of its component.

    Candidates on other chromosomes share the statistic ``s = 0`` and are
    stored as a count.
    """

    chosen_s: float  # s of the chosen partner (0 if it is on another chromosome)
    chosen_inter: bool  # the chosen partner is on another chromosome
    intra: tuple[float, ...]  # s of the candidates on the same chromosome
    num_inter: int  # candidates on other chromosomes


@dataclass(frozen=True)
class PatientChoices:
    patient_id: str
    choices: tuple[Choice, ...]
    skipped_components: int


@dataclass(frozen=True)
class Fit:
    alpha: float
    alpha_se: float
    inter_log_weight: float  # 0 if not fitted
    inter_log_weight_se: float | None  # None if not fitted
    converged: bool  # False if a parameter reached ±PARAMETER_LIMIT
    num_choices: int
    log_likelihood: float  # at the fit
    log_likelihood_uniform: float  # with every parameter 0

    def target_kwargs(self) -> dict:
        """Keyword arguments for :class:`ProximityTarget`."""
        return {"alpha": self.alpha, "inter_log_weight": self.inter_log_weight}


def split_patients(
    patient_ids: Iterable[str],
    test_fraction: float = 0.2,
    seed: int = 0,
    force_test: Sequence[str] = (),
) -> tuple[list[str], list[str]]:
    """Split patients into (train, test), deterministically for a given seed.

    ``round(test_fraction · N)`` patients go to the test set; ``force_test``
    (e.g. a case-study patient) are always among them.
    """
    ids = sorted(set(patient_ids))
    if not 0 < test_fraction < 1:
        raise ValueError(f"test_fraction must be in (0, 1), got {test_fraction}.")
    missing = [p for p in force_test if p not in ids]
    if missing:
        raise ValueError(f"Unknown patients in force_test: {missing}.")
    forced = set(force_test)
    size = max(round(test_fraction * len(ids)), len(forced))
    rest = [p for p in ids if p not in forced]
    random.Random(seed).shuffle(rest)
    test = sorted(forced | set(rest[: size - len(forced)]))
    train = sorted(p for p in ids if p not in test)
    return train, test


def patient_choices(patient: PatientRecord, **target_kwargs) -> PatientChoices:
    """Every observed junction of ``patient`` as two partner choices.

    ``target_kwargs`` (``scale_bp``, ``min_distance_bp``) define ``s`` exactly
    as the target does.
    """
    choices, skipped = [], 0
    for component in chromosome_components(patient):
        try:
            pc = convert(patient, component)
        except ConversionError:
            skipped += 1
            continue
        unit = ProximityTarget.from_patient(pc, alpha=1.0, **target_kwargs)
        chroms = unit.chromosomes
        for edge in pc.observed:
            for x, y in ((edge.end_a, edge.end_b), (edge.end_b, edge.end_a)):
                intra, inter = [], 0
                for z in range(len(chroms)):
                    if z in (x, x ^ 1):
                        continue
                    if chroms[z] == chroms[x]:
                        intra.append(unit.log_edge_weight(x, z))
                    else:
                        inter += 1
                chosen_inter = chroms[y] != chroms[x]
                chosen_s = 0.0 if chosen_inter else unit.log_edge_weight(x, y)
                choices.append(Choice(chosen_s, chosen_inter, tuple(intra), inter))
    return PatientChoices(patient.patient_id, tuple(choices), skipped)


def log_likelihood(
    choices: Sequence[Choice], alpha: float, inter_log_weight: float = 0.0
) -> float:
    """Composite log-likelihood ``Σ log P(chosen)``."""
    return sum(_terms(c, alpha, inter_log_weight)[0] for c in choices)


def fit(
    choices: Sequence[Choice],
    fit_inter: bool = True,
    tol: float = 1e-10,
    max_iter: int = 100,
) -> Fit:
    """Maximum composite-likelihood parameters.

    With ``fit_inter`` the inter-chromosomal log-weight is fitted with
    ``alpha``; otherwise it is fixed at 0. ``alpha`` may come out negative
    (distant ends favoured).
    """
    if not choices:
        raise ValueError("No choices to fit.")
    if fit_inter and not any(c.num_inter for c in choices):
        raise ValueError(
            "No choice has a candidate on another chromosome; the "
            "inter-chromosomal weight cannot be fitted."
        )
    k = 2 if fit_inter else 1
    theta = [0.0] * k

    def evaluate(params):
        a, b = params[0], params[1] if k == 2 else 0.0
        value, grad, hess = 0.0, [0.0, 0.0], [[0.0, 0.0], [0.0, 0.0]]
        for c in choices:
            lp, g, h = _terms(c, a, b)
            value += lp
            for i in range(2):
                grad[i] += g[i]
                for j in range(2):
                    hess[i][j] += h[i][j]
        return value, grad[:k], [row[:k] for row in hess[:k]]

    # Stop on the Newton step, not the gradient: under separation the gradient
    # vanishes exponentially while the step stays large and drives a parameter
    # to the limit.
    value, grad, hess = evaluate(theta)
    for _ in range(max_iter):
        step = _newton_step(grad, hess)
        if max(abs(s) for s in step) < tol:
            break
        scale = 1.0
        for _ in range(60):  # step halving keeps the ascent monotone
            trial = [
                max(-PARAMETER_LIMIT, min(PARAMETER_LIMIT, t + scale * s))
                for t, s in zip(theta, step)
            ]
            trial_value, trial_grad, trial_hess = evaluate(trial)
            if trial_value >= value:
                break
            scale /= 2
        else:
            break
        moved = max(abs(t - u) for t, u in zip(trial, theta))
        theta, value, grad, hess = trial, trial_value, trial_grad, trial_hess
        if moved < 1e-12:
            break
    errors = _standard_errors(hess)
    converged = all(abs(t) < PARAMETER_LIMIT for t in theta) and all(
        math.isfinite(e) for e in errors
    )
    return Fit(
        alpha=theta[0],
        alpha_se=errors[0],
        inter_log_weight=theta[1] if k == 2 else 0.0,
        inter_log_weight_se=errors[1] if k == 2 else None,
        converged=converged,
        num_choices=len(choices),
        log_likelihood=value,
        log_likelihood_uniform=log_likelihood(choices, 0.0, 0.0),
    )


def _terms(c: Choice, a: float, b: float):
    """``log P(chosen)``, its gradient and its Hessian in ``(alpha, beta)``.

    Features: ``(s, 0)`` for an intra-chromosomal candidate, ``(0, 1)`` for an
    inter-chromosomal one.
    """
    exponents = [a * s for s in c.intra]
    inter_exponent = b + math.log(c.num_inter) if c.num_inter else -math.inf
    top = max([*exponents, inter_exponent])
    weights = [math.exp(e - top) for e in exponents]
    inter = math.exp(inter_exponent - top) if c.num_inter else 0.0
    z = sum(weights) + inter
    mean_s = sum(w * s for w, s in zip(weights, c.intra)) / z
    mean_s2 = sum(w * s * s for w, s in zip(weights, c.intra)) / z
    mean_t = inter / z
    log_p = (b if c.chosen_inter else a * c.chosen_s) - top - math.log(z)
    grad = [c.chosen_s - mean_s, float(c.chosen_inter) - mean_t]
    var_s = mean_s2 - mean_s * mean_s
    var_t = mean_t - mean_t * mean_t
    cov = -mean_s * mean_t
    return log_p, grad, [[-var_s, -cov], [-cov, -var_t]]


def _newton_step(grad, hess):
    if len(grad) == 1:
        return [-grad[0] / hess[0][0]] if hess[0][0] else [0.0]
    (p, q), (r, s) = hess
    det = p * s - q * r
    if det == 0:
        return [0.0, 0.0]
    return [-(s * grad[0] - q * grad[1]) / det, -(-r * grad[0] + p * grad[1]) / det]


def _standard_errors(hess) -> list[float]:
    """Square roots of the diagonal of the inverse observed information."""
    if len(hess) == 1:
        info = -hess[0][0]
        return [1 / math.sqrt(info) if info > 0 else math.inf]
    (p, q), (r, s) = hess
    det = p * s - q * r
    if det <= 0:
        return [math.inf, math.inf]
    return [math.sqrt(max(-s / det, 0.0)), math.sqrt(max(-p / det, 0.0))]
