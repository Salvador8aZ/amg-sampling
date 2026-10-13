"""Exact weighted distributions by enumeration (the oracle for weighted targets).

Under a target with edge weights ``w(e)`` (:mod:`amg_sampling.samplers.targets`),
a state ``r`` of the state space has probability ``w(r) / Z`` with
``w(r) = ∏ w(e)`` over its rejoin edges and ``Z = Σ_r w(r)``. This module
enumerates every state (:func:`amg_sampling.exact.enumerate.iter_states`), so,
like the uniform enumerator, it is for small ``n`` only.

Two different questions have the answer "the most likely cycle structure":

* the **most probable cycle structure**: the structure ``C`` maximising
  ``P(C) = Σ_{r : C(r) = C} w(r) / Z``, which sums over all states with it;
* the structure of the **most likely state** (the MAP state): the cycle
  structure of the single state maximising ``w(r)``.

They can differ: a structure shared by many moderately likely states can
outweigh the structure of the single best one. Both are returned.

Weights are handled in log space (log-sum-exp), so large ``alpha`` or many
edges cannot overflow.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import CycleStructure, cycle_structure
from amg_sampling.core.matching import PartialMatching, RejoinMatching
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.enumerate import iter_states


@dataclass(frozen=True)
class WeightedDistribution:
    """Exact distribution of cycle structures under a weighted target."""

    probabilities: dict[CycleStructure, float]  # sorted by decreasing probability
    num_states: int
    log_partition: float  # log Z
    map_state: RejoinMatching  # the first most likely state in enumeration order
    map_probability: float  # P(map_state)
    num_map_states: int  # states tied for the maximum weight (all of them if uniform)

    @property
    def most_probable_structure(self) -> CycleStructure:
        return next(iter(self.probabilities))

    def map_structure(self, theta: InitialConfiguration) -> CycleStructure:
        return cycle_structure(theta, self.map_state)


def weighted_cycle_distribution(
    theta: InitialConfiguration,
    space: StateSpace,
    target,
    fixed: PartialMatching | None = None,
    tie_tolerance: float = 1e-9,
) -> WeightedDistribution:
    """Exact cycle-structure probabilities under ``target``.

    ``target`` has ``log_weight(partners)``; ``None`` is the uniform target.
    Ties for the MAP state go to the first state in enumeration order; states
    whose log-weights differ by less than ``tie_tolerance`` count as tied.
    """
    log_mass: dict[CycleStructure, list[float]] = {}
    num_states = 0
    best_state, best_log, num_best = None, -math.inf, 0
    for state in iter_states(theta, space, fixed):
        num_states += 1
        lw = 0.0 if target is None else target.log_weight(state.partners)
        log_mass.setdefault(cycle_structure(theta, state), []).append(lw)
        if lw > best_log + tie_tolerance:
            best_state, best_log, num_best = state, lw, 1
        elif lw >= best_log - tie_tolerance:
            num_best += 1
    if best_state is None:
        raise ValueError("Empty distribution: the state space has no states.")
    by_structure = {c: _logsumexp(v) for c, v in log_mass.items()}
    log_z = _logsumexp(list(by_structure.values()))
    probabilities = {c: math.exp(lw - log_z) for c, lw in by_structure.items()}
    ordered = dict(
        sorted(probabilities.items(), key=lambda item: (-item[1], item[0].parts))
    )
    return WeightedDistribution(
        ordered, num_states, log_z, best_state, math.exp(best_log - log_z), num_best
    )


def total_variation(p: dict, q: dict) -> float:
    """Total variation distance ``½ Σ |p(x) − q(x)|`` between two distributions."""
    keys = set(p) | set(q)
    return 0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in keys)


def _logsumexp(values: list[float]) -> float:
    top = max(values)
    if top == -math.inf:
        return top
    return top + math.log(sum(math.exp(v - top) for v in values))
