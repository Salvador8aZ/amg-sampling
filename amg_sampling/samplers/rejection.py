"""IID uniform sampling from a state space by rejection.

Draw ``r`` uniformly from ALL (:func:`sample_matching`); return it if it lies
in the requested state space, otherwise draw again independently. Accepted
states are IID and uniform on the state space (``docs/theory/sampling.md`` §2).

With ``fixed`` rejoins, proposals are drawn uniformly from the completions of
``fixed`` instead (:func:`sample_completion`), so accepted states are uniform
on ``{states of the space containing every fixed edge}`` (§9).
"""

from __future__ import annotations

import random
from collections.abc import Iterator
from dataclasses import dataclass
from fractions import Fraction

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import PartialMatching, RejoinMatching
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.formulas import (
    num_all_states,
    num_deranged_states,
    num_proper_states,
    num_single_cycle_states,
)
from amg_sampling.samplers.uniform import sample_completion, sample_matching


@dataclass
class SamplingStats:
    """Running counts of proposals drawn from ALL and proposals accepted."""

    proposals: int = 0
    accepted: int = 0

    @property
    def rejected(self) -> int:
        return self.proposals - self.accepted

    @property
    def acceptance_rate(self) -> float | None:
        """Empirical ``accepted / proposals``, or ``None`` before any proposal."""
        return self.accepted / self.proposals if self.proposals else None


def sample_state(
    theta: InitialConfiguration,
    space: StateSpace,
    rng: random.Random,
    *,
    stats: SamplingStats | None = None,
    max_proposals: int | None = None,
    fixed: PartialMatching | None = None,
) -> RejoinMatching:
    """One state drawn uniformly from ``space(Θ)`` by rejection from ALL.

    With ``fixed``, the state is uniform among states of ``space(Θ)`` that
    contain every fixed edge. ``stats``, if given, is updated in place.
    ``max_proposals`` bounds the number of proposals for this call;
    ``RuntimeError`` is raised if it is exhausted.
    """
    _check_nonempty(theta, space)
    if fixed is not None and fixed.num_ends != theta.num_ends:
        raise ValueError(
            f"Fixed rejoins are on {fixed.num_ends} ends but Θ{theta.breaks} "
            f"has {theta.num_ends}."
        )
    n = theta.num_dsbs
    proposals = 0
    while max_proposals is None or proposals < max_proposals:
        r = sample_matching(n, rng) if fixed is None else sample_completion(fixed, rng)
        proposals += 1
        accepted = space.contains(theta, r)
        if stats is not None:
            stats.proposals += 1
            stats.accepted += accepted
        if accepted:
            return r
    raise RuntimeError(
        f"No state of {space.value} accepted in {max_proposals} proposals."
    )


def iter_samples(
    theta: InitialConfiguration,
    space: StateSpace,
    rng: random.Random,
    count: int | None = None,
    *,
    stats: SamplingStats | None = None,
    max_proposals_per_sample: int | None = None,
    fixed: PartialMatching | None = None,
) -> Iterator[RejoinMatching]:
    """Yield ``count`` IID uniform states of ``space(Θ)`` (unbounded if
    ``count`` is None).

    With ``fixed``, states are uniform among those containing every fixed edge.
    """
    _check_nonempty(theta, space)
    produced = 0
    while count is None or produced < count:
        yield sample_state(
            theta,
            space,
            rng,
            stats=stats,
            max_proposals=max_proposals_per_sample,
            fixed=fixed,
        )
        produced += 1


def theoretical_acceptance(
    theta: InitialConfiguration, space: StateSpace, *, max_chromosomes: int = 12
) -> Fraction:
    """Exact acceptance probability ``|space(Θ)| / |ALL|`` of the rejection sampler.

    * ALL: 1.
    * DERANGED: ``κ'(n) / (2n−1)!!`` (theorem 5; independent of the layout).
    * PROPER with every ``b_i = 1``: ``2^{n−1}(n−1)! / (2n−1)!!`` (lemma 14; a
      proof is in ``docs/theory/sampling.md`` §4).
    * PROPER otherwise: :func:`num_proper_states`, whose cost is ``O(3^k)`` in
      the number of chromosomes ``k``. ``ValueError`` is raised when
      ``k > max_chromosomes`` rather than starting an infeasible computation.
    """
    n = theta.num_dsbs
    total = num_all_states(n)
    if space is StateSpace.ALL:
        return Fraction(1)
    if space is StateSpace.DERANGED:
        return Fraction(num_deranged_states(n), total)
    if all(b == 1 for b in theta.breaks):
        return Fraction(num_single_cycle_states(n) if n >= 2 else 0, total)
    if theta.num_chromosomes > max_chromosomes:
        raise ValueError(
            f"|PROPER| for k = {theta.num_chromosomes} chromosomes needs the "
            "O(3^k) recursion; "
            f"raise max_chromosomes (currently {max_chromosomes}) to attempt it."
        )
    return Fraction(num_proper_states(theta), total)


def _check_nonempty(theta: InitialConfiguration, space: StateSpace) -> None:
    # DERANGED and PROPER are empty iff n = 1: for n >= 2 every single-cycle
    # state C_n is proper (docs/theory/exact.md §5, P3).
    if space is not StateSpace.ALL and theta.num_dsbs == 1:
        raise ValueError(
            f"{space.value} is empty for Θ{theta.breaks}; rejection would never stop."
        )
