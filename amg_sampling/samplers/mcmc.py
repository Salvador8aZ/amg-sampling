"""Markov chain Monte Carlo on a state space, with 2-switch (reversal) moves.

Each step picks an end ``u`` uniformly among the movable ends and an end ``v``
uniformly among the movable ends other than ``u`` and its partner. It then
proposes the 2-switch that rejoins ``u`` to ``v`` and their old partners to
each other (:mod:`amg_sampling.core.moves`). The proposal is accepted if the
new state lies in the requested state space, and otherwise the chain stays
where it is.

For a uniform target this is the Metropolis–Hastings algorithm. The proposal
is symmetric, so the acceptance probability ``min(1, π(s)/π(r))`` is 1 inside
the state space and 0 outside it. The uniform distribution is therefore
stationary, and the chain converges to it whenever the switch graph restricted
to the state space is connected (``docs/theory/mcmc.md``). A non-uniform
(biologically weighted) target changes only the acceptance probability.

With ``fixed`` (observed) rejoins, only the free ends move, so every state of
the chain contains every fixed edge.

Unlike :mod:`amg_sampling.samplers.rejection`, consecutive states are
correlated. Estimates need standard errors that account for this; see
:func:`amg_sampling.analysis.estimates.batch_means_ess`.

By default the chain starts from one exact uniform draw of the rejection
sampler, so for a uniform target it is stationary from the first step. The
burn-in exists for starting states supplied by the caller and for future
non-uniform targets.
"""

from __future__ import annotations

import random
from collections.abc import Iterator
from fractions import Fraction

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import PartialMatching, RejoinMatching
from amg_sampling.core.moves import movable_ends, switch, switch_in_place
from amg_sampling.core.statespace import StateSpace
from amg_sampling.samplers.rejection import SamplingStats, sample_state


def iter_chain(
    theta: InitialConfiguration,
    space: StateSpace,
    rng: random.Random,
    count: int | None = None,
    *,
    burn_in: int = 0,
    thin: int = 1,
    fixed: PartialMatching | None = None,
    initial: RejoinMatching | None = None,
    stats: SamplingStats | None = None,
    max_proposals_for_start: int | None = None,
) -> Iterator[RejoinMatching]:
    """Yield ``count`` states of a Markov chain on ``space(Θ)`` (unbounded if
    ``count`` is None).

    The chain first takes ``burn_in`` steps, then yields the state after every
    ``thin`` further steps. ``initial`` defaults to one uniform draw from the
    rejection sampler, which may use up to ``max_proposals_for_start``
    proposals. ``stats``, if given, counts the chain's proposed and accepted
    moves (not the rejection sampler's proposals).
    """
    if not isinstance(rng, random.Random):
        raise TypeError("rng must be a random.Random instance.")
    if burn_in < 0:
        raise ValueError(f"burn_in must be non-negative, got {burn_in}.")
    if thin < 1:
        raise ValueError(f"thin must be at least 1, got {thin}.")
    if initial is None:
        initial = sample_state(
            theta, space, rng, fixed=fixed, max_proposals=max_proposals_for_start
        )
    elif not space.contains(theta, initial):
        raise ValueError(f"The initial state is not in {space.value}.")
    movable = movable_ends(initial, fixed)
    partners = list(initial.partners)
    step = _stepper(theta, space, rng, movable, partners, stats)

    for _ in range(burn_in):
        step()
    produced = 0
    while count is None or produced < count:
        for _ in range(thin):
            step()
        yield RejoinMatching(partners)
        produced += 1


def _stepper(theta, space, rng, movable, partners, stats):
    """One Metropolis–Hastings step on ``partners`` (in place), as a closure."""
    f = len(movable)
    if f < 4:
        # At most one movable edge: the chain has a single state and never moves.
        def stay():
            if stats is not None:
                stats.proposals += 1

        return stay

    def step():
        u = movable[rng.randrange(f)]
        a = partners[u]
        v = movable[rng.randrange(f)]
        while v == u or v == a:  # uniform among the f − 2 other ends
            v = movable[rng.randrange(f)]
        b = partners[v]
        switch_in_place(partners, u, v)  # now {u, v} and {a, b}
        accepted = _admissible(theta, space, partners, u, v, a, b)
        if not accepted:
            switch_in_place(partners, u, a)  # restores {u, a} and {v, b}
        if stats is not None:
            stats.proposals += 1
            stats.accepted += accepted

    return step


def _admissible(theta, space, partners, u, v, a, b) -> bool:
    # The current state is in the space, so only the two new edges can create a
    # repaired DSB (a C1 cycle). Connectivity is checked in full.
    if space is StateSpace.ALL:
        return True
    if u ^ 1 == v or a ^ 1 == b:
        return False
    return space is StateSpace.DERANGED or space.contains_partners(theta, partners)


def transition_probabilities(
    theta: InitialConfiguration,
    space: StateSpace,
    state: RejoinMatching,
    fixed: PartialMatching | None = None,
) -> dict[RejoinMatching, Fraction]:
    """Exact one-step transition probabilities of the chain from ``state``.

    This enumerates every ordered pair ``(u, v)`` that the chain can draw,
    each with probability ``1 / (f (f − 2))`` for ``f`` movable ends, and
    applies the same acceptance rule. It is for verifying the chain on small
    state spaces (``docs/theory/mcmc.md`` §3).
    """
    if not space.contains(theta, state):
        raise ValueError(f"The state is not in {space.value}.")
    movable = movable_ends(state, fixed)
    f = len(movable)
    if f < 4:
        return {state: Fraction(1)}
    r = state.partners
    weight = Fraction(1, f * (f - 2))
    out: dict[RejoinMatching, Fraction] = {}
    for u in movable:
        for v in movable:
            if v in (u, r[u]):
                continue
            proposal = switch(state, u, v)
            target = proposal if space.contains(theta, proposal) else state
            out[target] = out.get(target, Fraction(0)) + weight
    return out
