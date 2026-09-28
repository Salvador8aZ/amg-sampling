"""Exact enumeration of the states of a state space (the small-n oracle).

The enumerator lists every rejoin matching of Θ's free ends exactly once and
keeps those that belong to the requested :class:`StateSpace`. It is meant to be
trustworthy ground truth for small ``n``; it is not meant to be fast.

Algorithm
---------
Canonical decomposition of perfect matchings: the smallest unmatched end ``v``
is paired with each larger unmatched end ``w`` in turn, and the procedure
recurses on the remaining ends. Every perfect matching is produced exactly
once, because the partner of the smallest unmatched end determines the branch
(see ``docs/theory/exact.md``).

For ``DERANGED`` and ``PROPER`` the branch ``w == d(v)`` is skipped: any
matching containing a DSB edge is not deranged, so this removes no admissible
state. Membership is then checked explicitly with :meth:`StateSpace.contains`
before a state is yielded, so correctness never depends on the pruning.

With ``fixed`` rejoins, the fixed ends are pre-matched and only the free ends
are paired, so the procedure lists every completion of ``fixed`` exactly once.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator, Mapping
from fractions import Fraction

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import CycleStructure, cycle_structure
from amg_sampling.core.matching import PartialMatching, RejoinMatching
from amg_sampling.core.statespace import StateSpace


def iter_states(
    theta: InitialConfiguration, space: StateSpace, fixed: PartialMatching | None = None
) -> Iterator[RejoinMatching]:
    """Yield every state of ``space`` for ``theta`` exactly once, in a
    deterministic order.

    If ``fixed`` is given, only states containing every fixed rejoin edge are listed.
    """
    num_ends = theta.num_ends
    skip_dsb_edges = space is not StateSpace.ALL
    partners = [-1] * num_ends
    if fixed is not None:
        _check_fixed(theta, fixed)
        for u, v in fixed.pairs():
            partners[u], partners[v] = v, u

    def rec(v: int) -> Iterator[RejoinMatching]:
        while v < num_ends and partners[v] != -1:
            v += 1
        if v == num_ends:
            state = RejoinMatching(partners)
            if space.contains(theta, state):
                yield state
            return
        for w in range(v + 1, num_ends):
            if partners[w] != -1 or (skip_dsb_edges and w == v ^ 1):
                continue
            partners[v], partners[w] = w, v
            yield from rec(v + 1)
            partners[v] = partners[w] = -1

    return rec(0)


def count_states(
    theta: InitialConfiguration, space: StateSpace, fixed: PartialMatching | None = None
) -> int:
    """``|space(Θ)|`` (or the number of its states containing ``fixed``) by
    exhaustive enumeration.
    """
    return sum(1 for _ in iter_states(theta, space, fixed))


def exact_cycle_distribution(
    theta: InitialConfiguration, space: StateSpace, fixed: PartialMatching | None = None
) -> dict[CycleStructure, int]:
    """Number of states of ``space`` with each cycle structure, by exhaustive
    enumeration.

    Keys are in reverse lexicographic order of ``parts`` (``C4, C3+C1, 2C2, ...``).
    Cycle structures with no states are omitted.
    """
    counts = Counter(
        cycle_structure(theta, r) for r in iter_states(theta, space, fixed)
    )
    return dict(sorted(counts.items(), key=lambda item: item[0].parts, reverse=True))


def probabilities(
    counts: Mapping[CycleStructure, int],
) -> dict[CycleStructure, Fraction]:
    """Exact probabilities ``count(C) / total`` from integer counts."""
    total = sum(counts.values())
    if total == 0:
        raise ValueError("Empty distribution: the state space has no states.")
    return {c: Fraction(m, total) for c, m in counts.items()}


def _check_fixed(theta: InitialConfiguration, fixed: PartialMatching) -> None:
    if fixed.num_ends != theta.num_ends:
        raise ValueError(
            f"Fixed rejoins are on {fixed.num_ends} ends but Θ{theta.breaks} "
            f"has {theta.num_ends}."
        )
