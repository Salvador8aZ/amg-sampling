"""The 2-switch (reversal) move on rejoin matchings.

A 2-switch takes two rejoin edges ``{u, r(u)}`` and ``{v, r(v)}`` and reconnects
their four ends the other way: ``{u, v}`` and ``{r(u), r(v)}``. The two edges
can be reconnected in two ways (``{a, c}, {b, d}`` or ``{a, d}, {b, c}``); the
choice of which end is ``u`` and which is ``v`` selects one of them.

Biologically, re-pairing the ends of two junctions is the combinatorial
analogue of an inversion or a reciprocal exchange: when both junctions lie on
one chromosome, the segment between them is reversed. Whether this is exactly
the reversal move of Sheth et al. is not asserted here; it is the move used by
the Markov chain in :mod:`amg_sampling.samplers.mcmc`.

Moves can be restricted to a set of *movable* ends. With observed (fixed)
rejoins, the movable ends are the free ends of the
:class:`~amg_sampling.core.matching.PartialMatching`: in a completion they are
matched among themselves, so a switch never touches a fixed edge.

Properties (``docs/theory/mcmc.md``):

* A matching with ``m`` movable edges has exactly ``m(m−1)`` distinct switch
  neighbours, and each arises from exactly four ordered pairs ``(u, v)``.
  (PROVED)
* Choosing ``u`` uniformly among the movable ends and ``v`` uniformly among the
  movable ends other than ``u`` and ``r(u)`` therefore proposes each neighbour
  with probability ``1 / (m(m−1))``: the proposal is symmetric. (PROVED)
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence

from amg_sampling.core.matching import (
    InvalidMatchingError,
    PartialMatching,
    RejoinMatching,
)


def movable_ends(
    matching: RejoinMatching, fixed: PartialMatching | None = None
) -> tuple[int, ...]:
    """Ends a switch may touch: all ends, or the free ends of ``fixed``.

    ``matching`` must complete ``fixed``, so the free ends are matched among
    themselves and every switch keeps the fixed edges.
    """
    if fixed is None:
        return tuple(range(matching.num_ends))
    if not fixed.is_completed_by(matching):
        raise InvalidMatchingError("The matching does not contain every fixed edge.")
    return fixed.free_ends


def switch(matching: RejoinMatching, u: int, v: int) -> RejoinMatching:
    """Replace ``{u, r(u)}`` and ``{v, r(v)}`` by ``{u, v}`` and ``{r(u), r(v)}``.

    ``v`` must differ from ``u`` and from its partner ``r(u)``; otherwise the two
    edges are the same edge and there is nothing to switch.
    """
    partners = list(matching.partners)
    for x in (u, v):
        if not 0 <= x < len(partners):
            raise IndexError(f"End {x} is out of range for {len(partners)} ends.")
    if v in (u, partners[u]):
        raise ValueError(f"Ends {u} and {v} are on the same rejoin edge.")
    switch_in_place(partners, u, v)
    return RejoinMatching(partners)


def switch_in_place(partners: list[int], u: int, v: int) -> None:
    """:func:`switch` on a mutable partner list, without validation.

    For inner loops: the caller guarantees that ``partners`` is a perfect
    matching and that ``v`` is neither ``u`` nor ``partners[u]``.
    """
    a, b = partners[u], partners[v]
    partners[u], partners[v] = v, u
    partners[a], partners[b] = b, a


def switch_neighbours(
    matching: RejoinMatching, movable: Sequence[int] | None = None
) -> Iterator[RejoinMatching]:
    """Every distinct switch neighbour of ``matching``, each exactly once.

    ``movable`` defaults to every end. It must be a union of rejoin edges of
    ``matching`` (as the free ends of a partial matching it completes are).
    """
    r = matching.partners
    ends = tuple(range(len(r))) if movable is None else tuple(movable)
    _require_union_of_edges(r, ends)
    edges = sorted({(min(x, r[x]), max(x, r[x])) for x in ends})
    for i, (a, b) in enumerate(edges):
        for c, d in edges[i + 1 :]:
            yield switch(matching, a, c)  # {a, c}, {b, d}
            yield switch(matching, a, d)  # {a, d}, {b, c}


def num_switch_neighbours(num_movable_edges: int) -> int:
    """``m(m−1)``: the number of switch neighbours of a matching with ``m``
    movable edges.
    """
    if num_movable_edges < 0:
        raise ValueError("The number of edges must be non-negative.")
    return num_movable_edges * (num_movable_edges - 1)


def _require_union_of_edges(r: Sequence[int], ends: Sequence[int]) -> None:
    ends_set = set(ends)
    if len(ends_set) != len(ends):
        raise ValueError("Movable ends must be distinct.")
    for x in ends:
        if not 0 <= x < len(r):
            raise IndexError(f"End {x} is out of range for {len(r)} ends.")
        if r[x] not in ends_set:
            raise ValueError(
                f"End {x} is movable but its partner {r[x]} is not; movable ends "
                "must be a union of rejoin edges."
            )
