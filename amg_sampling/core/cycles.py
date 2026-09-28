"""Cycle structures of the exchange multigraph.

The exchange multigraph Ξ has the free DSB ends as vertices and two kinds of
edges: initial (DSB) edges and final (rejoin) edges (paper, section 2). Both
edge sets are perfect matchings, so every vertex has degree 2 and Ξ is a
disjoint union of cycles alternating between DSB and rejoin edges. A cycle
through ``l`` DSB edges has ``2l`` edges and is written ``C_l``.

``C_1`` is a DSB whose two ends are rejoined to each other (a correctly
repaired DSB): the DSB edge and the rejoin edge are parallel. Chromatin edges
are not part of Ξ, so a ring (a rejoin parallel to a chromatin edge) does not
affect the cycle structure.

The computation walks the two partner maps directly; parallel edges are never
merged because the DSB and rejoin relations are stored separately.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import (
    RejoinMatching,
    require_compatible,
    validate_partner_map,
)


@dataclass(frozen=True, slots=True)
class CycleStructure:
    """A cycle structure ``m_1 C_{l_1} + ... + m_c C_{l_c}``, i.e. a partition of ``n``.

    ``parts`` lists the cycle lengths ``l`` (number of DSBs per cycle) in
    non-increasing order, one entry per cycle; this is the row-length list of
    the Young diagram (paper, section 5).
    """

    parts: tuple[int, ...]

    def __init__(self, parts: Iterable[int]):
        parts = tuple(parts)
        for p in parts:
            if isinstance(p, bool) or not isinstance(p, int) or p < 1:
                raise ValueError(f"Cycle lengths must be positive integers, got {p!r}.")
        object.__setattr__(self, "parts", tuple(sorted(parts, reverse=True)))

    @classmethod
    def from_multiplicities(cls, multiplicities: Mapping[int, int]) -> CycleStructure:
        """Build from ``{l: m_l}``, the number ``m_l`` of cycles ``C_l``."""
        parts: list[int] = []
        for length, count in multiplicities.items():
            if isinstance(count, bool) or not isinstance(count, int) or count < 0:
                raise ValueError(
                    f"Multiplicities must be non-negative integers, got {count!r}."
                )
            parts += [length] * count
        return cls(parts)

    @property
    def n(self) -> int:
        """Number of DSBs, ``sum(l * m_l)``."""
        return sum(self.parts)

    @property
    def num_cycles(self) -> int:
        """``|C|``, the number of cycles."""
        return len(self.parts)

    @property
    def rank(self) -> int:
        """Rank ``n - |C|`` in the graded poset L(n) (paper, theorem 25)."""
        return self.n - len(self.parts)

    def multiplicities(self) -> dict[int, int]:
        """``{l: m_l}`` ordered by decreasing ``l``."""
        return dict(Counter(self.parts))

    def count(self, length: int) -> int:
        """Number of cycles ``C_length``."""
        return self.parts.count(length)

    def __str__(self) -> str:
        return (
            "+".join(
                f"{m}C{l}" if m > 1 else f"C{l}"
                for l, m in self.multiplicities().items()
            )
            or "∅"
        )


def _exchange_cycle_lengths(dsb: Sequence[int], rejoin: Sequence[int]) -> list[int]:
    seen = bytearray(len(dsb))
    lengths = []
    for start in range(len(dsb)):
        if seen[start]:
            continue
        length, v = 0, start
        while True:
            w = dsb[v]
            seen[v] = seen[w] = 1
            length += 1
            v = rejoin[w]
            if v == start:
                break
        lengths.append(length)
    return lengths


def exchange_cycles(
    dsb: Sequence[int], rejoin: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    """Vertex sequences of the cycles of Ξ.

    Each cycle starts at its smallest end ``v`` and follows
    ``v, dsb[v], rejoin[dsb[v]], dsb[rejoin[dsb[v]]], ...``, so consecutive
    vertices alternate DSB and rejoin edges. A ``C_l`` has ``2l`` vertices.
    """
    validate_partner_map(dsb)
    validate_partner_map(rejoin)
    if len(dsb) != len(rejoin):
        raise ValueError("DSB and rejoin matchings must be on the same ends.")
    seen = bytearray(len(dsb))
    cycles = []
    for start in range(len(dsb)):
        if seen[start]:
            continue
        cycle, v = [], start
        while True:
            w = dsb[v]
            seen[v] = seen[w] = 1
            cycle += [v, w]
            v = rejoin[w]
            if v == start:
                break
        cycles.append(tuple(cycle))
    return tuple(cycles)


def cycle_structure_from_partners(
    dsb: Sequence[int], rejoin: Sequence[int]
) -> CycleStructure:
    """Cycle structure of the union of two perfect matchings on the same ends."""
    validate_partner_map(dsb)
    validate_partner_map(rejoin)
    if len(dsb) != len(rejoin):
        raise ValueError("DSB and rejoin matchings must be on the same ends.")
    return CycleStructure(_exchange_cycle_lengths(dsb, rejoin))


def cycle_structure(
    theta: InitialConfiguration, matching: RejoinMatching
) -> CycleStructure:
    """Cycle structure C(Ω) of the AMG with initial configuration ``theta`` and
    rejoins ``matching``.

    Depends only on the DSB and rejoin matchings, not on the chromosome layout.
    """
    require_compatible(theta, matching)
    return CycleStructure(
        _exchange_cycle_lengths(theta.dsb_partners, matching.partners)
    )


def repaired_dsbs(
    theta: InitialConfiguration, matching: RejoinMatching
) -> tuple[int, ...]:
    """DSBs whose two ends are rejoined to each other; each is a ``C_1`` cycle."""
    require_compatible(theta, matching)
    r = matching.partners
    return tuple(j for j in range(theta.num_dsbs) if r[2 * j] == 2 * j + 1)


def ring_fragments(
    theta: InitialConfiguration, matching: RejoinMatching
) -> tuple[int, ...]:
    """Indices into ``theta.fragments`` of fragments closed into a ring.

    Fragment ``(a, b)`` is a ring when its two ends are rejoined to each other,
    so the rejoin edge is parallel to the chromatin edge (paper, figure 3(b)).
    """
    require_compatible(theta, matching)
    r = matching.partners
    return tuple(i for i, (a, b) in enumerate(theta.fragments) if r[a] == b)
