"""The initial configuration Θ(k, (b_1, ..., b_k)).

Θ is the fixed part of an aberration multigraph: chromosomes, DSBs, free DSB
ends, telomeres, chromatin edges and DSB edges. It does not depend on how the
free ends are rejoined; the rejoin matching (the state) is kept separately.

Labelling convention
--------------------
DSBs are numbered ``j = 0, ..., n-1`` consecutively along chromosome 0, then
chromosome 1, and so on. DSB ``j`` has two free ends::

    2*j      left end  (joined by chromatin towards the start of the chromosome)
    2*j + 1  right end (joined by chromatin towards the end of the chromosome)

so the DSB partner of end ``v`` is ``v ^ 1``. Chromosome ``i`` is the path::

    telomere(i, LEFT) - [2j, 2j+1] - [2j+2, 2j+3] - ... - telomere(i, RIGHT)

where ``-`` is a chromatin edge and ``[x, y]`` is a DSB edge. Chromatin edges
between two free ends ``(2j+1, 2j+2)`` are the internal chromatin edges, called
*fragments*. Telomeres are not free ends and carry no integer end label.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import IntEnum

TELOMERE = -1
"""Sentinel returned by :meth:`InitialConfiguration.chromatin_partner` when the
chromatin neighbour of a free end is a telomere."""


class Side(IntEnum):
    LEFT = 0
    RIGHT = 1


@dataclass(frozen=True, slots=True)
class Telomere:
    """A chromosome end. ``end`` is the free DSB end it is joined to by chromatin."""

    chromosome: int
    side: Side
    end: int


@dataclass(frozen=True, slots=True)
class InitialConfiguration:
    """Θ(k, (b_1, ..., b_k)): ``k`` chromosomes, chromosome ``i`` carrying ``b_i`` DSBs.

    Every chromosome must carry at least one DSB: a chromosome without a DSB
    has no free ends, so it can never be joined to the rest of an AMG.
    """

    breaks: tuple[int, ...]

    num_chromosomes: int = field(init=False, repr=False)
    num_dsbs: int = field(init=False, repr=False)
    num_ends: int = field(init=False, repr=False)
    dsb_chromosome: tuple[int, ...] = field(init=False, repr=False)
    chromatin_partners: tuple[int, ...] = field(init=False, repr=False)
    dsb_partners: tuple[int, ...] = field(init=False, repr=False)
    fragments: tuple[tuple[int, int], ...] = field(init=False, repr=False)
    telomeres: tuple[Telomere, ...] = field(init=False, repr=False)

    def __init__(self, breaks: Sequence[int]):
        breaks = tuple(breaks)
        if not breaks:
            raise ValueError("Θ must contain at least one chromosome.")
        for b in breaks:
            if isinstance(b, bool) or not isinstance(b, int):
                raise TypeError(f"DSB counts must be integers, got {b!r}.")
            if b < 1:
                raise ValueError(
                    f"Every chromosome must carry at least one DSB, got {b}."
                )

        dsb_chromosome = tuple(i for i, b in enumerate(breaks) for _ in range(b))
        n = len(dsb_chromosome)

        chromatin = [TELOMERE] * (2 * n)
        fragments = []
        telomeres = []
        first = 0
        for i, b in enumerate(breaks):
            last = first + b - 1
            telomeres.append(Telomere(i, Side.LEFT, 2 * first))
            telomeres.append(Telomere(i, Side.RIGHT, 2 * last + 1))
            for j in range(first, last):
                a, c = 2 * j + 1, 2 * j + 2
                chromatin[a], chromatin[c] = c, a
                fragments.append((a, c))
            first = last + 1

        set_ = object.__setattr__
        set_(self, "breaks", breaks)
        set_(self, "num_chromosomes", len(breaks))
        set_(self, "num_dsbs", n)
        set_(self, "num_ends", 2 * n)
        set_(self, "dsb_chromosome", dsb_chromosome)
        set_(self, "chromatin_partners", tuple(chromatin))
        set_(self, "dsb_partners", tuple(v ^ 1 for v in range(2 * n)))
        set_(self, "fragments", tuple(fragments))
        set_(self, "telomeres", tuple(telomeres))

    # -- free ends -----------------------------------------------------------

    def dsb_of(self, end: int) -> int:
        """Index of the DSB that created free end ``end``."""
        self._check_end(end)
        return end >> 1

    def dsb_partner(self, end: int) -> int:
        """The other free end of the same DSB (the initial-edge partner)."""
        self._check_end(end)
        return end ^ 1

    def chromatin_partner(self, end: int) -> int:
        """Chromatin neighbour of ``end``, or :data:`TELOMERE`."""
        self._check_end(end)
        return self.chromatin_partners[end]

    def chromosome_of(self, end: int) -> int:
        """Chromosome containing free end ``end``."""
        self._check_end(end)
        return self.dsb_chromosome[end >> 1]

    def ends_of_chromosome(self, chromosome: int) -> range:
        """Free ends of ``chromosome`` in order along the chromosome."""
        if not 0 <= chromosome < self.num_chromosomes:
            raise IndexError(f"No chromosome {chromosome} in {self!r}.")
        first = sum(self.breaks[:chromosome])
        return range(2 * first, 2 * (first + self.breaks[chromosome]))

    # -- edges ---------------------------------------------------------------

    def dsb_edges(self) -> tuple[tuple[int, int], ...]:
        """Initial (DSB) edges ``(2j, 2j+1)``."""
        return tuple((2 * j, 2 * j + 1) for j in range(self.num_dsbs))

    @property
    def num_fragments(self) -> int:
        """Number of internal chromatin edges, ``n - k``."""
        return len(self.fragments)

    def _check_end(self, end: int) -> None:
        if not 0 <= end < self.num_ends:
            raise IndexError(f"Free end {end} out of range for {self!r}.")
