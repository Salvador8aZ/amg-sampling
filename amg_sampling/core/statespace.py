"""State spaces of rejoin matchings for a fixed initial configuration Θ.

* ``ALL``: every perfect matching of the 2n free ends (every complete
  exchange process, paper theorem 1).
* ``DERANGED``: no rejoin edge equals a DSB edge, i.e. no ``C_1`` cycle.
* ``PROPER``: deranged and the whole AMG is connected (paper, section 3).

Rings (a rejoin parallel to a chromatin edge) are allowed in all three.

``PROPER`` is the paper's definition of a proper AMG. It is not necessarily the
state space that patient-level inference will use; that is a separate
modelling decision (see ``docs/theory/representation.md``).

Connectivity is decided on chromosomes rather than on the full AMG: each
chromosome is a path through its telomeres and free ends using chromatin and
DSB edges, so the AMG is connected iff the graph whose vertices are the
chromosomes, with an edge for every rejoin between two different chromosomes,
is connected.
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import Enum

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import RejoinMatching, require_compatible


class StateSpace(Enum):
    ALL = "all"
    DERANGED = "deranged"
    PROPER = "proper"

    def contains(self, theta: InitialConfiguration, matching: RejoinMatching) -> bool:
        """Whether ``matching`` is a state of this space for ``theta``."""
        require_compatible(theta, matching)
        if self is StateSpace.ALL:
            return True
        if self is StateSpace.DERANGED:
            return _is_deranged(matching.partners)
        return _is_deranged(matching.partners) and _is_connected(
            theta, matching.partners
        )

    def contains_partners(
        self, theta: InitialConfiguration, partners: Sequence[int]
    ) -> bool:
        """:meth:`contains` for a raw partner map, without validating it.

        For inner loops such as a Markov chain: the caller guarantees that
        ``partners`` is a perfect matching on the free ends of ``theta``.
        """
        if self is StateSpace.ALL:
            return True
        if not _is_deranged(partners):
            return False
        return self is StateSpace.DERANGED or _is_connected(theta, partners)


def is_deranged(theta: InitialConfiguration, matching: RejoinMatching) -> bool:
    """No DSB is rejoined to its own partner end (no ``C_1`` cycle)."""
    require_compatible(theta, matching)
    return _is_deranged(matching.partners)


def is_connected(theta: InitialConfiguration, matching: RejoinMatching) -> bool:
    """The AMG (chromatin + DSB + rejoin edges, with telomeres) is connected."""
    require_compatible(theta, matching)
    return _is_connected(theta, matching.partners)


def is_proper(theta: InitialConfiguration, matching: RejoinMatching) -> bool:
    """Proper AMG in the sense of the paper: deranged and connected."""
    return StateSpace.PROPER.contains(theta, matching)


def chromosome_components(
    theta: InitialConfiguration, matching: RejoinMatching
) -> tuple[frozenset[int], ...]:
    """Connected components of the chromosome graph, ordered by smallest chromosome."""
    require_compatible(theta, matching)
    parent = _chromosome_union_find(theta, matching.partners)
    groups: dict[int, set[int]] = {}
    for i in range(theta.num_chromosomes):
        groups.setdefault(_find(parent, i), set()).add(i)
    return tuple(sorted((frozenset(g) for g in groups.values()), key=min))


def _is_deranged(r: Sequence[int]) -> bool:
    return all(r[v] != v ^ 1 for v in range(0, len(r), 2))


def _is_connected(theta: InitialConfiguration, r: Sequence[int]) -> bool:
    parent = _chromosome_union_find(theta, r)
    return len({_find(parent, i) for i in range(theta.num_chromosomes)}) == 1


def _chromosome_union_find(theta: InitialConfiguration, r: Sequence[int]) -> list[int]:
    parent = list(range(theta.num_chromosomes))
    chrom = theta.dsb_chromosome
    for v, w in enumerate(r):
        if v < w:
            a, b = _find(parent, chrom[v >> 1]), _find(parent, chrom[w >> 1])
            if a != b:
                parent[a] = b
    return parent


def _find(parent: list[int], i: int) -> int:
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i
