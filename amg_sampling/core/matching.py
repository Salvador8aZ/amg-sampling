"""Rejoin matchings: the state of an aberration multigraph.

A complete exchange process pairs every free DSB end with exactly one other
free end (paper, section 2: the rejoin map ρ : V_F → V_F is a bijection with no
fixed points, and ρ∘ρ = id since rejoin edges are undirected). We store ρ as a
partner map ``r`` with ``r[r[v]] == v`` and ``r[v] != v``.

A rejoin matching is defined on the free ends ``0, ..., 2n-1`` only; it does
not know about chromosomes. Whether it is admissible for a particular Θ is
decided in :mod:`amg_sampling.core.statespace`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, Iterator, Sequence

if TYPE_CHECKING:
    from amg_sampling.core.configuration import InitialConfiguration


class InvalidMatchingError(ValueError):
    """Raised when a partner map is not a perfect matching."""


def validate_partner_map(partners: Sequence[int]) -> None:
    """Check that ``partners`` is a fixed-point-free involution on ``range(len(partners))``."""
    m = len(partners)
    if m % 2:
        raise InvalidMatchingError(f"A perfect matching needs an even number of ends, got {m}.")
    for v, w in enumerate(partners):
        if isinstance(w, bool) or not isinstance(w, int):
            raise InvalidMatchingError(f"Partner of end {v} is not an integer: {w!r}.")
        if not 0 <= w < m:
            raise InvalidMatchingError(f"Partner of end {v} is out of range: {w}.")
        if w == v:
            raise InvalidMatchingError(f"End {v} is matched to itself.")
        if partners[w] != v:
            raise InvalidMatchingError(
                f"Not an involution: {v} -> {w} but {w} -> {partners[w]}."
            )


class RejoinMatching:
    """Immutable perfect matching of free DSB ends, stored as a partner map."""

    __slots__ = ("_partners",)

    def __init__(self, partners: Iterable[int]):
        partners = tuple(partners)
        validate_partner_map(partners)
        self._partners = partners

    @classmethod
    def from_pairs(cls, num_ends: int, pairs: Iterable[tuple[int, int]]) -> RejoinMatching:
        """Build from rejoin edges. Every end in ``range(num_ends)`` must occur exactly once."""
        partners: list[int | None] = [None] * num_ends
        for u, v in pairs:
            for x in (u, v):
                if not 0 <= x < num_ends:
                    raise InvalidMatchingError(f"End {x} is out of range for {num_ends} ends.")
                if partners[x] is not None:
                    raise InvalidMatchingError(f"End {x} occurs in more than one rejoin edge.")
            if u == v:
                raise InvalidMatchingError(f"End {u} is matched to itself.")
            partners[u], partners[v] = v, u
        missing = [v for v, w in enumerate(partners) if w is None]
        if missing:
            raise InvalidMatchingError(f"Ends {missing} are not matched.")
        return cls(partners)  # type: ignore[arg-type]

    @property
    def partners(self) -> tuple[int, ...]:
        return self._partners

    @property
    def num_ends(self) -> int:
        return len(self._partners)

    @property
    def num_edges(self) -> int:
        return len(self._partners) // 2

    def partner(self, end: int) -> int:
        if not 0 <= end < len(self._partners):
            raise IndexError(f"Free end {end} out of range.")
        return self._partners[end]

    def pairs(self) -> tuple[tuple[int, int], ...]:
        """Rejoin edges ``(u, v)`` with ``u < v``, sorted by ``u``."""
        return tuple((v, w) for v, w in enumerate(self._partners) if v < w)

    def relabel(self, permutation: Sequence[int]) -> RejoinMatching:
        """The matching ``π r π⁻¹``: end ``v`` is renamed ``permutation[v]``."""
        if sorted(permutation) != list(range(len(self._partners))):
            raise ValueError("Not a permutation of the free ends.")
        new = [0] * len(self._partners)
        for v, w in enumerate(self._partners):
            new[permutation[v]] = permutation[w]
        return RejoinMatching(new)

    def __len__(self) -> int:
        return len(self._partners)

    def __getitem__(self, end: int) -> int:
        return self._partners[end]

    def __iter__(self) -> Iterator[int]:
        return iter(self._partners)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, RejoinMatching):
            return NotImplemented
        return self._partners == other._partners

    def __hash__(self) -> int:
        return hash(self._partners)

    def __repr__(self) -> str:
        return f"RejoinMatching.from_pairs({len(self._partners)}, {list(self.pairs())})"


def require_compatible(theta: InitialConfiguration, matching: RejoinMatching) -> None:
    """Raise ``ValueError`` unless ``matching`` is defined on the free ends of ``theta``."""
    if matching.num_ends != theta.num_ends:
        raise ValueError(
            f"Matching has {matching.num_ends} ends but Θ{theta.breaks} has {theta.num_ends}."
        )


class PartialMatching:
    """Immutable set of fixed rejoin edges covering some of the free ends.

    Used for observed rejoins: a *completion* is a perfect matching that
    contains every fixed edge. The remaining ``free_ends`` are matched among
    themselves.
    """

    __slots__ = ("_num_ends", "_pairs", "_partner")

    def __init__(self, num_ends: int, pairs: Iterable[tuple[int, int]] = ()):
        if isinstance(num_ends, bool) or not isinstance(num_ends, int) or num_ends < 0 or num_ends % 2:
            raise InvalidMatchingError(f"num_ends must be a non-negative even integer, got {num_ends!r}.")
        pairs = tuple(pairs)  # may be a one-shot iterator; it is traversed twice below
        partner: dict[int, int] = {}
        for u, v in pairs:
            for x in (u, v):
                if isinstance(x, bool) or not isinstance(x, int) or not 0 <= x < num_ends:
                    raise InvalidMatchingError(f"End {x!r} is out of range for {num_ends} ends.")
                if x in partner:
                    raise InvalidMatchingError(f"End {x} occurs in more than one fixed edge.")
            if u == v:
                raise InvalidMatchingError(f"End {u} is matched to itself.")
            partner[u], partner[v] = v, u
        self._num_ends = num_ends
        self._partner = partner
        self._pairs = tuple(sorted((min(u, v), max(u, v)) for u, v in pairs))

    @property
    def num_ends(self) -> int:
        return self._num_ends

    def pairs(self) -> tuple[tuple[int, int], ...]:
        """Fixed edges ``(u, v)`` with ``u < v``, sorted."""
        return self._pairs

    def partner(self, end: int) -> int | None:
        """Fixed partner of ``end``, or ``None`` if ``end`` is free."""
        return self._partner.get(end)

    @property
    def free_ends(self) -> tuple[int, ...]:
        return tuple(v for v in range(self._num_ends) if v not in self._partner)

    def is_completed_by(self, matching: RejoinMatching) -> bool:
        """Whether ``matching`` contains every fixed edge."""
        if matching.num_ends != self._num_ends:
            return False
        return all(matching[u] == v for u, v in self._pairs)

    def __len__(self) -> int:
        return len(self._pairs)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PartialMatching):
            return NotImplemented
        return (self._num_ends, self._pairs) == (other._num_ends, other._pairs)

    def __hash__(self) -> int:
        return hash((self._num_ends, self._pairs))

    def __repr__(self) -> str:
        return f"PartialMatching({self._num_ends}, {list(self._pairs)})"
