"""Exactly uniform random perfect matchings of the 2n free DSB ends.

Algorithm: shuffle the labels ``0, ..., 2n-1`` uniformly and pair positions
``(0, 1), (2, 3), ..., (2n-2, 2n-1)``.

Why this is uniform (``docs/theory/sampling.md`` §1): a perfect matching ``M``
arises from exactly ``2^n n!`` of the ``(2n)!`` orderings (its ``n`` pairs can be
placed in the ``n`` position slots in ``n!`` orders, and each pair in 2
orientations), so every matching has probability
``2^n n! / (2n)! = 1 / (2n-1)!!``.

The shuffle is :meth:`random.Random.shuffle` (Fisher-Yates with unbiased
``_randbelow``), so the ordering is uniform up to the quality of the
generator. Randomness comes only from the ``rng`` argument; no global random
state is used.
"""

from __future__ import annotations

import random
from typing import Sequence

from amg_sampling.core.matching import RejoinMatching


def matching_from_ordering(ends: Sequence[int]) -> RejoinMatching:
    """Pair consecutive positions of an ordering of the free ends: ``(ends[0], ends[1])``, ...

    ``ends`` must be a permutation of ``range(len(ends))`` with even length.
    """
    if sorted(ends) != list(range(len(ends))):
        raise ValueError("ends must be a permutation of range(len(ends)).")
    if len(ends) % 2:
        raise ValueError("An even number of ends is required.")
    partners = [0] * len(ends)
    for i in range(0, len(ends), 2):
        a, b = ends[i], ends[i + 1]
        partners[a], partners[b] = b, a
    return RejoinMatching(partners)


def sample_matching(num_dsbs: int, rng: random.Random) -> RejoinMatching:
    """A uniformly random perfect matching of the ``2 * num_dsbs`` free ends."""
    if not isinstance(rng, random.Random):
        raise TypeError("rng must be a random.Random instance.")
    if isinstance(num_dsbs, bool) or not isinstance(num_dsbs, int) or num_dsbs < 0:
        raise ValueError(f"num_dsbs must be a non-negative integer, got {num_dsbs!r}.")
    ends = list(range(2 * num_dsbs))
    rng.shuffle(ends)
    return matching_from_ordering(ends)
