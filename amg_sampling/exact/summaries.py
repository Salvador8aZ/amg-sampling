"""Exact distributions of scalar cycle summaries over ALL and DERANGED, for any n.

These avoid summing theorem 2 over all partitions of ``n`` (``p(80)`` is about
1.6·10^7). All are derived from one decomposition (``docs/theory/exact.md`` §8):
in any rejoin matching, the exchange cycle through DSB 0 has some length
``l``; its other ``l−1`` DSBs can be chosen in ``C(n−1, l−1)`` ways, and a
single cycle through ``l`` given DSBs can be formed in ``2^{l−1}(l−1)!`` ways
(corollary 3). The remaining ``n−l`` DSBs carry an arbitrary matching with the
same constraints. Hence, if ``A(m)`` counts matchings of ``m`` DSBs whose
cycle lengths all lie in a set ``L``,

    A(0) = 1,   A(m) = Σ_{l ∈ L, l ≤ m} (m−1)!/(m−l)! · 2^{l−1} · A(m−l).

ALL uses ``L = {1, 2, ...}``; DERANGED uses ``L = {2, 3, ...}``. For Θ(1,(n)),
PROPER = DERANGED. For other Θ these are **not** PROPER distributions.

Counts are integers; divide by ``num_all_states(n)`` or
``num_deranged_states(n)`` for probabilities.
"""

from __future__ import annotations

from math import comb, factorial
from typing import Callable

from amg_sampling.core.statespace import StateSpace


def _min_length(space: StateSpace) -> int:
    if space is StateSpace.ALL:
        return 1
    if space is StateSpace.DERANGED:
        return 2
    raise ValueError("Summary distributions are layout-free: use ALL or DERANGED.")


def _single_cycles(l: int) -> int:
    return 2 ** (l - 1) * factorial(l - 1)


def _count_with_lengths(n: int, allowed: Callable[[int], bool]) -> list[int]:
    """``A(m)`` for ``m = 0..n``: matchings of ``m`` DSBs using only allowed cycle lengths."""
    a = [1] + [0] * n
    for m in range(1, n + 1):
        a[m] = sum(
            factorial(m - 1) // factorial(m - l) * 2 ** (l - 1) * a[m - l]
            for l in range(1, m + 1)
            if allowed(l)
        )
    return a


def cycle_count_distribution(n: int, space: StateSpace) -> dict[int, int]:
    """``{c: number of states with exactly c exchange cycles}`` over ALL or DERANGED."""
    lo = _min_length(space)
    # table[m][c]: matchings of m DSBs with c cycles, all of length >= lo.
    table = [[0] * (n + 1) for _ in range(n + 1)]
    table[0][0] = 1
    for m in range(1, n + 1):
        for l in range(lo, m + 1):
            weight = factorial(m - 1) // factorial(m - l) * 2 ** (l - 1)
            prev = table[m - l]
            row = table[m]
            for c in range(1, m + 1):
                if prev[c - 1]:
                    row[c] += weight * prev[c - 1]
    return {c: k for c, k in enumerate(table[n]) if k}


def largest_cycle_distribution(n: int, space: StateSpace) -> dict[int, int]:
    """``{L: number of states whose longest exchange cycle is C_L}`` over ALL or DERANGED."""
    lo = _min_length(space)
    at_most = {}
    for longest in range(lo - 1, n + 1):
        at_most[longest] = _count_with_lengths(n, lambda l: lo <= l <= longest)[n]
    return {
        longest: at_most[longest] - at_most[longest - 1]
        for longest in range(lo, n + 1)
        if at_most[longest] - at_most[longest - 1]
    }


def cycle_length_multiplicity_distribution(
    n: int, length: int, space: StateSpace
) -> dict[int, int]:
    """``{m: number of states with exactly m cycles C_length}`` over ALL or DERANGED.

    Choose the ``m·length`` DSBs of those cycles, split them into ``m`` unordered
    blocks and form a single cycle on each block; the remaining DSBs carry a
    matching with no cycle of that length.
    """
    lo = _min_length(space)
    if length < lo:
        raise ValueError(f"No C_{length} cycles exist in {space.value}.")
    rest = _count_with_lengths(n, lambda l: l >= lo and l != length)
    out = {}
    for m in range(0, n // length + 1):
        k = m * length
        blocks = factorial(k) // (factorial(length) ** m * factorial(m))
        count = comb(n, k) * blocks * _single_cycles(length) ** m * rest[n - k]
        if count:
            out[m] = count
    return out
