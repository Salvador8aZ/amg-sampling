"""Closed-form and recursive exact counts.

All arithmetic is on Python integers (``fractions.Fraction`` where a formula has
rational terms); no floating point. Each function states its source:

* **paper**: a numbered result of Sheth, Arsuaga & Sazdanovic (2026);
* **derived**: derived in this project, with the argument in
  ``docs/theory/exact.md``.

Notation: ``n`` DSBs; ``C = Σ m_l C_l`` a cycle structure (a partition of
``n``); ``|C|`` its number of cycles.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator
from fractions import Fraction
from functools import cache
from math import comb, factorial, prod

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import CycleStructure


def integer_partitions(n: int) -> Iterator[CycleStructure]:
    """Every partition of ``n`` as a :class:`CycleStructure`, in reverse
    lexicographic order.

    For ``n = 0`` yields the empty partition once.
    """
    if n < 0:
        raise ValueError("n must be non-negative.")

    def rec(remaining: int, largest: int) -> Iterator[tuple[int, ...]]:
        if remaining == 0:
            yield ()
            return
        for first in range(min(remaining, largest), 0, -1):
            for rest in rec(remaining - first, first):
                yield (first,) + rest

    for parts in rec(n, n):
        yield CycleStructure(parts)


# -- ALL: every rejoin matching ------------------------------------------------


def num_all_states(n: int) -> int:
    """``|ALL| = κ(n) = (2n−1)!! = ∏_{i=1}^{n} (2i−1)``.

    Paper, theorem 1 (equation 1). Depends only on ``n``.
    """
    _check_n(n)
    return prod(range(1, 2 * n, 2))


def num_states_with_cycle_structure(c: CycleStructure) -> int:
    """Number of rejoin matchings of ``n = c.n`` DSBs whose cycle structure is ``c``.

    ``κ_C(n) = 2^{n−|C|} n! / ∏_l (l^{m_l} m_l!)``, the ``m_l ≥ 0`` form of the
    paper's theorem 2 (equation 2), which allows ``C_1`` parts. Counts states of
    ALL; depends only on ``c``, not on the chromosome layout.
    """
    n = c.n
    denominator = prod(l**m * factorial(m) for l, m in c.multiplicities().items())
    numerator = 2 ** (n - c.num_cycles) * factorial(n)
    count, remainder = divmod(numerator, denominator)
    assert remainder == 0
    return count


def num_completions(num_free_ends: int) -> int:
    """Number of completions of a partial matching with ``f`` free ends: ``(f−1)!!``.

    Derived: completions correspond one-to-one to perfect matchings of the free
    ends (theorem 1 applied to ``f`` points). ALL only; no state-space condition.
    """
    if (
        isinstance(num_free_ends, bool)
        or not isinstance(num_free_ends, int)
        or num_free_ends < 0
        or num_free_ends % 2
    ):
        raise ValueError(
            f"num_free_ends must be a non-negative even integer, got {num_free_ends!r}."
        )
    return prod(range(1, num_free_ends, 2))


def num_single_cycle_states(n: int) -> int:
    """``2^{n−1} (n−1)!`` rejoin matchings with cycle structure ``C_n`` (``n ≥ 1``).

    Paper, corollary 3 (theorem 2 with ``C = C_n``); also the count in lemma 14.
    """
    _check_n(n)
    if n == 0:
        raise ValueError("C_n needs n >= 1.")
    return 2 ** (n - 1) * factorial(n - 1)


def cycle_distribution_all(n: int) -> dict[CycleStructure, int]:
    """Exact distribution of cycle structures over ALL, from theorem 2 (paper).

    Keys in reverse lexicographic order; every partition of ``n`` appears.
    """
    _check_n(n)
    return {c: num_states_with_cycle_structure(c) for c in integer_partitions(n)}


def unsigned_stirling_first_kind(n: int, k: int) -> int:
    """Unsigned Stirling number of the first kind ``[n k]``: permutations of
    ``n`` elements with ``k`` cycles.
    """
    if n < 0 or k < 0:
        raise ValueError("n and k must be non-negative.")
    return _stirling(n, k)


@cache
def _stirling(n: int, k: int) -> int:
    if n == k:
        return 1
    if n == 0 or k == 0:
        return 0
    return _stirling(n - 1, k - 1) + (n - 1) * _stirling(n - 1, k)


def cycle_count_all(n: int, cycles: int) -> int:
    """Number of rejoin matchings in ALL with exactly ``cycles`` exchange
    cycles: ``2^{n−c} [n c]``.

    Derived (``docs/theory/exact.md`` §7.2). A matching corresponds bijectively to a
    permutation of the ``n`` DSBs together with one orientation bit for every
    DSB except one reference DSB per cycle. Consistent with theorems 1 and 2
    of the paper. Section 3.1 of the paper prints ``[n c]`` for this count; see
    the discrepancy note in the theory document.
    """
    _check_n(n)
    if cycles < 0 or cycles > n:
        return 0
    return 2 ** (n - cycles) * unsigned_stirling_first_kind(n, cycles)


# -- DERANGED: no C_1 ------------------------------------------------------------


def num_deranged_states(n: int) -> int:
    """``|DERANGED| = κ'(n)`` by the recurrence ``κ'(n) = 2(n−1)(κ'(n−1) + κ'(n−2))``.

    Paper, theorem 5 (equation 4), stated there for proper AMGs on one
    chromosome, Θ(1,(n)). Initial values ``κ'(0) = 1`` (the empty matching) and
    ``κ'(1) = 0``; the paper gives ``κ'(1) = 0, κ'(2) = 2``, which the recurrence
    reproduces. Derived: DERANGED does not depend on the chromosome layout, so
    this is ``|DERANGED(Θ)|`` for every Θ with ``n`` DSBs.
    """
    _check_n(n)
    a, b = 1, 0  # κ'(0), κ'(1)
    if n == 0:
        return a
    for m in range(2, n + 1):
        a, b = b, 2 * (m - 1) * (b + a)
    return b


def num_deranged_states_closed_form(n: int) -> int:
    """``κ'(n) = n! Σ_{i=0}^{n} (−1)^{n−i} / (n−i)! · C(2i, i) / 2^i``.

    Paper, theorem 5 (equation 5), evaluated exactly with ``Fraction``.
    """
    _check_n(n)
    total = sum(
        Fraction((-1) ** (n - i) * comb(2 * i, i), factorial(n - i) * 2**i)
        for i in range(n + 1)
    )
    value = total * factorial(n)
    assert value.denominator == 1
    return value.numerator


def cycle_distribution_deranged(n: int) -> dict[CycleStructure, int]:
    """Exact distribution over DERANGED: theorem 2 restricted to partitions
    with no part 1.

    Paper, equation 3 (for Θ(1,(n))). Derived: valid for every Θ with ``n`` DSBs.
    """
    _check_n(n)
    return {
        c: num_states_with_cycle_structure(c)
        for c in integer_partitions(n)
        if c.count(1) == 0
    }


def num_two_cycle_deranged_states(n: int) -> int:
    """Deranged rejoin matchings with exactly two cycles:
    ``2^{n−2} (n−1)! Σ_{l=2}^{n−2} 1/l``.

    Paper, proposition 8 (for Θ(1,(n))), evaluated exactly with ``Fraction``.
    """
    _check_n(n)
    if n < 4:
        return 0
    value = (
        2 ** (n - 2) * factorial(n - 1) * sum(Fraction(1, l) for l in range(2, n - 1))
    )
    assert value.denominator == 1
    return value.numerator


# -- PROPER: deranged and connected (depends on Θ) -----------------------------


def lemma15_two_cycle_count(num_dsbs: int) -> int:
    """Proper AMGs with exactly two exchange cycles for Θ(N−1, (2, 1, …, 1)),
    ``N = num_dsbs ≥ 3``.

    ``2^{N−2} (N−2)! (N−3)``, with cycle structures ``C_j + C_{N−j}``, ``2 ≤ j ≤ N−2``.

    Paper, lemma 15 (ii), restated with ``N`` = number of DSBs. The lemma's
    heading uses ``n`` = number of DSBs, but bullet (ii) and its proof use ``n``
    for a configuration with ``n + 1`` DSBs, so the printed expression at ``n``
    equals this function at ``n + 1`` (``docs/theory/exact.md`` §7.1).
    The single-cycle count of lemma 15 (i) is :func:`num_single_cycle_states`.
    """
    _check_n(num_dsbs)
    if num_dsbs < 3:
        raise ValueError("Θ(N−1, (2, 1, ..., 1)) needs N >= 3 DSBs.")
    n = num_dsbs
    return 2 ** (n - 2) * factorial(n - 2) * (n - 3)


def num_proper_states(theta: InitialConfiguration) -> int:
    """``|PROPER(Θ)|`` by recursion on the connected component of the first chromosome.

    Derived (``docs/theory/exact.md``); generalises the paper's lemma 10
    (two chromosomes), lemma 12 (three) and lemma 14 (all ``b_i = 1``).
    Cost is ``O(3^k)`` in the number of chromosomes ``k``, independent of ``n``.
    """
    return sum(cycle_distribution_proper(theta, _scalar=True).values())


def cycle_distribution_proper(
    theta: InitialConfiguration, *, _scalar: bool = False
) -> dict[CycleStructure, int]:
    """Exact distribution of cycle structures over PROPER(Θ), without enumeration.

    Derived (``docs/theory/exact.md``). For a set ``S`` of chromosomes let
    ``D_S`` be the deranged distribution for ``n_S = Σ_{i∈S} b_i`` DSBs and
    ``P_S`` the distribution of deranged matchings on ``S`` whose chromosome
    graph is connected. Splitting by the component ``T`` of the smallest
    chromosome of ``S``:

        D_S = Σ_{T ∋ min S, T ⊆ S} P_T ⊛ D_{S∖T}

    where ``⊛`` combines cycle structures by union of parts. Solving for the
    ``T = S`` term gives ``P_S``. The cost grows with the number of partitions
    of the block sizes, so this is intended for moderate ``n``.
    """
    k = theta.num_chromosomes
    sizes = theta.breaks
    full = (1 << k) - 1

    def block_n(mask: int) -> int:
        return sum(sizes[i] for i in range(k) if mask >> i & 1)

    def deranged(m: int) -> dict[CycleStructure, int]:
        if _scalar:
            return {CycleStructure(()): num_deranged_states(m)}
        return cycle_distribution_deranged(m) if m else {CycleStructure(()): 1}

    def combine(a: dict, b: dict) -> dict:
        out: dict[CycleStructure, int] = defaultdict(int)
        for ca, ma in a.items():
            for cb, mb in b.items():
                key = ca if _scalar else CycleStructure(ca.parts + cb.parts)
                out[key] += ma * mb
        return out

    proper: dict[int, dict[CycleStructure, int]] = {}
    for mask in sorted(range(1, full + 1), key=lambda s: bin(s).count("1")):
        low = mask & -mask
        dist = defaultdict(int, deranged(block_n(mask)))
        sub = (mask - 1) & mask
        while sub:
            if sub & low:
                for key, m in combine(
                    proper[sub], deranged(block_n(mask ^ sub))
                ).items():
                    dist[key] -= m
            sub = (sub - 1) & mask
        proper[mask] = {c: m for c, m in dist.items() if m}
        assert all(m > 0 for m in proper[mask].values())

    result = proper[full]
    return dict(sorted(result.items(), key=lambda item: item[0].parts, reverse=True))


def _check_n(n: int) -> None:
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError(f"n must be a non-negative integer, got {n!r}.")
