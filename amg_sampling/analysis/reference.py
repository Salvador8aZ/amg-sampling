"""Exact references for an analysis, computed only when they are feasible.

Every entry records either a value or the reason it is unavailable, so the
report never silently omits an exact result.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import CycleStructure
from amg_sampling.core.matching import PartialMatching
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.enumerate import exact_cycle_distribution
from amg_sampling.exact.formulas import (
    cycle_distribution_all,
    cycle_distribution_deranged,
    cycle_distribution_proper,
    num_all_states,
    num_completions,
    num_deranged_states,
    num_proper_states,
    num_single_cycle_states,
)
from amg_sampling.exact.summaries import (
    cycle_count_distribution,
    largest_cycle_distribution,
)

# Feasibility limits (see docs/theory/exact.md §5 and README "Limitations").
MAX_N_FOR_PARTITIONS = 40  # theorem 2 sums over p(n) partitions; p(40) = 37 338
MAX_K_FOR_PROPER_COUNT = 12  # O(3^k) recursion
MAX_K_FOR_PROPER_DISTRIBUTION = 6
MAX_N_FOR_PROPER_DISTRIBUTION = 30
MAX_COMPLETIONS_TO_ENUMERATE = 2_100_000  # (f−1)!! for f = 16 free ends is 2 027 025


@dataclass
class ExactReference:
    """Exact quantities for uniform ``space(Θ)`` (the combinatorial null)."""

    total_states: int | None = None
    cycle_distribution: dict[CycleStructure, int] | None = None
    cycle_counts: dict[int, int] | None = None
    largest_cycle: dict[int, int] | None = None
    # When summaries are only available for DERANGED, they bound PROPER within
    # `tv_bound`.
    summaries_space: str | None = None
    tv_bound: Fraction | None = None
    unavailable: dict[str, str] = field(default_factory=dict)

    def probability(self, counts: dict | None, key) -> Fraction | None:
        if counts is None:
            return None
        total = sum(counts.values())
        return Fraction(counts.get(key, 0), total) if total else None


def exact_reference(theta: InitialConfiguration, space: StateSpace) -> ExactReference:
    ref = ExactReference()
    n, k = theta.num_dsbs, theta.num_chromosomes
    all_ones = all(b == 1 for b in theta.breaks)

    # |space(Θ)|
    if space is StateSpace.ALL:
        ref.total_states = num_all_states(n)
    elif space is StateSpace.DERANGED:
        ref.total_states = num_deranged_states(n)
    elif all_ones:
        ref.total_states = num_single_cycle_states(n) if n >= 2 else 0
    elif k <= MAX_K_FOR_PROPER_COUNT:
        ref.total_states = num_proper_states(theta)
    else:
        ref.unavailable["total_states"] = (
            f"O(3^k) recursion with k = {k} > {MAX_K_FOR_PROPER_COUNT}"
        )

    # Full cycle-structure distribution.
    layout_free = space is not StateSpace.PROPER or k == 1
    if layout_free:
        if n <= MAX_N_FOR_PARTITIONS:
            deranged = space is not StateSpace.ALL
            ref.cycle_distribution = (
                cycle_distribution_deranged(n)
                if deranged
                else cycle_distribution_all(n)
            )
        else:
            ref.unavailable["cycle_distribution"] = (
                f"theorem 2 needs all partitions of n = {n} > {MAX_N_FOR_PARTITIONS}"
            )
    elif all_ones and n >= 2:
        ref.cycle_distribution = {CycleStructure([n]): num_single_cycle_states(n)}
    elif k <= MAX_K_FOR_PROPER_DISTRIBUTION and n <= MAX_N_FOR_PROPER_DISTRIBUTION:
        ref.cycle_distribution = cycle_distribution_proper(theta)
    else:
        ref.unavailable["cycle_distribution"] = (
            f"PROPER recursion limited to k <= {MAX_K_FOR_PROPER_DISTRIBUTION} and "
            f"n <= {MAX_N_FOR_PROPER_DISTRIBUTION} (here k = {k}, n = {n})"
        )
    if ref.cycle_distribution is not None:
        ref.cycle_distribution = {c: m for c, m in ref.cycle_distribution.items() if m}

    # Scalar summaries.
    if ref.cycle_distribution is not None:
        ref.cycle_counts = _aggregate(ref.cycle_distribution, lambda c: c.num_cycles)
        ref.largest_cycle = _aggregate(ref.cycle_distribution, lambda c: c.parts[0])
        ref.summaries_space = space.value
    elif layout_free:
        summary_space = (
            StateSpace.ALL if space is StateSpace.ALL else StateSpace.DERANGED
        )
        ref.cycle_counts = cycle_count_distribution(n, summary_space)
        ref.largest_cycle = largest_cycle_distribution(n, summary_space)
        ref.summaries_space = space.value
    elif ref.total_states is not None:
        # Uniform(PROPER) and Uniform(DERANGED) differ by total variation
        # 1 − |PROPER|/|DERANGED| (docs/theory/sampling.md §10).
        ref.cycle_counts = cycle_count_distribution(n, StateSpace.DERANGED)
        ref.largest_cycle = largest_cycle_distribution(n, StateSpace.DERANGED)
        ref.summaries_space = StateSpace.DERANGED.value
        ref.tv_bound = 1 - Fraction(ref.total_states, num_deranged_states(n))
    else:
        ref.unavailable["summaries"] = "no exact PROPER summary for this Θ"
    return ref


@dataclass
class CompletionReference:
    """Exact completions of observed rejoins within a state space, if enumerable."""

    num_completions_all: int
    cycle_distribution: dict[CycleStructure, int] | None
    unavailable: str | None


def completion_reference(
    theta: InitialConfiguration, space: StateSpace, fixed: PartialMatching
) -> CompletionReference:
    total = num_completions(len(fixed.free_ends))
    if total > MAX_COMPLETIONS_TO_ENUMERATE:
        return CompletionReference(
            total,
            None,
            f"{total:,} completions exceed the enumeration limit "
            f"{MAX_COMPLETIONS_TO_ENUMERATE:,}",
        )
    return CompletionReference(
        total, exact_cycle_distribution(theta, space, fixed), None
    )


def _aggregate(dist: dict[CycleStructure, int], key) -> dict[int, int]:
    out: dict[int, int] = {}
    for c, m in dist.items():
        out[key(c)] = out.get(key(c), 0) + m
    return dict(sorted(out.items()))
