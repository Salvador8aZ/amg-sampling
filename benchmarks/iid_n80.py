"""Throughput of IID rejection sampling at n = 80 for several chromosome layouts.

Run with ``uv run python benchmarks/iid_n80.py``. Timings are wall-clock and
depend on the machine; the sample sizes are chosen for stable timings, not
for statistical validation (see amg_sampling/tests/test_iid_n80.py for that).
"""

from __future__ import annotations

import math
import platform
import random
import time

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import cycle_structure
from amg_sampling.core.statespace import StateSpace
from amg_sampling.samplers.rejection import SamplingStats, iter_samples, theoretical_acceptance
from amg_sampling.samplers.uniform import sample_matching

N = 80
LAYOUTS = [
    ("Θ(1,(80))", (80,)),
    ("Θ(2,(40,40))", (40, 40)),
    ("Θ(4,(20,…))", (20,) * 4),
    ("Θ(8,(10,…))", (10,) * 8),
    ("Θ(20,(4,…))", (4,) * 20),
    ("Θ(40,(2,…))", (2,) * 40),
    ("Θ(80,(1,…,1))", (1,) * 80),
]
PROPOSAL_TIMING_DRAWS = 5_000
ACCEPTED_SAMPLES = 2_000


def per_call_seconds(func, calls: int) -> float:
    start = time.perf_counter()
    for _ in range(calls):
        func()
    return (time.perf_counter() - start) / calls


def mean_seconds_over(func, states) -> float:
    """Average time of ``func(state)`` over a pool of random proposals."""
    start = time.perf_counter()
    for state in states:
        func(state)
    return (time.perf_counter() - start) / len(states)


def main() -> None:
    print(f"Python {platform.python_version()} ({platform.machine()}), n = {N}\n")
    rng = random.Random(2026)
    proposal_s = per_call_seconds(lambda: sample_matching(N, rng), PROPOSAL_TIMING_DRAWS)
    print(f"sample_matching: {proposal_s * 1e6:.1f} µs/proposal ({1 / proposal_s:,.0f} proposals/s)\n")

    header = (
        "| layout | space | theoretical acceptance | empirical acceptance (± 1 s.e.) "
        "| proposals/s | accepted/s | cycle_structure µs | PROPER check µs |"
    )
    print(header)
    print("|" + "---|" * 8)
    pool = [sample_matching(N, rng) for _ in range(5_000)]
    for index, (label, breaks) in enumerate(LAYOUTS):
        theta = InitialConfiguration(breaks)
        spaces = [StateSpace.PROPER]
        if breaks == (80,):
            spaces = [StateSpace.ALL, StateSpace.DERANGED, StateSpace.PROPER]
        # Costs are averaged over random proposals: is_proper returns early
        # for the ~40% of proposals that contain a C1.
        cs_s = mean_seconds_over(lambda r: cycle_structure(theta, r), pool)
        pr_s = mean_seconds_over(lambda r: StateSpace.PROPER.contains(theta, r), pool)
        for space in spaces:
            try:
                theory = f"{float(theoretical_acceptance(theta, space)):.6f}"
            except ValueError:
                theory = "n/a (O(3^k) recursion, k > 12)"
            stats = SamplingStats()
            start = time.perf_counter()
            for _ in iter_samples(theta, space, random.Random(100 + index), ACCEPTED_SAMPLES, stats=stats):
                pass
            elapsed = time.perf_counter() - start
            p = stats.acceptance_rate
            # Negative-binomial stopping: s.e. of p̂ ≈ p * sqrt((1-p)/a).
            se = p * math.sqrt((1 - p) / stats.accepted)
            print(
                f"| {label} | {space.value} | {theory} | {p:.4f} ± {se:.4f} "
                f"| {stats.proposals / elapsed:,.0f} | {stats.accepted / elapsed:,.0f} "
                f"| {cs_s * 1e6:.1f} | {pr_s * 1e6:.1f} |"
            )


if __name__ == "__main__":
    main()
