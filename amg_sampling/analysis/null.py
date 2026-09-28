"""Draw IID samples from a combinatorial null model and tally what the analyses need.

The null model is ``Uniform(space(Θ))``, optionally conditioned on containing a
set of fixed rejoins (``Uniform(space(Θ) ∩ completions of F)``). It is a
combinatorial reference: every admissible AMG is equally likely. It is not a
model of biological rejoining.
"""

from __future__ import annotations

import random
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import CycleStructure, cycle_structure
from amg_sampling.core.matching import PartialMatching
from amg_sampling.core.statespace import StateSpace
from amg_sampling.samplers.rejection import SamplingStats, iter_samples


@dataclass
class NullSample:
    """Tallies from ``num_samples`` IID states."""

    num_samples: int
    stats: SamplingStats
    cycle_structures: Counter = field(default_factory=Counter)
    edge_hits: list[int] = field(default_factory=list)
    joint_hits: int = 0


def sample_null(
    theta: InitialConfiguration,
    space: StateSpace,
    rng: random.Random,
    num_samples: int,
    *,
    track_edges: Sequence[tuple[int, int]] = (),
    fixed: PartialMatching | None = None,
    max_proposals_per_sample: int | None = None,
) -> NullSample:
    """Sample and count cycle structures, hits of each tracked edge, and hits of all of
    them.
    """
    if num_samples <= 0:
        raise ValueError("num_samples must be positive.")
    stats = SamplingStats()
    result = NullSample(num_samples, stats, edge_hits=[0] * len(track_edges))
    for r in iter_samples(
        theta,
        space,
        rng,
        num_samples,
        stats=stats,
        fixed=fixed,
        max_proposals_per_sample=max_proposals_per_sample,
    ):
        result.cycle_structures[cycle_structure(theta, r)] += 1
        partners = r.partners
        all_present = True
        for i, (a, b) in enumerate(track_edges):
            if partners[a] == b:
                result.edge_hits[i] += 1
            else:
                all_present = False
        if track_edges and all_present:
            result.joint_hits += 1
    return result


def summary_counts(cycle_structures: Counter, key) -> dict[int, int]:
    out: Counter = Counter()
    for c, m in cycle_structures.items():
        out[key(c)] += m
    return dict(sorted(out.items()))


def num_cycles(c: CycleStructure) -> int:
    return c.num_cycles


def largest_cycle(c: CycleStructure) -> int:
    return c.parts[0]
