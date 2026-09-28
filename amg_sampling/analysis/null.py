"""Sample a combinatorial null model and tally what the analyses need.

Samples are either IID (rejection, ``method="iid"``) or the states of a Markov
chain with 2-switch moves (``method="mcmc"``). Chain samples are correlated, so
their tallies are also kept per batch of consecutive samples. The analyses use
these batches for effective sample sizes and standard errors
(:func:`amg_sampling.analysis.estimates.batch_means_ess`).

The null model is ``Uniform(space(Θ))``, optionally conditioned on containing a
set of fixed rejoins (``Uniform(space(Θ) ∩ completions of F)``). It is a
combinatorial reference: every admissible AMG is equally likely. It is not a
model of biological rejoining.
"""

from __future__ import annotations

import math
import random
from collections import Counter
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

from amg_sampling.analysis.estimates import (
    ProportionEstimate,
    batch_means_ess,
    batch_means_standard_error,
)
from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import CycleStructure, cycle_structure
from amg_sampling.core.matching import PartialMatching
from amg_sampling.core.statespace import StateSpace
from amg_sampling.samplers.mcmc import iter_chain
from amg_sampling.samplers.rejection import SamplingStats, iter_samples

METHODS = ("iid", "mcmc")


@dataclass
class NullSample:
    """Tallies from ``num_samples`` states (IID, or consecutive chain states).

    For chain samples ``batch_size > 0``, and the ``batch_*`` lists hold the same
    tallies for each full batch of ``batch_size`` consecutive samples.
    """

    num_samples: int
    stats: SamplingStats
    cycle_structures: Counter = field(default_factory=Counter)
    edge_hits: list[int] = field(default_factory=list)
    joint_hits: int = 0
    method: str = "iid"
    batch_size: int = 0
    batch_cycle_structures: list[Counter] = field(default_factory=list)
    batch_edge_hits: list[list[int]] = field(default_factory=list)
    batch_joint_hits: list[int] = field(default_factory=list)

    @property
    def correlated(self) -> bool:
        return self.batch_size > 0

    def reference_ess(self) -> float:
        """ESS of the number of cycles; used where an indicator's own ESS
        cannot be estimated (no hits, or every sample a hit).
        """
        if not self.correlated:
            return float(self.num_samples)
        means = [
            sum(c.num_cycles * m for c, m in batch.items()) / self.batch_size
            for batch in self.batch_cycle_structures
        ]
        se = batch_means_standard_error(means)
        n = self.num_samples
        mean = sum(c.num_cycles * m for c, m in self.cycle_structures.items()) / n
        var = sum(
            m * (c.num_cycles - mean) ** 2 for c, m in self.cycle_structures.items()
        ) / max(n - 1, 1)
        if not se or var == 0:
            return float(n)
        # Var(mean) = var / ESS, and se² estimates Var(mean) over full batches.
        full = self.batch_size * len(means)
        return min(float(n), var / (se * se) * n / full)

    def ess(self, hits: int, batch_hits: Sequence[int]) -> float | None:
        """Effective sample size for a hit count (``None`` for IID samples)."""
        if not self.correlated:
            return None
        ess = batch_means_ess(batch_hits, self.batch_size, hits, self.num_samples)
        return self.reference_ess() if ess is None else ess

    def estimate(
        self, hits: int, batch_hits: Sequence[int], confidence: float = 0.95
    ) -> ProportionEstimate:
        return ProportionEstimate(
            hits, self.num_samples, confidence, self.ess(hits, batch_hits)
        )

    def edge_estimate(self, i: int, confidence: float = 0.95) -> ProportionEstimate:
        return self.estimate(
            self.edge_hits[i], [b[i] for b in self.batch_edge_hits], confidence
        )

    def joint_estimate(self, confidence: float = 0.95) -> ProportionEstimate:
        return self.estimate(self.joint_hits, self.batch_joint_hits, confidence)

    def ess_by_key(
        self, key: Callable[[CycleStructure], object] | None = None
    ) -> dict | None:
        """ESS of each cycle structure (or of each ``key(structure)``), or
        ``None`` for IID samples. Keys never sampled get the reference ESS.
        """
        if not self.correlated:
            return None
        key = key or (lambda c: c)
        totals: Counter = Counter()
        batches = []
        for batch in self.batch_cycle_structures:
            grouped: Counter = Counter()
            for c, m in batch.items():
                grouped[key(c)] += m
            batches.append(grouped)
        for c, m in self.cycle_structures.items():
            totals[key(c)] += m
        return {
            k: self.ess(hits, [b.get(k, 0) for b in batches])
            for k, hits in totals.items()
        }

    def mean_standard_error(self, value: Callable[[CycleStructure], int]) -> float:
        """Standard error of the sample mean of ``value(structure)``."""
        n = self.num_samples
        mean = sum(value(c) * m for c, m in self.cycle_structures.items()) / n
        if self.correlated:
            se = batch_means_standard_error(
                [
                    sum(value(c) * m for c, m in batch.items()) / self.batch_size
                    for batch in self.batch_cycle_structures
                ]
            )
            if se is not None:
                return se
        var = sum(
            m * (value(c) - mean) ** 2 for c, m in self.cycle_structures.items()
        ) / max(n - 1, 1)
        return math.sqrt(var / n)


def sample_null(
    theta: InitialConfiguration,
    space: StateSpace,
    rng: random.Random,
    num_samples: int,
    *,
    track_edges: Sequence[tuple[int, int]] = (),
    fixed: PartialMatching | None = None,
    max_proposals_per_sample: int | None = None,
    method: str = "iid",
    burn_in: int = 0,
    thin: int = 1,
) -> NullSample:
    """Sample and count cycle structures, hits of each tracked edge, and hits of all of
    them.

    ``method="mcmc"`` runs one chain (``burn_in`` and ``thin`` as in
    :func:`~amg_sampling.samplers.mcmc.iter_chain`); its starting state is drawn
    by rejection with at most ``max_proposals_per_sample`` proposals.
    """
    if num_samples <= 0:
        raise ValueError("num_samples must be positive.")
    if method not in METHODS:
        raise ValueError(
            f"Unknown sampling method {method!r}; choose one of {METHODS}."
        )
    stats = SamplingStats()
    k = len(track_edges)
    result = NullSample(num_samples, stats, edge_hits=[0] * k, method=method)
    if method == "iid":
        states = iter_samples(
            theta,
            space,
            rng,
            num_samples,
            stats=stats,
            fixed=fixed,
            max_proposals_per_sample=max_proposals_per_sample,
        )
    else:
        states = iter_chain(
            theta,
            space,
            rng,
            num_samples,
            burn_in=burn_in,
            thin=thin,
            fixed=fixed,
            stats=stats,
            max_proposals_for_start=max_proposals_per_sample,
        )
        result.batch_size = max(1, math.isqrt(num_samples))
    size = result.batch_size
    num_full = num_samples // size if size else 0
    for index, r in enumerate(states):
        c = cycle_structure(theta, r)
        result.cycle_structures[c] += 1
        partners = r.partners
        all_present = True
        hits = []
        for i, (a, b) in enumerate(track_edges):
            hit = partners[a] == b
            hits.append(hit)
            if hit:
                result.edge_hits[i] += 1
            else:
                all_present = False
        joint = bool(track_edges) and all_present
        result.joint_hits += joint
        if size and index // size < num_full:
            if index % size == 0:
                result.batch_cycle_structures.append(Counter())
                result.batch_edge_hits.append([0] * k)
                result.batch_joint_hits.append(0)
            result.batch_cycle_structures[-1][c] += 1
            for i, hit in enumerate(hits):
                result.batch_edge_hits[-1][i] += hit
            result.batch_joint_hits[-1] += joint
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
