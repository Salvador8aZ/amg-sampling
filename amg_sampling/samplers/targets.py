"""Non-uniform target distributions for the Markov chain (``samplers/mcmc.py``).

A target gives every rejoin matching ``r`` a weight ``w(r) = ∏ w(e)`` over its
rejoin edges ``e``, so a 2-switch changes the weight by the ratio of two new
edge weights to two old ones. The chain's acceptance probability becomes
``min(1, w(s)/w(r))`` (``docs/theory/mcmc.md`` §5). Edges fixed by observed
rejoins appear in every state and cancel from that ratio.

The uniform target (every admissible state equally likely) is ``None``
throughout, so the uniform chain is exactly the one verified in §3.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from amg_sampling.data.sheth.conversion import PatientConfiguration

MB = 1_000_000


@dataclass(frozen=True)
class ProximityTarget:
    """Genomic proximity: ends close together on a chromosome rejoin more often.

    An edge joining ends at distance ``d`` (base pairs) on one chromosome has
    weight ``(d / scale_bp) ** (-alpha)``, with ``d`` floored at
    ``min_distance_bp``. An edge between two chromosomes has weight 1, the
    weight of an intra-chromosomal edge at distance ``scale_bp``. ``alpha = 0``
    is the uniform target.
    """

    chromosomes: tuple[int, ...]  # per end
    coordinates: tuple[int, ...]  # per end, base pairs
    alpha: float
    scale_bp: int = MB
    min_distance_bp: int = 1_000

    name = "proximity"

    def __post_init__(self):
        if len(self.chromosomes) != len(self.coordinates):
            raise ValueError("chromosomes and coordinates must have equal length.")
        if not math.isfinite(self.alpha) or self.alpha < 0:
            raise ValueError(f"alpha must be finite and >= 0, got {self.alpha}.")
        if self.scale_bp <= 0 or self.min_distance_bp <= 0:
            raise ValueError("scale_bp and min_distance_bp must be positive.")

    @classmethod
    def from_patient(cls, pc: PatientConfiguration, alpha: float, **kwargs):
        return cls(
            tuple(site.chromosome for site in pc.ends),
            tuple(site.coordinate for site in pc.ends),
            float(alpha),
            **kwargs,
        )

    def log_edge_weight(self, u: int, v: int) -> float:
        if self.chromosomes[u] != self.chromosomes[v]:
            return 0.0
        d = max(abs(self.coordinates[u] - self.coordinates[v]), self.min_distance_bp)
        return -self.alpha * math.log(d / self.scale_bp)

    def log_weight(self, partners: Sequence[int]) -> float:
        """``log w(r)`` of a full partner map."""
        return sum(self.log_edge_weight(u, v) for u, v in enumerate(partners) if u < v)

    def describe(self) -> str:
        return f"proximity, alpha = {self.alpha:g}"
