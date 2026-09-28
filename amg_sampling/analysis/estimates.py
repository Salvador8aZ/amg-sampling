"""Monte Carlo estimates of probabilities from IID samples.

A probability ``p`` is estimated by the hit fraction ``k/N`` of ``N`` IID samples.

* Standard error: ``sqrt(p̂(1−p̂)/N)``, the plug-in binomial standard error.
* Interval: the Wilson score interval, which stays inside [0, 1] and is
  informative when ``k`` is 0 or ``N``.
* Zero hits: the estimate is reported as "0 of N", never as probability 0,
  together with the exact one-sided Clopper–Pearson upper bound
  ``1 − (1 − c)^{1/N}`` at confidence ``c`` (about ``3/N`` for c = 0.95).
* Fewer than ``FEW_HITS`` hits: the count and the Wilson interval are reported
  instead of "estimate ± standard error", which is unreliable for rare events.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from statistics import NormalDist

FEW_HITS = 10


@dataclass(frozen=True)
class ProportionEstimate:
    hits: int
    trials: int
    confidence: float = 0.95

    def __post_init__(self):
        if self.trials <= 0 or not 0 <= self.hits <= self.trials:
            raise ValueError(
                f"Need 0 <= hits <= trials and trials > 0, "
                f"got {self.hits}/{self.trials}."
            )
        if not 0 < self.confidence < 1:
            raise ValueError("confidence must be in (0, 1).")

    @property
    def estimate(self) -> float:
        return self.hits / self.trials

    @property
    def standard_error(self) -> float:
        p = self.estimate
        return math.sqrt(p * (1 - p) / self.trials)

    @property
    def wilson_interval(self) -> tuple[float, float]:
        z = NormalDist().inv_cdf(0.5 + self.confidence / 2)
        n, p = self.trials, self.estimate
        centre = (p + z * z / (2 * n)) / (1 + z * z / n)
        half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
        return max(0.0, centre - half), min(1.0, centre + half)

    @property
    def zero_hit_upper_bound(self) -> float | None:
        """Exact one-sided upper confidence bound when no hit was observed."""
        if self.hits:
            return None
        return 1 - (1 - self.confidence) ** (1 / self.trials)

    def describe(self) -> str:
        if self.hits == 0:
            return (
                f"0 occurrences in {self.trials:,} samples "
                f"({self.confidence:.0%} upper bound "
                f"{format_probability(self.zero_hit_upper_bound)})"
            )
        if self.hits < FEW_HITS:
            low, high = self.wilson_interval
            return (
                f"{self.hits} occurrence{'s' if self.hits > 1 else ''} "
                f"in {self.trials:,} samples "
                f"({self.confidence:.0%} interval "
                f"{format_probability(low)}–{format_probability(high)})"
            )
        return (
            f"{format_probability(self.estimate)} ± "
            f"{format_probability(self.standard_error)}"
        )

    def to_dict(self) -> dict:
        low, high = self.wilson_interval
        return {
            "hits": self.hits,
            "samples": self.trials,
            "estimate": self.estimate,
            "standard_error": self.standard_error,
            "confidence": self.confidence,
            "wilson_interval": [low, high],
            "zero_hit_upper_bound": self.zero_hit_upper_bound,
        }


def format_probability(p: float | None) -> str:
    """Four significant digits, scientific notation for small values."""
    if p is None:
        return "n/a"
    if p == 0:
        return "0"
    if abs(p) < 1e-3:
        return f"{p:.3e}"
    return f"{p:.4g}"
