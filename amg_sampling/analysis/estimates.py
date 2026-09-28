"""Monte Carlo estimates of probabilities from IID or Markov chain samples.

A probability ``p`` is estimated by the hit fraction ``k/N`` of ``N`` samples.

* Standard error: ``sqrt(p̂(1−p̂)/N)``, the plug-in binomial standard error.
* Interval: the Wilson score interval, which stays inside [0, 1] and is
  informative when ``k`` is 0 or ``N``.
* Zero hits: the estimate is reported as "0 of N", never as probability 0,
  together with the exact one-sided Clopper–Pearson upper bound
  ``1 − (1 − c)^{1/N}`` at confidence ``c`` (about ``3/N`` for c = 0.95).
* Fewer than ``FEW_HITS`` hits: the count and the Wilson interval are reported
  instead of "estimate ± standard error", which is unreliable for rare events.

Markov chain samples are correlated, so ``N`` overstates the information they
carry. For them, every formula above uses the *effective sample size*
``N_eff`` in place of ``N``, estimated by non-overlapping batch means
(:func:`batch_means_ess`). The chain is cut into ``B`` batches of ``L``
consecutive samples, and the variance of the batch means estimates the
variance of the overall mean. The resulting ``N_eff = N p̂(1−p̂) / σ̂²_BM`` is
capped at ``N``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from statistics import NormalDist

FEW_HITS = 10


@dataclass(frozen=True)
class ProportionEstimate:
    hits: int
    trials: int
    confidence: float = 0.95
    effective_trials: float | None = None  # N_eff for correlated (MCMC) samples

    def __post_init__(self):
        if self.trials <= 0 or not 0 <= self.hits <= self.trials:
            raise ValueError(
                f"Need 0 <= hits <= trials and trials > 0, "
                f"got {self.hits}/{self.trials}."
            )
        if not 0 < self.confidence < 1:
            raise ValueError("confidence must be in (0, 1).")
        if self.effective_trials is not None and not (
            0 < self.effective_trials <= self.trials
        ):
            raise ValueError(
                f"effective_trials must be in (0, trials], got {self.effective_trials}."
            )

    @property
    def n_eff(self) -> float:
        """The sample size used for errors: ``effective_trials``, else ``trials``."""
        return self.trials if self.effective_trials is None else self.effective_trials

    @property
    def estimate(self) -> float:
        return self.hits / self.trials

    @property
    def standard_error(self) -> float:
        p = self.estimate
        return math.sqrt(p * (1 - p) / self.n_eff)

    @property
    def wilson_interval(self) -> tuple[float, float]:
        z = NormalDist().inv_cdf(0.5 + self.confidence / 2)
        n, p = self.n_eff, self.estimate
        centre = (p + z * z / (2 * n)) / (1 + z * z / n)
        half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
        return max(0.0, centre - half), min(1.0, centre + half)

    @property
    def zero_hit_upper_bound(self) -> float | None:
        """Exact one-sided upper confidence bound when no hit was observed."""
        if self.hits:
            return None
        return 1 - (1 - self.confidence) ** (1 / self.n_eff)

    def describe(self) -> str:
        samples = f"{self.trials:,} samples"
        if self.effective_trials is not None:
            samples += f", ESS {self.effective_trials:,.0f}"
        if self.hits == 0:
            return (
                f"0 occurrences in {samples} "
                f"({self.confidence:.0%} upper bound "
                f"{format_probability(self.zero_hit_upper_bound)})"
            )
        if self.hits < FEW_HITS:
            low, high = self.wilson_interval
            return (
                f"{self.hits} occurrence{'s' if self.hits > 1 else ''} "
                f"in {samples} "
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
            "effective_samples": self.effective_trials,
            "estimate": self.estimate,
            "standard_error": self.standard_error,
            "confidence": self.confidence,
            "wilson_interval": [low, high],
            "zero_hit_upper_bound": self.zero_hit_upper_bound,
        }


def batch_means_ess(
    batch_hits: Sequence[int], batch_size: int, hits: int, trials: int
) -> float | None:
    """Effective sample size of a hit fraction from non-overlapping batch means.

    ``batch_hits[b]`` is the number of hits in batch ``b`` of ``batch_size``
    consecutive samples. Only full batches are passed, while ``hits`` and
    ``trials`` count every sample. Returns ``None`` when the ESS cannot be
    estimated: fewer than two batches, or ``p̂`` equal to 0 or 1, or batch
    means with zero variance. The result is capped at ``trials``.
    """
    num_batches = len(batch_hits)
    p = hits / trials
    if num_batches < 2 or batch_size < 1 or p in (0.0, 1.0):
        return None
    means = [h / batch_size for h in batch_hits]
    centre = sum(means) / num_batches
    var_means = sum((m - centre) ** 2 for m in means) / (num_batches - 1)
    sigma2 = batch_size * var_means  # asymptotic variance of the chain average
    if sigma2 <= 0:
        return None
    return min(float(trials), trials * p * (1 - p) / sigma2)


def batch_means_standard_error(batch_means: Sequence[float]) -> float | None:
    """Standard error of an overall mean from the means of equal batches."""
    b = len(batch_means)
    if b < 2:
        return None
    centre = sum(batch_means) / b
    var = sum((m - centre) ** 2 for m in batch_means) / (b - 1)
    return math.sqrt(var / b)


def format_probability(p: float | None) -> str:
    """Four significant digits, scientific notation for small values."""
    if p is None:
        return "n/a"
    if p == 0:
        return "0"
    if abs(p) < 1e-3:
        return f"{p:.3e}"
    return f"{p:.4g}"
