"""Sampling algorithms. Every sampler takes an explicit ``random.Random`` instance."""

from amg_sampling.samplers.rejection import (
    SamplingStats,
    iter_samples,
    sample_state,
    theoretical_acceptance,
)
from amg_sampling.samplers.uniform import (
    matching_from_ordering,
    sample_completion,
    sample_matching,
)

__all__ = [
    "SamplingStats",
    "iter_samples",
    "matching_from_ordering",
    "sample_completion",
    "sample_matching",
    "sample_state",
    "theoretical_acceptance",
]
