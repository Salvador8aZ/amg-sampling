"""Sampling algorithms. Every sampler takes an explicit ``random.Random`` instance."""

from amg_sampling.samplers.uniform import matching_from_ordering, sample_matching

__all__ = ["matching_from_ordering", "sample_matching"]
