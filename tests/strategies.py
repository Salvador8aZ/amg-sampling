"""Hypothesis strategies shared by the test suite."""

from hypothesis import strategies as st


def breaks(max_chromosomes: int = 6, max_breaks: int = 6):
    """DSB counts per chromosome, each at least 1."""
    return st.lists(
        st.integers(min_value=1, max_value=max_breaks),
        min_size=1,
        max_size=max_chromosomes,
    )
