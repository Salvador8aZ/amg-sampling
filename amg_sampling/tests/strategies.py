"""Hypothesis strategies shared by the test suite."""

from hypothesis import strategies as st

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import RejoinMatching


def breaks(max_chromosomes: int = 6, max_breaks: int = 6):
    """DSB counts per chromosome, each at least 1."""
    return st.lists(
        st.integers(min_value=1, max_value=max_breaks),
        min_size=1,
        max_size=max_chromosomes,
    )


def matchings(num_ends: int):
    """Perfect matchings on ``range(num_ends)``: shuffle the ends, pair consecutive ones."""
    return st.permutations(range(num_ends)).map(
        lambda p: RejoinMatching.from_pairs(num_ends, zip(p[::2], p[1::2]))
    )


@st.composite
def states(draw, max_chromosomes: int = 6, max_breaks: int = 6):
    """A pair (Θ, r) with r a rejoin matching on the free ends of Θ."""
    theta = InitialConfiguration(draw(breaks(max_chromosomes, max_breaks)))
    return theta, draw(matchings(theta.num_ends))


@st.composite
def matching_and_permutation(draw, max_ends: int = 40):
    """A perfect matching together with a permutation of its ends."""
    num_ends = 2 * draw(st.integers(min_value=1, max_value=max_ends // 2))
    return draw(matchings(num_ends)), draw(st.permutations(range(num_ends)))
