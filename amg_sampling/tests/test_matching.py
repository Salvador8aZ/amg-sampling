import pytest
from hypothesis import given

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import (
    InvalidMatchingError,
    RejoinMatching,
    require_compatible,
    validate_partner_map,
)
from amg_sampling.tests.strategies import matching_and_permutation, states

# -- deterministic -----------------------------------------------------------


def test_from_pairs_builds_partner_map():
    # Paper, figure 1: Θ(2,(1,1)), free ends b,c,f,g = 0,1,2,3; rejoins (b,g),(c,f).
    r = RejoinMatching.from_pairs(4, [(0, 3), (2, 1)])
    assert r.partners == (3, 2, 1, 0)
    assert r.pairs() == ((0, 3), (1, 2))
    assert r.num_ends == 4
    assert r.num_edges == 2
    assert r.partner(2) == 1
    assert r[0] == 3
    assert list(r) == [3, 2, 1, 0]


def test_equality_and_hash_are_by_partner_map():
    a = RejoinMatching.from_pairs(4, [(0, 3), (1, 2)])
    b = RejoinMatching([3, 2, 1, 0])
    assert a == b
    assert hash(a) == hash(b)
    assert a != RejoinMatching([1, 0, 3, 2])
    assert len({a, b}) == 1


def test_repr_round_trips():
    r = RejoinMatching([2, 3, 0, 1])
    assert eval(repr(r), {"RejoinMatching": RejoinMatching}) == r


def test_matching_may_pair_dsb_partners():
    # Repaired DSBs are valid rejoin matchings (they form C1 cycles).
    # Whether they are *admissible* is a state-space question, not a matching one.
    r = RejoinMatching([1, 0, 3, 2])
    assert r.pairs() == ((0, 1), (2, 3))


@pytest.mark.parametrize(
    "partners",
    [
        [0, 1],  # fixed points
        [1, 0, 2],  # odd length
        [1, 2, 0, 3],  # not an involution
        [1, 0, 3, 4],  # out of range
        [-1, 0],  # out of range
        [1.0, 0],  # not an int
        [True, False],  # bools are not end labels
    ],
)
def test_invalid_partner_maps_are_rejected(partners):
    with pytest.raises(InvalidMatchingError):
        RejoinMatching(partners)


@pytest.mark.parametrize(
    "num_ends, pairs",
    [
        (4, [(0, 1)]),  # end 2, 3 unmatched
        (4, [(0, 1), (1, 2)]),  # end 1 used twice
        (4, [(0, 1), (2, 2)]),  # self-loop
        (4, [(0, 1), (2, 5)]),  # out of range
    ],
)
def test_invalid_pairs_are_rejected(num_ends, pairs):
    with pytest.raises(InvalidMatchingError):
        RejoinMatching.from_pairs(num_ends, pairs)


def test_is_immutable():
    r = RejoinMatching([1, 0])
    with pytest.raises(AttributeError):
        r.extra = 1  # type: ignore[attr-defined]
    with pytest.raises(TypeError):
        r.partners[0] = 0  # type: ignore[index]


def test_require_compatible():
    theta = InitialConfiguration([2])
    require_compatible(theta, RejoinMatching([2, 3, 0, 1]))
    with pytest.raises(ValueError):
        require_compatible(theta, RejoinMatching([1, 0]))


def test_relabel_rejects_non_permutations():
    with pytest.raises(ValueError):
        RejoinMatching([1, 0]).relabel([0, 0])


# -- properties --------------------------------------------------------------


@given(states())
def test_generated_matchings_are_fixed_point_free_involutions(state):
    theta, r = state
    assert r.num_ends == theta.num_ends
    for v in range(r.num_ends):
        assert r[v] != v
        assert r[r[v]] == v


@given(states())
def test_pairs_round_trip(state):
    _, r = state
    pairs = r.pairs()
    assert len(pairs) == r.num_edges
    assert all(u < v for u, v in pairs)
    assert RejoinMatching.from_pairs(r.num_ends, pairs) == r


@given(states())
def test_validate_accepts_every_generated_matching(state):
    validate_partner_map(state[1].partners)


@given(matching_and_permutation())
def test_relabel_is_conjugation(data):
    r, perm = data
    s = r.relabel(perm)
    for v in range(r.num_ends):
        assert s[perm[v]] == perm[r[v]]
    inverse = [0] * len(perm)
    for v, p in enumerate(perm):
        inverse[p] = v
    assert s.relabel(inverse) == r
