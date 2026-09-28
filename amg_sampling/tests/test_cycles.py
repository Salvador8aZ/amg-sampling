import pytest
from hypothesis import given
from hypothesis import strategies as st

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import (
    CycleStructure,
    cycle_structure,
    cycle_structure_from_partners,
    exchange_cycles,
    repaired_dsbs,
    ring_fragments,
)
from amg_sampling.core.matching import InvalidMatchingError, RejoinMatching
from amg_sampling.tests.reference import (
    all_perfect_matchings,
    amg_multigraph,
    reference_cycle_parts,
)
from amg_sampling.tests.strategies import breaks, matching_and_permutation, matchings, states

# Paper letter labels for Θ(2,(2,2)) (figures 2 and 4): chromosome 0 is
# a-b-c-d-e-f, chromosome 1 is g-h-i-j-k-l; telomeres a, f, g, l.
L22 = dict(b=0, c=1, d=2, e=3, h=4, i=5, j=6, k=7)


def letters(labels, *pairs):
    return RejoinMatching.from_pairs(len(labels), [(labels[p[0]], labels[p[1]]) for p in pairs])


# -- CycleStructure value type -----------------------------------------------


def test_cycle_structure_is_a_sorted_partition():
    c = CycleStructure([1, 3, 2, 2])
    assert c.parts == (3, 2, 2, 1)
    assert c == CycleStructure((2, 1, 3, 2))
    assert c.n == 8
    assert c.num_cycles == 4
    assert c.rank == 4
    assert c.multiplicities() == {3: 1, 2: 2, 1: 1}
    assert c.count(2) == 2 and c.count(5) == 0
    assert str(c) == "C3+2C2+C1"


@pytest.mark.parametrize(
    "parts, text",
    [((4,), "C4"), ((2, 2), "2C2"), ((3, 1), "C3+C1"), ((1, 1, 1), "3C1"), ((), "∅")],
)
def test_cycle_structure_string(parts, text):
    assert str(CycleStructure(parts)) == text


def test_from_multiplicities():
    assert CycleStructure.from_multiplicities({2: 3, 4: 1, 1: 0}) == CycleStructure((4, 2, 2, 2))
    with pytest.raises(ValueError):
        CycleStructure.from_multiplicities({2: -1})


@pytest.mark.parametrize("bad", [[0], [-1], [1.5], [True]])
def test_cycle_structure_rejects_non_positive_parts(bad):
    with pytest.raises(ValueError):
        CycleStructure(bad)


# -- paper examples ----------------------------------------------------------


def test_figure_2_c4_versus_two_c2():
    theta = InitialConfiguration((2, 2))
    # Exchange process b c h i d e k j: C4.
    c4 = letters(L22, "ch", "id", "ek", "jb")
    # Exchange processes b c h i + d e k j: C2 + C2.
    c22 = letters(L22, "ch", "ib", "ek", "jd")
    assert cycle_structure(theta, c4) == CycleStructure([4])
    assert cycle_structure(theta, c22) == CycleStructure([2, 2])


def test_figure_4_states_before_and_after_reversal_of_cd():
    # Only the two states are checked here; the reversal move itself is Stage C.
    theta = InitialConfiguration((2, 2))
    omega = letters(L22, "ci", "bh", "ek", "dj")  # b c i h + d e k j
    omega_prime = letters(L22, "cj", "di", "bh", "ek")  # b c j k e d i h
    assert str(cycle_structure(theta, omega)) == "2C2"
    assert str(cycle_structure(theta, omega_prime)) == "C4"


def test_example_2_inversion_and_ring_are_both_c2():
    theta = InitialConfiguration([2])  # ends b, c, d, e = 0, 1, 2, 3
    inversion = RejoinMatching.from_pairs(4, [(0, 2), (1, 3)])  # (b,d), (c,e)
    ring = RejoinMatching.from_pairs(4, [(0, 3), (1, 2)])  # (b,e), (c,d)
    assert cycle_structure(theta, inversion) == CycleStructure([2])
    assert cycle_structure(theta, ring) == CycleStructure([2])
    assert ring_fragments(theta, inversion) == ()
    assert ring_fragments(theta, ring) == (0,)  # fragment (c, d)


def test_figure_1_two_chromosomes_one_dsb_each():
    theta = InitialConfiguration((1, 1))  # ends b, c, f, g = 0, 1, 2, 3
    r = RejoinMatching.from_pairs(4, [(0, 3), (1, 2)])  # (b,g), (c,f)
    assert cycle_structure(theta, r) == CycleStructure([2])


# -- C1 regression tests (old repository dropped these) ----------------------


def test_every_dsb_repaired_is_n_c1():
    theta = InitialConfiguration([3])
    r = RejoinMatching.from_pairs(6, [(0, 1), (2, 3), (4, 5)])
    assert cycle_structure(theta, r) == CycleStructure([1, 1, 1])
    assert str(cycle_structure(theta, r)) == "3C1"
    assert repaired_dsbs(theta, r) == (0, 1, 2)


def test_one_repaired_dsb_gives_c2_plus_c1():
    # The old simple-graph implementation reported only {4: 1} here.
    theta = InitialConfiguration([3])
    r = RejoinMatching.from_pairs(6, [(0, 1), (2, 4), (3, 5)])
    assert cycle_structure(theta, r) == CycleStructure([2, 1])
    assert repaired_dsbs(theta, r) == (0,)
    assert exchange_cycles(theta.dsb_partners, r.partners) == ((0, 1), (2, 3, 5, 4))


def test_c3_plus_c1_from_section_7_2():
    # Θ(1,(4)): DSB 0 repaired, DSBs 1-3 in one cycle 2-3-4-5-6-7-2.
    theta = InitialConfiguration([4])
    r = RejoinMatching.from_pairs(8, [(0, 1), (3, 4), (5, 6), (7, 2)])
    assert str(cycle_structure(theta, r)) == "C3+C1"
    assert repaired_dsbs(theta, r) == (0,)


def test_ring_and_repaired_dsb_together_keep_every_edge():
    # Θ(1,(3)): DSB 0 repaired (C1); fragment (3,4) closed into a ring;
    # DSBs 1 and 2 form a C2 through the ring's rejoin edge.
    theta = InitialConfiguration([3])
    pairs = [(0, 1), (3, 4), (2, 5)]
    r = RejoinMatching.from_pairs(6, pairs)
    assert str(cycle_structure(theta, r)) == "C2+C1"
    assert repaired_dsbs(theta, r) == (0,)
    assert ring_fragments(theta, r) == (1,)

    # The chromatin, DSB and rejoin relations are all still present.
    assert theta.chromatin_partner(3) == 4 and r[3] == 4
    assert theta.dsb_partner(0) == 1 and r[0] == 1
    g = amg_multigraph(theta.breaks, pairs)
    assert set(g[3][4]) == {"chromatin", "rejoin"}
    assert set(g[0][1]) == {"dsb", "rejoin"}
    n, k = theta.num_dsbs, theta.num_chromosomes
    assert g.number_of_edges() == (n + k) + n + n


# -- errors ------------------------------------------------------------------


def test_incompatible_sizes_raise():
    with pytest.raises(ValueError):
        cycle_structure(InitialConfiguration([2]), RejoinMatching([1, 0]))
    with pytest.raises(ValueError):
        cycle_structure_from_partners([1, 0, 3, 2], [1, 0])
    with pytest.raises(InvalidMatchingError):
        exchange_cycles([0, 1], [1, 0])


# -- exhaustive and property checks ------------------------------------------


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_all_matchings_agree_with_networkx(n):
    theta = InitialConfiguration([n])
    for pairs in all_perfect_matchings(2 * n):
        r = RejoinMatching.from_pairs(2 * n, pairs)
        assert cycle_structure(theta, r).parts == reference_cycle_parts(n, pairs)


@given(states())
def test_agrees_with_networkx_multigraph(state):
    theta, r = state
    assert cycle_structure(theta, r).parts == reference_cycle_parts(theta.num_dsbs, r.pairs())


@given(states())
def test_parts_form_a_partition_of_n(state):
    theta, r = state
    c = cycle_structure(theta, r)
    assert c.n == theta.num_dsbs
    assert all(p >= 1 for p in c.parts)
    assert list(c.parts) == sorted(c.parts, reverse=True)


@given(states())
def test_c1_count_equals_repaired_dsbs(state):
    theta, r = state
    assert cycle_structure(theta, r).count(1) == len(repaired_dsbs(theta, r))


@given(states())
def test_exchange_cycles_alternate_and_cover_every_end(state):
    theta, r = state
    cycles = exchange_cycles(theta.dsb_partners, r.partners)
    assert sorted(v for c in cycles for v in c) == list(range(theta.num_ends))
    for c in cycles:
        assert c[0] == min(c)
        for i in range(0, len(c), 2):
            assert theta.dsb_partner(c[i]) == c[i + 1]
            assert r[c[i + 1]] == c[(i + 2) % len(c)]
    assert sorted((len(c) // 2 for c in cycles), reverse=True) == list(
        cycle_structure(theta, r).parts
    )


@given(matching_and_permutation())
def test_invariant_under_relabelling_both_matchings(data):
    # Renaming the ends (conjugating d and r by the same permutation) is a
    # graph isomorphism of Ξ, so it cannot change the cycle structure.
    r, perm = data
    d = RejoinMatching([v ^ 1 for v in range(r.num_ends)])
    before = cycle_structure_from_partners(d.partners, r.partners)
    after = cycle_structure_from_partners(d.relabel(perm).partners, r.relabel(perm).partners)
    assert before == after


@given(st.data())
def test_invariant_under_dsb_preserving_relabelling(data):
    # Permuting DSBs and swapping the two ends of some DSBs fixes d, so the
    # cycle structure of Θ with the relabelled rejoin matching is unchanged.
    n = data.draw(st.integers(1, 20))
    theta = InitialConfiguration([n])
    r = data.draw(matchings(2 * n))
    sigma = data.draw(st.permutations(range(n)))
    flips = data.draw(st.lists(st.booleans(), min_size=n, max_size=n))
    perm = [2 * sigma[v >> 1] + ((v & 1) ^ flips[v >> 1]) for v in range(2 * n)]
    assert cycle_structure(theta, r.relabel(perm)) == cycle_structure(theta, r)


@given(st.data())
def test_independent_of_chromosome_layout(data):
    # Ξ contains no chromatin edges, so moving DSBs between chromosomes
    # (same n, same labels) leaves the cycle structure unchanged.
    bs = data.draw(breaks())
    theta = InitialConfiguration(bs)
    single = InitialConfiguration([theta.num_dsbs])
    r = data.draw(matchings(theta.num_ends))
    assert cycle_structure(theta, r) == cycle_structure(single, r)
