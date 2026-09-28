import pytest
from hypothesis import given

from amg_sampling.core.configuration import TELOMERE, InitialConfiguration, Side, Telomere
from amg_sampling.tests.strategies import breaks


# -- deterministic -----------------------------------------------------------


def test_single_chromosome_two_dsbs_matches_paper_example_2():
    # Paper, example 2: Θ(1,(2)) with chromosome edges (a,b),(c,d),(e,f) and
    # DSB edges (b,c),(d,e). Free ends b,c,d,e are labelled 0,1,2,3.
    theta = InitialConfiguration([2])
    assert theta.num_chromosomes == 1
    assert theta.num_dsbs == 2
    assert theta.num_ends == 4
    assert theta.dsb_edges() == ((0, 1), (2, 3))
    assert theta.fragments == ((1, 2),)  # chromosome edge (c, d)
    assert theta.chromatin_partners == (TELOMERE, 2, 1, TELOMERE)
    assert theta.telomeres == (
        Telomere(0, Side.LEFT, 0),  # a - b
        Telomere(0, Side.RIGHT, 3),  # e - f
    )


def test_two_chromosomes_matches_paper_figure_2_layout():
    # Θ(2,(2,2)): chromosome 0 = a b c d e f, chromosome 1 = g h i j k l.
    theta = InitialConfiguration((2, 2))
    assert theta.dsb_chromosome == (0, 0, 1, 1)
    assert theta.dsb_edges() == ((0, 1), (2, 3), (4, 5), (6, 7))
    assert theta.fragments == ((1, 2), (5, 6))  # (c,d) and (i,j)
    assert [t.end for t in theta.telomeres] == [0, 3, 4, 7]
    assert list(theta.ends_of_chromosome(1)) == [4, 5, 6, 7]
    # No chromatin edge crosses chromosomes.
    assert theta.chromatin_partner(3) == TELOMERE
    assert theta.chromatin_partner(4) == TELOMERE


def test_one_dsb_per_chromosome_has_no_fragments():
    theta = InitialConfiguration((1, 1, 1))
    assert theta.fragments == ()
    assert theta.num_fragments == 0
    assert all(c == TELOMERE for c in theta.chromatin_partners)


def test_accessors():
    theta = InitialConfiguration((3, 1))
    assert theta.dsb_of(5) == 2
    assert theta.dsb_partner(5) == 4
    assert theta.chromatin_partner(3) == 4
    assert theta.chromosome_of(5) == 0
    assert theta.chromosome_of(6) == 1


@pytest.mark.parametrize("bad", [[], [0], [2, 0], [-1]])
def test_rejects_invalid_break_counts(bad):
    with pytest.raises(ValueError):
        InitialConfiguration(bad)


@pytest.mark.parametrize("bad", [[1.0], [True], ["2"]])
def test_rejects_non_integer_break_counts(bad):
    with pytest.raises(TypeError):
        InitialConfiguration(bad)


def test_out_of_range_end_raises():
    theta = InitialConfiguration([2])
    with pytest.raises(IndexError):
        theta.chromatin_partner(4)
    with pytest.raises(IndexError):
        theta.dsb_partner(-1)
    with pytest.raises(IndexError):
        theta.ends_of_chromosome(1)


def test_is_immutable_and_hashable():
    theta = InitialConfiguration([2, 1])
    with pytest.raises(AttributeError):
        theta.breaks = (1,)
    assert theta == InitialConfiguration((2, 1))
    assert hash(theta) == hash(InitialConfiguration((2, 1)))
    assert theta != InitialConfiguration((1, 2))


# -- properties --------------------------------------------------------------


@given(breaks())
def test_counts(bs):
    theta = InitialConfiguration(bs)
    k, n = len(bs), sum(bs)
    assert theta.num_dsbs == n
    assert theta.num_ends == 2 * n
    assert theta.num_fragments == n - k
    assert len(theta.telomeres) == 2 * k


@given(breaks())
def test_dsb_partner_is_fixed_point_free_involution(bs):
    theta = InitialConfiguration(bs)
    d = theta.dsb_partners
    for v in range(theta.num_ends):
        assert d[v] != v
        assert d[d[v]] == v
        assert theta.dsb_of(v) == theta.dsb_of(d[v])


@given(breaks())
def test_chromatin_partner_is_an_involution_on_non_telomeric_ends(bs):
    theta = InitialConfiguration(bs)
    c = theta.chromatin_partners
    for v in range(theta.num_ends):
        if c[v] != TELOMERE:
            assert c[v] != v
            assert c[c[v]] == v
            assert theta.chromosome_of(c[v]) == theta.chromosome_of(v)
            assert c[v] != theta.dsb_partner(v)


@given(breaks())
def test_every_end_has_exactly_one_chromatin_neighbour(bs):
    # Each free end is joined by chromatin either to one other free end
    # (a fragment) or to one telomere, never both.
    theta = InitialConfiguration(bs)
    fragment_ends = [v for f in theta.fragments for v in f]
    telomere_ends = [t.end for t in theta.telomeres]
    assert sorted(fragment_ends + telomere_ends) == list(range(theta.num_ends))
    for v in telomere_ends:
        assert theta.chromatin_partner(v) == TELOMERE


@given(breaks())
def test_chromosomes_are_paths_in_labelling_order(bs):
    # Walking telomere -> chromatin -> DSB -> chromatin -> ... -> telomere
    # visits the chromosome's ends in increasing label order.
    theta = InitialConfiguration(bs)
    for i in range(theta.num_chromosomes):
        left, right = theta.telomeres[2 * i], theta.telomeres[2 * i + 1]
        assert (left.chromosome, left.side) == (i, Side.LEFT)
        assert (right.chromosome, right.side) == (i, Side.RIGHT)
        walk, v = [], left.end
        while True:
            walk += [v, theta.dsb_partner(v)]
            nxt = theta.chromatin_partner(theta.dsb_partner(v))
            if nxt == TELOMERE:
                break
            v = nxt
        assert walk == list(theta.ends_of_chromosome(i))
        assert walk[-1] == right.end
