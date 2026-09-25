from collections import Counter
from fractions import Fraction
from math import comb, factorial

import pytest
from hypothesis import given
from hypothesis import strategies as st

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import CycleStructure
from amg_sampling.core.statespace import StateSpace
from amg_sampling.exact.enumerate import count_states, exact_cycle_distribution
from amg_sampling.exact.formulas import (
    cycle_count_all,
    cycle_distribution_all,
    cycle_distribution_deranged,
    cycle_distribution_proper,
    integer_partitions,
    lemma15_two_cycle_count,
    num_all_states,
    num_deranged_states,
    num_deranged_states_closed_form,
    num_proper_states,
    num_single_cycle_states,
    num_states_with_cycle_structure,
    num_two_cycle_deranged_states,
    unsigned_stirling_first_kind,
)
from tests.paper_tables import TABLE_1, TABLE_3
from tests.reference import compositions
from tests.strategies import breaks

SMALL_N = [1, 2, 3, 4, 5, 6]


def by_string(dist):
    return {str(c): m for c, m in dist.items()}


def small_compositions(max_n):
    """Every DSB distribution with 1 <= n <= max_n, as pytest params."""
    return [pytest.param(b, id=str(b)) for n in range(1, max_n + 1) for b in compositions(n)]


# -- partitions ----------------------------------------------------------------


@pytest.mark.parametrize("n, p", [(0, 1), (1, 1), (4, 5), (6, 11), (10, 42), (20, 627)])
def test_partition_counts(n, p):
    parts = list(integer_partitions(n))
    assert len(parts) == p == len(set(parts))
    assert all(c.n == n for c in parts)


# -- ALL: theorems 1 and 2 ----------------------------------------------------------


@pytest.mark.parametrize("n", SMALL_N)
def test_theorem_1_against_enumeration(n):
    assert num_all_states(n) == count_states(InitialConfiguration([n]), StateSpace.ALL)


@pytest.mark.parametrize("n", range(0, 30))
def test_theorem_1_double_factorial(n):
    assert num_all_states(n) == factorial(2 * n) // (2**n * factorial(n))


@pytest.mark.parametrize("n", SMALL_N)
def test_theorem_2_against_enumeration(n):
    exact = exact_cycle_distribution(InitialConfiguration([n]), StateSpace.ALL)
    formula = {c: m for c, m in cycle_distribution_all(n).items()}
    assert exact == formula  # every partition of n occurs in ALL


@pytest.mark.parametrize("n", range(1, 21))
def test_theorem_2_sums_to_theorem_1(n):
    assert sum(cycle_distribution_all(n).values()) == num_all_states(n)


@pytest.mark.parametrize("n", range(1, 21))
def test_corollary_3_is_theorem_2_for_single_cycle(n):
    assert num_single_cycle_states(n) == num_states_with_cycle_structure(CycleStructure([n]))


def test_theorem_2_examples():
    assert num_states_with_cycle_structure(CycleStructure([4])) == 48
    assert num_states_with_cycle_structure(CycleStructure([3, 1])) == 32
    assert num_states_with_cycle_structure(CycleStructure([2, 2])) == 12
    assert num_states_with_cycle_structure(CycleStructure([1, 1, 1, 1])) == 1


# -- number of cycles: 2^{n-c} [n c] (derived) --------------------------------------


def paper_section_3_1_printed_count(n, c):
    """The expression printed in section 3.1 of the paper, kept only for regression tests."""
    return unsigned_stirling_first_kind(n, c)


def enumerated_cycle_counts(n, space=StateSpace.ALL):
    counts = Counter()
    for c, m in exact_cycle_distribution(InitialConfiguration([n]), space).items():
        counts[c.num_cycles] += m
    return counts


@pytest.mark.parametrize("n", SMALL_N)
def test_cycle_count_all_against_enumeration(n):
    exact = enumerated_cycle_counts(n)
    assert dict(exact) == {c: cycle_count_all(n, c) for c in range(1, n + 1)}


@pytest.mark.parametrize("n", SMALL_N)
def test_printed_section_3_1_expression_differs_from_enumeration(n):
    # The printed [n c] matches neither ALL nor DERANGED (= PROPER for k = 1)
    # for n >= 2; the ratio to ALL is the orientation factor 2^{n-c}.
    all_counts = enumerated_cycle_counts(n)
    deranged = enumerated_cycle_counts(n, StateSpace.DERANGED)
    for c in range(1, n + 1):
        printed = paper_section_3_1_printed_count(n, c)
        assert all_counts[c] == 2 ** (n - c) * printed
        if n >= 2 and c < n:
            assert all_counts[c] != printed
            assert deranged[c] != printed


@pytest.mark.parametrize("n", range(0, 41))
def test_cycle_counts_normalise_to_theorem_1(n):
    # Σ_c 2^{n-c} [n c] = (2n-1)!!
    assert sum(cycle_count_all(n, c) for c in range(0, n + 1)) == num_all_states(n)


@pytest.mark.parametrize("n", range(1, 21))
def test_three_level_consistency(n):
    # κ_C  ->  2^{n-c} [n c]  ->  (2n-1)!!
    by_cycles = Counter()
    for c in integer_partitions(n):
        by_cycles[c.num_cycles] += num_states_with_cycle_structure(c)
    for c in range(1, n + 1):
        assert by_cycles[c] == cycle_count_all(n, c)
    assert sum(by_cycles.values()) == num_all_states(n)


def mean_cycles_ewens_half(n):
    """Derived (not used by the library): E[#cycles | uniform ALL] = Σ_{i<n} 1/(2i+1)."""
    return sum(Fraction(1, 2 * i + 1) for i in range(n))


@pytest.mark.parametrize("n", SMALL_N)
def test_mean_number_of_cycles_against_enumeration(n):
    counts = enumerated_cycle_counts(n)
    mean = Fraction(sum(c * m for c, m in counts.items()), sum(counts.values()))
    assert mean == mean_cycles_ewens_half(n)


@pytest.mark.parametrize("n", range(1, 31))
def test_mean_number_of_cycles_against_cycle_counts(n):
    total = num_all_states(n)
    mean = Fraction(sum(c * cycle_count_all(n, c) for c in range(1, n + 1)), total)
    assert mean == mean_cycles_ewens_half(n)


def test_stirling_values():
    assert [unsigned_stirling_first_kind(5, k) for k in range(6)] == [0, 24, 50, 35, 10, 1]
    assert sum(unsigned_stirling_first_kind(7, k) for k in range(8)) == factorial(7)


# -- DERANGED: theorem 5, equation 3, proposition 8, corollary 7 ---------------------


@pytest.mark.parametrize("n", range(0, 40))
def test_theorem_5_recurrence_equals_closed_form(n):
    assert num_deranged_states(n) == num_deranged_states_closed_form(n)


@pytest.mark.parametrize("n", range(0, 25))
def test_equation_3_sum_equals_theorem_5(n):
    assert sum(cycle_distribution_deranged(n).values()) == num_deranged_states(n)


def test_deranged_counts_oeis_a053871():
    assert [num_deranged_states(n) for n in range(9)] == [1, 0, 2, 8, 60, 544, 6040, 79008, 1190672]


@pytest.mark.parametrize("b", small_compositions(6))
def test_deranged_against_enumeration_for_every_layout(b):
    theta = InitialConfiguration(b)
    n = theta.num_dsbs
    assert count_states(theta, StateSpace.DERANGED) == num_deranged_states(n)
    assert exact_cycle_distribution(theta, StateSpace.DERANGED) == cycle_distribution_deranged(n)


@pytest.mark.parametrize("n", range(0, 22))
def test_proposition_8_against_theorem_2(n):
    two = sum(m for c, m in cycle_distribution_deranged(n).items() if c.num_cycles == 2)
    assert num_two_cycle_deranged_states(n) == two


@pytest.mark.parametrize("n", range(4, 22))
def test_corollary_7_against_theorem_2(n):
    for l in range(2, n // 2 + 1):
        c = CycleStructure([l, n - l])
        if 2 * l == n:
            expected = 2 ** (n - 1) * factorial(n) // n**2
        else:
            expected = 2 ** (n - 2) * factorial(n) // (l * (n - l))
        assert num_states_with_cycle_structure(c) == expected


# -- PROPER: derived recursion vs enumeration and the paper's lemmas --------------


@pytest.mark.parametrize("b", small_compositions(6))
def test_proper_distribution_against_enumeration(b):
    theta = InitialConfiguration(b)
    exact = exact_cycle_distribution(theta, StateSpace.PROPER)
    assert cycle_distribution_proper(theta) == exact
    assert num_proper_states(theta) == sum(exact.values())


@pytest.mark.slow
@pytest.mark.parametrize("b", [(7,), (4, 3), (3, 2, 2), (2, 2, 1, 1, 1), (1,) * 7])
def test_proper_distribution_against_enumeration_n7(b):
    theta = InitialConfiguration(b)
    for space, formula in [
        (StateSpace.ALL, cycle_distribution_all(7)),
        (StateSpace.DERANGED, cycle_distribution_deranged(7)),
        (StateSpace.PROPER, cycle_distribution_proper(theta)),
    ]:
        exact = exact_cycle_distribution(theta, space)
        assert exact == {c: m for c, m in formula.items() if m}


@pytest.mark.parametrize("b", list(TABLE_1) + list(TABLE_3))
def test_proper_distribution_reproduces_paper_tables(b):
    expected = {**TABLE_1, **TABLE_3}[b]
    assert by_string(cycle_distribution_proper(InitialConfiguration(b))) == expected


@pytest.mark.parametrize("n", range(2, 16))
def test_lemma_10_two_chromosomes(n):
    # κ'(2,(l, n-l)) = κ'(n) − κ'(l) κ'(n−l). The paper states it for
    # ⌊n/2⌋ <= l <= n−2; the derivation holds for every 1 <= l <= n−1.
    k = num_deranged_states
    for l in range(1, n):
        assert num_proper_states(InitialConfiguration((l, n - l))) == k(n) - k(l) * k(n - l)


@pytest.mark.parametrize("n", range(3, 13))
def test_lemma_12_three_chromosomes(n):
    # The paper states lemma 12 for l1 >= l2 >= l3 >= 2; checked here for all l_i >= 1.
    k = num_deranged_states
    for l1 in range(1, n - 1):
        for l2 in range(1, n - l1):
            l3 = n - l1 - l2
            ls = (l1, l2, l3)
            formula = k(n) - sum(k(li) * k(n - li) for li in ls) + 2 * k(l1) * k(l2) * k(l3)
            assert num_proper_states(InitialConfiguration(ls)) == formula


@pytest.mark.parametrize("n", range(1, 11))
def test_lemma_14_one_dsb_per_chromosome(n):
    theta = InitialConfiguration([1] * n)
    expected = {CycleStructure([n]): num_single_cycle_states(n)} if n >= 2 else {}
    assert cycle_distribution_proper(theta) == expected


def paper_lemma_15_ii_printed(n):
    """Lemma 15 (ii) exactly as printed, kept only for regression tests."""
    return 2 ** (n - 1) * factorial(n - 1) * (n - 2)


def lemma_15_theta(num_dsbs):
    return InitialConfiguration((2,) + (1,) * (num_dsbs - 2))


@pytest.mark.parametrize("num_dsbs", range(3, 11))
def test_lemma_15_against_proper_distribution(num_dsbs):
    dist = cycle_distribution_proper(lemma_15_theta(num_dsbs))
    two = {c: m for c, m in dist.items() if c.num_cycles == 2}
    assert dist[CycleStructure([num_dsbs])] == num_single_cycle_states(num_dsbs)  # lemma 15 (i)
    assert sum(two.values()) == lemma15_two_cycle_count(num_dsbs)
    assert all(c.num_cycles <= 2 for c in dist)
    # Structures C_j + C_{N-j}, 2 <= j <= N-2, each ordered choice contributing 2^{N-2}(N-2)!.
    unit = 2 ** (num_dsbs - 2) * factorial(num_dsbs - 2)
    for j in range(2, num_dsbs - 1):
        c = CycleStructure([j, num_dsbs - j])
        assert two[c] == (unit if 2 * j == num_dsbs else 2 * unit)


@pytest.mark.parametrize("num_dsbs", range(3, 7))
def test_lemma_15_against_enumeration(num_dsbs):
    dist = exact_cycle_distribution(lemma_15_theta(num_dsbs), StateSpace.PROPER)
    two = sum(m for c, m in dist.items() if c.num_cycles == 2)
    assert two == lemma15_two_cycle_count(num_dsbs)
    assert two != paper_lemma_15_ii_printed(num_dsbs)


@pytest.mark.parametrize("n", range(2, 30))
def test_lemma_15_printed_expression_is_the_restatement_shifted_by_one(n):
    assert paper_lemma_15_ii_printed(n) == lemma15_two_cycle_count(n + 1)


@pytest.mark.parametrize("bad", [0, 1, 2])
def test_lemma_15_domain(bad):
    with pytest.raises(ValueError):
        lemma15_two_cycle_count(bad)


@pytest.mark.parametrize("b", small_compositions(6))
def test_proposition_4_cycle_counts_are_realised(b):
    # Proposition 4: for every c in [1, min(max b_j, ⌊n/2⌋)] some proper AMG has c cycles.
    n = sum(b)
    realised = {c.num_cycles for c in cycle_distribution_proper(InitialConfiguration(b))}
    for c in range(1, min(max(b), n // 2) + 1):
        assert c in realised


# -- properties at sizes beyond enumeration --------------------------------------


@given(breaks(max_chromosomes=5, max_breaks=4))
def test_proper_distribution_is_consistent(bs):
    theta = InitialConfiguration(bs)
    n = theta.num_dsbs
    dist = cycle_distribution_proper(theta)
    assert sum(dist.values()) == num_proper_states(theta) <= num_deranged_states(n)
    deranged = cycle_distribution_deranged(n)
    for c, m in dist.items():
        assert c.n == n and c.count(1) == 0
        assert 0 < m <= deranged[c]
    if n >= 2:
        # Derived: a single exchange cycle is always deranged and connected.
        assert dist[CycleStructure([n])] == num_single_cycle_states(n)


@given(st.permutations([3, 2, 2, 1]))
def test_proper_counts_do_not_depend_on_chromosome_order(order):
    theta = InitialConfiguration(order)
    assert cycle_distribution_proper(theta) == cycle_distribution_proper(
        InitialConfiguration((3, 2, 2, 1))
    )


@pytest.mark.parametrize("bad", [-1, 1.5, True])
def test_invalid_n(bad):
    with pytest.raises(ValueError):
        num_all_states(bad)


def test_closed_form_uses_exact_arithmetic():
    # Large n would lose precision in floating point.
    n = 60
    assert num_deranged_states_closed_form(n) == num_deranged_states(n)
    assert num_all_states(n) == comb(2 * n, n) * factorial(n) // 2**n
