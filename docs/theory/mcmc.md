# The 2-switch move and Markov chain sampling

Covers `amg_sampling.core.moves` and `amg_sampling.samplers.mcmc`. Labels as in
[`representation.md`](representation.md): **PROVED**, **VERIFIED(n ≤ N)**,
**CITED**, **CONJECTURE**.

## 1. The 2-switch move

Let `r` be a rejoin matching, and let `u` and `v` be ends on two different
rejoin edges `{u, a}` and `{v, b}` (so `a = r(u)`, `b = r(v)`, and
`v ∉ {u, a}`). The **2-switch** `s = switch(r, u, v)` replaces those two edges
by `{u, v}` and `{a, b}` and keeps every other edge.

Two rejoin edges `{a, b}`, `{c, d}` can be reconnected in exactly two other
ways, `{a, c}, {b, d}` and `{a, d}, {b, c}`. Biologically, re-pairing the ends
of two junctions is the combinatorial analogue of an inversion (both junctions
on one chromosome: the segment between them is reversed) or of a reciprocal
exchange. This document does not assert that the move equals the reversal move
of Sheth et al.

**Movable ends.** With observed (fixed) rejoins `F`, only the free ends of `F`
may move. In a completion of `F` the free ends are matched among themselves,
so a switch among them never touches a fixed edge. Write `m` for the number of
movable rejoin edges (`m = n` with no fixed edges) and `f = 2m` for the number
of movable ends.

**Proposition 1 (PROVED).** A matching has exactly `m(m−1)` distinct switch
neighbours, and none of them is the matching itself.

*Proof.* Each unordered pair of movable edges (`C(m, 2)` choices) gives two
neighbours, one for each other reconnection of its four ends. A neighbour
differs from `r` in exactly the two chosen edges, so it determines the pair
(the edges of `r` it lacks) and the reconnection (the edges it has instead).
All `2 · C(m, 2) = m(m−1)` neighbours are therefore distinct, and none equals
`r`. ∎

**Proposition 2 (PROVED).** Each neighbour arises from exactly four ordered
pairs `(u, v)`.

*Proof.* The neighbour with new edges `{u, v}` and `{a, b}` is produced by
`(u, v)`, `(v, u)`, `(a, b)` and `(b, a)`, and by no other pair: the first end
chosen must lie on one of the new edges, and the second end is then its new
partner. ∎

**Corollary 3 (PROVED).** Choose `u` uniformly among the `f` movable ends,
then `v` uniformly among the other `f − 2` movable ends not on `u`'s edge. Each
neighbour is proposed with probability `4 / (f(f−2)) = 1 / (m(m−1))`. This is
the same for every state, so the proposal is symmetric:
`q(r → s) = q(s → r)`.

VERIFIED(n ≤ 5): every matching of `Θ(1,(n))` has `n(n−1)` distinct
neighbours (`test_every_matching_has_m_m_minus_1_distinct_neighbours`).
VERIFIED(n ≤ 4): each neighbour arises from exactly four ordered pairs,
also with a fixed edge (`test_each_neighbour_arises_from_exactly_four_ordered_pairs`,
`test_moves_with_fixed_edges_keep_them`).

## 2. Connectivity of the restricted switch graphs

The chain can reach every state of a state space `S` only if the switch graph
restricted to `S` is connected (irreducibility).

**Proposition 4 (PROVED).** The switch graph on `ALL`, and on the completions
of any fixed set `F`, is connected.

*Proof.* Let `r ≠ s` be matchings of the same movable ends. Their symmetric
difference is a union of alternating cycles. Take an end `u` with
`r(u) ≠ s(u)`, and set `v = s(u)`. Then `v ∉ {u, r(u)}`, and
`switch(r, u, v)` contains the edge `{u, v}` of `s`, so it agrees with `s` on
at least one more edge than `r` did and loses no common edge (both removed
edges were not in `s`). Repeating this reaches `s` in fewer than `m` steps. ∎

`DERANGED` and `PROPER` are not closed under switches, so Proposition 4 does
not settle them. The path in the proof can pass through a repaired DSB or a
disconnected AMG.

**Claim 5.**
- **VERIFIED(n ≤ 6), CONJECTURE for larger n:** the switch graph restricted to
  `DERANGED`, and to `PROPER`, is connected for every layout `Θ` with
  `n ≤ 6` (all 63 layouts; `test_every_state_space_is_connected_under_switches`,
  plus `…_n6`, marked slow).
- **VERIFIED(n ≤ 5), CONJECTURE for larger n:** the same holds for the
  completions of every fixed set `F` that leaves at least two movable edges.
  That is 90 650 (layout, `F`, space) cases for `n ≤ 5`
  (`test_completion_spaces_are_connected_under_switches` for `n ≤ 4`, plus
  `…_n5`, marked slow).
