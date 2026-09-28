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

## 3. The Markov chain (`samplers/mcmc.py`)

**Algorithm.** Start from a state `r₀ ∈ S`. By default this is one exact
uniform draw of the rejection sampler (see [`sampling.md`](sampling.md) §2 and §9).
At each step:

1. Draw `u` uniformly from the `f` movable ends, then `v` uniformly from the
   `f − 2` movable ends other than `u` and `r(u)`.
2. Propose `s = switch(r, u, v)`.
3. Move to `s` if `s ∈ S`; otherwise stay at `r`.

After `burn_in` steps, record the state after every `thin` steps.

**Proposition 6 (PROVED).** `Uniform(S)` is stationary. If the switch graph
restricted to `S` is connected (Claim 5), the chain is irreducible and
aperiodic, so its distribution converges to `Uniform(S)` from any start, and
averages along the chain converge to expectations under `Uniform(S)`.

*Proof.*
- **Stationarity.** For `r ≠ s` in `S`, `P(r → s) = q(r → s) = 1/(m(m−1))` if
  `s` is a neighbour of `r`, and 0 otherwise. This is symmetric in `r` and `s`
  (Corollary 3), so the uniform distribution satisfies detailed balance. This
  is Metropolis–Hastings with a symmetric proposal and target `π ∝ 1_S`, whose
  acceptance probability `min(1, π(s)/π(r))` is 1 inside `S` and 0 outside it.
- **Irreducibility.** This is exactly connectivity of the restricted switch
  graph.
- **Aperiodicity.** If some state of `S` has a neighbour outside `S`, that
  state has a self-loop. Otherwise `S` is closed under switches, and since
  the switch graph on the movable ends is connected (Proposition 4), `S` is
  the whole of `ALL` or of the completions of `F`. There, the two
  reconnections of one pair of edges and the original state form a triangle,
  an odd cycle (for `m ≥ 2`). Either way the chain is aperiodic. ∎

For a uniform target the chain adds nothing that IID rejection does not
already provide (§7 of [`sampling.md`](sampling.md)). It is the foundation for
non-uniform targets: a biologically weighted target `π(r) ∝ w(r)` on `S` needs
only the acceptance probability `min(1, w(s)/w(r))` in step 3. Rejection from
uniform proposals would no longer give that distribution.

**Implementation notes.**
- For `DERANGED` only the two new edges can create a `C₁`, so membership is
  checked in `O(1)`. For `PROPER`, connectivity is recomputed by union-find in
  `O(n α(k))` per step.
- A rejected proposal is undone by a second switch, `switch(s, u, r(u))`.
- With at most one movable edge, the state space has a single state and the
  chain never moves.

**Verification.**

| Level | Check | Test |
|---|---|---|
| 0: kernel | Rows and columns of the exact transition matrix sum to 1, and it is symmetric, for every layout with `n ≤ 4` and every space, and with fixed rejoins | `test_kernel_is_symmetric_and_preserves_uniform`, `test_kernel_with_fixed_rejoins` |
| 1: steps | The fast chain's one-step distribution matches the exact kernel (χ²) | `test_one_step_distribution_matches_the_kernel` |
| 2: distributions | Chain frequencies of cycle structures match exact distributions within 4 batch-means standard errors, for `PROPER`, `DERANGED` and `ALL` | `test_chain_reproduces_exact_cycle_distributions` |
| 2: convergence | From the atypical start `C₆`, the chain forgets it after burn-in | `test_chain_converges_from_an_atypical_start` |
| 2: large `n` | `Θ(1,(40))` `DERANGED`: number of cycles against the exact distribution | `test_chain_at_n40_matches_the_exact_number_of_cycles` |
| 2: patient | P05-1657, chromosomes 8 and 12: the chain on the 945 completions reproduces the paper's table 2 | `test_mcmc_reproduces_the_p05_1657_completions` (needs the dataset) |

## 4. Standard errors for correlated samples

Consecutive chain states are correlated, so `N` samples carry the information
of fewer, `N_eff`, independent ones. For each reported probability the
analyses estimate `N_eff` by non-overlapping batch means
(`analysis.estimates.batch_means_ess`):

- **Batches.** Cut the `N` samples into `B = ⌊N/L⌋` batches of
  `L = ⌊√N⌋` consecutive samples.
- **Variance.** `σ̂²_BM = L · Var(batch means)` estimates the asymptotic
  variance of the chain average.
- **Effective size.** `N_eff = N p̂(1−p̂) / σ̂²_BM`, capped at `N`.

`N_eff` then replaces `N` in the standard error, the Wilson interval and the
zero-hit upper bound (`analysis.estimates`). When an indicator's own `N_eff`
cannot be estimated (no hits, every sample a hit, or constant batches), the
`N_eff` of the number of cycles is used instead. The report prints that
`N_eff` and the integrated autocorrelation time `τ = N / N_eff` as a mixing
diagnostic.

VERIFIED: on IID Bernoulli samples `N_eff ≈ N`, and on a two-state chain with
flip probability `q` it matches the known `N q/(1−q)` within 25%
(`test_batch_means_ess_*`).

**Limitations.**
- Batch means with `L = √N` is consistent but can underestimate the variance
  for chains that mix slowly relative to `L`.
- A single chain cannot detect a region of `S` it never visits. Several chains
  from dispersed starts (R̂) are future work, and matter most once targets are
  non-uniform.
