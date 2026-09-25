# IID uniform sampling

Covers `amg_sampling.samplers`. Labels as in
[`representation.md`](representation.md): **PROVED**, **VERIFIED(n ≤ N)**,
**CITED**, **CONJECTURE**.

## 1. Uniform perfect matchings by shuffled pairing

**Algorithm** (`sample_matching`). Put the free-end labels `0, …, 2n−1` in a
list, shuffle it uniformly, and pair positions `(0,1), (2,3), …, (2n−2, 2n−1)`
(`matching_from_ordering`).

**Theorem (PROVED).** If the ordering is uniform over the `(2n)!` permutations,
the resulting matching is uniform over the `(2n−1)!!` perfect matchings.

*Proof.* Fix a perfect matching `M` with pairs `e_1, …, e_n`. An ordering
produces `M` iff each of the `n` position slots `(2t, 2t+1)` holds the two ends
of one pair of `M`. Such orderings are obtained by choosing which pair goes in
which slot (`n!` ways) and, independently for each pair, which of its two ends
comes first (`2^n` ways). Distinct choices give distinct orderings, and every
ordering producing `M` arises this way. So exactly `2^n n!` orderings produce
`M`, independently of `M`, and

    P(M) = 2^n n! / (2n)! = 1 / (2n−1)!!,

using `(2n)! = 2^n n! · (2n−1)!!`. ∎

VERIFIED(n ≤ 4): mapping all `(2n)!` orderings gives every matching exactly
`2^n n!` times (`test_every_matching_arises_from_exactly_2n_nfactorial_orderings`).

**Random number generator.**
- `random.Random.shuffle` is the Fisher–Yates shuffle. Each index comes from
  `_randbelow`, which draws `k` random bits and rejects out-of-range values,
  so no modulo bias is introduced. With ideal random bits the ordering is
  exactly uniform.
- The underlying generator is the Mersenne Twister (MT19937), with a
  19 937-bit state. At n = 80 the number of orderings is `160! ≈ 2^946`, far
  below the generator's state space. The standard library is therefore
  adequate, and no runtime dependency (NumPy) is needed.
- It is a pseudo-random generator for simulation, not a cryptographic one.

**Reproducibility.**
- Every sampler takes an explicit `random.Random` instance. There is no hidden
  global state, and tests check that the global `random` state is untouched.
- The same seed with the same Python version gives the same sequence, and two
  instances with the same seed give identical sequences (tested).
- Python guarantees reproducibility across versions only for
  `Random.random()`. `shuffle` depends on `getrandbits` and `_randbelow`,
  whose use may change between Python versions. **Sequences are therefore
  only guaranteed reproducible on the same Python implementation and version.**
