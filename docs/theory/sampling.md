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

## 2. Rejection preserves uniformity

**Algorithm** (`sample_state`, `iter_samples`). Repeatedly draw `r` uniformly
from ALL and return the first draw that lies in the state space `X`.

**Theorem (PROVED).** The returned state is uniform on `X`.

*Proof.* Let `Y_1, Y_2, …` be independent uniform draws from ALL and
`τ = min{t : Y_t ∈ X}` (finite with probability 1 when `X ≠ ∅`). For `x ∈ X`,

    P(Y_τ = x) = Σ_t P(Y_1 ∉ X)^{t−1} P(Y_t = x) = (1/|ALL|) / (|X|/|ALL|) = 1/|X|. ∎

## 3. Accepted samples are IID

**Theorem (PROVED).** If proposals `Y_1, Y_2, …` are IID and each is accepted
iff it lies in the fixed set `X`, then the accepted proposals
`Z_1, Z_2, …` are IID with law `P(Y = x | Y ∈ X)`.

*Proof.* The acceptance times `τ_1 < τ_2 < …` split the proposal sequence into
blocks `(Y_{τ_{j−1}+1}, …, Y_{τ_j})`. Each block is determined by the
proposals inside it. Because the proposals are IID and the stopping rule
("stop at the first member of `X`") is the same for every block, the blocks
are IID. `Z_j = Y_{τ_j}` is a fixed function of block `j`, so the `Z_j` are
IID. Each has the law computed in §2. ∎

With uniform proposals the law is uniform on `X`. There is no Markov chain:
no burn-in, no autocorrelation, no mixing time. The independence test in
`tests/test_rejection.py` (non-overlapping consecutive pairs are uniform on
`X × X`) is only a sanity check; the argument above is what establishes independence.

In practice "independent" means independent up to the quality of the
pseudo-random generator (§1).

## 4. Acceptance probability

    acceptance(X) = |X| / |ALL|,   |ALL| = (2n−1)!!

`theoretical_acceptance(theta, space)` returns this as a `Fraction`:

| X | Numerator | Source | Cost |
|---|---|---|---|
| ALL | `(2n−1)!!` | theorem 1 | O(n) |
| DERANGED | `κ'(n)` | theorem 5; layout-independent (exact.md §1) | O(n) |
| PROPER, all `b_i = 1` | `2^{n−1}(n−1)!` (n ≥ 2) | lemma 14; proof below | O(n) |
| PROPER, general Θ | `num_proper_states(Θ)` | derived recursion (exact.md §5) | O(3^k); refused above `max_chromosomes` (default 12) |

**PROPER with one DSB per chromosome (PROVED).** If every chromosome has one
DSB, each chromosome's two ends are the two ends of one DSB, so they lie on
the same exchange cycle. The rejoin edges of a cycle stay inside the cycle.
So the chromosomes of different cycles are never linked, and the AMG is
connected iff there is exactly one cycle. A single cycle `C_n` (n ≥ 2) is
deranged. Hence `PROPER = {C_n states}`, of size `2^{n−1}(n−1)!` (corollary 3).
The paper states the count (lemma 14). This argument also shows the
states are exactly the single-cycle ones. Checked against the recursion for n ≤ 10.

**Empty state spaces.** DERANGED and PROPER are empty iff `n = 1`, since for
`n ≥ 2` every `C_n` state is proper (exact.md §5, P3). The sampler refuses an
empty space instead of looping forever. `max_proposals` optionally bounds the
work per sample.

**Empirical acceptance** (`SamplingStats`: `proposals`, `accepted`, `rejected`,
`acceptance_rate`) is kept separate from the theoretical value. Because the
sampler stops after a fixed number `a` of acceptances, the number of
rejections is negative binomial with mean `a(1−p)/p` and variance
`a(1−p)/p²`. The acceptance tests use that standard error. (With a fixed
number `N` of proposals, the binomial standard error `√(p(1−p)/N)` would
apply instead.)

The expected number of proposals per accepted sample is `1/p`.

## 5. ALL, DERANGED and PROPER as sampling targets

The same proposal distribution (uniform on ALL) serves all three targets. Only
the acceptance test changes:

- **Uniform on ALL** includes repaired DSBs (`C₁`) and disconnected AMGs. Its
  cycle-structure law is theorem 2, and its cycle-count law is Ewens(θ = ½)
  (exact.md §7.3).
- **Uniform on DERANGED** is ALL conditioned on "no `C₁`". Its cycle-structure
  law is theorem 2 restricted to partitions without 1s, the same for every Θ.
- **Uniform on PROPER(Θ)** also conditions on connectivity. Its law depends on
  Θ and comes from the recursion (exact.md §5) or from enumeration.

## 6. Exact references for validation

| Level | Reference |
|---|---|
| individual states | enumeration (`iter_states`), for spaces up to a few hundred states |
| cycle structure | theorem 2 (ALL, DERANGED, any n); recursion (PROPER; moderate n and k); enumeration (n ≤ 7) |
| scalar summaries (number of cycles, largest cycle, number of `C_ℓ`), ALL and DERANGED, any n | `exact/summaries.py` (exact.md §8); Ewens(½) mean and variance for ALL (exact.md §7.3) |
| acceptance rate | `theoretical_acceptance` |

Statistical checks use Pearson's χ² goodness-of-fit, with categories of
expected count < 5 pooled. The upper tail is computed from the regularised
incomplete gamma function in `tests/stats.py`, which is itself tested against
known quantiles and a negative control. Each check uses a fixed seed and
significance level `α = 10⁻³`. The tests are deterministic, and for a random
seed each would fail with probability about `α` if the sampler were correct.

## 7. Why IID rejection is preferable to MCMC for uniform targets

When `acceptance(X)` is not tiny, rejection sampling gives exact, independent
draws with a cost of `1/acceptance` proposals each, with no burn-in, no
autocorrelation, and no need to show that a chain is irreducible or mixes
quickly. MCMC is only worth its complications when:

1. acceptance is too small, e.g. when many constraints (patient data) are imposed; or
2. the target is not uniform (biological weights), where plain rejection
   from uniform proposals no longer gives the right distribution.

For the uniform targets studied here, the IID sampler is the baseline that
any MCMC sampler must reproduce.

## 8. Results at n = 80

**Validation** (`tests/test_iid_n80.py`, 8 000 samples each, α = 10⁻³):
- **Θ(1,(80)), PROPER (= DERANGED):** χ² tests of the number of cycles, the
  largest cycle and the number of `C₂` against the exact distributions of §6
  all pass.
- **ALL:** the same tests plus the number of `C₁` pass. The sample mean and
  variance of the number of cycles agree with the Ewens(½) values
  `Σ 1/(2i+1)` and `Σ 2i/(2i+1)²` within 4 standard errors. The variance's
  standard error uses the exact fourth moment.

**PROPER vs DERANGED at n = 80** (exact, from the recursion):

| Θ | 1 − \|PROPER\|/\|DERANGED\| |
|---|---|
| Θ(2,(40,40)) | 7·10⁻²⁵ |
| Θ(4,(20,…)) | 1·10⁻¹⁹ |
| Θ(8,(10,…)) | 5·10⁻¹³ |
| Θ(10,(8,…)) | 4·10⁻¹¹ |
| Θ(80,(1,…,1)) | 0.836 (PROPER is the single-cycle states only) |

When every chromosome carries several DSBs, uniform PROPER and uniform
DERANGED are practically the same distribution at n = 80. The difference
matters when many chromosomes carry one or two DSBs.

**Throughput** (`benchmarks/iid_n80.py`, pure Python 3.11, Apple arm64): about 43 µs per
proposal (≈ 23 000/s). Accepted samples per second: ALL ≈ 23 000;
DERANGED/PROPER with several DSBs per chromosome ≈ 10 000–13 500; Θ(80,(1,…,1))
≈ 1 700 (acceptance ≈ 0.099). A proposal's cost is about 55% shuffle and 30%
matching validation, so a mutable, non-validating representation could at
most roughly double throughput. This has not been done.

## 9. Conditional completion: fixed (observed) rejoins

A `PartialMatching` `F` fixes some rejoin edges. Its free ends `V_free`, with
`|V_free| = f`, are the ends not covered by `F`. A **completion** of `F` is a
perfect matching of all `2n` ends that contains every edge of `F`.

**Lemma (PROVED).** Completions of `F` correspond one-to-one to perfect matchings
of `V_free`. *Proof.* Removing the edges of `F` from a completion leaves a
perfect matching of `V_free`. Adding `F` to a perfect matching of `V_free` gives
a completion. The two maps are inverse to each other. ∎ There are therefore
`(f−1)!!` completions (`num_completions`).

**Theorem (PROVED).** `sample_completion` (shuffle `V_free`, pair consecutive
positions, add `F`) returns every completion with probability `1/(f−1)!!`.
*Proof.* §1 applied to the `f` free ends gives a uniform perfect matching of
`V_free`; the lemma transfers uniformity to completions. ∎

No proposal is ever rejected for missing a fixed edge, so this is not
rejection from ALL. Rejection by state-space membership is then applied
exactly as in §2–§3. The accepted states are therefore IID and uniform on

    X_F(Θ) = { r ∈ X(Θ) : F ⊆ r },

with acceptance probability `|X_F(Θ)| / (f−1)!!`. For small `f` the target is
enumerable (`iter_states(theta, space, fixed)`). The tests compare state
frequencies against it and check the enumeration against filtering the
full enumeration, for every layout with n ≤ 4.

## 10. Using exact DERANGED results as a reference for PROPER

For large `n` and `k ≥ 2`, no exact PROPER cycle distribution is available
(exact.md §5 limits), but `|PROPER(Θ)|` often is. The exact scalar summaries
of DERANGED (exact.md §8) then give a reference with a guaranteed error.

**Lemma (PROVED).** Let `P ⊆ D` be finite and non-empty, and let `μ_P`, `μ_D`
be the uniform distributions on them. Then for every event `A`,

    |μ_P(A) − μ_D(A)| ≤ TV(μ_P, μ_D) = 1 − |P|/|D|.

*Proof.* `μ_P(x) − μ_D(x)` equals `1/|P| − 1/|D| ≥ 0` on `P` and `−1/|D|` on
`D∖P`. The total variation is the sum of the positive parts,
`|P|(1/|P| − 1/|D|) = 1 − |P|/|D|`. For any `A`, `μ_P(A) − μ_D(A)` lies between
minus the sum of the negative parts and the sum of the positive parts. Both
sums equal the total variation. ∎

With `P = PROPER(Θ)` and `D = DERANGED`, the report prints the DERANGED
summaries together with this bound. For Θ(5,(20,20,20,20,20)) the bound is
about `8·10⁻²²`. When the bound is large (e.g. Θ(80,(1,…,1)), 0.84) the DERANGED
reference is still printed with its bound, but it is not informative.
