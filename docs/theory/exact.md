# Exact counts and cycle-structure distributions

Covers `amg_sampling.exact` (enumeration and formulas). Labels as in
[`representation.md`](representation.md): **PROVED**, **VERIFIED(n ≤ N)**,
**CITED**, **CONJECTURE**. "Paper" is Sheth, Arsuaga & Sazdanovic (2026).

Throughout, Θ = Θ(k, (b₁, …, b_k)) with `n = Σ bᵢ` DSBs, `d` the DSB matching,
`r` the rejoin matching.

## 1. Three probability spaces

The uniform distribution on each state space induces a cycle-structure
distribution:

    P(C | S, Θ) = #{r ∈ S(Θ) : C(r) = C} / |S(Θ)|,   S ∈ {ALL, DERANGED, PROPER}

| | Set of matchings | Depends on | Cycle-structure distribution |
|---|---|---|---|
| ALL | all perfect matchings of the 2n ends | `n` only | `κ_C(n) / (2n−1)!!` over every partition `C ⊢ n` |
| DERANGED | ALL with no `r(v) = d(v)` | `n` only | `κ_C(n) / κ'(n)` over partitions with no part 1 |
| PROPER | DERANGED ∩ connected | `n` **and** Θ | no closed form in the paper; computed by §5 |

**ALL and DERANGED do not depend on Θ (PROVED).** Both are defined using only
`d` and `r`, and `d(v) = v XOR 1` for every Θ with `n` DSBs. So for fixed `n`
every Θ has literally the same set of states. The cycle structure also
depends only on `(d, r)` (representation.md §3). Hence the distributions agree
for every layout. The test suite also checks this for the layouts of n = 4 and
n = 5 by comparing state sets and distributions (VERIFIED(n ≤ 5), in addition
to the proof).

**DERANGED is ALL conditioned on "no C₁" (PROVED).** DERANGED = {r ∈ ALL :
C(r) has no part 1}, so

    P(C | DERANGED) = P(C | ALL) / P(no C₁ | ALL)   for C with no part 1, else 0.

**PROPER depends on Θ (VERIFIED, counterexamples).** Θ(1,(4)), Θ(2,(2,2)) and
Θ(4,(1,1,1,1)) all have n = 4, but give P(2C₂ | PROPER) = 12/60, 8/56 and 0.
Connectivity is a property of how rejoins link *chromosomes*, so it cannot be
read from the cycle structure alone.

### What probability space does theorem 2 describe?

Theorem 2 counts perfect matchings of the `2n` free ends by cycle structure,
with no connectivity condition and (in the `m_l ≥ 0` form) with `C₁` allowed.
It therefore describes **ALL**: `κ_C(n)/(2n−1)!!` is `P(C | ALL)`, the same for
every Θ with `n` DSBs. Restricted to partitions without 1s (paper eq. 3) it
describes **DERANGED**, and for Θ(1,(n)) only, where DERANGED = PROPER, it
also describes PROPER.

For any Θ with `k ≥ 2`, theorem 2 does **not** give `P(C | PROPER, Θ)`.
What it does give is that the count of states with cycle structure C is at
most `κ_C(n)`. For `C = C_n` the count equals `κ_C(n)` (§5, fact P3).

Consequence for sampling: a sampler targeting uniform PROPER(Θ) must be
validated against `P(C | PROPER, Θ)` (§5 or enumeration), not against theorem 2.

## 2. |ALL|

`|ALL| = (2n−1)!! = ∏_{i=1}^n (2i−1)` (CITED theorem 1; PROVED by the
canonical decomposition below; VERIFIED(n ≤ 7)).

## 3. |DERANGED|

`|DERANGED| = κ'(n)` = number of perfect matchings of `2n` points that avoid a
fixed perfect matching (OEIS A053871: 1, 0, 2, 8, 60, 544, 6040, 79008, …).
The paper derives it for Θ(1,(n)) (theorem 5: the recurrence
`κ'(n) = 2(n−1)(κ'(n−1)+κ'(n−2))`, the closed form eq. 5, and the EGF
`e^{−x}(1−2x)^{−1/2}`). Because DERANGED does not depend on Θ (§1), the same
numbers hold for every Θ (PROVED). The recurrence, the closed form, and the sum
of theorem 2 over partitions without 1s agree for `n < 40`. They also match
enumeration on every layout with n ≤ 6, and on five layouts with n = 7
(VERIFIED).

## 4. |PROPER| and the counting problem

`PROPER(Θ)` = deranged matchings whose chromosome graph is connected
(representation.md §4). This does depend on Θ, and the paper gives formulas
only for special cases:

| Paper result | Setting | Status here |
|---|---|---|
| Theorem 5 / eq. 3 | Θ(1,(n)) | PROVED equal to DERANGED; VERIFIED |
| Lemma 10: `κ'(n) − κ'(l)κ'(n−l)` | Θ(2,(l,n−l)), stated for `⌊n/2⌋ ≤ l ≤ n−2` | Special case of §5; checked for **all** `1 ≤ l ≤ n−1`, n ≤ 15 |
| Lemma 12 | Θ(3,(ℓ₁,ℓ₂,ℓ₃)), stated for `ℓᵢ ≥ 2` | Special case of §5; checked for all `ℓᵢ ≥ 1`, n ≤ 12 |
| Lemma 14: `2^{n−1}(n−1)!`, all `C_n` | Θ(n,(1,…,1)) | Special case of §5; checked n ≤ 10 |
| Lemma 15, bullet 1: `2^{n−1}(n−1)!` with `C_n` | Θ(n−1,(2,1,…,1)) | agrees (n ≤ 10) |
| Lemma 15 (ii): `2^{n−1}(n−1)!(n−2)` two-cycle AMGs, `C_j + C_{n−j+1}` | same | Uses `n` for `n+1` DSBs; restated with `N` DSBs as `2^{N−2}(N−2)!(N−3)`. See §7.1 |

## 5. Exact PROPER distribution without enumeration (derived)

For a set `S` of chromosomes, let `n_S = Σ_{i∈S} bᵢ`. Let `D_S` be the
cycle-structure distribution of deranged matchings on the ends of `S` (equal
to theorem 2 without 1-parts for `n_S`, by §1). Let `P_S` be the same
distribution restricted to matchings whose chromosome graph on `S` is
connected.

**Recursion (PROVED).** Let `s = min S`. Every deranged matching on `S`
decomposes uniquely as follows. `T` is the set of chromosomes in the connected
component of `s`, so `s ∈ T ⊆ S`. The matching restricted to `T` is deranged
and connected. The matching restricted to `S∖T` is deranged and has no rejoin
into `T`. Each matching with this decomposition has cycle structure `C_T ⊎ C_{S∖T}`,
because exchange cycles never cross between `T` and `S∖T`. Conversely, every
such pair gives a deranged matching on `S` whose component of `s` is `T`.
Hence

    D_S = Σ_{s ∈ T ⊆ S} P_T ⊛ D_{S∖T},     D_∅ = {∅ : 1}

where `⊛` multiplies counts and unites parts. The `T = S` term is `P_S`, so
`P_S` is obtained by subtracting the other terms. `cycle_distribution_proper`
implements this, and `num_proper_states` does the same with scalar counts.
The number of subset pairs is `O(3^k)`. For distributions, the cost also grows
with the number of partitions of the block sizes.

Specialising to `k = 2` gives lemma 10. For `k = 3`, expanding the
recursion gives lemma 12, i.e. Möbius inversion over set partitions with
`μ = (−1)^{|π|−1}(|π|−1)!` (PROVED by expansion; also checked numerically).

VERIFIED: agrees with enumeration on every layout with n ≤ 6, on five layouts
with n = 7, and on every row of tables 1 and 3.

**Facts derived from the recursion or definitions:**

- **P1 (PROVED).** `|PROPER(Θ)| ≤ |DERANGED| = κ'(n)`. Equality holds iff no
  deranged matching splits the chromosomes into two groups with no rejoin
  between them. Examples: `k = 1`, and Θ(2,(n−1,1)) because `κ'(1) = 0`
  (table 3, row (5,1)).
- **P2 (PROVED).** A one-DSB chromosome can never be a connected component on
  its own, because its two ends would form a `C₁`.
- **P3 (PROVED).** For `n ≥ 2`, every state of ALL with cycle structure `C_n`
  is PROPER, for every Θ. *Proof:* a single cycle through all `n` DSBs has no
  `C₁`. Its rejoin edges link consecutive DSBs of the cycle, so all
  chromosomes lie in one component of the chromosome graph. Hence
  `#{C_n in PROPER(Θ)} = 2^{n−1}(n−1)!` for every Θ. This is why the `C_n`
  column of tables 1 and 3 is constant.

## 6. How cycle-structure counts are obtained

1. **Enumeration** (`iter_states`, `exact_cycle_distribution`). The
   **canonical decomposition** pairs the smallest unmatched end `v` with each
   larger unmatched `w` and recurses. Each perfect matching is produced
   exactly once, since the partner of the smallest unmatched end identifies
   the branch at each level (PROVED; VERIFIED against brute force as
   *sets*, every layout n ≤ 6). For DERANGED and PROPER the branch
   `w = d(v)` is skipped, which is safe because such matchings are not
   deranged. Membership is then checked explicitly. Integer counts are the
   authoritative output; `probabilities()` gives `Fraction`s.
2. **Formulas** (`exact/formulas.py`). Theorem 2 for ALL and DERANGED, and §5
   for PROPER. Integer and `Fraction` arithmetic only.

## 7. Differences from the paper found in Stage B

Each item says exactly what disagrees and what evidence supports the
expression used here. Where both a derivation and independent enumeration
agree, the code uses the derived expression. The printed one is kept only
inside the tests.

### 7.1 Lemma 15 (ii): use of `n`

**PAPER STATEMENT.** The heading says Θ(n−1, (2,1,…,1)) is "n−2 chromosomes with 1 DSB each and 1
chromosome with 2 DSBs", i.e. `n` DSBs. (i): `2^{n−1}(n−1)!` AMGs with `C_n`.
(ii): `2^{n−1}(n−1)!(n−2)` AMGs with two cycles `C_j + C_{n−j+1}`, `2 ≤ j ≤ n−1`.
The proof of (ii) selects "i of the n−1 chromosomes with a single DSB"
(`1 ≤ i ≤ n−2`) and forms cycles of lengths `2(i+1)` and `2(n−i)`.

**WHAT IT COUNTS.** PROPER AMGs of that Θ with exactly two exchange cycles.

**OBSERVATION.** The heading and (i) use `n` = number of DSBs. Bullet (ii) and
its proof describe one 2-DSB chromosome plus `n−1` single-DSB chromosomes,
i.e. `n+1` DSBs. The two cycle lengths also sum to `2(n+1)` edges.

**DERIVATION** (with `N` = number of DSBs; 2-DSB chromosome with DSBs α, β;
single-DSB chromosomes `s_1…s_{N−2}`).
1. α and β lie in different cycles. Otherwise the other cycle holds only
   single-DSB chromosomes, whose ends all rejoin inside that cycle, so they
   are disconnected from chromosome A.
2. Each cycle holds at least one single-DSB chromosome. A cycle through α
   (or β) alone is a `C₁`.
3. Conversely, any two-cycle matching satisfying 1 and 2 is deranged (both
   cycles have length ≥ 2) and connected (both cycles contain chromosome A).
4. Choose the set `I` of single DSBs in α's cycle, `|I| = i`, `1 ≤ i ≤ N−3`.
   By corollary 3 there are `2^{l−1}(l−1)!` single cycles through `l` given
   DSBs, so the count is `C(N−2,i) · 2^i i! · 2^{N−2−i}(N−2−i)! = 2^{N−2}(N−2)!`
   for each `i`. Total: **`2^{N−2}(N−2)!(N−3)`**, structures `C_j + C_{N−j}` with `2 ≤ j ≤ N−2`.

This equals the printed expression with `n = N−1`: printed`(n)` = restated`(n+1)`.

**COMPUTATIONAL VERIFICATION.**

| N | printed at n = N | restated | production enumerator | brute-force reference | reference repository |
|---|---|---|---|---|---|
| 3 | 8 | 0 | 0 | 0 | 0 |
| 4 | 96 | 8 | 8 | 8 | 8 |
| 5 | 1152 | 96 | 96 | 96 | 96 |
| 6 | 15360 | 1152 | 1152 | 1152 | 1152 |
| 7 | 230400 | 15360 | 15360 | — | 15360 |

The printed value at `n = N` also differs from the two-cycle counts in ALL
and DERANGED (N = 4: 44 and 12). The ratio printed/actual, `2(N−1)(N−2)/(N−3)`,
is not constant, so no fixed symmetry or ordering factor explains it. Tests:
`test_lemma_15_*`; §5 recursion checked for `N ≤ 10`.

**CONCLUSION.** This is an indexing inconsistency within the lemma, not an
arithmetic error. After restating everything with `N` = number of DSBs, the
derivation, the paper's own proof and enumeration agree. The API exposes only
`lemma15_two_cycle_count(num_dsbs)`.

### 7.2 Section 3.1: number of AMGs with `c` cycles

**PAPER STATEMENT.** "In the case of AMGs on a single chromosome with n DSBs the number of AMGs with
exactly c cycles equals the Stirling number of the first kind" `[n c]`. The
section also says that the probability of `c = k` is `[n k]/n!` and that a typical
AMG has about `log n` cycles.

**DERIVATION** (a bijection).
- Take an exchange cycle through a set `B` of `l` DSBs and let `j₀ = min B`.
- Walk it from end `2j₀`, crossing DSB `j₀` first. Then follow `r` into DSB
  `j₁`, cross it, follow `r` into `j₂`, and so on back to `2j₀`.
- This records a cyclic order `(j₀, j₁, …, j_{l−1})` of `B` and, for each
  `t ≥ 1`, a bit saying which end of DSB `j_t` was entered.
- Starting at the minimum DSB, through its left end, fixes both the rotation
  and the direction. So each cycle yields exactly one (cyclic order, bits) pair.
- Conversely, every cyclic order with `l−1` bits rebuilds a unique set of
  rejoin edges.
- Hence single cycles on `B` correspond to `(l−1)!` cyclic permutations ×
  `2^{l−1}` orientations. For `l = 1` there are no bits and the result is a
  repaired DSB (`C₁`).
- Over all cycles, ALL is in bijection with pairs `(σ, ε)`: σ a permutation of
  the `n` DSBs with `c(σ)` cycles, and `ε ∈ {0,1}^{n−c(σ)}`. Therefore

      #{r ∈ ALL : c cycles} = 2^{n−c} [n c].

The factor `2^{n−c}` is one free orientation per DSB, except the reference DSB
of each cycle, whose orientation fixes the direction of traversal.

**CONSISTENCY WITH THE PAPER'S THEOREMS.**
- Grouping theorem 2 by `|C|` gives `Σ_{|C|=c} 2^{n−c} n!/∏ l^{m_l} m_l!`,
  which is `2^{n−c}` times the number of permutations with `c` cycles.
- The proof of theorem 2 itself mentions "2 orientations" per initial edge
  and "one rejoin edge per cycle as a reference".
- Theorem 1 follows: `Σ_c 2^{n−c}[n c] = 2^n ∏_{i=0}^{n−1}(1/2 + i) = ∏_{i=0}^{n−1}(2i+1) = (2n−1)!!`.

**COMPUTATIONAL VERIFICATION.**
- `cycle_count_all(n, c)` equals enumeration of ALL for every `c`, n ≤ 6.
- The printed `[n c]` equals neither ALL nor DERANGED (= PROPER for `k = 1`)
  for `n ≥ 2, c < n`. Example: n = 4, c = 1 gives `[4 1] = 6`, but ALL = DERANGED = 48.
- Counting ALL up to swapping the two ends of each DSB also does not give
  `[n c]` (n = 4, c = 1: 3 orbits).
- Normalisation to `(2n−1)!!` checked for n ≤ 40. The three-level
  consistency `κ_C → 2^{n−c}[n c] → (2n−1)!!` checked for n ≤ 20.

**CONCLUSION.** The count of AMGs with `c` cycles implied by theorems 1 and 2
and by enumeration is `2^{n−c}[n c]`. Section 3.1 prints `[n c]`, which counts
the permutation σ alone, without orientations. The library implements
`cycle_count_all(n, c) = 2^{n−c}[n c]`. The printed expression appears only in tests.

### 7.3 DERIVED: the cycle count under uniform ALL is Ewens with θ = 1/2

Not used by any code; recorded as a derived result.

From the rising-factorial identity `Σ_c [n c] y^c = y(y+1)⋯(y+n−1)` with `y = x/2`:

    G(x) := Σ_c 2^{n−c}[n c] x^c = 2^n ∏_{i=0}^{n−1}(x/2 + i) = ∏_{i=0}^{n−1}(x + 2i).

- **Normalisation:** `G(1) = ∏(2i+1) = (2n−1)!! = |ALL|` (theorem 1).
- **Distribution:** dividing numerator and denominator by `2^n`,

      P(c | ALL) = 2^{n−c}[n c]/(2n−1)!! = θ^c [n c] / (θ(θ+1)⋯(θ+n−1)),  θ = 1/2,

  which is the Ewens distribution of the number of cycles with θ = 1/2.
- **Independent-sum form:** `G(x)/G(1) = ∏_i (x + 2i)/(1 + 2i)`, so the number
  of cycles is a sum of independent Bernoulli variables with success
  probabilities `1/(2i+1)`, `i = 0,…,n−1`.
- **Mean:**

      E[c | ALL] = Σ_{i=0}^{n−1} 1/(2i+1) = H_{2n} − H_n/2 ~ (1/2) log n + log 2 + γ/2.

**COMPUTATIONAL VERIFICATION.** The exact mean equals the enumerated mean as a
`Fraction` for n ≤ 6 (e.g. n = 6: 1.8782… vs `H_6` = 2.45). It also equals the
mean from `cycle_count_all` for n ≤ 30.

**CONCLUSION (derived).** Under uniform ALL the mean number of cycles grows
like `(1/2) log n`, not `log n`. The paper's section 3.1 heuristics correspond
to θ = 1, i.e. uniform random permutations. This concerns ALL only. DERANGED
and PROPER have different cycle-count distributions (no `C₁`), and no
analogous closed form is claimed here.

### 7.4 Minor points

1. **Corollary 7 remark.** "The number of cycle structures `C_l + C_{n−l}` with
   `l = n−l` is two times bigger than `2C_l`" is unclear. The corollary's two
   formulas themselves agree with theorem 2 (checked n ≤ 21).
2. **Lemmas 10 and 12 ranges.** Their formulas also hold outside the stated
   ranges, for all block sizes ≥ 1 (they are special cases of §5; checked
   n ≤ 15 and n ≤ 12).
3. **Theorem 5 initial conditions.** Stated "for n ⩾ 1, with κ'(1)=0 and κ'(2)=2".
   The recurrence applies from `n ≥ 3`, or from `n ≥ 2` with `κ'(0) = 1`.
4. **Lemma 14** writes the cycle structure as "n1". Context indicates `C_n`.
5. **Proposition 4** asserts existence of proper AMGs with `c` cycles for
   `c ≤ min(max bⱼ, ⌊n/2⌋)`. It is not an upper bound: Θ(4,(2,2,1,1)) has
   proper `3C₂` states. The existence claim is VERIFIED for every layout n ≤ 6.

## 8. Summary: what depends on what

| Quantity | Depends only on n | Depends on Θ |
|---|---|---|
| ALL, its size, its cycle-structure distribution | ✓ | |
| DERANGED, its size, its distribution | ✓ | |
| number of `C_n` states in PROPER | ✓ (P3) | |
| PROPER, its size, its distribution | | ✓ |
| which states are PROPER | | ✓ (chromosome membership only, not DSB order within a chromosome) |

The last row is PROVED: connectivity uses only `chrom(v)`. Reordering the
chromosomes also changes nothing (PROVED): the relabelling that moves each
chromosome's block of DSB ends to its new position preserves `d` and maps
chromosome membership consistently. So conjugating by it is a bijection
between the two PROPER sets that preserves cycle structure. A property test
checks this as well.
