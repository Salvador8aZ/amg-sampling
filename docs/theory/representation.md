# Mathematical representation

This note covers only what is implemented in `amg_sampling.core`. Every claim
carries one of these labels:

- **PROVED**: an argument is given here.
- **VERIFIED(n ≤ N)**: checked exhaustively by the test suite up to `n = N`.
  This is evidence, not a proof.
- **CITED**: stated in Sheth, Arsuaga & Sazdanovic (2026), not re-derived here.
- **CONJECTURE / UNVERIFIED**: not established; not used by any code.

Paper references are to Sheth, Arsuaga, Sazdanovic, *J. Phys. A* **59** 115601 (2026).

## 1. Initial configuration Θ (`core/configuration.py`)

Θ(k, (b₁, …, b_k)) has `k` chromosomes, chromosome `i` carrying `bᵢ ≥ 1` DSBs,
and `n = Σ bᵢ` DSBs in total.

**Labelling.** DSBs are numbered `j = 0, …, n−1` consecutively along chromosome
0, then chromosome 1, etc. DSB `j` has free ends `2j` (left) and `2j+1` (right).
The set of free ends is `V_F = {0, …, 2n−1}`.

Θ is stored as three fixed maps on `V_F`:

| Map | Definition | Kind |
|---|---|---|
| DSB partner `d` | `d(v) = v XOR 1` | fixed-point-free involution |
| chromatin partner `c` | `c(2j+1) = 2j+2` and `c(2j+2) = 2j+1` if DSBs `j`, `j+1` are on the same chromosome; otherwise `c(v) = TELOMERE` | partial involution |
| chromosome membership | `chrom(v) = chromosome of DSB ⌊v/2⌋` | |

Each chromosome `i` is the path

    T_i^L —c— 2f —d— 2f+1 —c— 2f+2 —d— 2f+3 —c— … —d— 2l+1 —c— T_i^R

where `f` and `l` are its first and last DSB. The internal chromatin edges
`(2j+1, 2j+2)` are called **fragments**. There are `n − k` fragments and `2k`
telomeres. Telomeres are recorded explicitly (`Telomere(chromosome, side, end)`)
but are not free ends and carry no integer label.

Facts (**PROVED**, immediate from the construction; also property-tested):
each free end has exactly one chromatin neighbour, either one other free end
or one telomere; `c(v) ≠ d(v)`; `c` never crosses chromosomes.

**Relation to other labellings.** The paper labels vertices with letters,
telomeres included (e.g. Θ(1,(2)): `a b c d e f`, free ends `b c d e` = `0 1 2 3`).
The reference repository numbers telomeres and free ends together along
each chromosome. The PennyLane script in the reference repository uses DSB
ends `(2i, 2i+1)` like this project, but omits telomeres and the chromatin
edges to them.

## 2. State: rejoin matching (`core/matching.py`)

A state is a rejoin matching `r`: a fixed-point-free involution of `V_F`
(`r(r(v)) = v`, `r(v) ≠ v`), stored as an immutable tuple `r[v]`. This is the
paper's rejoin map ρ with rejoin edges `E_R = {(v, ρ(v))}`.

An AMG is the pair (Θ, r). The chromatin and DSB edges belong to Θ and are
never stored per state.

**Why partner maps, not a graph.** The three relations `d`, `c`, `r` stay
separate, so a rejoin edge parallel to a DSB edge (repaired DSB) or to a
chromatin edge (ring) never overwrites anything. The reference repository
stored all edges in a simple `networkx.Graph`, which merges parallel edges; as
a result it lost `C₁` cycles. NetworkX appears in this project only as an
independent reference in tests (a `MultiGraph` keyed by edge type).

## 3. Cycle structure (`core/cycles.py`)

The exchange multigraph Ξ has vertex set `V_F` and edge multiset
`{DSB edges} ⊔ {rejoin edges}` (paper §2).

**Fact (PROVED).** Ξ is a disjoint union of even cycles that alternate DSB and
rejoin edges. *Proof:* `d` and `r` are perfect matchings of `V_F`, so every
vertex has degree exactly 2 (counting parallel edges), and the two edges at a
vertex are of different types.

A cycle through `l` DSB edges has `2l` edges and is written `C_l`. The cycle
structure `C(Ω) = Σ m_l C_l` is a partition of `n` (`Σ l·m_l = n`, CITED §2).
It is stored as `CycleStructure(parts)`, the parts in non-increasing order,
i.e. the rows of its Young diagram. The string form lists cycles by decreasing
length with multiplicities: `C4`, `2C2`, `C3+C1`. The paper is not consistent
about the order (`C2+C3` in table 3, `C8+C2` in table 2); only the partition matters.

**C₁.** `r(2j) = 2j+1` means DSB `j` is correctly repaired. Its DSB and rejoin
edges are parallel and form a `C₁`. `repaired_dsbs(Θ, r)` lists these DSBs.

**Rings.** `r(a) = b` for a fragment `(a, b)` means the rejoin edge is parallel
to a chromatin edge (paper figure 3(b)). `ring_fragments(Θ, r)` lists them.
Chromatin edges are not in Ξ, so rings do not affect the cycle structure.

**Computation.** Starting from any unvisited end `v`, follow `v → d(v) → r(d(v)) → …`
until returning to `v`; the number of `d`-steps is `l`. Time `O(n)`.

**Fact (PROVED).** The cycle structure depends only on `(d, r)`, not on the
chromosome layout (Ξ contains no chromatin edges). Hence it is also
invariant under any relabelling that conjugates both `d` and `r`
(an isomorphism of Ξ). Both are property-tested.

## 4. State spaces (`core/statespace.py`)

For fixed Θ:

| Space | Definition | Size |
|---|---|---|
| `ALL` | every rejoin matching of `V_F` | `(2n−1)!!` (CITED thm 1; VERIFIED(n ≤ 5)) |
| `DERANGED` | `r(v) ≠ d(v)` for all `v`, i.e. no `C₁` | OEIS A053871 (VERIFIED(n ≤ 5)); independent of the chromosome layout (PROVED: the definition involves only `d` and `r`) |
| `PROPER` | deranged and the AMG is connected | tables 1 and 3 (VERIFIED for every row, n ≤ 6) |

`PROPER` follows the paper (§3): "connected, and its exchange multigraph does
not contain parallel edges". The only possible parallel edges in Ξ are a
rejoin edge parallel to a DSB edge, so "no parallel edges" = `DERANGED`.
Rings are allowed: example 2 counts the ring AMG as proper.

**Connectivity criterion (PROVED).** Let `H` be the multigraph whose vertices
are the chromosomes, with one edge `{chrom(v), chrom(r(v))}` per rejoin edge
between different chromosomes. The AMG is connected iff `H` is connected.
*Proof:* each chromosome's telomeres and free ends are connected by the
path in §1, which uses only chromatin and DSB edges. So contracting each
chromosome to a point preserves connectivity. After contraction the only
remaining edges are rejoin edges between different chromosomes, which are
exactly the edges of `H`. ∎ Evaluated by union-find in `O(n α(k))`.
It also agrees with `networkx.is_connected` on the explicit AMG multigraph
for every DSB distribution and every matching with n ≤ 6, and on random
states (Hypothesis).

Consequence (PROVED): for `k = 1`, `PROPER = DERANGED`.

### PROPER is not necessarily the patient-level state space

`PROPER` implements the paper's combinatorial definition. The state space for
patient inference is a separate modelling decision and may differ, e.g. by
fixing observed rejoins, by allowing some repaired DSBs, or by relaxing
connectivity when only a subset of chromosomes is analysed. Nothing in the
core assumes that sampling will happen on `PROPER`. Each future sampler must
name its state space explicitly.

## 5. Claims deferred to later stages (UNVERIFIED; not used by any code)

These appeared in planning discussions. They must be derived, checked
exhaustively for small `n`, and labelled before any code depends on them.

1. Every state has exactly `n(n−1)` distinct 2-switch neighbours.
2. Every 2-switch neighbour arises from exactly four ordered endpoint pairs `(u, v)`.
3. The resulting switch proposal kernel is symmetric on `ALL`.
4. The switch graph restricted to `DERANGED` or `PROPER` is connected. (CONJECTURE)
5. The set of ring fragments is invariant under reversals (an argument was
   sketched during planning; it will be written up with the reversal move).
6. IID rejection acceptance rates: `|DERANGED| / |ALL| → e^{−1/2}` for large `n`,
   and `|PROPER| / |ALL| ≈ √(π/(4n))` for Θ(n, (1, …, 1)). Still unproved; exact
   values are computed instead (`theoretical_acceptance`), e.g. 0.6046 and
   0.0992 at n = 80.
7. Uniform IID sampling from `ALL` by shuffling the `2n` ends and pairing
   consecutive ones. **Now PROVED** in [`sampling.md`](sampling.md) §1.
