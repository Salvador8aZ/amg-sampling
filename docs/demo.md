# Demonstration script (5–10 minutes)

Before the meeting: run `uv sync --extra plots`, import the data
(`docs/data.md`), and run each command once so nothing is downloaded live.

## 1. The question (1 min)

> "Cancer genomes contain rearrangements: chromosomes break (double-strand
> breaks) and the broken ends are rejoined, sometimes to the wrong partner.
> Sheth, Arsuaga and Sazdanovic model this as an aberration multigraph: the
> breaks, plus a pairing of the broken ends. Following break → rejoin → break
> gives exchange cycles, and the cycle structure describes the rearrangement
> pathway. Sequencing shows us some of the rejoins but not all. So we ask:
> which pathways are compatible with what we observe, and how unusual is it
> compared with a model where every compatible rearrangement is equally likely?"

## 2. Why exact enumeration fails, and why sampling works (1–2 min)

- With `n` breaks there are `(2n−1)!!` possible pairings of the `2n` ends:
  945 for n = 5, 6.5·10⁸ for n = 10, about 10¹⁸⁶ for n = 100. Listing them
  stops working around n = 8.
- For the uniform question we don't need to list them. Shuffling the ends and
  pairing neighbours gives an **exactly uniform** random pairing (each pairing
  arises from exactly `2ⁿn!` orderings). We keep the draws that satisfy the
  constraints and discard the rest. The kept draws are **independent and
  exactly uniform**: no Markov chain, no burn-in, no convergence diagnostics.
- About 60% of draws are accepted, so at n = 100 we still get thousands of
  exact samples per second.

## 3. Synthetic run (2 min)

```bash
uv run amg-sampling problem=configuration problem.breaks='[20,20,20,20,20]' \
    statespace=proper sampler=iid sampler.num_samples=100000 analysis=distribution seed=42
```

Point at:
- **Header:** Θ(5,(20,…)), 100 DSBs, 200 free ends, the seed, 100 000 samples
  from about 165 000 proposals.
- **Acceptance rate:** the empirical rate equals the exact value (≈ 0.605). The
  sampler's efficiency is known in advance, not guessed.
- **|PROPER(Θ)| ≈ 4·10¹⁸⁶:** computed exactly by a recursion, even though the
  states cannot be listed.
- **Number of cycles:** sampled vs exact columns agree within the standard
  errors. The exact column comes from formulas (Uniform(DERANGED), which the
  report shows differs from Uniform(PROPER) by at most ~10⁻²¹).
- **Typical pathway:** one long cycle (C100) has probability ≈ 0.15, and the
  mean number of cycles is ≈ 2.8.

**What PROPER means:** no break is simply repaired to its own partner (no
`C1`), and the whole rearrangement is connected. These are the "proper AMGs"
of the paper.

## 4. The patient (3 min)

```bash
uv run amg-sampling patients --patient P05-1657
uv run amg-sampling problem=patient problem.patient_id=P05-1657 problem.chromosomes='[8,12]' \
    statespace=proper sampler=iid analysis=observed seed=42
```

- **Why this patient:** it is the paper's own case study (§7.1, figure 13,
  table 2), and the data file is identical to the one the authors used (the
  checksum is printed).
- **Components:** junctions link chromosomes 8 and 12, and separately
  chromosomes 7, 4 and 21. Following the paper, each component is analysed on
  its own.
- **Chromosomes 8 + 12:** 10 breaks, 5 observed rejoins, 10 unmatched ends, so
  **945 possible completions**. Their cycle structures reproduce the paper's
  table 2 exactly (C10 = 384, C8+C2 = 240, …, 5C2 = 1). The software computes
  this both by exact enumeration and by uniform conditional sampling.
- **Share vs null columns:** the table compares the pathways compatible with
  the data against the uniform null for the same Θ. Example: C10 is 40.6% of
  the compatible AMGs but 48.0% of all proper AMGs for Θ(2,(9,1)). The
  observed junctions shift weight towards pathways with small 2-cycles.
- **Per-edge probabilities:** each observed junction has null probability
  ≈ 0.055, just above 1/19 (any given pair under uniform pairing).
- **All five together:** exact probability 945 / 387 099 936 ≈ 2.4·10⁻⁶. With
  100 000 samples we expect about 0.24 hits. The report shows the count with a
  confidence interval, or with an upper bound when there are 0 hits, and never
  reports "probability 0".

**What the patient probability means:** it is a statement about the
**combinatorial null**: "if every proper AMG for this Θ were equally likely,
this pattern would appear with probability p". It is **not** the biological
probability that the cell rearranged this way.

Show `outputs/<date>/<time>/observed_cycle_structures.png` and
`rejoin_probabilities.png`.

## 5. Limitations and next steps (1–2 min)

- **Main limitation:** the null is uniform. Real repair favours nearby ends
  and has biological structure. The numbers are a reference point, not a model
  of the cell. Also:
  - the strand → end convention follows the original code;
  - exact PROPER results are limited to about 12 chromosomes (an `O(3^k)`
    recursion), after which only sampling is used.
- **Why MCMC is future work, not needed now:** for the uniform target,
  independent rejection sampling is already exact and fast, even at 80–100
  breaks, with acceptance ≥ 0.1 in every case tested. MCMC becomes necessary
  once the target is **non-uniform** (weights for biologically plausible
  rejoins) or heavily constrained. Rejection from uniform proposals is then
  inefficient or no longer correct. This IID sampler will then be the baseline
  that any MCMC must reproduce.
- **Next:**
  - the paper's reversal moves;
  - a weighted rejoining model;
  - MCMC and ABC-SMC for inferring model parameters from patients.
