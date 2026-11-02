# amg-sampling

Exact combinatorics and uniform sampling of **aberration multigraphs (AMGs)**,
the graph model of chromosome rearrangements used by

> S. Sheth, J. Arsuaga, R. Sazdanovic (2026). *Characterizing cancer chromosome
> aberration pathways using multigraphs.* J. Phys. A: Math. Theor. **59** 115601.
> https://doi.org/10.1088/1751-8121/ae4d8f

The long-term goal is to **enumerate AMGs with quantum computing**. The
project gets there in three parts, built on one shared exact-combinatorics
core:

| Part | Role | Entry point | Status |
|---|---|---|---|
| **I. IID sampling** | The base case: exactly uniform, independent AMGs (the combinatorial null), checked against exact counts | `sampler=iid` | Implemented; validated up to n = 80 |
| **II. MCMC with proximity weights** | The proposed sampler when biological restrictions make the target non-uniform and rejection sampling stops working | `sampler=mcmc target=proximity` | Implemented; validated against exact weighted enumeration |
| **III. Quantum annealing** | Counting PROPER AMGs as a Montanaro (2015) partition-function problem: the route to a quantum speedup | `scripts/quantum_annealing.py` (on the `quantum` branch) | Exact small-n emulation |

Each part is checked against the one before it: the exact results check IID,
IID checks MCMC at α = 0, and the exact energy histogram checks the
annealing schedule.

The shared core reproduces the paper's published counts and extends them from
small cases, where every configuration can be listed, to realistic sizes
(80+ double-strand breaks) by exact formulas. It also turns a patient's
observed rearrangements into Θ plus observed rejoins, so every part can be run
on real data.

## The scientific question

- A **double-strand break (DSB)** cuts a chromosome and leaves two **free ends**.
- With `n` DSBs on `k` chromosomes there are `2n` free ends. Repair pairs them
  up: a **rejoin matching**. An end rejoined to its original partner is a
  correct repair; any other pairing is an exchange.
- The **initial configuration** `Θ(k, (b₁,…,b_k))` records how many DSBs each
  chromosome has. An **AMG** is Θ together with one rejoin matching.
- Following DSB → rejoin → DSB → … around an AMG gives closed loops, its
  **exchange cycles**. The multiset of their lengths is the **cycle structure**,
  e.g. `C10` (one cycle through all 10 DSBs) or `C8+C2`. The cycle structure
  describes the exchange pathway: many small cycles suggest independent
  events, and one long cycle a coordinated one (e.g. chromoplexy).
- A **proper** AMG (paper §3) has no correctly repaired DSB and is connected.

Sequencing observes some rejoins (junctions) but usually not all. The question
is then: *which exchange pathways are compatible with the observations, and
how unusual is what we see compared with a model in which every compatible
AMG is equally likely?*

The number of AMGs is `(2n−1)!!`: about 6.5·10⁸ for n = 10 and 7·10¹⁸⁶ for
n = 100. Listing them all works only for n ≲ 8, so larger problems need
formulas or sampling.

## What the software does

### Shared core

| Capability | How |
|---|---|
| **Exact combinatorial analysis** | Enumerates every AMG for small n, and applies the paper's theorems plus derived formulas and recursions for any n. Reproduces every row of the paper's tables 1–3. |
| **Synthetic configurations** | Any Θ, e.g. five chromosomes with 20 DSBs each. |
| **Sheth patient reconstruction** | Loads the dataset used in the paper, splits each patient into chromosome components, and builds Θ and the observed rejoins with provenance. Reproduces the paper's case study (patient P05-1657: 945 pathways, table 2). |
| **Patient vs null comparison** | Compares the patient's observed rejoins and the cycle structures compatible with them against the uniform combinatorial null. Exact values are used where feasible, sampled estimates with standard errors otherwise. |

### Part I: IID sampling (base case)

Draws uniform random rejoin matchings and keeps those in the chosen state space
(ALL, DERANGED or PROPER). The samples are exactly uniform and independent: no
Markov chain, burn-in or autocorrelation. With observed rejoins, only the
unmatched ends are drawn. Validated against the exact results at the state,
distribution and acceptance-rate levels, including at n = 80
([docs/theory/sampling.md](docs/theory/sampling.md)).

### Part II: MCMC with proximity weights (biological restrictions)

Uniform rejection sampling is exact, but it cannot sample a **non-uniform**
target. Biology gives such a target, because repair favours nearby ends. Part
II is a Markov chain whose moves are 2-switches: two junctions re-pair their
ends, the combinatorial analogue of an inversion.

- **Target.** `target=proximity` weights each rejoin by
  `(distance / scale_bp)^(−α)` within a chromosome, and by
  `exp(inter_log_weight)` between chromosomes. `α = 0` gives back the
  uniform null, and the chain then reproduces Part I.
- **Calibration.** α and the inter-chromosomal weight are fitted to the Baca
  et al. junctions on a train/test split of patients
  (`analysis/calibration.py`, `scripts/biological_toy.py`).
- **Validation.** Detailed balance of the exact kernel, and agreement with
  exact weighted enumeration (`exact/weighted.py`). Standard errors use batch
  means and an effective sample size
  ([docs/theory/mcmc.md](docs/theory/mcmc.md)).

### Part III: quantum annealing (towards quantum enumeration)

`|PROPER(Θ)|` is the ground-state count `Z(∞)` of a defect Hamiltonian on all
AMGs. It is estimated by a telescoping product of ratios along a
B-Chebyshev cooling schedule, with each ratio estimated by amplitude
estimation (Montanaro 2015, §3). The walk operators are the Part II 2-switch
chain at each temperature. The code emulates the ideal algorithm exactly for
small Θ: partition function, schedule, overlaps, relaxation times, and
classical-vs-quantum error scaling. It lives on the `quantum` branch
(`amg_sampling/quantum.py`, `scripts/quantum_annealing.py`,
`docs/theory/quantum.md`); `git checkout quantum` to use it.

## Installation

Requires Python 3.11–3.13 (Hydra does not yet run on 3.14) and
[Poetry](https://python-poetry.org). Create an isolated environment first
(venv, Conda or pyenv), then install Poetry into it:

```bash
git clone https://github.com/Salvador8aZ/amg-sampling.git
cd amg-sampling
python3 -m venv .venv
source .venv/bin/activate
pip install poetry
poetry install --extras plots   # runtime + dev dependencies (+ matplotlib for figures)
```

The project defines two dependency groups:

- **dev**: testing and code quality tools (pytest, hypothesis, networkx, ruff,
  coverage, vulture, deptry). Installed by default, so a fresh clone can run
  `make check` immediately.
- **notebook**: Jupyter, JupyterLab, Matplotlib and Seaborn. Marked
  `optional = true`, so it is installed only with `poetry install --with notebook`.

Install runtime dependencies alone, without dev tooling, with
`poetry install --only main`. The `plots` extra is what consumers of the
package install for figures (`pip install "amg_sampling[plots]"`).

The patient data are **not** included; see [docs/data.md](docs/data.md).
They are Table S3C of the supplementary material of Baca et al. (2013),
downloaded from the publisher and converted to `data/external/nihms.csv` on
first use (falling back to the copy distributed by Sheth et al.), or
explicitly with:

```bash
make data   # = poetry run amg-sampling fetch-data
```

## Quick start

```bash
# Shared: run the test suite (about 2 minutes; patient tests skip if the data are absent)
make test

# Shared: list the patients in the Sheth dataset, then inspect the demonstration patient
poetry run amg-sampling patients
poetry run amg-sampling patients --patient P05-1657

# Shared: the demonstration patient, exact results only (no sampling)
poetry run amg-sampling problem=patient analysis=observed sampler=exact
```

**Part I: IID (base case)**

```bash
# A synthetic configuration: Θ(5,(20,20,20,20,20)), 100 DSBs, 200 free ends
poetry run amg-sampling problem=configuration problem.breaks='[20,20,20,20,20]' \
    statespace=proper sampler=iid sampler.num_samples=100000 analysis=distribution seed=42

# Chromosomes 8 and 12 of P05-1657 vs the uniform null
poetry run amg-sampling problem=patient problem.patient_id=P05-1657 problem.chromosomes='[8,12]' \
    statespace=proper sampler=iid analysis=observed seed=42

# Sampling throughput at n = 80
poetry run python benchmarks/iid_n80.py
```

**Part II: MCMC with proximity weights**

```bash
# The 2-switch chain with the uniform target: reproduces Part I
poetry run amg-sampling problem=patient problem.patient_id=P05-1657 problem.chromosomes='[8,12]' \
    statespace=proper sampler=mcmc analysis=observed seed=42

# A proximity-weighted target (alpha = 0 is uniform)
poetry run amg-sampling problem=patient problem.patient_id=P05-1657 problem.chromosomes='[8,12]' \
    statespace=proper sampler=mcmc target=proximity target.alpha=1.0 analysis=observed seed=42

# Fit alpha and the inter-chromosomal weight to the Baca et al. junctions
poetry run python scripts/biological_toy.py
```

**Part III: quantum annealing** (on the `quantum` branch)

```bash
git checkout quantum
poetry install --extras "plots quantum"
poetry run python scripts/quantum_annealing.py
```

The synthetic IID run takes about 15 s and the patient IID run about 3 s on a laptop.
[docs/demo.md](docs/demo.md) is a 5–10 minute walkthrough of these commands.

## Configuration

Experiments are configured with [Hydra](https://hydra.cc). Four independent
groups live in `amg_sampling/conf/`, and any key can be overridden on the
command line.

| Group | Options | Main keys |
|---|---|---|
| `problem` | `configuration` (default), `patient` | `problem.breaks=[…]`; `problem.patient_id`, `problem.chromosomes=[…]` |
| `statespace` | `proper` (default), `deranged`, `all` | — |
| `sampler` | `iid` (default), `mcmc`, `exact` | `sampler.num_samples`; `sampler.burn_in`, `sampler.thin` (mcmc) |
| `analysis` | `distribution` (default), `observed`, `rejoin_probability` | `analysis.confidence` |
| top level | — | `seed`, `data.sheth_csv`, `output.plots`, `output.top` |

- **`iid`**: exactly uniform, independent samples by rejection.
- **`mcmc`**: a Markov chain whose moves are 2-switches (re-pairing the ends
  of two junctions, the combinatorial analogue of an inversion). With observed
  rejoins only the unmatched ends move. For uniform targets it reproduces
  `iid` and exists as the basis for biologically weighted targets. Its
  standard errors use batch means and an effective sample size, printed in the
  report ([docs/theory/mcmc.md](docs/theory/mcmc.md)).
- **`exact`**: exact results only; fails when none is feasible.

- **`distribution`**: cycle-structure, number-of-cycles and largest-cycle
  distributions under the null, with exact references where available.
- **`observed`** (patients): what the observed rejoins determine, the cycle
  structures of AMGs consistent with them, and their probability under the null.
- **`rejoin_probability`** (patients): the null probability of each observed
  rejoin, and of all of them together.

## Output

Every run prints a report and writes to `outputs/<date>/<time>/`:

| File | Content |
|---|---|
| `config.yaml` | the fully resolved configuration (Hydra also keeps `.hydra/`) |
| `results.json` | every reported number, exact values as fractions, and data provenance |
| `cycle_distribution.csv` / `null_cycle_distribution.csv` | cycle structures: exact and sampled probabilities, standard errors |
| `cycle_count_distribution.csv`, `largest_cycle_distribution.csv` | scalar summaries (analysis=distribution) |
| `completion_distribution.csv`, `observed_cycle_structures.csv` | completions of the observed rejoins vs the null (analysis=observed) |
| `rejoin_probabilities.csv` | each observed rejoin: genomic ends, source CSV line, null probability |
| `*.png` | figures (only with the `plots` extra) |

## Interpretation: combinatorial null vs biological probability

All probabilities are computed under a **uniform combinatorial null model**.
Every AMG in the chosen state space and compatible with Θ is taken as equally
likely. They answer: *"if every admissible rejoining were equally likely, how
often would we see this?"*

They are **not biological probabilities**. Real repair favours nearby ends,
depends on nuclear architecture, and is shaped by selection, and none of this
is modelled. A small null probability means the observation is unusual
*relative to uniform rejoining*. It does not measure how likely the
rearrangement is in a cell.

## Project status

**Shared core (implemented)**
- Exact mathematics: the AMG representation, exact enumeration, the paper's
  theorems, and derived formulas (exact PROPER distributions, cycle summaries
  for any n).
- Patient data: Sheth/Baca dataset import, conversion with provenance and the
  corrected strand convention, reproduction of the paper's case study, and
  patient-vs-null comparison. Junction homology and foreign-sequence lengths
  are reported with each observed rejoin.
- Hydra experiments with saved JSON/CSV results and optional figures.

**Part I: IID (implemented)**
- Uniform rejection sampling, validated at the state, distribution and
  acceptance-rate levels, including at n = 80.

**Part II: MCMC (implemented)**
- The 2-switch (reversal-type) move, with its neighbour count and proposal
  symmetry proved, and connectivity of DERANGED and PROPER verified
  exhaustively for n ≤ 6.
- The chain (`sampler=mcmc`), validated against the exact kernel and exact
  distributions, including the paper's P05-1657 table 2.
- The proximity target (`target=proximity`: α and an inter-chromosomal
  weight), fitted to the Baca et al. junctions, validated by detailed balance
  and against exact weighted enumeration (`docs/theory/mcmc.md` §5).

**Part III: quantum annealing (`quantum` branch, emulation)**
- The defect Hamiltonian, B-Chebyshev schedule, overlaps and relaxation times,
  computed exactly for n ≤ 7. The ideal amplitude-estimation estimator is
  emulated and compared with classical Monte Carlo.

**Planned (not implemented)**
- Part II: per-patient (hierarchical) proximity models, junction sequence
  features, several chains per run with R̂, and ABC-SMC for parameter inference.
- Part III: the walk reflections and their errors (Montanaro theorems 3.4 and
  3.5), and circuit-level experiments (e.g. PennyLane).
- The paper's reversal networks (the graph of AMGs linked by reversals).
- MLflow experiment tracking.

## Limitations

- The exact PROPER recursion costs `O(3^k)` in the number of chromosomes `k`.
  Exact PROPER totals are computed for k ≤ 12, and full distributions for
  k ≤ 6 and n ≤ 30. Beyond that, sampling (always available) is used.
- Exact enumeration of the completions of observed rejoins is limited to
  about 2·10⁶ completions (16 unmatched ends). Beyond that they are sampled.
- Strand signs are mapped to DSB ends with the dRanger semantics of the source
  data (Baca et al. 2013): `+` joins the DNA left of the break (the DSB's left
  end), `-` the DNA right of it. Sheth et al. (2026) used the reverse mapping;
  `problem.strand_convention=sheth2026` reproduces their figures. Cycle-structure
  counts are the same under both; which ends are observed differs. See
  [docs/data.md](docs/data.md#strand-convention).
- Patients are analysed one chromosome component at a time, assuming no
  exchange between chromosomes that no observed junction links (paper §7.1).

## Documentation

- [docs/demo.md](docs/demo.md): supervisor demonstration script.
- [docs/data.md](docs/data.md): dataset provenance, licensing and conversion.
- [docs/theory/representation.md](docs/theory/representation.md): data model and state spaces.
- [docs/theory/exact.md](docs/theory/exact.md): exact counts, the paper's
  theorems, and documented differences from the paper.
- [docs/theory/sampling.md](docs/theory/sampling.md) (Part I): uniform sampling,
  rejection, IID proofs, and fixed-rejoin completion.
- [docs/theory/mcmc.md](docs/theory/mcmc.md) (Part II): the 2-switch move, the
  Markov chain, standard errors for correlated samples, and the proximity
  target and its calibration.
- `docs/theory/quantum.md` (Part III, on the `quantum` branch): PROPER counting
  as a Montanaro (2015) partition-function problem, with an exact small-n
  emulation.

Results in the theory notes are labelled PROVED, VERIFIED(n ≤ N) or CONJECTURE.

## Repository layout

```
amg_sampling/
  __init__.py, _version.py   package version, read from the installed metadata
  directories.py             absolute paths to package, repository, data and test directories
  py.typed                   marks the package as typed for consumers
  core/        shared: Θ, rejoin matchings, cycle structures, state spaces, 2-switch moves
  exact/       shared: enumeration, formulas, summaries; weighted.py (II): exact weighted references
  samplers/    uniform.py, rejection.py (I): IID; mcmc.py, targets.py (II): chain, proximity target
  data/sheth/  shared: dataset records, download, loader, conversion to Θ + observed rejoins
  analysis/    shared: null-model analyses, report, JSON/CSV output, figures;
               calibration.py (II): fitting the proximity target
  conf/        Hydra configuration groups
  app.py       Hydra application;  cli.py  command-line entry point
  tests/       pytest suite (exact oracles, statistical tests, CLI smoke tests)
    test_data/ synthetic fixtures (no real patient data)
data/          repository-level data; external/ holds the dataset and is git-ignored
notebooks/     Jupyter notebooks for exploration or documentation
scripts/       import_boundaries.py (tooling); biological_toy.py (II): proximity fit
benchmarks/    iid_n80.py (I): sampling throughput at n = 80
docs/          demo, data and theory notes
.github/       continuous integration and Dependabot
pyproject.toml project metadata, dependencies and tool configuration
poetry.lock    dependency lockfile
Makefile       formatting, linting, testing and coverage commands
```

The package directory sits at the repository root rather than under `src/`,
and the tests live inside it, so they ship with the wheel. `make test-wheel`
runs them against the installed package.

## Development

The Makefile runs every tool inside the Poetry environment:

| Command                  | Description                                                        |
|--------------------------|--------------------------------------------------------------------|
| `make test`              | Run the test suite.                                                |
| `make format`            | Apply safe lint fixes, then format, with Ruff.                     |
| `make format-check`      | Verify formatting without rewriting files. Used by CI.             |
| `make lint`              | Run Ruff lint checks.                                              |
| `make check`             | Formatting check, lint and tests. Does not modify files.           |
| `make coverage`          | Run tests with coverage enforcement.                               |
| `make coverage-html`     | Create an HTML coverage report.                                    |
| `make import-boundaries` | Verify optional dependencies stay isolated behind one module.      |
| `make deps-check`        | Verify imported packages are declared, and in the right group.     |
| `make deadcode`          | Report unused code. Advisory, not a gate.                          |
| `make deadcode-baseline` | Baseline existing dead code so only new dead code surfaces.        |
| `make test-wheel`        | Build a wheel and run the shipped tests against the installed one. |

`make check` deliberately does not reformat. Run `make format` to fix what is
fixable, then `make check` to verify.

Other useful commands:

```bash
poetry run pytest -m "not slow"                  # fast subset (~20 s)
HYPOTHESIS_PROFILE=ci poetry run pytest          # as in CI (derandomised)
poetry run python benchmarks/iid_n80.py          # IID sampling throughput at n = 80
```

matplotlib is imported only by `amg_sampling.analysis.plots`, which is the
boundary that `make import-boundaries` checks for the `plots` extra.

## License

MIT; see [LICENSE](LICENSE). The Sheth / Baca et al. dataset is not covered by
this license and is not distributed here.
