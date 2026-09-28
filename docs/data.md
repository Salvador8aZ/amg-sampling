# The Sheth patient dataset

## What it is

Sheth, Arsuaga & Sazdanovic (2026) analyse prostate-cancer rearrangements
from Baca et al., *Punctuated evolution of prostate cancer genomes*,
Cell 153:666–677 (2013). Their code uses a CSV of that study's
structural-variant calls, distributed as `data/nihms.csv` in
[`siddharthsheth/aberration_multigraph`](https://github.com/siddharthsheth/aberration_multigraph)
(also in `Salvador8aZ/aberration_multigraph`).

- 5 710 rows (one rearrangement junction each), 57 patients.
- Reference copy SHA-256:
  `13a70c63d09af08b423fe372fd09a4d6c1299bb9ad5a1424e9cc490af1a6a3bc`.

## Why it is not in this repository

The reference repository's MIT license covers its code. No licence statement
covers the data file itself, and the underlying calls come from a published
study's supplementary material. Redistribution rights are therefore unclear,
so **this repository does not include the data**. `data/external/` is in
`.gitignore`.

## How to obtain it

1. Clone the reference repository (or download only `data/nihms.csv` from it):

   ```bash
   git clone https://github.com/siddharthsheth/aberration_multigraph
   ```

2. Import the file. This checks the format and the checksum, then copies it
   to `data/external/nihms.csv`:

   ```bash
   uv run amg-sampling import-data aberration_multigraph/data/nihms.csv
   ```

   Alternatively set `AMG_SHETH_DATA=/path/to/nihms.csv`, or pass
   `data.sheth_csv=/path/to/nihms.csv` to an experiment.

Every analysis records the file name, its SHA-256 and whether it is the
reference copy in `results.json`.

## Fields used

See `amg_sampling/data/sheth/models.py`. From each row:

| Column | Use |
|---|---|
| `Individual` | patient identifier |
| `Number` | junction number within the patient (provenance only) |
| `Breakpoint 1/2 chromosome` | chromosome (integers 1–24; 23/24 are not relabelled) |
| `Breakpoint 1/2 position` | breakpoint base-pair position |
| `Breakpoint 1/2 strand` | `+` or `-` |
| `Class` | reported rearrangement class (provenance only) |

The CSV line number of every junction is kept and reported with each edge.

## From junctions to Θ and observed rejoins

See `amg_sampling/data/sheth/conversion.py` for details.

| Step | Status |
|---|---|
| Junction breakpoints (chromosome, position, strand) | **observed** (in the data) |
| Patient split into chromosome components linked by junctions | Sheth methodology (paper §7.1): no exchange assumed between unlinked chromosomes |
| Every distinct breakpoint is one DSB, ordered along its chromosome | reconstructed |
| Strand `+` → the DSB's right end, `-` → its left end | reference-repository convention; reproduces the paper's table 2, but the biological strand semantics are not independently verified |
| Each junction is a rejoin edge between those two ends | observed junction, with a reconstructed end assignment |
| Ends not in any junction | **unmatched**: their partners are unknown, and completions are enumerated or sampled |
| The only possible completion when exactly two ends are unmatched | reconstructed (labelled "unique completion", never "observed") |

Conversion refuses, rather than silently alters, components where a DSB end
would be used by two junctions. This affects 6 components in 6 patients
(listed by `amg-sampling patients`).
