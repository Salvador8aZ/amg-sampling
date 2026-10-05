# The Sheth patient dataset

## What it is

Sheth, Arsuaga & Sazdanovic (2026) analyse prostate-cancer rearrangements
from Baca et al., *Punctuated evolution of prostate cancer genomes*,
Cell 153:666–677 (2013). Their code uses a CSV of that study's
structural-variant calls, distributed as `data/nihms.csv` in
[`siddharthsheth/aberration_multigraph`](https://github.com/siddharthsheth/aberration_multigraph)
(also in `Salvador8aZ/aberration_multigraph`).

- 5 710 rows (one rearrangement junction each), 57 patients.
- Original publication: Baca et al., Cell 153:666–677 (2013),
  doi:[10.1016/j.cell.2013.03.027](https://doi.org/10.1016/j.cell.2013.03.027).
- Reference copy SHA-256:
  `13a70c63d09af08b423fe372fd09a4d6c1299bb9ad5a1424e9cc490af1a6a3bc`.

## Why it is not in this repository

The reference repository's MIT license covers its code. No licence statement
covers the data file itself, and the underlying calls come from a published
study's supplementary material. Redistribution rights are therefore unclear,
so **this repository does not include the data**. `data/external/` is in
`.gitignore`.

## How to obtain it

The file is downloaded from its source, never from this repository, in the
style of [Pooch](https://www.fatiando.org/pooch/) (as `python-minecraft-data`
does for Minecraft data):

- the URL is pinned to commit `5cc6a8b` of `siddharthsheth/aberration_multigraph`,
  so it cannot change;
- the download is kept only if its SHA-256 matches the reference copy;
- it is cached in `data/external/nihms.csv`, with a provenance and citation
  note in `data/external/nihms.csv.SOURCE.txt`. Both are git-ignored.

Either download it explicitly:

```bash
make data   # = poetry run amg-sampling fetch-data
```

or just run a patient experiment: if `data.sheth_csv` does not exist, it is
fetched first. The source is set in `amg_sampling/conf/config.yaml`
(`data.url`, `data.sha256`); pass `data.fetch=false` to work offline.

To use a local copy instead, run
`poetry run amg-sampling import-data /path/to/nihms.csv`, set
`AMG_SHETH_DATA=/path/to/nihms.csv`, or pass `data.sheth_csv=/path/to/nihms.csv`.
An existing file is never overwritten by a download.

If the download fails with `CERTIFICATE_VERIFY_FAILED` on macOS, the
python.org installer's certificates are missing: run
`/Applications/Python 3.x/Install Certificates.command`, or prefix the command
with `SSL_CERT_FILE=/etc/ssl/cert.pem`.

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
| `Homology length` | base pairs of microhomology at the junction (optional; reported, not used by the model) |
| `Foreign sequence length` | base pairs of non-templated sequence inserted at the junction (optional; reported, not used by the model) |

The two sequence-feature columns may be absent (the file still loads). In the
reference copy, 124 junctions have `-1` in both columns, with `Foreign
sequence` = `failed`: the junction sequence could not be assembled. One
junction (line 5573) has `NaN` in both. All of these are stored as "not
measured", never as 0. The features hint at the repair mechanism (for
example, microhomology of a few base pairs is typical of
microhomology-mediated end joining). They are carried with every observed
junction, in `results.json`, `rejoin_probabilities.csv` and the terminal
report, for future biologically weighted models.

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
