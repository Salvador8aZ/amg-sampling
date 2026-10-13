# The Sheth patient dataset

## What it is

Sheth, Arsuaga & Sazdanovic (2026) analyse prostate-cancer rearrangements
from Baca et al., *Punctuated evolution of prostate cancer genomes*,
Cell 153:666–677 (2013),
doi:[10.1016/j.cell.2013.03.021](https://doi.org/10.1016/j.cell.2013.03.021).
The data are that study's structural-variant calls, **Table S3C** of its
supplementary workbook `mmc3.xlsx`. **The data belong to Baca et al.**

- 5 710 rows (one rearrangement junction each), 57 patients.
- Sheth et al. distribute the same table as `data/nihms.csv` in
  [`siddharthsheth/aberration_multigraph`](https://github.com/siddharthsheth/aberration_multigraph)
  (an Excel "Save as CSV" of Table S3C; "nihms" is the NIH manuscript
  submission of the paper, NIHMS472000).

Two copies are recognised, by SHA-256:

| Copy | SHA-256 | Row order |
|---|---|---|
| Converted from the Baca et al. supplement (default) | `8b128017090d6ccd70fd62be13646467a32b92b3871c55804ae870dd96a75f35` | patient, then junction `Number` |
| Sheth et al. `data/nihms.csv` (fallback) | `13a70c63d09af08b423fe372fd09a4d6c1299bb9ad5a1424e9cc490af1a6a3bc` | the same, except patient P05-1657, listed by breakpoint 1 (chromosome, then position) |

The two have identical rows and values; only the order of P05-1657's nine
rows (lines 312–320) differs. Results do not depend on row order: DSBs are
ordered by position along each chromosome. Only the order in which observed
junctions are listed follows the file.

Neither CSV keeps the workbook's own row order, so a CSV line number does not
point at a row of Table S3C. A junction is identified everywhere (workbook and
both copies) by patient and junction `Number`, which is unique within each
patient (5 710 distinct pairs).

## Why it is not in this repository

No licence statement covers the data. They come from a published study's
supplementary material (and the reference repository's MIT licence covers
its code only). Redistribution rights are therefore unclear, so **this
repository does not include the data**. `data/external/` is in `.gitignore`.

## How to obtain it

The file is downloaded from its source, never from this repository, in the
style of [Pooch](https://www.fatiando.org/pooch/) (as `python-minecraft-data`
does for Minecraft data). Sources are tried in order; the first that verifies
is kept:

1. **The Baca et al. supplement.** `mmc3.xlsx` is downloaded from Elsevier
   (`https://ars.els-cdn.com/content/image/1-s2.0-S0092867413003437-mmc3.xlsx`)
   and Table S3C is converted to `nihms.csv` by
   `amg_sampling/data/sheth/baca.py`. The conversion writes cells as Excel
   displays them (Quality and Score in `0.00E+00` form, e.g. `9.81E-01`),
   drops the sheet's trailing blank rows, sorts rows by patient and junction
   number, and uses a UTF-8 byte-order mark and CRLF line endings, like Excel's
   CSV export. The result must have the SHA-256 above.
2. **The same workbook, found from the DOI.** If the URL no longer works, the
   DOI is resolved to the publisher's current landing page and the workbook
   is looked for there (from the article's PII, and from `.xlsx` links on the
   page). Any workbook found must again convert to the SHA-256 above.
3. **The Sheth et al. CSV**, at pinned commit `5cc6a8b`, which must have its
   own SHA-256 above.

Elsevier's server sits behind Cloudflare, which answers Python's HTTP client
with 403 Forbidden. The download is then retried with `curl` (installed on
macOS, most Linux distributions and Windows 10+).

The file is saved as `data/external/nihms.csv`, with a note in
`data/external/nihms.csv.SOURCE.txt` giving the URL it came from, its SHA-256,
the download date and the citation. After a download the same note is printed.
Both files are git-ignored. If no source verifies, nothing is written and the
error lists each attempt.

Either download it explicitly:

```bash
make data   # = poetry run amg-sampling fetch-data
```

or just run a patient experiment: if `data.sheth_csv` does not exist, it is
fetched first. The sources are set in `amg_sampling/conf/config.yaml`
(`data.baca_url`, `data.doi`, `data.sha256`, `data.sheth_url`,
`data.sheth_sha256`); set a source to `null` to skip it, or pass
`data.fetch=false` to work offline.

To use a local copy instead, run
`poetry run amg-sampling import-data /path/to/nihms.csv`, set
`AMG_SHETH_DATA=/path/to/nihms.csv`, or pass `data.sheth_csv=/path/to/nihms.csv`.
An existing file is never overwritten by a download; delete
`data/external/nihms.csv` to download it again.

If the python.org installer's certificates are missing on macOS, the system
store `/etc/ssl/cert.pem` is used. If the download still fails with
`CERTIFICATE_VERIFY_FAILED`, run
`/Applications/Python 3.x/Install Certificates.command`, or prefix the command
with `SSL_CERT_FILE=/etc/ssl/cert.pem`.

Every analysis records the file name, its SHA-256 and which known copy it is
(if any) in `results.json`.

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

The two sequence-feature columns may be absent (the file still loads). In both
known copies, 124 junctions have `-1` in both columns, with `Foreign
sequence` = `failed`: the junction sequence could not be assembled. One
junction (line 5573) has `NaN` in both. All of these are stored as "not
measured", never as 0. The features hint at the repair mechanism (for
example, microhomology of a few base pairs is typical of
microhomology-mediated end joining). They are carried with every observed
junction, in `results.json`, `rejoin_probabilities.csv` and the terminal
report, for future biologically weighted models.

Every observed edge is reported with its junction `Number` (with the
patient, the identifier shared by Table S3C and both copies) and the line of
the local CSV it was read from, which depends on the copy.

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
