# amg-sampling

Exact enumeration and sampling of aberration multigraph (AMG) rejoin
configurations.

The mathematical framework follows

> S. Sheth, J. Arsuaga, R. Sazdanovic (2026). *Characterizing cancer chromosome
> aberration pathways using multigraphs.* J. Phys. A: Math. Theor. **59** 115601.
> https://doi.org/10.1088/1751-8121/ae4d8f

The repository `Salvador8aZ/aberration_multigraph` (a copy of
`siddharthsheth/aberration_multigraph`) is used only as a reference
implementation: its validated numerical results serve as ground truth, but its
API and data model are not reused.

## Scientific goal

Given `k` chromosomes carrying a total of `n` double-strand breaks (DSBs), the
`2n` free DSB ends are paired by a *rejoin matching*. The project aims to

1. represent AMG configurations efficiently,
2. enumerate them exactly for small `n` (ground truth),
3. sample them for large `n` (patients may have 80+ DSBs),
4. estimate the induced distribution over cycle structures.

Theory notes: [`representation.md`](docs/theory/representation.md) (data model,
state spaces) and [`exact.md`](docs/theory/exact.md) (exact counts and
cycle-structure distributions).

## Development

Requires [uv](https://docs.astral.sh/uv/).

```bash
uv sync
uv run pytest
```

Hypothesis profiles: `dev` (default) and `ci` (derandomized). Select with
`HYPOTHESIS_PROFILE=ci uv run pytest`.

## Roadmap

| Stage | Content | Status |
|---|---|---|
| A | Project skeleton; initial configuration Θ; rejoin matchings; cycle structures; state spaces | done |
| B | Exact enumeration (oracle); exact formulas (PROPER recursion is O(3^k) in the number of chromosomes; see `exact.md` §5) | done |
| C | IID uniform and rejection sampling (first sampling baseline) | in progress |
| D | Reversal and switch moves | not started |
| E | Small-n transition-graph laboratory | not started |
| F | MCMC (only after the state space and proposal kernel are established) | not started |
| G | Hydra configuration, MLflow tracking | not started |
| later | Patient-specific constraints; biological weights; ABC-SMC; quantum experiments | not started |

Validation hierarchy: exact enumeration → exact formulas → IID rejection
sampler → MCMC. For small `n` all applicable methods must agree on the
cycle-structure distribution.

## License

MIT — see [`LICENSE`](LICENSE).
