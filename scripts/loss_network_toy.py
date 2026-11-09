"""Per-junction rejoin probabilities, after the telecom loss network of Wurtz et al.

Wurtz et al. (2022, arXiv:2205.08500, §XI, pp. 22–23) ask of a phone network
"what is the chance a call on route r gets through?" and answer with a ratio of
two partition functions over independent sets. For AMGs the routes are possible
rejoins (u, v), two rejoins clash when they share a free end, and the weights
are Part II's proximity weights ``w(r) = ∏ w(e)``. The same question becomes
"what is the chance end u rejoins end v?":

    P(u–v) = Σ_{r ∋ u–v} w(r) / Σ_r w(r) = Z(u–v forced) / Z

computed two ways that must agree: as a weighted average of the indicator
``[u–v ∈ r]`` over one enumeration (the paper's ⟨F⟩), and as the ratio of two
partition functions (the paper's loss-network formula, and the kind of quantity
Part III's Montanaro estimator targets).

1. Toy: 4 free ends at 0, 4, 1, 5 Mb (a, b, c, d), α = 1, every matching (ALL).
2. P05-1657, chromosomes 8 and 12 (PROPER, 5 observed rejoins, 10 free ends):
   the partner table of every free end under uniform and fitted proximity.
3. Leave-one-out: hide each observed junction in turn, recompute, and report
   how likely the model finds it.

Fitted parameters are read from ``outputs/biological-toy/results.json`` (run
``scripts/biological_toy.py`` first) unless ``--alpha``/``--inter`` are given.
Exact enumeration only: fine here (≤ 10 395 completions), not for large n.

Usage:  poetry run python scripts/loss_network_toy.py [--out DIR]
Writes ``results.json`` to DIR (default ``outputs/loss-network``).
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from amg_sampling import directories
from amg_sampling.analysis.calibration import fit, patient_choices, split_patients
from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import PartialMatching
from amg_sampling.core.statespace import StateSpace
from amg_sampling.data.sheth import convert, default_data_path, load_dataset
from amg_sampling.exact.enumerate import iter_states
from amg_sampling.samplers.targets import MB, ProximityTarget

PATIENT, COMPONENT = "P05-1657", (8, 12)
FITS = directories.base("outputs/biological-toy/results.json")
SEED, TEST_FRACTION = 2026, 0.2  # the split scripts/biological_toy.py uses


def logsumexp(values: list[float]) -> float:
    if not values:
        return -math.inf
    top = max(values)
    if top == -math.inf:
        return top
    return top + math.log(sum(math.exp(v - top) for v in values))


def log_weight(target, partners) -> float:
    return 0.0 if target is None else target.log_weight(partners)


def log_partition(theta, space, target, fixed) -> float:
    """``log Z``: log of the total weight of every state containing ``fixed``."""
    return logsumexp(
        [log_weight(target, s.partners) for s in iter_states(theta, space, fixed)]
    )


def junction_table(theta, space, target, fixed):
    """P(u–v) for every pair of free ends, as a weighted average (one pass)."""
    free = set(fixed.free_ends)
    mass: dict[tuple[int, int], list[float]] = defaultdict(list)
    every: list[float] = []
    for state in iter_states(theta, space, fixed):
        lw = log_weight(target, state.partners)
        every.append(lw)
        for u, v in enumerate(state.partners):
            if u < v and u in free:
                mass[(u, v)].append(lw)
    log_z = logsumexp(every)
    return (
        {e: math.exp(logsumexp(m) - log_z) for e, m in mass.items()},
        log_z,
        len(every),
    )


def ratio_check(theta, space, target, fixed, table, log_z) -> float:
    """Recompute every P(u–v) as Z(u–v forced) / Z; return the largest gap."""
    free = fixed.free_ends
    worst = 0.0
    for i, u in enumerate(free):
        for v in free[i + 1 :]:
            forced = PartialMatching(fixed.num_ends, (*fixed.pairs(), (u, v)))
            p = math.exp(log_partition(theta, space, target, forced) - log_z)
            worst = max(worst, abs(p - table.get((u, v), 0.0)))
    return worst


def partners_of(end: int, table) -> list[tuple[int, float]]:
    """Every partner of ``end`` with its probability, most likely first."""
    out = [(v if u == end else u, p) for (u, v), p in table.items() if end in (u, v)]
    return sorted(out, key=lambda item: -item[1])


# -- 1. toy ---------------------------------------------------------------------------


def toy() -> dict:
    names = "abcd"
    theta = InitialConfiguration((2,))  # 2 DSBs on one chromosome: 4 free ends
    coords = (0, 4 * MB, 1 * MB, 5 * MB)  # a, b, c, d
    target = ProximityTarget((0, 0, 0, 0), coords, alpha=1.0)
    fixed = PartialMatching(theta.num_ends)
    print("1. TOY: ends a, b, c, d at 0, 4, 1, 5 Mb; w = 1/distance (α = 1); ALL")
    out = {}
    for name, t in (("uniform", None), ("alpha=1", target)):
        table, log_z, num = junction_table(theta, StateSpace.ALL, t, fixed)
        gap = ratio_check(theta, StateSpace.ALL, t, fixed, table, log_z)
        print(
            f"   [{name}]  {num} matchings, Z = {math.exp(log_z):.4f}, "
            f"ratio-form gap {gap:.1e}"
        )
        for v, p in partners_of(0, table):
            print(f"      P(a–{names[v]}) = {p:.3f}")
        out[name] = {f"{names[u]}-{names[v]}": p for (u, v), p in table.items()}
    return out


# -- 2 and 3. patient -----------------------------------------------------------------


def models_from(args) -> dict:
    """Uniform plus the proximity fits: ``--alpha``, else the saved fits of
    ``scripts/biological_toy.py``, else the same fits recomputed here.
    """
    if args.alpha is not None:
        return {"uniform": None, "fitted": (args.alpha, args.inter or 0.0)}
    if FITS.exists():
        cal = json.loads(FITS.read_text())["calibration"]
        print(f"Fitted parameters from {FITS}")
        return {
            "uniform": None,
            "alpha only": (cal["alpha only"]["alpha"], 0.0),
            "alpha + inter": (
                cal["alpha + inter"]["alpha"],
                cal["alpha + inter"]["inter_log_weight"],
            ),
        }
    print(f"{FITS} not found: fitting the proximity target (seed {SEED}) ...")
    patients = load_dataset(default_data_path()).patients
    train, _ = split_patients(patients, TEST_FRACTION, SEED, force_test=[PATIENT])
    choices = [c for p in train for c in patient_choices(patients[p]).choices]
    alpha_only = fit(choices, fit_inter=False)
    both = fit(choices, fit_inter=True)
    return {
        "uniform": None,
        "alpha only": (alpha_only.alpha, 0.0),
        "alpha + inter": (both.alpha, both.inter_log_weight),
    }


def patient(models: dict) -> dict:
    pc = convert(load_dataset(default_data_path()).patients[PATIENT], COMPONENT)
    theta, fixed, space = pc.theta, pc.fixed, StateSpace.PROPER
    label = lambda e: pc.ends[e].label()  # noqa: E731
    targets = {
        name: None
        if m is None
        else ProximityTarget.from_patient(pc, alpha=m[0], inter_log_weight=m[1])
        for name, m in models.items()
    }
    print(
        f"\n2. {PATIENT}, chromosomes {COMPONENT}: {pc.describe_theta()}, "
        f"{len(pc.observed)} observed rejoins, {len(pc.free_ends)} free ends, PROPER"
    )
    tables, out = {}, {"partners": {}, "leave_one_out": []}
    for name, t in targets.items():
        table, log_z, num = junction_table(theta, space, t, fixed)
        gap = ratio_check(theta, space, t, fixed, table, log_z)
        tables[name] = table
        print(f"   [{name}] {num} completions; ratio-form gap {gap:.1e}")
    print("\n   Most likely partner of each free end (probability per model):")
    head = "".join(f"{n:>15}" for n in targets)
    print(f"   {'free end':<28}{'top partner (fitted)':<30}{head}")
    last = list(targets)[-1]
    for u in pc.free_ends:
        v, _ = partners_of(u, tables[last])[0]
        e = (min(u, v), max(u, v))
        probs = "".join(f"{tables[n].get(e, 0.0):>15.3f}" for n in targets)
        print(f"   {label(u):<28}{label(v):<30}{probs}")
        out["partners"][label(u)] = {
            n: {label(w): p for w, p in partners_of(u, tables[n])} for n in targets
        }

    print("\n3. Leave-one-out: hide one observed junction, how likely is it?")
    print(f"   {'observed junction':<56}{head}")
    for edge in pc.observed:
        rest = PartialMatching(
            theta.num_ends,
            [(e.end_a, e.end_b) for e in pc.observed if e is not edge],
        )
        e = (min(edge.end_a, edge.end_b), max(edge.end_a, edge.end_b))
        row = {}
        for name, t in targets.items():
            table, _, _ = junction_table(theta, space, t, rest)
            p = table.get(e, 0.0)
            others = [q for _, q in partners_of(edge.end_a, table)]
            rank = 1 + sum(q > p + 1e-12 for q in others)
            ties = sum(abs(q - p) <= 1e-12 for q in others)  # includes the true one
            row[name] = {
                "probability": p,
                "rank_among_partners": rank,
                "tied_with": ties - 1,
                "num_candidates": len(others),
            }
        text = f"{label(edge.end_a)} — {label(edge.end_b)}"
        cells = ""
        for n in targets:
            r = row[n]
            k, n_ = r["rank_among_partners"], r["num_candidates"]
            where = (
                f"#{k}-{k + r['tied_with']}/{n_}" if r["tied_with"] else f"#{k}/{n_}"
            )
            cells += f"{r['probability']:>8.3f} {where:<8}"
        print(f"   {text:<56}{cells}")
        out["leave_one_out"].append(
            {
                "junction": text,
                "sv_class": edge.variant.sv_class if edge.variant else None,
                **row,
            }
        )
    print("   (#k/N: rank of the true partner among the N candidate partners of the")
    print(
        "    first end; #k-m/N: tied for ranks k to m. Uniform is nearly flat, tilted"
    )
    print("    only by PROPER's constraints.)")
    out["models"] = {
        n: None if m is None else {"alpha": m[0], "inter_log_weight": m[1]}
        for n, m in models.items()
    }
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--alpha", type=float, help="proximity exponent (default: fitted)"
    )
    parser.add_argument("--inter", type=float, help="inter-chromosomal log-weight")
    parser.add_argument(
        "--out", type=Path, default=directories.base("outputs/loss-network")
    )
    args = parser.parse_args()
    results = {"toy": toy(), "patient": patient(models_from(args))}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "results.json").write_text(json.dumps(results, indent=2))
    print(f"\nWrote {args.out / 'results.json'}")


if __name__ == "__main__":
    main()
