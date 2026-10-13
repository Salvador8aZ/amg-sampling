"""Biologically weighted null for P05-1657: calibration and exact toy example.

1. Split the Baca et al. patients 80/20 by patient (P05-1657 always in the
   test set), fit the proximity target on the training patients, with and
   without an inter-chromosomal weight, and test both on the held-out
   patients (``amg_sampling.analysis.calibration``).
2. For P05-1657, chromosomes 8 and 12 (945 completions), compute the exact
   cycle-structure distribution under each fitted target, the most probable
   cycle structure and the most likely single completion
   (``amg_sampling.exact.weighted``).
3. Sample the same distributions with the 2-switch Markov chain and compare
   with the exact ones (total variation distance).

Usage:  poetry run python scripts/biological_toy.py [--seed 2026] [--out DIR]
Writes ``results.json`` to DIR (default ``outputs/biological-toy``).
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from amg_sampling.analysis.calibration import (
    fit,
    log_likelihood,
    patient_choices,
    split_patients,
)
from amg_sampling.analysis.null import sample_null
from amg_sampling.core.statespace import StateSpace
from amg_sampling.data.sheth import convert, default_data_path, load_dataset
from amg_sampling.exact.weighted import total_variation, weighted_cycle_distribution
from amg_sampling.samplers.targets import ProximityTarget

PATIENT, COMPONENT = "P05-1657", (8, 12)
SPACE = StateSpace.PROPER


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--samples", type=int, default=50_000)
    parser.add_argument("--out", type=Path, default=Path("outputs/biological-toy"))
    args = parser.parse_args()

    dataset = load_dataset(default_data_path())
    results = {
        "data": {"file": dataset.source.file_name, "sha256": dataset.source.sha256}
    }

    # -- 1. calibration -------------------------------------------------------------
    train, test = split_patients(
        dataset.patients, args.test_fraction, args.seed, force_test=[PATIENT]
    )
    choices = {p: patient_choices(dataset.patients[p]) for p in dataset.patients}
    train_choices = [c for p in train for c in choices[p].choices]
    test_choices = [c for p in test for c in choices[p].choices]
    print(f"Patients: {len(train)} training, {len(test)} test ({', '.join(test)})")
    print(
        f"Choices:  {len(train_choices)} training, {len(test_choices)} test; "
        f"{sum(c.skipped_components for c in choices.values())} components skipped "
        "(an end used by two junctions)\n"
    )
    models = {"uniform": None}
    calibration = {}
    for name, fit_inter in (("alpha only", False), ("alpha + inter", True)):
        f = fit(train_choices, fit_inter=fit_inter)
        gain = (
            log_likelihood(test_choices, f.alpha, f.inter_log_weight)
            - log_likelihood(test_choices, 0.0, 0.0)
        ) / len(test_choices)
        per_patient = {}
        for p in test:
            ch = choices[p].choices
            if ch:
                per_patient[p] = (
                    log_likelihood(ch, f.alpha, f.inter_log_weight)
                    - log_likelihood(ch, 0.0, 0.0)
                ) / len(ch)
        improved = sum(v > 0 for v in per_patient.values())
        calibration[name] = {
            "alpha": f.alpha,
            "alpha_se": f.alpha_se,
            "inter_log_weight": f.inter_log_weight,
            "inter_log_weight_se": f.inter_log_weight_se,
            "converged": f.converged,
            "train_gain_per_choice": (f.log_likelihood - f.log_likelihood_uniform)
            / f.num_choices,
            "test_gain_per_choice": gain,
            "test_patients_improved": f"{improved}/{len(per_patient)}",
            "test_gain_by_patient": per_patient,
        }
        inter = (
            f", inter log-weight = {f.inter_log_weight:.3f}"
            f" ± {f.inter_log_weight_se:.3f}"
            if fit_inter
            else ""
        )
        print(f"Model '{name}': alpha = {f.alpha:.3f} ± {f.alpha_se:.3f}{inter}")
        print(
            f"  log pseudo-likelihood gain over uniform per choice: "
            f"train {calibration[name]['train_gain_per_choice']:+.3f}, "
            f"test {gain:+.3f} ({improved}/{len(per_patient)} test patients improved)"
        )
        models[name] = f
    results["split"] = {"seed": args.seed, "train": train, "test": test}
    results["calibration"] = calibration

    # -- 2 and 3. exact toy example and the chain --------------------------------------
    pc = convert(dataset.patients[PATIENT], COMPONENT)
    fixed = pc.fixed
    print(
        f"\n{PATIENT}, chromosomes {COMPONENT}: {pc.describe_theta()}, "
        f"{len(pc.observed)} observed rejoins, {len(pc.free_ends)} free ends"
    )
    toy = {}
    for name, f in models.items():
        target = (
            None if f is None else ProximityTarget.from_patient(pc, **f.target_kwargs())
        )
        exact = weighted_cycle_distribution(pc.theta, SPACE, target, fixed)
        free = set(pc.free_ends)
        map_edges = [
            (u, v)
            for u, v in enumerate(exact.map_state.partners)
            if u < v and u in free
        ]
        entry = {
            "num_completions": exact.num_states,
            "distribution": {str(c): p for c, p in exact.probabilities.items()},
            "most_probable_structure": str(exact.most_probable_structure),
            "map_structure": str(exact.map_structure(pc.theta)),
            "map_probability": exact.map_probability,
            "num_map_states": exact.num_map_states,
            "map_completion_edges": [
                [pc.ends[u].label(), pc.ends[v].label()] for u, v in map_edges
            ],
        }
        print(f"\n  [{name}]")
        print(
            "   exact: "
            + ", ".join(
                f"{c} {p:.3f}" for c, p in list(exact.probabilities.items())[:5]
            )
        )
        print(
            f"   most probable structure: {exact.most_probable_structure} "
            f"({exact.probabilities[exact.most_probable_structure]:.3f})"
        )
        if exact.num_map_states > 1:
            print(
                f"   most likely completion:  {exact.num_map_states} completions tie "
                "(no single most likely one)"
            )
            entry["map_structure"] = entry["map_completion_edges"] = None
        else:
            print(
                f"   most likely completion:  {entry['map_structure']} "
                f"(probability {exact.map_probability:.4f}; "
                f"uniform 1/{exact.num_states})"
            )
            for a, b in entry["map_completion_edges"]:
                print(f"      {a}  —  {b}")
        if target is not None:
            sample = sample_null(
                pc.theta,
                SPACE,
                random.Random(args.seed),
                args.samples,
                fixed=fixed,
                method="mcmc",
                burn_in=5_000,
                thin=5,
                target=target,
            )
            sampled = {
                str(c): m / sample.num_samples
                for c, m in sample.cycle_structures.items()
            }
            tv = total_variation(sampled, entry["distribution"])
            entry["mcmc"] = {
                "samples": sample.num_samples,
                "acceptance": sample.stats.accepted / sample.stats.proposals,
                "distribution": sampled,
                "total_variation_to_exact": tv,
            }
            print(
                f"   MCMC ({sample.num_samples} states, acceptance "
                f"{entry['mcmc']['acceptance']:.2f}): total variation to exact {tv:.4f}"
            )
        toy[name] = entry
    results["toy"] = {
        "patient": PATIENT,
        "chromosomes": list(COMPONENT),
        "space": SPACE.value,
        "models": toy,
    }

    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / "results.json"
    path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nWritten {path}")


if __name__ == "__main__":
    main()
