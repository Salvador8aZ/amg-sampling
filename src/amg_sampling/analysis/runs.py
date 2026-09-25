"""The three analyses offered by the application, as plain functions.

Each returns a JSON-serialisable ``dict``. Printing, CSV export and plotting
all read this dictionary, so every reported number is also saved.

Terminology: the *combinatorial null model* is the uniform distribution on a
state space of AMGs compatible with Θ. Probabilities under it answer "if every
admissible AMG were equally likely, how often would we see this?". They are
not biological rejoining probabilities.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from fractions import Fraction

from amg_sampling.analysis.estimates import ProportionEstimate
from amg_sampling.analysis.null import largest_cycle, num_cycles, sample_null, summary_counts
from amg_sampling.analysis.reference import completion_reference, exact_reference
from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import CycleStructure, cycle_structure
from amg_sampling.core.statespace import StateSpace
from amg_sampling.data.sheth.conversion import PatientConfiguration, RejoinEdge
from amg_sampling.data.sheth.models import SourceInfo
from amg_sampling.samplers.rejection import theoretical_acceptance

SCHEMA_VERSION = 1

INTERPRETATION = (
    "Probabilities are computed under a uniform combinatorial null model: every AMG "
    "in the chosen state space and compatible with the initial configuration is "
    "equally likely. They are not biological rejoining probabilities."
)


class AnalysisError(ValueError):
    """The requested analysis does not apply to the requested problem."""


@dataclass(frozen=True)
class Problem:
    kind: str  # "configuration" or "patient"
    theta: InitialConfiguration
    label: str
    patient: PatientConfiguration | None = None
    source: SourceInfo | None = None


@dataclass(frozen=True)
class SamplerSettings:
    name: str  # "iid" or "exact"
    num_samples: int = 100_000
    max_proposals_per_sample: int | None = 1_000_000

    @property
    def samples(self) -> bool:
        return self.name == "iid"


# -- shared pieces ----------------------------------------------------------------


MAX_EXACT_DIGITS = 40  # longer exact fractions and counts are written as strings / omitted


def fraction_dict(value: Fraction | None) -> dict | None:
    """``{"fraction": "p/q", "value": float}``; the exact fraction is omitted when very long."""
    if value is None:
        return None
    text = f"{value.numerator}/{value.denominator}"
    return {"fraction": text if len(text) <= 2 * MAX_EXACT_DIGITS else None, "value": float(value)}


def json_int(value: int | None) -> int | str | None:
    """Integers beyond 2^53 are written as decimal strings so every JSON reader keeps them exact."""
    if value is None or abs(value) < 2**53:
        return value
    return str(value)


def theta_label(theta: InitialConfiguration) -> str:
    return f"Θ({theta.num_chromosomes},({','.join(map(str, theta.breaks))}))"


def problem_dict(problem: Problem) -> dict:
    theta = problem.theta
    out = {
        "kind": problem.kind,
        "label": problem.label,
        "theta": {
            "notation": theta_label(theta),
            "num_chromosomes": theta.num_chromosomes,
            "breaks": list(theta.breaks),
            "num_dsbs": theta.num_dsbs,
            "num_free_ends": theta.num_ends,
        },
    }
    if problem.patient is not None:
        pc = problem.patient
        out["patient_id"] = pc.patient_id
        out["chromosomes"] = list(pc.chromosomes)
        out["observed_rejoins"] = [edge_dict(pc, e, i) for i, e in enumerate(pc.observed, 1)]
        out["num_unmatched_ends"] = len(pc.free_ends)
        out["notes"] = list(pc.notes)
    if problem.source is not None:
        s = problem.source
        out["data_source"] = {"file_name": s.file_name, "sha256": s.sha256, "is_reference_copy": s.is_reference_copy}
    return out


def edge_dict(pc: PatientConfiguration, edge: RejoinEdge, index: int) -> dict:
    v = edge.variant
    return {
        "edge": f"e{index}",
        "end_a": edge.end_a,
        "end_b": edge.end_b,
        "site_a": pc.ends[edge.end_a].label(),
        "site_b": pc.ends[edge.end_b].label(),
        "status": edge.status.value,
        "source_line": v.source_line if v else None,
        "sv_class": v.sv_class if v else None,
    }


def header(analysis: str, problem: Problem, space: StateSpace, sampler: SamplerSettings, seed: int) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "analysis": analysis,
        "problem": problem_dict(problem),
        "state_space": space.value,
        "null_model": f"Uniform({space.value.upper()}(Θ))",
        "sampler": {
            "name": sampler.name,
            "num_samples": sampler.num_samples if sampler.samples else None,
            "seed": seed,
        },
        "interpretation": INTERPRETATION,
    }


def sampling_dict(sample, theoretical: Fraction | None) -> dict:
    stats = sample.stats
    return {
        "samples": sample.num_samples,
        "proposals": stats.proposals,
        "acceptance_rate": stats.acceptance_rate,
        "theoretical_acceptance": fraction_dict(theoretical),
    }


def safe_theoretical_acceptance(theta: InitialConfiguration, space: StateSpace) -> Fraction | None:
    try:
        return theoretical_acceptance(theta, space)
    except ValueError:
        return None


def distribution_rows(exact: dict | None, sampled: dict | None, key_name: str) -> list[dict]:
    """Merge exact and sampled counts into rows (by probability for cycle structures, else by key)."""
    keys = set(exact or {}) | set(sampled or {})
    exact_total = sum((exact or {}).values())
    sampled_total = sum((sampled or {}).values())
    rows = []
    for key in keys:
        row = {key_name: str(key) if isinstance(key, CycleStructure) else key}
        if isinstance(key, CycleStructure):
            row["num_cycles"] = key.num_cycles
        if exact is not None:
            count = exact.get(key, 0)
            row["exact_count"] = json_int(count)
            row["exact_probability"] = fraction_dict(Fraction(count, exact_total)) if exact_total else None
        if sampled is not None:
            est = ProportionEstimate(sampled.get(key, 0), sampled_total)
            row["sampled_count"] = est.hits
            row["sampled_probability"] = est.estimate
            row["standard_error"] = est.standard_error
        rows.append(row)

    def order(row):
        p = row.get("exact_probability")
        return (-(p["value"] if p else 0.0), -row.get("sampled_probability", 0.0), str(row[key_name]))

    if key_name == "cycle_structure":
        rows.sort(key=order)
    else:
        rows.sort(key=lambda row: row[key_name])
    return rows


def exact_dict(ref) -> dict:
    return {
        "total_states": json_int(ref.total_states),
        "cycle_distribution_available": ref.cycle_distribution is not None,
        "summaries_space": ref.summaries_space,
        "total_variation_bound": fraction_dict(ref.tv_bound),
        "unavailable": ref.unavailable,
    }


# -- analysis=distribution ------------------------------------------------------------


def run_distribution(problem: Problem, space: StateSpace, sampler: SamplerSettings, seed: int) -> dict:
    # For a patient this describes the unconditional null Uniform(space(Θ));
    # observed rejoins are used by analysis=observed.
    theta = problem.theta
    ref = exact_reference(theta, space)
    results = header("distribution", problem, space, sampler, seed)
    results["exact"] = exact_dict(ref)
    sampled = None
    if sampler.samples:
        sample = sample_null(
            theta, space, random.Random(seed), sampler.num_samples,
            max_proposals_per_sample=sampler.max_proposals_per_sample,
        )
        results["sampling"] = sampling_dict(sample, safe_theoretical_acceptance(theta, space))
        sampled = sample.cycle_structures
    elif ref.cycle_distribution is None and ref.cycle_counts is None:
        raise AnalysisError(f"No exact result is available for {theta_label(theta)}: {ref.unavailable}. Use sampler=iid.")

    results["cycle_distribution"] = distribution_rows(ref.cycle_distribution, sampled, "cycle_structure")
    results["cycle_count_distribution"] = distribution_rows(
        ref.cycle_counts, summary_counts(sampled, num_cycles) if sampled else None, "num_cycles"
    )
    results["largest_cycle_distribution"] = distribution_rows(
        ref.largest_cycle, summary_counts(sampled, largest_cycle) if sampled else None, "largest_cycle"
    )
    results["means"] = {
        "num_cycles": mean_summary(ref.cycle_counts, summary_counts(sampled, num_cycles) if sampled else None),
        "largest_cycle": mean_summary(ref.largest_cycle, summary_counts(sampled, largest_cycle) if sampled else None),
    }
    return results


def mean_summary(exact: dict | None, sampled: dict | None) -> dict:
    """Exact mean, and sample mean with standard error sd/sqrt(N)."""
    out = {"exact_mean": None, "sampled_mean": None, "sampled_standard_error": None}
    if exact:
        total = sum(exact.values())
        out["exact_mean"] = float(Fraction(sum(k * m for k, m in exact.items()), total))
    if sampled:
        n = sum(sampled.values())
        mean = sum(k * m for k, m in sampled.items()) / n
        var = sum(m * (k - mean) ** 2 for k, m in sampled.items()) / max(n - 1, 1)
        out["sampled_mean"] = mean
        out["sampled_standard_error"] = (var / n) ** 0.5
    return out


# -- analysis=observed and analysis=rejoin_probability -------------------------------------


def _require_patient(problem: Problem, analysis: str) -> PatientConfiguration:
    if problem.patient is None:
        raise AnalysisError(f"analysis={analysis} needs observed rejoins; use problem=patient.")
    return problem.patient


def _null_and_edges(problem, space, sampler, seed, confidence, results) -> tuple:
    pc = problem.patient
    theta = problem.theta
    ref = exact_reference(theta, space)
    results["null"] = {"exact": exact_dict(ref)}
    edges = [(e.end_a, e.end_b) for e in pc.observed]
    sample = None
    if sampler.samples:
        sample = sample_null(
            theta, space, random.Random(seed), sampler.num_samples, track_edges=edges,
            max_proposals_per_sample=sampler.max_proposals_per_sample,
        )
        results["null"]["sampling"] = sampling_dict(sample, safe_theoretical_acceptance(theta, space))

    rejoins = []
    for i, edge in enumerate(pc.observed):
        row = edge_dict(pc, edge, i + 1)
        if sample is not None:
            est = ProportionEstimate(sample.edge_hits[i], sample.num_samples, confidence)
            row["null_probability"] = est.to_dict()
            row["null_probability_text"] = est.describe()
        rejoins.append(row)
    results["rejoins"] = rejoins
    return ref, sample


def _joint(problem, space, ref, completions, sample, confidence) -> dict:
    joint = {"event": "every observed rejoin is present", "num_observed_rejoins": len(problem.patient.observed)}
    if completions is not None and completions.cycle_distribution is not None and ref.total_states:
        numerator = sum(completions.cycle_distribution.values())
        joint["exact"] = fraction_dict(Fraction(numerator, ref.total_states))
        joint["exact_numerator"] = json_int(numerator)
        joint["exact_denominator"] = json_int(ref.total_states)
    else:
        joint["exact"] = None
    if sample is not None and problem.patient.observed:
        est = ProportionEstimate(sample.joint_hits, sample.num_samples, confidence)
        joint["sampled"] = est.to_dict()
        joint["sampled_text"] = est.describe()
    return joint


def run_rejoin_probability(
    problem: Problem, space: StateSpace, sampler: SamplerSettings, seed: int, confidence: float = 0.95
) -> dict:
    pc = _require_patient(problem, "rejoin_probability")
    if not sampler.samples:
        raise AnalysisError("analysis=rejoin_probability estimates probabilities by sampling; use sampler=iid.")
    results = header("rejoin_probability", problem, space, sampler, seed)
    ref, sample = _null_and_edges(problem, space, sampler, seed, confidence, results)
    completions = completion_reference(problem.theta, space, pc.fixed)
    results["joint_observed_rejoins"] = _joint(problem, space, ref, completions, sample, confidence)
    return results


def run_observed(
    problem: Problem, space: StateSpace, sampler: SamplerSettings, seed: int, confidence: float = 0.95
) -> dict:
    pc = _require_patient(problem, "observed")
    theta = problem.theta
    results = header("observed", problem, space, sampler, seed)
    ref, sample = _null_and_edges(problem, space, sampler, seed, confidence, results)

    # 1. What the data determine.
    completions = completion_reference(theta, space, pc.fixed)
    observed = {
        "num_observed_rejoins": len(pc.observed),
        "num_unmatched_ends": len(pc.free_ends),
        "num_completions_all": completions.num_completions_all,
    }
    reconstructed = pc.reconstructed_matching()
    c_obs = None
    if reconstructed is not None:
        matching, inferred = reconstructed
        c_obs = cycle_structure(theta, matching)
        edges = list(pc.observed) + list(inferred)
        observed["reconstructed"] = {
            "rejoins": [edge_dict(pc, e, i) for i, e in enumerate(edges, 1)],
            "cycle_structure": str(c_obs),
            "membership": {s.value: s.contains(theta, matching) for s in StateSpace},
        }
    else:
        observed["reconstructed"] = None

    # 2. Completions of the observed rejoins within the state space.
    comp = {"space": space.value, "exact_count": None, "unavailable": completions.unavailable}
    comp_exact = completions.cycle_distribution
    if comp_exact is not None:
        comp["exact_count"] = sum(comp_exact.values())
    comp_sampled = None
    if sampler.samples and reconstructed is None and (comp_exact is None or sum(comp_exact.values()) > 0):
        try:
            cond = sample_null(
                theta, space, random.Random(seed + 1_000_003), sampler.num_samples, fixed=pc.fixed,
                max_proposals_per_sample=sampler.max_proposals_per_sample,
            )
        except RuntimeError as error:
            comp["sampling_error"] = str(error)
        else:
            comp_sampled = cond.cycle_structures
            comp["sampling"] = sampling_dict(cond, None)
    comp["distribution"] = distribution_rows(comp_exact, comp_sampled, "cycle_structure")
    observed["completions"] = comp
    results["observed"] = observed

    # 3. The observed (or each completion-compatible) cycle structure under the null.
    if c_obs is not None:
        structures = [c_obs]
    else:
        structures = list(comp_exact or comp_sampled or {})
    null_total = sum(sample.cycle_structures.values()) if sample else 0
    rows = []
    for c in structures:
        row = {"cycle_structure": str(c), "num_cycles": c.num_cycles}
        if comp_exact:
            row["share_among_completions"] = fraction_dict(Fraction(comp_exact.get(c, 0), sum(comp_exact.values())))
        if comp_sampled:
            share = ProportionEstimate(comp_sampled.get(c, 0), sum(comp_sampled.values()), confidence)
            row["share_among_completions_sampled"] = share.to_dict()
            row["share_among_completions_sampled_text"] = share.describe()
        row["null_exact"] = fraction_dict(ref.probability(ref.cycle_distribution, c))
        if sample is not None:
            est = ProportionEstimate(sample.cycle_structures.get(c, 0), null_total, confidence)
            row["null_sampled"] = est.to_dict()
            row["null_sampled_text"] = est.describe()
        rows.append(row)
    rows.sort(key=lambda r: -(r["share_among_completions"]["value"] if r.get("share_among_completions") else 0))
    results["observed_cycle_structures_under_null"] = rows
    results["observed_cycle_structure"] = str(c_obs) if c_obs is not None else None

    results["null_cycle_distribution"] = distribution_rows(
        ref.cycle_distribution, sample.cycle_structures if sample else None, "cycle_structure"
    )
    results["joint_observed_rejoins"] = _joint(problem, space, ref, completions, sample, confidence)
    return results


ANALYSES = {
    "distribution": run_distribution,
    "observed": run_observed,
    "rejoin_probability": run_rejoin_probability,
}
