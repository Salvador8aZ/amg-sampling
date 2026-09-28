"""Analyses: estimates, exact references, observed/rejoin probabilities, result files."""

import csv
import json
import math
from fractions import Fraction

import pytest

from amg_sampling import directories
from amg_sampling.analysis.estimates import ProportionEstimate, format_probability
from amg_sampling.analysis.output import write_results
from amg_sampling.analysis.plots import make_plots
from amg_sampling.analysis.reference import exact_reference
from amg_sampling.analysis.report import format_report
from amg_sampling.analysis.runs import (
    AnalysisError,
    Problem,
    SamplerSettings,
    run_distribution,
    run_observed,
    run_rejoin_probability,
)
from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.statespace import StateSpace
from amg_sampling.data.sheth import convert, load_dataset
from amg_sampling.exact.enumerate import exact_cycle_distribution, iter_states

SYNTHETIC = directories.test_data("synthetic_sheth.csv")
PROPER = StateSpace.PROPER


def patient_problem(patient_id="SYN-1", chromosomes=(1, 2)):
    ds = load_dataset(SYNTHETIC)
    pc = convert(ds.patients[patient_id], chromosomes)
    return Problem("patient", pc.theta, patient_id, pc, ds.source)


# -- estimates ----------------------------------------------------------------------


def test_proportion_estimate():
    est = ProportionEstimate(250, 1000)
    assert est.estimate == 0.25
    assert est.standard_error == pytest.approx(math.sqrt(0.25 * 0.75 / 1000))
    low, high = est.wilson_interval
    assert low < 0.25 < high
    assert est.zero_hit_upper_bound is None
    assert est.describe() == "0.25 ± 0.01369"


def test_zero_hits_are_never_reported_as_probability_zero():
    est = ProportionEstimate(0, 100_000)
    assert est.zero_hit_upper_bound == pytest.approx(1 - 0.05 ** (1 / 100_000))
    assert est.zero_hit_upper_bound == pytest.approx(3 / 100_000, rel=0.01)  # "rule of three"
    text = est.describe()
    assert text.startswith("0 occurrences in 100,000 samples")
    assert "upper bound 2.996e-05" in text
    assert est.wilson_interval[0] == 0.0


def test_few_hits_report_count_and_interval():
    text = ProportionEstimate(3, 10_000).describe()
    assert text.startswith("3 occurrences in 10,000 samples (95% interval")


@pytest.mark.parametrize("bad", [(1, 0), (-1, 5), (6, 5)])
def test_invalid_estimates(bad):
    with pytest.raises(ValueError):
        ProportionEstimate(*bad)


def test_format_probability():
    assert format_probability(None) == "n/a"
    assert format_probability(0) == "0"
    assert format_probability(0.123456) == "0.1235"
    assert format_probability(2.441e-6) == "2.441e-06"


# -- exact references -----------------------------------------------------------------


def test_exact_reference_small_proper():
    theta = InitialConfiguration((3, 1))
    ref = exact_reference(theta, PROPER)
    assert ref.cycle_distribution == exact_cycle_distribution(theta, PROPER)
    assert ref.total_states == sum(ref.cycle_distribution.values())


def test_exact_reference_large_proper_uses_bounded_deranged_summaries():
    ref = exact_reference(InitialConfiguration((20,) * 5), PROPER)
    assert ref.cycle_distribution is None and "cycle_distribution" in ref.unavailable
    assert ref.summaries_space == "deranged"
    assert 0 < ref.tv_bound < Fraction(1, 10**20)


def test_exact_reference_refuses_infeasible_proper_count():
    ref = exact_reference(InitialConfiguration((2,) * 20), PROPER)
    assert ref.total_states is None and "total_states" in ref.unavailable


# -- analyses -------------------------------------------------------------------------


def test_distribution_analysis_matches_exact_and_is_deterministic():
    problem = Problem("configuration", InitialConfiguration((3, 2)), "Θ(2,(3,2))")
    a = run_distribution(problem, PROPER, SamplerSettings("iid", 3000), seed=5)
    b = run_distribution(problem, PROPER, SamplerSettings("iid", 3000), seed=5)
    assert a == b
    assert a["sampling"]["samples"] == 3000
    rows = {r["cycle_structure"]: r for r in a["cycle_distribution"]}
    assert set(rows) == {"C5", "C3+C2"}
    assert rows["C5"]["exact_probability"]["fraction"] == "8/11"  # 384 / 528 (paper, table 3)
    for r in rows.values():
        assert abs(r["sampled_probability"] - r["exact_probability"]["value"]) < 5 * r["standard_error"]


def test_exact_sampler_needs_an_exact_result():
    problem = Problem("configuration", InitialConfiguration((2,) * 20), "big")
    with pytest.raises(AnalysisError):
        run_distribution(problem, PROPER, SamplerSettings("exact"), seed=0)


def test_observed_analysis_against_enumeration():
    problem = patient_problem()  # Θ(2,(3,1)), 2 observed rejoins, 4 unmatched ends
    pc = problem.patient
    results = run_observed(problem, PROPER, SamplerSettings("iid", 20_000), seed=11)
    proper = list(iter_states(pc.theta, PROPER))
    consistent = [r for r in proper if pc.fixed.is_completed_by(r)]

    obs = results["observed"]
    assert obs["num_unmatched_ends"] == 4 and obs["num_completions_all"] == 3
    assert obs["completions"]["exact_count"] == len(consistent)
    assert obs["reconstructed"] is None

    joint = results["joint_observed_rejoins"]
    assert Fraction(joint["exact"]["fraction"]) == Fraction(len(consistent), len(proper))

    # Rejoin probabilities: sampled vs exact P(e in R) by enumeration.
    for row in results["rejoins"]:
        exact = sum(r[row["end_a"]] == row["end_b"] for r in proper) / len(proper)
        est = row["null_probability"]
        assert abs(est["estimate"] - exact) < 5 * math.sqrt(exact * (1 - exact) / est["samples"])

    # Observed cycle structures: exact null probability agrees with enumeration.
    null = exact_cycle_distribution(pc.theta, PROPER)
    for row in results["observed_cycle_structures_under_null"]:
        (c,) = [c for c in null if str(c) == row["cycle_structure"]]
        assert Fraction(row["null_exact"]["fraction"]) == Fraction(null[c], len(proper))


def test_observed_analysis_with_unique_completion():
    results = run_observed(patient_problem("SYN-2", [5]), PROPER, SamplerSettings("iid", 2000), seed=1)
    rec = results["observed"]["reconstructed"]
    assert rec["cycle_structure"] == "C2"
    assert rec["membership"] == {"all": True, "deranged": True, "proper": True}
    assert [e["status"] for e in rec["rejoins"]] == ["observed junction", "reconstructed (unique completion)"]
    assert results["observed_cycle_structure"] == "C2"


def test_joint_zero_hit_handling():
    # P(both observed rejoins | uniform ALL) = 3/105; with 5 samples and this seed there is no hit.
    problem = patient_problem()
    results = run_rejoin_probability(problem, StateSpace.ALL, SamplerSettings("iid", 5), seed=3)
    joint = results["joint_observed_rejoins"]
    assert Fraction(joint["exact"]["fraction"]) == Fraction(3, 105)
    assert joint["sampled"]["hits"] == 0
    assert joint["sampled"]["estimate"] == 0
    assert joint["sampled"]["zero_hit_upper_bound"] == pytest.approx(1 - 0.05 ** (1 / 5))
    assert joint["sampled_text"].startswith("0 occurrences in 5 samples")


def test_observed_needs_a_patient():
    problem = Problem("configuration", InitialConfiguration((2, 2)), "x")
    with pytest.raises(AnalysisError, match="problem=patient"):
        run_observed(problem, PROPER, SamplerSettings("iid", 10), seed=0)
    with pytest.raises(AnalysisError):
        run_rejoin_probability(patient_problem(), PROPER, SamplerSettings("exact"), seed=0)


# -- results files, report, plots --------------------------------------------------------


@pytest.fixture(scope="module")
def observed_results():
    return run_observed(patient_problem(), PROPER, SamplerSettings("iid", 2000), seed=2)


def test_results_json_schema(tmp_path, observed_results):
    write_results(observed_results, tmp_path)
    data = json.loads((tmp_path / "results.json").read_text(encoding="utf-8"))
    assert data["schema_version"] == 1
    assert data["analysis"] == "observed"
    assert data["state_space"] == "proper"
    assert data["null_model"] == "Uniform(PROPER(Θ))"
    assert "not biological" in data["interpretation"]
    for key in ("problem", "sampler", "null", "observed", "rejoins", "joint_observed_rejoins",
                "observed_cycle_structures_under_null", "null_cycle_distribution"):
        assert key in data
    assert data["problem"]["theta"] == {
        "notation": "Θ(2,(3,1))", "num_chromosomes": 2, "breaks": [3, 1], "num_dsbs": 4, "num_free_ends": 8,
    }
    assert data["problem"]["data_source"]["file_name"] == "synthetic_sheth.csv"


def test_csv_tables(tmp_path, observed_results):
    files = {p.name for p in write_results(observed_results, tmp_path)}
    assert {"results.json", "rejoin_probabilities.csv", "observed_cycle_structures.csv",
            "null_cycle_distribution.csv", "completion_distribution.csv"} <= files
    with open(tmp_path / "rejoin_probabilities.csv", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    for column in ("edge", "site_a", "site_b", "status", "source_line", "null_probability_estimate",
                   "null_probability_standard_error", "null_probability_zero_hit_upper_bound"):
        assert column in rows[0]


def test_report_mentions_the_null_and_its_interpretation(observed_results):
    text = format_report(observed_results)
    assert "Patient:          SYN-1" in text
    assert "uniform combinatorial null" in text
    assert "not biological" in text


def test_plots(tmp_path, observed_results):
    pytest.importorskip("matplotlib")
    written = {p.name for p in make_plots(observed_results, tmp_path)}
    assert written == {"observed_cycle_structures.png", "rejoin_probabilities.png"}
    assert all((tmp_path / name).stat().st_size > 0 for name in written)
