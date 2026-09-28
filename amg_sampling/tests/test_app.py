"""Hydra configuration, the application runner, and command-line smoke tests."""

import json
import subprocess
import sys
from pathlib import Path

import pytest
from hydra import compose, initialize_config_module

from amg_sampling.app import build_problem, execute
from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.data.sheth import default_data_path

SYNTHETIC = Path(__file__).parent / "test_data" / "synthetic_sheth.csv"
REAL = default_data_path()


def config(*overrides):
    with initialize_config_module(version_base=None, config_module="amg_sampling.conf"):
        return compose(config_name="config", overrides=list(overrides))


# -- configuration composition ---------------------------------------------------------


def test_defaults():
    cfg = config()
    assert cfg.problem.kind == "configuration"
    assert list(cfg.problem.breaks) == [4, 4, 4]
    assert cfg.statespace.name == "proper"
    assert cfg.sampler.name == "iid" and cfg.sampler.num_samples == 100_000
    assert cfg.analysis.name == "distribution"
    assert cfg.seed == 42


@pytest.mark.parametrize(
    "group, option, key, value",
    [
        ("statespace", "all", "name", "all"),
        ("statespace", "deranged", "name", "deranged"),
        ("sampler", "exact", "name", "exact"),
        ("analysis", "observed", "name", "observed"),
        ("analysis", "rejoin_probability", "confidence", 0.95),
        ("problem", "patient", "patient_id", "P05-1657"),
    ],
)
def test_groups_are_independent(group, option, key, value):
    cfg = config(f"{group}={option}")
    assert cfg[group][key] == value
    # Other groups keep their defaults.
    assert cfg.seed == 42
    if group != "sampler":
        assert cfg.sampler.name == "iid"


def test_synthetic_configuration_parsing():
    cfg = config("problem.breaks=[20,20,20,20,20]")
    problem = build_problem(cfg)
    assert problem.theta == InitialConfiguration((20, 20, 20, 20, 20))
    assert problem.theta.num_dsbs == 100 and problem.theta.num_ends == 200
    with pytest.raises(ValueError):
        build_problem(config("problem.breaks=[3,0]"))


def test_patient_problem_from_synthetic_data():
    cfg = config("problem=patient", "problem.patient_id=SYN-1", "problem.chromosomes=[1,2]",
                 f"data.sheth_csv={SYNTHETIC}")
    problem = build_problem(cfg)
    assert problem.kind == "patient" and problem.theta.breaks == (3, 1)
    cfg = config("problem=patient", "problem.patient_id=SYN-1", "problem.chromosomes=null",
                 f"data.sheth_csv={SYNTHETIC}")
    assert build_problem(cfg).patient.chromosomes == (1, 2)  # only one component


# -- execute (no Hydra runtime needed) -------------------------------------------------------


def test_execute_writes_files_and_is_deterministic(tmp_path):
    cfg = config("problem.breaks=[3,2]", "sampler.num_samples=500", "output.plots=false")
    results_a, files = execute(cfg, tmp_path / "a")
    results_b, _ = execute(cfg, tmp_path / "b")
    assert results_a == results_b
    names = {p.name for p in files}
    assert {"config.yaml", "results.json", "cycle_distribution.csv", "cycle_count_distribution.csv",
            "largest_cycle_distribution.csv"} <= names
    assert json.loads((tmp_path / "a" / "results.json").read_text(encoding="utf-8")) == json.loads(
        (tmp_path / "b" / "results.json").read_text(encoding="utf-8")
    )
    assert "breaks:\n  - 3\n  - 2" in (tmp_path / "a" / "config.yaml").read_text(encoding="utf-8")


def test_execute_rejects_bad_seed(tmp_path):
    with pytest.raises(ValueError, match="seed"):
        execute(config("seed=abc", "sampler.num_samples=10"), tmp_path)


@pytest.mark.skipif(not REAL.is_file(), reason="Sheth dataset not available (see docs/data.md)")
def test_demonstration_patient_exact(tmp_path):
    cfg = config("problem=patient", "analysis=observed", "sampler=exact", "output.plots=false",
                 f"data.sheth_csv={REAL.resolve()}")
    results, _ = execute(cfg, tmp_path)
    assert results["observed"]["completions"]["exact_count"] == 945
    joint = results["joint_observed_rejoins"]
    assert (joint["exact_numerator"], joint["exact_denominator"]) == (945, 387_099_936)


# -- command line ---------------------------------------------------------------------------


def run_cli(*args, cwd):
    return subprocess.run(
        [sys.executable, "-m", "amg_sampling", *args], cwd=cwd, capture_output=True, text=True, timeout=300
    )


def test_cli_configuration_run(tmp_path):
    run_dir = tmp_path / "run"
    out = run_cli("problem.breaks=[3,3]", "sampler.num_samples=300", "output.plots=false",
                  f"hydra.run.dir={run_dir}", cwd=tmp_path)
    assert out.returncode == 0, out.stderr
    assert "AMG Sampling — Configuration analysis: distribution" in out.stdout
    assert "Θ(2,(3,3))" in out.stdout
    assert "not biological" in out.stdout
    assert (run_dir / "results.json").is_file() and (run_dir / "config.yaml").is_file()


def test_cli_patient_run_on_synthetic_data(tmp_path):
    run_dir = tmp_path / "run"
    out = run_cli("problem=patient", "problem.patient_id=SYN-1", "problem.chromosomes=[1,2]",
                  f"data.sheth_csv={SYNTHETIC}", "analysis=observed", "sampler.num_samples=300",
                  f"hydra.run.dir={run_dir}", cwd=tmp_path)
    assert out.returncode == 0, out.stderr
    assert "Patient:          SYN-1" in out.stdout
    assert (run_dir / "rejoin_probabilities.csv").is_file()
    assert (run_dir / "rejoin_probabilities.png").is_file()


def test_cli_reports_user_errors(tmp_path):
    out = run_cli("analysis=observed", "sampler.num_samples=10", f"hydra.run.dir={tmp_path / 'r'}", cwd=tmp_path)
    assert out.returncode == 2
    assert "needs observed rejoins" in out.stderr


def test_cli_patients(tmp_path):
    out = run_cli("patients", "--data", str(SYNTHETIC), cwd=tmp_path)
    assert out.returncode == 0, out.stderr
    assert "4 patients" in out.stdout and "SYN-1" in out.stdout
    out = run_cli("patients", "--data", str(SYNTHETIC), "--patient", "SYN-3", cwd=tmp_path)
    assert out.returncode == 0 and "not convertible" in out.stdout
    out = run_cli("patients", "--data", str(tmp_path / "missing.csv"), cwd=tmp_path)
    assert out.returncode == 2 and "docs/data.md" in out.stderr


def test_cli_import_data(tmp_path):
    out = run_cli("import-data", str(SYNTHETIC), "--dest", str(tmp_path / "copy.csv"), cwd=tmp_path)
    assert out.returncode == 0, out.stderr
    assert "does NOT match the reference copy" in out.stdout
    assert (tmp_path / "copy.csv").read_bytes() == SYNTHETIC.read_bytes()
