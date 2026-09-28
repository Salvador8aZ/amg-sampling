"""Hydra application: configuration → problem → analysis → report and files.

Hydra is used only here. ``execute`` takes a composed configuration and an
output directory, so it can be called without Hydra's runtime (as the tests do).
"""

from __future__ import annotations

import sys
from pathlib import Path

import hydra
from hydra.core.hydra_config import HydraConfig
from hydra.utils import to_absolute_path
from omegaconf import DictConfig, OmegaConf

from amg_sampling.analysis.output import write_results
from amg_sampling.analysis.plots import make_plots, matplotlib_available
from amg_sampling.analysis.report import format_report
from amg_sampling.analysis.runs import (
    ANALYSES,
    AnalysisError,
    Problem,
    SamplerSettings,
    theta_label,
)
from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.statespace import StateSpace
from amg_sampling.data.sheth import (
    ConversionError,
    DatasetFormatError,
    chromosome_components,
    convert,
    load_dataset,
)

USER_ERRORS = (
    AnalysisError,
    ConversionError,
    DatasetFormatError,
    FileNotFoundError,
    ValueError,
    TypeError,
)


def build_problem(cfg: DictConfig) -> Problem:
    kind = cfg.problem.kind
    if kind == "configuration":
        breaks = list(cfg.problem.breaks)
        theta = InitialConfiguration([int(b) for b in breaks])
        return Problem(
            "configuration",
            theta,
            f"{theta_label(theta)}: {theta.num_dsbs} DSBs on "
            f"{theta.num_chromosomes} chromosomes",
        )
    if kind == "patient":
        path = Path(to_absolute_path(str(cfg.data.sheth_csv)))
        dataset = load_dataset(path)
        patient_id = str(cfg.problem.patient_id)
        if patient_id not in dataset.patients:
            raise ValueError(
                f"Unknown patient {patient_id!r}; run 'amg-sampling patients' "
                "to list them."
            )
        patient = dataset.patients[patient_id]
        chromosomes = cfg.problem.chromosomes
        if chromosomes is None:
            components = chromosome_components(patient)
            if len(components) > 1:
                raise ValueError(
                    f"{patient_id} has several components {list(components)}; "
                    "set problem.chromosomes "
                    "(e.g. problem.chromosomes="
                    f"'[{','.join(map(str, components[0]))}]')."
                )
            chromosomes = components[0]
        pc = convert(patient, [int(c) for c in chromosomes])
        label = f"{patient_id} chromosomes {', '.join(map(str, pc.chromosomes))}"
        return Problem("patient", pc.theta, label, pc, dataset.source)
    raise ValueError(f"Unknown problem kind {kind!r}.")


def execute(cfg: DictConfig, output_dir: Path) -> tuple[dict, list[Path]]:
    """Run the configured analysis, write result files, and return (results, files)."""
    seed = cfg.seed
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError(f"seed must be an integer, got {seed!r}.")
    problem = build_problem(cfg)
    space = StateSpace(cfg.statespace.name)
    sampler = SamplerSettings(
        name=cfg.sampler.name,
        num_samples=int(cfg.sampler.get("num_samples", 0) or 0),
        max_proposals_per_sample=cfg.sampler.get("max_proposals_per_sample"),
        burn_in=int(cfg.sampler.get("burn_in", 0) or 0),
        thin=int(cfg.sampler.get("thin", 1) or 1),
    )
    if sampler.samples and sampler.num_samples <= 0:
        raise ValueError("sampler.num_samples must be positive.")
    analysis = cfg.analysis.name
    if analysis not in ANALYSES:
        raise ValueError(
            f"Unknown analysis {analysis!r}; choose one of {sorted(ANALYSES)}."
        )
    kwargs = (
        {"confidence": float(cfg.analysis.confidence)}
        if "confidence" in cfg.analysis
        else {}
    )
    results = ANALYSES[analysis](problem, space, sampler, seed, **kwargs)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "config.yaml").write_text(
        OmegaConf.to_yaml(cfg, resolve=True), encoding="utf-8"
    )
    files = [output_dir / "config.yaml", *write_results(results, output_dir)]
    if cfg.output.plots:
        files += make_plots(results, output_dir)
    return results, files


@hydra.main(version_base=None, config_path="conf", config_name="config")
def hydra_main(cfg: DictConfig) -> None:
    output_dir = Path(HydraConfig.get().runtime.output_dir)
    try:
        results, files = execute(cfg, output_dir)
    except USER_ERRORS as error:
        print(f"amg-sampling: error: {error}", file=sys.stderr)
        sys.exit(2)
    print(format_report(results, top=int(cfg.output.top)))
    print()
    print(f"Results written to {_display_path(output_dir)}/")
    for path in files:
        print(f"  {path.name}")
    if cfg.output.plots and not matplotlib_available():
        print(
            "  (figures skipped: install the 'plots' extra, "
            "e.g. `poetry install --extras plots`)"
        )


def _display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(Path.cwd().resolve()))
    except ValueError:
        return str(path)
