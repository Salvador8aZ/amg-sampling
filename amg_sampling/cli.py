"""Command-line entry point ``amg-sampling``.

* ``amg-sampling [hydra overrides...]`` runs an experiment (see ``conf/``).
* ``amg-sampling patients [--data PATH] [--patient ID]`` lists the Sheth patients.
* ``amg-sampling import-data SOURCE`` copies a local copy of the dataset into
  ``data/external/`` and checks it against the reference checksum.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

SUBCOMMANDS = ("patients", "import-data")


def main(argv: list[str] | None = None) -> int | None:
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] in SUBCOMMANDS:
        return _subcommand(argv)
    from amg_sampling.app import hydra_main

    sys.argv = [sys.argv[0], *argv]
    hydra_main()
    return 0


def _subcommand(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="amg-sampling")
    commands = parser.add_subparsers(dest="command", required=True)
    patients = commands.add_parser("patients", help="list patients in the Sheth dataset")
    patients.add_argument("--data", type=Path, default=None, help="path to the dataset CSV")
    patients.add_argument("--patient", default=None, help="show the chromosome components of one patient")
    importer = commands.add_parser("import-data", help="copy a local dataset file into data/external/")
    importer.add_argument("source", type=Path)
    importer.add_argument("--dest", type=Path, default=Path("data/external/nihms.csv"))
    args = parser.parse_args(argv)
    try:
        if args.command == "patients":
            return _patients(args.data, args.patient)
        return _import(args.source, args.dest)
    except (FileNotFoundError, ValueError) as error:
        print(f"amg-sampling: error: {error}", file=sys.stderr)
        return 2


def _patients(data: Path | None, patient_id: str | None) -> int:
    from amg_sampling.data.sheth import ConversionError, chromosome_components, convert, default_data_path, load_dataset
    from amg_sampling.exact.formulas import num_completions

    path = data or default_data_path()
    dataset = load_dataset(path)
    source = dataset.source
    tag = "reference copy" if source.is_reference_copy else "differs from the reference copy"
    print(f"{source.file_name}: {source.rows:,} junctions, {len(dataset.patients)} patients ({tag})")

    if patient_id is not None:
        if patient_id not in dataset.patients:
            raise ValueError(f"Unknown patient {patient_id!r}.")
        patient = dataset.patients[patient_id]
        print(f"\n{patient_id}: {len(patient.variants)} junctions on chromosomes {', '.join(map(str, patient.chromosomes))}")
        print(f"{'chromosomes':<18}{'Θ':<22}{'observed':>9}{'unmatched':>10}{'completions':>14}  status")
        for component in chromosome_components(patient):
            chroms = ",".join(map(str, component))
            try:
                pc = convert(patient, component)
            except ConversionError as error:
                print(f"{chroms:<18}{'—':<22}{'':>9}{'':>10}{'':>14}  not convertible: {error}")
                continue
            completions = num_completions(len(pc.free_ends))
            status = "unique completion" if len(pc.free_ends) <= 2 else "ok"
            if pc.notes:
                status += f" ({len(pc.notes)} note(s))"
            print(
                f"{chroms:<18}{pc.describe_theta():<22}{len(pc.observed):>9}{len(pc.free_ends):>10}"
                f"{completions:>14,}  {status}"
            )
        first = ",".join(map(str, chromosome_components(patient)[0]))
        print(f"\nExample: amg-sampling problem=patient problem.patient_id={patient_id} "
              f"problem.chromosomes='[{first}]' analysis=observed")
        return 0

    print(f"\n{'patient':<12}{'junctions':>10}{'DSBs':>7}{'components':>11}  largest component (chromosomes: DSBs, unmatched ends)")
    for pid, patient in dataset.patients.items():
        components = chromosome_components(patient)
        dsbs = len({(b.chromosome, b.position) for v in patient.variants for b in v.breakpoints})
        largest = components[0]
        try:
            pc = convert(patient, largest)
            detail = f"{','.join(map(str, largest))}: {pc.theta.num_dsbs}, {len(pc.free_ends)}"
        except ConversionError:
            detail = f"{','.join(map(str, largest))}: not convertible (an end is used twice)"
        if len(detail) > 60:
            detail = detail[:57] + "..."
        print(f"{pid:<12}{len(patient.variants):>10}{dsbs:>7}{len(components):>11}  {detail}")
    print("\nDetails of one patient: amg-sampling patients --patient P05-1657")
    return 0


def _import(source: Path, dest: Path) -> int:
    from amg_sampling.data.sheth import REFERENCE_SHA256, load_dataset
    from amg_sampling.data.sheth.loader import file_sha256

    if not source.is_file():
        raise FileNotFoundError(f"{source} does not exist.")
    load_dataset(source)  # validates the format before copying
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, dest)
    sha = file_sha256(dest)
    match = "matches" if sha == REFERENCE_SHA256 else "does NOT match"
    print(f"Copied to {dest} (sha256 {sha[:12]}…, {match} the reference copy used by Sheth et al.)")
    return 0
