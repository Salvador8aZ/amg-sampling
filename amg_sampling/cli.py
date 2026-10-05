"""Command-line entry point ``amg-sampling``.

* ``amg-sampling [hydra overrides...]`` runs an experiment (see ``conf/``).
* ``amg-sampling patients [--data PATH] [--patient ID]`` lists the Sheth patients.
* ``amg-sampling import-data SOURCE`` copies a local copy of the dataset into
  ``data/external/`` and checks it against the reference checksum.
* ``amg-sampling fetch-data`` downloads the dataset from its source into
  ``data/external/`` (experiments also do this on first use; see ``conf/``).
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from amg_sampling import directories

SUBCOMMANDS = ("patients", "import-data", "fetch-data")


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
    patients = commands.add_parser(
        "patients", help="list patients in the Sheth dataset"
    )
    patients.add_argument(
        "--data", type=Path, default=None, help="path to the dataset CSV"
    )
    patients.add_argument(
        "--patient", default=None, help="show the chromosome components of one patient"
    )
    importer = commands.add_parser(
        "import-data", help="copy a local dataset file into data/external/"
    )
    importer.add_argument("source", type=Path)
    importer.add_argument(
        "--dest", type=Path, default=directories.data("external/nihms.csv")
    )
    fetcher = commands.add_parser(
        "fetch-data", help="download the dataset from its source into data/external/"
    )
    fetcher.add_argument(
        "--dest", type=Path, default=directories.data("external/nihms.csv")
    )
    args = parser.parse_args(argv)
    try:
        if args.command == "patients":
            return _patients(args.data, args.patient)
        if args.command == "fetch-data":
            return _fetch(args.dest)
        return _import(args.source, args.dest)
    except (OSError, ValueError) as error:
        print(f"amg-sampling: error: {error}", file=sys.stderr)
        return 2


def _patients(data: Path | None, patient_id: str | None) -> int:
    from amg_sampling.data.sheth import (
        ConversionError,
        chromosome_components,
        convert,
        default_data_path,
        load_dataset,
    )
    from amg_sampling.exact.formulas import num_completions

    path = data or default_data_path()
    dataset = load_dataset(path)
    source = dataset.source
    tag = (
        "reference copy"
        if source.is_reference_copy
        else "differs from the reference copy"
    )
    print(
        f"{source.file_name}: {source.rows:,} junctions, "
        f"{len(dataset.patients)} patients ({tag})"
    )

    if patient_id is not None:
        if patient_id not in dataset.patients:
            raise ValueError(f"Unknown patient {patient_id!r}.")
        patient = dataset.patients[patient_id]
        print(
            f"\n{patient_id}: {len(patient.variants)} junctions on chromosomes "
            f"{', '.join(map(str, patient.chromosomes))}"
        )
        rows = []
        for component in chromosome_components(patient):
            chroms = ",".join(map(str, component))
            try:
                pc = convert(patient, component)
            except ConversionError as error:
                rows.append((chroms, "—", "", "", "", f"not convertible: {error}"))
                continue
            status = "unique completion" if len(pc.free_ends) <= 2 else "ok"
            if pc.notes:
                status += f" ({len(pc.notes)} note(s))"
            rows.append(
                (
                    chroms,
                    pc.describe_theta(),
                    f"{len(pc.observed):,}",
                    f"{len(pc.free_ends):,}",
                    _count(num_completions(len(pc.free_ends))),
                    status,
                )
            )
        # Columns widen to fit long chromosome lists and Θ labels.
        w_chroms = max(len("chromosomes"), *(len(r[0]) for r in rows)) + 2
        w_theta = max(20, *(len(r[1]) for r in rows)) + 2
        print(
            f"{'chromosomes':<{w_chroms}}{'Θ':<{w_theta}}{'observed':>9}"
            f"{'unmatched':>10}{'completions':>14}  status"
        )
        for chroms, theta, observed, free, completions, status in rows:
            print(
                f"{chroms:<{w_chroms}}{theta:<{w_theta}}{observed:>9}{free:>10}"
                f"{completions:>14}  {status}"
            )
        first = ",".join(map(str, chromosome_components(patient)[0]))
        print(
            f"\nExample: amg-sampling problem=patient problem.patient_id={patient_id} "
            f"problem.chromosomes='[{first}]' analysis=observed"
        )
        return 0

    print(
        f"\n{'patient':<12}{'junctions':>10}{'DSBs':>7}{'components':>11}  "
        "largest component (chromosomes: DSBs, unmatched ends)"
    )
    for pid, patient in dataset.patients.items():
        components = chromosome_components(patient)
        dsbs = len(
            {
                (b.chromosome, b.position)
                for v in patient.variants
                for b in v.breakpoints
            }
        )
        largest = components[0]
        try:
            pc = convert(patient, largest)
            detail = (
                f"{','.join(map(str, largest))}: {pc.theta.num_dsbs}, "
                f"{len(pc.free_ends)}"
            )
        except ConversionError:
            detail = (
                f"{','.join(map(str, largest))}: not convertible (an end is used twice)"
            )
        if len(detail) > 60:
            detail = detail[:57] + "..."
        print(
            f"{pid:<12}{len(patient.variants):>10}{dsbs:>7}"
            f"{len(components):>11}  {detail}"
        )
    print("\nDetails of one patient: amg-sampling patients --patient P05-1657")
    return 0


def _count(value: int) -> str:
    """Exact with thousands separators up to 10^12, else ``1.234e+701``.

    Completion counts ``(f−1)!!`` can have hundreds of digits, too large for a
    float, so the scientific form is built from the decimal digits.
    """
    if value < 10**12:
        return f"{value:,}"
    digits = str(value)
    return f"{digits[0]}.{digits[1:4]}e+{len(digits) - 1}"


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
    print(
        f"Copied to {dest} (sha256 {sha[:12]}…, {match} the reference copy "
        "used by Sheth et al.)"
    )
    return 0


def _fetch(dest: Path) -> int:
    from amg_sampling.data.sheth import load_dataset
    from amg_sampling.data.sheth.fetch import SOURCE_URL, fetch

    existed = dest.is_file()
    fetch(dest)
    load_dataset(dest)  # validates the format
    if existed:
        print(f"{dest} already exists; not downloaded again.")
    else:
        print(f"Downloaded {SOURCE_URL}\n  to {dest} (SHA-256 verified)")
    return 0
