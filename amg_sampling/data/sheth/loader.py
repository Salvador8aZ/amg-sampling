"""Read the Sheth / Baca et al. structural-variant CSV into :mod:`models` records."""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path

from amg_sampling import directories
from amg_sampling.data.sheth.models import (
    Breakpoint,
    PatientRecord,
    ShethDataset,
    SourceInfo,
    StructuralVariant,
)

# SHA-256 of data/nihms.csv in siddharthsheth/aberration_multigraph
# (commit 188e107, "tidied data"), the file used for the paper's case study.
REFERENCE_SHA256 = "13a70c63d09af08b423fe372fd09a4d6c1299bb9ad5a1424e9cc490af1a6a3bc"

DATA_ENV_VAR = "AMG_SHETH_DATA"
DEFAULT_RELATIVE_PATH = Path("data/external/nihms.csv")

REQUIRED_COLUMNS = (
    "Individual",
    "Number",
    "Breakpoint 1 chromosome",
    "Breakpoint 1 strand",
    "Breakpoint 1 position",
    "Breakpoint 2 chromosome",
    "Breakpoint 2 strand",
    "Breakpoint 2 position",
    "Class",
)

# Optional junction-sequence features. -1 in the source means the junction
# sequence could not be assembled; it, NaN and empty cells are read as None
# ("not measured").
OPTIONAL_LENGTH_COLUMNS = {
    "homology_length": "Homology length",
    "foreign_sequence_length": "Foreign sequence length",
}


class DatasetFormatError(ValueError):
    """The file does not have the expected columns or values."""


def default_data_path() -> Path:
    """``$AMG_SHETH_DATA`` if set, else ``data/external/nihms.csv``, resolved
    by :func:`resolve_data_path`.
    """
    return resolve_data_path(os.environ.get(DATA_ENV_VAR, DEFAULT_RELATIVE_PATH))


def resolve_data_path(
    path: str | os.PathLike, start: str | os.PathLike | None = None
) -> Path:
    """A relative ``path`` is taken from ``start`` (default: the working
    directory) if it exists there, and otherwise from the repository root.
    Runs started in a subdirectory (an IDE run configuration, a notebook) then
    still find ``data/external/``.
    """
    path = Path(path)
    if path.is_absolute():
        return path
    local = Path(start) / path if start is not None else path
    if local.exists():
        return local
    return directories.base(path)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 16), b""):
            digest.update(block)
    return digest.hexdigest()


def load_dataset(path: str | os.PathLike) -> ShethDataset:
    """Parse the CSV. Raises ``FileNotFoundError`` or :class:`DatasetFormatError`."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"Sheth dataset not found at {path}. See docs/data.md for how to obtain it "
            f"(or set {DATA_ENV_VAR})."
        )
    patients: dict[str, list[StructuralVariant]] = {}
    rows = 0
    # utf-8-sig strips the byte-order mark present in the reference copy.
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise DatasetFormatError(f"{path.name}: missing columns {missing}.")
        for line, row in enumerate(reader, start=2):
            rows += 1
            variant = StructuralVariant(
                source_line=line,
                number=_integer(row, "Number", line),
                breakpoint_1=_breakpoint(row, 1, line),
                breakpoint_2=_breakpoint(row, 2, line),
                sv_class=row["Class"].strip(),
                **{
                    field: _optional_length(row, column, line)
                    for field, column in OPTIONAL_LENGTH_COLUMNS.items()
                },
            )
            patients.setdefault(row["Individual"].strip(), []).append(variant)
    sha = file_sha256(path)
    return ShethDataset(
        source=SourceInfo(path.name, sha, rows, sha == REFERENCE_SHA256),
        patients={pid: PatientRecord(pid, tuple(v)) for pid, v in patients.items()},
    )


def _integer(row: dict, column: str, line: int) -> int:
    value = (row.get(column) or "").strip()
    try:
        return int(value)
    except ValueError:
        raise DatasetFormatError(
            f"line {line}: column {column!r} is not an integer: {value!r}."
        ) from None


def _optional_length(row: dict, column: str, line: int) -> int | None:
    """A non-negative length, or ``None`` if the column is absent, empty, NaN or -1."""
    value = (row.get(column) or "").strip()
    if value in ("", "-1") or value.lower() == "nan":
        return None
    try:
        length = int(value)
    except ValueError:
        raise DatasetFormatError(
            f"line {line}: column {column!r} is not an integer: {value!r}."
        ) from None
    if length < 0:
        raise DatasetFormatError(
            f"line {line}: column {column!r} must be -1 or non-negative, got {length}."
        )
    return length


def _breakpoint(row: dict, which: int, line: int) -> Breakpoint:
    prefix = f"Breakpoint {which}"
    strand = (row.get(f"{prefix} strand") or "").strip()
    if strand not in ("+", "-"):
        raise DatasetFormatError(
            f"line {line}: {prefix} strand must be '+' or '-', got {strand!r}."
        )
    return Breakpoint(
        chromosome=_integer(row, f"{prefix} chromosome", line),
        position=_integer(row, f"{prefix} position", line),
        strand=strand,
    )
