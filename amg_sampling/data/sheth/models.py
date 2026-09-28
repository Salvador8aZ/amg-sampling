"""Typed records for the structural-variant table used by Sheth et al. (2026).

The table is the prostate-cancer rearrangement data of Baca et al., *Punctuated
evolution of prostate cancer genomes*, Cell 153:666–677 (2013), as distributed
in the reference repository ``siddharthsheth/aberration_multigraph`` under
``data/nihms.csv``. It is **not** included in this repository (see
``docs/data.md``).

Each CSV row is one rearrangement junction joining two breakpoints. Only
fields that are present in the source are stored here:

============================  ===========================================
Field                         Source column
============================  ===========================================
``PatientRecord.patient_id``  ``Individual``
``source_line``               line number in the CSV file (header = line 1)
``number``                    ``Number``
``breakpoint_1.chromosome``   ``Breakpoint 1 chromosome`` (integer; the source
                              uses 1–24, with 23/24 presumably X/Y, which is
                              not asserted here)
``breakpoint_1.strand``       ``Breakpoint 1 strand`` (``+`` or ``-``)
``breakpoint_1.position``     ``Breakpoint 1 position`` (integer base pair)
``breakpoint_2.*``            the corresponding ``Breakpoint 2`` columns
``sv_class``                  ``Class`` (e.g. ``inter_chr``, ``inversion``)
============================  ===========================================

Other columns (read support, annotations, validation results) are not used by
the combinatorial model and are not stored.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Breakpoint:
    """One side of a junction as reported in the source table."""

    chromosome: int
    position: int
    strand: str  # "+" or "-", exactly as in the source


@dataclass(frozen=True, slots=True)
class StructuralVariant:
    """One row of the source table: a junction between two breakpoints."""

    source_line: int
    number: int
    breakpoint_1: Breakpoint
    breakpoint_2: Breakpoint
    sv_class: str

    @property
    def breakpoints(self) -> tuple[Breakpoint, Breakpoint]:
        return (self.breakpoint_1, self.breakpoint_2)


@dataclass(frozen=True, slots=True)
class PatientRecord:
    """All junctions reported for one individual, in source order."""

    patient_id: str
    variants: tuple[StructuralVariant, ...]

    @property
    def chromosomes(self) -> tuple[int, ...]:
        return tuple(sorted({b.chromosome for v in self.variants for b in v.breakpoints}))


@dataclass(frozen=True, slots=True)
class SourceInfo:
    """Where a dataset was read from. ``file_name`` omits the directory on purpose."""

    file_name: str
    sha256: str
    rows: int
    is_reference_copy: bool  # sha256 equals the copy used by the Sheth repository


@dataclass(frozen=True, slots=True)
class ShethDataset:
    source: SourceInfo
    patients: dict[str, PatientRecord]  # in order of first appearance
