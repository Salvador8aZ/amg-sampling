"""Convert a patient's junctions into an initial configuration and observed rejoins.

What is observed and what is reconstructed
------------------------------------------
Observed in the source table (``models``): each junction's two breakpoints,
given as (chromosome, position, strand).

Reconstructed here, following the Sheth methodology as implemented in the
reference repository (``examples/patient_analysis/nihms_patient.py``):

1. **Chromosome set.** A patient is analysed one *component* at a time: a set
   of chromosomes linked by observed junctions. No exchange is assumed between
   chromosomes that no observed junction links (paper, section 7.1).
2. **DSBs.** Every distinct breakpoint (chromosome, position) in the component
   is a DSB. DSBs are ordered by chromosome, then position. DSB ``j`` has a left end
   ``2j`` (joined by chromatin towards lower coordinates, labelled ``position``)
   and a right end ``2j+1`` (labelled ``position + 1``).
3. **Rejoin ends.** The strand sign says which side of the break is fused into
   the junction. Default (``convention="baca"``, the dRanger semantics of the
   source data, Baca et al. 2013): strand ``+`` means the DNA to the *left* of
   the break (lower coordinates) is joined, so the junction uses that DSB's
   left end; strand ``-`` uses its right end. This is why Baca's "potential
   deletion" class is always ``(+, -)``: the outer flanks are joined and the
   middle is lost, as the Table S3B copy number confirms (e.g. the ~3 Mb
   TMPRSS2-ERG deletion in P05-1657, junction #26).

   ``convention="sheth2026"`` is the reverse mapping (``+`` -> right end)
   used by Sheth et al. 2026 and their reference repository
   (``nihms_patient.py::_rejoin_edges``: ``if loc1 == '+': bp1 += 1``). It is
   kept only to reproduce the paper's figures. Flipping the convention moves
   every rejoin to the other end of its own DSB: cycle-structure counts are
   unchanged, but which ends are observed and which are free, and hence
   proximity weights and rejoin probabilities, change.
4. **Missing rejoins.** Ends not used by any junction are *free*. Their
   partners are unknown; completions are enumerated or sampled, never invented.
   When exactly two ends are free, the completion is unique. That edge is
   reported as *reconstructed (unique completion)*, never as observed.

Differences from the reference implementation (all make failures explicit):

* Breakpoints with identical (chromosome, position) are one DSB. The reference
  code would create duplicate DSB edges.
* A junction with one breakpoint outside the chosen chromosomes is an error.
  The reference code silently dropped such junctions.
* An end used by two junctions is an error: it cannot be a rejoin matching.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.matching import PartialMatching, RejoinMatching
from amg_sampling.data.sheth.models import Breakpoint, PatientRecord, StructuralVariant

# Which end of a DSB a strand sign selects: 0 = left end (2j, labelled
# ``position``), 1 = right end (2j+1, labelled ``position + 1``).
STRAND_CONVENTIONS: dict[str, dict[str, int]] = {
    "baca": {"+": 0, "-": 1},
    "sheth2026": {"+": 1, "-": 0},
}
DEFAULT_CONVENTION = "baca"


class ConversionError(ValueError):
    """The requested chromosomes cannot be represented as Θ plus a partial
    rejoin matching.
    """


class EdgeStatus(Enum):
    OBSERVED = "observed junction"
    UNIQUE_COMPLETION = "reconstructed (unique completion)"


@dataclass(frozen=True, slots=True)
class EndSite:
    """Genomic meaning of a free DSB end."""

    end: int
    chromosome: int
    position: int  # breakpoint position of the DSB
    side: str  # "left" (end 2j) or "right" (end 2j+1)

    @property
    def coordinate(self) -> int:
        """Vertex coordinate: position (left end) or position + 1 (right end)."""
        return self.position if self.side == "left" else self.position + 1

    def label(self) -> str:
        """Breakpoint position and side, e.g. ``chr8:40,932,539 (right)``."""
        return f"chr{self.chromosome}:{self.position:,} ({self.side})"


@dataclass(frozen=True, slots=True)
class RejoinEdge:
    end_a: int
    end_b: int
    status: EdgeStatus
    variant: StructuralVariant | None  # the source row, for observed junctions


@dataclass(frozen=True)
class PatientConfiguration:
    """Θ and the observed rejoins of one patient component."""

    patient_id: str
    chromosomes: tuple[int, ...]
    theta: InitialConfiguration
    ends: tuple[EndSite, ...]
    observed: tuple[RejoinEdge, ...]
    notes: tuple[str, ...]
    convention: str = DEFAULT_CONVENTION

    @property
    def fixed(self) -> PartialMatching:
        return PartialMatching(
            self.theta.num_ends, [(e.end_a, e.end_b) for e in self.observed]
        )

    @property
    def free_ends(self) -> tuple[int, ...]:
        return self.fixed.free_ends

    def reconstructed_matching(
        self,
    ) -> tuple[RejoinMatching, tuple[RejoinEdge, ...]] | None:
        """The unique completion when at most two ends are free, with its
        inferred edges.
        """
        free = self.free_ends
        if len(free) > 2:
            return None
        inferred = ()
        pairs = [(e.end_a, e.end_b) for e in self.observed]
        if free:
            inferred = (
                RejoinEdge(free[0], free[1], EdgeStatus.UNIQUE_COMPLETION, None),
            )
            pairs.append(free)
        return RejoinMatching.from_pairs(self.theta.num_ends, pairs), inferred

    def describe_theta(self) -> str:
        return (
            f"Θ({self.theta.num_chromosomes},({','.join(map(str, self.theta.breaks))}))"
        )


def chromosome_components(patient: PatientRecord) -> tuple[tuple[int, ...], ...]:
    """Sets of chromosomes linked by observed junctions, largest (by
    breakpoints) first.
    """
    parent: dict[int, int] = {}

    def find(x: int) -> int:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for v in patient.variants:
        a, b = find(v.breakpoint_1.chromosome), find(v.breakpoint_2.chromosome)
        parent[a] = b
    groups: dict[int, set[int]] = {}
    for c in patient.chromosomes:
        groups.setdefault(find(c), set()).add(c)
    sizes = {
        root: len(_distinct_breakpoints(patient, chroms))
        for root, chroms in groups.items()
    }
    ordered = sorted(groups, key=lambda root: (-sizes[root], min(groups[root])))
    return tuple(tuple(sorted(groups[root])) for root in ordered)


def convert(
    patient: PatientRecord, chromosomes, convention: str = DEFAULT_CONVENTION
) -> PatientConfiguration:
    """Θ, end sites and observed rejoins for ``chromosomes`` (a union of components).

    ``convention`` maps strand signs to DSB ends: ``"baca"`` (default, correct)
    or ``"sheth2026"`` (the paper's reversed mapping); see the module docstring.
    """
    if convention not in STRAND_CONVENTIONS:
        raise ConversionError(
            f"Unknown strand convention {convention!r}; "
            f"choose one of {sorted(STRAND_CONVENTIONS)}."
        )
    side_of = STRAND_CONVENTIONS[convention]
    chosen = tuple(sorted(set(int(c) for c in chromosomes)))
    if not chosen:
        raise ConversionError("No chromosomes selected.")
    unknown = [c for c in chosen if c not in patient.chromosomes]
    if unknown:
        raise ConversionError(
            f"{patient.patient_id} has no junction on chromosome(s) {unknown}."
        )
    inside = set(chosen)
    variants, crossing = [], []
    for v in patient.variants:
        n_in = sum(b.chromosome in inside for b in v.breakpoints)
        if n_in == 2:
            variants.append(v)
        elif n_in == 1:
            crossing.append(v)
    if crossing:
        lines = ", ".join(_ref(v) for v in crossing)
        raise ConversionError(
            f"Junctions {lines} join the selected chromosomes "
            "to others; "
            f"select a union of components: {chromosome_components(patient)}."
        )

    dsbs = _distinct_breakpoints(patient, inside)
    breaks = tuple(sum(1 for c, _ in dsbs if c == chrom) for chrom in chosen)
    theta = InitialConfiguration(breaks)
    index = {site: j for j, site in enumerate(dsbs)}
    ends = tuple(
        EndSite(2 * j + side, chrom, pos, "left" if side == 0 else "right")
        for j, (chrom, pos) in enumerate(dsbs)
        for side in (0, 1)
    )

    used: dict[int, StructuralVariant] = {}
    observed = []
    for v in variants:
        pair = tuple(_end_of(index, b, side_of) for b in v.breakpoints)
        for end in pair:
            if end in used:
                raise ConversionError(
                    f"End {ends[end].label()} is used by junctions "
                    f"{_ref(used[end])} and {_ref(v)}; "
                    "not a rejoin matching."
                )
            used[end] = v
        if pair[0] == pair[1]:
            raise ConversionError(f"Junction {_ref(v)} joins an end to itself.")
        observed.append(RejoinEdge(pair[0], pair[1], EdgeStatus.OBSERVED, v))

    notes = []
    positions = set(dsbs)
    for chrom, pos in dsbs:
        if (chrom, pos + 1) in positions:
            notes.append(
                f"Breakpoints chr{chrom}:{pos} and chr{chrom}:{pos + 1} are "
                "adjacent; they are "
                "kept as separate DSBs (the reference code would merge their vertices)."
            )
    return PatientConfiguration(
        patient.patient_id,
        chosen,
        theta,
        ends,
        tuple(observed),
        tuple(notes),
        convention,
    )


def _ref(v: StructuralVariant) -> str:
    """``11 (CSV line 315)``: junction Number, then the line of the local copy."""
    return f"{v.number} (CSV line {v.source_line})"


def _distinct_breakpoints(patient: PatientRecord, chromosomes) -> list[tuple[int, int]]:
    return sorted(
        {
            (b.chromosome, b.position)
            for v in patient.variants
            for b in v.breakpoints
            if b.chromosome in chromosomes
        }
    )


def _end_of(
    index: dict[tuple[int, int], int], b: Breakpoint, side_of: dict[str, int]
) -> int:
    j = index[(b.chromosome, b.position)]
    return 2 * j + side_of[b.strand]
