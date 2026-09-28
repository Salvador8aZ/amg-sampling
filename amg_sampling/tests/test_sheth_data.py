from pathlib import Path

import pytest

from amg_sampling.core.configuration import InitialConfiguration
from amg_sampling.core.cycles import cycle_structure
from amg_sampling.core.matching import PartialMatching
from amg_sampling.core.statespace import StateSpace
from amg_sampling.data.sheth import (
    REFERENCE_SHA256,
    ConversionError,
    DatasetFormatError,
    EdgeStatus,
    chromosome_components,
    convert,
    default_data_path,
    load_dataset,
)
from amg_sampling.exact.enumerate import count_states, exact_cycle_distribution
from amg_sampling.tests.paper_patient import P05_1657, TABLE_2

SYNTHETIC = Path(__file__).parent / "test_data" / "synthetic_sheth.csv"
REAL = default_data_path()
needs_real_data = pytest.mark.skipif(
    not REAL.is_file(), reason=f"Sheth dataset not available at {REAL} (see docs/data.md)"
)


def by_string(dist):
    return {str(c): m for c, m in dist.items()}


# -- loading --------------------------------------------------------------------


def test_load_synthetic_dataset():
    ds = load_dataset(SYNTHETIC)
    assert list(ds.patients) == ["SYN-1", "SYN-2", "SYN-3", "SYN-4"]
    assert ds.source.rows == 7
    assert ds.source.file_name == "synthetic_sheth.csv"
    assert not ds.source.is_reference_copy
    v = ds.patients["SYN-1"].variants[0]
    assert v.source_line == 2 and v.number == 1 and v.sv_class == "inter_chr"
    assert (v.breakpoint_1.chromosome, v.breakpoint_1.position, v.breakpoint_1.strand) == (1, 100, "+")
    assert ds.patients["SYN-1"].chromosomes == (1, 2)


def test_missing_file_explains_where_to_get_data(tmp_path):
    with pytest.raises(FileNotFoundError, match="docs/data.md"):
        load_dataset(tmp_path / "absent.csv")


@pytest.mark.parametrize(
    "content, message",
    [
        ("Individual,Number\nA,1\n", "missing columns"),
        (SYNTHETIC.read_text(encoding="utf-8-sig").replace("1,+,100", "1,x,100"), "strand"),
        (SYNTHETIC.read_text(encoding="utf-8-sig").replace("1,+,100", "1,+,1e2"), "not an integer"),
    ],
)
def test_malformed_files_are_rejected(tmp_path, content, message):
    path = tmp_path / "bad.csv"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(DatasetFormatError, match=message):
        load_dataset(path)


def test_default_path_uses_environment_variable(monkeypatch, tmp_path):
    monkeypatch.setenv("AMG_SHETH_DATA", str(tmp_path / "x.csv"))
    assert default_data_path() == tmp_path / "x.csv"
    monkeypatch.delenv("AMG_SHETH_DATA")
    assert default_data_path() == Path("data/external/nihms.csv")


# -- conversion (synthetic) --------------------------------------------------------


@pytest.fixture(scope="module")
def synthetic():
    return load_dataset(SYNTHETIC).patients


def test_conversion_of_two_linked_chromosomes(synthetic):
    patient = synthetic["SYN-1"]
    assert chromosome_components(patient) == ((1, 2),)
    pc = convert(patient, (1, 2))
    # chromosome 1: DSBs at 100, 300, 700; chromosome 2: DSB at 500.
    assert pc.theta == InitialConfiguration((3, 1))
    assert [(e.chromosome, e.position, e.side) for e in pc.ends[:2]] == [(1, 100, "left"), (1, 100, "right")]
    # '+' uses the right end, '-' the left end.
    pairs = [(e.end_a, e.end_b) for e in pc.observed]
    assert pairs == [(1, 6), (2, 5)]
    assert all(e.status is EdgeStatus.OBSERVED for e in pc.observed)
    assert pc.observed[0].variant.source_line == 2
    assert pc.free_ends == (0, 3, 4, 7)
    assert pc.reconstructed_matching() is None
    assert pc.ends[5].label() == "chr1:700 (right)"


def test_unique_completion_is_marked_as_reconstructed(synthetic):
    pc = convert(synthetic["SYN-2"], [5])
    matching, inferred = pc.reconstructed_matching()
    assert matching.pairs() == ((0, 3), (1, 2))
    assert [(e.end_a, e.end_b, e.status) for e in inferred] == [(0, 3, EdgeStatus.UNIQUE_COMPLETION)]
    assert str(cycle_structure(pc.theta, matching)) == "C2"


def test_end_used_twice_is_an_error(synthetic):
    with pytest.raises(ConversionError, match="used by junctions on source lines 5 and 6"):
        convert(synthetic["SYN-3"], [3])


def test_same_position_on_both_strands_is_one_dsb_with_both_ends_used(synthetic):
    pc = convert(synthetic["SYN-4"], [6, 7])
    assert pc.theta.breaks == (2, 1)  # chr6: 50, 51; chr7: 80
    assert sorted(v for e in pc.observed for v in (e.end_a, e.end_b)) == [0, 1, 2, 5]
    assert any("adjacent" in note for note in pc.notes)


def test_selection_must_be_a_union_of_components(synthetic):
    with pytest.raises(ConversionError, match="union of components"):
        convert(synthetic["SYN-1"], [1])
    with pytest.raises(ConversionError, match="no junction"):
        convert(synthetic["SYN-1"], [9])


# -- regression: the paper's case study (no data file needed) ---------------------


def fixture_state(component):
    spec = P05_1657[component]
    theta = InitialConfiguration(spec["breaks"])
    return theta, PartialMatching(theta.num_ends, spec["observed"])


def test_table_2_chromosomes_8_and_12():
    theta, fixed = fixture_state((8, 12))
    assert theta.num_dsbs == 10 and len(fixed.free_ends) == 10
    assert count_states(theta, StateSpace.ALL, fixed) == 945
    assert by_string(exact_cycle_distribution(theta, StateSpace.PROPER, fixed)) == TABLE_2


def test_chromosome_7_has_three_completions():
    # Paper, section 7.1: two C4 realisations and one 2C2.
    theta, fixed = fixture_state((7,))
    assert by_string(exact_cycle_distribution(theta, StateSpace.PROPER, fixed)) == {"C4": 2, "2C2": 1}


@pytest.mark.parametrize("component", [(4,), (21,)])
def test_rings_on_chromosomes_4_and_21(component):
    # Paper: vertices 2 and 3 form a ring, so 1 and 4 must be joined; one 2-cycle.
    theta, fixed = fixture_state(component)
    assert by_string(exact_cycle_distribution(theta, StateSpace.PROPER, fixed)) == {"C2": 1}


# -- the real dataset, when available ---------------------------------------------


@needs_real_data
def test_real_dataset_is_the_reference_copy():
    ds = load_dataset(REAL)
    assert ds.source.sha256 == REFERENCE_SHA256
    assert ds.source.rows == 5710
    assert len(ds.patients) == 57


@needs_real_data
def test_real_conversion_of_p05_1657_matches_fixture():
    patient = load_dataset(REAL).patients["P05-1657"]
    assert chromosome_components(patient) == tuple(P05_1657)
    for component, spec in P05_1657.items():
        pc = convert(patient, component)
        assert pc.theta.breaks == spec["breaks"]
        assert tuple((e.end_a, e.end_b) for e in pc.observed) == spec["observed"]
