"""Pooch-style download of the dataset, tested offline through file:// URLs."""

import hashlib

import pytest
from omegaconf import OmegaConf

from amg_sampling import directories
from amg_sampling.data.sheth import REFERENCE_SHA256, SHETH_SHA256, load_dataset
from amg_sampling.data.sheth.baca import WorkbookError, read_sheet, to_csv
from amg_sampling.data.sheth.fetch import (
    BACA_DOI,
    BACA_URL,
    SHETH_URL,
    DownloadError,
    fetch,
    supplement_candidates,
)
from amg_sampling.data.sheth.loader import file_sha256
from amg_sampling.tests.workbooks import build_xlsx

SYNTHETIC = directories.test_data("synthetic_sheth.csv")

HEADER = [
    "Individual", "Number",
    "Breakpoint 1 chromosome", "Breakpoint 1 strand", "Breakpoint 1 position",
    "Breakpoint 2 chromosome", "Breakpoint 2 strand", "Breakpoint 2 position",
    "Class", "Span", "Quality", "Foreign sequence",
]  # fmt: skip
QUALITY = HEADER.index("Quality")
# Out of junction order, as in the Sheth copy of P05-1657, and followed by the
# blank rows found at the end of Table S3C.
ROWS = [
    ["SYN-2", 1, 5, "+", 1000, 5, "-", 2000, "long_range", 1000, 1, None],
    ["SYN-1", 10, 1, "-", 300, 1, "+", 700, "inversion", 400, 0.98069139999999999, "A"],
    ["SYN-1", 2, 1, "+", 100, 2, "-", 500, "inter_chr", "NaN", 0.5894924, None],
    [None] * len(HEADER),
]
EXPECTED_CSV = (
    "﻿"
    + ",".join(HEADER)
    + "\r\nSYN-1,2,1,+,100,2,-,500,inter_chr,NaN,5.89E-01,"
    + "\r\nSYN-1,10,1,-,300,1,+,700,inversion,400,9.81E-01,A"
    + "\r\nSYN-2,1,5,+,1000,5,-,2000,long_range,1000,1.00E+00,"
).encode("utf-8")


def workbook(rows=ROWS, sheet="Table S3C") -> bytes:
    return build_xlsx(
        {"Table S3A": [["other"]], sheet: [HEADER, *rows]},
        scientific=frozenset({QUALITY}),
    )


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@pytest.fixture
def supplement(tmp_path):
    path = tmp_path / "source" / "mmc3.xlsx"
    path.parent.mkdir()
    path.write_bytes(workbook())
    return path


# -- converting the Baca et al. workbook -------------------------------------------


def test_to_csv_formats_sorts_and_drops_blank_rows():
    assert to_csv(workbook()) == EXPECTED_CSV


def test_converted_csv_loads(tmp_path):
    path = tmp_path / "nihms.csv"
    path.write_bytes(to_csv(workbook()))
    ds = load_dataset(path)
    assert ds.source.rows == 3
    assert [v.number for v in ds.patients["SYN-1"].variants] == [2, 10]
    assert not ds.source.is_reference_copy


def test_read_sheet_pads_rows_to_the_header():
    rows = read_sheet(workbook())
    assert {len(r) for r in rows} == {len(HEADER)}
    assert rows[-1][-1] == ""


def test_missing_sheet_or_non_workbook_is_rejected():
    with pytest.raises(WorkbookError, match="Table S3C"):
        to_csv(workbook(sheet="Table 1"))
    with pytest.raises(WorkbookError, match="zip"):
        to_csv(b"<html>moved</html>")


# -- downloading -------------------------------------------------------------------


def test_fetch_converts_the_baca_supplement(tmp_path, supplement):
    dest = tmp_path / "external" / "nihms.csv"
    messages = []
    got = fetch(
        dest,
        baca_url=supplement.as_uri(),
        doi=None,
        sheth_url=None,
        sha256=sha(EXPECTED_CSV),
        notify=messages.append,
    )
    assert got == dest
    assert dest.read_bytes() == EXPECTED_CSV
    note = (tmp_path / "external" / "nihms.csv.SOURCE.txt").read_text()
    assert supplement.as_uri() in note
    assert "Baca" in note and "belong to Baca et al." in note
    assert messages and "Baca" in messages[0]
    assert list(dest.parent.glob("*.part")) == []


def test_fetch_falls_back_to_the_sheth_csv(tmp_path):
    dest = tmp_path / "nihms.csv"
    fetch(
        dest,
        baca_url=(tmp_path / "moved.xlsx").as_uri(),
        doi=None,
        sheth_url=SYNTHETIC.as_uri(),
        sheth_sha256=file_sha256(SYNTHETIC),
    )
    assert dest.read_bytes() == SYNTHETIC.read_bytes()
    note = (tmp_path / "nihms.csv.SOURCE.txt").read_text()
    assert SYNTHETIC.as_uri() in note
    assert "Sheth" in note and "belong to Baca et al." in note


def test_fetch_rejects_a_converted_checksum_mismatch(tmp_path, supplement):
    dest = tmp_path / "nihms.csv"
    with pytest.raises(DownloadError, match="SHA-256"):
        fetch(
            dest,
            baca_url=supplement.as_uri(),
            doi=None,
            sheth_url=None,
            sha256="0" * 64,
        )
    assert not dest.exists()


def test_fetch_reports_every_failed_source(tmp_path):
    dest = tmp_path / "nihms.csv"
    with pytest.raises(DownloadError) as error:
        fetch(
            dest,
            baca_url=(tmp_path / "a.xlsx").as_uri(),
            doi=None,
            sheth_url=SYNTHETIC.as_uri(),
            sheth_sha256="0" * 64,
        )
    assert "Baca et al. supplement" in str(error.value)
    assert "Sheth et al. CSV" in str(error.value)
    assert list(tmp_path.iterdir()) == []


def test_fetch_keeps_an_existing_file(tmp_path):
    dest = tmp_path / "nihms.csv"
    dest.write_text("local copy")
    assert fetch(dest, baca_url="file:///does/not/exist", sha256="0" * 64) == dest
    assert dest.read_text() == "local copy"


# -- locating the workbook from the DOI ---------------------------------------------


def test_candidates_from_an_elsevier_landing_address():
    landing = "https://linkinghub.elsevier.com/retrieve/pii/S0092867413003437"
    assert supplement_candidates(landing, "") == [BACA_URL]


def test_candidates_from_page_links_put_table_s3_first():
    page = (
        '<a href="/files/mmc1.xlsx">S1</a> <a href="https://cdn.example/x/mmc3.xlsx">'
        "S3</a> cited: pii/S0092867413009999"
    )
    assert supplement_candidates("https://journal.example/article/1", page) == [
        "https://cdn.example/x/mmc3.xlsx",
        "https://journal.example/files/mmc1.xlsx",
    ]


def test_config_matches_the_pinned_sources():
    data = OmegaConf.load(directories.code("conf/config.yaml")).data
    assert data.baca_url == BACA_URL
    assert data.doi == BACA_DOI
    assert data.sha256 == REFERENCE_SHA256
    assert data.sheth_url == SHETH_URL
    assert data.sheth_sha256 == SHETH_SHA256
