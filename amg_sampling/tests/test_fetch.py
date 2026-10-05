"""Pooch-style download of the dataset, tested offline through file:// URLs."""

import pytest
from omegaconf import OmegaConf

from amg_sampling import directories
from amg_sampling.data.sheth import REFERENCE_SHA256
from amg_sampling.data.sheth.fetch import SOURCE_URL, fetch
from amg_sampling.data.sheth.loader import file_sha256

SYNTHETIC = directories.test_data("synthetic_sheth.csv")


def test_fetch_downloads_verifies_and_records_provenance(tmp_path):
    dest = tmp_path / "external" / "nihms.csv"
    got = fetch(dest, SYNTHETIC.as_uri(), file_sha256(SYNTHETIC))
    assert got == dest
    assert dest.read_bytes() == SYNTHETIC.read_bytes()
    note = (tmp_path / "external" / "nihms.csv.SOURCE.txt").read_text()
    assert SYNTHETIC.as_uri() in note and "Baca" in note
    assert list(dest.parent.glob("*.part")) == []


def test_fetch_rejects_a_checksum_mismatch(tmp_path):
    dest = tmp_path / "nihms.csv"
    with pytest.raises(ValueError, match="SHA-256"):
        fetch(dest, SYNTHETIC.as_uri(), "0" * 64)
    assert list(tmp_path.iterdir()) == []


def test_fetch_keeps_an_existing_file(tmp_path):
    dest = tmp_path / "nihms.csv"
    dest.write_text("local copy")
    assert fetch(dest, "file:///does/not/exist", "0" * 64) == dest
    assert dest.read_text() == "local copy"


def test_config_matches_the_pinned_source():
    data = OmegaConf.load(directories.code("conf/config.yaml")).data
    assert data.url == SOURCE_URL
    assert data.sha256 == REFERENCE_SHA256
