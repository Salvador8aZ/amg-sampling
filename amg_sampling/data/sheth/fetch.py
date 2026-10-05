"""Download the Sheth dataset on first use, in the style of Pooch.

As in ``python-minecraft-data``, the file is never redistributed: it is fetched
from its source at a pinned (immutable) commit, checked against a known SHA-256
and kept in a local cache, ``data/external/``. A file already there is used as
is, so ``import-data`` copies and ``$AMG_SHETH_DATA`` are never overwritten.
"""

from __future__ import annotations

import os
import tempfile
import urllib.request
from datetime import date
from pathlib import Path

from amg_sampling.data.sheth.loader import REFERENCE_SHA256, file_sha256

SOURCE_COMMIT = "5cc6a8bb23d148096f80c7998207281527a40f0d"
SOURCE_URL = (
    "https://raw.githubusercontent.com/siddharthsheth/aberration_multigraph/"
    f"{SOURCE_COMMIT}/data/nihms.csv"
)

CITATION = """\
Data: structural-variant calls from Baca S.C. et al., "Punctuated evolution of
prostate cancer genomes", Cell 153:666-677 (2013), doi:10.1016/j.cell.2013.03.027.
Tidied as data/nihms.csv by Sheth, Arsuaga & Sazdanovic in
https://github.com/siddharthsheth/aberration_multigraph.
This copy was downloaded from its source and is not redistributed by amg_sampling.
Cite the original publications when using it.
"""


def fetch(
    dest: str | os.PathLike, url: str = SOURCE_URL, sha256: str = REFERENCE_SHA256
) -> Path:
    """Return ``dest``, downloading it from ``url`` first if it does not exist.

    The download goes to a temporary file and is moved into place only if its
    SHA-256 equals ``sha256``. A provenance note is written next to it.
    """
    dest = Path(dest)
    if dest.is_file():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    handle, tmp = tempfile.mkstemp(dir=dest.parent, suffix=".part")
    try:
        with os.fdopen(handle, "wb") as out, urllib.request.urlopen(url) as response:
            for block in iter(lambda: response.read(1 << 16), b""):
                out.write(block)
        actual = file_sha256(Path(tmp))
        if actual != sha256:
            raise ValueError(
                f"Download from {url} has SHA-256 {actual}, expected {sha256}; "
                "file discarded."
            )
        os.replace(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    dest.with_name(dest.name + ".SOURCE.txt").write_text(
        f"Source: {url}\nSHA-256: {sha256}\nDownloaded: {date.today()}\n\n" + CITATION,
        encoding="utf-8",
    )
    return dest
