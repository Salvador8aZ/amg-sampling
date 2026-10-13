"""Download the patient dataset on first use, in the style of Pooch.

As in ``python-minecraft-data``, the file is never redistributed: it is fetched
from its source, checked against a known SHA-256 and kept in a local cache,
``data/external/``. A file already there is used as is, so ``import-data``
copies and ``$AMG_SHETH_DATA`` are never overwritten.

Sources are tried in order, and the first one that verifies is kept:

1. the supplementary workbook of Baca et al. (2013) at its publisher URL,
   converted to ``nihms.csv`` by :mod:`baca`;
2. the same workbook found again from the paper's DOI, in case the publisher
   has moved it;
3. the CSV distributed by Sheth et al. at a pinned commit.
"""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import ssl
import subprocess
import tempfile
import urllib.error
import urllib.request
from collections.abc import Callable
from datetime import date
from html import unescape
from pathlib import Path
from urllib.parse import urljoin

from amg_sampling.data.sheth import baca
from amg_sampling.data.sheth.loader import REFERENCE_SHA256, SHETH_SHA256

BACA_DOI = "10.1016/j.cell.2013.03.021"
BACA_PII = "S0092867413003437"
BACA_URL = f"https://ars.els-cdn.com/content/image/1-s2.0-{BACA_PII}-mmc3.xlsx"
DOI_RESOLVER = "https://doi.org/"

SHETH_COMMIT = "5cc6a8bb23d148096f80c7998207281527a40f0d"
SHETH_URL = (
    "https://raw.githubusercontent.com/siddharthsheth/aberration_multigraph/"
    f"{SHETH_COMMIT}/data/nihms.csv"
)

# Certificate store of macOS and most Linux distributions.
_SYSTEM_CA_FILE = "/etc/ssl/cert.pem"
# Some publisher servers refuse urllib's default user agent.
_HEADERS = {"User-Agent": "Mozilla/5.0 (amg_sampling dataset download)"}

BACA_CITATION = """\
Data: structural-variant calls (Table S3C) from the supplementary material of
Baca S.C. et al., "Punctuated evolution of prostate cancer genomes",
Cell 153:666-677 (2013), doi:10.1016/j.cell.2013.03.021.
The data belong to Baca et al. They were downloaded from the publisher and
converted to CSV by amg_sampling; the file is not redistributed by amg_sampling.
Cite the original publication when using it, and Sheth, Arsuaga & Sazdanovic
(2026), doi:10.1088/1751-8121/ae4d8f, for the analysis methodology.
"""

SHETH_CITATION = """\
Data: structural-variant calls (Table S3C) from the supplementary material of
Baca S.C. et al., "Punctuated evolution of prostate cancer genomes",
Cell 153:666-677 (2013), doi:10.1016/j.cell.2013.03.021.
The data belong to Baca et al. The Baca supplement could not be downloaded, so
this is the copy tidied as data/nihms.csv by Sheth, Arsuaga & Sazdanovic in
https://github.com/siddharthsheth/aberration_multigraph. It has the same rows
as the supplement; patient P05-1657 is listed by breakpoint 1 (chromosome,
then position) instead of junction number. The file is not redistributed by
amg_sampling.
Cite the original publications when using it.
"""


class DownloadError(OSError):
    """No source produced a file with the expected SHA-256."""


def fetch(
    dest: str | os.PathLike,
    *,
    baca_url: str | None = BACA_URL,
    doi: str | None = BACA_DOI,
    sheth_url: str | None = SHETH_URL,
    sha256: str = REFERENCE_SHA256,
    sheth_sha256: str = SHETH_SHA256,
    notify: Callable[[str], None] | None = None,
) -> Path:
    """Return ``dest``, downloading it first if it does not exist.

    A source is skipped by passing ``None``. Converted Baca data must have
    SHA-256 ``sha256``; the Sheth CSV must have ``sheth_sha256``. Nothing is
    written unless a source verifies. A provenance note is written next to
    the file, and passed to ``notify`` (e.g. ``print``) after a download.
    """
    dest = Path(dest)
    if dest.is_file():
        return dest
    failures = []
    tried = set()

    def attempt(label: str, url: str, get: Callable[[str], bytes], expected: str):
        tried.add(url)
        try:
            data = get(url)
        except (OSError, ValueError) as error:
            failures.append(f"{label} ({url}): {error}")
            return None
        actual = _sha256(data)
        if actual != expected:
            failures.append(f"{label} ({url}): SHA-256 {actual}, expected {expected}")
            return None
        return data

    result = None
    if baca_url:
        data = attempt("Baca et al. supplement", baca_url, _baca_csv, sha256)
        if data is not None:
            result = (data, baca_url, BACA_CITATION)
    if result is None and doi:
        try:
            candidates = [url for url in locate_by_doi(doi) if url not in tried]
        except (OSError, ValueError) as error:
            failures.append(f"DOI {doi}: {error}")
            candidates = []
        else:
            if not candidates:
                failures.append(f"DOI {doi}: no other location of the workbook found")
        for url in candidates:
            label = f"Baca et al. supplement via DOI {doi}"
            data = attempt(label, url, _baca_csv, sha256)
            if data is not None:
                result = (data, url, BACA_CITATION)
                break
    if result is None and sheth_url:
        data = attempt("Sheth et al. CSV", sheth_url, _download, sheth_sha256)
        if data is not None:
            result = (data, sheth_url, SHETH_CITATION)
    if result is None:
        hint = ""
        if any("CERTIFICATE_VERIFY_FAILED" in f for f in failures):
            hint = (
                "\nPython has no certificate store: run 'Install Certificates.command' "
                "in /Applications/Python 3.x/, or set SSL_CERT_FILE (see docs/data.md)."
            )
        raise DownloadError(
            "Could not download the dataset; file not written. Tried:\n  "
            + "\n  ".join(failures or ["no sources"])
            + hint
        )

    data, url, citation = result
    _write_atomically(dest, data)
    note = (
        f"Source: {url}\nSHA-256: {_sha256(data)}\nDownloaded: {date.today()}\n\n"
        + citation
    )
    dest.with_name(dest.name + ".SOURCE.txt").write_text(note, encoding="utf-8")
    if notify is not None:
        notify(f"Downloaded {dest}\n{note}")
    return dest


def locate_by_doi(doi: str, resolver: str = DOI_RESOLVER) -> list[str]:
    """Candidate URLs of the supplementary workbook, found from the DOI.

    The DOI is resolved to the publisher's landing page. If its address names
    an Elsevier article (PII), the workbook's address on Elsevier's content
    server is derived from it; ``.xlsx`` links on the page follow. Every
    candidate is still checked against the expected SHA-256.
    """
    request = urllib.request.Request(resolver + doi, headers=_HEADERS)
    with urllib.request.urlopen(request, context=_ssl_context()) as response:
        landing = response.url
        page = response.read(1 << 22).decode("utf-8", errors="replace")
    return supplement_candidates(landing, page)


def supplement_candidates(landing_url: str, page: str) -> list[str]:
    """Workbook URLs derived from a landing page address and its HTML.

    The PII is read from the address only: the page also cites other articles.
    """
    candidates = []
    for pii in re.findall(r"pii/(S\d{16}|S\d{15}X)", landing_url):
        candidates.append(
            f"https://ars.els-cdn.com/content/image/1-s2.0-{pii}-mmc3.xlsx"
        )
    for link in re.findall(
        r"""(?:href|src)\s*=\s*["']([^"']+\.xlsx)["']""", page, re.I
    ):
        candidates.append(urljoin(landing_url, unescape(link)))
    candidates.sort(key=lambda url: "mmc3" not in url)  # Table S3 first
    return list(dict.fromkeys(candidates))


def _download(url: str) -> bytes:
    """The content at ``url``.

    Elsevier's content server sits behind Cloudflare, which answers Python's
    HTTP client with 403 Forbidden but serves the same public file to curl.
    On a 403 the download is therefore retried with curl, if installed.
    """
    request = urllib.request.Request(url, headers=_HEADERS)
    try:
        with urllib.request.urlopen(request, context=_ssl_context()) as response:
            return response.read()
    except urllib.error.HTTPError as error:
        curl = shutil.which("curl")
        if error.code != 403 or curl is None:
            raise
    result = subprocess.run(
        [curl, "--fail", "--silent", "--show-error", "--location",
         "--max-time", "300", "--user-agent", _HEADERS["User-Agent"], url],
        capture_output=True,
        check=False,
    )  # fmt: skip
    if result.returncode != 0:
        message = result.stderr.decode(errors="replace").strip()
        raise OSError(f"HTTP 403 Forbidden; retry with curl failed: {message}")
    return result.stdout


def _ssl_context() -> ssl.SSLContext:
    """Default verification; certificates are always checked.

    The python.org installer for macOS ships without a certificate store until
    'Install Certificates.command' is run. Then, unless SSL_CERT_FILE or
    SSL_CERT_DIR is set, the operating system's store is used instead.
    """
    context = ssl.create_default_context()
    paths = ssl.get_default_verify_paths()
    if (
        not context.get_ca_certs()
        and paths.cafile is None
        and paths.capath is None
        and os.path.isfile(_SYSTEM_CA_FILE)
    ):
        context.load_verify_locations(_SYSTEM_CA_FILE)
    return context


def _baca_csv(url: str) -> bytes:
    return baca.to_csv(_download(url))


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_atomically(dest: Path, data: bytes) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    handle, tmp = tempfile.mkstemp(dir=dest.parent, suffix=".part")
    try:
        with os.fdopen(handle, "wb") as out:
            out.write(data)
        os.replace(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
