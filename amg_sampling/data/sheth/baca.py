"""Convert the Baca et al. (2013) supplementary workbook into ``nihms.csv``.

The structural-variant calls are Table S3C of ``mmc3.xlsx``, the supplementary
material of Baca et al., Cell 153:666-677 (2013). Sheth, Arsuaga & Sazdanovic
distribute the same table as ``data/nihms.csv`` (an Excel "Save as CSV" of it).
This module rebuilds that CSV from the workbook itself, using only the
standard library (an ``.xlsx`` file is a zip archive of XML parts):

- cells are written as Excel displays them: integers without a decimal point,
  cells formatted ``0.00E+00`` (Quality, Score) as e.g. ``9.81E-01``;
- blank rows at the end of the sheet are dropped;
- rows are sorted by patient (``Individual``) and junction ``Number``. The
  workbook lists them in another order of its own, so a CSV line number does
  not point at a workbook row: patient and ``Number`` identify a junction in
  the workbook and in both CSV copies;
- the file has a UTF-8 byte-order mark and CRLF line endings, like the
  Excel export.

The result has the same columns and values as the Sheth copy. The two differ
only in the row order of patient P05-1657, which the Sheth copy lists by
breakpoint 1 (chromosome, then position) instead of junction number.
"""

from __future__ import annotations

import csv
import io
import re
import zipfile
from xml.etree import ElementTree

SHEET = "Table S3C"

_MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_PACKAGE_REL = "{http://schemas.openxmlformats.org/package/2006/relationships}"
# Built-in number format 11 is "0.00E+00".
_SCIENTIFIC = 11


class WorkbookError(ValueError):
    """The workbook is not readable or lacks the expected sheet."""


def read_sheet(xlsx: bytes, name: str = SHEET) -> list[list[str]]:
    """The cells of sheet ``name`` as displayed text, one list per row.

    Rows are padded to the header width; rows with no value are dropped.
    """
    with _open(xlsx) as book:
        paths = _sheet_paths(book)
        if name not in paths:
            raise WorkbookError(f"No sheet {name!r}; sheets: {list(paths)}.")
        strings = _shared_strings(book)
        formats = _cell_formats(book)
        root = ElementTree.fromstring(book.read(paths[name]))
    rows = []
    for row in root.iter(f"{_MAIN}row"):
        cells = {}
        for cell in row.iter(f"{_MAIN}c"):
            text = _cell_text(cell, strings, formats)
            if text != "":
                cells[_column(cell.get("r"))] = text
        if cells:
            rows.append(cells)
    if not rows:
        raise WorkbookError(f"Sheet {name!r} is empty.")
    width = max(rows[0]) + 1
    return [[cells.get(j, "") for j in range(width)] for cells in rows]


def to_csv(xlsx: bytes, name: str = SHEET) -> bytes:
    """``nihms.csv`` built from sheet ``name`` of the Baca et al. workbook."""
    header, *rows = read_sheet(xlsx, name)
    try:
        individual, number = header.index("Individual"), header.index("Number")
        rows.sort(key=lambda r: (r[individual], int(r[number])))
    except ValueError as error:
        raise WorkbookError(f"Sheet {name!r}: {error}.") from None
    buffer = io.StringIO()
    csv.writer(buffer, lineterminator="\r\n").writerows([header, *rows])
    return ("﻿" + buffer.getvalue().removesuffix("\r\n")).encode("utf-8")


def _open(xlsx: bytes) -> zipfile.ZipFile:
    try:
        return zipfile.ZipFile(io.BytesIO(xlsx))
    except zipfile.BadZipFile:
        raise WorkbookError("Not an .xlsx workbook (not a zip archive).") from None


def _sheet_paths(book: zipfile.ZipFile) -> dict[str, str]:
    rels = ElementTree.fromstring(book.read("xl/_rels/workbook.xml.rels"))
    targets = {
        r.get("Id"): r.get("Target") for r in rels.iter(f"{_PACKAGE_REL}Relationship")
    }
    workbook = ElementTree.fromstring(book.read("xl/workbook.xml"))
    paths = {}
    for sheet in workbook.iter(f"{_MAIN}sheet"):
        target = targets[sheet.get(f"{_REL}id")]
        paths[sheet.get("name")] = (
            target.lstrip("/") if target.startswith("/") else f"xl/{target}"
        )
    return paths


def _shared_strings(book: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in book.namelist():
        return []
    root = ElementTree.fromstring(book.read("xl/sharedStrings.xml"))
    return ["".join(t.text or "" for t in si.iter(f"{_MAIN}t")) for si in root]


def _cell_formats(book: zipfile.ZipFile) -> list[int]:
    if "xl/styles.xml" not in book.namelist():
        return []
    xfs = ElementTree.fromstring(book.read("xl/styles.xml")).find(f"{_MAIN}cellXfs")
    return [] if xfs is None else [int(xf.get("numFmtId", 0)) for xf in xfs]


def _cell_text(cell, strings: list[str], formats: list[int]) -> str:
    kind = cell.get("t")
    if kind == "inlineStr":
        return "".join(t.text or "" for t in cell.iter(f"{_MAIN}t"))
    value = cell.find(f"{_MAIN}v")
    if value is None or value.text is None:
        return ""
    if kind == "s":
        return strings[int(value.text)]
    if kind in ("str", "b", "e"):
        return value.text
    number = float(value.text)
    style = int(cell.get("s", 0))
    if style < len(formats) and formats[style] == _SCIENTIFIC:
        return f"{number:.2E}"
    if number.is_integer():
        return str(int(number))
    # Excel's "General" format shows at most 15 significant digits.
    return f"{number:.15g}"


def _column(reference: str) -> int:
    index = 0
    for letter in re.match(r"[A-Z]+", reference).group():
        index = index * 26 + ord(letter) - ord("A") + 1
    return index - 1
