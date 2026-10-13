"""Build small .xlsx workbooks for tests, with the parts Excel writes.

Cells are given as Python values: ``str`` becomes a shared string, ``int`` and
``float`` a number, ``None`` an empty cell. Numbers in ``scientific`` columns
get Excel's built-in ``0.00E+00`` format (as Quality and Score in Baca et al.).
"""

import io
import zipfile
from xml.sax.saxutils import escape

_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_PACKAGE_REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def build_xlsx(
    sheets: dict[str, list[list]], scientific: frozenset[int] = frozenset()
) -> bytes:
    strings: list[str] = []
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as book:
        names = list(sheets)
        for i, name in enumerate(names, 1):
            book.writestr(
                f"xl/worksheets/sheet{i}.xml", _sheet(sheets[name], strings, scientific)
            )
        book.writestr(
            "xl/workbook.xml",
            f'<workbook xmlns="{_MAIN}" xmlns:r="{_REL}"><sheets>'
            + "".join(
                f'<sheet name="{escape(n)}" sheetId="{i}" r:id="rId{i}"/>'
                for i, n in enumerate(names, 1)
            )
            + "</sheets></workbook>",
        )
        book.writestr(
            "xl/_rels/workbook.xml.rels",
            f'<Relationships xmlns="{_PACKAGE_REL}">'
            + "".join(
                f'<Relationship Id="rId{i}" Type="{_REL}/worksheet" '
                f'Target="worksheets/sheet{i}.xml"/>'
                for i in range(1, len(names) + 1)
            )
            + "</Relationships>",
        )
        book.writestr(
            "xl/sharedStrings.xml",
            f'<sst xmlns="{_MAIN}">'
            + "".join(f"<si><t>{escape(s)}</t></si>" for s in strings)
            + "</sst>",
        )
        book.writestr(
            "xl/styles.xml",
            f'<styleSheet xmlns="{_MAIN}"><cellXfs count="2">'
            '<xf numFmtId="0"/><xf numFmtId="11" applyNumberFormat="1"/>'
            "</cellXfs></styleSheet>",
        )
    return buffer.getvalue()


def _sheet(rows: list[list], strings: list[str], scientific: frozenset[int]) -> str:
    out = []
    for r, row in enumerate(rows, 1):
        cells = []
        for c, value in enumerate(row):
            ref = f"{_letters(c)}{r}"
            if value is None:
                continue
            if isinstance(value, str):
                strings.append(value)
                cells.append(f'<c r="{ref}" t="s"><v>{len(strings) - 1}</v></c>')
            else:
                style = ' s="1"' if c in scientific else ""
                cells.append(f'<c r="{ref}"{style}><v>{value!r}</v></c>')
        out.append(f'<row r="{r}">{"".join(cells)}</row>')
    return (
        f'<worksheet xmlns="{_MAIN}"><sheetData>{"".join(out)}</sheetData></worksheet>'
    )


def _letters(index: int) -> str:
    letters = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        letters = chr(ord("A") + rem) + letters
    return letters
