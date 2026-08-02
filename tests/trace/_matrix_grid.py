"""Read matrices back out of the rendered report.

Marco (2026-07-29): matrices are displayed as INDEXED GRIDS -- an HTML table with a
header row of column indices, a header column of row indices, and a light rule --
rather than as a typeset ``pmatrix``. These helpers let the renderer tests keep
asserting on the ENTRIES (and on "shown in full") without hard-coding either
presentation.
"""
import html as _html
import re

_TABLE = re.compile(r'<table class="ql-matrix">(?P<body>.*?)</table>', re.S)
_CAYLEY = re.compile(r'<table class="ql-matrix ql-cayley">(?P<body>.*?)</table>', re.S)
_ROW = re.compile(r"<tr>(.*?)</tr>", re.S)
_CELL = re.compile(r"<td>(.*?)</td>", re.S)
_CELLX = re.compile(r"<td[^>]*>(.*?)</td>", re.S)     # tolerant of a degree-boundary class
_HEAD = re.compile(r"<th[^>]*>(.*?)</th>", re.S)
_ANNOT = re.compile(r"<annotation[^>]*>(.*?)</annotation>", re.S)


def grids(html):
    """Every matrix grid on the page, as a list of rows of cell strings (the index
    headers are stripped -- what comes back is the matrix itself)."""
    out = []
    for m in _TABLE.finditer(html):
        rows = []
        for tr in _ROW.findall(m.group("body")):
            cells = _CELL.findall(tr)
            if cells:                       # the header row has no <td> at all
                rows.append([c.strip() for c in cells])
        out.append(rows)
    return out


def has_grid(html, matrix):
    """True when some grid on the page is exactly ``matrix`` (entries compared as
    strings, so 1 and "1" match)."""
    want = [[str(x) for x in row] for row in matrix]
    return want in grids(html)


def _cell_source(cell):
    """The TeX source of a Cayley cell -- the ``<annotation>`` body of the typeset math
    (HTML-unescaped), or the raw cell text for a plain ``0``."""
    a = _ANNOT.search(cell)
    return _html.unescape(a.group(1).strip()) if a else cell.strip()


def cayley_cells(html):
    """Every Cayley product grid on the page, as rows of cell TeX strings (the index
    header row/column are ``<th>`` and drop out). A nonzero cell yields its TeX
    (e.g. ``\\alpha^{2}_{1}``); a zero cell yields ``"0"``."""
    out = []
    for m in _CAYLEY.finditer(html):
        rows = []
        for tr in _ROW.findall(m.group("body")):
            cells = _CELLX.findall(tr)
            if cells:
                rows.append([_cell_source(c) for c in cells])
        out.append(rows)
    return out


def cayley_headers(html):
    """``(col_labels, row_labels)`` TeX of the FIRST Cayley grid -- the class-name
    headers (corner cell dropped), so a test can pin they carry α/z class names."""
    m = _CAYLEY.search(html)
    if not m:
        return [], []
    trs = _ROW.findall(m.group("body"))
    cols = [_cell_source(h) for h in _HEAD.findall(trs[0])][1:]     # drop the corner
    rows = [_cell_source(_HEAD.findall(tr)[0]) for tr in trs[1:] if _HEAD.findall(tr)]
    return cols, rows


def grid_indices(html):
    """``(column_headers, row_headers)`` of the FIRST grid -- the extra index row and
    column themselves, so a test can pin that they are present and 1-based."""
    m = _TABLE.search(html)
    if not m:
        return [], []
    trs = _ROW.findall(m.group("body"))
    cols = [h.strip() for h in _HEAD.findall(trs[0])][1:]     # drop the corner cell
    rows = [_HEAD.findall(tr)[0].strip() for tr in trs[1:] if _HEAD.findall(tr)]
    return cols, rows
