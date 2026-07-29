"""Read the report's Result dimensions back out of the rendered HTML.

Marco (2026-07-29): the Result section used to be one long equation inside a
horizontally-scrolling box; it is now a degree TABLE, the same shape the GUI uses
for Ext/Tor. These helpers let the renderer tests keep asserting on the VALUES
without hard-coding either presentation.
"""
import re

_TABLE = re.compile(
    r'<table class="ql-dims">\s*<tr>(?P<head>.*?)</tr>\s*<tr>(?P<body>.*?)</tr>\s*'
    r'</table>', re.S)
_CELL = re.compile(r"<t[dh]>(.*?)</t[dh]>", re.S)


def result_tables(html):
    """Every Result table in the page, as ``(row_label, [dim, ...])`` pairs."""
    out = []
    for m in _TABLE.finditer(html):
        cells = _CELL.findall(m.group("body"))
        if not cells:
            continue
        out.append((cells[0].strip(), [int(c) for c in cells[1:]]))
    return out


def result_dims(html, label=None):
    """The dims of the (first matching) Result table; ``[]`` when there is none.

    ``label`` filters on the row label, e.g. ``"dim HH^n"`` / ``"dim Ext^n"``.
    """
    for row_label, dims in result_tables(html):
        if label is None or row_label == label:
            return dims
    return []


def has_result(html, kind, dims):
    """True when the page shows exactly ``dims`` for the ``HH^`` / ``HH_`` (or
    ``Ext^`` / ``Tor_``) row -- the table-form replacement for the old
    ``"HH^{0} = 3" in html`` assertion."""
    return result_dims(html, "dim %sn" % kind) == list(dims)
