"""Marco 2026-08-01 -- double zebra striping on every matrix grid: alternate DATA rows
(even rows light grey) AND DATA columns (even cols a translucent overlay), so an entry is
easy to locate by its (row, column) index. Structure-safe: striping is pure CSS
nth-child on the existing <td>, so the grid parser still reads entries verbatim.

This pins the CSS is present on BOTH surfaces (report <style> + the shared gui.css) and
-- the load-bearing constraint -- that a rendered grid still parses byte-for-byte."""
import pathlib

import pytest

from quiverlab.trace import render_html as rh
from quiverlab.trace.render_html import matrix_grid

_SRC = pathlib.Path(rh.__file__).read_text(encoding="utf-8")
_ROOT = pathlib.Path(__file__).resolve().parents[2]


@pytest.mark.oracle_selfcert
def test_report_striping_css_present():
    for needle in ("table.ql-matrix tr:nth-child(even) td",
                   "table.ql-matrix td:nth-child(even)",
                   "print-color-adjust:exact"):
        assert needle in _SRC, needle


@pytest.mark.oracle_selfcert
def test_gui_striping_css_present_both_copies():
    a = (_ROOT / "docs/gui/gui.css").read_text(encoding="utf-8")
    b = (_ROOT / "webapp/static/gui/gui.css").read_text(encoding="utf-8")
    assert a == b, "gui.css copies must be byte-identical"
    for needle in ("table.qlgui-matrix tr:nth-child(even) td",
                   "table.qlgui-matrix td:nth-child(even)",
                   'data-md-color-scheme="slate"'):
        assert needle in a, needle


@pytest.mark.oracle_selfcert
def test_striping_does_not_change_grid_structure():
    """The grid parser reads entries from literal <tr>/<td> -- the striping must not add
    attributes to those tags (it uses nth-child only), so entries still parse."""
    import sys
    sys.path.insert(0, str(_ROOT / "tests" / "trace"))
    from _matrix_grid import grids
    html = matrix_grid([["1", "0"], ["2", "3"]], label="d")
    assert grids(html) == [[["1", "0"], ["2", "3"]]]
    assert "<tr>" in html and "<td>1</td>" in html      # literal tags, no attributes
