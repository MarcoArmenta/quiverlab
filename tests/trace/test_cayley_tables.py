"""Marco 2026-08-01: HH product bidegrees render as CAYLEY MULTIPLICATION TABLES.

Each nonzero product family bidegree is a grid -- rows the left-factor classes
(``alpha^p_i``), columns the right-factor classes (``alpha^q_j``, ``z^n_j`` for cap),
each cell the product written DIRECTLY in the target basis (``0`` or a signed
combination). Zeros are SHOWN inside a table that is nonzero somewhere; an all-zero
bidegree keeps the one-line vanish statement. GF(p) coefficients display as balanced
representatives (``c > p/2`` -> ``c-p``); the JSON keeps the raw residues. Honest
structural notes ("all squares are 0", graded (anti)symmetry) are DERIVED from the
constants, never asserted blindly.

Render-only: the block data (structure constants) is unchanged.
"""
import json

import pytest

import quiverlab as ql
from quiverlab.trace.products import (
    balanced_coeff, balanced_rep_note, cayley_table, cell_tex,
    prime_from_basis, structural_notes, structural_note_line)
from quiverlab.trace.render_html import cayley_grid_html, render_html
from quiverlab.trace.render_json import render_json
from quiverlab.trace.results_html import results_section
from quiverlab.trace.products import products_chapter
from tests.trace._matrix_grid import cayley_cells, cayley_headers

pytestmark = [pytest.mark.oracle_selfcert]


# --------------------------------------------------------------------------- #
# Balanced representatives (display only; the JSON keeps raw residues).
# --------------------------------------------------------------------------- #
def test_balanced_coeff_gf7():
    assert balanced_coeff("6", 7) == "-1"          # p-1 -> -1
    assert balanced_coeff("5", 7) == "-2"
    assert balanced_coeff("4", 7) == "-3"
    assert balanced_coeff("3", 7) == "3"           # 3 <= 7/2 -> unchanged
    assert balanced_coeff("0", 7) == "0"
    assert balanced_coeff("1", 7) == "1"


def test_balanced_coeff_gf2_and_verbatim():
    assert balanced_coeff("1", 2) == "1"           # 1 is not > 2/2
    assert balanced_coeff("6", None) == "6"        # no prime -> verbatim
    assert balanced_coeff("1/2", 7) == "1/2"       # non-residue (fraction) -> verbatim
    assert balanced_coeff("x+1", 4) == "x+1"       # extension-field repr -> verbatim


def test_prime_from_basis():
    assert prime_from_basis("bar/GF(7)") == 7
    assert prime_from_basis("cs/GF(2)") == 2
    assert prime_from_basis("cs/QQ") is None
    assert prime_from_basis(None) is None


def test_balanced_rep_legend_wording():
    assert "balanced representatives mod 7" in balanced_rep_note(7)
    assert balanced_rep_note(None) == ""


# --------------------------------------------------------------------------- #
# Cell text: the product in the target basis, balanced + signed.
# --------------------------------------------------------------------------- #
def test_cell_tex_combination_and_signs():
    # 1 * a^2_1 + 5 * a^2_3 over GF(7): 5 balances to -2
    tex = cell_tex("cup", 2, ["1", "0", "5"], 7)
    assert tex == r"\alpha^{2}_{1} - 2\,\alpha^{2}_{3}"
    assert cell_tex("cup", 2, ["0", "0"], 7) == "0"          # all-zero -> "0"
    assert cell_tex("cap", 1, ["1"], 2) == r"z^{1}_{1}"      # cap outputs z


# --------------------------------------------------------------------------- #
# Structural notes: derived, true AND false.
# --------------------------------------------------------------------------- #
def test_structural_notes_true_antisymmetric_zero_diagonal():
    # cup p=q=1 over GF(3): zero diagonal, antisymmetric off-diagonal (sign -1).
    # k=0 block [[0,1],[2,0]] (2 == -1 mod 3), single output generator.
    constants = ((("0", "1"), ("2", "0")),)
    notes = structural_notes("cup", (1, 1), (2, 2, 1), constants, 3)
    assert "all squares are 0" in notes
    assert "the table is graded-antisymmetric" in notes


def test_structural_notes_false_when_sabotaged():
    # break antisymmetry: [[0,1],[1,0]] is symmetric, not -1-mirrored over GF(3).
    constants = ((("0", "1"), ("1", "0")),)
    notes = structural_notes("cup", (1, 1), (2, 2, 1), constants, 3)
    assert "all squares are 0" in notes                       # diagonal still zero
    assert not any("antisymmetric" in n for n in notes)       # but NOT antisymmetric
    # break the diagonal too: no squares note.
    constants2 = ((("1", "1"), ("1", "0")),)
    notes2 = structural_notes("cup", (1, 1), (2, 2, 1), constants2, 3)
    assert "all squares are 0" not in notes2


def test_structural_notes_only_for_square_cohomology_bidegrees():
    assert structural_notes("cap", (1, 1), (1, 1, 1), ((("1",),),), 7) == []
    assert structural_notes("cup", (0, 1), (1, 1, 1), ((("1",),),), 7) == []   # p != q


def test_structural_note_line_capitalizes_and_joins():
    constants = ((("0", "1"), ("2", "0")),)
    line = structural_note_line("cup", (1, 1), (2, 2, 1), constants, 3)
    assert line == "All squares are 0; the table is graded-antisymmetric."


# --------------------------------------------------------------------------- #
# cayley_table struct: labels + cells; cap uses z columns.
# --------------------------------------------------------------------------- #
def test_cayley_table_cells_and_labels():
    A = ql.truncated_polynomial(2, field=ql.GF(2))
    t = A.cup_products(2).tables[(1, 1)]
    tbl = cayley_table("cup", t.degrees, t.out_degree, list(t.dims), t.constants, 2)
    assert tbl["corner"] == r"\cup"
    assert tbl["row_labels"] == [r"\alpha^{1}_{1}", r"\alpha^{1}_{2}"]
    assert tbl["col_labels"] == [r"\alpha^{1}_{1}", r"\alpha^{1}_{2}"]
    assert tbl["cells"] == [[r"\alpha^{2}_{1}", r"\alpha^{2}_{2}"],
                            [r"\alpha^{2}_{2}", "0"]]


def test_cap_table_rows_alpha_cols_z():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    obj = A.cap_products(2)
    t = next(t for t in obj.tables.values() if t.dims[0] and t.dims[1])
    tbl = cayley_table("cap", t.degrees, t.out_degree, list(t.dims), t.constants,
                       prime_from_basis(obj.basis))
    assert tbl["corner"] == r"\cap"
    assert all(lab.startswith(r"\alpha") for lab in tbl["row_labels"])
    assert all(lab.startswith("z") for lab in tbl["col_labels"])


# --------------------------------------------------------------------------- #
# Display cap: a >=50 x/y table points at the JSON, no grid.
# --------------------------------------------------------------------------- #
def test_cayley_grid_display_cap():
    struct = {"corner": r"\cup", "row_labels": [], "col_labels": [], "cells": [],
              "dl": 50, "dr": 2, "note": ""}
    html = cayley_grid_html(struct)
    assert "ql-cayley" not in html
    assert "display cap" in html and "JSON record" in html
    assert "50×2 product table" in html


# --------------------------------------------------------------------------- #
# The two report surfaces + the JSON.
# --------------------------------------------------------------------------- #
def test_results_surface_renders_grid_and_balanced_legend():
    A = ql.truncated_polynomial(2, field=ql.GF(2))
    block = A.cup_products(2).blocks()
    html = "".join(results_section({"cup": block}))
    assert "ql-cayley" in html
    assert "balanced representatives mod 2" in html
    tables = cayley_cells(html)
    assert len(tables) == 1                          # one big degree-major table
    big = tables[0]
    assert big[2][2] == r"\alpha^{2}_{1}"            # a known nonzero product
    assert big[3][3] == "0"                          # a SHOWN computed zero
    assert big[2][4] == "—"                          # a beyond-window em dash


def test_json_keeps_raw_residues_not_balanced():
    # The bracket over GF(7) has a raw residue 5 (displayed -2). The machine record
    # (trace.json) keeps the RAW residue; the balanced -2 is display-only.
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    events = products_chapter(A, "bracket", A.gerstenhaber_brackets(2))
    obj = json.loads(render_json(events, title="bracket"))
    dumped = json.dumps(obj)
    assert '"5"' in dumped                          # the raw residue survives
    # the balanced display string never leaks into the machine record
    assert "-2\\," not in dumped and "- 2 \\," not in dumped


def test_em_dash_marks_beyond_window_not_zero():
    # The honesty mark: a cell whose target degree exceeds the computed top is an em
    # dash (NOT computed), a computed vanishing product is 0 -- both appear and are
    # distinct, and the legend explains the em dash.
    from tests.trace._matrix_grid import cayley_cells as _cc
    A = ql.truncated_polynomial(2, field=ql.GF(2))
    html = "".join(results_section({"cup": A.cup_products(2).blocks()}))
    big = _cc(html)[0]
    flat = [c for row in big for c in row]
    assert "—" in flat and "0" in flat               # both present, distinct
    assert "beyond the computed window" in html
    assert "shown as 0" in html                       # legend states computed zeros stay 0


def test_all_zero_family_keeps_one_line_vanish(monkeypatch):
    # A whole family that vanishes keeps the one-line statement, never a grid of 0s.
    from quiverlab.trace import products as P
    fake_tables = [{"degrees": (p, q), "out_degree": p + q, "dims": [1, 1, 1],
                    "constants": [[["0"]]]}
                   for p in range(2) for q in range(2) if p + q <= 1]
    from quiverlab.trace.render_html import family_cayley_html
    html = "".join(family_cayley_html("cup", fake_tables, 2))
    assert "All cup products in the served bidegrees vanish." in html
    assert "ql-cayley" not in html


def test_over_axis_cap_falls_back_to_per_bidegree(monkeypatch):
    # When a family's per-axis class count exceeds the cap, the big table is dropped
    # for per-bidegree grids with a stated note + JSON pointer.
    from quiverlab.trace.render_html import family_cayley_html
    # two bidegrees, each 30 classes per axis -> 60 > 50 total on each axis
    def mk(p, q, dl, dr):
        constants = [[["1" if (i == 0 and j == 0) else "0" for j in range(dr)]
                      for i in range(dl)]]
        return {"degrees": (p, q), "out_degree": p + q, "dims": [dl, dr, 1],
                "constants": constants}
    tables = [mk(0, 0, 30, 30), mk(0, 1, 30, 30), mk(1, 0, 30, 30)]
    html = "".join(family_cayley_html("cup", tables, 2))
    assert "exceeding the 50-class display cap" in html
    assert "per-bidegree tables follow" in html
    assert "JSON record" in html
    assert "ql-cayley" in html                        # the per-bidegree grids still render


def test_degree_boundary_separators_present():
    A = ql.truncated_polynomial(2, field=ql.GF(2))
    html = "".join(results_section({"cup": A.cup_products(2).blocks()}))
    assert "ql-degrow" in html and "ql-degcol" in html    # heavier rule at boundaries


def test_combined_note_true_and_sabotaged():
    from quiverlab.trace.products import combined_cayley
    # a graded-antisymmetric, zero-diagonal cup family over GF(3): degrees {0,1}, with
    # HH^0 dim 1 (unit) and HH^1 dim 2; only the (1,1)->HH^2 block matters for squares.
    tables = [
        {"degrees": (0, 0), "out_degree": 0, "dims": [1, 1, 1],
         "constants": [[["1"]]]},                                  # unit * unit = unit
        {"degrees": (0, 1), "out_degree": 1, "dims": [1, 2, 2],
         "constants": [[["1", "0"]], [["0", "1"]]]},               # unit acts as identity
        {"degrees": (1, 0), "out_degree": 1, "dims": [2, 1, 2],
         "constants": [[["1"], ["0"]], [["0"], ["1"]]]},
        {"degrees": (1, 1), "out_degree": 2, "dims": [2, 2, 1],
         "constants": [[["0", "1"], ["2", "0"]]]},                 # antisym, zero diag
    ]
    good = combined_cayley("cup", tables, 3)
    assert "graded-commutative" in good["note"]
    # sabotage the (1,1) block -> no graded law
    sab = [dict(t) for t in tables]
    sab[3] = {"degrees": (1, 1), "out_degree": 2, "dims": [2, 2, 1],
              "constants": [[["0", "1"], ["1", "0"]]]}             # now symmetric, not -1
    bad = combined_cayley("cup", sab, 3)
    assert "graded-commutative" not in bad["note"]
