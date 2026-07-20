"""Plan 07 acceptance (spec §3.7, §3.8, §3.9, D9): one algebra, drawn, TikZ'd, and
computed with a worked-steps document whose claims equal the computed values.

Every claim below is INDEPENDENT, not a snapshot: the drawn figure's structure is
checked against the exact layout, the TikZ node coordinates are rebuilt from that
same layout, the worked-steps dimensions are re-derived from the recorded events
(`derive_dims`) and required to equal the engine's own `.dims`, and the References
resolve through Plan 06's `bibliography()`. Reuses the patterns pinned by the
Task 5/6/9/11 tests (`test_draw`, `test_tikz`, `test_render_text`, `test_references`)."""
from fractions import Fraction

import quiverlab
from quiverlab import Quiver, truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace.render_text import derive_dims
from quiverlab.trace.provenance import references_for, resolve_references
from quiverlab.viz.layout import layout


def _tikz_coord(z):
    """Independent oracle for TikZ's exact coordinate formatting (mirrors
    viz.tikz._coord): an integer prints as itself, a Fraction p/q as `{p/q}`."""
    z = Fraction(z)
    if z.denominator == 1:
        return str(z.numerator)
    return "{%d/%d}" % (z.numerator, z.denominator)


def _assert_draw_matches_layout(fig, quiver, relations):
    """The drawn figure has one axes and one text label per vertex, each sitting at
    that vertex's exact layout coordinate (compared via Fraction -- the halves are
    dyadic, so the matplotlib float position is exact; no float literal is used)."""
    assert fig.__class__.__name__ == "Figure"
    assert len(fig.axes) == 1
    ax = fig.axes[0]
    texts = {t.get_text(): t.get_position() for t in ax.texts}
    L = layout(quiver, relations=relations)
    for v, (x, y) in L.positions.items():
        assert str(v) in texts, "no label drawn for vertex %r" % (v,)
        px, py = texts[str(v)]
        assert (Fraction(px), Fraction(py)) == (Fraction(x), Fraction(y))
    return L


def _assert_tikz_shares_layout(src, quiver, relations, L):
    """The TikZ source is one picture whose every vertex node sits at that vertex's
    exact layout coordinate -- i.e. draw() and tikz() are driven by ONE layout."""
    assert r"\begin{tikzpicture}" in src and r"\end{tikzpicture}" in src
    for v, (x, y) in L.positions.items():
        node = "(v%s) at (%s, %s)" % (v, _tikz_coord(x), _tikz_coord(y))
        assert node in src, "TikZ node for vertex %r not at its layout coord: %s" % (v, node)


def test_end_to_end_draw_tikz_and_worked_steps(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    A = truncated_polynomial(2, field=CC)
    relations = A.relations or []

    # draw + tikz are driven by the SAME exact layout (spec §3.7).
    fig = A.draw(file="A.svg")
    assert (tmp_path / "A.svg").exists()
    L = _assert_draw_matches_layout(fig, A.quiver, relations)
    _assert_tikz_shares_layout(A.tikz(), A.quiver, relations, L)

    # Verbose default (D9) writes exactly one worked-steps file (.pdf or .html per
    # the toolchain) to ./quiverlab_traces/.
    quiverlab.verbose = True
    try:
        table = A.hochschild_cohomology(2)
    finally:
        quiverlab.verbose = False
    out = tmp_path / "quiverlab_traces"
    files = list(out.glob("HHc_*"))
    assert files, "verbose run must write a worked-steps file (D9)"
    assert len(files) == 1, "exactly one worked-steps file per computation"
    assert files[0].suffix in (".pdf", ".html")

    # The document's claim is re-derived independently from the recorded events and
    # required to equal the engine's own dimensions (the binding golden discipline):
    # HH^n = collapsed_dim(n) - rank_n - rank_{n-1}.
    tr = Trace()
    again = A.hochschild_cohomology(2, trace=tr)          # explicit trace => NO file
    assert derive_dims(list(tr)) == table.dims == again.dims == [2, 1, 1]
    assert len(list(out.glob("HHc_*"))) == 1              # explicit trace wrote no file

    # The result carries the merged family+engine citation union; here, no family, so
    # exactly ("bar",). Task 11 must NOT overwrite it engine-only.
    assert table.references == A.citations() == ("bar",)

    # The References section resolves those keys through Plan 06's bibliography()
    # (independent: build the expected entry straight from bibliography()).
    bib = {e.key: e for e in quiverlab.bibliography()}
    pairs = resolve_references(references_for(list(tr)))
    assert pairs == ((bib["bar"].bibtex_key, bib["bar"].formatted),)
    assert pairs[0][0] == "Hochschild1945" and "Hochschild" in pairs[0][1]


def test_end_to_end_commuting_square_draw_tikz_shares_layout(tmp_path, monkeypatch):
    """The README's worked example: a non-trivial (multi-vertex, half-integer-row)
    layout drawn and TikZ'd from one shared layout, with the relation shown."""
    monkeypatch.chdir(tmp_path)
    Q = Quiver(vertices=[1, 2, 3, 4],
               arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    relations = ["a*b - c*d"]

    fig = A.draw(file="square.png")
    assert (tmp_path / "square.png").exists() and (tmp_path / "square.png").stat().st_size > 0
    L = _assert_draw_matches_layout(fig, A.quiver, relations)

    src = A.tikz()
    _assert_tikz_shares_layout(src, A.quiver, relations, L)
    # the half-integer rows are the load-bearing case (matplotlib float vs pgf {p/q}).
    assert "(v2) at (1, {1/2})" in src and "(v3) at (1, {-1/2})" in src
    assert "relations: $a*b - c*d$" in src


def test_per_call_verbose_false_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    quiverlab.verbose = True
    try:
        truncated_polynomial(2, field=CC).hochschild_cohomology(2, verbose=False)
    finally:
        quiverlab.verbose = False
    assert not (tmp_path / "quiverlab_traces").exists()
