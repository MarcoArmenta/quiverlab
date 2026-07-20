"""Plain-text worked-steps renderer + the HH*(k[x]/x^2) golden trace (spec §3.8, §8).
The binding discipline: every claim in the rendered text equals what the engine
computed (recomputed independently here)."""
import pathlib

import quiverlab
from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace.render_text import render_text, derive_dims
from quiverlab.trace.events import RankStep

GOLDEN = pathlib.Path(__file__).parent / "golden" / "hh_kx2.txt"


def _record():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    return A, tr, table


def _record_homology():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = A.hochschild_homology(2, trace=tr)
    return A, tr, table


def _bar_reference():
    """The bar engine's citation, taken from Plan 06's registry (never hardcoded):
    the `bar` registry key's (bibtex_key, formatted)."""
    bib = quiverlab.bibliography()
    entry = {e.key: e for e in bib}["bar"]
    return entry.bibtex_key, entry.formatted


def test_derive_dims_matches_computed():
    A, tr, table = _record()
    assert derive_dims(list(tr)) == table.dims == [2, 1, 1]


def test_golden_text_trace_matches_with_registry_reference():
    A, tr, table = _record()
    title = "HH^ of " + repr(A).splitlines()[0]
    bibtex_key, formatted = _bar_reference()
    assert bibtex_key == "Hochschild1945"          # stable BibTeX id Plan 06 backs `bar` with
    refs = ((bibtex_key, formatted),)
    text = render_text(list(tr), title=title, references=refs)
    # every line but the registry-sourced citation is golden-fixed; splice the
    # formatted line in from the registry so Plan 06's exact wording is not duplicated.
    expected = GOLDEN.read_text().replace("__BAR_FORMATTED__", formatted)
    assert text == expected
    assert "Hochschild" in formatted


def test_claims_equal_computed_ranks():
    A, tr, table = _record()
    ranks = [e.rank for e in tr if isinstance(e, RankStep)]
    assert ranks == [0, 1, 0]
    # the "rank = k" lines in the text are exactly these
    text = render_text(list(tr))
    assert text.count("rank = 0") == 2 and text.count("rank = 1") == 1


def test_elision_note_rendered(tmp_path):
    from quiverlab.trace.events import RankStep as RS
    ev = [RS(degree=5, side="cochain", nrows=30, ncols=30, rank=3, field="CC",
             matrix=None, elided=True, note="30x30 matrix over CC elided (> 400 cells); rank 3 recorded")]
    text = render_text(ev)
    assert "elided" in text and "rank = 3" in text


def test_homology_derive_dims_uses_the_homology_pairing():
    """Regression for the cohomology-vs-homology dim swap (review C1).

    The engine emits RankStep(degree=n) <-> b_n for HOMOLOGY, so the correct
    pairing is HH_n = C_n - rank(b_n) - rank(b_{n+1}) (rk[n+1], not rk[n-1]).
    k[x]/(x^2) distinguishes the two: the WRONG cohomology formula
    (C_n - rk[n] - rk[n-1]) yields [2, 2, 1] for this trace, while the engine
    (and the correct homology formula) yields [2, 1, 1]."""
    A, tr, table = _record_homology()
    ev = list(tr)
    assert table.dims == [2, 1, 1]
    # derive_dims must agree with the engine (and NOT reproduce the [2, 2, 1] the
    # old cohomology formula produced on these same homology events).
    assert derive_dims(ev) == table.dims == A.hochschild_homology(2).dims
    assert derive_dims(ev) != [2, 2, 1]
    # the rendered Result line must show the CORRECT homology dims, not the swap.
    text = render_text(ev)
    assert "Result: HH_0 = 2   HH_1 = 1   HH_2 = 1" in text
    assert "HH_1 = 2" not in text


def test_verbose_homology_written_document_result_matches_engine(tmp_path, monkeypatch):
    """End-to-end: a verbose hochschild_homology run WRITES a worked-steps
    document whose headline Result dims equal the engine's .dims. Force the no-JS
    HTML fallback (toolchain-agnostic) so the written Result is readable without a
    LaTeX toolchain, then assert its dims are the engine's [2, 1, 1] -- NOT the
    [2, 2, 1] the pre-fix cohomology formula produced (review C1)."""
    import quiverlab.trace.writer as writer

    monkeypatch.setattr(writer, "have_latex", lambda: None)  # always HTML fallback
    monkeypatch.chdir(tmp_path)
    A = truncated_polynomial(2, field=CC)
    expected = A.hochschild_homology(2).dims
    assert expected == [2, 1, 1]
    A.hochschild_homology(2, verbose=True)  # trace is None -> write path fires
    written = list((tmp_path / "quiverlab_traces").glob("HHh_*.html"))
    assert len(written) == 1
    doc = written[0].read_text()
    for i, d in enumerate(expected):
        assert "HH_{%d} = %d" % (i, d) in doc
    assert "HH_{1} = 2" not in doc  # the old-formula value must not appear
