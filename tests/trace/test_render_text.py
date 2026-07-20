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
