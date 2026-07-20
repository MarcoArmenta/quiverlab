"""Renderer parity (spec §8): the SAME events produce CONSISTENT claims in all
three renderers -- identical resulting dims and identical per-degree ranks. The
binding discipline holds regardless of output format."""
import re

from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace.events import RankStep
from quiverlab.trace.render_text import render_text, derive_dims
from quiverlab.trace.render_html import render_html
from quiverlab.trace.render_latex import render_latex

# A generic (bibtex_key, formatted) fixture: parity is about the renderers agreeing,
# not about any specific citation, so this stays decoupled from Plan 06's registry.
REFS = (("Refkey2020", "A. Author, A Journal 1 (2020), 1-2."),)


def _events():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    return list(tr), table


def _homology_events():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = A.hochschild_homology(2, trace=tr)
    return list(tr), table


def _dims_in(s):
    return [int(m) for m in re.findall(r"HH[\^_]?\{?\d+\}?\s*=\s*(\d+)", s)]


def test_all_three_renderers_agree_on_dims():
    ev, table = _events()
    txt = render_text(ev, title="t", references=REFS)
    html = render_html(ev, title="t", references=REFS)
    tex = render_latex(ev, title="t", references=REFS)
    assert derive_dims(ev) == table.dims == [2, 1, 1]
    assert _dims_in(txt) == [2, 1, 1]
    assert _dims_in(html) == [2, 1, 1]
    assert _dims_in(tex) == [2, 1, 1]


def test_all_three_renderers_agree_on_homology_dims():
    """Homology mirror of the cohomology parity test (review C1): all three
    renderers must derive the engine's HOMOLOGY dims. For k[x]/(x^2) the correct
    dims are [2, 1, 1]; the pre-fix cohomology formula (rk[n-1]) yielded [2, 2, 1]
    on these homology events, so this test fails under the old derive_dims."""
    ev, table = _homology_events()
    txt = render_text(ev, title="t", references=REFS)
    html = render_html(ev, title="t", references=REFS)
    tex = render_latex(ev, title="t", references=REFS)
    assert derive_dims(ev) == table.dims == [2, 1, 1]
    assert _dims_in(txt) == [2, 1, 1]
    assert _dims_in(html) == [2, 1, 1]
    assert _dims_in(tex) == [2, 1, 1]


def test_all_three_render_the_same_ranks():
    ev, table = _events()
    ranks = [e.rank for e in ev if isinstance(e, RankStep)]
    assert ranks == [0, 1, 0]
    for render in (render_text, render_html, render_latex):
        s = render(ev, title="t", references=REFS)
        # every recorded rank appears as a "rank = k" / "rank} = k" claim
        found = [int(m) for m in re.findall(r"rank[^\d=]*=\s*(\d+)", s)]
        assert found == ranks


def test_all_three_carry_the_reference():
    ev, table = _events()
    for render in (render_text, render_html, render_latex):
        assert "Refkey2020" in render(ev, title="t", references=REFS)
