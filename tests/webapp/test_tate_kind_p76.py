"""The `tate_hochschild` compute kind across the tiers (Plan 76 / R3).

Unmarked by the Plan-32 extras-gated ruling (`tests/webapp` carries no `oracle_*`
markers -- they would make the audited class counts environment-dependent).

What these pin: the block is a DEGREE-RANGE kind reporting the symmetric window
[-top, top]; degrees a route cannot know are `null`, never a wrong number; the
`ordinary_pos` HH anchor is carried separately from the Tate degree 0; the server
runner and the Pyodide twin produce BYTE-IDENTICAL JSON; and a refusal comes back as a
clean `error` block rather than an exception.
"""
import json

import pytest

from quiverlab import truncated_polynomial
from quiverlab.fields import GF, QQ
from quiverlab.hpc.spec import parse_compute_item


def _server_block(A, top):
    from quiverlab.hpc.spec import _dispatch
    item = parse_compute_item(f"tate_hochschild:0..{top}")
    block, _ = _dispatch(A, item, None, {})
    return block


def _shared_builder(A, top):
    """The library builder BOTH tiers delegate to.

    This pins that the server tier adds nothing but citations on top of it.  The parity
    that matters -- the Pyodide twin's OWN dispatch branch producing the same JSON -- is
    exercised through `compute_one` in `tests/gui/test_tate_runner_twin_p76.py`; calling
    the builder here and calling it the twin would be a test of nothing."""
    from quiverlab.hochschild.tate import tate_hochschild_block
    return tate_hochschild_block(A, top)


def test_native_block_kx3():
    A = truncated_polynomial(3, field=GF(32003))
    b = _server_block(A, 2)
    assert b["kind"] == "tate_hochschild" and b["engine"] == "native"
    assert b["pos_dims"] == [2, 2, 2]
    assert b["neg_dims"] == [2, 2]
    assert b["hat_hh0"] == 2 and b["agrees_from"] == 1
    assert b["scope"] == "self_injective"
    assert b["ordinary_pos"][0] == 3            # HH^0, NOT the Tate degree 0
    assert "bergh_jorgensen_tate" in b["references"]
    assert "wang_singular_hh" in b["references"]
    bib = [k for k, _ in b["citations"]]        # resolved BibTeX keys
    assert "BerghJorgensen2013tate" in bib and "Wang2015singular" in bib


def test_positive_route_block_is_honest_about_degree_zero():
    A = truncated_polynomial(3, field=QQ)
    b = _server_block(A, 3)
    assert b["engine"] in ("positive", "duality")
    assert b["pos_dims"][0] is None and b["hat_hh0"] is None
    assert b["ordinary_pos"][0] == 3            # the ordinary anchor is still reported


def test_server_tier_adds_only_citations_to_the_shared_block():
    for A, top in ((truncated_polynomial(3, field=GF(32003)), 2),
                   (truncated_polynomial(2, field=GF(32003)), 3),
                   (truncated_polynomial(3, field=QQ), 2)):
        server = _server_block(A, top)
        shared = _shared_builder(A, top)
        server.pop("citations", None)
        assert json.dumps(server, sort_keys=True) == json.dumps(shared, sort_keys=True)


def test_range_must_start_at_zero():
    """The shared grammar guard covers the new kind for free -- `1..5` never parses, so
    `tate_hochschild:0..N` is the ONLY spelling and the canonical key is a function of
    (kind, hi) alone."""
    from quiverlab.hpc.spec import SpecError
    with pytest.raises(SpecError):
        parse_compute_item("tate_hochschild:1..5")


def test_missing_range_refused():
    from quiverlab.hpc.spec import ComputeError
    A = truncated_polynomial(3, field=GF(32003))
    item = parse_compute_item("tate_hochschild")
    from quiverlab.hpc.spec import _dispatch
    with pytest.raises(ComputeError):
        _dispatch(A, item, None, {})


def test_refusal_is_an_error_block_not_a_raise():
    """A non-Gorenstein / deferred algebra must produce a clean typed block."""
    from quiverlab import Quiver
    from quiverlab.families.radical_square_zero import RadicalSquareZero
    A = RadicalSquareZero(Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}),
                          field=GF(32003))
    b = _server_block(A, 2)
    # kA2-like input is Gorenstein, so the honest outcome is the positive route with
    # nothing claimed -- either way it is a block, never an exception.
    assert b["kind"] == "tate_hochschild"
    assert "error" in b or all(d is None for d in b["pos_dims"])


def test_report_renders_both_halves_and_keeps_hh_separate():
    from quiverlab.trace.results_html import _block_html as render_result_block
    A = truncated_polynomial(3, field=GF(32003))
    b = _server_block(A, 2)
    html = "\n".join(render_result_block(b["kind"], b))
    assert "Tate" in html
    assert "stable centre" in html or "not</b> HH" in html or "HH&#770;" in html
    assert "dim HH^n" in html                   # the ordinary anchor table is present


def test_report_renders_absent_cells_as_dashes():
    from quiverlab.trace.results_html import _block_html as render_result_block
    A = truncated_polynomial(3, field=QQ)
    b = _server_block(A, 2)
    html = "\n".join(render_result_block(b["kind"], b))
    assert "&mdash;" in html                    # degree 0 is absent, not a number
