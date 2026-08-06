"""Regression gate for the 2026-08-06 derived_fingerprint report crash.

``_derived_fingerprint_html`` built its chunk list but never returned it, so
``results_section`` did ``out.extend(None)`` and the whole worked-steps HTML
bundle died for ANY request containing a derived_fingerprint block (the JSON
record survived; the report did not).  Two gates:

* the live render: a derived_fingerprint block flows through
  ``results_section`` and lands in the page;
* the class gate: every per-kind ``*_html`` helper in ``results_html.py``
  ends in an explicit ``return`` — a fall-through renderer returns ``None``
  and kills the report, so it must fail here first.
"""

import ast
import inspect

import pytest

from quiverlab.trace import results_html

pytestmark = pytest.mark.oracle_selfcert


def _fingerprint_results():
    return {
        "derived_fingerprint": {
            "fingerprint": {
                "coxeter_polynomial": "t^3 + t^2 + t + 1",
                "cartan_det": 1,
                "cartan_smith": [1, 1, 1],
                "hh_cohomology_dims": [1, 0, 0],
                "hh_homology_dims": [3, 0, 0],
                "cyclic_dims": {"error": "max_cells exceeded"},
                "center_dim": 1,
                "gl_dim": 1,
            },
        },
    }


def test_derived_fingerprint_block_returns_chunks():
    chunks = results_html._block_html(
        "derived_fingerprint", _fingerprint_results()["derived_fingerprint"])
    assert isinstance(chunks, list) and chunks
    page = "".join(chunks)
    assert "Coxeter polynomial" in page
    assert "t^3 + t^2 + t + 1" in page
    # the captured-error cell renders its message instead of vanishing
    assert "max_cells exceeded" in page


def test_derived_fingerprint_renders_through_results_section():
    page = "".join(results_html.results_section(_fingerprint_results()))
    assert "derived" in page and "Coxeter polynomial" in page


def test_every_renderer_helper_ends_with_return():
    tree = ast.parse(inspect.getsource(results_html))
    helpers = [node for node in tree.body
               if isinstance(node, ast.FunctionDef) and node.name.endswith("_html")]
    assert len(helpers) >= 20, "renderer helpers went missing — wrong module?"
    fallthrough = [(fn.name, fn.body[-1].lineno) for fn in helpers
                   if not isinstance(fn.body[-1], ast.Return)]
    assert not fallthrough, (
        "renderer helpers whose last statement is not a return (a fall-through "
        "returns None and results_section dies on out.extend(None)): %r"
        % fallthrough)
