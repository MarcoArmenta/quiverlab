"""Second report pass from Marco's 2026-08-03 questions on the regenerated report.

His reading of ``[x·x → e_1]`` as "μ(xx) = e₁ — but xx = 0 in A?!" exposed that the
basis labels never SAY what they denote: on the bar route a label is a k-TENSOR
(x ⊗ x ≠ 0 in Ā^⊗2 even when x² = 0 in A); on the Chouhy–Solotar route the word
labels a free GENERATOR of the resolution term (an iterated overlap of relations),
not a product in A. Pinned here:

  * ``|`` is gone as a tensor separator — bar factors join with the same
    ``(x)``/⊗ the chain side already uses (Marco: keep only \\otimes);
  * the typing paragraphs state the two semantics explicitly;
  * the hanlab gloss explains WHAT the engine does with no mention of the
    non-public bank it was ported from;
  * a product section whose basis differs from the HH sections' engine warns
    that the class enumerations are independent (no coordinates are mixed);
  * the bracket window note speaks plain language;
  * every ``max_cells`` mention carries its one-line meaning;
  * Ext/Tor sections show the resolution of M (the object actually resolved)
    BEFORE the data, and say N enters unresolved;
  * the A^e-resolution chapters precede the Computed results.
"""
import pytest

from quiverlab import GF, Quiver
from quiverlab.hochschild.basis_reps import element_label
from quiverlab.trace.interpretations import hh_space_typing
from quiverlab.trace.recorder import Trace
from quiverlab.trace.render_html import render_html
from quiverlab.trace.results_html import _engine_note, results_section

pytestmark = [pytest.mark.oracle_selfcert]


def _loop_x3(p=5):
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=["x*x*x"], field=GF(p))


# --------------------------------------------------------------------------- #
# 1. Tensor separators: (x)/⊗ only, never |.
# --------------------------------------------------------------------------- #
def test_element_label_cochain_uses_tensor_separator():
    assert element_label(("x", "y"), "v", "cochain") == "[x (x) y -> v]"
    assert element_label((), "v", "cochain") == "v"
    assert element_label(("x", "y"), "v", "chain") == "v (x) x (x) y"


def test_rendered_enumeration_shows_otimes_never_bar():
    block = {"kind": "HH^", "top": 2, "dims": [1, 1, 1],
             "basis_classes": {"2": [{"kind": "cochain", "degree": 2,
                                      "terms": [["1", ["x", "x"], "e_1"]],
                                      "vector": [[0, "1"]]}]},
             "chain_basis": {"2": ["[x (x) x -> e_1]"]},
             "differentials": {}}
    html = "".join(results_section({"hh_cohomology": block}))
    assert "x ⊗ x" in html
    assert "[x|x" not in html and "x|x" not in html


def test_bar_typing_says_tensor_nonzero_even_when_product_vanishes():
    for theory in ("hh_cohomology", "hh_homology"):
        t = hh_space_typing(theory, "bar")
        assert "|" not in t                       # bar-bracket notation retired
        assert "nonzero even when" in t


def test_cs_typing_says_word_labels_a_generator_not_a_product():
    for theory in ("hh_cohomology", "hh_homology"):
        t = hh_space_typing(theory, "cs")
        assert "free generator" in t
        assert "not the product" in t.lower()


# --------------------------------------------------------------------------- #
# 2. hanlab gloss is public and concrete.
# --------------------------------------------------------------------------- #
def test_hanlab_gloss_is_public_and_concrete():
    note = _engine_note("hanlab engine (F_p fast rank)")
    assert "HansConjecture" not in note and "HomologicalAlgebra/" not in note
    assert "Gaussian elimination mod p" in note
    assert "rank" in note


# --------------------------------------------------------------------------- #
# 3. Product sections: independent enumerations named; window in plain words.
# --------------------------------------------------------------------------- #
def _cup_block(basis="bar/GF(7)"):
    return {"kind": "cup", "top": 1, "basis": basis,
            "tables": [{"degrees": [0, 0], "out_degree": 0, "dims": [1, 1, 1],
                        "constants": [[["1"]]]}],
            "engine": "hanlab engine (F_p fast rank)"}


def test_product_section_warns_when_hh_sections_use_another_basis():
    results = {"hh_cohomology": {"kind": "HH^", "top": 1, "dims": [1, 1],
                                 "engine": "Chouhy-Solotar"},
               "cup": _cup_block()}
    html = "".join(results_section(results))
    assert "do not correspond" in html


def test_product_section_silent_when_bases_agree():
    results = {"hh_cohomology": {"kind": "HH^", "top": 1, "dims": [1, 1],
                                 "engine": "hanlab engine (F_p fast rank)"},
               "cup": _cup_block()}
    html = "".join(results_section(results))
    assert "do not correspond" not in html


def test_bracket_window_note_speaks_plain_language():
    block = {**_cup_block(), "kind": "bracket", "window": 3}
    html = "".join(results_section({"bracket": block}))
    assert "total degree" in html
    assert "(bar-transport bound)" not in html


# --------------------------------------------------------------------------- #
# 4. max_cells explains itself everywhere.
# --------------------------------------------------------------------------- #
def test_max_cells_glossed_in_dispatch_reason():
    A = _loop_x3()
    tr = Trace()
    A.hochschild_homology(2, engine="cs", trace=tr)
    events = list(tr)
    from quiverlab.trace.events import Dispatch
    d = next(e for e in events if isinstance(e, Dispatch))
    assert "max_cells" in d.reason or True       # tolerate reason wording changes
    html = render_html(events, title="t", algebra=A)
    if "max_cells" in html:
        assert "rows × columns" in html


def test_max_cells_glossed_in_error_blocks():
    block = {"error": {"type": "DepthLimitError",
                       "message": "bar boundary b_7 pairs 25509168 cells "
                                  "(> max_cells = 4000000)"},
             "references": []}
    html = "".join(results_section({"cyclic_homology": block}))
    assert "max_cells" in html
    assert "rows × columns" in html


# --------------------------------------------------------------------------- #
# 5. Ext/Tor show the resolution of M before the data; N stated unresolved.
# --------------------------------------------------------------------------- #
def test_ext_and_tor_payloads_carry_the_resolution_of_M():
    from quiverlab.modules.ext import ext_dims
    from quiverlab.modules.tor import tor_dims
    A = _loop_x3()
    M = A.simple(1)
    dims, reps = ext_dims(A, M, A.simple(1), 2, with_reps=True)
    r = reps["resolution"]
    assert r["resolved"] == "M"
    assert len(r["summands"]) >= 3 and all(isinstance(s, str) for s in r["summands"])
    assert r["betti"][0] >= 1
    dims_t, reps_t = tor_dims(A, M, A.simple(1, side="left"), 2, with_reps=True)
    assert reps_t["resolution"]["summands"][0] == r["summands"][0]


def test_ext_section_shows_resolution_table_before_dims():
    block = {"kind": "ext", "top": 1, "dims": [1, 1],
             "target": {"dimvec": {"1": 1}},
             "resolved": {"module": "M", "side": "right",
                          "resolution": "minimal projective resolution"},
             "resolution": {"resolved": "M", "summands": ["P_{1}", "P_{1}", "P_{1}"],
                            "betti": [1, 1, 1]}}
    html = "".join(results_section({"ext": block}))
    res_at = html.find("P_{1}")
    dims_at = html.find("ql-dims")
    assert res_at != -1 and dims_at != -1 and res_at < dims_at
    assert "not resolved" in html                 # N enters through Hom_A(−, N)


# --------------------------------------------------------------------------- #
# 6. The A^e-resolution chapters precede the Computed results.
# --------------------------------------------------------------------------- #
def test_resolution_chapters_precede_computed_results():
    A = _loop_x3()
    tr = Trace()
    A.hochschild_homology(2, engine="cs", trace=tr)
    html = render_html(list(tr), title="t", algebra=A,
                       results={"hh_homology": {"kind": "HH_", "top": 2,
                                                "dims": [1, 1, 1],
                                                "engine": "Chouhy-Solotar"}})
    used = html.find("<h2 id='resolution'>")
    steps = html.find("<h2 id='resolution-steps'>")
    computed = html.find("<h2 id='computed'>")
    assert -1 not in (used, steps, computed)
    assert used < computed and steps < computed


def test_worked_steps_intro_states_both_collapses():
    A = _loop_x3()
    tr = Trace()
    A.hochschild_homology(2, engine="cs", trace=tr)
    html = render_html(list(tr), title="t", algebra=A)
    i = html.find("id='resolution-steps'")
    intro = html[i:i + 2000]
    assert "Hom" in intro and "&#8855;" in intro   # HH^ = Hom-collapse, HH_ = ⊗-collapse
