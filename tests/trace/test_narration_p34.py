"""Plan 34 (post-critique) -- the worked-steps NARRATION is correct and self-checking:

  * MAJOR-4: a LEFT module was narrated with right-sided prose ("the right ... module",
    rad M = M J, m |-> m.a); the narration is now side-aware.
  * MINOR-5: the trace re-derives ranks/nullities in parallel with the engine; a drift
    between the shown pivot-rank/nullity and the engine's reported dimension now raises a
    loud QuiverlabError (the honest move) rather than printing a wrong number.
  * MINOR-7: the tau trace now shows coker(d_1^*) = Tr M (dim/vector) and the application
    of D, and the tau^- branch shows its transposed presentation matrix symmetrically.
  * MINOR-8: every displayed reduced matrix carries a one-line "column-reduce: pivots at
    column(s) ..., rank ..." summary (the homework-acceptable middle ground).
  * MINOR-9: a char-p decomposition the certificate cannot settle refuses out loud, and
    the trace narrates the honest refusal (never a silent verdict)."""
import pytest

from quiverlab import Quiver, GF, truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.trace import modules as TM
from quiverlab.trace.modules import (
    trace_module_report, trace_radical, trace_socle, trace_top, trace_tau,
    trace_decompose, trace_ext,
)
from quiverlab.trace.render_text import render_text


def _kA2():
    return Quiver(vertices=[1, 2], arrows={"a": (1, 2)}).algebra(relations=[], field=GF(2))


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
                  ).algebra(relations=["a*b - c*d"], field=GF(2))


# --------------------------------------------------------------------------- #
# MAJOR-4: left modules are narrated left-sided.
# --------------------------------------------------------------------------- #
def test_left_module_narration_is_side_aware():
    A = _kA2()
    M = A.simple(1, side="left")
    M.name = "M"
    assert M.side == "left"
    ev = trace_module_report(A, M, top=2, with_tau=False, with_decompose=False)
    txt = render_text(ev, title="t")
    intro = txt.split("compute")[0]
    assert "left" in intro and "right" not in intro         # "the left A-module M ..."
    assert "rad M = J M" in txt                             # left radical (J M, not M J)
    assert "m |-> a.m" in txt                               # left action a.m, not m.a
    assert "J m = 0" in txt                                 # left socle annihilator
    assert "rad M = M J" not in txt and "m |-> m.a" not in txt


def test_right_module_narration_unchanged():
    A = _kA2()
    M = A.simple(1)                                          # right (default)
    M.name = "M"
    ev = trace_module_report(A, M, top=2, with_tau=False, with_decompose=False)
    txt = render_text(ev, title="t")
    assert "the right " in txt.split("compute")[0]
    assert "rad M = M J" in txt and "m |-> m.a" in txt


# --------------------------------------------------------------------------- #
# MINOR-6: the radical justification is A/J semisimple, not "Nakayama's lemma".
# --------------------------------------------------------------------------- #
def test_radical_justification_is_semisimple_quotient_not_nakayama():
    A = _square()
    M = A.projective(1).radical()
    M.name = "M"
    txt = render_text(trace_radical(M)[0], title="t")
    assert "A/J is semisimple" in txt and "annihilates every simple" in txt
    assert "Nakayama" not in txt


# --------------------------------------------------------------------------- #
# MINOR-5: the shown rank/nullity is asserted against the engine (loud on drift).
# --------------------------------------------------------------------------- #
# The guard compares the trace's SHOWN pivot-rank/nullity against the engine's reported
# dimension. We force drift by corrupting the DIMENSION the engine reports (returning a
# module whose .dim is off by one) -- NOT by patching the shared linalg (which the engine
# itself uses, so that would break the computation rather than exercise the guard).
def test_radical_drift_guard_raises_loudly(monkeypatch):
    A = _square()
    M = A.projective(1).radical()
    M.name = "M"
    trace_radical(M)                                        # no drift -> fine
    real_radical = type(M).radical

    def wrong_dim(self):
        R = real_radical(self)
        R.dim += 1                                          # engine now reports a wrong dim
        return R
    monkeypatch.setattr(type(M), "radical", wrong_dim)
    with pytest.raises(QuiverlabError, match="drift"):
        trace_radical(M)


def test_socle_drift_guard_raises_loudly(monkeypatch):
    A = _square()
    M = A.projective(1).radical()
    M.name = "M"
    trace_socle(M)
    real_socle = type(M).socle

    def wrong_dim(self):
        S = real_socle(self)
        S.dim += 1
        return S
    monkeypatch.setattr(type(M), "socle", wrong_dim)
    with pytest.raises(QuiverlabError, match="disagrees"):
        trace_socle(M)


def test_top_drift_guard_raises_loudly(monkeypatch):
    A = _square()
    M = A.projective(1).radical()
    M.name = "M"
    trace_top(M)
    real_top = type(M).top

    def wrong_dim(self):
        T = real_top(self)
        T.dim += 1                                          # T.dim != M.dim - dim rad now
        return T
    monkeypatch.setattr(type(M), "top", wrong_dim)
    with pytest.raises(QuiverlabError, match="disagrees"):
        trace_top(M)


def test_ext_drift_guard_matches_engine():
    A = _kA2()
    from quiverlab.modules.ext import ext_dims
    ev, dims = trace_ext(A, A.simple(1), A.simple(2), 3)    # asserts internally vs ext_dims
    assert dims == ext_dims(A, A.simple(1), A.simple(2), 3)


# --------------------------------------------------------------------------- #
# MINOR-7: tau shows coker(d_1^*) = Tr M + D; tau^- shows a transposed matrix.
# --------------------------------------------------------------------------- #
def test_tau_shows_transpose_cokernel_and_duality():
    A = _square()
    M = A.projective(1).radical()
    M.name = "M"
    txt = render_text(trace_tau(M, "tau")[0], title="t")
    assert "Tr M = coker(d_1^*) has dimension" in txt       # the cokernel step (dim shown)
    assert "duality: D = Hom_k(-,k)" in txt                 # the application of D
    assert "d_{1}^{*}" in txt                               # the transposed matrix appears


def test_tau_minus_shows_presentation_and_transpose_symmetrically():
    A = _square()
    M = A.projective(1).radical()                           # tau^- M nonzero here (dim 3)
    M.name = "M"
    ev, t = trace_tau(M, "tau_minus")
    assert t.dim > 0
    txt = render_text(ev, title="t")
    assert "inverse AR translate" in txt                    # the phrase is kept
    assert "Projective presentation" in txt                 # D M's presentation is shown
    assert "d_{1}^{*}" in txt                               # ...and its transpose matrix
    assert "Tr(D" in txt                                    # the cokernel = tau^- directly


# --------------------------------------------------------------------------- #
# MINOR-8: every reduced matrix carries a pivots/rank one-liner, same phrasing.
# --------------------------------------------------------------------------- #
def test_reduction_summary_phrasing_is_consistent():
    A = _square()
    M = A.projective(1).radical()
    M.name = "M"
    rad = render_text(trace_radical(M)[0], title="t")
    soc = render_text(trace_socle(M)[0], title="t")
    tau = render_text(trace_tau(M, "tau")[0], title="t")
    assert "Column-reducing G over GF(2): pivots at column(s)" in rad
    assert "Column-reducing K over GF(2): pivots at column(s)" in soc
    assert "Column-reducing d_1* over GF(2): pivots at column(s)" in tau


# --------------------------------------------------------------------------- #
# MINOR-9: a char-p decomposition the certificate can't settle refuses out loud.
# --------------------------------------------------------------------------- #
def test_decompose_char_p_refusal_is_narrated():
    # k[x]/(x^2) over GF(2): the regular module is indecomposable but char 2 <= dim 2, so
    # the trace-form certificate is unreliable and no Fitting split exists -> honest refusal.
    A = truncated_polynomial(2, field=GF(2))
    M = A.projective(1)
    M.name = "M"
    ev, result = trace_decompose(M)
    assert result is None                                   # the engine refused loudly
    txt = render_text(ev, title="t")
    assert "honest refusal" in txt
    assert "could not be certified within budget" in txt
