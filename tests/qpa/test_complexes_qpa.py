"""QPA (GAP) as an EXTERNAL oracle for the Plan-39 complex surface.

QPA 1.37 SHIPS the Ch.10 complex package (probed live below): ``StalkComplex``,
``FiniteComplex``, ``HomologyOfComplex``, ``Shift``, ``FiniteChainMap``,
``MappingCone``. We build the SAME bounded complexes in both systems and compare
homology dimensions computed by QPA's ``HomologyOfComplex`` against our
``ChainComplex.homology_dims``:

  * a stalk complex (homology = the module);
  * a bounded 2-term complex ``[M --f--> N]`` -- which IS our mapping cone of the
    degree-0 chain map ``stalk(M) -> stalk(N)`` -- via QPA ``FiniteComplex``; the
    canonical (``dim Hom = 1``) map is used, so the two systems build the same
    complex up to a scalar and the homology dims agree exactly;
  * ``Shift`` bookkeeping, with the documented CONVENTION FLIP: QPA's
    ``Shift(C, k)`` sends homology at degree 0 to degree ``-k`` (verified live),
    the OPPOSITE of our ``X[k]`` (degree ``+k``), so our ``shift(k)`` is
    cross-checked against QPA ``Shift(C, -k)``.

HONEST FALLBACK (the plan anticipated it): QPA's ``MappingCone`` returns a complex
object that ``HomologyOfComplex`` refuses through libgap ("no method found on 2
arguments" -- a category/stream representation ``HomologyOfComplex`` does not
accept when scripted one statement at a time). So the mapping cone is cross-checked
via the EQUIVALENT explicit ``FiniteComplex`` (the 2-term complex our cone equals),
recorded verbatim in ``test_qpa_complex_surface_probe``. This is genuine coverage of
our cone/homology, never a silent skip. Kept test-only (read-only script/module
helpers), mirroring ``test_tor_qpa.py``.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex, ChainMap
from quiverlab.modules.qpa_module import graded_form
from quiverlab.qpa import scripts, session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


# --------------------------------------------------------------------------- #
# Helpers: declare modules in a GAP base script, run homology-list queries.
# --------------------------------------------------------------------------- #
def _base(A, named_modules):
    """A GAP base script declaring the algebra (as ``A``), the category
    ``cat := CatOfRightAlgebraModules(A)``, and each named module."""
    s = scripts.quiver_and_algebra_script(A)
    for name, M in named_modules.items():
        dv, arr = graded_form(M)
        s += "\n" + scripts.module_decl(A, dv, arr, name)
    s += "\ncat := CatOfRightAlgebraModules(A);;"
    return s


def _homology_list(base, build_lines, degs):
    """``{i: dim H_i}`` from a GAP script that ends by binding the complex ``C``
    (``build_lines``) then reading ``Dimension(HomologyOfComplex(C, i))`` for each
    ``i`` in ``degs``."""
    listexpr = "[" + ",".join(
        "Dimension(HomologyOfComplex(C,%d))" % i for i in degs) + "];"
    out = session.run(base + "\n" + build_lines + "\n" + listexpr)
    return {i: int(v) for i, v in zip(degs, out)}


def _kA3():
    return linear_path_algebra(3, field=QQ)


# --------------------------------------------------------------------------- #
# The probe: record exactly what QPA exposes and what forced the fallback.
# --------------------------------------------------------------------------- #
def test_qpa_complex_surface_probe():
    lg = session.libgap_handle()
    for name in ("Complex", "FiniteComplex", "StalkComplex", "HomologyOfComplex",
                 "Shift", "FiniteChainMap", "MappingCone"):
        assert bool(lg.IsBoundGlobal(name)), "QPA is missing %s" % name

    # HomologyOfComplex accepts a StalkComplex (2-arg form works) ...
    A = _kA3()
    base = _base(A, {"M1": A.simple(1)})
    ok = session.run(base + "\nsc := StalkComplex(cat, M1, 0);;\n"
                     "Dimension(HomologyOfComplex(sc, 0));")
    assert int(ok) == 1

    # ... but it does NOT accept QPA's MappingCone object through libgap (the
    # documented fallback reason). Assert the fallback still holds so this note is
    # revisited if a future QPA makes the cone directly homology-scriptable.
    from quiverlab.errors import QuiverlabError
    cone_build = ("hh := HomOverAlgebra(M1, M1);;\nf := hh[1];;\n"
                  "cP := StalkComplex(cat, M1, 0);;\ncS := StalkComplex(cat, M1, 0);;\n"
                  "phi := FiniteChainMap(cP, cS, 0, [f]);;\n"
                  "cone := MappingCone(phi);;\nDimension(HomologyOfComplex(cone, 0));")
    scriptable = True
    try:
        session.run(base + "\n" + cone_build)
    except Exception:                    # noqa: BLE001 -- any GAP error = not scriptable
        scriptable = False
    assert not scriptable, (
        "QPA's MappingCone became homology-scriptable -- wire a DIRECT cone "
        "crosscheck (HomologyOfComplex on the cone) and drop the FiniteComplex "
        "fallback used by test_cone_homology_vs_qpa.")


# --------------------------------------------------------------------------- #
# (a) stalk homology.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("v", [1, 2, 3])
def test_stalk_homology_vs_qpa(v):
    A = _kA3()
    M = A.simple(v)
    base = _base(A, {"M": M})
    degs = [-1, 0, 1]
    qpa = _homology_list(base, "C := StalkComplex(cat, M, 0);;", degs)
    ours = ChainComplex.stalk(M, 0).homology_dims()
    for i in degs:
        assert qpa[i] == ours.get(i, 0), (v, i, qpa, ours)


# --------------------------------------------------------------------------- #
# (b) mapping cone == bounded 2-term complex [M --f--> N] (degrees 1, 0),
#     cross-checked against QPA FiniteComplex(cat, 1, [f]).
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("src_kind,tgt_kind,label", [
    ("projective", "simple", "P1 ->> S1 (surjection: H_1 = rad P1)"),
    ("simple", "projective", "S3 >-> P1 (injection: H_0 = coker)"),
    ("projective", "projective", "P1 -> P1 (iso: acyclic cone)"),
])
def test_cone_homology_vs_qpa(src_kind, tgt_kind, label):
    A = _kA3()
    # pick vertices so dim Hom(src, tgt) == 1 (canonical map): P1/S1, S3/P1, P1/P1
    verts = {("projective", "simple"): (1, 1),
             ("simple", "projective"): (3, 1),
             ("projective", "projective"): (1, 1)}[(src_kind, tgt_kind)]
    M = getattr(A, src_kind)(verts[0])
    N = getattr(A, tgt_kind)(verts[1])
    hb = A.hom_basis(M, N)
    assert len(hb) == 1, (label, "expected a canonical 1-dim Hom", len(hb))
    f = hb[0].matrix
    # our mapping cone of the degree-0 chain map stalk(M) -> stalk(N)
    cone = ChainMap(ChainComplex.stalk(M, 0), ChainComplex.stalk(N, 0), {0: f}).cone()
    ours = cone.homology_dims()
    # the SAME complex in QPA: FiniteComplex(cat, 1, [f]) = [M(deg1) --f--> N(deg0)]
    base = _base(A, {"MM": M, "NN": N})
    degs = [-1, 0, 1, 2]
    qpa = _homology_list(
        base, "hh := HomOverAlgebra(MM, NN);;\nf := hh[1];;\n"
              "C := FiniteComplex(cat, 1, [f]);;", degs)
    for i in degs:
        assert qpa[i] == ours.get(i, 0), (label, i, qpa, ours)


# --------------------------------------------------------------------------- #
# (c) Shift bookkeeping -- QPA Shift(C, -k) == our X[k] (opposite convention).
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("k", [1, 2, -1])
def test_shift_bookkeeping_vs_qpa(k):
    A = _kA3()
    M = A.simple(2)
    base = _base(A, {"M": M})
    degs = [-3, -2, -1, 0, 1, 2, 3]
    # our shift by k moves homology from degree 0 to degree +k
    ours = ChainComplex.stalk(M, 0).shift(k).homology_dims()
    # QPA's Shift(C, k) moves it to degree -k (verified live), so use -k here
    qpa = _homology_list(base,
                         "sc := StalkComplex(cat, M, 0);;\nC := Shift(sc, %d);;" % (-k),
                         degs)
    for i in degs:
        assert qpa[i] == ours.get(i, 0), (k, i, qpa, ours)
    assert ours.get(k) == M.dim          # our convention: homology lands at +k
