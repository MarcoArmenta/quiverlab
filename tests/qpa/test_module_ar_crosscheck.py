"""QPA (GAP) as the oracle for the Plan-23 module surface (Tier 1b item 4).

Crosschecks against QPA's DTr / TrD / ProjectiveResolution / DualOfModule /
InjDimensionOfModule across a zoo incl. the Plan-18 multi-vertex line algebra:
tau/tau^- dimension vectors + isomorphism class, projective & injective resolution
term dimension vectors, and injective dimension. qpa-marked: skips locally,
mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import Quiver, GF, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.qpa import session
from quiverlab.qpa.scripts import module_decl, quiver_and_algebra_script

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _square(field=QQ):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _cn(n, field=QQ):
    verts = list(range(1, n + 1))
    arrows = {f"a{i}": (i, i % n + 1) for i in range(1, n + 1)}
    rels = [f"a{i}*a{i % n + 1}" for i in range(1, n + 1)]
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


def _line_abc_cde(field=QQ):
    # Plan-18 multi-vertex record: 1->2->...->6 with rels a*b*c and c*d*e.
    Q = Quiver([1, 2, 3, 4, 5, 6],
               {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 5), "e": (5, 6)})
    return Q.algebra(relations=["a*b*c", "c*d*e"], field=field)


def _indec_modules(A):
    return ([A.simple(v) for v in A.quiver.vertices]
            + [A.projective(v) for v in A.quiver.vertices]
            + [A.injective(v) for v in A.quiver.vertices])


# --- graded-form emitter sanity (our module -> QPA -> IsomorphicModules) ----
def test_graded_form_roundtrips_projective():
    A = linear_path_algebra(3, field=QQ)
    from quiverlab.modules.qpa_module import graded_form
    dv, arr = graded_form(A.projective(1))
    base = quiver_and_algebra_script(A) + "\n" + module_decl(A, dv, arr, "M")
    base += "\nP := IndecProjectiveModules(A)[1];;"
    assert bool(session.run(base + "\nIsomorphicModules(M, P);"))


# --- tau / tau^- : dimension vectors + isomorphism class --------------------
@pytest.mark.parametrize("A", [linear_path_algebra(2, field=QQ),
                               linear_path_algebra(3, field=QQ),
                               _square(), _cn(3), _line_abc_cde()])
def test_tau_crosscheck(A):
    for M in _indec_modules(A):
        A.crosscheck("tau", M).assert_agree()
        A.crosscheck("tau_minus", M).assert_agree()


@pytest.mark.parametrize("field", [GF(2), GF(3)])
def test_tau_dimvec_crosscheck_over_gfp(field):
    A = linear_path_algebra(3, field=field)
    for M in _indec_modules(A):
        rep = A.crosscheck("tau", M)
        assert rep.ours == rep.qpa, rep       # dimension vectors match over GF(p)


# --- projective / injective resolution terms --------------------------------
@pytest.mark.parametrize("A", [linear_path_algebra(3, field=QQ), _square(),
                               _line_abc_cde()])
def test_proj_resolution_crosscheck(A):
    for v in A.quiver.vertices:
        A.crosscheck("proj_resolution", A.simple(v), 4).assert_agree()


@pytest.mark.parametrize("A", [linear_path_algebra(3, field=QQ), _square(),
                               _line_abc_cde()])
def test_inj_resolution_crosscheck(A):
    for v in A.quiver.vertices:
        A.crosscheck("inj_resolution", A.simple(v), 4).assert_agree()


# --- injective dimension (finite AND infinite/self-injective) ---------------
@pytest.mark.parametrize("A", [linear_path_algebra(3, field=QQ), _square(),
                               _line_abc_cde()])
def test_inj_dimension_crosscheck(A):
    for v in A.quiver.vertices:
        A.crosscheck("inj_dimension", A.simple(v), 8).assert_agree()


def test_inj_dimension_infinite_selfinjective():
    A = _cn(3)                                # self-injective => non-proj has inf inj.dim
    rep = A.crosscheck("inj_dimension", A.simple(1), 6)
    assert rep.ours is None and rep.qpa is None    # both report infinite
    rep.assert_agree()
