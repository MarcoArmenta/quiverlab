"""QPA (GAP) as the oracle for the Plan-44 tilting + approximation surface (C7).

Probe outcome (verified live 2026-08-05): QPA 1.37 DOES ship a computational tilting and
approximation surface, so this is a genuine crosscheck, NOT the honest skip the plan
anticipated:

  * ``IsTiltingModule`` is a stored PROPERTY (Has/Set), not a computation -- calling it on
    an arbitrary module raises "no method found". The COMPUTATIONAL check is
    ``TiltingModule(M, n)``: it returns ``false`` when ``M`` is not an ``n``-tilting module
    and a ``[true, coresolution]`` list when it is. We crosscheck
    ``is_tilting_module(T).is_tilting == (TiltingModule(T, n) <> false)``.
  * ``MinimalRightAddMApproximation(M, C)`` / ``MinimalLeftAddMApproximation(M, C)`` ARE
    computational; we crosscheck our ``right_add_approximation`` / ``left_add_approximation``
    against the dimension vector of the QPA approximation's source / range.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.approximations import (left_add_approximation,
                                              right_add_approximation)
from quiverlab.modules.morphism import direct_sum
from quiverlab.modules.qpa_module import graded_form
from quiverlab.modules.tilting import is_tilting_module
from quiverlab.qpa import session
from quiverlab.qpa.scripts import module_decl, quiver_and_algebra_script

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def _kA3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=QQ)


def _dv_list(A, M):
    dv = M.dimension_vector()
    return [dv[v] for v in A.quiver.vertices]


def _qpa_tilting(A, T, n, var):
    dv, arr = graded_form(T)
    scr = quiver_and_algebra_script(A) + "\n" + module_decl(A, dv, arr, var)
    scr += "\nt := TiltingModule(%s, %d);" % (var, n)
    return bool(session.run(scr))                         # false -> not tilting; list -> tilting


def _qpa_right_approx_source_dv(A, M, C):
    scr = quiver_and_algebra_script(A) + "\n"
    dvm, arrm = graded_form(M)
    dvc, arrc = graded_form(C)
    scr += module_decl(A, dvm, arrm, "MM") + "\n" + module_decl(A, dvc, arrc, "CC")
    scr += "\nf := MinimalRightAddMApproximation(MM, CC);"
    scr += "\nDimensionVector(Source(f));"
    return list(session.run(scr))


def _qpa_left_approx_range_dv(A, M, C):
    scr = quiver_and_algebra_script(A) + "\n"
    dvm, arrm = graded_form(M)
    dvc, arrc = graded_form(C)
    scr += module_decl(A, dvm, arrm, "MM") + "\n" + module_decl(A, dvc, arrc, "CC")
    # QPA arg-order quirk: RIGHT is (M, C) but LEFT is (C, M) -- verified live 2026-08-05.
    scr += "\ng := MinimalLeftAddMApproximation(CC, MM);"
    scr += "\nDimensionVector(Range(g));"
    return list(session.run(scr))


# --------------------------------------------------------------------------- #
# is_tilting_module vs QPA TiltingModule(T, n)
# --------------------------------------------------------------------------- #
def _tilting_zoo():
    A2 = _kA2()
    A3 = _kA3()
    reg2, _, _ = direct_sum(A2.projective(1), A2.projective(2))
    apr2, _, _ = direct_sum(A2.projective(1), A2.simple(1))
    reg3, _, _ = direct_sum(*[A3.projective(v) for v in (1, 2, 3)])
    return [
        (A2, reg2, True), (A2, A2.projective(1), False), (A2, apr2, True),
        (A3, reg3, True), (A3, A3.projective(1), False),
    ]


@pytest.mark.parametrize("idx", range(5))
def test_is_tilting_crosscheck(idx):
    A, T, expected = _tilting_zoo()[idx]
    ours = is_tilting_module(T, n=1).is_tilting
    theirs = _qpa_tilting(A, T, 1, "TT%d" % idx)
    assert ours == theirs == expected


# --------------------------------------------------------------------------- #
# minimal right / left add(M)-approximation vs QPA
# --------------------------------------------------------------------------- #
def test_right_add_approximation_crosscheck():
    A = _kA2()
    P1, S1 = A.projective(1), A.simple(1)
    f = right_add_approximation(P1, S1)                   # P1^k ->> S1
    ours = _dv_list(A, f.src)
    theirs = _qpa_right_approx_source_dv(A, P1, S1)
    assert ours == theirs                                 # source dimension vector agrees


def test_left_add_approximation_crosscheck():
    A = _kA2()
    P1, S2 = A.projective(1), A.simple(2)
    g = left_add_approximation(P1, S2)                    # S2 >-> P1^k
    ours = _dv_list(A, g.tgt)
    theirs = _qpa_left_approx_range_dv(A, P1, S2)
    assert ours == theirs                                 # target (M^k) dimension vector agrees
