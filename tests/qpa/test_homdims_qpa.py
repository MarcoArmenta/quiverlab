"""QPA (GAP) as an EXTERNAL oracle for the Plan-40 C6 homological-dimension family.

QPA 1.37 ships bound-parametrised ``GlobalDimensionOfAlgebra(A, n)`` /
``DominantDimensionOfAlgebra(A, n)`` / ``GorensteinDimensionOfAlgebra(A, n)``, each
returning an int or GAP ``infinity`` -- so global / dominant / Gorenstein dimensions
get a genuine external certificate (``crosscheck_global_dimension`` /
``crosscheck_dominant_dimension`` / ``crosscheck_gorenstein``, all computing the QPA
side with ITS OWN homological machinery).

QPA has NO Igusa-Todorov surface: a live ``NamesGVars()`` probe finds no ``Igusa`` /
``Todorov`` / ``phiDimension`` / ``psiDimension`` name (verified 2026-08-05). So phi/psi
have no QPA oracle -- their coverage is the Task-B literature battery
(``tests/modules/test_igusa_todorov.py``: phi = pd for finite pd, the Barrios-Mata
self-injective all-zero closed form, projective additivity). This file's IT probe SKIPS
that comparison honestly and FAILS only if QPA ever ships an IT surface.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import GF, Quiver, truncated_polynomial
from quiverlab.qpa import session
from quiverlab.qpa.crosscheck import (crosscheck_dominant_dimension,
                                      crosscheck_global_dimension,
                                      crosscheck_gorenstein)

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

F = GF(7)


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=F)


def _kA3_rel():
    # kA3 / (a*b): Nakayama, gl.dim 2, dominant dimension 2.
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=F)


def _line_abc_cde():
    # Plan-18 multi-vertex diversity record: 1 -> ... -> 6, rels a*b*c, c*d*e.
    return Quiver([1, 2, 3, 4, 5, 6],
                  {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 5), "e": (5, 6)}
                  ).algebra(relations=["a*b*c", "c*d*e"], field=F)


def _kx3():
    # k[x]/(x^3): self-injective -- gl.dim and dominant dimension both INFINITE
    # (QPA `infinity` <-> our unresolved/infinite marker, both -> None), Gorenstein 0.
    return truncated_polynomial(3, field=F)


_ZOO = [("kA2", _kA2()), ("kA3_rel", _kA3_rel()),
        ("line_abc_cde", _line_abc_cde()), ("kx3", _kx3())]


# --- IT probe: honest skip, FAILS if QPA ever ships an Igusa-Todorov surface ---
_IT_NAMES = ("IgusaTodorovFunction", "PhiDimension", "PsiDimension",
             "IgusaTodorov", "phiDim", "psiDim")


def test_qpa_exposes_no_igusa_todorov_surface():
    lg = session.libgap_handle()
    # (1) scan the global name table FIRST (IsBoundGlobal REGISTERS the queried name
    #     into NamesGVars(), so scanning after would echo the queries back).
    names = [str(n) for n in lg.eval("NamesGVars()")]
    it_like = sorted(n for n in names
                     if "igusa" in n.lower() or "todorov" in n.lower()
                     or "phidimension" in n.lower() or "psidimension" in n.lower())
    assert it_like == [], (
        "QPA now exposes Igusa-Todorov-like names %s -- wire a crosscheck for phi/psi "
        "and drop this honest-skip guard" % it_like)
    # (2) the specific entry points we would use are unbound (callable check).
    present = [nm for nm in _IT_NAMES if bool(lg.eval(f'IsBoundGlobal("{nm}")'))]
    assert present == [], (
        "QPA now binds Igusa-Todorov entry point(s) %s -- add the phi/psi crosscheck "
        "(the literature battery is no longer the only oracle)" % present)


@pytest.mark.parametrize("name,A", _ZOO)
def test_global_dimension_matches_qpa(name, A):
    crosscheck_global_dimension(A).assert_agree()


@pytest.mark.parametrize("name,A", _ZOO)
def test_dominant_dimension_matches_qpa(name, A):
    crosscheck_dominant_dimension(A).assert_agree()


@pytest.mark.parametrize("name,A", _ZOO)
def test_gorenstein_matches_qpa(name, A):
    crosscheck_gorenstein(A).assert_agree()
