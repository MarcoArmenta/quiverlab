"""QPA (GAP) as the oracle for the Plan-27 Yoneda / Ext-algebra surface (Tier 2).

Five crosschecks against QPA 1.37, all < 100 ms per call, over a Koszul-discriminating
battery. The Yoneda algebra E(A) = Ext^*_A(A/J, A/J) is realised in QPA as the Ext-algebra
of ``M := DirectSumOfQPAModules(SimpleModules(A))`` (simples in quiver-vertex order):

  * ``ext_algebra_dims``      -- total ``dim E^n`` vs ``ExtAlgebraGenerators(M, n)[1]``.
  * ``ext_generator_degrees`` -- new-generator counts per degree vs ``...[2]`` (degree 0
    maps to ``|Q_0|``, the vertex idempotents QPA counts as the base R = k^{Q_0}).
  * ``ext_quiver``            -- ``dim Ext^1(S_i, S_j)`` corner matrix vs
    ``Length(ExtOverAlgebra(S_i, S_j)[2])`` (pins the corner/direction convention).
  * ``quadratic``             -- ``is_quadratic(A)`` vs ``IsQuadraticIdeal(rels)``.
  * ``koszul_derived``        -- our three-valued ``koszul`` vs the QPA-derivable verdict
    (quadratic AND no degree->=2 generator). QPA ships NO ``IsKoszul``/``KoszulDual``;
    the certifier (G-quadratic / Priddy PBW) is our theory oracle, and QPA validates its
    INPUTS (dims, generator degrees, quadraticity).

The Koszulity discriminator is explicit in the battery: rad^2 = 0 A_3 and the cubic A_4
BOTH have ``dim Ext^2 = 1``, yet A_3's degree-2 class is decomposable (Koszul, gens[2]=0)
while A_4's is a genuine new generator (not Koszul, gens[2]=1). QuantumCI is checked at a
clean GF(32003) point (dim E^n = n+1 regardless of the commutation scalar, matching the
Chouhy-Solotar chain count). Fields: GF(32003) primarily + a QQ spot-check.

qpa-marked: skips locally (no GAP), mandatory under QUIVERLAB_REQUIRE_QPA=1 in CI.
"""
import pytest

from quiverlab import GF, QuantumCI, Quiver, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

F = GF(32003)
_TOP = 5


# --- battery constructors ---------------------------------------------------
def _rad2_a3(field=F):
    # 1 -> 2 -> 3, rad^2 = 0: E = kQ (free); dim Ext^2(S_1,S_3) = 1 is DECOMPOSABLE.
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=field)


def _kronecker2(field=F):
    # the 2-Kronecker path algebra (hereditary): two parallel arrows 1 => 2.
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=field)


def _cubic_a4(field=F):
    # 1 -> 2 -> 3 -> 4 with the cubic relation a*b*c: dim Ext^2(S_1,S_4) = 1 is a genuine
    # NEW generator (the non-Koszul discriminator against rad^2 = 0 A_3).
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=["a*b*c"], field=field)


def _square(field=F):
    # commutative square a*b - c*d: gl.dim 2, Koszul self-dual, E ~= A.
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _line_abc_cde(field=F):
    # Plan-18 multi-vertex diversity record: 1 -> ... -> 6 with rels a*b*c and c*d*e.
    return Quiver([1, 2, 3, 4, 5, 6],
                  {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 5), "e": (5, 6)}
                  ).algebra(relations=["a*b*c", "c*d*e"], field=field)


_BATTERY = [
    ("kx2", truncated_polynomial(2, field=F)),            # E = k[y], Koszul
    ("kx3", truncated_polynomial(3, field=F)),            # E = k[y,z]/(y^2), not Koszul
    ("rad2_a3", _rad2_a3()),                              # Koszul; Ext^2 decomposable
    ("kronecker2", _kronecker2()),                        # hereditary, Koszul
    ("cubic_a4", _cubic_a4()),                            # not Koszul; Ext^2 generator
    ("quantum_ci", QuantumCI("2", field=F)),             # quantum plane, dim E^n = n+1
    ("square", _square()),                               # gl.dim 2, Koszul, E ~= A
    ("line_abc_cde", _line_abc_cde()),                   # multi-vertex, not Koszul
]


# --- the five crosschecks over the battery ----------------------------------
@pytest.mark.parametrize("name,A", _BATTERY)
def test_ext_algebra_dims_crosscheck(name, A):
    A.crosscheck("ext_algebra_dims", _TOP).assert_agree()


@pytest.mark.parametrize("name,A", _BATTERY)
def test_ext_generator_degrees_crosscheck(name, A):
    A.crosscheck("ext_generator_degrees", _TOP).assert_agree()


@pytest.mark.parametrize("name,A", _BATTERY)
def test_ext_quiver_crosscheck(name, A):
    A.crosscheck("ext_quiver").assert_agree()


@pytest.mark.parametrize("name,A", _BATTERY)
def test_quadratic_crosscheck(name, A):
    A.crosscheck("quadratic").assert_agree()


@pytest.mark.parametrize("name,A", _BATTERY)
def test_koszul_derived_crosscheck(name, A):
    A.crosscheck("koszul_derived", _TOP).assert_agree()


# --- the Koszulity discriminator, made explicit -----------------------------
def test_koszul_discriminator_rad2_vs_cubic():
    # Same dim Ext^2 = 1, opposite generator status: rad^2 = 0 A_3 is Koszul (the class
    # is decomposable), the cubic A_4 is not (the class is a new generator).
    rad2_gen = _rad2_a3().crosscheck("ext_generator_degrees", 3).assert_agree()
    cubic_gen = _cubic_a4().crosscheck("ext_generator_degrees", 3).assert_agree()
    assert rad2_gen.ours[2] == 0 and rad2_gen.qpa[2] == 0     # decomposable
    assert cubic_gen.ours[2] == 1 and cubic_gen.qpa[2] == 1   # new generator
    # ... yet dim Ext^2 = 1 for BOTH (QPA side).
    assert _rad2_a3().crosscheck("ext_algebra_dims", 3).qpa[2] == 1
    assert _cubic_a4().crosscheck("ext_algebra_dims", 3).qpa[2] == 1
    # and the derived Koszul verdict splits accordingly.
    assert _rad2_a3().crosscheck("koszul_derived", 3).assert_agree().qpa is True
    assert _cubic_a4().crosscheck("koszul_derived", 3).assert_agree().qpa is False


# --- QQ spot-check (the fields QPA supports exactly are QQ + prime GF(p)) ----
@pytest.mark.parametrize("A", [truncated_polynomial(3, field=QQ), _square(field=QQ)])
def test_ext_algebra_crosscheck_over_QQ(A):
    A.crosscheck("ext_algebra_dims", 4).assert_agree()
    A.crosscheck("ext_generator_degrees", 4).assert_agree()
    A.crosscheck("ext_quiver").assert_agree()
    A.crosscheck("quadratic").assert_agree()
    A.crosscheck("koszul_derived", 4).assert_agree()
