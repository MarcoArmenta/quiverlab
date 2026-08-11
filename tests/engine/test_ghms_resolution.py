"""The GHMS comultiplicative minimal A^e-resolution of a Koszul algebra: terms
P_n = A (x)_S K_n (x)_S A, ranks = dim K_n, d.d = 0 (Plan 75 / R10).

Koszul-gated on the FULL three-valued ext_algebra verdict: loud on False (naming the Ext
obstruction) and loud on None (honestly inconclusive) -- never on `g_quadratic` alone,
which would mislabel the genuinely-not-Koszul preprojective A_3 as merely "uncertified".
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.families import IncidenceAlgebra
from quiverlab.families.exterior import ExteriorAlgebra
from quiverlab.families.preprojective import PreprojectiveAlgebra
from quiverlab.fields import GF, QQ
from quiverlab.hochschild import koszul_ghms as G
from quiverlab.hochschild.koszul_ghms import GHMSResolution

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _kz3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}).algebra(
        relations=["a*b", "b*c", "c*a"], field=QQ)


def _diamond():
    return IncidenceAlgebra([("0", "a"), ("0", "b"), ("a", "1"), ("b", "1")], field=QQ)


# ------------------------------------------------------------------ terms + d.d = 0
@selfcert
@pytest.mark.parametrize("dom", [GF(32003), QQ])
def test_ghms_terms_and_dd_zero(dom):
    res = GHMSResolution(ExteriorAlgebra(2, field=dom), top=5)
    assert [res.rank(n) for n in range(6)] == [1, 2, 3, 4, 5, 6]      # dim K_n
    assert res.assert_dd_zero(5) is True


@selfcert
@pytest.mark.parametrize("build,ranks", [
    (_diamond, [4, 4, 1, 0, 0]),          # gl.dim 2 -- the resolution TERMINATES
    (_kz3, [3, 3, 3, 3, 3]),              # self-injective -- PERIODIC, never terminates
])
def test_ghms_multi_vertex(build, ranks):
    """Corner-typed multi-vertex terms: K_0 is indexed by VERTICES, so rank(0) = |Q_0|."""
    A = build()
    res = GHMSResolution(A, top=4)
    assert [res.rank(n) for n in range(5)] == ranks
    assert res.rank(0) == len(A.quiver.vertices)
    assert res.assert_dd_zero(4) is True
    # every generator carries a corner (v, w)
    for n in range(3):
        assert len(res.term(n)) == res.rank(n)
        assert all(isinstance(c, tuple) and len(c) == 2 for c in res.term(n))


# ---------------------------------------------- what d.d = 0 DOES and DOES NOT pin
@selfcert
@pytest.mark.parametrize("mode", ["constant_minus", "constant_plus"])
def test_dd_zero_catches_a_NON_ALTERNATING_sign(mode):
    """NON-VACUITY of the self-cert. Writing d_n = L_n + s_n R_n, and using L.L = R.R = 0
    plus the comultiplicative identity, composition leaves (s_n + s_{n-1})*L_{n-1}R_n. So
    d.d = 0 holds IFF the signs ALTERNATE. Replacing the shipped -(-1)^n by a CONSTANT
    sign must therefore be caught -- if it were not, `assert_dd_zero` would be decoration.
    """
    orig = G.GHMSResolution.differential
    try:
        G.GHMSResolution.differential = _patched(orig, mode)
        with pytest.raises(QuiverlabError, match=r"d_\d+ . d_\d+ != 0|sign"):
            GHMSResolution(ExteriorAlgebra(2, field=QQ), top=4).assert_dd_zero(4)
    finally:
        G.GHMSResolution.differential = orig


@selfcert
def test_dd_zero_is_BLIND_to_a_global_sign_flip():
    """THE HONEST LIMIT, recorded as a test rather than left implicit.

    Flipping the sign of the right-splitting term in EVERY degree keeps the signs
    alternating, so d.d = 0 still holds -- the self-cert cannot see it. That is not a
    defect: a global flip is the generator rescaling w_n |-> (-1)^n w_n, an ISOMORPHISM of
    complexes, so it cannot change homology either, and the cross-engine HH anchors cannot
    pin it any more than this can. The overall sign is a CONVENTION (chosen to match the
    published GHMS formula); what is genuinely arbitrated is the ALTERNATION, above.
    """
    orig = G.GHMSResolution.differential
    try:
        G.GHMSResolution.differential = _patched(orig, "global_flip")
        res = GHMSResolution(ExteriorAlgebra(2, field=QQ), top=4)
        assert res.assert_dd_zero(4) is True          # undetected, and provably harmless
    finally:
        G.GHMSResolution.differential = orig


def _patched(orig, mode):
    def f(self, n):
        D = orig(self, n)
        one = list(self.algebra.unit)
        shipped = -1 if (n % 2 == 0) else 1           # the shipped -(-1)^n
        want = {"global_flip": -shipped, "constant_minus": -1, "constant_plus": 1}[mode]
        r = want // shipped
        return [[[(c * r if u == one else c, u, v) for (c, u, v) in cell] for cell in row]
                for row in D]
    return f


# ------------------------------------------------------------------ canonicality
@selfcert
def test_ghms_canonical_differential():
    """Byte-reproducible: the splittings are READ OFF the path coordinates and the only
    solve is against an independent basis, so the differential is unique by construction --
    there is nothing to canonicalize away (contrast the CS correction solve, Plan 17)."""
    a = GHMSResolution(ExteriorAlgebra(3, field=GF(2)), top=4)
    b = GHMSResolution(ExteriorAlgebra(3, field=GF(2)), top=4)
    assert [a.differential(n) for n in range(1, 5)] == [b.differential(n) for n in range(1, 5)]


# ------------------------------------------------------------------ the Koszul gate
@selfcert
def test_non_quadratic_refuses():
    from quiverlab.families.basic import truncated_polynomial
    with pytest.raises(QuiverlabError, match="[Kk]oszul|quadratic|QUADRATIC"):
        GHMSResolution(truncated_polynomial(3, field=QQ))


@lit
def test_not_koszul_refuses_naming_the_obstruction():
    """Preprojective A_3 is NOT Koszul (Brenner-Butler-King almost-Koszul breaks at
    h - 1 = 3): ext_algebra reports a NEW GENERATOR in degree 3. GHMS must refuse, and the
    refusal must name that -- not merely say "uncertified"."""
    with pytest.raises(QuiverlabError, match="NOT Koszul|obstruction"):
        GHMSResolution(PreprojectiveAlgebra("A3", field=QQ))
    E = PreprojectiveAlgebra("A3", field=QQ).ext_algebra(3)
    assert E.koszul is False and E.koszul_obstruction[0] == 3


@selfcert
def test_gate_reports_why_it_accepted():
    res = GHMSResolution(ExteriorAlgebra(2, field=QQ), top=3)
    assert "Koszul" in res.koszul_reason
