"""Coxeter *spectral* literature-oracle battery (Plan 58 / R20).

The record's oracle set plus the m-Kronecker ladder and the ADE / affine cyclotomic
sweep, all as EXACT certified values (Φ_n labels, certified algebraic ρ/M, finite
Coxeter order Φ^m = I). Every pin is a fetched/recomputed literature value; a mismatch
means fix the quiver / transcription, never the pin (P38 rule).

Sources (registry keys owned by the citation registry; cited here as prose):
  * de la Peña, "On the Mahler measure of the Coxeter polynomial of an algebra",
    Adv. Math. 2014, arXiv:1310.1910 (Lehmer / E10).       [dlPena2014mahler]
  * de la Peña, "Algebras whose Coxeter polynomials are products of cyclotomic
    polynomials", 2013, arXiv:1310.1557 (periodic ⊊ cyclotomic-type). [dlPena2013cyclotomic]
  * Lenzing-de la Peña, "Spectral analysis...", ICRA XII, arXiv:0805.1018 (ADE/affine
    tables).                                               [lenzing_delapena_spectral]
"""
import sympy as sp
import pytest

from quiverlab import CC, GF, Quiver, linear_path_algebra
from quiverlab.families import dynkin_quiver

pytestmark = [pytest.mark.oracle_literature]

t = sp.Symbol("t")


# --------------------------------------------------------------------------- #
# builders (mirrors tests/invariants/test_coxeter_literature.py)
# --------------------------------------------------------------------------- #
def _kron(m, field=None):
    """The m-Kronecker: m parallel arrows 1 -> 2 (hereditary)."""
    field = field if field is not None else CC
    return Quiver([1, 2], {f"a{i}": (1, 2) for i in range(m)}).algebra(field=field)


def _star(arm_edges, field=CC):
    """The star tree: center (vertex 0) with arms of the given edge-counts, all
    oriented toward the center. [p1,..,pt] in de la Peña's weight notation has arm i
    of p_i - 1 edges, e.g. the wild star [2,3,7] = _star([1, 2, 6]) (10 vertices)."""
    verts, arrows, nxt = [0], {}, 1
    for ai, L in enumerate(arm_edges):
        prev = 0
        for e in range(L):
            verts.append(nxt)
            arrows[f"s{ai}_{e}"] = (nxt, prev)
            prev, nxt = nxt, nxt + 1
    return Quiver(sorted(verts), arrows).algebra(relations=[], field=field)


def _canonical(*weights, field=CC):
    """Canonical algebra C(p_1,..,p_t): a source (0) and a sink joined by t arms,
    arm i a directed path of p_i arrows. Two arms (t = 2) is the tame hereditary
    affine algebra A~_{p_1,p_2}."""
    sink = 1 + sum(w - 1 for w in weights)
    arrows, arm_full, nxt = {}, [], 1
    for ai, p in enumerate(weights):
        seq = [0]
        for _ in range(p - 1):
            seq.append(nxt)
            nxt += 1
        seq.append(sink)
        names = []
        for i in range(p):
            nm = f"a{ai}_{i}"
            arrows[nm] = (seq[i], seq[i + 1])
            names.append(nm)
        arm_full.append(names)
    rels = []
    if len(weights) >= 3:
        a1, a2 = "*".join(arm_full[0]), "*".join(arm_full[1])
        for i in range(2, len(weights)):
            rels.append(f"{a1} - {a2} + {'*'.join(arm_full[i])}")
    return Quiver(list(range(sink + 1)), arrows).algebra(relations=rels, field=field)


def _coxeter_number(A, cap=200):
    """min m >= 1 with Phi^m = I (exact matrix powers), or None -- the independent
    cross-check of coxeter_spectral's lcm-based coxeter_order."""
    Phi = sp.Matrix(A.coxeter_matrix())
    ident = sp.eye(Phi.rows)
    power = Phi.copy()
    for m in range(1, cap + 1):
        if power == ident:
            return m
        power = power * Phi
    return None


_LEHMER = t**10 + t**9 - t**7 - t**6 - t**5 - t**4 - t**3 + t + 1


# --------------------------------------------------------------------------- #
# 1. kA_n: v_{n+1} = product of cyclotomics, order = Coxeter number n+1
# --------------------------------------------------------------------------- #
def test_ka2_phi3():
    r = linear_path_algebra(2, field=CC).coxeter_spectral()
    assert [f["cyclotomic_index"] for f in r["factorization"]] == [3]
    assert r["cyclotomic"] and r["coxeter_order"] == 3
    assert r["outside_unit_circle_count"] == 0
    assert r["spectral_radius"]["value"] == "1" and r["mahler_measure"]["value"] == "1"


def test_ka4_single_phi5():
    """A_4 : v_5 = Phi_5 (5 prime => a single cyclotomic factor); order 5."""
    r = linear_path_algebra(4, field=CC).coxeter_spectral()
    assert [f["cyclotomic_index"] for f in r["factorization"]] == [5]
    assert r["coxeter_order"] == 5 and r["outside_unit_circle_count"] == 0


def test_ka5_three_factors():
    """A_5 : v_6 = Phi_2 * Phi_3 * Phi_6 (three cyclotomic factors); order lcm = 6."""
    r = linear_path_algebra(5, field=CC).coxeter_spectral()
    idx = sorted(f["cyclotomic_index"] for f in r["factorization"])
    assert idx == [2, 3, 6]
    assert all(f["multiplicity"] == 1 for f in r["factorization"])
    assert r["cyclotomic"] and r["coxeter_order"] == 6


# --------------------------------------------------------------------------- #
# 2. 3-Kronecker: the record's wild anchor
# --------------------------------------------------------------------------- #
def test_three_kronecker_record():
    r = _kron(3).coxeter_spectral()
    assert not r["cyclotomic"] and r["quasi_unipotent"] is False
    assert r["coxeter_order"] is None
    assert [f["cyclotomic_index"] for f in r["factorization"]] == [None]
    assert r["outside_unit_circle_count"] == 1
    c = r["spectral_radius"]
    assert c["minpoly"] == [1, -7, 1] and c["degree"] == 2
    assert (sp.Rational(c["interval"][0]), sp.Rational(c["interval"][1])) == (6, 7)
    assert sp.simplify(sp.sympify(c["value"]) - (7 + 3*sp.sqrt(5))/2) == 0
    # rho = M for the reciprocal Coxeter polynomial with a single outside root
    assert sp.simplify(sp.sympify(r["mahler_measure"]["value"]) - sp.sympify(c["value"])) == 0


# --------------------------------------------------------------------------- #
# 3. m-Kronecker ladder m in {1,2,3,4,5}
# --------------------------------------------------------------------------- #
def test_m_kronecker_ladder():
    # m = 1: A_2 = Phi_3 (cyclotomic, order 3)
    r1 = _kron(1).coxeter_spectral()
    assert [f["cyclotomic_index"] for f in r1["factorization"]] == [3]
    assert r1["coxeter_order"] == 3
    # m = 2: (t-1)^2 = Phi_1^2 (cyclotomic type, order None -- Jordan/affine A~_1)
    r2 = _kron(2).coxeter_spectral()
    assert [f["cyclotomic_index"] for f in r2["factorization"]] == [1]
    assert r2["factorization"][0]["multiplicity"] == 2
    assert r2["cyclotomic"] and r2["coxeter_order"] is None
    # m >= 3: wild x^2 - (m^2-2)x + 1, one root outside, rho = M = (a + sqrt(a^2-4))/2
    for m in (3, 4, 5):
        a = m * m - 2
        r = _kron(m).coxeter_spectral()
        assert not r["cyclotomic"] and r["coxeter_order"] is None
        assert r["outside_unit_circle_count"] == 1
        assert r["spectral_radius"]["minpoly"] == [1, -a, 1]
        expected = (a + sp.sqrt(a*a - 4)) / 2
        assert sp.simplify(sp.sympify(r["spectral_radius"]["value"]) - expected) == 0
        assert sp.simplify(sp.sympify(r["mahler_measure"]["value"]) - expected) == 0


# --------------------------------------------------------------------------- #
# 4. T_{2,3,7} = E10 = Lehmer  --  de la Peña arXiv:1310.1910
# --------------------------------------------------------------------------- #
def test_lehmer_e10_star():
    """The wild hereditary star [2,3,7] (10 vertices) has Coxeter polynomial =
    Lehmer's polynomial; rho = M is the smallest known Salem number (minpoly =
    Lehmer, interval (1,2)), one root outside the circle."""
    A = _star([1, 2, 6])
    assert len(A.quiver.vertices) == 10
    r = A.coxeter_spectral()
    assert not r["cyclotomic"] and r["coxeter_order"] is None
    assert [f["cyclotomic_index"] for f in r["factorization"]] == [None]
    assert r["outside_unit_circle_count"] == 1
    c = r["spectral_radius"]
    assert sp.Poly(c["minpoly"], sp.Symbol("x")) == sp.Poly(_LEHMER.subs(t, sp.Symbol("x")),
                                                            sp.Symbol("x"))
    assert (sp.Rational(c["interval"][0]), sp.Rational(c["interval"][1])) == (1, 2)
    assert sp.simplify(sp.sympify(r["mahler_measure"]["value"]) - sp.sympify(c["value"])) == 0


# --------------------------------------------------------------------------- #
# 5. Dynkin / affine cyclotomic sweep  --  Lenzing-de la Peña arXiv:0805.1018
# --------------------------------------------------------------------------- #
def test_dynkin_finite_orders():
    """A4/D5/E6/E7/E8 : all cyclotomic, finite Coxeter order = Coxeter number,
    rho = M = 1; the order agrees with the independent brute-force _coxeter_number."""
    stated = {
        linear_path_algebra(4, field=CC): 5,               # A4 : h = 5
        dynkin_quiver("D5").algebra(relations=[]): 8,       # D5 : h = 8
        dynkin_quiver("E6").algebra(relations=[]): 12,      # E6 : h = 12
        dynkin_quiver("E7").algebra(relations=[]): 18,      # E7 : h = 18
        dynkin_quiver("E8").algebra(relations=[]): 30,      # E8 : h = 30
    }
    for A, h in stated.items():
        r = A.coxeter_spectral()
        assert r["cyclotomic"] and r["quasi_unipotent"]
        assert r["coxeter_order"] == h == _coxeter_number(A)
        assert r["outside_unit_circle_count"] == 0
        assert r["spectral_radius"]["value"] == "1" and r["mahler_measure"]["value"] == "1"


def test_affine_cyclotomic_but_infinite_order():
    """Affine A~_{1,3} / D~4 / E~6 : all cyclotomic == True but coxeter_order is None
    (defect / Jordan block), rho = M = 1 -- the periodic ⊊ cyclotomic-type separation
    (de la Peña arXiv:1310.1557)."""
    for A in (_canonical(1, 3),           # A~_{1,3}
              _star([1, 1, 1, 1]),        # D~4
              _star([2, 2, 2])):          # E~6
        r = A.coxeter_spectral()
        assert r["cyclotomic"] and r["quasi_unipotent"]
        assert r["coxeter_order"] is None and _coxeter_number(A) is None
        assert r["outside_unit_circle_count"] == 0
        assert r["spectral_radius"]["value"] == "1" and r["mahler_measure"]["value"] == "1"


# --------------------------------------------------------------------------- #
# 6. Lehmer-class note is honest documentation only
# --------------------------------------------------------------------------- #
def test_lehmer_note_is_documentation_only():
    from quiverlab.invariants.coxeter_spectral import coxeter_spectral_block
    block = coxeter_spectral_block(_kron(3))
    note = block["lehmer_class_note"]
    assert "Lehmer" in note and "μ" in note
    # references mu_0 by its minimal polynomial (Lehmer's degree-10 polynomial)
    assert "t^10" in note
    # NEVER a general verdict on the open Lehmer problem
    assert "open" in note.lower()
    # mu_0 is given by its minimal polynomial, NEVER as a decimal float (the record's
    # 1.17628... is deliberately absent; the arXiv id 1310.1910 is not a float value)
    assert "1.17628" not in note and "1.176" not in note
