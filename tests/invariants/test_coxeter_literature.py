"""Coxeter / spectral literature-oracle battery (Plan 29 Part 1).

Every value below is an EXACT polynomial / integer pinned against a fetched,
recomputed source (see docs/plans/2026-07-25-literature-oracles-deep-research.md,
de la Pena / Lenzing / Marcos cluster). The convention is mapped and verified:
quiverlab's `coxeter_polynomial` = charpoly(-C^{-T} C) equals the papers'
charpoly(-C^{-1} C^t) exactly (variable rename only) for every triangular Cartan,
which covers all examples here.

Sources (registry keys owned by a sibling agent -- cited here as prose only):
  * Nakayama polynomials + Coxeter numbers:
    Lenzing-Meltzer-Ruan, "Nakayama algebras and Fuchsian singularities",
    arXiv:2112.15587, Section 2.7 and Prop 6.1.   [lenzing_meltzer_ruan]
  * Dynkin / affine / canonical tables, B3 parametric, perfect-square (Prop 1.5):
    Lenzing-de la Pena, "Spectral analysis of finite dimensional algebras and
    singularities", ICRA XII, EMS Congr. Rep. 2008, arXiv:0805.1018.
    [lenzing_delapena_spectral]
  * Lehmer's polynomial + Mahler measure ordering:
    de la Pena, "On the Mahler measure of the Coxeter polynomial of an algebra",
    Adv. Math. 2014, arXiv:1310.1910, Section 2.5.   [delapena_mahler]

Fast bucket: the Coxeter number search (Phi^m = I) and the Lehmer spectral radius
are all sub-second (the degree-14 Coxeter-number-126 case runs in ~0.01 s).
"""
import math

import pytest
import sympy as sp

from quiverlab import CC, GF, Quiver, linear_path_algebra
from quiverlab.families import dynkin_quiver, TruncatedPathAlgebra
from quiverlab.invariants.spectral import spectral_radius, mahler_measure

pytestmark = [pytest.mark.oracle_literature]

t = sp.Symbol("t")


# --------------------------------------------------------------------------- #
# builders + small helpers
# --------------------------------------------------------------------------- #
def _v(n):
    """v_n = (x^n - 1)/(x - 1) = 1 + x + ... + x^{n-1}."""
    return sum(t**k for k in range(n))


def _cox(A):
    """The library's Coxeter polynomial as an expanded sympy expression in t."""
    return sp.expand(A.coxeter_polynomial().as_expr())


def _eq(p, q):
    """Exact polynomial equality on integer/rational coefficients."""
    return sp.expand(p - q) == 0


def _cyc(k):
    return sp.cyclotomic_poly(k, t)


def _nakayama(n, r):
    """N_n(r): the uniserial Nakayama algebra = linear equioriented A_n with
    rad^r = 0 (kill every path of length r).  n vertices, Loewy length r."""
    return TruncatedPathAlgebra(dynkin_quiver(f"A{n}"), r)


def _star(arm_edges, field=CC):
    """The star tree: a center (vertex 0) with arms of the given edge-counts, all
    arms oriented toward the center.  Hereditary => Coxeter polynomial is
    orientation-free.  [p1,..,pt] in de la Pena's weight notation has arm i of
    p_i - 1 edges, e.g. the wild star [2,3,7] = _star([1, 2, 6]) (10 vertices)."""
    verts, arrows, nxt = [0], {}, 1
    for ai, L in enumerate(arm_edges):
        prev = 0
        for e in range(L):
            verts.append(nxt)
            arrows[f"s{ai}_{e}"] = (nxt, prev)     # inner -> ... -> center
            prev, nxt = nxt, nxt + 1
    return Quiver(sorted(verts), arrows).algebra(relations=[], field=field)


def _canonical(*weights, field=CC):
    """Canonical algebra C(p_1,..,p_t): a source (0) and a sink joined by t arms,
    arm i a directed path of p_i arrows (0 -> ... -> sink).  For t >= 2 there are
    2 + sum(p_i - 1) vertices.  For t >= 3 the t - 2 canonical relations make the
    parallel full-arm paths alpha_i linearly dependent (here alpha_1 - alpha_2 +
    alpha_i = 0, i.e. lambda = 1, generic for t = 3).  Two arms (t = 2, no
    relation) is exactly the tame hereditary affine algebra A~_{p_1,p_2}."""
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
    """min m >= 1 with Phi^m = I (exact integer/rational matrix powers), or None."""
    Phi = sp.Matrix(A.coxeter_matrix())
    ident = sp.eye(Phi.rows)
    power = Phi.copy()
    for m in range(1, cap + 1):
        if power == ident:
            return m
        power = power * Phi
    return None


def _perfect_square(value):
    v = int(value)
    r = math.isqrt(v) if v >= 0 else -1
    return v >= 0 and r * r == v


# --------------------------------------------------------------------------- #
# 1. Nakayama N_n(r)  --  LMR arXiv:2112.15587 Section 2.7 + Prop 6.1
# --------------------------------------------------------------------------- #
# The six verbatim polynomials (Coxeter number in the comment), all recomputed
# under quiverlab's convention in the research doc and re-pinned live here.
_NAKAYAMA = {
    (17, 8): (t + 1) * (t**16 + t**8 + 1),                   # Coxeter number 24
    (16, 3): (t + 1) * (t**6 - t**3 + 1) * (t**9 + 1),       # 18
    (15, 6): (t + 1) * (t**8 + t**4 + 1) * (t**6 + 1),       # 12
    (15, 4): (t + 1) * (t**8 + t**4 + 1) * (t**6 + 1),       # 12 (= chi(15,6))
    (15, 5): (t + 1) * (t**4 + 1) * (t**5 + 1)**2,           # 40
    (14, 7): (t + 1) * (t**6 - t**3 + 1) * (t**7 + 1),       # 126
}


def test_nakayama_six_polynomials():
    """LMR Section 2.7: the six exact Nakayama Coxeter polynomials, incl. the
    chi(15,6) = chi(15,4) coincidence and the (t^5+1)^2 repeated-factor case."""
    for (n, r), expected in _NAKAYAMA.items():
        chi = _cox(_nakayama(n, r))
        assert sp.Poly(chi, t).degree() == n            # degree == #vertices
        assert _eq(chi, expected), f"N_{n}({r}) mismatch: {chi}"
    # the stated coincidence, at the polynomial level
    assert _eq(_cox(_nakayama(15, 6)), _cox(_nakayama(15, 4)))


def test_nakayama_prop61_family():
    """LMR Prop 6.1: r >= 9  =>  chi(r+7, r) = (t+1)(t^6 - t^3 + 1)(t^r + 1).
    Two instances r = 9, 10 (n = 16, 17)."""
    for r in (9, 10):
        n = r + 7
        expected = (t + 1) * (t**6 - t**3 + 1) * (t**r + 1)
        assert _eq(_cox(_nakayama(n, r)), expected)
    # r = 9 lands exactly on chi(16,3) (both Coxeter number 18) -- a consistency
    # check between the Section-2.7 table and the Prop-6.1 family.
    assert _eq(_cox(_nakayama(16, 9)), _cox(_nakayama(16, 3)))


def test_nakayama_coxeter_numbers():
    """LMR: Coxeter number = min m with Phi^m = I, for the five stated cases.
    The degree-14 case (126) is included and remains sub-second."""
    stated = {(17, 8): 24, (16, 3): 18, (15, 6): 12, (15, 5): 40, (14, 7): 126}
    for (n, r), h in stated.items():
        assert _coxeter_number(_nakayama(n, r)) == h, f"N_{n}({r}) Coxeter number"


# --------------------------------------------------------------------------- #
# 2. Dynkin table  --  Lenzing-de la Pena arXiv:0805.1018
# --------------------------------------------------------------------------- #
def test_dynkin_A():
    """A_n : chi = v_{n+1} (Coxeter number n+1)."""
    for n in range(2, 7):
        assert _eq(_cox(linear_path_algebra(n)), _v(n + 1))


def test_dynkin_D_both_forms():
    """D_n : chi = Phi_2 * (x^{n-1} + 1) -- SHIP the v-factorization, not the
    mis-transcribed cyclotomic-condition column (research-doc caveat).  Assert it
    agrees with the cyclotomic products D4 = Phi_2^2 Phi_6, D5 = Phi_2 Phi_8,
    D6 = Phi_2^2 Phi_10, and that both forms coincide."""
    cyc = {4: _cyc(2)**2 * _cyc(6), 5: _cyc(2) * _cyc(8), 6: _cyc(2)**2 * _cyc(10)}
    for n in (4, 5, 6):
        chi = _cox(dynkin_quiver(f"D{n}").algebra(relations=[]))
        vform = _cyc(2) * (t**(n - 1) + 1)
        assert _eq(cyc[n], vform)          # the two literature forms agree
        assert _eq(chi, vform)             # and the library reproduces them


def test_dynkin_E():
    """E6 = Phi_3 Phi_12, E7 = Phi_2 Phi_18, E8 = Phi_30."""
    forms = {6: _cyc(3) * _cyc(12), 7: _cyc(2) * _cyc(18), 8: _cyc(30)}
    for n, form in forms.items():
        assert _eq(_cox(dynkin_quiver(f"E{n}").algebra(relations=[])), form)


def test_dynkin_orientation_independence():
    """The Coxeter polynomial of a Dynkin path algebra is orientation-independent
    (derived invariance).  Two genuinely different orientations of A4 and of D4 --
    including a non-monotone zigzag whose Cartan matrix differs -- agree."""
    a_lin = linear_path_algebra(4)
    a_zig = dynkin_quiver("A4", {"e12": (1, 2), "e23": (3, 2), "e34": (3, 4)}).algebra(relations=[])
    a_rev = dynkin_quiver("A4", "reverse").algebra(relations=[])
    assert a_lin.cartan_matrix() != a_zig.cartan_matrix()      # different orientation
    assert _eq(_cox(a_lin), _cox(a_zig))
    assert _eq(_cox(a_lin), _cox(a_rev))

    d_lin = dynkin_quiver("D4", "linear").algebra(relations=[])
    d_out = dynkin_quiver("D4", {"e12": (2, 1), "e23": (2, 3), "e24": (2, 4)}).algebra(relations=[])
    assert _eq(_cox(d_lin), _cox(d_out))


# --------------------------------------------------------------------------- #
# 3. Affine + canonical  --  Lenzing-de la Pena arXiv:0805.1018
# --------------------------------------------------------------------------- #
def test_affine_A_pq_and_orientation_dependence():
    """A~_{p,q} = C(p,q) (two hereditary arms) : chi = (x-1)^2 v_p v_q.  This is
    the ONLY orientation-dependent Euclidean case.  (1,2) and (2,1) agree by
    symmetry; on the 4-cycle (1,3) and (2,2) are different acyclic orientations
    with different Coxeter polynomials."""
    for (p, q) in [(1, 2), (2, 1), (1, 3), (2, 2)]:
        assert _eq(_cox(_canonical(p, q)), (t - 1)**2 * _v(p) * _v(q))
    assert _eq(_cox(_canonical(1, 2)), _cox(_canonical(2, 1)))          # symmetry
    assert not _eq(_cox(_canonical(1, 3)), _cox(_canonical(2, 2)))      # orientation-dependent


def test_affine_D4_E6():
    """D~4 = (x-1)^2 v_2^3  (star with four 1-edge arms, 5 vertices);
    E~6 = (x-1)^2 v_2 v_3^2 (star with three 2-edge arms, 7 vertices)."""
    assert _eq(_cox(_star([1, 1, 1, 1])), (t - 1)**2 * _v(2)**3)
    assert _eq(_cox(_star([2, 2, 2])), (t - 1)**2 * _v(2) * _v(3)**2)


def test_canonical_algebras():
    """Canonical algebras C(p_1,..,p_t), t = 3 : chi = (x-1)^2 prod v_{p_i}
    (depends only on the weights, not the parameter lambda).  Built as a
    quiver-with-relations over both CC and GF(32003); 2 + sum(p_i - 1) vertices."""
    cases = {
        (2, 2, 2): (t - 1)**2 * _v(2)**3,
        (3, 3, 3): (t - 1)**2 * _v(3)**3,
        (2, 3, 5): (t - 1)**2 * _v(2) * _v(3) * _v(5),
    }
    for weights, expected in cases.items():
        for field in (CC, GF(32003)):
            A = _canonical(*weights, field=field)
            assert len(A.quiver.vertices) == 2 + sum(p - 1 for p in weights)
            assert _eq(_cox(A), expected), f"C{weights} over {field}"


# --------------------------------------------------------------------------- #
# 4. Lehmer  --  de la Pena arXiv:1310.1910 Section 2.5
# --------------------------------------------------------------------------- #
_LEHMER = t**10 + t**9 - t**7 - t**6 - t**5 - t**4 - t**3 + t + 1


def test_lehmer_polynomial():
    """The wild hereditary star [2,3,7] (10 vertices) has Coxeter polynomial =
    Lehmer's polynomial, and its Mahler measure equals its spectral radius (the
    smallest known Salem number), EXACTLY."""
    chi = _cox(_star([1, 2, 6]))
    assert _eq(chi, _LEHMER)
    rho = spectral_radius(chi)
    mm = mahler_measure(chi)
    assert sp.simplify(rho - mm) == 0            # exact algebraic equality
    # rho is the Salem root: its minimal polynomial is Lehmer (irreducible monic)
    # and it exceeds 1 -- both exact.
    assert sp.minimal_polynomial(rho, t) == _LEHMER
    assert (rho - 1).is_positive


def test_lehmer_spectral_ordering():
    """de la Pena Section 2.5: strict ordering rho[2,4,5] > rho[2,3,8] > rho[2,3,7],
    decided by exact sympy sign of the differences (no float thresholds)."""
    r245 = spectral_radius(_cox(_star([1, 3, 4])))     # [2,4,5]
    r238 = spectral_radius(_cox(_star([1, 2, 7])))     # [2,3,8]
    r237 = spectral_radius(_cox(_star([1, 2, 6])))     # [2,3,7] = Lehmer
    assert (r245 - r238).is_positive
    assert (r238 - r237).is_positive


def test_mahler_measure_one_iff_cyclotomic_type():
    """Mahler measure == 1 exactly for the cyclotomic-type (Dynkin / affine)
    members -- the clean cyclotomic-type predicate."""
    cyclotomic = [
        linear_path_algebra(4),
        dynkin_quiver("D4").algebra(relations=[]),
        dynkin_quiver("E6").algebra(relations=[]),
        _canonical(1, 3),                # affine A~_{1,3}
        _star([1, 1, 1, 1]),             # affine D~4
        _star([2, 2, 2]),                # affine E~6
        _canonical(2, 2, 2),             # canonical C(2,2,2) (cyclotomic)
    ]
    for A in cyclotomic:
        assert mahler_measure(_cox(A)) == sp.Integer(1)
    # contrast: Lehmer's wild star is NOT cyclotomic type
    assert mahler_measure(_LEHMER) != sp.Integer(1)


# --------------------------------------------------------------------------- #
# 5. B3 parametric  --  Lenzing-de la Pena arXiv:0805.1018
# --------------------------------------------------------------------------- #
def _cox_from_cartan(C):
    """Matrix-level Coxeter polynomial, quiverlab's convention Phi = -C^{-T} C.
    Used for the general (a,b,c) triangular Cartan, which no single hereditary
    kQ realizes independently -- the SOURCE statement is itself matrix-level."""
    C = sp.Matrix(C)
    Phi = -C.inv().T * C
    return sp.expand(Phi.charpoly(t).as_expr())


def test_b3_parametric_formula():
    """Cartan [[1,a,b],[0,1,c],[0,0,1]] => chi = x^3 + alpha x^2 + alpha x + 1
    with alpha = abc - a^2 - b^2 - c^2 + 3 (matrix-level; the paper's statement is
    itself over arbitrary triangular 3x3 Cartans)."""
    for (a, b, c) in [(1, 1, 1), (2, 1, 1), (1, 2, 3), (2, 3, 4)]:
        alpha = a * b * c - a**2 - b**2 - c**2 + 3
        expected = t**3 + alpha * t**2 + alpha * t + 1
        assert _eq(_cox_from_cartan([[1, a, b], [0, 1, c], [0, 0, 1]]), expected)


def test_b3_parametric_a3_instance():
    """(a,b,c) = (1,1,1) is realized by the hereditary linear A3 whose Cartan is
    exactly [[1,1,1],[0,1,1],[0,0,1]]: the library's coxeter_polynomial agrees
    with the matrix-level formula and with A3 = v_4 (alpha = 1)."""
    A3 = linear_path_algebra(3)
    assert A3.cartan_matrix() == [[1, 1, 1], [0, 1, 1], [0, 0, 1]]
    assert _eq(_cox(A3), _cox_from_cartan([[1, 1, 1], [0, 1, 1], [0, 0, 1]]))
    assert _eq(_cox(A3), t**3 + t**2 + t + 1)      # alpha = 1*1*1 - 3 + 3 = 1
    assert _eq(_cox(A3), _v(4))


# --------------------------------------------------------------------------- #
# 6. Happel / Lenzing-de la Pena Prop 1.5 : chi_A(-1) is a perfect square
# --------------------------------------------------------------------------- #
def test_perfect_square_sweep():
    """chi_A(-1) is a perfect square for every (triangular) battery member
    (Happel; Lenzing-de la Pena Prop 1.5)."""
    members = [
        linear_path_algebra(4), linear_path_algebra(6),
        dynkin_quiver("D5").algebra(relations=[]),
        dynkin_quiver("E7").algebra(relations=[]),
        dynkin_quiver("E8").algebra(relations=[]),
        _nakayama(17, 8), _nakayama(15, 5), _nakayama(14, 7),
        _canonical(1, 2), _canonical(1, 3), _canonical(2, 2),
        _star([1, 1, 1, 1]), _star([2, 2, 2]),
        _canonical(2, 2, 2), _canonical(3, 3, 3), _canonical(2, 3, 5),
        _star([1, 2, 6]),                          # Lehmer star
    ]
    for A in members:
        value = A.coxeter_polynomial().as_expr().subs(t, -1)
        assert _perfect_square(value), f"chi(-1) = {value} is not a perfect square"
