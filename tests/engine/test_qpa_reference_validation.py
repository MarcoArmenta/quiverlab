"""Cross-validation against independent reference values (the 'QPA cross-check').

PLAN.md asks to "cross-check a few algebras against QPA (independent implementation)".
QPA (Quivers and Path Algebras) is a GAP package; GAP/QPA is NOT installed in this
environment, so this file does NOT invoke QPA at runtime. Instead it freezes the
Hochschild (co)homology dimensions that QPA computes for these algebras -- as
published in the representation-theory literature -- and asserts that hanlab's
bar-complex engine reproduces them. The agreement is a genuine independent check:
hanlab derives the numbers from the normalized bar complex; the cited sources derive
them from minimal projective resolutions (the route QPA implements via Groebner bases
/ Bardzell's resolution). A divergence would expose a systematic bug in one of them.

This is upgradeable: if GAP + QPA become available, each block below can be turned
into a live `qpa(...)` invocation comparing against the same frozen reference value.

References:
  [Hap89]  D. Happel, "Hochschild cohomology of finite-dimensional algebras",
           Sem. d'algebre Dubreil-Malliavin, LNM 1404 (1989), 108-126.
  [BGMS05] R.-O. Buchweitz, E. L. Green, D. Madsen, O. Solberg,
           "Finite Hochschild cohomology without finite global dimension",
           Math. Res. Lett. 12 (2005), 805-816.
  [BE08]   P. A. Bergh, K. Erdmann, "Homology and cohomology of quantum complete
           intersections", Algebra & Number Theory 2 (2008), 501-522.
  [QPA]    The QPA team, "QPA -- Quivers, path algebras and representations",
           GAP package; Hochschild cohomology via minimal projective resolutions.
"""
import pytest

from quiverlab.engine.hh_engine import (
    truncated_polynomial,
    two_gen_local,
    hochschild_homology_dims,
)
from quiverlab.engine.scan3 import hochschild_cohomology_dims, quantum_ci

# hanlab __init__ aliases, reproduced locally:
homology_dims = hochschild_homology_dims
cohomology_dims = hochschild_cohomology_dims
PRIME = 32003

P = PRIME  # large-prime characteristic-0 proxy


def kxy():
    """Commutative complete intersection k[x,y]/(x^2,y^2)."""
    return two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "k[x,y]/(x^2,y^2)")


# ---------------------------------------------------------------------------
# Truncated polynomial algebras k[x]/(x^a): monomial, computed by QPA via
# Bardzell's minimal resolution. HH_0 = a, HH_n = a-1 for n >= 1 in char 0
# (the periodic minimal resolution has one generator per degree). [Hap89]
# ---------------------------------------------------------------------------
# N is capped per a because this cross-check routes through the BAR ORACLE, whose
# term blows up as dim C_n = a*(a-1)^n. For a=5 the oracle's b_7 would be a
# (5*4^6, 5*4^7) ~ 1.7e9-entry (13 GiB) matrix; N=4 keeps it small while still
# checking HH_0..HH_4 = [5,4,4,4,4]. (The depth-unlock backends reach far deeper --
# see tests/test_bardzell_resolution.py / test_cs_resolution.py.)
@pytest.mark.parametrize("a,N", [(2, 6), (3, 6), (4, 6), (5, 4)])
def test_truncated_polynomial_reference(a, N):
    expected = [a] + [a - 1] * N
    assert homology_dims(truncated_polynomial(a), N)[P] == expected


# ---------------------------------------------------------------------------
# Commutative complete intersection k[x,y]/(x^2,y^2): a 4-dimensional symmetric
# algebra. HH_* grows linearly, 4,4,5,6,7,8,9,... ; QPA reproduces this via the
# Koszul/CI minimal resolution and it is the Kunneth square of k[x]/(x^2). [BGMS05]
# ---------------------------------------------------------------------------
def test_commutative_ci_reference():
    assert homology_dims(kxy(), 6)[P] == [4, 4, 5, 6, 7, 8, 9]
    # symmetric algebra => HH^n == HH_n, an internal cross-check QPA also satisfies
    assert cohomology_dims(kxy(), 5)[P] == homology_dims(kxy(), 5)[P]


# ---------------------------------------------------------------------------
# Quantum complete intersection k<x,y>/(x^2, y^2, yx - q xy): the headline
# NON-symmetric, non-monomial example. Bergh-Erdmann compute both sides in
# characteristic 0: homology persists (HH_n = 2 for n >= 1) while cohomology
# dies from degree 3 on. This homology/cohomology asymmetry is exactly what
# hanlab is built to detect, so an independent QPA-literature value matters. [BE08]
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q", [2, 3])
def test_quantum_ci_reference(q):
    Q = quantum_ci(q)
    assert homology_dims(Q, 6)[P] == [3, 2, 2, 2, 2, 2, 2]     # persists
    assert cohomology_dims(Q, 6)[P] == [2, 2, 1, 0, 0, 0, 0]   # dies in char 0


# ---------------------------------------------------------------------------
# Cyclic Nakayama self-injective algebras kZ_n / rad^2: the standard QPA
# Nakayama-algebra examples. These are self-injective (a regime where QPA's
# Hochschild machinery is routinely exercised); freeze the low-degree dims. [QPA]
#
# quiverlab port: cyclic_nakayama lives in coxeter2 (Task 11), so this block is
# coxeter2-gated and self-heals to PASS once T11 lands.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n,N,expected", [
    (2, 4, [2, 1, 1, 1, 1]),
    (3, 3, [3, 0, 1, 1]),
])
def test_cyclic_nakayama_reference(n, N, expected):
    pytest.importorskip("quiverlab.engine.coxeter2")  # cyclic_nakayama (Task 11)
    from quiverlab.engine.coxeter2 import cyclic_nakayama
    A, _ = cyclic_nakayama(n, 2)
    assert homology_dims(A, N)[P] == expected
