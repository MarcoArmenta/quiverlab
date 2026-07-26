"""Plan 31 — the TrivialExtension double-quiver PRESENTATION battery.

``TrivialExtension(A)`` of a quiver-presented ``A`` now returns a genuine
``kQ_T/I_T``-presented ``Algebra`` (built via ``Quiver.algebra``), so the
path-basis invariants (``is_symmetric`` / ``is_frobenius`` /
``is_weakly_symmetric`` / ``is_selfinjective`` / ``cartan_matrix`` /
``loewy_length`` / ``center`` / ``engine="cs"``) all serve ``T(A)``. The
construction is certified PER INSTANCE by ``dim T(A) == 2 * dim A`` (design
decision D2) and QPA-oracled natively in
``tests/qpa/test_trivial_extension_qpa.py`` — no transcribed theorem is
load-bearing (Happel 1988, LMS LNS 119, is the classical prose reference).

Oracle classes pinned here (registry keys in parentheses):

* structural certificate — ``dim T(A) = 2*dim A`` on every base and
  ``_validate()`` (associativity + two-sided unit); the number of new dual
  arrows = ``dim soc_{A^e} A`` (bimodule socle), new arrows named ``te0, te1, …``.
* Nakayama oracle (``nakayama``, ``assem_book``) — ``T(kA_n) ≅ NakayamaAlgebra(n,
  n+1, cyclic=True)`` = the symmetric cyclic Nakayama / Brauer star
  ``kZ_n/J^{n+1}`` (``n | (L-1)`` with ``L = n+1``): equal dim, Cartan, Loewy
  length ``n+1``, and all four self-injectivity booleans.
* commutative-CI oracle (``assem_book``; ``bar``) — ``T(k[x]/(x^a)) ≅
  k⟨x,y⟩/(x^a, y^2, xy − yx)`` (the PLAIN commutator, every characteristic; D(A)
  is an honest bimodule with no Koszul sign): equal dim ``2a``, Cartan ``[[2a]]``,
  Loewy ``a+1``, symmetry; for ``a=2`` this is ``k[x,y]/(x²,y²)`` and reproduces
  the ``HH_• = [4,4,5,6]`` oracle (see the iso-invariance test).
* Cartan identity (``assem_book``) — ``C_{T(A)} = C_A + C_A^T`` elementwise (repo
  convention ``C[i][j] = dim e_i A e_j``), incl. the multi-vertex zoo record
  ``line_abc_cde`` (Plan 18).
* symmetry (``skowronski_yamagata``) — every ``T(A)`` is symmetric / weakly
  symmetric / Frobenius / self-injective over GF(32003) and QQ, and over GF(2)
  for the Nakayama and dual-number images (the trace-form positive branch is a
  field-agnostic full-rank Gram witness).
* iso-invariance (``bar``) — the presented build and the RETAINED ⋉
  structure-constant build (``_trivial_extension_structure_constants``, kept as
  the oracle) agree degreewise on HH^• and HH_•, and HH^n = HH_n (T(A) is
  symmetric). Over GF(32003) the auto engine is the fast bar-rank engine — the
  same normalized-bar complex as the pure oracle, identical dims.
* cross-engine (``chouhy_solotar``; ``bar``) — the CS engine, REFUSED on the old
  presentation-free build (``ValueError``: CS needs a quiver), now computes on the
  presented ``T(kA_2)`` and agrees with the bar oracle.
* fallback (``assem_book``) — a presentation-less base (raw structure constants)
  returns the UNCHANGED ⋉ build (``quiver is None``); ``is_symmetric`` /
  ``is_frobenius`` refuse loudly (``FieldError``) while ``dim`` and the bar HH
  oracle still serve (design decision D3).

``tests/families/`` auto-assigns to the deep CI bucket (``tests/conftest.py``).
"""
import pytest

from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.families import NakayamaAlgebra, TrivialExtension
from quiverlab.families import trivial_extension as _te_mod  # ⋉ oracle lives here
from quiverlab.families.zoo import build_from_record, load_catalog
from quiverlab.errors import FieldError
from quiverlab.fields import QQ

F = GF(32003)      # characteristic-0 proxy (large prime; fast int64 rank engine)
F2 = GF(2)


# ----------------------------------------------------------------- base builders
def _kA(n, field):
    return linear_path_algebra(n, field=field)


def _kron(field):
    """The 2-Kronecker quiver 1 ==> 2 (two parallel arrows), hereditary."""
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=field)


def _comm_square(field):
    """The commutative square with the relation a*c - b*d (dim 9)."""
    return Quiver([1, 2, 3, 4],
                  {"a": (1, 2), "b": (1, 3), "c": (2, 4), "d": (3, 4)}
                  ).algebra(relations=["a*c - b*d"], field=field)


def _dual(a, field):
    return truncated_polynomial(a, field=field)


_BASES = {
    "kA2": lambda fld: _kA(2, fld),
    "kA3": lambda fld: _kA(3, fld),
    "kA4": lambda fld: _kA(4, fld),
    "kron": lambda fld: _kron(fld),
    "comm": lambda fld: _comm_square(fld),
    "dual2": lambda fld: _dual(2, fld),
    "dual3": lambda fld: _dual(3, fld),
}

# frozen dims (research 2026-07-26, machine-verified): dim T(A) = 2 * dim A.
_DIM = {"kA2": 6, "kA3": 12, "kA4": 20, "kron": 8, "comm": 18, "dual2": 4, "dual3": 6}
# frozen new-arrow counts = dim of the bimodule socle of A.
_NEW_ARROWS = {"kA2": 1, "kA3": 1, "kron": 2, "dual3": 1, "comm": 1}

# module-level memo caches (Algebra objects are immutable value objects; reuse is
# safe and keeps the presented builds — the expensive step — off the hot path).
_BASE_CACHE = {}
_TE_CACHE = {}


def _base(name, field):
    key = (name, repr(field))
    A = _BASE_CACHE.get(key)
    if A is None:
        A = _BASES[name](field)
        _BASE_CACHE[key] = A
    return A


def _TE(name, field):
    """Presented ``TrivialExtension`` of the named base over ``field`` (memoized)."""
    key = (name, repr(field))
    T = _TE_CACHE.get(key)
    if T is None:
        T = TrivialExtension(_base(name, field))
        _TE_CACHE[key] = T
    return T


def _plus_transpose(C):
    n = len(C)
    return [[C[i][j] + C[j][i] for j in range(n)] for i in range(n)]


# =============================================================== structural / D2
@pytest.mark.parametrize("name,expected_dim", list(_DIM.items()))
def test_dim_is_twice_base_and_validates(name, expected_dim):
    """D2 certificate: ``dim T(A) = 2*dim A`` on every base; ``_validate()`` re-checks
    that the block multiplication stays associative with a two-sided unit."""
    A = _base(name, F)
    T = _TE(name, F)
    assert T.quiver is not None                      # presented route (not the ⋉ fallback)
    assert T.dim == 2 * A.dim == expected_dim
    T._validate()                                    # raises on any structure-constant defect


@pytest.mark.parametrize("name,expected_new", list(_NEW_ARROWS.items()))
def test_new_arrow_count_is_bimodule_socle_dim(name, expected_new):
    """#new dual arrows = ``dim soc_{A^e} A`` (frozen: kA_2→1, kA_3→1, 2-Kronecker→2,
    k[x]/(x^3)→1, comm-square→1). New arrows are named ``te0, te1, …`` (D4;
    underscore-prefixed only on a name collision, which none of these bases hit)."""
    A = _base(name, F)
    T = _TE(name, F)
    te_named = [nm for nm in T.quiver.arrows if nm not in A.quiver.arrows]
    assert len(T.quiver.arrows) - len(A.quiver.arrows) == expected_new
    assert len(te_named) == expected_new
    assert all(nm.lstrip("_").startswith("te") for nm in te_named)


# ================================================================ Nakayama oracle
@pytest.mark.parametrize("n", [2, 3, 4])
def test_TkAn_is_symmetric_cyclic_nakayama(n):
    """``T(kA_n) ≅ NakayamaAlgebra(n, n+1, cyclic=True)`` = the Brauer star
    ``kZ_n/J^{n+1}`` (``nakayama``; ``assem_book``; ``skowronski_yamagata``): equal
    dim, Cartan, Loewy length ``n+1``, and all four self-injectivity booleans."""
    T = _TE(f"kA{n}", F)
    N = NakayamaAlgebra(n=n, l=n + 1, cyclic=True, field=F)
    assert T.dim == N.dim == n * (n + 1)
    assert T.cartan_matrix() == N.cartan_matrix()
    assert T.loewy_length() == N.loewy_length() == n + 1
    for A in (T, N):
        assert A.is_symmetric() is True
        assert A.is_weakly_symmetric() is True
        assert A.is_frobenius() is True
        assert A.is_selfinjective() is True


# =========================================================== commutative-CI oracle
@pytest.mark.parametrize("a", [2, 3])
def test_Tdual_is_two_loop_commutator(a):
    """``T(k[x]/(x^a)) ≅ k⟨x,y⟩/(x^a, y^2, xy − yx)`` — the PLAIN commutator in every
    characteristic (``assem_book``; ``skowronski_yamagata``): equal dim ``2a``,
    Cartan ``[[2a]]``, Loewy ``a+1``, symmetry. For ``a=2`` this is
    ``k[x,y]/(x²,y²)`` (``HH_• = [4,4,5,6]`` pinned in the iso-invariance test)."""
    T = _TE(f"dual{a}", F)
    two_loop = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=[f"x^{a}", "y^2", "x*y - y*x"], field=F)
    assert T.dim == two_loop.dim == 2 * a
    assert T.cartan_matrix() == two_loop.cartan_matrix() == [[2 * a]]
    assert T.loewy_length() == two_loop.loewy_length() == a + 1
    for A in (T, two_loop):
        assert A.is_symmetric() is True
        assert A.is_weakly_symmetric() is True
        assert A.is_frobenius() is True
        assert A.is_selfinjective() is True


# =============================================================== Cartan identity
@pytest.mark.parametrize("name", ["kA2", "kA3", "kron", "comm"])
def test_cartan_identity(name):
    """``C_{T(A)} = C_A + C_A^T`` elementwise (``assem_book``); always symmetric.
    2-Kronecker gives ``[[2,2],[2,2]]``."""
    A = _base(name, F)
    T = _TE(name, F)
    assert T.cartan_matrix() == _plus_transpose(A.cartan_matrix())


def test_cartan_identity_zoo_line_abc_cde():
    """``C_{T(A)} = C_A + C_A^T`` on the multi-vertex zoo record ``line_abc_cde``
    (Plan 18 diversity record; base ``han_conjecture``/``chouhy_solotar``,
    identity ``assem_book``)."""
    rec = next(r for r in load_catalog() if r.get("name") == "line_abc_cde")
    A = build_from_record(rec, field=F)
    T = TrivialExtension(A)
    assert T.quiver is not None
    assert T.dim == 2 * A.dim
    assert T.cartan_matrix() == _plus_transpose(A.cartan_matrix())


# ========================================================================= symmetry
@pytest.mark.parametrize("field", [F, QQ])
@pytest.mark.parametrize("name", list(_BASES))
def test_every_TA_is_symmetric_over_gf_and_qq(name, field):
    """Every ``T(A)`` is symmetric / weakly symmetric / Frobenius / self-injective
    over GF(32003) and QQ (``skowronski_yamagata`` trace-form certifier;
    ``assem_book``)."""
    T = _TE(name, field)
    assert T.quiver is not None
    assert T.is_frobenius() is True
    assert T.is_symmetric() is True
    assert T.is_weakly_symmetric() is True
    assert T.is_selfinjective() is True


@pytest.mark.parametrize("name", ["kA2", "kA3", "kA4", "dual2", "dual3"])
def test_TA_is_symmetric_over_gf2(name):
    """GF(2) is safe on the Nakayama and dual-number images: the trace-form positive
    branch is a field-agnostic full-rank Gram witness (research finding 5). (Kronecker
    / comm-square are cut from GF(2) per the research cut-list, rank 1.)"""
    T = _TE(name, F2)
    assert T.is_frobenius() is True
    assert T.is_symmetric() is True
    assert T.is_weakly_symmetric() is True
    assert T.is_selfinjective() is True


# =================================================================== iso-invariance
# (name, top, HH^• == HH_• dims): the presented build and the retained ⋉ build agree
# degreewise AND HH^n == HH_n (T(A) symmetric). Windows kept shallow — the bar complex
# blows up on self-injectives; dim >= 12 is pinned at top 1 only.
_ISO = [
    ("dual2", 3, [4, 4, 5, 6]),
    ("kA2", 3, [3, 1, 1, 1]),
    ("kron", 2, [3, 4, 6]),
    ("kA3", 1, [4, 1]),
    ("kA4", 1, [5, 1]),
    ("comm", 1, [5, 1]),
]


@pytest.mark.parametrize("name,top,expected", _ISO)
def test_iso_invariance_presented_vs_lagrange_build(name, top, expected):
    """HH^• and HH_• of the presented ``T(A)`` equal those of the retained ⋉ build
    ``_trivial_extension_structure_constants(A)`` degreewise, and equal each other
    (HH^n == HH_n) — the iso-invariance oracle for the presentation (``bar``). Over
    GF(32003) the auto engine is the fast bar-rank engine (same normalized-bar
    complex; identical dims)."""
    A = _base(name, F)
    T_pres = _TE(name, F)
    T_lag = _te_mod._trivial_extension_structure_constants(A)
    assert T_lag.quiver is None                          # the retained structure-constant oracle
    co_pres = T_pres.hochschild_cohomology(top).dims
    co_lag = T_lag.hochschild_cohomology(top).dims
    ho_pres = T_pres.hochschild_homology(top).dims
    ho_lag = T_lag.hochschild_homology(top).dims
    assert co_pres == co_lag == expected                 # presented == ⋉ oracle, cohomology
    assert ho_pres == ho_lag == expected                 # presented == ⋉ oracle, homology (= HH^n)


def test_dualnum_matches_explicit_kxy():
    """The a=2 dual-number cross-link: ``T(k[x]/(x^2))`` (presented), its ⋉ build, AND
    the explicit ``k[x,y]/(x²,y²)`` presentation all give ``HH_• = [4,4,5,6]`` at
    top 3 (``bar``; the verification-page CI oracle)."""
    kxy = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x^2", "y^2", "x*y - y*x"], field=F)
    assert kxy.hochschild_homology(3).dims == [4, 4, 5, 6]
    assert _TE("dual2", F).hochschild_homology(3).dims == [4, 4, 5, 6]


# ===================================================================== cross-engine
def test_cs_now_computes_on_presented_TkA2():
    """The Chouhy-Solotar engine was REFUSED on the old presentation-free ``T(A)``
    (needs a quiver + relations). The presented ``T(kA_2)`` carries both, so
    ``engine="cs"`` now COMPUTES and agrees with the bar oracle: ``[3,1,1]`` at
    top 2 (``chouhy_solotar``; ``bar``)."""
    T = _TE("kA2", F)
    cs = T.hochschild_cohomology(2, engine="cs").dims
    bar = T.hochschild_cohomology(2).dims                # auto = fast bar-rank over GF(p)
    assert cs == bar == [3, 1, 1]


# ========================================================================= fallback
def test_presentationless_base_falls_back_to_lagrange_build():
    """D3 fallback: a base with NO quiver (raw structure constants) returns the
    UNCHANGED ⋉ build — ``quiver is None``, ``is_symmetric`` / ``is_frobenius`` refuse
    loudly (``FieldError``, exactly as on any presentation-less algebra), while ``dim``
    and the bar HH oracle still serve (``assem_book``)."""
    from quiverlab.core.algebra import Algebra
    # k[x]/(x^2) as raw structure constants, no quiver (mirrors test_symmetric_regression).
    sc = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(sc, [1, 0], field=CC)
    assert A.quiver is None
    T = TrivialExtension(A)
    assert T.quiver is None                              # fallback: ⋉ structure-constant build
    assert T.dim == 2 * A.dim == 4
    with pytest.raises(FieldError):
        T.is_symmetric()
    with pytest.raises(FieldError):
        T.is_frobenius()
    assert T.hochschild_cohomology(1, engine="bar").dims == [4, 4]   # bar still serves


def test_cc_algebraic_base_falls_back_per_d3():
    """D3 second fallback branch: a base that HAS a quiver but lives over CC
    (coefficients not string-representable for Quiver.algebra) also returns the ⋉
    build — the certificate never lets a wrong presented algebra through
    (``assem_book``). quiver is None; is_symmetric refuses loudly; dim and the bar
    HH oracle still serve. Distinct from the presentation-less trigger above."""
    A = linear_path_algebra(2, field=CC)                 # quiver present, CC coefficients
    assert A.quiver is not None
    T = TrivialExtension(A)
    assert T.quiver is None                              # CC-algebraic -> ⋉ fallback (D3)
    assert T.dim == 2 * A.dim == 6
    with pytest.raises(FieldError):
        T.is_symmetric()
    assert T.hochschild_cohomology(1, engine="bar").dims[1] >= 1     # bar still serves


# ================================================================== Loewy / center
@pytest.mark.parametrize("name,expected_loewy", [
    ("kA2", 3), ("kA3", 4), ("kA4", 5),      # T(kA_n): Loewy = n+1
    ("dual2", 3), ("dual3", 4),              # T(k[x]/x^a): Loewy = a+1
    ("kron", 3),                             # T(2-Kronecker): Loewy = 3
])
def test_loewy_length(name, expected_loewy):
    """``loewy_length`` (nilpotency index of rad): T(kA_n)=n+1, T(k[x]/x^a)=a+1,
    T(2-Kronecker)=3 (``assem_book``)."""
    assert _TE(name, F).loewy_length() == expected_loewy


@pytest.mark.parametrize("n", [2, 3, 4])
def test_center_equals_hh0_for_TkAn(n):
    """``dim Z(T(kA_n)) = dim HH^0(T(kA_n)) = n+1`` (Brauer star ``kZ_n/J^{n+1}``;
    ``nakayama``; ``bar``)."""
    T = _TE(f"kA{n}", F)
    assert T.center()[0] == n + 1
    assert T.hochschild_cohomology(0).dims == [n + 1]
