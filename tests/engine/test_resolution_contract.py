"""The 'Hard contract' every resolution backend must honour (CLAUDE.md / resolutions.py).

`test_resolution_protocol.py` checks that BarResolution reproduces the default engine.
This file is the *generic* guard: it runs the structural contract against EVERY backend
(the reference bar complex, the small/periodic backends, and both Chouhy-Solotar
presentations) so a newly added backend cannot silently violate it. The contract:

  * differential matrices are integer (int64), NEVER pre-reduced mod p;
  * `differential_matrix` returns shape `(dim_{n-1}, dim_n)` -- i.e. exactly
    `(len(term_basis(n-1)), len(term_basis(n)))` -- so the consumers' rank identities
    are unchanged;
  * `term_basis` returns an ordered list of hashable generators (it indexes cleanly);
  * the differentials form a complex: `d_{n-1} . d_n == 0` over the integers.

The `d.d == 0` check is computed in Python-`object` dtype so it tests the genuine
integer maps (not a mod-p reduction and not a silently-overflowing int64 product).
"""
import numpy as np
import pytest

from quiverlab.engine.hh_engine import truncated_polynomial, hochschild_homology_dims
from quiverlab.engine.resolutions import BarResolution, TruncatedPolynomialResolution
from quiverlab.engine.resolutions_periodic import CyclicNakayamaResolution

# hanlab __init__ alias, reproduced locally:
homology_dims = hochschild_homology_dims


def _backends():
    """(label, resolution, algebra) for every AVAILABLE backend, with a matching algebra.

    bar/trunc3 and truncpoly/trunc3 are always present (ported in resolutions.py). The
    Bardzell and periodic-Nakayama backends self-heal once resolutions_bardzell (Task 10)
    and coxeter2.cyclic_nakayama (Task 11) land. Three bank cases are NOT ported here:
    cs/trunc3, cs/qci2 and periodic/qci2 all depend on resolutions_cs (the Chouhy-Solotar
    closed form), which is EXCLUDED from Plan 02 -- so QuantumCIResolution is dormant.
    """
    A3 = truncated_polynomial(3)
    cases = [
        ("bar/trunc3", BarResolution(), A3),
        ("truncpoly/trunc3", TruncatedPolynomialResolution(), A3),
    ]
    try:                                                  # Bardzell minimal resolution (Task 10)
        from quiverlab.engine.resolutions_bardzell import (
            BardzellResolution, MonomialPresentation)
    except ImportError:
        return cases
    cases.append(("bardzell/trunc3",
                  BardzellResolution(MonomialPresentation.truncated_polynomial(3)), A3))
    try:                                                  # cyclic_nakayama algebras (Task 11)
        from quiverlab.engine.coxeter2 import cyclic_nakayama
    except ImportError:
        return cases
    Acn, _ = cyclic_nakayama(3, 2)
    cases.append(("periodic/cyclicNak(3,2)", CyclicNakayamaResolution(3, 2), Acn))
    return cases


_CASES = _backends()
_IDS = [c[0] for c in _CASES]
_DEPTH = 5  # term sizes stay small for the periodic/CS backends at this depth


@pytest.mark.parametrize("label,res,alg", _CASES, ids=_IDS)
def test_differential_dtype_and_shape(label, res, alg):
    bases = {n: res.term_basis(alg, n) for n in range(_DEPTH + 1)}
    index = {n: {g: i for i, g in enumerate(bases[n])} for n in range(_DEPTH + 1)}
    # term_basis must be a list of hashable generators with no collisions
    for n, b in bases.items():
        assert isinstance(b, list), f"{label}: term_basis({n}) is not a list"
        assert len(index[n]) == len(b), f"{label}: term_basis({n}) has duplicate generators"
    for n in range(1, _DEPTH + 1):
        M = res.differential_matrix(alg, n, bases[n], index[n - 1])
        assert M.dtype == np.int64, f"{label}: d_{n} dtype {M.dtype} != int64"
        expected = (len(bases[n - 1]), len(bases[n]))
        assert M.shape == expected, f"{label}: d_{n} shape {M.shape} != {expected}"


@pytest.mark.parametrize("label,res,alg", _CASES, ids=_IDS)
def test_is_a_complex(label, res, alg):
    """d_{n-1} . d_n == 0 over the integers (object dtype: no mod-p, no overflow)."""
    bases = {n: res.term_basis(alg, n) for n in range(_DEPTH + 1)}
    index = {n: {g: i for i, g in enumerate(bases[n])} for n in range(_DEPTH + 1)}
    d = {n: res.differential_matrix(alg, n, bases[n], index[n - 1]) for n in range(1, _DEPTH + 1)}
    for n in range(2, _DEPTH + 1):
        prod = d[n - 1].astype(object) @ d[n].astype(object)
        assert np.all(prod == 0), f"{label}: d_{n-1} . d_{n} != 0"


def test_not_pre_reduced_mod_p():
    """Backends must keep integer entries un-reduced so small primes carry torsion
    signal. The truncated-polynomial backend's even differential is multiply-by-`a`
    (= 3 for k[x]/(x^3)); a backend that pre-reduced mod 3 would zero it out and
    report the wrong char-3 homology. Cross-check the raw entry and the consequence.
    """
    A = truncated_polynomial(3)
    res = TruncatedPolynomialResolution()
    basis2 = res.term_basis(A, 2)
    basis1 = res.term_basis(A, 1)
    idx1 = {g: i for i, g in enumerate(basis1)}
    M = res.differential_matrix(A, 2, basis2, idx1)
    assert 3 in np.unique(M), "even differential should carry the raw entry a=3, un-reduced"
    # and the consequence: char-3 homology does NOT collapse the way a pre-reduced
    # matrix would. HH_n(k[x]/(x^3)) stays full-rank 3 in characteristic 3.
    dims = homology_dims(A, 4, resolution=res)
    assert dims[3] == [3, 3, 3, 3, 3]
