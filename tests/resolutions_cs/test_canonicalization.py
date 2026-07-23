"""Plan 17: the CS correction solve is canonicalized modulo its nullspace.

reduce_mod_nullspace(x, A, dom) is the unique element of x + Null(A) with zero
coordinates at every free (non-pivot) column of A's RREF.  solve() already
returns that representative (free variables set to 0), so canonicalization is a
NO-OP today -- the point is that it is now an explicit, tested guarantee instead
of a solver-convention accident (see the WARNING block that used to live in
test_battery_bank_oracle.py).  The adversarial test below proves the CS
differential bytes no longer depend on WHICH solution the solver returns."""
import pytest

from quiverlab import Quiver, GF
from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace, reduce_mod_nullspace, solve


def _dom_cases():
    return (GF(5), QQ)


def test_coset_invariance_and_idempotence():
    """Every solution of A y = b canonicalizes to the SAME vector; applying the
    reduction twice equals applying it once."""
    for dom in _dom_cases():
        i = dom.coerce
        # rank-2 system with a 2-dim nullspace (4 unknowns)
        A = [[i(1), i(2), i(0), i(1)],
             [i(0), i(0), i(1), i(3)],
             [i(1), i(2), i(1), i(4)]]          # row3 = row1 + row2 (dependent)
        b = [i(1), i(2), i(3)]
        x0 = solve(A, b, dom)
        assert x0 is not None
        canon = reduce_mod_nullspace(x0, A, dom)
        assert canon == reduce_mod_nullspace(canon, A, dom)      # idempotent
        for v in nullspace(A, dom):
            shifted = [dom.add(a_, b_) for a_, b_ in zip(x0, v)]
            assert reduce_mod_nullspace(shifted, A, dom) == canon  # coset-invariant


def test_solver_output_is_already_canonical():
    """solve()'s free-variables-zero particular solution IS the canonical
    representative -- the no-op property that keeps every byte pin passing."""
    for dom in _dom_cases():
        i = dom.coerce
        A = [[i(1), i(2), i(0), i(1)],
             [i(0), i(0), i(1), i(3)]]
        b = [i(4), i(1)]
        x0 = solve(A, b, dom)
        assert reduce_mod_nullspace(x0, A, dom) == x0


def test_full_rank_is_untouched():
    """No nullspace -> the reduction returns the input unchanged."""
    dom = GF(7)
    i = dom.coerce
    A = [[i(1), i(1)], [i(0), i(1)]]
    b = [i(3), i(2)]
    x0 = solve(A, b, dom)
    assert reduce_mod_nullspace(x0, A, dom) == x0


def _qci(field=None):
    from quiverlab.groebner import build_reduction_system
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    f = field or GF(5)
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", "y*x - 2*x*y"]
    return ChouhySolotarResolution(Q.algebra(relations=rels, field=f),
                                   build_reduction_system(Q, rels, f), max_degree=7)


def _correction_nullities(res, upto):
    """Nullities of the correction systems actually solved while building d_1..d_upto."""
    from quiverlab.resolutions_cs.pelt import terms_to_pelt, apply_lower
    dom = res.dom
    out = [0]
    for n in range(1, upto + 1):
        for chain in res.ss.S(n):
            gens = res._lower_generators(n, chain)
            if not gens:
                continue
            delta = res.delta_terms(n, chain)
            rhs_pe = apply_lower(res, n, terms_to_pelt(res, delta))
            cols = [apply_lower(res, n, terms_to_pelt(res, [g])) for g in gens]
            keys = sorted(set(rhs_pe) | {k for col in cols for k in col})
            M = [[col.get(k, dom.zero()) for col in cols] for k in keys]
            out.append(len(nullspace(M, dom)) if M else len(gens))
    return out


def test_correction_nullity_is_nonzero_somewhere():
    """The adversarial battery is NOT vacuous: on the quantum CI some degree has a
    correction system with genuine nullspace freedom (the bank docstring measured
    nullity growing with degree)."""
    res = _qci()
    assert max(_correction_nullities(res, upto=6)) > 0


def test_cs_bytes_survive_an_adversarial_solver(monkeypatch):
    """THE Plan-17 theorem: shift the solver's answer by a nullspace vector (a
    DIFFERENT valid solution) -- the built differentials must be byte-identical
    anyway, because _d_general canonicalizes.  Before this plan the bytes moved
    (the bank-oracle WARNING); now they cannot."""
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    baseline = _qci()
    ref = {n: [sorted((baseline.to_int(c), a, t, cc)
                      for (c, a, t, cc) in baseline.d_terms(n, q))
               for q in baseline.ss.S(n)] for n in range(1, 7)}

    real_solve = ChouhySolotarResolution._solve
    shifted = {"count": 0}

    def adversarial_solve(self, M, rhs, ncols):
        sol = real_solve(self, M, rhs, ncols)
        if sol is None or not M:
            return sol
        basis = nullspace(M, self.dom)
        if not basis:
            return sol
        shifted["count"] += 1
        return [self.dom.add(s, v) for s, v in zip(sol, basis[0])]

    monkeypatch.setattr(ChouhySolotarResolution, "_solve", adversarial_solve)
    perturbed = _qci()
    got = {n: [sorted((perturbed.to_int(c), a, t, cc)
                      for (c, a, t, cc) in perturbed.d_terms(n, q))
               for q in perturbed.ss.S(n)] for n in range(1, 7)}
    assert shifted["count"] > 0            # the adversary actually fired
    assert got == ref                      # ...and the bytes did not move
    perturbed.assert_dd_zero(upto=6, side="hom")
    perturbed.assert_order_condition(upto=6)
