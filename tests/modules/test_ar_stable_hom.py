"""Stable Hom (mod projectives). Self-cert: underline Hom(P, -) = 0 (source
projective), and a cover composite factors through a projective. Cross-engine:
the Auslander-Reiten formula dim Ext^1(M,N) = dim underline Hom(tau^- N, M)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.ar import hom_factors_through_projective, stable_hom_dim
from quiverlab.modules.morphism import hom_basis

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _a4():
    return linear_path_algebra(4, field=QQ)


@selfcert
def test_stable_hom_from_projective_vanishes():
    A = _a4()
    P1 = A.projective(1)
    for v in (1, 2, 3, 4):
        assert stable_hom_dim(P1, A.simple(v)) == 0     # source projective


@selfcert
def test_cover_composite_factors_through_projective():
    A = _a4()
    S1 = A.simple(1)
    pi = S1.projective_cover_hom()                       # P(S1) ->> S1, factors trivially
    # any g: S1 -> P(S1) composed with pi factors through the projective P(S1)
    for g in hom_basis(S1, pi.src):
        assert hom_factors_through_projective(g.then(pi))


@xeng
def test_auslander_reiten_formula():
    A = _a4()
    simples = [A.simple(v) for v in (1, 2, 3, 4)]
    for M in simples:
        for N in simples:
            lhs = A.ext(M, N, 1)
            rhs = stable_hom_dim(N.tau_minus(), M)
            assert lhs == rhs, (M.name, N.name, lhs, rhs)
