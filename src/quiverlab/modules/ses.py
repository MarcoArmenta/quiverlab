"""Short exact sequences, split test, pushout/pullback (Plan 37 / C1).

Thin exact-linear-algebra layer over :class:`quiverlab.modules.morphism.ModuleHom`:
a :class:`ShortExactSequence` certifies its own exactness at construction (the
rank identity that ``im f = ker g``), :meth:`ShortExactSequence.is_split` solves for
a section over the ``Hom(N, M)`` basis, and :func:`pushout` / :func:`pullback`
realize the two universal squares as quotient / kernel of a biproduct. Float-free.
"""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.morphism import ModuleHom, direct_sum, hom_basis


class ShortExactSequence:
    """``0 -> L --f--> M --g--> N -> 0``, exactness certified at construction:
    ``f`` mono, ``g`` epi, ``g.f = 0``, and ``rank f + rank g = dim M`` (together
    these force ``im f = ker g``)."""

    def __init__(self, f: ModuleHom, g: ModuleHom, check=True):
        if f.tgt is not g.src:
            raise QuiverlabError("not composable: f.tgt is not g.src")
        self.f, self.g = f, g
        self.L, self.M, self.N = f.src, f.tgt, g.tgt
        if check:
            ok = (f.is_mono() and g.is_epi() and f.then(g).is_zero()
                  and f.rank() + g.rank() == self.M.dim)
            if not ok:
                raise QuiverlabError(
                    "sequence is not exact "
                    "(mono/epi/g.f=0/rank identity failed)")

    def is_split(self) -> bool:
        """True iff a section ``s: N -> M`` with ``s.then(g) = id_N`` exists.
        The composite ``s.then(g)`` is linear in ``s``, so expand ``s`` over the
        ``Hom(N, M)`` basis and solve for coefficients hitting ``vec(id_N)``."""
        dom = self.M.domain
        n = self.N.dim
        if n == 0:
            return True                       # N = 0: 0 -> L -> M -> 0 splits (L ~ M)
        basis = hom_basis(self.N, self.M)
        if not basis:
            return False                      # N != 0 but no section can exist
        cols = []
        for s in basis:
            comp = s.then(self.g).matrix      # N -> N, n x n
            cols.append([comp[i][j] for j in range(n) for i in range(n)])
        ident = lm.identity(n, dom)
        target = [[ident[i][j]] for j in range(n) for i in range(n)]   # vec(id_N) column
        B = lm.cols_to_matrix(cols)           # (n*n) x len(basis)
        return lm.solve_columns(B, target, dom) is not None


def pushout(f: ModuleHom, g: ModuleHom):
    """Pushout of ``B <--f-- A --g--> C`` (shared source ``A``):
    ``P = (B (+) C) / <(f(a), -g(a))>``. Returns ``(P, inB, inC)`` with
    ``f.then(inB) == g.then(inC)`` (the pushout square commutes)."""
    if f.src is not g.src:
        raise QuiverlabError("pushout needs a shared source")
    from quiverlab.modules.yoneda import _quotient_with_maps
    D, (iB, iC), _ = direct_sum(f.tgt, g.tgt)
    dom = D.domain
    A = f.src
    diag_cols = []
    for j in range(A.dim):
        fb = [f.matrix[i][j] for i in range(f.tgt.dim)]
        gc = [dom.neg(g.matrix[i][j]) for i in range(g.tgt.dim)]
        diag_cols.append(fb + gc)             # psi(a) = (f(a), -g(a)); a module map
    P, proj_mat, _ = _quotient_with_maps(D, diag_cols, dom, name="pushout")
    proj = ModuleHom(D, P, proj_mat, check=False)
    return P, iB.then(proj), iC.then(proj)


def pullback(f: ModuleHom, g: ModuleHom):
    """Pullback of ``B --f--> D <--g-- C`` (shared target ``D``):
    ``P = ker(B (+) C --[f,-g]--> D)``. Returns ``(P, prB, prC)`` with
    ``prB.then(f) == prC.then(g)`` (the pullback square commutes)."""
    if f.tgt is not g.tgt:
        raise QuiverlabError("pullback needs a shared target")
    D, (iB, iC), (pB, pC) = direct_sum(f.src, g.src)
    dom = D.domain
    diff = ModuleHom(
        D, f.tgt,
        [[*f.matrix[i], *[dom.neg(x) for x in g.matrix[i]]]
         for i in range(f.tgt.dim)],
        check=False)                          # [f | -g] : B (+) C -> D
    P, iota = diff.kernel()
    return P, iota.then(pB), iota.then(pC)
