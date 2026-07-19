"""Hom and End spaces of right A-modules over any Domain (spec §3.6, §5 c.7).

phi: M -> N is a right-module map iff N.action[b] @ phi = phi @ M.action[b] for every
generator b (arrows and idempotents). Column-stacking vec: the constraint per b is
(I_{dimM} (x) N.action[b] - M.action[b]^T (x) I_{dimN}) vec(phi) = 0. dim Hom = dim ker."""
from quiverlab.modules import linalg_mod as lm


def _generators(M):
    labels = [f"e_{v}" for v in M.algebra.quiver.vertices]
    labels += list(M.algebra.quiver.arrows)
    return labels


def hom_space(M, N):
    dom = M.domain
    dm, dn = M.dim, N.dim
    Im = lm.identity(dm, dom)
    In = lm.identity(dn, dom)
    blocks = []
    for b in _generators(M):
        Nb, Mb = N.action[b], M.action[b]
        left = lm.kron(Im, Nb, dom)                     # I_dm (x) N.action[b]
        right = lm.kron(lm.transpose(Mb), In, dom)      # M.action[b]^T (x) I_dn
        blocks.append([[dom.sub(left[i][j], right[i][j]) for j in range(dm * dn)]
                       for i in range(dm * dn)])
    stacked = lm.vstack(blocks) if blocks else lm.zeros(0, dm * dn, dom)
    ker = lm.kernel_columns(stacked, dom)
    # reshape each kernel vector (length dm*dn, column-stacked) into a dn x dm matrix
    homs = []
    for z in ker:
        phi = [[z[j * dn + i] for j in range(dm)] for i in range(dn)]
        homs.append(phi)
    return homs


def hom_dim(M, N):
    return len(hom_space(M, N))


def end_dim(M):
    return hom_dim(M, M)
