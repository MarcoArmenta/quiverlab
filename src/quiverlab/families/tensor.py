"""Tensor product A (x)_k B as a structure-constant algebra (spec §3.4).
Basis a_i (x) b_j in row-major order i*dim(B)+j; (a(x)b)(a'(x)b') = aa' (x) bb'."""
from quiverlab.core.algebra import Algebra
from quiverlab.errors import FieldError


def TensorProduct(A, B):
    if A.domain.name != B.domain.name:
        raise FieldError(
            f"TensorProduct needs one shared field; got {A.domain.name} and {B.domain.name}",
            hint="build both factors over the same field")
    dom = A.domain
    da, db = A.dim, B.dim
    m = da * db
    zero = dom.zero()

    def idx(i, j):
        return i * db + j

    T = [[[zero] * m for _ in range(m)] for _ in range(m)]
    for i in range(da):
        for k in range(da):
            aik = A.T[i][k]
            for j in range(db):
                for l in range(db):
                    bjl = B.T[j][l]
                    row = idx(i, j)
                    col = idx(k, l)
                    vec = T[row][col]
                    for p in range(da):
                        ap = aik[p]
                        if dom.is_zero(ap):
                            continue
                        for r in range(db):
                            br = bjl[r]
                            if dom.is_zero(br):
                                continue
                            vec[idx(p, r)] = dom.add(vec[idx(p, r)], dom.mul(ap, br))
    unit = [zero] * m
    for i in range(da):
        if dom.is_zero(A.unit[i]):
            continue
        for j in range(db):
            if dom.is_zero(B.unit[j]):
                continue
            unit[idx(i, j)] = dom.mul(A.unit[i], B.unit[j])
    la = A.basis_labels or [f"a{i}" for i in range(da)]
    lb = B.basis_labels or [f"b{j}" for j in range(db)]
    labels = [f"{la[i]}(x){lb[j]}" for i in range(da) for j in range(db)]
    C = Algebra(dom, T, unit, basis_labels=labels)
    C._family_citations = ("tensor_product", "hodge")
    return C
