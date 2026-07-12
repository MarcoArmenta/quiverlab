"""Exact Gaussian elimination over any Domain. Correctness first; fast GF(p)
kernels replace this behind the same names in Plan 02."""


def rref(rows, dom):
    A = [list(r) for r in rows]
    nr = len(A)
    nc = len(A[0]) if nr else 0
    pivots = []
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if not dom.is_zero(A[i][c])), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = dom.inv(A[r][c])
        A[r] = [dom.mul(inv, x) for x in A[r]]
        for i in range(nr):
            if i != r and not dom.is_zero(A[i][c]):
                f = A[i][c]
                A[i] = [dom.sub(x, dom.mul(f, y)) for x, y in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
        if r == nr:
            break
    return A, pivots


def rank(rows, dom):
    return len(rref(rows, dom)[1])


def nullspace(rows, dom):
    if not rows:
        return []
    R, pivots = rref(rows, dom)
    nc = len(rows[0])
    free = [c for c in range(nc) if c not in pivots]
    basis = []
    for fc in free:
        v = [dom.zero()] * nc
        v[fc] = dom.one()
        for r, pc in enumerate(pivots):
            v[pc] = dom.neg(R[r][fc])
        basis.append(v)
    return basis


def solve(A, b, dom):
    if not A:
        return [] if all(dom.is_zero(x) for x in b) else None
    nc = len(A[0])
    aug = [list(row) + [x] for row, x in zip(A, b)]
    R, pivots = rref(aug, dom)
    if nc in pivots:  # pivot in the augmented column: inconsistent
        return None
    x = [dom.zero()] * nc
    for r, pc in enumerate(pivots):
        x[pc] = R[r][nc]
    return x
