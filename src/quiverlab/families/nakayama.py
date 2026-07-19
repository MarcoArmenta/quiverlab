"""Nakayama (serial) algebras (spec §3.4). Kupisch series K=[c_1..c_n], c_i = dim P_i.
Linear (A_n) if min(K)=1; cyclic (Z_n) if min(K)>=2. Paths read left-to-right."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import AdmissibilityError


def _validate_linear(K):
    n = len(K)
    if K[-1] != 1:
        raise AdmissibilityError(
            f"linear Kupisch series must end in 1 (sink projective is simple), got {K}",
            hint="use min(series) >= 2 for a cyclic Nakayama algebra")
    for i in range(n - 1):
        if K[i] < 2:
            raise AdmissibilityError(f"c_{i+1} = {K[i]} < 2 in the interior of {K}",
                                     hint="interior projectives have length >= 2")
        if K[i] > K[i + 1] + 1:
            raise AdmissibilityError(
                f"Kupisch admissibility fails at index {i+1}: {K[i]} > {K[i+1]} + 1",
                hint="need c_i <= c_{i+1} + 1")


def _validate_cyclic(K):
    n = len(K)
    for i in range(n):
        if K[i] < 2:
            raise AdmissibilityError(f"cyclic Kupisch entry c_{i+1} = {K[i]} < 2 in {K}",
                                     hint="every cyclic projective has length >= 2")
        if K[i] > K[(i + 1) % n] + 1:
            raise AdmissibilityError(
                f"cyclic admissibility fails at index {i+1}: {K[i]} > {K[(i+1)%n]} + 1",
                hint="need c_i <= c_{i+1} + 1 cyclically")


def NakayamaAlgebra(kupisch=None, *, n=None, l=None, cyclic=False, field=None):
    if kupisch is None:
        if n is None or l is None:
            raise AdmissibilityError(
                "give a Kupisch series, or n and l",
                hint="NakayamaAlgebra([3,2,2]) or NakayamaAlgebra(n=4, l=3, cyclic=True)")
        kupisch = [l] * n if cyclic else [min(l, n - i) for i in range(n)]  # e.g. l=3,n=4 -> [3,3,2,1]
    K = [int(c) for c in kupisch]
    if len(K) < 1 or any(c < 1 for c in K):
        raise AdmissibilityError(f"Kupisch entries must be >= 1, got {K}", hint="e.g. [3,2,2]")
    is_cyclic = min(K) >= 2
    m = len(K)
    if is_cyclic:
        _validate_cyclic(K)
        verts = list(range(1, m + 1))
        arrows = {f"a{i}": (i, i % m + 1) for i in range(1, m + 1)}   # i -> (i mod m)+1
    else:
        _validate_linear(K)
        verts = list(range(1, m + 1))
        arrows = {f"a{i}": (i, i + 1) for i in range(1, m)}           # 1->2->...->n
    Q = Quiver(verts, arrows)
    order = list(arrows)                                              # a1, a2, ...
    rels = []
    for i in range(m):                                               # length-c_i path from vertex i+1
        length = K[i]
        if is_cyclic:
            path = [order[(i + k) % m] for k in range(length)]
        else:
            if i + length > m - 1:                                   # runs off the sink (arrows a1..a_{m-1}): already zero
                continue
            path = [order[i + k] for k in range(length)]             # a_{i+1} .. (length arrows)
            if len(path) < length:
                continue
        if len(path) >= 2:
            rels.append("*".join(path))
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("nakayama", "assem_book")
    return A
