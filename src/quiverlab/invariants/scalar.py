"""Scalar algebra invariants: Loewy length, center, complexity (spec section 3.5).

Loewy length and center are exact over ANY Domain (radical-ideal powers / commutant
linear algebra). Complexity routes through the fast GF(p) engine (the minimal A^e
resolution's term growth), so it is prime-field-only and fails loudly otherwise, exactly
like the other engine-backed invariants."""
from quiverlab.fields import linalg


def _radical_basis_indices(A):
    """Indices of the basis labels lying in rad A = the arrow ideal (all non-idempotent
    basis paths, i.e. every label that is not 'e_v')."""
    return [i for i, lab in enumerate(A.basis_labels) if not lab.startswith("e_")]


def _ideal_product_span(A, gens_vectors, rad_idx):
    """Span (as a list of coordinate vectors) of gens * rad A = { g * r : g in gens,
    r a radical basis element }, reduced to an independent set via rref."""
    dom = A.domain
    rows = []
    for g in gens_vectors:
        for ri in rad_idx:
            prod = A.multiply(g, A._basis_vec(ri))
            rows.append(prod)
    if not rows:
        return []
    R, piv = linalg.rref(rows, dom)
    return [R[i] for i in range(len(piv))]


def loewy_length(A):
    """Least n with rad^n A = 0. rad^1 = arrow ideal; rad^{k+1} = rad^k * rad."""
    dom = A.domain
    rad_idx = _radical_basis_indices(A)
    if not rad_idx:
        return 1                                    # semisimple: rad = 0, rad^1 = 0
    current = [A._basis_vec(i) for i in rad_idx]     # rad^1 spanning set
    n = 1
    while current:
        nxt = _ideal_product_span(A, current, rad_idx)
        if not nxt:
            return n + 1                             # rad^{n+1} = 0
        current = nxt
        n += 1
        if n > A.dim + 1:                            # safety: nilpotent within dim steps
            return n
    return n


def center(A):
    """(dim Z(A), basis). Z(A) = { z : z*b = b*z for every basis element b }. Solve the
    stacked commutator system over the Domain; nullspace = the center."""
    dom = A.domain
    m = A.dim
    rows = []
    # unknown z = (z_0..z_{m-1}); for each basis b, (z*b - b*z) = 0 gives m linear rows.
    # coefficient of z_k in coordinate-c of (z*b - b*z): mult(e_k,b)[c] - mult(b,e_k)[c].
    for b in range(m):
        eb = A._basis_vec(b)
        for c in range(m):
            row = []
            for k in range(m):
                ek = A._basis_vec(k)
                left = A.multiply(ek, eb)[c]         # (z*b) coord c, coefficient of z_k
                right = A.multiply(eb, ek)[c]        # (b*z) coord c
                row.append(dom.sub(left, right))
            rows.append(row)
    basis = linalg.nullspace(rows, dom)
    return len(basis), basis


def complexity(A, n):
    """Apparent complexity of A from the minimal A^e (bimodule) resolution's term-
    dimension growth up to degree n (fast GF(p) engine). Returns complexity_of's honest
    label (int / None / '>=2').

    CAVEAT — this can UNDER-REPORT; treat the number as a lower-bound estimate,
    trustworthy as EXACT only on local / single-vertex inputs:

    (a) LOCAL-ONLY radical. The growth sequence is read from
        ``engine.resolutions_minimal``, whose radical logic is LOCAL-ONLY. On a
        multi-vertex algebra of genuinely INFINITE complexity this can still return a
        FINITE value. So the result is exact only for local / single-vertex algebras;
        elsewhere read it as a lower bound. (No vertex-count guard is imposed on
        purpose — one would wrongly break the correct ``linear_path_algebra(2) -> 0``
        pin.)

    (b) SILENT TRUNCATION PREFIX. A memory-truncated resolution build contributes a
        silent prefix to the growth sequence: the build's truncation marker
        (``truncated_at`` / the discarded fourth return value below) is NOT consulted
        here, so a run that stopped early for memory reasons is read as if complete.
    """
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import minimal_resolution
    from quiverlab.engine.scan3 import complexity_of
    A._require_prime_field("complexity")             # loud FieldError off a prime field
    eng = to_engine(A)
    p = A.domain.p
    rks, cols, _e, _trunc = minimal_resolution(eng, n, p)
    # rks[k] = number of A^e generators of P_k; term k-dimension is m^2 * rks[k].
    # MINIMAL FIX vs brief: minimal_resolution returns `rks` as a dict keyed by degree,
    # so `for r in rks` iterates the KEYS (0,1,2,..) not the ranks. Iterate the rank
    # VALUES in degree order so complexity_of reads the term-dimension growth sequence.
    seq = [max(0, rks[k]) for k in sorted(rks)]
    return complexity_of(seq)
