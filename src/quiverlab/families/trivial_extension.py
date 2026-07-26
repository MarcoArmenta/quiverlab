"""Trivial extension T(A) = A |x D(A), D(A) = Hom_k(A, k) (spec §3.4) -- symmetric.

Two builds (Plan 31):

* PRESENTED route (default when A carries a path-type basis over QQ or a prime
  field GF(p)): a genuine ``kQ_T / I_T`` algebra via ``Quiver.algebra``, so every
  path-basis invariant (``is_symmetric``, ``cartan_matrix``, modules,
  ``engine="cs"``, ...) serves T(A). ``Q_T`` has the vertices Q_0 and the arrows
  of Q plus one dual arrow per corner-homogeneous basis vector of the bimodule
  socle ``soc_{A^e}A = {x : (rad A)·x = 0 = x·(rad A)}``, with REVERSED direction
  (``w in e_i A e_j`` yields ``te : j -> i``). ``I_T`` is extracted as
  ``ker(pi: kQ_T -> A|xD(A))`` by a length-lex mini-Groebner and CERTIFIED per
  instance by ``dim(kQ_T / I_T) == 2 * dim A`` (D1/D2). The construction is not a
  transcribed theorem: every instance is certified, and QPA's native
  ``TrivialExtensionOfQuiverAlgebra`` cross-checks it.

* FALLBACK (presentation-less base, or a domain whose relation coefficients are
  not string-representable for ``Quiver.algebra`` -- CC, GF(p^n), QQ(i)): the
  structure-constant build ``_trivial_extension_structure_constants``,
  byte-identical to the pre-Plan-31 algebra. It also doubles as the
  iso-invariance oracle for the presented route.

Structure constants (fallback AND the pi target): basis ``a_0..a_{n-1}`` then
``a_0^*..a_{n-1}^*``; ``(a, f)(b, g) = (ab, a.g + f.b)``; the D(A) bimodule action
``(a.f.b)(c) = f(b.c.a)``; D(A) is a square-zero ideal.
"""
from fractions import Fraction

from quiverlab.core.algebra import Algebra
from quiverlab.errors import FieldError, NotFiniteDimensionalError, QuiverlabError

_PRESENTED_CITATIONS = ("assem_book", "skowronski_yamagata",
                        "happel_trivial_extension", "cmrs_split")


def _trivial_extension_structure_constants(A):
    """A |x D(A) as raw structure constants (dim 2*dim A). The pre-Plan-31 build,
    kept verbatim: the presented route's pi target AND the fallback for
    presentation-less / non-string-representable bases."""
    dom = A.domain
    n = A.dim
    m = 2 * n
    zero, one = dom.zero(), dom.one()

    def basis(i):
        v = [zero] * n
        v[i] = one
        return v

    # left/right action of A on itself in coordinates: prod[i][k] = a_i * a_k
    prod = [[A.multiply(basis(i), basis(k)) for k in range(n)] for i in range(n)]
    T = [[[zero] * m for _ in range(m)] for _ in range(m)]
    # (a_i, 0)(a_k, 0) = (a_i a_k, 0)
    for i in range(n):
        for k in range(n):
            for t in range(n):
                T[i][k][t] = prod[i][k][t]
    # (a_i,0)(0, a_k^*) = (0, a_i . a_k^*) ; (a_i . g)(c) = g(c a_i) => coeff on a_l^* is [a_l a_i]_k
    for i in range(n):
        for k in range(n):
            for l in range(n):
                T[i][n + k][n + l] = prod[l][i][k]
    # (0, a_i^*)(a_k, 0) = (0, a_i^* . a_k) ; (f . b)(c) = f(b c) => coeff on a_l^* is [a_k a_l]_i
    for i in range(n):
        for k in range(n):
            for l in range(n):
                T[n + i][k][n + l] = prod[k][l][i]
    # (0, *)(0, *) = 0  (D(A) is a square-zero ideal) -- already zero.
    unit = [zero] * m
    for t in range(n):
        unit[t] = A.unit[t]
    la = A.basis_labels or [f"a{i}" for i in range(n)]
    labels = list(la) + [f"{lbl}*" for lbl in la]
    Text = Algebra(dom, T, unit, basis_labels=labels)
    Text._family_citations = ("assem_book",)
    return Text


def _is_string_representable_domain(dom):
    """True iff relation coefficients over ``dom`` render as ``Quiver.algebra``
    strings: QQ (Fraction -> 'p/q') or a prime field GF(p) (int in 0..p-1). The
    GF(p^n) tuple / CC-sympy / QQ(i) domains are not, so they fall back (D3)."""
    from quiverlab.fields.primefield import PrimeField
    from quiverlab.fields.rationals import RationalField
    return isinstance(dom, (RationalField, PrimeField))


def _coeff_sign_mag(coeff):
    """(sign, magnitude_string) of an exact coefficient. QQ coefficients are
    ``Fraction`` (may be negative); prime-field coefficients are ints in 0..p-1
    (never negative -- the field already folds any minus)."""
    if isinstance(coeff, Fraction):
        sign = -1 if coeff < 0 else 1
        m = abs(coeff)
        mag = str(m.numerator) if m.denominator == 1 else f"{m.numerator}/{m.denominator}"
        return sign, mag
    return 1, str(int(coeff))


def _relation_string(terms):
    """A parallel relation (list of ``(coeff, word)`` with coeff an exact domain
    element) as a ``combinat.relations`` grammar string: coefficient BEFORE the
    arrows, signs folded into ' + ' / ' - ' term separators (a leading '-' only
    on a negative first term, which the grammar accepts)."""
    parts = []
    for idx, (coeff, word) in enumerate(terms):
        path = "*".join(word)
        sign, mag = _coeff_sign_mag(coeff)
        if idx == 0:
            sep = "" if sign > 0 else "-"
        else:
            sep = " + " if sign > 0 else " - "
        parts.append(f"{sep}{path}" if mag == "1" else f"{sep}{mag}*{path}")
    return "".join(parts)


def _solve_combo(cols, image, dom):
    """The unique ``c`` with ``sum_k c[k]*cols[k] == image`` (``cols`` are
    linearly independent by construction), or None when ``image`` is not in their
    span. An empty span with a zero ``image`` returns ``[]`` (the empty
    combination -> a monomial relation ``word``)."""
    from quiverlab.fields import linalg
    if not cols:
        return [] if all(dom.is_zero(x) for x in image) else None
    ncols = len(cols)
    mat = [[cols[k][t] for k in range(ncols)] for t in range(len(image))]
    return linalg.solve(mat, image, dom)


def _bimodule_socle(A, rad, dom):
    """A basis of ``soc_{A^e}A = {x : (rad A)·x = 0 = x·(rad A)}`` as coordinate
    vectors, via one exact nullspace over A's domain. ``rad`` empty (A
    semisimple) => the socle is all of A."""
    from quiverlab.fields import linalg
    n = A.dim
    rows = []
    for ri in rad:
        for s in range(n):
            rows.append([A.T[ri][j][s] for j in range(n)])     # (r . x)[s] = 0
        for s in range(n):
            rows.append([A.T[i][ri][s] for i in range(n)])      # (x . r)[s] = 0
    if not rows:
        return [A._basis_vec(t) for t in range(n)]
    return linalg.nullspace(rows, dom)


def _corner_maps(A, idem, rad, src, tgt):
    """Per basis index t, the (source vertex, target vertex) corner it lives in,
    read off the certified path-type structure (idem index -> vertex via the
    'e_<v>' label, the ``invariants.cartan`` idiom)."""
    verts = list(A.quiver.vertices)
    idem_vertex = {}
    for e in idem:
        lab = A.basis_labels[e]
        idem_vertex[e] = next(w for w in verts if f"e_{w}" == lab)
    vsrc, vtgt = {}, {}
    for e in idem:
        vsrc[e] = vtgt[e] = idem_vertex[e]
    for i in rad:
        vsrc[i] = idem_vertex[src[i]]
        vtgt[i] = idem_vertex[tgt[i]]
    return vsrc, vtgt


def _socle_arrows(A, socle, vsrc, vtgt, dom):
    """Split the bimodule socle into a corner-homogeneous basis and its dual
    covectors. Returns ``[(i_vertex, j_vertex, phi)]`` -- one entry per new arrow,
    where ``w_s in e_i A e_j`` (a path i -> j) with dual ``phi_s`` (``phi_s(w_t) =
    delta_st``, supported on the (i, j) corner so pi(te_s) lands in the right
    corner). Deterministic corner order (verts x verts)."""
    from quiverlab.fields import linalg
    n = A.dim
    verts = list(A.quiver.vertices)
    out = []
    for i in verts:
        for j in verts:
            cols = [t for t in range(n) if vsrc[t] == i and vtgt[t] == j]
            if not cols:
                continue
            projs = []
            for w in socle:
                if any(not dom.is_zero(w[t]) for t in cols):
                    wij = [dom.zero()] * n
                    for t in cols:
                        wij[t] = w[t]
                    projs.append(wij)
            if not projs:
                continue
            R, pivots = linalg.rref(projs, dom)
            basis_ij = [list(R[r]) for r in range(len(pivots))]
            wsub = [[b[t] for t in cols] for b in basis_ij]     # k x |cols|
            k = len(basis_ij)
            for kk in range(k):
                e_kk = [dom.one() if r == kk else dom.zero() for r in range(k)]
                phi_sub = linalg.solve(wsub, e_kk, dom)
                if phi_sub is None:                             # rows independent: cannot happen
                    raise QuiverlabError(
                        "TrivialExtension: bimodule-socle dual solve is inconsistent",
                        hint="the corner-homogeneous socle basis is not independent; "
                             "please report this presentation")
                phi = [dom.zero()] * n
                for idx, t in enumerate(cols):
                    phi[t] = phi_sub[idx]
                out.append((i, j, phi))
    return out


def _presented_trivial_extension(A, dom, idem, rad, src, tgt):
    from quiverlab.combinat.quiver import Quiver

    n = A.dim
    verts = list(A.quiver.vertices)

    # 1. bimodule socle, corner-split into new arrows + dual covectors.
    socle = _bimodule_socle(A, rad, dom)
    vsrc, vtgt = _corner_maps(A, idem, rad, src, tgt)
    arrows_info = _socle_arrows(A, socle, vsrc, vtgt, dom)

    # 2. Q_T: the arrows of Q plus one dual arrow te{s}: j -> i per socle vector
    #    (reversed), with the D4 disjointness guard.
    existing = set(A.quiver.arrows)
    prefix = ""
    while any((prefix + f"te{s}") in existing for s in range(len(arrows_info))):
        prefix += "_"
    new_names = [prefix + f"te{s}" for s in range(len(arrows_info))]
    arrows_T = dict(A.quiver.arrows)                            # originals first (insertion order)
    for s, (i, j, _phi) in enumerate(arrows_info):
        arrows_T[new_names[s]] = (j, i)                        # reversed direction
    Q_T = Quiver(verts, arrows_T)

    # 3. pi: kQ_T -> A|xD(A) on generators (arrows only -- relations live in J).
    Tsc = _trivial_extension_structure_constants(A)
    label_index = {lab: t for t, lab in enumerate(A.basis_labels)}
    img = {}
    for a in A.quiver.arrows:                                   # original arrow -> (alpha, 0)
        t = label_index[a]                                     # A-basis index of the length-1 path a
        vec = [dom.zero()] * (2 * n)
        vec[t] = dom.one()
        img[a] = vec
    for s, (_i, _j, phi) in enumerate(arrows_info):            # te{s} -> (0, phi_s)
        vec = [dom.zero()] * (2 * n)
        for t in range(n):
            vec[n + t] = phi[t]
        img[new_names[s]] = vec

    def pi_word(word):
        acc = img[word[0]]
        for a in word[1:]:
            acc = Tsc.multiply(acc, img[a])
        return acc

    def extract_relations(max_len):
        """Length-lex mini-Groebner: process words by increasing length; keep a
        per-corner reduced-echelon set of images of normal words; extend only
        normal words (on the right); a word whose image reduces against prior
        normal forms of its corner emits the parallel relation (a zero image
        emits a monomial). Deterministic (length, then discovery order)."""
        corner_basis = {}                                      # (s_vertex, t_vertex) -> [(word, image)]
        rels = []
        level = [(a,) for a in Q_T.arrows]                     # length-1 words, arrow order
        length = 1
        while level and length <= max_len:
            nxt = []
            for word in level:
                image = pi_word(word)
                corner = (Q_T.word_source(word), Q_T.word_target(word))
                basis = corner_basis.setdefault(corner, [])
                c = _solve_combo([bimg for _bw, bimg in basis], image, dom)
                if c is None:                                  # image independent -> new normal
                    basis.append((word, image))
                    tv = Q_T.word_target(word)
                    for a in Q_T.arrows:
                        if Q_T.source(a) == tv:
                            nxt.append(word + (a,))
                else:                                          # reducible -> emit the relation
                    terms = [(dom.one(), word)]
                    for kk, ck in enumerate(c):
                        if not dom.is_zero(ck):
                            terms.append((dom.neg(ck), basis[kk][0]))
                    rels.append(_relation_string(terms))
            level = nxt
            length += 1
        return rels

    # 4. I_T = ker(pi); certify dim(kQ_T / I_T) == 2*dim A (D2). The window is
    #    L(A) + 2 (comfortably past the top relation length L(A) + 1); a missed
    #    relation shows up as dim > 2n (or an infinite quotient) -- widen once.
    base_bound = A.loewy_length() + 2
    B = None
    got = None
    for bound in (base_bound, base_bound + 1):
        rel_strings = extract_relations(bound)
        try:
            cand = Q_T.algebra(relations=rel_strings, field=dom)
        except NotFiniteDimensionalError:
            got = "infinite-dimensional"
            continue
        if cand.dim == 2 * n:
            B = cand
            break
        got = cand.dim
    if B is None:
        raise QuiverlabError(
            f"TrivialExtension certificate failed: dim(kQ_T / I_T) must be "
            f"2*dim A = {2 * n}, got {got} even after widening the kernel-"
            f"enumeration window to length {base_bound + 1}",
            hint="the length-lex kernel enumeration did not capture every "
                 "relation of the trivial extension; please report this presentation")
    B._family_citations = _PRESENTED_CITATIONS
    return B


def TrivialExtension(A):
    """The trivial extension T(A) = A |x D(A) (symmetric for every f.d. A).

    Presented as a genuine ``kQ_T / I_T`` algebra when A carries a path-type basis
    over QQ or a prime field GF(p) (so all path-basis invariants serve T(A));
    otherwise the byte-identical structure-constant fallback (D3)."""
    dom = A.domain
    if not _is_string_representable_domain(dom):
        return _trivial_extension_structure_constants(A)
    from quiverlab.invariants.pathbasis import path_type_basis
    try:
        idem, rad, src, tgt = path_type_basis(A, "TrivialExtension")
    except FieldError:
        return _trivial_extension_structure_constants(A)
    return _presented_trivial_extension(A, dom, idem, rad, src, tgt)
