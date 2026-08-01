"""Module Tor_n^A(M, N) for a RIGHT A-module M and a LEFT A-module N, over any exact
Domain (Plan 29 Part 4). Tor is the homological sibling of ``modules/ext.py::ext_dims``.

  Tor_n^A(M, N) = H_n(P_* (x)_A N),   P_* -> M the minimal projective resolution of the
RIGHT module M (``modules/resolution.py``). Each projective term P_i = (+)_k e_{v_k} A
collapses summand-wise under ``- (x)_A N``:

  e_v A (x)_A N  ~=  e_v N   (the vertex-v component of the LEFT module N),

whose k-dimension is the dimension-vector entry of N at v. The induced differential
d_i (x) 1 sends the summand-k generator g_k = e_{v_k} to
d_i(g_k) = sum_beta c_beta * beta, where beta ranges over the path-basis elements of
P_{i-1} (beta in summand l, path label p_beta). Tensoring and using x*a (x) n = x (x) a*n,
  g_k (x) y  |->  sum_beta c_beta (g_l (x) p_beta . y)
collapses (component l) to  sum_beta c_beta (p_beta . y) in e_{w_l} N, where p_beta . y is
the LEFT action of the path p_beta on N. So the whole complex is exact linear algebra.

Sides (Plan 24). A left A-module N is stored as a RIGHT A^op-module: N.algebra = A^op,
N.side = "left", N.base_algebra = A. The left action of an A-path with A-label L on N is
the right A^op-action of the SAME underlying element, whose A^op-label is
``reverse_label(L)`` (``opposite.py``: index set preserved, labels reversed; the
anti-homomorphism order matches the column convention action[x*y] = action[y] @ action[x]).
Hence  p . y = N.action[reverse_label(label(p))] @ y. In particular e_v N = im N.action["e_v"].

Well-definedness is certified degreewise against the EXISTING Ext engine by the classical
duality  dim Tor_n^A(M, N) = dim Ext_A^n(M, DN)  (D side-aware, Plan 24: DN of a left
module is a RIGHT module), and by the balance of Tor (``resolve="second"`` resolves N over
A^op) -- see ``tests/modules/test_tor.py``.

Sources: Cartan-Eilenberg, "Homological Algebra", Princeton Univ. Press (1956)
(citations key ``tensor_product``); Assem-Simson-Skowronski, "Elements of the
Representation Theory of Associative Algebras, Vol. 1", Cambridge Univ. Press (2006)
(key ``assem_book``).
"""
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.builders import projective
from quiverlab.modules.opposite import reverse_label
from quiverlab.modules.resolution import minimal_resolution


def _assert_tor_compatible(A, M, N):
    """Refuse loudly unless M is a RIGHT A-module and N is a LEFT A-module over the SAME
    base algebra A. Mirrors ``hom.py::_assert_comparable``, but Tor pairs OPPOSITE sides:
    M lives in mod-A, N in A-mod = mod-A^op. A wrong side or a cross-algebra pair is a
    category error -- never a silent 0."""
    if M.side != "right":
        raise QuiverlabError(
            f"Tor: the first argument must be a RIGHT A-module, got a {M.side} module",
            hint="Tor_n^A(M, N) pairs a right module M with a left module N; dualize "
                 "(D exchanges sides) or swap the arguments")
    if N.side != "left":
        raise QuiverlabError(
            f"Tor: the second argument must be a LEFT A-module, got a {N.side} module",
            hint='Tor_n^A(M, N) pairs a right module M with a left module N; build N with '
                 'side="left" (or dualize a right module)')
    if M.base_algebra is not A:
        raise QuiverlabError(
            f"Tor: the right module is over a different algebra ({M.base_algebra} vs {A})",
            hint="Tor pairs a right and a left module over ONE fixed base algebra A")
    if N.base_algebra is not A:
        raise QuiverlabError(
            f"Tor: the left module is over a different algebra ({N.base_algebra} vs {A})",
            hint="Tor pairs a right and a left module over ONE fixed base algebra A")


def _left_action(N, a_label):
    """Matrix of the LEFT action on N of the A-path with A-label ``a_label``: the right
    A^op-action of the same element = ``N.action`` at the reversed label. Falls back to
    composing arrow actions if a composite label was not materialized on N (every builder
    / ``from_arrow_action`` materializes every basis label, so the fallback is defensive)."""
    key = reverse_label(a_label)
    mat = N.action.get(key)
    if mat is not None:
        return mat
    return N._action_of_word(tuple(key.split("*")))


def _vertex_basis(N, v, dom, cache):
    """Columns spanning e_v N = im(N.action["e_v"]) (an idempotent projection); its length
    is the dimension-vector entry of N at v. Cached per vertex."""
    if v not in cache:
        E = N.action[f"e_{v}"]
        piv = lm.column_space_pivots(E, dom)
        cache[v] = [lm.col(E, j) for j in piv]
    return cache[v]


def _proj(A, v, cache):
    if v not in cache:
        cache[v] = projective(A, v)
    return cache[v]


def _induced(A, N, term_src, term_tgt, dmat, pcache, vbcache, dom):
    """Matrix of the induced differential  T_src = P_src (x)_A N  ->  T_tgt = P_tgt (x)_A N
    (rank-honest). Its columns are the images of a basis of each summand block e_{v_k} N;
    the codomain is padded to (+)_l (full N) -- an injective embedding of T_tgt, so the
    column rank is exactly rank(T_src -> T_tgt)."""
    Vsrc = term_src.vertices                        # summand vertices of P_src (dmat columns)
    Wtgt = term_tgt.vertices                         # summand vertices of P_tgt (dmat rows)
    # source offsets + generator column of each summand within P_src
    src_off, o = [], 0
    for v in Vsrc:
        src_off.append(o)
        o += _proj(A, v, pcache).dim
    # target offsets + per-summand path-basis labels within P_tgt
    tgt_off, tgt_labels, o = [], [], 0
    for w in Wtgt:
        Pw = _proj(A, w, pcache)
        tgt_off.append(o)
        tgt_labels.append(Pw._pv_basis_labels)
        o += Pw.dim
    cols = []
    for k, vk in enumerate(Vsrc):
        Pk = _proj(A, vk, pcache)
        gen_col = src_off[k] + Pk._pv_basis_labels.index(f"e_{vk}")
        for y in _vertex_basis(N, vk, dom, vbcache):        # basis of e_{vk} N
            full = []
            for l, wl in enumerate(Wtgt):
                acc = [dom.zero()] * N.dim
                for t, lab in enumerate(tgt_labels[l]):
                    coeff = dmat[tgt_off[l] + t][gen_col]
                    if dom.is_zero(coeff):
                        continue
                    contrib = lm.matvec(_left_action(N, lab), y, dom)   # (lab) . y
                    acc = [dom.add(acc[j], dom.mul(coeff, contrib[j])) for j in range(N.dim)]
                full.extend(acc)
            cols.append(full)
    return lm.cols_to_matrix(cols)


def tor_dims(A, M, N, top, resolve="first", max_term_dim=200000, with_reps=False):
    """[dim Tor_0^A(M, N), ..., dim Tor_top^A(M, N)] for a RIGHT A-module M and a LEFT
    A-module N, over any exact Domain (Plan 29 Part 4).

    ``resolve="first"`` (default) resolves M and tensors with N. ``resolve="second"``
    computes the SAME groups by resolving N over A^op (the balance of Tor):
    Tor_n^A(M, N) = Tor_n^{A^op}(N_as_right_over_Aop, M_as_left_over_Aop). ``max_term_dim``
    guards the syzygy blow-up, passed straight to the minimal resolution.

    With ``with_reps=True`` returns ``(dims, payload)`` where ``payload`` carries the
    explicit cycle representatives (``basis_classes`` / ``chain_basis`` /
    ``differentials`` per degree, Plan 35 wave 3a; Tor_0 = M (x)_A N as the cokernel)
    captured from the SAME tensor complex -- see ``modules.complex_reps.tor_reps``.
    Only the default ``resolve="first"`` supports rep capture."""
    if with_reps:
        if resolve != "first":
            raise QuiverlabError(
                'tor_dims: explicit representatives are captured only for '
                'resolve="first" (resolving the right module M)')
        from quiverlab.modules.complex_reps import tor_reps
        return tor_reps(A, M, N, top)
    _assert_tor_compatible(A, M, N)
    if resolve == "second":
        return tor_dims(A.opposite(), N.with_side("right"), M.with_side("left"),
                        top, resolve="first", max_term_dim=max_term_dim)
    if resolve != "first":
        raise QuiverlabError(
            f'tor_dims: resolve must be "first" or "second", got {resolve!r}')
    dom = A.domain
    terms, dmats = minimal_resolution(M, top + 1, max_term_dim=max_term_dim)
    pcache, vbcache = {}, {}
    # dim T_i = sum over the summand vertices of P_i of dim e_v N
    Tdim = [sum(len(_vertex_basis(N, v, dom, vbcache)) for v in term.vertices)
            for term in terms]
    # parts[n] = the induced d_{n+1}: T_{n+1} -> T_n (from dmats[n+1]: P_{n+1} -> P_n)
    parts = [_induced(A, N, terms[n + 1], terms[n], dmats[n + 1], pcache, vbcache, dom)
             for n in range(len(terms) - 1)]
    out = []
    for n in range(top + 1):
        tn = Tdim[n] if n < len(Tdim) else 0
        r_out = lm.mat_rank(parts[n], dom) if n < len(parts) else 0            # rank d_{n+1}
        r_in = lm.mat_rank(parts[n - 1], dom) if 0 <= n - 1 < len(parts) else 0  # rank d_n
        out.append(tn - r_in - r_out)
    return out


def tor(A, M, N, n):
    """dim Tor_n^A(M, N) for a RIGHT A-module M and a LEFT A-module N."""
    return tor_dims(A, M, N, n)[n]
