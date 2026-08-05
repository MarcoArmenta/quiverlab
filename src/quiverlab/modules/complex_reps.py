"""Plan 35 wave 3a -- explicit Ext / Tor representatives + self-cert data.

The module counterpart of ``quiverlab.hochschild.basis_reps`` (which captures the HH
product (co)cycle representatives): here we capture, at dims-computation time, the
actual class representatives of ``Ext^n_A(M, N)`` and ``Tor_n^A(M, N)`` -- as
coordinate vectors over an ordered, LABELLED basis of the ambient Hom / tensor space
-- so a reader can read off each class AND verify it is a genuine (co)cycle.

Two coherent views of every class (mirroring ``basis_reps``):

  * ``terms``  -- the labelled term-sum. Ext: ``Sigma c * [g -> v]`` (the explicit
    hom sends the P_n-summand generator ``g`` to the N-basis vector ``v``); Tor:
    ``Sigma c * (g (x) v)`` (a tensor of a P_n generator and an N-basis vector);
  * ``vector`` -- SPARSE coordinates ``[[index, coeff_str], ...]`` over the ordered
    basis the ``chain_basis`` enumeration lists.

The ambient bases are the ADJUNCTION bases, which are the interpretable ones:

  * Ext: ``Hom_A(P_n, N)`` with ``P_n = (+)_k e_{v_k} A``; a hom out of ``e_v A`` is
    determined by the image of its generator ``e_v`` in ``N e_v``, so a basis is
    ``{phi_{k,y} : phi_{k,y}(g_k) = y, y in a basis of N e_{v_k}}``. The Hom-complex
    differential is ``delta^n(phi) = phi . d_{n+1}`` (``modules.ext._delta_matrix``);
  * Tor: ``P_n (x)_A N`` collapses summand-wise to ``(+)_k e_{v_k} N`` (Plan 29's
    ``modules.tor``); a basis is ``{g_k (x) y : y in a basis of e_{v_k} N}`` and the
    boundary is ``d_n (x) 1`` expressed in these TRUE tensor coordinates. Tor_0 is the
    cokernel ``M (x)_A N`` -- every 0-chain is a cycle, classes are its coset reps.

Because the dims are rank quantities they are basis-independent, so routing through
these adjunction bases returns EXACTLY the ``ext_dims`` / ``tor_dims`` list (asserted
by the self-check in each ``*_reps`` entry): the reps and the dims come from the same
Hom / tensor complex, never recomputed against a different basis.

Float-free: every coefficient and matrix entry is an exact Domain-element string (the
AST gate scans this file).
"""
from quiverlab.errors import QuiverlabError
from quiverlab.hochschild.basis_reps import MATRIX_CELL_CAP, serialize_differential
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.builders import projective
from quiverlab.modules.opposite import reverse_label
from quiverlab.modules.resolution import minimal_resolution


# --------------------------------------------------------------------------- #
# Labels
# --------------------------------------------------------------------------- #
def _summand_generator_labels(vertices, sym):
    """One display label per resolution-term summand: its generator, ``P_v`` (or
    ``P_v#k`` when the vertex repeats, so distinct summands never collide)."""
    counts = {}
    for v in vertices:
        counts[v] = counts.get(v, 0) + 1
    seen, out = {}, []
    for v in vertices:
        if counts[v] == 1:
            out.append("%s_%s" % (sym, v))
        else:
            seen[v] = seen.get(v, 0) + 1
            out.append("%s_%s#%d" % (sym, v, seen[v]))
    return out


def _value_label(v, j):
    """The j-th basis vector (1-based) of the vertex-v component of N."""
    return "n_%s,%d" % (v, j + 1)


def element_label(gen, val, kind):
    """One ordered-basis element as a display string. Ext (Hom): ``[g -> v]`` (the hom
    sending generator ``g`` to ``v``); Tor (tensor): ``g (x) v``."""
    if kind == "ext":
        return "[%s -> %s]" % (gen, val)
    return "%s (x) %s" % (gen, val)


def enumeration_labels(elements, kind):
    """The ordered basis enumeration as display labels (what the class vectors index
    into). Capped: an over-long enumeration ships an elided marker with the length and a
    rebuild pointer -- never silently truncated."""
    if len(elements) > MATRIX_CELL_CAP:
        return {"elided": True, "length": len(elements),
                "note": "enumeration exceeds the 250000-cell display cap; "
                        "rebuild from the ordered Hom / tensor basis"}
    return [element_label(g, v, kind) for (g, v) in elements]


# --------------------------------------------------------------------------- #
# One class = one column vector over an ordered, labelled basis enumeration.
# `elements[idx]` = (generator_label, value_label); `kind` in {"ext","tor"}.
# --------------------------------------------------------------------------- #
def serialize_class(column, elements, n, kind, dom):
    terms, vector = [], []
    for idx, coeff in enumerate(column):
        if dom.is_zero(coeff):
            continue
        cs = str(coeff)
        gen_label, val_label = elements[idx]
        terms.append([cs, gen_label, val_label])
        vector.append([idx, cs])
    return {"kind": kind, "degree": n, "terms": terms, "vector": vector}


def _classes_from_columns(columns, elements, n, kind, dom):
    return [serialize_class(col, elements, n, kind, dom) for col in columns]


def _identity_cols(dim, dom):
    return [lm.col(lm.identity(dim, dom), j) for j in range(dim)]


def _cols_of(mat):
    """The columns of a (possibly empty) matrix as column vectors."""
    if not mat or not mat[0]:
        return []
    return [lm.col(mat, j) for j in range(len(mat[0]))]


def _reps_from_complex(diff_here, diff_prev, space_dim, dom):
    """Basis of ``ker(diff_here) / im(diff_prev)`` as coordinate columns over the
    ambient degree basis. ``diff_here`` annihilates the classes (``None``/empty => the
    whole space is (co)cycles); ``diff_prev`` (``None``/empty => none) supplies the
    (co)boundaries. Since ``im diff_prev`` lies in ``ker diff_here`` (d.d = 0), the
    greedy independent-modulo pick returns exactly ``space_dim - rank_here -
    rank_prev`` classes -- the (co)homology dimension."""
    if diff_here and diff_here[0]:
        cycles = lm.kernel_columns(diff_here, dom)
    else:
        cycles = _identity_cols(space_dim, dom)
    bdys = _cols_of(diff_prev) if diff_prev is not None else []
    chosen = lm.independent_modulo(cycles, bdys, dom)
    return [cycles[i] for i in chosen]


# --------------------------------------------------------------------------- #
# Term reconstruction (shared by Ext and Tor). A resolution term ``P_n`` is the direct
# sum of the projectives P_{v_k}; rebuilding ``projective(A, v)`` reproduces the SAME
# ordered path basis the resolution's differentials use (deterministic factories), so
# the coordinates line up with ``dmats``.
# --------------------------------------------------------------------------- #
def _resolution_payload(terms):
    r"""The resolution of M the Ext/Tor engines actually walk, as display data
    (Marco 2026-08-03: the report shows the resolution BEFORE the Ext/Tor
    numbers). One latex ⊕-decomposition per degree + the Betti numbers; the same
    ``P_{v}^{k} \oplus ...`` shape the projective_resolution block uses."""
    summands, betti = [], []
    for t in terms:
        vs = list(t.vertices)
        betti.append(len(vs))
        if not vs:
            summands.append("0")
            continue
        counts = {}
        for v in vs:
            counts[v] = counts.get(v, 0) + 1
        parts = []
        for v in sorted(counts, key=lambda x: (str(type(x)), str(x))):
            base = "P_{%s}" % v
            parts.append(base if counts[v] == 1 else "%s^{%d}" % (base, counts[v]))
        summands.append(r" \oplus ".join(parts))
    return {"resolved": "M", "summands": summands, "betti": betti}


def _terms_info(A, terms):
    cache = {}

    def P(v):
        if v not in cache:
            cache[v] = projective(A, v)
        return cache[v]

    info = []
    for t in terms:
        vs = list(t.vertices)
        offs, gen_col, paths, o = [], [], [], 0
        for v in vs:
            labs = list(P(v)._pv_basis_labels)
            offs.append(o)
            paths.append(labs)
            gen_col.append(o + labs.index("e_%s" % v))
            o += len(labs)
        info.append({"vertices": vs, "offsets": offs, "gen_col": gen_col,
                     "paths": paths, "dim": o})
    return info


def _vertex_basis(N, v, dom, cache):
    """Columns spanning the vertex-v component ``N e_v = im N.action['e_v']`` (an
    idempotent projection). Cached per vertex; the SAME basis Plan 29's ``modules.tor``
    uses, so Tor coordinates match its rank bookkeeping."""
    if v not in cache:
        E = N.action["e_%s" % v]
        piv = lm.column_space_pivots(E, dom)
        cache[v] = [lm.col(E, j) for j in piv]
    return cache[v]


# --------------------------------------------------------------------------- #
# Ext
# --------------------------------------------------------------------------- #
def _hom_adjunction_basis(N, tinfo, dom, vbcache):
    """(elements, homs) for ``Hom_A(P_n, N)``: ``elements[i]=(gen_label,val_label)``,
    ``homs[i]`` = the dim(N) x dim(P_n) matrix of ``phi_{k,y}`` (generator g_k -> y,
    g_k.p -> y.p on the path basis, 0 on the other summands)."""
    gens = _summand_generator_labels(tinfo["vertices"], "P")
    elements, homs = [], []
    for k, v in enumerate(tinfo["vertices"]):
        off, paths = tinfo["offsets"][k], tinfo["paths"][k]
        for j, y in enumerate(_vertex_basis(N, v, dom, vbcache)):
            phi = lm.zeros(N.dim, tinfo["dim"], dom)
            for t, p in enumerate(paths):
                colvec = lm.matvec(N.action[p], y, dom)      # y . p (right action)
                for i in range(N.dim):
                    phi[i][off + t] = colvec[i]
            elements.append((gens[k], _value_label(v, j)))
            homs.append(phi)
    return elements, homs


def _ext_complex(A, M, N, top):
    """Shared core of :func:`ext_reps` and the Plan-41 accessor
    :func:`ext_cocycle_data`: the ``Hom_A(P_*, N)`` cochain complex of the minimal
    resolution of ``M`` through degree ``top``. Returns
    ``(terms, dmats, tinfo, elements, homs, deltas, cols_by_deg, dims)``. ONE
    construction, so the serialized reps, the ``dims``, and the Plan-41 End-action all
    read the very same ambient Hom basis (no basis mismatch)."""
    from quiverlab.modules.ext import _delta_matrix
    from quiverlab.modules.hom import _assert_comparable
    _assert_comparable(M, N, "Ext")
    dom = A.domain
    terms, dmats = minimal_resolution(M, top + 1)
    tinfo = _terms_info(M.algebra, terms)
    vbcache = {}
    elements, homs = [], []
    for ti in tinfo:
        e, h = _hom_adjunction_basis(N, ti, dom, vbcache)
        elements.append(e)
        homs.append(h)
    deltas = []
    for n in range(len(homs) - 1):
        dn1 = dmats[n + 1]
        deltas.append(_delta_matrix(homs[n], homs[n + 1], dn1, dom)
                      if (dn1 and dn1[0]) else
                      lm.zeros(len(homs[n + 1]), len(homs[n]), dom))
    cols_by_deg, dims = {}, []
    for n in range(top + 1):
        space = len(homs[n]) if n < len(homs) else 0
        here = deltas[n] if n < len(deltas) else None
        prev = deltas[n - 1] if 0 <= n - 1 < len(deltas) else None
        cols = _reps_from_complex(here, prev, space, dom) if space else []
        cols_by_deg[n] = cols
        dims.append(len(cols))
    return terms, dmats, tinfo, elements, homs, deltas, cols_by_deg, dims


def ext_cocycle_data(A, M, N, top):
    """Plan-41 internal accessor: the raw ``Ext^0..top_A(M, N)`` cohomology data WITHOUT
    serialization. Returns ``(terms, dmats, homs, cols_by_deg)`` where ``homs[n]`` is the
    ambient ``Hom(P_n, N)`` basis (matrices ``N.dim x P_n.dim``) and ``cols_by_deg[n]``
    the chosen cohomology-class coordinate columns over that basis -- the SAME complex
    :func:`ext_reps` walks (shared :func:`_ext_complex` core), so a class reconstructed
    from these columns coincides byte-for-byte with the serialized reps."""
    terms, dmats, _tinfo, _elements, homs, _deltas, cols_by_deg, _dims = _ext_complex(
        A, M, N, top)
    return terms, dmats, homs, cols_by_deg


def ext_reps(A, M, N, top, interpret=False):
    """``(dims, payload)`` for ``Ext^0..top_A(M, N)``. ``payload`` carries
    ``basis_classes`` / ``chain_basis`` / ``differentials`` keyed by ``str(degree)``
    (single side -- Ext is cohomological). ``dims`` is asserted equal to
    ``ext_dims(A, M, N, top)`` (basis independence).

    With ``interpret=True`` the payload also carries ``interpretation`` (Plan 35 wave
    3c): for each degree ``n >= 1`` class, the explicit Yoneda ``n``-fold exact
    sequence ``0 -> N -> Q -> P_{n-2} -> ... -> P_0 -> M -> 0`` realizing it -- the
    middle module ``Q`` as a full representation, every connecting map as a matrix, and
    the SELF-CERTIFIED exactness facts (a class whose sequence fails certification ships
    an honest ``error`` entry, never a wrong sequence)."""
    dom = A.domain
    terms, dmats, _tinfo, elements, homs, deltas, cols_by_deg, dims = _ext_complex(
        A, M, N, top)
    bc, cb, diffs = {}, {}, {}
    for n in range(top + 1):
        space = len(homs[n]) if n < len(homs) else 0
        here = deltas[n] if n < len(deltas) else None
        cols = cols_by_deg[n]
        elems = elements[n] if n < len(elements) else []
        bc[str(n)] = _classes_from_columns(cols, elems, n, "ext", dom)
        cb[str(n)] = enumeration_labels(elems, "ext")
        diffs[str(n)] = _ext_differential(here, space,
                                          len(homs[n + 1]) if n + 1 < len(homs) else 0, n)
    _cross_check(dims, A, M, N, top, "Ext")
    payload = {"basis_classes": bc, "chain_basis": cb, "differentials": diffs,
               "resolution": _resolution_payload(terms)}
    if interpret:
        payload["interpretation"] = _ext_interpretation(
            A, M, N, top, terms, dmats, homs, cols_by_deg)
    return dims, payload


def _ext_interpretation(A, M, N, top, terms, dmats, homs, cols_by_deg):
    """The Yoneda-extension interpretation of every degree ``n >= 1`` Ext class: build
    the ``n``-fold exact sequence, self-certify it, and serialize it. Guarded per class
    -- a construction/certification failure yields an honest ``error`` entry."""
    dom = A.domain
    sequences = {}
    for n in range(1, top + 1):
        cols = cols_by_deg.get(n) or []
        if not cols:
            continue
        width = len(homs[n][0][0]) if (n < len(homs) and homs[n]) else 0
        out = []
        for i, col in enumerate(cols):
            name = "\\alpha^{%d}_{%d}" % (n, i + 1)
            out.append(_one_ext_sequence(A, M, N, n, col, homs[n], width, terms,
                                         dmats, name))
        sequences[str(n)] = out
    return {"theory": "ext", "sequences": sequences}


def _reconstruct_cocycle(col, homs_n, n_dim, width, dom):
    """The class's actual cocycle ``f: P_n -> N`` (``n_dim x width`` matrix) as the
    ``col``-combination of the ambient Hom-basis matrices ``homs_n``."""
    f = lm.zeros(n_dim, width, dom)
    for idx, c in enumerate(col):
        if dom.is_zero(c):
            continue
        h = homs_n[idx]
        for r in range(n_dim):
            fr, hr = f[r], h[r]
            for j in range(width):
                fr[j] = dom.add(fr[j], dom.mul(c, hr[j]))
    return f


def _summand_display(mod):
    """A resolution term module's display label ``P_{v} (+) P_{w}`` from its summand
    vertices, or its dimension-vector shape when the summands are not recorded."""
    verts = getattr(mod, "_summand_vertices", None)
    if not verts:
        return "P"
    counts = {}
    for v in verts:
        counts[v] = counts.get(v, 0) + 1
    parts = []
    for v in sorted(counts, key=lambda x: (0, x) if isinstance(x, int)
                    and not isinstance(x, bool) else (1, str(x))):
        c = counts[v]
        parts.append("P_{%s}^{%d}" % (v, c) if c > 1 else "P_{%s}" % v)
    return " \\oplus ".join(parts)


def _safe_identify(mod):
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.hom import identify_standard
    try:
        return identify_standard(mod)
    except QuiverlabError:
        return None


def _one_ext_sequence(A, M, N, n, col, homs_n, width, terms, dmats, name):
    """Build + certify + serialize ONE class's Yoneda sequence. Never raises -- a loud
    library refusal (bad cocycle, oversize, undecidable) becomes an honest error entry."""
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.qpa_module import module_blocks
    from quiverlab.modules.yoneda import yoneda_sequence
    dom = A.domain
    try:
        f = _reconstruct_cocycle(col, homs_n, N.dim, width, dom)
        seq = yoneda_sequence(M, N, f, n, terms, dmats)
        ok, info = seq.check_exact()
    except QuiverlabError as exc:
        return {"class_name": name, "degree": n, "certified": False,
                "error": str(exc)}
    mods_out = []
    for mod, role in zip(seq.modules, seq.roles):
        dv = {str(v): int(d) for v, d
              in sorted(mod.dimension_vector().items(), key=lambda kv: str(kv[0]))}
        entry = {"role": role, "dim": mod.dim, "dimvec": dv}
        if role == "sub":
            entry["label"] = "N"
        elif role == "quotient":
            entry["label"] = "M"
        elif role == "middle":
            entry["label"] = "E" if n == 1 else "Q"
            entry["module"] = module_blocks(mod)
            std = _safe_identify(mod)
            if std is not None:
                entry["standard"] = {"kind": std[0], "vertex": str(std[1])}
        else:                                            # resolution_term
            entry["label"] = _summand_display(mod)
        mods_out.append(entry)
    maps_out = []
    for i, fmat in enumerate(seq.maps):
        rows = len(fmat)
        cols = len(fmat[0]) if (fmat and fmat[0]) else 0
        ser = serialize_differential((rows, cols), (lambda fm=fmat: fm),
                                     "connecting map", _Str)
        ser["from"] = mods_out[i]["label"]
        ser["to"] = mods_out[i + 1]["label"]
        maps_out.append(ser)
    return {"class_name": name, "degree": n,
            "kind": "baer" if n == 1 else "spliced",
            "certified": bool(ok),
            "modules": mods_out, "maps": maps_out,
            "facts": info if ok else [], **({} if ok else {"error": info})}


def _ext_differential(here, cols, rows, n):
    """Serialize ``delta^n : Hom(P_n, N) -> Hom(P_{n+1}, N)`` (rows x cols) -- the
    coboundary that annihilates a degree-n Ext cocycle. Capped at 250k cells."""
    mat = here if (here and here[0]) else []
    shape = (rows if mat else 0, cols)
    note = ("delta^%d = Hom(d_{%d}, N); rebuild via "
            "quiverlab.modules.complex_reps.ext_reps(A, M, N, top)" % (n, n + 1))
    return serialize_differential(shape, (lambda: mat), note, _Str)


# --------------------------------------------------------------------------- #
# Tor
# --------------------------------------------------------------------------- #
def _left_action(N, a_label):
    """Matrix of the LEFT action on the left module N of the A-path ``a_label`` (Plan
    24: N stored as a right A^op-module; left A-action = right A^op-action of the same
    element, whose A^op-label is the reversed label)."""
    key = reverse_label(a_label)
    mat = N.action.get(key)
    if mat is not None:
        return mat
    return N._action_of_word(tuple(key.split("*")))


def _tor_tensor_basis(N, tinfo, dom, vbcache):
    """(elements, offsets, dim): ordered basis of ``P_n (x)_A N = (+)_k e_{v_k} N``.
    ``elements[i]=(gen_label, val_label)`` for ``g_k (x) y``; ``offsets[k]`` is the block
    start of summand k; ``dim`` the total."""
    gens = _summand_generator_labels(tinfo["vertices"], "P")
    elements, offsets, o = [], [], 0
    for k, v in enumerate(tinfo["vertices"]):
        yb = _vertex_basis(N, v, dom, vbcache)
        offsets.append(o)
        for j in range(len(yb)):
            elements.append((gens[k], _value_label(v, j)))
        o += len(yb)
    return elements, offsets, o


def _tor_boundary(M, N, tinfo, dmats, n, tens, dom, vbcache):
    """``d_n (x) 1 : T_n -> T_{n-1}`` in TRUE tensor coordinates (n >= 1). Column of the
    source basis element ``g_k (x) y``: ``d_n(g_k) (x) y``, and ``(g_l . p) (x) y =
    g_l (x) (p . y)`` collapses each contribution into ``e_{w_l} N`` (solved back into
    that summand's ordered basis)."""
    src, tgt = tinfo[n], tinfo[n - 1]
    dmat = dmats[n]
    tgt_elems, tgt_off, tgt_dim = tens[n - 1]
    ybasis = {}                                          # target summand -> basis matrix

    def ycols(l, w):
        if l not in ybasis:
            ybasis[l] = lm.cols_to_matrix(_vertex_basis(N, w, dom, vbcache))
        return ybasis[l]

    cols = []
    for k, vk in enumerate(src["vertices"]):
        gcol = src["gen_col"][k]
        for y in _vertex_basis(N, vk, dom, vbcache):
            full = [dom.zero()] * tgt_dim
            for l, wl in enumerate(tgt["vertices"]):
                acc = [dom.zero()] * N.dim
                nonzero = False
                for t, p in enumerate(tgt["paths"][l]):
                    coeff = dmat[tgt["offsets"][l] + t][gcol]
                    if dom.is_zero(coeff):
                        continue
                    nonzero = True
                    py = lm.matvec(_left_action(N, p), y, dom)
                    acc = [dom.add(acc[i], dom.mul(coeff, py[i])) for i in range(N.dim)]
                if not nonzero:
                    continue
                coords = lm.solve_columns(ycols(l, wl),
                                          lm.cols_to_matrix([acc]), dom)
                if coords is None:
                    raise QuiverlabError(
                        "tor_reps: a boundary contribution left the vertex-%s "
                        "component of N -- the tensor collapse is inconsistent" % wl)
                base = tgt_off[l]
                for i, c in enumerate(coords[0]):
                    full[base + i] = c
            cols.append(full)
    return lm.cols_to_matrix(cols) if cols else lm.zeros(tgt_dim, 0, dom)


def tor_reps(A, M, N, top):
    """``(dims, payload)`` for ``Tor_0..top^A(M, N)`` (M right, N left). ``payload``
    carries ``basis_classes`` / ``chain_basis`` / ``differentials`` keyed by
    ``str(degree)`` (single side -- Tor is homological). Tor_0 is the cokernel
    ``M (x)_A N`` (every 0-chain is a cycle). ``dims`` is asserted equal to
    ``tor_dims(A, M, N, top)``."""
    from quiverlab.modules.tor import _assert_tor_compatible
    _assert_tor_compatible(A, M, N)
    dom = A.domain
    terms, dmats = minimal_resolution(M, top + 1)
    tinfo = _terms_info(M.algebra, terms)
    vbcache = {}
    tens = [_tor_tensor_basis(N, ti, dom, vbcache) for ti in tinfo]   # (elems, off, dim)
    bnds = {}                                            # n -> d_n (n >= 1)
    for n in range(1, len(tinfo)):
        bnds[n] = _tor_boundary(M, N, tinfo, dmats, n, tens, dom, vbcache)
    bc, cb, diffs, dims = {}, {}, {}, []
    for n in range(top + 1):
        elems, _off, tdim = tens[n] if n < len(tens) else ([], [], 0)
        here = bnds.get(n)                               # d_n : T_n -> T_{n-1}
        prev = bnds.get(n + 1)                           # d_{n+1} : T_{n+1} -> T_n
        cols = _reps_from_complex(here, prev, tdim, dom) if tdim else []
        dims.append(len(cols))
        bc[str(n)] = _classes_from_columns(cols, elems, n, "tor", dom)
        cb[str(n)] = enumeration_labels(elems, "tor")
        diffs[str(n)] = _tor_differential(here, tdim, n)
    _cross_check(dims, A, M, N, top, "Tor")
    return dims, {"basis_classes": bc, "chain_basis": cb, "differentials": diffs,
                  "resolution": _resolution_payload(terms)}


def _tor_differential(here, cols, n):
    """Serialize ``d_n : T_n -> T_{n-1}`` (the boundary annihilating a degree-n Tor
    cycle). ``d_0 = 0`` (every 0-chain is a cycle; Tor_0 = M (x)_A N is the cokernel)."""
    if n == 0 or here is None:
        note = ("d_0 = 0 (every 0-chain is a cycle; Tor_0 = M (x)_A N = coker d_1)"
                if n == 0 else "d_%d not built (beyond the resolved length)" % n)
        # a zero-row map carries no entries; keep the explanatory note visible so the
        # renderer can state "Tor_0 = M (x)_A N = coker d_1" (not present on the
        # non-elided serialize_differential shape otherwise).
        return {"shape": [0, cols], "rows": [], "note": note}
    rows = len(here)
    note = ("d_%d = d_%d (x) 1 on P_%d (x)_A N; rebuild via "
            "quiverlab.modules.complex_reps.tor_reps(A, M, N, top)" % (n, n, n))
    return serialize_differential((rows, cols), (lambda: here), note, _Str)


# --------------------------------------------------------------------------- #
# serialize_differential wants a coefficient stringifier; the module matrices already
# hold Domain elements, so `str` is exact. basis_reps._coeff_str calls dom.coerce when
# dom is not None; a tiny shim exposes just the `.coerce` used, keeping this exact.
# --------------------------------------------------------------------------- #
class _StrDom:
    @staticmethod
    def coerce(x):
        return x

    @staticmethod
    def is_zero(x):                                      # unused by serialize_differential
        return False


_Str = _StrDom()


def _cross_check(dims, A, M, N, top, op):
    """The rep pass must return EXACTLY the engine dims (basis independence). A drift is
    a bug in the capture -- refuse loudly rather than ship reps that disagree with the
    number the report prints."""
    if op == "Ext":
        from quiverlab.modules.ext import ext_dims
        engine = ext_dims(A, M, N, top)
    else:
        from quiverlab.modules.tor import tor_dims
        engine = tor_dims(A, M, N, top)
    if list(dims) != list(engine):
        raise QuiverlabError(
            "%s reps: captured dims %s disagree with the engine %s -- the explicit "
            "representatives drifted from the computed dimensions" % (op, dims, engine))
