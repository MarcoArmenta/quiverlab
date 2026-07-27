"""Module-computation worked-step builders (Plan 30 Part C; Plan 34 homework depth).

These are the sanctioned ``quiverlab.trace`` hooks for the module surface: given a
module computation they RUN it through the public surface and emit the typed step
events (``ModuleTerm`` / ``ModuleDifferential`` / ``ExtDegree`` / ``StepNote``, see
``trace.events``) that the renderers turn into an exhaustive, replay-by-hand
pdf/tex/html bundle.

Plan 34 (Marco's feedback -- "say how every object is computed as I would demand an
undergraduate student in his homework, writing everything and with justifications"):
every traced object now states, in full sentences, WHAT is computed, the DEFINITION
used, the ACTUAL computation (matrices shown), and WHY the conclusion follows, with a
one-line literature justification (routed through the citation registry:
``assem_book`` for the module theory, ``minimal_resolution`` / ``module_ext`` for
resolutions and Ext, ``tensor_product`` for Tor). A ``StepNote`` with ``heading=True``
opens a numbered worked step; the LaTeX renderer numbers them ``Step N.``.

Design (mirrors how ``hochschild.bar`` hooks the HH engine):
  * ``trace_radical`` / ``trace_top`` / ``trace_socle`` show the arrow-action matrices
    and the span-of-images / quotient / annihilator linear algebra behind rad/top/soc;
  * a *projective resolution* trace reads the PUBLIC ``ProjectiveResolution`` object --
    its ``.term(n)`` summand vertices and ``.differential(n)`` matrices -- and narrates
    the projective cover + syzygy at each step;
  * an *Ext* / *Tor* trace re-runs the same pipeline (minimal resolution + Hom / tensor
    collapse), emitting the per-degree rank bookkeeping with the dimension count spelled
    out; cross-check tests pin the dims to ``ext_dims`` / ``tor_dims``;
  * a ``tau`` trace shows the projective presentation P_1 -> P_0, the transpose (Hom
    over A^op), the duality D and the result; a ``decompose`` trace shows the
    endomorphism-ring split (or the locality certificate, honestly labelled).

Matrix bodies are elided past ``recorder.MATRIX_ELISION_CELLS`` (shape + rank kept) --
every computational step still appears, never silently omitted. Float-free: all data
are ints and domain-element strings.
"""
from quiverlab.errors import QuiverlabError
from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm
from quiverlab.trace.events import ModuleTerm, StepNote
from quiverlab.trace.recorder import module_differential, ext_degree


# --------------------------------------------------------------------------- #
# Small shared helpers (all float-free: ints + domain-element strings).
# --------------------------------------------------------------------------- #

def _side_words(M):
    """Side-aware narration fragments (Plan 24 sides; Plan 34 MAJOR-4). A RIGHT
    A-module has ``rad M = M J`` and arrows acting on the right (``m |-> m.a``); a LEFT
    A-module has ``rad M = J M`` and arrows acting on the left (``m |-> a.m``), and the
    socle annihilator is ``J m = 0`` rather than ``m J = 0``. The stored matrices are the
    same representation either way -- only the prose changes, so the algorithms stay
    side-blind while the homework report reads correctly for the user's chosen side."""
    if getattr(M, "side", "right") == "left":
        return {"side": "left", "rad": "rad M = J M", "act": "m |-> a.m",
                "annih": "J m = 0", "otherside": "right"}
    return {"side": "right", "rad": "rad M = M J", "act": "m |-> m.a",
            "annih": "m J = 0", "otherside": "left"}


def _reduction_note(symbol, piv, rank, dom):
    """The one-line column-reduction summary shown after every displayed matrix we reduce
    (Plan 34 MINOR-8): the pivot columns and the resulting rank -- the homework-acceptable
    middle ground (no full RREF display). ``symbol`` is PLAIN text (StepNotes are TeX-
    escaped, so no LaTeX macros here)."""
    return StepNote(
        "Column-reducing %s over %s: pivots at column(s) %s, so rank %s = %d."
        % (symbol, dom.name, list(piv), symbol, rank))


def _dv(dimvec):
    """A dimension vector as a str-keyed, vertex-sorted dict (golden-stable)."""
    return {str(v): int(n) for v, n in sorted(dimvec.items(), key=lambda kv: str(kv[0]))}


def _mat_shape(D):
    return (len(D), len(D[0]) if (D and D[0]) else 0)


def _rank(D, dom):
    return lm.mat_rank(D, dom) if (D and D[0]) else 0


def _fmt_dimvec(dv):
    """A plain-text dimension vector ``(1:2, 3:1)`` (shared narrative formatting)."""
    return "(" + ", ".join("%s:%s" % (k, v) for k, v in dv.items()) + ")"


def _factor_stack_str(dv):
    """A semisimple module as a composition-factor sum ``S_1 (+) S_2^2`` (prose)."""
    parts = ["S_%s" % v if m == 1 else "S_%s^%d" % (v, m)
             for v, m in dv.items() if m > 0]
    return " (+) ".join(parts) if parts else "0"


def _m_vertices(M):
    """The vertices where M is supported (dim M e_v > 0), in the quiver order, each
    listed once: the k-vector-space decomposition M = (+)_v M e_v as a summand list for
    ``oplus_tex`` (sym ``Me`` -> ``Me_v``)."""
    dv = M.dimension_vector()
    return [v for v in M.algebra.quiver.vertices if dv.get(v, 0) > 0]


def _vertex_basis(M, v):
    """A k-basis of the vertex component M e_v as column vectors of M (the pivot columns
    of the idempotent projection e_v)."""
    dom = M.domain
    P = M.vertex_projection(v)
    piv = lm.column_space_pivots(P, dom)
    return [lm.col(P, j) for j in piv]


def _arrow_symbol(a, name="M"):
    """The arrow-action map symbol used in the matrix header (raw TeX). Arrow names in
    quiverlab quivers are alphanumeric, so they are TeX-safe as written; the action of
    the arrow ``a`` on the traced module is written ``\\rho_M(a)`` (right
    multiplication m |-> m.a), with the module's actual name in the subscript."""
    sub = name if len(name) == 1 else "{%s}" % name
    return r"\rho_%s(%s)" % (sub, a)


def _arrow_action_events(M):
    """One ``ModuleDifferential`` per arrow a: the matrix of the right action
    ``rho_M(a): M -> M`` (m |-> m.a) in the fixed k-basis of M. These are the raw data
    every rad/top/soc computation is read off."""
    dom = M.domain
    verts = _m_vertices(M)
    out = []
    for a in M.algebra.quiver.arrows:
        Aa = M.action[a]
        out.append(module_differential(
            degree=0, kind="module", sym="Me", symbol=_arrow_symbol(a, M.name),
            dom_summands=verts, cod_summands=verts, D=Aa,
            nrows=M.dim, ncols=M.dim, dom=dom, cod_is_module=True,
            dom_is_module=True, mod_name=M.name, rank=_rank(Aa, dom)))
    return out


def _poly_str(coeffs, dom):
    """An ascending Domain coefficient list ``[a0, a1, ..., 1]`` as a readable monic
    polynomial ``x^d + ... + a0`` (plain text; escaped-caret renders as a literal
    caret in the PDF, which is fine for a polynomial)."""
    terms = []
    for i in range(len(coeffs) - 1, -1, -1):
        c = coeffs[i]
        if dom.is_zero(c):
            continue
        cs = str(c)
        if i == 0:
            terms.append(cs)
        elif i == 1:
            terms.append("x" if cs == "1" else "%s*x" % cs)
        else:
            terms.append("x^%d" % i if cs == "1" else "%s*x^%d" % (cs, i))
    return " + ".join(terms) if terms else "0"


# --------------------------------------------------------------------------- #
# Radical, top, socle: the arrow-action linear algebra, shown in full.
# --------------------------------------------------------------------------- #

def trace_radical(M):
    """Worked steps for rad M = M J: the arrow-action matrices, the assembled matrix of
    their image columns, its column reduction, and the resulting basis / dimension
    vector. Returns ``(events, radM)``."""
    dom = M.domain
    R = M.radical()
    sw = _side_words(M)
    ev = [StepNote(
        "Radical of %s" % M.name,
        "By definition the radical %s, where J = rad A is the arrow ideal (the Jacobson "
        "radical) of A. Since A/J is semisimple, J annihilates every simple module, so "
        "rad M is the smallest submodule with semisimple quotient; for a bound quiver "
        "algebra this is exactly the sum of the images of the arrow actions, rad M = sum "
        "over the arrows a of image(%s) (Assem-Simson-Skowronski, Elements of the "
        "Representation Theory of Associative Algebras I, Cambridge 2006, I.3). We first "
        "record how each arrow acts on M." % (sw["rad"], sw["act"]),
        heading=True)]
    ev += _arrow_action_events(M)
    arrows = list(M.algebra.quiver.arrows)
    cols, dom_summands = [], []
    for a in arrows:
        s = M.algebra.quiver.source(a)
        for u in _vertex_basis(M, s):
            cols.append(lm.matvec(M.action[a], u, dom))
        if any(True for _ in _vertex_basis(M, s)):
            dom_summands.append(s)
    G = lm.cols_to_matrix(cols) if cols else lm.zeros(M.dim, 0, dom)
    ncols = len(cols)
    piv = lm.column_space_pivots(G, dom) if (G and G[0]) else []
    r = len(piv)
    if r != R.dim:
        raise QuiverlabError(
            "trace_radical: the shown column-rank %d of the arrow-image matrix disagrees "
            "with the engine's dim rad %s = %d -- the worked steps drifted from the "
            "computation" % (r, M.name, R.dim))
    if ncols:
        ev.append(StepNote(
            "The images of the arrows are the (nonzero) columns above. Assembling them "
            "into one matrix G, the radical is the column span rad M = colspace(G)."))
        ev.append(module_differential(
            degree=0, kind="module", sym="Me", symbol="G",
            dom_summands=dom_summands, cod_summands=[], D=G,
            nrows=M.dim, ncols=ncols, dom=dom, cod_is_module=True, rank=r))
        ev.append(_reduction_note("G", piv, r, dom))
        ev.append(StepNote(
            "Therefore rad %s has dimension %d and dimension vector %s; as a submodule "
            "of M it is %s." % (M.name, R.dim, _fmt_dimvec(_dv(R.dimension_vector())),
                                _radical_series_hint(R))))
    else:
        ev.append(StepNote(
            "There are no arrows, so every arrow action is trivial and rad %s = 0 "
            "(M is semisimple)." % M.name))
    return ev, R


def _radical_series_hint(R):
    dv = _dv(R.dimension_vector())
    return "the submodule with composition factors %s" % _factor_stack_str(dv) \
        if R.dim else "0"


def trace_top(M):
    """Worked steps for top M = M / rad M: the per-vertex multiplicity arithmetic and
    the resulting semisimple quotient. Returns ``(events, topM)``."""
    R = M.radical()
    T = M.top()
    if T.dim != M.dim - R.dim:
        raise QuiverlabError(
            "trace_top: the shown quotient dimension %d = dim %s - dim rad %s disagrees "
            "with the engine's dim top %s = %d" % (M.dim - R.dim, M.name, M.name,
                                                   M.name, T.dim))
    dv_M, dv_R, dv_T = (_dv(M.dimension_vector()), _dv(R.dimension_vector()),
                        _dv(T.dimension_vector()))
    arithmetic = "; ".join(
        "at vertex %s: %d - %d = %d" % (v, dv_M.get(v, 0), dv_R.get(v, 0), dv_T.get(v, 0))
        for v in sorted(dv_M))
    ev = [StepNote(
        "Top of %s" % M.name,
        "The top top M = M / rad M is the largest semisimple quotient of M "
        "(Assem-Simson-Skowronski I, 2006). Its multiplicity of the simple S_v equals "
        "dim (M e_v) - dim (rad M e_v); coset representatives are the standard basis "
        "vectors of M that are linearly independent modulo rad M.", heading=True),
        StepNote(
        "Per-vertex multiplicities of the simples in top M -- %s." % arithmetic,
        "Hence top %s = %s, a semisimple module of dimension %d with dimension vector "
        "%s." % (M.name, _factor_stack_str(dv_T), T.dim, _fmt_dimvec(dv_T)))]
    return ev, T


def trace_socle(M):
    """Worked steps for soc M = {m : m J = 0} = intersection of the arrow kernels: the
    stacked homogeneous system (each arrow restricted to its target component) and its
    nullspace. Returns ``(events, socM)``."""
    dom = M.domain
    S = M.socle()
    sw = _side_words(M)
    ev = [StepNote(
        "Socle of %s" % M.name,
        "The socle soc M = { m in M : %s } = the intersection over the arrows a of "
        "ker(%s); it is the largest semisimple submodule of M "
        "(Assem-Simson-Skowronski I, 2006). Concretely it is the simultaneous nullspace "
        "of the arrow-action matrices: stack them into one system and solve."
        % (sw["annih"], sw["act"]),
        heading=True)]
    arrows = list(M.algebra.quiver.arrows)
    blocks, cod_summands = [], []
    for a in arrows:
        t = M.algebra.quiver.target(a)
        Bt = _vertex_basis(M, t)
        if not Bt:
            continue
        Bt_mat = lm.cols_to_matrix(Bt)
        # rows = coordinates of (m.a) in the M e_t basis (m.a lives in M e_t)
        coeffs = lm.solve_columns(Bt_mat, M.action[a], dom)
        blocks.append(lm.cols_to_matrix(coeffs))
        cod_summands.append(t)
    if blocks:
        K = lm.vstack(blocks)
        nrows = len(K)
        rank_K = _rank(K, dom)
        if M.dim - rank_K != S.dim:
            raise QuiverlabError(
                "trace_socle: the shown nullity %d = dim %s - rank K disagrees with the "
                "engine's dim soc %s = %d" % (M.dim - rank_K, M.name, M.name, S.dim))
        pivK = lm.column_space_pivots(K, dom) if (K and K[0]) else []
        ev.append(module_differential(
            degree=0, kind="module", sym="Me", symbol="K",
            dom_summands=_m_vertices(M), cod_summands=cod_summands, D=K,
            nrows=nrows, ncols=M.dim, dom=dom, cod_is_module=False,
            dom_is_module=True, mod_name=M.name, rank=rank_K))
        ev.append(_reduction_note("K", pivK, rank_K, dom))
        ev.append(StepNote(
            "The nullspace of this stacked system (dimension dim %s - rank K = %d) is "
            "exactly soc %s." % (M.name, M.dim - rank_K, M.name)))
    else:
        # No arrow acts into a nonzero vertex-component of M, so m.a = 0 holds
        # automatically for every m: the annihilator condition is vacuous.
        ev.append(StepNote(
            "Every arrow maps %s into a vertex where it vanishes, so m.a = 0 holds "
            "automatically for all m; the annihilator condition is vacuous and "
            "soc %s = %s." % (M.name, M.name, M.name)))
    dv_S = _dv(S.dimension_vector())
    ev.append(StepNote(
        "Solving, the socle has dimension %d and dimension vector %s; as a semisimple "
        "module soc %s = %s." % (S.dim, _fmt_dimvec(dv_S), M.name, _factor_stack_str(dv_S))))
    return ev, S


# --------------------------------------------------------------------------- #
# The projectives and injectives of A (Marco #4): the P/I section data.
# --------------------------------------------------------------------------- #

def _radical_layers(M):
    """The Loewy (radical) layers of a module, top to bottom, each a str-keyed
    composition-factor multiplicity dict: [top(M), top(rad M), top(rad^2 M), ...].
    Computed through the PUBLIC ``radical()`` / ``top()`` surface only."""
    layers = []
    cur = M
    while cur.dim > 0:
        layers.append(_dv(cur.top().dimension_vector()))
        r = cur.radical()
        if r.dim >= cur.dim:            # radical must strictly shrink for f.d. modules
            break
        cur = r
    return layers


def algebra_objects(A):
    """For each vertex v: the indecomposable projective P_v and injective I_v with
    dimension vector, Loewy (radical) layers, top and socle. Uses only the public
    module surface (``A.projective`` / ``A.injective`` + radical/top/socle). The
    simples S_v are deliberately omitted (Marco: obvious). Returns a list of dicts,
    one per vertex, all data str/int (float-free)."""
    out = []
    for v in A.quiver.vertices:
        row = {"vertex": str(v)}
        for sym, mod in (("P", A.projective(v)), ("I", A.injective(v))):
            row[sym] = {
                "dimvec": _dv(mod.dimension_vector()),
                "dim": mod.dim,
                "layers": _radical_layers(mod),
                "top": _dv(mod.top().dimension_vector()),
                "socle": _dv(mod.socle().dimension_vector()),
            }
        out.append(row)
    return out


# --------------------------------------------------------------------------- #
# Projective resolution trace (the acceptance bar: every differential verbatim).
# --------------------------------------------------------------------------- #

def trace_projective_resolution(M, top):
    """Worked-step events for the minimal projective resolution of ``M`` to length
    ``top``: the projective-cover overview, then each term Q_n = (+) P_v with its
    dimension vector, interleaved with the differential OUT of it (the augmentation
    epsilon: Q_0 -> M, then d_n: Q_n -> Q_{n-1}) rendered as a matrix, and the syzygy
    that each term covers. Returns ``(events, res)`` where ``res`` is the public
    ``ProjectiveResolution``."""
    res = M.projective_resolution(top)
    dom = M.domain
    events = [StepNote(
        "Minimal projective resolution of %s" % M.name,
        "We build ... -> Q_1 -> Q_0 -> M -> 0 by iterated projective covers "
        "(Green-Solberg-Zacharia, Minimal projective resolutions, Trans. AMS 353 "
        "(2001) 2915-2939): cover M by Q_0 = (+) P_v over a basis of top M, take the "
        "syzygy Omega_1 = ker(Q_0 -> M), cover that, and repeat. Choosing the "
        "generators modulo the radical keeps the resolution minimal, so the number of "
        "P_v summands in Q_n is the true n-th Betti number.", heading=True),
        StepNote(
        "Projective cover: top %s has composition factors %s, so its projective cover "
        "is Q_0 = %s (one P_v per simple in the top)."
        % (M.name, _fmt_dimvec(_dv(M.top().dimension_vector())),
           _oplus_str(res.term(0))))]
    n_terms = res.length
    for n in range(n_terms):
        t = res.terms[n]
        verts = res.term(n)
        dimvec = _dv(t.module.dimension_vector()) if t.module is not None else {}
        # narrate the syzygy this term covers (n >= 1): Omega_n = ker(d_{n-1})
        if n >= 1 and (n - 1) < len(res.dmats):
            dprev = res.differential(n - 1)
            ncols_prev = _mat_shape(dprev)[1]
            nullity = ncols_prev - _rank(dprev, dom)
            if verts:
                events.append(StepNote(
                    "Syzygy: Omega_%d = ker(%s) has dimension %d (the nullity of the "
                    "previous differential). Its top needs %d generator(s), so it is "
                    "covered by Q_%d = %s."
                    % (n, (r"epsilon" if n == 1 else "d_%d" % (n - 1)), nullity,
                       len(verts), n, _oplus_str(verts))))
            else:
                events.append(StepNote(
                    "Syzygy: Omega_%d = ker(%s) has dimension %d; it is zero (or "
                    "projective), so the resolution stops here."
                    % (n, (r"epsilon" if n == 1 else "d_%d" % (n - 1)), nullity)))
        events.append(ModuleTerm(degree=n, kind="projective", sym="P",
                                 summands=list(verts), dim=t.dim, dimvec=dimvec))
        # the differential out of Q_n, if recorded: d_0 = epsilon (Q_0 -> M),
        # d_n = (Q_n -> Q_{n-1}).
        if n < len(res.dmats):
            D = res.differential(n)
            nrows, ncols = _mat_shape(D)
            cod_is_module = (n == 0)
            symbol = r"\varepsilon" if n == 0 else "d_{%d}" % n
            rankD = _rank(D, dom)
            events.append(module_differential(
                degree=n, kind="projective", sym="P", symbol=symbol,
                dom_summands=verts,
                cod_summands=[] if cod_is_module else res.term(n - 1),
                D=D, nrows=nrows, ncols=ncols, dom=dom,
                cod_is_module=cod_is_module, rank=rankD))
            pivD = lm.column_space_pivots(D, dom) if (D and D[0]) else []
            events.append(_reduction_note("epsilon" if n == 0 else "d_%d" % n,
                                          pivD, rankD, dom))
    pd = res.pd()
    events.append(StepNote(
        "The resolution terminates once a syzygy is zero, so the projective dimension is "
        "pd(%s) = %s." % (M.name, ("infinite (unresolved within the requested length)"
                                   if pd is None else pd))))
    return events, res


def trace_injective_resolution(M, top):
    """Worked-step events for the minimal injective coresolution of ``M`` to length
    ``top``: each term E^n = (+) I_v with its dimension vector, plus a narration of the
    duality E^n = D(P_n over A^op). Returns ``(events, res)`` with the public
    ``InjectiveResolution``."""
    res = M.injective_resolution(top)
    events = [StepNote(
        "Minimal injective coresolution of %s" % M.name,
        "Dually to the projective case, 0 -> M -> E^0 -> E^1 -> ... is obtained by "
        "dualizing the minimal projective resolution of D M over A^op: E^n = D(P_n), "
        "P_* the minimal projective resolution of D%s over A^op "
        "(Assem-Simson-Skowronski I, 2006). D preserves dimension vectors, so the "
        "injective summands are read off the P_n." % M.name, heading=True)]
    for n in range(res.length):
        t = res.terms[n]
        verts = res.term(n)
        dim = t.dim if t is not None else 0
        dimvec = _dv(t.dimension_vector()) if t is not None else {}
        events.append(ModuleTerm(degree=n, kind="injective", sym="I",
                                 summands=list(verts), dim=dim, dimvec=dimvec))
    idim = res.injective_dimension()
    events.append(StepNote(
        "The coresolution terminates once a term is zero, so the injective dimension is "
        "id(%s) = %s." % (M.name, ("infinite (unresolved within the requested length)"
                                   if idim is None else idim))))
    return events, res


# --------------------------------------------------------------------------- #
# Ext trace: the Hom-collapse matrices + per-degree rank bookkeeping.
# --------------------------------------------------------------------------- #

def trace_ext(A, M, N, top):
    """Worked-step events for Ext^0..top_A(M, N): the minimal resolution of M, the
    Hom(Q_n, N) space dimensions, the connecting maps delta^n rendered as matrices, and
    Ext^n = dim Hom(Q_n,N) - rank(delta^n) - rank(delta^{n-1}) with the arithmetic
    spelled out. Re-runs the ``ext_dims`` pipeline (a cross-check test pins the emitted
    dims to it). Returns ``(events, dims)``."""
    from quiverlab.modules.ext import _delta_matrix
    from quiverlab.modules.hom import _assert_comparable, hom_space
    from quiverlab.modules.resolution import minimal_resolution
    _assert_comparable(M, N, "Ext")
    dom = A.domain
    terms, dmats = minimal_resolution(M, top + 1)
    Qs = [t.module for t in terms]
    Homs = [hom_space(Q, N) if (Q is not None and Q.dim) else [] for Q in Qs]
    deltas = []
    for n in range(len(Qs) - 1):
        dn1 = dmats[n + 1]
        deltas.append(_delta_matrix(Homs[n], Homs[n + 1], dn1, dom)
                      if (dn1 and dn1[0]) else
                      lm.zeros(len(Homs[n + 1]), len(Homs[n]), dom))
    events = [StepNote(
        "The Ext groups Ext(%s, %s) in each degree" % (M.name, N.name),
        "Ext^n_A(M, N) = H^n(Hom_A(Q_*, N)), where Q_* -> M is the minimal projective "
        "resolution above (Green-Solberg-Zacharia 2001). Applying Hom_A(-, N) and using "
        "Hom_A(P_v, N) = N e_v (a map out of P_v is determined by the image of its "
        "generator e_v) collapses the complex to finite-dimensional vector spaces; the "
        "n-th cohomology is Ext^n.", heading=True)]
    dims = []
    for n in range(top + 1):
        space_dim = len(Homs[n]) if n < len(Homs) else 0
        Dn = deltas[n] if n < len(deltas) else None
        r_n = _rank(Dn, dom) if Dn is not None else 0
        r_nm1 = _rank(deltas[n - 1], dom) if 0 <= n - 1 < len(deltas) else 0
        res_dim = space_dim - r_n - r_nm1
        dims.append(res_dim)
        D = Dn if Dn is not None else lm.zeros(0, space_dim, dom)
        nrows, ncols = _mat_shape(D)
        events.append(ext_degree(
            degree=n, op="Ext", space_dim=space_dim, rank_here=r_n,
            rank_prev=r_nm1, result_dim=res_dim, D=D, nrows=nrows, ncols=ncols,
            dom=dom))
    from quiverlab.modules.ext import ext_dims
    engine_dims = ext_dims(A, M, N, top)
    if dims != engine_dims:
        raise QuiverlabError(
            "trace_ext: the shown Ext dimensions %s disagree with the engine's ext_dims "
            "%s for Ext(%s, %s) -- the worked steps drifted from the computation"
            % (dims, engine_dims, M.name, N.name))
    return events, dims


# --------------------------------------------------------------------------- #
# Tor trace: the tensor-collapse matrices + per-degree rank bookkeeping.
# --------------------------------------------------------------------------- #

def trace_tor(A, M, N, top):
    """Worked-step events for Tor_0..top^A(M, N) with M a RIGHT and N a LEFT A-module:
    the minimal resolution of M, the tensor spaces P_n (x)_A N (= (+) N e_v over the
    summands of Q_n), the induced differentials d_n rendered as matrices, and Tor_n =
    dim(P_n (x) N) - rank(d_n) - rank(d_{n+1}) with the arithmetic spelled out. Re-runs
    the ``tor_dims`` pipeline (a cross-check test pins the emitted dims to it). Returns
    ``(events, dims)``."""
    from quiverlab.modules.tor import _assert_tor_compatible, _induced, _vertex_basis as _tvb
    from quiverlab.modules.resolution import minimal_resolution
    _assert_tor_compatible(A, M, N)
    dom = A.domain
    terms, dmats = minimal_resolution(M, top + 1)
    pcache, vbcache = {}, {}
    Tdim = [sum(len(_tvb(N, v, dom, vbcache)) for v in t.vertices) for t in terms]
    parts = [_induced(A, N, terms[n + 1], terms[n], dmats[n + 1], pcache, vbcache, dom)
             for n in range(len(terms) - 1)]
    events = [StepNote(
        "The Tor groups Tor(%s, %s) in each degree" % (M.name, N.name),
        "Tor_n^A(M, N) = H_n(Q_* (x)_A N), where Q_* -> M is the minimal projective "
        "resolution above (Cartan-Eilenberg, Homological Algebra, Princeton 1956). "
        "Applying - (x)_A N and using P_v (x)_A N = N e_v (the vertex-v component of "
        "the left module N) collapses the complex to finite-dimensional vector spaces; "
        "the n-th homology is Tor_n.", heading=True)]
    dims = []
    for n in range(top + 1):
        tn = Tdim[n] if n < len(Tdim) else 0
        Dout = parts[n] if n < len(parts) else None          # d_{n+1}: T_{n+1} -> T_n
        r_out = _rank(Dout, dom) if Dout is not None else 0
        r_in = _rank(parts[n - 1], dom) if 0 <= n - 1 < len(parts) else 0
        res_dim = tn - r_out - r_in
        dims.append(res_dim)
        D = Dout if Dout is not None else lm.zeros(0, tn, dom)
        nrows, ncols = _mat_shape(D)
        events.append(ext_degree(
            degree=n, op="Tor", space_dim=tn, rank_here=r_out, rank_prev=r_in,
            result_dim=res_dim, D=D, nrows=nrows, ncols=ncols, dom=dom))
    from quiverlab.modules.tor import tor_dims
    engine_dims = tor_dims(A, M, N, top)
    if dims != engine_dims:
        raise QuiverlabError(
            "trace_tor: the shown Tor dimensions %s disagree with the engine's tor_dims "
            "%s for Tor(%s, %s) -- the worked steps drifted from the computation"
            % (dims, engine_dims, M.name, N.name))
    return events, dims


# --------------------------------------------------------------------------- #
# AR translate trace: the projective presentation + the D/Tr steps.
# --------------------------------------------------------------------------- #

def _presentation_events(res, dom, name="M"):
    """The minimal projective PRESENTATION Q_1 --d_1--> Q_0 --> ``name`` --> 0 as
    worked-step events (Q_0, the augmentation epsilon, Q_1, d_1) -- the input the
    transpose needs, shown WITHOUT the full-resolution narration / projective-dimension
    line (which the resolution step, if present, already gives). ``name`` is the module
    being presented (``M`` for tau, ``DM`` for the tau^- branch)."""
    out = [StepNote(
        "Projective presentation: Q_0 = %s covers %s via the augmentation epsilon, and "
        "Q_1 = %s covers the first syzygy, giving Q_1 --d_1--> Q_0 --> %s --> 0."
        % (_oplus_str(res.term(0)), name, _oplus_str(res.term(1)), name))]
    t0 = res.terms[0]
    out.append(ModuleTerm(
        degree=0, kind="projective", sym="P", summands=list(res.term(0)), dim=t0.dim,
        dimvec=_dv(t0.module.dimension_vector()) if t0.module is not None else {}))
    if res.dmats:
        D0 = res.differential(0)
        nr, nc = _mat_shape(D0)
        out.append(module_differential(
            degree=0, kind="projective", sym="P", symbol=r"\varepsilon",
            dom_summands=res.term(0), cod_summands=[], D=D0, nrows=nr, ncols=nc,
            dom=dom, cod_is_module=True, rank=_rank(D0, dom)))
    if res.length > 1:
        t1 = res.terms[1]
        out.append(ModuleTerm(
            degree=1, kind="projective", sym="P", summands=list(res.term(1)),
            dim=t1.dim,
            dimvec=_dv(t1.module.dimension_vector()) if t1.module is not None else {}))
        D1 = res.differential(1)
        nr, nc = _mat_shape(D1)
        if nc:
            out.append(module_differential(
                degree=1, kind="projective", sym="P", symbol="d_{1}",
                dom_summands=res.term(1), cod_summands=res.term(0), D=D1,
                nrows=nr, ncols=nc, dom=dom, cod_is_module=False, rank=_rank(D1, dom)))
    return out


def trace_tau(M, kind="tau", show_presentation=True):
    """Worked-step events for the AR translate tau M (kind="tau") or tau^- M
    (kind="tau_minus"): the projective presentation Q_1 -> Q_0 of M, the transpose
    Tr = coker(Hom(-,A)) shown as the transposed differential over A^op, the duality D,
    and the resulting translate's dimension vector (or ``tau = 0`` for a projective /
    ``tau^- = 0`` for an injective). ``show_presentation=False`` skips re-showing the
    presentation (used inside a combined report where the resolution step already did)."""
    dom = M.domain
    res = M.projective_resolution(1)
    d1 = res.differential(1) if res.length > 1 else lm.zeros(res.terms[0].dim, 0, dom)
    if kind == "tau":
        events = [StepNote(
            "Auslander-Reiten translate tau %s" % M.name,
            "From the minimal projective presentation Q_1 --d_1--> Q_0 --> M --> 0 the "
            "translate is tau M = D(Tr M), where Tr M = coker(Hom_A(d_1, A)) is the "
            "transpose (Assem-Simson-Skowronski IV.2, 2006). Applying the contravariant "
            "Hom_A(-, A) turns d_1 into the transposed map d_1^* : Hom(Q_0,A) -> "
            "Hom(Q_1,A) over A^op (at the level of the summand generators this is the "
            "transpose of d_1), Tr M is its cokernel over A^op, and D = Hom_k(-,k) sends "
            "it back to a right A-module.", heading=True)]
        if show_presentation:
            events += _presentation_events(res, dom)
        else:
            events.append(StepNote(
                "We reuse the presentation Q_1 --d_1--> Q_0 --> M --> 0 computed above."))
        events.append(StepNote(
            "transpose: Tr M = coker(Hom(d_1, A) : Hom(Q_0,A) -> Hom(Q_1,A))"))
        dstar = lm.transpose(d1)
        nr, nc = _mat_shape(dstar)
        if nc:
            rank_ds = _rank(dstar, dom)
            events.append(module_differential(
                degree=1, kind="projective", sym="P", symbol="d_{1}^{*}",
                dom_summands=res.term(0), cod_summands=res.term(1), D=dstar,
                nrows=nr, ncols=nc, dom=dom, cod_is_module=False, rank=rank_ds))
            pivds = lm.column_space_pivots(dstar, dom) if (dstar and dstar[0]) else []
            events.append(_reduction_note("d_1*", pivds, rank_ds, dom))
        # coker(d_1^*) = Tr M, then D: the two homework steps the old trace skipped.
        TrM = M.transpose()
        events.append(StepNote(
            "cokernel: Tr %s = coker(d_1^*) has dimension %d and dimension vector %s "
            "(a module over A^op)."
            % (M.name, TrM.dim, _fmt_dimvec(_dv(TrM.dimension_vector())))))
        events.append(StepNote(
            "duality: D = Hom_k(-,k) sends Tr %s back to a right A-module and preserves "
            "dimension vectors, so tau %s = D(Tr %s)." % (M.name, M.name, M.name)))
        t = M.tau()
    else:
        DM = M.dualize()
        DM.name = "D%s" % M.name
        events = [StepNote(
            "inverse Auslander-Reiten translate tau^- %s" % M.name,
            "Dually, tau^- M = Tr(D M): first dualize M to D M = Hom_k(M,k) (a module "
            "over A^op), take ITS minimal projective presentation, and transpose it back "
            "with Hom_{A^op}(-, A^op); the cokernel is a %s A-module "
            "(Assem-Simson-Skowronski IV.2, 2006). tau^-(injective) = 0."
            % _side_words(M)["side"], heading=True)]
        events.append(StepNote(
            "inverse AR translate tau^- M = Tr(D M): the transpose of the projective "
            "presentation of D M (tau^-(injective) = 0)."))
        resd = DM.projective_resolution(1)
        events += _presentation_events(resd, dom, name=DM.name)
        d1d = (resd.differential(1) if resd.length > 1
               else lm.zeros(resd.terms[0].dim, 0, dom))
        events.append(StepNote(
            "transpose: Tr(D%s) = coker(Hom(d_1, A^op) : Hom(Q_0,A^op) -> Hom(Q_1,A^op))"
            % M.name))
        dstar = lm.transpose(d1d)
        nr, nc = _mat_shape(dstar)
        if nc:
            rank_ds = _rank(dstar, dom)
            events.append(module_differential(
                degree=1, kind="projective", sym="P", symbol="d_{1}^{*}",
                dom_summands=resd.term(0), cod_summands=resd.term(1), D=dstar,
                nrows=nr, ncols=nc, dom=dom, cod_is_module=False, rank=rank_ds))
            pivds = lm.column_space_pivots(dstar, dom) if (dstar and dstar[0]) else []
            events.append(_reduction_note("d_1*", pivds, rank_ds, dom))
        t = M.tau_minus()
        events.append(StepNote(
            "cokernel: tau^- %s = Tr(D%s) = coker(d_1^*) directly (no further duality)."
            % (M.name, M.name)))
    sym = "tau" if kind == "tau" else "tau^-"
    if t.dim == 0:
        events.append(StepNote(
            "%s %s = 0 (M is %s, and the translate of such a module vanishes)."
            % (sym, M.name, "projective" if kind == "tau" else "injective")))
    else:
        events.append(StepNote(
            "%s %s has dimension vector %s (dim %d)."
            % (sym, M.name, _fmt_dimvec(_dv(t.dimension_vector())), t.dim)))
    return events, t


# --------------------------------------------------------------------------- #
# Krull-Schmidt decomposition trace: the End-ring split or locality certificate.
# --------------------------------------------------------------------------- #

def trace_decompose(M):
    """Worked-step events for the Krull-Schmidt decomposition of ``M``: the endomorphism
    ring End_A(M) = Hom_A(M,M), a Fitting split witness (a non-scalar endomorphism whose
    minimal polynomial factors coprimely, with its two kernels) OR the local-ring
    certificate that M is indecomposable -- honestly labelled when the certificate is
    refused within budget. Returns ``(events, result)`` where ``result`` is the
    ``decompose()`` list, or ``None`` when the engine refused loudly."""
    from quiverlab.modules.decompose import (
        _min_poly_coeffs, _factor_min_poly, _candidate_endomorphisms, _factoring_supported,
        _DEFAULT_BUDGET,
    )
    from quiverlab.modules.hom import hom_space
    from quiverlab.errors import QuiverlabError
    dom = M.domain
    H = hom_space(M, M)
    ev = [StepNote(
        "Krull-Schmidt decomposition of %s" % M.name,
        "By the Krull-Schmidt theorem M decomposes uniquely (up to isomorphism and "
        "order) as a direct sum of indecomposable summands. We test indecomposability "
        "through the endomorphism ring End_A(M) = Hom_A(M, M): if some endomorphism has "
        "a minimal polynomial that factors into coprime pieces, Fitting's lemma splits M "
        "as a direct sum of the corresponding kernels; if no such split exists and "
        "End_A(M) is local, M is indecomposable (Assem-Simson-Skowronski I, 2006).",
        heading=True),
        StepNote("The endomorphism ring End_A(%s) = Hom_A(%s, %s) has dimension %d over "
                 "%s." % (M.name, M.name, M.name, len(H), dom.name))]
    # search for a Fitting split witness (deterministic candidate order)
    witness = None
    if _factoring_supported(dom) and H:
        for phi in _candidate_endomorphisms(H, dom, _DEFAULT_BUDGET):
            m = _min_poly_coeffs(phi, dom)
            if len(m) <= 2:
                continue
            factors = _factor_min_poly(m, dom)
            if len(factors) >= 2:
                witness = (phi, m, factors)
                break
    if witness is not None:
        phi, m, factors = witness
        nr, nc = _mat_shape(phi)
        ev.append(StepNote(
            "A splitting endomorphism phi in End_A(%s) has minimal polynomial "
            "m(x) = %s, which factors into coprime pieces -- so Fitting's lemma applies."
            % (M.name, _poly_str(m, dom))))
        rank_phi = _rank(phi, dom)
        ev.append(module_differential(
            degree=0, kind="module", sym="Me", symbol=r"\varphi",
            dom_summands=_m_vertices(M), cod_summands=_m_vertices(M), D=phi,
            nrows=nr, ncols=nc, dom=dom, cod_is_module=True, dom_is_module=True,
            mod_name=M.name, rank=rank_phi))
        pivphi = lm.column_space_pivots(phi, dom) if (phi and phi[0]) else []
        ev.append(_reduction_note("phi", pivphi, rank_phi, dom))
        ev.append(StepNote(
            "Writing m = f * g with f, g coprime, M = ker f(phi) (+) ker g(phi) "
            "(a genuine module direct sum, since polynomials in the module map phi are "
            "again module maps). We recurse on the two summands."))
    try:
        result = M.decompose()
    except QuiverlabError as exc:
        ev.append(StepNote(
            "The decomposition could not be certified within budget (honest refusal): "
            "%s" % exc))
        return ev, None
    if witness is None:
        cert = ("dim End = 1, so End_A(M) = k*id is a field, hence local"
                if len(H) <= 1 else
                "no Fitting split exists and End_A(M) is local (certified by the exact "
                "endomorphism-ring criterion)")
        ev.append(StepNote(
            "No coprime split was found; %s. Therefore %s is indecomposable."
            % (cert, M.name)))
    trivial = (len(result) == 1 and result[0][1] == 1)
    if trivial:
        ev.append(StepNote(
            "Hence the Krull-Schmidt decomposition of %s is trivial: %s is itself "
            "indecomposable (%s)." % (M.name, M.name, _summand_desc(result[0][0]))))
    else:
        parts = " (+) ".join(
            "(%s)%s" % (_summand_desc(s), ("^%d" % k if k > 1 else ""))
            for s, k in result)
        ev.append(StepNote(
            "Krull-Schmidt decomposition: %s = %s." % (M.name, parts),
            "Each summand is certified indecomposable and the summands reassemble %s."
            % M.name))
    return ev, result


def _summand_desc(s):
    """A short description of an indecomposable summand: its dimension vector and its
    top (a hint at which indecomposable it is), e.g.
    ``dimension vector (2:1, 4:1), top S_2``."""
    return "dimension vector %s, top %s" % (
        _fmt_dimvec(_dv(s.dimension_vector())),
        _factor_stack_str(_dv(s.top().dimension_vector())))


# --------------------------------------------------------------------------- #
# The combined module report: rad/top/soc + resolution + tau + decompose (+ Ext).
# --------------------------------------------------------------------------- #

def trace_module_report(A, M, N=None, top=3, with_tau=True, with_decompose=True):
    """A full, homework-standard worked-steps report for the module ``M``: its radical,
    top and socle; the minimal projective resolution; (optionally) the AR translate
    tau M and the Krull-Schmidt decomposition; and, when a second module ``N`` is
    supplied, Ext^*(M, N). Returns the flat event list (feed straight to a renderer)."""
    ev = [StepNote(
        "We study the %s module %s (dimension %d, dimension vector %s) over the %s and compute "
        "its radical, top, socle, minimal projective resolution%s below. Each step gives "
        "the definition, the explicit matrices, and the justification."
        % (_side_words(M)["side"], M.name, M.dim,
           _fmt_dimvec(_dv(M.dimension_vector())), repr(A).splitlines()[0],
           (", Auslander-Reiten translate and Krull-Schmidt decomposition"
            if (with_tau or with_decompose) else "")))]
    ev += trace_radical(M)[0]
    ev += trace_top(M)[0]
    ev += trace_socle(M)[0]
    ev += trace_projective_resolution(M, top)[0]
    if N is not None:
        ev += trace_ext(A, M, N, top)[0]
    if with_tau:
        ev += trace_tau(M, "tau", show_presentation=False)[0]
    if with_decompose:
        ev += trace_decompose(M)[0]
    return ev


def _oplus_str(summands):
    """A plain-text direct sum of projectives ``P_1^2 (+) P_3`` from a vertex list with
    repetition (narrative use; the renderers have their own typeset version)."""
    if not summands:
        return "0"
    groups = []
    for v in summands:
        for g in groups:
            if g[0] == v:
                g[1] += 1
                break
        else:
            groups.append([v, 1])
    return " (+) ".join("P_%s" % v if c == 1 else "P_%s^%d" % (v, c) for v, c in groups)
