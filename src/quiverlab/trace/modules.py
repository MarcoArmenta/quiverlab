"""Module-computation worked-step builders (Plan 30, Part C).

These are the sanctioned ``quiverlab.trace`` hooks for the module surface: given a
module computation they RUN it through the public surface and emit the typed
step events (``ModuleTerm`` / ``ModuleDifferential`` / ``ExtDegree`` / ``StepNote``,
see ``trace.events``) that the renderers turn into an exhaustive, replay-by-hand
pdf/tex/html bundle.

Design (mirrors how ``hochschild.bar`` hooks the HH engine):
  * a *projective resolution* trace reads the PUBLIC ``ProjectiveResolution``
    object -- its ``.term(n)`` summand vertices and ``.differential(n)`` matrices;
  * an *Ext* trace re-runs the same ``ext_dims`` pipeline (minimal resolution +
    ``hom_space`` + the Hom-collapse ``delta`` matrices), emitting the per-degree
    rank bookkeeping; a cross-check test pins its dims to ``ext_dims``;
  * a ``tau`` trace shows the projective presentation P_1 -> P_0 and narrates the
    D/Tr steps of the AR translate.

Matrix bodies are elided past ``recorder.MATRIX_ELISION_CELLS`` (shape + rank kept)
-- every computational step still appears, never silently omitted. Float-free:
all data are ints and domain-element strings.
"""
from quiverlab.modules import linalg_mod as lm
from quiverlab.trace.events import ModuleTerm, StepNote
from quiverlab.trace.recorder import module_differential, ext_degree


def _dv(dimvec):
    """A dimension vector as a str-keyed, vertex-sorted dict (golden-stable)."""
    return {str(v): int(n) for v, n in sorted(dimvec.items(), key=lambda kv: str(kv[0]))}


def _mat_shape(D):
    return (len(D), len(D[0]) if (D and D[0]) else 0)


def _rank(D, dom):
    return lm.mat_rank(D, dom) if (D and D[0]) else 0


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
    epsilon: Q_0 -> M, then d_n: Q_n -> Q_{n-1}) rendered as a matrix. Returns
    ``(events, res)`` where ``res`` is the public ``ProjectiveResolution``."""
    res = M.projective_resolution(top)
    dom = M.domain
    events = [StepNote(
        "Minimal projective resolution of %s (top generators chosen modulo rad M)"
        % M.name,
        "projective cover of M: top M has composition factors %s"
        % _fmt_dimvec(_dv(M.top().dimension_vector())))]
    n_terms = res.length
    for n in range(n_terms):
        t = res.terms[n]
        verts = res.term(n)
        dimvec = _dv(t.module.dimension_vector()) if t.module is not None else {}
        events.append(ModuleTerm(degree=n, kind="projective", sym="P",
                                 summands=list(verts), dim=t.dim, dimvec=dimvec))
        # the differential out of Q_n, if recorded: d_0 = epsilon (Q_0 -> M),
        # d_n = (Q_n -> Q_{n-1}).
        if n < len(res.dmats):
            D = res.differential(n)
            nrows, ncols = _mat_shape(D)
            cod_is_module = (n == 0)
            symbol = r"\varepsilon" if n == 0 else "d_{%d}" % n
            events.append(module_differential(
                degree=n, kind="projective", sym="P", symbol=symbol,
                dom_summands=verts,
                cod_summands=[] if cod_is_module else res.term(n - 1),
                D=D, nrows=nrows, ncols=ncols, dom=dom,
                cod_is_module=cod_is_module, rank=_rank(D, dom)))
    pd = res.pd()
    events.append(StepNote(
        "projective dimension pd(M) = %s" % ("infinite (unresolved within top)"
                                             if pd is None else pd)))
    return events, res


def trace_injective_resolution(M, top):
    """Worked-step events for the minimal injective coresolution of ``M`` to length
    ``top``: each term E^n = (+) I_v with its dimension vector, plus a narration of
    the duality E^n = D(P_n over A^op). The injective differentials are the k-duals
    of the projective differentials of DM (not re-rendered here). Returns
    ``(events, res)`` with the public ``InjectiveResolution``."""
    res = M.injective_resolution(top)
    events = [StepNote(
        "Minimal injective coresolution of %s" % M.name,
        "E^n = D(P_n), P_* the minimal projective resolution of D%s over A^op"
        % M.name)]
    for n in range(res.length):
        t = res.terms[n]
        verts = res.term(n)
        dim = t.dim if t is not None else 0
        dimvec = _dv(t.dimension_vector()) if t is not None else {}
        events.append(ModuleTerm(degree=n, kind="injective", sym="I",
                                 summands=list(verts), dim=dim, dimvec=dimvec))
    idim = res.injective_dimension()
    events.append(StepNote(
        "injective dimension id(M) = %s" % ("infinite (unresolved within top)"
                                            if idim is None else idim)))
    return events, res


# --------------------------------------------------------------------------- #
# Ext trace: the Hom-collapse matrices + per-degree rank bookkeeping.
# --------------------------------------------------------------------------- #

def trace_ext(A, M, N, top):
    """Worked-step events for Ext^0..top_A(M, N): the minimal resolution of M, the
    Hom(Q_n, N) space dimensions, the connecting maps delta^n rendered as matrices,
    and Ext^n = dim Hom(Q_n,N) - rank(delta^n) - rank(delta^{n-1}). Re-runs the
    ``ext_dims`` pipeline (a cross-check test pins the emitted dims to it). Returns
    ``(events, dims)``."""
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
        "Ext^*(%s, %s) = H^*(Hom(Q_*, %s)) over the minimal resolution Q_* of %s"
        % (M.name, N.name, N.name, M.name))]
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
    return events, dims


# --------------------------------------------------------------------------- #
# AR translate trace: the projective presentation + the D/Tr steps.
# --------------------------------------------------------------------------- #

def trace_tau(M, kind="tau"):
    """Worked-step events for the AR translate tau M (kind="tau") or tau^- M
    (kind="tau_minus"): the projective presentation Q_1 -> Q_0 of M, the transpose
    Tr / duality D narration, and the resulting translate's dimension vector (or
    ``tau = 0`` for a projective, ``tau^- = 0`` for an injective)."""
    events, _ = trace_projective_resolution(M, 1)
    if kind == "tau":
        events.append(StepNote(
            "transpose: Tr M = coker(Hom(d_1, A) : Hom(Q_0,A) -> Hom(Q_1,A))",
            "AR translate tau M = D(Tr M)  (tau(projective) = 0)"))
        t = M.tau()
    else:
        events.append(StepNote(
            "inverse AR translate tau^- M = Tr(D M)",
            "transpose of the projective presentation of D M "
            "(tau^-(injective) = 0)"))
        t = M.tau_minus()
    sym = "tau" if kind == "tau" else "tau^-"
    if t.dim == 0:
        events.append(StepNote("%s %s = 0" % (sym, M.name)))
    else:
        events.append(StepNote(
            "%s %s has dimension vector %s (dim %d)"
            % (sym, M.name, _fmt_dimvec(_dv(t.dimension_vector())), t.dim)))
    return events, t


def _fmt_dimvec(dv):
    """A plain-text dimension vector ``(1:2, 3:1)`` (shared narrative formatting)."""
    return "(" + ", ".join("%s:%s" % (k, v) for k, v in dv.items()) + ")"
