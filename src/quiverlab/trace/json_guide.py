"""Marco 2026-07-31 (ADDENDUM 2): the ``json_guide`` -- a per-computation index of HOW
to recover every object this example produced from the machine-readable ``result.json``.

``build_json_guide(results)`` returns a list of ``{"object", "path", "note"}`` recipes,
GENERATED from what was actually computed: for every result kind present it emits path
recipes with the REAL keys of THIS result (the Ext basis at a degree, the cup product
table, a resolution's differentials, ...), and NOTHING generic that does not exist here.

Two invariants make the guide unable to lie:

  * every recipe is emitted through ``_emit``, which first WALKS the relative key path
    against the block and skips the recipe unless it resolves -- so a path can never
    point at an absent object (the guide-cannot-lie test re-walks every path);
  * paths use ONLY dot and ``["key"]`` / ``[int]`` steps rooted at ``results.<kind>``,
    matching the ``result.json`` envelope shape (``{kind: block}``), so they are
    trivially parseable.

Accepts either runner's ``results`` shape via ``results_html.normalize``. Float-free.
"""
from quiverlab.trace.results_html import normalize


def build_json_guide(results):
    """The list of ``{object, path, note}`` recipes for the computed ``results``."""
    guide = []
    for kind, block in normalize(results):
        if not isinstance(block, dict) or block.get("error"):
            continue
        _recipes_for(kind, block, guide)
    return guide


# --------------------------------------------------------------------------- #
# Path model: a step is (key, dotted). ``dotted`` renders ``.key`` (a static field
# name); otherwise ``["key"]`` (a dynamic string key) or ``[int]`` (a list index).
# _walk applies obj[key] for every step, so an emitted path always resolves.
# --------------------------------------------------------------------------- #
def _render_path(kind, steps):
    out = "results.%s" % kind
    for key, dotted in steps:
        if dotted:
            out += ".%s" % key
        elif isinstance(key, int):
            out += "[%d]" % key
        else:
            out += '["%s"]' % key
    return out


def _walk(block, steps):
    obj = block
    for key, _dotted in steps:
        obj = obj[key]
    return obj


def _emit(guide, block, kind, obj_name, steps, note):
    """Append the recipe iff its path resolves against ``block`` (else skip silently)."""
    try:
        _walk(block, steps)
    except (KeyError, IndexError, TypeError):
        return
    guide.append({"object": obj_name, "path": _render_path(kind, steps), "note": note})


def _first_degree(d):
    """The smallest degree key of a ``{str(degree): ...}`` dict, or None."""
    if not isinstance(d, dict) or not d:
        return None
    try:
        return str(min(int(k) for k in d))
    except (ValueError, TypeError):
        return next(iter(d))


# --------------------------------------------------------------------------- #
# Per-kind recipes -- each key referenced is guarded by _emit's walk.
# --------------------------------------------------------------------------- #
def _hh_reps(kind, b, guide, cls_word, sided=False):
    """The shared basis_classes / chain_basis / differentials reps recipes. HH / Ext /
    Tor / cyclic carry the ``{str(degree): ...}`` payload (``sided=False``); the product
    blocks (cup/cap/bracket/connes_b) carry ``{side: {str(degree): ...}}``
    (``sided=True``), so the recipe descends through the side key first."""
    def _pick(field):
        d = b.get(field)
        if not isinstance(d, dict) or not d:
            return None
        if sided:
            side = next((s for s in ("coh", "hom") if d.get(s)), None)
            deg = _first_degree(d[side]) if side is not None else None
            if side is None or deg is None:
                return None
            return [(field, True), (side, False), (deg, False)], ("side %s, " % side), deg
        deg = _first_degree(d)
        if deg is None:
            return None
        return [(field, True), (deg, False)], "", deg

    p = _pick("basis_classes")
    if p:
        steps, where, deg = p
        _emit(guide, b, kind, "%s in degree %s" % (cls_word, deg), steps,
              "%slist; entry j is class j -- .terms is the labelled term-sum, .vector "
              "is the sparse coordinate vector [[index, coeff], ...] over the ordered "
              "basis in the matching chain_basis entry. Sibling keys are the other "
              "degrees." % where)
    p = _pick("chain_basis")
    if p:
        steps, _where, deg = p
        _emit(guide, b, kind, "ordered basis of degree %s" % deg, steps,
              "the ordered basis element labels the .vector indices point into.")
    p = _pick("differentials")
    if p:
        steps, _where, deg = p
        _emit(guide, b, kind, "annihilating differential in degree %s" % deg, steps,
              "the (co)boundary whose kernel/quotient gives this degree; .rows is the "
              "exact matrix (or elided with a rebuild note past the recorder cap).")


def _recipes_for(kind, b, guide):
    if kind in ("hh_cohomology", "hh_homology"):
        word = "cohomology" if kind == "hh_cohomology" else "homology"
        _emit(guide, b, kind, "Hochschild %s dimensions" % word, [("dims", True)],
              "list; dim HH in degree n is dims[n].")
        _emit(guide, b, kind, "inner-derivation dimensions", [("inner_dims", True)],
              "rank of the coboundary per degree (HH^1's inner derivations).")
        _hh_reps(kind, b, guide, "Hochschild %s classes" % word)
        return
    if kind == "cyclic_homology":
        _emit(guide, b, kind, "cyclic homology dimensions", [("dims", True)],
              "list; dim HC_n is dims[n].")
        _hh_reps(kind, b, guide, "cyclic homology classes")
        cs = b.get("column_structure")
        m = _first_degree(cs)
        if m is not None:
            _emit(guide, b, kind, "total-complex columns of Tot_%s" % m,
                  [("column_structure", True), (m, False)],
                  "the Tot_n = C_n (+) C_{n-2} (+) ... column layout (degree, offset, "
                  "dim) the coordinate vectors are sliced by.")
        return
    if kind in ("ext", "tor"):
        _emit(guide, b, kind, "%s dimensions" % kind.title(), [("dims", True)],
              "list; dim in degree n is dims[n].")
        _emit(guide, b, kind, "the second module N", [("target", True), ("dimvec", True)],
              "N's dimension vector (the module Ext/Tor was taken against).")
        _hh_reps(kind, b, guide, "%s classes" % kind.title())
        if kind == "ext":
            interp = b.get("interpretation")
            seqs = interp.get("sequences") if isinstance(interp, dict) else None
            n = _first_degree(seqs)
            if n is not None:
                _emit(guide, b, kind, "Yoneda exact sequences in degree %s" % n,
                      [("interpretation", True), ("sequences", True), (n, False)],
                      "each Ext^n class as the constructed + certified n-fold exact "
                      "sequence 0 -> N -> Q -> ... -> M -> 0 (modules, maps, facts).")
        return
    if kind in ("cup", "cap", "bracket"):
        op = {"cup": "∪", "cap": "∩", "bracket": "[-,-]"}[kind]
        _emit(guide, b, kind, "%s product tables" % kind, [("tables", True)],
              "list of tables. For α_i^p %s α_j^q find the table whose degrees == [p, q]; "
              "constants[k][i][j] is the coefficient of the k-th output class; "
              "dims == [dl, dr, dout]." % op)
        _hh_reps(kind, b, guide, "the operand/output classes", sided=True)
        return
    if kind == "connes_b":
        mats = b.get("matrices")
        n = _first_degree(mats)
        if n is not None:
            _emit(guide, b, kind, "induced Connes differential B_%s" % n,
                  [("matrices", True), (n, False)],
                  "the matrix of B_n: HH_n -> HH_{n+1} on the recorded homology bases "
                  "(rows index HH_{n+1}, columns HH_n).")
            _emit(guide, b, kind, "rank of B_%s" % n, [("ranks", True), (n, False)],
                  "the induced rank of B_n.")
        _hh_reps(kind, b, guide, "the cycle classes", sided=True)
        return
    if kind in ("projective_resolution", "injective_resolution"):
        _emit(guide, b, kind, "resolution terms", [("terms", True)],
              "list; terms[n] is the degree-n term's dimension vector {vertex: mult}.")
        _emit(guide, b, kind, "the degree-0 differential", [("differentials", True),
              (0, False), ("matrix", True)],
              "differentials[n].matrix is the exact matrix of d_n (rows: target basis, "
              "columns: source basis); an elided body carries a rebuild note.")
        _emit(guide, b, kind, "ordered basis of term 0", [("term_basis", True), (0, False)],
              "term_basis[n] lists the ordered path basis of term n -- the index order "
              "the differential grids use for their columns (proj) / rows (inj).")
        for key in ("pd", "injective_dimension"):
            _emit(guide, b, kind, "the homological dimension", [(key, True)],
                  "the projective/injective dimension (null = beyond the probed length).")
        return
    # Single-object kinds: point at whatever key this block actually carries.
    _SIMPLE = {
        "cartan": ("matrix", "the Cartan matrix C = (dim e_i A e_j)."),
        "coxeter_polynomial": ("latex", "the Coxeter polynomial (LaTeX source)."),
        "global_dimension": ("text", "the global dimension, stated in prose."),
        "center": ("dim", "dim Z(A), the dimension of the centre."),
        "dimension": ("value", "dim_k A."),
        "dimension_vector": ("latex", "the module's dimension vector (LaTeX)."),
        "projective_dimension": ("value", "pd M (null = beyond the probed length)."),
        "injective_dimension": ("value", "id M (null = beyond the probed length)."),
    }
    if kind in _SIMPLE:
        field, note = _SIMPLE[kind]
        _emit(guide, b, kind, kind.replace("_", " "), [(field, True)], note)
        return
    if kind == "rad_top_soc":
        for field in ("radical", "top", "socle"):
            _emit(guide, b, kind, "%s of M" % field, [(field, True)],
                  "a full representation {dims: dimension vector, maps: per-arrow "
                  "matrices} of rad/top/soc M.")
        return
    if kind in ("tau", "tau_minus"):
        _emit(guide, b, kind, "the AR translate", [("repr", True)],
              "the translate's full representation (dims + per-arrow maps); "
              "targets[] carries the second module's translate when one was named.")
        return
    if kind == "decompose":
        _emit(guide, b, kind, "Krull-Schmidt summands", [("summands", True)],
              "list of indecomposable summands, each with multiplicity, dim_vector, and "
              "maps (or a standard-indecomposable name).")
        return
