"""Render the COMPUTED RESULT BLOCKS into the worked-steps HTML report.

Marco (2026-07-29): "the idea of the report html is that it has saved everything
that we have computed with the gui, and they can keep that report." The report used
to carry only the worked STEPS of one traced computation; every other answer the
page displayed (the Cartan matrix, the Coxeter polynomial, rad/top/soc, the AR
translates, Ext/Tor tables, the resolutions...) lived only in the browser and was
lost when the tab closed.

This module renders exactly the blocks the two runners produce
(``quiverlab.hpc.spec._dispatch`` / ``_dispatch_module`` and their Pyodide twin in
``docs/gui/runner.py``), in the SAME shapes the GUI shows them, so the saved report
is a faithful record of the session.

Presentation rules carried over from the GUI (all Marco, 2026-07-29):
  * matrices are shown COMPLETE -- no scrollbar, no clipping (a wide one is typeset
    a size down by :func:`quiverlab.trace.render_html._math`);
  * an arrow acting as the exact ZERO map is named, not printed as a zero block;
  * a differential equal to one already shown is REFERENCED, not repeated.

Float-free: every number is copied from block data (ints / exact strings).
"""
from quiverlab.trace.render_html import (
    _dims_table, _esc, _math, _math_inline, matrix_grid)

# Human headings per compute kind. A kind missing here still renders (the key is
# shown verbatim), so a newly added invariant degrades to an honest label rather
# than vanishing from the report.
_HEADINGS = {
    "hh_cohomology": "Hochschild cohomology",
    "hh_homology": "Hochschild homology",
    "cyclic_homology": "Cyclic homology",
    # Plan 35 HH product surface -- the gui.js PRODUCT_TITLE i18n titles.
    "cup": "Cup product tables",
    "cap": "Cap product tables",
    "bracket": "Gerstenhaber bracket tables",
    "connes_b": "Connes differentials",
    "cartan": "Cartan matrix",
    "coxeter_polynomial": "Coxeter polynomial",
    "global_dimension": "Global dimension",
    "center": "Centre",
    "dimension": "Dimension",
    "dimension_vector": "Dimension vector of M",
    "rad_top_soc": "Radical, top and socle of M",
    "tau": "AR translate τM",
    "tau_minus": "Inverse AR translate τ⁻M",
    "decompose": "Krull–Schmidt decomposition of M",
    "ext": "Ext",
    "tor": "Tor",
    "projective_resolution": "Projective resolution of M",
    "injective_resolution": "Injective resolution of M",
    "projective_dimension": "Projective dimension of M",
    "injective_dimension": "Injective dimension of M",
}

_TARGET_ROLE = {"ext_target": "the Ext target", "tor_target": "the Tor target"}

_TAU_NOTES = {"mod.tau_additive": "τ computed summand-wise (τ is additive)"}


def normalize(results):
    """``results`` in either runner's shape -> ``[(kind, block), ...]``.

    The server/HPC tier returns ``{kind: block}``; the Pyodide runner accumulates a
    LIST of blocks each tagged with its ``invariant`` string (``"ext:0..4"``), and a
    failed computation contributes ``{"invariant": ..., "error": {...}}``. Both are
    accepted so one renderer serves both tiers. Anything else yields ``[]``."""
    out = []
    if isinstance(results, dict):
        for kind, block in results.items():
            if isinstance(block, dict):
                out.append((str(kind), block))
        return out
    if isinstance(results, (list, tuple)):
        for block in results:
            if not isinstance(block, dict):
                continue
            spec = block.get("invariant") or block.get("kind") or ""
            out.append((str(spec).split(":")[0], block))
    return out


def results_section(results):
    """HTML chunks for every computed block, or ``[]`` when there is nothing."""
    items = normalize(results)
    if not items:
        return []
    out = []
    for kind, block in items:
        heading = _HEADINGS.get(kind, kind.replace("_", " "))
        out.append("<h3>%s</h3>" % _esc(heading))
        err = block.get("error")
        if err:
            out.append("<p class='ql-note'>not computed — %s</p>"
                       % _esc(_error_text(err)))
            continue
        out.extend(_block_html(kind, block))
        out.extend(_citations_html(block))
    return out


# --------------------------------------------------------------------------- #
# Per-kind rendering
# --------------------------------------------------------------------------- #

def _block_html(kind, b):
    if kind in ("hh_cohomology", "hh_homology"):
        sup = kind == "hh_cohomology"
        chunks = [_dims_table("dim HH%sn" % ("^" if sup else "_"), b.get("dims") or [])]
        if b.get("engine"):
            chunks.append("<p class='ql-note'>engine: %s</p>" % _esc(str(b["engine"])))
        return chunks
    if kind == "cyclic_homology":
        # HC is a homology-style subscript table HC_n (Plan-35 follow-up), rendered
        # exactly like the HH dims tables above.
        chunks = [_dims_table("dim HC_n", b.get("dims") or [])]
        if b.get("engine"):
            chunks.append("<p class='ql-note'>engine: %s</p>" % _esc(str(b["engine"])))
        return chunks
    if kind == "cartan":
        if b.get("matrix"):
            return [matrix_grid(b["matrix"], label="C")]
        return [_math("C = " + b["latex"])] if b.get("latex") else []
    if kind == "coxeter_polynomial":
        return [_math(r"\chi(t) = " + b["latex"])] if b.get("latex") else []
    if kind == "global_dimension":
        return ["<p>%s</p>" % _esc(str(b.get("text", "")))]
    if kind == "center":
        return [_math(r"\dim Z(A) = %s" % _num(b.get("dim")))]
    if kind == "dimension":
        return [_math(r"\dim_k A = %s" % _num(b.get("value")))]
    if kind == "dimension_vector":
        return [_math(b["latex"])] if b.get("latex") else []
    if kind == "rad_top_soc":
        return _rad_top_soc_html(b)
    if kind in ("tau", "tau_minus"):
        return _tau_html(kind, b)
    if kind == "decompose":
        return _decompose_html(b)
    if kind in ("ext", "tor"):
        op = "Ext^" if kind == "ext" else "Tor_"
        chunks = []
        target = (b.get("target") or {}).get("dimvec")
        if target:
            chunks.append("<p>against N, dimension vector %s.</p>" % _esc(_dv(target)))
        chunks.append(_dims_table("dim %sn" % op, b.get("dims") or []))
        # Plan 35 wave 3a: the per-degree EXPLICIT REPRESENTATIVES (ordered basis ->
        # classes -> differential + verification), when the block carries them.
        from quiverlab.trace.render_html import module_reps_sections
        secs = module_reps_sections(b.get("basis_classes"), b.get("chain_basis"),
                                    b.get("differentials"), kind,
                                    anchor_prefix="cr")
        if secs:
            chunks.append("<p><i>Explicit representatives by degree — each class as a "
                          "term-sum and a coordinate vector over the ordered basis, with "
                          "the differential that annihilates it.</i></p>")
            chunks.extend(secs)
        return chunks
    if kind in ("cup", "cap", "bracket"):
        return _product_tables_html(kind, b)
    if kind == "connes_b":
        return _connes_b_html(b)
    if kind in ("projective_resolution", "injective_resolution"):
        return _resolution_html(kind, b)
    if kind in ("projective_dimension", "injective_dimension"):
        chunks = []
        if b.get("latex"):
            chunks.append(_math(b["latex"]))
        if b.get("note"):
            chunks.append("<p class='ql-note'>%s</p>" % _esc(str(b["note"])))
        return chunks
    # An unknown kind still leaves a trace of what was asked for.
    return ["<p class='ql-note'>computed; see the JSON record for its data.</p>"]


def _rad_top_soc_html(b):
    trio = [("rad M", b.get("radical")), ("top M", b.get("top")),
            ("soc M", b.get("socle"))]
    if any(v is None or v.get("dims") is None for _, v in trio):
        return ["<p class='ql-note'>this result predates the full-representation "
                "format; recompute it to record the per-arrow matrices.</p>"]
    rows = ["<tr><th></th><th>dim vector</th></tr>"]
    for label, view in trio:
        rows.append("<tr><th>%s</th><td>%s</td></tr>"
                    % (_esc(label), _esc(_dv(view.get("dims")))))
    out = ['<table class="ql-table">%s</table>' % "".join(rows)]
    if any(v.get("display_only") for _, v in trio):
        out.append("<p class='ql-note'>display only — entries lie outside the "
                   "integer/fraction input grammar (e.g. GF(p^n) elements).</p>")
    for label, view in trio:
        out.extend(_maps_html(label, view))
    return out


def _tau_html(kind, b):
    name = "M"
    out = _translate_html(kind, b, name)
    for t in b.get("targets") or []:
        role = _TARGET_ROLE.get(t.get("role"))
        out.append("<p>and for N%s:</p>" % (", %s" % _esc(role) if role else ""))
        out.extend(_translate_html(kind, t, "N"))
    return out


def _translate_html(kind, t, name):
    sym = ("τ" if kind == "tau" else "τ⁻") + name
    if t.get("error"):
        return ["<p class='ql-note'>%s is unavailable: %s</p>"
                % (_esc(sym), _esc(str(t["error"])))]
    out = []
    if t.get("latex"):
        out.append(_math(t["latex"]))
    rep = t.get("repr")
    if rep:
        out.extend(_maps_html(sym, rep))
    out.extend(_input_certificate_html(t, name))
    return out


def _input_certificate_html(t, name):
    if t.get("indecomposable") is True:
        return ["<p class='ql-note'>input %s is indecomposable.</p>" % _esc(name)]
    if t.get("decomposition"):
        parts = "  ⊕  ".join(
            _dv(s.get("dim_vector")) + ("^%s" % _num(s.get("multiplicity"))
                                        if (s.get("multiplicity") or 1) > 1 else "")
            for s in t["decomposition"])
        note = _TAU_NOTES.get(t.get("note_key"), "")
        return ["<p class='ql-note'>input %s ≅ %s%s</p>"
                % (_esc(name), _esc(parts), " — " + _esc(note) if note else "")]
    return []


_STD_SYM = {"simple": "S", "projective": "P", "injective": "I"}


def _summand_name(s, i):
    """``S_2`` / ``P_1`` / ``I_3`` for a summand recognised as a standard
    indecomposable, else the positional ``M_i`` (Marco 2026-07-29)."""
    std = s.get("standard")
    if isinstance(std, dict) and std.get("kind") in _STD_SYM:
        return "%s_%s" % (_STD_SYM[std["kind"]], std.get("vertex"))
    return "M_%d" % i


def _decompose_html(b):
    summands = b.get("summands") or []
    rows = ["<tr><th>summand</th><th>multiplicity</th><th>dim vector</th></tr>"]
    for i, s in enumerate(summands):
        rows.append("<tr><td>%s</td><td>%s</td><td>%s</td></tr>"
                    % (_esc(_summand_name(s, i + 1)), _num(s.get("multiplicity")),
                       _esc(_dv(s.get("dim_vector")))))
    out = ["<p>%s indecomposable summand(s):</p>"
           % _num(b.get("iso_classes", len(summands))),
           '<table class="ql-table">%s</table>' % "".join(rows)]
    # A summand named S_v / P_v / I_v needs no matrices -- the name IS the module.
    # Every other one carries its full action, since its dimension vector does not
    # determine it (Marco 2026-07-29).
    for i, s in enumerate(summands):
        if s.get("standard") or not s.get("maps"):
            continue
        out.extend(_maps_html(_summand_name(s, i + 1), s))
    if any(s.get("display_only") for s in summands):
        out.append("<p class='ql-note'>display only — entries lie outside the "
                   "integer/fraction input grammar (e.g. GF(p^n) elements).</p>")
    return out


def _resolution_html(kind, b):
    proj = kind == "projective_resolution"
    summ = b.get("summands") or []
    terms = b.get("terms") or []
    # term | (+)-decomposition, the table the GUI shows (Plan 30, Marco #3).
    rows = ["<tr><th>term</th><th>&#8853;-decomposition</th></tr>"]
    for n in range(len(terms)):
        tex = summ[n] if n < len(summ) and summ[n] is not None else "0"
        rows.append("<tr><td>%d</td><td>%s</td></tr>" % (n, _math_inline(tex)))
    out = ['<table class="ql-table">%s</table>' % "".join(rows)]
    d = b.get("pd") if proj else b.get("injective_dimension")
    out.append("<p>%s = %s</p>" % ("pd M" if proj else "id M",
                                   _num(d) if d is not None
                                   else "&#8734; (beyond the probed length)"))
    out.extend(_term_basis_html(b, proj))
    out.extend(_differentials_html(b, proj))
    return out


# ordered term-basis entries shown inline before a machine-record pointer.
_TERM_BASIS_DISPLAY = 64


def _term_basis_html(b, proj):
    """The ordered k-basis of each resolution term (Plan 35 UNIT 2), as a numbered
    (1-based) list per term -- the SAME index order the differential grids use for
    their columns (projective) / rows (injective), so a matrix entry can be traced to
    the basis vector it acts on. Omitted when the block carries no ``term_basis`` (a
    structure-constants algebra, or an older cached result) -- tolerance."""
    tb = b.get("term_basis")
    if not tb:
        return []
    out = ["<p><i>Ordered basis of each resolution term — the 1-based index order the "
           "differential grids below use for their %s.</i></p>"
           % ("columns" if proj else "rows")]
    for n, labels in enumerate(tb):
        name = ("Q_{%d}" % n) if proj else ("E^{%d}" % n)
        if not labels:
            out.append("<p>%s = 0.</p>" % _math_inline(name))
            continue
        items = "".join("<li>%s</li>" % _esc(str(x))
                        for x in labels[:_TERM_BASIS_DISPLAY])
        out.append("<p>%s basis:</p><ol class='ql-enum'>%s</ol>"
                   % (_math_inline(name), items))
        if len(labels) > _TERM_BASIS_DISPLAY:
            out.append("<p class='ql-note'>… %d more (full list in the machine "
                       "record).</p>" % (len(labels) - _TERM_BASIS_DISPLAY))
    return out


def _differentials_html(b, proj):
    diffs = b.get("differentials") or []
    if not diffs:
        return []
    out = ["<p><i>Differentials (rows: target basis, columns: source basis; "
           "%s).</i></p>" % ("d<sub>0</sub> = &#949;: Q<sub>0</sub> &#8594; M" if proj
                             else "d<sup>0</sup> = &#953;: M &#8594; E<sup>0</sup>")]
    seen = {}
    for n, d in enumerate(diffs):
        sym = ("d_{%d}" if proj else "d^{%d}") % n
        label = ("d_%d" if proj else "d^%d") % n
        matrix = d.get("matrix")
        if d.get("elided") or matrix is None:
            out.append("<p class='ql-note'>%s: %s&#215;%s matrix (body not recorded "
                       "— it exceeded the recorder's memory backstop).</p>"
                       % (_esc(label), _num(d.get("rows")), _num(d.get("cols"))))
            continue
        key = _matrix_key(matrix)
        if key in seen:
            out.append(_math("%s = %s" % (sym, seen[key])))
            out.append("<p class='ql-note'>(the same matrix as above; not "
                       "repeated)</p>")
            continue
        seen[key] = sym
        out.append(matrix_grid(matrix, label=sym))
    return out


# --------------------------------------------------------------------------- #
# Plan 35 HH product surface: cup / cap / bracket tables + the Connes differential.
# The block shapes are the frozen result objects' ``.blocks()`` (see
# quiverlab.hochschild.products): cup/cap/bracket carry ``tables`` (each with
# ``degrees`` / ``out_degree`` / ``dims=[dl,dr,dout]`` / ``constants[k][i][j]``
# exact strings); connes_b carries per-n ``matrices`` + ``ranks``. Rendered the way
# the GUI shows them (webapp/static/gui/gui.js), sharing the equation builder with
# the worked-steps chapter so the report and the chapter never drift.
# --------------------------------------------------------------------------- #

def _product_heading(kind, degrees, out_degree):
    """The map label above one table: ``HH^p ∪ HH^q → HH^{p+q}`` (cup),
    ``HH^p ∩ HH_n → HH_{n-p}`` (cap), ``[HH^p, HH^q] → HH^{p+q-1}`` (bracket)."""
    p = degrees[0] if degrees else 0
    q = degrees[1] if len(degrees) > 1 else 0
    if kind == "cap":
        return r"HH^{%s} \cap HH_{%s} \to HH_{%s}" % (p, q, out_degree)
    if kind == "bracket":
        return r"[HH^{%s}, HH^{%s}] \to HH^{%s}" % (p, q, out_degree)
    return r"HH^{%s} \cup HH^{%s} \to HH^{%s}" % (p, q, out_degree)


def _product_tables_html(kind, b):
    """cup / cap / bracket: the notation legend, the per-degree EXPLICIT
    REPRESENTATIVES (ordered basis -> classes -> differential + verification, Plan 35
    UNIT 2), then one map heading plus its nonzero-product equation lines per bidegree
    (an all-vanishing table states so), the bracket's served-window note, and the
    engine provenance. The equation lines come from the SAME builder the worked-steps
    chapter uses (``quiverlab.trace.products.equation_lines``); the degree sections
    come from the SAME renderer the chapter uses
    (``quiverlab.trace.render_html.product_degree_sections``) -- one implementation,
    no drift. ``b["differentials"]`` here is the product ``{side:{degree:...}}`` shape,
    read ONLY inside this kind-scoped function (the module-resolution block ships a
    LIST under the same key -- never read it shape-blind)."""
    from quiverlab.trace.products import equation_lines, notation_legend
    from quiverlab.trace.render_html import (
        product_degree_sections, product_table_reference)
    out = ["<p class='ql-note'>%s</p>"
           % _esc(notation_legend(kind, "", b.get("basis")))]
    secs = product_degree_sections(b.get("basis_classes"), b.get("chain_basis"),
                                   b.get("differentials"), anchor_prefix="cr-" + kind)
    have_reps = bool(secs)
    if have_reps:
        out.append("<p><b>Explicit representatives by degree</b></p>")
        out.extend(secs)
        out.append("<p><b>Structure-constant tables</b> (in the explicit classes "
                   "above; each table links to its operands' and output's degree "
                   "sections):</p>")
    for t in (b.get("tables") or []):
        degrees = list(t.get("degrees") or [])
        out_degree = t.get("out_degree")
        dims = list(t.get("dims") or [0, 0, 0])
        constants = t.get("constants") or []
        out.append('<p class="ql-mlabel">%s</p>'
                   % _math_inline(_product_heading(kind, degrees, out_degree)))
        if have_reps and degrees:               # MINOR 1: link to the degree sections
            out.append(product_table_reference(kind, degrees, out_degree, "cr-" + kind))
        lines = equation_lines(kind, degrees, out_degree, dims, constants)
        if lines:
            out.extend(_math(line) for line in lines)
        else:
            out.append("<p class='ql-note'>every product in this bidegree "
                       "vanishes.</p>")
    if kind == "bracket" and b.get("window") is not None:
        out.append("<p class='ql-note'>bracket structure constants served to the "
                   "degree window %s (bar-transport bound).</p>"
                   % _num(b.get("window")))
    if b.get("engine"):
        out.append("<p class='ql-note'>engine: %s</p>" % _esc(str(b["engine"])))
    return out


def _connes_b_html(b):
    """connes_b: the legend, the per-degree explicit homology-cycle representatives
    (Plan 35 UNIT 2), then one induced Connes differential grid
    ``B_n : HH_n → HH_{n+1}`` per degree with its induced rank, and the engine
    provenance. ``b["differentials"]`` is the product ``{side:{degree:...}}`` shape,
    read ONLY here (kind-scoped)."""
    from quiverlab.trace.products import notation_legend
    from quiverlab.trace.render_html import (
        product_degree_sections, product_table_reference)
    out = ["<p class='ql-note'>%s</p>"
           % _esc(notation_legend("connes_b", "", None))]
    secs = product_degree_sections(b.get("basis_classes"), b.get("chain_basis"),
                                   b.get("differentials"), anchor_prefix="cr-connes_b")
    have_reps = bool(secs)
    if have_reps:
        out.append("<p><b>Explicit representatives by degree</b></p>")
        out.extend(secs)
        out.append("<p><b>Induced Connes differentials</b> (in the explicit cycle "
                   "classes above; each links to its source/target degree "
                   "sections):</p>")
    matrices = b.get("matrices") or {}
    ranks = b.get("ranks") or {}
    for key in sorted(matrices, key=lambda s: int(s)):
        n = int(key)
        out.append('<p class="ql-mlabel">%s</p>'
                   % _math_inline(r"B_{%d} : HH_{%d} \to HH_{%d}" % (n, n, n + 1)))
        if have_reps:                           # MINOR 1: link to the degree sections
            out.append(product_table_reference("connes_b", (n,), n, "cr-connes_b"))
        out.append(matrix_grid(matrices[key], label="B_{%d}" % n))
        out.append("<p class='ql-note'>induced rank B_%d = %s</p>"
                   % (n, _num(ranks.get(key))))
    if b.get("engine"):
        out.append("<p class='ql-note'>engine: %s</p>" % _esc(str(b["engine"])))
    return out


# --------------------------------------------------------------------------- #
# Small shared helpers
# --------------------------------------------------------------------------- #

def _maps_html(label, view):
    """One matrix line per arrow acting NON-trivially; the zero arrows are named in
    a single line (Marco 2026-07-29: an exactly-zero block carries no information,
    but omitting it silently would hide that the arrow exists). ONE implementation,
    shared with the modules section."""
    from quiverlab.trace.render_html import _arrow_maps_html
    return _arrow_maps_html(label, (view or {}).get("maps") or {})


def _is_zero(matrix):
    return all(str(x) == "0" for row in (matrix or []) for x in (row or []))


def _matrix_key(matrix):
    return tuple(tuple(str(x) for x in row) for row in matrix)


def _cols(matrix):
    rows = matrix or []
    return len(rows[0]) if rows and rows[0] is not None else 0


def _pmatrix(matrix):
    rows = matrix or []
    if not rows or not rows[0]:
        return "0"
    body = r" \\ ".join(" & ".join(str(x) for x in row) for row in rows)
    return r"\begin{pmatrix} %s \end{pmatrix}" % body


def _dv(dimvec):
    if not isinstance(dimvec, dict):
        return "{}"
    return "{" + ", ".join("%s: %s" % (k, dimvec[k]) for k in dimvec) + "}"


def _num(x):
    return "?" if x is None else str(x)


def _error_text(err):
    if isinstance(err, dict):
        return "%s: %s" % (err.get("type", "error"), err.get("message", ""))
    return str(err)


def _citations_html(b):
    cites = b.get("citations") or []
    if not cites:
        return []
    return ["<p class='ql-note'>%s</p>"
            % _esc(" · ".join(str(c[1]) for c in cites if len(c) > 1))]
