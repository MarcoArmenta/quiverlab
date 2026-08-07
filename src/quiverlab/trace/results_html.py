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
    _dims_table, _esc, _math, _math_inline, gloss_max_cells, matrix_grid)

# Human headings per compute kind. A kind missing here still renders (the key is
# shown verbatim), so a newly added invariant degrades to an honest label rather
# than vanishing from the report.
_HEADINGS = {
    "hh_cohomology": "Hochschild cohomology",
    "hh_homology": "Hochschild homology",
    "cyclic_homology": "Cyclic homology",
    "ss_hochschild": "Hochschild (b,B) spectral sequence",
    "radical_filtration_ss": "Radical-filtration spectral sequence",
    "ar_quiver": "Auslander–Reiten quiver",
    "derived_compare": "Derived fingerprint comparison",
    # Plan 35 HH product surface -- the gui.js PRODUCT_TITLE i18n titles.
    "cup": "Cup product tables",
    "cap": "Cap product tables",
    "bracket": "Gerstenhaber bracket tables",
    "connes_b": "Connes differentials",
    "cartan": "Cartan matrix",
    "coxeter_polynomial": "Coxeter polynomial",
    "global_dimension": "Global dimension",
    "homological_profile": "Homological dimensions",
    "center": "Centre",
    "dimension": "Dimension",
    "ext_algebra": "Yoneda Ext-algebra and Koszulity",
    "recognizers": "Structural recognizers and type",
    "derived_fingerprint": "Derived fingerprint",
    "strings": "Strings and bands",
    "quasi_hereditary": "Quasi-hereditary structure",
    "dimension_vector": "Dimension vector of M",
    "rad_top_soc": "Radical, top and socle of M",
    "tau": "AR translate τM",
    "tau_minus": "Inverse AR translate τ⁻M",
    "decompose": "Krull–Schmidt decomposition of M",
    "almost_split": "Almost-split sequence of M",
    "ext": "Ext",
    "tor": "Tor",
    "projective_resolution": "Projective resolution of M",
    "injective_resolution": "Injective resolution of M",
    "projective_dimension": "Projective dimension of M",
    "injective_dimension": "Injective dimension of M",
    "tilting_check": "Tilting test",
    "orbit_geometry": "Orbit geometry",
    "tau_tilting": "τ-tilting: support τ-tilting pairs, exchange graph and fan",
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
    # Marco 2026-08-03: the product sections must SAY when their recorded basis
    # differs from the route the HH sections were computed on (independent class
    # enumerations, indices do not correspond). The HH engines are the context.
    ctx = {"hh_engines": [b.get("engine") for k, b in items
                          if k in ("hh_cohomology", "hh_homology") and b.get("engine")]}
    for kind, block in items:
        heading = _HEADINGS.get(kind, kind.replace("_", " "))
        out.append("<h3 id='cr-%s'>%s</h3>" % (_esc(kind), _esc(heading)))
        err = block.get("error")
        if err:
            out.append("<p class='ql-note'>not computed — %s</p>"
                       % _esc(_error_text(err)))
            continue
        out.extend(_block_html(kind, block, ctx))
        out.extend(_citations_html(block))
    return out


# --------------------------------------------------------------------------- #
# Per-kind rendering
# --------------------------------------------------------------------------- #

# Marco 2026-08-03: an engine provenance line must SAY what the engine is, not
# just drop an internal codename ("hanlab"). Substring-keyed glosses; an engine
# string matching none renders verbatim, unchanged.
_ENGINE_GLOSS = (
    ("hanlab", "quiverlab's exact GF(p) linear-algebra engine: it assembles the "
               "boundary/coboundary matrices of the chosen (co)chain complex with "
               "integer entries mod p and computes their exact rank by Gaussian "
               "elimination mod p; every dimension follows by rank-nullity -- "
               "nothing numerical, no floating point"),
    ("chouhy", "the Chouhy–Solotar projective bimodule resolution built from "
               "the admissible presentation, certified per instance "
               "(d∘d = 0 + the order gate)"),
    ("solotar", "the Chouhy–Solotar projective bimodule resolution built from "
                "the admissible presentation, certified per instance "
                "(d∘d = 0 + the order gate)"),
    ("(b,b)", "the exact (b, B) mixed-complex engine on the normalized bar "
              "complex: b is the Hochschild boundary, B the Connes boundary"),
)


def _engine_note(engine):
    """The ``engine: ...`` provenance line, glossed when the name is one of ours."""
    s = str(engine)
    low = s.lower()
    for key, gloss in _ENGINE_GLOSS:
        if key in low:
            return ("<p class='ql-note'>engine: %s — %s.</p>"
                    % (_esc(s), _esc(gloss)))
    return "<p class='ql-note'>engine: %s</p>" % _esc(s)


def _resolved_note(kind, b):
    """Ext/Tor: WHICH module was resolved and BY WHICH resolution, stated before
    any number (Marco 2026-08-03 -- the products precedent: name the objects
    first). Rendered only when the block carries the runner's ``resolved``
    provenance (an older cached block renders as before -- tolerance)."""
    res = b.get("resolved")
    if not res:
        return []
    side = res.get("side", "right")
    resolution = res.get("resolution", "minimal projective resolution")
    if kind == "ext":
        formula = _math_inline(
            r"\operatorname{Ext}^{n}(M, N) = "
            r"H^{n}(\operatorname{Hom}_{A}(P_{\bullet}, N))")
    else:
        formula = _math_inline(
            r"\operatorname{Tor}_{n}(M, N) = H_{n}(P_{\bullet} \otimes_{A} N)")
    out = ["<p>Object resolved: the %s A-module M, by its %s "
           "%s; then %s.</p>"
           % (_esc(side), _esc(resolution),
              _math_inline(r"P_{\bullet} \to M"), formula)]
    # Marco 2026-08-03: SHOW the resolution before the data. The runner ships the
    # terms of the resolution of M it actually used (payload key `resolution`).
    rt = b.get("resolution")
    if rt and rt.get("summands"):
        rows = ["<tr><th>n</th><th>P<sub>n</sub></th></tr>"]
        for n, tex in enumerate(rt["summands"]):
            rows.append("<tr><td>%d</td><td>%s</td></tr>"
                        % (n, _math_inline(tex)))
        out.append("<p><i>The resolution of M used (terms shown to the depth "
                   "this computation needed):</i></p>")
        out.append('<table class="ql-table">%s</table>' % "".join(rows))
    if kind == "ext":
        out.append("<p class='ql-note'>N is not resolved: it enters through "
                   "Hom_A(−, N) applied to this resolution. (Equivalently one "
                   "could coresolve N injectively — Ext is balanced — but the "
                   "engine resolves M.)</p>")
    else:
        out.append("<p class='ql-note'>N is not resolved: it enters through "
                   "− ⊗_A N applied to this resolution. (Equivalently one could "
                   "resolve the left module N over A^op — Tor is balanced — but "
                   "the engine resolves M.)</p>")
    return out


_RECOGNIZER_LABELS = {
    "is_semisimple": "semisimple", "is_radical_square_zero": "radical square zero",
    "is_hereditary": "hereditary", "is_basic": "basic", "is_nakayama": "Nakayama",
    "is_special_biserial": "special biserial", "is_string": "string",
    "is_gentle": "gentle", "is_selfinjective": "self-injective",
    "is_symmetric": "symmetric",
}
_FORM_TYPE_GLOSS = {
    "finite": "positive definite Tits form (finite representation type when hereditary)",
    "tame": "positive semidefinite, not definite (tame type when hereditary)",
    "wild": "indefinite Tits form (wild type when hereditary)",
}


def _ext_algebra_html(b):
    """The Yoneda Ext-algebra block: the three-valued Koszul verdict, the graded
    (Betti) dimensions of E(A), and its minimal generators/relations by degree."""
    koszul, reason = b.get("koszul"), b.get("koszul_reason")
    obstruction, cdeg = b.get("obstruction"), b.get("certified_through_degree")
    if koszul is True:
        verdict = "A is <b>Koszul</b>"
    elif koszul is False:
        if obstruction:
            verdict = ("A is <b>not Koszul</b> — obstruction at degree %s (%s)"
                       % (_esc(str(obstruction[0])), _esc(str(obstruction[1]))))
        else:
            verdict = "A is <b>not Koszul</b>"
    else:
        verdict = "Koszulity is <b>undecided</b> through degree %s" % _esc(str(cdeg))
    if reason and koszul is not False:
        verdict += " — %s" % _esc(str(reason))
    out = ["<p>%s. Graded (Betti) data for the Yoneda algebra %s through degree "
           "%s:</p>" % (verdict,
                        _math_inline(r"E(A) = \operatorname{Ext}^{*}_{A}(A/J, A/J)"),
                        _esc(str(cdeg)))]
    dims = b.get("graded_dims") or []
    if dims:
        out.append(_dims_table("dim E^n", dims))
    if b.get("latex"):
        out.append(_math(b["latex"]))
    gens, rels = b.get("generators_by_degree") or {}, b.get("relations_by_degree") or {}

    def _by_degree(d):
        if not d:
            return "none"
        return ", ".join("degree %s: %s" % (_esc(k), _esc(str(v)))
                         for k, v in sorted(d.items(), key=lambda kv: int(kv[0])))
    out.append("<p>Minimal generators of E(A): %s. Minimal relations: %s.</p>"
               % (_by_degree(gens), _by_degree(rels)))
    return out


def _tau_tilting_html(b):
    """The C4 tau-tilting block (Plan 45): the support tau-tilting pairs (label +
    g-matrix + support), the brick-labelled Hasse edges, the AIR four-way counts +
    maximal-green-sequence count, and the wall-and-chamber fan (chambers + walls),
    with the honest complete-iff-tau-tilting-finite status."""
    n = b.get("n")
    out = ["<p>Support τ-tilting pairs (M, P) of A, their g-matrices and the mutation "
           "exchange graph (Adachi–Iyama–Reiten). Each pair is a maximal cone of the "
           "g-vector fan; each exchange edge crosses a wall labelled by a brick "
           "(King θ-stability). The BFS from (A, 0) is complete iff A is "
           "τ-tilting-finite.</p>"]
    if not b.get("complete"):
        out.append("<p class='ql-note'>The exchange graph did not close "
                   "(status: <b>%s</b>) — A is τ-tilting-infinite or the pair budget was "
                   "hit. %d pairs were found before the cap; the fan, the four-way counts "
                   "and the green-sequence count are omitted (a partial value would "
                   "mislead).</p>" % (_esc(str(b.get("status"))), b.get("num_pairs", 0)))
    counts = b.get("counts")
    if counts:
        out.append("<p>The AIR four-way count identity holds on this run: "
                   "#support τ-tilting pairs = #functorially-finite torsion classes = "
                   "#2-term silting = #semibricks = <b>%d = %d = %d = %d</b>.</p>"
                   % (counts.get("pairs"), counts.get("torsion"),
                      counts.get("silting"), counts.get("semibricks")))
    if b.get("green_count") is not None:
        out.append("<p>Maximal green sequences (directed maximal chains "
                   "(A,0) → (0,A) of left mutations): <b>%d</b>.</p>"
                   % b["green_count"])
    # pairs table
    rows = ["<tr><th>id</th><th>pair (M, P)</th><th>support</th>"
            "<th>g-matrix (columns = g-vectors)</th></tr>"]
    for p in (b.get("pairs") or []):
        gm = "; ".join("(" + ", ".join(str(p["g_matrix"][r][c])
                                       for r in range(len(p["g_matrix"])))
                       + ")" for c in range(len(p["g_matrix"][0]) if p["g_matrix"] else 0))
        star = " (initial)" if p.get("is_initial") else ""
        rows.append("<tr><td>%d</td><td>%s%s</td><td>%s</td><td>%s</td></tr>"
                    % (p["id"], _esc(str(p["label"])), star,
                       _esc(str(p.get("support") or [])), _esc(gm)))
    out.append("<table class='ql-tt-pairs'>%s</table>" % "".join(rows))
    # Hasse edges
    hasse = b.get("hasse") or []
    if hasse:
        lis = []
        for e in hasse:
            bd = e.get("brick_dimvec")
            bstr = (" — brick " + _dv(bd)) if bd else ""
            lis.append("<li>%d → %d%s</li>" % (e["from"], e["to"], _esc(bstr)))
        out.append("<p>Hasse quiver (downward = left mutation, from (A,0) to (0,A)):</p>")
        out.append("<ul class='ql-tt-hasse'>%s</ul>" % "".join(lis))
    # the fan (chambers + walls)
    fan = b.get("fan")
    if fan and fan.get("chambers"):
        out.append("<p>Wall-and-chamber fan (n = %s%s): one chamber per pair (the "
                   "g-vector cone), one wall per exchange edge (normal = the brick "
                   "dim-vector). Exact rational; the GUI draws it live.</p>"
                   % (_esc(str(fan.get("n"))),
                      ", L1/octahedron projection" if fan.get("projection") == "L1"
                      else ""))
        wl = []
        for w in (fan.get("walls") or []):
            bd = w.get("brick_dimvec")
            wl.append("<li>wall between chambers %s — normal (brick dim-vector) %s</li>"
                      % (_esc(str(w.get("between"))),
                         _dv(bd) if bd else "?"))
        if wl:
            out.append("<ul class='ql-tt-walls'>%s</ul>" % "".join(wl))
    return out


def _recognizers_html(b):
    """The recognizer batch + type block: each flag as yes/no (or an honest
    per-flag 'not decided' when it refused), the Dynkin/Euclidean type and the
    finite/tame/wild form type."""
    flags = b.get("flags") or {}
    lis = []
    for key in ("is_semisimple", "is_radical_square_zero", "is_hereditary",
                "is_basic", "is_nakayama", "is_special_biserial", "is_string",
                "is_gentle", "is_selfinjective", "is_symmetric"):
        if key not in flags:
            continue
        v = flags[key]
        label = _esc(_RECOGNIZER_LABELS.get(key, key))
        if isinstance(v, dict) and "error" in v:
            lis.append("<li>%s: not decided — %s</li>" % (label, _esc(str(v["error"]))))
        elif v is True:
            lis.append("<li>%s: <b>yes</b></li>" % label)
        else:
            lis.append("<li>%s: no</li>" % label)
    out = ["<ul class='ql-flags'>%s</ul>" % "".join(lis)] if lis else []
    dt = b.get("dynkin_type")
    out.append("<p>Underlying diagram type: <b>%s</b>.</p>" % _esc(str(dt)) if dt
               else "<p>Underlying diagram type: not a Dynkin/Euclidean diagram.</p>")
    ft = b.get("form_type")
    if ft:
        out.append("<p>Form type: <b>%s</b> — %s.</p>"
                   % (_esc(ft), _esc(_FORM_TYPE_GLOSS.get(ft, ""))))
    else:
        out.append("<p>Form type: undefined (the Cartan matrix is not unimodular, "
                   "so the Euler/Tits form has no integral matrix here).</p>")
    return out


def _quasi_hereditary_html(b):
    """The quasi-hereditary structure block (Plan 47): the verdict + the order (with the
    honest order-dependence note), the per-index brick / Delta-filtration certificates, and
    the standard-module dim-vectors table."""
    out = []
    verdict = "yes" if b.get("is_quasi_hereditary") else "no"
    order = ", ".join(str(v) for v in (b.get("order") or []))
    out.append("<p>Quasi-hereditary in this order: <b>%s</b>. Order (lowest to highest): "
               "%s.</p>" % (verdict, _esc(order)))
    if b.get("order_note"):
        out.append("<p><em>%s.</em></p>" % _esc(str(b["order_note"])))
    gl = b.get("gl_dim") or {}
    if gl:
        gtxt = ("%s (exact)" % gl.get("value") if gl.get("exact")
                else ">= %s (certified lower bound)" % gl.get("value"))
        out.append("<p>Global dimension: <b>%s</b> (quasi-heredity requires it finite).</p>"
                   % _esc(gtxt))
    if not b.get("is_quasi_hereditary") and b.get("note"):
        out.append("<p>Failing clause: %s.</p>" % _esc(str(b["note"])))
    # per-index certificates
    per = b.get("per_index") or {}
    if per:
        rows = []
        for v in (b.get("order") or list(per)):
            info = per.get(str(v)) or per.get(v) or {}
            brick = "yes" if info.get("brick") else "no"
            filt = "yes" if info.get("delta_filters_P") else "no"
            rows.append("<tr><th>%s</th><td>%s</td><td>%s</td></tr>"
                        % (_esc(str(v)), brick, filt))
        out.append("<table class='ql-qh'><thead><tr><th>i</th>"
                   "<th>End &Delta;(i)=k</th><th>P(i) &Delta;-filtered</th></tr></thead>"
                   "<tbody>%s</tbody></table>" % "".join(rows))
    # standard-module dim vectors
    sd = b.get("standard_dims") or {}
    if sd:
        rows = []
        for v in (b.get("order") or list(sd)):
            info = sd.get(str(v)) or sd.get(v) or {}
            dv = info.get("dimvec") or {}
            dvtxt = ", ".join("%s:%s" % (w, n) for w, n in dv.items())
            rows.append("<tr><th>&Delta;(%s)</th><td>%s</td><td>%s</td></tr>"
                        % (_esc(str(v)), _esc(str(info.get("dim"))), _esc(dvtxt)))
        out.append("<table class='ql-qh-delta'><thead><tr><th>standard</th><th>dim</th>"
                   "<th>dim vector</th></tr></thead><tbody>%s</tbody></table>"
                   % "".join(rows))
    # wave 2 enrichment (present only when quasi-hereditary): the characteristic
    # tilting module (summand dims) and the Ringel dual (dimension + Cartan matrix).
    ct = b.get("characteristic_tilting")
    if isinstance(ct, dict):
        if ct.get("error"):
            out.append("<p class='ql-note'>Characteristic tilting module: unavailable — "
                       "%s</p>" % _esc(str(ct["error"])))
        else:
            line = ("Characteristic tilting module T: dim %s, dim vector %s."
                    % (_num(ct.get("dim")), _esc(_dv(ct.get("dimvec") or {}))))
            summ = ct.get("summands")
            if summ:
                parts = []
                for s in summ:
                    m = s.get("mult", 1)
                    d = "%s" % _dv(s.get("dimvec") or {})
                    parts.append(d if m == 1 else "%s (&times;%s)" % (d, _num(m)))
                line += " Indecomposable summands (dim vectors): %s." % "; ".join(
                    _esc(p) for p in parts)
            elif ct.get("summands_error"):
                line += (" Summand split unavailable — %s."
                         % _esc(str(ct["summands_error"])))
            out.append("<p>%s</p>" % line)
    rd = b.get("ringel_dual")
    if isinstance(rd, dict):
        if rd.get("error"):
            out.append("<p class='ql-note'>Ringel dual R(A): unavailable — %s</p>"
                       % _esc(str(rd["error"])))
        else:
            out.append("<p>Ringel dual R(A) = End<sub>A</sub>(T)<sup>op</sup>: "
                       "dim %s.</p>" % _num(rd.get("dim")))
            if rd.get("cartan"):
                out.append(matrix_grid(rd["cartan"], label="C_{R(A)}"))
            elif rd.get("cartan_error"):
                out.append("<p class='ql-note'>Ringel dual Cartan matrix unavailable — "
                           "%s</p>" % _esc(str(rd["cartan_error"])))
            if rd.get("note"):
                out.append("<p class='ql-note'>%s</p>" % _esc(str(rd["note"])))
    return out


def _derived_fingerprint_html(b):
    """The derived fingerprint: the invariant tuple rendered as a table, plus the
    binding necessary-condition scope line. A field captured as an error prints its
    message; no field is silently dropped."""
    fp = b.get("fingerprint") or {}

    def _cell(val):
        if isinstance(val, dict) and "error" in val:
            return "unavailable — %s" % _esc(str(val["error"]))
        if isinstance(val, list):
            return _esc("[" + ", ".join(str(x) for x in val) + "]")
        return _esc(str(val))

    rows = [
        ("Coxeter polynomial", fp.get("coxeter_polynomial")),
        ("det C", fp.get("cartan_det")),
        ("Cartan Smith factors", fp.get("cartan_smith")),
        ("dim HH^• (cohomology)", fp.get("hh_cohomology_dims")),
        ("dim HH_• (homology)", fp.get("hh_homology_dims")),
        ("dim HC_• (cyclic)", fp.get("cyclic_dims")),
        ("dim Z(A)", fp.get("center_dim")),
        ("global dimension", fp.get("gl_dim")),
    ]
    body = "".join("<tr><th>%s</th><td>%s</td></tr>" % (_esc(k), _cell(v))
                   for k, v in rows)
    out = ["<table class='ql-fingerprint'><tbody>%s</tbody></table>" % body]
    out.append("<p><em>%s.</em></p>" % _esc(
        b.get("scope", "a derived-invariant fingerprint; equal values are a "
                       "necessary condition for derived equivalence, not a proof")))
    return out


def _fp_cell(val):
    """One fingerprint value as an HTML cell: a captured error prints its message, a
    list prints ``[a, b, ...]``, everything else prints verbatim (escaped)."""
    if isinstance(val, dict) and "error" in val:
        return "unavailable — %s" % _esc(str(val["error"]))
    if isinstance(val, list):
        return _esc("[" + ", ".join(str(x) for x in val) + "]")
    return _esc(str(val))


# The fingerprint fields, in the order both fingerprint renderers show them.
_FP_ROWS = (
    ("Coxeter polynomial", "coxeter_polynomial"),
    ("det C", "cartan_det"),
    ("Cartan Smith factors", "cartan_smith"),
    ("dim HH^• (cohomology)", "hh_cohomology_dims"),
    ("dim HH_• (homology)", "hh_homology_dims"),
    ("dim HC_• (cyclic)", "cyclic_dims"),
    ("dim Z(A)", "center_dim"),
    ("global dimension", "gl_dim"),
)


def _incomparable_note(fields, fa, fb):
    """Word the incomparable-fields note by HONESTLY counting the side(s) on which
    each field could not be computed -- a field can error on A only, B only, or BOTH
    (``compare_fingerprints`` collects all three into ``incomparable_fields``), so
    the flat 'errored on one side' phrasing was wrong for the both-sides case. Each
    field is annotated with where it failed: ``coxeter_polynomial (A and B)``."""
    parts = []
    for key in fields:
        sides = [s for s, fp in (("A", fa), ("B", fb))
                 if isinstance(fp.get(key), dict)]
        where = " and ".join(sides) if sides else "?"
        parts.append("%s (%s)" % (key, where))
    return ("<p class='ql-note'>Incomparable — a field could not be computed on the "
            "annotated algebra(s): %s.</p>" % _esc(", ".join(parts)))


def _derived_compare_html(b):
    """The derived_compare block (wave 2): the derived fingerprints of A and B side by
    side, the honest verdict (never an equivalence claim), and the incomparable fields
    (each annotated with the side(s) it errored on). A captured error prints its
    message per cell."""
    fa = b.get("fingerprint_a") or {}
    fb = b.get("fingerprint_b") or {}
    body = "".join(
        "<tr><th>%s</th><td>%s</td><td>%s</td></tr>"
        % (_esc(label), _fp_cell(fa.get(key)), _fp_cell(fb.get(key)))
        for label, key in _FP_ROWS)
    out = ["<table class='ql-fingerprint'><thead><tr><th>invariant</th><th>A</th>"
           "<th>B</th></tr></thead><tbody>%s</tbody></table>" % body]
    verdict = b.get("verdict_text") or b.get("verdict") or ""
    out.append("<p><b>Verdict:</b> %s.</p>" % _esc(str(verdict)))
    if b.get("incomparable_fields"):
        out.append(_incomparable_note(b["incomparable_fields"], fa, fb))
    out.append("<p><em>%s.</em></p>" % _esc(
        b.get("scope", "equal fingerprints are a necessary condition for derived "
                       "equivalence, not a proof")))
    return out


def _radical_filtration_ss_html(b):
    """The radical-filtration spectral-sequence block (wave 2): the abutment table (the
    homology of the resolution complex per degree), the E_inf netPage grid, and the
    convergence prose -- the ss_hochschild presentation, on the associated-graded SS of
    the minimal resolution of the sum of the simple modules."""
    # Defensive-only (the established ss_hochschild style): `results_section` catches
    # any block carrying an `error` field and renders the generic "not computed — …"
    # path BEFORE this per-kind renderer, so a refusal never reaches here in practice.
    # Kept total so the renderer is safe if a caller ever dispatches it directly.
    if b.get("error"):
        return ["<p class='ql-note'>%s</p>" % gloss_max_cells(_esc(str(b["error"])))]
    chunks = [_dims_table("dim E_inf total (= H_n of the complex)", b.get("abutment") or [])]
    grid = (b.get("grid") or "").replace("```", "").strip()
    if grid:
        chunks.append("<pre class=\"ql-ss-grid\">%s</pre>" % _esc(grid))
    if b.get("prose"):
        chunks.append("<p>%s</p>" % _esc(str(b["prose"])))
    return chunks


def _ar_quiver_html(b):
    """The Auslander–Reiten quiver block (wave 2): the completeness status, the
    indecomposables (name + dimension vector), the irreducible-map arrows with
    multiplicities and the tau-orbits. An honest budget cap ships a partial list; a
    self-injective / uncertifiable refusal is the library's loud message."""
    # Defensive-only (the established ss_hochschild style): `results_section` catches
    # any block carrying an `error` field and renders the generic "not computed — …"
    # path BEFORE this per-kind renderer (see test_ar_quiver_self_injective_renders),
    # so a refusal never reaches here in practice. Kept so the renderer stays total.
    if b.get("error"):
        return ["<p class='ql-note'>AR quiver not knitted: %s</p>" % _esc(str(b["error"]))]
    chunks = []
    status = b.get("status", "?")
    if b.get("complete"):
        chunks.append("<p>The algebra is representation-finite: the AR quiver has "
                      "%s indecomposables and %s irreducible-map arrows.</p>"
                      % (_num(b.get("num_vertices")), _num(b.get("num_arrows"))))
    else:
        chunks.append("<p class='ql-note'>Incomplete (status: %s) — %s</p>"
                      % (_esc(status), _esc(str(b.get("partial_note")
                         or "the knit did not close within the budget"))))
    verts = b.get("vertices") or []
    if verts:
        rows = "".join(
            "<tr><th>%s</th><td>%s</td><td>%s</td></tr>"
            % (_num(v.get("id")), _esc(str(v.get("name") or "—")), _esc(_dv(v.get("dimvec") or {})))
            for v in verts)
        chunks.append("<table class='ql-table'><thead><tr><th>id</th><th>module</th>"
                      "<th>dim vector</th></tr></thead><tbody>%s</tbody></table>" % rows)
    arrows = b.get("arrows") or []
    if arrows:
        parts = []
        for a in arrows:
            m = a.get("mult", 1)
            label = "%s &rarr; %s" % (_num(a.get("from")), _num(a.get("to")))
            parts.append(label if m == 1 else "%s (&times;%s)" % (label, _num(m)))
        chunks.append("<p>Irreducible maps: %s.</p>" % "; ".join(parts))
    orbits = b.get("tau_orbits") or []
    if orbits:
        otxt = "; ".join("{" + ", ".join(_num(x) for x in orbit) + "}" for orbit in orbits)
        chunks.append("<p>&tau;-orbits: %s.</p>" % otxt)
    return chunks


def _strings_html(b):
    """The gentle / string subsystem block (Plan 46): recognizer verdicts + string
    census + band presence + honest rep-type + (gentle) AG invariant."""
    rc = b.get("recognizers") or {}
    out = []
    lis = []
    for key, label in (("is_special_biserial", "special biserial"),
                       ("is_string", "string"), ("is_gentle", "gentle")):
        if key not in rc:
            continue
        v = rc[key]
        if isinstance(v, dict) and "error" in v:
            lis.append("<li>%s: not decided — %s</li>" % (label, _esc(str(v["error"]))))
        elif v is True:
            lis.append("<li>%s: <b>yes</b></li>" % label)
        else:
            lis.append("<li>%s: no</li>" % label)
    if lis:
        out.append("<ul class='ql-flags'>%s</ul>" % "".join(lis))
    s = b.get("strings")
    if s:
        if s.get("status") == "complete":
            out.append("<p>String modules: <b>%s</b> (a complete list up to length %s — "
                       "all indecomposable string modules).</p>"
                       % (_num(s.get("count")), _num(s.get("max_length"))))
        else:
            out.append("<p>String modules: <b>%s</b> (a length-capped sample up to length "
                       "%s; not claimed complete — the algebra has bands).</p>"
                       % (_num(s.get("count")), _num(s.get("max_length"))))
        if s.get("sample"):
            out.append("<p>Sample walks: %s.</p>"
                       % _esc(", ".join(str(x) for x in s["sample"])))
    bands = b.get("bands")
    if bands:
        if bands.get("exist"):
            out.append("<p>Bands: <b>yes</b> — %s (the algebra is representation-infinite)."
                       "</p>" % _esc(", ".join(str(x) for x in bands.get("sample", []))))
        else:
            out.append("<p>Bands: none.</p>")
    rt = {"finite": "representation-finite", "infinite": "representation-infinite",
          "unknown": "undetermined (budget/length cut)"}
    out.append("<p>Representation type: <b>%s</b>.</p>"
               % _esc(rt.get(b.get("rep_type"), str(b.get("rep_type")))))
    ag = b.get("ag_invariant")
    if isinstance(ag, list):
        pairs = ", ".join("(%s, %s)" % (p[0], p[1]) for p in ag)
        out.append("<p>Avella-Alaminos–Geiss invariant (a DERIVED invariant, honestly "
                   "NOT complete): {%s}.</p>" % _esc(pairs))
    elif isinstance(ag, dict) and "error" in ag:
        out.append("<p>AG invariant: not decided — %s.</p>" % _esc(str(ag["error"])))
    if b.get("note"):
        out.append("<p>%s.</p>" % _esc(str(b["note"])))
    return out


def _block_html(kind, b, ctx=None):
    if kind in ("hh_cohomology", "hh_homology"):
        sup = kind == "hh_cohomology"
        from quiverlab.trace.render_html import (
            hh_element_interpretation, hh_reps_sections, hh_typing_html,
            route_of_engine)
        # Typing statement at the TOP of the section (Marco 2026-07-31): exactly what
        # the engine computes and what the bar-bracket / tensor notation means.
        chunks = [hh_typing_html(kind, route_of_engine(b.get("engine"))),
                  _dims_table("dim HH%sn" % ("^" if sup else "_"), b.get("dims") or [])]
        if b.get("engine"):
            chunks.append(_engine_note(b["engine"]))
        chunks.extend(_dictionary_framing_html(kind, b.get("dims") or []))
        # Plan 35 wave 3d: the element-wise dictionary read-offs (central elements /
        # derivations / deformation cochain / commutator residues) and the per-degree
        # EXPLICIT REPRESENTATIVES, when the block carries them (hochschild.hh_reps).
        chunks.extend(hh_element_interpretation(kind, b.get("basis_classes"),
                                                b.get("inner_dims")))
        secs = hh_reps_sections(kind, b.get("basis_classes"), b.get("chain_basis"),
                                b.get("differentials"), anchor_prefix="cr-" + kind)
        if secs:
            chunks.append("<p><i>Explicit representatives by degree — each class "
                          "written over the ordered (co)chain basis, with the "
                          "differential that annihilates it. Coordinate vectors are "
                          "in the JSON.</i></p>")
            chunks.extend(secs)
        return chunks
    if kind == "cyclic_homology":
        # HC is a homology-style subscript table HC_n (Plan-35 follow-up), rendered
        # like the HH dims tables above, THEN the Plan 35 wave-3b per-degree EXPLICIT
        # REPRESENTATIVES over the total complex Tot_n = C_n (+) C_{n-2} (+) ...
        chunks = [_dims_table("dim HC_n", b.get("dims") or [])]
        if b.get("engine"):
            chunks.append(_engine_note(b["engine"]))
        chunks.extend(_dictionary_framing_html(kind, b.get("dims") or []))
        from quiverlab.trace.render_html import cyclic_degree_sections
        secs = cyclic_degree_sections(b.get("basis_classes"), b.get("chain_basis"),
                                      b.get("differentials"), b.get("column_structure"),
                                      anchor_prefix="cr")
        if secs:
            chunks.append("<p><i>Explicit representatives by degree — each class as its "
                          "representative, with the total differential b+B that "
                          "annihilates it. The total complex Tot_n and its ordered basis "
                          "(into which the coordinate vectors index) are in the "
                          "JSON.</i></p>")
            chunks.extend(secs)
        return chunks
    if kind == "ss_hochschild":
        # The Hochschild (b, B) spectral sequence (Plan 42): the abutment table
        # (E_inf totals == HC_n), the netPage E_inf grid, and the convergence prose.
        if b.get("error"):
            return ["<p>%s</p>" % gloss_max_cells(_esc(str(b["error"])))]
        chunks = [_dims_table("dim E_inf total (= HC_n)", b.get("abutment") or [])]
        grid = (b.get("grid") or "").replace("```", "").strip()
        if grid:
            chunks.append("<pre class=\"ql-ss-grid\">%s</pre>" % _esc(grid))
        if b.get("prose"):
            chunks.append("<p>%s</p>" % _esc(str(b["prose"])))
        return chunks
    if kind == "radical_filtration_ss":
        return _radical_filtration_ss_html(b)
    if kind == "ar_quiver":
        return _ar_quiver_html(b)
    if kind == "cartan":
        if b.get("matrix"):
            return [matrix_grid(b["matrix"], label="C")]
        return [_math("C = " + b["latex"])] if b.get("latex") else []
    if kind == "coxeter_polynomial":
        return [_math(r"\chi(t) = " + b["latex"])] if b.get("latex") else []
    if kind == "global_dimension":
        return ["<p>%s</p>" % _esc(str(b.get("text", "")))]
    if kind == "homological_profile":
        return _homological_profile_html(b)
    if kind == "center":
        return [_math(r"\dim Z(A) = %s" % _num(b.get("dim")))]
    if kind == "dimension":
        return [_math(r"\dim_k A = %s" % _num(b.get("value")))]
    if kind == "ext_algebra":
        return _ext_algebra_html(b)
    if kind == "tau_tilting":
        return _tau_tilting_html(b)
    if kind == "recognizers":
        return _recognizers_html(b)
    if kind == "derived_fingerprint":
        return _derived_fingerprint_html(b)
    if kind == "derived_compare":
        return _derived_compare_html(b)
    if kind == "strings":
        return _strings_html(b)
    if kind == "quasi_hereditary":
        return _quasi_hereditary_html(b)
    if kind == "dimension_vector":
        return [_math(b["latex"])] if b.get("latex") else []
    if kind == "rad_top_soc":
        return _rad_top_soc_html(b)
    if kind in ("tau", "tau_minus"):
        return _tau_html(kind, b)
    if kind == "decompose":
        return _decompose_html(b)
    if kind == "almost_split":
        return _almost_split_html(b)
    if kind in ("ext", "tor"):
        op = "Ext^" if kind == "ext" else "Tor_"
        chunks = []
        target = (b.get("target") or {}).get("dimvec")
        if target:
            chunks.append("<p>against N, dimension vector %s.</p>" % _esc(_dv(target)))
        # Marco 2026-08-03: name the objects BEFORE the numbers -- which module was
        # resolved, and by which resolution (the products precedent).
        chunks.extend(_resolved_note(kind, b))
        chunks.append(_dims_table("dim %sn" % op, b.get("dims") or []))
        # Plan 35 wave 3c: the classical dictionary framing per degree (what each class
        # MEANS), keyed only off the theory + the number of degrees computed.
        chunks.extend(_dictionary_framing_html(kind, b.get("dims") or []))
        # Plan 35 wave 3a: the per-degree EXPLICIT REPRESENTATIVES (ordered basis ->
        # classes -> differential + verification), when the block carries them.
        from quiverlab.trace.render_html import module_reps_sections
        secs = module_reps_sections(b.get("basis_classes"), b.get("chain_basis"),
                                    b.get("differentials"), kind,
                                    anchor_prefix="cr")
        if secs:
            chunks.append("<p><i>Explicit representatives by degree — each class as its "
                          "representative, with the differential that annihilates it. The "
                          "ordered basis (into which the coordinate vectors index) is in "
                          "the JSON.</i></p>")
            chunks.extend(secs)
        # Plan 35 wave 3c: the Yoneda exact-sequence interpretation (Ext only) -- each
        # class as the constructed + certified exact sequence 0 -> N -> Q -> ... -> M -> 0.
        if kind == "ext":
            from quiverlab.trace.render_html import ext_interpretation_sections
            chunks.extend(ext_interpretation_sections(b.get("interpretation"),
                                                      anchor_prefix="cr"))
        return chunks
    if kind in ("cup", "cap", "bracket"):
        return _product_tables_html(kind, b, ctx)
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
    if kind == "tilting_check":
        if b.get("error"):
            return ["<p class='ql-note'>%s</p>" % _esc(str(b["error"]))]
        verdict = ("M is a tilting module."
                   if b.get("is_tilting")
                   else "M is not a tilting module: %s." % _esc(str(b.get("note", ""))))
        pd = b.get("pd")
        rows = [
            ("tilting", "yes" if b.get("is_tilting") else "no"),
            ("n (pd bound tested)", _num(b.get("n"))),
            ("pd M", _num(pd) if pd is not None else "&gt; bound"),
            ("Ext<sup>i</sup>(M,M)=0 (1&#8804;i&#8804;n)",
             "yes" if b.get("self_ext_vanishes") else "no"),
            ("# non-iso indec. summands", _num(b.get("num_summands"))),
            ("# vertices (rank K<sub>0</sub>)", _num(b.get("num_vertices"))),
        ]
        tbl = "".join("<tr><th>%s</th><td>%s</td></tr>" % (k, v) for k, v in rows)
        return ["<p>%s</p>" % verdict,
                '<table class="ql-table">%s</table>' % tbl]
    if kind == "orbit_geometry":
        return _orbit_geometry_html(b)
    # An unknown kind still leaves a trace of what was asked for.
    return ["<p class='ql-note'>computed; see the JSON record for its data.</p>"]


def _orbit_geometry_html(b):
    """The orbit_geometry block (Plan 49 / C8): orbit dim + rigidity + honest codim +
    (hereditary Dynkin) the Kac canonical decomposition."""
    if b.get("error"):
        return ["<p class='ql-note'>%s</p>" % _esc(str(b["error"]))]
    chunks = []
    if b.get("latex"):
        chunks.append(_math(b["latex"]))
    rows = [
        ("dimension vector d", _esc(_dv(b.get("dim_vector") or {}))),
        ("dim GL(d) = &sum;<sub>v</sub> d<sub>v</sub><sup>2</sup>", _num(b.get("group_dim"))),
        ("dim Rep(Q,d) (ambient)", _num(b.get("rep_variety_dim"))),
        ("dim End<sub>A</sub>(M)", _num(b.get("end_dim"))),
        ("dim O<sub>M</sub> (orbit)", _num(b.get("orbit_dim"))),
        ("dim Ext<sup>1</sup>(M,M)", _num(b.get("ext1_self"))),
    ]
    tbl = "".join("<tr><th>%s</th><td>%s</td></tr>" % (k, v) for k, v in rows)
    chunks.append('<table class="ql-table">%s</table>' % tbl)
    # rigidity verdict + honest codim gloss (hereditary = codim; general = upper bound)
    rigid = bool(b.get("rigid"))
    verdict = ("M is rigid: Ext<sup>1</sup>(M,M) = 0, so the orbit O<sub>M</sub> is open (Voigt)."
               if rigid else
               "M is not rigid: Ext<sup>1</sup>(M,M) &gt; 0.")
    if b.get("codim_semantics") == "hereditary":
        gloss = ("A is hereditary, so dim Ext<sup>1</sup>(M,M) IS the codimension of the "
                 "orbit closure in Rep(Q,d) (Voigt, Rep smooth).")
    else:
        gloss = ("A = kQ/I is not hereditary, so dim Ext<sup>1</sup>(M,M) is only an UPPER "
                 "BOUND on the codimension (the module variety is cut by the relations).")
    chunks.append("<p>%s %s</p>" % (verdict, gloss))
    # the Kac canonical decomposition (hereditary Dynkin) or the honest refusal note
    cd = b.get("canonical_decomposition")
    if cd:
        parts = []
        for c in cd:
            raw = c.get("name") or ("(%s)" % ", ".join(str(x) for x in c.get("root", ())))
            name = _esc(raw)
            m = c.get("multiplicity", 1)
            parts.append(name if m == 1 else "%s<sup>%d</sup>" % (name, m))
        chunks.append("<p>Kac canonical decomposition: d = %s (each component a positive "
                      "root; the generic module is rigid).</p>" % " &oplus; ".join(parts))
    elif b.get("canonical_note"):
        chunks.append("<p class='ql-note'>Canonical decomposition not computed: %s</p>"
                      % _esc(str(b["canonical_note"])))
    # wave 2 enrichment: the P49 degeneration (= hom) order of M's dimension vector,
    # present only when it is informative (a genuine multi-class poset) or an honest
    # library refusal. Absent (byte-identical) for a unique module of that dim vector.
    deg = b.get("degeneration_order")
    if isinstance(deg, dict):
        if not deg.get("complete"):
            chunks.append("<p class='ql-note'>Degeneration order not computed "
                          "(status: %s) — %s</p>"
                          % (_esc(str(deg.get("status", "?"))),
                             _esc(str(deg.get("note", "")))))
        else:
            rows = []
            for v in deg.get("vertices", []):
                summ = " &oplus; ".join(
                    (s.get("name") or _dv(s.get("dimvec") or {}))
                    if s.get("mult", 1) == 1
                    else "%s(&times;%s)" % (s.get("name") or "?", s.get("mult"))
                    for s in v.get("summands", []))
                marks = []
                if v.get("is_generic"):
                    marks.append("generic")
                if deg.get("m_index") == v.get("index"):
                    marks.append("= M")
                tag = (" (%s)" % ", ".join(marks)) if marks else ""
                rows.append("<tr><th>%s</th><td>%s</td><td>%s</td></tr>"
                            % (_num(v.get("index")), _esc(summ + tag),
                               _num(v.get("orbit_dim"))))
            chunks.append("<p>Degeneration (hom) order of the iso-classes of dimension "
                          "vector d (a &le;<sub>deg</sub> b means a degenerates from b):</p>")
            chunks.append("<table class='ql-table'><thead><tr><th>class</th>"
                          "<th>module</th><th>dim O</th></tr></thead><tbody>%s"
                          "</tbody></table>" % "".join(rows))
            covers = deg.get("covers") or []
            if covers:
                ctxt = "; ".join("%s &lt; %s" % (_num(a), _num(bb)) for a, bb in covers)
                chunks.append("<p>Hasse covers: %s.</p>" % ctxt)
    return chunks


def _dictionary_framing_html(theory, dims):
    """Plan 35 wave 3c -- the classical dictionary framing per degree: what each class
    of this space MEANS. Keyed only off the theory + the number of degrees computed
    (a shared prose builder both surfaces use). The higher-degree framing sentence,
    identical across degrees, is shown once."""
    from quiverlab.trace.interpretations import sentence
    lis, seen = [], set()
    for n in range(len(dims)):
        s = sentence(theory, n)
        if s is None or s in seen:
            continue
        seen.add(s)
        lis.append("<li>%s</li>" % _esc(s))
    if not lis:
        return []
    return ["<p><i>Interpretation of the spaces (the classical dictionary):</i></p>",
            "<ul class='ql-interp'>%s</ul>" % "".join(lis)]


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
    out.extend(_loewy_series_html(b.get("series") or []))
    for label, view in trio:
        out.extend(_maps_html(label, view))
    return out


def _loewy_series_html(series):
    """The Loewy (radical) series as a stacked diagram, one row per layer top to
    bottom, factors as S_v^m (Plan 37). Empty series => nothing rendered."""
    if not series:
        return []
    rows = []
    for i, layer in enumerate(series):
        factors = " ⊕ ".join(
            ("S_%s" % v if m == 1 else "S_%s^%d" % (v, m))
            for v, m in sorted(layer.items()) if m)
        rows.append("<tr><th>layer %d</th><td>%s</td></tr>"
                    % (i + 1, _esc(factors or "0")))
    return ["<p class='ql-note'>Loewy (radical) series, top to bottom:</p>",
            '<table class="ql-loewy">%s</table>' % "".join(rows)]


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


def _findim_text(f):
    if f.get("exact"):
        return "findim = %s  (%s)" % (_num(f.get("lower")), f.get("note", ""))
    if f.get("upper") is not None:
        return "findim in [%s, %s]  (%s)" % (
            _num(f.get("lower")), _num(f.get("upper")), f.get("note", ""))
    return "findim >= %s  (%s)" % (_num(f.get("lower")), f.get("note", ""))


def _homological_profile_html(b):
    """The C6 homological-dimensions family (Plan 40) as a labelled table -- each row
    its own honest marker (exact value / certified bound / infinity / undecided /
    per-entry error), never a bare number the engine did not resolve."""
    it = b.get("igusa_todorov") or {}
    if it.get("error"):
        it_text = "not computed: %s" % it["error"]
    else:
        it_text = "of %s: φ = %s, ψ = %s" % (
            it.get("module", ""), _num(it.get("phi")), _num(it.get("psi")))
    rows = [
        ("Global dimension", (b.get("global_dimension") or {}).get("text", "")),
        ("Finitistic dimension", _findim_text(b.get("finitistic") or {})),
        ("Dominant dimension", (b.get("dominant") or {}).get("text", "")),
        ("Gorenstein", (b.get("gorenstein") or {}).get("text", "")),
        ("Igusa–Todorov φ/ψ", it_text),
    ]
    body = "".join("<tr><th>%s</th><td>%s</td></tr>" % (_esc(lab), _esc(str(val)))
                   for lab, val in rows)
    return ['<table class="ql-table">%s</table>' % body]
def _almost_split_html(b):
    """The almost-split sequence 0 → τM → E → M → 0 (Plan 41): τM as a full
    representation, then E's Krull–Schmidt summands (standard summands named, others
    with full matrices). An ``exists: false`` block renders the honest refusal."""
    if b.get("exists") is False:
        return ["<p>No almost-split sequence: %s</p>"
                % _esc(b.get("reason", "input not eligible"))]
    out = []
    if b.get("latex"):
        out.append(_math(b["latex"]))
    out.append("<p>M is indecomposable and non-projective; τM as a full "
               "representation:</p>")
    if b.get("tau"):
        out.extend(_maps_html("τM", b["tau"]))
    summands = (b.get("middle") or {}).get("summands") or []
    rows = ["<tr><th>summand</th><th>multiplicity</th><th>dim vector</th></tr>"]
    for i, s in enumerate(summands):
        rows.append("<tr><td>%s</td><td>%s</td><td>%s</td></tr>"
                    % (_esc(_summand_name(s, i + 1)), _num(s.get("multiplicity")),
                       _esc(_dv(s.get("dim_vector")))))
    out.append("<p>middle term E — %d Krull–Schmidt summand(s):</p>" % len(summands))
    out.append('<table class="ql-table">%s</table>' % "".join(rows))
    for i, s in enumerate(summands):
        if s.get("standard") or not s.get("maps"):
            continue
        out.extend(_maps_html(_summand_name(s, i + 1), s))
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


# ordered term-basis entries shown inline before a machine-record pointer (Marco
# 2026-08-02: at most the first 20 basis elements of a listed space are displayed).
_TERM_BASIS_DISPLAY = 20


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
        if _is_zero(matrix):
            # Marco 2026-08-03: a zero map is STATED, never drawn -- and never
            # cross-referenced ("d_3 = d_1" between two zero maps hides the fact).
            out.append(matrix_grid(matrix, label=sym))
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


def _product_tables_html(kind, b, ctx=None):
    """cup / cap / bracket: the notation legend, ONE flat list of ALL (co)homology basis
    classes across degrees (Marco 2026-08-02 -- the chain bases live in the HH sections
    above, so the products just remind the classes then show the table), then the whole
    family's Cayley multiplication table, the bracket's served-window note, and the engine
    provenance. The flat class list comes from the SAME builder the worked-steps chapter
    uses (``quiverlab.trace.render_html.product_flat_classes_html``) and the Cayley table
    from ``family_cayley_html`` -- one implementation, no drift."""
    from quiverlab.trace.products import (
        notation_legend, balanced_rep_note, prime_from_basis)
    from quiverlab.trace.render_html import (
        product_flat_classes_html, family_cayley_html)
    prime = prime_from_basis(b.get("basis"))
    out = ["<p class='ql-note'>%s</p>"
           % _esc(notation_legend(kind, "", b.get("basis")))]
    if prime is not None:                             # the balanced-rep legend, once
        out.append("<p class='ql-note'>%s</p>" % _esc(balanced_rep_note(prime)))
    # Marco 2026-08-03: structure constants are basis-dependent. When the HH
    # sections were computed on a DIFFERENT route (say Chouhy-Solotar) than this
    # table's recorded basis (say bar), the two class enumerations are independent
    # -- say so, loudly, so nobody cross-reads indices between them. No coordinates
    # from different resolutions are ever mixed inside one table.
    from quiverlab.trace.render_html import route_of_engine
    hh_routes = {route_of_engine(e) for e in (ctx or {}).get("hh_engines", [])}
    if hh_routes and route_of_engine(b.get("basis")) not in hh_routes:
        out.append("<p class='ql-note'>Note: this table's classes are enumerated "
                   "over its recorded basis (%s); the Hochschild sections above "
                   "were computed on a different route (%s). The two class lists "
                   "are independent enumerations — their indices do not correspond "
                   "— and no coordinates from different resolutions are mixed "
                   "inside this table; the dimension tables agree because "
                   "dimensions are basis-independent.</p>"
                   % (_esc(str(b.get("basis"))),
                      _esc(", ".join(sorted(str(e) for e in (ctx or {}).get("hh_engines", []))))))
    # Marco 2026-08-02: for the products, one flat list of ALL (co)homology basis classes
    # across degrees, then the multiplication table right away -- no per-degree
    # sub-sections (those live in the HH cohomology/homology sections above). The flat
    # list's own intro line heads it (symmetric with the worked-steps chapter).
    out.extend(product_flat_classes_html(b.get("basis_classes")))
    # ONE big degree-major Cayley table for the whole family (Marco 2026-08-01), with
    # the em-dash beyond-window mark and the >cap per-bidegree fallback -- the SAME
    # builder + renderer the worked-steps chapter uses (no drift).
    out.extend(family_cayley_html(kind, list(b.get("tables") or []), prime))
    if kind == "bracket" and b.get("window") is not None:
        out.append("<p class='ql-note'>bracket structure constants are served for "
                   "arguments of total degree ≤ %s — the largest window this "
                   "bar-route computation certifies; a cell beyond it is marked "
                   "—. The bracket is computed entirely on the bar (co)chain "
                   "route.</p>" % _num(b.get("window")))
    if b.get("engine"):
        out.append(_engine_note(b["engine"]))
    return out


def _connes_b_html(b):
    """connes_b (Marco 2026-08-02): the legend, ONE flat homology-cycle class list
    (z^n_j, degree-major, all degrees at once -- the SAME builder the products use, no
    per-degree sub-sections and no chain enumerations), then one induced Connes
    differential grid ``B_n : HH_n → HH_{n+1}`` per degree with its induced rank, and the
    engine provenance."""
    from quiverlab.trace.products import notation_legend
    from quiverlab.trace.render_html import product_flat_classes_html
    out = ["<p class='ql-note'>%s</p>"
           % _esc(notation_legend("connes_b", "", None))]
    out.extend(product_flat_classes_html(b.get("basis_classes")))
    out.append("<p class='ql-note'>Each cycle class is written over the chain basis "
               "enumerated in the Hochschild homology sections / JSON; the induced "
               "matrices below act on these classes (rows index HH_{n+1}, columns "
               "HH_n).</p>")
    matrices = b.get("matrices") or {}
    ranks = b.get("ranks") or {}
    for key in sorted(matrices, key=lambda s: int(s)):
        n = int(key)
        out.append('<p class="ql-mlabel">%s</p>'
                   % _math_inline(r"B_{%d} : HH_{%d} \to HH_{%d}" % (n, n, n + 1)))
        out.append(matrix_grid(matrices[key], label="B_{%d}" % n))
        out.append("<p class='ql-note'>induced rank B_%d = %s</p>"
                   % (n, _num(ranks.get(key))))
    if b.get("engine"):
        out.append(_engine_note(b["engine"]))
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
    from quiverlab.trace.render_html import gloss_max_cells
    if isinstance(err, dict):
        return gloss_max_cells("%s: %s" % (err.get("type", "error"),
                                           err.get("message", "")))
    return gloss_max_cells(str(err))


def _citations_html(b):
    cites = b.get("citations") or []
    if not cites:
        return []
    return ["<p class='ql-note'>%s</p>"
            % _esc(" · ".join(str(c[1]) for c in cites if len(c) > 1))]
