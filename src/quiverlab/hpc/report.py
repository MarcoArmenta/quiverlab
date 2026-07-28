"""Plan 28 -- render a ``result.json`` into a human report (HTML / plain text).

The report is a self-contained, no-JS HTML page or plain text; ``--format txt``
and ``--format html`` always work with no toolchain, and ``--format json`` emits
the worked-steps event stream (trace.json). PDF/TeX report output has been
removed -- ``fmt="pdf"``/``"tex"`` is refused loudly.

The HTML shows RENDERED math -- matrices as tables, sub/superscripted summand
notation ``P_1^2 + P_3`` -- never LaTeX source (Marco, 2026-07-28). Oversized
matrices keep the Plan-34 contract: past 25 rows/columns they are STATED-elided
(shown in full in the JSON result), so one huge action matrix cannot dominate
the page.

The envelope carries ``result_schema`` (an int; absent == legacy 1). ``render``
refuses a NEWER ``result_schema`` loudly (the CLI maps that to exit 65) and warns
on library-version skew between the result and the running quiverlab."""
from __future__ import annotations

import json
import pathlib
import re

import quiverlab as ql
from quiverlab.hpc.spec import RESULT_SCHEMA


class ReportError(Exception):
    """The report could not be produced (unsupported result_schema, unwritable
    output, or a removed output format that was explicitly requested)."""


class ResultSchemaError(ReportError):
    """The result was written by a NEWER envelope schema than this tool speaks
    (the CLI maps this to exit 65)."""


class ReportWriteError(ReportError):
    """The report output path could not be written (the CLI maps this to 73)."""


# --------------------------------------------------------------------------- #
# Loading + version gating
# --------------------------------------------------------------------------- #

def load_result(path) -> dict:
    p = pathlib.Path(path)
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ReportError(f"cannot read result file {p}: {exc}")
    except json.JSONDecodeError as exc:
        raise ReportError(f"{p} is not valid JSON: {exc}")


def check_result_schema(result: dict) -> None:
    schema = result.get("result_schema", 1)
    if not isinstance(schema, int) or isinstance(schema, bool):
        raise ReportError(f"result_schema {schema!r} is not an integer")
    if schema > RESULT_SCHEMA:
        raise ResultSchemaError(
            f"result_schema {schema} is newer than this quiverlab-hpc supports "
            f"({RESULT_SCHEMA}); upgrade quiverlab to render it")


def _version_warning(result: dict) -> str | None:
    got = result.get("quiverlab_version")
    cur = getattr(ql, "__version__", "unknown")
    if got and got != cur:
        return (f"result was produced by quiverlab {got}, rendering with {cur}; "
                "numbers are shown verbatim from the file")
    return None


# --------------------------------------------------------------------------- #
# Small rendering helpers (shared by the txt and html builders)
# --------------------------------------------------------------------------- #

def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _dimvec_str(dv: dict) -> str:
    return "(" + ", ".join(f"{k}: {v}" for k, v in dv.items()) + ")"


# One summand of a resolution term as produced by the spec core:
# ``P_{1}``, ``P_{1}^{2}``, ``I_{v}^{3}`` (and ``S_{v}`` in Loewy stacks).
_SUMMAND_RE = re.compile(r"([PIS])_\{([^{}]*)\}(?:\^\{([0-9]+)\})?")


def _summand_txt(s: str) -> str:
    """``"P_{1}^{2} \\oplus P_{3}"`` -> ``"P_1^2 + P_3"`` (Marco's notation);
    anything unrecognized is passed through verbatim."""
    if not s or s == "0":
        return "0"
    parts = []
    for p in s.split(r"\oplus"):
        m = _SUMMAND_RE.fullmatch(p.strip())
        if m is None:
            return s
        letter, v, e = m.groups()
        parts.append(f"{letter}_{v}" + (f"^{e}" if e else ""))
    return " + ".join(parts)


def _summand_html(s: str) -> str:
    """Same summand string rendered as real HTML sub/superscripts."""
    if not s or s == "0":
        return "0"
    parts = []
    for p in s.split(r"\oplus"):
        m = _SUMMAND_RE.fullmatch(p.strip())
        if m is None:
            return _esc(s)
        letter, v, e = m.groups()
        h = f"{letter}<sub>{_esc(v)}</sub>"
        if e:
            h += f"<sup>{_esc(e)}</sup>"
        parts.append(h)
    return " + ".join(parts)


def _matrix_html(mat) -> str:
    """A matrix as a rendered HTML table with 1-based row/column index headers
    (a header row on top, an index column on the left) over a white/grey striped
    grid, so an entry's position reads off directly. Past 25 rows or columns it
    is STATED-elided (the full matrix lives in the JSON result), keeping the
    Plan-34 wide-matrix contract."""
    nrows = len(mat)
    ncols = max((len(row) for row in mat), default=0)
    if nrows > 25 or ncols > 25:
        return ("<em>[%d&times;%d matrix -- shown in full in the JSON result]</em>"
                % (nrows, ncols))
    if nrows == 0 or ncols == 0:
        return "<em>(empty matrix)</em>"
    head = ("<tr><th></th>"
            + "".join(f"<th>{j + 1}</th>" for j in range(ncols)) + "</tr>")
    body = "".join(
        f"<tr><th>{i + 1}</th>"
        + "".join(f"<td>{_esc(str(x))}</td>" for x in row) + "</tr>"
        for i, row in enumerate(mat))
    return f"<table class='matrix'>{head}{body}</table>"


def _expr_html(text: str) -> str:
    """A relation / polynomial string rendered lightly for HTML: products drop the
    ``*``, powers become superscripts (handles both ``^n`` and sympy's ``**n``)."""
    t = _esc(text)
    t = re.sub(r"\*\*(\d+)", r"<sup>\1</sup>", t)
    t = re.sub(r"\^(\d+)", r"<sup>\1</sup>", t)
    t = t.replace("*", "&middot;")
    return t


def _hh_txt(kind_label: str, i: int, d) -> str:
    return f"{kind_label}{i} = {d}"


def _hh_html(kind_label: str, i: int, d) -> str:
    tag = "sup" if kind_label.endswith("^") else "sub"
    return f"HH<{tag}>{i}</{tag}> = {d}"


# --------------------------------------------------------------------------- #
# Intermediate: result.json -> a list of neutral sections
# --------------------------------------------------------------------------- #

class _Section:
    __slots__ = ("title", "rows", "html")

    def __init__(self, title, rows=None, html=None):
        self.title = title
        self.rows = rows or []          # plain-text lines (the txt report)
        self.html = html or []          # raw-HTML blocks; when present, the HTML
        #                                 report shows these INSTEAD of ``rows``


def _hh_section(kind_label, block) -> _Section:
    dims = block.get("dims", [])
    rows = [f"engine: {block.get('engine', '?')}",
            "dimensions: " + ", ".join(_hh_txt(kind_label, i, d)
                                       for i, d in enumerate(dims))]
    html = [f"<p>engine: {_esc(str(block.get('engine', '?')))}</p>",
            "<p>" + ",&ensp;".join(_hh_html(kind_label, i, d)
                                   for i, d in enumerate(dims)) + "</p>"]
    return _Section(f"Hochschild {'cohomology' if kind_label == 'HH^' else 'homology'}",
                    rows, html)


def _module_view_row(label, view) -> str:
    return f"{label}: dim {view.get('dim')}  {_dimvec_str(view.get('dimvec', {}))}"


def _repr_rows_html(label, view) -> tuple:
    """(txt rows, html blocks) for a full-representation module block
    ``{"dims": ..., "maps": ...}`` (Plan 34): the dimension VECTOR and the exact
    per-arrow action matrices, each labeled by its arrow (referenceable against
    the arrow list in the Computation section)."""
    rows = [f"{label}: dimension vector {_dimvec_str(view.get('dims', {}))}"]
    html = [f"<p><b>{_esc(label)}</b> &mdash; dimension vector "
            f"{_esc(_dimvec_str(view.get('dims', {})))}</p>"]
    maps = view.get("maps") or {}
    if not maps:
        rows.append("    (all arrow actions vanish)")
        html.append("<p><em>(all arrow actions vanish)</em></p>")
    for a in sorted(maps):
        rows.append(f"    action of arrow {a}: {maps[a]}")
        html.append("<div class='arrowmap'><span>action of arrow "
                    f"<code>{_esc(a)}</code>:</span> {_matrix_html(maps[a])}</div>")
    return rows, html


def _section_for(name: str, block: dict) -> _Section:
    if name == "hh_cohomology":
        return _hh_section("HH^", block)
    if name == "hh_homology":
        return _hh_section("HH_", block)
    if name == "cartan":
        mat = block.get("matrix", [])
        return _Section("Cartan matrix",
                        ["C = " + "; ".join(str(r) for r in mat)],
                        [f"<div class='arrowmap'><span>C =</span> {_matrix_html(mat)}</div>"])
    if name == "coxeter_polynomial":
        text = block.get("text", "")
        return _Section("Coxeter polynomial", [text],
                        [f"<p>{_expr_html(text)}</p>"])
    if name == "global_dimension":
        return _Section("Global dimension",
                        [block.get("text", ""), f"value: {block.get('value')} "
                         f"(exact: {block.get('exact')})"])
    if name == "center":
        rows = [f"dimension: {block.get('dim')}"]
        for row in block.get("basis", []):
            rows.append("  basis: [" + ", ".join(str(x) for x in row) + "]")
        return _Section("Center Z(A)", rows)
    if name == "dimension":
        return _Section("Algebra dimension", [f"dim A = {block.get('value')}"])
    if name == "dimension_vector":
        return _Section("Module dimension vector",
                        [f"side: {block.get('side')}",
                         _module_view_row("module", block)])
    if name == "rad_top_soc":
        rows = [f"side: {block.get('side')}"]
        html = [f"<p>side: {_esc(str(block.get('side')))}</p>"]
        for label, key in (("rad M", "radical"), ("top M", "top"), ("soc M", "socle")):
            r, h = _repr_rows_html(label, block.get(key, {}))
            rows += r
            html += h
        return _Section("Radical / top / socle", rows, html)
    if name in ("tau", "tau_minus"):
        sym = "tau M" if name == "tau" else "tau^- M"
        rows = [f"side: {block.get('side')}"]
        html = [f"<p>side: {_esc(str(block.get('side')))}</p>"]
        if block.get("is_zero"):
            rows.append(f"{sym} = 0")
            html.append(f"<p>{_esc(sym)} = 0</p>")
        elif isinstance(block.get("repr"), dict):
            r, h = _repr_rows_html(sym, block["repr"])
            rows += r
            html += h
        else:                                   # older result: dim vector only
            rows.append(_module_view_row(sym, block))
            html.append(f"<p>{_esc(_module_view_row(sym, block))}</p>")
        # The Plan-30 certificate is about the INPUT M, not the translate:
        # True means M is indecomposable (hence so is a nonzero translate);
        # the additive note means M decomposes and tau distributes over (+).
        if block.get("note_key") == "mod.tau_additive":
            n_sum = len(block.get("decomposition") or [])
            line = (f"note: the input M is decomposable ({n_sum} indecomposable "
                    f"summands, see the Krull-Schmidt section); the AR translate "
                    f"is additive, so {sym} is the direct sum of the summands' "
                    "translates (projective/injective summands are killed)")
            rows.append(line)
            html.append(f"<p><em>{_esc(line)}</em></p>")
        elif block.get("indecomposable") is True:
            line = (f"indecomposable: True (the input M is indecomposable, hence "
                    f"so is {sym} when nonzero)")
            rows.append(line)
            html.append(f"<p>{_esc(line)}</p>")
        return _Section("Auslander-Reiten translate", rows, html)
    if name == "ext":
        dims = block.get("dims", [])
        rows = ["target " + _module_view_row("N", block.get("target", {})),
                "Ext dims: " + ", ".join(f"Ext^{i} = {d}" for i, d in enumerate(dims))]
        html = ["<p>target " + _esc(_module_view_row("N", block.get("target", {})))
                + "</p>",
                "<p>" + ",&ensp;".join(f"Ext<sup>{i}</sup> = {d}"
                                       for i, d in enumerate(dims)) + "</p>"]
        return _Section("Ext^*(M, N)", rows, html)
    if name == "tor":
        dims = block.get("dims", [])
        rows = ["target " + _module_view_row("N", block.get("target", {})),
                "Tor dims: " + ", ".join(f"Tor_{i} = {d}" for i, d in enumerate(dims))]
        html = ["<p>target " + _esc(_module_view_row("N", block.get("target", {})))
                + "</p>",
                "<p>" + ",&ensp;".join(f"Tor<sub>{i}</sub> = {d}"
                                       for i, d in enumerate(dims)) + "</p>"]
        return _Section("Tor_*(M, N)", rows, html)
    if name == "decompose":
        summands = block.get("summands", [])
        rows = [f"side: {block.get('side')}",
                f"indecomposable summands: {block.get('iso_classes')} iso-class(es)"]
        html = [f"<p>side: {_esc(str(block.get('side')))}</p>",
                f"<p>indecomposable summands: "
                f"{_esc(str(block.get('iso_classes')))} iso-class(es)</p>"]
        for s in summands:
            line = (f"  dim vector {_dimvec_str(s.get('dim_vector', {}))}"
                    f" x {s.get('multiplicity')}"
                    f"  (indecomposable: {s.get('indecomposable')})")
            rows.append(line)
            html.append(f"<p class='resterm'>{_esc(line.strip())}</p>")
        return _Section("Krull-Schmidt decomposition", rows, html)
    if name in ("projective_resolution", "injective_resolution"):
        proj = name == "projective_resolution"
        letter = "P" if proj else "I"
        terms = block.get("terms", [])
        betti = block.get("betti", [])
        summands = block.get("summands", [])
        diffs = block.get("differentials", [])
        top = block.get("top")
        rows = [f"computed window: degrees 0..{top}"]
        html = [f"<p>computed window: degrees 0..{_esc(str(top))}</p>"]
        if diffs:
            conv = ("differential matrices: rows index the target basis, columns "
                    "the source basis (vertex-ordered)")
            rows.append(conv)
            html.append(f"<p class='dv'>{_esc(conv)}</p>")
        for i, dv in enumerate(terms):
            # Preferred: the summand notation P_j^n + P_i^m per degree; fall back
            # to Betti counts when an older result carries no summand strings.
            if i < len(summands):
                s_txt, s_html = _summand_txt(summands[i]), _summand_html(summands[i])
            else:
                b = betti[i] if i < len(betti) else "?"
                s_txt = s_html = f"{letter}^{b}"
            rows.append(f"  deg {i}: {s_txt}   (dim vector {_dimvec_str(dv)})")
            html.append(f"<p class='resterm'>deg {i}:&ensp;{s_html}&ensp;"
                        f"<span class='dv'>(dim vector "
                        f"{_esc(_dimvec_str(dv))})</span></p>")
            if i < len(diffs) and isinstance(diffs[i], dict):
                d = diffs[i]
                if proj:      # maps OUT of degree i: eps: Q_0 -> M, d_i: Q_i -> Q_{i-1}
                    label_txt = ("eps: deg 0 -> M" if i == 0
                                 else f"d_{i}: deg {i} -> deg {i - 1}")
                    label_html = ("&epsilon;: deg 0 &rarr; M" if i == 0
                                  else f"d<sub>{i}</sub>: deg {i} &rarr; deg {i - 1}")
                else:         # maps INTO degree i: iota: M -> E^0, d^i: E^{i-1} -> E^i
                    label_txt = ("iota: M -> deg 0" if i == 0
                                 else f"d^{i}: deg {i - 1} -> deg {i}")
                    label_html = ("&iota;: M &rarr; deg 0" if i == 0
                                  else f"d<sup>{i}</sup>: deg {i - 1} &rarr; deg {i}")
                if d.get("elided"):
                    note = (f"[{d.get('rows')}x{d.get('cols')} matrix -- elided "
                            "past the size cap; rerun with a worked-steps trace "
                            "for the full map]")
                    rows.append(f"    {label_txt}: {note}")
                    html.append(f"<p class='resterm'>&ensp;{label_html}: "
                                f"<em>{_esc(note)}</em></p>")
                else:
                    mat = d.get("matrix", [])
                    rows.append(f"    {label_txt}: {mat}")
                    html.append("<div class='arrowmap'><span>"
                                f"{label_html}:</span> {_matrix_html(mat)}</div>")
        def _dim_line(label, v):
            if v is not None:
                return f"{label}: {v}"
            return (f"{label}: not resolved within the computed window "
                    f"(> {top}; possibly infinite)")
        if "pd" in block:
            line = _dim_line("projective dimension", block["pd"])
            rows.append(line)
            html.append(f"<p>{_esc(line)}</p>")
        if "injective_dimension" in block:
            line = _dim_line("injective dimension", block["injective_dimension"])
            rows.append(line)
            html.append(f"<p>{_esc(line)}</p>")
        label = ("Projective resolution" if name == "projective_resolution"
                 else "Injective resolution")
        return _Section(label, rows, html)
    if name == "projective_dimension":
        v = block.get("value")
        return _Section("Projective dimension",
                        [f"pd M = {v if v is not None else 'infinite'} "
                         f"(finite: {block.get('finite')})"])
    if name == "injective_dimension":
        v = block.get("value")
        return _Section("Injective dimension",
                        [f"id M = {v if v is not None else 'infinite'} "
                         f"(finite: {block.get('finite')})"])
    # Unknown block: dump its keys honestly rather than dropping it -- but never
    # the inline citation texts (the full references live at the bottom of the
    # report); only the reference keys are mentioned.
    slim = {k: v for k, v in block.items()
            if k not in ("citations", "references", "latex")}
    rows = [json.dumps(slim, default=str, sort_keys=True)]
    keys = block.get("references") or []
    if keys:
        rows.append("references: " + ", ".join(str(k) for k in keys)
                    + "  (full citations at the bottom)")
    return _Section(name, rows)


def _rebuild_algebra(result: dict):
    """Best-effort algebra rebuild from the result's algebra spec (for the
    descriptive sections). Never raises -- returns None on any failure."""
    algebra = result.get("algebra")
    if not isinstance(algebra, dict):
        return None
    try:
        from quiverlab.hpc.spec import build_algebra
        return build_algebra(algebra)
    except Exception:
        return None


def _computation_section(result: dict) -> _Section:
    """The report header: version, the raw algebra spec, and -- whenever the
    algebra can be rebuilt -- its quiver presentation: the vertex list, the LABELED
    arrow list (so later sections can reference arrows by name), and the
    relations, rendered (not LaTeX source)."""
    algebra = result.get("algebra", {})
    rows = [f"quiverlab version: {result.get('quiverlab_version', '?')}",
            "algebra: " + json.dumps(algebra, default=str, sort_keys=True)]
    html = [f"<p>quiverlab version: {_esc(str(result.get('quiverlab_version', '?')))}</p>",
            "<p>algebra: <code>" + _esc(json.dumps(algebra, default=str, sort_keys=True))
            + "</code></p>"]
    A = _rebuild_algebra(result)
    Q = getattr(A, "quiver", None) if A is not None else None
    if Q is not None:
        verts = ", ".join(str(v) for v in Q.vertices)
        arrows_txt = "; ".join(f"{n}: {s} -> {t}" for n, (s, t) in Q.arrows.items())
        rows.append("vertices: " + verts)
        rows.append("arrows: " + arrows_txt)
        html.append(f"<p><b>vertices:</b> {_esc(verts)}</p>")
        html.append("<p><b>arrows:</b> " + ";&ensp;".join(
            f"<code>{_esc(str(n))}</code>: {_esc(str(s))} &rarr; {_esc(str(t))}"
            for n, (s, t) in Q.arrows.items()) + "</p>")
        rels = [repr(r) for r in (getattr(A, "relations", None) or [])]
        if rels:
            rows.append("relations: " + "; ".join(rels))
            html.append("<p><b>relations:</b> " + ",&ensp;".join(
                _expr_html(r) for r in rels) + "</p>")
    return _Section("Computation", rows, html)


def _module_input_sections(result: dict) -> list:
    """The computed-on modules as full representations (per-arrow matrices),
    from the CLI-envelope echo written beside ``result_schema``. Absent on webapp
    results and pre-change files -- then no section is emitted."""
    out = []
    for key, title in (("module", "The module M"),
                       ("ext_target", "The Ext target N"),
                       ("tor_target", "The Tor target N (a left module)")):
        blk = result.get(key)
        if not (isinstance(blk, dict) and "dims" in blk):
            continue
        rows = [f"side: {blk.get('side')}"]
        html = [f"<p>side: {_esc(str(blk.get('side')))}</p>"]
        r, h = _repr_rows_html("M" if key == "module" else "N", blk)
        rows += r
        html += h
        out.append(_Section(title, rows, html))
    return out


def _projectives_injectives_section(result: dict):
    """The "projectives and injectives of A" section (Plan 30, Marco #4): rebuild
    the algebra from the result's algebra spec and describe each P_v / I_v by its
    dimension vector and Loewy (radical) layers. Descriptive-only: any failure to
    rebuild (an exotic spec, a version skew) skips the section rather than sinking
    the whole report. The simples S_v are omitted (Marco: obvious)."""
    A = _rebuild_algebra(result)
    if A is None:
        return None
    try:
        from quiverlab.trace.modules import algebra_objects
        from quiverlab.trace.render_text import factor_stack_text
        objects = algebra_objects(A)
    except Exception:
        return None
    if not objects:
        return None
    legend = ("Loewy layers are the semisimple slices rad^k X / rad^(k+1) X, "
              "listed top to bottom and separated by '|'; S_v^m means the simple "
              "S_v repeated m times. The simples S_v are omitted.")
    rows = [legend]
    for row in objects:
        v = row["vertex"]
        for sym in ("P", "I"):
            d = row[sym]
            stack = " | ".join(factor_stack_text(L) for L in d["layers"]) or "0"
            rows.append(f"{sym}_{v}: dim {d['dim']}  {_dimvec_str(d['dimvec'])}  "
                        f"Loewy: {stack}  (top {factor_stack_text(d['top'])}, "
                        f"soc {factor_stack_text(d['socle'])})")
    return _Section("The projectives and injectives of A", rows)


def build_sections(result: dict) -> list:
    sections = [_computation_section(result)]
    sections += _module_input_sections(result)
    pi = _projectives_injectives_section(result)
    if pi is not None:
        sections.append(pi)
    per_kind = (result.get("resources") or {}).get("per_kind") or {}
    for name, block in result.get("results", {}).items():
        if isinstance(block, dict):
            sec = _section_for(name, block)
            pk = per_kind.get(name)
            if isinstance(pk, dict) and pk.get("wall_ms") is not None:
                line = "resources: wall " + _fmt_ms(pk["wall_ms"])
                if pk.get("peak_rss_mib") is not None:
                    line += (f"; process peak RSS {pk['peak_rss_mib']} MiB "
                             "by the end of this step")
                sec.rows.append(line)
                if sec.html:
                    sec.html.append(f"<p class='dv'>{_esc(line)}</p>")
            sections.append(sec)
    refs = result.get("references", [])
    if refs:
        rows = []
        for e in refs:
            key = e.get("key") if isinstance(e, dict) else str(e)
            formatted = e.get("formatted", key) if isinstance(e, dict) else str(e)
            rows.append(f"[{key}] {formatted}")
        sections.append(_Section("References", rows))
    repro = result.get("reproduce")
    if repro:
        sections.append(_Section("Reproduce locally", repro.splitlines()))
    res = result.get("resources")
    if isinstance(res, dict) and res:
        sections.append(_resources_section(res))
    return sections


def _fmt_ms(ms: int) -> str:
    return f"{ms // 1000}.{(ms % 1000) // 100} s  ({ms} ms)"


def _resources_section(res: dict) -> _Section:
    """The resources footer (CLI-envelope ``resources``): wall time, peak memory,
    detected cores/RAM -- exact ints from the run, formatted without floats. The
    ``per_kind`` breakdown is rendered inside each section, not here."""
    rows = []
    ms = res.get("wall_ms")
    if ms is not None:
        rows.append("wall time: " + _fmt_ms(ms))
    if res.get("peak_rss_mib") is not None:
        rows.append(f"peak memory (RSS): {res['peak_rss_mib']} MiB")
    if res.get("cores_detected") is not None:
        rows.append(f"cpu cores available: {res['cores_detected']}")
    if res.get("ram_mib_detected") is not None:
        rows.append(f"ram detected: {res['ram_mib_detected']} MiB")
    for k in sorted(res):
        if k not in ("wall_ms", "peak_rss_mib", "cores_detected",
                     "ram_mib_detected", "per_kind"):
            rows.append(f"{k}: {res[k]}")
    return _Section("Resources used", rows)


# --------------------------------------------------------------------------- #
# Format renderers
# --------------------------------------------------------------------------- #

def render_text(result: dict) -> str:
    out = ["quiverlab report", "=" * 40, ""]
    for sec in build_sections(result):
        out.append(sec.title)
        out.append("-" * len(sec.title))
        out.extend(sec.rows)
        out.append("")
    return "\n".join(out).rstrip("\n") + "\n"


_HTML_STYLE = (
    "<style>body{font-family:sans-serif;max-width:52rem;margin:2rem auto}"
    "pre{background:#f4f4f4;padding:6px;overflow-x:auto}"
    "h2{border-bottom:1px solid #ccc}"
    "table.matrix{display:inline-table;border-collapse:collapse;margin:2px 6px;"
    "border-left:1px solid #444;border-right:1px solid #444;border-radius:6px}"
    "table.matrix td{padding:1px 9px;text-align:right}"
    # index headers + a white/grey grid: odd rows white, even rows grey, with a
    # faint column stripe on even columns so both indices trace easily.
    "table.matrix th{padding:1px 7px;text-align:right;font-weight:normal;"
    "font-size:85%;color:#888;background:#fafafa;border-bottom:1px solid #ddd}"
    "table.matrix tr>th:first-child{border-right:1px solid #ddd;border-bottom:0}"
    "table.matrix tr:nth-child(even) td{background:#f0f0f0}"
    "table.matrix td:nth-child(odd){box-shadow:inset 0 0 0 99px rgba(0,0,0,0.025)}"
    ".arrowmap{margin:4px 0;display:flex;align-items:center;gap:6px;flex-wrap:wrap}"
    "p.resterm{margin:2px 0}.dv{color:#666}"
    "</style>")


def render_html(result: dict) -> str:
    body = ["<!doctype html><html><head><meta charset='utf-8'>", _HTML_STYLE,
            "<title>quiverlab report</title></head><body>",
            "<h1>quiverlab report</h1>",
            "<p><i>This report is print-ready and self-contained (no JavaScript); "
            "use your browser's Print &rarr; Save as PDF for a typeset copy.</i></p>"]
    for sec in build_sections(result):
        body.append("<h2>%s</h2>" % _esc(sec.title))
        if sec.html:
            body.extend(sec.html)
        elif sec.rows:
            body.append("<pre>" + _esc("\n".join(sec.rows)) + "</pre>")
    body.append("</body></html>")
    return "\n".join(body) + "\n"


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

# PDF/TeX report output has been removed; only these formats are produced.
_REMOVED_FORMATS = frozenset({"pdf", "tex"})


def default_out_name(fmt: str) -> str:
    # ``json`` is the worked-steps EVENT STREAM (trace.json), not a rendered report,
    # so it keeps that name rather than ``report.json``.
    if fmt == "json":
        return "trace.json"
    return "report." + {"html": "html", "txt": "txt"}.get(fmt, "txt")


def _sibling_trace_json(src_path) -> str:
    """The worked-steps machine record (``trace.json``) that the writer produced
    beside ``result.json`` (Plan 34). ``--format json`` emits it VERBATIM -- the same
    bytes the writer wrote, promoted by the hpc spec -- rather than re-rendering the
    result envelope, since the JSON deliverable IS the complete event stream. Raises
    ``ReportError`` when the source was an in-memory dict (no directory to look in) or
    no trace.json exists (no worked-steps computation ran)."""
    if src_path is None:
        raise ReportError(
            "--format json emits the worked-steps event stream (trace.json), which "
            "lives beside result.json; pass the result.json path, not an in-memory dict")
    tj = pathlib.Path(src_path).parent / "trace.json"
    try:
        return tj.read_text(encoding="utf-8")
    except OSError:
        raise ReportError(
            f"no trace.json beside {src_path} (the worked-steps event stream is "
            "produced only when a traced computation ran); use --format "
            "html/txt for the result report")


def render(result, out_path=None, fmt: str = "auto", on_warn=None):
    """Render ``result`` (a dict or a path to result.json) to ``out_path``.

    ``fmt`` in {auto, html, txt, json}. ``auto`` writes the self-contained HTML
    report; ``txt`` writes the plain-text report; ``json`` emits the worked-steps
    event stream (trace.json) that lives beside the result (Plan 34). PDF/TeX report
    output has been removed -- ``fmt="pdf"``/``"tex"`` raises ``ReportError``. Returns
    ``(pathlib.Path, actual_fmt)``."""
    # Remember the source PATH before loading (``json`` reads the sibling trace.json).
    src_path = None if isinstance(result, dict) else result
    if not isinstance(result, dict):
        result = load_result(result)
    check_result_schema(result)
    warn = _version_warning(result)
    if warn and on_warn is not None:
        on_warn(warn)

    if fmt in _REMOVED_FORMATS:
        raise ReportError("PDF/TeX output has been removed; use html or json "
                          f"(requested {fmt!r})")
    actual = "html" if fmt == "auto" else fmt

    out = pathlib.Path(out_path) if out_path is not None else pathlib.Path(default_out_name(actual))
    try:
        if actual == "txt":
            out.write_text(render_text(result), encoding="utf-8")
        elif actual == "html":
            out.write_text(render_html(result), encoding="utf-8")
        elif actual == "json":
            out.write_text(_sibling_trace_json(src_path), encoding="utf-8")
        else:
            raise ReportError(
                f"unknown format {fmt!r} (expected auto/html/txt/json)")
    except OSError as exc:
        raise ReportWriteError(f"cannot write report to {out}: {exc}")
    return out, actual
