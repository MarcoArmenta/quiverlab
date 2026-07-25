"""Plan 28 -- render a ``result.json`` into a human report (LaTeX PDF / HTML /
plain text), reusing ``quiverlab.trace``'s LaTeX escaping + engine detection.

The engine ladder mirrors ``quiverlab.trace.writer``: tectonic/pdflatex on PATH
-> compile a standalone article to PDF; otherwise a self-contained no-JS HTML
(math shown as TeX source). ``--format txt`` always works with no toolchain.

The envelope carries ``result_schema`` (an int; absent == legacy 1). ``render``
refuses a NEWER ``result_schema`` loudly (the CLI maps that to exit 65) and warns
on library-version skew between the result and the running quiverlab."""
from __future__ import annotations

import json
import pathlib
import shutil
import subprocess
import tempfile

import quiverlab as ql
from quiverlab.hpc.spec import RESULT_SCHEMA
from quiverlab.trace.render_latex import _tex_escape
from quiverlab.trace.writer import have_latex


class ReportError(Exception):
    """The report could not be produced (unsupported result_schema, unwritable
    output, or a required LaTeX toolchain that was explicitly requested)."""


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
# Intermediate: result.json -> a list of neutral sections
# --------------------------------------------------------------------------- #

class _Section:
    __slots__ = ("title", "rows", "math")

    def __init__(self, title, rows=None, math=None):
        self.title = title
        self.rows = rows or []          # plain-text lines
        self.math = math or []          # TeX-source strings (displayed math)


def _dimvec_str(dv: dict) -> str:
    return "(" + ", ".join(f"{k}: {v}" for k, v in dv.items()) + ")"


def _hh_section(kind_label, block) -> _Section:
    dims = block.get("dims", [])
    rows = [f"engine: {block.get('engine', '?')}",
            "dimensions: " + ", ".join(f"{kind_label}{i} = {d}"
                                       for i, d in enumerate(dims))]
    math = [r",\quad ".join(r"%s{%d} = %d" % (kind_label, i, d)
                            for i, d in enumerate(dims))] if dims else []
    return _Section(f"Hochschild {'cohomology' if kind_label == 'HH^' else 'homology'}",
                    rows, math)


def _module_view_row(label, view) -> str:
    return f"{label}: dim {view.get('dim')}  {_dimvec_str(view.get('dimvec', {}))}"


def _section_for(name: str, block: dict) -> _Section:
    if name == "hh_cohomology":
        return _hh_section("HH^", block)
    if name == "hh_homology":
        return _hh_section("HH_", block)
    if name == "cartan":
        return _Section("Cartan matrix",
                        ["rows: " + "; ".join(str(r) for r in block.get("matrix", []))],
                        [block["latex"]] if block.get("latex") else [])
    if name == "coxeter_polynomial":
        return _Section("Coxeter polynomial",
                        [block.get("text", "")],
                        [block["latex"]] if block.get("latex") else [])
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
        return _Section("Radical / top / socle",
                        [f"side: {block.get('side')}",
                         _module_view_row("rad M", block.get("radical", {})),
                         _module_view_row("top M", block.get("top", {})),
                         _module_view_row("soc M", block.get("socle", {}))])
    if name in ("tau", "tau_minus"):
        sym = "tau M" if name == "tau" else "tau^- M"
        rows = [f"side: {block.get('side')}"]
        rows.append(f"{sym} = 0" if block.get("is_zero")
                    else _module_view_row(sym, block))
        return _Section("Auslander-Reiten translate", rows)
    if name == "ext":
        dims = block.get("dims", [])
        rows = ["target " + _module_view_row("N", block.get("target", {})),
                "Ext dims: " + ", ".join(f"Ext^{i} = {d}" for i, d in enumerate(dims))]
        return _Section("Ext^*(M, N)", rows)
    if name in ("projective_resolution", "injective_resolution"):
        terms = block.get("terms", [])
        rows = [f"top: {block.get('top')}"]
        for i, dv in enumerate(terms):
            rows.append(f"  P_{i}: {_dimvec_str(dv)}  (betti {block.get('betti', [None])[i] if i < len(block.get('betti', [])) else '?'})")
        if "pd" in block:
            rows.append(f"projective dimension: {block['pd']}")
        if "injective_dimension" in block:
            rows.append(f"injective dimension: {block['injective_dimension']}")
        label = ("Projective resolution" if name == "projective_resolution"
                 else "Injective resolution")
        return _Section(label, rows)
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
    # Unknown block: dump its keys honestly rather than dropping it.
    return _Section(name, [json.dumps(block, default=str, sort_keys=True)])


def build_sections(result: dict) -> list:
    algebra = result.get("algebra", {})
    head_rows = [f"quiverlab version: {result.get('quiverlab_version', '?')}",
                 "algebra: " + json.dumps(algebra, default=str, sort_keys=True)]
    sections = [_Section("Computation", head_rows)]
    for name, block in result.get("results", {}).items():
        if isinstance(block, dict):
            sections.append(_section_for(name, block))
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
    return sections


# --------------------------------------------------------------------------- #
# Format renderers
# --------------------------------------------------------------------------- #

def render_text(result: dict) -> str:
    out = ["quiverlab report", "=" * 40, ""]
    for sec in build_sections(result):
        out.append(sec.title)
        out.append("-" * len(sec.title))
        out.extend(sec.rows)
        for m in sec.math:
            out.append("    " + m)
        out.append("")
    return "\n".join(out).rstrip("\n") + "\n"


_HTML_STYLE = ("<style>body{font-family:sans-serif;max-width:52rem;margin:2rem auto}"
               "pre{background:#f4f4f4;padding:6px;overflow-x:auto}"
               "h2{border-bottom:1px solid #ccc}</style>")


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def render_html(result: dict) -> str:
    body = ["<!doctype html><html><head><meta charset='utf-8'>", _HTML_STYLE,
            "<title>quiverlab report</title></head><body>",
            "<h1>quiverlab report</h1>",
            "<p><i>Math is shown as TeX source (no JavaScript); compile a PDF with "
            "pdflatex/tectonic for typeset output.</i></p>"]
    for sec in build_sections(result):
        body.append("<h2>%s</h2>" % _esc(sec.title))
        if sec.rows:
            body.append("<pre>" + _esc("\n".join(sec.rows)) + "</pre>")
        for m in sec.math:
            body.append("<pre><code>%s</code></pre>" % _esc(m))
    body.append("</body></html>")
    return "\n".join(body) + "\n"


def render_latex(result: dict) -> str:
    out = [r"\documentclass{article}", r"\usepackage{amsmath}",
           r"\usepackage[T1]{fontenc}", r"\begin{document}",
           r"\section*{quiverlab report}"]
    for sec in build_sections(result):
        out.append(r"\subsection*{%s}" % _tex_escape(sec.title))
        for row in sec.rows:
            out.append(r"\noindent %s\\" % _tex_escape(row))
        for m in sec.math:
            out.append(r"\[ %s \]" % m)
    out.append(r"\end{document}")
    return "\n".join(out) + "\n"


def _compile_pdf(tex: str, out_pdf: pathlib.Path, engine: str) -> None:
    with tempfile.TemporaryDirectory() as d:
        src = pathlib.Path(d) / "report.tex"
        src.write_text(tex, encoding="utf-8")
        if engine == "tectonic":
            cmd = ["tectonic", "-o", d, str(src)]
        else:
            cmd = ["pdflatex", "-interaction=nonstopmode", "-output-directory", d, str(src)]
        subprocess.run(cmd, cwd=d, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=True)
        built = pathlib.Path(d) / "report.pdf"
        shutil.copyfile(built, out_pdf)


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def default_out_name(fmt: str) -> str:
    return "report." + {"pdf": "pdf", "html": "html", "txt": "txt"}.get(fmt, "txt")


def render(result, out_path=None, fmt: str = "auto", on_warn=None):
    """Render ``result`` (a dict or a path to result.json) to ``out_path``.

    ``fmt`` in {auto, pdf, html, txt}. ``auto`` compiles a PDF when tectonic/
    pdflatex is on PATH, else writes HTML. ``pdf`` without a toolchain raises
    ``ReportError``. Returns ``(pathlib.Path, actual_fmt)``."""
    if not isinstance(result, dict):
        result = load_result(result)
    check_result_schema(result)
    warn = _version_warning(result)
    if warn and on_warn is not None:
        on_warn(warn)

    engine = have_latex()
    actual = fmt
    if fmt == "auto":
        actual = "pdf" if engine else "html"
    if actual == "pdf" and engine is None:
        raise ReportError("no LaTeX toolchain (tectonic/pdflatex) on PATH; use "
                          "--format html or --format txt")

    out = pathlib.Path(out_path) if out_path is not None else pathlib.Path(default_out_name(actual))
    try:
        if actual == "txt":
            out.write_text(render_text(result), encoding="utf-8")
        elif actual == "html":
            out.write_text(render_html(result), encoding="utf-8")
        elif actual == "pdf":
            try:
                _compile_pdf(render_latex(result), out, engine)
            except (subprocess.SubprocessError, OSError) as exc:
                raise ReportError(f"LaTeX compilation failed ({engine}): {exc}; "
                                  "use --format html or --format txt")
        else:
            raise ReportError(f"unknown format {fmt!r} (expected auto/pdf/html/txt)")
    except OSError as exc:
        raise ReportWriteError(f"cannot write report to {out}: {exc}")
    return out, actual
