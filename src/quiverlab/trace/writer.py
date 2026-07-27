"""Renderer selection + output-path contract + the printed one-liner (spec §3.8).

Selection: pdflatex or tectonic on PATH -> compile LaTeX to PDF; otherwise (or if a
found toolchain fails to compile) write the self-contained no-JS HTML (TeX source)
and print a LOUD, HONEST one-liner distinguishing "no toolchain found" from
"compilation failed". Output:
./quiverlab_traces/HHc_<hash>.<ext> (cohomology) / HHh_<hash>.<ext> (homology)
(Plan 09 collects the newest *.pdf, else *.html, from this directory -- the glob is
extension-based, so the safe stem does not affect it).

The filename hash is 12 hex chars (48 bits): the plan's original 4-hex (16-bit)
stem would collide by the birthday bound at ~256 distinct traces and silently
overwrite; 12 hex pushes the collision horizon out of practical reach while staying
fully deterministic (no floats)."""
import hashlib
import pathlib
import shutil
import subprocess
import tempfile

from quiverlab.trace.render_latex import render_latex
from quiverlab.trace.render_html import render_html
from quiverlab.trace.render_json import render_json

# Filesystem-safe filename stems for the caret-bearing kinds (no "^" in a filename).
_SAFE_STEM = {"HH^": "HHc", "HH_": "HHh"}


def have_latex():
    for engine in ("tectonic", "pdflatex"):
        if shutil.which(engine):
            return engine
    return None


def _hash(algebra, kind, top):
    h = hashlib.sha1(("%s|%s|%s" % (repr(algebra), kind, top)).encode("utf-8"))
    return h.hexdigest()[:12]


def _needed_matrix_cols(events):
    """The widest rendered ``pmatrix`` (a differential) across the events, used to
    lift amsmath's column ceiling. Elided matrices render as ``\\text{...}`` (no
    columns), so only non-elided matrix events contribute."""
    best = 0
    for e in events:
        if getattr(e, "matrix", None) is not None and not getattr(e, "elided", False):
            best = max(best, getattr(e, "ncols", 0) or 0)
    return best


def _widen_matrix_cols(tex, ncols):
    r"""amsmath's ``pmatrix`` caps columns at ``MaxMatrixCols`` (default 10); a
    wider differential makes pdflatex/tectonic abort with "Extra alignment tab has
    been changed to \cr", so the PDF silently degrades to the HTML fallback *even
    with a toolchain present* (this was the "print report gives no PDF" bug). Raise
    the ceiling to the widest matrix actually rendered so BOTH the compiled PDF and
    the persisted, compile-it-yourself ``.tex`` build. Idempotent (skips if the
    counter is already set upstream) and anchored on ``\begin{document}`` -- present
    in every LaTeX document, so this survives preamble edits."""
    if ncols <= 10 or r"\setcounter{MaxMatrixCols}" in tex:
        return tex
    return tex.replace(
        r"\begin{document}",
        r"\setcounter{MaxMatrixCols}{%d}" % ncols + "\n" + r"\begin{document}", 1)


def _compile_pdf(tex, out_pdf, engine):
    """Compile `tex` to `out_pdf` with `engine`; return the page count (best effort)."""
    with tempfile.TemporaryDirectory() as d:
        src = pathlib.Path(d) / "trace.tex"
        src.write_text(tex)
        if engine == "tectonic":
            cmd = ["tectonic", "-o", d, str(src)]
        else:
            cmd = ["pdflatex", "-interaction=nonstopmode", "-output-directory", d, str(src)]
        subprocess.run(cmd, cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        built = pathlib.Path(d) / "trace.pdf"
        shutil.copyfile(built, out_pdf)
        data = built.read_bytes()
        return max(data.count(b"/Type /Page") - data.count(b"/Type /Pages"), 1)


def write_trace(events, table, algebra, kind, top, references=(), out_dir=None,
                meta_out=None):
    """Render the worked steps to PDF (a LaTeX toolchain on PATH) or the print-ready
    HTML fallback, persist the .tex alongside, print the one-liner, and return the
    produced path (str). ``meta_out`` (optional dict) is the split-deployment fix
    (Plan 34 MAJOR-3): the WORKER records HERE why a PDF was not produced --
    ``"no_toolchain"`` or ``"compile"`` -- so ``pages.py`` reports the RECORDED reason
    instead of re-probing ``have_latex()`` at page-render time on a different tier.
    Set only when the HTML fallback is taken; untouched (so absent) when a PDF is
    produced -- keeping the PDF-present result dict byte-identical."""
    events = list(events)
    out = pathlib.Path(out_dir) if out_dir is not None else (pathlib.Path.cwd() / "quiverlab_traces")
    out.mkdir(parents=True, exist_ok=True)
    stem = "%s_%s" % (_SAFE_STEM.get(kind, kind), _hash(algebra, kind, top))
    title = "%s of %s" % (kind, repr(algebra).splitlines()[0])
    # Build the LaTeX source ONCE and persist it next to the produced pdf/html:
    # the .tex is the exhaustive, downloadable-everywhere source (Plan 30 C1), so it
    # is written whether or not a toolchain is present -- persisted, not
    # temp-discarded by the compile step.
    tex = render_latex(events, title=title, references=references, algebra=algebra)
    # Lift amsmath's 10-column pmatrix ceiling to the widest matrix actually
    # rendered, or a >10-column differential aborts the compile (PDF -> HTML
    # fallback even with a toolchain present). Applied to the persisted .tex too,
    # so the "compile it yourself" source builds standalone.
    tex = _widen_matrix_cols(tex, _needed_matrix_cols(events))
    (out / (stem + ".tex")).write_text(tex)
    # The JSON machine record (Plan 34): the complete event stream, deterministic and
    # schema-versioned, written beside the .tex whether or not a toolchain is present
    # -- the third mandated artifact (PDF/HTML + this). A pure function of the events,
    # so it is byte-identical for identical input regardless of the PDF/HTML branch.
    (out / (stem + ".json")).write_text(
        render_json(events, title=title, references=references, algebra=algebra))
    engine = have_latex()
    html_note = "no LaTeX toolchain found -- install pdflatex or tectonic for a PDF"
    reason = "no_toolchain"                       # why no PDF (recorded for pages.py)
    if engine is not None:
        pdf = out / (stem + ".pdf")
        try:
            pages = _compile_pdf(tex, str(pdf), engine)
            print("Worked steps: %s (%d pp)" % (_rel(pdf), pages))
            return str(pdf)
        except Exception:
            # a toolchain WAS found but the compile failed: say so honestly; never
            # claim "no toolchain found" when one is on PATH.
            html_note = "LaTeX compilation failed (%s); wrote HTML fallback" % engine
            reason = "compile"
    html = out / (stem + ".html")
    html.write_text(render_html(events, title=title, references=references, algebra=algebra))
    if meta_out is not None:
        meta_out["pdf_fallback_reason"] = reason
    print("Worked steps: %s (HTML, no JavaScript; %s)" % (_rel(html), html_note))
    return str(html)


def _rel(p):
    try:
        return str(p.relative_to(pathlib.Path.cwd()))
    except ValueError:
        return str(p)
