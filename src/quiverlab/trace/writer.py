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


def write_trace(events, table, algebra, kind, top, references=(), out_dir=None):
    out = pathlib.Path(out_dir) if out_dir is not None else (pathlib.Path.cwd() / "quiverlab_traces")
    out.mkdir(parents=True, exist_ok=True)
    stem = "%s_%s" % (_SAFE_STEM.get(kind, kind), _hash(algebra, kind, top))
    title = "%s of %s" % (kind, repr(algebra).splitlines()[0])
    engine = have_latex()
    html_note = "no LaTeX toolchain found -- install pdflatex or tectonic for a PDF"
    if engine is not None:
        pdf = out / (stem + ".pdf")
        try:
            pages = _compile_pdf(render_latex(events, title=title, references=references),
                                 str(pdf), engine)
            print("Worked steps: %s (%d pp)" % (_rel(pdf), pages))
            return str(pdf)
        except Exception:
            # a toolchain WAS found but the compile failed: say so honestly; never
            # claim "no toolchain found" when one is on PATH.
            html_note = "LaTeX compilation failed (%s); wrote HTML fallback" % engine
    html = out / (stem + ".html")
    html.write_text(render_html(events, title=title, references=references))
    print("Worked steps: %s (HTML, no JavaScript; %s)" % (_rel(html), html_note))
    return str(html)


def _rel(p):
    try:
        return str(p.relative_to(pathlib.Path.cwd()))
    except ValueError:
        return str(p)
