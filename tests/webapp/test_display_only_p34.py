"""Plan 34 (post-critique) -- MAJOR-4 (GF(p^n) entries render readably + a
``display_only`` flag, never the internal tuple ``(0, 0)``), MINOR-5 (the GUI JS
elides oversized matrices), MINOR-6 (an old-shape cached rad/top/soc says "recompute"
instead of a fabricated zero action) and MINOR-7 wiring.

The Python serializer (``module_blocks``) is tested directly over GF(4). The JS
helpers (pure) are extracted textually and evaluated under node (skipped when node is
unavailable, as the sibling JS tests do); the DOM-side notes are pinned by static
scan + i18n/template wiring, in the idiom of tests/webapp/test_pages.py.
"""
import json
import pathlib
import re
import shutil
import subprocess

import pytest

from quiverlab import Quiver, GF
from quiverlab.modules.qpa_module import module_blocks, _json_entry, _reenterable

_ROOT = pathlib.Path(__file__).resolve().parents[2]
_APP_JS = _ROOT / "webapp" / "static" / "app.js"
_GUI_JS = _ROOT / "docs" / "gui" / "gui.js"
_INDEX = _ROOT / "webapp" / "templates" / "index.html"
_EN = _ROOT / "webapp" / "server" / "i18n" / "en.json"
_ES = _ROOT / "webapp" / "server" / "i18n" / "es.json"


# --------------------------------------------------------------------------- #
# MAJOR-4: module_blocks over GF(4) -- readable entries + the display_only flag.
# --------------------------------------------------------------------------- #
def _gf4_module():
    """A kA_2 representation over GF(4) whose arrow acts by the field generator x, so
    the serialized action matrix carries an extension-field element (not int/fraction)."""
    F = GF(4)
    x = F.gen()                                   # (0, 1), i.e. the generator x
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    A = Q.algebra(relations=[], field=F)
    return A.module({1: 1, 2: 1}, {"a": [[F.zero(), F.zero()], [x, F.zero()]]}), F


def test_gf4_entries_render_readably_not_as_internal_tuple():
    M, F = _gf4_module()
    blk = module_blocks(M)
    # the arrow entry is the domain's own readable rendering, NOT "(0, 1)":
    assert blk["maps"]["a"] == [["x^1"]], blk["maps"]["a"]
    assert "(0, 1)" not in json.dumps(blk), "internal coefficient tuple leaked into the block"
    # and the block is flagged display-only (the entry is not re-enterable):
    assert blk.get("display_only") is True


def test_json_entry_uses_domain_to_str_for_extension_elements():
    F = GF(4)
    assert _json_entry(F.gen(), F) == "x^1"        # readable
    assert _json_entry((1, 1), F) == "1 + x^1"     # 1 + x
    # and re-enterability is judged correctly:
    assert not _reenterable("x^1") and not _reenterable("(0, 1)")
    assert _reenterable(3) and _reenterable("1/2") and _reenterable("-2")


def test_gfp_and_prime_field_modules_stay_byte_clean_no_flag():
    """GF(p) / integer modules keep int entries and gain NO display_only key (the
    common case stays byte-identical)."""
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    A = Q.algebra(relations=[], field=GF(2))
    M = A.module({1: 1, 2: 1}, {"a": [[0, 0], [1, 0]]})
    blk = module_blocks(M)
    assert blk["maps"]["a"] == [[1]]
    assert "display_only" not in blk


# --------------------------------------------------------------------------- #
# The DOM-side notes exist and are wired for i18n (webapp) / hardcoded (gui docs).
# --------------------------------------------------------------------------- #
def test_i18n_has_display_only_and_stale_strings_en_es():
    en = json.loads(_EN.read_text(encoding="utf-8"))
    es = json.loads(_ES.read_text(encoding="utf-8"))
    for key in ("mod.display_only", "mod.stale_recompute"):
        assert key in en and en[key], "missing EN %s" % key
        assert key in es and es[key], "missing ES %s" % key
        assert en[key] != es[key], "EN/ES not translated for %s" % key


def test_index_template_wires_the_new_data_attributes():
    html = _INDEX.read_text(encoding="utf-8")
    assert "data-mod-display-only=" in html and "mod.display_only" in html
    assert "data-mod-stale=" in html and "mod.stale_recompute" in html


def test_app_js_renders_display_only_and_stale_notes():
    src = _APP_JS.read_text(encoding="utf-8")
    assert "d.modDisplayOnly" in src              # MAJOR-4 note (i18n)
    assert "d.modStale" in src                    # MINOR-6 note (i18n)
    assert "radTopSocStale" in src and "radTopSocDisplayOnly" in src


def test_gui_js_renders_notes_and_handles_blocked_popup():
    src = _GUI_JS.read_text(encoding="utf-8")
    assert "radTopSocStale" in src and "radTopSocDisplayOnly" in src
    assert "recompute" in src                     # MINOR-6 honest notice
    assert "display only" in src                  # MAJOR-4 note
    assert "popup blocked" in src                 # MINOR-7 visible message


# --------------------------------------------------------------------------- #
# The pure JS helpers under node: elision threshold, stale + display-only guards.
# --------------------------------------------------------------------------- #
def _grab_const(src, name):
    m = re.search(r"(?:var|const)\s+%s\s*=\s*[^;]+;" % re.escape(name), src)
    assert m, "constant %s not found" % name
    return m.group(0)


def _grab_fn(src, name):
    start = src.index("function %s" % name)
    depth, started, i = 0, False, start
    while i < len(src):
        c = src[i]
        if c == "{":
            depth += 1
            started = True
        elif c == "}":
            depth -= 1
            if started and depth == 0:
                return src[start:i + 1]
        i += 1
    raise AssertionError("function %s not found / unbalanced" % name)


def _node(harness):
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not available")
    proc = subprocess.run([node, "-e", harness], capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


@pytest.mark.parametrize("path", [_APP_JS, _GUI_JS])
def test_js_matlatex_renders_full_matrices_with_only_a_backstop(path):
    """Plan 34 (Marco, mid-flight): matrices are NOT elided at a small size -- they are
    shown IN FULL (and scroll-wrapped, see below). Only a very large sanity backstop
    trips the shape-only note, so a pathological payload cannot hang the browser."""
    src = path.read_text(encoding="utf-8")
    harness = "\n".join([
        _grab_const(src, "MAT_BACKSTOP_CELLS"),
        _grab_fn(src, "matTooBig"),
        _grab_fn(src, "matLatex"),
        "const m30 = Array.from({length: 30}, () => Array.from({length: 30}, () => 1));",
        # 600x600 = 360000 cells > the 250000 backstop (mirrors recorder MATRIX_ELISION_CELLS)
        "const huge = Array.from({length: 600}, () => Array.from({length: 600}, () => 1));",
        "const ok = [[1,0],[0,1]];",
        "process.stdout.write(JSON.stringify({"
        "  m30big: matTooBig(m30), hugeBig: matTooBig(huge), okBig: matTooBig(ok),"
        "  m30full: matLatex(m30).indexOf('\\\\begin{pmatrix}') === 0,"
        "  hugeBackstop: matLatex(huge).indexOf('backstop') >= 0,"
        "  okFull: matLatex(ok).indexOf('backstop') < 0"
        "}));",
    ])
    got = _node(harness)
    assert got["m30big"] is False                  # 30x30 (900 cells) renders in full
    assert got["m30full"] is True                  # ...as a real pmatrix, not elided
    assert got["hugeBig"] is True                  # 400x400 (160k cells) trips the backstop
    assert got["hugeBackstop"] is True             # ...stated by shape
    assert got["okBig"] is False and got["okFull"] is True


@pytest.mark.parametrize("path", [_APP_JS, _GUI_JS])
def test_js_wraps_matrices_in_a_horizontal_scroll_box(path):
    """MINOR-5 (revised): matrices are wrapped in an overflow-x:auto box so a wide
    differential scrolls INSIDE its box and the page body never scrolls sideways."""
    src = path.read_text(encoding="utf-8")
    fn = _grab_fn(src, "mathScroll")
    assert "overflowX" in fn and "auto" in fn
    assert "maxWidth" in fn                         # box is capped to the container width
    assert "mathScroll(matLatex(" in src           # used by the rad/top/soc renderer


@pytest.mark.parametrize("path", [_APP_JS, _GUI_JS])
def test_js_stale_and_display_only_guards(path):
    src = path.read_text(encoding="utf-8")
    harness = "\n".join([
        _grab_fn(src, "radTopSocStale"),
        _grab_fn(src, "radTopSocDisplayOnly"),
        # old-shape cached result: views carry a dimvec but no {dims, maps}
        "const old = {radical:{dimvec:{}}, top:{dimvec:{}}, socle:{dimvec:{}}};",
        "const fresh = {radical:{dims:{},maps:{}}, top:{dims:{},maps:{}}, socle:{dims:{},maps:{}}};",
        "const disp = {radical:{dims:{},maps:{},display_only:true},"
        "             top:{dims:{},maps:{}}, socle:{dims:{},maps:{}}};",
        "process.stdout.write(JSON.stringify({"
        "  staleOld: radTopSocStale(old), staleFresh: radTopSocStale(fresh),"
        "  doDisp: radTopSocDisplayOnly(disp), doFresh: radTopSocDisplayOnly(fresh)"
        "}));",
    ])
    got = _node(harness)
    assert got["staleOld"] is True                 # MINOR-6: old shape -> recompute
    assert got["staleFresh"] is False
    assert got["doDisp"] is True                   # MAJOR-4: display-only detected
    assert got["doFresh"] is False
