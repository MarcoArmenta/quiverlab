"""Plan 34 (post-critique) -- MAJOR-4 (GF(p^n) entries render readably + a
``display_only`` flag, never the internal tuple ``(0, 0)``), MINOR-5 (matrix
display), MINOR-6 (an old-shape cached rad/top/soc says "recompute" instead of a
fabricated zero action) and MINOR-7 wiring.

MINOR-5 has since been REVERSED twice by Marco: first from "elide oversized
matrices" to "show them in full inside a scroll box" (Plan 34, mid-flight), then
on 2026-07-29 to "show them in full with NO scrollbar" -- an over-wide matrix is
shrunk to the column width instead. The tests at the bottom of this file pin the
current contract plus the rest of that 2026-07-29 pass.

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
    proc = subprocess.run([node, "-e", harness], capture_output=True, text=True,
                          encoding="utf-8", timeout=30)
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
def test_js_shows_matrices_complete_without_a_scrollbar(path):
    """Marco 2026-07-29 (supersedes the MINOR-5 scroll box): a matrix must be visible
    COMPLETE, never behind a scrollbar. `mathFit` clips nothing; an over-wide matrix
    is SHRUNK to the column width by the post-typeset `fitMath` pass instead."""
    src = path.read_text(encoding="utf-8")
    assert "mathScroll" not in src                  # the scroll box is gone
    assert "overflowX" not in src                   # ...and with it every clip
    fit = _grab_fn(src, "fitMath")
    assert "scale(" in fit and "scrollWidth" in fit  # shrink-to-fit, measured
    # Matrices themselves are INDEXED GRIDS (Marco 2026-07-29), not typeset
    # pmatrices, so nothing about them can clip either.
    assert "matrixGrid(maps[a])" in src            # used by the rad/top/soc renderer


@pytest.mark.parametrize("path", [_APP_JS, _GUI_JS])
def test_js_hides_arrows_acting_as_zero(path):
    """Marco 2026-07-29: an arrow acting as the exact zero map carries no
    information, so its matrix is not printed (the socle of his example showed a
    2x2 zero block for arrow d). The predicate is exact-string, never numeric."""
    src = path.read_text(encoding="utf-8")
    harness = "\n".join([
        _grab_fn(src, "matIsZero"),
        "const z = [['0','0'],['0','0']];",
        "const nz = [['0','0'],['0','1']];",
        "const frac = [['0','1/2']];",
        "process.stdout.write(JSON.stringify({"
        "  zero: matIsZero(z), nonzero: matIsZero(nz), frac: matIsZero(frac),"
        "  empty: matIsZero([])"
        "}));",
    ])
    got = _node(harness)
    assert got["zero"] is True and got["nonzero"] is False
    assert got["frac"] is False                     # a rational entry is not zero
    assert got["empty"] is True                     # no arrows -> "acts as zero"


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


# --------------------------------------------------------------------------- #
# Marco 2026-07-29: the two renderers show the AR translates as full modules, the
# resolutions' differentials, and never repeat an identical differential.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("path", [_APP_JS, _GUI_JS])
def test_js_renders_the_ar_translate_as_a_full_module(path):
    """The translate IS a module: its per-arrow action matrices are shown, and the
    SECOND module N's translate is shown too when the request named one."""
    src = path.read_text(encoding="utf-8")
    assert "t.repr" in src                       # the translate's {dims, maps}
    assert "b.targets" in src or "block.targets" in src
    assert "ext_target" in src and "tor_target" in src


@pytest.mark.parametrize("path", [_APP_JS, _GUI_JS])
def test_js_renders_resolution_differentials_without_repeats(path):
    src = path.read_text(encoding="utf-8")
    assert "differentials" in src                # the block field is rendered...
    assert "not repeated" in src                 # ...and repeats are referenced
    assert "rows: target basis, columns: source basis" in src


def test_gui_js_never_typesets_a_missing_pd_id_latex():
    """An older cached pd/id block has no `latex`; composing it from `value` is what
    stops the literal "undefined" Marco saw twice in example-a."""
    src = _GUI_JS.read_text(encoding="utf-8")
    assert "homdimLatex" in src
    harness = "\n".join([
        _grab_fn(src, "homdimLatex"),
        "process.stdout.write(JSON.stringify({"
        "  finite: homdimLatex('projective_dimension', {value: 3}),"
        "  unresolved: homdimLatex('injective_dimension', {value: null, bound: 32}),"
        "  legacy: homdimLatex('projective_dimension', {})"
        "}));",
    ])
    got = _node(harness)
    assert got["finite"] == "\\operatorname{pd} M = 3"
    assert got["unresolved"] == "\\operatorname{id} M > 32"
    assert got["legacy"] == "\\operatorname{pd} M > 32"     # never "undefined"


def test_the_2026_07_29_strings_are_bilingual_and_wired():
    """The webapp pages are EN/ES; every string the 2026-07-29 pass added to the
    module blocks has both translations and a template data-attribute, so the
    Spanish page does not silently degrade to English."""
    en = json.loads(_EN.read_text(encoding="utf-8"))
    es = json.loads(_ES.read_text(encoding="utf-8"))
    html = _INDEX.read_text(encoding="utf-8")
    src = _APP_JS.read_text(encoding="utf-8")
    for key, attr, prop in (
            ("mod.arrows_acting_zero", "data-mod-arrows-acting-zero", "modArrowsActingZero"),
            ("mod.differentials_proj", "data-mod-differentials-proj", "modDifferentialsProj"),
            ("mod.differentials_inj", "data-mod-differentials-inj", "modDifferentialsInj"),
            ("mod.same_matrix", "data-mod-same-matrix", "modSameMatrix"),
            ("mod.matrix_too_large", "data-mod-matrix-too-large", "modMatrixTooLarge"),
            ("mod.tau_target_ext", "data-mod-tau-target-ext", "modTauTargetExt"),
            ("mod.tau_target_tor", "data-mod-tau-target-tor", "modTauTargetTor"),
            ("mod.tau_unavailable", "data-mod-tau-unavailable", "modTauUnavailable"),
            ("mod.arrows_zero", "data-mod-arrows-zero", "modArrowsZero"),
            ("mod.rad_top_soc", "data-mod-rad-top-soc", "modRadTopSoc")):
        assert en.get(key), "missing EN %s" % key
        assert es.get(key), "missing ES %s" % key
        assert en[key] != es[key], "EN/ES not translated for %s" % key
        assert attr in html and key in html, "template does not wire %s" % key
        assert prop in src, "app.js does not read %s" % prop


def test_js_names_standard_summands_and_shows_the_rest(tmp_path):
    """Marco 2026-07-29: a Krull-Schmidt summand isomorphic to S_v / P_v / I_v is
    NAMED (its matrices would be noise); every other summand is shown in full."""
    for path in (_APP_JS, _GUI_JS):
        src = path.read_text(encoding="utf-8")
        harness = "\n".join([
            _grab_const(src, "STD_SYM"),
            _grab_fn(src, "summandName"),
            "process.stdout.write(JSON.stringify({"
            "  simple: summandName({standard:{kind:'simple',vertex:'2'}}, 1),"
            "  proj:   summandName({standard:{kind:'projective',vertex:'1'}}, 2),"
            "  inj:    summandName({standard:{kind:'injective',vertex:'3'}}, 3),"
            "  plain:  summandName({maps:{}}, 4),"
            "  unknown: summandName({standard:{kind:'weird',vertex:'9'}}, 5)"
            "}));",
        ])
        got = _node(harness)
        assert got["simple"] == "S_2" and got["proj"] == "P_1" and got["inj"] == "I_3"
        assert got["plain"] == "M_4"
        assert got["unknown"] == "M_5", "an unrecognised kind must not be named"
        # the non-standard summands' matrices are rendered
        assert "appendSummandMaps" in src or "summandName(s, i + 1), s, d" in src


@pytest.mark.parametrize("path", [_APP_JS, _GUI_JS])
def test_js_matrices_are_indexed_grids(path):
    """Marco 2026-07-29: every matrix gets an extra row and column of indices and a
    light grid, so an entry can be read off by position. 1-based."""
    src = path.read_text(encoding="utf-8")
    assert "ql-matrix" in src or "qlgui-matrix" in src        # the grid class
    assert "ql-corner" in src or "qlgui-corner" in src        # the header corner
    fn = _grab_fn(src, "matrixGrid")
    assert "j + 1" in fn and "i + 1" in fn                    # 1-based indices
    assert "matTooBig" in fn                                  # the sanity backstop holds
