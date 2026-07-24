"""Task 11 -- server-rendered bilingual pages + artifact downloads.

Adapted from the brief against the REAL in-repo interfaces (recorded in
task-11-report.md):
  * the catalog family names are CapitalCase library classes (``NakayamaAlgebra``
    etc.); ``truncated_polynomial`` is NOT a family (same finding as Task 9's
    test_api.py), so the "catalog injected" witness asserts a real family name.
  * the runner's HTML worked-steps fallback is ``trace_steps.html`` (NOT
    ``trace.html``); the download whitelist matches the runner's real names.
  * job ids are validated with Task 9's ``valid_ulid`` helper.

Beyond the brief's page tests this file also pins the adjudicated hardening:
download traversal + whitelist + done-gate, no untranslated key leak in /es
(brief keys AND a key-shaped regex over visible text), and a no-inline-script
assertion on every rendered page (strict CSP: every <script> must carry src=).
"""
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.store import JobStore


def _client(tmp_path):
    return TestClient(create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})))


def _done_job(cfg, result_json: str):
    """Create -> claim -> mark_done a job with the given result.json bytes."""
    store = JobStore(cfg.db_path)
    store.init_schema()
    jid = store.create_job({}, ip="h")
    store.claim_next()
    art = cfg.artifacts_dir / jid
    art.mkdir(parents=True, exist_ok=True)
    (art / "result.json").write_text(result_json)
    store.mark_done(jid, str(art))
    return jid


# --------------------------------------------------------------------------- #
# Brief's page tests (family-name witness adapted to the real catalog).
# --------------------------------------------------------------------------- #
def test_index_renders_form(tmp_path):
    r = _client(tmp_path).get("/")
    assert r.status_code == 200
    assert "quiverlab" in r.text.lower()
    assert "NakayamaAlgebra" in r.text     # catalog injected (real family name)


def test_about_renders(tmp_path):
    r = _client(tmp_path).get("/about")
    assert r.status_code == 200
    assert "pip install quiverlab" in r.text


def test_job_page_404_for_malformed_id(tmp_path):
    r = _client(tmp_path).get("/job/NOSUCHJOB")     # not a ULID -> 404 before store
    assert r.status_code == 404


def test_job_page_404_for_wellformed_but_absent(tmp_path):
    r = _client(tmp_path).get("/job/01AN4Z07BY79KA1307SR9X4MV3")  # valid ULID, absent
    assert r.status_code == 404


def test_download_404_for_unknown(tmp_path):
    r = _client(tmp_path).get("/download/NOSUCHJOB/result.json")
    assert r.status_code == 404


def test_done_job_page_renders_downloads_and_reproduce(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid = _done_job(cfg, '{"reproduce": "import quiverlab as ql", '
                         '"quiverlab_version": "9.9.9"}')
    r = TestClient(create_app(cfg)).get("/job/" + jid)
    assert r.status_code == 200
    assert "Reproduce locally" in r.text
    assert "result.json" in r.text
    assert "9.9.9" in r.text


def test_es_index_is_spanish(tmp_path):
    r = _client(tmp_path).get("/es")
    assert r.status_code == 200
    assert "Polinomio de Coxeter" in r.text     # translated chrome
    assert "Coxeter polynomial" not in r.text   # English string does not leak


def test_no_untranslated_key_leaks_in_es(tmp_path):
    from webapp.server.i18n import all_keys
    r = _client(tmp_path).get("/es")
    for key in all_keys():
        assert key not in r.text, f"untranslated key leaked into /es HTML: {key}"


def test_language_toggle_hrefs(tmp_path):
    c = _client(tmp_path)
    assert 'href="/es"' in c.get("/").text
    assert 'href="/"' in c.get("/es").text
    assert 'href="/es/about"' in c.get("/about").text
    assert 'href="/about"' in c.get("/es/about").text


def test_about_es_translates(tmp_path):
    r = _client(tmp_path).get("/es/about")
    assert r.status_code == 200
    assert "Acerca de quiverlab-web" in r.text


def test_literature_page_renders(tmp_path):
    r = _client(tmp_path).get("/literature")
    assert r.status_code == 200
    assert "Literature" in r.text                 # translated chrome (lit.h1)


def test_literature_page_es(tmp_path):
    r = _client(tmp_path).get("/es/literature")
    assert r.status_code == 200
    assert "Bibliografía" in r.text


def test_no_untranslated_key_leaks_in_es_literature(tmp_path):
    from webapp.server.i18n import all_keys
    r = _client(tmp_path).get("/es/literature")
    for key in all_keys():
        assert key not in r.text, f"untranslated key leaked into /es/literature: {key}"


def test_done_job_page_shows_references(tmp_path):
    import json as _json
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid = _done_job(cfg, _json.dumps({
        "reproduce": "import quiverlab as ql", "quiverlab_version": "9.9.9",
        "references": [{"key": "ChouhySolotar2015",
                        "formatted": "Chouhy, Solotar. Projective resolutions... (2015).",
                        "doi_url": "https://doi.org/10.1016/j.jalgebra.2015.02.019",
                        "arxiv_url": "https://arxiv.org/abs/1406.2300"}]}))
    r = TestClient(create_app(cfg)).get("/job/" + jid)
    assert r.status_code == 200
    assert "References" in r.text
    assert "Chouhy, Solotar" in r.text
    assert "https://arxiv.org/abs/1406.2300" in r.text


def test_docs_link_absent_by_default(tmp_path):
    r = _client(tmp_path).get("/")
    assert ">Docs<" not in r.text        # no docs_url configured -> no link


def test_docs_link_present_when_configured(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_DOCS_URL": "https://docs.example"})
    c = TestClient(create_app(cfg))
    assert 'href="https://docs.example"' in c.get("/").text
    assert ">Docs<" in c.get("/").text
    assert ">Documentación<" in c.get("/es").text


def test_big_warning_templates_localized(tmp_path):
    en = _client(tmp_path).get("/").text
    assert "data-big-warn=" in en
    assert "builds ~{cells}-cell" in en            # EN template carries the placeholder
    es = _client(tmp_path).get("/es").text
    assert "matrices de ~{cells} celdas" in es     # ES template, same placeholder


# --------------------------------------------------------------------------- #
# Added hardening: downloads (whitelist + done-gate + traversal), key-shaped
# regex leak scan, and no-inline-script on every page.
# --------------------------------------------------------------------------- #
def test_download_serves_only_whitelisted_names_for_done_job(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid = _done_job(cfg, '{"ok": true}')
    art = cfg.artifacts_dir / jid
    (art / "secret.txt").write_text("token=hunter2")   # present but NOT whitelisted
    c = TestClient(create_app(cfg))

    ok = c.get(f"/download/{jid}/result.json")
    assert ok.status_code == 200
    assert ok.headers["content-type"].startswith("application/json")
    assert 'attachment' in ok.headers.get("content-disposition", "")

    bad = c.get(f"/download/{jid}/secret.txt")         # non-whitelisted -> 404
    assert bad.status_code == 404
    assert "hunter2" not in bad.text


def test_download_rejects_traversal(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid = _done_job(cfg, '{"ok": true}')
    # A secret sitting one level up from the job's artifact dir.
    (cfg.artifacts_dir / "secret.txt").write_text("token=hunter2")
    c = TestClient(create_app(cfg))
    for payload in ["../../etc/passwd", "..%2f..%2fsecret.txt",
                    "%2e%2e%2fsecret.txt", "result.json.bak"]:
        r = c.get(f"/download/{jid}/{payload}")
        assert r.status_code in (400, 404), payload
        assert "root:" not in r.text and "hunter2" not in r.text


def test_download_404_for_non_done_job(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    store = JobStore(cfg.db_path)
    store.init_schema()
    jid = store.create_job({}, ip="h")                 # pending, never marked done
    art = cfg.artifacts_dir / jid
    art.mkdir(parents=True, exist_ok=True)
    (art / "result.json").write_text('{"ok": true}')
    r = TestClient(create_app(cfg)).get(f"/download/{jid}/result.json")
    assert r.status_code == 404


def test_download_trace_html_fallback_served(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid = _done_job(cfg, '{"ok": true}')
    (cfg.artifacts_dir / jid / "trace_steps.html").write_text("<p>steps</p>")
    r = TestClient(create_app(cfg)).get(f"/download/{jid}/trace_steps.html")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/html")


_KEY_SHAPE = re.compile(r"[a-z][a-z_]*\.[a-z_]+")
_TAG = re.compile(r"<[^>]+>")
_SCRIPT_OPEN = re.compile(r"<script\b[^>]*>", re.IGNORECASE)


def _visible_text(html: str) -> str:
    return _TAG.sub(" ", html)


def test_no_key_shaped_token_leaks_in_es_visible_text(tmp_path):
    """Stricter than the exact-key scan: any catalog-key-shaped token appearing
    in the /es *visible* text (tags stripped) that is an actual catalog key is a
    leak. Filenames/paths live in attributes and are stripped, so a real family
    name or asset path never trips this."""
    from webapp.server.i18n import all_keys
    keys = all_keys()
    for page in ("/es", "/es/about", "/es/literature"):
        visible = _visible_text(_client(tmp_path).get(page).text)
        leaked = [m.group(0) for m in _KEY_SHAPE.finditer(visible) if m.group(0) in keys]
        assert not leaked, f"{page}: leaked keys {leaked}"


def test_no_inline_script_on_any_page(tmp_path):
    """Strict CSP is script-src 'self': every <script> MUST carry a src= (no
    inline JS anywhere)."""
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    jid = _done_job(cfg, '{"reproduce": "import quiverlab as ql", '
                         '"quiverlab_version": "9.9.9"}')
    c = TestClient(create_app(cfg))
    for page in ("/", "/es", "/about", "/es/about", "/literature",
                 "/es/literature", f"/job/{jid}", f"/es/job/{jid}"):
        html = c.get(page).text
        for tag in _SCRIPT_OPEN.findall(html):
            assert "src=" in tag, f"{page}: inline <script> found: {tag}"


_APP_JS = Path(__file__).resolve().parents[2] / "webapp" / "static" / "app.js"

# The whole HTML-injection class -- not just innerHTML. Any of these can parse an
# attacker-controlled string as markup; app.js must use none of them.
_HTML_INJECTION_SINKS = (
    "innerHTML", "outerHTML", "insertAdjacentHTML",
    "document.write", "setHTMLUnsafe", "createContextualFragment",
)


def _strip_js_comments_and_strings(src: str) -> str:
    """Best-effort removal of ``//`` line comments and string literals
    (``"..."`` / ``'...'`` / ``` `...` ```) from JS source so the sink scan sees
    only executable code -- honest comments and URL literals may now name the
    forbidden tokens without tripping the guard.

    This is a GUARD, not a JS parser. Known limits (all absent from app.js, which
    is why the approximation is adequate): it does not handle ``/* block
    comments */``, regex literals, or template-literal ``${...}`` interpolation.
    Strings are stripped BEFORE comments so a ``//`` inside a string (e.g.
    ``"https://"``) is never mistaken for a comment; string patterns are
    newline-bounded so an unbalanced quote inside a ``//`` comment cannot run
    away across lines."""
    src = re.sub(r'"(?:[^"\\\n]|\\.)*"', "", src)   # double-quoted
    src = re.sub(r"'(?:[^'\\\n]|\\.)*'", "", src)   # single-quoted
    src = re.sub(r"`(?:[^`\\]|\\.)*`", "", src)     # template literal (best-effort)
    src = re.sub(r"//[^\n]*", "", src)              # // line comments
    return src


def test_app_js_has_no_html_injection_sinks():
    """The strict CSP blocks script *execution* but NOT HTML/link injection, so a
    regression to any HTML-parsing sink would be a live XSS/link-injection hole.
    app.js is a served static file (not exercised by TestClient's non-JS render),
    so this static scan is the honest, checkable fence. Comments and string
    literals are stripped first (see ``_strip_js_comments_and_strings``) so the
    assertion fires only on the tokens in *executable* code."""
    assert _APP_JS.is_file(), _APP_JS
    code = _strip_js_comments_and_strings(_APP_JS.read_text())
    for token in _HTML_INJECTION_SINKS:
        assert token not in code, (
            f"app.js contains HTML-injection sink '{token}' in executable code; "
            "render server data via textContent / DOM nodes instead"
        )


def test_is_http_url_rejects_non_http_schemes():
    """Table-test app.js's ``isHttpUrl`` gate under node: extract the function
    source textually and evaluate it against a case table in one ``node -e`` run.
    Skipped (never failed) when node is unavailable -- a missing runtime is not a
    regression."""
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not available")

    src = _APP_JS.read_text()
    m = re.search(r"function isHttpUrl\s*\([^)]*\)\s*\{.*?\n\}", src, re.DOTALL)
    assert m, "isHttpUrl function not found in app.js"
    fn_src = m.group(0)

    rejected = ["//evil.com", " https://x", "HTTP://x", "https:@evil",
                "javascript:alert(1)", "data:text/html,x", "", None]
    accepted = ["https://doi.org/10.1/x", "http://example.org"]

    harness = fn_src + "\n" + (
        f"const rejected = {json.dumps(rejected)};\n"
        f"const accepted = {json.dumps(accepted)};\n"
        "process.stdout.write(JSON.stringify({\n"
        "  rejected: rejected.map((u) => isHttpUrl(u)),\n"
        "  accepted: accepted.map((u) => isHttpUrl(u)),\n"
        "}));\n"
    )
    proc = subprocess.run([node, "-e", harness],
                          capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stderr
    got = json.loads(proc.stdout)
    assert got["rejected"] == [False] * len(rejected), \
        f"isHttpUrl accepted a bad href: {list(zip(rejected, got['rejected']))}"
    assert got["accepted"] == [True] * len(accepted), \
        f"isHttpUrl rejected a good href: {list(zip(accepted, got['accepted']))}"
