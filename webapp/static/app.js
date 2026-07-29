// External static script (no inline JS — the CSP forbids it).
//
// The strict CSP (script-src 'self') blocks script *execution*, but it does NOT
// block HTML/link injection. So every server-provided string below is placed via
// textContent / DOM node APIs (never innerHTML / insertAdjacentHTML / document.write)
// and reference links are built with document.createElement("a"), accepting
// http(s) hrefs only.

// Interpolate {cells}/{minutes}/{mem}/{maxcells} into a t()-sourced template
// string. {mem} is the server's human-readable memory estimate (estimate.mem_human).
function interp(tmpl, est, maxcells) {
  return (tmpl || "")
    .replace("{cells}", (est && est.cells != null) ? est.cells.toLocaleString() : "?")
    .replace("{minutes}", (est && est.minutes != null) ? est.minutes : "?")
    .replace("{mem}", (est && est.mem_human != null) ? est.mem_human : "?")
    .replace("{maxcells}", maxcells ? Number(maxcells).toLocaleString() : "?");
}

// A styled error box carrying plain text (no markup interpretation).
function errDiv(text) {
  const div = document.createElement("div");
  div.className = "err";
  div.textContent = text;
  return div;
}

// Honesty note: this gate intentionally drops every non-http(s):// href — bare
// DOIs (10.1/x), protocol-relative (//host), uppercase schemes (HTTP://) — and the
// server currently emits only https://… or null (see webapp/server/references.py).
// True only for http:// or https:// hrefs — everything else (javascript:, data:,
// relative, …) is rejected so we never build an attacker-controlled link.
function isHttpUrl(url) {
  return typeof url === "string"
    && (url.startsWith("https://") || url.startsWith("http://"));
}

// Append " · <a href=url>label</a>" to li, but only for http(s) URLs.
function appendLink(li, url, label) {
  if (!isHttpUrl(url)) return;
  li.append(" · ");
  const a = document.createElement("a");
  a.href = url;               // property assignment, not markup parsing
  a.textContent = label;
  li.appendChild(a);
}

// Format a dimension-vector object like {"1": 2, "2": 0} as "{1: 2, 2: 0}".
function dvText(dv) {
  return "{" + Object.keys(dv || {}).map(function (k) {
    return k + ": " + dv[k];
  }).join(", ") + "}";
}

// Build a small <table> with a header row of textContent strings.
function tableWithHeader(headers) {
  const tbl = document.createElement("table");
  const hr = document.createElement("tr");
  for (const label of headers) {
    const th = document.createElement("th");
    th.textContent = label;
    hr.appendChild(th);
  }
  tbl.appendChild(hr);
  return tbl;
}

function cellText(text) {
  const td = document.createElement("td");
  td.textContent = text;
  return td;
}

// Plan 30 (Marco #3): a resolution's rendered table shows term | ⊕-decomposition.
// The "# summands" and dim-vector columns are gone from the RENDERING; the raw
// betti/terms fields stay in the JSON dump above. The decomposition cell carries
// the block's LaTeX (P_1^2 ⊕ P_3), typeset later by renderMath.
function resolutionTable(block, d) {
  const tbl = tableWithHeader([d.colTerm || "term",
                               d.colDecomp || "⊕-decomposition"]);
  const summands = block.summands || [];
  const terms = block.terms || [];
  for (let n = 0; n < terms.length; n++) {
    const tr = document.createElement("tr");
    tr.appendChild(cellText(String(n)));
    const td = document.createElement("td");
    const span = document.createElement("span");
    span.className = "arithmatex";
    span.textContent = "\\(" + (summands[n] != null ? summands[n] : "0") + "\\)";
    td.appendChild(span);
    tr.appendChild(td);
    tbl.appendChild(tr);
  }
  return tbl;
}

// Krull–Schmidt summand table: summand | multiplicity | dim vector (Plan 30).
function decomposeBlock(block, d) {
  const wrap = document.createElement("div");
  const p = document.createElement("p");
  p.textContent = (d.modDecompHeading || "Krull–Schmidt decomposition")
    + " (" + block.iso_classes + ")";
  wrap.appendChild(p);
  const tbl = tableWithHeader([d.modSummand || "summand",
                               d.modMult || "multiplicity",
                               d.modDimvec || "dim vector"]);
  (block.summands || []).forEach(function (s, i) {
    const tr = document.createElement("tr");
    tr.appendChild(cellText("M_" + (i + 1)));
    tr.appendChild(cellText(String(s.multiplicity)));
    tr.appendChild(cellText(dvText(s.dim_vector)));
    tbl.appendChild(tr);
  });
  wrap.appendChild(tbl);
  return wrap;
}

// Plan 34 (Marco): rad/top/soc as FULL representations -- the dim VECTOR per object
// (the redundant total-dim column is gone) plus each arrow's exact action matrix,
// typeset by renderMath (KaTeX). The LaTeX is built from the block's {dims, maps}.

// Matrices are the COMPLETE human record (Plan 34, Marco): shown IN FULL, wrapped in
// a horizontally-scrollable container (mathScroll) so the page body never scrolls
// sideways -- NOT elided at a small size. MAT_BACKSTOP_CELLS is only a SANITY cap
// mirroring the trace recorder's record-time memory backstop: a pathological/corrupt
// payload past it is stated by shape instead of hanging the browser. One constant per
// file (gui.js has the same one, same comment).
const MAT_BACKSTOP_CELLS = 250000;
function matTooBig(mat) {
  const rows = (mat || []).length;
  const cols = rows ? (mat[0] || []).length : 0;
  return rows * cols > MAT_BACKSTOP_CELLS;
}
function matLatex(mat) {                    // [[..],[..]] -> \begin{pmatrix}..\end{pmatrix}
  mat = mat || [];
  if (matTooBig(mat)) {                     // sanity backstop only (never normal use)
    const cols = mat.length ? (mat[0] || []).length : 0;
    return "\\text{[" + mat.length + "\\times" + cols + " matrix beyond the display backstop]}";
  }
  const body = mat.map(function (row) {
    return row.map(String).join(" & ");
  }).join(" \\\\ ");
  return "\\begin{pmatrix} " + body + " \\end{pmatrix}";
}

// A matrix wrapped in a horizontally-scrollable inline box (Plan 34, Marco): the
// full matrix scrolls INSIDE this box, so a wide differential never makes the page
// body scroll sideways. Inline-block keeps it inside the surrounding <p>.
function mathScroll(latex) {
  const box = document.createElement("span");
  box.style.display = "inline-block";
  box.style.maxWidth = "100%";
  box.style.overflowX = "auto";
  box.style.verticalAlign = "middle";
  const span = document.createElement("span");
  span.className = "arithmatex";
  span.textContent = "\\(" + latex + "\\)";
  box.appendChild(span);
  return box;
}

// A pre-Plan-34 cached rad/top/soc lacked the per-view {dims, maps}; guard so an old
// shape is called out honestly (MINOR-6) rather than silently rendering a fabricated
// "every arrow acts as zero".
function radTopSocStale(block) {
  return [block.radical, block.top, block.socle].some(function (v) {
    return !v || v.dims == null || v.maps == null;
  });
}
function radTopSocDisplayOnly(block) {      // any view carries non-re-enterable entries
  return [block.radical, block.top, block.socle].some(function (v) {
    return v && v.display_only === true;
  });
}

function radTopSocBlock(block, d) {
  const wrap = document.createElement("div");
  const p = document.createElement("p");
  p.textContent = d.modRadTopSoc || "radical / top / socle";
  wrap.appendChild(p);
  if (radTopSocStale(block)) {              // MINOR-6: honest "recompute", never a fake zero
    wrap.appendChild(errDiv(d.modStale ||
      "this result was computed by an older version — recompute to see the full representation."));
    return wrap;
  }
  if (radTopSocDisplayOnly(block)) {        // MAJOR-4: extension-field entries are not re-enterable
    const note = document.createElement("p");
    note.className = "pdf-note";
    note.textContent = d.modDisplayOnly || "display only — entries are not re-enterable in the module panel.";
    wrap.appendChild(note);
  }
  const trio = [["rad M", block.radical], ["top M", block.top], ["soc M", block.socle]];
  const tbl = tableWithHeader(["", d.modDimvec || "dim vector"]);
  for (const pair of trio) {
    const tr = document.createElement("tr");
    const th = document.createElement("th");
    th.textContent = pair[0];
    tr.appendChild(th);
    tr.appendChild(cellText(dvText((pair[1] || {}).dims)));
    tbl.appendChild(tr);
  }
  wrap.appendChild(tbl);
  for (const pair of trio) {
    const label = pair[0], maps = (pair[1] || {}).maps || {};
    const arrows = Object.keys(maps);
    if (!arrows.length) {
      const q = document.createElement("p");
      q.textContent = label + ": " + (d.modArrowsZero || "every arrow acts as zero");
      wrap.appendChild(q);
      continue;
    }
    for (const a of arrows) {
      const q = document.createElement("p");
      q.appendChild(document.createTextNode(label + ", arrow " + a + ": "));
      q.appendChild(mathScroll(matLatex(maps[a])));   // full matrix, scrollable box
      wrap.appendChild(q);
    }
  }
  return wrap;
}

// The AR-translate input certificate (Marco #1): indecomposable, or the input's
// decomposition + the additivity note. null when the block carries no certificate
// (the decompose engine was unavailable), so the note never lies.
function tauCertNote(block, d) {
  const p = document.createElement("p");
  if (block.indecomposable === true) {
    p.textContent = d.modIndecomposable || "input M is indecomposable";
    return p;
  }
  if (block.decomposition) {
    const parts = block.decomposition.map(function (s) {
      return dvText(s.dim_vector) + (s.multiplicity > 1 ? "^" + s.multiplicity : "");
    }).join("  ⊕  ");
    p.textContent = "M ≅ " + parts + " — "
      + (d.modTauAdditive || "τ computed summand-wise (τ is additive)");
    return p;
  }
  return null;
}

// Render structured tables for any module blocks in the result (reachable via the
// /api/compute API). Additive: it complements the raw JSON dump, never replaces it.
function renderModuleBlocks(out, res) {
  const results = res.results || {};
  const d = form ? form.dataset : {};
  for (const kind of Object.keys(results)) {
    const b = results[kind];
    if (!b || typeof b !== "object") continue;
    if (kind === "projective_resolution" || kind === "injective_resolution") {
      out.appendChild(resolutionTable(b, d));
    } else if (kind === "rad_top_soc") {
      out.appendChild(radTopSocBlock(b, d));
    } else if (kind === "decompose") {
      out.appendChild(decomposeBlock(b, d));
    } else if (kind === "tau" || kind === "tau_minus") {
      const note = tauCertNote(b, d);
      if (note) out.appendChild(note);
    }
  }
}

// Render a compute result into `out` using DOM nodes (server strings via
// textContent). Keeps the same visual structure as before (h2/pre/p/ul/li/h3).
function renderResult(out, res) {
  out.textContent = "";       // clear prior render

  const h2 = document.createElement("h2");
  h2.textContent = "Result";
  out.appendChild(h2);

  const pre = document.createElement("pre");
  pre.textContent = JSON.stringify(res.results, null, 2);
  out.appendChild(pre);

  renderModuleBlocks(out, res);   // Plan 30: structured module tables

  if (res.quiverlab_version) {
    const p = document.createElement("p");
    p.append("quiverlab ");
    const code = document.createElement("code");
    code.textContent = res.quiverlab_version;
    p.appendChild(code);
    out.appendChild(p);
  }

  if (res.references && res.references.length) {
    const h3 = document.createElement("h3");
    h3.textContent = form.dataset.refsLabel || "References";
    out.appendChild(h3);
    const ul = document.createElement("ul");
    for (const ref of res.references) {
      const li = document.createElement("li");
      li.textContent = ref.formatted;
      appendLink(li, ref.doi_url, "DOI");
      appendLink(li, ref.arxiv_url, "arXiv");
      ul.appendChild(li);
    }
    out.appendChild(ul);
  }

  if (res.reproduce) {
    const h3 = document.createElement("h3");
    h3.textContent = "Reproduce locally";
    out.appendChild(h3);
    const pre2 = document.createElement("pre");
    pre2.textContent = res.reproduce;
    out.appendChild(pre2);
  }
}

// Index page: POST the form to /api/compute, render results or redirect to a job.
const form = document.getElementById("compute-form");

// ---- per-parameter inputs (Marco 2026-07-28: explain every parameter, prefill
// defaults, never demand JSON from scratch). Data source: GET /api/catalog --
// per-family params carry kind/default plus curated example + bilingual help.
let CATALOG = null;

function paramPrefill(p) {
  const v = p.default !== null && p.default !== undefined ? p.default
          : (p.example !== null && p.example !== undefined ? p.example : "");
  if (v === "") return "";
  return typeof v === "string" ? v : JSON.stringify(v);
}

function renderParamFields() {
  const box = document.getElementById("param-fields");
  const summary = document.getElementById("family-summary");
  if (!box || !CATALOG) return;
  const lang = (form.dataset.lang === "es") ? "es" : "en";
  const fam = CATALOG.families.find((f) => f.name === form.elements.family.value);
  box.textContent = "";
  if (summary) summary.textContent = (fam && fam.summary) ? (fam.summary[lang] || "") : "";
  if (!fam) return;
  for (const p of fam.params) {
    const row = document.createElement("div");
    row.className = "param-row";
    const label = document.createElement("label");
    label.textContent = p.name;
    label.setAttribute("for", "param-" + p.name);
    row.appendChild(label);
    const input = document.createElement("input");
    input.id = "param-" + p.name;
    input.dataset.param = p.name;
    if (p.kind === "bool") {
      input.type = "checkbox";
      input.checked = p.default === true;
    } else {
      input.type = "text";
      input.value = paramPrefill(p);
    }
    row.appendChild(input);
    if (p.help && p.help[lang]) {
      const help = document.createElement("div");
      help.className = "param-help";
      help.textContent = p.help[lang];
      row.appendChild(help);
    }
    box.appendChild(row);
  }
}

function collectParams() {
  const out = {};
  for (const input of document.querySelectorAll("#param-fields [data-param]")) {
    const name = input.dataset.param;
    if (input.type === "checkbox") {
      if (input.checked) out[name] = true;
      continue;
    }
    const raw = input.value.trim();
    if (!raw) continue;                    // empty -> omit (builder default wins)
    try {
      out[name] = JSON.parse(raw);         // numbers, lists, true/false
    } catch (err) {
      out[name] = raw;                     // exact strings: "A3", "1/2", ...
    }
  }
  return out;
}

async function initParamForm() {
  if (!form || !document.getElementById("param-fields")) return;
  try {
    const r = await fetch("/api/catalog");
    CATALOG = await r.json();
  } catch (err) {
    return;                                // form still submits with {} params
  }
  renderParamFields();
  form.elements.family.addEventListener("change", renderParamFields);
}
initParamForm();

function readComputeBody() {
  const fd = new FormData(form);
  const fieldRaw = fd.get("field").trim();
  const field = fieldRaw === "CC"
    ? {kind: "CC"}
    : {kind: "GF", p: parseInt(fieldRaw.replace(/[^0-9]/g, ""), 10), n: 1};
  return {
    schema: 1,
    algebra: {kind: "family", family: fd.get("family"),
              params: collectParams(), field},
    compute: fd.getAll("compute"),
    artifacts: {pdf: fd.get("pdf") === "1", tikz: false},
  };
}

if (form) form.addEventListener("submit", async (e) => {
  e.preventDefault();
  const body = readComputeBody();
  const r = await fetch("/api/compute", {method: "POST",
    headers: {"content-type": "application/json"}, body: JSON.stringify(body)});
  const out = document.getElementById("out");
  const data = await r.json();
  // Over the anonymous cap but runnable as a big job: reveal the email field.
  if (r.status === 202 && data.tier === "big") {
    document.getElementById("big-warn").textContent = interp(form.dataset.bigWarn, data.estimate);
    document.getElementById("big-warn").style.display = "block";
    document.getElementById("big-email-row").style.display = "block";
    return;
  }
  if (r.status === 202 && data.tier === "queued") { window.location = "/job/" + data.job_id; return; }
  // Cache hit (Plan 25): this exact computation was done before -- go straight to
  // its result page (downloads + references + the "previously computed" note).
  if (r.status === 200 && data.tier === "cached") { window.location = "/job/" + data.job_id; return; }
  // Beyond big caps, or big tier disabled: honest message + local hint.
  if (r.status === 422 && data.reason) {
    const tmpl = data.reason === "beyond_big_cap" ? form.dataset.bigReject : form.dataset.bigDisabled;
    out.replaceChildren(errDiv(interp(tmpl, data.estimate, form.dataset.bigCap)));
    return;
  }
  if (!r.ok) { out.replaceChildren(errDiv(data.error_type + ": " + data.message)); return; }
  renderResult(out, data.result);
  renderMath(out);
});

// Big-job email submit: send the spec + email to /api/jobs/big for a magic link.
const bigSend = document.getElementById("big-send");
if (bigSend && form) bigSend.addEventListener("click", async () => {
  const email = document.getElementById("big-email").value.trim();
  if (!email) return;
  const body = readComputeBody();
  body.email = email;
  body.lang = form.dataset.lang || "en";
  const r = await fetch("/api/jobs/big", {method: "POST",
    headers: {"content-type": "application/json"}, body: JSON.stringify(body)});
  const sent = document.getElementById("big-sent");
  // Cache hit (Plan 25): already computed by someone -- served with NO email, NO
  // token. Skip the inbox message and go straight to the cached result page.
  if (r.status === 200) { const d = await r.json();
    if (d.status === "cached") { window.location = "/job/" + d.job_id; return; } }
  if (r.status === 202) { sent.textContent = form.dataset.bigSent || "Check your inbox."; }
  else { const d = await r.json(); sent.replaceChildren(errDiv(d.message || "error")); }
});

// Job page: while the job runs, poll its status every 2s; reload on completion
// so the server-rendered downloads + reproduce snippet appear.
const meta = document.getElementById("job-meta");
if (meta && ["pending", "running"].includes(meta.dataset.status)) {
  const jobId = meta.dataset.jobId;
  const timer = setInterval(async () => {
    const r = await fetch("/api/jobs/" + jobId);
    if (!r.ok) return;
    const s = await r.json();
    const statusEl = document.getElementById("status");
    if (statusEl) statusEl.textContent = s.status;
    const prog = document.getElementById("progress");
    if (prog && s.progress) prog.textContent = JSON.stringify(s.progress);
    if (s.status === "done" || s.status === "failed") {
      clearInterval(timer);
      window.location.reload();
    }
  }, 2000);
}

// Feedback page: POST the form to /api/feedback, show the reference id.
// CSP-safe: every server string goes through textContent / errDiv, never
// innerHTML. Reuses the errDiv() helper defined above.
const fbForm = document.getElementById("feedback-form");
if (fbForm) {
  // Reveal the literature-only fields when that category is selected.
  const catSel = document.getElementById("category");
  const litFields = document.getElementById("lit-fields");
  const syncLit = () => {
    if (litFields) {
      litFields.style.display =
        (catSel && catSel.value === "literature") ? "block" : "none";
    }
  };
  if (catSel) catSel.addEventListener("change", syncLit);
  syncLit();

  fbForm.addEventListener("submit", async (e) => {
    e.preventDefault();
    const fd = new FormData(fbForm);
    const body = {
      category: fd.get("category"),
      message: fd.get("message"),
      contact: fd.get("contact") || null,
      job_ref: fd.get("job_ref") || null,
      reference: fd.get("reference") || null,
      why_relevant: fd.get("why_relevant") || null,
      website: fd.get("website") || "",
    };
    const r = await fetch("/api/feedback", {method: "POST",
      headers: {"content-type": "application/json"}, body: JSON.stringify(body)});
    const out = document.getElementById("fb-out");
    const data = await r.json();
    if (r.status === 201) {
      out.textContent =
        (fbForm.dataset.thanks || "Thank you. Your reference is") + " " + data.reference;
      fbForm.reset();
      syncLit();
    } else {
      out.replaceChildren(errDiv(data.message || "error"));
    }
  });
}

// KaTeX (vendored, no CDN): typeset $…$ / \(…\) / $$…$$ math. Guarded so the page
// still works if the optional contrib script is unavailable. app.js is the last
// script in <body>, so the DOM is fully parsed here.
function renderMath(root) {
  if (typeof renderMathInElement !== "function") return;
  renderMathInElement(root || document.body, {
    delimiters: [
      {left: "$$", right: "$$", display: true},
      {left: "\\(", right: "\\)", display: false},
      {left: "$", right: "$", display: false},
    ],
    throwOnError: false,
  });
}
renderMath(document.body);
