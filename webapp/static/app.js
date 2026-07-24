// External static script (no inline JS — the CSP forbids it).
//
// The strict CSP (script-src 'self') blocks script *execution*, but it does NOT
// block HTML/link injection. So every server-provided string below is placed via
// textContent / DOM node APIs (never innerHTML / insertAdjacentHTML / document.write)
// and reference links are built with document.createElement("a"), accepting
// http(s) hrefs only.

// Interpolate {cells}/{minutes}/{maxcells} into a t()-sourced template string.
function interp(tmpl, est, maxcells) {
  return (tmpl || "")
    .replace("{cells}", (est && est.cells != null) ? est.cells.toLocaleString() : "?")
    .replace("{minutes}", (est && est.minutes != null) ? est.minutes : "?")
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

function readComputeBody() {
  const fd = new FormData(form);
  const fieldRaw = fd.get("field").trim();
  const field = fieldRaw === "CC"
    ? {kind: "CC"}
    : {kind: "GF", p: parseInt(fieldRaw.replace(/[^0-9]/g, ""), 10), n: 1};
  return {
    schema: 1,
    algebra: {kind: "family", family: fd.get("family"),
              params: JSON.parse(fd.get("params") || "{}"), field},
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
