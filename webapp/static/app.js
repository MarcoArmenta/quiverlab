// External static script (no inline JS — the CSP forbids it).

// Interpolate {cells}/{minutes}/{maxcells} into a t()-sourced template string.
function interp(tmpl, est, maxcells) {
  return (tmpl || "")
    .replace("{cells}", (est && est.cells != null) ? est.cells.toLocaleString() : "?")
    .replace("{minutes}", (est && est.minutes != null) ? est.minutes : "?")
    .replace("{maxcells}", maxcells ? Number(maxcells).toLocaleString() : "?");
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
    out.innerHTML = '<div class="err">' + interp(tmpl, data.estimate, form.dataset.bigCap) + "</div>";
    return;
  }
  if (!r.ok) { out.innerHTML = '<div class="err">' + data.error_type + ": " + data.message + "</div>"; return; }
  const res = data.result;
  let html = "<h2>Result</h2><pre>" + JSON.stringify(res.results, null, 2) + "</pre>";
  if (res.quiverlab_version) html += "<p>quiverlab <code>" + res.quiverlab_version + "</code></p>";
  if (res.references && res.references.length) {
    html += "<h3>" + (form.dataset.refsLabel || "References") + "</h3><ul>";
    for (const r of res.references) {
      html += "<li>" + r.formatted;
      if (r.doi_url) html += ' · <a href="' + r.doi_url + '">DOI</a>';
      if (r.arxiv_url) html += ' · <a href="' + r.arxiv_url + '">arXiv</a>';
      html += "</li>";
    }
    html += "</ul>";
  }
  if (res.reproduce) html += "<h3>Reproduce locally</h3><pre>" + res.reproduce + "</pre>";
  out.innerHTML = html;
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
  else { const d = await r.json(); sent.innerHTML = '<div class="err">' + (d.message || "error") + "</div>"; }
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
