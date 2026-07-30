/* quiverlab draw-page Web Worker -- the SERVER-BACKED twin of docs/gui/worker.js.
   Same message protocol as the Pyodide worker (gui.js is vendored byte-identical
   and cannot tell them apart); compute runs on this webapp's own guarded API
   instead of in-browser Pyodide:

     init      -> ready            (no engine download: the server IS the engine)
     calibrate -> calibrated       (nominal factor; wall times are server-side)
     probe     -> probe            (POST /api/gui/probe: build-only, dim + label)
     run       -> built, result*, [trace], done   (POST /api/compute + job poll)

   Absolute API paths: a worker resolves relative URLs against ITS OWN script URL
   (/gui/worker.js), which would mangle them. The Spanish mount serves this same
   file at /es/gui/worker.js; the API is unprefixed either way. */
"use strict";
var VERSION = "";

self.onmessage = function (e) {
  var m = e.data;
  var job = m.cmd === "init" ? init(m.manifest)
          : m.cmd === "run" ? run(m.request)
          : m.cmd === "probe" ? probe(m)
          : m.cmd === "calibrate" ? calibrate()
          : Promise.reject(new Error("unknown cmd " + m.cmd));
  job.catch(function (err) {
    self.postMessage({ type: "fatal", message: String(err && err.message || err) });
  });
};

async function init(manifest) {
  VERSION = (manifest && manifest.quiverlab_version) || "";
  self.postMessage({ type: "ready", version: VERSION });
}

async function calibrate() {
  // gui.js probes only once factor !== null. Wall time is measured server-side,
  // so a nominal factor just switches the live dim/ETA hints on.
  self.postMessage({ type: "calibrated", factor: 1, seconds: 0 });
}

async function postJSON(url, body) {
  var r = await fetch(url, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  var data = null;
  try { data = await r.json(); } catch (err) { /* non-JSON body */ }
  return { status: r.status, data: data };
}

function protocolError(status, data, fallback) {
  if (data && data.error) return { ok: false, error: data.error };
  if (data && data.error_type) {
    return { ok: false, error: { type: data.error_type, message: data.message || fallback } };
  }
  // FastAPI validation failure (422): surface the actual field errors instead
  // of a bare status code -- the user typed something the schema refuses.
  if (data && data.detail) {
    var d = data.detail;
    var msg = typeof d === "string" ? d
      : (Array.isArray(d) ? d.map(function (e) {
          return ((e.loc || []).join(".") + ": " + (e.msg || "")).trim();
        }).join("; ") : JSON.stringify(d));
    return { ok: false, error: { type: "ValidationError", message: msg } };
  }
  return { ok: false, error: { type: "HTTP " + status, message: fallback } };
}

async function probe(m) {
  var r = await postJSON("/api/gui/probe", m.request);
  var data = (r.data && typeof r.data.ok === "boolean")
    ? r.data
    : protocolError(r.status, r.data, "probe failed");
  self.postMessage({ type: "probe", seq: m.seq, data: data });
}

async function fetchArtifact(jobId, name) {
  try {
    var r = await fetch("/download/" + jobId + "/" + name);
    if (!r.ok) return "";
    return await r.text();
  } catch (err) {
    return "";
  }
}

async function pollJob(jobId) {
  // Poll until the job reaches a TERMINAL state. There is deliberately no client
  // deadline: the offline app has no wall cap, so a real computation may run for
  // hours -- the page giving up at 30 minutes on a job that is still going was
  // just a lie about a live job (Marco 2026-07-30). The deployed server does cap
  // its jobs, so there the poll still ends promptly, via the failed status.
  //
  // The loop only ever exits on: done, failed/error, or the job vanishing (a
  // purge or a server restart) -- the last after a few consecutive misses, so a
  // transient blip does not abandon a running job.
  var waited = 0, misses = 0;
  for (;;) {
    var delay = waited < 20000 ? 1000 : (waited < 600000 ? 5000 : 15000);
    await new Promise(function (res) { setTimeout(res, delay); });
    waited += delay;
    var r;
    try {
      r = await fetch("/api/jobs/" + jobId);
    } catch (err) {                       // network blip: keep waiting
      if (++misses >= 20) return "lost contact with the local server while the " +
        "job was running (its permalink keeps working: /job/" + jobId + ")";
      continue;
    }
    if (!r.ok) {
      if (++misses >= 20) return "job " + jobId + ": HTTP " + r.status;
      continue;
    }
    misses = 0;
    var job = await r.json();
    if (job.status === "done") return null;
    if (job.status === "failed" || job.status === "error") {
      return job.error || "job failed";
    }
  }
}

async function run(request) {
  // 1. "built": the same bounded build the probe endpoint does (label + dim).
  var b = await postJSON("/api/gui/probe", request);
  var built = (b.data && typeof b.data.ok === "boolean")
    ? b.data
    : protocolError(b.status, b.data, "build failed");
  self.postMessage({ type: "built", data: built, eta: null });
  if (!built.ok) { self.postMessage({ type: "done" }); return; }

  // 2. compute: instant answers inline; queued/cached resolve via the job's
  //    result.json artifact (the job status payload carries no result body).
  var t0 = Date.now();
  var r = await postJSON("/api/compute", request);
  var body = r.data || {};
  var result = null, jobId = null, failure = null;

  if (body.tier === "instant" && body.result) {
    result = body.result;
  } else if (body.tier === "queued" || body.tier === "cached") {
    jobId = body.job_id;
    if (body.tier === "queued") {
      self.postMessage({ type: "probe", seq: -1, data: built }); // keep the hint fresh
      failure = await pollJob(jobId);
    }
    if (!failure) {
      var raw = await fetchArtifact(jobId, "result.json");
      if (raw) result = JSON.parse(raw);
      else failure = "finished job exposed no result.json";
    }
  } else if (body.tier === "big") {
    failure = "this example is beyond the local instant/queued tiers" +
      (body.estimate ? " (estimate: " + JSON.stringify(body.estimate) + ")" : "");
  } else {
    var pe = protocolError(r.status, body, "compute failed");
    failure = pe.error.type + ": " + pe.error.message;
  }

  if (failure || !result || !result.results) {
    self.postMessage({ type: "result",
      data: { ok: false, error: { type: "ComputeError",
                                  message: failure || "no results returned" } } });
    self.postMessage({ type: "done" });
    return;
  }

  // 3. one "result" message per requested invariant, in request order -- the
  //    exact shape the Pyodide runner's compute_one returns: {ok, invariant, block}.
  var elapsed = Date.now() - t0;
  for (var i = 0; i < request.compute.length; i++) {
    var inv = request.compute[i];
    var kind = inv.split(":")[0];
    var block = result.results[kind];
    if (block === undefined) {
      self.postMessage({ type: "result",
        data: { ok: false, error: { type: "MissingBlock",
                                    message: "no result block for " + inv } } });
      continue;
    }
    self.postMessage({ type: "result",
      data: { ok: true, invariant: inv, block: block },
      elapsed_ms: i === 0 ? elapsed : 0, eta: null });
  }

  // 4. worked-steps trace, when the job persisted it (queued/cached jobs only:
  //    the instant tier returns no artifact directory).
  var tikz = "";
  if (jobId) {
    var html = await fetchArtifact(jobId, "trace_steps.html");
    var json = await fetchArtifact(jobId, "trace.json");
    if (html || json) self.postMessage({ type: "trace", html: html, json: json });
    tikz = await fetchArtifact(jobId, "tikz.tex");
  }
  // bundle = the full result document; snippet = runnable library code for the
  // same request (the Pyodide runner generates these in-engine; here the result
  // is already at hand and the snippet is derivable from the request alone).
  self.postMessage({ type: "artifacts",
    tikz: tikz,
    snippet: pythonSnippet(request),
    bundle: JSON.stringify(result, null, 1) });

  self.postMessage({ type: "done" });
}

function pythonSnippet(request) {
  var alg = request.algebra || {};
  var f = alg.field || {};
  var field = f.kind === "GF"
    ? "ql.GF(" + f.p + (f.n && f.n > 1 ? ", " + f.n : "") + ")"
    : "ql.CC";
  var lines = ["import quiverlab as ql", ""];
  if (alg.kind === "family") {
    lines.push("A = ql." + alg.family + "(" +
      Object.keys(alg.params || {}).map(function (k) {
        return k + "=" + JSON.stringify(alg.params[k]);
      }).join(", ") + (Object.keys(alg.params || {}).length ? ", " : "") +
      "field=" + field + ")");
  } else {
    lines.push("Q = ql.Quiver(vertices=" + JSON.stringify(alg.vertices) + ",");
    lines.push("              arrows=" + JSON.stringify(alg.arrows).replace(/\[/g, "(").replace(/\]/g, ")") + ")");
    lines.push("A = Q.algebra(relations=" + JSON.stringify(alg.relations || []) +
               ", field=" + field + ")");
  }
  (request.compute || []).forEach(function (inv) {
    var parts = inv.split(":"), kind = parts[0];
    var top = parts[1] ? parts[1].split("..").pop() : null;
    var call = { hh_cohomology: "hochschild_cohomology", hh_homology: "hochschild_homology" }[kind];
    if (call) lines.push("print(A." + call + "(" + top + "))");
    else lines.push("print(A." + kind + "(" + (top !== null ? top : "") + "))  # see docs for the exact call");
  });
  return lines.join("\n") + "\n";
}
