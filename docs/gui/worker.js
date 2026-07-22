/* quiverlab GUI Web Worker (Plans 10+11): Pyodide + the site's own wheel.
   PINNED Pyodide version — never a floating tag (supply-chain posture of
   mkdocs.yml). Protocol: see gui.js. All runner I/O is JSON strings. */
"use strict";
var PYODIDE_BASE = "https://cdn.jsdelivr.net/pyodide/v0.28.3/full/";
var pyodide = null, runner = null;

self.onmessage = function (e) {
  var m = e.data;
  var job = m.cmd === "init" ? init(m.manifest)
          : m.cmd === "run" ? run(m.request, m.factor)
          : m.cmd === "probe" ? probe(m)
          : m.cmd === "calibrate" ? calibrate()
          : Promise.reject(new Error("unknown cmd " + m.cmd));
  job.catch(function (err) {
    self.postMessage({ type: "fatal", message: String(err && err.message || err) });
  });
};

async function init(manifest) {
  importScripts(PYODIDE_BASE + "pyodide.js");
  pyodide = await loadPyodide({ indexURL: PYODIDE_BASE });
  await pyodide.loadPackage(["micropip", "numpy", "sympy", "matplotlib"]);
  var micropip = pyodide.pyimport("micropip");
  await micropip.install(new URL(manifest.wheel, self.location.href).href);
  var resp = await fetch(new URL("runner.py", self.location.href));
  pyodide.FS.writeFile("runner.py", await resp.text());
  runner = pyodide.pyimport("runner");
  self.postMessage({ type: "ready", version: manifest.quiverlab_version });
}

async function calibrate() {
  var out = JSON.parse(runner.calibrate());
  self.postMessage({ type: "calibrated", factor: out.factor, seconds: out.seconds });
}

async function probe(m) {
  // A probe is a REAL build: true dimension + verbatim early errors.
  var built = JSON.parse(runner.run_build(JSON.stringify(m.request)));
  var data = built;
  if (built.ok) {
    var est = JSON.parse(runner.estimate(m.factor));
    if (est.ok) data = Object.assign({}, built, { eta: est });
  }
  self.postMessage({ type: "probe", seq: m.seq, data: data });
}

function remainingEta(breakdown, doneCount, factor) {
  var units = 0;
  for (var i = doneCount; i < breakdown.length; i++) units += breakdown[i].units;
  var b = JSON.parse(runner.bucket_for_seconds(units * factor));
  return { bucket: b.bucket, label: b.label, units: units };
}

async function run(request, factor) {
  // Pass the FULL request: python_snippet/result_bundle read request.compute
  // from the stored request (run_build itself only validates schema+algebra).
  var built = JSON.parse(runner.run_build(JSON.stringify(request)));
  var est = null;
  if (built.ok && typeof factor === "number") {
    var e0 = JSON.parse(runner.estimate(factor));
    if (e0.ok) est = e0;
  }
  self.postMessage({ type: "built", data: built, eta: est });
  if (built.ok) {
    for (var i = 0; i < request.compute.length; i++) {
      var t0 = performance.now();
      var res = JSON.parse(runner.compute_one(request.compute[i]));
      var elapsed = performance.now() - t0;
      var eta = null;
      if (est) {
        // Re-scale from measured speed (EMA) before bucketing the remainder.
        var u = est.breakdown[i].units;
        if (u >= 0.05) factor = 0.5 * factor + 0.5 * (elapsed / 1000 / u);
        eta = remainingEta(est.breakdown, i + 1, factor);
      }
      self.postMessage({ type: "result", data: res, elapsed_ms: elapsed, eta: eta });
    }
    self.postMessage({ type: "trace",
                       html: request.artifacts && request.artifacts.pdf
                             ? runner.trace_html() : "" });
    self.postMessage({ type: "artifacts", tikz: runner.tikz(),
                       snippet: runner.python_snippet(),
                       bundle: runner.result_bundle() });
  }
  self.postMessage({ type: "done" });
}
