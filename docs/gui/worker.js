/* quiverlab GUI Web Worker (Plan 10): Pyodide + the site's own wheel.
   PINNED Pyodide version — never a floating tag (supply-chain posture of
   mkdocs.yml). Protocol: see gui.js. All runner I/O is JSON strings. */
"use strict";
var PYODIDE_BASE = "https://cdn.jsdelivr.net/pyodide/v0.28.3/full/";
var pyodide = null, runner = null;

self.onmessage = function (e) {
  var m = e.data;
  var job = m.cmd === "init" ? init(m.manifest)
          : m.cmd === "run" ? run(m.request)
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

async function run(request) {
  // Pass the FULL request: python_snippet/result_bundle read request.compute
  // from the stored request (run_build itself only validates schema+algebra).
  var built = JSON.parse(runner.run_build(JSON.stringify(request)));
  self.postMessage({ type: "built", data: built });
  if (built.ok) {
    for (var i = 0; i < request.compute.length; i++) {
      var res = JSON.parse(runner.compute_one(request.compute[i]));
      self.postMessage({ type: "result", data: res });
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
