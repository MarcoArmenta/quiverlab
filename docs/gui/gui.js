/* quiverlab landing-page GUI (Plan 10). Vanilla JS + SVG, zero dependencies.
   Loaded site-wide via extra_javascript: exits immediately unless #qlgui exists
   (the landing page only). Compute lives in gui/worker.js (Pyodide); this file
   is the editor, the request builder, and the results renderer. */
(function () {
  "use strict";
  var root = document.getElementById("qlgui");
  if (!root) return;

  var SVGNS = "http://www.w3.org/2000/svg";
  var R = 16;                       // vertex radius (px in SVG user units)
  var S = { vertices: [], arrows: [], nextId: 1, selected: null, dragFrom: null,
            dragMoved: false, dragOrigin: null,
            worker: null, engineReady: false, manifest: null, busy: false,
            artifacts: { tikz: "", snippet: "", bundle: "", traceHtml: "" } };

  // ---------- tiny DOM helpers ----------
  function h(tag, attrs) {
    var el = document.createElement(tag);
    for (var k in (attrs || {})) {
      if (k === "text") el.textContent = attrs[k];
      else el.setAttribute(k, attrs[k]);
    }
    for (var i = 2; i < arguments.length; i++) el.appendChild(arguments[i]);
    return el;
  }
  function sv(tag, attrs) {
    var el = document.createElementNS(SVGNS, tag);
    for (var k in (attrs || {})) el.setAttribute(k, attrs[k]);
    return el;
  }

  // ---------- static shell ----------
  root.innerHTML =
    '<div class="qlgui-row">' +
    '  <label>Preset <select id="qlgui-preset"><option value="">— build your own —</option></select></label>' +
    '  <label>Field <select id="qlgui-field"><option value="CC">CC</option><option value="GF">GF(p^n)</option></select></label>' +
    '  <label id="qlgui-p-wrap" style="display:none">p <input type="number" id="qlgui-p" value="2" min="2"></label>' +
    '  <label id="qlgui-n-wrap" style="display:none">n <input type="number" id="qlgui-n" value="1" min="1"></label>' +
    '  <button id="qlgui-clear" class="qlgui-secondary" type="button">Clear</button>' +
    '  <span id="qlgui-status">engine loads on first use</span>' +
    '</div>' +
    '<div id="qlgui-canvas-wrap">' +
    '  <svg id="qlgui-canvas" viewBox="0 0 800 340" preserveAspectRatio="xMidYMid meet"></svg>' +
    '  <input type="text" id="qlgui-rename">' +
    '</div>' +
    '<p class="qlgui-hint">Click empty space: add a vertex. Drag vertex → vertex: add an arrow ' +
    '(onto itself: a loop). Click: select. Double-click an arrow label: rename. Delete key: remove.</p>' +
    '<div class="qlgui-row"><label style="flex:1 1 260px">Relations ' +
    '<input type="text" id="qlgui-relations" placeholder="e.g. a*b - c, x*x*x"></label></div>' +
    '<div class="qlgui-row" id="qlgui-invariants">' +
    '  <label><input type="checkbox" id="qlgui-hhc" checked> HH^0..<select id="qlgui-hhc-top"></select></label>' +
    '  <label><input type="checkbox" id="qlgui-hhh"> HH_0..<select id="qlgui-hhh-top"></select></label>' +
    '  <label><input type="checkbox" id="qlgui-cartan" checked> Cartan matrix</label>' +
    '  <label><input type="checkbox" id="qlgui-coxeter_polynomial"> Coxeter polynomial</label>' +
    '  <label><input type="checkbox" id="qlgui-global_dimension"> gl.dim</label>' +
    '  <label><input type="checkbox" id="qlgui-center"> center</label>' +
    '  <label><input type="checkbox" id="qlgui-trace" checked> worked-steps report</label>' +
    '</div>' +
    '<div class="qlgui-row">' +
    '  <button id="qlgui-compute" type="button" disabled>Compute</button>' +
    '  <button id="qlgui-cancel" class="qlgui-secondary" type="button" disabled>Cancel</button>' +
    '  <button id="qlgui-print" class="qlgui-secondary" type="button" disabled>Print report (PDF)</button>' +
    '  <button id="qlgui-tikz" class="qlgui-secondary" type="button" disabled>TikZ</button>' +
    '  <button id="qlgui-json" class="qlgui-secondary" type="button" disabled>JSON</button>' +
    '  <button id="qlgui-snippet" class="qlgui-secondary" type="button" disabled>Copy Python</button>' +
    '</div>' +
    '<div id="qlgui-results"></div>';

  var el = {};
  ["preset", "field", "p-wrap", "n-wrap", "p", "n", "clear", "status", "canvas",
   "rename", "relations", "hhc", "hhc-top", "hhh", "hhh-top", "cartan",
   "coxeter_polynomial", "global_dimension", "center", "trace", "compute",
   "cancel", "print", "tikz", "json", "snippet", "results"]
    .forEach(function (id) { el[id] = document.getElementById("qlgui-" + id); });
  [el["hhc-top"], el["hhh-top"]].forEach(function (sel) {
    for (var i = 0; i <= 10; i++) sel.appendChild(h("option", { text: String(i) }));
    sel.value = "4";
  });

  function setStatus(text, cls) {
    el.status.textContent = text;
    el.status.className = cls || "";
  }

  // ---------- editor state ops ----------
  function nextArrowName() {
    var used = {};
    S.arrows.forEach(function (a) { used[a.name] = 1; });
    for (var suffix = 0; ; suffix++) {
      for (var c = 97; c <= 122; c++) {
        var n = String.fromCharCode(c) + (suffix ? String(suffix) : "");
        if (!used[n]) return n;
      }
    }
  }
  function vertexAt(id) {
    return S.vertices.filter(function (v) { return v.id === id; })[0];
  }
  function removeSelected() {
    if (!S.selected) return;
    if (S.selected.type === "vertex") {
      S.vertices = S.vertices.filter(function (v) { return v.id !== S.selected.key; });
      S.arrows = S.arrows.filter(function (a) {
        return a.s !== S.selected.key && a.t !== S.selected.key; });
    } else {
      S.arrows = S.arrows.filter(function (a) { return a.name !== S.selected.key; });
    }
    S.selected = null;
    render();
  }

  // ---------- geometry + render ----------
  function siblings(a) {   // arrows sharing the same UNORDERED vertex pair
    return S.arrows.filter(function (b) {
      return (b.s === a.s && b.t === a.t) || (b.s === a.t && b.t === a.s);
    });
  }
  function arrowPath(a) {
    var p = vertexAt(a.s), q = vertexAt(a.t);
    if (a.s === a.t) {                       // loop(s), stacked above the vertex
      var loops = S.arrows.filter(function (b) { return b.s === a.s && b.t === a.s; });
      var i = loops.indexOf(a), off = 54 + 34 * i;
      return { d: "M " + (p.x - 11) + " " + (p.y - 12) +
                  " C " + (p.x - 42) + " " + (p.y - off) + ", " +
                          (p.x + 42) + " " + (p.y - off) + ", " +
                          (p.x + 11) + " " + (p.y - 12),
               lx: p.x, ly: p.y - off + 6 };
    }
    var sib = siblings(a), i = sib.indexOf(a), n = sib.length;
    // Perpendicular from the CANONICAL (min,max) order so opposite-direction
    // arrows never collapse onto the same curve.
    var c0 = vertexAt(Math.min(a.s, a.t)), c1 = vertexAt(Math.max(a.s, a.t));
    var dx = c1.x - c0.x, dy = c1.y - c0.y, len = Math.sqrt(dx * dx + dy * dy) || 1;
    var px = -dy / len, py = dx / len;
    var k = 36 * (i - (n - 1) / 2);
    var mx = (p.x + q.x) / 2 + k * px, my = (p.y + q.y) / 2 + k * py;
    var ux = (q.x - p.x), uy = (q.y - p.y), ul = Math.sqrt(ux * ux + uy * uy) || 1;
    ux /= ul; uy /= ul;
    var sx = p.x + ux * (R + 2), sy = p.y + uy * (R + 2);
    var ex = q.x - ux * (R + 6), ey = q.y - uy * (R + 6);
    return { d: "M " + sx + " " + sy + " Q " + mx + " " + my + " " + ex + " " + ey,
             lx: 0.25 * sx + 0.5 * mx + 0.25 * ex,
             ly: 0.25 * sy + 0.5 * my + 0.25 * ey - 5 };
  }
  function render() {
    var svg = el.canvas;
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    var defs = sv("defs");
    var marker = sv("marker", { id: "qlgui-arrowhead", viewBox: "0 0 10 10",
      refX: "9", refY: "5", markerWidth: "7", markerHeight: "7", orient: "auto-start-reverse" });
    var head = sv("path", { d: "M 0 0 L 10 5 L 0 10 z" });
    head.style.fill = "currentColor";
    marker.appendChild(head); defs.appendChild(marker); svg.appendChild(defs);
    S.arrows.forEach(function (a) {
      var g = sv("g", { "class": "qla" +
        (S.selected && S.selected.type === "arrow" && S.selected.key === a.name ? " sel" : "") });
      var geo = arrowPath(a);
      var path = sv("path", { d: geo.d, "marker-end": "url(#qlgui-arrowhead)" });
      var label = sv("text", { x: geo.lx, y: geo.ly, "text-anchor": "middle" });
      label.textContent = a.name;
      [path, label].forEach(function (node) {
        node.addEventListener("click", function (e) {
          e.stopPropagation(); S.selected = { type: "arrow", key: a.name }; render();
        });
        node.addEventListener("dblclick", function (e) {
          e.stopPropagation(); startRename(a, geo);
        });
      });
      g.appendChild(path); g.appendChild(label); svg.appendChild(g);
    });
    S.vertices.forEach(function (v) {
      var g = sv("g", { "class": "qlv" +
        (S.selected && S.selected.type === "vertex" && S.selected.key === v.id ? " sel" : "") });
      var c = sv("circle", { cx: v.x, cy: v.y, r: R });
      c.addEventListener("mousedown", function (e) {
        e.preventDefault(); e.stopPropagation();
        S.dragFrom = v.id; S.dragMoved = false;
        S.dragOrigin = [e.clientX, e.clientY];
      });
      c.addEventListener("mouseup", function (e) {
        if (S.dragFrom === null) return;
        e.stopPropagation();
        // A press-release with no movement on the SAME vertex is a click
        // (select), not a self-loop; loops need a real drag gesture.
        if (S.dragFrom !== v.id || S.dragMoved) {
          S.arrows.push({ name: nextArrowName(), s: S.dragFrom, t: v.id });
          S.dragFrom = null; S.dragMoved = false; render();
        } else {
          S.dragFrom = null; S.dragMoved = false;
        }
      });
      c.addEventListener("click", function (e) {
        e.stopPropagation(); S.selected = { type: "vertex", key: v.id }; render();
      });
      var t = sv("text", { x: v.x, y: v.y + 5, "text-anchor": "middle" });
      t.textContent = String(v.id);
      g.appendChild(c); g.appendChild(t); svg.appendChild(g);
    });
    el.compute.disabled = !(S.engineReady && S.vertices.length && !S.busy);
  }

  function canvasPoint(e) {
    var pt = el.canvas.createSVGPoint();
    pt.x = e.clientX; pt.y = e.clientY;
    return pt.matrixTransform(el.canvas.getScreenCTM().inverse());
  }
  el.canvas.addEventListener("click", function (e) {
    if (e.target !== el.canvas) return;
    var p = canvasPoint(e);
    S.vertices.push({ id: S.nextId++, x: p.x, y: p.y });
    S.selected = null; render();
  });
  el.canvas.addEventListener("mousemove", function (e) {
    var old = el.canvas.querySelector(".qlgui-rubber");
    if (old) old.remove();
    if (S.dragFrom === null) return;
    if (S.dragOrigin && Math.hypot(e.clientX - S.dragOrigin[0],
                                   e.clientY - S.dragOrigin[1]) > 8) {
      S.dragMoved = true;
    }
    var v = vertexAt(S.dragFrom), p = canvasPoint(e);
    el.canvas.appendChild(sv("line", { "class": "qlgui-rubber",
      x1: v.x, y1: v.y, x2: p.x, y2: p.y }));
  });
  document.addEventListener("mouseup", function () {
    S.dragFrom = null; S.dragMoved = false; S.dragOrigin = null;
    var old = el.canvas.querySelector(".qlgui-rubber");
    if (old) old.remove();
  });
  document.addEventListener("keydown", function (e) {
    var tag = document.activeElement.tagName;
    if ((e.key === "Delete" || e.key === "Backspace") &&
        tag !== "INPUT" && tag !== "SELECT" && tag !== "TEXTAREA") { removeSelected(); }
  });

  // ---------- arrow rename ----------
  function startRename(arrow, geo) {
    var box = el.canvas.getBoundingClientRect();
    var scaleX = box.width / 800, scaleY = box.height / 340;
    el.rename.style.display = "block";
    el.rename.style.left = (geo.lx * scaleX - 24) + "px";
    el.rename.style.top = (geo.ly * scaleY - 12) + "px";
    el.rename.value = arrow.name;
    el.rename.focus(); el.rename.select();
    el.rename.onkeydown = function (e) {
      if (e.key === "Enter") commit();
      if (e.key === "Escape") {
        // Explicit blur: a display:none input silently keeps keyboard focus,
        // which makes the document-level Delete guard eat later deletions.
        el.rename.onblur = null;
        el.rename.style.display = "none";
        el.rename.blur();
      }
    };
    el.rename.onblur = commit;
    function commit() {
      var name = el.rename.value.trim();
      var taken = S.arrows.some(function (b) { return b !== arrow && b.name === name; });
      if (/^[A-Za-z][A-Za-z0-9_]*$/.test(name) && !taken) {
        arrow.name = name;
        el.rename.onblur = null;      // avoid blur() re-entering commit
        el.rename.style.display = "none";
        el.rename.blur();
        render();
      } else { el.rename.style.borderColor = "#c62828"; }
    }
  }

  el.clear.addEventListener("click", function () {
    S.vertices = []; S.arrows = []; S.nextId = 1; S.selected = null;
    el.relations.value = ""; el.results.innerHTML = ""; render();
  });
  el.field.addEventListener("change", function () {
    var gf = el.field.value === "GF";
    el["p-wrap"].style.display = gf ? "" : "none";
    el["n-wrap"].style.display = gf ? "" : "none";
  });

  // ---------- presets ----------
  fetch("gui/presets.json").then(function (r) { return r.ok ? r.json() : []; })
    .then(function (presets) {
      presets.forEach(function (p, i) {
        el.preset.appendChild(h("option", { value: String(i), text: p.label }));
      });
      el.preset.addEventListener("change", function () {
        if (el.preset.value === "") return;
        var p = presets[parseInt(el.preset.value, 10)];
        S.vertices = p.vertices.map(function (id, i) {
          var angle = 2 * Math.PI * i / p.vertices.length - Math.PI / 2;
          var rad = p.vertices.length === 1 ? 0 : 110;
          return { id: id, x: 400 + rad * Math.cos(angle), y: 185 + rad * Math.sin(angle) };
        });
        S.nextId = Math.max.apply(null, p.vertices.concat([0])) + 1;
        S.arrows = Object.keys(p.arrows).map(function (name) {
          return { name: name, s: p.arrows[name][0], t: p.arrows[name][1] };
        });
        el.relations.value = p.relations.join(", ");
        el.field.value = p.field.kind === "CC" ? "CC" : "GF";
        el.field.dispatchEvent(new Event("change"));
        if (p.field.kind === "GF") { el.p.value = p.field.p; el.n.value = p.field.n || 1; }
        S.selected = null; el.results.innerHTML = ""; render();
      });
    }).catch(function () { /* presets are a convenience; the editor still works */ });

  // ---------- request ----------
  function buildRequest() {
    var arrows = {};
    S.arrows.forEach(function (a) { arrows[a.name] = [a.s, a.t]; });
    var relations = el.relations.value.split(",")
      .map(function (s) { return s.trim(); }).filter(Boolean);
    var field = el.field.value === "CC" ? { kind: "CC" }
      : { kind: "GF", p: parseInt(el.p.value, 10) || 0, n: parseInt(el.n.value, 10) || 1 };
    var compute = [];
    if (el.hhc.checked) compute.push("hh_cohomology:0.." + el["hhc-top"].value);
    if (el.hhh.checked) compute.push("hh_homology:0.." + el["hhh-top"].value);
    ["cartan", "coxeter_polynomial", "global_dimension", "center"].forEach(function (k) {
      if (el[k].checked) compute.push(k);
    });
    return { schema: 1,
             algebra: { kind: "quiver",
                        vertices: S.vertices.map(function (v) { return v.id; }),
                        arrows: arrows, relations: relations, field: field },
             compute: compute,
             artifacts: { pdf: el.trace.checked, tikz: true } };
  }

  // ---------- engine (Pyodide worker) ----------
  function startWorker() {
    if (S.worker) S.worker.terminate();
    S.engineReady = false;
    S.worker = new Worker("gui/worker.js");
    S.worker.onmessage = onWorkerMessage;
    S.worker.onerror = function (e) {
      setStatus("engine error — see console", "err");
      console.error("qlgui worker:", e);
    };
    if (S.manifest) {                    // restart path (after Cancel)
      setStatus("engine reloading…");
      S.worker.postMessage({ cmd: "init", manifest: S.manifest });
      return;
    }
    fetch("gui/manifest.json")
      .then(function (r) { if (!r.ok) throw new Error("no manifest"); return r.json(); })
      .then(function (m) {
        if (!m.wheel) { setStatus("engine payload not built (QLGUI_SKIP_WHEEL)", "err"); return; }
        S.manifest = m;
        setStatus("engine loading… (~60 MB once, then cached)");
        S.worker.postMessage({ cmd: "init", manifest: m });
      })
      .catch(function () { setStatus("engine manifest missing — editor-only preview", "err"); });
  }

  function setBusy(b) {
    S.busy = b;
    el.cancel.disabled = !b;
    render();
  }

  function onWorkerMessage(e) {
    var m = e.data;
    if (m.type === "ready") {
      S.engineReady = true;
      setStatus("engine ready — quiverlab " + m.version, "ok");
      render();
    } else if (m.type === "built") {
      if (m.data.ok) {
        el.results.appendChild(h("div", { "class": "qlgui-block",
          text: m.data.algebra }));
      } else { renderError(m.data); }
    } else if (m.type === "result") {
      if (m.data.ok) renderBlock(m.data); else renderError(m.data);
    } else if (m.type === "trace") {
      S.artifacts.traceHtml = m.html;
      el.print.disabled = !m.html;
    } else if (m.type === "artifacts") {
      S.artifacts.tikz = m.tikz; S.artifacts.snippet = m.snippet;
      S.artifacts.bundle = m.bundle;
      el.tikz.disabled = el.json.disabled = el.snippet.disabled = false;
    } else if (m.type === "done") {
      setBusy(false);
      setStatus("engine ready — quiverlab " + (S.manifest ? S.manifest.quiverlab_version : ""), "ok");
    } else if (m.type === "fatal") {
      console.error("qlgui engine:", m.message);
      setStatus("engine failed — see console", "err");
      setBusy(false);
    }
  }

  function renderError(res) {
    if (res.detail) console.error(res.detail);
    var div = h("div", { "class": "qlgui-block qlgui-error",
      text: res.error.type + ": " + res.error.message });
    el.results.appendChild(div);
  }

  function citesLine(block) {
    return h("div", { "class": "qlgui-cites",
      text: (block.citations || []).map(function (c) { return c[1]; }).join(" · ") });
  }

  function renderBlock(res) {
    var b = res.block, name = res.invariant.split(":")[0];
    var div = h("div", { "class": "qlgui-block" });
    if (name === "hh_cohomology" || name === "hh_homology") {
      var sup = name === "hh_cohomology";
      var head = h("tr"), row = h("tr");
      head.appendChild(h("th", { text: "n" }));
      row.appendChild(h("th", { text: sup ? "dim HH^n" : "dim HH_n" }));
      b.dims.forEach(function (d, n) {
        head.appendChild(h("td", { text: String(n) }));
        row.appendChild(h("td", { text: String(d) }));
      });
      div.appendChild(h("p", { text: sup ? "Hochschild cohomology" : "Hochschild homology" }));
      div.appendChild(h("table", {}, head, row));
      div.appendChild(h("div", { "class": "qlgui-cites", text: b.engine }));
    } else if (name === "cartan") {
      div.appendChild(h("p", { text: "Cartan matrix:" }));
      // className EXACTLY "arithmatex": the site's MathJax config matches
      // class patterns against the full className string, so a combined
      // "qlgui-block arithmatex" is silently skipped (found live).
      div.appendChild(h("p", { "class": "arithmatex", text: "\\[ C = " + b.latex + " \\]" }));
    } else if (name === "coxeter_polynomial") {
      div.appendChild(h("p", { "class": "arithmatex", text: "\\[ \\chi(t) = " + b.latex + " \\]" }));
    } else if (name === "global_dimension") {
      div.appendChild(h("p", { text: b.text }));
    } else if (name === "center") {
      div.appendChild(h("p", { "class": "arithmatex", text: "\\( \\dim Z(A) = " + b.dim + " \\)" }));
    }
    div.appendChild(citesLine(b));
    el.results.appendChild(div);
    if (window.MathJax && window.MathJax.typesetPromise) {
      // Full-page sweep, NOT typesetPromise([div]): with explicit roots the
      // walker never consults the root's own class, so the site's
      // ignoreHtmlClass ".*|" config silently skips the block (found live).
      window.MathJax.typesetPromise();
    }
  }

  // ---------- buttons ----------
  el.compute.addEventListener("click", function () {
    if (S.busy || !S.engineReady) return;
    el.results.innerHTML = "";
    S.artifacts = { tikz: "", snippet: "", bundle: "", traceHtml: "" };
    el.print.disabled = el.tikz.disabled = el.json.disabled = el.snippet.disabled = true;
    setBusy(true);
    setStatus("computing…");
    S.worker.postMessage({ cmd: "run", request: buildRequest() });
  });
  el.cancel.addEventListener("click", function () {
    if (!S.busy) return;
    setBusy(false);
    startWorker();                             // sets its own transient status…
    setStatus("cancelled — engine reloading…"); // …so the acknowledgment must land last
    render();                                   // engineReady just went false — disable Compute
  });
  function download(name, text, type) {
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], { type: type }));
    a.download = name;
    a.click();
    URL.revokeObjectURL(a.href);
  }
  el.print.addEventListener("click", function () {
    var url = URL.createObjectURL(new Blob([S.artifacts.traceHtml], { type: "text/html" }));
    var w = window.open(url, "_blank");
    if (w) w.addEventListener("load", function () { w.print(); });
  });
  el.tikz.addEventListener("click", function () {
    download("quiver.tex", S.artifacts.tikz, "text/x-tex");
  });
  el.json.addEventListener("click", function () {
    download("quiverlab-result.json", S.artifacts.bundle, "application/json");
  });
  el.snippet.addEventListener("click", function () {
    navigator.clipboard.writeText(S.artifacts.snippet).then(function () {
      setStatus("Python snippet copied", "ok");
    });
  });

  window.QLGUI = { S: S, buildRequest: buildRequest };
  render();
  // Engine loads on FIRST INTENT (whole-branch review decision): pure readers
  // never pay the ~60 MB download; the first GUI touch starts it.
  var engineStarted = false;
  function ensureEngine() {
    if (engineStarted) return;
    engineStarted = true;
    startWorker();
  }
  el.canvas.addEventListener("mousedown", ensureEngine, true); // capture: circle handlers stopPropagation
  el.preset.addEventListener("change", ensureEngine);
  el.relations.addEventListener("focus", ensureEngine);
})();
