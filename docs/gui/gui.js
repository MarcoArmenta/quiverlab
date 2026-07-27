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
            dragMoved: false, dragOrigin: null, pressOnEmpty: false,
            worker: null, engineReady: false, manifest: null, busy: false,
            artifacts: { tikz: "", snippet: "", bundle: "", traceHtml: "", traceTex: "", traceJson: "" },
            // Plan 26 no-code module panel: dims/maps hold the explicit module the
            // user types (matrix entries are exact strings); side/vertex track the
            // toggle + builtin pick-list. All read back in buildRequest().
            module: { enabled: false, side: "right", vertex: null,
                      dims: {}, maps: {} },
            // Plan 30: the SECOND argument N (Ext/Tor target). Same no-code editor
            // as the main module -- an explicit dims+matrices module OR an S/P/I
            // pick-list, over the same quiver. Tor forces N to a LEFT A-module.
            target: { mode: "simple", side: "right", vertex: null,
                      dims: {}, maps: {} } };

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
    // ---- Plan 26: no-code module panel ----
    '<fieldset id="qlgui-module" class="qlgui-fieldset">' +
    '  <legend><label><input type="checkbox" id="qlgui-mod-enable"> Module (no code)</label></legend>' +
    '  <div class="qlgui-row">' +
    '    <label>build <select id="qlgui-mod-mode">' +
    '      <option value="explicit">explicit (dims + matrices)</option>' +
    '      <option value="simple">S(v) simple</option>' +
    '      <option value="projective">P(v) projective</option>' +
    '      <option value="injective">I(v) injective</option>' +
    '    </select></label>' +
    '    <label>side <select id="qlgui-mod-side">' +
    '      <option value="right">right</option><option value="left">left</option>' +
    '    </select></label>' +
    '  </div>' +
    '  <div id="qlgui-mod-body"></div>' +
    '  <div class="qlgui-row" id="qlgui-mod-kinds">' +
    '    <label><input type="checkbox" id="qlgui-dimension_vector" checked> dim vector</label>' +
    '    <label><input type="checkbox" id="qlgui-rad_top_soc"> rad/top/soc</label>' +
    '    <label><input type="checkbox" id="qlgui-tau"> τ</label>' +
    '    <label><input type="checkbox" id="qlgui-tau_minus"> τ⁻</label>' +
    '    <label><input type="checkbox" id="qlgui-projective_dimension"> proj.dim</label>' +
    '    <label><input type="checkbox" id="qlgui-injective_dimension"> inj.dim</label>' +
    '    <label><input type="checkbox" id="qlgui-projective_resolution"> proj.res 0..<select id="qlgui-pr-top"></select></label>' +
    '    <label><input type="checkbox" id="qlgui-injective_resolution"> inj.res 0..<select id="qlgui-ir-top"></select></label>' +
    '    <label><input type="checkbox" id="qlgui-decompose"> decompose</label>' +
    '    <label><input type="checkbox" id="qlgui-ext"> Ext 0..<select id="qlgui-ext-top"></select></label>' +
    '    <label><input type="checkbox" id="qlgui-tor"> Tor 0..<select id="qlgui-tor-top"></select></label>' +
    '  </div>' +
    // ---- Plan 30: second-argument N editor (Ext/Tor target) ----
    '  <fieldset id="qlgui-target" class="qlgui-fieldset" style="display:none">' +
    '    <legend>second argument N (Ext / Tor target)</legend>' +
    '    <div class="qlgui-row">' +
    '      <label>build <select id="qlgui-target-mode">' +
    '        <option value="explicit">explicit (dims + matrices)</option>' +
    '        <option value="simple">S(v) simple</option>' +
    '        <option value="projective">P(v) projective</option>' +
    '        <option value="injective">I(v) injective</option>' +
    '      </select></label>' +
    '      <label>side <select id="qlgui-target-side">' +
    '        <option value="right">right</option><option value="left">left</option>' +
    '      </select></label>' +
    '    </div>' +
    '    <div id="qlgui-target-body"></div>' +
    '    <p class="qlgui-hint" id="qlgui-target-note"></p>' +
    '  </fieldset>' +
    '</fieldset>' +
    '<div id="qlgui-eta" class="qlgui-hint"></div>' +
    '<div class="qlgui-row">' +
    '  <button id="qlgui-compute" type="button" disabled>Compute</button>' +
    '  <button id="qlgui-cancel" class="qlgui-secondary" type="button" disabled>Cancel</button>' +
    '  <button id="qlgui-print" class="qlgui-secondary" type="button" disabled title="Open the typeset report and print it to PDF from your browser">Print report</button>' +
    '  <button id="qlgui-report-html" class="qlgui-secondary" type="button" disabled title="Download the print-ready HTML report (print to PDF from your browser)">Report (HTML)</button>' +
    '  <button id="qlgui-worked-tex" class="qlgui-secondary" type="button" disabled title="Download the LaTeX source and compile it yourself">Report (TeX)</button>' +
    '  <button id="qlgui-report-json" class="qlgui-secondary" type="button" disabled title="Download the complete worked-steps event stream (exact, machine-readable JSON)">Report data (JSON)</button>' +
    '  <button id="qlgui-tikz" class="qlgui-secondary" type="button" disabled>TikZ</button>' +
    '  <button id="qlgui-json" class="qlgui-secondary" type="button" disabled>JSON</button>' +
    '  <button id="qlgui-snippet" class="qlgui-secondary" type="button" disabled>Copy Python</button>' +
    '  <button id="qlgui-config" class="qlgui-secondary" type="button" disabled>Config (YAML)</button>' +
    '</div>' +
    '<div id="qlgui-results"></div>';

  var el = {};
  ["preset", "field", "p-wrap", "n-wrap", "p", "n", "clear", "status", "canvas",
   "rename", "relations", "hhc", "hhc-top", "hhh", "hhh-top", "cartan",
   "coxeter_polynomial", "global_dimension", "center", "trace", "compute",
   "cancel", "print", "report-html", "worked-tex", "report-json", "tikz", "json", "snippet", "config", "results", "eta",
   // Plan 26 module panel + Plan 30 (tor / decompose / second-argument editor)
   "module", "mod-enable", "mod-mode", "mod-side", "mod-body", "mod-kinds",
   "dimension_vector", "rad_top_soc", "tau", "tau_minus",
   "projective_dimension", "injective_dimension",
   "projective_resolution", "pr-top", "injective_resolution", "ir-top",
   "decompose", "ext", "ext-top", "tor", "tor-top",
   "target", "target-mode", "target-side", "target-body", "target-note"]
    .forEach(function (id) { el[id] = document.getElementById("qlgui-" + id); });
  [el["hhc-top"], el["hhh-top"], el["pr-top"], el["ir-top"], el["ext-top"],
   el["tor-top"]]
    .forEach(function (sel) {
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
        // A gesture that begins on a vertex must never spawn a vertex on the
        // trailing synthesized click (stopPropagation keeps the canvas
        // mousedown below from resetting this flag for us).
        S.pressOnEmpty = false;
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
    // Config export needs only the request builder (pure, client-side): usable as
    // soon as there is a quiver, WITHOUT the engine -- design here, run on a cluster.
    el.config.disabled = !(S.vertices.length && !S.busy);
    renderModulePanel();
    scheduleProbe();
  }

  // ---------- Plan 26: no-code module panel ----------
  var MOD_KIND_IDS = ["dimension_vector", "rad_top_soc", "tau", "tau_minus",
    "projective_dimension", "injective_dimension",
    "projective_resolution", "injective_resolution", "decompose", "ext", "tor"];

  // Generic matrix-editor helpers over a module-state {dims, maps} + a side. Used
  // by BOTH the main module panel (S.module) and the second-argument editor
  // (S.target), so the target gets the SAME explicit dims+matrices editor (Plan 30).
  function stDim(mstate, id) { var n = mstate.dims[id]; return (n == null) ? 0 : n; }

  function syncDims(mstate) {    // one dim per current vertex (default 1); drop stale
    var d = {};
    S.vertices.forEach(function (v) {
      d[v.id] = (mstate.dims[v.id] != null) ? mstate.dims[v.id] : 1;
    });
    mstate.dims = d;
  }

  function matrixDimsFor(a, mstate, side) {   // [rows, cols] of arrow a's block, per side
    var ds = stDim(mstate, a.s), dt = stDim(mstate, a.t);
    // right: a:s→t acts M_s→M_t, block is dim[t]×dim[s]; left is over A^op, so the
    // block is transposed (dim[s]×dim[t]) — the runner places it by the rep quiver.
    return (side === "left") ? [ds, dt] : [dt, ds];
  }

  function syncMaps(mstate, side) {   // rows×cols string grid per arrow, preserving overlap
    var maps = {};
    S.arrows.forEach(function (a) {
      var rc = matrixDimsFor(a, mstate, side), rows = rc[0], cols = rc[1];
      var old = mstate.maps[a.name] || [];
      var grid = [];
      for (var i = 0; i < rows; i++) {
        var row = [];
        for (var j = 0; j < cols; j++) {
          row.push((old[i] && old[i][j] != null) ? old[i][j] : "0");
        }
        grid.push(row);
      }
      maps[a.name] = grid;
    });
    mstate.maps = maps;
  }

  // Thin S.module-bound wrappers (keep the main panel's existing call sites terse).
  function modDim(id) { return stDim(S.module, id); }
  function syncModuleDims() { syncDims(S.module); }
  function matrixDims(a) { return matrixDimsFor(a, S.module, S.module.side); }
  function syncModuleMaps() { syncMaps(S.module, S.module.side); }

  // Render the explicit dims-picker + per-arrow matrix grids for `mstate`/`side`
  // into `bodyEl`. `onEdit` fires on every dim change (rebuild) or cell edit.
  function renderExplicitBody(bodyEl, mstate, side, onDimChange) {
    syncDims(mstate); syncMaps(mstate, side);
    var dimRow = h("div", { "class": "qlgui-row" });
    S.vertices.forEach(function (v) {
      var inp = h("input", { type: "number", min: "0",
                             value: String(stDim(mstate, v.id)) });
      inp.addEventListener("input", function () {
        mstate.dims[v.id] = Math.max(0, parseInt(inp.value, 10) || 0);
      });
      inp.addEventListener("change", onDimChange);
      dimRow.appendChild(h("label", { text: "dim " + v.id + " " }, inp));
    });
    bodyEl.appendChild(h("div", { "class": "qlgui-mrow" },
      h("span", { "class": "qlgui-mlabel", text: "dimension vector" }), dimRow));
    S.arrows.forEach(function (a) {
      var rc = matrixDimsFor(a, mstate, side), rows = rc[0], cols = rc[1];
      var grid = h("div", { "class": "qlgui-mgrid" });
      grid.style.gridTemplateColumns = "repeat(" + Math.max(cols, 1) + ", 3.2em)";
      if (rows === 0 || cols === 0) {
        grid.appendChild(h("span", { "class": "qlgui-hint",
          text: rows + "×" + cols + " (empty block)" }));
      }
      for (var i = 0; i < rows; i++) {
        for (var j = 0; j < cols; j++) {
          (function (ii, jj) {
            var cell = h("input", { type: "text", value: mstate.maps[a.name][ii][jj] });
            cell.addEventListener("input", function () {
              mstate.maps[a.name][ii][jj] = cell.value; scheduleProbe();
            });
            grid.appendChild(cell);
          })(i, j);
        }
      }
      bodyEl.appendChild(h("div", { "class": "qlgui-mrow" },
        h("span", { "class": "qlgui-mlabel",
          text: a.name + ": " + a.s + "→" + a.t + "  [" + rows + "×" + cols + "]" }), grid));
    });
  }

  function fillVertexOptions(sel, current) {
    sel.innerHTML = "";
    S.vertices.forEach(function (v) {
      sel.appendChild(h("option", { value: String(v.id), text: String(v.id) }));
    });
    if (current != null && S.vertices.some(function (v) { return v.id === current; })) {
      sel.value = String(current);
    } else if (S.vertices.length) { sel.value = String(S.vertices[0].id); }
  }

  function renderModulePanel() {
    var on = el["mod-enable"].checked;
    S.module.enabled = on;
    S.module.side = el["mod-side"].value;
    el["mod-mode"].disabled = el["mod-side"].disabled = !on;
    MOD_KIND_IDS.forEach(function (k) { if (el[k]) el[k].disabled = !on; });
    ["pr-top", "ir-top", "ext-top", "tor-top", "target-mode", "target-side"]
      .forEach(function (k) { if (el[k]) el[k].disabled = !on; });
    var body = el["mod-body"];
    body.innerHTML = "";
    if (on) {
      var mode = el["mod-mode"].value;
      if (mode !== "explicit") {                  // S(v)/P(v)/I(v) pick-list
        var vsel = h("select", {});
        fillVertexOptions(vsel, S.module.vertex);
        S.module.vertex = parseInt(vsel.value, 10);
        vsel.addEventListener("change", function () {
          S.module.vertex = parseInt(vsel.value, 10); scheduleProbe();
        });
        body.appendChild(h("div", { "class": "qlgui-row" },
          h("label", { text: mode + " at vertex " }, vsel)));
      } else if (!S.vertices.length) {
        body.appendChild(h("p", { "class": "qlgui-hint",
          text: "add vertices on the canvas to define the module" }));
      } else {
        renderExplicitBody(body, S.module, S.module.side, function () {
          renderModulePanel(); scheduleProbe();
        });
      }
    }
    renderTargetPanel(on);
  }

  // ---------- Plan 30: the second-argument N editor (Ext/Tor target) ----------
  function renderTargetPanel(moduleOn) {
    // Shown only when the module panel is on AND Ext or Tor is requested.
    var wantExt = moduleOn && el.ext.checked, wantTor = moduleOn && el.tor.checked;
    var show = wantExt || wantTor;
    el.target.style.display = show ? "" : "none";
    var note = el["target-note"];
    // Tor's N must be a LEFT A-module: force + lock the side toggle when Tor is on
    // (and Ext is not competing for a right N). The note states it honestly.
    if (wantTor && !wantExt) {
      el["target-side"].value = "left";
      el["target-side"].disabled = true;
      S.target.side = "left";
      note.textContent = "Tor's second argument N is a LEFT A-module (side forced left).";
    } else {
      el["target-side"].disabled = !moduleOn;
      S.target.side = el["target-side"].value;
      note.textContent = wantExt && wantTor
        ? "Ext uses N as typed; Tor reads the same N as a LEFT module."
        : "";
    }
    var body = el["target-body"];
    body.innerHTML = "";
    if (!show) return;
    var mode = el["target-mode"].value;
    if (mode !== "explicit") {                    // S(v)/P(v)/I(v) pick-list
      var vsel = h("select", {});
      fillVertexOptions(vsel, S.target.vertex);
      S.target.vertex = parseInt(vsel.value, 10);
      vsel.addEventListener("change", function () {
        S.target.vertex = parseInt(vsel.value, 10); scheduleProbe();
      });
      body.appendChild(h("div", { "class": "qlgui-row" },
        h("label", { text: mode + " at vertex " }, vsel)));
      return;
    }
    if (!S.vertices.length) {
      body.appendChild(h("p", { "class": "qlgui-hint",
        text: "add vertices on the canvas to define N" }));
      return;
    }
    renderExplicitBody(body, S.target, S.target.side, function () {
      renderModulePanel(); scheduleProbe();
    });
  }

  function normGrid(g) {        // trim; an empty cell means 0
    return (g || []).map(function (row) {
      return row.map(function (x) {
        var s = (x == null ? "" : String(x)).trim();
        return s === "" ? "0" : s;
      });
    });
  }

  function moduleSpec() {
    S.module.side = el["mod-side"].value;
    var mode = el["mod-mode"].value, side = S.module.side;
    if (mode !== "explicit") {
      return { builtin: { kind: mode, vertex: S.module.vertex }, side: side };
    }
    syncModuleDims(); syncModuleMaps();
    var dims = {}, maps = {};
    S.vertices.forEach(function (v) { dims[String(v.id)] = modDim(v.id); });
    S.arrows.forEach(function (a) { maps[a.name] = normGrid(S.module.maps[a.name]); });
    return { dims: dims, maps: maps, side: side };
  }

  // The second argument N, read from the S.target editor at the requested side.
  // Same no-code surface as the main module: explicit dims+matrices OR an S/P/I
  // pick-list (Plan 30). Ext reads N at S.target.side; Tor forces N to a LEFT module.
  function targetSpecWith(side) {
    var mode = el["target-mode"].value;
    var fallbackV = S.vertices.length ? S.vertices[0].id : null;
    if (mode !== "explicit") {
      return { builtin: { kind: mode,
                          vertex: (S.target.vertex != null ? S.target.vertex : fallbackV) },
               side: side };
    }
    syncDims(S.target); syncMaps(S.target, side);
    var dims = {}, maps = {};
    S.vertices.forEach(function (v) { dims[String(v.id)] = stDim(S.target, v.id); });
    S.arrows.forEach(function (a) { maps[a.name] = normGrid(S.target.maps[a.name]); });
    return { dims: dims, maps: maps, side: side };
  }
  function extTargetSpec() { return targetSpecWith(S.target.side); }
  function torTargetSpec() { return targetSpecWith("left"); }

  function canvasPoint(e) {
    var pt = el.canvas.createSVGPoint();
    pt.x = e.clientX; pt.y = e.clientY;
    return pt.matrixTransform(el.canvas.getScreenCTM().inverse());
  }
  // Record whether a press began on genuinely empty canvas. Dragging vertex →
  // vertex ends the vertex mouseup with a render() that detaches the pressed
  // circles, so the browser synthesizes a click on their nearest surviving
  // ancestor — the <svg>, i.e. el.canvas — which would otherwise pass the guard
  // below and spawn a phantom vertex under the release point. Gating on this
  // flag means only a real empty-canvas press adds a vertex.
  el.canvas.addEventListener("mousedown", function (e) {
    S.pressOnEmpty = (e.target === el.canvas);
  });
  el.canvas.addEventListener("click", function (e) {
    if (e.target !== el.canvas || !S.pressOnEmpty) return;
    S.pressOnEmpty = false;
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
    S.module = { enabled: false, side: "right", vertex: null, dims: {}, maps: {} };
    S.target = { mode: "simple", side: "right", vertex: null, dims: {}, maps: {} };
    el["mod-enable"].checked = false; el["mod-mode"].value = "explicit";
    el["mod-side"].value = "right";
    el["target-mode"].value = "explicit"; el["target-side"].value = "right";
    el.ext.checked = false; el.tor.checked = false;
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
    var module = null, extTarget = null, torTarget = null;
    if (el["mod-enable"].checked) {          // read live, independent of render timing
      module = moduleSpec();
      ["dimension_vector", "rad_top_soc", "tau", "tau_minus",
       "projective_dimension", "injective_dimension", "decompose"].forEach(function (k) {
        if (el[k].checked) compute.push(k);
      });
      if (el.projective_resolution.checked)
        compute.push("projective_resolution:0.." + el["pr-top"].value);
      if (el.injective_resolution.checked)
        compute.push("injective_resolution:0.." + el["ir-top"].value);
      if (el.ext.checked) {
        compute.push("ext:0.." + el["ext-top"].value);
        extTarget = extTargetSpec();
      }
      if (el.tor.checked) {
        compute.push("tor:0.." + el["tor-top"].value);
        torTarget = torTargetSpec();
      }
    }
    var req = { schema: 1,
                algebra: { kind: "quiver",
                           vertices: S.vertices.map(function (v) { return v.id; }),
                           arrows: arrows, relations: relations, field: field },
                compute: compute,
                artifacts: { pdf: el.trace.checked, tikz: true } };
    if (module) req.module = module;
    if (extTarget) req.ext_target = extTarget;
    if (torTarget) req.tor_target = torTarget;
    return req;
  }

  // Hand-rolled YAML emitter for the exported cluster config. MIRRORS
  // webapp/server/clusterconfig.py: both must round-trip through a YAML parser
  // back to the same compute request (the export test pins the round-trip, not
  // byte-identity). Self-contained + pure so it can be extracted and node-tested.
  // Handles the shallow request shape: scalars, block lists, nested dicts, and
  // matrix (list-of-list) module maps.
  function configYaml(obj) {
    function scalar(v) {
      if (v === null || v === undefined) return "null";
      if (typeof v === "boolean") return v ? "true" : "false";
      if (typeof v === "number") return String(v);
      var s = String(v);
      // Quote unless it is an unambiguous plain scalar. A pure integer string
      // ("0", "-3") stays plain (round-trips to an int, which the schema accepts);
      // everything else (empty, YAML-special chars, leading space, reserved words,
      // fractions like "1/2") is double-quoted (JSON strings are valid YAML).
      var plainInt = /^-?[0-9]+$/.test(s);
      if (!plainInt && (s === "" || /[:#\-?\[\]{}&*!|>'"%@`,]/.test(s) ||
          /^\s|\s$/.test(s) || /^(true|false|null|yes|no|on|off|~)$/i.test(s) ||
          /^[0-9.+]/.test(s))) {
        return JSON.stringify(s);
      }
      return s;
    }
    function keyScalar(k) {
      // Keys must keep their STRING type through the round-trip: quote anything a
      // YAML parser would otherwise read as a non-string. Notably module dims keys
      // ("1","2") are numeric strings and MUST be quoted (the schema's dims keys
      // are strings), unlike matrix VALUES where a plain int is acceptable.
      var s = String(k);
      if (s === "" || /^-?[0-9]+$/.test(s) || /[:#\-?\[\]{}&*!|>'"%@`,]/.test(s) ||
          /^\s|\s$/.test(s) || /^(true|false|null|yes|no|on|off|~)$/i.test(s) ||
          /^[0-9.+]/.test(s)) {
        return JSON.stringify(s);
      }
      return s;
    }
    function emit(o, indent) {
      var pad = new Array(indent + 1).join("  "), out = "";
      if (Array.isArray(o)) {
        if (o.length === 0) return pad + "[]\n";
        o.forEach(function (item) {
          if (item !== null && typeof item === "object") {
            out += pad + "-\n" + emit(item, indent + 1);   // nested list/dict
          } else {
            out += pad + "- " + scalar(item) + "\n";
          }
        });
        return out;
      }
      if (o !== null && typeof o === "object") {
        Object.keys(o).sort().forEach(function (k) {
          var v = o[k], kk = keyScalar(k),
            nonEmptyObj = v !== null && typeof v === "object" &&
            (Array.isArray(v) ? v.length : Object.keys(v).length);
          if (nonEmptyObj) {
            out += pad + kk + ":\n" + emit(v, indent + 1);
          } else if (v !== null && typeof v === "object") {
            out += pad + kk + ": " + (Array.isArray(v) ? "[]" : "{}") + "\n";
          } else {
            out += pad + kk + ": " + scalar(v) + "\n";
          }
        });
        return out;
      }
      return pad + scalar(o) + "\n";
    }
    return emit(obj, 0);
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

  // ---------- wait estimates (Plan 11) ----------
  S.factor = null; S.probeSeq = 0; S.eta = null;
  S.tickerId = null; S.computeT0 = 0;
  var probeTimer = null;

  function setEta(text, isErr) {
    el.eta.textContent = text;
    el.eta.className = "qlgui-hint" + (isErr ? " err" : "");
  }

  function scheduleProbe() {
    if (probeTimer) clearTimeout(probeTimer);
    probeTimer = setTimeout(function () {
      probeTimer = null;
      if (!S.engineReady || S.busy || S.factor === null || !S.vertices.length) return;
      S.probeSeq++;
      S.worker.postMessage({ cmd: "probe", request: buildRequest(),
                             seq: S.probeSeq, factor: S.factor });
    }, 600);
  }

  function startTicker() {
    stopTicker();
    S.computeT0 = Date.now();
    S.tickerId = setInterval(function () {
      var secs = Math.round((Date.now() - S.computeT0) / 1000);
      setStatus("computing… · " + secs + " s elapsed" +
                (S.eta ? " · " + S.eta.label : ""));
    }, 1000);
  }

  function stopTicker() {
    if (S.tickerId) { clearInterval(S.tickerId); S.tickerId = null; }
  }

  function onWorkerMessage(e) {
    var m = e.data;
    if (m.type === "ready") {
      S.engineReady = true;
      setStatus("engine ready — quiverlab " + m.version, "ok");
      render();
      S.worker.postMessage({ cmd: "calibrate" });
    } else if (m.type === "calibrated") {
      S.factor = m.factor;
      scheduleProbe();
    } else if (m.type === "probe") {
      if (m.seq === S.probeSeq && !S.busy) {
        if (m.data.ok) {
          setEta("dim = " + m.data.dim +
                 (m.data.eta ? " · " + m.data.eta.label : ""), false);
        } else {
          if (m.data.detail) console.error(m.data.detail);
          setEta(m.data.error.type + ": " + m.data.error.message, true);
        }
      }
    } else if (m.type === "built") {
      if (m.data.ok) {
        el.results.appendChild(h("div", { "class": "qlgui-block",
          text: m.data.algebra }));
        if (m.eta) S.eta = m.eta;
      } else { renderError(m.data); }
    } else if (m.type === "result") {
      if (m.data.ok) renderBlock(m.data); else renderError(m.data);
      if (m.eta) S.eta = m.eta;
    } else if (m.type === "trace") {
      S.artifacts.traceHtml = m.html;
      S.artifacts.traceTex = m.tex || "";
      S.artifacts.traceJson = m.json || "";
      el.print.disabled = el["report-html"].disabled = !m.html;  // typeset MathML report
      el["worked-tex"].disabled = !S.artifacts.traceTex;   // .tex download (Plan 30 C1)
      el["report-json"].disabled = !S.artifacts.traceJson; // .json machine record (Plan 34)
    } else if (m.type === "artifacts") {
      S.artifacts.tikz = m.tikz; S.artifacts.snippet = m.snippet;
      S.artifacts.bundle = m.bundle;
      el.tikz.disabled = el.json.disabled = el.snippet.disabled = false;
    } else if (m.type === "done") {
      stopTicker();
      setBusy(false);
      setStatus("engine ready — quiverlab " + (S.manifest ? S.manifest.quiverlab_version : ""), "ok");
      scheduleProbe();
    } else if (m.type === "fatal") {
      console.error("qlgui engine:", m.message);
      stopTicker();
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

  // ---- module-block rendering helpers (Plan 26) ----
  function dvText(dv) {
    return "{" + Object.keys(dv).map(function (k) { return k + ": " + dv[k]; }).join(", ") + "}";
  }
  function degreeTable(rowLabel, dims) {   // n | value_n, like the HH table
    var head = h("tr"), row = h("tr");
    head.appendChild(h("th", { text: "n" }));
    row.appendChild(h("th", { text: rowLabel }));
    dims.forEach(function (d, n) {
      head.appendChild(h("td", { text: String(n) }));
      row.appendChild(h("td", { text: String(d) }));
    });
    return h("table", {}, head, row);
  }
  // Plan 34 (Marco): rad/top/soc are shown as FULL representations -- the dim
  // VECTOR per object (the redundant total-dim column is gone) PLUS each arrow's
  // exact action matrix, typeset like the other matrix blocks.

  // Matrices are the COMPLETE human record (Plan 34, Marco): shown IN FULL, wrapped in
  // a horizontally-scrollable container (mathScroll) so the page body never scrolls
  // sideways -- NOT elided at a small size. MAT_BACKSTOP_CELLS is only a SANITY cap
  // mirroring the trace recorder's record-time memory backstop: a pathological/corrupt
  // payload past it is stated by shape instead of hanging the browser. One constant per
  // file (app.js has the same one, same comment).
  var MAT_BACKSTOP_CELLS = 250000;
  function matTooBig(mat) {
    var rows = (mat || []).length, cols = rows ? (mat[0] || []).length : 0;
    return rows * cols > MAT_BACKSTOP_CELLS;
  }
  function matLatex(mat) {                  // [[..],[..]] -> \begin{pmatrix}..\end{pmatrix}
    mat = mat || [];
    if (matTooBig(mat)) {                   // sanity backstop only (never normal use)
      var c = mat.length ? (mat[0] || []).length : 0;
      return "\\text{[" + mat.length + "\\times" + c + " matrix beyond the display backstop]}";
    }
    var body = mat.map(function (row) {
      return row.map(String).join(" & ");
    }).join(" \\\\ ");
    return "\\begin{pmatrix} " + body + " \\end{pmatrix}";
  }
  // A matrix wrapped in a horizontally-scrollable inline box (Plan 34, Marco): the
  // full matrix scrolls INSIDE this box, so a wide differential never makes the page
  // body scroll sideways.
  function mathScroll(latex) {
    var box = h("span");
    box.style.display = "inline-block";
    box.style.maxWidth = "100%";
    box.style.overflowX = "auto";
    box.style.verticalAlign = "middle";
    box.appendChild(h("span", { "class": "arithmatex", text: "\\(" + latex + "\\)" }));
    return box;
  }
  // A pre-Plan-34 cached rad/top/soc lacked the per-view {dims, maps} (MINOR-6); an
  // extension-field module carries non-re-enterable entries (MAJOR-4). Pure predicates.
  function radTopSocStale(block) {
    return [block.radical, block.top, block.socle].some(function (v) {
      return !v || v.dims == null || v.maps == null;
    });
  }
  function radTopSocDisplayOnly(block) {
    return [block.radical, block.top, block.socle].some(function (v) {
      return v && v.display_only === true;
    });
  }
  function repDimTable(pairs) {             // label | dim vector (NO total-dim column)
    var head = h("tr");
    ["", "dim vector"].forEach(function (t) { head.appendChild(h("th", { text: t })); });
    var tbl = h("table", {}, head);
    pairs.forEach(function (p) {
      var r = h("tr");
      r.appendChild(h("th", { text: p[0] }));
      r.appendChild(h("td", { text: dvText(p[1].dims) }));
      tbl.appendChild(r);
    });
    return tbl;
  }
  function appendRepMaps(div, label, view) {   // one arrow-matrix line per arrow
    var maps = (view && view.maps) || {}, arrows = Object.keys(maps);
    if (!arrows.length) {
      div.appendChild(h("p", { "class": "qlgui-cites",
        text: label + ": every arrow acts as zero" }));
      return;
    }
    arrows.forEach(function (a) {
      div.appendChild(h("p", {}, document.createTextNode(label + ", arrow " + a + ": "),
        mathScroll(matLatex(maps[a]))));   // full matrix, horizontally-scrollable box
    });
  }
  // Plan 30 (Marco #3): term | ⊕-decomposition. The "# summands" and dim-vector
  // columns are gone from the RENDERING; the raw `betti`/`terms` fields stay in the
  // block JSON for backward-compat. `b.summands[n]` is the LaTeX P_1^2 ⊕ P_3.
  function resTable(b) {
    var head = h("tr");
    ["term", "⊕-decomposition"].forEach(function (t) {
      head.appendChild(h("th", { text: t }));
    });
    var tbl = h("table", {}, head);
    var summ = b.summands || [];
    (b.terms || []).forEach(function (dv, n) {
      var r = h("tr");
      r.appendChild(h("td", { text: String(n) }));
      var latex = (summ[n] != null) ? summ[n] : "0";
      // className EXACTLY "arithmatex" so the site MathJax config typesets it.
      r.appendChild(h("td", {}, h("span", { "class": "arithmatex",
        text: "\\(" + latex + "\\)" })));
      tbl.appendChild(r);
    });
    return tbl;
  }

  // Krull–Schmidt summand table: summand | multiplicity | dim vector (Plan 30).
  function decompTable(summands) {
    var head = h("tr");
    ["summand", "multiplicity", "dim vector"].forEach(function (t) {
      head.appendChild(h("th", { text: t }));
    });
    var tbl = h("table", {}, head);
    (summands || []).forEach(function (s, i) {
      var r = h("tr");
      r.appendChild(h("td", { text: "M_" + (i + 1) }));
      r.appendChild(h("td", { text: String(s.multiplicity) }));
      r.appendChild(h("td", { text: dvText(s.dim_vector) }));
      tbl.appendChild(r);
    });
    return tbl;
  }

  // The AR-translate input certificate (Marco #1): indecomposable, or M's
  // Krull–Schmidt decomposition + the additivity note. No-op when the block
  // carries no certificate (the decompose engine was unavailable).
  var TAU_NOTES = { "mod.tau_additive": "τ computed summand-wise (τ is additive)" };
  function appendInputCert(div, b) {
    if (b.indecomposable === true) {
      div.appendChild(h("p", { "class": "qlgui-hint", text: "input M is indecomposable" }));
    } else if (b.decomposition) {
      var parts = b.decomposition.map(function (s) {
        return dvText(s.dim_vector) + (s.multiplicity > 1 ? "^" + s.multiplicity : "");
      }).join("  ⊕  ");
      div.appendChild(h("p", { "class": "qlgui-hint",
        text: "input M ≅ " + parts + " — " + (TAU_NOTES[b.note_key] || "") }));
    }
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
    } else if (name === "dimension_vector" || name === "tau" || name === "tau_minus" ||
               name === "projective_dimension" || name === "injective_dimension") {
      div.appendChild(h("p", { "class": "arithmatex", text: "\\[ " + b.latex + " \\]" }));
      if (name === "tau" || name === "tau_minus") appendInputCert(div, b);
    } else if (name === "rad_top_soc") {
      div.appendChild(h("p", { text: "radical / top / socle:" }));
      if (radTopSocStale(b)) {              // MINOR-6: honest "recompute", never a fake zero
        div.appendChild(h("p", { "class": "qlgui-error",
          text: "this result was computed by an older version — recompute to see the "
                + "full representation." }));
      } else {
        if (radTopSocDisplayOnly(b)) {      // MAJOR-4: extension-field entries not re-enterable
          div.appendChild(h("p", { "class": "qlgui-hint",
            text: "display only — entries lie outside the integer/fraction input grammar "
                  + "(e.g. GF(p^n) elements) and cannot be re-entered in the module panel." }));
        }
        var trio = [["rad M", b.radical], ["top M", b.top], ["soc M", b.socle]];
        div.appendChild(repDimTable(trio));
        trio.forEach(function (p) { appendRepMaps(div, p[0], p[1]); });
      }
    } else if (name === "decompose") {
      div.appendChild(h("p", { text: "Krull–Schmidt decomposition — " + b.iso_classes +
        " indecomposable summand(s):" }));
      div.appendChild(decompTable(b.summands));
    } else if (name === "ext") {
      div.appendChild(h("p", { text: "Ext to the target module — dim vector " + dvText(b.target.dimvec) + ":" }));
      div.appendChild(degreeTable("dim Ext^n", b.dims));
    } else if (name === "tor") {
      div.appendChild(h("p", { text: "Tor with the target left module — dim vector " + dvText(b.target.dimvec) + ":" }));
      div.appendChild(degreeTable("dim Tor_n", b.dims));
    } else if (name === "projective_resolution" || name === "injective_resolution") {
      var proj = name === "projective_resolution";
      div.appendChild(h("p", { text: proj ? "projective resolution" : "injective resolution" }));
      div.appendChild(resTable(b));
      var d = proj ? b.pd : b.injective_dimension;
      div.appendChild(h("p", { text: (proj ? "pd = " : "id = ") +
        (d == null ? "∞ (beyond the probed length)" : String(d)) }));
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
    S.artifacts = { tikz: "", snippet: "", bundle: "", traceHtml: "", traceTex: "" };
    el.print.disabled = el.tikz.disabled = el.json.disabled = el.snippet.disabled = true;
    el["report-html"].disabled = el["worked-tex"].disabled = el["report-json"].disabled = true;
    S.eta = null;
    setBusy(true);
    setStatus("computing…");
    startTicker();
    S.worker.postMessage({ cmd: "run", request: buildRequest(), factor: S.factor });
  });
  el.cancel.addEventListener("click", function () {
    if (!S.busy) return;
    stopTicker();
    S.factor = null;                            // worker died: recalibrate on ready
    setBusy(false);
    startWorker();                             // sets its own transient status…
    setStatus("cancelled — engine reloading…"); // …so the acknowledgment must land last
    render();                                   // engineReady just went false — disable Compute
  });
  el.relations.addEventListener("input", scheduleProbe);
  [el.field, el.p, el.n, el.hhc, el["hhc-top"], el.hhh, el["hhh-top"], el.cartan,
   el.coxeter_polynomial, el.global_dimension, el.center]
    .forEach(function (x) { x.addEventListener("change", scheduleProbe); });
  // Module panel: enable/mode/side rebuild the dynamic body; the kind controls
  // just re-probe. The panel itself refreshes on every render() (vertex/arrow ops).
  el["mod-enable"].addEventListener("change", function () {
    renderModulePanel(); scheduleProbe();
  });
  [el["mod-mode"], el["mod-side"]].forEach(function (x) {
    x.addEventListener("change", function () { renderModulePanel(); scheduleProbe(); });
  });
  [el.dimension_vector, el.rad_top_soc, el.tau, el.tau_minus,
   el.projective_dimension, el.injective_dimension, el.decompose,
   el.projective_resolution, el["pr-top"], el.injective_resolution, el["ir-top"],
   el["ext-top"], el["tor-top"]]
    .forEach(function (x) { x.addEventListener("change", scheduleProbe); });
  // Ext/Tor toggles + the second-argument editor's mode/side rebuild the target
  // panel (show/hide, S/P/I vs explicit, side) then re-probe.
  [el.ext, el.tor, el["target-mode"], el["target-side"]].forEach(function (x) {
    x.addEventListener("change", function () { renderModulePanel(); scheduleProbe(); });
  });
  function download(name, text, type) {
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], { type: type }));
    a.download = name;
    a.click();
    URL.revokeObjectURL(a.href);
  }
  el.print.addEventListener("click", function () {
    // The report is self-contained typeset MathML (render_html) -- open it in a
    // new tab and print it; the browser's "Save as PDF" produces the PDF. MathML
    // typesets during layout (no async), so printing on load is safe.
    var url = URL.createObjectURL(new Blob([S.artifacts.traceHtml], { type: "text/html" }));
    var w = window.open(url, "_blank");
    if (w) {
      w.addEventListener("load", function () { w.print(); });
    } else {
      // MINOR-7: window.open returns null when the browser blocks the popup -- say so
      // visibly (the existing notice channel) instead of silently doing nothing.
      setStatus("popup blocked — allow popups to print, or use “Report (HTML)” to download it", "err");
    }
  });
  el["report-html"].addEventListener("click", function () {  // print-ready typeset report
    download("report.html", S.artifacts.traceHtml, "text/html");
  });
  el["worked-tex"].addEventListener("click", function () {   // Plan 30 C1: the .tex source
    download("worked-steps.tex", S.artifacts.traceTex, "text/x-tex");
  });
  el["report-json"].addEventListener("click", function () {  // Plan 34: the JSON machine record
    download("trace.json", S.artifacts.traceJson, "application/json");
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
  el.config.addEventListener("click", function () {
    // The exported request as a runnable cluster config: a one-line header naming
    // the command, then the YAML request body.
    var text = "# run on a cluster: quiverlab-hpc run this-file.yaml -o result.json\n" +
               configYaml(buildRequest());
    download("quiverlab-config.yaml", text, "text/yaml");
  });

  window.QLGUI = { S: S, buildRequest: buildRequest, configYaml: configYaml };
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
  el["mod-enable"].addEventListener("change", ensureEngine);   // enabling the module panel is intent
})();
