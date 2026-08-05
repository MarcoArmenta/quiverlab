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
            artifacts: { tikz: "", snippet: "", bundle: "", traceHtml: "", traceJson: "" },
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
    // ---- Plan 35: HH product surface (cup / cap / bracket / connes_b) ----
    '  <label><input type="checkbox" id="qlgui-cup"> cup 0..<input type="number" id="qlgui-cup-top" value="2" min="0"></label>' +
    '  <label><input type="checkbox" id="qlgui-cap"> cap 0..<input type="number" id="qlgui-cap-top" value="2" min="0"></label>' +
    '  <label><input type="checkbox" id="qlgui-bracket"> bracket 0..<input type="number" id="qlgui-bracket-top" value="2" min="0"></label>' +
    '  <label><input type="checkbox" id="qlgui-connes_b"> Connes B 0..<input type="number" id="qlgui-connes_b-top" value="2" min="0"></label>' +
    // Plan-35 follow-up: cyclic homology HC_0..HC_n (default 6), right after Connes B.
    '  <label><input type="checkbox" id="qlgui-cyclic_homology"> cyclic homology 0..<input type="number" id="qlgui-cyclic_homology-top" value="6" min="0"></label>' +
    // Plan-42: the Hochschild (b, B) spectral sequence (abutting to HC), right after cyclic homology.
    '  <label><input type="checkbox" id="qlgui-ss_hochschild"> (b,B) spectral sequence 0..<input type="number" id="qlgui-ss_hochschild-top" value="4" min="0"></label>' +
    '  <label><input type="checkbox" id="qlgui-cartan" checked> Cartan matrix</label>' +
    '  <label><input type="checkbox" id="qlgui-coxeter_polynomial"> Coxeter polynomial</label>' +
    '  <label><input type="checkbox" id="qlgui-global_dimension"> gl.dim</label>' +
    '  <label><input type="checkbox" id="qlgui-center"> center</label>' +
    // ---- Plan 38: Yoneda Ext-algebra + Koszulity, and the recognizer batch ----
    '  <label><input type="checkbox" id="qlgui-ext_algebra"> Ext-algebra / Koszul 0..<input type="number" id="qlgui-ext_algebra-top" value="6" min="0"></label>' +
    '  <label><input type="checkbox" id="qlgui-recognizers"> recognizers + type</label>' +
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
    '  <button id="qlgui-report-json" class="qlgui-secondary" type="button" disabled title="Download the complete worked-steps event stream (exact, machine-readable JSON)">Report data (JSON)</button>' +
    '  <button id="qlgui-tikz" class="qlgui-secondary" type="button" disabled>TikZ</button>' +
    '  <button id="qlgui-json" class="qlgui-secondary" type="button" disabled>JSON</button>' +
    '  <button id="qlgui-snippet" class="qlgui-secondary" type="button" disabled>Copy Python</button>' +
    '  <button id="qlgui-config" class="qlgui-secondary" type="button" disabled>Config (YAML)</button>' +
    '</div>' +
    '<div id="qlgui-results"></div>';

  var el = {};
  ["preset", "field", "p-wrap", "n-wrap", "p", "n", "clear", "status", "canvas",
   "rename", "relations", "hhc", "hhc-top", "hhh", "hhh-top",
   // Plan 35 HH product surface: cup / cap / bracket / connes_b + degree pickers
   "cup", "cup-top", "cap", "cap-top", "bracket", "bracket-top",
   "connes_b", "connes_b-top", "cyclic_homology", "cyclic_homology-top",
   // Plan 42: the Hochschild (b, B) spectral sequence + its degree picker
   "ss_hochschild", "ss_hochschild-top", "cartan",
   "coxeter_polynomial", "global_dimension", "center",
   // Plan 38: Ext-algebra/Koszul (with a degree picker) + the recognizer batch
   "ext_algebra", "ext_algebra-top", "recognizers",
   "trace", "compute",
   "cancel", "print", "report-html", "report-json", "tikz", "json", "snippet", "config", "results", "eta",
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
    // Plan 35 HH product surface, in the Task-12 curated-request order (cup, cap,
    // bracket, connes_b) immediately after hh_homology -- do not reorder. Cyclic
    // homology (HC) follows connes_b (Plan-35 follow-up).
    if (el.cup.checked) compute.push("cup:0.." + el["cup-top"].value);
    if (el.cap.checked) compute.push("cap:0.." + el["cap-top"].value);
    if (el.bracket.checked) compute.push("bracket:0.." + el["bracket-top"].value);
    if (el.connes_b.checked) compute.push("connes_b:0.." + el["connes_b-top"].value);
    if (el.cyclic_homology.checked)
      compute.push("cyclic_homology:0.." + el["cyclic_homology-top"].value);
    // Plan 42: the (b, B) spectral sequence, right after cyclic homology.
    if (el.ss_hochschild.checked)
      compute.push("ss_hochschild:0.." + el["ss_hochschild-top"].value);
    if (el.ext_algebra.checked)
      compute.push("ext_algebra:0.." + el["ext_algebra-top"].value);
    ["cartan", "coxeter_polynomial", "global_dimension", "center",
     "recognizers"].forEach(function (k) {
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
    // Schema tag is HONEST: module / Ext / Tor blocks are the schema-2 surface
    // (the webapp validator refuses them under schema 1 -- found live when the
    // ported canvas 422'd on every module request); plain algebra requests keep
    // the schema-1 tag so their cache keys stay byte-stable.
    var req = { schema: (module || extTarget || torTarget) ? 2 : 1,
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
      S.artifacts.traceJson = m.json || "";
      el.print.disabled = el["report-html"].disabled = !m.html;  // typeset MathML report
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

  // Matrices are the COMPLETE human record (Plan 34, Marco): shown IN FULL and never
  // elided at a small size. MAT_BACKSTOP_CELLS is only a SANITY cap mirroring the
  // trace recorder's record-time memory backstop: a pathological/corrupt payload past
  // it is stated by shape instead of hanging the browser. One constant per file
  // (app.js has the same one, same comment).
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
  // A matrix shown COMPLETE, with NO scrollbar (Marco, 2026-07-29): the box carries
  // the class the post-typeset `fitMath` pass measures, and a matrix wider than the
  // column is SHRUNK to fit (transform: scale) rather than clipped behind a scroll
  // bar. Nothing is ever hidden -- worst case the glyphs get smaller.
  function mathFit(latex) {
    var box = h("span", { "class": "qlgui-fit" });
    box.appendChild(h("span", { "class": "arithmatex", text: "\\(" + latex + "\\)" }));
    return box;
  }
  // A matrix as an INDEXED GRID (Marco 2026-07-29): a header row of column indices,
  // a header column of row indices, and a light rule between cells, so an entry can
  // be read off by position. 1-based, the mathematician's convention. Entries are
  // shown verbatim (exact ints / fraction strings), never reformatted.
  function matrixGrid(mat) {
    mat = mat || [];
    var ncols = mat.length ? (mat[0] || []).length : 0;
    // Marco 2026-08-03: a zero MAP is stated ("0"), never drawn as a grid of 0s
    // (mirrors quiverlab.trace.render_html.matrix_grid).
    if (!mat.length || !ncols || matIsZero(mat)) {
      return h("p", { "class": "arithmatex", text: "\\( 0 \\)" });
    }
    // Marco 2026-08-02: only matrices with fewer than 20 rows AND columns are shown;
    // a larger one states its size and points at the JSON (the complete matrix is
    // always in the accompanying record). Single chokepoint for every ORDINARY grid.
    // Product Cayley tables do NOT route through here (cayleyBigGrid) -- uncapped.
    if (mat.length >= 20 || ncols >= 20) {
      return h("p", { "class": "qlgui-cites", text: mat.length + "×" + ncols
        + " matrix (exceeds the 20-line display cap); the complete matrix is in the "
        + "accompanying JSON record." });
    }
    if (matTooBig(mat)) {                   // sanity backstop only (never normal use)
      return h("p", { "class": "qlgui-cites",
        text: mat.length + "\u00d7" + ncols + " matrix beyond the display backstop" });
    }
    var head = h("tr");
    head.appendChild(h("th", { "class": "qlgui-corner", text: "" }));
    for (var j = 0; j < ncols; j++) head.appendChild(h("th", { text: String(j + 1) }));
    var tbl = h("table", { "class": "qlgui-matrix" }, head);
    mat.forEach(function (row, i) {
      var r = h("tr");
      r.appendChild(h("th", { text: String(i + 1) }));
      (row || []).forEach(function (x) { r.appendChild(h("td", { text: String(x) })); });
      tbl.appendChild(r);
    });
    return tbl;
  }
  // Shrink-to-fit every .qlgui-fit box whose typeset content overflows its column.
  // Called after MathJax/KaTeX finishes; scale is shrink-ONLY (never magnifies), and
  // the wrapper's height is corrected so the scaled box does not overlap its
  // neighbours. Idempotent: the measurement reads the UNSCALED width each time.
  function fitMath() {
    var boxes = el.results.querySelectorAll(".qlgui-fit");
    for (var i = 0; i < boxes.length; i++) {
      var box = boxes[i], inner = box.firstChild;
      if (!inner) continue;
      inner.style.transform = "";
      inner.style.display = "inline-block";
      inner.style.transformOrigin = "left top";
      var avail = box.parentNode ? box.parentNode.clientWidth : 0;
      var want = inner.scrollWidth;
      if (!avail || !want || want <= avail) { box.style.height = ""; continue; }
      var k = avail / want;
      inner.style.transform = "scale(" + k + ")";
      box.style.height = Math.ceil(inner.offsetHeight * k) + "px";
    }
  }
  // A matrix is EXACTLY zero -- its arrow carries no information, so the display
  // omits it and names it in one line instead (Marco, 2026-07-29).
  // Marco 2026-08-03: an engine provenance line says what the engine IS (mirrors
  // quiverlab.trace.results_html._ENGINE_GLOSS).
  function engineNote(engine) {
    var s = String(engine || "");
    var low = s.toLowerCase();
    if (low.indexOf("hanlab") !== -1)
      return s + " \u2014 quiverlab's exact GF(p) linear-algebra engine: it assembles "
        + "the boundary/coboundary matrices of the chosen (co)chain complex with "
        + "integer entries mod p and computes their exact rank by Gaussian "
        + "elimination mod p; every dimension follows by rank-nullity \u2014 nothing "
        + "numerical, no floating point.";
    if (low.indexOf("chouhy") !== -1 || low.indexOf("solotar") !== -1)
      return s + " \u2014 the Chouhy\u2013Solotar projective bimodule resolution built "
        + "from the admissible presentation, certified per instance "
        + "(d\u2218d = 0 + the order gate).";
    if (low.indexOf("(b,b)") !== -1)
      return s + " \u2014 the exact (b, B) mixed-complex engine on the normalized bar "
        + "complex: b is the Hochschild boundary, B the Connes boundary.";
    return s;
  }
  function matIsZero(mat) {
    return (mat || []).every(function (row) {
      return (row || []).every(function (x) { return String(x) === "0"; });
    });
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
  function loewyFactors(layer) {           // a Loewy layer as S_v ⊕ S_w^m (Plan 37)
    var parts = [];
    Object.keys(layer).sort().forEach(function (v) {
      var m = layer[v];
      if (m) parts.push(m === 1 ? "S_" + v : "S_" + v + "^" + m);
    });
    return parts.length ? parts.join(" ⊕ ") : "0";
  }
  function loewySeriesTable(series) {       // layer | factors, top to bottom
    var head = h("tr");
    ["layer", "factors"].forEach(function (t) { head.appendChild(h("th", { text: t })); });
    var tbl = h("table", {}, head);
    series.forEach(function (layer, i) {
      var r = h("tr");
      r.appendChild(h("td", { text: String(i + 1) }));
      r.appendChild(h("td", { text: loewyFactors(layer) }));
      tbl.appendChild(r);
    });
    return tbl;
  }
  // One arrow-matrix line per arrow that acts NON-trivially. An arrow acting as the
  // zero map carries no information, so its matrix is not printed (Marco,
  // 2026-07-29) -- the arrows are named in a single trailing line so the reader can
  // still tell "acts as zero" apart from "not an arrow of the quiver".
  function appendRepMaps(div, label, view) {
    var maps = (view && view.maps) || {}, arrows = Object.keys(maps);
    var zero = arrows.filter(function (a) { return matIsZero(maps[a]); });
    var live = arrows.filter(function (a) { return !matIsZero(maps[a]); });
    if (!arrows.length || !live.length) {
      div.appendChild(h("p", { "class": "qlgui-cites",
        text: label + ": every arrow acts as zero" }));
      return;
    }
    live.forEach(function (a) {
      div.appendChild(h("p", { text: label + ", arrow " + a + ":" }));
      div.appendChild(matrixGrid(maps[a]));       // indexed grid, shown complete
    });
    if (zero.length) {
      div.appendChild(h("p", { "class": "qlgui-cites",
        text: label + ": arrow" + (zero.length > 1 ? "s " : " ") + zero.join(", ")
              + " act" + (zero.length > 1 ? "" : "s") + " as zero" }));
    }
  }
  // pd/id display for a block that predates the runners' `latex` key (an older
  // cached result). Mirrors quiverlab.hpc.spec._homdim_latex: an unresolved probe
  // states the certified lower bound, never a bare ∞ it did not prove.
  function homdimLatex(name, b) {
    var op = name === "projective_dimension" ? "pd" : "id";
    if (b.value == null) {
      return "\\operatorname{" + op + "} M > " + (b.bound != null ? b.bound : 32);
    }
    return "\\operatorname{" + op + "} M = " + b.value;
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
  // A summand recognised as a STANDARD indecomposable is NAMED S_v / P_v / I_v and
  // needs no matrices; any other one is shown in full below the table, since its
  // dimension vector does not determine it (Marco 2026-07-29).
  var STD_SYM = { simple: "S", projective: "P", injective: "I" };
  function summandName(s, i) {
    var std = s && s.standard;
    if (std && STD_SYM[std.kind]) return STD_SYM[std.kind] + "_" + std.vertex;
    return "M_" + i;
  }
  function decompTable(summands) {
    var head = h("tr");
    ["summand", "multiplicity", "dim vector"].forEach(function (t) {
      head.appendChild(h("th", { text: t }));
    });
    var tbl = h("table", {}, head);
    (summands || []).forEach(function (s, i) {
      var r = h("tr");
      r.appendChild(h("td", { text: summandName(s, i + 1) }));
      r.appendChild(h("td", { text: String(s.multiplicity) }));
      r.appendChild(h("td", { text: dvText(s.dim_vector) }));
      tbl.appendChild(r);
    });
    return tbl;
  }
  function appendSummandMaps(div, summands) {
    (summands || []).forEach(function (s, i) {
      if (s.standard || !s.maps) return;
      appendRepMaps(div, summandName(s, i + 1), s);
    });
  }

  // The differentials of a projective/injective resolution, as full matrices
  // (Marco, 2026-07-29). A differential EQUAL to one already shown is not repeated:
  // the line references the earlier degree instead, which is the whole point in a
  // periodic resolution where every second map is the same. `d.elided` (the
  // recorder's memory backstop) states the shape honestly.
  function appendDifferentials(div, b, proj) {
    var diffs = b.differentials;
    if (!diffs || !diffs.length) return;
    div.appendChild(h("p", { text: proj
      ? "differentials (rows: target basis, columns: source basis; d_0 = ε: Q_0 → M)"
      : "differentials (rows: target basis, columns: source basis; d^0 = ι: M → E^0)" }));
    var seen = [];                     // [{key, label}] in first-seen order
    diffs.forEach(function (d, n) {
      var sym = proj ? ("d_{" + n + "}") : ("d^{" + n + "}");
      var label = (proj ? "d_" : "d^") + n;
      if (d.elided || !d.matrix) {
        div.appendChild(h("p", { text: label + ": " + d.rows + "×" + d.cols
          + " matrix (too large to display; it is complete in the report data)" }));
        return;
      }
      if (matIsZero(d.matrix)) {
        // Marco 2026-08-03: a zero map is stated, never drawn or cross-referenced.
        div.appendChild(h("p", { "class": "arithmatex", text: "\\(" + sym + " = 0\\)" }));
        return;
      }
      var key = JSON.stringify(d.matrix);
      var prior = null;
      for (var i = 0; i < seen.length; i++) {
        if (seen[i].key === key) { prior = seen[i].label; break; }
      }
      if (prior !== null) {
        div.appendChild(h("p", { "class": "qlgui-hint",
          text: label + " = " + prior + " (same matrix; not repeated)" }));
        return;
      }
      seen.push({ key: key, label: label });
      div.appendChild(h("p", { "class": "arithmatex", text: "\\(" + sym + " =\\)" }));
      div.appendChild(matrixGrid(d.matrix));
    });
  }

  // The ordered k-basis of each resolution term (Plan 35 UNIT 2), as a numbered
  // (1-based) list per term -- the SAME index order the differential grids use for
  // their columns (projective) / rows (injective). Tolerant of a block WITHOUT
  // `term_basis` (an older cached result / a structure-constants algebra): renders
  // nothing. Mirrors quiverlab.trace.results_html._term_basis_html.
  var TERM_BASIS_DISPLAY = 20;
  function appendTermBasis(div, b, proj) {
    var tb = b.term_basis;
    if (!tb || !tb.length) return;
    div.appendChild(h("p", { text: "Ordered basis of each resolution term — the "
      + "1-based index order the differential grids below use for their "
      + (proj ? "columns." : "rows.") }));
    tb.forEach(function (labels, n) {
      var name = proj ? ("Q_" + n) : ("E^" + n);
      if (!labels || !labels.length) {
        div.appendChild(h("p", { text: name + " = 0" }));
        return;
      }
      div.appendChild(h("p", { text: name + " basis:" }));
      var ol = h("ol");
      labels.slice(0, TERM_BASIS_DISPLAY).forEach(function (x) {
        ol.appendChild(h("li", { text: String(x) }));
      });
      div.appendChild(ol);
      if (labels.length > TERM_BASIS_DISPLAY) {
        div.appendChild(h("p", { "class": "qlgui-cites", text: "… "
          + (labels.length - TERM_BASIS_DISPLAY) + " more (full list in the report data)" }));
      }
    });
  }

  // The AR-translate input certificate (Marco #1): indecomposable, or the input's
  // Krull–Schmidt decomposition + the additivity note. No-op when the block
  // carries no certificate (the decompose engine was unavailable).
  var TAU_NOTES = { "mod.tau_additive": "τ computed summand-wise (τ is additive)" };
  function appendInputCert(div, b, name) {
    name = name || "M";
    if (b.indecomposable === true) {
      div.appendChild(h("p", { "class": "qlgui-hint",
        text: "input " + name + " is indecomposable" }));
    } else if (b.decomposition) {
      var parts = b.decomposition.map(function (s) {
        return dvText(s.dim_vector) + (s.multiplicity > 1 ? "^" + s.multiplicity : "");
      }).join("  ⊕  ");
      div.appendChild(h("p", { "class": "qlgui-hint",
        text: "input " + name + " ≅ " + parts + " — " + (TAU_NOTES[b.note_key] || "") }));
    }
  }

  // One AR translate: its dimension-vector line, the FULL per-arrow action matrices
  // of the translate itself, and the input's indecomposability certificate. Used for
  // M and, when the request names a second module, for N too (Marco, 2026-07-29).
  var TARGET_ROLE = { ext_target: "the Ext target", tor_target: "the Tor target" };
  function appendTranslate(div, t, kind, name) {
    var sym = (kind === "tau" ? "τ" : "τ⁻") + name;
    if (t.error) {
      div.appendChild(h("p", { "class": "qlgui-error",
        text: sym + " is unavailable: " + t.error }));
      return;
    }
    div.appendChild(h("p", { "class": "arithmatex", text: "\\[ " + t.latex + " \\]" }));
    if (t.repr) appendRepMaps(div, sym, t.repr);
    appendInputCert(div, t, name);
  }

  // ---- HH product surface rendering (Plan 35) ----
  // .blocks() shapes: cup/cap/bracket carry `tables`, each with degrees /
  // out_degree / dims=[dl,dr,dout] / constants[k][i][j] (exact strings);
  // connes_b carries per-n `matrices` + `ranks`. gui.js hardcodes English (the
  // webapp families page carries the i18n twins block.*.title / products.zero).
  var PRODUCT_TITLE = { cup: "Cup product tables", cap: "Cap product tables",
                        bracket: "Gerstenhaber bracket tables",
                        connes_b: "Connes differentials" };
  var PRODUCT_OP = { cup: "\\cup", cap: "\\cap" };
  function coeffTerm(c, sym) {              // exact-string coeff -> "c \, sym" LaTeX
    return c === "1" ? sym : c + " \\, " + sym;   // matches trace.products._term
  }
  // The nonzero structure-constant equations of ONE table, in the report's SHARED
  // notation (quiverlab.trace.products.equation_lines): alpha^p_i / beta^q_j (z^n_j
  // for cap) on the left, gamma^{p+q}_k (w^{n-p}_k for cap) on the right; zero terms
  // skipped and a fully-zero equation omitted (the whole-table-zero case is the
  // caller's "vanish" line). One notation everywhere -- report, results, GUI.
  function productEquations(name, t) {
    var K = t.constants || [], dims = t.dims || [0, 0, 0];
    var degs = t.degrees || [0, 0], out = t.out_degree;
    var dl = dims[0], dr = dims[1], dout = dims[2], lines = [];
    var p = degs[0], q = degs[1];
    // Marco 2026-07-31: UNIFORM notation -- every cohomology class is α^deg_index and
    // every homology class z^deg_index; β/γ/w are gone. cup/bracket act on cohomology
    // (all α); cap outputs homology (z).
    var rightSym = name === "cap" ? "z" : "\\alpha";
    var outSym = name === "cap" ? "z" : "\\alpha";
    for (var i = 0; i < dl; i++) {
      for (var j = 0; j < dr; j++) {
        var terms = [];
        for (var k = 0; k < dout; k++) {
          var c = String(((K[k] || [])[i] || [])[j]);
          if (c === "0" || c === "undefined") continue;
          terms.push(coeffTerm(c, outSym + "^{" + out + "}_{" + (k + 1) + "}"));
        }
        if (!terms.length) continue;
        var L = "\\alpha^{" + p + "}_{" + (i + 1) + "}";
        var R = rightSym + "^{" + q + "}_{" + (j + 1) + "}";
        var lhs = (name === "bracket") ? "[" + L + ", " + R + "]"
          : L + " " + PRODUCT_OP[name] + " " + R;
        lines.push("\\( " + lhs + " = " + terms.join(" + ") + " \\)");
      }
    }
    return lines;
  }
  // Plan-35 follow-up (Marco): the notation legend defining the product-table
  // symbols (alpha/beta/gamma/z/w), shown above each family's tables. Hardcoded
  // English like PRODUCT_TITLE; the report surfaces build the same wording in
  // quiverlab.trace.products.notation_legend. The concrete recorded basis (b.basis,
  // e.g. "bar/GF(2)") is named so the reader knows what the constants refer to.
  function productLegend(name, b) {
    var onBasis = b.basis ? "relative to the recorded basis " + b.basis
                          : "relative to the recorded class basis";
    // Plan-35 UNIT 2: the legend points at the explicit per-degree listings, where
    // every α/β/γ/z/w is printed as its (co)cycle term-sum + coordinate vector with
    // the annihilating differential (mirrors trace.products.notation_legend).
    var explicit = " Each class is listed explicitly by degree below as a combination "
      + "of the ordered basis elements; the coordinate vectors are recorded in the JSON.";
    if (name === "cup")
      return "α^n_j denotes the j-th basis class of HH^n (superscript = degree, "
        + "subscript = index). Every table line states "
        + "α^p_i ∪ α^q_j = Σ_k c·α^{p+q}_k, " + onBasis + "; the constants c are "
        + "basis-dependent." + explicit;
    if (name === "bracket")
      return "α^n_j denotes the j-th basis class of HH^n (superscript = degree, "
        + "subscript = index). Every table line states "
        + "[α^p_i, α^q_j] = Σ_k c·α^{p+q-1}_k in degree p+q−1, " + onBasis
        + "; the constants c are basis-dependent." + explicit;
    return "α^p_j denotes the j-th basis class of HH^p (cohomology) and z^n_j the j-th "
      + "of HH_n (homology). Every table line states α^p_i ∩ z^n_j = Σ_k c·z^{n-p}_k, "
      + onBasis + "; the constants c are basis-dependent." + explicit;
  }
  // ---- Cayley multiplication tables (Marco 2026-08-01) ----
  // Every product bidegree renders as a GRID (not equation lines): rows = left classes,
  // cols = right classes, each cell the product written DIRECTLY in the target basis
  // (0 / a signed combination). Zeros are SHOWN inside a table that is nonzero
  // somewhere. Mirrors quiverlab.trace.products.cayley_table byte-for-byte in logic.
  var PRODUCT_CORNER = { cup: "\\cup", cap: "\\cap", bracket: "[-,-]" };
  function primeFromBasis(basis) {           // "bar/GF(7)" -> 7, "cs/QQ" -> null
    var m = /GF\((\d+)\)/.exec(String(basis || ""));
    return m ? parseInt(m[1], 10) : null;
  }
  function balancedCoeff(c, prime) {         // residue c > p/2 shown as c-p (display only)
    c = String(c);
    if (prime == null) return c;
    if (!/^\d+$/.test(c)) return c;          // fraction / non-residue -> verbatim
    var v = parseInt(c, 10);
    return (v >= 0 && v < prime && 2 * v > prime) ? String(v - prime) : c;
  }
  function balancedRepNote(prime) {
    return "Coefficients are shown as balanced representatives mod " + prime
      + " (a residue c > " + prime + "/2 is written c-" + prime + "); the JSON record "
      + "keeps the raw residues.";
  }
  function signedJoinTex(pieces) {           // (coeff, gen) terms -> signed TeX sum
    if (!pieces.length) return "0";          // spacing matches trace.products._signed_join
    return pieces.map(function (p, i) {
      var neg = p.c.charAt(0) === "-", mag = neg ? p.c.slice(1) : p.c;
      var term = (mag === "1") ? p.g : (mag + "\\," + p.g);
      if (i === 0) return neg ? ("-" + term) : term;
      return (neg ? " - " : " + ") + term;
    }).join("");
  }
  function cellTex(name, out, coeffs, prime) {   // one cell: Σ_k c_k·g_k in the target basis
    var outSym = name === "cap" ? "z" : "\\alpha", pieces = [];
    for (var k = 0; k < coeffs.length; k++) {
      var disp = balancedCoeff(coeffs[k], prime);
      if (String(disp) === "0") continue;
      pieces.push({ c: disp, g: outSym + "^{" + out + "}_{" + (k + 1) + "}" });
    }
    return signedJoinTex(pieces);
  }
  function mirrorSign(name, p, q) {          // graded transpose sign (+1 / -1)
    if (name === "cup") return ((p * q) % 2) ? -1 : 1;
    return (((p - 1) * (q - 1)) % 2 === 0) ? -1 : 1;   // bracket
  }
  function isIntStr(s) { return /^-?\d+$/.test(String(s)); }
  // Honest structural notes DERIVED from the constants (cup/bracket, square bidegree).
  function cayleyStructuralNotes(name, degrees, dims, constants, prime) {
    if (name !== "cup" && name !== "bracket") return [];
    var p = degrees[0], q = degrees[1], dl = dims[0], dr = dims[1], dout = dims[2];
    if (p !== q || dl !== dr || !dl) return [];
    var n = dl, K = constants || [], notes = [], i, j, k;
    var squares = true;
    for (i = 0; i < n && squares; i++)
      for (k = 0; k < dout; k++)
        if (String(((K[k] || [])[i] || [])[i]) !== "0") { squares = false; break; }
    if (squares) notes.push("all squares are 0");
    if (prime != null) {
      var allInt = true, mirrored = true, sign = mirrorSign(name, p, q);
      for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < dout; k++) {
        var a = String(((K[k] || [])[i] || [])[j]), b = String(((K[k] || [])[j] || [])[i]);
        if (!isIntStr(a) || !isIntStr(b)) { allInt = false; }
        else if (((parseInt(a, 10) - sign * parseInt(b, 10)) % prime + prime) % prime !== 0)
          mirrored = false;
      }
      if (allInt && mirrored)
        notes.push(sign === -1 ? "the table is graded-antisymmetric"
                               : "the table is graded-commutative (symmetric)");
    }
    return notes;
  }
  function cayleyNoteLine(name, degrees, dims, constants, prime) {
    var notes = cayleyStructuralNotes(name, degrees, dims, constants, prime);
    if (!notes.length) return "";
    var s = notes.join("; ");
    return s.charAt(0).toUpperCase() + s.slice(1) + ".";
  }
  // A Cayley grid as an indexed table (reuses qlgui-matrix: double zebra + corner):
  // corner = the product operator, header row = right classes, header column = left
  // classes, each cell the product math in the target basis.
  function mathCell(tex) {
    return h("span", { "class": "arithmatex", text: "\\(" + tex + "\\)" });
  }
  // ---- ONE big degree-major Cayley table per family (Marco 2026-08-01 addendum) ----
  // Rows/columns run over ALL (co)homology classes, degree-major; a cell whose target
  // degree is beyond the computed window is an em dash (not computed), a computed
  // vanishing product is 0. Mirrors quiverlab.trace.products.combined_cayley.
  // Marco 2026-08-02: the product Cayley table is UNCAPPED -- the one big degree-graded
  // table ALWAYS renders (product tables can be big), so there is no per-axis cap and no
  // per-bidegree fallback (mirrors quiverlab.trace.products.combined_cayley).
  var EM_DASH = "—";
  var FAMILY_HEADING = { cup: "HH^{*} \\cup HH^{*} \\to HH^{*}",
    cap: "HH^{*} \\cap HH_{*} \\to HH_{*}", bracket: "[HH^{*}, HH^{*}] \\to HH^{*}" };
  var FAMILY_AXIS_NOTE = "One row per left class and one column per right class, "
    + "ordered degree-major (the degree is the class' superscript); a heavier rule "
    + "marks each degree boundary.";
  function beyondWindowNote() {
    return EM_DASH + " marks a cell whose target degree lies beyond the computed "
      + "window (not computed); a computed vanishing product is shown as 0.";
  }
  function combinedOutDegree(name, p, q) {
    if (name === "cup") return p + q;
    if (name === "bracket") return p + q - 1;
    return q - p;                            // cap: (p, n) -> n - p
  }
  function combinedNote(name, tbl, rowMeta, prime) {
    if (name !== "cup" && name !== "bracket") return "";
    var notes = [], diagSeen = 0, diagZero = 0;
    rowMeta.forEach(function (rm) {
      var p = rm[0], ii = rm[1], t = tbl[p + "," + p];
      if (!t) return;
      diagSeen++;
      var dout = t.dims[2], z = true;
      for (var k = 0; k < dout; k++)
        if (String(((t.constants[k] || [])[ii] || [])[ii]) !== "0") z = false;
      if (z) diagZero++;
    });
    if (diagSeen && diagZero === diagSeen) notes.push("all squares are 0");
    if (prime != null) {
      var seen = 0, ok = 0, allInt = true;
      rowMeta.forEach(function (rm) {
        var p = rm[0], ii = rm[1];
        rowMeta.forEach(function (cm) {
          var q = cm[0], jj = cm[1], t = tbl[p + "," + q], tT = tbl[q + "," + p];
          if (!t || !tT || !allInt) return;
          var dout = t.dims[2], sign = mirrorSign(name, p, q), good = true;
          for (var k = 0; k < dout; k++) {
            var av = ((t.constants[k] || [])[ii] || [])[jj];
            var bv = ((tT.constants[k] || [])[jj] || [])[ii];
            if (!isIntStr(av) || !isIntStr(bv)) { allInt = false; return; }
            if ((((parseInt(av, 10) - sign * parseInt(bv, 10)) % prime) + prime)
                % prime !== 0) good = false;
          }
          seen++; if (good) ok++;
        });
      });
      if (allInt && seen && ok === seen)
        notes.push(name === "cup" ? "the cup product is graded-commutative"
                                  : "the Gerstenhaber bracket is graded-antisymmetric");
    }
    if (!notes.length) return "";
    var s = notes.join("; ");
    return s.charAt(0).toUpperCase() + s.slice(1) + ".";
  }
  function combinedCayley(name, tables, prime) {
    var tbl = {}, leftDims = {}, rightDims = {};
    (tables || []).forEach(function (t) {
      var d = t.degrees || [0, 0], dims = t.dims || [0, 0, 0];
      tbl[d[0] + "," + d[1]] = { dims: dims, constants: t.constants };
      leftDims[d[0]] = dims[0]; rightDims[d[1]] = dims[1];
    });
    var num = function (a, b) { return a - b; };
    var leftDegs = Object.keys(leftDims).map(Number).sort(num);
    var rightDegs = Object.keys(rightDims).map(Number).sort(num);
    var totalRows = leftDegs.reduce(function (s, p) { return s + leftDims[p]; }, 0);
    var totalCols = rightDegs.reduce(function (s, q) { return s + rightDims[q]; }, 0);
    // The block's computed top = max recorded out_degree; classify cells against it
    // EXPLICITLY (target > top -> beyond window / em dash), not by a missing key (an
    // in-window bidegree that is absent means the recording broke -> throw). Mirrors
    // quiverlab.trace.products.combined_cayley (Marco 2026-08-02).
    var top = 0;
    (tables || []).forEach(function (t) { if (t.out_degree > top) top = t.out_degree; });
    var rowLabels = [], rowDegsep = [], rowMeta = [];
    leftDegs.forEach(function (p, bi) {
      for (var i = 0; i < leftDims[p]; i++) {
        rowLabels.push("\\alpha^{" + p + "}_{" + (i + 1) + "}");
        rowDegsep.push(i === 0 && bi > 0); rowMeta.push([p, i]);
      }
    });
    var rightSym = name === "cap" ? "z" : "\\alpha";
    var colLabels = [], colDegsep = [], colMeta = [];
    rightDegs.forEach(function (q, bj) {
      for (var j = 0; j < rightDims[q]; j++) {
        colLabels.push(rightSym + "^{" + q + "}_{" + (j + 1) + "}");
        colDegsep.push(j === 0 && bj > 0); colMeta.push([q, j]);
      }
    });
    var cells = [], hasBeyond = false;
    rowMeta.forEach(function (rm) {
      var p = rm[0], i = rm[1], row = [];
      colMeta.forEach(function (cm) {
        var q = cm[0], j = cm[1], target = combinedOutDegree(name, p, q);
        if (target < 0) { row.push("0"); return; }   // cap below degree 0: structural zero
        if (target > top) {                           // beyond the computed window
          if (name === "cap")
            throw new Error("combined Cayley cap table: cell (" + p + ", " + q
              + ") target degree " + target + " > top " + top
              + ", impossible for a cap (n-p <= n <= top)");
          hasBeyond = true; row.push(EM_DASH); return;
        }
        var t = tbl[p + "," + q];
        if (!t)                                        // in-window but not recorded: a bug
          throw new Error("combined Cayley table: in-window bidegree (" + p + ", " + q
            + ") (target " + target + " <= top " + top
            + ") has no recorded structure constants -- recording is incomplete");
        var dout = t.dims[2], coeffs = [];
        for (var k = 0; k < dout; k++) coeffs.push(((t.constants[k] || [])[i] || [])[j]);
        row.push(cellTex(name, target, coeffs, prime));
      });
      cells.push(row);
    });
    return { corner: PRODUCT_CORNER[name], rowLabels: rowLabels,
      colLabels: colLabels, rowDegsep: rowDegsep, colDegsep: colDegsep, cells: cells,
      dl: totalRows, dr: totalCols, hasBeyond: hasBeyond,
      note: combinedNote(name, tbl, rowMeta, prime) };
  }
  function cayleyBigGrid(c) {
    var head = h("tr");
    head.appendChild(h("th", { "class": "qlgui-corner" }, mathCell(c.corner)));
    c.colLabels.forEach(function (lbl, j) {
      head.appendChild(h("th", c.colDegsep[j] ? { "class": "qlgui-degcol" } : {},
                         mathCell(lbl)));
    });
    var tbl = h("table", { "class": "qlgui-matrix qlgui-cayley" }, head);
    c.cells.forEach(function (row, i) {
      var r = h("tr");
      r.appendChild(h("th", c.rowDegsep[i] ? { "class": "qlgui-degrow" } : {},
                      mathCell(c.rowLabels[i])));
      row.forEach(function (cell, j) {
        var cls = [];
        if (c.rowDegsep[i]) cls.push("qlgui-degrow");
        if (c.colDegsep[j]) cls.push("qlgui-degcol");
        var td = h("td", cls.length ? { "class": cls.join(" ") } : {});
        td.appendChild((cell === "0" || cell === EM_DASH)
                       ? document.createTextNode(cell) : mathCell(cell));
        r.appendChild(td);
      });
      tbl.appendChild(r);
    });
    return tbl;
  }
  // ---- Plan 35 UNIT 2: per-degree explicit representatives (product / Connes) ----
  // The block carries basis_classes / chain_basis / differentials as {side: {degree:
  // ...}} (quiverlab.hochschild.basis_reps). Render one sub-section per (side, degree):
  // the ordered (co)chain enumeration, the explicit classes as term-sum + coordinate
  // vector over that enumeration, and the annihilating differential + a one-line
  // verification sentence -- mirroring quiverlab.trace.render_html.product_degree_sections.
  // `b.differentials` here is the product {side:{degree:...}} shape, read ONLY inside
  // these product renderers (the module-resolution block ships a LIST under the same
  // key -- never read it shape-blind). Tolerant of a block WITHOUT these fields.
  // Marco 2026-07-31: typing statements + Ext/Tor label notes, mirroring the Python
  // single source quiverlab.trace.interpretations (hh_space_typing / module_reps_label_note).
  var NOTATION_TAIL = " Here ⊗ is the tensor product over k, and · inside a word "
    + "is composition of arrows along a path (left to right) -- not a scalar.";
  function hhTyping(theory, route) {
    var coh = (theory === "hh_cohomology" || theory === "HH^");
    if (route === "cs") {
      if (coh)
        return "What the engine computes: Hochschild cohomology as the cohomology of "
          + "Hom_{A^e}(P_•, A), where P_• → A is the Chouhy–Solotar projective bimodule "
          + "resolution with P_n = ⊕_σ A e_{o(σ)} ⊗ e_{t(σ)} A (σ over the degree-n "
          + "ambiguity chains). Degree n collapses to the corner space "
          + "C^n = ⊕_σ e_{o(σ)} A e_{t(σ)}. A basis element v ⊗ p pairs a path v ∈ A "
          + "with a chain word p = a_1·a_2·…. The word p NAMES the free generator of "
          + "P_n attached to that degree-n chain (an iterated overlap of the "
          + "relations): it is not the product of its letters in A -- that product is "
          + "typically zero there, the generator never is. A cochain [p → v] sends "
          + "THAT generator to v ∈ A." + NOTATION_TAIL;
      return "What the engine computes: Hochschild homology as the homology of "
        + "A ⊗_{A^e} P_•, where P_• → A is the Chouhy–Solotar projective bimodule "
        + "resolution with P_n = ⊕_σ A e_{o(σ)} ⊗ e_{t(σ)} A. Degree n collapses to "
        + "the corner space C_n = ⊕_σ e_{t(σ)} A e_{o(σ)}. A basis element v ⊗ p pairs "
        + "a path v ∈ A with a chain word p = a_1·a_2·…. The word p NAMES the free "
        + "generator of P_n attached to that degree-n chain (an iterated overlap of "
        + "the relations): it is not the product of its letters in A -- that product "
        + "is typically zero there, the generator never is." + NOTATION_TAIL;
    }
    if (coh)
      return "What the engine computes: Hochschild cohomology as the cohomology of the "
        + "normalized bar cochain complex. A degree-n cochain is a k-linear map "
        + "C^n = Hom_k(Ā^⊗n, A), where Ā = A/(k·1) is the algebra modulo its unit "
        + "(spanned by the arrows and longer paths); by the tensor–hom adjunction this "
        + "is Hom_{A^e}(A ⊗ Ā^⊗n ⊗ A, A), a bimodule map out of the n-th term of the "
        + "bar resolution of A. A basis functional written [w_1 ⊗ … ⊗ w_n ↦ v] sends "
        + "the single basis tensor w_1 ⊗ … ⊗ w_n (each w_i ∈ Ā) to v ∈ A and every "
        + "other basis tensor to 0. The argument is a TENSOR over k, not a product in "
        + "A: w ⊗ w is a basis element of Ā^⊗2 and stays nonzero even when w·w = 0 "
        + "in A." + NOTATION_TAIL;
    return "What the engine computes: Hochschild homology as the homology of the "
      + "normalized bar chain complex C_n = A ⊗_k Ā^⊗n (equivalently "
      + "A ⊗_{A^e} (A ⊗ Ā^⊗n ⊗ A)). A basis chain written v ⊗ w_1 ⊗ … ⊗ w_n has v ∈ A "
      + "and each w_i ∈ Ā. It is a TENSOR over k, not a product in A: v ⊗ w ⊗ w stays "
      + "nonzero even when w·w = 0 in A." + NOTATION_TAIL;
  }
  function moduleRepsLabelNote(kind) {
    var common = "Notation. In a class term-sum, P_v is the generator of the projective "
      + "summand A e_v of the resolution term P_n, and P_v#k the k-th copy of P_v when "
      + "the vertex v repeats among the summands of P_n; ";
    if (kind === "ext")
      return common + "n_{v,j} is the j-th basis vector of the vertex-v part N e_v of N. "
        + "A basis functional [P_v#k → n_{w,j}] is the A-linear map sending that "
        + "generator to n_{w,j} and every other generator to 0.";
    return common + "n_{v,j} is the j-th basis vector of the vertex-v part e_v N of N. "
      + "A basis element P_v#k ⊗ n_{w,j} is the tensor of that generator with n_{w,j} "
      + "in P_n ⊗_A N.";
  }
  var REPS_ENUM_DISPLAY = 20;
  var REPS_SIDE = {
    coh: { longName: "cohomology", dsym: "\\delta", isCoh: true, letter: "\\alpha",
           cyc: "cocycle" },
    hom: { longName: "homology", dsym: "b", isCoh: false, letter: "z", cyc: "cycle" }
  };
  function prettyLabel(s) {                 // "-> " / "(x)" -> "→" / "⊗" (display)
    return String(s).replace(/->/g, "→").replace(/\(x\)/g, "⊗");
  }
  function coeffSplit(c) {
    c = String(c);
    var neg = c.charAt(0) === "-";
    return { neg: neg, mag: neg ? c.slice(1) : c };
  }
  function signedJoin(pieces) {
    if (!pieces.length) return "0";
    return pieces.map(function (p, i) {
      if (i === 0) return p.neg ? ("-" + p.mag) : p.mag;
      return (p.neg ? " - " : " + ") + p.mag;
    }).join("");
  }
  function termSumText(vector, enumLabels) {  // reuse the UNIT-1 enumeration labels
    if (!Array.isArray(enumLabels) || !enumLabels.length) return null;  // elided -> coord only
    return signedJoin((vector || []).map(function (pair) {
      var cs = coeffSplit(pair[1]), lab = prettyLabel(enumLabels[pair[0]]);
      return { neg: cs.neg, mag: cs.mag === "1" ? lab : cs.mag + " " + lab };
    }));
  }
  function appendRepsEnumeration(div, enumLabels, S, n) {
    var amb = (S.isCoh ? "C^" : "C_") + n;
    div.appendChild(h("p", { text: "Ordered basis of the degree-" + n + " " + S.longName
      + " space " + amb + " (entry k is the k-th basis element; the JSON coordinate vectors index into it):" }));
    if (enumLabels && enumLabels.elided) {
      div.appendChild(h("p", { "class": "qlgui-cites", text: enumLabels.length
        + " elements; the full enumeration is in the report data" }));
      return;
    }
    if (!enumLabels || !enumLabels.length) {
      div.appendChild(h("p", { "class": "qlgui-cites", text: "the space is zero-dimensional" }));
      return;
    }
    var ol = h("ol");
    enumLabels.slice(0, REPS_ENUM_DISPLAY).forEach(function (lbl) {
      ol.appendChild(h("li", { text: prettyLabel(lbl) }));
    });
    div.appendChild(ol);
    if (enumLabels.length > REPS_ENUM_DISPLAY) {
      div.appendChild(h("p", { "class": "qlgui-cites", text: "… "
        + (enumLabels.length - REPS_ENUM_DISPLAY) + " more (full enumeration in the report data)" }));
    }
  }
  function appendRepsClasses(div, classes, enumLabels, S, n) {
    if (!classes || !classes.length) {
      div.appendChild(h("p", { "class": "qlgui-cites",
        text: "no classes (the space is zero in this degree)" }));
      return;
    }
    // Marco 2026-08-02: the (co)homology CLASS list is UNCAPPED ("no limit for the bases
    // of (co)homology") -- only the chain-space enumeration above it is capped.
    div.appendChild(h("p", { text: "Basis classes, each written over the ordered basis above:" }));
    classes.forEach(function (cl, i) {
      var nm = S.letter + "^{" + n + "}_{" + (i + 1) + "}";
      var term = termSumText(cl.vector, enumLabels);
      var p = h("p");
      p.appendChild(h("span", { "class": "arithmatex", text: "\\(" + nm + "\\)" }));
      p.appendChild(document.createTextNode(
        " = " + (term != null ? term : "(recorded in the report data)")));
      div.appendChild(p);
    });
  }
  function appendRepsDifferential(div, diff, S, n, nClasses) {
    if (!diff) return;
    var sym = S.isCoh ? (S.dsym + "^{" + n + "}") : (S.dsym + "_{" + n + "}");
    var lo = S.isCoh ? (n + 1) : Math.max(n - 1, 0);
    var arrow = S.isCoh ? (sym + " : C^{" + n + "} \\to C^{" + (n + 1) + "}")
                        : (sym + " : C_{" + n + "} \\to C_{" + lo + "}");
    div.appendChild(h("p", { "class": "arithmatex", text: "\\(" + arrow + "\\)" }));
    if (diff.elided) {
      var sh = diff.shape || [0, 0];
      div.appendChild(h("p", { "class": "qlgui-cites", text: sh[0] + "×" + sh[1]
        + " matrix (body in the report data; rebuild: " + (diff.note || "") + ")" }));
    } else {
      div.appendChild(matrixGrid(diff.rows || []));
    }
    if (!nClasses) return;
    var sentence = (!S.isCoh && n === 0)
      ? ("every 0-chain is a " + S.cyc + " (" + sym + " vanishes)")
      : ("each " + S.letter + "^{" + n + "}_i is a " + S.cyc + ": applying " + sym
         + " to its coordinate vector gives 0");
    div.appendChild(h("p", { "class": "qlgui-hint", text: "Verification: " + sentence }));
  }
  function appendProductReps(div, b, name, showDiff) {
    var bc = b.basis_classes;
    if (!bc) return false;
    var cb = b.chain_basis || {}, diffs = b.differentials || {};
    var rendered = false;
    ["coh", "hom"].forEach(function (side) {
      var byDeg = bc[side];
      if (!byDeg) return;
      var S = REPS_SIDE[side];
      Object.keys(byDeg).map(Number).sort(function (a, c) { return a - c; })
        .forEach(function (n) {
          rendered = true;
          var key = String(n);
          div.appendChild(h("p", { id: "gui-" + name + "-hh-" + side + "-deg-" + n },
            h("b", { text: "Hochschild " + S.longName + " in degree " + n })));
          if (!(byDeg[key] || []).length) {     // zero space: one line, keep the anchor
            div.appendChild(h("p", { "class": "arithmatex",
              text: "\\(" + (S.isCoh ? "HH^{" : "HH_{") + n + "} = 0\\)" }));
            return;
          }
          appendRepsEnumeration(div, (cb[side] || {})[key], S, n);
          appendRepsClasses(div, byDeg[key], (cb[side] || {})[key], S, n);
          // Marco 2026-07-31: product sections drop the annihilating differential
          // (it lives in the plain HH degree sections); showDiff picks the surface.
          if (showDiff)
            appendRepsDifferential(div, (diffs[side] || {})[key], S, n,
              (byDeg[key] || []).length);
        });
    });
    return rendered;
  }

  // Marco 2026-08-02: for the products, ONE flat list of ALL (co)homology basis classes
  // across degrees (cohomology α^n_i then, for the cap, homology z^n_i), degree-major,
  // each written over the chain basis enumerated in the HH sections above -- no per-degree
  // sub-sections here (that is Marco's point: just remind the classes, all at once, then
  // show the table). Respects the 50-class display cap. Mirrors
  // quiverlab.trace.render_html.product_flat_classes_html.
  function appendProductFlatClasses(div, b) {
    var bc = b.basis_classes;
    if (!bc) return false;
    // Marco 2026-08-02: the flat (co)homology class list is UNCAPPED ("no limit for the
    // bases of (co)homology").
    var cb = b.chain_basis || {}, items = [];
    ["coh", "hom"].forEach(function (side) {
      var byDeg = bc[side];
      if (!byDeg) return;
      var S = REPS_SIDE[side];
      Object.keys(byDeg).map(Number).sort(function (a, c) { return a - c; })
        .forEach(function (n) {
          var key = String(n), enumLabels = (cb[side] || {})[key];
          (byDeg[key] || []).forEach(function (cl, i) {
            items.push({ nm: S.letter + "^{" + n + "}_{" + (i + 1) + "}",
                         term: termSumText(cl.vector, enumLabels) });
          });
        });
    });
    if (!items.length) return false;
    div.appendChild(h("p", { text: "The Hochschild (co)homology basis classes, over the "
      + "chain bases enumerated in the sections above:" }));
    items.forEach(function (it) {
      var p = h("p");
      p.appendChild(h("span", { "class": "arithmatex", text: "\\(" + it.nm + "\\)" }));
      p.appendChild(document.createTextNode(
        " = " + (it.term != null ? it.term : "(recorded in the report data)")));
      div.appendChild(p);
    });
    return true;
  }

  // Plan 35 wave 3d: the PLAIN hh_cohomology / hh_homology blocks. The reps payload is
  // SINGLE-side {str(degree): ...} (like Ext / Tor / HC); we reuse appendProductReps by
  // wrapping the single side, and read off the element-wise classical dictionary
  // (central elements / derivations / deformation cochain / commutator residues) from
  // the captured term-sums -- mirrors quiverlab.trace.interpretations +
  // render_html.hh_element_interpretation. Element-wise ONLY where reps are present.
  function magnitude(coeff, value) {
    return String(coeff) === "1" ? value : coeff + " " + value;
  }
  function arrowWord(word) {
    return Array.isArray(word) ? word.join("·") : String(word);
  }
  function elementFromTerms(terms) {              // degree-0 (co)chain -> element of A
    var parts = (terms || []).map(function (t) { return magnitude(t[0], t[2]); });
    return parts.length ? parts.join(" + ") : "0";
  }
  function groupedReadoff(terms, keyOf, label, sortKeys) {
    var by = {}, order = [];
    (terms || []).forEach(function (t) {
      var k = keyOf(t);
      if (!by[k]) { by[k] = []; order.push(k); }
      by[k].push(magnitude(t[0], t[2]));
    });
    if (sortKeys) order.sort();
    return order.map(function (k) { return label + "(" + k + ") = " + by[k].join(" + "); });
  }
  function derivationValues(terms) {
    return groupedReadoff(terms, function (t) { return arrowWord(t[1]); }, "D", true);
  }
  function deformationCochain(terms) {
    return groupedReadoff(terms, function (t) { return arrowWord(t[1]); }, "μ", false);
  }
  var HH_INTERP = {
    hh_cohomology: {
      0: { head: "HH⁰ = Z(A): each class is a central element z (za = az for all a in A).",
           readoff: function (t) { return [elementFromTerms(t)]; } },
      1: { head: "HH¹ = Der(A) / Inn(A): each class is an outer derivation, read off as "
             + "D(arrow) = value.", readoff: derivationValues },
      2: { head: "HH² = infinitesimal deformations: each class is the 2-cocycle μ(a, b) "
             + "of a first-order deformation a * b = ab + t·μ(a, b).", readoff: deformationCochain }
    },
    hh_homology: {
      0: { head: "HH₀ = A / [A, A]: each class is the residue of an element modulo the "
             + "commutators ab − ba.", readoff: function (t) { return [elementFromTerms(t)]; } }
    }
  };
  function appendHHInterpretation(div, name, b) {
    var bc = b.basis_classes;
    if (!bc) return;
    var letter = name === "hh_cohomology" ? "\\alpha" : "z", byN = HH_INTERP[name] || {};
    Object.keys(bc).map(Number).sort(function (a, c) { return a - c; }).forEach(function (n) {
      var cfg = byN[n], classes = bc[String(n)] || [];
      if (!cfg || !classes.length) return;
      div.appendChild(h("p", {}, h("i", { text: cfg.head })));
      var ul = h("ul", { "class": "qlgui-interp" });
      classes.forEach(function (cl, i) {
        var nm = letter + "^{" + n + "}_{" + (i + 1) + "}";
        var lines = cfg.readoff(cl.terms || []);
        var body = lines.map(function (s) { return String(s).replace("->", "→"); }).join("; ") || "0";
        var li = h("li");
        li.appendChild(h("span", { "class": "arithmatex", text: "\\(" + nm + "\\)" }));
        li.appendChild(document.createTextNode(": " + body));
        ul.appendChild(li);
      });
      div.appendChild(ul);
      if (name === "hh_cohomology" && n === 1 && b.inner_dims) {
        div.appendChild(h("p", { "class": "qlgui-cites",
          text: "inner derivations (the coboundaries a ↦ ax − xa): dimension "
            + (b.inner_dims["1"] != null ? b.inner_dims["1"] : "?") + " = rank δ⁰." }));
      }
    });
  }
  function appendHHReps(div, name, b) {
    if (!b.basis_classes) return;
    var side = name === "hh_cohomology" ? "coh" : "hom", w = {
      basis_classes: {}, chain_basis: {}, differentials: {} };
    w.basis_classes[side] = b.basis_classes;
    w.chain_basis[side] = b.chain_basis || {};
    w.differentials[side] = b.differentials || {};
    appendProductReps(div, w, name, true);      // plain HH keeps the differential
  }

  // Plan 35 wave 3a: module Ext / Tor explicit representatives. The block carries a
  // SINGLE-side {str(degree): ...} payload (Ext cohomological, Tor homological) --
  // distinct from the products {side:{degree}} shape and the module-resolution LIST
  // `differentials`, so this reader is kind-scoped. Mirrors
  // quiverlab.trace.render_html.module_reps_sections; reuses the UNIT-1 term-sum /
  // coordinate-vector / matrix-grid helpers. Tolerant of a block WITHOUT these fields.
  var MODULE_REPS = {
    ext: { letter: "\\alpha", cyc: "cocycle", head: "Ext^",
           ambient: function (n) { return "\\mathrm{Hom}_A(P_{" + n + "}, N)"; },
           arrow: function (n) { return "\\delta^{" + n + "} : \\mathrm{Hom}(P_{" + n
             + "},N) \\to \\mathrm{Hom}(P_{" + (n + 1) + "},N)"; },
           dsym: function (n) { return "\\delta^{" + n + "}"; } },
    tor: { letter: "z", cyc: "cycle", head: "Tor_",
           ambient: function (n) { return "P_{" + n + "} \\otimes_A N"; },
           arrow: function (n) { return "d_{" + n + "} : P_{" + n
             + "}\\otimes_A N \\to P_{" + Math.max(n - 1, 0) + "}\\otimes_A N"; },
           dsym: function (n) { return "d_{" + n + "}"; } }
  };
  function appendModuleRepsDiff(div, diff, cfg, n, nClasses) {
    if (!diff) return;
    div.appendChild(h("p", { "class": "arithmatex", text: "\\(" + cfg.arrow(n) + "\\)" }));
    if (diff.elided) {
      var sh = diff.shape || [0, 0];
      div.appendChild(h("p", { "class": "qlgui-cites", text: sh[0] + "×" + sh[1]
        + " matrix (body in the report data; rebuild: " + (diff.note || "") + ")" }));
    } else if (diff.shape && diff.shape[0] === 0) {
      if (diff.note) div.appendChild(h("p", { "class": "qlgui-cites", text: diff.note }));
      return;
    } else {
      div.appendChild(matrixGrid(diff.rows || []));
    }
    if (!nClasses) return;
    div.appendChild(h("p", { "class": "qlgui-hint", text: "Verification: each "
      + cfg.letter + "^{" + n + "}_i is a " + cfg.cyc + ": applying " + prettyLabel(cfg.dsym(n))
      + " to its coordinate vector gives 0" }));
  }
  function appendModuleReps(div, b, kind) {
    var bc = b.basis_classes;
    if (!bc) return false;
    var cb = b.chain_basis || {}, diffs = b.differentials || {}, cfg = MODULE_REPS[kind];
    var rendered = false;
    div.appendChild(h("p", { "class": "qlgui-cites", text: moduleRepsLabelNote(kind) }));
    // Marco 2026-08-02: the ordered Hom/tensor basis is NOT enumerated per degree here --
    // one pointer states it lives in the report data; the class list stays (capped at 20).
    div.appendChild(h("p", { "class": "qlgui-cites", text: "The ordered basis of each "
      + "Hom/tensor space (into which the coordinate vectors index) is recorded in the "
      + "report data." }));
    Object.keys(bc).map(Number).sort(function (a, c) { return a - c; }).forEach(function (n) {
      rendered = true;
      var key = String(n);
      div.appendChild(h("p", { id: "gui-" + kind + "-deg-" + n }, h("b", {}, [
        h("span", { "class": "arithmatex", text: "\\(" + cfg.head + "{" + n + "}\\)" }),
        document.createTextNode(" in degree " + n)
      ])));
      var classes = bc[key] || [], enumLabels = cb[key];
      if (!classes.length) {                    // zero group: one line, keep the anchor
        div.appendChild(h("p", { "class": "arithmatex",
          text: "\\(" + cfg.head + "{" + n + "} = 0\\)" }));
        return;
      }
      div.appendChild(h("p", { text: "Basis classes, each as its labelled representative:" }));
      classes.slice(0, 20).forEach(function (cl, i) {
        var nm = cfg.letter + "^{" + n + "}_{" + (i + 1) + "}";
        var term = termSumText(cl.vector, enumLabels);
        var p = h("p");
        p.appendChild(h("span", { "class": "arithmatex", text: "\\(" + nm + "\\)" }));
        p.appendChild(document.createTextNode(
          " = " + (term != null ? term : "(recorded in the report data)")));
        div.appendChild(p);
      });
      if (classes.length > 20)
        div.appendChild(h("p", { "class": "qlgui-cites", text: "… and "
          + (classes.length - 20) + " more classes (see the report data)" }));
      appendModuleRepsDiff(div, diffs[key], cfg, n, classes.length);
    });
    return rendered;
  }

  // Plan 35 wave 3c: the classical DICTIONARY of what each (co)homology space MEANS.
  // Mirrors quiverlab.trace.interpretations (the SAME prose, one shared source of truth).
  var DICTIONARY = {
    ext: function (n) {
      if (n === 0) return "Ext⁰(M, N) = Homₐ(M, N): its basis classes are the "
        + "A-module homomorphisms M → N.";
      if (n === 1) return "Each basis class of Ext¹(M, N) is a short exact sequence "
        + "0 → N → E → M → 0 (a Baer extension of M by N), up to "
        + "equivalence. Below, each extension module E is constructed explicitly as a "
        + "pushout and its exactness is verified.";
      return "Each basis class of Extⁿ(M, N) is an n-fold exact sequence "
        + "0 → N → Q → P_{n-2} → ⋯ → P_0 → M → 0 "
        + "(Yoneda), up to equivalence. Below, each is spliced explicitly from the pushout "
        + "module Q and the minimal resolution of M, and its exactness is verified at "
        + "every joint.";
    },
    tor: function (n) {
      if (n === 0) return "Tor₀(M, N) = M ⊗ₐ N, the tensor product itself: "
        + "the coequalizer of the two actions M ⊗ A ⊗ N ⇉ M ⊗ N. Its "
        + "classes are the cosets m ⊗ n.";
      if (n === 1) return "Tor₁(M, N) measures the failure of M (equivalently N) to "
        + "be flat: a nonzero class is a syzygy relation among the generators that "
        + "− ⊗ N does not see, i.e. an obstruction to flatness.";
      return "Torₙ(M, N) is the n-th derived functor of − ⊗ₐ N at M: a "
        + "higher syzygy / flatness obstruction (homological framing).";
    },
    hh_cohomology: function (n) {
      if (n === 0) return "HH⁰(A) = Z(A), the CENTRE of A: the classes are the "
        + "central elements z (those with za = az for all a in A).";
      if (n === 1) return "HH¹(A) = Der(A) / Inn(A), the OUTER DERIVATIONS: each "
        + "class is a derivation D : A → A (a k-linear map with the Leibniz rule "
        + "D(ab) = D(a) b + a D(b)), determined by its values on the arrow generators, "
        + "taken modulo the inner derivations a ↦ ax − xa.";
      if (n === 2) return "HH²(A) classifies the INFINITESIMAL DEFORMATIONS of A: "
        + "each class is a 2-cocycle μ(a, b) giving a first-order (square-zero) "
        + "deformation of the multiplication a * b = ab + t·μ(a, b); the "
        + "coboundaries are the trivial (gauge) deformations.";
      return "HHⁿ(A) controls the higher obstructions to deforming A (its Yoneda "
        + "product with HH² carries the obstruction cocycles); homological framing.";
    },
    hh_homology: function (n) {
      if (n === 0) return "HH₀(A) = A / [A, A], the COMMUTATOR QUOTIENT: the classes "
        + "are the residues of A modulo the subspace spanned by the commutators ab − ba.";
      return "HHₙ(A) is the n-th Hochschild homology, Tor^{A^e}_n(A, A) -- a "
        + "derived-functor / cyclic-theory invariant (homological framing).";
    },
    cyclic_homology: function (n) {
      if (n === 0) return "HC₀(A) = A / [A, A]: the same space as HH₀, read as "
        + "the TRACE FUNCTIONALS on A (a trace τ with τ(ab) = τ(ba) is "
        + "exactly a linear form on A / [A, A]).";
      return "HCₙ(A) is cyclic homology in degree n, the homology of Connes' (b, B) "
        + "total complex -- it packages the S, B, I periodicity of the Hochschild theory "
        + "(homological framing).";
    }
  };
  var DICTIONARY_ALIAS = { "HH^": "hh_cohomology", "HH_": "hh_homology", "HC_": "cyclic_homology" };

  function appendDictionaryFraming(div, theory, dims) {
    var fn = DICTIONARY[DICTIONARY_ALIAS[theory] || theory];
    if (!fn || !dims) return;
    var lis = [], seen = {};
    for (var n = 0; n < dims.length; n++) {
      var s = fn(n);
      if (!s || seen[s]) continue;
      seen[s] = true;
      lis.push(h("li", { text: s }));
    }
    if (!lis.length) return;
    div.appendChild(h("p", {}, h("i", { text: "Interpretation of the spaces (the classical dictionary):" })));
    div.appendChild(h("ul", { "class": "qlgui-interp" }, lis));
  }

  // Plan 35 wave 3c: the Yoneda exact-sequence interpretation of every Ext class -- the
  // constructed + self-certified exact sequence 0 -> N -> Q -> ... -> M -> 0. Mirrors
  // quiverlab.trace.render_html.ext_interpretation_sections. Tolerant of an old-cache
  // block with no `interpretation` field.
  function dvInline(dv) {
    if (!dv) return "()";
    return "(" + Object.keys(dv).sort().map(function (k) { return dv[k]; }).join(", ") + ")";
  }
  function yonedaMiddle(div, mid) {
    var label = mid.label || "Q", dv = dvInline(mid.dimvec);
    if (mid.standard) {
      var sym = { simple: "S", projective: "P", injective: "I" }[mid.standard.kind] || "?";
      var p = h("p");
      p.appendChild(h("span", { "class": "arithmatex", text: "\\(" + label + "\\)" }));
      p.appendChild(document.createTextNode(" is the standard indecomposable "));
      p.appendChild(h("span", { "class": "arithmatex", text: "\\(" + sym + "_{" + mid.standard.vertex + "}\\)" }));
      p.appendChild(document.createTextNode(" (dimension vector " + dv + ")."));
      div.appendChild(p);
      return;
    }
    div.appendChild(h("p", { text: label + ", the extension module, dimension vector " + dv + " — its action:" }));
    appendRepMaps(div, label, mid.module || {});
    if (mid.module && mid.module.display_only)
      div.appendChild(h("p", { "class": "qlgui-cites",
        text: "display only — entries lie outside the integer/fraction grammar" }));
  }
  function yonedaFacts(facts) {
    return (facts || []).map(function (f) {
      if (f.fact === "injective") return "injective at " + f.node + " (rank " + f.rank + " = dim " + f.dim + ")";
      if (f.fact === "surjective") return "surjective at " + f.node + " (rank " + f.rank + " = dim " + f.dim + ")";
      if (f.fact === "im=ker") return "image = kernel at " + f.node + " (" + f.rank_in + " + " + f.rank_out + " = " + f.dim + ")";
      return "";
    }).filter(Boolean).join("; ");
  }
  function appendExtInterpretation(div, interp) {
    if (!interp || !interp.sequences) return;
    var seqs = interp.sequences;
    Object.keys(seqs).map(Number).sort(function (a, c) { return a - c; }).forEach(function (n) {
      var list = seqs[String(n)] || [];
      if (!list.length) return;
      var head = h("h4", { id: "gui-ext-yoneda-deg-" + n });
      head.appendChild(document.createTextNode("Interpretation: "));
      head.appendChild(h("span", { "class": "arithmatex", text: "\\(\\mathrm{Ext}^{" + n + "}(M, N)\\)" }));
      head.appendChild(document.createTextNode(" as exact sequences"));
      div.appendChild(head);
      div.appendChild(h("p", {}, h("i", { text: DICTIONARY.ext(n) })));
      list.forEach(function (seq) {
        div.appendChild(h("h5", {}, h("span", { "class": "arithmatex", text: "\\(" + seq.class_name + "\\)" })));
        if (!seq.certified) {
          div.appendChild(h("p", { "class": "qlgui-cites",
            text: "this class's exact sequence could not be certified (" + (seq.error || "unknown")
              + "); it is omitted rather than shown wrongly." }));
          return;
        }
        var mods = seq.modules || [];
        var line = "0 \\to " + mods.map(function (m) { return m.label || "?"; }).join(" \\to ") + " \\to 0";
        div.appendChild(h("p", {}, h("span", { "class": "arithmatex", text: "\\(" + line + "\\)" })));
        div.appendChild(h("p", { "class": "qlgui-cites",
          text: "dimension vectors: " + mods.map(function (m) { return (m.label || "?") + " " + dvInline(m.dimvec); }).join(", ") + "." }));
        var mid = mods.filter(function (m) { return m.role === "middle"; })[0];
        if (mid) yonedaMiddle(div, mid);
        (seq.maps || []).forEach(function (mp) {
          div.appendChild(h("p", {}, h("span", { "class": "arithmatex",
            text: "\\(" + (mp.from || "?") + " \\to " + (mp.to || "?") + "\\)" })));
          if (mp.elided) {
            var sh = mp.shape || [0, 0];
            div.appendChild(h("p", { "class": "qlgui-cites",
              text: sh[0] + "×" + sh[1] + " matrix (body in the machine record)" }));
          } else {
            div.appendChild(matrixGrid(mp.rows || []));
          }
        });
        var fs = yonedaFacts(seq.facts);
        if (fs) div.appendChild(h("p", { "class": "qlgui-cites", text: "Exactness verified: " + fs + "." }));
      });
    });
  }

  // Plan 35 wave 3b: cyclic homology HC explicit representatives. The block carries a
  // SINGLE-side {str(degree): ...} payload (HC is homological) PLUS a column_structure.
  // Marco 2026-08-02: the Tot column enumerations (the Tot_n = C_n (+) C_{n-2} (+) ...
  // decomposition heading and the ordered Tot_n basis) are NOT re-listed -- one pointer
  // states they live in the report data; each degree keeps its class list + total
  // differential. Mirrors quiverlab.trace.render_html.cyclic_degree_sections.
  function appendCyclicClasses(div, classes, enumLabels, n) {
    if (!classes.length) {
      div.appendChild(h("p", { "class": "qlgui-cites", text: "no classes (HC_" + n + " is zero)" }));
      return;
    }
    div.appendChild(h("p", { text: "Basis classes, each as its labelled representative:" }));
    classes.slice(0, 20).forEach(function (cl, i) {
      var nm = "z^{" + n + "}_{" + (i + 1) + "}";
      var term = termSumText(cl.vector, enumLabels);
      var p = h("p");
      p.appendChild(h("span", { "class": "arithmatex", text: "\\(" + nm + "\\)" }));
      p.appendChild(document.createTextNode(
        " = " + (term != null ? term : "(recorded in the report data)")));
      div.appendChild(p);
    });
    if (classes.length > 20)
      div.appendChild(h("p", { "class": "qlgui-cites", text: "… and "
        + (classes.length - 20) + " more classes (see the report data)" }));
  }
  function appendCyclicDiff(div, diff, n, nClasses) {
    if (!diff) return;
    var arrow = "D_{" + n + "} : \\mathrm{Tot}_{" + n + "} \\to \\mathrm{Tot}_{"
      + Math.max(n - 1, 0) + "}";
    div.appendChild(h("p", { "class": "arithmatex", text: "\\(" + arrow + "\\)" }));
    if (diff.elided) {
      var sh = diff.shape || [0, 0];
      div.appendChild(h("p", { "class": "qlgui-cites", text: sh[0] + "×" + sh[1]
        + " matrix (body in the report data; rebuild: " + (diff.note || "") + ")" }));
    } else if (diff.shape && diff.shape[0] === 0) {
      if (diff.note) div.appendChild(h("p", { "class": "qlgui-cites", text: diff.note }));
      return;
    } else {
      div.appendChild(matrixGrid(diff.rows || []));
    }
    if (!nClasses) return;
    div.appendChild(h("p", { "class": "qlgui-hint", text: "Verification: each z^{" + n
      + "}_i is a cycle of the total complex: applying D_{" + n + "} = b + B to its "
      + "coordinate vector gives 0" }));
  }
  function appendCyclicReps(div, b) {
    var bc = b.basis_classes;
    if (!bc) return false;
    var cb = b.chain_basis || {}, diffs = b.differentials || {};
    var rendered = false;
    // Marco 2026-08-02: one pointer for the whole section (the Tot column structure +
    // ordered basis are in the report data; the coordinate vectors index into it).
    div.appendChild(h("p", { "class": "qlgui-cites", text: "The total complex "
      + "Tot_n = C_n ⊕ C_{n-2} ⊕ … and the ordered basis of each Tot_n (into which the "
      + "coordinate vectors index) are recorded in the report data." }));
    Object.keys(bc).map(Number).sort(function (a, c) { return a - c; }).forEach(function (n) {
      rendered = true;
      var key = String(n);
      div.appendChild(h("p", { id: "gui-cyclic-hc-deg-" + n }, h("b", {}, [
        h("span", { "class": "arithmatex", text: "\\(HC_{" + n + "}\\)" }),
        document.createTextNode(" in degree " + n)
      ])));
      if (!(bc[key] || []).length) {            // zero space: one line, keep the anchor
        div.appendChild(h("p", { "class": "arithmatex", text: "\\(HC_{" + n + "} = 0\\)" }));
        return;
      }
      appendCyclicClasses(div, bc[key] || [], cb[key], n);
      appendCyclicDiff(div, diffs[key], n, (bc[key] || []).length);
    });
    return rendered;
  }

  function renderProductTables(div, name, b) {
    var prime = primeFromBasis(b.basis);
    div.appendChild(h("p", { text: PRODUCT_TITLE[name] }));
    div.appendChild(h("p", { "class": "qlgui-hint", text: productLegend(name, b) }));
    if (prime != null)
      div.appendChild(h("p", { "class": "qlgui-hint", text: balancedRepNote(prime) }));
    appendProductFlatClasses(div, b);      // flat class list, then the table right away
    var tables = b.tables || [];
    // Marco 2026-07-31: a product FAMILY whose every bidegree vanishes collapses to
    // one section-level line -- no empty tables.
    var allZero = tables.length && tables.every(function (t) {
      var d = t.dims || [0, 0, 0];
      return !d[0] || !d[1] || (t.constants || []).every(matIsZero);
    });
    if (allZero) {
      var fam = { cup: "cup products", cap: "cap products",
                  bracket: "Gerstenhaber brackets" }[name] || "products";
      div.appendChild(h("p", { "class": "qlgui-hint",
        text: "All " + fam + " in the served bidegrees vanish." }));
    } else {
      // ONE big degree-major Cayley table for the family (Marco 2026-08-01), UNCAPPED
      // (Marco 2026-08-02: product tables can be big; the one table always renders).
      div.appendChild(h("p", { "class": "arithmatex",
        text: "\\(" + FAMILY_HEADING[name] + "\\)" }));
      div.appendChild(h("p", { "class": "qlgui-hint", text: FAMILY_AXIS_NOTE }));
      var c = combinedCayley(name, tables, prime);
      if (c.note) div.appendChild(h("p", { "class": "qlgui-hint", text: c.note }));
      if (c.hasBeyond)
        div.appendChild(h("p", { "class": "qlgui-hint", text: beyondWindowNote() }));
      div.appendChild(cayleyBigGrid(c));
    }
    if (name === "bracket" && b.window != null) {
      div.appendChild(h("p", { "class": "qlgui-hint",
        text: "served to degree window " + b.window + " (bar-transport bound)" }));
    }
    div.appendChild(h("div", { "class": "qlgui-cites", text: engineNote(b.engine) }));
  }
  function renderConnesB(div, b) {
    div.appendChild(h("p", { text: PRODUCT_TITLE.connes_b }));
    div.appendChild(h("p", { "class": "qlgui-hint",
      text: "each induced Connes differential B_n: HH_n → HH_{n+1} is written on "
          + "the recorded homology bases — rows index HH_{n+1}, columns index HH_n. "
          + "The cycle classes z^n_j are listed below, all degrees at once." }));
    // Marco 2026-08-02: ONE flat homology (z^n_j) class list -- the SAME builder the
    // products use -- then the induced B matrices per degree with rank lines. No
    // per-degree sub-sections and no chain enumerations (those live in the HH sections).
    appendProductFlatClasses(div, b);
    div.appendChild(h("p", { "class": "qlgui-cites", text: "Each cycle class is written "
      + "over the chain basis enumerated in the Hochschild homology sections / report "
      + "data; the induced matrices below act on these classes." }));
    var keys = Object.keys(b.matrices || {})
      .map(Number).sort(function (a, c) { return a - c; });
    keys.forEach(function (n) {
      div.appendChild(h("p", { "class": "arithmatex",
        text: "\\( B_{" + n + "} : HH_{" + n + "} \\to HH_{" + (n + 1) + "} \\)" }));
      div.appendChild(matrixGrid(b.matrices[String(n)]));
      div.appendChild(h("p", { text: "rank B_" + n + " = " + b.ranks[String(n)] }));
    });
    div.appendChild(h("div", { "class": "qlgui-cites", text: engineNote(b.engine) }));
  }

  function renderBlock(res) {
    var b = res.block, name = res.invariant.split(":")[0];
    var div = h("div", { "class": "qlgui-block" });
    if (name === "hh_cohomology" || name === "hh_homology") {
      var sup = name === "hh_cohomology";
      // Typing statement at the top (Marco 2026-07-31): exactly what the engine
      // computes and what the bar-bracket / tensor notation means.
      var route = /cs|chouhy|solotar/.test(String(b.engine || "").toLowerCase())
        ? "cs" : "bar";
      div.appendChild(h("p", { "class": "qlgui-hint", text: hhTyping(name, route) }));
      var head = h("tr"), row = h("tr");
      head.appendChild(h("th", { text: "n" }));
      row.appendChild(h("th", { text: sup ? "dim HH^n" : "dim HH_n" }));
      b.dims.forEach(function (d, n) {
        head.appendChild(h("td", { text: String(n) }));
        row.appendChild(h("td", { text: String(d) }));
      });
      div.appendChild(h("p", { text: sup ? "Hochschild cohomology" : "Hochschild homology" }));
      div.appendChild(h("table", {}, head, row));
      div.appendChild(h("div", { "class": "qlgui-cites", text: engineNote(b.engine) }));
      appendDictionaryFraming(div, name, b.dims);
      // Plan 35 wave 3d: the element-wise dictionary read-offs + per-degree explicit
      // representatives (when the block carries them; tolerant of an old-cache block).
      appendHHInterpretation(div, name, b);
      if (b.basis_classes)
        div.appendChild(h("p", {}, h("b", { text: "Explicit representatives by degree:" })));
      appendHHReps(div, name, b);
    } else if (name === "cup" || name === "cap" || name === "bracket") {
      renderProductTables(div, name, b);
    } else if (name === "connes_b") {
      renderConnesB(div, b);
    } else if (name === "cyclic_homology") {
      // Plan-35 follow-up: HC is a homology-style subscript table HC_n, rendered
      // exactly like the HH dims tables above.
      var chead = h("tr"), crow = h("tr");
      chead.appendChild(h("th", { text: "n" }));
      crow.appendChild(h("th", { text: "dim HC_n" }));
      b.dims.forEach(function (d, n) {
        chead.appendChild(h("td", { text: String(n) }));
        crow.appendChild(h("td", { text: String(d) }));
      });
      div.appendChild(h("p", { text: "Cyclic homology" }));
      div.appendChild(h("table", {}, chead, crow));
      div.appendChild(h("div", { "class": "qlgui-cites", text: engineNote(b.engine) }));
      appendDictionaryFraming(div, name, b.dims);
      // Plan 35 wave 3b: the per-degree explicit representatives over the total complex.
      if (b.basis_classes)
        div.appendChild(h("p", {}, h("b", { text: "Explicit representatives by degree:" })));
      appendCyclicReps(div, b);
    } else if (name === "ss_hochschild") {
      // Plan 42: the Hochschild (b, B) spectral sequence. Abutment table (E_inf
      // totals == HC_n), the netPage E_inf grid, and the convergence prose. A loud
      // DepthLimit guard is reported as an honest error line, never a crash.
      if (b.error) {
        div.appendChild(h("p", { text: b.error }));
      } else {
        var shead = h("tr"), srow = h("tr");
        shead.appendChild(h("th", { text: "n" }));
        srow.appendChild(h("th", { text: "dim E_inf total (= HC_n)" }));
        (b.abutment || []).forEach(function (d, n) {
          shead.appendChild(h("td", { text: String(n) }));
          srow.appendChild(h("td", { text: String(d) }));
        });
        div.appendChild(h("p", { text: "Hochschild (b,B) spectral sequence" }));
        div.appendChild(h("table", {}, shead, srow));
        if (b.grid)
          div.appendChild(h("pre", { text: b.grid.replace(/```/g, "").trim() }));
        if (b.prose) div.appendChild(h("p", { text: b.prose }));
      }
    } else if (name === "cartan") {
      div.appendChild(h("p", { text: "Cartan matrix:" }));
      div.appendChild(matrixGrid(b.matrix));
    } else if (name === "coxeter_polynomial") {
      div.appendChild(h("p", { "class": "arithmatex", text: "\\[ \\chi(t) = " + b.latex + " \\]" }));
    } else if (name === "global_dimension") {
      div.appendChild(h("p", { text: b.text }));
    } else if (name === "center") {
      div.appendChild(h("p", { "class": "arithmatex", text: "\\( \\dim Z(A) = " + b.dim + " \\)" }));
    } else if (name === "tau" || name === "tau_minus") {
      // M's translate WITH its full per-arrow matrices, then the same for the
      // second module N when the request named one (Marco, 2026-07-29).
      appendTranslate(div, b, name, "M");
      (b.targets || []).forEach(function (t) {
        var role = TARGET_ROLE[t.role];
        div.appendChild(h("p", { text: role ? "and for N, " + role + ":" : "and for N:" }));
        appendTranslate(div, t, name, "N");
      });
    } else if (name === "dimension_vector" || name === "projective_dimension" ||
               name === "injective_dimension") {
      // A pre-2026-07-29 cached pd/id block carries no `latex` (that missing key is
      // exactly what typeset as a literal "undefined"); compose it from the value.
      div.appendChild(h("p", { "class": "arithmatex",
        text: "\\[ " + (b.latex || homdimLatex(name, b)) + " \\]" }));
      if (b.note) div.appendChild(h("p", { "class": "qlgui-hint", text: b.note }));
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
        if (b.series && b.series.length) {   // the Loewy (radical) series (Plan 37)
          div.appendChild(h("p", { text: "Loewy (radical) series, top to bottom:" }));
          div.appendChild(loewySeriesTable(b.series));
        }
        trio.forEach(function (p) { appendRepMaps(div, p[0], p[1]); });
      }
    } else if (name === "decompose") {
      div.appendChild(h("p", { text: "Krull–Schmidt decomposition — " + b.iso_classes +
        " indecomposable summand(s):" }));
      div.appendChild(decompTable(b.summands));
      appendSummandMaps(div, b.summands);
    } else if (name === "ext" || name === "tor") {
      var isExt = name === "ext";
      div.appendChild(h("p", { text: (isExt ? "Ext to the target module — dim vector "
        : "Tor with the target left module — dim vector ") + dvText(b.target.dimvec) + ":" }));
      if (b.resolved) {
        // Marco 2026-08-03: name the objects BEFORE the numbers -- which module
        // was resolved, and by which resolution.
        div.appendChild(h("p", { text: "Object resolved: the " + b.resolved.side
          + " A-module M, by its " + b.resolved.resolution + "; "
          + (isExt ? "Ext^n(M, N) = H^n(Hom_A(P_\u2022, N))."
                   : "Tor_n(M, N) = H_n(P_\u2022 \u2297_A N).") }));
      }
      if (b.resolution && b.resolution.summands && b.resolution.summands.length) {
        // ... and SHOW that resolution before the data (terms to the depth used).
        div.appendChild(h("p", { text: "The resolution of M used:" }));
        var rtbl = h("table", { "class": "qlgui-table" });
        var rhead = h("tr"); rhead.appendChild(h("th", { text: "n" }));
        rhead.appendChild(h("th", { text: "P_n" })); rtbl.appendChild(rhead);
        b.resolution.summands.forEach(function (tex, n) {
          var tr = h("tr"); tr.appendChild(h("td", { text: String(n) }));
          var td = h("td", { "class": "arithmatex", text: "\\(" + tex + "\\)" });
          tr.appendChild(td); rtbl.appendChild(tr);
        });
        div.appendChild(rtbl);
        div.appendChild(h("p", { "class": "qlgui-hint", text: isExt
          ? "N is not resolved: it enters through Hom_A(\u2212, N) applied to this resolution."
          : "N is not resolved: it enters through \u2212 \u2297_A N applied to this resolution." }));
      }
      div.appendChild(degreeTable(isExt ? "dim Ext^n" : "dim Tor_n", b.dims));
      appendDictionaryFraming(div, name, b.dims);
      if (b.basis_classes)
        div.appendChild(h("p", {}, h("b", { text: "Explicit representatives by degree:" })));
      appendModuleReps(div, b, name);
      if (isExt) appendExtInterpretation(div, b.interpretation);
    } else if (name === "projective_resolution" || name === "injective_resolution") {
      var proj = name === "projective_resolution";
      div.appendChild(h("p", { text: proj ? "projective resolution" : "injective resolution" }));
      div.appendChild(resTable(b));
      var d = proj ? b.pd : b.injective_dimension;
      div.appendChild(h("p", { text: (proj ? "pd = " : "id = ") +
        (d == null ? "∞ (beyond the probed length)" : String(d)) }));
      appendTermBasis(div, b, proj);
      appendDifferentials(div, b, proj);
    } else if (name === "ext_algebra") {
      // Plan 38: the three-valued Koszul verdict + graded (Betti) data of E(A).
      var kv;
      if (b.koszul === true) kv = "A is Koszul" + (b.koszul_reason ? " — " + b.koszul_reason : "");
      else if (b.koszul === false) kv = b.obstruction
        ? "A is not Koszul — obstruction at degree " + b.obstruction[0] + " (" + b.obstruction[1] + ")"
        : "A is not Koszul";
      else kv = "Koszulity undecided through degree " + b.certified_through_degree
        + (b.koszul_reason ? " — " + b.koszul_reason : "");
      div.appendChild(h("p", { text: kv + "." }));
      div.appendChild(degreeTable("dim E^n", b.graded_dims || []));
      if (b.latex)
        div.appendChild(h("p", { "class": "arithmatex", text: "\\[ " + b.latex + " \\]" }));
      var byDeg = function (d) {
        d = d || {};
        var ks = Object.keys(d).sort(function (a, c) { return (+a) - (+c); });
        return ks.length ? ks.map(function (k) { return "degree " + k + ": " + d[k]; }).join(", ") : "none";
      };
      div.appendChild(h("p", { text: "Minimal generators of E(A): " + byDeg(b.generators_by_degree)
        + "; minimal relations: " + byDeg(b.relations_by_degree) + "." }));
    } else if (name === "recognizers") {
      // Plan 38: the recognizer flags + Dynkin/Euclidean type + form type.
      var RLBL = { is_semisimple: "semisimple", is_radical_square_zero: "radical square zero",
        is_hereditary: "hereditary", is_basic: "basic", is_nakayama: "Nakayama",
        is_special_biserial: "special biserial", is_string: "string", is_gentle: "gentle",
        is_selfinjective: "self-injective", is_symmetric: "symmetric" };
      var RORD = ["is_semisimple", "is_radical_square_zero", "is_hereditary", "is_basic",
        "is_nakayama", "is_special_biserial", "is_string", "is_gentle",
        "is_selfinjective", "is_symmetric"];
      var ul = h("ul");
      RORD.forEach(function (k) {
        if (!b.flags || !(k in b.flags)) return;
        var v = b.flags[k], txt;
        if (v && typeof v === "object" && "error" in v) txt = RLBL[k] + ": not decided — " + v.error;
        else txt = RLBL[k] + ": " + (v === true ? "yes" : "no");
        ul.appendChild(h("li", { text: txt }));
      });
      div.appendChild(ul);
      div.appendChild(h("p", { text: "Diagram type: "
        + (b.dynkin_type || "not a Dynkin/Euclidean diagram") }));
      div.appendChild(h("p", { text: "Form type: "
        + (b.form_type || "undefined (Cartan not unimodular)") }));
    }
    div.appendChild(citesLine(b));
    el.results.appendChild(div);
    if (window.MathJax && window.MathJax.typesetPromise) {
      // Full-page sweep, NOT typesetPromise([div]): with explicit roots the
      // walker never consults the root's own class, so the site's
      // ignoreHtmlClass ".*|" config silently skips the block (found live).
      // fitMath runs AFTER typesetting -- it measures the rendered width.
      var done = window.MathJax.typesetPromise();
      if (done && done.then) done.then(fitMath); else fitMath();
    } else {
      fitMath();
    }
  }

  // ---------- buttons ----------
  el.compute.addEventListener("click", function () {
    if (S.busy || !S.engineReady) return;
    el.results.innerHTML = "";
    S.artifacts = { tikz: "", snippet: "", bundle: "", traceHtml: "", traceJson: "" };
    el.print.disabled = el.tikz.disabled = el.json.disabled = el.snippet.disabled = true;
    el["report-html"].disabled = el["report-json"].disabled = true;
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
  // The shrink-to-fit factor depends on the column width, so re-measure on resize
  // (debounced) -- a widened window must give the matrices their full size back.
  var fitTimer = null;
  window.addEventListener("resize", function () {
    if (fitTimer) clearTimeout(fitTimer);
    fitTimer = setTimeout(fitMath, 150);
  });
  el.relations.addEventListener("input", scheduleProbe);
  [el.field, el.p, el.n, el.hhc, el["hhc-top"], el.hhh, el["hhh-top"],
   el.cup, el["cup-top"], el.cap, el["cap-top"], el.bracket, el["bracket-top"],
   el.connes_b, el["connes_b-top"],
   el.cyclic_homology, el["cyclic_homology-top"],
   el.ss_hochschild, el["ss_hochschild-top"], el.cartan,
   el.coxeter_polynomial, el.global_dimension, el.center,
   el.ext_algebra, el["ext_algebra-top"], el.recognizers]
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
