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

// A summand recognised as a STANDARD indecomposable is NAMED S_v / P_v / I_v and
// needs no matrices; any other one is shown in full, since its dimension vector
// does not determine it (Marco 2026-07-29).
const STD_SYM = { simple: "S", projective: "P", injective: "I" };
function summandName(s, i) {
  const std = s && s.standard;
  if (std && STD_SYM[std.kind]) return STD_SYM[std.kind] + "_" + std.vertex;
  return "M_" + i;
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
    tr.appendChild(cellText(summandName(s, i + 1)));
    tr.appendChild(cellText(String(s.multiplicity)));
    tr.appendChild(cellText(dvText(s.dim_vector)));
    tbl.appendChild(tr);
  });
  wrap.appendChild(tbl);
  (block.summands || []).forEach(function (s, i) {
    if (s.standard || !s.maps) return;
    appendRepMaps(wrap, summandName(s, i + 1), s, d);
  });
  return wrap;
}

// Plan 41 (C3): the almost-split (Auslander–Reiten) sequence 0 → τM → E → M → 0,
// with τM as a full representation and E's Krull–Schmidt summands; an honest refusal
// block for a projective / decomposable / undecidable input.
function almostSplitBlock(block, d) {
  const wrap = document.createElement("div");
  if (block.exists === false) {
    const pr = document.createElement("p");
    pr.textContent = (d.modAlmostSplitRefused || "No almost-split sequence")
      + " — " + (block.reason || "input not eligible") + ".";
    wrap.appendChild(pr);
    return wrap;
  }
  const p = document.createElement("p");
  p.textContent = (d.modAlmostSplitSeq || "Almost-split sequence")
    + ":  0 → τM → E → M → 0";
  wrap.appendChild(p);
  appendRepMaps(wrap, "τM", block.tau, d);
  const p2 = document.createElement("p");
  p2.textContent = (d.modAlmostSplitMiddle || "middle term E — Krull–Schmidt summands");
  wrap.appendChild(p2);
  const summ = (block.middle && block.middle.summands) || [];
  const tbl = tableWithHeader([d.modSummand || "summand",
                               d.modMult || "multiplicity",
                               d.modDimvec || "dim vector"]);
  summ.forEach(function (s, i) {
    const tr = document.createElement("tr");
    tr.appendChild(cellText(summandName(s, i + 1)));
    tr.appendChild(cellText(String(s.multiplicity)));
    tr.appendChild(cellText(dvText(s.dim_vector)));
    tbl.appendChild(tr);
  });
  wrap.appendChild(tbl);
  summ.forEach(function (s, i) {
    if (s.standard || !s.maps) return;
    appendRepMaps(wrap, summandName(s, i + 1), s, d);
  });
  return wrap;
}

// Plan 44 (C7): the tilting_check verdict as a small key/value table + a one-line
// verdict sentence. Labels are i18n via the form dataset (fallback English).
function tiltingBlock(block, d) {
  const wrap = document.createElement("div");
  if (block.error) {
    const e = document.createElement("p");
    e.className = "error";
    e.textContent = block.error;
    wrap.appendChild(e);
    return wrap;
  }
  const p = document.createElement("p");
  p.textContent = block.is_tilting
    ? (d.modTiltingYes || "M is a tilting module.")
    : (d.modTiltingNo || "M is not a tilting module") + ": " + block.note + ".";
  wrap.appendChild(p);
  const rows = [
    [d.modTiltingRow || "tilting", block.is_tilting ? "yes" : "no"],
    [d.modTiltingN || "n (pd bound tested)", String(block.n)],
    ["pd M", (block.pd === null || block.pd === undefined) ? "> bound" : String(block.pd)],
    [d.modTiltingSelfExt || "Ext^i(M,M)=0 (1≤i≤n)",
     block.self_ext_vanishes ? "yes" : "no"],
    [d.modTiltingSummands || "# non-iso indec. summands", String(block.num_summands)],
    [d.modTiltingVertices || "# vertices (rank K_0)", String(block.num_vertices)]
  ];
  const tbl = document.createElement("table");
  rows.forEach(function (r) {
    const tr = document.createElement("tr");
    const th = document.createElement("th");
    th.textContent = r[0];
    const td = document.createElement("td");
    td.textContent = r[1];
    tr.appendChild(th); tr.appendChild(td);
    tbl.appendChild(tr);
  });
  wrap.appendChild(tbl);
  return wrap;
}

// Plan 49 (C8): orbit geometry -- orbit dim + Voigt rigidity + HONEST codim +
// (hereditary Dynkin) the Kac canonical decomposition. Labels i18n via the form
// dataset (fallback English).
function orbitGeometryBlock(block, d) {
  const wrap = document.createElement("div");
  if (block.error) {
    const e = document.createElement("p");
    e.className = "error";
    e.textContent = block.error;
    wrap.appendChild(e);
    return wrap;
  }
  const dvText = function (dv) {
    return "{" + Object.keys(dv || {}).map(function (k) { return k + ": " + dv[k]; }).join(", ") + "}";
  };
  const rows = [
    [d.ogDimVec || "dimension vector d", dvText(block.dim_vector)],
    [d.ogGroupDim || "dim GL(d) = Σ d_v²", String(block.group_dim)],
    [d.ogRepDim || "dim Rep(Q,d) (ambient)", String(block.rep_variety_dim)],
    [d.ogEndDim || "dim End_A(M)", String(block.end_dim)],
    [d.ogOrbitDim || "dim O_M (orbit)", String(block.orbit_dim)],
    [d.ogExt1 || "dim Ext¹(M,M)", String(block.ext1_self)]
  ];
  const tbl = document.createElement("table");
  rows.forEach(function (r) {
    const tr = document.createElement("tr");
    const th = document.createElement("th");
    th.textContent = r[0];
    const td = document.createElement("td");
    td.textContent = r[1];
    tr.appendChild(th); tr.appendChild(td);
    tbl.appendChild(tr);
  });
  wrap.appendChild(tbl);
  const verdict = block.rigid
    ? (d.ogRigid || "M is rigid: Ext¹(M,M) = 0, so the orbit O_M is open (Voigt).")
    : (d.ogNotRigid || "M is not rigid: Ext¹(M,M) > 0.");
  const gloss = block.codim_semantics === "hereditary"
    ? (d.ogCodimHered || " A is hereditary, so dim Ext¹(M,M) IS the codimension of the "
        + "orbit closure in Rep(Q,d) (Voigt; Rep smooth).")
    : (d.ogCodimGeneral || " A = kQ/I is not hereditary, so dim Ext¹(M,M) is only an "
        + "UPPER BOUND on the codimension (the module variety is cut by the relations).");
  const pv = document.createElement("p");
  pv.textContent = verdict + gloss;
  wrap.appendChild(pv);
  if (block.canonical_decomposition && block.canonical_decomposition.length) {
    const parts = block.canonical_decomposition.map(function (c) {
      const nm = c.name || ("(" + (c.root || []).join(", ") + ")");
      return c.multiplicity === 1 ? nm : nm + "^" + c.multiplicity;
    });
    const pc = document.createElement("p");
    pc.textContent = (d.ogCanonical || "Kac canonical decomposition") + ": d = "
      + parts.join(" ⊕ ") + " (each component a positive root; the generic module is rigid).";
    wrap.appendChild(pc);
  } else if (block.canonical_note) {
    const pn = document.createElement("p");
    pn.className = "hint";
    pn.textContent = (d.ogCanonicalNone || "Canonical decomposition not computed") + ": "
      + block.canonical_note;
    wrap.appendChild(pn);
  }
  return wrap;
}

// Plan 34 (Marco): rad/top/soc as FULL representations -- the dim VECTOR per object
// (the redundant total-dim column is gone) plus each arrow's exact action matrix,
// typeset by renderMath (KaTeX). The LaTeX is built from the block's {dims, maps}.

// Matrices are the COMPLETE human record (Plan 34, Marco): shown IN FULL and never
// elided at a small size. MAT_BACKSTOP_CELLS is only a SANITY cap mirroring the
// trace recorder's record-time memory backstop: a pathological/corrupt payload past
// it is stated by shape instead of hanging the browser. One constant per file
// (gui.js has the same one, same comment).
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

// A matrix shown COMPLETE, with NO scrollbar (Marco, 2026-07-29): the box clips
// nothing -- fitMath() shrinks an over-wide matrix to the column width after
// typesetting, so the page body never scrolls sideways either. Same contract as
// gui.js's mathFit.
function mathFit(latex) {
  const box = document.createElement("span");
  box.className = "ql-fit";
  box.style.display = "block";
  box.style.maxWidth = "100%";
  const span = document.createElement("span");
  span.className = "arithmatex";
  span.style.display = "inline-block";
  span.style.transformOrigin = "left top";
  span.textContent = "\\(" + latex + "\\)";
  box.appendChild(span);
  return box;
}

// Shrink-to-fit every .ql-fit box whose typeset content overflows its column.
// Shrink-ONLY (never magnifies); the wrapper height is corrected so a scaled box
// does not overlap its neighbours. Idempotent (it re-measures unscaled each time).
function fitMath(root) {
  const boxes = (root || document).querySelectorAll(".ql-fit");
  for (const box of boxes) {
    const inner = box.firstChild;
    if (!inner) continue;
    inner.style.transform = "";
    const avail = box.parentNode ? box.parentNode.clientWidth : 0;
    const want = inner.scrollWidth;
    if (!avail || !want || want <= avail) { box.style.height = ""; continue; }
    const k = avail / want;
    inner.style.transform = "scale(" + k + ")";
    box.style.height = Math.ceil(inner.offsetHeight * k) + "px";
  }
}

// A matrix as an INDEXED GRID (Marco 2026-07-29): a header row of column indices,
// a header column of row indices, and a light rule between cells, so an entry can
// be read off by position. 1-based, the mathematician's convention.
function matrixGrid(mat) {
  mat = mat || [];
  const ncols = mat.length ? (mat[0] || []).length : 0;
  // Marco 2026-08-03: a zero MAP is stated ("0"), never drawn as a grid of 0s.
  if (!mat.length || !ncols || matIsZero(mat)) {
    const p = document.createElement("p");
    p.className = "arithmatex";
    p.textContent = "\\( 0 \\)";
    return p;
  }
  if (matTooBig(mat)) {                     // sanity backstop only (never normal use)
    const p = document.createElement("p");
    p.className = "pdf-note";
    p.textContent = mat.length + "\u00d7" + ncols + " matrix beyond the display backstop";
    return p;
  }
  const tbl = document.createElement("table");
  tbl.className = "ql-matrix";
  const head = document.createElement("tr");
  const corner = document.createElement("th");
  corner.className = "ql-corner";
  head.appendChild(corner);
  for (let j = 0; j < ncols; j++) {
    const th = document.createElement("th");
    th.textContent = String(j + 1);
    head.appendChild(th);
  }
  tbl.appendChild(head);
  mat.forEach(function (row, i) {
    const tr = document.createElement("tr");
    const th = document.createElement("th");
    th.textContent = String(i + 1);
    tr.appendChild(th);
    (row || []).forEach(function (x) { tr.appendChild(cellText(String(x))); });
    tbl.appendChild(tr);
  });
  return tbl;
}

// An arrow acting as the EXACT zero map carries no information, so its matrix is
// not printed (Marco, 2026-07-29) -- the arrows are named in one line instead.
function matIsZero(mat) {
  return (mat || []).every((row) => (row || []).every((x) => String(x) === "0"));
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
  if (block.series && block.series.length) {   // the Loewy (radical) series (Plan 37)
    const p2 = document.createElement("p");
    p2.textContent = d.modLoewySeries || "Loewy (radical) series, top to bottom:";
    wrap.appendChild(p2);
    wrap.appendChild(loewySeriesTable(block.series));
  }
  for (const pair of trio) {
    appendRepMaps(wrap, pair[0], pair[1], d);
  }
  return wrap;
}

// A Loewy layer as S_v ⊕ S_w^m; the series as a layer|factors table (Plan 37).
function loewyFactors(layer) {
  const parts = [];
  for (const v of Object.keys(layer).sort()) {
    const m = layer[v];
    if (m) parts.push(m === 1 ? "S_" + v : "S_" + v + "^" + m);
  }
  return parts.length ? parts.join(" ⊕ ") : "0";
}

function loewySeriesTable(series) {
  const tbl = tableWithHeader(["layer", "factors"]);
  series.forEach((layer, i) => {
    const tr = document.createElement("tr");
    tr.appendChild(cellText(String(i + 1)));
    tr.appendChild(cellText(loewyFactors(layer)));
    tbl.appendChild(tr);
  });
  return tbl;
}

// One arrow-matrix line per arrow that acts NON-trivially; the zero arrows are
// named in a single trailing line, so "acts as zero" stays distinguishable from
// "not an arrow of the quiver" without printing a zero matrix (Marco 2026-07-29).
function appendRepMaps(wrap, label, view, d) {
  const maps = (view || {}).maps || {};
  const arrows = Object.keys(maps);
  const live = arrows.filter((a) => !matIsZero(maps[a]));
  const zero = arrows.filter((a) => matIsZero(maps[a]));
  if (!live.length) {
    const q = document.createElement("p");
    q.textContent = label + ": " + ((d || {}).modArrowsZero || "every arrow acts as zero");
    wrap.appendChild(q);
    return;
  }
  for (const a of live) {
    const q = document.createElement("p");
    q.textContent = label + ", arrow " + a + ":";
    wrap.appendChild(q);
    wrap.appendChild(matrixGrid(maps[a]));          // indexed grid, shown complete
  }
  if (zero.length) {
    const q = document.createElement("p");
    q.className = "pdf-note";
    q.textContent = label + ": "
      + ((d || {}).modArrowsActingZero || "arrows acting as zero:")
      + " " + zero.join(", ");
    wrap.appendChild(q);
  }
}

// The AR-translate input certificate (Marco #1): indecomposable, or the input's
// decomposition + the additivity note. null when the block carries no certificate
// (the decompose engine was unavailable), so the note never lies.
function tauCertNote(block, d, name) {
  name = name || "M";
  const p = document.createElement("p");
  if (block.indecomposable === true) {
    p.textContent = (d.modIndecomposable || "input M is indecomposable")
      .replace(" M ", " " + name + " ");
    return p;
  }
  if (block.decomposition) {
    const parts = block.decomposition.map(function (s) {
      return dvText(s.dim_vector) + (s.multiplicity > 1 ? "^" + s.multiplicity : "");
    }).join("  ⊕  ");
    p.textContent = name + " ≅ " + parts + " — "
      + (d.modTauAdditive || "τ computed summand-wise (τ is additive)");
    return p;
  }
  return null;
}

// A tau / tau^- block: the dimension-vector line, the translate's FULL per-arrow
// matrices, the input certificate -- and the same again for the second module N
// when the request named one (Marco, 2026-07-29).
// The second module N is either the Ext argument or the Tor argument; name which.
function targetRoleText(role, d) {
  if (role === "tor_target") {
    return (d || {}).modTauTargetTor || "and for N, the Tor target:";
  }
  if (role === "ext_target") {
    return (d || {}).modTauTargetExt || "and for N, the Ext target:";
  }
  return "and for N:";
}
function tauBlock(block, d, kind) {
  const wrap = document.createElement("div");
  appendTranslate(wrap, block, d, kind, "M");
  for (const t of block.targets || []) {
    const p = document.createElement("p");
    p.textContent = targetRoleText(t.role, d);
    wrap.appendChild(p);
    appendTranslate(wrap, t, d, kind, "N");
  }
  return wrap;
}

function appendTranslate(wrap, t, d, kind, name) {
  const sym = (kind === "tau" ? "τ" : "τ⁻") + name;
  if (t.error) {
    wrap.appendChild(errDiv(sym + " " + ((d || {}).modTauUnavailable || "is unavailable")
                            + ": " + t.error));
    return;
  }
  if (t.latex) {
    const eq = document.createElement("p");
    eq.className = "arithmatex";
    eq.textContent = "\\[ " + t.latex + " \\]";
    wrap.appendChild(eq);
  }
  if (t.repr) appendRepMaps(wrap, sym, t.repr, d);
  const note = tauCertNote(t, d, name);
  if (note) wrap.appendChild(note);
}

// The differentials of a projective/injective resolution, as full matrices. A
// differential EQUAL to one already shown references the earlier degree instead of
// repeating the matrix (Marco, 2026-07-29) -- decisive for periodic resolutions.
function differentialsBlock(block, proj, d) {
  const wrap = document.createElement("div");
  const diffs = block.differentials || [];
  if (!diffs.length) return wrap;
  const head = document.createElement("p");
  head.textContent = proj
    ? ((d || {}).modDifferentialsProj
       || "differentials (rows: target basis, columns: source basis; d_0 = ε: Q_0 → M)")
    : ((d || {}).modDifferentialsInj
       || "differentials (rows: target basis, columns: source basis; d^0 = ι: M → E^0)");
  wrap.appendChild(head);
  const seen = new Map();
  diffs.forEach(function (df, n) {
    const label = (proj ? "d_" : "d^") + n;
    const sym = proj ? "d_{" + n + "}" : "d^{" + n + "}";
    const p = document.createElement("p");
    if (df.elided || !df.matrix) {
      p.textContent = label + ": " + df.rows + "×" + df.cols + " — "
        + ((d || {}).modMatrixTooLarge
           || "matrix too large to display; it is complete in the report data");
      wrap.appendChild(p);
      return;
    }
    if (matIsZero(df.matrix)) {
      // Marco 2026-08-03: a zero map is stated, never drawn or cross-referenced.
      p.className = "arithmatex";
      p.textContent = "\\(" + sym + " = 0\\)";
      wrap.appendChild(p);
      return;
    }
    const key = JSON.stringify(df.matrix);
    if (seen.has(key)) {
      p.className = "pdf-note";
      p.textContent = label + " = " + seen.get(key) + " ("
        + ((d || {}).modSameMatrix || "the same matrix as above; not repeated") + ")";
      wrap.appendChild(p);
      return;
    }
    seen.set(key, label);
    p.className = "arithmatex";
    p.textContent = "\\(" + sym + " =\\)";
    wrap.appendChild(p);
    wrap.appendChild(matrixGrid(df.matrix));
  });
  return wrap;
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
      out.appendChild(differentialsBlock(b, kind === "projective_resolution", d));
    } else if (kind === "rad_top_soc") {
      out.appendChild(radTopSocBlock(b, d));
    } else if (kind === "decompose") {
      out.appendChild(decomposeBlock(b, d));
    } else if (kind === "almost_split") {
      out.appendChild(almostSplitBlock(b, d));
    } else if (kind === "tilting_check") {
      out.appendChild(tiltingBlock(b, d));
    } else if (kind === "orbit_geometry") {
      out.appendChild(orbitGeometryBlock(b, d));
    } else if (kind === "tau" || kind === "tau_minus") {
      out.appendChild(tauBlock(b, d, kind));
    } else if (kind === "homological_profile") {
      out.appendChild(homologicalProfileBlock(b, d));
    } else if (kind === "derived_fingerprint") {
      out.appendChild(derivedFingerprintBlock(b, d));
    }
  }
}

// The derived fingerprint (Plan 43) as a labelled table + necessary-condition scope.
// A field captured as an error prints its message; no field is silently dropped.
function derivedFingerprintBlock(block, d) {
  const wrap = document.createElement("div");
  const p = document.createElement("p");
  p.textContent = d.dfTitle || "Derived fingerprint";
  wrap.appendChild(p);
  const fp = block.fingerprint || {};
  const cell = function (v) {
    if (v && typeof v === "object" && "error" in v) return "unavailable — " + v.error;
    if (Array.isArray(v)) return "[" + v.join(", ") + "]";
    return String(v);
  };
  const rows = [
    ["Coxeter polynomial", fp.coxeter_polynomial],
    ["det C", fp.cartan_det],
    ["Cartan Smith factors", fp.cartan_smith],
    ["dim HH^• (cohomology)", fp.hh_cohomology_dims],
    ["dim HH_• (homology)", fp.hh_homology_dims],
    ["dim HC_• (cyclic)", fp.cyclic_dims],
    ["dim Z(A)", fp.center_dim],
    ["global dimension", fp.gl_dim],
  ];
  const tbl = document.createElement("table");
  for (const r of rows) {
    const tr = document.createElement("tr");
    const th = document.createElement("th");
    th.textContent = r[0];
    const td = document.createElement("td");
    td.textContent = cell(r[1]);
    tr.appendChild(th);
    tr.appendChild(td);
    tbl.appendChild(tr);
  }
  wrap.appendChild(tbl);
  const scope = document.createElement("p");
  scope.className = "hint";
  scope.textContent = block.scope || (d.dfScope
    || "a derived-invariant fingerprint; equal values are a necessary condition for derived equivalence, not a proof");
  wrap.appendChild(scope);
  return wrap;
}

// The C6 homological-dimensions family (Plan 40) as a labelled table. Values are
// the block's own honest text fields (exact value / certified bound / infinity /
// undecided / per-entry error); labels are i18n via the form dataset.
function homProfileFindim(f) {
  if (f.exact) return "findim = " + f.lower + "  (" + f.note + ")";
  if (f.upper !== null && f.upper !== undefined)
    return "findim ∈ [" + f.lower + ", " + f.upper + "]  (" + f.note + ")";
  return "findim ≥ " + f.lower + "  (" + f.note + ")";
}

function homologicalProfileBlock(block, d) {
  const wrap = document.createElement("div");
  const p = document.createElement("p");
  p.textContent = d.hpTitle || "Homological dimensions";
  wrap.appendChild(p);
  const rows = [
    [d.hpGldim || "Global dimension", block.global_dimension.text],
    [d.hpFindim || "Finitistic dimension", homProfileFindim(block.finitistic)],
    [d.hpDomdim || "Dominant dimension", block.dominant.text],
    [d.hpGoren || "Gorenstein", block.gorenstein.text],
  ];
  const it = block.igusa_todorov;
  rows.push([(d.hpIgusa || "Igusa–Todorov φ/ψ"),
             it.error ? ("not computed: " + it.error)
                      : ("of " + it.module + ": φ = " + it.phi + ", ψ = " + it.psi)]);
  const tbl = document.createElement("table");
  for (const r of rows) {
    const tr = document.createElement("tr");
    const th = document.createElement("th");
    th.textContent = r[0];
    const td = document.createElement("td");
    td.textContent = r[1];
    tr.appendChild(th);
    tr.appendChild(td);
    tbl.appendChild(tr);
  }
  wrap.appendChild(tbl);
  return wrap;
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
  const lang = ["en", "es", "fr", "zh"].includes(form.dataset.lang)
    ? form.dataset.lang : "en";
  const fam = CATALOG.families.find((f) => f.name === form.elements.family.value);
  box.textContent = "";
  if (summary) summary.textContent =
    (fam && fam.summary) ? (fam.summary[lang] || fam.summary.en || "") : "";
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
    const helpText = p.help && (p.help[lang] || p.help.en);
    if (helpText) {
      const help = document.createElement("div");
      help.className = "param-help";
      help.textContent = helpText;
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
  // The fields the library accepts come from the catalog (CATALOG.fields), NOT a
  // hardcoded list -- a scalar field is one that `needs` no extra input (CC, QQ);
  // a field that needs `p` (GF) reads the prime/degree. A bare "CC"/"QQ" matches
  // its catalog entry; anything else falls through to GF(p).
  const fieldsMeta = (CATALOG && CATALOG.fields) || {};
  const meta = fieldsMeta[fieldRaw];
  const field = (meta && !(meta.needs || []).length)
    ? {kind: fieldRaw}
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
  if (typeof renderMathInElement === "function") {
    renderMathInElement(root || document.body, {
      delimiters: [
        {left: "$$", right: "$$", display: true},
        {left: "\\(", right: "\\)", display: false},
        {left: "$", right: "$", display: false},
      ],
      throwOnError: false,
    });
  }
  fitMath(root);       // measure AFTER typesetting: shrink over-wide matrices
}
renderMath(document.body);
// The fit factor depends on the column width, so re-measure on resize (debounced):
// a widened window must give the matrices their full size back.
let fitTimer = null;
window.addEventListener("resize", function () {
  if (fitTimer) clearTimeout(fitTimer);
  fitTimer = setTimeout(function () { fitMath(); }, 150);
});
