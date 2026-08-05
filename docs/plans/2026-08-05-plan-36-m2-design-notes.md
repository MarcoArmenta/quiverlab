# Plan 36 — Macaulay2 design-mining note (for P39 / P42)

While building the M2 oracle bridge (Plan 36) we drove Macaulay2's `Complexes`
and `SpectralSequences` packages directly. This note records the three API
shapes worth **mirroring** when P39 (Complexes / chain maps / cones / hyper-Ext)
and P42 (spectral-sequence engine + presets) design their surfaces. It is a
design reference, not a spec — the goal is ergonomic parity with a mature,
battle-tested homological-algebra UI, not a port (M2 is commutative-only and has
no quiver/idempotent type, so none of its algebra is reusable here; see the
Plan-36 honest-scope entry on the verification page).

Sources (Macaulay2 ≥ 1.24 packaged documentation; the versions verified live at
plan time were `AssociativeAlgebras`, `Complexes`, and `SpectralSequences` under
M2 1.26.06):

- `Complexes` — https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/Complexes/html/
- `SpectralSequences` — https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/SpectralSequences/html/
- `AssociativeAlgebras` — https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/AssociativeAlgebras/html/
  (flagged "interface will change" upstream — hence the version-pin policy in the
  Plan-36 bridge).

---

## 1. `Complexes` ergonomics to mirror (feeds P39)

A `Complex` in M2 exposes three things P39's complex type should offer under the
same shape:

- **A single indexable differential object.** `C.dd` is one `ComplexMap` (verified
  `class C.dd == ComplexMap`); the degree-`i` differential is read as `dd^C_i`
  (returns the `Matrix`). Contrast our current per-degree `differential(n)` calls:
  a single `C.dd` handle that you *index* (rather than a family of methods) reads
  cleaner and composes (`C.dd * C.dd == 0` is one expression). P39 should expose
  the differential as one indexable object, not N separate accessors.
- **Terms as `C_i`.** `C_i` is the degree-`i` term (verified `C_2` prints `R^3`).
  Subscript access to terms + `dd^C_i` for maps is the whole navigation surface —
  small and memorable. Map ours to `C[i]` / `C.d[i]` or the Pythonic `C.term(i)` /
  `C.diff(i)`, but keep the *pairing* (term `i`, differential into term `i-1`).
- **`LengthLimit => n` as first-class truncation.** `freeResolution(M, LengthLimit
  => 5)` computes exactly the first 5 terms — truncation is a *build parameter*,
  not a post-hoc slice, and over Artinian rings it is mandatory (the resolution is
  infinite). This is **exactly our `top`**: P39/P42 should thread a single
  `top`/`length` bound through every complex constructor (resolution, cone, Hom
  double complex), refusing loudly when a caller forgets it on a
  non-perfect input — the same discipline the Plan-36 bridge already applies
  (`LengthLimit => {top}` in `commutative_ext_script`).

## 2. The ASCII `betti` grid as the resolution-table display target (feeds P39)

`betti C` prints a compact grid that is the de-facto standard for reading a
resolution's shape. Verified output for `k[x,y]/(x²,y²)`'s residue-field
resolution:

```
       0 1 2 3 4 5
total: 1 2 3 4 5 6
    0: 1 2 3 4 5 6
```

The conventions to copy in the worked-steps resolution tables:

- **Columns = homological degree** (`0 1 2 3 …`), the axis our resolution tables
  already order by.
- **A `total:` row** summing each column — the graded Betti numbers `dim C_n`,
  which is precisely the number our M2 crosscheck compares (`rank C_n`).
- **Internal-degree rows** (`0:`, `1:`, …) splitting each total by internal degree,
  with a **`.` in every zero slot** so the grid stays sparse and scannable (this
  example has a single internal-degree row, so no dots appear; a
  `k[x,y]/(x²,y³)`-style example fills several rows with dots between them).

Our resolution tables already render `C_n = ⊕ P(v,w)` term-by-term; adding the
`betti`-style total/internal-degree grid (columns = degree, dotted zeros) as an
*at-a-glance* header would match what every commutative algebraist reads first.
P39's resolution rendering should offer this grid as the summary line above the
per-term detail.

## 3. The `SpectralSequences` page API as the SS-surface shape (feeds P42)

M2's `SpectralSequences` package is the closest prior art to P42's engine, and
its surface is the one to mirror:

- **`E^r` is a page object; `E^r_{p,q}` returns the actual module** at position
  `(p,q)` (not just its dimension) — so a caller can inspect representatives, not
  only Betti numbers. P42's `reduce_mod_nullspace`-canonical page representatives
  (byte-reproducible, CS precedent) should be reachable the same way:
  `page(r)[p, q]` returns the subquotient with its chosen basis, `page(r).dim(p,q)`
  the number.
- **`spots`** enumerates the non-zero positions of a page — the support, so display
  and convergence code iterate only over live cells. P42 should expose the same
  (`page(r).spots`), which also makes the "single row/column ⇒ E₂ collapse"
  shortcut a one-line check.
- **`netPage`** prints the page as a 2-D grid (p across, q up) with the module
  dimensions in the cells. This is the display target for the worked-steps
  page-grid + differential-matrix rendering P42 already plans (metaplan P42:
  "page-grid display, M2 `netPage` style"). Copy the grid layout (p→columns,
  q→rows, `.`/blank for empty cells, differential arrows `d_r` annotated between
  cells) for the report renderer.

Net: P42's public shape should be `SS(filtered_or_double_complex)` → a page
sequence where each page is indexable by `(p,q)` returning a canonical module,
carries a `spots` support set, and renders as a `netPage`-style grid — with the
convergence report (E_∞ at page `max(width,height)+1`, degeneration decidable by
rank) computed off `fields.linalg` exactly as the metaplan §2 spectral-sequence
brief specifies.

---

*Referenced by the metaplan v0.2.0 cards P39 and P42 (see
`docs/plans/2026-08-05-metaplan-v0.2.0.md`, "Ergonomics follow the P36 mining
note").*
