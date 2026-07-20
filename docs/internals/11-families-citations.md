# 11 — Families and citations

## The mathematics
A *family* is a recipe that turns a few numbers or a diagram name into a
finite-dimensional algebra `kQ/I`. quiverlab ships the standard catalogue of
representation theory: Nakayama (serial) algebras by Kupisch series, hereditary
path algebras of Dynkin/Euclidean quivers, truncations `kQ/rad^r`, incidence
algebras of posets, quantum complete intersections, exterior algebras,
preprojective algebras, and two constructions that build a new algebra from old
ones (tensor product and trivial extension). A *citation registry* records, for
every algorithm and every family, the paper it comes from.

## How it is represented
Each family is a plain Python function returning a Plan-03 `Algebra`. There are
three construction *routes*:
- **monomial** — the relations are single forbidden paths; the algebra is built
  by the Plan-01 monomial route (`Quiver.algebra`, the forbidden-word automaton).
  Nakayama, PathAlgebra, TruncatedPathAlgebra, RadicalSquareZero.
- **general** — at least one relation is a genuine linear combination (e.g.
  `x*y + q*y*x`); the algebra is completed by the Plan-03 Gröbner engine.
  QuantumCI, ExteriorAlgebra, PreprojectiveAlgebra, IncidenceAlgebra.
- **structure-constant** — the multiplication table `T` is written directly from
  the factors, with no quiver. TensorProduct, TrivialExtension.
A family stamps the algebra with `_family_citations`, a tuple of registry keys.
The registry itself is a dict `key -> Reference(key, bibtex_key, kind, title,
annotation, tags)`; the annotations are the ground truth the web /literature page
and `quiverlab.bibliography()` render. `references.bib` holds the verified BibTeX.

## How the computation runs
1. `NakayamaAlgebra([3,2,2])` reads the Kupisch series, decides linear vs cyclic
   (`min >= 2` -> cyclic), lays out the quiver, generates the length-`c_i`
   forbidden path from each vertex, and calls the monomial route. `dim = sum c_i`.
2. `QuantumCI(q="i")` writes the three relation strings, `x*y + i*y*x` among
   them; the relation parser accepts the exact token `i` (this chapter's one new
   grammar rule); the Gröbner engine rewrites `y*x -> i*x*y` and certifies dim 4.
3. `TensorProduct(A, B)` fills `T[i*db+j][k*db+l]` with the outer product of the
   two multiplication tables; `dim = dim A * dim B`.
4. `A.hochschild_cohomology(n)` attaches, to the returned `HHTable`, the citation
   keys of the engine paths it used (`.references`): `bar` always, plus `bardzell`
   or `chouhy_solotar` by dispatch. `A.citations()` unions those with the family
   keys. `bibliography(keys)` groups them by kind and prints the annotations.
5. `zoo(dim_max)` loads a bundled JSON catalogue (lifted from hanlab's open_zoo),
   rebuilds each confluent reduction system into an `Algebra`, and yields those
   with `dim <= dim_max`.

## A worked micro-example
`NakayamaAlgebra([3,2,2])`: cyclic `Z_3`, arrows `a1:1->2, a2:2->3, a3:3->1`,
forbidden paths `a1*a2*a3` (len 3 from 1), `a2*a3` (len 2 from 2), `a3*a1` (len 2
from 3). Irreducible basis `e_1,e_2,e_3,a1,a2,a3,a1*a2` -> **dim 7**. Cartan
`[[1,1,1],[0,1,1],[1,0,1]]` (row `i` = composition factors of `P_i`), `det 1`,
`sum = 7`. Centre = scalars, so `HH^0 = 1`. `A.citations()` -> `('nakayama',
'assem_book', 'bar')`; `bibliography(A.citations())` prints the ASS textbook and
Happel's bar-complex reference with their annotations.

## Where to look in the code
| concept | file | function/class |
|---|---|---|
| Kupisch series | `families/nakayama.py` | `NakayamaAlgebra` |
| Dynkin diagram -> quiver | `families/dynkin.py` | `dynkin_quiver` |
| hereditary path algebra | `families/path_algebra.py` | `PathAlgebra` |
| `kQ/rad^r` | `families/truncated.py` | `TruncatedPathAlgebra` |
| poset -> incidence algebra | `families/poset.py`, `families/incidence.py` | `Poset`, `IncidenceAlgebra` |
| quantum CI / exterior / preprojective | `families/{quantum,exterior,preprojective}.py` | `QuantumCI`, `ExteriorAlgebra`, `PreprojectiveAlgebra` |
| tensor / trivial extension | `families/{tensor,trivial_extension}.py` | `TensorProduct`, `TrivialExtension` |
| discoverability | `families/discover.py` | `families`, `CATALOG` |
| curated zoo | `families/zoo.py`, `families/zoo_catalog.json` | `zoo`, `build_from_record` |
| citation registry | `citations/registry.py`, `citations/references.bib` | `REGISTRY`, `reference`, `bibtex` |
| bibliography | `citations/bibliography.py` | `bibliography`, `Bibliography` |
| result references | `hochschild/table.py`, `core/algebra.py` | `HHTable.references`, `Algebra.citations` |
| batch persistence | `batch/db.py`, `batch/scan.py` | `ResultsDB`, `analyze`, `run_scan` |

This chapter is the **Plan 06** checkout. `sweep` (Plan 05) consumes this catalogue
(`families()` is the hook); the trace subsystem (Plan 07) will render the citation
keys these functions stamp.

<!-- Chapter-number note: on disk chapters 01-08 exist (08 = Groebner, Plan 03).
The Plan-04 renumber proposal reserves 09 for the Chouhy-Solotar resolution and 10
for the Plan-05 modules chapter; this families+citations chapter therefore takes the
next free index, 11, per the Plan-06 dispatch. Plan 06 owns the families + citations
chapter content regardless of its final number; if the 09/10 chapters land in a
different order at integration, renumber this to the next free index. -->
