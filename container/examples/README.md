# Container examples: one nontrivial algebra, computed two ways

The worked example is the **quantum complete intersection**

```
Lambda_q = k<x,y> / (x^2, y^2, x*y + 2*y*x),   k = GF(32003)      (QuantumCI q=2)
```

dim 4, non-monomial, noncommutative, **self-injective but not symmetric**
(q != 1, so the Nakayama automorphism is nontrivial), quadratic and **Koszul**.
Small enough to run in seconds, rich enough that every part of the library has
something nontrivial to say about it.

## The two ways

`qci-q2.yaml` is a single Plan-28 spec document. Run it **(1) in the container**
and **(2) with the code** (the wheel CLI, or the library directly):

```bash
# (1) container
mkdir -p out/container && chmod 777 out/container
docker run --rm -v "$PWD/container/examples:/cfg:ro" -v "$PWD/out/container:/out" \
    quiverlab:local run /cfg/qci-q2.yaml -o /out/result.json

# (2) code (wheel CLI in any venv with quiverlab[fast,hpc] installed)
quiverlab-hpc run container/examples/qci-q2.yaml -o out/local/result.json
```

The two `result.json` (and `tikz.tex`) files are **byte-identical** across the
container's Python 3.12 and a host Python 3.11 venv — modulo the per-run
`resources` telemetry footer (wall time / peak RSS), which is timing and thus
never reproducible. Check it like this:

```bash
python -c "
import json
strip = lambda p: json.dumps({k: v for k, v in json.load(open(p)).items()
                              if k != 'resources'}, sort_keys=True)
a, b = strip('out/local/result.json'), strip('out/container/result.json')
assert a == b, 'mathematical payload drifted'
print('payload byte-identical')"
```

Everything mathematical is byte-stable by construction (exact arithmetic,
canonical CS differentials, sorted-key JSON).

`compute_everything.py` is the direct-library companion: it recomputes
everything the YAML computes through `import quiverlab` and adds the surface
the spec does not expose (cyclic homology, complexity, the Frobenius/symmetry
certificates, the Nakayama automorphism, the Yoneda Ext-algebra, Koszulity,
`A^op`, the trivial extension, duality/transpose/Hom, cross-engine HH checks).
It also runs identically in both worlds:

```bash
python container/examples/compute_everything.py                    # code
docker run --rm -v "$PWD:/work:ro" --entrypoint python \
    quiverlab:local /work/container/examples/compute_everything.py # container
```

The two outputs diff clean.

## What comes out (all exact, all cross-checked)

| invariant | value |
|---|---|
| dim, basis | 4; `e_1, x, y, x*y`; Loewy length 3 |
| Cartan / Coxeter poly | `[[4]]`; `t + 1` |
| center | dim 2 = ⟨1, x·y⟩ |
| global dimension | infinite (certified `>= 32`) |
| self-injective / Frobenius | True / True |
| symmetric / weakly symm. | **False** / True (Nakayama `nu(x) = -2x`, `nu(y) = -y/2`) |
| HH^0..8 | `[2, 2, 1, 0, 0, 0, 0, 0, 0]` (Bergh–Erdmann; CS ≡ bar ≡ auto) |
| HH_0..8 | `[3, 2, 2, 2, 2, 2, 2, 2, 2]` |
| HC_0..4 | `[3, 0, 3, 0, 3]` |
| complexity | 2 (= codimension of the quantum CI) |
| Ext-algebra E(Λ) | 2 generators, 1 relation (quantum plane), **Koszul** |
| Ext^n(S,S) = Tor_n(S,S) | `n + 1` — the Koszul-dual graded dimensions |
| proj. resolution of S | term dims `4, 8, 12, …` (Betti `1, 2, 3, …`), pd = ∞ |
| tau(S), tau^-(S) | dim 5, indecomposable; `tau^- tau S ≅ S` |
| decompose (dim-3 test module) | 2-dim indec ⊕ simple (Krull–Schmidt certified) |
| T(A) = A ⋉ DA | dim 8, symmetric True |

Render a human-readable report from either result file:

```bash
docker run --rm -v "$PWD/out/container:/out" quiverlab:local render \
    /out/result.json -o /out/report.html --format html
```
