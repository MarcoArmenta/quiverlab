"""Plan 35 -- the HH product surface: structure-constant tables for the cup
product, the cap module action, the Gerstenhaber bracket, and the induced
Connes differential, as frozen result objects with one canonical block
serialization (consumed identically by hpc/spec.py and docs/gui/runner.py).

Constants are ALWAYS exact strings at the boundary (`str(entry)`): ints mod p
on the GF(p) routes, Domain reprs on the CS route. No floats can appear (the
AST gate scans this file)."""
from dataclasses import dataclass


@dataclass(frozen=True)
class ProductTable:
    kind: str            # "cup" | "cap" | "bracket"
    degrees: tuple       # (p, q) -- for cap, (p, n): HH^p (x) HH_n -> HH_{n-p}
    out_degree: int
    dims: tuple          # (dim_left, dim_right, dim_out)
    constants: tuple     # [k][i][j] -> str: left_i * right_j = sum_k c^k_ij out_k

    def as_dict(self):
        return {"degrees": list(self.degrees), "out_degree": self.out_degree,
                "dims": list(self.dims),
                "constants": [[[c for c in row] for row in mat]
                              for mat in self.constants]}


class HHProducts:
    """A family of product tables up to `top`. kind in {"cup","cap","bracket"}."""

    def __init__(self, kind, top, tables, engine, basis, window, references):
        self.kind = kind
        self.top = top
        self.tables = dict(tables)     # {(p, q): ProductTable}
        self.engine = engine
        self.basis = basis
        self.window = window           # int for bracket (served window), else None
        self.references = list(references)

    def blocks(self):
        out = {"kind": self.kind, "top": self.top, "engine": self.engine,
               "basis": self.basis,
               "tables": [self.tables[k].as_dict()
                          for k in sorted(self.tables)],
               "references": list(self.references)}
        if self.window is not None:
            out["window"] = self.window
        return out

    def __repr__(self):
        return (f"<HHProducts {self.kind} top={self.top} "
                f"tables={len(self.tables)} basis={self.basis!r}>")


class ConnesB:
    """Induced Connes differentials B: HH_n -> HH_{n+1} for 0 <= n < top."""

    def __init__(self, top, hh_dims, matrices, ranks, engine, references):
        self.top = top
        self.hh_dims = list(hh_dims)   # dim HH_0..HH_top
        self.matrices = dict(matrices) # {n: rows of str, shape hh_dims[n+1] x hh_dims[n]}
        self.ranks = dict(ranks)       # {n: int}
        self.engine = engine
        self.references = list(references)

    def blocks(self):
        return {"kind": "connes_b", "top": self.top,
                "hh_dims": list(self.hh_dims),
                "matrices": {str(n): self.matrices[n] for n in sorted(self.matrices)},
                "ranks": {str(n): self.ranks[n] for n in sorted(self.ranks)},
                "engine": self.engine, "references": list(self.references)}

    def __repr__(self):
        return f"<ConnesB top={self.top} ranks={self.ranks}>"


def _pairs(kind, top):
    """The degree pairs a `kind` table family covers up to `top`."""
    if kind == "cup":
        return [(p, q) for p in range(top + 1) for q in range(top + 1 - p)]
    if kind == "cap":
        return [(p, n) for n in range(top + 1) for p in range(n + 1)]
    if kind == "bracket":
        return [(p, q) for p in range(1, top + 2) for q in range(1, top + 2 - p + 1)
                if p + q - 1 <= top]
    raise ValueError(f"unknown product kind {kind!r}")


def gfp_product_tables(A, kind, top, max_cells):
    """The GF(p) bar-route table family: tt_calculus structure constants on the
    bar HH basis. A must be over a prime field (the caller routes)."""
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine import tt_calculus as TT
    from quiverlab.engine.scan3 import cochain_basis
    from quiverlab.errors import DepthLimitError
    prime = A.domain.p
    E = to_engine(A.unit_adapted())
    fn = {"cup": TT.cup_product_matrix, "cap": TT.cap_product_matrix,
          "bracket": TT.gerstenhaber_bracket_matrix}[kind]
    out_deg = {"cup": lambda p, q: p + q, "cap": lambda p, n: n - p,
               "bracket": lambda p, q: p + q - 1}[kind]
    tables = {}
    for (p, q) in _pairs(kind, top):
        cells = len(cochain_basis(E, p)) * len(cochain_basis(E, q))
        if cells > max_cells:
            raise DepthLimitError(
                f"{kind} product table ({p}, {q}): the two cochain bases pair "
                f"{cells} cells (> max_cells = {max_cells})",
                hint="raise max_cells or lower top")
        C, dl, dr, dout = fn(E, p, q, prime)
        tables[(p, q)] = ProductTable(
            kind=kind, degrees=(p, q), out_degree=out_deg(p, q),
            dims=(dl, dr, dout),
            constants=tuple(tuple(tuple(str(int(C[k, i, j])) for j in range(dr))
                                  for i in range(dl)) for k in range(dout)))
    window = top if kind == "bracket" else None
    return HHProducts(kind=kind, top=top, tables=tables,
                      engine="hanlab engine (F_p fast rank)",
                      basis=f"bar/GF({prime})", window=window,
                      references=_REFERENCES[kind])


_REFERENCES = {"cup": ["cup", "gerstenhaber"],
               "cap": ["cup", "gerstenhaber"],
               "bracket": ["bracket", "gerstenhaber"],
               "connes_b": ["cyclic"]}
