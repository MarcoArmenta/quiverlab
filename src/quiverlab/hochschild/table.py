"""Human-readable dimension tables (spec §3.5)."""


class HHTable:
    def __init__(self, dims, kind, algebra_repr, engine="normalized bar complex",
                 references=()):
        self.dims = list(dims)
        self.kind = kind
        self.top = len(dims) - 1
        self.algebra_repr = algebra_repr
        self.engine = engine
        self.references = tuple(references)

    def __getitem__(self, n):
        return self.dims[n]

    def __iter__(self):
        return iter(self.dims)

    def __eq__(self, other):
        if isinstance(other, HHTable):
            return self.kind == other.kind and self.dims == other.dims
        return NotImplemented

    def __repr__(self):
        head = f"{self.kind}n dimensions for {self.algebra_repr} (engine: {self.engine})"
        cells = "  ".join(f"{self.kind}{n} = {d}" for n, d in enumerate(self.dims))
        out = head + "\n" + cells
        if self.references:
            out += "\nrefs: " + ", ".join(self.references)
        return out
