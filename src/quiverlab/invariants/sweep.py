"""sweep(): the 'all the moving variables' table -- an algebra rebuilt over several
fields, with chosen invariants tabulated per field (spec section 3.9). The character-
dependence view. A cell whose invariant is unavailable over a field (e.g. an engine-
backed GF(p)-only invariant over CC) records the reason instead of crashing."""


def _default_invariants():
    return {
        "dimension": lambda A: A.dim,
        "loewy_length": lambda A: A.loewy_length(),
        "cartan_det": lambda A: _det(A.cartan_matrix()),
        "coxeter_polynomial": lambda A: A.coxeter_polynomial().as_expr(),
    }


def _det(mat):
    import sympy
    return int(sympy.Matrix(mat).det())


class SweepTable:
    def __init__(self, fields, invariants, data):
        self.fields = fields
        self.invariant_names = list(invariants)
        self._data = data                       # {(inv_name, field_name): value}

    def cell(self, inv, field):
        return self._data[(inv, getattr(field, "name", str(field)))]

    @property
    def rows(self):
        return [[inv] + [self._data[(inv, getattr(f, "name", str(f)))] for f in self.fields]
                for inv in self.invariant_names]

    def __repr__(self):
        head = ["invariant"] + [getattr(f, "name", str(f)) for f in self.fields]
        lines = [head] + [[str(x) for x in row] for row in self.rows]
        widths = [max(len(lines[r][c]) for r in range(len(lines))) for c in range(len(head))]
        out = []
        for r, line in enumerate(lines):
            out.append("  ".join(cell.ljust(widths[c]) for c, cell in enumerate(line)))
            if r == 0:
                out.append("  ".join("-" * widths[c] for c in range(len(head))))
        return "\n".join(out)

    def latex(self):
        cols = "l" + "c" * len(self.fields)
        head = " & ".join(["invariant"] + [getattr(f, "name", str(f)) for f in self.fields])
        body = " \\\\\n".join(
            " & ".join([row[0]] + [str(x) for x in row[1:]]) for row in self.rows)
        return f"\\begin{{tabular}}{{{cols}}}\n{head} \\\\\\hline\n{body}\n\\end{{tabular}}"


def sweep(builder, *args, fields, invariants=None, **kwargs):
    invs = invariants if invariants is not None else _default_invariants()
    data = {}
    for f in fields:
        A = builder(*args, field=f, **kwargs)
        fname = getattr(f, "name", str(f))
        for name, fn in invs.items():
            try:
                data[(name, fname)] = fn(A)
            except Exception as exc:              # record, do not crash the sweep
                data[(name, fname)] = f"n/a: {type(exc).__name__}"
    return SweepTable(fields, invs, data)
