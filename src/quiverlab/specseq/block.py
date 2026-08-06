"""The shared ``ss_hochschild`` no-code block builder (Plan 42, Task 9).

Both compute runners -- the server/HPC ``quiverlab.hpc.spec`` and the Pyodide twin
``docs/gui/runner.py`` -- call THIS builder, so the block is byte-identical across
tiers (the ``ext_algebra_block`` / ``recognizers_block`` precedent). It returns
``references`` (not ``citations``); each runner resolves the citation pairs the
same way, exactly like the other Plan-38 blocks."""


def specseq_block(A, top, max_cells=4_000_000):
    """The Hochschild ``(b, B)`` spectral-sequence block: the ``E_inf`` page support,
    the ``netPage`` grid, the abutment (== ``HC_*`` dims), the degeneration data and
    the convergence prose.

    The exponential bar basis is guarded: a ``DepthLimitError`` is caught and
    reported as ``{"kind": "ss_hochschild", "error": ...}`` (the Plan-30 per-entry
    honesty precedent), never a 500. ``references`` is always present so the
    per-runner ``block["citations"] = _citation_pairs(block["references"])`` wiring is
    uniform for both the success and error shapes."""
    from quiverlab.errors import DepthLimitError
    from quiverlab.specseq.presets import hochschild_bB_ss
    refs = ["cyclic", "weibel_homological"]
    try:
        ss = hochschild_bB_ss(A, top, max_cells=max_cells)
    except DepthLimitError as exc:
        return {"kind": "ss_hochschild", "top": top, "error": str(exc),
                "references": refs}
    conv = ss.convergence
    Einf = ss.page(conv.e_infinity_page)
    # The reported abutment + grid are the CERTIFIED HC range 0..top only; the
    # bicomplex is built one degree deeper (top+1) so D_{top+1} pins rank(HC_top),
    # but that boundary degree's homology is a truncation artifact, not HC_{top+1},
    # so it is trimmed here (the self-cert still runs on the full construction).
    einf = [[p, q, Einf.dim(p, q)] for (p, q) in Einf.spots if p + q <= top]
    abutment = [conv.abutment.get(n, 0) for n in range(top + 1)]
    if conv.degenerates_at == 1:
        deg = "it degenerates at E_1 (E_1 = E_inf; the filtration is trivial)"
    elif conv.degenerates_at is not None:
        deg = (f"it degenerates at E_{conv.degenerates_at} "
               "(all higher differentials vanish)")
    else:
        deg = (f"it stabilizes at the E_{conv.e_infinity_page} page "
               "(the bounded-filtration bound)")
    abut_str = ", ".join(f"HC_{n}={abutment[n]}" for n in range(top + 1))
    prose = ("The Hochschild (b, B) spectral sequence converges to cyclic homology "
             f"({abut_str}); {deg}.")
    return {"kind": "ss_hochschild", "top": top,
            "einf": einf,
            "grid": _trim_grid(Einf, top),
            "abutment": abutment,
            "degenerates_at": conv.degenerates_at,
            "collapse": conv.collapse(),
            "prose": prose,
            "references": refs}


def _trim_grid(page, top):
    """The ``netPage`` grid restricted to cells of total degree ``p + q <= top``
    (the certified HC range) -- the same layout as ``Page.grid`` but with the
    boundary degree ``top+1`` (a truncation artifact) dropped."""
    cells = {pq: sq for pq, sq in page._cells.items()
             if pq[0] + pq[1] <= top and sq.dim > 0}
    if not cells:
        return "```\nE_inf: (empty)\n```"
    ps = sorted({p for (p, q) in cells})
    qs = sorted({q for (p, q) in cells})
    w = max(2, max(len(str(sq.dim)) for sq in cells.values()))
    lines = ["E_inf page  (p across, q up):"]
    for q in reversed(qs):
        row = [f"q={q:>2} |"]
        for p in ps:
            sq = cells.get((p, q))
            row.append((str(sq.dim) if sq else ".").rjust(w))
        lines.append(" ".join(row))
    lines.append("      " + " ".join(f"p={p}".rjust(w) for p in ps))
    return "```\n" + "\n".join(lines) + "\n```"


# --------------------------------------------------------------------------- #
# radical_filtration_ss: the associated-graded / radical-filtration spectral
# sequence as an ALGEBRA-only no-code block (wave 2). It runs the P42
# ``radical_filtration_ss`` preset on the minimal projective resolution of the
# direct sum of all simple modules -- the canonical, vertex-symmetric algebra
# object whose radical filtration encodes the Koszul staircase -- to length
# ``top``. Same shared-builder discipline as ``specseq_block`` (byte-identical
# server / Pyodide twins); ``references`` present in BOTH the success and error
# shapes so the ``block["citations"] = _citation_pairs(...)`` wiring is uniform.
# --------------------------------------------------------------------------- #

_RADICAL_REFS = ["weibel_homological", "priddy", "froberg_koszul"]


def radical_filtration_ss_block(A, top, budget_dim=200000):
    """The ``radical_filtration_ss`` block for algebra ``A`` through resolution
    length ``top``: the ``E_inf`` page support, the netPage grid, the abutment
    (= the homology of the resolution complex per degree), the degeneration data
    and the convergence prose -- mirroring :func:`specseq_block`.

    The object filtered is the minimal projective resolution ``Q_top -> ... -> Q_0``
    of ``(+)_v S_v`` (the direct sum of the simple modules), so the sequence is an
    algebra invariant that needs no module input. A presentation-less algebra (no
    quiver / no simple modules) or a resolution that blows up is reported as
    ``{"kind": "radical_filtration_ss", "error": ...}`` -- the ``specseq_block``
    honesty precedent, never a crash."""
    from quiverlab.errors import DepthLimitError, QuiverlabError
    from quiverlab.modules.complexes import ChainComplex
    from quiverlab.modules.morphism import direct_sum
    from quiverlab.specseq.presets import radical_filtration_ss
    refs = list(_RADICAL_REFS)
    if getattr(A, "quiver", None) is None:
        return {"kind": "radical_filtration_ss", "top": top,
                "error": "radical_filtration_ss needs a quiver-presented algebra "
                         "(no simple modules for a structure-constant algebra)",
                "references": refs}
    try:
        simples = [A.simple(v) for v in A.quiver.vertices]
        M = direct_sum(*simples)[0] if len(simples) > 1 else simples[0]
        X = ChainComplex.from_projective_resolution(M, top)
        ss = radical_filtration_ss(X)
    except (DepthLimitError, QuiverlabError) as exc:
        return {"kind": "radical_filtration_ss", "top": top, "error": str(exc),
                "references": refs}
    conv = ss.convergence
    Einf = ss.page(conv.e_infinity_page)
    einf = [[p, q, Einf.dim(p, q)] for (p, q) in Einf.spots if Einf.dim(p, q) > 0]
    abutment = [int(conv.abutment.get(n, 0)) for n in range(top + 1)]
    prose = ("The radical-filtration (associated-graded) spectral sequence of the "
             "minimal projective resolution of the sum of the simple modules. "
             + conv.prose())
    return {"kind": "radical_filtration_ss", "top": top,
            "resolved": "minimal projective resolution of (+)_v S(v)",
            "einf": einf,
            "grid": Einf.grid(),
            "abutment": abutment,
            "degenerates_at": conv.degenerates_at,
            "collapse": conv.collapse(),
            "prose": prose,
            "references": refs}
