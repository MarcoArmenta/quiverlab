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
