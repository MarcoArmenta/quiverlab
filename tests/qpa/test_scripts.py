"""The QPA GAP-source builders are PURE string generation from an Algebra's
public presentation -- no GAP needed to exercise them. These `fast` tests run
everywhere (they escape the `qpa` bucket via the explicit marker) and pin the
Python-side correctness that the GAP-only path would otherwise only hit in CI:
reading `algebra.quiver`/`algebra.relations` (public attrs, NOT `_quiver`/
`_relations`) and `domain.characteristic` (an int attribute, NOT a method)."""
import pytest

from quiverlab import Quiver, GF
from quiverlab.errors import QuiverlabError
from quiverlab.qpa.scripts import quiver_and_algebra_script

pytestmark = pytest.mark.fast


def test_ka2_gf2_script_no_relations():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=GF(2))
    src = quiver_and_algebra_script(A)
    assert 'Quiver(2, [[1, 2, "a"]])' in src
    assert "PathAlgebra(GF(2), Q)" in src
    assert "A := kQ;;" in src            # no relations -> A is the free path algebra


def test_gfp_monomial_script_includes_relations():
    # a relational presentation must reach the `rels` branch -- the `_relations`
    # bug silently dropped it (getattr(algebra, "_relations", None) is always None).
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x*x"], field=GF(3))
    src = quiver_and_algebra_script(A)
    assert "PathAlgebra(GF(3), Q)" in src
    assert "rels :=" in src and "A := kQ/rels;;" in src


def test_builder_reads_public_attrs_not_private():
    # regression guard: a swap back to algebra._quiver / _relations / .characteristic()
    # would AttributeError/TypeError here, since Algebra has no _quiver/_relations attr
    # and Domain.characteristic is an int, not a callable.
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=GF(5))
    quiver_and_algebra_script(A)  # must not raise


def test_module_decl_shape_gfp_and_qq():
    # Plan 23: RightModuleOverPathAlgebra emitter (graded form -> QPA module). Runs
    # without GAP; pins the row-convention matrix + GF(p) element syntax.
    from quiverlab.fields import QQ
    from quiverlab.modules.qpa_module import graded_form
    from quiverlab.qpa.scripts import module_decl
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=GF(5))
    dv, arr = graded_form(A.projective(1))            # P_1 = [1,1], a acts nonzero
    line = module_decl(A, dv, arr, "M")
    assert line.startswith("M := RightModuleOverPathAlgebra(A, [1, 1], [")
    assert '"a"' in line and "One(GF(5))" in line
    Aq = Quiver([1, 2], {"a": (1, 2)}).algebra(field=QQ)
    dvq, arrq = graded_form(Aq.projective(1))
    lineq = module_decl(Aq, dvq, arrq, "M")
    assert "One(" not in lineq                          # QQ entries are bare rationals


@pytest.mark.parametrize("bad", ["in", "mod", "not", "rec", "do", "od", "fi",
                                 "then", "else", "function", "return", "local"])
def test_gap_reserved_word_arrow_name_refused_cleanly(bad):
    # A GAP reserved word used as an arrow name renders `kQ.<name>` into GAP source
    # unescaped -> a confusing GAP syntax error at crosscheck. The bridge now refuses it
    # LOUDLY (typed QuiverlabError naming the arrow) BEFORE emitting any GAP source.
    # No GAP needed: pure string-builder guard.
    A = Quiver([1], {bad: (1, 1)}).algebra(relations=[f"{bad}*{bad}"], field=GF(3))
    with pytest.raises(QuiverlabError, match=bad):
        quiver_and_algebra_script(A)


def test_validate_arrow_name_guard_branches():
    # Direct unit of the guard. Quiver already enforces the identifier rule
    # (r"[A-Za-z_][A-Za-z0-9_]*") at construction, so a non-identifier name cannot reach
    # the bridge through the public API -- but the guard defends that branch too (a name
    # built by a lower-level path). Reserved-word rejection is the reachable, value-adding
    # case; a keyword IS a valid identifier so Quiver lets it through but GAP would not.
    from quiverlab.qpa.scripts import _validate_arrow_name
    for bad in ("in", "mod", "not", "rec", "function", "return"):     # GAP keywords
        with pytest.raises(QuiverlabError, match="reserved word"):
            _validate_arrow_name(bad)
    for bad in ("a-b", "1x", "a.b", "a b", ""):                       # non-identifiers
        with pytest.raises(QuiverlabError, match="valid GAP identifier"):
            _validate_arrow_name(bad)
    for ok in ("a", "x", "a1", "alpha", "f_2", "In", "Mod"):          # valid + non-keyword
        _validate_arrow_name(ok)                                       # must not raise


def test_valid_arrow_names_pass_through_unchanged():
    # The guard must NOT touch valid names: single letters, digits, underscores all
    # emit exactly as before (byte-identical crosscheck source).
    for name in ("a", "x", "a1", "alpha", "f_2", "In", "Mod"):     # "In"/"Mod" != keywords
        A = Quiver([1, 2], {name: (1, 2)}).algebra(field=GF(5))
        src = quiver_and_algebra_script(A)
        assert f'[1, 2, "{name}"]' in src
