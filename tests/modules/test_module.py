"""Right A-modules over any Domain (spec §3.6, §5 component 7). Column-vector
convention: m*b = action[b] @ m; action[x*y] = action[y] @ action[x]."""
from quiverlab import linear_path_algebra, GF, Quiver
from quiverlab.modules.module import Module
from quiverlab.modules import linalg_mod as lm
from quiverlab.fields import QQ


def _ka3_p1():
    """Regular projective P_1 = e_1 A over kA_3 (1 --x--> 2 --y--> 3, no relations, so the
    length-2 path x*y is NONZERO). Basis columns [e_1, x, x*y] (paths starting at 1)."""
    A = Quiver([1, 2, 3], {"x": (1, 2), "y": (2, 3)}).algebra(field=QQ)
    o, z = QQ.coerce(1), QQ.coerce(0)
    action = {
        "e_1": [[o, z, z], [z, z, z], [z, z, z]],
        "e_2": [[z, z, z], [z, o, z], [z, z, z]],
        "e_3": [[z, z, z], [z, z, z], [z, z, o]],
        "x":   [[z, z, z], [o, z, z], [z, z, z]],   # e_1*x = x    : col0 -> col1 (1 -> 2)
        "y":   [[z, z, z], [z, z, z], [z, o, z]],   # x*y = x*y    : col1 -> col2 (2 -> 3)
        "x*y": [[z, z, z], [z, z, z], [o, z, z]],   # e_1*(x*y)=x*y: col0 -> col2 (1 -> 3)
    }
    return A, action


def test_matrix_helpers_over_qq():
    I = lm.identity(2, QQ)
    assert I == [[QQ.coerce(1), QQ.coerce(0)], [QQ.coerce(0), QQ.coerce(1)]]
    A = [[QQ.coerce(1), QQ.coerce(2)]]
    B = [[QQ.coerce(1)], [QQ.coerce(1)]]
    assert lm.matmul(A, B, QQ) == [[QQ.coerce(3)]]
    K = lm.kernel_columns([[QQ.coerce(1), QQ.coerce(-1)]], QQ)
    assert len(K) == 1 and K[0][0] == K[0][1]


def test_column_space_pivots_selects_the_nonzero_column():
    # regression: column_space_pivots must return pivot COLUMNS of M (not of M^T).
    # M = [[0, 0], [1, 0]] has column 0 nonzero, column 1 zero -> pivots == [0].
    o, z = QQ.coerce(1), QQ.coerce(0)
    M = [[z, z], [o, z]]
    assert lm.column_space_pivots(M, QQ) == [0]     # the rad-image basis, not the zero col


def test_regular_projective_p1_of_a2_is_a_module():
    # P_1 = e_1 A for kA_2 built directly (Task 3 builder); here we hand-assemble it to
    # pin the representation. basis of P_1: [e_1, a] (paths from vertex 1).
    A = linear_path_algebra(2)     # arrow auto-named "a1"; basis [e_1, e_2, a1]
    d = QQ
    o, z = d.coerce(1), d.coerce(0)
    # right action on column [e_1, a1]:  e_1*e_1=e_1, e_1*a1=a1, a1*e_2=a1
    action = {
        "e_1": [[o, z], [z, z]],   # projects onto e_1 component
        "e_2": [[z, z], [z, o]],   # projects onto a1 component (a1 = a1*e_2)
        "a1":  [[z, z], [o, z]],   # e_1*a1 = a1 : sends e_1-coord to a1-coord
    }
    P1 = Module(A, 2, action, name="P_1")
    ok, why = P1.check_module()
    assert ok, why
    assert P1.dimension_vector() == {1: 1, 2: 1}


def test_dimension_vector_sums_to_dim_over_gfp():
    A = linear_path_algebra(2)
    p = GF(5)
    o, z = p.coerce(1), p.coerce(0)
    S1 = Module(A, 1, {"e_1": [[o]], "e_2": [[z]], "a1": [[z]]}, name="S_1")
    assert S1.dimension_vector() == {1: 1, 2: 0}
    assert sum(S1.dimension_vector().values()) == S1.dim


def test_action_of_length2_path_is_antihomomorphism_not_homomorphism():
    # SWAP REGRESSION. kA_2 has no length-2 path, so every existing test is blind to a
    # left/right (homo vs anti-homo) swap in _action_of_word. Over kA_3 the length-2 path
    # x*y exercises the composition order. With the correct anti-homomorphism convention
    #     action[x*y] = action[y] @ action[x]
    # we get, on the basis [e_1, x, x*y]:
    #     action[y] @ action[x] = [[0,0,0],[0,0,0],[1,0,0]]   (NONZERO, = stored action[x*y])
    # whereas the swapped HOMOMORPHISM order would give
    #     action[x] @ action[y] = [[0,0,0],[0,0,0],[0,0,0]]   (ALL ZERO).
    # So the equality asserted below holds ONLY for the anti-homomorphism order; flipping
    # _action_of_word to `action[x] @ action[y]` turns the RHS into the zero matrix and the
    # assertion fails. (Confirmed empirically by monkeypatching _action_of_word to the homo
    # order: the composed value became all-zero and no longer equalled action[x*y].)
    A, action = _ka3_p1()
    M = Module(A, 3, action, name="P_1")
    ok, why = M.check_module()
    assert ok, why
    assert M.dimension_vector() == {1: 1, 2: 1, 3: 1}

    composed = M._action_of_word(("x", "y"))
    anti = lm.matmul(action["y"], action["x"], QQ)   # anti-homomorphism: action[y] @ action[x]
    homo = lm.matmul(action["x"], action["y"], QQ)    # the swap: action[x] @ action[y]
    assert composed == anti                            # must be the anti-homomorphism order
    assert composed == action["x*y"]                   # nonzero, matches the stored composite
    assert composed != homo                            # and DIFFERS from the homomorphism order
    assert any(not QQ.is_zero(x) for row in composed for x in row)  # RHS is genuinely nonzero


def test_check_module_relation_branch_rejects_relation_violation():
    # Exercises check_module's relation branch (untested before: kA_2 is hereditary). Over
    # B = kA_3 / (x*y) the relation forces action[y] @ action[x] == 0. We feed a representation
    # whose idempotents and per-arrow grading are all valid, but whose length-2 composite is
    # nonzero -- a genuine right kA_3-module that is NOT a B-module. check_module must reach the
    # relation branch and reject it, naming the offending relation.
    B = Quiver([1, 2, 3], {"x": (1, 2), "y": (2, 3)}).algebra(relations=["x*y"], field=QQ)
    assert repr(B.relations[0]) == "x*y"
    o, z = QQ.coerce(1), QQ.coerce(0)
    action = {
        "e_1": [[o, z, z], [z, z, z], [z, z, z]],
        "e_2": [[z, z, z], [z, o, z], [z, z, z]],
        "e_3": [[z, z, z], [z, z, z], [z, z, o]],
        "x":   [[z, z, z], [o, z, z], [z, z, z]],   # 1 -> 2, graded, valid on its own
        "y":   [[z, z, z], [z, z, z], [z, o, z]],   # 2 -> 3, graded, valid on its own
    }
    M = Module(B, 3, action, name="not-a-B-module")
    ok, why = M.check_module()
    assert not ok
    assert "relation not satisfied" in why and "x*y" in why


def test_check_module_rejects_grading_violation_over_hereditary_algebra():
    # The Fix-1 finding: over a hereditary algebra (no relations) the old check_module blessed
    # any arrow action because it only checked sum(P_v)==I and the (empty) relation set. A
    # grading-violating arrow action must now be caught by the grading branch. Here action[x]
    # sends the vertex-1 slot back to the vertex-1 slot (should land in vertex 2), so
    # action[x] != P_2 @ action[x] @ P_1.
    A, action = _ka3_p1()
    bad = dict(action)
    o, z = QQ.coerce(1), QQ.coerce(0)
    bad["x"] = [[o, z, z], [z, z, z], [z, z, z]]   # 1 -> 1, violates s(x)=1 -> t(x)=2
    del bad["x*y"]                                   # drop the stale composite; test the arrow
    M = Module(A, 3, bad, name="bad-grading")
    ok, why = M.check_module()
    assert not ok
    assert "grading" in why and "x" in why


def test_from_arrow_action_extends_to_all_basis_labels():
    A = linear_path_algebra(2)   # basis labels e_1, e_2, a1
    d = QQ
    o, z = d.coerce(1), d.coerce(0)
    P1 = Module.from_arrow_action(
        A, dimension_vector={1: 1, 2: 1},
        arrow_action={"a1": [[z, z], [o, z]]},
        name="P_1")
    # extension must have filled e_1, e_2, a1 and satisfy the module axioms
    assert set(P1.action) == {"e_1", "e_2", "a1"}
    ok, _ = P1.check_module()
    assert ok
