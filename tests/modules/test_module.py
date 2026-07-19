"""Right A-modules over any Domain (spec §3.6, §5 component 7). Column-vector
convention: m*b = action[b] @ m; action[x*y] = action[y] @ action[x]."""
from quiverlab import linear_path_algebra, GF
from quiverlab.modules.module import Module
from quiverlab.modules import linalg_mod as lm
from quiverlab.fields import QQ


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
