import pytest
import sympy

import detray_sympy.checks
import detray_sympy.matrices


# Check that matrix multiplication with matrix D, computed with the gradient,
# is idempotent in terms of the shape of the output. This is important because
# the transport Jacobian matrix $J_T$ is only ever updated in two ways:
#
# 1. $J_T = I$
# 2. $J_T = D J_T$
#
# This gives us a natural induction that the Jacobian matrix is always a
# product of D-shaped matrices, and if the shape of $D_1 D_2$ is the same as
# the shape of $D_1$, then the transport Jacobian matrix always has the exact
# same known substructure.
def test_shape_idempotence_D_with_gradient():
    dFdr = sympy.MatrixSymbol("dFdr", 3, 3).as_explicit().as_mutable()
    dGdr = sympy.MatrixSymbol("dGdr", 3, 3).as_explicit().as_mutable()
    dFdt = sympy.MatrixSymbol("dFdt", 3, 3).as_explicit().as_mutable()
    dGdt = sympy.MatrixSymbol("dGdt", 3, 3).as_explicit().as_mutable()
    dFdqop = sympy.MatrixSymbol("dFdqop", 3, 1).as_explicit().as_mutable()
    dGdqop = sympy.MatrixSymbol("dGdqop", 3, 1).as_explicit().as_mutable()
    dqopqop = sympy.Symbol("dqopqop")

    D = detray_sympy.matrices.get_matrix_D(
        dFdr, dGdr, dFdqop, dGdqop, dFdt, dGdt, dqopqop, gradient=True
    )
    assert detray_sympy.checks.has_same_known_substructure(D, D * D)
    assert detray_sympy.checks.has_same_known_substructure(D, D * D * D)
    assert detray_sympy.checks.has_same_known_substructure(D, D * D * D * D)


# Test that multiplication with D is idempotent even with the gradient
# disabled.
def test_shape_idempotence_D_without_gradient():
    dFdr = sympy.MatrixSymbol("dFdr", 3, 3).as_explicit().as_mutable()
    dGdr = sympy.MatrixSymbol("dGdr", 3, 3).as_explicit().as_mutable()
    dFdt = sympy.MatrixSymbol("dFdt", 3, 3).as_explicit().as_mutable()
    dGdt = sympy.MatrixSymbol("dGdt", 3, 3).as_explicit().as_mutable()
    dFdqop = sympy.MatrixSymbol("dFdqop", 3, 1).as_explicit().as_mutable()
    dGdqop = sympy.MatrixSymbol("dGdqop", 3, 1).as_explicit().as_mutable()
    dqopqop = sympy.Symbol("dqopqop")

    D = detray_sympy.matrices.get_matrix_D(
        dFdr, dGdr, dFdqop, dGdqop, dFdt, dGdt, dqopqop, gradient=False
    )
    assert detray_sympy.checks.has_same_known_substructure(D, D * D)
    assert detray_sympy.checks.has_same_known_substructure(D, D * D * D)
    assert detray_sympy.checks.has_same_known_substructure(D, D * D * D * D)
