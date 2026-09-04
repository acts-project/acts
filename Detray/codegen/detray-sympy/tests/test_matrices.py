import pytest
import sympy

import detray_sympy.checks
import detray_sympy.matrices


# Test that the D-matrix getters with gradients turned on give matrices with
# the same known substructure.
def test_generic_D_matrix_with_gradient():
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
    Dg = detray_sympy.matrices.get_generic_matrix_D(gradient=True)
    assert detray_sympy.checks.has_same_known_substructure(D, Dg)


# Test that the D-matrix getters with gradients turned off give matrices with
# the same known substructure.
def test_generic_D_matrix_without_gradient():
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
    Dg = detray_sympy.matrices.get_generic_matrix_D(gradient=False)
    assert detray_sympy.checks.has_same_known_substructure(D, Dg)
