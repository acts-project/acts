import sympy


def get_matrix_D(dFdr, dGdr, dFdqop, dGdqop, dFdt, dGdt, dqopqop, gradient=True):
    D = sympy.eye(8)

    assert dFdr.shape == (3, 3)
    assert dGdr.shape == (3, 3)
    assert dFdt.shape == (3, 3)
    assert dGdt.shape == (3, 3)
    assert dFdqop.shape == (3, 1)
    assert dGdqop.shape == (3, 1)
    assert type(dqopqop) == sympy.Symbol

    if gradient:
        D[0:3, 0:3] = dFdr
        D[4:7, 0:3] = dGdr

    D[0:3, 7:8] = dFdqop
    D[4:7, 7:8] = dGdqop
    D[0:3, 4:7] = dFdt
    D[4:7, 4:7] = dGdt
    D[7, 7] = dqopqop

    return D


def add_transport_jacobian_substructure(J, gradient=True):
    tmp = sympy.eye(8)

    if gradient:
        tmp[0:3, 0:3] = J[0:3, 0:3]
        tmp[4:7, 0:3] = J[4:7, 0:3]

    tmp[0:3, 7:8] = J[0:3, 7:8]
    tmp[4:7, 7:8] = J[4:7, 7:8]
    tmp[0:3, 4:7] = J[0:3, 4:7]
    tmp[4:7, 4:7] = J[4:7, 4:7]
    tmp[7, 7] = J[7, 7]

    return tmp


def get_generic_matrix_D(gradient=True):
    D = sympy.MatrixSymbol("D", 8, 8).as_explicit().as_mutable()
    return add_transport_jacobian_substructure(D, gradient=gradient)
