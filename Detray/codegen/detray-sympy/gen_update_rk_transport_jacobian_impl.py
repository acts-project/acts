import sys
import sympy

from detray_sympy.common import name_expr, cxx_printer_wo_known
from detray_sympy.output import write_out_file
from detray_sympy.codegen import gen_cxx_code
import detray_sympy.matrices


def gen_code(gradient=True):
    dFdr = sympy.MatrixSymbol("dFdr", 3, 3).as_explicit().as_mutable()
    dGdr = sympy.MatrixSymbol("dGdr", 3, 3).as_explicit().as_mutable()

    dFdt = sympy.MatrixSymbol("dFdt", 3, 3).as_explicit().as_mutable()
    dGdt = sympy.MatrixSymbol("dGdt", 3, 3).as_explicit().as_mutable()

    dFdqop = sympy.MatrixSymbol("dFdqop", 3, 1).as_explicit().as_mutable()
    dGdqop = sympy.MatrixSymbol("dGdqop", 3, 1).as_explicit().as_mutable()

    dqopqop = sympy.Symbol("dqopqop")

    D = detray_sympy.matrices.get_matrix_D(
        dFdr, dGdr, dFdqop, dGdqop, dFdt, dGdt, dqopqop, gradient=gradient
    )

    J_transport = sympy.MatrixSymbol("J_transport", 8, 8).as_explicit().as_mutable()
    J_transport = detray_sympy.matrices.add_transport_jacobian_substructure(
        J_transport, gradient=gradient
    )

    new_transport_jacobian = D * J_transport

    input_name_exprs = []
    input_name_exprs.append(name_expr("J_transport", J_transport))
    input_name_exprs.append(name_expr("dFdt", dFdt))
    input_name_exprs.append(name_expr("dGdt", dGdt))
    if gradient:
        input_name_exprs.append(name_expr("dFdr", dFdr))
        input_name_exprs.append(name_expr("dGdr", dGdr))
    input_name_exprs.append(name_expr("dFdqop", dFdqop))
    input_name_exprs.append(name_expr("dGdqop", dGdqop))
    input_name_exprs.append(name_expr("dqopqop", dqopqop))
    output_name_exprs = [name_expr("new_J", new_transport_jacobian)]
    code = gen_cxx_code(
        "update_transport_jacobian_"
        + ("with_gradient" if gradient else "without_gradient")
        + "_impl",
        input_name_exprs,
        output_name_exprs,
        run_cse=True,
        printer=cxx_printer_wo_known,
    )
    return code


if __name__ == "__main__":
    if len(sys.argv) > 1:
        output = sys.argv[1]
    else:
        output = None

    c1 = gen_code(False)
    c2 = gen_code(True)

    write_out_file(c1 + c2, output)
