import sys
import sympy

from detray_sympy.common import name_expr
from detray_sympy.output import write_out_file
from detray_sympy.codegen import gen_cxx_code

if __name__ == "__main__":
    if len(sys.argv) > 1:
        output = sys.argv[1]
    else:
        output = None

    C = sympy.MatrixSymbol("C", 6, 6).as_explicit().as_mutable()

    J_full = sympy.MatrixSymbol("J_full", 6, 6).as_explicit().as_mutable()
    tmp = sympy.eye(6)
    tmp[0:4, 0:5] = J_full[0:4, 0:5]
    tmp[5:6, 0:5] = J_full[5:6, 0:5]
    J_full = tmp

    input_name_exprs = [name_expr("C", C), name_expr("J_full", J_full)]
    output_name_exprs = [name_expr("new_C", J_full * C * J_full.T)]
    code = gen_cxx_code(
        "transport_covariance_to_bound_impl",
        input_name_exprs,
        output_name_exprs,
        run_cse=True,
    )

    write_out_file(code, output)
