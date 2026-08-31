import argparse
import contextlib
import sys

import sympy as sym
from sympy import MatrixSymbol

from codegen.sympy_common import (
    NamedExpr,
    name_expr,
    find_by_name,
    cxx_printer,
    my_expression_print,
)

step_path_derivatives = (
    MatrixSymbol("step_path_derivatives", 8, 1).as_explicit().as_mutable()
)
step_path_derivatives[7, 0] = 0  # qop

surface_path_derivatives = (
    MatrixSymbol("surface_path_derivatives", 1, 8).as_explicit().as_mutable()
)
surface_path_derivatives[0, 3] = 0
surface_path_derivatives[0, 7] = 0

# M is the bound-to-free jacobian transported to the current point.  loc0 and
# loc1 stay position only, phi and theta pick up a direction part (position too
# on a line surface), q/p picks up all three, and time stays exactly e_time.
J_bf = MatrixSymbol("M", 8, 6).as_explicit().as_mutable()
tmp = sym.zeros(8, 6)
tmp[0:3, 0:2] = J_bf[0:3, 0:2]
tmp[0:3, 2:5] = J_bf[0:3, 2:5]
tmp[3, 4] = J_bf[3, 4]
tmp[4:7, 2:5] = J_bf[4:7, 2:5]
tmp[7, 4] = J_bf[7, 4]
tmp[3, 5] = 1
J_bf = tmp

J_fb = MatrixSymbol("J_fb", 6, 8).as_explicit().as_mutable()
tmp = sym.zeros(6, 8)
tmp[0:2, 0:3] = J_fb[0:2, 0:3]
tmp[2:4, 4:7] = J_fb[2:4, 4:7]
tmp[5, 3] = 1
tmp[4, 7] = 1
J_fb = tmp


def full_transport_jacobian_generic() -> list[NamedExpr]:
    J_full = name_expr(
        "J_full",
        J_fb * (sym.eye(8) + step_path_derivatives * surface_path_derivatives) * J_bf,
    )

    return [J_full]


def full_transport_jacobian_curvilinear(direction: MatrixSymbol) -> list[NamedExpr]:
    surface_path_derivatives = (
        MatrixSymbol("surface_path_derivatives", 1, 8).as_explicit().as_mutable()
    )
    surface_path_derivatives[0, 0:3] = -direction.as_explicit().transpose()
    surface_path_derivatives[0, 3:8] = sym.zeros(1, 5)

    J_full = name_expr(
        "J_full",
        J_fb * (sym.eye(8) + step_path_derivatives * surface_path_derivatives) * J_bf,
    )

    return [J_full]


def my_full_transport_jacobian_generic_function_print(name_exprs, run_cse=True):
    printer = cxx_printer
    outputs = [find_by_name(name_exprs, name)[0] for name in ["J_full"]]

    lines = []

    head = (
        "template <typename T> void boundToBoundTransportJacobianImpl("
        "std::span<const T, 48> J_fb, std::span<const T, 48> M,"
        " std::span<const T, 8> step_path_derivatives,"
        " std::span<const T, 8> surface_path_derivatives,"
        " std::span<T, 36> J_full) {"
    )
    lines.append(head)

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)


def my_full_transport_jacobian_curvilinear_function_print(name_exprs, run_cse=True):
    printer = cxx_printer
    outputs = [find_by_name(name_exprs, name)[0] for name in ["J_full"]]

    lines = []

    head = (
        "template <typename T> void boundToCurvilinearTransportJacobianImpl("
        "std::span<const T, 48> J_fb, std::span<const T, 48> M,"
        " std::span<const T, 8> step_path_derivatives,"
        " std::span<const T, 3> dir, std::span<T, 36> J_full) {"
    )
    lines.append(head)

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)


def check_curvilinear_is_generic_specialised() -> None:
    """Assert the curvilinear jacobian is the generic one at a curvilinear surface.

    They are printed as two functions, so nothing else keeps them in step. A
    curvilinear surface has path derivatives -direction over the position part
    and zero elsewhere.
    """
    direction = MatrixSymbol("dir", 3, 1)
    generic = full_transport_jacobian_generic()[0].expr
    curvilinear = full_transport_jacobian_curvilinear(direction)[0].expr

    spd = MatrixSymbol("surface_path_derivatives", 1, 8).as_explicit()
    at_curvilinear = {spd[0, i]: -direction.as_explicit()[i, 0] for i in range(3)}
    at_curvilinear.update({spd[0, i]: 0 for i in range(3, 8)})

    diff = sym.expand(generic.subs(at_curvilinear) - curvilinear)
    if any(e != 0 for e in diff):
        bad = [
            (i, j)
            for i in range(diff.rows)
            for j in range(diff.cols)
            if diff[i, j] != 0
        ]
        raise AssertionError(
            f"curvilinear jacobian is not the generic one specialised, at {bad}"
        )


def check_covariance_transport_sparsity() -> None:
    """Assert the shape the covariance transport masks this jacobian down to.

    Live entries outside it would be dropped there silently, so check here.

    Raises AssertionError if the jacobian reaches outside that shape.
    """
    J = sym.expand(full_transport_jacobian_generic()[0].expr)
    qop, time = 4, 5
    # q/p depends on nothing but itself, its own diagonal free; nothing depends
    # on time.
    leaks = [(qop, j) for j in range(6) if j != qop and J[qop, j] != 0]
    leaks += [(i, time) for i in range(6) if J[i, time] != (1 if i == time else 0)]
    if leaks:
        raise AssertionError(
            "bound-to-bound jacobian is not of the shape the covariance "
            f"transport assumes; live entries outside it: {leaks}"
        )


HEADER = """// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Note: This file is generated by generate_sympy_jac.py
//       Do not modify it manually.

#pragma once

#include <cmath>
#include <span>
"""


def main(argv: list[str]) -> None:
    """Generate the transport jacobians.

    @param argv is the command line, argv[0] being the program name
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        help="file to write the generated jacobians to; stdout if omitted",
    )
    parser.add_argument(
        "--no-check",
        action="store_true",
        help="skip the symbolic assertions",
    )
    args = parser.parse_args(argv[1:])

    if not args.no_check:
        # If one of these fires the generated code would be wrong, so they
        # guard the generator rather than a test that might not be run.
        check_curvilinear_is_generic_specialised()
        check_covariance_transport_sparsity()

    with (
        open(args.output, "w") if args.output else contextlib.nullcontext(sys.stdout)
    ) as out:
        out.write(HEADER)
        out.write(
            my_full_transport_jacobian_generic_function_print(
                full_transport_jacobian_generic(), run_cse=True
            )
        )
        out.write("\n")
        out.write(
            my_full_transport_jacobian_curvilinear_function_print(
                full_transport_jacobian_curvilinear(MatrixSymbol("dir", 3, 1)),
                run_cse=True,
            )
        )
        out.write("\n")


if __name__ == "__main__":
    main(sys.argv)
