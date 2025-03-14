import sys

import numpy as np

import sympy as sym
from sympy import Symbol, Matrix, ImmutableMatrix, MatrixSymbol
from sympy.codegen.ast import Assignment

from codegen.sympy_common import (
    make_vector,
    NamedExpr,
    name_expr,
    find_by_name,
    my_subs,
    cxx_printer,
    my_expression_print,
)


output = sys.stdout
if len(sys.argv) > 1:
    output = open(sys.argv[1], "w")


# q/p
l = Symbol("l", real=True)

# step length
h = Symbol("h", real=True)

# path length
s = Symbol("s", real=True)

# position
p = make_vector("p", 3, real=True)

# direction
d = make_vector("d", 3, real=True)

# time
t = Symbol("t", real=True)

# mass
m = Symbol("m", real=True, positive=True)

# absolute momentum
p_abs = Symbol("p_abs", real=True, positive=True)

# energy loss per distance
g = Symbol("g", real=True)

# charge
q = Symbol("q", real=True)

# magnetic field
B = make_vector("B", 3, real=True)

# specific magnetic field values
B1 = make_vector("B1", 3, real=True)
B2 = make_vector("B2", 3, real=True)
B3 = make_vector("B3", 3, real=True)

# specific energy loss per distance values
g1 = Symbol("g1", real=True)
g2 = Symbol("g2", real=True)
g3 = Symbol("g3", real=True)
g4 = Symbol("g4", real=True)


def rk4_subexpr(f, x, y, ydot, h):
    k1 = name_expr("k1", f(1, x, y, ydot))
    x2 = name_expr("x2", x + h / 2)
    y2 = name_expr("y2", y + h / 2 * ydot + h**2 / 8 * k1.name)
    ydot2 = name_expr("ydot2", ydot + h / 2 * k1.name)

    k2 = name_expr("k2", f(2, x2.expr, y2.expr.as_explicit(), ydot2.expr.as_explicit()))
    ydot3 = name_expr("ydot3", ydot + h / 2 * k2.name)

    k3 = name_expr("k3", f(3, x2.expr, y2.expr.as_explicit(), ydot3.expr.as_explicit()))
    x3 = name_expr("x3", x + h)
    y3 = name_expr("y3", y + h * ydot + h**2 / 2 * k3.name)
    ydot4 = name_expr("ydot4", ydot + h * k3.name)

    k4 = name_expr("k4", f(4, x3.expr, y3.expr.as_explicit(), ydot4.expr.as_explicit()))

    new_y = name_expr("new_y", y + h * ydot + h**2 / 6 * (k1.name + k2.name + k3.name))
    new_ydot = name_expr(
        "new_ydot", ydot + h / 6 * (k1.name + 2 * (k2.name + k3.name) + k4.name)
    )

    dk1dyydot = name_expr("dk1dyydot", k1.expr.jacobian([y, ydot]))
    dk2dyydot = name_expr(
        "dk2dyydot",
        k2.expr.jacobian([y, ydot]) + k2.expr.jacobian(k1.name) * dk1dyydot.name,
    )
    dk3dyydot = name_expr(
        "dk3dyydot",
        k3.expr.jacobian([y, ydot]) + k3.expr.jacobian(k2.name) * dk2dyydot.name,
    )
    dk4dyydot = name_expr(
        "dk4dyydot",
        k4.expr.jacobian([y, ydot]) + k4.expr.jacobian(k3.name) * dk3dyydot.name,
    )

    dydyydot = name_expr(
        "dydyydot",
        new_y.expr.as_explicit().jacobian([y, ydot])
        + new_y.expr.as_explicit().jacobian(k1.name) * dk1dyydot.name
        + new_y.expr.as_explicit().jacobian(k2.name) * dk2dyydot.name
        + new_y.expr.as_explicit().jacobian(k3.name) * dk3dyydot.name,
    )
    dydotdyydot = name_expr(
        "dydotdyydot",
        new_ydot.expr.as_explicit().jacobian([y, ydot])
        + new_ydot.expr.as_explicit().jacobian(k1.name) * dk1dyydot.name
        + new_ydot.expr.as_explicit().jacobian(k2.name) * dk2dyydot.name
        + new_ydot.expr.as_explicit().jacobian(k3.name) * dk3dyydot.name
        + new_ydot.expr.as_explicit().jacobian(k4.name) * dk4dyydot.name,
    )

    return (
        ((new_y, new_ydot), (k1, k2, k3, k4)),
        ((dydyydot, dydotdyydot), (dk1dyydot, dk2dyydot, dk3dyydot, dk4dyydot)),
        (x2, y2, ydot2, ydot3, x3, y3, ydot4),
    )


def rk4_fullexpr(f, x, y, ydot, h):
    k1 = name_expr("k1", f(1, x, y, ydot))
    x2 = name_expr("x2", x + h / 2)
    y2 = name_expr("y2", y + h / 2 * ydot + h**2 / 8 * k1.expr)
    ydot2 = name_expr("ydot2", ydot + h / 2 * k1.expr)

    k2 = name_expr("k2", f(2, x2.expr, y2.expr, ydot2.expr))
    ydot3 = name_expr("ydot3", ydot + h / 2 * k2.expr)

    k3 = name_expr("k3", f(3, x2.expr, y2.expr, ydot3.expr))
    x3 = name_expr("x3", x + h)
    y3 = name_expr("y3", y + h * ydot + h**2 / 2 * k3.expr)
    ydot4 = name_expr("ydot4", ydot + h * k3.expr)

    k4 = name_expr("k4", f(4, x3.expr, y3.expr, ydot4.expr))

    new_y = name_expr("new_y", y + h * ydot + h**2 / 6 * (k1.expr + k2.expr + k3.expr))
    new_ydot = name_expr(
        "new_ydot", ydot + h / 6 * (k1.expr + 2 * (k2.expr + k3.expr) + k4.expr)
    )

    return (new_y, new_ydot), (k1, k2, k3, k4), (x2, y2, ydot2, ydot3, x3, y3, ydot4)


def rk4_vacuum_fullexpr2():
    def f(x, y, ydot):
        d = ydot[0:3, 0]
        return d.cross(l * B)

    def decorator(i, ydotdot):
        if i == 1:
            return ydotdot.subs(B, B1)
        if i in [2, 3]:
            return ydotdot.subs(B, B2)
        if i == 4:
            return ydotdot.subs(B, B3)

    (
        (new_y, new_ydot),
        (dydyydot, dydotdyydot),
        (k1, k2, k3, k4),
        (x2, y2, ydot2, ydot3, x3, y3, ydot4),
    ) = rk4_fullexpr(f, s, p, d, h, decorator)

    p2 = name_expr("p2", y2.expr[0:3, 0])
    p3 = name_expr("p3", y3.expr[0:3, 0])
    new_p = name_expr("new_p", new_y.expr[0:3, 0])
    new_d_tmp = name_expr("new_d_tmp", new_ydot.expr[0:3, 0])
    new_d = name_expr("new_d", new_d_tmp.expr / new_d_tmp.expr.norm())
    err = name_expr(
        "err",
        h**2
        * (k1.expr[0:3, 0] - k2.expr[0:3, 0] - k3.expr[0:3, 0] + k4.expr[0:3, 0]).norm(
            1
        ),
    )

    dtds = name_expr("dtds", sym.sqrt(1 + m**2 / p_abs**2))
    new_t = name_expr("new_t", t + h * dtds.expr)

    path_derivatives = name_expr("path_derivatives", sym.zeros(8, 1))
    path_derivatives.expr[0:3, 0] = new_d.expr
    path_derivatives.expr[3, 0] = dtds.expr
    path_derivatives.expr[4:7, 0] = k4.expr

    D = sym.eye(8)
    D[0:3, :] = new_p.expr.jacobian([p, t, d, l])
    D[4:7, :] = new_d_tmp.expr.jacobian([p, t, d, l])
    D[3, 7] = h * m**2 * l / dtds.expr

    J = MatrixSymbol("J", 8, 8).as_explicit().as_mutable()
    for indices in np.ndindex(J.shape):
        if D[indices] in [0, 1]:
            J[indices] = D[indices]
    J = ImmutableMatrix(J)

    new_J = name_expr("new_J", J * D)

    return [p2, p3, err, new_p, new_t, new_d, path_derivatives, new_J]


def rk4_vacuum_fullexpr():
    k1 = name_expr("k1", d.cross(l * B1))
    p2 = name_expr("p2", p + h / 2 * d + h**2 / 8 * k1.expr)

    k2 = name_expr("k2", (d + h / 2 * k1.expr).cross(l * B2))
    k3 = name_expr("k3", (d + h / 2 * k2.expr).cross(l * B2))
    p3 = name_expr("p3", p + h * d + h**2 / 2 * k3.expr)

    k4 = name_expr("k4", (d + h * k3.expr).cross(l * B3))

    err = name_expr("err", h**2 * (k1.expr - k2.expr - k3.expr + k4.expr).norm(1))

    new_p = name_expr("new_p", p + h * d + h**2 / 6 * (k1.expr + k2.expr + k3.expr))
    new_d_tmp = name_expr(
        "new_d_tmp", d + h / 6 * (k1.expr + 2 * (k2.expr + k3.expr) + k4.expr)
    )
    new_d = name_expr("new_d", new_d_tmp.expr / new_d_tmp.expr.norm())

    dtds = name_expr("dtds", sym.sqrt(1 + m**2 / p_abs**2))
    new_t = name_expr("new_t", t + h * dtds.expr)

    path_derivatives = name_expr("path_derivatives", sym.zeros(8, 1))
    path_derivatives.expr[0:3, 0] = new_d.expr
    path_derivatives.expr[3, 0] = dtds.expr
    path_derivatives.expr[4:7, 0] = k4.expr

    D = sym.eye(8)
    D[0:3, :] = new_p.expr.jacobian([p, t, d, l])
    D[4:7, :] = new_d_tmp.expr.jacobian([p, t, d, l])
    D[3, 7] = h * m**2 * l / dtds.expr

    J = MatrixSymbol("J", 8, 8).as_explicit().as_mutable()
    for indices in np.ndindex(J.shape):
        if D[indices] in [0, 1]:
            J[indices] = D[indices]
    J = ImmutableMatrix(J)

    new_J = name_expr("new_J", J * D)

    return [p2, p3, err, new_p, new_t, new_d, path_derivatives, new_J]


def rk4_vacuum_tunedexpr():
    k1 = name_expr("k1", d.cross(l * B1))
    p2 = name_expr("p2", p + h / 2 * d + h**2 / 8 * k1.name)

    k2 = name_expr("k2", (d + h / 2 * k1.name).as_explicit().cross(l * B2))
    k3 = name_expr("k3", (d + h / 2 * k2.name).as_explicit().cross(l * B2))
    p3 = name_expr("p3", p + h * d + h**2 / 2 * k3.name)

    k4 = name_expr("k4", (d + h * k3.name).as_explicit().cross(l * B3))

    err = name_expr(
        "err", h**2 * (k1.name - k2.name - k3.name + k4.name).as_explicit().norm(1)
    )

    new_p = name_expr("new_p", p + h * d + h**2 / 6 * (k1.name + k2.name + k3.name))
    new_d_tmp = name_expr(
        "new_d_tmp", d + h / 6 * (k1.name + 2 * (k2.name + k3.name) + k4.name)
    )
    new_d = name_expr("new_d", new_d_tmp.name / new_d_tmp.name.as_explicit().norm())

    dtds = name_expr("dtds", sym.sqrt(1 + m**2 / p_abs**2))
    new_t = name_expr("new_t", t + h * dtds.name)

    path_derivatives = name_expr("path_derivatives", sym.zeros(8, 1))
    path_derivatives.expr[0:3, 0] = new_d.name.as_explicit()
    path_derivatives.expr[3, 0] = dtds.name
    path_derivatives.expr[4:7, 0] = k4.name.as_explicit()

    dk1dTL = name_expr("dk1dTL", k1.expr.jacobian([d, l]))
    dk2dTL = name_expr(
        "dk2dTL", k2.expr.jacobian([d, l]) + k2.expr.jacobian(k1.name) * dk1dTL.expr
    )
    dk3dTL = name_expr(
        "dk3dTL",
        k3.expr.jacobian([d, l]) + k3.expr.jacobian(k2.name) * dk2dTL.name,
    )
    dk4dTL = name_expr(
        "dk4dTL",
        k4.expr.jacobian([d, l]) + k4.expr.jacobian(k3.name) * dk3dTL.name,
    )

    dFdTL = name_expr(
        "dFdTL",
        new_p.expr.as_explicit().jacobian([d, l])
        + new_p.expr.as_explicit().jacobian(k1.name) * dk1dTL.expr
        + new_p.expr.as_explicit().jacobian(k2.name) * dk2dTL.name
        + new_p.expr.as_explicit().jacobian(k3.name) * dk3dTL.name,
    )
    dGdTL = name_expr(
        "dGdTL",
        new_d_tmp.expr.as_explicit().jacobian([d, l])
        + new_d_tmp.expr.as_explicit().jacobian(k1.name) * dk1dTL.expr
        + new_d_tmp.expr.as_explicit().jacobian(k2.name) * dk2dTL.name
        + new_d_tmp.expr.as_explicit().jacobian(k3.name) * dk3dTL.name
        + new_d_tmp.expr.as_explicit().jacobian(k4.name) * dk4dTL.name,
    )

    D = sym.eye(8)
    D[0:3, 4:8] = dFdTL.name.as_explicit()
    D[4:7, 4:8] = dGdTL.name.as_explicit()
    D[3, 7] = h * m**2 * l / dtds.name

    J = Matrix(MatrixSymbol("J", 8, 8).as_explicit())
    for indices in np.ndindex(J.shape):
        if D[indices] in [0, 1]:
            J[indices] = D[indices]
    J = ImmutableMatrix(J)
    new_J = name_expr("new_J", J * D)

    return [
        k1,
        p2,
        k2,
        k3,
        p3,
        k4,
        err,
        new_p,
        dtds,
        new_t,
        new_d_tmp,
        new_d,
        path_derivatives,
        dk2dTL,
        dk3dTL,
        dk4dTL,
        dFdTL,
        dGdTL,
        new_J,
    ]


def rk4_dense_tunedexpr():
    def f(i, x, y, ydot):
        B = [B1, B2, B2, B3][i - 1]
        g = [g1, g2, g3, g4][i - 1]

        d = ydot[0:3, 0]
        dtds = ydot[3, 0]
        l = ydot[4, 0]
        return Matrix.vstack(
            d.cross(l * B),
            Matrix([g * m**2 * l**3 / q**3]),
            Matrix([dtds * l**2 * g / q]),
        )

    big_l = Symbol("big_l", real=True)
    dtds = name_expr("dtds", sym.sqrt(1 + m**2 / p_abs**2))

    (
        ((new_y, new_ydot), (k1, k2, k3, k4)),
        ((dydyydot, dydotdyydot), (dk1dyydot, dk2dyydot, dk3dyydot, dk4dyydot)),
        (x2, y2, ydot2, ydot3, x3, y3, ydot4),
    ) = rk4_subexpr(
        f,
        s,
        Matrix.vstack(p, Matrix([t, big_l])),
        Matrix.vstack(d, Matrix([dtds.name, l])),
        h,
    )

    p2 = name_expr("p2", y2.expr[0:3, 0])
    l2 = name_expr("l2", ydot2.expr[4, 0])
    p3 = name_expr("p3", y3.expr[0:3, 0])
    l3 = name_expr("l3", ydot3.expr[4, 0])
    l4 = name_expr("l4", ydot4.expr[4, 0])
    new_p = name_expr("new_p", new_y.expr[0:3, 0])
    new_t = name_expr("new_t", new_y.expr[3, 0])
    new_d_tmp = name_expr("new_d_tmp", new_ydot.expr[0:3, 0])
    new_d = name_expr("new_d", new_d_tmp.name / new_d_tmp.name.as_explicit().norm())
    new_l = name_expr("new_l", new_ydot.expr[4, 0])
    err = name_expr(
        "err",
        h**2
        * (k1.name[0:3, 0] - k2.name[0:3, 0] - k3.name[0:3, 0] + k4.name[0:3, 0])
        .as_explicit()
        .norm(1),
    )

    path_derivatives = name_expr("path_derivatives", sym.zeros(8, 1))
    path_derivatives.expr[0:3, 0] = new_d.name.as_explicit()
    path_derivatives.expr[3, 0] = new_ydot.name[3, 0]
    path_derivatives.expr[4:7, 0] = k4.name[0:3, 0].as_explicit()
    path_derivatives.expr[7, 0] = new_ydot.name[4, 0]

    dk1dTL = name_expr("dk1dTL", k1.expr.jacobian([t, d, l]))
    dk2dTL = name_expr(
        "dk2dTL", k2.expr.jacobian([t, d, l]) + k2.expr.jacobian(k1.name) * dk1dTL.expr
    )
    dk3dTL = name_expr(
        "dk3dTL",
        k3.expr.jacobian([t, d, l]) + k3.expr.jacobian(k2.name) * dk2dTL.name,
    )
    dk4dTL = name_expr(
        "dk4dTL",
        k4.expr.jacobian([t, d, l]) + k4.expr.jacobian(k3.name) * dk3dTL.name,
    )

    F = Matrix.vstack(new_p.expr.as_explicit(), Matrix([new_t.expr]))
    dFdTL = name_expr(
        "dFdTL",
        F.jacobian([t, d, l])
        + F.jacobian(k1.name) * dk1dTL.expr
        + F.jacobian(k2.name) * dk2dTL.name
        + F.jacobian(k3.name) * dk3dTL.name,
    )
    G = Matrix.vstack(new_d_tmp.expr.as_explicit(), Matrix([new_l.expr]))
    dGdTL = name_expr(
        "dGdTL",
        G.jacobian([t, d, l])
        + G.jacobian(k1.name) * dk1dTL.expr
        + G.jacobian(k2.name) * dk2dTL.name
        + G.jacobian(k3.name) * dk3dTL.name
        + G.jacobian(k4.name) * dk4dTL.name,
    )

    D = sym.eye(8)
    D[0:4, 3:8] = dFdTL.name.as_explicit()
    D[4:8, 3:8] = dGdTL.name.as_explicit()

    J = Matrix(MatrixSymbol("J", 8, 8).as_explicit())
    for indices in np.ndindex(J.shape):
        if D[indices] in [0, 1]:
            J[indices] = D[indices]
    J = ImmutableMatrix(J)
    new_J = name_expr("new_J", J * D)

    return [
        dtds,
        k1,
        y2,
        ydot2,
        ydot3,
        p2,
        l2,
        k2,
        l3,
        k3,
        y3,
        ydot4,
        p3,
        l4,
        k4,
        err,
        new_y,
        new_ydot,
        new_p,
        new_t,
        new_d_tmp,
        new_d,
        new_l,
        path_derivatives,
        dk2dTL,
        dk3dTL,
        dk4dTL,
        dFdTL,
        dGdTL,
        new_J,
    ]


def print_rk4_vacuum(name_exprs, run_cse=True):
    printer = cxx_printer
    outputs = [
        find_by_name(name_exprs, name)[0]
        for name in [
            "p2",
            "p3",
            "err",
            "new_p",
            "new_t",
            "new_d",
            "path_derivatives",
            "new_J",
        ]
    ]

    lines = []

    head = "template <typename T, typename GetB> Acts::Result<bool> rk4_vacuum(const T* p, const T* d, const T t, const T h, const T l, const T m, const T p_abs, GetB getB, T* err, const T errTol, T* new_p, T* new_t, T* new_d, T* path_derivatives, T* J) {"
    lines.append(head)

    lines.append("  const auto B1res = getB(p);")
    lines.append(
        "  if (!B1res.ok()) {\n    return Acts::Result<bool>::failure(B1res.error());\n  }"
    )
    lines.append("  const auto B1 = *B1res;")

    def pre_expr_hook(var):
        if str(var) == "p2":
            return "T p2[3];"
        if str(var) == "p3":
            return "T p3[3];"
        if str(var) == "new_J":
            return "T new_J[64];"
        return None

    def post_expr_hook(var):
        if str(var) == "p2":
            return "const auto B2res = getB(p2);\n  if (!B2res.ok()) {\n    return Acts::Result<bool>::failure(B2res.error());\n  }\n  const auto B2 = *B2res;"
        if str(var) == "p3":
            return "const auto B3res = getB(p3);\n  if (!B3res.ok()) {\n    return Acts::Result<bool>::failure(B3res.error());\n  }\n  const auto B3 = *B3res;"
        if str(var) == "err":
            return (
                "if (*err > errTol) {\n  return Acts::Result<bool>::success(false);\n}"
            )
        if str(var) == "new_d":
            return "if (J == nullptr) {\n  return Acts::Result<bool>::success(true);\n}"
        if str(var) == "new_J":
            return printer.doprint(Assignment(MatrixSymbol("J", 8, 8), var))
        return None

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
        pre_expr_hook=pre_expr_hook,
        post_expr_hook=post_expr_hook,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("  return Acts::Result<bool>::success(true);")

    lines.append("}")

    return "\n".join(lines)


def print_rk4_dense(name_exprs, run_cse=True):
    printer = cxx_printer
    outputs = [
        find_by_name(name_exprs, name)[0]
        for name in [
            "p2",
            "l2",
            "l3",
            "p3",
            "l4",
            "err",
            "new_p",
            "new_t",
            "new_d",
            "new_l",
            "path_derivatives",
            "new_J",
        ]
    ]

    lines = []

    head = "template <typename T, typename GetB, typename GetG> Acts::Result<bool> rk4_dense(const T* p, const T* d, const T t, const T h, const T l, const T m, const T q, const T p_abs, GetB getB, GetG getG, T* err, const T errTol, T* new_p, T* new_t, T* new_d, T* new_l, T* path_derivatives, T* J) {"
    lines.append(head)

    lines.append("  const auto B1res = getB(p);")
    lines.append(
        "  if (!B1res.ok()) {\n    return Acts::Result<bool>::failure(B1res.error());\n  }"
    )
    lines.append("  const auto B1 = *B1res;")
    lines.append("  const auto g1 = getG(p, l);")

    def pre_expr_hook(var):
        if str(var) == "p2":
            return "T p2[3];"
        if str(var) == "p3":
            return "T p3[3];"
        if str(var) == "l2":
            return "T l2[1];"
        if str(var) == "l3":
            return "T l3[1];"
        if str(var) == "l4":
            return "T l4[1];"
        if str(var) == "new_J":
            return "T new_J[64];"
        return None

    def post_expr_hook(var):
        if str(var) == "p2":
            return "const auto B2res = getB(p2);\n  if (!B2res.ok()) {\n    return Acts::Result<bool>::failure(B2res.error());\n  }\n  const auto B2 = *B2res;"
        if str(var) == "p3":
            return "const auto B3res = getB(p3);\n  if (!B3res.ok()) {\n    return Acts::Result<bool>::failure(B3res.error());\n  }\n  const auto B3 = *B3res;"
        if str(var) == "l2":
            return "const auto g2 = getG(p2, *l2);"
        if str(var) == "l3":
            return "const auto g3 = getG(p2, *l3);"
        if str(var) == "l4":
            return "const auto g4 = getG(p3, *l4);"
        if str(var) == "err":
            return (
                "if (*err > errTol) {\n  return Acts::Result<bool>::success(false);\n}"
            )
        if str(var) == "new_d":
            return "if (J == nullptr) {\n  return Acts::Result<bool>::success(true);\n}"
        if str(var) == "new_J":
            return printer.doprint(Assignment(MatrixSymbol("J", 8, 8), var))
        return None

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
        pre_expr_hook=pre_expr_hook,
        post_expr_hook=post_expr_hook,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("  return Acts::Result<bool>::success(true);")

    lines.append("}")

    return "\n".join(lines)


output.write(
    """
// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Note: This file is generated by generate_sympy_stepper.py
//       Do not modify it manually.

#pragma once

#include <Acts/Utilities/Result.hpp>

#include <cmath>
""".strip()
)

output.write("\n\n")

all_name_exprs = rk4_vacuum_tunedexpr()
# all_name_exprs = rk4_vacuum_fullexpr()
# all_name_exprs = rk4_vacuum_fullexpr2()

# attempted manual CSE which turns out a bit slower
#
# sub_name_exprs = [
#     name_expr("hlB1", h * l * B1),
#     name_expr("hlB2", h * l * B2),
#     name_expr("hlB3", h * l * B3),
#     name_expr("lB1", l * B1),
#     name_expr("lB2", l * B2),
#     name_expr("lB3", l * B3),
#     name_expr("h2_2", h**2 / 2),
#     name_expr("h_8", h / 8),
#     name_expr("h_6", h / 6),
#     name_expr("h_2", h / 2),
# ]
# all_name_exprs = [
#     NamedExpr(name, my_subs(expr, sub_name_exprs)) for name, expr in all_name_exprs
# ]
# all_name_exprs.extend(sub_name_exprs)

code = print_rk4_vacuum(
    all_name_exprs,
    run_cse=True,
)
output.write(code + "\n")

output.write("\n")

all_name_exprs = rk4_dense_tunedexpr()

code = print_rk4_dense(
    all_name_exprs,
    run_cse=True,
)
output.write(code + "\n")

if output is not sys.stdout:
    output.close()
