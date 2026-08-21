import sys

import numpy as np

import sympy as sym
from sympy import Symbol, Matrix, ImmutableMatrix
from sympy.codegen.ast import Assignment

from codegen.sympy_common import (
    make_vector,
    NamedExpr,
    Derivation,
    StructuredMatrix,
    explicit,
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


def dt_dqop(dtds):
    """d(t)/d(q/p) over a step of length h.

    h m**2 l / (q**2 dtds), with the q**2 folded into p_abs = |q| / |l|.
    """
    return h * m**2 / (p_abs**2 * l * dtds)


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


class AtlasStages:
    """The named quantities of one ATLAS-form RK4 step."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def _atlas_rk4_stages(b, taylor_norm):
    """Build the value path of an ATLAS-form RK4 vacuum step.

    ATLAS carries the half-step scaled field H = (h*l/2) * B rather than the
    plain slopes k_i, so every stage slope comes out as (h/2)*k_i directly and
    nothing is rescaled by h or l again.
    """
    # (h*l/2) is ATLAS' PS2.
    hl2 = b.add("hl2", h * l / 2)
    S3 = b.add("S3", h / 3)
    S4 = b.add("S4", h / 4)

    H0 = b.add("H0", hl2.name * B1)
    A0 = b.add("A0", d.cross(explicit(H0.name)))  # h/2 * k1
    A2 = b.add("A2", explicit(A0.name) + d)  # d + h/2 * k1
    A1 = b.add("A1", explicit(A2.name) + d)  # 2d + h/2 * k1
    b.add("p2", p + S4.name * explicit(A1.name))

    H1 = b.add("H1", hl2.name * B2)
    A3 = b.add("A3", d + explicit(A2.name).cross(explicit(H1.name)))  # d + h/2 * k2
    A4 = b.add("A4", d + explicit(A3.name).cross(explicit(H1.name)))  # d + h/2 * k3
    A5 = b.add("A5", 2 * explicit(A4.name) - d)  # d + h * k3
    b.add("p3", p + h * explicit(A4.name))

    H2 = b.add("H2", hl2.name * B3)
    A6 = b.add("A6", explicit(A5.name).cross(explicit(H2.name)))  # h/2 * k4

    # (A1+A6)-(A3+A4) is h/2 * (k1-k2-k3+k4), hence the leading 2*|h|.
    err_vec = (explicit(A1.name) + explicit(A6.name)) - (
        explicit(A3.name) + explicit(A4.name)
    )
    b.add("err", 2 * sym.Abs(h) * err_vec.norm(1))

    b.add(
        "new_p",
        p + S3.name * (explicit(A2.name) + explicit(A3.name) + explicit(A4.name)),
    )

    # An = 3 * (d + h/6*(k1 + 2k2 + 2k3 + k4)), the unnormalised new direction.
    An = b.add(
        "An",
        2 * explicit(A3.name)
        + (explicit(A0.name) + explicit(A5.name) + explicit(A6.name)),
    )

    if taylor_norm:
        # ATLAS replaces 1/|An| by its Taylor expansion around |An| = 3, which
        # RK4 stays within ~1e-7 of, so the truncation error is below rounding.
        Dv = b.add("Dv", (An.name[0] ** 2 + An.name[1] ** 2) + (An.name[2] ** 2 - 9))
        Dfac = b.add(
            "Dfac",
            sym.Rational(1, 3) - sym.Rational(1, 648) * Dv.name * (12 - Dv.name),
        )
        new_d = b.add("new_d", Dfac.name * explicit(An.name))
    else:
        inv_norm = b.add("inv_norm", 1 / explicit(An.name).norm())
        new_d = b.add("new_d", inv_norm.name * explicit(An.name))

    dtds = b.add("dtds", sym.sqrt(1 + m**2 / p_abs**2))
    b.add("new_t", t + h * dtds.name)

    Sl = b.add("Sl", 2 / h)  # undoes the h/2 scaling of A6 -> k4
    b.add(
        "path_derivatives",
        Matrix.vstack(
            explicit(new_d.name),
            Matrix([dtds.name]),
            Sl.name * explicit(A6.name),
            Matrix([0]),
        ),
    )

    return AtlasStages(
        S3=S3,
        H0=H0,
        H1=H1,
        H2=H2,
        A0=A0,
        A2=A2,
        A3=A3,
        A4=A4,
        A5=A5,
        A6=A6,
        dtds=dtds,
    )


def _field_contrib(b, what, stage, H, same_as, seed):
    """The term a tangent picks up from H's own dependence on l.

    H is linear in l, so `l * dH/dl == H`, making this the H-linear part of
    the stage, which is already named.  `same_as` encodes that identity and is
    checked against the plain chain rule before use.
    """
    contrib = stage.expr.jacobian(H.name) * seed
    b.check_same(what, contrib[:, seed.cols - 1], same_as)
    return Matrix.hstack(contrib[:, 0 : seed.cols - 1], same_as)


# The bound-to-free jacobian M, stored column major.  "hold" is never written,
# "step" is written by every step, "dense" only by a dense step.
# fmt: off
_B2F = StructuredMatrix("M", [
    #  loc0    loc1    phi     theta   qop      time
    ["hold", "hold", "step", "step", "step",     0],  # pos0
    ["hold", "hold", "step", "step", "step",     0],  # pos1
    ["hold", "hold", "step", "step", "step",     0],  # pos2
    [     0,      0,      0,      0, "step",     1],  # time
    [     0,      0, "step", "step", "step",     0],  # dir0
    [     0,      0, "step", "step", "step",     0],  # dir1
    [     0,      0, "step", "step", "step",     0],  # dir2
    [     0,      0,      0,      0, "dense",    0],  # qop
])
# fmt: on

_B2F_LIVE = _B2F.entries("step")
_B2F_DENSE_LIVE = _B2F.entries("step", "dense")


def b2f_step_update(D, live):
    """Apply a free-to-free step jacobian D to the bound-to-free jacobian M.

    Only the live columns are touched, and the structural zeros of _B2F keep
    the products sparse.  The vacuum kernel folds this into the RK recursion
    instead, and never builds D at all.
    """
    out = []
    for c in sorted({col for _, col in live}):
        new_v = D * _B2F.matrix[:, c]
        out.extend([new_v[i, 0]] for i, col in live if col == c)
    return Matrix(out)


def rk4_vacuum_b2f_atlasexpr(taylor_norm=False):
    """ATLAS-form RK4 vacuum step transporting the bound-to-free jacobian.

    The step jacobian D is never built: the columns of M go through the same
    RK recursion as the track state, which is why the code below mirrors the
    value path line for line.
    """
    b = Derivation()
    st = _atlas_rk4_stages(b, taylor_norm)
    H0, H1, H2 = st.H0, st.H1, st.H2
    A0, A3, A4, A6 = st.A0, st.A3, st.A4, st.A6
    S3, dtds = st.S3, st.dtds

    M = _B2F.matrix

    def col(c, rows):
        return Matrix([[M[i, c]] for i in rows])

    def tangent(name, stage, contribs, field=None):
        T = sym.zeros(3, 1) if field is None else field
        for var, Tvar in contribs:
            T = T + stage.expr.jacobian(var) * Tvar
        return b.add(name, T)

    def propagate(tag, c, scale):
        """Push column `c` of M through the step.

        `scale` is None for a column with no q/p component.  The q/p column
        carries its direction part scaled by l, ATLAS' pVector[40] convention.
        """
        if scale is None:
            seed = col(c, (4, 5, 6))
            fields = [None] * 4
            pos_fac, dir_fac = S3.name, sym.Rational(1, 3)
        else:
            vl = M[7, c]
            seed = explicit(b.add(f"u{tag}", l * col(c, (4, 5, 6))).name)
            fields = [
                _field_contrib(
                    b, f"{tag}0", A0, H0, explicit(A0.name) * vl, explicit(H0.name) * vl
                ),
                _field_contrib(
                    b,
                    f"{tag}3",
                    A3,
                    H1,
                    (explicit(A3.name) - d) * vl,
                    explicit(H1.name) * vl,
                ),
                _field_contrib(
                    b,
                    f"{tag}4",
                    A4,
                    H1,
                    (explicit(A4.name) - d) * vl,
                    explicit(H1.name) * vl,
                ),
                _field_contrib(
                    b, f"{tag}6", A6, H2, explicit(A6.name) * vl, explicit(H2.name) * vl
                ),
            ]
            pos_fac, dir_fac = scale[0], scale[1]

        v0 = tangent(f"{tag}0", A0, [(d, seed)], fields[0])
        v2 = b.add(f"{tag}2", explicit(v0.name) + seed)
        v3 = tangent(
            f"{tag}3", A3, [(d, seed), (st.A2.name, explicit(v2.name))], fields[1]
        )
        v4 = tangent(
            f"{tag}4", A4, [(d, seed), (A3.name, explicit(v3.name))], fields[2]
        )
        v5 = b.add(f"{tag}5", 2 * explicit(v4.name) - seed)
        v6 = tangent(f"{tag}6", A6, [(st.A5.name, explicit(v5.name))], fields[3])

        new_pos = col(c, (0, 1, 2)) + pos_fac * (
            explicit(v2.name) + explicit(v3.name) + explicit(v4.name)
        )
        new_dir = dir_fac * (
            explicit(v0.name)
            + 2 * explicit(v3.name)
            + explicit(v5.name)
            + explicit(v6.name)
        )
        return new_pos, new_dir

    inv_l = b.add("inv_l", 1 / l)
    S3_l = b.add("S3_l", S3.name * inv_l.name)
    third_l = b.add("third_l", inv_l.name / 3)

    phi_pos, phi_dir = propagate("Tp", 2, None)
    the_pos, the_dir = propagate("Tt", 3, None)
    qop_pos, qop_dir = propagate("Tq", 4, (S3_l.name, third_l.name))

    # dt/ds depends on q/p only
    dtdl = b.add("dtdl", dt_dqop(dtds.name))
    new_time = M[3, 4] + dtdl.name * M[7, 4]

    b.add(
        "new_M",
        Matrix.vstack(
            phi_pos,
            phi_dir,
            the_pos,
            the_dir,
            qop_pos,
            Matrix([new_time]),
            qop_dir,
        ),
    )

    return b.name_exprs


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

    # Unlike the vacuum kernel this one builds D, because energy loss couples
    # time and q/p into every stage.  D still applies to M, not to an 8x8.
    new_M = name_expr("new_M", b2f_step_update(D, _B2F_DENSE_LIVE))

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
        new_M,
    ]


def print_rk4_vacuum_b2f(name_exprs, run_cse=False):
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
            "new_M",
        ]
    ]

    lines = []

    lines.append(
        "template <typename T, typename GetB>\n"
        "Acts::Result<bool> rk4_vacuum(std::span<const T, 3> p,"
        " std::span<const T, 3> d, const T t, const T h, const T l, const T m,"
        " const T p_abs, GetB getB, T* err, const T errTol,"
        " std::span<T, 3> new_p, T* new_t, std::span<T, 3> new_d,"
        " std::span<T, 8> path_derivatives, std::span<T> M) {"
    )
    lines.append(f"  assert(M.empty() || M.size() == {_B2F.size});")

    lines.append("  const auto B1res = getB(p);")
    lines.append(
        "  if (!B1res.ok()) {\n    return Acts::Result<bool>::failure(B1res.error());\n  }"
    )
    lines.append("  const auto B1 = *B1res;")

    def pre_expr_hook(var):
        if str(var) == "p2":
            return "std::array<T, 3> p2{};"
        if str(var) == "p3":
            return "std::array<T, 3> p3{};"
        if str(var) == "new_M":
            return f"std::array<T, {len(_B2F_LIVE)}> new_M{{}};"
        return None

    def post_expr_hook(var):
        if str(var) == "p2":
            return "const auto B2res = getB(std::span<const T, 3>(p2));\n  if (!B2res.ok()) {\n    return Acts::Result<bool>::failure(B2res.error());\n  }\n  const auto B2 = *B2res;"
        if str(var) == "p3":
            return "const auto B3res = getB(std::span<const T, 3>(p3));\n  if (!B3res.ok()) {\n    return Acts::Result<bool>::failure(B3res.error());\n  }\n  const auto B3 = *B3res;"
        if str(var) == "err":
            return (
                "if (*err > errTol) {\n  return Acts::Result<bool>::success(false);\n}"
            )
        if str(var) == "new_d":
            return "if (M.empty()) {\n  return Acts::Result<bool>::success(true);\n}"
        if str(var) == "new_M":
            return "\n".join(
                f"M[{_B2F.flat_index(i, j)}] = new_M[{k}];"
                for k, (i, j) in enumerate(_B2F_LIVE)
            )
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
            "new_M",
        ]
    ]

    lines = []

    lines.append(
        "template <typename T, typename GetB, typename GetG>\n"
        "Acts::Result<bool> rk4_dense(std::span<const T, 3> p,"
        " std::span<const T, 3> d, const T t, const T h, const T l, const T m,"
        " const T q, const T p_abs, GetB getB, GetG getG, T* err,"
        " const T errTol, std::span<T, 3> new_p, T* new_t,"
        " std::span<T, 3> new_d, T* new_l, std::span<T, 8> path_derivatives,"
        " std::span<T> M) {"
    )
    lines.append(f"  assert(M.empty() || M.size() == {_B2F.size});")

    lines.append("  const auto B1res = getB(p);")
    lines.append(
        "  if (!B1res.ok()) {\n    return Acts::Result<bool>::failure(B1res.error());\n  }"
    )
    lines.append("  const auto B1 = *B1res;")
    lines.append("  const auto g1 = getG(p, l);")

    def pre_expr_hook(var):
        if str(var) == "p2":
            return "std::array<T, 3> p2{};"
        if str(var) == "p3":
            return "std::array<T, 3> p3{};"
        if str(var) == "l2":
            return "T l2[1];"
        if str(var) == "l3":
            return "T l3[1];"
        if str(var) == "l4":
            return "T l4[1];"
        if str(var) == "new_M":
            return f"std::array<T, {len(_B2F_DENSE_LIVE)}> new_M{{}};"
        return None

    def post_expr_hook(var):
        if str(var) == "p2":
            return "const auto B2res = getB(std::span<const T, 3>(p2));\n  if (!B2res.ok()) {\n    return Acts::Result<bool>::failure(B2res.error());\n  }\n  const auto B2 = *B2res;"
        if str(var) == "p3":
            return "const auto B3res = getB(std::span<const T, 3>(p3));\n  if (!B3res.ok()) {\n    return Acts::Result<bool>::failure(B3res.error());\n  }\n  const auto B3 = *B3res;"
        if str(var) == "l2":
            return "const auto g2 = getG(std::span<const T, 3>(p2), *l2);"
        if str(var) == "l3":
            return "const auto g3 = getG(std::span<const T, 3>(p2), *l3);"
        if str(var) == "l4":
            return "const auto g4 = getG(std::span<const T, 3>(p3), *l4);"
        if str(var) == "err":
            return (
                "if (*err > errTol) {\n  return Acts::Result<bool>::success(false);\n}"
            )
        if str(var) == "new_d":
            return "if (M.empty()) {\n  return Acts::Result<bool>::success(true);\n}"
        if str(var) == "new_M":
            return "\n".join(
                f"M[{_B2F.flat_index(i, j)}] = new_M[{k}];"
                for k, (i, j) in enumerate(_B2F_DENSE_LIVE)
            )
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


output.write("""
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

#include "Acts/Utilities/Result.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <span>
""".strip())

output.write("\n\n")

# taylor_norm and run_cse are off: neither is faster, and both obscure the
# correspondence to the ATLAS code.
all_name_exprs = rk4_vacuum_b2f_atlasexpr(taylor_norm=False)

code = print_rk4_vacuum_b2f(all_name_exprs)
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
