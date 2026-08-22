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
qop = Symbol("qop", real=True)

# step length
h = Symbol("h", real=True)

# path length
s = Symbol("s", real=True)

# position
position = make_vector("pos", 3, real=True)

# direction
direction = make_vector("dir", 3, real=True)

# time
time = Symbol("time", real=True)

# mass
mass = Symbol("mass", real=True, positive=True)

# absolute momentum
p_abs = Symbol("p_abs", real=True, positive=True)

# energy loss per distance, dE/ds
dEds = Symbol("dEds", real=True)

# charge
charge = Symbol("charge", real=True)

# magnetic field
B = make_vector("B", 3, real=True)

# specific magnetic field values
B1 = make_vector("B1", 3, real=True)
B2 = make_vector("B2", 3, real=True)
B3 = make_vector("B3", 3, real=True)

# energy loss per distance at each of the four stages
dEds1 = Symbol("dEds1", real=True)
dEds2 = Symbol("dEds2", real=True)
dEds3 = Symbol("dEds3", real=True)
dEds4 = Symbol("dEds4", real=True)


def dtime_dqop(dtds):
    """d(time)/d(q/p) over a step of length h.

    h mass^2 qop / (charge^2 dtds), with the charge^2 folded into
    p_abs = |charge| / |qop|.
    """
    return h * mass**2 / (p_abs**2 * qop * dtds)


def rk4_subexpr(eom, x, y, ydot, h):
    k1 = name_expr("k1", eom(1, x, y, ydot))
    x2 = name_expr("x2", x + h / 2)
    y2 = name_expr("y2", y + h / 2 * ydot + h**2 / 8 * k1.name)
    ydot2 = name_expr("ydot2", ydot + h / 2 * k1.name)

    k2 = name_expr(
        "k2", eom(2, x2.expr, y2.expr.as_explicit(), ydot2.expr.as_explicit())
    )
    ydot3 = name_expr("ydot3", ydot + h / 2 * k2.name)

    k3 = name_expr(
        "k3", eom(3, x2.expr, y2.expr.as_explicit(), ydot3.expr.as_explicit())
    )
    x3 = name_expr("x3", x + h)
    y3 = name_expr("y3", y + h * ydot + h**2 / 2 * k3.name)
    ydot4 = name_expr("ydot4", ydot + h * k3.name)

    k4 = name_expr(
        "k4", eom(4, x3.expr, y3.expr.as_explicit(), ydot4.expr.as_explicit())
    )

    new_y = name_expr("new_y", y + h * ydot + h**2 / 6 * (k1.name + k2.name + k3.name))
    new_ydot = name_expr(
        "new_ydot", ydot + h / 6 * (k1.name + 2 * (k2.name + k3.name) + k4.name)
    )

    dk1_dstate = name_expr("dk1_dstate", k1.expr.jacobian([y, ydot]))
    dk2_dstate = name_expr(
        "dk2_dstate",
        k2.expr.jacobian([y, ydot]) + k2.expr.jacobian(k1.name) * dk1_dstate.name,
    )
    dk3_dstate = name_expr(
        "dk3_dstate",
        k3.expr.jacobian([y, ydot]) + k3.expr.jacobian(k2.name) * dk2_dstate.name,
    )
    dk4_dstate = name_expr(
        "dk4_dstate",
        k4.expr.jacobian([y, ydot]) + k4.expr.jacobian(k3.name) * dk3_dstate.name,
    )

    dy_dstate = name_expr(
        "dy_dstate",
        new_y.expr.as_explicit().jacobian([y, ydot])
        + new_y.expr.as_explicit().jacobian(k1.name) * dk1_dstate.name
        + new_y.expr.as_explicit().jacobian(k2.name) * dk2_dstate.name
        + new_y.expr.as_explicit().jacobian(k3.name) * dk3_dstate.name,
    )
    dydot_dstate = name_expr(
        "dydot_dstate",
        new_ydot.expr.as_explicit().jacobian([y, ydot])
        + new_ydot.expr.as_explicit().jacobian(k1.name) * dk1_dstate.name
        + new_ydot.expr.as_explicit().jacobian(k2.name) * dk2_dstate.name
        + new_ydot.expr.as_explicit().jacobian(k3.name) * dk3_dstate.name
        + new_ydot.expr.as_explicit().jacobian(k4.name) * dk4_dstate.name,
    )

    return (
        ((new_y, new_ydot), (k1, k2, k3, k4)),
        ((dy_dstate, dydot_dstate), (dk1_dstate, dk2_dstate, dk3_dstate, dk4_dstate)),
        (x2, y2, ydot2, ydot3, x3, y3, ydot4),
    )


class AtlasStages:
    """The named quantities of one ATLAS-form RK4 step.

    `bend1..3` are the half-step bend vectors, `kick1`/`kick4` the first and
    last stage slopes (already scaled by h/2) and `dir2`/`dir3`/`dir4`/
    `dir_end` the direction estimates the intermediate stages are taken at.
    """

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def _atlas_rk4_stages(deriv, taylor_norm):
    """Build the value path of an ATLAS-form RK4 vacuum step.

    ATLAS carries a half-step bend vector, bend_i = (h*qop/2) * B_i, rather
    than the plain slopes k_i.  Every stage slope then comes out as (h/2)*k_i
    directly, as a bare cross product, and neither h nor qop appears again
    anywhere in the recursion.

    The stage names say which stage a quantity belongs to.  ATLAS' own names
    are positional; `docs/groups/sympy_codegen.md` maps the two onto each
    other (bend1..3 are H0..H2, and kick1/dir2/dir_half_sum/dir3/dir4/
    dir_end/kick4 are A0..A6).
    """
    # ATLAS calls this PS2.
    half_h_qop = deriv.add("half_h_qop", h * qop / 2)
    h_third = deriv.add("h_third", h / 3)
    h_quarter = deriv.add("h_quarter", h / 4)

    bend1 = deriv.add("bend1", half_h_qop.name * B1)
    kick1 = deriv.add("kick1", direction.cross(explicit(bend1.name)))  # h/2 * k1
    dir2 = deriv.add("dir2", explicit(kick1.name) + direction)  # dir + h/2 * k1
    dir_half_sum = deriv.add(
        "dir_half_sum", explicit(dir2.name) + direction
    )  # 2*dir + h/2 * k1
    deriv.add("pos2", position + h_quarter.name * explicit(dir_half_sum.name))

    bend2 = deriv.add("bend2", half_h_qop.name * B2)
    dir3 = deriv.add(
        "dir3", direction + explicit(dir2.name).cross(explicit(bend2.name))
    )  # dir + h/2 * k2
    dir4 = deriv.add(
        "dir4", direction + explicit(dir3.name).cross(explicit(bend2.name))
    )  # dir + h/2 * k3
    dir_end = deriv.add("dir_end", 2 * explicit(dir4.name) - direction)  # dir + h * k3
    deriv.add("pos3", position + h * explicit(dir4.name))

    bend3 = deriv.add("bend3", half_h_qop.name * B3)
    kick4 = deriv.add(
        "kick4", explicit(dir_end.name).cross(explicit(bend3.name))
    )  # h/2 * k4

    # (dir_half_sum+kick4)-(dir3+dir4) is h/2 * (k1-k2-k3+k4), hence the leading 2*|h|.
    err_vec = (explicit(dir_half_sum.name) + explicit(kick4.name)) - (
        explicit(dir3.name) + explicit(dir4.name)
    )
    deriv.add("err", 2 * sym.Abs(h) * err_vec.norm(1))

    deriv.add(
        "new_pos",
        position
        + h_third.name
        * (explicit(dir2.name) + explicit(dir3.name) + explicit(dir4.name)),
    )

    # new_dir_x3 = 3 * (dir + h/6*(k1 + 2k2 + 2k3 + k4)), the unnormalised new direction.
    new_dir_x3 = deriv.add(
        "new_dir_x3",
        2 * explicit(dir3.name)
        + (explicit(kick1.name) + explicit(dir_end.name) + explicit(kick4.name)),
    )

    if taylor_norm:
        # ATLAS replaces 1/|new_dir_x3| by its Taylor expansion around |new_dir_x3| = 3, which
        # RK4 stays within ~1e-7 of, so the truncation error is below rounding.
        norm_dev = deriv.add(
            "norm_dev",
            (new_dir_x3.name[0] ** 2 + new_dir_x3.name[1] ** 2)
            + (new_dir_x3.name[2] ** 2 - 9),
        )
        norm_fac = deriv.add(
            "norm_fac",
            sym.Rational(1, 3)
            - sym.Rational(1, 648) * norm_dev.name * (12 - norm_dev.name),
        )
        new_dir = deriv.add("new_dir", norm_fac.name * explicit(new_dir_x3.name))
    else:
        inv_norm = deriv.add("inv_norm", 1 / explicit(new_dir_x3.name).norm())
        new_dir = deriv.add("new_dir", inv_norm.name * explicit(new_dir_x3.name))

    dtds = deriv.add("dtds", sym.sqrt(1 + mass**2 / p_abs**2))
    deriv.add("new_time", time + h * dtds.name)

    # B3 is evaluated at pos + h*dir4, which differs from the step end point by
    # O(h^3).  Handing it back lets the caller reuse it as the next step's
    # first field sample and save one lookup per step, as ATLAS does.
    deriv.add("new_B", B3)

    two_over_h = deriv.add("two_over_h", 2 / h)  # undoes the h/2 scaling of kick4 -> k4
    deriv.add(
        "path_derivatives",
        Matrix.vstack(
            explicit(new_dir.name),
            Matrix([dtds.name]),
            two_over_h.name * explicit(kick4.name),
            Matrix([0]),
        ),
    )

    return AtlasStages(
        h_third=h_third,
        bend1=bend1,
        bend2=bend2,
        bend3=bend3,
        kick1=kick1,
        dir2=dir2,
        dir3=dir3,
        dir4=dir4,
        dir_end=dir_end,
        kick4=kick4,
        dtds=dtds,
    )


def _field_contrib(deriv, what, stage, bend, same_as, seed):
    """The term a tangent picks up from the bend vector's dependence on q/p.

    A bend vector is linear in q/p, so `qop * d(bend)/d(qop) == bend`, making
    this the bend-linear part of the stage, which is already named.  `same_as`
    encodes that identity and is checked against the plain chain rule before
    use.
    """
    contrib = stage.expr.jacobian(bend.name) * seed
    deriv.check_same(what, contrib[:, seed.cols - 1], same_as)
    return Matrix.hstack(contrib[:, 0 : seed.cols - 1], same_as)


# The bound-to-free jacobian M, stored column major.  "hold" is never written,
# "step" is written by every step, "dense" only by a dense step.
# fmt: off
_B2F = StructuredMatrix("M", [
    #  loc0    loc1    phi     theta   qop      time
    ["hold", "hold", "step", "step", "step",     0],  # pos[0]
    ["hold", "hold", "step", "step", "step",     0],  # pos[1]
    ["hold", "hold", "step", "step", "step",     0],  # pos[2]
    [     0,      0,      0,      0, "step",     1],  # time
    [     0,      0, "step", "step", "step",     0],  # dir[0]
    [     0,      0, "step", "step", "step",     0],  # dir[1]
    [     0,      0, "step", "step", "step",     0],  # dir[2]
    [     0,      0,      0,      0, "dense",    0],  # qop
])
# fmt: on

_B2F_LIVE = _B2F.entries("step")
_B2F_DENSE_LIVE = _B2F.entries("step", "dense")


def _b2f_column_name(col):
    """Name of the staging array holding live column `col` of the jacobian."""
    return f"new_M{col}"


def _by_column(entries):
    """`entries` grouped by column, as (staging array name, entries) pairs.

    One staging array per live column lets each column be written back as soon
    as it is finished, so only six or seven values are ever live at once
    instead of all nineteen -- the read-modify-write-per-column shape the
    ATLAS stepper has.
    """
    columns = {}
    for entry in entries:
        columns.setdefault(entry[1], []).append(entry)
    return [(_b2f_column_name(c), group) for c, group in sorted(columns.items())]


_B2F_LIVE_COLUMNS = _by_column(_B2F_LIVE)


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
    deriv = Derivation()
    st = _atlas_rk4_stages(deriv, taylor_norm)
    bend1, bend2, bend3 = st.bend1, st.bend2, st.bend3
    kick1, dir3, dir4, kick4 = st.kick1, st.dir3, st.dir4, st.kick4
    h_third, dtds = st.h_third, st.dtds

    M = _B2F.matrix

    def col(c, rows):
        return Matrix([[M[i, c]] for i in rows])

    def tangent(name, stage, contribs, field=None):
        acc = sym.zeros(3, 1) if field is None else field
        for var, tangent_of_var in contribs:
            acc = acc + stage.expr.jacobian(var) * tangent_of_var
        return deriv.add(name, acc)

    def propagate(tag, c, scale):
        """Push column `c` of M through the step.

        `scale` is None for a column with no q/p component.  The q/p column
        carries its direction part scaled by q/p, ATLAS' pVector[40]
        convention.
        """
        if scale is None:
            seed = col(c, (4, 5, 6))
            fields = [None] * 4
            pos_fac, dir_fac = h_third.name, third.name
        else:
            seed_qop_row = M[7, c]
            seed = explicit(deriv.add(f"{tag}_seed", qop * col(c, (4, 5, 6))).name)
            fields = [
                _field_contrib(
                    deriv,
                    f"{tag}_kick1",
                    kick1,
                    bend1,
                    explicit(kick1.name) * seed_qop_row,
                    explicit(bend1.name) * seed_qop_row,
                ),
                _field_contrib(
                    deriv,
                    f"{tag}_dir3",
                    dir3,
                    bend2,
                    (explicit(dir3.name) - direction) * seed_qop_row,
                    explicit(bend2.name) * seed_qop_row,
                ),
                _field_contrib(
                    deriv,
                    f"{tag}_dir4",
                    dir4,
                    bend2,
                    (explicit(dir4.name) - direction) * seed_qop_row,
                    explicit(bend2.name) * seed_qop_row,
                ),
                _field_contrib(
                    deriv,
                    f"{tag}_kick4",
                    kick4,
                    bend3,
                    explicit(kick4.name) * seed_qop_row,
                    explicit(bend3.name) * seed_qop_row,
                ),
            ]
            pos_fac, dir_fac = scale[0], scale[1]

        tan_kick1 = tangent(f"{tag}_kick1", kick1, [(direction, seed)], fields[0])
        tan_dir2 = deriv.add(f"{tag}_dir2", explicit(tan_kick1.name) + seed)
        tan_dir3 = tangent(
            f"{tag}_dir3",
            dir3,
            [(direction, seed), (st.dir2.name, explicit(tan_dir2.name))],
            fields[1],
        )
        tan_dir4 = tangent(
            f"{tag}_dir4",
            dir4,
            [(direction, seed), (dir3.name, explicit(tan_dir3.name))],
            fields[2],
        )
        tan_dir_end = deriv.add(f"{tag}_dir_end", 2 * explicit(tan_dir4.name) - seed)
        tan_kick4 = tangent(
            f"{tag}_kick4",
            kick4,
            [(st.dir_end.name, explicit(tan_dir_end.name))],
            fields[3],
        )

        new_pos = col(c, (0, 1, 2)) + pos_fac * (
            explicit(tan_dir2.name) + explicit(tan_dir3.name) + explicit(tan_dir4.name)
        )
        new_dir = dir_fac * (
            explicit(tan_kick1.name)
            + 2 * explicit(tan_dir3.name)
            + explicit(tan_dir_end.name)
            + explicit(tan_kick4.name)
        )
        return new_pos, new_dir

    # NOTE: dir_fac must be a Symbol, never a Number. sympy distributes a
    # Number over an Add, so Rational(1,3) * (v0 + 2*v3 + v5 + v6) emits four
    # multiplications per component where ATLAS' equivalent
    # ((d2A0 + 2*d2A3) + (d2A5 + d2A6)) * (1./3.) uses two -- eighteen wasted
    # multiplies once the three columns are counted. A Symbol does not
    # distribute, so the scaling stays factored and no temporary is needed.
    third = deriv.add("third", sym.Rational(1, 3))

    inv_qop = deriv.add("inv_qop", 1 / qop)
    qop_pos_weight = deriv.add("qop_pos_weight", h_third.name * inv_qop.name)
    qop_dir_weight = deriv.add("qop_dir_weight", inv_qop.name / 3)

    phi_pos, phi_dir = propagate("dphi", 2, None)
    the_pos, the_dir = propagate("dtheta", 3, None)
    qop_pos, qop_dir = propagate("dqop", 4, (qop_pos_weight.name, qop_dir_weight.name))

    # dt/ds depends on q/p only
    dt_dqop = deriv.add("dt_dqop", dtime_dqop(dtds.name))
    new_time = M[3, 4] + dt_dqop.name * M[7, 4]

    deriv.add(_b2f_column_name(2), Matrix.vstack(phi_pos, phi_dir))
    deriv.add(_b2f_column_name(3), Matrix.vstack(the_pos, the_dir))
    deriv.add(_b2f_column_name(4), Matrix.vstack(qop_pos, Matrix([new_time]), qop_dir))

    return deriv.name_exprs


def rk4_dense_tunedexpr():
    def eom(i, x, y, ydot):
        B = [B1, B2, B2, B3][i - 1]
        dEds = [dEds1, dEds2, dEds3, dEds4][i - 1]

        direction = ydot[0:3, 0]
        dtds = ydot[3, 0]
        qop = ydot[4, 0]
        return Matrix.vstack(
            direction.cross(qop * B),
            Matrix([dEds * mass**2 * qop**3 / charge**3]),
            Matrix([dtds * qop**2 * dEds / charge]),
        )

    qop_integral = Symbol("qop_integral", real=True)
    dtds = name_expr("dtds", sym.sqrt(1 + mass**2 / p_abs**2))

    (
        ((new_y, new_ydot), (k1, k2, k3, k4)),
        ((dy_dstate, dydot_dstate), (dk1_dstate, dk2_dstate, dk3_dstate, dk4_dstate)),
        (x2, y2, ydot2, ydot3, x3, y3, ydot4),
    ) = rk4_subexpr(
        eom,
        s,
        Matrix.vstack(position, Matrix([time, qop_integral])),
        Matrix.vstack(direction, Matrix([dtds.name, qop])),
        h,
    )

    pos2 = name_expr("pos2", y2.expr[0:3, 0])
    qop2 = name_expr("qop2", ydot2.expr[4, 0])
    pos3 = name_expr("pos3", y3.expr[0:3, 0])
    qop3 = name_expr("qop3", ydot3.expr[4, 0])
    qop4 = name_expr("qop4", ydot4.expr[4, 0])
    new_pos = name_expr("new_pos", new_y.expr[0:3, 0])
    new_time = name_expr("new_time", new_y.expr[3, 0])
    new_dir_unnorm = name_expr("new_dir_unnorm", new_ydot.expr[0:3, 0])
    new_dir = name_expr(
        "new_dir", new_dir_unnorm.name / new_dir_unnorm.name.as_explicit().norm()
    )
    new_qop = name_expr("new_qop", new_ydot.expr[4, 0])
    err = name_expr(
        "err",
        h**2
        * (k1.name[0:3, 0] - k2.name[0:3, 0] - k3.name[0:3, 0] + k4.name[0:3, 0])
        .as_explicit()
        .norm(1),
    )

    new_B = name_expr("new_B", B3)

    path_derivatives = name_expr("path_derivatives", sym.zeros(8, 1))
    path_derivatives.expr[0:3, 0] = new_dir.name.as_explicit()
    path_derivatives.expr[3, 0] = new_ydot.name[3, 0]
    path_derivatives.expr[4:7, 0] = k4.name[0:3, 0].as_explicit()
    path_derivatives.expr[7, 0] = new_ydot.name[4, 0]

    dk1_dfree = name_expr("dk1_dfree", k1.expr.jacobian([time, direction, qop]))
    dk2_dfree = name_expr(
        "dk2_dfree",
        k2.expr.jacobian([time, direction, qop])
        + k2.expr.jacobian(k1.name) * dk1_dfree.expr,
    )
    dk3_dfree = name_expr(
        "dk3_dfree",
        k3.expr.jacobian([time, direction, qop])
        + k3.expr.jacobian(k2.name) * dk2_dfree.name,
    )
    dk4_dfree = name_expr(
        "dk4_dfree",
        k4.expr.jacobian([time, direction, qop])
        + k4.expr.jacobian(k3.name) * dk3_dfree.name,
    )

    new_y_rows = Matrix.vstack(new_pos.expr.as_explicit(), Matrix([new_time.expr]))
    dy_dfree = name_expr(
        "dy_dfree",
        new_y_rows.jacobian([time, direction, qop])
        + new_y_rows.jacobian(k1.name) * dk1_dfree.expr
        + new_y_rows.jacobian(k2.name) * dk2_dfree.name
        + new_y_rows.jacobian(k3.name) * dk3_dfree.name,
    )
    new_ydot_rows = Matrix.vstack(
        new_dir_unnorm.expr.as_explicit(), Matrix([new_qop.expr])
    )
    dydot_dfree = name_expr(
        "dydot_dfree",
        new_ydot_rows.jacobian([time, direction, qop])
        + new_ydot_rows.jacobian(k1.name) * dk1_dfree.expr
        + new_ydot_rows.jacobian(k2.name) * dk2_dfree.name
        + new_ydot_rows.jacobian(k3.name) * dk3_dfree.name
        + new_ydot_rows.jacobian(k4.name) * dk4_dfree.name,
    )

    D = sym.eye(8)
    D[0:4, 3:8] = dy_dfree.name.as_explicit()
    D[4:8, 3:8] = dydot_dfree.name.as_explicit()

    # Unlike the vacuum kernel this one builds D, because energy loss couples
    # time and q/p into every stage.  D still applies to M, not to an 8x8.
    new_M = name_expr("new_M", b2f_step_update(D, _B2F_DENSE_LIVE))

    return [
        dtds,
        k1,
        y2,
        ydot2,
        ydot3,
        pos2,
        qop2,
        k2,
        qop3,
        k3,
        y3,
        ydot4,
        pos3,
        qop4,
        k4,
        err,
        new_B,
        new_y,
        new_ydot,
        new_pos,
        new_time,
        new_dir_unnorm,
        new_dir,
        new_qop,
        path_derivatives,
        dk2_dfree,
        dk3_dfree,
        dk4_dfree,
        dy_dfree,
        dydot_dfree,
        new_M,
    ]


def print_rk4_vacuum_b2f(name_exprs, run_cse=False):
    printer = cxx_printer
    outputs = [
        find_by_name(name_exprs, name)[0]
        for name in [
            "pos2",
            "pos3",
            "err",
            "new_B",
            "new_pos",
            "new_time",
            "new_dir",
            "path_derivatives",
            *[name for name, _ in _B2F_LIVE_COLUMNS],
        ]
    ]

    lines = []

    lines.append(
        "template <typename T, typename GetB>\n"
        "Acts::Result<bool> rk4_vacuum(std::span<const T, 3> pos,"
        " std::span<const T, 3> dir, const T time, const T h, const T qop, const T mass,"
        " const T p_abs, std::span<const T, 3> B1, GetB getB, T& err,"
        " const T errTol, std::span<T, 3> new_pos, T& new_time,"
        " std::span<T, 3> new_dir, std::span<T, 3> new_B,"
        " std::span<T, 8> path_derivatives, std::span<T> M) {"
    )
    lines.append(f"  assert(M.empty() || M.size() == {_B2F.size});")

    def pre_expr_hook(var):
        if str(var) == "pos2":
            return "std::array<T, 3> pos2{};"
        if str(var) == "pos3":
            return "std::array<T, 3> pos3{};"
        for name, group in _B2F_LIVE_COLUMNS:
            if str(var) == name:
                return f"std::array<T, {len(group)}> {name}{{}};"
        return None

    def post_expr_hook(var):
        if str(var) == "pos2":
            return "const auto B2res = getB(std::span<const T, 3>(pos2));\n  if (!B2res.ok()) {\n    return Acts::Result<bool>::failure(B2res.error());\n  }\n  const auto B2 = *B2res;"
        if str(var) == "pos3":
            return "const auto B3res = getB(std::span<const T, 3>(pos3));\n  if (!B3res.ok()) {\n    return Acts::Result<bool>::failure(B3res.error());\n  }\n  const auto B3 = *B3res;"
        if str(var) == "err":
            return (
                "if (err > errTol) {\n  return Acts::Result<bool>::success(false);\n}"
            )
        if str(var) == "new_dir":
            return "if (M.empty()) {\n  return Acts::Result<bool>::success(true);\n}"
        for name, group in _B2F_LIVE_COLUMNS:
            if str(var) == name:
                return "\n".join(
                    f"M[{_B2F.flat_index(i, j)}] = {name}[{k}];"
                    for k, (i, j) in enumerate(group)
                )
        return None

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
        pre_expr_hook=pre_expr_hook,
        post_expr_hook=post_expr_hook,
        scalar_outputs_by_pointer=False,
    )
    lines.extend([f"  {line}" for line in code.split("\n")])

    lines.append("  return Acts::Result<bool>::success(true);")

    lines.append("}")

    return "\n".join(lines)


def print_rk4_dense(name_exprs, run_cse=True):
    printer = cxx_printer
    outputs = [
        find_by_name(name_exprs, name)[0]
        for name in [
            "pos2",
            "qop2",
            "qop3",
            "pos3",
            "qop4",
            "err",
            "new_B",
            "new_pos",
            "new_time",
            "new_dir",
            "new_qop",
            "path_derivatives",
            "new_M",
        ]
    ]

    lines = []

    lines.append(
        "template <typename T, typename GetB, typename GetG>\n"
        "Acts::Result<bool> rk4_dense(std::span<const T, 3> pos,"
        " std::span<const T, 3> dir, const T time, const T h, const T qop, const T mass,"
        " const T charge, const T p_abs, std::span<const T, 3> B1, GetB getB,"
        " GetG getG, T& err, const T errTol, std::span<T, 3> new_pos, T& new_time,"
        " std::span<T, 3> new_dir, T& new_qop, std::span<T, 3> new_B,"
        " std::span<T, 8> path_derivatives, std::span<T> M) {"
    )
    lines.append(f"  assert(M.empty() || M.size() == {_B2F.size});")

    lines.append("  const auto dEds1 = getG(pos, qop);")

    def pre_expr_hook(var):
        if str(var) == "pos2":
            return "std::array<T, 3> pos2{};"
        if str(var) == "pos3":
            return "std::array<T, 3> pos3{};"
        if str(var) in ("qop2", "qop3", "qop4"):
            return f"T {var}{{}};"
        if str(var) == "new_M":
            return f"std::array<T, {len(_B2F_DENSE_LIVE)}> new_M{{}};"
        return None

    def post_expr_hook(var):
        if str(var) == "pos2":
            return "const auto B2res = getB(std::span<const T, 3>(pos2));\n  if (!B2res.ok()) {\n    return Acts::Result<bool>::failure(B2res.error());\n  }\n  const auto B2 = *B2res;"
        if str(var) == "pos3":
            return "const auto B3res = getB(std::span<const T, 3>(pos3));\n  if (!B3res.ok()) {\n    return Acts::Result<bool>::failure(B3res.error());\n  }\n  const auto B3 = *B3res;"
        if str(var) == "qop2":
            return "const auto dEds2 = getG(std::span<const T, 3>(pos2), qop2);"
        if str(var) == "qop3":
            return "const auto dEds3 = getG(std::span<const T, 3>(pos2), qop3);"
        if str(var) == "qop4":
            return "const auto dEds4 = getG(std::span<const T, 3>(pos3), qop4);"
        if str(var) == "err":
            return (
                "if (err > errTol) {\n  return Acts::Result<bool>::success(false);\n}"
            )
        if str(var) == "new_dir":
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
        scalar_outputs_by_pointer=False,
    )
    lines.extend([f"  {line}" for line in code.split("\n")])

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
