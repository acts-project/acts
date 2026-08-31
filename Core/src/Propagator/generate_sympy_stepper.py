import argparse
import contextlib
import functools
import random
import sys
from dataclasses import dataclass

import sympy as sym
from sympy import Symbol, Matrix

from codegen.sympy_common import (
    make_vector,
    NamedExpr,
    Derivation,
    StructuredMatrix,
    explicit,
    name_expr,
    find_by_name,
    cxx_printer,
    my_expression_print,
)

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


@dataclass
class Rk4Step:
    """Everything one plain RK4 step of a second order ODE produces."""

    #: the integrated state after the step
    new_y: NamedExpr
    #: its derivative after the step
    new_ydot: NamedExpr
    #: the four stage slopes
    k: tuple[NamedExpr, NamedExpr, NamedExpr, NamedExpr]
    #: d(new_y) / d(y, ydot)
    dy: NamedExpr
    #: d(new_ydot) / d(y, ydot)
    dydot: NamedExpr
    #: d(k_i) / d(y, ydot)
    dk: tuple[NamedExpr, NamedExpr, NamedExpr, NamedExpr]
    #: the stage points, kept because the kernels sample the field at them
    x2: NamedExpr
    y2: NamedExpr
    ydot2: NamedExpr
    ydot3: NamedExpr
    x3: NamedExpr
    y3: NamedExpr
    ydot4: NamedExpr

    def named_exprs(self) -> list[NamedExpr]:
        """Every named quantity the step introduced, in dependency order.

        `Derivation.resolve` relies on that order, so a name has to appear
        before the ones whose bodies mention it.
        """
        return [
            *self.k,
            self.x2,
            self.y2,
            self.ydot2,
            self.ydot3,
            self.x3,
            self.y3,
            self.ydot4,
            self.new_y,
            self.new_ydot,
            *self.dk,
            self.dy,
            self.dydot,
        ]


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

    return Rk4Step(
        new_y=new_y,
        new_ydot=new_ydot,
        k=(k1, k2, k3, k4),
        dy=dy_dstate,
        dydot=dydot_dstate,
        dk=(dk1_dstate, dk2_dstate, dk3_dstate, dk4_dstate),
        x2=x2,
        y2=y2,
        ydot2=ydot2,
        ydot3=ydot3,
        x3=x3,
        y3=y3,
        ydot4=ydot4,
    )


@dataclass
class AtlasStages:
    """The named quantities of one ATLAS-form RK4 step."""

    #: h/3, the position update weight
    h_third: NamedExpr
    #: the half-step bend vector at each of the three sample points
    bend1: NamedExpr
    bend2: NamedExpr
    bend3: NamedExpr
    #: the stage slopes and direction estimates, already scaled by h/2
    kick1: NamedExpr
    dir2: NamedExpr
    dir3: NamedExpr
    dir4: NamedExpr
    dir_end: NamedExpr
    kick4: NamedExpr
    #: dt/ds
    dtds: NamedExpr


def _atlas_rk4_stages(deriv: Derivation, taylor_norm: bool) -> AtlasStages:
    """Build the value path of an ATLAS-form RK4 vacuum step.

    ATLAS carries a half-step bend vector, bend_i = (h*qop/2) * B_i, rather
    than the plain slopes k_i, so every stage slope is a bare cross product
    and neither h nor qop appears again in the recursion.

    `docs/groups/sympy_codegen.md` maps these names onto ATLAS' positional
    ones: bend1..3 are H0..H2, kick1/dir2/dir_half_sum/dir3/dir4/dir_end/
    kick4 are A0..A6.
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

    # B3 sits O(h^3) from the step end point, so the caller can hand it back as
    # the next step's first sample and save one lookup.
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


def _field_contrib(
    deriv: Derivation,
    what: str,
    stage: NamedExpr,
    bend: NamedExpr,
    same_as: Matrix,
) -> Matrix:
    """The term a tangent picks up from the bend vector's dependence on q/p.

    A bend vector is linear in q/p, so `d(bend)/dlog|qop| == bend` and the
    term is the stage's bend-linear part, which the value path already named.
    @p same_as is that name; it is checked against the chain rule before use.
    """
    chain_rule = stage.expr.jacobian(bend.name) * explicit(bend.name)
    deriv.check_same(what, chain_rule, same_as)
    return same_as


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

    One staging array per column lets each be written back as soon as it is
    finished, instead of keeping every entry live until the end of the kernel.
    """
    columns = {}
    for entry in entries:
        columns.setdefault(entry[1], []).append(entry)
    return [(_b2f_column_name(c), group) for c, group in sorted(columns.items())]


_B2F_LIVE_COLUMNS = _by_column(_B2F_LIVE)

# The q/p column of M differentiates by the log of the current q/p, not by the
# starting q/p:
#
#     M[i, 4] = dFree_i/dlog|qop| = qop * (dFree_i/dqop_0) / (dqop/dqop_0)
#     M[7, 4] = dqop/dqop_0, kept plain -- it converts between the two
#
# The plain row makes the conversion exact both ways, and in this form each
# field term is a stage's bend-linear part (see _field_contrib). A vacuum step
# preserves it; rk4_dense moves q/p and converts around its M update.
# Derivation in docs/groups/sympy_codegen.md.
_B2F_QOP_COLUMN = 4
_B2F_QOP_ROW = 7


def _scale_qop_column(v, factor):
    """Scale the free rows of a q/p column, leaving its q/p row alone.

    That row is the same number in both conventions (see _B2F_QOP_COLUMN).
    """
    out = v.copy()
    out[:_B2F_QOP_ROW, 0] = v[:_B2F_QOP_ROW, 0] * factor
    return out


def b2f_step_update(
    D: Matrix, live: list[tuple[int, int]], qop_in=None, qop_out=None
) -> Matrix:
    """Apply a free-to-free step jacobian D to the bound-to-free jacobian M.

    Only the live columns are touched, and the structural zeros of _B2F keep
    the products sparse.  The vacuum kernel folds this into the RK recursion
    instead, and never builds D at all.

    qop_in and qop_out convert the q/p column to plain and back, which a step
    that moves q/p needs (see _B2F_QOP_COLUMN). Without them every column is
    treated as plain.
    """
    assert (qop_in is None) == (qop_out is None)
    out = []
    for c in sorted({col for _, col in live}):
        v = _B2F.matrix[:, c]
        scaled = c == _B2F_QOP_COLUMN and qop_in is not None
        if scaled:
            # to plain: undo the log factor qop_in, redo the chain rule row
            v = _scale_qop_column(v, v[_B2F_QOP_ROW, 0] / qop_in)
        new_v = D * v
        if scaled:
            # back again, against the row and the q/p the step leaves behind
            new_v = _scale_qop_column(new_v, qop_out / new_v[_B2F_QOP_ROW, 0])
        out.extend([new_v[i, 0]] for i, col in live if col == c)
    return Matrix(out)


def rk4_vacuum_b2f_atlasexpr(taylor_norm: bool = False) -> list[NamedExpr]:
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

    def propagate(tag, c, with_field):
        """Push column `c` of M through the step.

        `with_field` is False for a column with no q/p component, which picks
        up nothing from the bend vectors' own q/p dependence. The q/p column is
        read in its scaled form (see _B2F_QOP_COLUMN), where each field term is
        a stage quantity the value path already computed.
        """
        tan_in = col(c, (4, 5, 6))
        pos_fac, dir_fac = h_third.name, third.name
        if not with_field:
            fields = [None] * 4
        else:
            fields = [
                _field_contrib(
                    deriv, f"{tag}_kick1", kick1, bend1, explicit(kick1.name)
                ),
                _field_contrib(
                    deriv, f"{tag}_dir3", dir3, bend2, explicit(dir3.name) - direction
                ),
                _field_contrib(
                    deriv, f"{tag}_dir4", dir4, bend2, explicit(dir4.name) - direction
                ),
                _field_contrib(
                    deriv, f"{tag}_kick4", kick4, bend3, explicit(kick4.name)
                ),
            ]

        tan_kick1 = tangent(f"{tag}_kick1", kick1, [(direction, tan_in)], fields[0])
        tan_dir2 = deriv.add(f"{tag}_dir2", explicit(tan_kick1.name) + tan_in)
        tan_dir3 = tangent(
            f"{tag}_dir3",
            dir3,
            [(direction, tan_in), (st.dir2.name, explicit(tan_dir2.name))],
            fields[1],
        )
        tan_dir4 = tangent(
            f"{tag}_dir4",
            dir4,
            [(direction, tan_in), (dir3.name, explicit(tan_dir3.name))],
            fields[2],
        )
        tan_dir_end = deriv.add(f"{tag}_dir_end", 2 * explicit(tan_dir4.name) - tan_in)
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

    # Must be a Symbol, not a Number: sympy distributes a Number over an Add,
    # which emits one multiplication per term instead of one per component.
    third = deriv.add("third", sym.Rational(1, 3))

    phi_pos, phi_dir = propagate("dphi", 2, False)
    the_pos, the_dir = propagate("dtheta", 3, False)
    qop_pos, qop_dir = propagate("dqop", 4, True)

    # dt/ds depends on q/p only: the time row moves for the q/p column alone,
    # by d(h dt/ds)/dlog|qop|.
    dt_dqop = deriv.add("dt_dqop", qop * dtime_dqop(dtds.name))
    new_time = M[3, 4] + dt_dqop.name

    deriv.add(_b2f_column_name(2), Matrix.vstack(phi_pos, phi_dir))
    deriv.add(_b2f_column_name(3), Matrix.vstack(the_pos, the_dir))
    deriv.add(_b2f_column_name(4), Matrix.vstack(qop_pos, Matrix([new_time]), qop_dir))

    return deriv.name_exprs


# The naive reference: equations of motion, plain RK4, chain-rule derivatives.
# Both kernels are the same derivation with a different right hand side.
#
# The state is split the way a second order ODE wants it, laid out to match
# FreeVector so the generated code takes the parameter vector as it stands:
#
#   y    = (pos, time)          -> free indices 0,1,2 and 3
#   ydot = (dir, dt/ds, q/p)    -> free indices 4,5,6 and 7, plus dt/ds
#
# dt/ds rides in `ydot` because the dense right hand side evolves it, but it is
# not independent: dtds = sqrt(1 + mass^2 qop^2 / charge^2). Every derivative
# with respect to q/p therefore picks up the path through it, which is what
# `free_param_seed` carries.


def eom_vacuum(i, x, y, ydot):
    """d(ydot)/ds in vacuum: the field bends the direction, nothing else moves.

    @param i is the Runge-Kutta stage, 1 to 4, which selects the field sample
    """
    del x, y
    B = [B1, B2, B2, B3][i - 1]
    direction = ydot[0:3, 0]
    qop = ydot[4, 0]
    return Matrix.vstack(
        direction.cross(qop * B),
        Matrix([0]),  # no energy loss, so dt/ds is constant
        Matrix([0]),  # ... and so is q/p
    )


def eom_dense(i, x, y, ydot):
    """d(ydot)/ds through material: the same, plus energy loss on q/p.

    @param i is the Runge-Kutta stage, 1 to 4, which selects the samples
    """
    del x, y
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


def propagate_tangent(
    deriv: Derivation, tag: str, step: Rk4Step, y, ydot, u: Matrix
) -> tuple[Matrix, Matrix]:
    """Push one tangent vector through the Runge-Kutta recursion.

    The chain rule `rk4_subexpr` uses for the full derivative blocks,
    contracted with a vector at every stage instead of carried as a matrix. A
    column of the bound-to-free jacobian is such a vector, so this transports
    it without ever forming the 8x8 step jacobian.

    @param tag prefixes the stage tangent names, one set per column
    @param u is the seed tangent, in (y, ydot) space
    @return the tangents of the new state and of its new derivative
    """
    k = step.k
    ts = []
    for i in range(4):
        contrib = k[i].expr.jacobian([y, ydot]) * u
        if i > 0:
            contrib = contrib + k[i].expr.jacobian(k[i - 1].name) * ts[i - 1]
        ts.append(explicit(deriv.add(f"{tag}k{i + 1}", contrib).name))

    def combine(out):
        total = out.expr.as_explicit().jacobian([y, ydot]) * u
        for i in range(4):
            block = out.expr.as_explicit().jacobian(k[i].name)
            if any(e != 0 for e in block):
                total = total + block * ts[i]
        return total

    return combine(step.new_y), combine(step.new_ydot)


@dataclass
class DenseDerivation:
    """The dense derivation, with the pieces its self check needs."""

    #: every named quantity, in the order the printer wants them
    name_exprs: list[NamedExpr]
    #: the derivation the stage tangents were collected in
    deriv: Derivation
    #: q/p after the step, which the column rescaling needs
    new_qop: NamedExpr
    #: the transported bound-to-free columns
    new_M: NamedExpr


@functools.cache
def rk4_dense_derivation() -> DenseDerivation:
    qop_integral = Symbol("qop_integral", real=True)
    dtds = name_expr("dtds", sym.sqrt(1 + mass**2 / p_abs**2))

    step = rk4_subexpr(
        eom_dense,
        s,
        Matrix.vstack(position, Matrix([time, qop_integral])),
        Matrix.vstack(direction, Matrix([dtds.name, qop])),
        h,
    )

    k1, k2, k3, k4 = step.k
    new_y, new_ydot = step.new_y, step.new_ydot

    pos2 = name_expr("pos2", step.y2.expr[0:3, 0])
    qop2 = name_expr("qop2", step.ydot2.expr[4, 0])
    pos3 = name_expr("pos3", step.y3.expr[0:3, 0])
    qop3 = name_expr("qop3", step.ydot3.expr[4, 0])
    qop4 = name_expr("qop4", step.ydot4.expr[4, 0])
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

    # Transport the bound-to-free jacobian the way the vacuum kernel does:
    # push each live column through the same Runge-Kutta recursion as the
    # state, rather than assembling the 8x8 step jacobian and multiplying.
    deriv = Derivation()
    carried = [dtds, *step.named_exprs()]
    for ne in carried:
        deriv.record(ne)

    y = Matrix.vstack(position, Matrix([time, qop_integral]))
    ydot = Matrix.vstack(direction, Matrix([dtds.name, qop]))
    seed = free_param_seed()

    columns = []
    for tag, c in (("dphi", 2), ("dtheta", 3), ("dqop", _B2F_QOP_COLUMN)):
        v = _B2F.matrix[:, c]
        scaled = c == _B2F_QOP_COLUMN
        if scaled:
            # to plain: undo the log factor qop, redo the chain rule row
            v = _scale_qop_column(v, v[_B2F_QOP_ROW, 0] / qop)

        tan_y, tan_ydot = propagate_tangent(deriv, tag, step, y, ydot, seed * v)
        new_v = Matrix([tan_y[0:4, 0], tan_ydot[0:3, 0], Matrix([tan_ydot[4, 0]])])
        if scaled:
            # back again, against the row and the q/p the step leaves behind
            new_v = _scale_qop_column(new_v, new_qop.expr / new_v[_B2F_QOP_ROW, 0])
        columns.extend([new_v[i, 0]] for i, col in _B2F_DENSE_LIVE if col == c)

    tangent_exprs = deriv.name_exprs[len(carried) :]
    new_M = name_expr("new_M", Matrix(columns))

    name_exprs = [
        dtds,
        k1,
        step.y2,
        step.ydot2,
        step.ydot3,
        pos2,
        qop2,
        k2,
        qop3,
        k3,
        step.y3,
        step.ydot4,
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
        # the per-column stage tangents, in the order they were derived
        *tangent_exprs,
        new_M,
    ]
    return DenseDerivation(name_exprs, deriv, new_qop, new_M)


def rk4_dense_tangentexpr() -> list[NamedExpr]:
    """RK4 step through material, transporting the bound-to-free jacobian.

    The value path is the plain arrangement, not the ATLAS one: energy loss
    moves q/p every stage, so the half-step bend vector cannot be hoisted out
    of the recursion.
    """
    return rk4_dense_derivation().name_exprs


# Preconditions the derivations rely on but the signatures cannot express.
# p_abs is zero for a neutral, which makes dt/ds infinite.
INPUT_ASSERTS = """\
  assert(std::abs(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] - 1) <
         1e-8 && "direction must be a unit vector");
  assert(p_abs > 0 && "absolute momentum must be positive");
  assert(h != 0 && "step length must be non-zero");\
"""


def print_rk4_vacuum_b2f(name_exprs: list[NamedExpr], run_cse: bool = False) -> str:
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
        "Rk4Status rk4_vacuum(std::span<const T, 3> pos,"
        " std::span<const T, 3> dir, const T time, const T h, const T qop, const T mass,"
        " const T p_abs, std::span<const T, 3> B1, GetB getB, T& err,"
        " const T errTol, std::error_code& fieldErr,"
        " std::span<T, 3> new_pos, T& new_time,"
        " std::span<T, 3> new_dir, std::span<T, 3> new_B,"
        " std::span<T, 8> path_derivatives, std::span<T> M) {"
    )
    lines.append(INPUT_ASSERTS)
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
            return "const auto B2res = getB(std::span<const T, 3>(pos2));\n  if (!B2res.ok()) {\n    fieldErr = B2res.error();\n    return Rk4Status::FieldError;\n  }\n  const auto B2 = *B2res;"
        if str(var) == "pos3":
            return "const auto B3res = getB(std::span<const T, 3>(pos3));\n  if (!B3res.ok()) {\n    fieldErr = B3res.error();\n    return Rk4Status::FieldError;\n  }\n  const auto B3 = *B3res;"
        if str(var) == "err":
            return "if (err > errTol) {\n  return Rk4Status::Rejected;\n}"
        if str(var) == "new_dir":
            return "if (M.empty()) {\n  return Rk4Status::Accepted;\n}"
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
    )
    lines.extend([f"  {line}" for line in code.split("\n")])

    lines.append("  return Rk4Status::Accepted;")

    lines.append("}")

    return "\n".join(lines)


def print_rk4_dense(name_exprs: list[NamedExpr], run_cse: bool = True) -> str:
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
    lines.append(INPUT_ASSERTS)
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
        if str(var) == "new_qop":
            # new_qop carries the energy loss, so it belongs to the step and
            # has to be written before the jacobian-only part is skipped
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
    lines.extend([f"  {line}" for line in code.split("\n")])

    lines.append("  return Acts::Result<bool>::success(true);")

    lines.append("}")

    return "\n".join(lines)


def naive_rk4(eom) -> tuple[Derivation, Rk4Step]:
    """One plain RK4 step of the equations of motion, nothing rearranged.

    The reference the optimised forms are checked against: no named stage
    shortcuts, no reuse identities, no scaled field.

    @param eom is the right hand side, `f(stage, x, y, ydot) -> d(ydot)/ds`
    @return the Derivation every name resolves through, and the step
    """
    dtds_sym = Symbol("dtds", real=True, positive=True)
    # `rk4_subexpr` integrates y with y' = ydot, so the two have to be the same
    # length; `qop_integral` is the unused integral of q/p that pads y to match
    qop_integral = Symbol("qop_integral", real=True)
    y = Matrix.vstack(position, Matrix([time, qop_integral]))
    ydot = Matrix.vstack(direction, Matrix([dtds_sym, qop]))
    step = rk4_subexpr(eom, s, y, ydot, h)

    deriv = Derivation()
    # dt/ds is a component of the state but not an independent one; binding it
    # here is what lets a derivative with respect to q/p find its way through
    deriv.record(NamedExpr(dtds_sym, sym.sqrt(1 + mass**2 / p_abs**2)))
    for ne in step.named_exprs():
        deriv.record(ne)
    return deriv, step


def check_atlas_matches_naive() -> None:
    """Assert the ATLAS-form vacuum step computes the plain RK4 step.

    The arrangement cannot be *derived* from the naive form by substitution:
    once the naive expression is expanded, compound stage quantities like
    dir2 = kick1 + dir no longer appear syntactically. It is proved instead,
    by resolving every ATLAS name back down to the inputs and comparing.
    """
    naive, naive_step = naive_rk4(eom_vacuum)
    naive_y = explicit(naive.resolve(naive_step.new_y.expr))
    naive_ydot = explicit(naive.resolve(naive_step.new_ydot.expr))

    atlas = Derivation()
    _atlas_rk4_stages(atlas, taylor_norm=False)

    def resolved(name):
        return explicit(atlas.resolve(find_by_name(atlas.name_exprs, name)[1]))

    # the naive stage slopes, for the quantities ATLAS expresses over them
    naive_k = [explicit(naive.resolve(k.expr))[0:3, 0] for k in naive_step.k]

    checks = [
        ("new_pos", resolved("new_pos"), naive_y[0:3, 0]),
        ("new_time", Matrix([resolved("new_time")]), Matrix([naive_y[3, 0]])),
        # new_dir_x3 is three times the unnormalised new direction, by construction
        ("new_dir_x3", resolved("new_dir_x3"), 3 * naive_ydot[0:3, 0]),
        # Compared before the norm is taken: expand cannot reason about the
        # absolute values the two forms scale by.
        (
            "err vector",
            (resolved("dir_half_sum") + resolved("kick4"))
            - (resolved("dir3") + resolved("dir4")),
            (h / 2) * (naive_k[0] - naive_k[1] - naive_k[2] + naive_k[3]),
        ),
        # the path derivatives ATLAS reaches through two_over_h*kick4 are just k4
        (
            "path_derivatives (dir rows)",
            resolved("path_derivatives")[4:7, 0],
            naive_k[3],
        ),
    ]
    for what, atlas, naive_expr in checks:
        diff = sym.expand(Matrix(atlas) - Matrix(naive_expr))
        if any(e != 0 for e in diff):
            raise AssertionError(
                f"ATLAS form disagrees with plain RK4 on {what}:\n{diff}"
            )


def free_param_seed(at: dict | None = None) -> Matrix:
    """Map derivatives against (y, ydot) onto derivatives against free params.

    `rk4_subexpr` differentiates against (pos, time, qop_integral) and
    (dir, dt/ds, q/p). Two things have to happen to reach the eight free
    parameters:

    - `qop_integral`, the dummy that pads y to ydot's length, is not one, so its
      column goes;
    - `dt/ds` is not independent of q/p, so its column folds into the q/p one
      weighted by d(dt/ds)/d(q/p).

    @param at optionally evaluates the seed at a point
    @return the 10 x 8 map
    """
    dtds_sym = Symbol("dtds", real=True, positive=True)
    seed = sym.zeros(10, 8)
    for i in range(4):  # pos, time
        seed[i, i] = 1
    for i in range(3):  # dir
        seed[5 + i, 4 + i] = 1
    seed[8, 7] = mass**2 * qop / (charge**2 * dtds_sym)  # d(dt/ds)/d(q/p)
    seed[9, 7] = 1
    return seed if at is None else seed.subs(at)


#: the point every check evaluates at
SAMPLE_SEED = 20250814


@functools.cache
def sample_point(seed: int = SAMPLE_SEED) -> dict:
    """A rational point to evaluate a symbolic identity at.

    An identity between rational functions either holds everywhere or fails on
    a set of measure zero, so one arbitrary rational point settles it, and
    unlike the closed forms it stays small. Exact rationals, never floats, so
    a difference of zero means zero.

    @param seed fixes the point, so a failure is reproducible
    @return a substitution for every input the kernels take
    """
    rng = random.Random(seed)

    def r():
        return sym.Rational(rng.randint(2, 40), rng.randint(2, 40))

    # mass and momentum off a Pythagorean triple, so that
    # dt/ds = sqrt(1 + mass^2/p_abs^2) is rational too. A surd here propagates
    # into every stage and the nested radicals stop collapsing.
    u, v = rng.randint(3, 12), rng.randint(1, 2)

    point = {
        h: r(),
        time: r(),
        qop: r(),
        mass: sym.Integer(u**2 - v**2),
        p_abs: sym.Integer(2 * u * v),
        s: r(),
    }
    for vec in (position, direction, B1, B2, B3):
        for i in range(3):
            point[vec[i]] = r()
    for gi in (dEds1, dEds2, dEds3, dEds4):
        point[gi] = r()
    point[Symbol("dtds", real=True, positive=True)] = sym.Rational(
        u**2 + v**2, 2 * u * v
    )
    # p = |q| / |q/p| ties the three momentum inputs together
    point[charge] = point[p_abs] * point[qop]
    return point


@functools.cache
def naive_free_step_jacobian(eom, at_key: int | None = None) -> Matrix:
    """The free-to-free step jacobian D, straight from the chain rule.

    `free_param_seed` maps the (y, ydot) derivatives `rk4_subexpr` produces
    onto the eight free parameters. The direction rows are the *unnormalised*
    direction's, matching what every stepper here transports.

    Cached, since several checks want the same matrix and it is the expensive
    part of all of them.

    @param eom is the right hand side, as `eom_vacuum` or `eom_dense`
    @param at_key is a `sample_point` seed, or None to stay symbolic
    @return the 8x8 step jacobian over (pos, time, dir, q/p)
    """
    at = None if at_key is None else sample_point(at_key)
    deriv, step = naive_rk4(eom)
    dy = explicit(deriv.resolve(step.dy.expr, at))
    dydot = explicit(deriv.resolve(step.dydot.expr, at))

    seed = free_param_seed(at)

    D = sym.zeros(8, 8)
    D[0:4, :] = dy[0:4, :] * seed  # position and time
    D[4:7, :] = dydot[0:3, :] * seed  # direction, unnormalised
    D[7, :] = dydot[4, :] * seed  # q/p
    return D


def _check_transport_matches_naive(what, got, eom, live, qop_out) -> None:
    """Assert transported bound-to-free columns equal D times M.

    Neither kernel builds D; both contract the columns through the same
    Runge-Kutta recursion as the state, which is not obviously the same
    operation as the matrix product. Settled at a rational point, since the
    closed forms of the derivative blocks run to megabytes.

    @param what names the kernel, for the error message
    @param got the transported columns the kernel emits, at the point
    @param eom the right hand side, to build the naive D from
    @param live the entries the kernel writes back
    @param qop_out q/p after the step, which the column rescaling undoes with
    """
    at = sample_point()
    D = naive_free_step_jacobian(eom, SAMPLE_SEED)
    expected = b2f_step_update(D, live, at[qop], qop_out).subs(at)

    diff = Matrix(got) - Matrix(expected)
    bad = [i for i, e in enumerate(diff) if sym.simplify(e) != 0]
    if bad:
        # one entry is enough to debug from; these expressions are huge
        raise AssertionError(
            f"{what}: transported bound-to-free columns disagree with D * M at "
            f"live entries {[live[i] for i in bad]}; first difference:\n"
            f"{str(diff[bad[0]])[:2000]}"
        )


def check_jacobian_matches_naive() -> None:
    """Assert the vacuum kernel's transported columns equal D times M."""
    at = sample_point()
    deriv = Derivation()
    name_exprs = rk4_vacuum_b2f_atlasexpr()
    for ne in name_exprs:
        deriv.record(ne)
    got = Matrix.vstack(
        *[
            explicit(deriv.resolve(find_by_name(name_exprs, n)[1], at))
            for n in [name for name, _ in _B2F_LIVE_COLUMNS]
        ]
    )
    # q/p does not change in vacuum, so the column's scaling is the same on
    # both sides of the step
    _check_transport_matches_naive("vacuum", got, eom_vacuum, _B2F_LIVE, at[qop])


def check_dense_jacobian_matches_naive() -> None:
    """Assert the dense kernel's transported columns equal D times M.

    The only place the q/p column's rescaling is exercised against a q/p the
    step itself moved. SympyStepperTests runs the dense kernel on vacuum
    material, where the rescaling is the identity.
    """
    at = sample_point()
    dense = rk4_dense_derivation()
    got = explicit(dense.deriv.resolve(dense.new_M.expr, at))
    _check_transport_matches_naive(
        "dense",
        got,
        eom_dense,
        _B2F_DENSE_LIVE,
        dense.deriv.resolve(dense.new_qop.expr, at),
    )


def check_dropped_rows_stay_zero() -> None:
    """Assert the rows the bound-to-free update drops can never be populated.

    _B2F declares the time and q/p rows of the phi and theta columns zero and
    the update reads them as literal zeros. With those rows zero on the way
    in, (D M)[3, c] is the sum over j not in {3, 7} of D[3, j] M[j, c], so it
    vanishes for arbitrary M only if the time and q/p rows of D have support
    nowhere but the time and q/p columns.
    """
    # the free rows _B2F holds at zero in the phi and theta columns
    zero_rows = sorted({row for row, col in _B2F.entries(0) if col in (2, 3)})
    at = sample_point()
    seed = free_param_seed(at)
    for what, eom in (("vacuum", eom_vacuum), ("dense", eom_dense)):
        deriv, step = naive_rk4(eom)
        # only the time and q/p rows are needed; the direction rows are the
        # expensive part of the jacobian and have nothing to do with this
        rows = {
            3: explicit(deriv.resolve(step.dy.expr[3, :], at)) * seed,
            7: explicit(deriv.resolve(step.dydot.expr[4, :], at)) * seed,
        }
        assert sorted(rows) == zero_rows
        for row, values in rows.items():
            leaks = [
                col
                for col in range(8)
                if col not in zero_rows and sym.simplify(values[0, col]) != 0
            ]
            if leaks:
                raise AssertionError(
                    f"{what}: row {row} of the step jacobian reaches columns "
                    f"{leaks}, so _B2F would drop a live term"
                )


def check_dense_reduces_to_vacuum() -> None:
    """Assert the dense equations of motion become the vacuum ones without material.

    The two kernels alternate freely along a trajectory, so at zero energy
    loss they have to agree exactly.
    """
    vac, vac_step = naive_rk4(eom_vacuum)
    den, den_step = naive_rk4(eom_dense)

    no_loss = [(gi, 0) for gi in (dEds1, dEds2, dEds3, dEds4)]
    for what, a, b in [
        ("state", vac.resolve(vac_step.new_y.expr), den.resolve(den_step.new_y.expr)),
        (
            "derivative",
            vac.resolve(vac_step.new_ydot.expr),
            den.resolve(den_step.new_ydot.expr),
        ),
    ]:
        diff = sym.expand(explicit(a) - explicit(b).subs(no_loss))
        if any(e != 0 for e in diff):
            raise AssertionError(
                f"dense with zero energy loss differs from vacuum on {what}:\n{diff}"
            )


HEADER = """\
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
#include <system_error>

/// Outcome of one attempted Runge-Kutta step.
enum class Rk4Status {
  /// the error estimate exceeded the tolerance; nothing was written
  Rejected,
  /// the step was taken and the outputs written
  Accepted,
  /// a field lookup failed; the error is in the `fieldErr` output
  FieldError,
};
"""


def main(argv: list[str]) -> None:
    """Generate the Runge-Kutta kernels.

    @param argv is the command line, argv[0] being the program name
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        help="file to write the generated kernels to; stdout if omitted",
    )
    parser.add_argument(
        "--no-check",
        action="store_true",
        help="skip the symbolic equivalence assertions (they dominate the run "
        "time; only for iterating on the printers, never for a real build)",
    )
    args = parser.parse_args(argv[1:])

    if not args.no_check:
        # If one of these fires the generated code would be wrong, so they
        # guard the generator rather than a test that might not be run.
        check_atlas_matches_naive()
        check_dense_reduces_to_vacuum()
        check_jacobian_matches_naive()
        check_dense_jacobian_matches_naive()
        check_dropped_rows_stay_zero()

    with (
        open(args.output, "w") if args.output else contextlib.nullcontext(sys.stdout)
    ) as out:
        out.write(HEADER)
        out.write("\n")
        # taylor_norm and run_cse are off for the vacuum kernel: neither is
        # faster, and both obscure the correspondence to the ATLAS code.
        out.write(print_rk4_vacuum_b2f(rk4_vacuum_b2f_atlasexpr(taylor_norm=False)))
        out.write("\n\n")
        out.write(print_rk4_dense(rk4_dense_tangentexpr(), run_cse=True))
        out.write("\n")


if __name__ == "__main__":
    main(sys.argv)
