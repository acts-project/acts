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


def dt_dqop(dtds):
    """d(t)/d(q/p) over a step of length h.

    Equals h m**2 l / (q**2 dtds), with the q**2 folded into p_abs via
    p_abs = |q| / |l|.
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


def _ex(x):
    """Explicit matrix view of a MatrixSymbol (or pass through)."""
    return x.as_explicit() if hasattr(x, "as_explicit") else x


class _Builder:
    """Collects named expressions and can resolve them back to closed form.

    The closed form is only needed for the self checks: the ATLAS recursion
    below reuses a handful of already computed quantities as derivative terms,
    and every one of those shortcuts is verified against what the plain chain
    rule produces before it is allowed into the generated code.
    """

    def __init__(self):
        self.name_exprs = []
        self._by_name = {}

    def add(self, name, expr):
        if isinstance(expr, Matrix):
            expr = ImmutableMatrix(expr)
        ne = name_expr(name, expr)
        self.name_exprs.append(ne)
        self._by_name[ne.name] = ne.expr
        return ne

    def resolve(self, expr):
        subs = list(self._by_name.items())
        for _ in range(64):
            new = expr.subs(subs)
            if new == expr:
                return sym.expand(new)
            expr = new
        raise RuntimeError("named expressions did not resolve")

    def check_same(self, what, expr_a, expr_b):
        diff = sym.simplify(sym.expand(self.resolve(expr_a) - self.resolve(expr_b)))
        if any(e != 0 for e in diff):
            raise AssertionError(f"{what}: shortcut does not match chain rule\n{diff}")


# tangent seed for the direction: [ dd/dd | l * dd/dl ] = [ I | 0 ]
class _Stages:
    """The named quantities of one ATLAS-form RK4 step."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def _atlas_rk4_stages(b, taylor_norm):
    """Build the value path of an ATLAS-form RK4 vacuum step.

    ATLAS never evaluates the plain slopes k_i.  It carries the half-step
    scaled field H = (h*l/2) * B, so every stage slope comes out as (h/2)*k_i
    directly and no further scaling by h or l appears anywhere in the
    recursion.
    """
    # (h*l/2) is ATLAS' PS2.
    hl2 = b.add("hl2", h * l / 2)
    S3 = b.add("S3", h / 3)
    S4 = b.add("S4", h / 4)

    H0 = b.add("H0", hl2.name * B1)
    A0 = b.add("A0", d.cross(_ex(H0.name)))  # h/2 * k1
    A2 = b.add("A2", _ex(A0.name) + d)  # d + h/2 * k1
    A1 = b.add("A1", _ex(A2.name) + d)  # 2d + h/2 * k1
    b.add("p2", p + S4.name * _ex(A1.name))

    H1 = b.add("H1", hl2.name * B2)
    A3 = b.add("A3", d + _ex(A2.name).cross(_ex(H1.name)))  # d + h/2 * k2
    A4 = b.add("A4", d + _ex(A3.name).cross(_ex(H1.name)))  # d + h/2 * k3
    A5 = b.add("A5", 2 * _ex(A4.name) - d)  # d + h * k3
    b.add("p3", p + h * _ex(A4.name))

    H2 = b.add("H2", hl2.name * B3)
    A6 = b.add("A6", _ex(A5.name).cross(_ex(H2.name)))  # h/2 * k4

    # (A1+A6)-(A3+A4) is h/2 * (k1-k2-k3+k4), hence the leading 2*|h|.
    err_vec = (_ex(A1.name) + _ex(A6.name)) - (_ex(A3.name) + _ex(A4.name))
    b.add("err", 2 * sym.Abs(h) * err_vec.norm(1))

    b.add("new_p", p + S3.name * (_ex(A2.name) + _ex(A3.name) + _ex(A4.name)))

    # An = 3 * (d + h/6*(k1 + 2k2 + 2k3 + k4)), three times the unnormalised
    # new direction.
    An = b.add("An", 2 * _ex(A3.name) + (_ex(A0.name) + _ex(A5.name) + _ex(A6.name)))

    if taylor_norm:
        # ATLAS replaces 1/|An| by its second order Taylor expansion around
        # |An| = 3, which needs neither a square root nor a division.  RK4
        # keeps |An| within ~1e-7 of 3, so the O(u^3) truncation error is far
        # below double precision.
        Dv = b.add("Dv", (An.name[0] ** 2 + An.name[1] ** 2) + (An.name[2] ** 2 - 9))
        Dfac = b.add(
            "Dfac",
            sym.Rational(1, 3) - sym.Rational(1, 648) * Dv.name * (12 - Dv.name),
        )
        new_d = b.add("new_d", Dfac.name * _ex(An.name))
    else:
        inv_norm = b.add("inv_norm", 1 / _ex(An.name).norm())
        new_d = b.add("new_d", inv_norm.name * _ex(An.name))

    dtds = b.add("dtds", sym.sqrt(1 + m**2 / p_abs**2))
    b.add("new_t", t + h * dtds.name)

    # B3 is evaluated at p + h*A4, which differs from the step end point by
    # O(h^3).  Handing it back lets the caller reuse it as the next step's
    # first field sample and save one lookup per step, as ATLAS does.
    b.add("new_B", B3)

    Sl = b.add("Sl", 2 / h)  # undoes the h/2 scaling of A6 -> k4
    b.add(
        "path_derivatives",
        Matrix.vstack(
            _ex(new_d.name),
            Matrix([dtds.name]),
            Sl.name * _ex(A6.name),
            Matrix([0]),
        ),
    )

    return _Stages(
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

    H is linear in l, so `l * dH/dl == H` exactly; with the tangent's l part
    scaled by l, this term is the H-linear part of the stage, which is already
    available as a named quantity.  `same_as` encodes that identity and is
    verified against the plain chain rule before use -- it is what keeps the
    recursion free of any extra cross product.
    """
    contrib = stage.expr.jacobian(H.name) * seed
    b.check_same(what, contrib[:, seed.cols - 1], same_as)
    return Matrix.hstack(contrib[:, 0 : seed.cols - 1], same_as)


# Entries of the bound-to-free jacobian that a vacuum step can change, as
# (free row, bound column) with the ACTS orderings
# rows    (pos0, pos1, pos2, time, dir0, dir1, dir2, qop)
# columns (loc0, loc1, phi, theta, qop, time).
#
# loc0 and loc1 have no direction and no q/p component, and the time column is
# exactly e_time, so those three columns cannot move: with zero direction and
# zero q/p perturbation the step leaves the position rows alone as well.  That
# is the same reduction the ATLAS stepper makes when it propagates only the
# blocks at pVector[24], [32] and [40].  Within the live columns the q/p row
# never changes because l' = l in vacuum, and the time row only moves for the
# q/p column because dt/ds depends on q/p alone.
# One group per live column. Keeping them separate lets each column be written
# back as soon as it is finished, so only six or seven values are ever live at
# once instead of all nineteen -- the same read-modify-write-per-column shape
# the ATLAS stepper has.
_B2F_LIVE_COLUMNS = [
    ("new_Mp", [(i, 2) for i in (0, 1, 2, 4, 5, 6)]),
    ("new_Mt", [(i, 3) for i in (0, 1, 2, 4, 5, 6)]),
    ("new_Mq", [(i, 4) for i in (0, 1, 2, 3, 4, 5, 6)]),
]

_B2F_LIVE = [e for _, g in _B2F_LIVE_COLUMNS for e in g]

# A dense step additionally changes q/p itself, so the q/p row of the q/p
# column moves too.
_B2F_DENSE_LIVE = (
    [(i, 2) for i in (0, 1, 2, 4, 5, 6)]
    + [(i, 3) for i in (0, 1, 2, 4, 5, 6)]
    + [(i, 4) for i in (0, 1, 2, 3, 4, 5, 6, 7)]
)

# Rows that are structurally zero on input, for the columns that are live at
# all: neither a phi nor a theta perturbation changes the time or the q/p
# coordinate, in vacuum or in matter.
_B2F_ZERO_ROWS = {2: (3, 7), 3: (3, 7)}

# bound parameter count -- M is 8 x _B2F_COLS, stored column major
_B2F_COLS = 6

# The q/p column of M is not stored as the plain bound-to-free jacobian.  Its
# free rows are multiplied by q/p and divided by the column's own q/p row:
#
#     M[i, 4] = l * dFree_i/dl_0 / (dl/dl_0)     for i in 0..6
#     M[7, 4] =     dl/dl_0                      unchanged
#
# Both factors are storage convention only, and together they are what makes
# this column cost exactly what the ATLAS stepper's pVector[40] block costs:
# the l makes the term each RK stage picks up from H's own l dependence equal
# to the stage itself instead of a fresh cross product, and the normalisation
# makes the coefficient of that term one instead of a load of the q/p row.
# In vacuum both factors are constant across a step, so the invariant holds by
# itself; rk4_dense changes both and restores it explicitly.  The stepper
# converts to and from the plain jacobian where the covariance engine wants it.
_B2F_QOP_COLUMN = 4


def b2f_step_update(D, live, l_in=None, l_out=None):
    """Apply a free-to-free step jacobian D to the bound-to-free jacobian M.

    Only the live columns are touched, and the structurally zero rows are
    dropped from the input so the products stay sparse.  This is the generic
    fallback; the vacuum kernel instead folds the same operation into the RK
    recursion, which is cheaper still because it never builds D.

    The q/p column is stored scaled (see _B2F_QOP_SCALING).  A step that
    changes q/p changes both factors of that scaling, so it is undone on the
    way in, using the q/p the column was scaled with, and redone on the way
    out with the new one.  Pass l_in and l_out to do that; leaving them None
    treats every column as plain.
    """
    M = MatrixSymbol("M", 8, 6)
    qop_col = _B2F_QOP_COLUMN
    out = []
    for c in sorted({col for _, col in live}):
        zero = _B2F_ZERO_ROWS.get(c, ())
        v = Matrix([[0 if i in zero else M[i, c]] for i in range(8)])
        scaled = c == qop_col and l_in is not None
        if scaled:
            # stored = l_in * plain / plain[7], so plain = stored * stored[7] / l_in
            f = v[7, 0] / l_in
            v = Matrix([[v[i, 0] * f] if i < 7 else [v[i, 0]] for i in range(8)])
        new_v = D * v
        if scaled:
            g = l_out / new_v[7, 0]
            new_v = Matrix(
                [[new_v[i, 0] * g] if i < 7 else [new_v[i, 0]] for i in range(8)]
            )
        out.extend([new_v[i, 0]] for i, col in live if col == c)
    return Matrix(out)


def rk4_vacuum_b2f_atlasexpr(taylor_norm=False):
    """ATLAS-form RK4 vacuum step transporting the bound-to-free jacobian.

    The free-to-free step jacobian D is never built.  Instead the columns of
    the bound-to-free jacobian M are pushed through the very same RK recursion
    as the track state itself -- which is why the code below mirrors the value
    path line for line, exactly as ATLAS' d2A/d3A/d4A block mirrors its A0..A6
    block.  That removes both the fourth tangent direction (a bound
    parametrisation has two direction degrees of freedom, a free-to-free D
    needs three) and the 8x8 composition.
    """
    b = _Builder()
    st = _atlas_rk4_stages(b, taylor_norm)
    H0, H1, H2 = st.H0, st.H1, st.H2
    A0, A3, A4, A6 = st.A0, st.A3, st.A4, st.A6
    S3, dtds = st.S3, st.dtds

    M = MatrixSymbol("M", 8, 6)

    def col(c, rows):
        return Matrix([[M[i, c]] for i in rows])

    def tangent(name, stage, contribs, field=None):
        T = sym.zeros(3, 1) if field is None else field
        for var, Tvar in contribs:
            T = T + stage.expr.jacobian(var) * Tvar
        return b.add(name, T)

    def propagate(tag, c, with_field):
        """Push column `c` of M through the step.

        `with_field` is False for a column with no q/p component, which then
        picks up nothing from H's own l dependence.

        The q/p column is stored in the scaled form described at
        _B2F_QOP_SCALING: its free rows carry a factor l and are divided by the
        column's own q/p row.  Both factors are pure storage convention, but
        together they are what makes this column cost exactly what the ATLAS
        stepper's pVector[40] block costs.  The l turns the term each stage
        picks up from H's l dependence into the stage itself rather than a
        fresh cross product, and the normalisation makes the coefficient of
        that term a literal one instead of a load of the q/p row.  Nothing is
        assumed about the q/p row -- rk4_dense restores the invariant whenever
        it changes it.
        """
        seed = col(c, (4, 5, 6))
        pos_fac, dir_fac = S3.name, third.name
        if not with_field:
            fields = [None] * 4
        else:
            one = sym.Integer(1)
            fields = [
                _field_contrib(
                    b, f"{tag}0", A0, H0, _ex(A0.name) * one, _ex(H0.name) * one
                ),
                _field_contrib(
                    b, f"{tag}3", A3, H1, (_ex(A3.name) - d) * one, _ex(H1.name) * one
                ),
                _field_contrib(
                    b, f"{tag}4", A4, H1, (_ex(A4.name) - d) * one, _ex(H1.name) * one
                ),
                _field_contrib(
                    b, f"{tag}6", A6, H2, _ex(A6.name) * one, _ex(H2.name) * one
                ),
            ]

        v0 = tangent(f"{tag}0", A0, [(d, seed)], fields[0])
        v2 = b.add(f"{tag}2", _ex(v0.name) + seed)
        v3 = tangent(f"{tag}3", A3, [(d, seed), (st.A2.name, _ex(v2.name))], fields[1])
        v4 = tangent(f"{tag}4", A4, [(d, seed), (A3.name, _ex(v3.name))], fields[2])
        v5 = b.add(f"{tag}5", 2 * _ex(v4.name) - seed)
        v6 = tangent(f"{tag}6", A6, [(st.A5.name, _ex(v5.name))], fields[3])

        new_pos = col(c, (0, 1, 2)) + pos_fac * (
            _ex(v2.name) + _ex(v3.name) + _ex(v4.name)
        )
        new_dir = dir_fac * (
            _ex(v0.name) + 2 * _ex(v3.name) + _ex(v5.name) + _ex(v6.name)
        )
        return new_pos, new_dir

    # NOTE: dir_fac must be a Symbol, never a Number. sympy distributes a
    # Number over an Add, so Rational(1,3) * (v0 + 2*v3 + v5 + v6) emits four
    # multiplications per component where ATLAS' equivalent
    # ((d2A0 + 2*d2A3) + (d2A5 + d2A6)) * (1./3.) uses two -- eighteen wasted
    # multiplies once the three columns are counted. A Symbol does not
    # distribute, so the scaling stays factored and no temporary is needed.

    third = b.add("third", sym.Rational(1, 3))

    phi_pos, phi_dir = propagate("Tp", 2, False)
    the_pos, the_dir = propagate("Tt", 3, False)
    qop_pos, qop_dir = propagate("Tq", 4, True)

    # dt/ds depends on q/p only, so the time row moves for the q/p column alone.
    # The l is the scaling the stored column carries; the division by the q/p
    # row that goes with it has already cancelled, since the plain entry is
    # d(t)/d(q/p) times that row.
    dtdl = b.add("dtdl", l * dt_dqop(dtds.name))
    new_time = M[3, 4] + dtdl.name

    b.add("new_Mp", Matrix.vstack(phi_pos, phi_dir))
    b.add("new_Mt", Matrix.vstack(the_pos, the_dir))
    b.add("new_Mq", Matrix.vstack(qop_pos, Matrix([new_time]), qop_dir))

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

    new_B = name_expr("new_B", B3)

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

    # Unlike the vacuum kernel this one does build D, because energy loss
    # couples time and q/p into every stage and folding that into the RK
    # recursion buys little on a path that is not hot.  It still applies D to
    # the bound-to-free jacobian directly rather than to an 8x8 transport,
    # which is also what removes the composition-order question entirely.
    new_M = name_expr("new_M", b2f_step_update(D, _B2F_DENSE_LIVE, l, new_l.expr))

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
        new_B,
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
            "new_B",
            "new_p",
            "new_t",
            "new_d",
            "path_derivatives",
            "new_Mp",
            "new_Mt",
            "new_Mq",
        ]
    ]

    lines = []

    lines.append(
        "template <typename T, typename GetB>\n"
        "Acts::Result<bool> rk4_vacuum(std::span<const T, 3> p,"
        " std::span<const T, 3> d, const T t, const T h, const T l, const T m,"
        " const T p_abs, std::span<const T, 3> B1, GetB getB, T* err,"
        " const T errTol, std::span<T, 3> new_p, T* new_t,"
        " std::span<T, 3> new_d, std::span<T, 3> new_B,"
        " std::span<T, 8> path_derivatives, std::span<T> M) {"
    )
    lines.append(f"  assert(M.empty() || M.size() == {_B2F_COLS * 8});")

    def pre_expr_hook(var):
        if str(var) == "p2":
            return "std::array<T, 3> p2;"
        if str(var) == "p3":
            return "std::array<T, 3> p3;"
        for name, group in _B2F_LIVE_COLUMNS:
            if str(var) == name:
                return f"std::array<T, {len(group)}> {name};"
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
        for name, group in _B2F_LIVE_COLUMNS:
            if str(var) == name:
                return "\n".join(
                    f"M[{i + 8 * j}] = {name}[{k}];" for k, (i, j) in enumerate(group)
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
            "new_B",
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
        " const T q, const T p_abs, std::span<const T, 3> B1, GetB getB,"
        " GetG getG, T* err, const T errTol, std::span<T, 3> new_p, T* new_t,"
        " std::span<T, 3> new_d, T* new_l, std::span<T, 3> new_B,"
        " std::span<T, 8> path_derivatives, std::span<T> M) {"
    )
    lines.append(f"  assert(M.empty() || M.size() == {_B2F_COLS * 8});")

    lines.append("  const auto g1 = getG(p, l);")

    def pre_expr_hook(var):
        if str(var) == "p2":
            return "std::array<T, 3> p2;"
        if str(var) == "p3":
            return "std::array<T, 3> p3;"
        if str(var) == "l2":
            return "T l2[1];"
        if str(var) == "l3":
            return "T l3[1];"
        if str(var) == "l4":
            return "T l4[1];"
        if str(var) == "new_M":
            return f"std::array<T, {len(_B2F_DENSE_LIVE)}> new_M;"
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
                f"M[{i + 8 * j}] = new_M[{k}];"
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

# Against the free-to-free kernel this replaces (clang 17, -O2, arm64, timing
# SympyStepper::step with chained steps): ~3-5% faster with covariance
# transport, ~3-5% slower without.  The vacuum path is latency bound on the
# serial chain between the three field lookups, so pre-scaling the field adds
# one multiply to exactly that chain, which is what the no-covariance case pays.
#
# Two knobs measured and left off: taylor_norm=True, ATLAS' sqrt-free direction
# normalisation, trades a sqrt and a division for more multiplications at the
# same speed; run_cse=True is within noise and obscures the correspondence to
# the ATLAS code.
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
