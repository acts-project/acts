# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Transport Jacobian code generation
#
# This file has two lives. The build runs it as a plain Python script to
# generate `codegen/sympy_jac_math.hpp`, the header behind
# `Acts::detail::sympy::boundToBoundTransportJacobian` and
# `boundToCurvilinearTransportJacobian`. Opened in a jupytext-capable viewer
# (JupyterLab or VS Code with the jupytext extension, `jupytext --to ipynb`
# for anything else) the same file is a notebook that derives, shows and
# explains every matrix that goes into that header.
#
# The two never drift apart because there is only one file: the `# %%` markers
# are Python comments, so what the notebook shows is exactly what the build
# runs. The rendered form of this notebook is part of the ACTS documentation;
# see `docs/render_notebook.py` for how it gets there.
#
# ## How the build invokes this
#
# `codegen/manifest.json` is the single source of truth for every generated
# file. Its entry for this one says:
#
# ```json
# "Core/src/Propagator/codegen/sympy_jac_math.hpp": {
#     "script": "Core/src/Propagator/detail/generate_sympy_jac.py",
#     "output": "codegen/sympy_jac_math.hpp",
#     "with_requirements": ["codegen/requirements.txt"],
#     "with": ["codegen"],
#     "isolated": true,
#     "python_version": "3.13"
# }
# ```
#
# `cmake/ActsCodegen.cmake` reads that entry and runs the script in a pinned,
# isolated `uv` environment with the shared `codegen` package available,
# passing the destination path as `argv[1]`. `CI/pregenerate_codegen.py` reads
# the same entry to produce the file ahead of a build, which is how release
# source archives ship generated code and build without Python or network.
#
# To run it by hand exactly the way the build does:
#
# ```console
# $ uv run --isolated --python 3.13 --no-project \
#       --with-requirements codegen/requirements.txt --with ./codegen \
#       Core/src/Propagator/detail/generate_sympy_jac.py /dev/stdout
# ```
#
# To open it as a notebook, the only extra requirement is that the `codegen`
# package is importable:
#
# ```console
# $ uv run --with ./codegen --with jupyterlab --with jupytext \
#       jupyter lab Core/src/Propagator/detail/generate_sympy_jac.py
# ```

# %% [markdown]
# ## The problem
#
# A track is stored on a surface as six *bound* parameters with a $6\times6$
# covariance. Propagating that covariance from a surface $A$ to a surface $B$
# needs the Jacobian of the bound parameters at $B$ with respect to the bound
# parameters at $A$,
#
# $$ J_\text{full} = \frac{\partial q_B}{\partial q_A}, \qquad C_B = J_\text{full}\, C_A\, J_\text{full}^{T}. $$
#
# The stepper does not work in bound parameters, so $J_\text{full}$ is
# assembled out of four pieces (this is the standard result, spelled out in the
# [transport Jacobian correction white
# paper](https://acts.readthedocs.io/en/latest/white_papers/correction-for-transport-jacobian.html)):
#
# $$ J_\text{full} = J_{fb}\,\bigl(\mathbb{1}_8 + \tfrac{\partial r}{\partial s} \tfrac{\partial s}{\partial r}\bigr)\, J_t\, J_{bf} $$
#
# read right to left:
#
# | factor | shape | meaning |
# | --- | --- | --- |
# | $J_{bf}$ | $8\times6$ | bound parameters at $A$ to free parameters at $A$ |
# | $J_t$ | $8\times8$ | free-to-free transport from $A$ to $B$ |
# | $\mathbb{1} + \frac{\partial r}{\partial s}\frac{\partial s}{\partial r}$ | $8\times8$ | path-length correction at $B$ (below) |
# | $J_{fb}$ | $6\times8$ | free parameters at $B$ to bound parameters at $B$ |
#
# Doing that product with dense matrices costs $6\cdot8\cdot8 + \ldots$
# multiply-adds, and the overwhelming majority of them multiply a structural
# zero. All four factors are extremely sparse, and *which* entries are zero is
# known at compile time — it follows from the parametrisation, not from the
# track. That is what makes this worth generating: sympy is told the sparsity
# pattern of $M = J_t J_{bf}$ and the remaining factors; the printed C++
# contains only the operations that survive.
#
# The equivalent dense implementation is
# `Acts::detail::boundToBoundTransportJacobian` in
# `Core/src/Propagator/detail/JacobianEngine.cpp` — a useful reference when
# reading what follows, and the thing the generated code has to agree with.

# %% [markdown]
# ## Setup
#
# The imports, and the one piece of dual-mode machinery: whether we are a
# notebook or a build step. `ipykernel` is only ever imported by a Jupyter
# kernel, so its presence in `sys.modules` separates the two cases without any
# environment variable or flag.

# %%
import argparse
import sys
from pathlib import Path

import sympy as sym
from sympy import MatrixSymbol

# Notebook or build step? In a kernel we render; in the build we write the
# header cmake asked for, whose path arrives as argv[1].
IN_NOTEBOOK = "ipykernel" in sys.modules
if IN_NOTEBOOK:
    OUTPUT_PATH = None
    RUN_CHECKS = True
else:
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
    args = parser.parse_args()
    OUTPUT_PATH = Path(args.output) if args.output else None
    RUN_CHECKS = not args.no_check

try:
    import codegen.sympy_common  # noqa: F401
except ModuleNotFoundError:
    # Opened as a notebook without the shared package installed: fall back to
    # the copy in this checkout.
    for _candidate in [Path.cwd(), *Path.cwd().parents]:
        if (_candidate / "codegen" / "src" / "codegen").is_dir():
            sys.path.insert(0, str(_candidate / "codegen" / "src"))
            break

from codegen.sympy_common import (  # noqa: E402
    cxx_printer,
    my_cse,
    my_expression_print,
    name_expr,
)


def structure(matrix):
    """Sparsity pattern of `matrix`: 1 where an entry is not a literal zero.

    The interesting property of every matrix below is which entries are
    structurally zero, and that is much easier to see as a pattern than as a
    page of symbols.
    """
    return sym.Matrix(*matrix.shape, lambda i, j: 0 if matrix[i, j] == 0 else 1)


# %% [markdown]
# ## Parameter conventions
#
# Both parameter vectors are fixed by ACTS and used verbatim below.
#
# The **free** vector $r$ has eight components, indexed as in
# `Acts::FreeIndices`:
#
# | index | 0, 1, 2 | 3 | 4, 5, 6 | 7 |
# | --- | --- | --- | --- | --- |
# | | $x,\ y,\ z$ | $t$ | $d_x,\ d_y,\ d_z$ | $q/p$ |
#
# The **bound** vector $q$ has six, indexed as in `Acts::BoundIndices`:
#
# | index | 0, 1 | 2 | 3 | 4 | 5 |
# | --- | --- | --- | --- | --- | --- |
# | | $l_0,\ l_1$ | $\phi$ | $\theta$ | $q/p$ | $t$ |
#
# Note that the two orderings do not agree: free puts time third, bound puts
# it last. Several of the "1" entries in the matrices below are exactly that
# reshuffling.

# %% [markdown]
# ## $M$ — the transported bound-to-free Jacobian
#
# The stepper now transports the bound-to-free Jacobian directly:
# $M = J_t J_{bf}$. The generator receives this $8\times6$ matrix instead
# of multiplying the two factors again. Position depends on the initial local
# coordinates, angles and $q/p$; direction depends on angles and $q/p$;
# time depends on $q/p$ and the initial time. The $q/p$ diagonal is retained
# as an input, while the time column is exactly the free-time unit vector.
#
# The structural zeros are shared with the stepper and allow sympy to discard
# operations that cannot contribute.

# %%
M_dense = MatrixSymbol("M", 8, 6).as_explicit()
M = sym.zeros(8, 6)
M[0:3, 0:5] = M_dense[0:3, 0:5]
M[3, 4] = M_dense[3, 4]
M[4:7, 2:5] = M_dense[4:7, 2:5]
M[7, 4] = M_dense[7, 4]
M[3, 5] = 1

structure(M)

# %% [markdown]
# ## $J_{fb}$ — free to bound at the target surface
#
# The projection back to bound parameters, built by
# `Acts::Surface::freeToBoundJacobian`: local coordinates from the global
# position, $\phi$ and $\theta$ from the direction, and $q/p$ and $t$ passed
# through into their (differently numbered) bound slots.

# %%
J_fb_dense = MatrixSymbol("J_fb", 6, 8).as_explicit()

J_fb = sym.zeros(6, 8)
J_fb[0:2, 0:3] = J_fb_dense[0:2, 0:3]  # local coordinates vs. position
J_fb[2:4, 4:7] = J_fb_dense[2:4, 4:7]  # angles vs. direction
J_fb[5, 3] = 1  # bound time   <- free time
J_fb[4, 7] = 1  # bound q/p    <- free q/p

J_fb

# %% [markdown]
# ## The path-length correction
#
# The stepper stops at the target surface, not after a fixed path length. If
# the start parameters are varied, the varied trajectory reaches the surface
# after a *different* path length, and the free parameters have to be walked
# forward (or back) along the trajectory to the surface before being projected
# onto it. To first order that is
#
# $$ r \mapsto r + \frac{\partial r}{\partial s}\, \delta s, \qquad \delta s = \frac{\partial s}{\partial r}\, \delta r, $$
#
# which is the rank-one correction $\mathbb{1}_8 + \frac{\partial r}{\partial
# s}\frac{\partial s}{\partial r}$ sandwiched between $J_t$ and $J_{fb}$ in the
# product.
#
# **$\partial r / \partial s$**, `step_path_derivatives`, comes from the
# stepper (`SympyStepper::State::derivative`): position changes along the
# direction, time along $\mathrm{d}t/\mathrm{d}s$, direction along the
# curvature term.

# %%
step_path_derivatives = (
    MatrixSymbol("step_path_derivatives", 8, 1).as_explicit().as_mutable()
)
step_path_derivatives[7, 0] = 0  # q/p is constant along a vacuum step

step_path_derivatives.T

# %% [markdown]
# That last zero keeps the q/p row of the final Jacobian diagonal when
# combined with the sparsity of $M$. Its diagonal is not forced to one:
# the transported q/p derivative is preserved. A nonzero q/p path derivative
# would introduce further terms, outside the shape assumed by this generator
# and the covariance transport (asserted below).

# %% [markdown]
# **$\partial s / \partial r$**, `surface_path_derivatives`, is the derivative
# of the path length at which the trajectory intersects the target surface
# with respect to the free parameters. It comes from
# `Acts::Surface::freeToPathDerivative`. Neither the time nor $q/p$ moves the
# intersection point, so two of its eight entries are structurally zero.

# %%
surface_path_derivatives = (
    MatrixSymbol("surface_path_derivatives", 1, 8).as_explicit().as_mutable()
)
surface_path_derivatives[0, 3] = 0  # time does not move the intersection
surface_path_derivatives[0, 7] = 0  # neither does q/p

surface_path_derivatives

# %% [markdown]
# ## The generic product
#
# With $J_t J_{bf}$ already transported as $M$, the assembly is one line. `name_expr` pairs an
# expression with the symbol the generated code should assign it to — here a
# $6\times6$ `MatrixSymbol` called `J_full`, which becomes the output span
# in the C++ signature.

# %%
J_full_generic = name_expr(
    "J_full",
    J_fb * (sym.eye(8) + step_path_derivatives * surface_path_derivatives) * M,
)

structure(J_full_generic.expr)

# %% [markdown]
# Two pieces of structure survive the product. The time column is empty apart
# from its diagonal: time is carried along by the propagation and never fed
# back, so nothing but time itself depends on where the clock started. And the
# q/p row depends only on its own initial value. Its diagonal entry is
# taken from the transported matrix $M$, rather than fixed to `1`.
#
# Everywhere else the result is dense: each of the remaining bound parameters
# at the target ends up depending on all five live parameters at the start.
# The saving is not in the shape of the answer but in the cost of reaching it,
# and counting operations makes that point:

# %%
# Only computed when read as a notebook; the build has no use for it.
if IN_NOTEBOOK:
    print(f"raw expression: {sym.count_ops(J_full_generic.expr):>6} operations")
    print(
        "after CSE:      "
        f"{sum(sym.count_ops(e) for _, e in my_cse([J_full_generic])):>6} operations"
    )

# %% [markdown]
# ## The curvilinear special case
#
# Propagating to a curvilinear frame — the plane through the current position
# perpendicular to the direction — is common enough to be worth its own
# function. There the intersection condition is explicit: the plane is
# perpendicular to $\vec d$, so moving the free position by $\delta\vec r$
# changes the path length to the plane by $-\vec d \cdot \delta \vec r$, and
# nothing else contributes. $\partial s / \partial r$ collapses from six
# unknown entries to $(-\vec d^{\,T}, 0, 0, 0, 0, 0)$, which sympy then
# multiplies out.
#
# Everything else in the product is unchanged, so only this factor is
# rebuilt.

# %%
direction = MatrixSymbol("dir", 3, 1).as_explicit()

surface_path_derivatives_curvilinear = sym.Matrix.hstack(
    -direction.transpose(), sym.zeros(1, 5)
)

J_full_curvilinear = name_expr(
    "J_full",
    J_fb
    * (sym.eye(8) + step_path_derivatives * surface_path_derivatives_curvilinear)
    * M,
)

surface_path_derivatives_curvilinear

# %% [markdown]
# ## The shape the covariance transport assumes
#
# That structure is not just an observation: the *next* generator
# in the chain is built on them. `generate_sympy_cov.py`, which prints
# $C \mapsto J C J^{T}$ for the bound covariance, masks its `J_full` input down
# to exactly this shape before multiplying, so that the zeros cost nothing
# there either.
#
# Nothing links the two files, so if a change here made an entry outside that
# shape live, the covariance transport would silently drop it. Rather than
# leave that to chance, assert the shape at generation time: a violation stops
# the build with an error instead of quietly producing a wrong covariance.


# %%
def check_covariance_transport_sparsity(J_full):
    """Assert the shape the covariance transport masks this jacobian down to.

    Live entries outside it would be dropped there silently, so check here.

    Raises AssertionError if the jacobian reaches outside that shape.
    """
    J = sym.expand(J_full)
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


def check_curvilinear_is_generic_specialised():
    """Check the curvilinear expression against the specialised generic one."""
    spd = MatrixSymbol("surface_path_derivatives", 1, 8).as_explicit()
    at_curvilinear = {spd[0, i]: -direction[i, 0] for i in range(3)}
    at_curvilinear.update({spd[0, i]: 0 for i in range(3, 8)})
    diff = sym.expand(
        J_full_generic.expr.subs(at_curvilinear) - J_full_curvilinear.expr
    )
    if any(e != 0 for e in diff):
        raise AssertionError("curvilinear jacobian is not the generic one specialised")


if RUN_CHECKS:
    check_curvilinear_is_generic_specialised()
    check_covariance_transport_sparsity(J_full_generic.expr)

# %% [markdown]
# ## From expressions to C++
#
# The printing side is shared with every other ACTS generator and lives in
# `codegen/src/codegen/sympy_common.py`. `my_expression_print` does three
# things:
#
# 1. **common subexpression elimination** (`my_cse`): sympy's `cse` on the
#    *inflated* expression list, i.e. with matrices broken into scalars so
#    subexpressions can be shared across entries, then reassembled. This is
#    where the operation count above drops.
# 2. **ordering** (`order_exprs_by_output`): topologically sorts the
#    assignments so every temporary is defined before it is used, grouped by
#    the output that needs it.
# 3. **printing** (`cxx_printer`): a `CXX17CodePrinter` subclass that renders
#    matrix elements as flat column-major indices `var[i + j*rows]` — matching
#    Eigen's default storage, so the C++ side can pass spans over its data — and reorders sums so plain terms come first, which lets the
#    compiler contract them into fused multiply-adds.
#
# Temporaries are emitted as `const auto`, outputs are written through the
# span parameters. All that is left for a generator to supply is the
# function signature.


# %%
def function_source(signature, named_expr, run_cse=True):
    """Render one `NamedExpr` as the body of a C++ function."""
    body = my_expression_print(
        cxx_printer, [named_expr], [named_expr.name], run_cse=run_cse
    )
    indented = "\n".join(f"  {line}" for line in body.split("\n"))
    return f"{signature} {{\n{indented}\n}}"


GENERIC_SIGNATURE = (
    "template <typename T> void boundToBoundTransportJacobianImpl("
    "std::span<const T, 48> J_fb, std::span<const T, 48> M,"
    " std::span<const T, 8> step_path_derivatives,"
    " std::span<const T, 8> surface_path_derivatives,"
    " std::span<T, 36> J_full)"
)

CURVILINEAR_SIGNATURE = (
    "template <typename T> void boundToCurvilinearTransportJacobianImpl("
    "std::span<const T, 48> J_fb, std::span<const T, 48> M,"
    " std::span<const T, 8> step_path_derivatives,"
    " std::span<const T, 3> dir, std::span<T, 36> J_full)"
)

# %% [markdown]
# ## The generated header
#
# Assemble the file: licence header, include guard, and the two functions.

# %%
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

source = (
    HEADER
    + function_source(GENERIC_SIGNATURE, J_full_generic)
    + "\n"
    + function_source(CURVILINEAR_SIGNATURE, J_full_curvilinear)
    + "\n"
)

# %% [markdown]
# In the build this is written to the path cmake passed in. In a notebook
# there is nothing to write, so the header is simply shown — which is the
# whole point of reading this as a notebook: the code below is what the
# compiler actually sees.

# %%
if OUTPUT_PATH is not None:
    OUTPUT_PATH.write_text(source)
else:
    print(source, end="")
