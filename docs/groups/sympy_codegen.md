@defgroup sympy_codegen Symbolic code generation
@ingroup propagation
@brief How the @ref Acts::SympyStepper kernels are derived with sympy, and the ATLAS-form arrangement they use

> [!tip]
> This page is about how the stepping kernels are *produced*. For the stepper
> interface and the rest of the propagation machinery see @ref propagation.

The @ref Acts::SympyStepper does not contain a hand-written Runge-Kutta step.
Its two inner kernels — `rk4_vacuum` and `rk4_dense` — are emitted at build
time by `Core/src/Propagator/generate_sympy_stepper.py`, which derives them
symbolically with [sympy](https://www.sympy.org) and prints them as C++.

## Why generate them

A Runge-Kutta stepper has to produce two things per step: the new track state,
and the Jacobian of that state with respect to the old one. The second is
simply the derivative of the first, but written by hand the two drift apart —
a term added to the value path is easy to forget in the derivative path, and
the resulting bug is a covariance that is quietly a few percent wrong rather
than a test that fails.

Deriving the kernel symbolically removes that failure mode. The generator
states the equations of motion once, applies the Runge-Kutta recursion to
them, and lets sympy differentiate *the expressions that are actually
evaluated*. Structural zeros, constant folding and — where it pays —
common-subexpression elimination then happen on the symbolic form, before a
single line of C++ exists.

The same machinery generates the bound/free Jacobian helpers
(`generate_sympy_jac.py`) and the covariance transport
(`generate_sympy_cov.py`); the printer and the symbolic helpers they share
live in the `codegen` package at the repository root, which is also used by
the [detray](https://github.com/acts-project/detray) backend.

## How the build wires it up

`cmake/ActsCodegen.cmake` runs each generator through `uv` into a throw-away
environment pinned by `codegen/requirements.txt`, so a build never depends on
the ambient Python. `codegen/manifest.json` is the single place that records
which generator produces which header; releases ship the generated files
pre-built, and a build fully covered by them never invokes Python at all.

## The equations of motion

With @f$\vec r@f$ the position, @f$\vec T@f$ the unit direction,
@f$\lambda = q/p@f$, and @f$s@f$ the path length, the vacuum kernel integrates

@f[
  \frac{d\vec r}{ds} = \vec T, \qquad
  \frac{d\vec T}{ds} = \lambda\, \vec T \times \vec B, \qquad
  \frac{dt}{ds} = \frac{E}{p} = \sqrt{1 + \frac{m^2}{p^2}}, \qquad
  \frac{d\lambda}{ds} = 0 .
@f]

This is a second-order system in @f$\vec r@f$, which is what makes the
Runge-Kutta-Nyström form @cite Bugge:1981 the natural integrator. Stages 2 and
3 share the midpoint, so the field is sampled at three points per step — at
@f$\vec r@f$, at the midpoint, and near the end point — rather than four.

The dense kernel adds continuous energy loss, so @f$\lambda@f$ — and with it
@f$dt/ds@f$ — evolve along the step:

@f[
  \frac{d\lambda}{ds} = \frac{\lambda^2}{q}\,\frac{dt}{ds}\,g ,
  \qquad
  \frac{d}{ds}\!\left(\frac{dt}{ds}\right) = \frac{m^2\lambda^3}{q^3}\, g ,
@f]

with @f$g@f$ the mean energy loss per unit path length @cite Lund:2008ad.
Because @f$\lambda@f$ now enters every stage, the dense kernel cannot use the
pre-scaled arrangement below and builds the step Jacobian explicitly instead.

## The ATLAS arrangement

The vacuum kernel does not evaluate the plain Runge-Kutta slopes
@f$\vec k_i@f$. It follows the arrangement of the ATLAS `RungeKuttaPropagator`
— which @ref Acts::AtlasStepper transcribes directly — and carries a
**half-step bend vector** at each of the three field samples:

@f[
  \vec H_i = \frac{h\lambda}{2}\,\vec B_i .
@f]

Every stage slope then comes out already scaled by @f$h/2@f$, as a bare cross
product, and neither @f$h@f$ nor @f$\lambda@f$ appears again anywhere in the
recursion:

@f[
\begin{aligned}
  \texttt{kick1}         &= \vec T \times \vec H_1                  &&= \tfrac{h}{2}\vec k_1 \\
  \texttt{dir2}          &= \vec T + \texttt{kick1}                 &&= \vec T + \tfrac{h}{2}\vec k_1 \\
  \texttt{dir\_half\_sum} &= \vec T + \texttt{dir2}                  &&= 2\vec T + \tfrac{h}{2}\vec k_1 \\
  \texttt{dir3}          &= \vec T + \texttt{dir2} \times \vec H_2  &&= \vec T + \tfrac{h}{2}\vec k_2 \\
  \texttt{dir4}          &= \vec T + \texttt{dir3} \times \vec H_2  &&= \vec T + \tfrac{h}{2}\vec k_3 \\
  \texttt{dir\_end}      &= 2\,\texttt{dir4} - \vec T               &&= \vec T + h\,\vec k_3 \\
  \texttt{kick4}         &= \texttt{dir\_end} \times \vec H_3       &&= \tfrac{h}{2}\vec k_4
\end{aligned}
@f]

The two intermediate field sample points, the step result and the embedded
error estimate all fall out of these:

@f[
\begin{aligned}
  \vec r_2      &= \vec r + \tfrac{h}{4}\,\texttt{dir\_half\_sum}, \qquad
  \vec r_3      = \vec r + h\,\texttt{dir4}, \\
  \vec r'       &= \vec r + \tfrac{h}{3}\left(\texttt{dir2} + \texttt{dir3} + \texttt{dir4}\right), \\
  3\,\vec T'_{\text{un}} &= \texttt{kick1} + 2\,\texttt{dir3} + \texttt{dir\_end} + \texttt{kick4}, \\
  \varepsilon   &= 2|h| \left\|\,\texttt{dir\_half\_sum} + \texttt{kick4}
                   - \texttt{dir3} - \texttt{dir4}\,\right\|_1
                 = h^2 \left\|\vec k_1 - \vec k_2 - \vec k_3 + \vec k_4\right\|_1 .
\end{aligned}
@f]

> [!note]
> ATLAS also replaces the direction normalisation @f$1/\|3\vec T'_{\text{un}}\|@f$
> by a Taylor expansion around @f$3@f$, avoiding a square root and a division.
> The generator can emit that form (`taylor_norm=True`), but it is off: it
> trades the root for several multiplications on a chain that is not the
> binding one, and measured neutral to slower.

## Transporting the Jacobian

The kernel is handed the bound-to-free Jacobian @f$M@f$ and updates it in place,
never forming the 8×8 free-to-free step Jacobian @f$D@f$. Each live column is
pushed through *the same recursion as the state*, which is why the tangent block
mirrors the value block line for line — as ATLAS' `d2A`/`d3A`/`d4A` block
mirrors its `A0`..`A6`; see @ref sympy_codegen_naming for the two sets of
names.

Rows are the eight free parameters @f$x_i@f$ **now**, columns the six bound
parameters @f$b_j@f$ **at the start surface**,
@f$M_{ij} = \partial x_i / \partial b_j@f$:

|  | @f$l_0@f$ | @f$l_1@f$ | @f$\phi@f$ | @f$\theta@f$ | @f$\lambda@f$ | @f$t@f$ |
|---|---|---|---|---|---|---|
| @f$\vec r@f$ *(3 rows)* | hold | hold | step | step | step | · |
| @f$t@f$ | · | · | · | · | step | 1 |
| @f$\vec T@f$ *(3 rows)* | · | · | step | step | step | · |
| @f$\lambda@f$ | · | · | · | · | dense | · |

`·` is a structural zero, `1` a constant one, `hold` an entry no step writes,
`step` one every step writes, `dense` one only a dense step writes. The
generator declares the sparsity in this one table and reads its index sets back
off it.

### The @f$q/p@f$ column

@f$\lambda@f$ enters the recursion only through the bend vector
@f$\vec H = (h\lambda/2)\,\vec B@f$, so the @f$\lambda@f$ column is the only one
with a term from the field's own @f$\lambda@f$ dependence: one 3-vector at each
of the four stages that use a bend vector. Stored plainly, each of the four
carries a factor @f$M_{\lambda\lambda}@f$, and the column has to be scaled by
@f$\lambda@f$ into the recursion and unscaled out of it.

Storing it differentiated by @f$\log|\lambda|@f$ of the *current*
@f$\lambda@f$, with the @f$\lambda@f$ row kept plain, removes both:

@f[
  M_{i\lambda} \;\equiv\; \frac{\partial x_i}{\partial \log|\lambda|}
    \;=\; \lambda \, \frac{\partial x_i / \partial \lambda_0}
                          {\partial \lambda / \partial \lambda_0} ,
    \qquad i < 7 .
@f]

This is exact because nothing but @f$\lambda_0@f$ can change @f$\lambda@f$, so
@f$M_{\lambda\lambda}@f$ is the whole chain rule from
@f$\partial/\partial\lambda_0@f$ to @f$\partial/\partial\lambda@f$: the division
applies it, the factor carries on to the log, and the plain row inverts both.

It is cheaper because @f$\vec H@f$ is homogeneous of degree one in
@f$\lambda@f$, so @f$\partial\vec H/\partial\log|\lambda| = \vec H@f$. Each
stage's field term is then that stage's bend-linear part, with no
@f$M_{\lambda\lambda}@f$ factor and no scaling around the recursion: nineteen
multiplications and a division fewer per step, 395 floating-point operations
with covariance transport instead of 415, against ATLAS' 394.

That identity is the one place the recursion departs from the plain chain rule,
so the generator forms the chain-rule product as well and checks the two agree
(`Derivation.check_same`).

@f$\lambda@f$ and @f$M_{\lambda\lambda}@f$ are constant across a vacuum step, so
it carries the scaled form into itself; `rk4_dense` moves @f$\lambda@f$ and
converts explicitly. The stepper state holds the scaled form, and
`detail::sympy::toScaledBoundToFree` and its inverse convert where the
covariance engine wants the plain Jacobian. The convention is singular at
@f$\lambda = 0@f$, where the plain column already is.

ATLAS' `pVector[40]` block is the @f$M_{\lambda\lambda} = 1@f$ case: without
dense material the row stays one, and the block stays permanently scaled by
@f$\lambda@f$.

## Naming {#sympy_codegen_naming}

The kernels were originally transcribed with ATLAS' variable names, which are
positional rather than descriptive and are not defined in any published note —
they are Athena source convention. The generator now uses names that say what
the quantity is. The correspondence, for anyone reading the two side by side:

| generated name | ATLAS name | quantity |
|---|---|---|
| `pos`, `dir`, `qop` | `R`, `A`, `P[7]` (`CM`) | position, unit direction, @f$q/p@f$ |
| `half_h_qop` | `PS2` (@f$=@f$ `Pi`@f$\cdot h@f$, with `Pi`@f$=\lambda/2@f$) | @f$h\lambda/2@f$ |
| `bend1`, `bend2`, `bend3` | `H0`, `H1`, `H2` | @f$(h\lambda/2)\vec B_i@f$ at the three field samples |
| `kick1` | `A0`, `B0`, `C0` | @f$\tfrac{h}{2}\vec k_1@f$ |
| `dir2` | `A2`, `B2`, `C2` | @f$\vec T + \tfrac{h}{2}\vec k_1@f$ |
| `dir_half_sum` | `A1`, `B1`, `C1` | @f$2\vec T + \tfrac{h}{2}\vec k_1@f$ |
| `dir3` | `A3`, `B3`, `C3` | @f$\vec T + \tfrac{h}{2}\vec k_2@f$ |
| `dir4` | `A4`, `B4`, `C4` | @f$\vec T + \tfrac{h}{2}\vec k_3@f$ |
| `dir_end` | `A5`, `B5`, `C5` | @f$\vec T + h\vec k_3@f$ |
| `kick4` | `A6`, `B6`, `C6` | @f$\tfrac{h}{2}\vec k_4@f$ |
| `new_dir_x3` | — | @f$3\times@f$ the unnormalised new direction |
| `h_third`, `h_quarter`, `two_over_h` | `S3`, `S4`, `Sl` | @f$h/3@f$, @f$h/4@f$, @f$2/h@f$ |
| `dphi_*`, `dtheta_*`, `dqop_*` | `d2A*`, `d3A*`, `d4A*` | the tangent of the correspondingly named stage, one set per live column |
| `M` | `pVector[8..55]` | bound-to-free Jacobian, column major |
| `dEds`, `dEds1..4` | — | energy loss per unit path, per stage (@ref Acts::AtlasStepper has no material) |

Two things the rename fixes rather than preserves. ATLAS numbers its bend
vectors `H0`..`H2` against field samples `B1`..`B3`; here both are numbered
`1..3`. And the `A`-family index is not a stage index — `A1` is used *after*
`A2` — whereas `kick1`/`dir2`/`dir3`/`dir4`/`dir_end`/`kick4` say which stage
each belongs to.

## References

The Runge-Kutta-Nyström track model and the semi-analytic transport of the
derivatives alongside the trajectory come from @cite Bugge:1981. Adaptive
Runge-Kutta-Nyström step control and error estimation in ATLAS are studied in
@cite Lund:2009, and the extrapolation package the `RungeKuttaPropagator`
belongs to is described in @cite Salzburger:2007.

None of these notes define the variable names the Athena implementation uses —
they are source convention, and the table above is the only written record of
what they mean.
