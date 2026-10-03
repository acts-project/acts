@page alignment_composite_jacobians Composite alignment Jacobians

# Composite alignment Jacobians

This page documents the composite (hierarchy) Jacobians, which relate the
alignment parameters of a component (e.g. a sensor) to those of a rigid
structure it belongs to (e.g. a stave, a layer or a half-barrel). It covers the
parameter conventions they rely on, the derivation, the implementation, and what
is still missing.

[TOC]

## Motivation

ACTS computes derivatives of the residuals w.r.t. the alignment parameters of
single surfaces. Aligning every sensor individually does not scale to a real
tracker: the first steps of a real alignment are done on larger rigid bodies.

Aligning a rigid structure requires the chain rule from the sensor alignment
parameters, for which ACTS computes derivatives, to the structure parameters:

@f[
\frac{\partial r}{\partial \mathbf{A}_\mathrm{structure}}
  = \frac{\partial r}{\partial \mathbf{a}_\mathrm{sensor}}
    \cdot \frac{\partial \mathbf{a}_\mathrm{sensor}}{\partial \mathbf{A}_\mathrm{structure}}
@f]

The second factor is the *composite* or *hierarchy Jacobian*. The derivatives
of a structure are the sum of this product over all its sensors on the track.

## ACTS alignment parameter convention

The derivation depends on exactly what the ACTS surface alignment parameters
mean. This was checked in the code rather than assumed:

- `Acts::AlignmentIndices` (`Core/include/Acts/Definitions/Alignment.hpp`):
  `eAlignmentCenter0..2` are the "center of geometry object in global 3D
  cartesian coordinates", and `eAlignmentRotation0..2` are the "rotation angle
  around local x/y/z axis".
- `Surface::alignmentToBoundDerivativeWithoutCorrection`
  (`Core/src/Surfaces/Surface.cpp`) builds the translation block of
  @f$\partial(\text{local 3D position})/\partial\mathbf{a}@f$ from
  @f$-\hat{x}_\mathrm{local}^T, -\hat{y}_\mathrm{local}^T, -\hat{z}_\mathrm{local}^T@f$
  expressed in global coordinates. The translations are therefore **global**
  shifts of the surface center.
- The rotation block uses @f$(\mathbf{p} - \mathbf{c})@f$, the position relative
  to the surface center, so the rotation pivot is the surface center, i.e. the
  origin of its local frame.
- `detail::rotationToLocalAxesDerivative`
  (`Core/src/Surfaces/detail/AlignmentHelper.cpp`) documents the rotated frame
  as `compositeRotation * deltaRotation`, with @f$\Delta R@f$ built from small
  rotations about the **local** x, y, z axes. At first order
  @f$\Delta R = 1 + [\boldsymbol{\omega}]_\times@f$ with the standard
  cross-product matrix, so a local rotation vector @f$\boldsymbol{\omega}_l@f$
  corresponds to the global rotation vector
  @f$\boldsymbol{\omega}_g = R\,\boldsymbol{\omega}_l@f$.

In short, **current ACTS surface alignment parameters are mixed**: global-frame
translations of the center, and local-frame rotations about the center.

ACTS is expected to move to fully local-frame parameters. The implementation is
structured so that only one function changes when that happens (see
@ref alignment_composite_convention_adapter).

## Derivation

### Local-frame parameters

For an object with local-to-global transform @f$(R, \mathbf{c})@f$, the
*local-frame* alignment parameters are @f$(\delta\mathbf{t}, \delta\boldsymbol{\omega})@f$:
a translation along the object's own local axes, and small rotation angles about
its own local axes, pivoting on its local origin. The moved object is

@f[
R' = R\,(1 + [\delta\boldsymbol{\omega}]_\times), \qquad
\mathbf{c}' = \mathbf{c} + R\,\delta\mathbf{t},
@f]

i.e. @f$T' = T \cdot \Delta@f$ with @f$\Delta = (1 + [\delta\boldsymbol{\omega}]_\times,\ \delta\mathbf{t})@f$.
This is the convention used by CMS and Belle II for hierarchical alignment.

### Composite to component

Let the composite (structure) have transform @f$(R_c, \mathbf{c}_c)@f$ and
parameters @f$(\delta\mathbf{T}, \delta\mathbf{W})@f$, and the component
(sensor) have @f$(R_s, \mathbf{c}_s)@f$. A rigid motion of the composite moves
the global point @f$\mathbf{x}@f$ by
@f$R_c\,\delta\mathbf{T} + (R_c\,\delta\mathbf{W}) \times (\mathbf{x} - \mathbf{c}_c)@f$
and rotates by the global rotation vector @f$R_c\,\delta\mathbf{W}@f$.
With the lever arm @f$\mathbf{d} = \mathbf{c}_s - \mathbf{c}_c@f$ (global
coordinates), the component therefore moves by

@f[
\delta\mathbf{c}_s = R_c\,\delta\mathbf{T} - [\mathbf{d}]_\times R_c\,\delta\mathbf{W},
@f]
@f[
R_s' = (1 + [R_c\,\delta\mathbf{W}]_\times)\,R_s
     = R_s\,(1 + [R_s^T R_c\,\delta\mathbf{W}]_\times).
@f]

Expressing @f$\delta\mathbf{c}_s@f$ along the component axes
(@f$\delta\mathbf{t}_s = R_s^T\,\delta\mathbf{c}_s@f$) and reading off
@f$\delta\boldsymbol{\omega}_s = R_s^T R_c\,\delta\mathbf{W}@f$ gives

@f[
\frac{\partial(\delta\mathbf{t}, \delta\boldsymbol{\omega})_s}
     {\partial(\delta\mathbf{T}, \delta\mathbf{W})_c}
= \begin{pmatrix}
    R_\mathrm{rel}^T & -R_s^T\,[\mathbf{d}]_\times\,R_c \\
    0                & R_\mathrm{rel}^T
  \end{pmatrix},
\qquad R_\mathrm{rel} = R_c^T R_s .
@f]

### Relation to the CMS form

The CMS form of this matrix is

@f[
\begin{pmatrix} R^{-1} & (\mathbf{T} \times \mathbf{a}_s)\cdot\mathbf{b} \\ 0 & R^{-1} \end{pmatrix},
@f]

with @f$R@f$ the component rotation relative to the composite, @f$\mathbf{T}@f$
the component position relative to the composite, @f$\mathbf{a}_s@f$ the
component axes and @f$\mathbf{b}@f$ the composite axes. The two agree:

- the diagonal blocks are both @f$R_\mathrm{rel}^T = R_s^T R_c@f$;
- element @f$(i, j)@f$ of the top-right block is
  @f$-\hat{a}_i \cdot (\mathbf{d} \times \hat{b}_j) = \hat{b}_j \cdot (\mathbf{d} \times \hat{a}_i)
  = (\mathbf{d} \times \hat{a}_i)\cdot\hat{b}_j@f$, the CMS expression with
  @f$\mathbf{T} = \mathbf{d}@f$, including the sign.

This form is only correct when **both** levels use local-frame parameters. It
does not apply directly to current ACTS surface parameters, whose translations
are global. Porting CMS or Belle II code without accounting for this gives a
wrong top-left block. On a pure translation test with @f$R_c \approx 1@f$ that
mistake can look correct.

@anchor alignment_composite_convention_adapter
### Convention adapter

The ACTS parameters of a surface relate to its local-frame parameters by
@f$\delta\mathbf{c} = R_s\,\delta\mathbf{t}@f$, with the rotations unchanged:

@f[
\frac{\partial\,\mathbf{a}_\mathrm{ACTS}}{\partial(\delta\mathbf{t}, \delta\boldsymbol{\omega})}
= \begin{pmatrix} R_s & 0 \\ 0 & 1 \end{pmatrix}.
@f]

The Jacobian used on the derivatives ACTS provides is therefore

@f[
\frac{\partial\,\mathbf{a}_\mathrm{ACTS}}{\partial(\delta\mathbf{T}, \delta\mathbf{W})}
= \begin{pmatrix} R_s & 0 \\ 0 & 1 \end{pmatrix}
  \begin{pmatrix} R_\mathrm{rel}^T & -R_s^T [\mathbf{d}]_\times R_c \\ 0 & R_\mathrm{rel}^T \end{pmatrix}
= \begin{pmatrix} R_c & -[\mathbf{d}]_\times R_c \\ 0 & R_s^T R_c \end{pmatrix}.
@f]

The adapter is only needed when chaining to structure parameters:

@f[
\frac{\partial r}{\partial(\delta\mathbf{T}, \delta\mathbf{W})}
= \frac{\partial r}{\partial \mathbf{a}_\mathrm{ACTS}}
  \begin{pmatrix} R_s & 0 \\ 0 & 1 \end{pmatrix}
  \frac{\partial(\delta\mathbf{t}, \delta\boldsymbol{\omega})_s}
       {\partial(\delta\mathbf{T}, \delta\mathbf{W})_c}.
@f]

An alignment of single surfaces uses the derivatives
@f$\partial r/\partial \mathbf{a}_\mathrm{ACTS}@f$ as they are: the solver then
returns ACTS parameters, which is the convention corrections are applied in, and
no adapter enters. Applying the adapter to the surface derivatives alone (the
surface as its own structure) would instead give parameters in the surface's
local frame.

When ACTS moves to local-frame parameters, the adapter becomes the identity and
the composite Jacobian is used unchanged.

The alternative of carrying two composite Jacobians, one per convention, was
rejected. The structure-level maths is identical in both conventions; only the
interpretation of the sensor derivatives changes. With two versions, every
alignment parameter, every solver result and every correction applied back to
the geometry would have two possible meanings, and mixing them up gives silent
sign or frame errors.

## Implementation

`Acts/Surfaces/detail/AlignmentHelper.hpp` provides two free functions next to
`rotationToLocalAxesDerivative`:

```cpp
namespace Acts::detail {

/// d(dt, dw)_component / d(dT, dW)_composite, both in local-frame parameters
AlignmentMatrix compositeToComponentJacobian(
    const Transform3& compositeTransform, const Transform3& componentTransform);

/// d(ACTS alignment parameters) / d(local-frame parameters) = diag(R, 1)
AlignmentMatrix localFrameToAlignmentParametersJacobian(
    const Transform3& surfaceTransform);

}  // namespace Acts::detail
```

Both take local-to-global transforms. `localFrameToAlignmentParametersJacobian`
is the only place that encodes the ACTS convention.

## Limitations and open items

- **Structure frames.** The Jacobian needs the frame of the composite: its
  origin (the pivot of the rotations, entering through the lever arm) and its
  orientation @f$R_c@f$. The choice is up to the user; a common one is the
  global orientation centred on the structure, so that the translations are
  global shifts.
- **Several levels.** The Jacobians chain: for a sensor in a stave in a layer,
  the Jacobian w.r.t. the layer is the product of the sensor-to-stave and
  stave-to-layer Jacobians, both in local-frame parameters. No helper for
  multi-level chains is provided yet.
- **Degrees of freedom at several levels.** When both a structure and its
  components float, the structure motion is degenerate with the common motion
  of its components. The usual way out (CMS, Belle II) is 6 linear constraints
  per structure, @f$\sum_{s \in \mathrm{structure}} J_s^T\, \mathbf{a}_s = 0@f$,
  so that the component corrections carry no net rigid-body motion of their
  structure. Since the component and structure parameters are in different
  frames, a component degree of freedom maps to a linear combination of
  structure ones; the degeneracy is not visible in the alignment masks alone.
- **Linearisation point.** The Jacobians depend on the transforms they are
  evaluated with. Evaluating them on the nominal geometry is fine for small
  misalignments; an iterative alignment with large corrections would have to
  re-evaluate them on the updated alignment.
