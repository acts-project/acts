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

ACTS surface alignment parameters are local-frame parameters: translations
along the surface's own local axes, and rotations about those axes.

The adapter from these parameters to the local-frame parameters used below is
the identity. It is kept in one function so that a change of convention stays
local to it (see @ref alignment_composite_convention_adapter).

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

Both levels use local-frame parameters, which is also the ACTS surface
convention, so this matrix applies directly to ACTS alignment derivatives.

@anchor alignment_composite_convention_adapter
### Convention adapter

ACTS alignment parameters are the local-frame parameters
@f$(\delta\mathbf{t}, \delta\boldsymbol{\omega})@f$. They match the parameters
of the composite Jacobian, so

@f[
\frac{\partial\,\mathbf{a}_\mathrm{ACTS}}{\partial(\delta\mathbf{t}, \delta\boldsymbol{\omega})}
= I.
@f]

The derivative with respect to the composite parameters is therefore

@f[
\frac{\partial r}{\partial(\delta\mathbf{T}, \delta\mathbf{W})}
= \frac{\partial r}{\partial \mathbf{a}_\mathrm{ACTS}}
  \frac{\partial(\delta\mathbf{t}, \delta\boldsymbol{\omega})_s}
       {\partial(\delta\mathbf{T}, \delta\mathbf{W})_c}.
@f]

`localFrameToAlignmentParametersJacobian` returns that identity. It is the
single place that records the convention. While the center parameters were a
shift in global coordinates, it returned @f$\mathrm{diag}(R_s, 1)@f$.

## Implementation

`Acts/Surfaces/detail/AlignmentHelper.hpp` provides two free functions next to
`rotationToLocalAxesDerivative`:

```cpp
namespace Acts::detail {

/// d(dt, dw)_component / d(dT, dW)_composite, both in local-frame parameters
AlignmentMatrix compositeToComponentJacobian(
    const Transform3& compositeTransform, const Transform3& componentTransform);

/// d(ACTS alignment parameters) / d(local-frame parameters) = Identity
AlignmentMatrix localFrameToAlignmentParametersJacobian(
    const Transform3& surfaceTransform);

}  // namespace Acts::detail
```

`compositeToComponentJacobian` takes the local-to-global transforms.
`localFrameToAlignmentParametersJacobian` still takes the surface transform so
a later convention change can use it; the argument is unused while the
conventions agree. It is the only place that encodes the ACTS convention.

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
