@defgroup kf_material_effects Material effects in the Kalman filters
@ingroup track_fitting
@brief Pointwise material formalism shared by the KF and CKF

> [!tip]
> This page documents the concrete formalism used to apply material effects
> during Kalman filtering. For the high-level conceptual picture see
> @ref material-eff, and for the material description itself see @ref material.

Both the @ref Acts::KalmanFitter "Kalman Filter" (KF) and the
@ref Acts::CombinatorialKalmanFilter "Combinatorial Kalman Filter" (CKF) share
the exact same material-update machinery. Rather than integrating material
continuously along the trajectory, they use a **pointwise ("thin scatterer")
model**: all material assigned to a surface is collapsed onto a single point,
where it is applied as a deterministic correction to the momentum plus an
additive Gaussian process noise on the diagonal of the covariance matrix. This
matches the @f$\vec w_{k-1}@f$ process-noise term @f$\mathbf Q_{k-1}@f$ of the
@ref kalman-formalism "Kalman formalism" @cite Fruhwirth:1987fm. The
parametrizations used follow the ATLAS treatment of energy loss and multiple
scattering @cite Lund:2008ad.

The same code is used by the standalone @ref Acts::MaterialInteractor
propagator actor (used e.g. by plain propagation and the
@ref Acts::GaussianSumFitter "GSF"); the KF and CKF, however, do not install
that actor. Instead they call the underlying routine directly from their Kalman
actors so the material noise can be interleaved with the Kalman update at
precisely the right point on each surface.

## Where the update is applied {#kf-material-where}

On every surface the fitters apply the pointwise material update, using a
@ref Acts::MaterialUpdateMode that avoids double counting:

- On a **measurement surface** the interaction is split into a `PreUpdate`
  (applied *before* the Kalman filter update, after the covariance has been
  transported to the surface) and a `PostUpdate` (applied *after* the update),
  so that noise is added on both approach and exit.
- On a **passive / material-only surface** (including holes) a single
  `FullUpdate` is applied.
- At the start and target surface the mode is automatically restricted to
  `PostUpdate` / `PreUpdate` respectively.

Whether scattering and/or energy loss are applied is controlled by the
`multipleScattering` and `energyLoss` flags on the fitter options
(`Acts::KalmanFitterOptions`, `Acts::CombinatorialKalmanFilterOptions`), both
enabled by default.

Before the material is evaluated, the slab thickness is scaled by the surface
path correction @f$1/\cos\alpha@f$ to account for the incidence angle @f$\alpha@f$
of the trajectory.

## Multiple Coulomb scattering {#kf-material-scattering}

Scattering does not change the momentum magnitude but widens the direction
uncertainty. A single scattering standard deviation @f$\theta_0@f$ is computed
per surface. For all particles except electrons the **Highland formula**
@cite Highland:1975pq, in the parametrization of the Particle Data Group
@cite ParticleDataGroup:2018ovx (eq. 33.15), is used:

@f[
  \theta_0 = \frac{13.6\,\mathrm{MeV}}{\beta c\, p}\, q\,
             \sqrt{\frac{x}{X_0}}
             \left( 1 + 0.038\, \ln\!\left(\frac{x}{X_0}\frac{q^2}{\beta^2}\right) \right),
@f]

with the path length in radiation lengths @f$x/X_0@f$, the momentum @f$p@f$, the
velocity @f$\beta c@f$ and the charge @f$q@f$. For **electrons and positrons** the
Rossi–Greisen form (with a @f$17.5\,\mathrm{MeV}@f$ prefactor) is used instead.
The Highland term is evaluated as:

@snippet{trimleft} Interactions.cpp Highland theta0

The single angle @f$\theta_0@f$ is projected onto the two bound angular
parameters, giving the process-noise variances that are added to the covariance:

@f[
  \sigma^2(\theta) = \theta_0^2,
  \qquad
  \sigma^2(\phi) = \left(\frac{\theta_0}{\sin\theta}\right)^2 .
@f]

The @f$1/\sin\theta@f$ factor on @f$\phi@f$ accounts for the metric of the polar
parametrization. Note that only the diagonal @f$\sigma^2(\phi)@f$ and
@f$\sigma^2(\theta)@f$ entries are modified; no @f$\phi@f$–@f$\theta@f$ correlation
is introduced:

@snippet{trimleft} PointwiseMaterialInteraction.cpp scattering variance

## Energy loss {#kf-material-eloss}

Energy loss enters in two ways: a deterministic shift of the mean @f$q/p@f$ and an
additional variance on @f$q/p@f$.

**Mean energy loss (state update).** The mean ionization loss is computed from
the **Bethe formula** @cite ParticleDataGroup:2018ovx (eq. 33.5), including the
density-effect correction (eq. 33.6) and the maximum single-collision energy
transfer @f$W_\text{max}@f$ (eq. 33.4). The resulting energy loss @f$\Delta E@f$ is
applied to the particle energy,

@f[
  E' = \sqrt{m^2 + p^2} - \Delta E \cdot s,
  \qquad
  p' = \sqrt{E'^2 - m^2},
@f]

where @f$s = \pm 1@f$ is the propagation direction (energy decreases in the
forward direction, increases in the backward/smoothing direction). The updated
@f$q/p'@f$ is written back to the track state. A floor of @f$p' \geq 10\,\mathrm{MeV}@f$
is applied so that a too-large loss does not push the particle to negative
momentum:

@snippet{trimleft} PointwiseMaterialInteraction.hpp energy loss update

> [!note]
> The mean correction inside the KF/CKF uses the **ionization (Bethe) term
> only**. Radiative (bremsstrahlung / Bethe–Heitler) losses
> are *not* added to the mean here. Because bremsstrahlung is strongly
> non-Gaussian, it cannot be modelled adequately by this pointwise Gaussian
> update; electron fitting should therefore use the
> @ref Acts::GaussianSumFitter "GSF", which models the Bethe–Heitler
> distribution as a Gaussian mixture.

**Energy-loss straggling (variance).** The fluctuation of the ionization loss is
described by the **Landau–Vavilov** distribution. Its full width at half maximum
(@f$4\varepsilon@f$, @cite ParticleDataGroup:2018ovx fig. 33.7) is converted to an equivalent Gaussian
standard deviation @f$\sigma_E@f$ via
@f$\sigma_E = \mathrm{fwhm}/(2\sqrt{2\ln 2})@f$, which is then propagated to a
variance on @f$q/p@f$ through the Jacobian @f$\mathrm{d}(q/p)/\mathrm{d}E@f$:

@f[
  \sigma^2(q/p) = \left(\frac{\mathrm{d}(q/p)}{\mathrm{d}E}\right)^2 \sigma_E^2 .
@f]

This variance is added to the diagonal @f$\sigma^2(q/p)@f$ covariance entry:

@snippet{trimleft} PointwiseMaterialInteraction.cpp energy loss variance

## Covariance update {#kf-material-covariance}

The three variances computed above are added onto the corresponding **diagonal**
entries of the bound covariance, @f$\sigma^2(\phi)@f$, @f$\sigma^2(\theta)@f$ and
@f$\sigma^2(q/p)@f$. The sign is set by the
@ref Acts::NoiseUpdateMode -- noise is
*added* during the forward filtering pass. Each variance is floored at zero to
protect against numerical underflow. No off-diagonal correlations are created,
and the local position and time entries are untouched:

@snippet{trimleft} PointwiseMaterialInteraction.hpp covariance update

## Summary {#kf-material-summary}

| Effect | Formalism | Applied to |
| --- | --- | --- |
| Multiple scattering | Highland (Rossi–Greisen for @f$e^\pm@f$) | @f$\sigma^2(\phi)@f$, @f$\sigma^2(\theta)@f$ |
| Mean energy loss | Bethe (ionization) | mean @f$q/p@f$ |
| Energy-loss straggling | Landau width @f$\to@f$ Gaussian @f$\sigma@f$ | @f$\sigma^2(q/p)@f$ |
| Bremsstrahlung mean | — (use @ref Acts::GaussianSumFitter "GSF") | — |
| Correlations | none (diagonal noise only) | — |

## Implementation pointers {#kf-material-implementation}

- Physics formulas: @ref Acts::computeMultipleScatteringTheta0 (Highland /
  Rossi–Greisen), @ref Acts::computeEnergyLossBethe (mean ionization),
  @ref Acts::computeEnergyLossLandauSigmaQOverP (straggling) and
  @ref Acts::computeEnergyLossRadiative (radiative).
- Track fitters: @ref Acts::KalmanFitter, @ref Acts::CombinatorialKalmanFilter.
- Standalone propagator actor: @ref Acts::MaterialInteractor.
- Bethe–Heitler (GSF only): @ref Acts::AtlasBetheHeitlerApprox.

The pointwise application and covariance update themselves live in internal
(non-public) `detail` code; the snippets shown in the sections above are
extracted directly from it.
