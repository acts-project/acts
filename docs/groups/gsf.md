@defgroup gsf The Gaussian Sum Filter
@ingroup track_fitting
@brief Multi-component fitter for electrons with non-Gaussian energy loss

> [!tip]
> This page documents the @ref Acts::GaussianSumFitter "GSF" as implemented in
> ACTS today. For the pointwise multiple-scattering and ionization formalism that
> the GSF shares with the ordinary Kalman filters, see @ref kf_material_effects;
> for the high-level conceptual picture see @ref material-eff. The reference
> implementation and the tuning studies quoted below are described in detail in
> @cite Huth:2024.

## Why a mixture of Gaussians?

The @ref Acts::KalmanFitter "Kalman Filter" (KF) is the optimal estimator only as
long as every involved distribution is Gaussian. For electrons this breaks down:
their dominant energy loss is **bremsstrahlung**, whose Bethe–Heitler
distribution is strongly non-Gaussian and heavily tailed. Feeding that through
the single-Gaussian pointwise update of the KF (see the
@ref kf-material-eloss "energy-loss note") biases the momentum estimate and
wrecks its error estimate.

The **Gaussian Sum Filter** (GSF) @cite Fruhwirth:1997 addresses this by
modelling the track state as a *weighted mixture* of Gaussians instead of a
single one,

@f[
  p(\vec x) = \sum_i^{N_c} w_i\, \mathcal N(\vec x \mid \vec\mu_i, \mathbf\Sigma_i),
  \qquad \sum_i^{N_c} w_i = 1 ,
@f]

and running, in effect, one Kalman filter per component. Each component carries a
weight, a bound parameter vector and a bound covariance:

@snippet{trimleft} GsfComponent.hpp gsf component

> [!note]
> The GSF is substantially more expensive than the KF and is therefore only run
> when a track is likely to be an electron — typically to re-fit a silicon track
> that has been associated with an electromagnetic-calorimeter cluster
> @cite Huth:2024. It is **not** a track-finding algorithm: the measurement
> sequence is taken as given.

## Bethe–Heitler energy loss as a mixture {#gsf-bethe-heitler}

On a material surface the GSF replaces the KF's deterministic ionization loss and
Landau straggling on @f$q/p@f$ with an explicit mixture model of the
Bremsstrahlung loss. The Bethe–Heitler probability density of the energy
retention @f$z = E_f/E_i@f$ depends only on the traversed thickness in radiation
lengths @f$x/X_0@f$,

@f[
  f(z) = \frac{[-\ln z]^{c-1}}{\Gamma(c)}, \quad 0 \le z \le 1,
  \qquad c = \frac{x/X_0}{\ln 2} .
@f]

This is approximated by a 1D Gaussian mixture in @f$z@f$,
@f$f(z) \approx \sum_n^{N_{bh}} \pi_n\,\mathcal N(z \mid \mu_{z,n}, \sigma_{z,n})@f$.
For a **single** component the natural choice keeps the first two moments of the
exact distribution @cite Huth:2024,

@f[
  \mu_z = e^{-t}, \qquad \sigma_z^2 = 3^{-c} - 4^{-c}, \qquad t = x/X_0 ,
@f]

which is what @ref Acts::BetheHeitlerApproxSingleCmp evaluates:

@snippet{trimleft} BetheHeitlerApprox.hpp single component moments

A single Gaussian reflects the true, tailed distribution very poorly, so in
practice a multi-component approximation is used. Because the mixture cannot be
derived in closed form, its weights, means and variances are pre-fitted
(minimising either the Kullback–Leibler divergence or the CDF distance to
@f$f(z)@f$) and stored as polynomials in @f$x/X_0@f$ so they can be interpolated at
run time. Any approximation is accessed through the abstract interface:

@snippet{trimleft} BetheHeitlerApprox.hpp bethe heitler interface

`mixture()` writes @f$N_{bh}@f$ one-dimensional components (weight, mean, variance
in @f$z@f$) into a caller-provided span:

@snippet{trimleft} BetheHeitlerApprox.hpp gaussian component

The concrete @ref Acts::PolynomialBetheHeitlerApprox implements the polynomial
form; the default parametrisation shipped in the source
(@ref Acts::makeDefaultBetheHeitlerApprox) and the reference JSON configuration
used by the examples (`betheHeitler_geantSim_cdf_nC6_O5.json`) both use a
**six-component, fifth-order** CDF fit, split into a low- and a high-thickness
range at @f$x/X_0 = 0.1@f$ @cite Huth:2024. When a surface exceeds the valid
@f$x/X_0@f$ range the fitter counts the occurrence and emits a warning; the
deprecated alias `Acts::AtlasBetheHeitlerApprox` remains for the ATLAS/Athena
`.par` data format.

Applying the loss convolves every track-state component with every Bethe–Heitler
component, so a mixture of @f$N_c@f$ components becomes @f$N_c \cdot N_{bh}@f$. In
the backward pass an effective energy *gain* is applied instead.

## The algorithm on a surface {#gsf-algorithm}

The GSF actor (`Acts::detail::Gsf::GsfActor`) drives the fit as a propagator
actor. When
the multi-stepper reports the state on a surface, the actor executes the
following, in this order:

1. Transport each component's covariance onto the surface.
2. If the surface has material, apply **multiple scattering** as `PreUpdate`
   (measurement surface) or `FullUpdate` (passive surface) — this reuses exactly
   the pointwise machinery documented in @ref kf-material-scattering.
3. Perform the **Kalman update**: run the measurement update for every component
   and re-weight them. Up to a normalisation, a component's weight is scaled by
   the likelihood of the measurement given that component,
   @f[
     w_{k|k}^i \;\propto\; w_{k|k-1}^i \,
       \mathcal N\!\bigl(m_k \mid \mathbf H_k \vec x_{k|k-1}^i,\;
                          \mathbf V_k + \mathbf H_k \mathbf\Sigma_{k|k-1}^i \mathbf H_k^{\mathsf T}\bigr),
   @f]
   i.e. components incompatible with the measurement are exponentially
   suppressed. On a passive surface a no-measurement update is done instead
   (which may flag a hole).
4. Apply the **Bethe–Heitler** convolution of @ref gsf-bethe-heitler to every
   component, expanding the mixture.
5. **Reduce** the mixture back down (see @ref gsf-reduction) and drop components
   below the weight cutoff.
6. Push the reduced mixture back into the stepper, then apply the `PostUpdate`
   scattering on measurement surfaces.

> [!note]
> The measurement update (step&nbsp;3) always happens **before** the
> Bremsstrahlung loss (step&nbsp;4): the update must use the state estimated from
> the material seen between the *previous* and the *current* surface. The mixture
> reduction (step&nbsp;5) runs only *after* the loss, because the Kalman update
> cannot increase the component count — only the convolution can.

Concretely, the re-weighting multiplies each component's weight by
@f$\sqrt{1/\det R}\,\exp(-\tfrac12\chi^2)@f$ (the smallest @f$\chi^2@f$ over the
components is factored out for numerical stability, and the weights are
normalised afterwards):

@snippet{trimleft} GsfUtils.hpp posterior weights

## Mixture reduction {#gsf-reduction}

Left unchecked, the component count would grow by a factor @f$N_{bh}@f$ per
material surface. To keep it bounded, a reducer is invoked after each convolution
to bring the mixture down to `min(stepper.maxComponents, maxComponents)`
components. The reducer is a delegate, so it can be swapped out:

@snippet{trimleft} GsfOptions.hpp mixture reducer

Two production reducers are provided (a third, `…Naive`, is a
reference/benchmark baseline):

- @ref Acts::reduceMixtureLargestWeights — simply discards the lowest-weight
  components. Fast, but loses information.
- @ref Acts::reduceMixtureWithKLDistance — greedily merges the pair of components
  with the smallest symmetric Kullback–Leibler distance (evaluated on the
  @f$q/p@f$ dimension only) until the target count is reached. Slower, but
  markedly better, and the recommended choice @cite Huth:2024.

The pairwise distance driving that greedy merge is the symmetric KL divergence,
restricted to the @f$q/p@f$ dimension:

@snippet{trimleft} GsfComponentMerging.cpp kl divergence

## Reducing a mixture to one estimate {#gsf-merging}

Several steps — storing an intermediate state, and producing the final fitted
parameters — need a *single* parameter vector and covariance from a mixture. The
method is selectable via the @ref Acts::ComponentMergeMethod enum:

@snippet{trimleft} GsfOptions.hpp component merge method

- `eMean` keeps the first two moments (weighted mean and covariance of the
  mixture). The mean is a poor summary of a tailed distribution and can bias the
  result.
- `eMaxWeight` (the default) takes the parameters of the highest-weight
  component as the point estimate while still reporting the full mixture
  covariance. When one component dominates, this approximates the mode well and
  avoids the bias.

Merging must respect cyclic bound coordinates. Which coordinates are cyclic
depends on the surface type, encoded as compile-time angle descriptions (note
that on a cylinder the local @f$R\phi@f$ coordinate is cyclic, scaled by the
radius):

@snippet{trimleft} GsfComponentMerging.hpp angle description

The mean itself is then formed with complex-phase arithmetic — each cyclic
coordinate is mapped onto the unit circle, averaged as a complex number, and
converted back with `std::arg` — so that angles wrap correctly:

@snippet{trimleft} GsfComponentMerging.hpp circular mean

## Multi-component transport {#gsf-multistepper}

Each component must be transported individually, so the GSF runs on the
@ref Acts::MultiEigenStepperLoop rather than the single-component stepper. The
navigator, however, must see a *single* trajectory. The stepper therefore
presents a reduced representation to the navigation, configurable through the
reducer type; the default is the highest-weight component
(`Acts::MaxWeightReducerLoop`, with `Acts::MaxMomentumReducerLoop` as an
alternative), which keeps the navigation stream close to the bulk of the mixture.

Determining when the whole multi-component state has "reached" a surface is
handled by @ref Acts::MultiStepperSurfaceReached, which by default treats the
state as on-surface once its *average* is within tolerance. This guards against a
pathology described in @cite Huth:2024 — low-momentum components approaching a
cylinder on a straight-line intersection can spiral indefinitely while always
reporting *reachable*. A step limit that engages once the first component lands
on the surface (`stepLimitAfterFirstComponentOnSurface`, default 50) forces the
remaining stragglers to *unreachable* and removes them, after which the weights
are renormalised:

@snippet{trimleft} MultiStepperLoop.ipp step limit

## Forward/backward passes and output {#gsf-passes}

@ref Acts::GaussianSumFitter is constructed from a propagator, a shared
Bethe–Heitler approximation and a logger, and exposes two `fit` overloads: one
for the standard @ref Acts::Navigator, and one taking an explicit surface
sequence for use with the @ref Acts::DirectNavigator — the latter is the
re-fitting entry point used in the electron workflow above.

A fit runs a forward pass from the start parameters, then a backward pass that
starts from the last measurement with its covariance inflated by
`reverseFilteringCovarianceScaling` (default 100) and targets the reference
surface. Measurement surfaces that were seen going forward but not on the way
back are flagged as outliers. The multi-component state is merged
(@ref gsf-merging) into the single set of parameters that downstream algorithms
expect; the full final mixture can optionally be attached to the track.

> [!note]
> The current implementation stores only the **means** of the per-surface states
> in the @ref Acts::MultiTrajectory, so it performs **no** dedicated component
> smoothing of the kind originally described for the GSF — the backward pass
> plays the role of the smoother.

## Configuration and tuning {#gsf-configuration}

The knobs on @ref Acts::GsfOptions trade physics performance against runtime. The
values below summarise the scan in @cite Huth:2024; the ACTS example chain uses
12 components, KL-distance reduction, `eMaxWeight` merging and a weight cutoff of
@f$10^{-4}@f$.

| Option | Effect | Guidance @cite Huth:2024 |
| --- | --- | --- |
| `maxComponents` | mixture size after each reduction | runtime grows @f$\approx@f$ quadratically; physics plateaus beyond @f$\sim 12@f$ (library default 4, example default 12) |
| `weightCutoff` | discard components below this weight | @f$10^{-4}@f$ is a good default; @f$0.1@f$ is too aggressive (fit failures spike) |
| `mixtureReducer` | reduction algorithm | KL-distance clearly beats the weight cut at modest extra cost |
| `componentMergeMethod` | mixture → single estimate | `eMaxWeight` avoids the @f$q/p@f$ bias seen with `eMean` |
| Bethe–Heitler approx | mixture model of the loss | 6-component CDF polynomials, split at @f$x/X_0=0.1@f$ |
| `reverseFilteringCovarianceScaling` | covariance inflation for the backward pass | default 100 (not tuned for all setups) |
| `disableAllMaterialHandling` | switch off convolution and scattering | debugging only |

The payoff: against the KF, the 12-component GSF turns a heavily one-sided
@f$q/p@f$ residual into a near-symmetric one and shrinks its width, while the KF's
@f$q/p@f$ pull — its error estimate — is badly distorted by the non-Gaussian loss
@cite Huth:2024. A single-component GSF (equivalent to a KF using the Bethe–Heitler
mean and variance) is visibly biased, which is what motivates the mixture in the
first place.

## Implementation pointers {#gsf-implementation}

- Fitter and options: @ref Acts::GaussianSumFitter, @ref Acts::GsfOptions,
  @ref Acts::GsfComponent, @ref Acts::ComponentMergeMethod.
- Bethe–Heitler approximation: @ref Acts::BetheHeitlerApprox,
  @ref Acts::PolynomialBetheHeitlerApprox,
  @ref Acts::makeDefaultBetheHeitlerApprox.
- Mixture reduction: @ref Acts::reduceMixtureWithKLDistance,
  @ref Acts::reduceMixtureLargestWeights.
- Multi-component transport: @ref Acts::MultiEigenStepperLoop,
  @ref Acts::MultiStepperSurfaceReached.
- Shared material formalism (scattering / ionization): @ref kf_material_effects.

The per-surface algorithm itself lives in the internal `Acts::detail::Gsf`
code; the snippets above are extracted directly from the corresponding headers.
