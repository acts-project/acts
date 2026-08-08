@defgroup track_fitting Track Fitting
@brief Algorithms and data structures for track fitting

The track fitters estimate the track parameters from a fixed set of
measurements. Material effects (multiple scattering and energy loss) are applied
on each surface using a shared pointwise formalism; see
@ref kf_material_effects for the details of that formalism and its
implementation.

Currently, there are implementations for three different fitters:

- Kalman Filter (KF)
- Gaussian Sum Filter (GSF)
- Global Chi-Square Fitter (GX2F)

Even though all of them are least-squares fits, the concepts are quite
different, so identical results should not be expected from all of them.

## Kalman Filter (KF)

The Kalman Filter is an iterative fitter. It successively combines
measurements to obtain an estimate of the track parameters. The KF needs an
estimate as a starting point, and the procedure alternates between two steps:
extrapolating the current state to the next surface, and updating that
extrapolation using the measurement on the new surface. Depending on which
data is used, the KF gives three different interpretations of the track
parameters at a given surface:

- *predicted*: an extrapolation from the previous surfaces onto the current one.
- *filtered*: the predicted parameters updated with the measurement on the
  current surface, effectively a weighted mean.
- *smoothed*: an estimate that also incorporates information from later
  surfaces, only available once the full propagation has completed once.

See @ref kalman-formalism "the Kalman formalism" for the full mathematical
treatment of prediction, filtering and smoothing, and @ref Acts::KalmanFitter
for the concrete implementation.

## Gaussian Sum Filter (GSF)

The GSF is an extension of the Kalman Filter that allows handling non-Gaussian
errors by modelling the track state as a Gaussian mixture:

@f[
p(\vec{x}) = \sum_i w_i \varphi(\vec{x}; \mu_i, \Sigma_i), \quad \sum_i w_i = 1
@f]

A common use case for this is electron fitting. The energy loss from
Bremsstrahlung for electrons in matter is highly non-Gaussian, and thus cannot
be modelled accurately by the pointwise Gaussian material formalism used by
the KF (see @ref kf_material_effects). Instead, the Bremsstrahlung is modelled
as a Bethe-Heitler distribution, where @f$z@f$ is the fraction of the energy
remaining after the interaction (@f$E_f/E_i@f$), and @f$t@f$ is the material
thickness in terms of the radiation length:

@f[
f(z) = \frac{(- \ln z)^{c-1}}{\Gamma(c)}, \quad c = t/\ln 2
@f]

@anchor fig_bethe_heitler

![The true Bethe-Heitler distribution compared with a Gaussian mixture approximation (in thin lines the individual components are drawn) at t = 0.1 (corresponds to ~10mm Silicon).](track_fitting/gsf_bethe_heitler_approx.svg) {width=450px}

To be able to handle this with the Kalman filter mechanics, this distribution
is approximated by a Gaussian mixture as well (see @ref fig_bethe_heitler
"the figure above"). The GSF algorithm then works as follows (see also
@ref fig_gsf_overview "the figure below"):

- On a surface with material, the Bethe-Heitler energy-loss distribution is
  approximated with a fixed number of Gaussian components for each existing
  component. Since the number of components would otherwise grow
  exponentially with each material interaction, components that are close in
  terms of their *Kullback-Leibler divergence* are merged to limit the
  computational cost.
- On a measurement surface, a Kalman update is performed for each component.
  Afterwards, the component weights are corrected according to each
  component's compatibility with the measurement.

@anchor fig_gsf_overview

![Simplified overview of the GSF algorithm.](track_fitting/gsf_overview.svg) {width=450px}

### The Multi-Stepper

To implement the GSF, a special stepper is needed that can handle a
multi-component state internally: the class template @ref
Acts::MultiStepperLoop can extend any valid single-component stepper (e.g.,
@ref Acts::EigenStepper or @ref Acts::SympyStepper) to a multi-component
stepper. It interfaces to the navigation as one aggregate state to limit the
navigation overhead, but internally processes a multi-component state. How
this aggregation is performed can be configured via a template parameter; by
default the maximum weight is used (@ref Acts::MaxWeightReducerLoop).

Even though the multi-stepper interface exposes only one aggregate state and
thus is compatible with most standard tools, a special aborter is required to
stop the navigation when the surface is reached: @ref
Acts::MultiStepperSurfaceReached. It checks if all components have already
reached the target surface and updates their state accordingly. Optionally,
it can also stop the propagation when the aggregate state reaches the
surface.

### Using the GSF

The GSF is implemented in @ref Acts::GaussianSumFitter. The interface of its
`fit(...)` functions is very similar to the one of @ref Acts::KalmanFitter
(one overload for the standard @ref Acts::Navigator and one for the @ref
Acts::DirectNavigator that takes an additional `std::vector<const
Acts::Surface *>` argument).

The fit can be customized with several options. Important ones are:

- *maximum components*: how many components at most should be kept.
- *weight cut*: when to drop components.
- *component merging*: how a multi-component state is reduced to a single set
  of parameters and covariance, chosen via the @ref Acts::ComponentMergeMethod
  enum. Two methods are currently supported:
    - the *mean* computes the mean and the covariance of the mean.
    - *max weight* takes the parameters of the component with the maximum
      weight and computes the variance around them. This is a cheap
      approximation of the mode, which is not implemented currently.
- *mixture reduction*: how the number of components is reduced to the maximum
  allowed number, configured via an @ref Acts::Delegate:
    - *weight cut*: keep only the N components with the largest weights,
      implemented in @ref Acts::reduceMixtureLargestWeights.
    - *KL distance*: merge the closest components until the required amount
      is reached. The distance measure is the *Kullback-Leibler distance* in
      the *q/p* component, implemented in @ref
      Acts::reduceMixtureWithKLDistance.

> [!note]
> A good starting configuration is to use 12 components, the *max weight*
> merging and the *KL distance* reduction.

All options can be found in @ref Acts::GsfOptions.

If the GSF finds the column with the string identifier
`"gsf-final-multi-component-state"` (defined in
`Acts::GsfConstants::kFinalMultiComponentStateColumn`) in the track container,
it adds the final multi-component state to the track as a
`std::optional<Acts::MultiComponentBoundTrackParameters<SinglyCharged>>`
object.

A GSF example can be found in the ACTS Examples Framework
[here](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/truth_tracking_gsf.py).

### Customising the Bethe-Heitler approximation

The GSF needs an approximation of the Bethe-Heitler distribution as a
Gaussian mixture on each material interaction (see above). This task is
delegated to a separate class, which can be provided by a template parameter
to @ref Acts::GaussianSumFitter, so in principle it can be implemented in
different ways.

However, ACTS ships with @ref Acts::AtlasBetheHeitlerApprox, which implements
the ATLAS strategy for this task: to be able to evaluate the approximation of
the Bethe-Heitler distribution for different materials and thicknesses, the
individual Gaussian components (weight, mean, variance of the ratio
@f$E_f/E_i@f$) are parametrised as polynomials in @f$x/x_0@f$. This class can
load files in the ATLAS format found
[here](https://gitlab.cern.ch/atlas/athena/-/tree/main/Tracking/TrkFitter/TrkGaussianSumFilter/Data).
A default parameterization can be created with @ref
Acts::makeDefaultBetheHeitlerApprox.

@ref Acts::AtlasBetheHeitlerApprox is constructed with two parameterizations,
allowing different parameterizations to be used for different @f$x/x_0@f$. In
particular, it has this behaviour:

- @f$x/x_0 < 0.0001@f$: return no change.
- @f$x/x_0 < 0.002@f$: return a single Gaussian approximation.
- @f$x/x_0 < 0.1@f$: return the approximation for low @f$x/x_0@f$.
- @f$x/x_0 \geq 0.1@f$: return the approximation for high @f$x/x_0@f$. The
  maximum possible value is @f$x/x_0 = 0.2@f$; for higher values it is
  clipped to 0.2 and the GSF emits a warning.

### Further reading

- Thomas Atkinson, *Electron reconstruction with the ATLAS inner detector*,
  2006, see [here](https://cds.cern.ch/record/1448253).
- R. Frühwirth, *Track fitting with non-Gaussian noise*, 1997, see
  [here](https://doi.org/10.1016/S0010-4655(96)00155-5).
- R. Frühwirth, *A Gaussian-mixture approximation of the Bethe-Heitler model
  of electron energy loss by bremsstrahlung*, 2003, see
  [here](https://doi.org/10.1016/S0010-4655(03)00292-3).

## Global Chi-Square Fitter (GX2F)

In general, the GX2F is a weighted least-squares fit, minimising the

@f[
\chi^2 = \sum_i \frac{r_i^2}{\sigma_i^2}
@f]

of a track. Here, @f$r_i@f$ are the residuals, weighted with @f$\sigma_i^2@f$,
the covariance of the measurement (a detector property). Unlike the KF and
the GSF, the GX2F looks at all measurements at the same time and iteratively
minimises the starting parameters.

With the GX2F we obtain the final parameters @f$\vec\alpha_n@f$ from starting
parameters @f$\vec\alpha_0@f$. We set @f$\chi^2 = \chi^2(\vec\alpha)@f$ as a
function of the track parameters, though the @f$\chi^2@f$ minimisation could
be used for many other problems. Even in the context of track fitting, there
is a lot of freedom in how the GX2F is used: the residuals @f$r_i@f$ can have
many interpretations. Most of the time they are the distance between a
measurement and the prediction, but scattering angles or energy loss can also
be used as residuals. The subscript @f$i@f$ therefore mostly stands for a
measurement surface, since the sum runs over all of them.

This section covers:

- mathematical description of the base algorithm.
- mathematical description of multiple scattering.
- (coming soon) mathematical description of energy loss.
- implementation in ACTS.
- pros/cons.

### Mathematical description of the base algorithm

> [!note]
> The mathematical derivation is shortened at some places. A publication
> including the full derivation is planned.

To begin with, here is a short overview of the algorithm, with each step
described in more detail below:

1. Minimise the @f$\chi^2@f$ function.
2. Update the initial parameters (iteratively).
3. Calculate the covariance for the final parameters.

Before going into detail, a few symbols need to be introduced. As already
mentioned, there are track parameters @f$\vec\alpha@f$ to be fitted. To fit
them, the residuals are calculated as

@f[
r_i = m_i - f_i^m(\vec\alpha)
@f]

where @f$f^m(\vec\alpha)@f$ is the projection of the propagation function
@f$f(\vec\alpha)@f$ into the measurement dimension. Basically, if there is a
pixel measurement, the projection is onto the surface, discarding all angular
information. This projection can be different for each measurement surface.

#### 1. Minimise the @f$\chi^2@f$ function

The minimum of the @f$\chi^2@f$ function is expected at

@f[
\frac{\partial\chi^2(\vec\alpha)}{\partial\vec\alpha} = 0.
@f]

To find the zero(s) of this function, any method could be used, but ACTS
uses a modified [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method),
since it requires just another derivative of the @f$\chi^2@f$ function.

#### 2. Update the initial parameters (iteratively)

Since the Newton-Raphson method is used to find the minimum of the
@f$\chi^2@f$ function, iteration is needed. Each iteration (should) give
improved parameters @f$\vec\alpha@f$. While iterating, a system is updated in
this form:

@f[
\vec\alpha_{n+i} = \vec\alpha_n + \vec{\delta\alpha}_n.
@f]

After some derivations of the @f$\chi^2@f$ function and the Newton-Raphson
method, a matrix equation to calculate @f$\vec{\delta\alpha}_n@f$ is found:

@f[
[a_{kl}] \vec{\delta\alpha}_n = \vec b
@f]

with

@f[
a_{kl} = \sum_{i=1}^N \frac{1}{\sigma_i^2} \frac{\partial f_i^m(\vec\alpha)}{\partial \alpha_k}\frac{\partial f_i^m(\vec\alpha)}{\partial \alpha_l}
@f]

(where second order derivatives are omitted) and

@f[
b_k = \sum_{i=1}^N \frac{r_i}{\sigma_i^2} \frac{\partial f_i^m(\vec\alpha)}{\partial \alpha_k}.
@f]

At first sight, these expressions might seem intimidating and hard to
compute. But looking closer, these derivatives already exist in the
framework: all derivatives are elements of the Jacobian

@f[
\mathbf{J} = \begin{pmatrix}
                 \cdot & \dots & \cdot\\
                 \vdots & \frac{\partial f^m(\vec\alpha)}{\partial \alpha_k} & \vdots\\
                 \cdot & \dots & \cdot
             \end{pmatrix}.
@f]

At this point all information needed to perform a parameter update is
available; this is repeated until the parameters @f$\vec\alpha@f$ converge.

#### 3. Calculate the covariance for the final parameters

The calculation of the covariance of the final parameters is quite simple
compared to the steps before:

@f[
cov_{\vec\alpha} = [a_{kl}]^{-1}
@f]

Since it only depends on @f$[a_{kl}]@f$ of the last iteration, the GX2F does
not need an initial estimate for the covariance.

### Mathematical description of multiple scattering

To describe multiple scattering, the GX2F can fit the scattering angles as if
they were normal parameters. Of course, fitting more parameters increases the
dimensions of all matrices, making it computationally more expensive.

To recap multiple scattering: to describe scattering on a surface, only the
two angles @f$\theta@f$ and @f$\phi@f$ are needed, where:

- @f$\theta@f$ is the angle between the extrapolation of the incoming
  trajectory and the scattered trajectory.
- @f$\phi@f$ is the rotation around the extrapolation of the incoming
  trajectory.

This description is only valid for thin materials. To model thicker
materials, multiple thin materials could in theory be added together. It can
be shown that it is enough to use two sets of @f$\theta@f$ and @f$\phi@f$ on
both sides of the material, named @f$\theta_{in}@f$, @f$\theta_{out}@f$,
@f$\phi_{in}@f$, and @f$\phi_{out}@f$ -- but in the end they are just
multiple parameters in the fit, so only @f$\theta@f$ and @f$\phi@f$ (like for
thin materials) are considered here.

By defining residuals and covariances for the scattering angles, they can be
put into the @f$\chi^2@f$ function. For the residuals, since the expected
angle is 0:

@f[
r_s = \begin{cases}
         0 - \theta_s(\vec\alpha) \\
         0 - \sin(\theta_{loc})\phi_s(\vec\alpha)
      \end{cases}
@f]

with @f$\theta_{loc}@f$ the angle between the incoming trajectory and the
normal of the surface (there is no angle information @f$\phi@f$ if the
trajectory is perpendicular). For the covariances, the Highland form
(@cite ParticleDataGroup:2018ovx) is used:

@f[
\sigma_{scat} = \frac{13.6 \text{MeV}}{\beta c p} Z\prime \sqrt{\frac{x}{X_0}} \left( 1 + 0.038 \ln{\frac{x}{X_0}} \right)
@f]

with

- @f$x@f$: material layer thickness.
- @f$X_0@f$: radiation length.
- @f$p@f$: particle momentum.
- @f$Z\prime@f$: charge number.
- @f$\beta c@f$: velocity.

Combining these terms, the @f$\chi^2@f$ function including multiple
scattering can be written as

@f[
\chi^2 = \sum_{i=1}^N \frac{r_i^2}{\sigma_i^2} + \sum_{s}^S \left(\frac{\theta_s^2}{\sigma_s^2} + \frac{\sin^2{(\theta_{loc})}\phi_s^2}{\sigma_s^2}\right)
@f]

Note that both scattering angles have the same covariance.

### (coming soon) Mathematical description of energy loss

> [!note]
> The mathematical description of the energy loss for the GX2F is not yet
> written -- the corresponding development work has not finished yet.

### Implementation in ACTS

The implementation is in some points similar to the KF, since the KF
interface was chosen as a starting point. This makes it easier to replace
both fitters with each other. The structure of the GX2F implementation
follows coarsely the mathematical outline given above. It is best to start
reading the implementation from `fit()`:

1. Set up the fitter:
    - Actor
    - Aborter
    - Propagator
    - Variables needed for longer than one iteration
2. Iterate:
    1. Update parameters.
    2. Propagate through the geometry.
    3. Loop over the track and calculate and sum over:
        - @f$\chi^2@f$
        - @f$[a_{kl}]@f$
        - @f$\vec b@f$
    4. Solve @f$[a_{kl}] \vec{\delta\alpha}_n = \vec b@f$.
    5. Check for convergence.
3. Calculate the covariance of the final parameters.
4. Prepare and return the final track.

#### Configuration

Here is a simplified example of the configuration of the fitter.

```cpp
template <typename traj_t>
struct Gx2FitterOptions {
Gx2FitterOptions( ... ) : ... {}

Gx2FitterOptions() = delete;

...
//common options:
// geoContext, magFieldContext, calibrationContext, extensions,
// propagatorPlainOptions, referenceSurface, multipleScattering,
// energyLoss, freeToBoundCorrection

/// Max number of iterations during the fit (abort condition)
size_t nUpdateMax = 5;

/// Check for convergence (abort condition). Set to 0 to skip.
double relChi2changeCutOff = 1e-7;
};
```

Common options like the geometry context or toggling of the energy loss are
similar to the other fitters. For now there are three GX2F-specific options:

1. `nUpdateMax` sets an abort condition for the parameter update as a
   maximum number of allowed iterations. This condition is not meant to be
   hit in practice, but it stops the fit in case of poor convergence.
2. `relChi2changeCutOff` is the desired convergence criterion: at each step
   of the iteration, the current @f$\chi^2@f$ is compared to the previous
   one, and if the relative change is small enough, the fit finishes.

### Pros/Cons

There are some reasons for and against the GX2F. The biggest issue of the
GX2F is its performance: currently, the most expensive part is the
propagation, and since a full propagation is needed each iteration, at least
4-5 full propagations are needed. This is a lot compared to the 2
propagations of the KF. However, since the GX2F is a global fitter, it can
more easily resolve left-right-ambiguous measurements, like in the TRT
(Transition Radiation Tracker -- straw tubes).

Further reading:

- Kalman formalism reference: <https://twiki.cern.ch/twiki/pub/LHCb/ParametrizedKalman/paramKalmanV01.pdf>
- Highland form / multiple scattering derivation: <https://cds.cern.ch/record/1005181/files/thesis-2006-072.pdf#page=80>
