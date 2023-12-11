# Track Fitting

The track fitting algorithms estimate the track parameters.
It is part of the pattern recognition/track  reconstruction/tracking.
We can run the track fitting algorithms, after we allocated all hits to single tracks with the help of a track finding algorithm.
It is not necessary, that all points of a track are present.

Currently, we have implementations for three different fitters:
* Kalman Filter
* GSF
* Global Chi-Square Fitter (GX2F) [wip]
Even though all of them are least-squares fits, the concepts are quite different.
Therefore, we should not expect identical results from all of them.

(kf_core)=
## Kalman Filter (KF) [wip]
The Kalman Filter is an iterative fitter.
It successively combines measurements to obtain an estimate of the track parameters.
The KF needs an estimate as a starting point. The procedure alternates between two methods:
1. Extrapolate the current state to the next surface.
2. Update the extrapolation using the measurement of the new surface.[^billoir]
The meaning of "this surface" and "the next surface" changes with the context.
There are three different interpretations for this.
The KF can give us those three interpretations as sets of track parameters:
    * predicted: Uses "older" data (i.e. from the last surfaces) to make the prediction. This prediction is an extrapolation from the old data onto the current surface.
    * filtered: Uses the "current" data (i.e. the predicted data updated with the measurement on the current surface). It is some kind of weighted mean.
    * smoothed: Uses the "future" data to predict the current parameters. This can only be evaluated if the whole propagation is finished once. This can be done in to ways: one uses backwards-propagation and one does not.

:::{todo}
Complete Kalman Filter description
:::

(gsf_core)=
## Gaussian Sum Filter (GSF)

The GSF is an extension of the Kalman-Filter that allows to handle non-gaussian errors by modelling the track state as a gaussian mixture:

$$
p(\vec{x}) = \sum_i w_i \varphi(\vec{x}; \mu_i, \Sigma_i), \quad \sum_i w_i = 1
$$

A common use case of this is electron fitting. The energy-loss of Bremsstrahlung for electrons in matter are highly non-Gaussian, and thus cannot be modelled accurately by the default material interactions in the Kalman Filter. Instead, the Bremsstrahlung is modelled as a Bethe-Heitler distribution, where $z$ is the fraction of the energy remaining after the interaction ($E_f/E_i$), and $t$ is the material thickness in terms of the radiation length:

$$
f(z) = \frac{(- \ln z)^{c-1}}{\Gamma(c)}, \quad c = t/\ln 2
$$

(figBetheHeitler)=
:::{figure} figures/gsf_bethe_heitler_approx.svg
:width: 450px
:align: center
The true Bethe-Heitler distribution compared with a gaussian mixture approximation (in thin lines the individual components are drawn) at t = 0.1 (corresponds to ~ 10mm Silicon).
:::

To be able to handle this with the Kalman filter mechanics, this distribution is approximated by a gaussian mixture as well (see {numref}`figBetheHeitler`). The GSF Algorithm works then as follows (see also {numref}`figGsf`)

* On a surface with material, the Bethe-Heitler energy-loss distribution is approximated with a fixed number of gaussian components for each component. Since this way the number of components would grow exponentially with each material interaction, components that are close in terms of their *Kullback–Leibler divergence* are merged to limit the computational cost.
* On a measurement surface, for each component a Kalman update is performed. Afterwards, the component weights are corrected according to each component's compatibility with the measurement.

(figGsf)=
:::{figure} figures/gsf_overview.svg
:width: 450px
:align: center
Simplified overview of the GSF algorithm.
:::

### The Multi-Stepper
To implement the GSF, a special stepper is needed, that can handle a multi-component state internally: The {class}`Acts::MultiEigenStepperLoop`, which is based on the {class}`Acts::EigenStepper` and thus shares a lot of code with it. It interfaces to the navigation as one aggregate state to limit the navigation overhead, but internally processes a multi-component state. How this aggregation is performed can be configured via a template parameter, by default weighted average is used ({struct}`Acts::WeightedComponentReducerLoop`).

Even though the multi-stepper interface exposes only one aggregate state and thus is compatible with most standard tools, there is a special aborter is required to stop the navigation when the surface is reached, the {struct}`Acts::MultiStepperSurfaceReached`. It checks if all components have reached the target surface already and updates their state accordingly. Optionally, it also can stop the propagation when the aggregate state reaches the surface.


### Using the GSF

The GSF is implemented in the class {struct}`Acts::GaussianSumFitter`. The interface of its `fit(...)`-functions is very similar to the one of the {class}`Acts::KalmanFitter` (one for the standard {class}`Acts::Navigator` and one for the {class}`Acts::DirectNavigator` that takes an additional `std::vector<const Acts::Surface *>` as an argument):

```{doxygenstruct} Acts::GaussianSumFitter
---
members: fit
outline:
---
```

The fit can be customized with several options. Important ones are:
* *maximum components*: How many components at maximum should be kept.
* *weight cut*: When to drop components.
* *component merging*: How a multi-component state is reduced to a single set of parameters and covariance. The method can be chosen with the enum {enum}`Acts::ComponentMergeMethod`. Two methods are supported currently:
    * The *mean* computes the mean and the covariance of the mean.
    * *max weight* takes the parameters of component with the maximum weight and computes the variance around these. This is a cheap approximation of the mode, which is not implemented currently.
* *mixture reduction*: How the number of components is reduced to the maximum allowed number. Can be configured via a {class}`Acts::Delegate`:
    * *Weight cut*: Keep only the N components with the largest weights. Implemented in {func}`Acts::reduceMixtureLargestWeights`.
    * *KL distance*: Merge the closest components until the required amount is reached. The distance measure is the *Kullback-Leibler distance* in the *q/p* component. Implemented in {func}`Acts::reduceMixtureWithKLDistance`.

:::{note}
A good starting configuration is to use 12 components, the *max weight* merging and the *KL distance* reduction.
:::

All options can be found in the {struct}`Acts::GsfOptions`:

```{doxygenstruct} Acts::GsfOptions
---
outline:
---
```

If the GSF finds the column with the string identifier *"gsf-final-multi-component-state"* (defined in `Acts::GsfConstants::kFinalMultiComponentStateColumn`) in the track container, it adds the final multi-component state to the track as a `std::optional<Acts::MultiComponentBoundTrackParameters<SinglyCharged>>` object.

A GSF example can be found in the ACTS Examples Framework [here](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/truth_tracking_gsf.py).

### Customising the Bethe-Heitler approximation

The GSF needs an approximation of the Bethe-Heitler distribution as a Gaussian mixture on each material interaction (see above). This task is delegated to a separate class, that can be provided by a template parameter to {struct}`Acts::GaussianSumFitter`, so in principle it can be implemented in different ways.

However, ACTS ships with the class {class}`Acts::AtlasBetheHeitlerApprox` that implements the ATLAS strategy for this task: To be able to evaluate the approximation of the Bethe-Heitler distribution for different materials and thicknesses, the individual Gaussian components (weight, mean, variance of the ratio $E_f/E_i$) are parametrised as polynomials in $x/x_0$. This class can load files in the ATLAS format that can be found [here](https://gitlab.cern.ch/atlas/athena/-/tree/main/Tracking/TrkFitter/TrkGaussianSumFilter/Data). A default parameterization can be created with {func}`Acts::makeDefaultBetheHeitlerApprox`.

The {class}`Acts::AtlasBetheHeitlerApprox` is constructed with two parameterizations, allowing to use different parameterizations for different $x/x_0$. In particular, it has this behaviour:
* $x/x_0 < 0.0001$: Return no change
* $x/x_0 < 0.002$: Return a single gaussian approximation
* $x/x_0 < 0.1$: Return the approximation for low $x/x_0$.
* $x/x_0 \geq 0.1$: Return the approximation for high $x/x_0$. The maximum possible value is $x/x_0 = 0.2$, for higher values it is clipped to 0.2 and the GSF emits a warning.

### Further reading

* *Thomas Atkinson*, Electron reconstruction with the ATLAS inner detector, 2006, see [here](https://cds.cern.ch/record/1448253)
* *R Frühwirth*, Track fitting with non-Gaussian noise, 1997, see [here](https://doi.org/10.1016/S0010-4655(96)00155-5)
* *R Frühwirth*, A Gaussian-mixture approximation of the Bethe–Heitler model of electron energy loss by bremsstrahlung, 2003, see [here](https://doi.org/10.1016/S0010-4655(03)00292-3)

(gx2f_core)=
## Global Chi-Square Fitter (GX2F)

In general the GX2F is a weighted least squares fit, minimising the $\chi^2$
$$
\chi^2 = \sum_i \frac{r_i^2}{\sigma_i^2}
$$
of a track.
Here, $r_i$ are our residuals that we weight with $\sigma_i^2$, the covariance of the measurement (a detector property).
Unlike the KF and the GSF, the GX2F looks at all measurements at the same time and iteratively minimises the starting parameters.

With the GX2F we can obtain the final parameters $\vec\alpha_n$ from starting parameters $\vec\alpha_0$.
We set the $\chi^2 = \chi^2(\vec\alpha)$ as a function of the track parameters, but the $\chi^2$-minimisation could be used for many other problems.
Even in the context of track fitting, we are quite free on how to use the GX2F.
Especially the residuals $r_i$ can have many interpretations.
Most of the time we will see them as the distance between a measurement and our prediction.
But we can also use scattering angles, energy loss, ... as residuals.

This chapter on the GX2F guides through:
- Mathematical description of the base algorithm
- Mathematical description of the multiple scattering [wip]
- (coming soon) Mathematical description of the energy loss
- Implementation in ACTS[wip]
- Pros/Cons [wip]

### Mathematical description of the algorithm

To begin with, there will be a short overview on the algorithm.
Later in this section, each step is described in more detail.
1. Minimise the $\chi^2$ function
2. Update the initial parameters (iteratively)
3. Calculate the covariance for the final parameters

#### 1. Minimise the $\chi^2$ function

We expect the minimum of the $\chi^2$ function at
$$
\frac{\partial\chi^2(\vec\alpha)}{\partial\vec\alpha} = 0.
$$
To find the zero(s) of this function we could use any method but we will stick to a modified [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method),
since it requires just another derivative of the $\chi^2$ function.

TODO:


#### 2. Update the initial parameters (iteratively)

Since we are using the Newton-Raphson method to find the minimum of the $\chi^2$ function, we need to iterate.


#### 3. Calculate the covariance for the final parameters

The calculation of the covariance of the final parameters is quite simple compared to the steps before:
$$
cov_{\vec\alpha} = [a_{kl}]^{-1}
$$
Since it only depends on the $[a_{kl}]$ of the last iteration, the GX2F does not need an initial estimate for the covariance.


### Mathematical description of the multiple scattering [wip]

:::{todo}
Write GX2F: Mathematical description of the multiple scattering
:::

### (coming soon) Mathematical description of the energy loss [wip]

:::{todo}
Write GX2F: Mathematical description of the energy loss

The development work on the energy loss has not finished yet.
:::


### Implementation in ACTS [wip]

:::{todo}
Write GX2F: Implementation in ACTS
:::

### Pros/Cons [wip]
:::{todo}
Write GX2F: Pros/Cons
:::

[^billoir]: https://twiki.cern.ch/twiki/pub/LHCb/ParametrizedKalman/paramKalmanV01.pdf
