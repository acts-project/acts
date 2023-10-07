# Track Fitting

The track fitting algorithms estimate the track parameters.
It is part of the pattern recognition/track  reconstruction/tracking.
We can run the track fitting algorithms, after we allocated all hits to single tracks with the help of a track finding algorithm.
It is not necessary, that all points of a track are present.

Currently we have implementations for three different fitters:
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

:::{note}
This chapter will be extended in the future.
:::

(gsf_core)=
## Gaussian Sum Filter (GSF)

:::{note}
The GSF is not considered as production ready yet, therefore it is located in the namespace `Acts::Experimental`.
:::

The GSF is an extension of the Kalman-Filter that allows to handle non-gaussian errors by modelling the track state as a gaussian mixture:

$$
p(\vec{x}) = \sum_i w_i \varphi(\vec{x}; \mu_i, \Sigma_i), \quad \sum_i w_i = 1
$$

A common use case of this is electron fitting. The energy-loss of Bremsstrahlung for electrons in matter are highly non-Gaussian, and thus cannot be modelled accurately by the default material interactions in the Kalman Filter. Instead, the Bremsstrahlung is modelled as a Bethe-Heitler distribution, where $z$ is the fraction of the energy remaining after the interaction ($E_f/E_i$), and $t$ is the material thickness in terms of the radiation length:

$$
f(z) = \frac{(- \ln z)^{c-1}}{\Gamma(c)}, \quad c = t/\ln 2
$$

(figBetheHeitler)=
:::{figure} ../figures/gsf_bethe_heitler_approx.svg
:width: 450px
:align: center
The true Bethe-Heitler distribution compared with a gaussian mixture approximation (in thin lines the individual components are drawn) at t = 0.1 (corresponds to ~ 10mm Silicon).
:::

To be able to handle this with the Kalman filter mechanics, this distribution is approximated by a gaussian mixture as well (see {numref}`figBetheHeitler`). The GSF Algorithm works then as follows (see also {numref}`figGsf`)

* On a surface with material, the Bethe-Heitler energy-loss distribution is approximated with a fixed number of gaussian components for each component. Since this way the number of components would grow exponentially with each material interaction, components that are close in terms of their *Kullback–Leibler divergence* are merged to limit the computational cost.
* On a measurement surface, for each component a Kalman update is performed. Afterwards, the component weights are corrected according to each components compatibility with the measurement.

(figGsf)=
:::{figure} ../figures/gsf_overview.svg
:width: 450px
:align: center
Simplified overview of the GSF algorithm.
:::

### The Multi-Stepper
To implement the GSF, a special stepper is needed, that can handle a multi-component state internally: The {class}`Acts::MultiEigenStepperLoop`, which is based on the {class}`Acts::EigenStepper` and thus shares a lot of code with it. It interfaces to the navigation as one aggregate state to limit the navigation overhead, but internally processes a multi-component state. How this aggregation is performed can be configured via a template parameter, by default weighted average is used ({struct}`WeightedComponentReducerLoop`).

At the end of the fit the multi-component state must be reduced to a single set of parameters with a corresponding covariance matrix. This is supported by the {class}`Acts::MultiEigenStepperLoop` in two different ways currently: The *mean* method computes the mean and the covariance matrix of the multi-component state, whereas the *maximum weight* method just returns the component with the maximum weight. This can be configured in the constructor of the {class}`Acts::MultiEigenStepperLoop` with the {enum}`Acts::FinalReductionMethod`. In the future there is planned to add a *mode* finding method as well.

:::{note}
In practice it turned out that the *maximum weight* method leads to better results so far.
:::

Even though the multi-stepper interface exposes only one aggregate state and thus is compatible with most standard tools, there is a special aborter is required to stop the navigation when the surface is reached, the {class}`Acts::MultiStepperSurfaceReached`. It checks if all components have reached the target surface already and updates their state accordingly. Optionally, it also can stop the propagation when the aggregate state reaches the surface.


### Using the GSF

The GSF is implemented in the class {class}`Acts::Experimental::GaussianSumFitter`. The interface of its `fit(...)`-functions is very similar to the one of the {class}`Acts::KalmanFitter` (one for the standard {class}`Acts::Navigator` and one for the {class}`Acts::DirectNavigator` that takes an additional `std::vector<const Acts::Surface *>` as an argument):

```{doxygenstruct} Acts::Experimental::GaussianSumFitter
---
members: fit
outline:
---
```

The fit can be customized with several options, e.g., the maximum number of components. All options can be found in the {struct}`Acts::GsfOptions`.

To simplify integration, the GSF returns a {class}`Acts::KalmanFitterResult` object, the same as the {class}`Acts::KalmanFitter`. This allows to use the same analysis tools for both fitters. Currently, the states of the individual components are not returned by the fitter.

A GSF example can be found in the Acts Examples Framework [here](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/truth_tracking_gsf.py).

### Customising the Bethe-Heitler approximation

The GSF needs an approximation of the Bethe-Heitler distribution as a Gaussian mixture on each material interaction (see above). This task is delegated to a separate class, that can be provided by a template parameter to {class}`Acts::Experimental::GaussianSumFitter`, so in principle it can be implemented in different ways.

However, ACTS ships with the class {class}`Acts::Experimental::AtlasBetheHeitlerApprox` that implements the ATLAS strategy for this task: To be able to evaluate the approximation of the Bethe-Heitler distribution for different materials and thicknesses, the individual Gaussian components (weight, mean, variance of the ratio $E_f/E_i$) are parametrised as polynomials in $x/x_0$. This class can load files in the ATLAS format that can be found [here](https://gitlab.cern.ch/atlas/athena/-/tree/master/Tracking/TrkFitter/TrkGaussianSumFilter/Data). A default parameterization can be created with {func}`Acts::Experimental::makeDefaultBetheHeitlerApprox`.

The {class}`Acts::Experimental::AtlasBetheHeitlerApprox` is constructed with two parameterizations, allowing to use different parameterizations for different $x/x_0$. In particular, it has this behaviour:
* $x/x_0 < 0.0001$: Return no change
* $x/x_0 < 0.002$: Return a single gaussian approximation
* $x/x_0 < 0.1$: Return the approximation for low $x/x_0$.
* $x/x_0 \geq 0.1$: Return the approximation for high $x/x_0$. The maximum possible value is $x/x_0 = 0.2$, for higher values it is clipped to 0.2 and the GSF emits a warning.

### Further reading
* *Thomas Atkinson*, Electron reconstruction with the ATLAS inner detector, 2006, see [here](https://cds.cern.ch/record/1448253)
* *R Frühwirth*, Track fitting with non-Gaussian noise, 1997, see [here](https://doi.org/10.1016/S0010-4655(96)00155-5)
* *R Frühwirth*, A Gaussian-mixture approximation of the Bethe–Heitler model of electron energy loss by bremsstrahlung, 2003, see [here](https://doi.org/10.1016/S0010-4655(03)00292-3)

(gx2f_core)=
## Global Chi-Square Fitter (GX2F) [wip]
:::{note}
This chapter will be added soon.
:::

[^billoir]: https://twiki.cern.ch/twiki/pub/LHCb/ParametrizedKalman/paramKalmanV01.pdf
