# Track Fitting

(gsf_core)=
## Gaussian Sum Filter

The GSF is an extension of the Kalman-Filter that allows to handle non-gaussian errors by modelling the track state as a gaussian mixture:

$$
p(\vec{x}) = \sum_i w_i \varphi(\vec{x}; \mu_i, \Sigma_i), \quad \sum_i w_i = 1
$$

A common use case of this is electron fitting. The energy-loss of Bremsstrahlung for electrons in matter are highly non-gaussian, and thus are not modeled accurately by the default material interactions in the Kalman Filter. Instead, the Bremsstrahlung is modeled as a Bethe-Heitler distribution, which is approximated as a gaussian mixture.

### Implementation

To implement the GSF, a special stepper is needed, that can handle a multi-component state internally: The {class}`Acts::MultiEigenStepperLoop`. On a surface with material, the Bethe-Heitler energy-loss distribution is approximated with a fixed number of gaussian distributions for each component. Since the number of components would grow exponentially with each material interaction, components that are close in terms of their *Kullback–Leibler divergence* are merged to limit the computational cost. The Kalman update mechanism is based on the code for the {class}`Acts::KalmanFitter`.

### Using the GSF

The GSF is implemented in the class {class}`Acts::GaussianSumFitter`. The interface of its `fit(...)`-functions is very similar to the one of the {class}`Acts::KalmanFitter` (one for the standard {class}`Acts::Navigator` and one for the {class}`Acts::DirectNavigator` that takes an additional `std::vector<const Acts::Surface *>` as an argument):

```{doxygenstruct} Acts::GaussianSumFitter
---
members: fit
outline:
---
```

The fit can be customized with several options, e.g., the maximum number of components. All options can be found in the {struct}`Acts::GsfOptions`.

To simplify integration, the GSF returns a {class}`Acts::KalmanFitterResult` object, the same as the {class}`Acts::KalmanFitter`. This allows to use the same analysis tools for both fitters. Currently, the states of the individual components are not returned by the fitter.

A GSF example can be found in the Acts Examples Framework [here](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/truth_tracking_gsf.py).

### Customizing the Bethe-Heitler approximation

To be able to evaluate the approximation of the Bethe-Heitler distribution for different materials and thicknesses, the individual gaussian components (weight, mean, variance of the ratio $E_f/E_i$) are parameterized as polynomials in $x/x_0$. The default parametrization uses 6 components and 5th order polynomials.

This approximation of the Bethe-Heitler distribution is described in {class}`Acts::BetheHeitlerApprox`. The class is templated on the number of components and the degree of the polynomial, and is designed to be used with the [parameterization files from ATLAS](https://gitlab.cern.ch/atlas/athena/-/tree/master/Tracking/TrkFitter/TrkGaussianSumFilter/Data). However, in principle the GSF could be constructed with custom classes with the same interface as {class}`Acts::BetheHeitlerApprox`.

For small $x/x_0$ the {class}`Acts::BetheHeitlerApprox` only returns a one-component mixture or no change at all. When loading a custom parametrization, it is possible to specify different parameterizations for high and for low $x/x_0$. The thresholds are currently not configurable.

### Further reading
* *Thomas Atkinson*, Electron reconstruction with the ATLAS inner detector, 2006, see [here](https://cds.cern.ch/record/1448253)
* *R Frühwirth*, Track fitting with non-Gaussian noise, 1997, see [here](https://doi.org/10.1016/S0010-4655(96)00155-5)

(kf_core)=
## Kalman Filter

:::{note}
This is a stub!
:::
