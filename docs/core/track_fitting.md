# Track Fitting

## Gaussian Sum Filter

The GSF is an extension of the Kalman-Filter that allows to handle non-gaussian errors by modelling the track state as a gaussian mixture:

$$
p(\vec{x}) = \sum_i w_i \varphi(\vec{x}; \mu_i, \Sigma_i), \quad \sum_i w_i = 1
$$

In particular this is used in Acts for electron refitting. The energy-loss of Brehmsstrahlung for electrons in matter are highly non-gaussian, and thus are not modeled accurately by the default material interactions in the Kalman Filter. Instead, the Brehmsstrahlung is modeled as a Bethe-Heitler distribution, which is approximated as a gaussian mixture. For further reading see [here](https://cds.cern.ch/record/1448253).

### Implementation

To implement the GSF, a special stepper is needed, that can handle a multi-component state internally: The {class}`Acts::MultiEigenStepperLoop`. On a surface with material, for each component the material loss is modeled with a approximated Bethe-Heitler distribution with a fixed number of components. Since the number of components would grow exponentially with each material interaction, close components are merged afterwards to limit the computational cost.

The kalman update mechanism is based on the code for the {class}`Acts::KalmanFitter`. For the smoothing, the GSF implements the smoothing algorithm described [here](https://doi.org/10.1016/S0010-4655(96)00155-5).

### Using the GSF

The GSF is implemented in the class {class}`Acts::GaussianSumFitter`. The interface of its `fit(...)`-functions is very similar to the one of the {class}`Acts::KalmanFitter`. The fit can be customized with several options, e.g., the maximum number of components. All options can be found in the {struct}`GsfOptions`.

To simplify integration, the GSF also returns a {class}`Acts::KalmanFitter` object. This allows to use the same analysis tools for both fitters. Currently, the states of the individual components are not returned by the fitter.

For an usage example in the Acts Examples Framework see [here](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/truth_tracking_gsf.py).

### Customizing the Bethe-Heitler approximation

The approximation of the Bethe-Heitler distribution is described in {class}`Acts::BetheHeitlerApprox`. The classis  templated on the number of components and the degree of the polynomial, and is designed to be used with the [parameterization files from ATLAS](https://gitlab.cern.ch/atlas/athena/-/tree/master/Tracking/TrkFitter/TrkGaussianSumFilter/Data). However, in principle the GSF could be constructed with custom classes with the same interface of {class}`Acts::BetheHeitlerApprox`.

To be able to evaluate the approximation of the Bethe-Heitler distribution for different materials and thicknesses, the individuall gaussian components (weight, mean, variance of the ratio $E_f/E_i$) are parameterized as polynomials in $x/x_0$. The default parameterization uses 6 components and 5th order polynomials.

For small $x/x_0$ the {class}`Acts::BetheHeitlerApprox` only returns a one-component mixture or even no change at all. When loading custom parameterizations, it is possible to specify one for high momentum and one for low $x/x_0$. The thresholds cannot be customized currently.
