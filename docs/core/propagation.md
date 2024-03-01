(propagation_impl)=
# Propagation and extrapolation

The track propagation is an essential part of track reconstruction. This section describes the high-level classes and concepts used for this task in ACTS.

## Overview: Steppers, Navigators and Actors

The propagation through a geometry is based on the interaction of two different components:

* The **Stepper** provides the implementation of the solution of the equation of motion (either by analytical means or through numerical integration).
* The **Navigator** keeps track of the current position in the geometry and adjusts the step size so that the stepper does not step through a surface.

Following the general ACTS design, these classes do not manage their internal state via member variables, but provide an internal `State` struct which contains all relevant data and is managed by the propagator.

The interaction of these two components is handled by the {class}`Acts::Propagator` class template that takes the stepper and the navigator as template parameters:

```cpp
Propagator<Navigator, Stepper>
```

Additional to these mandatory components, the propagator can be equipped with **Actors** and **Aborters** to allow for custom behaviour. These are function objects that are hooked in the propagation loop. Actors just perform some action on the propagator state (e.g. the {class}`Acts::KalmanFitter` is an actor), aborts can abort propagation (e.g., the {struct}`Acts::PathLimitReached`).

The propagator exposes its state to the actors and aborters as arguments to `operator()`. Actors must define a default-constructable `result_type`, which can be modified in each call:

```cpp
template<typename propagator_state_t, typename stepper_t>
auto operator(propagator_state_t &state, const stepper_t &stepper, result_type &result) const {
  const auto &navigatorState = state.navigation;
  const auto &stepperState = state.stepping;
  const auto &options = state.options;
}
```

The result of a propagation consists of the track parameters at the endpoint of the propagation as well as the results of all actors.

## Initialization and running

The {class}`Acts::Propagator` is initialized with the helper class
{struct}`Acts::PropagatorOptions`, which is templated on the list of actors and
aborters. This is done with the classes {struct}`Acts::ActionList` and
{struct}`Acts::AbortList` (which are in fact small wrappers around
`std::tuple`):

```cpp
using MyOptions = Acts::PropagatorOptions<
                    Acts::ActionList<MyActor1, MyActor2>,
                    Acts::AbortList<MyAborter1, MyAborter2>
                  >;
```

The actors and aborters are instantiated with the options and can be accessed with the `get`-method that expects the corresponding actor type as template parameter. Besides this, the {struct}`Acts::PropagatorOptions` also contains a lot of general options like the `maxStepSize`:

```cpp
auto options = MyOptions();
options.actionList.get<MyActor1>().foo = bar;
options.maxStepSize = 100;
```

All available options can be found in the {struct}`Acts::PropagatorPlainOptions`, from which {struct}`Acts::PropagatorOptions` inherits.

:::{tip}
The propagator also contains a loop-protection mechanism. It estimates a circle perimeter from the momentum and the magnetic field, and aborts the propagation when a certain fraction (default: 0.5) of the circle has been propagated. This behaviour can be changed in the {struct}`Acts::PropagatorOptions` via the boolean `loopProtection` and the float `loopFraction`.
:::

To run the propagation, we must call the member function `propagate(...)` with the initial track parameters and the propagator options. There are several overloads to the `propagate(...)` function, which allow further customization:
* With/without a target surface: The overload with a target surface automatically adds an aborter for the passed `Surface` to the `AbortList`.
* With/without a prepared result object: Without a result object, a suitable result object is default-constructed internally.

The result is an instance of {class}`Acts::Result`. It contains the actual result, or an error code in case something went wrong. In the actual result, the results of the different actors can again be accessed via a `get` method:

```cpp
auto res = propagator.propagate(myParams, options);

if( res.ok() ) {
  res.value().get<MyActor1::result_type>();
}
```

## Navigators

ACTS comes with two navigators: The standard navigator {class}`Acts::Navigator` that performs the full navigation in the volume/layer/surface hierarchy, and the {class}`Acts::DirectNavigator` that takes a sequence of surfaces and just navigates to one after the other. This sequence must be initialized with a special actor, the {struct}`Acts::DirectNavigator::Initializer`.

The navigators provide information about the current position inside the geometry in their state variable ({struct}`Acts::Navigator::State` and {struct}`Acts::DirectNavigator::State`), e.g. pointers to the `currentSurface` and the `currentVolume`.

:::{tip}
The {class}`Acts::Navigator` by default does a straight-line extrapolation to resolve layer candidates. In certain geometries (e.g., telescope) this can lead to the effect that bent tracks miss some layers. This can be mitigated by disabling the bound-check for layer resolution with the `boundaryCheckLayerResolving` option in {struct}`Acts::Navigator::Config`.
:::

## Steppers

ACTS also provides a variety of stepper implementations. Since these in general can work very differently internally, the state itself is not the main interface to the steppers. Instead, all steppers provide a common API, to that we can pass instances of the stepper state. This allows a generic and template-based design even for very different steppers:

```cpp
template<typename propagator_state_t, typename stepper_t>
auto operator(propagator_state_t &state, const stepper_t &stepper) const {
  stepper.foo(state.stepping);
}
```

### AtlasStepper

The {class}`Acts::AtlasStepper` is a pure transcript from the ATLAS `RungeKuttaPropagator` and `RungeKuttaUtils`.

### StraightLineStepper

The {class}`Acts::StraightLineStepper` is a very stripped down stepper that just performs a linear propagation without magnetic field. It can be used for debugging, validation or other simple tasks.

### EigenStepper

The {class}`Acts::EigenStepper` implements the same functionality as the ATLAS stepper, however, the stepping code is rewritten by using `Eigen` primitives. Thus, it also uses a 4th-order Runge-Kutta algorithm for the integration of the EOM. Additionally, the {class}`Acts::EigenStepper` allows to customize the concrete integration step via **extensions**.

The extensions encapsulate the relevant equations for different environments. There exists a {struct}`Acts::DefaultExtension` that is suited for propagation in a vacuum, and the {struct}`Acts::DenseEnvironmentExtension`, that contains additional code to handle the propagation inside materials. Which extension is used is selected by a bidding-system.

The extension can be configured via the {struct}`Acts::StepperExtensionList`:

```cpp
using Stepper = Acts::EigenStepper<
                  Acts::StepperExtensionList<
                    Acts::DefaultExtension,
                    Acts::DenseEnvironmentExtension
                  >
                >;
```

By default, the {class}`Acts::EigenStepper` only uses the {struct}`Acts::DefaultExtension`.

### MultiEigenStepperLoop

The {class}`Acts::MultiEigenStepperLoop` is an extension of the {class}`Acts::EigenStepper` and is designed to internally handle a multi-component state, while interfacing as a single component to the navigator. It is mainly used for the {struct}`Acts::GaussianSumFitter`.
