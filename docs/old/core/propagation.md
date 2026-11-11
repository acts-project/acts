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

Additional to these mandatory components, the propagator runs exactly one **actor** to allow for custom behaviour. The actor is a function object that is hooked in the propagation loop. It performs actions on the propagator state (e.g. the {class}`Acts::KalmanFitter` uses an actor), and it decides when the propagation stops.

The propagator exposes its state to the actor as arguments to `act()` and `checkAbort()`. An actor can define a default-constructable `state_type`. The propagator holds one instance of it in `state.actorState` and hands it to each call:

```cpp
template <typename propagator_state_t, typename stepper_t, typename navigator_t>
Result<void> act(propagator_state_t &state, const stepper_t &stepper,
                 const navigator_t &navigator, state_type &actorState,
                 const Logger &logger) const {
  const auto &navigatorState = state.navigation;
  const auto &stepperState = state.stepping;
  const auto &options = state.options;
  return Result<void>::success();
}
```

The actor is fully responsible for the end of the propagation. The free functions in `Acts/Propagator/detail/AbortConditions.hpp` check the usual conditions, for example the path limit and the target surface.

The result of a propagation holds the track parameters at the endpoint of the propagation, the transport Jacobian, the number of steps, the path length and the statistics. The actor state stays in the propagator state.

## Initialization and running

The {class}`Acts::Propagator` is initialized with the helper class
{struct}`Acts::PropagatorOptions`, which is templated on the actor. If no actor is given, the propagator uses {struct}`Acts::DefaultActor`. It stops at the target surface, at the path limit (lowered by the loop protection), outside the volumes that the options allow, and, if `stopAtEndOfWorld` is set, at the end of the world.

To run several actors, combine them with the composite actor {struct}`Acts::ActorList`. Add {struct}`Acts::DefaultActor` to keep its abort conditions. Put it first, so the loop protection sees the start momentum before other actors change it:

```cpp
using MyOptions = Propagator::Options<
                    Acts::ActorList<Acts::DefaultActor, MyActor1, MyActor2>
                  >;
```

The actors of the list can be accessed with the `get`-method that expects the corresponding actor type as template parameter. Besides this, the {struct}`Acts::PropagatorOptions` also contains a lot of general options like the `maxStepSize`:

```cpp
auto options = MyOptions(gctx, mctx);
options.actor.get<MyActor1>().foo = bar;
options.stepping.maxStepSize = 100;
```

All available options can be found in the {struct}`Acts::PropagatorPlainOptions`, from which {struct}`Acts::PropagatorOptions` inherits.

:::{tip}
The propagator also contains a loop-protection mechanism. It estimates a circle perimeter from the momentum and the magnetic field, and aborts the propagation when a certain fraction (default: 0.5) of the circle has been propagated. This behaviour can be changed in the {struct}`Acts::PropagatorOptions` via the boolean `loopProtection` and the float `loopFraction`.
:::

To run the propagation, we must call the member function `propagate(...)` with the initial track parameters and the propagator options. With a target surface, the navigator steers to the surface, and {struct}`Acts::DefaultActor` stops there.

The result is an instance of {class}`Acts::Result`. It contains the actual result, or an error code in case something went wrong. To read the actor state after the propagation, make the propagator state yourself and propagate with it:

```cpp
auto state = propagator.makeState(options);
auto res = propagator.propagate(state, myParams);

if (res.ok()) {
  state.actorState.get<MyActor1::state_type>();
}
```

## Navigators

ACTS comes with a couple of different navigator implementations:
- The standard navigator {class}`Acts::Navigator` which performs the full navigation in the volume/layer/surface hierarchy
- The {class}`Acts::DirectNavigator` which takes a sequence of surfaces and just navigates to one after the other.
- The {class}`Acts::TryAllNavigator` which, as the name suggests, tries to intersect all available surfaces without acceleration structure and special assumptions. This navigator is meant for validation rather than production.
- The {class}`Acts::TryAllOverstepNavigator` which is similar to the {class}`Acts::TryAllNavigator`, but deliberately oversteps and then intersects surfaces which might have been missed by that step and targets them. This navigator is meant for validation rather than production.

The navigators provide information about the current position inside the geometry in their state variable (e.g. {struct}`Acts::Navigator::State` and {struct}`Acts::DirectNavigator::State`), e.g. pointers to the `currentSurface` and the `currentVolume`.

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

The {class}`Acts::EigenStepper` implements the same functionality as the ATLAS stepper, however, the stepping code is rewritten by using `Eigen` primitives. Thus, it also uses a 4th-order Runge-Kutta algorithm for the integration of the EOM. Additionally, the {class}`Acts::EigenStepper` allows to customize the concrete integration step via **extension**.

The extension encapsulate the relevant equations for different environments. There exists a {struct}`Acts::EigenStepperDefaultExtension` that is suited for propagation in a vacuum, and the {struct}`Acts::EigenStepperDenseExtension`, that contains additional code to handle the propagation inside materials. Which extension is used is decided by the user.

The extension can be configured via the {class}`Acts::EigenStepper`:

```c++
using Stepper = Acts::EigenStepper<Acts::EigenStepperDenseExtension>;
```

By default, the {class}`Acts::EigenStepper` only uses the {struct}`Acts::EigenStepperDenseExtension`.

### MultiStepperLoop

The {class}`Acts::MultiStepperLoop` can extend any single-component stepper and is designed to internally handle a multi-component state, while interfacing as a single component to the navigator. It is mainly used for the {struct}`Acts::GaussianSumFitter`.
