@defgroup propagation Propagation
@brief Particle propagation

## Steppers

### GenericRungeKuttaStepper

@ref Acts::GenericRungeKuttaStepper is a reference stepper, similar to what
@ref Acts::TryAllNavigator is for navigation. It is not optimised for speed.
It integrates the equation of motion with any explicit Runge-Kutta tableau and
an exact transport jacobian, so other steppers can be validated against it.
