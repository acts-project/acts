@defgroup propagation Propagation
@brief Particle propagation

The transport Jacobian used for covariance transport is not written by hand: it
is derived symbolically with sympy and printed as C++ at build time. See
@ref codegen_sympy_jac for the derivation, the sparsity patterns it exploits
and the code it produces.

## Steppers

### GenericRungeKuttaStepper

@ref Acts::GenericRungeKuttaStepper is a reference stepper, similar to what
@ref Acts::Experimental::TryAllNavigator is for navigation. It is not optimised
for speed. It integrates the equation of motion with any explicit Runge-Kutta
tableau and an exact transport jacobian, so other steppers can be validated
against it.
