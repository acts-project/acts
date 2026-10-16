@defgroup propagation Propagation
@brief Particle propagation

The transport Jacobian used for covariance transport is not written by hand: it
is derived symbolically with sympy and printed as C++ at build time. See
@ref codegen_sympy_jac for the derivation, the sparsity patterns it exploits
and the code it produces.
