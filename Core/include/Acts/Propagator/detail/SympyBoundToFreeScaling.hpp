// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

namespace Acts::detail::sympy {

// Storage convention for the q/p column of the bound-to-free jacobian.
//
// The plain column differentiates the free parameters by the *starting* q/p.
// While stepping, the sympy stepper differentiates them by the log of the
// *current* q/p:
//
//     stored(i, eBoundQOverP) = d(free_i) / d(log|q/p|)
//                             = qOverP * plain(i, eBoundQOverP)
//                                      / plain(eFreeQOverP, eBoundQOverP)
//
// The division is the chain rule to the current q/p -- no other bound parameter
// reaches q/p, so that one row is all of it -- and the factor carries it on to
// the log. The row is stored plain, so the conversion is exact both ways.
//
// It also makes the column cheap: a bend vector is linear in q/p, so it is its
// own derivative by log q/p and each Runge-Kutta stage's field term is the
// stage itself. Singular at q/p == 0, as log q/p is. See @ref sympy_codegen.

/// Convert a plain bound-to-free jacobian to the stored form.
///
/// @param [in,out] jacobian jacobian to convert in place
/// @param [in] qOverP the current q/p
inline void toScaledBoundToFree(BoundToFreeMatrix& jacobian, double qOverP) {
  jacobian.block<eFreeQOverP, 1>(0, eBoundQOverP) *=
      qOverP / jacobian(eFreeQOverP, eBoundQOverP);
}

/// Recover the plain bound-to-free jacobian from the stored form.
///
/// @param [in,out] jacobian jacobian to convert in place
/// @param [in] qOverP the q/p the column is stored against
inline void fromScaledBoundToFree(BoundToFreeMatrix& jacobian, double qOverP) {
  jacobian.block<eFreeQOverP, 1>(0, eBoundQOverP) *=
      jacobian(eFreeQOverP, eBoundQOverP) / qOverP;
}

/// Move a stored column onto a new q/p, leaving the plain jacobian it stands
/// for unchanged. Needed whenever q/p changes outside the step kernel.
///
/// @param [in,out] jacobian jacobian to rescale in place
/// @param [in] oldQOverP the q/p the column is stored against
/// @param [in] newQOverP the q/p to store it against instead
inline void rescaleBoundToFree(BoundToFreeMatrix& jacobian, double oldQOverP,
                               double newQOverP) {
  jacobian.block<eFreeQOverP, 1>(0, eBoundQOverP) *= newQOverP / oldQOverP;
}

}  // namespace Acts::detail::sympy
