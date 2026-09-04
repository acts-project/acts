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

// Scaling convention for the q/p column of the bound-to-free jacobian.
//
// While stepping, the sympy stepper differentiates the free parameters by the
// log of the *current* q/p rather than by the starting q/p:
//
//     scaled(i, eBoundQOverP) = d(free_i) / d(log|q/p|)
//                             = qOverP * plain(i, eBoundQOverP)
//                                      / plain(eFreeQOverP, eBoundQOverP)
//
// The q/p row itself stays plain, which is what makes the conversion exact both
// ways. Why the stepper wants this form: @ref sympy_codegen.
//
// A freshly initialized surface jacobian has e_qop as its q/p column, so
// scaling it is a no-op; code that writes one into the state directly
// (LoopComponentProxy::update) needs no conversion. Singular at q/p == 0, as
// the plain column already is.

/// Convert a plain bound-to-free jacobian to the scaled form.
///
/// @param [in,out] jacobian jacobian to convert in place
/// @param [in] qOverP the current q/p
inline void toScaledBoundToFree(BoundToFreeMatrix& jacobian, double qOverP) {
  jacobian.block<eFreeQOverP, 1>(0, eBoundQOverP) *=
      qOverP / jacobian(eFreeQOverP, eBoundQOverP);
}

/// Recover the plain bound-to-free jacobian from the scaled form.
///
/// @param [in,out] jacobian jacobian to convert in place
/// @param [in] qOverP the q/p the column is scaled against
inline void fromScaledBoundToFree(BoundToFreeMatrix& jacobian, double qOverP) {
  jacobian.block<eFreeQOverP, 1>(0, eBoundQOverP) *=
      jacobian(eFreeQOverP, eBoundQOverP) / qOverP;
}

/// Move a scaled column onto a new q/p, leaving the plain jacobian it stands
/// for unchanged. Needed whenever q/p changes outside the step kernel.
///
/// @param [in,out] jacobian jacobian to rescale in place
/// @param [in] oldQOverP the q/p the column is scaled against
/// @param [in] newQOverP the q/p to scale it against instead
inline void rescaleBoundToFree(BoundToFreeMatrix& jacobian, double oldQOverP,
                               double newQOverP) {
  jacobian.block<eFreeQOverP, 1>(0, eBoundQOverP) *= newQOverP / oldQOverP;
}

}  // namespace Acts::detail::sympy
