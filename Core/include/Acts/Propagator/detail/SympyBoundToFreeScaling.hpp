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

/// @brief Storage convention for the q/p column of the transported
///        bound-to-free jacobian in the sympy stepper
///
/// While a track is being stepped, the sympy stepper does not hold the plain
/// bound-to-free jacobian.  The free rows of its q/p column are multiplied by
/// q/p and divided by the column's own q/p row:
///
///     stored(i, eBoundQOverP) = qOverP * plain(i, eBoundQOverP)
///                                      / plain(eFreeQOverP, eBoundQOverP)
///
/// for the seven rows above eFreeQOverP, while the q/p row itself is stored
/// unchanged.  Both factors are storage convention only and carry no
/// approximation.  Together they are what lets the step kernel match the ATLAS
/// stepper's operation count on this column: the q/p factor turns the term
/// each Runge-Kutta stage picks up from the field's own q/p dependence into
/// the stage itself rather than a fresh cross product, and the normalisation
/// makes the coefficient of that term a literal one instead of a load of the
/// q/p row.
///
/// In vacuum neither factor changes across a step, so the convention is
/// preserved for free.  A dense step changes both and undoes and redoes the
/// scaling explicitly -- nothing here assumes the q/p row stays one.
///
/// Note that this is singular for q/p == 0.  So is the plain jacobian's own
/// q/p column, which the covariance transport divides by elsewhere, so the
/// convention adds no restriction that was not already there.
///
/// @param [in,out] jacobian bound-to-free jacobian to convert in place
/// @param [in] qOverP the q/p the column is scaled with
inline void toScaledBoundToFree(BoundToFreeMatrix& jacobian, double qOverP) {
  jacobian.block<eFreeQOverP, 1>(0, eBoundQOverP) *=
      qOverP / jacobian(eFreeQOverP, eBoundQOverP);
}

/// Undo @ref toScaledBoundToFree, recovering the plain bound-to-free jacobian.
///
/// @param [in,out] jacobian scaled bound-to-free jacobian to convert in place
/// @param [in] qOverP the q/p the column was scaled with
inline void fromScaledBoundToFree(BoundToFreeMatrix& jacobian, double qOverP) {
  jacobian.block<eFreeQOverP, 1>(0, eBoundQOverP) *=
      jacobian(eFreeQOverP, eBoundQOverP) / qOverP;
}

}  // namespace Acts::detail::sympy
