// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <tuple>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/StepperState.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class Surface;

namespace detail {

/// @brief These functions perform the transport of a covariance matrix using
/// given Jacobians. The required data is provided by a @p StepperState object
/// with some additional data. Since this is a purely algebraic problem the
/// calculations are identical for @c StraightLineStepper and @c EigenStepper.
/// As a consequence the methods can be located in a seperate file.

/// Construct bound parameters at the current position.
///
/// @param [in] freeParams free parameters vector
/// @param [in] boundCov bound covariance matrix
/// @param [in] covIsValid whether the covariance matrix contains valid entries
/// @param [in] surface surface on which the parameters are bound to
/// @param [in] geoCtx geometry context
/// @return bound parameters on the surface
///
/// @note Assumes that the position is already on the surface and covariance
///       (optionally) has already been transported.
BoundParameters makeBoundParameters(const FreeVector& freeParams,
                                    const BoundSymMatrix& boundCov,
                                    bool covIsValid, const Surface& surface,
                                    const GeometryContext& geoCtx);

/// Construct a curvilinear state at the current position
///
/// @param [in] freeParams free parameters vector
/// @param [in] curvilinearCov curvilinear covariance matrix
/// @param [in] covIsValid whether the covariance matrix contains valid entries
/// @return curvilinear parameters at the current position
///
/// @note Assumes that the covariance (optionally) has already been transported.
CurvilinearParameters makeCurvilinearParameters(
    const FreeVector& freeParams, const BoundSymMatrix& curvilinearCov,
    bool covIsValid);

/// Transport covariance to the bound state on an arbitrary surface.
///
/// @param [in,out] state The stepper state
/// @param [in] surface is the surface to which the covariance is forwarded to
///
/// @note This assumes that the current global parameter state has already been
///       propagated to be on the surface.
///
/// The Jacobians are reset such that they represent the propagation starting at
/// the surface.
void transportCovarianceToBound(StepperState& state, const Surface& surface);

/// Transport covariance to the bound state on the current curvilinear frame.
///
/// @param [in,out] state The stepper state
///
/// The Jacobians are reset such that they represent the propagation starting at
/// the current curvilinear frame.
void transportCovarianceToCurvilinear(StepperState& state);

}  // namespace detail
}  // namespace Acts
