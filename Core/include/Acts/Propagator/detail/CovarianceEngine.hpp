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

/// Construct a bound state at the current position.
///
/// @param [in] state State that will be presented as @c BoundState
/// @param [in] surface The surface to which we bind the state
/// @return A bound state:
///   - the bound parameters at the surface
///   - the stepwise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
///
/// @note Assumes that the position is already on the surface and covariance
///       (optionally) has already been transported.
std::tuple<BoundParameters, BoundMatrix, double> boundState(
    const StepperState& state, const Surface& surface);

/// Construct a curvilinear state at the current position
///
/// @param [in] state State that will be presented as @c CurvilinearState
/// @return A curvilinear state:
///   - the curvilinear parameters at given position
///   - the stepwise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
///
/// @note Assumes that the covariance (optionally) has already been transported.
std::tuple<CurvilinearParameters, BoundMatrix, double> curvilinearState(
    const StepperState& state);

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
