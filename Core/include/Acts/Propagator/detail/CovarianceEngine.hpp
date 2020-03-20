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

/// Create and return the bound state at the current position.
///
/// @param [in] state State that will be presented as @c BoundState
/// @param [in] surface The surface to which we bind the state
/// @return A bound state:
///   - the parameters at the surface
///   - the stepwise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
///
/// @warning Calls the covariance transport internally and modifies the state
///          accordingly on every call, but only if the covariance transport is
///          enabled. Multiple calls at the same position can thus lead to
///          ill-defined outputs.
std::tuple<BoundParameters, BoundMatrix, double> boundState(
    StepperState& state, const Surface& surface);

/// Create and return a curvilinear state at the current position
///
/// @param [in] state State that will be presented as @c CurvilinearState
/// @return A curvilinear state:
///   - the curvilinear parameters at given position
///   - the stepweise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
///
/// @warning Calls the covariance transport internally and modifies the state
///          accordingly on every call, but only if the covariance transport is
///          enabled. Multiple calls at the same position can thus lead to
///          ill-defined outputs.
std::tuple<CurvilinearParameters, BoundMatrix, double> curvilinearState(
    StepperState& state);

/// @brief Method for on-demand transport of the covariance to a new frame at
/// current position in parameter space
///
/// @param [in,out] state The stepper state
/// @param [in] surface is the surface to which the covariance is
///        forwarded to
/// @note No check is done if the position is actually on the surface
///
/// @return Projection jacobian from global to bound parameters
void covarianceTransport(StepperState& state, const Surface* surface = nullptr);

}  // namespace detail
}  // namespace Acts