// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <sstream>
#include <string>
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

/// @brief TargetOptions struct for geometry interface
struct TargetOptions {
  /// Navigation direction
  NavigationDirection navDir = forward;

  /// Target Boundary check directive - always false here
  BoundaryCheck boundaryCheck = false;

  /// Object to check against - always nullptr here
  const Surface* startObject = nullptr;

  /// The path limit
  double pathLimit = std::numeric_limits<double>::max();

  /// create target options
  TargetOptions(NavigationDirection ndir) : navDir(ndir) {}
};

/// The debug logging for standard aborters
///
/// It needs to be fed by a lambda function that returns a string,
/// that guarantees that the lambda is only called in the
/// state.options.debug == true case in order not to spend time
/// when not needed.
///
/// @param state the propagator cache for the debug flag, prefix/stream
/// @param logAction is a callable function that returns a streamable object
template <typename propagator_state_t>
void targetDebugLog(propagator_state_t& state, const std::string& status,
                    const std::function<std::string()>& logAction) {
  if (state.options.debug) {
    std::stringstream dstream;
    dstream << " " << status << " ";
    dstream << std::setw(state.options.debugPfxWidth);
    dstream << " Target "
            << " | ";
    dstream << std::setw(state.options.debugMsgWidth);
    dstream << logAction() << '\n';
    state.options.debugString += dstream.str();
  }
}

/// This is the condition that the pathLimit has been reached
struct PathLimitReached {
  /// Boolean switch for Loop protection
  double internalLimit = std::numeric_limits<double>::max();

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in,out] state The propagation state object
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& state,
                  const stepper_t& /*unused*/) const {
    if (state.navigation.targetReached) {
      return true;
    }
    // Check if the maximum allowed step size has to be updated
    double distance = state.stepping.navDir * std::abs(internalLimit) -
                      state.stepping.pathAccumulated;
    double tolerance = state.options.targetTolerance;
    state.stepping.stepSize.update(distance, ConstrainedStep::aborter);
    bool limitReached = (distance * distance < tolerance * tolerance);
    if (limitReached) {
      targetDebugLog(state, "x", [&] {
        std::stringstream dstream;
        dstream << "Path limit reached at distance " << distance;
        return dstream.str();
      });
      // reaching the target means navigation break
      state.navigation.targetReached = true;
    } else {
      targetDebugLog(state, "o", [&] {
        std::stringstream dstream;
        dstream << "Target stepSize (path limit) updated to ";
        dstream << state.stepping.stepSize.toString();
        return dstream.str();
      });
    }
    // path limit check
    return limitReached;
  }
};

/// This is the condition that the Surface has been reached
/// it then triggers an propagation abort of the propagation
struct SurfaceReached {
  SurfaceReached() = default;

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for propagation
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper) const {
    return (*this)(state, stepper, *state.navigation.targetSurface);
  }

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for the progation
  /// @param [in] targetSurface The target surface
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const Surface& targetSurface) const {
    if (state.navigation.targetReached) {
      return true;
    }
    // Check if the cache filled the currentSurface - or if we are on the
    // surface
    // @todo - do not apply the isOnSurface check here, but handle by the
    // intersectionEstimate
    if ((state.navigation.currentSurface &&
         state.navigation.currentSurface == &targetSurface)) {
      targetDebugLog(state, "x", [&] {
        std::string ds("Target surface reached.");
        return ds;
      });
      // reaching the target calls a navigation break
      state.navigation.targetReached = true;
      return true;
    }
    // Calculate the distance to the surface
    const double tolerance = state.options.targetTolerance;
    const auto sIntersection = targetSurface.intersect(
        state.geoContext, stepper.position(state.stepping),
        state.stepping.navDir * stepper.direction(state.stepping), true);

    // The target is reached
    bool targetReached =
        (sIntersection.intersection.status == Intersection::Status::onSurface);
    double distance = sIntersection.intersection.pathLength;

    // Return true if you fall below tolerance
    if (targetReached) {
      targetDebugLog(state, "x", [&] {
        std::stringstream dstream;
        dstream << "Target surface reached at distance (tolerance) ";
        dstream << distance << " (" << tolerance << ")";
        return dstream.str();
      });
      // assigning the currentSurface
      state.navigation.currentSurface = &targetSurface;
      targetDebugLog(state, "x", [&] {
        std::stringstream dstream;
        dstream << "Current surface set to target surface  ";
        dstream << state.navigation.currentSurface->geoID();
        return dstream.str();
      });
      // reaching the target calls a navigation break
      state.navigation.targetReached = true;
    } else {
      // Target is not reached, update the step size
      const double overstepLimit = stepper.overstepLimit(state.stepping);
      // Check the alternative solution
      if (distance < overstepLimit and sIntersection.alternative) {
        // Update the distance to the alternative solution
        distance = sIntersection.alternative.pathLength;
      }
      state.stepping.stepSize.update(state.stepping.navDir * distance,
                                     ConstrainedStep::aborter);
      targetDebugLog(state, "o", [&] {
        std::stringstream dstream;
        dstream << "Target stepSize (surface) updated to ";
        dstream << state.stepping.stepSize.toString();
        return dstream.str();
      });
    }
    // path limit check
    return targetReached;
  }
};

/// This is the condition that the end of World has been reached
/// it then triggers an propagation abort
struct EndOfWorldReached {
  EndOfWorldReached() = default;

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param [in,out] state The propagation state object
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& state,
                  const stepper_t& /*unused*/) const {
    if (state.navigation.currentVolume != nullptr) {
      return false;
    }
    state.navigation.targetReached = true;
    return true;
  }
};

}  // namespace Acts
