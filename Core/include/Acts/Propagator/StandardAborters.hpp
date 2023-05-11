// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <limits>
#include <sstream>
#include <string>

namespace Acts {

/// @brief TargetOptions struct for geometry interface
struct TargetOptions {
  /// Navigation direction
  Direction navDir = Direction::Forward;

  /// Target Boundary check directive - always false here
  BoundaryCheck boundaryCheck = false;

  /// Object to check against - always nullptr here
  const Surface* startObject = nullptr;

  /// The path limit
  double pathLimit = std::numeric_limits<double>::max();

  /// create target options
  TargetOptions(Direction ndir) : navDir(ndir) {}
};

/// This is the condition that the pathLimit has been reached
struct PathLimitReached {
  /// Boolean switch for Loop protection
  double internalLimit = std::numeric_limits<double>::max();

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for propagation
  /// @param [in] navigator Navigator used for propagation
  /// @param logger a logger instance
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, const Logger& logger) const {
    if (navigator.targetReached(state.navigation)) {
      return true;
    }
    // Check if the maximum allowed step size has to be updated
    double distance = state.stepping.navDir * std::abs(internalLimit) -
                      state.stepping.pathAccumulated;
    double tolerance = state.options.targetTolerance;
    stepper.setStepSize(state.stepping, distance, ConstrainedStep::aborter,
                        false);
    bool limitReached = (std::abs(distance) < std::abs(tolerance));
    if (limitReached) {
      ACTS_VERBOSE("Target: x | "
                   << "Path limit reached at distance " << distance);
      // reaching the target means navigation break
      navigator.targetReached(state.navigation, true);
    } else {
      ACTS_VERBOSE("Target: 0 | "
                   << "Target stepSize (path limit) updated to "
                   << stepper.outputStepSize(state.stepping));
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
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for propagation
  /// @param [in] navigator Navigator used for propagation
  /// @param logger a logger instance
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, const Logger& logger) const {
    return (*this)(state, stepper, navigator,
                   *navigator.targetSurface(state.navigation), logger);
  }

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for the progation
  /// @param [in] navigator Navigator used for the progation
  /// @param [in] targetSurface The target surface
  /// @param logger a logger instance
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, const Surface& targetSurface,
                  const Logger& logger) const {
    if (navigator.targetReached(state.navigation)) {
      return true;
    }
    // Check if the cache filled the currentSurface - or if we are on the
    // surface
    if ((navigator.currentSurface(state.navigation) &&
         navigator.currentSurface(state.navigation) == &targetSurface)) {
      ACTS_VERBOSE("Target: x | "
                   << "Target surface reached.");
      // reaching the target calls a navigation break
      navigator.targetReached(state.navigation, true);
      return true;
    }
    // Calculate the distance to the surface
    const double tolerance = state.options.targetTolerance;
    const auto sIntersection = targetSurface.intersect(
        state.geoContext, stepper.position(state.stepping),
        state.stepping.navDir * stepper.direction(state.stepping), true);

    // The target is reached
    bool targetReached = (sIntersection.intersection.status ==
                          Intersection3D::Status::onSurface);
    double distance = sIntersection.intersection.pathLength;

    // Return true if you fall below tolerance
    if (targetReached) {
      ACTS_VERBOSE("Target: x | "
                   << "Target surface reached at distance (tolerance) "
                   << distance << " (" << tolerance << ")");
      // assigning the currentSurface
      navigator.currentSurface(state.navigation, &targetSurface);
      ACTS_VERBOSE("Target: x | "
                   << "Current surface set to target surface  "
                   << navigator.currentSurface(state.navigation)->geometryId());

      // reaching the target calls a navigation break
      navigator.targetReached(state.navigation, true);
    } else {
      // Target is not reached, update the step size
      const double overstepLimit = stepper.overstepLimit(state.stepping);
      // Check the alternative solution
      if (distance < overstepLimit and sIntersection.alternative) {
        // Update the distance to the alternative solution
        distance = sIntersection.alternative.pathLength;
      }
      stepper.setStepSize(state.stepping, state.stepping.navDir * distance,
                          ConstrainedStep::aborter, false);

      ACTS_VERBOSE("Target: 0 | "
                   << "Target stepSize (surface) updated to "
                   << stepper.outputStepSize(state.stepping));
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
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] navigator The navigator object
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool operator()(propagator_state_t& state, const stepper_t& /*stepper*/,
                  const navigator_t& navigator,
                  const Logger& /*logger*/) const {
    if (navigator.currentVolume(state.navigation) != nullptr) {
      return false;
    }
    navigator.targetReached(state.navigation, true);
    return true;
  }
};

/// If the particle stopped (p=0) abort the propagation
struct ParticleStopped {
  ParticleStopped() = default;

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper The stepper object
  /// @param [in] navigator The navigator object
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator,
                  const Logger& /*logger*/) const {
    if (stepper.momentum(state.stepping) > 0) {
      return false;
    }
    navigator.targetReached(state.navigation, true);
    return true;
  }
};

}  // namespace Acts
