// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <sstream>
#include <string>
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

namespace detail {

  /// @brief TargetOptions struct for geometry interface
  struct TargetOptions
  {

    /// the navigation direction
    NavigationDirection navDir = forward;

    /// the boundary check directive - always false here
    BoundaryCheck boundaryCheck = false;

    /// object to check against - always nullptr here
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
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_state_t>
  void
  targetDebugLog(propagator_state_t&          state,
                 std::string                  status,
                 std::function<std::string()> logAction)
  {
    if (state.options.debug) {
      std::stringstream dstream;
      dstream << " " << status << " ";
      dstream << std::setw(state.options.debugPfxWidth);
      dstream << " target aborter "
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth);
      dstream << logAction() << '\n';
      state.options.debugString += dstream.str();
    }
  }

  /// This is the condition that the pathLimit has been reached
  struct PathLimitReached
  {

    /// boolean operator for abort condition using the result
    template <typename propagator_state_t, typename result_t>
    bool
    operator()(const result_t& /*r*/, propagator_state_t& state) const
    {
      return operator()(state);
    }

    /// boolean operator for abort condition without using the result
    ///
    /// @tparam propagator_state_t Type of the propagator state
    ///
    /// @param[in,out] state The propagation state object
    template <typename propagator_state_t>
    bool
    operator()(propagator_state_t& state) const
    {
      // Check if the maximum allowed step size has to be updated
      double distance
          = state.options.pathLimit - state.stepping.pathAccumulated;
      double tolerance = state.options.targetTolerance;
      state.stepping.stepSize.update(distance, ConstrainedStep::aborter);
      bool limitReached = (distance * distance < tolerance * tolerance);
      if (limitReached) {
        targetDebugLog(state, "x", [&] {
          std::stringstream dstream;
          dstream << "Path limit reached at distance " << distance;
          return dstream.str();
        });
        // reaching the target means navigaiton break
        state.navigation.targetReached = true;
      } else
        targetDebugLog(state, "o", [&] {
          std::stringstream dstream;
          dstream << "Target stepSize (path limit) updated to ";
          dstream << state.stepping.stepSize.toString();
          return dstream.str();
        });
      // path limit check
      return limitReached;
    }
  };

  /// This is the condition that the Surface has been reached
  /// it then triggers an propagation abort of the propagation
  struct SurfaceReached
  {
    /// Default Constructor
    SurfaceReached() {}

    /// boolean operator for abort condition using the result (ignored)
    template <typename propagator_state_t, typename result_t>
    bool
    operator()(const result_t&, propagator_state_t& state) const
    {
      return operator()(state);
    }

    /// boolean operator for abort condition without using the result
    ///
    /// @tparam propagator_state_t Type of the propagator state
    ///
    /// @param[in,out] state The propagation state object
    template <typename propagator_state_t>
    bool
    operator()(propagator_state_t& state) const
    {
      if (state.navigation.targetReached) return true;

      // check if the cache filled the currentSurface
      if (state.navigation.currentSurface
          && state.navigation.currentSurface
              == state.navigation.targetSurface) {
        targetDebugLog(state, "x", [&] {
          std::string ds("Target surface reached.");
          return ds;
        });
        // reaching the target calls a navigation break
        state.navigation.targetReached = true;
        return true;
      }
      // calculate the distance to the surface
      // @todo: add corrector
      const double tolerance = state.options.targetTolerance;
      const auto   iestimate
          = state.navigation.targetSurface->intersectionEstimate(
              state.stepping, TargetOptions(state.options.direction), nullptr);
      const double distance = iestimate.intersection.pathLength;
      // Adjust the step size so that we cannot cross the target surface
      state.stepping.stepSize.update(distance, ConstrainedStep::aborter);
      // return true if you fall below tolerance
      bool targetReached = (distance * distance <= tolerance * tolerance);
      if (targetReached) {
        targetDebugLog(state, "x", [&] {
          std::stringstream dstream;
          dstream << "Target surface reached at distance (tolerance) ";
          dstream << distance << " (" << tolerance << ")";
          return dstream.str();
        });
        // assigning the currentSurface
        state.navigation.currentSurface = state.navigation.targetSurface;
        targetDebugLog(state, "x", [&] {
          std::stringstream dstream;
          dstream << "Current surface set to target surface  ";
          dstream << state.navigation.currentSurface->geoID().toString();
          return dstream.str();
        });
        // reaching the target calls a navigation break
        state.navigation.targetReached = true;
      } else {
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

  /// This is the condition that the Surface has been reached
  /// it then triggers an propagation abort of the propagation
  struct EndOfWorldReached
  {
    /// Default Constructor
    EndOfWorldReached() {}

    /// boolean operator for abort condition using the result (ignored)
    template <typename propagator_state_t, typename result_t>
    bool
    operator()(const result_t&, propagator_state_t& state) const
    {
      return operator()(state);
    }

    /// boolean operator for abort condition without using the result
    ///
    /// @tparam propagator_state_t Type of the propagator state
    ///
    /// @param[in,out] state The propagation state object
    template <typename propagator_state_t>
    bool
    operator()(propagator_state_t& state) const
    {
      if (state.navigation.currentVolume) return false;
      state.navigation.targetReached = true;
      return true;
    }
  };

}  // namespace detail
}  // namespace Acts
