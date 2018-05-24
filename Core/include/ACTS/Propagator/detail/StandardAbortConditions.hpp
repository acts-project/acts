// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_STANDARD_ABORT_CONDITIONS_HPP
#define ACTS_STANDARD_ABORT_CONDITIONS_HPP

#include <limits>
#include <sstream>
#include <string>
#include "ACTS/Propagator/detail/ConstrainedStep.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

namespace detail {

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

    /// this is direction * absolute path limit
    double signedPathLimit = std::numeric_limits<double>::max();

    /// the tolerance used to defined "reached"
    double tolerance = 0.;

    /// write out debug information
    double debug = false;

    /// constructor
    ///
    /// @param tlimit is the signed path limit for this propagation
    /// @param ttolerance The tolerance to declare "reached"
    PathLimitReached(double tlimit     = std::numeric_limits<double>::max(),
                     double ttolerance = 0.001)
      : signedPathLimit(tlimit), tolerance(std::abs(ttolerance))
    {
    }

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
      double diffToLimit = signedPathLimit - state.stepping.accumulatedPath;
      state.stepping.stepSize.update(diffToLimit, ConstrainedStep::aborter);
      bool limitReached = (std::abs(diffToLimit) < tolerance);
      if (limitReached) {
        targetDebugLog(state, "x", [&] {
          std::stringstream dstream;
          dstream << "Path limit reached at distance " << diffToLimit;
          return dstream.str();
        });
        // reaching the target means navigaiton break
        state.targetReached = true;
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
  template <typename Surface>
  struct SurfaceReached
  {
    /// the plain pointer to the surface
    /// - safe as the condition lives shorter than the surface
    const Surface* surface = nullptr;

    /// the direction
    NavigationDirection direction = forward;
    /// the tolerance to be defined on surface
    double tolerance = 0.;
    /// output debug
    bool debug = false;

    /// Default Constructor
    ///
    SurfaceReached()
      : surface(nullptr), tolerance(std::numeric_limits<double>::max())
    {
    }

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
      if (!surface) return false;

      // check if the cache filled the currentSurface
      if (state.currentSurface == surface) {
        targetDebugLog(state, "x", [&] {
          std::string ds("Target surface reached.");
          return ds;
        });
        // reaching the target calls a navigation break
        state.targetReached = true;
        return true;
      }

      // calculate the distance to the surface
      const double distance
          = surface
                ->intersectionEstimate(state.stepping.position(),
                                       state.stepping.direction(),
                                       direction,
                                       false)
                .pathLength;
      // Adjust the step size so that we cannot cross the target surface
      state.stepping.stepSize.update(distance, ConstrainedStep::aborter);
      // return true if you fall below tolerance
      bool targetReached = (std::abs(distance) <= tolerance);
      if (targetReached) {
        targetDebugLog(state, "x", [&] {
          std::stringstream dstream;
          dstream << "Target surface reached at distance (tolerance) ";
          dstream << distance << " (" << tolerance << ")";
          return dstream.str();
        });
        // assigning the currentSurface
        state.currentSurface = surface;
        targetDebugLog(state, "x", [&] {
          std::stringstream dstream;
          dstream << "Current surface set to target surface  ";
          dstream << state.currentSurface->geoID().toString();
          return dstream.str();
        });
        // reaching the target calls a navigation break
        state.targetReached = true;
      } else
        targetDebugLog(state, "o", [&] {
          std::stringstream dstream;
          dstream << "Target stepSize (surface) updated to ";
          dstream << state.stepping.stepSize.toString();
          return dstream.str();
        });
      // path limit check
      return targetReached;
    }
  };

}  // namespace detail
}  // namespace Acts

#endif
