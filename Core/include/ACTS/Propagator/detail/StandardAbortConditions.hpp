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
  /// that guarantees that the lambda is only called in the cache.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @param cache the stepper cache for the debug flag, prefix and length
  /// @param logAction is a callable function that returns a stremable object
  template <typename cache_t>
  void
  targetDebugLog(cache_t&                     cache,
                 std::string                  status,
                 std::function<std::string()> logAction)
  {
    if (cache.debug) {
      std::stringstream dstream;
      dstream << " " << status << " " << std::setw(cache.debugPfxWidth);
      dstream << " target aborter "
              << " | ";
      dstream << std::setw(cache.debugMsgWidth) << logAction() << '\n';
      cache.debugString += dstream.str();
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
    template <typename cache_t, typename result_t>
    bool
    operator()(const result_t& /*r*/, cache_t& cache) const
    {
      return operator()(cache);
    }

    /// boolean operator for abort condition without using the result
    /// @param cache The propagation cache
    /// @param stepMax Maximum step for the propagation cache it might
    ///        be adapted to the remainim path length
    template <typename cache_t>
    bool
    operator()(cache_t& cache) const
    {
      // Check if the maximum allowed step size has to be updated
      double diffToLimit = signedPathLimit - cache.accumulatedPath;
      cache.stepSize.update(diffToLimit, ConstrainedStep::aborter);
      bool limitReached = (std::abs(diffToLimit) < tolerance);
      if (limitReached) {
        targetDebugLog(cache, "x", [&] {
          std::stringstream dstream;
          dstream << "Path limit reached at distance " << diffToLimit;
          return dstream.str();
        });
        // reaching the target means navigaiton break
        cache.targetReached = true;
      } else
        targetDebugLog(cache, "o", [&] {
          std::stringstream dstream;
          dstream << "Target stepSize (path limit) updated to ";
          dstream << cache.stepSize.toString();
          return dstream.str();
        });
      // path limit check
      return limitReached;
    }
  };

  /// This is the condition that the Surface has been reached
  /// it then triggers an propagation abort of the propgation
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

    /// constructor
    ///
    /// @tparam Surface Type of the surface
    ///
    /// @param tsurface The target surface
    /// @param ttolerance The tolerance to declare "reached"
    SurfaceReached()
      : surface(nullptr), tolerance(std::numeric_limits<double>::max())
    {
    }

    /// boolean operator for abort condition using the result (ignored)
    template <typename cache_t, typename result_t>
    bool
    operator()(const result_t&, cache_t& cache) const
    {
      return operator()(cache);
    }

    /// boolean operator for abort condition without using the result
    /// @param cache The propagation cache
    /// @param stepMax Maximum step for the propagation cache it might
    ///        be adapted to the remainim path length
    template <typename cache_t>
    bool
    operator()(cache_t& cache) const
    {
      if (!surface) return false;

      // check if the cache filled the currentSurface
      if (cache.currentSurface == surface) {
        targetDebugLog(cache, "x", [&] {
          std::string ds("Target surface reached.");
          return ds;
        });
        // reaching the target calls a navigation break
        cache.targetReached = true;
        return true;
      }

      // calculate the distance to the surface
      // @todo that might cause problems with a cylinder
      const double distance
          = surface
                ->intersectionEstimate(cache.position(),
                                       direction * cache.direction(),
                                       true,
                                       false)
                .pathLength;
      // Adjust the step size so that we cannot cross the target surface
      cache.stepSize.update(cache.navDir * distance, ConstrainedStep::aborter);
      // return true if you fall below tolerance
      bool targetReached = (std::abs(distance) <= tolerance);
      if (targetReached) {
        targetDebugLog(cache, "x", [&] {
          std::stringstream dstream;
          dstream << "Target surface reached at distance (tolerance) ";
          dstream << distance << " (" << tolerance << ")";
          return dstream.str();
        });
        // assigning the currentSurface
        cache.currentSurface = surface;
        targetDebugLog(cache, "x", [&] {
          std::stringstream dstream;
          dstream << "Current surface set to target surface  ";
          dstream << cache.currentSurface->geoID().toString();
          return dstream.str();
        });
        // reaching the target calls a navigation break
        cache.targetReached = true;
      } else
        targetDebugLog(cache, "o", [&] {
          std::stringstream dstream;
          dstream << "Target stepSize (surface) updated to ";
          dstream << cache.stepSize.toString();
          return dstream.str();
        });
      // path limit check
      return targetReached;
    }
  };

}  // namespace detail
}  // namespace Acts

#endif
