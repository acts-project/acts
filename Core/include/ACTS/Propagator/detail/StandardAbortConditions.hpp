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
#include "ACTS/Propagator/detail/ConstrainedStep.hpp"
#include "ACTS/Utilities/Definitions.hpp"

#ifndef ABORTER_DEBUG_OUTPUTS
#define ABORTER_DEBUG_OUTPUTS
#define TARGETLOG(cache, status, message)                                      \
  if (debug) {                                                                 \
    std::stringstream dstream;                                                 \
    dstream << " " << status << " " << std::setw(cache.debugPfxWidth);         \
    dstream << " target aborter "                                              \
            << " | ";                                                          \
    dstream << std::setw(cache.debugMsgWidth) << message << '\n';              \
    cache.debugString += dstream.str();                                        \
  }
#endif

namespace Acts {

namespace detail {

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
        TARGETLOG(cache, "x", "Path limit reached.");
        // reaching the target means navigaiton break
        cache.targetReached = true;
      } else
        TARGETLOG(cache,
                  "o",
                  "Target stepSize (path limit) updated to "
                      << cache.stepSize.toString());
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
        TARGETLOG(cache, "x", "Target surface reached.");
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
        TARGETLOG(cache, "x", "Target surface reached.");
        // assigning the currentSurface
        cache.currentSurface = surface;
        TARGETLOG(cache,
                  "x",
                  "Current surface set to target surface "
                      << cache.currentSurface->geoID().toString());
        // reaching the target calls a navigation break
        cache.targetReached = true;
      } else
        TARGETLOG(cache,
                  "o",
                  "Target stepSize (surface) updated to "
                      << cache.stepSize.toString());
      // path limit check
      return targetReached;
    }
  };

}  // namespace detail
}  // namespace Acts

#endif
