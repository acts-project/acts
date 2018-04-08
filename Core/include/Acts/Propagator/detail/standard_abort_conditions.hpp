// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <limits>
#include "ACTS/Propagator/detail/constrained_step.hpp"
#include "ACTS/Utilities/Definitions.hpp"

#ifndef ABORTER_DEBUG_OUTPUTS
#define ABORTER_DEBUG_OUTPUTS
#define TARGETLOG(cache, status, message)                                      \
  if (debug) {                                                                 \
    std::stringstream dstream;                                                 \
    dstream << " " << status << " " << std::setw(cache.debug_pfx_width);       \
    dstream << " target aborter "                                              \
            << " | ";                                                          \
    dstream << std::setw(cache.debug_msg_width) << message << '\n';            \
    cache.debug_string += dstream.str();                                       \
  }
#endif

namespace Acts {

namespace detail {

  /// This is the condition that the pathLimit has been reached
  struct path_limit_reached
  {

    /// this is direction * absolute path limit
    double signed_path_limit = std::numeric_limits<double>::max();

    /// the tolerance used to defined "reached"
    double tolerance = 0.;

    /// write out debug information
    double debug = false;

    /// constructor
    ///
    /// @param tlimit is the signed path limit for this propagation
    /// @param ttolerance The tolerance to declare "reached"
    path_limit_reached(double tlimit     = std::numeric_limits<double>::max(),
                       double ttolerance = 0.001)
      : signed_path_limit(tlimit), tolerance(std::abs(ttolerance))
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
      double diff_to_limit = signed_path_limit - cache.accumulated_path;
      cache.step_size.update(diff_to_limit, constrained_step::aborter);
      bool limit_reached = (std::abs(diff_to_limit) < tolerance);
      if (limit_reached) {
        TARGETLOG(cache, "x", "Path limit reached.");
        // reaching the target means navigaiton break
        cache.target_reached = true;
      } else
        TARGETLOG(cache,
                  "o",
                  "Target step_size (path limit) updated to "
                      << cache.step_size.toString());
      // path limit check
      return limit_reached;
    }
  };

  /// This is the condition that the Surface has been reached
  /// it then triggers an propagation abort of the propgation
  template <typename Surface>
  struct surface_reached
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
    surface_reached()
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

      // check if the cache filled the current_surface
      if (cache.current_surface == surface) {
        TARGETLOG(cache, "x", "Target surface reached.");
        // reaching the target calls a navigation break
        cache.target_reached = true;
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
      cache.step_size.update(cache.nav_dir * distance,
                             constrained_step::aborter);
      // return true if you fall below tolerance
      bool targed_reached = (std::abs(distance) <= tolerance);
      if (targed_reached) {
        TARGETLOG(cache, "x", "Target surface reached.");
        // assigning the current_surface
        cache.current_surface = surface;
        TARGETLOG(cache,
                  "x",
                  "Current surface set to target surface "
                      << cache.current_surface->geoID().toString());
        // reaching the target calls a navigation break
        cache.target_reached = true;
      } else
        TARGETLOG(cache,
                  "o",
                  "Target step_size (surface) updated to "
                      << cache.step_size.toString());
      // path limit check
      return targed_reached;
    }
  };

}  // namespace detail
}  // namespace Acts