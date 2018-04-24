// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <limits>

namespace Acts {

namespace detail {

  /// This is the condition that the pathLimit has been reached
  struct path_limit_reached
  {

    /// this is direction * absolute path limit
    double signed_path_limit = std::numeric_limits<double>::max();

    /// the tolerance used to defined "reached"
    double tolerance = 0.;

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
      if (std::abs(cache.step_size)
          > std::abs(signed_path_limit - cache.accumulated_path))
        cache.step_size = signed_path_limit - cache.accumulated_path;

      // path limit check
      return (std::abs(signed_path_limit - cache.accumulated_path) < tolerance);
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
    int direction = 1;
    /// the tolerance to be defined on surface
    double tolerance = 0.;

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
      if (std::abs(cache.step_size) > std::abs(distance))
        cache.step_size = distance;
      // return true if you fall below tolerance
      return (std::abs(distance) <= tolerance);
    }
  };

}  // namespace detail
}  // namespace Acts