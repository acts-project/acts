// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_STANDARD_ABORT_CONDITIONS_HPP
#define ACTS_STANDARD_ABORT_CONDITIONS_HPP

#include <limits>

namespace Acts {

namespace detail {

  /// This is a default non-abort condition
  struct just_continue
  {

    /// boolean operator for abort condition using the result
    /// always return false
    template <typename input, typename result_t>
    bool
    operator()(const result_t&, input&, double&) const
    {
      return false;
    }

    /// boolean operator for abort condition without using the result
    /// always return false
    template <typename input>
    bool
    operator()(input&, double&) const
    {
      return false;
    }
  };

  /// This is the condition that the pathLimit has been reached
  struct path_limit_reached
  {

    /// this is direction * absolute path limit
    double singed_path_limit = std::numeric_limits<double>::max();
    /// the tolerance used to defined "reached"
    double tolerance = 0.;

    /// constructor
    ///
    /// @param tlimit is the signed path limit for this propagation
    /// @param ttolerance The tolerance to declare "reached"
    path_limit_reached(double tlimit     = std::numeric_limits<double>::max(),
                       double ttolerance = 0.001)
      : singed_path_limit(tlimit), tolerance(std::abs(ttolerance))
    {
    }

    /// constructor
    ///
    /// @param direction_sign is the propagation direction
    /// @param abs_limit is the absolute path limit for this propagation
    /// @param ttolerance The tolerance to declare "reached"
    path_limit_reached(int direction_sign, double abs_limit, double ttolerance)
      : singed_path_limit(direction_sign * std::abs(abs_limit))
      , tolerance(std::abs(ttolerance))
    {
    }

    /// boolean operator for abort condition using the result
    template <typename input, typename result_t>
    bool
    operator()(const result_t& r, input& cache, double& stepMax) const
    {
      return operator()(cache, stepMax);
    }

    /// boolean operator for abort condition without using the result
    /// @param cache The propagation cache
    /// @param stepMax Maximum step for the propagation cache it might
    ///        be adapted to the remainim path length
    template <typename input>
    bool
    operator()(input& cache, double& stepMax) const
    {
      // Check if the maximum allowed step size has to be updated
      if (std::abs(stepMax)
          > std::abs(singed_path_limit - cache.accumulated_path))
        stepMax = singed_path_limit - cache.accumulated_path;
      // path limit check
      return (std::abs(singed_path_limit - cache.accumulated_path) < tolerance);
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
    /// the tolerance to be defined on surface
    double tolerance = 0.;

    /// constructor
    ///
    /// @tparam Surface Type of the surface
    ///
    /// @param tsurface The target surface
    /// @param ttolerance The tolerance to declare "reached"
    surface_reached(const Surface& tsurface, double ttolerance = 0.)
      : surface(&tsurface), tolerance(ttolerance)
    {
    }

    /// boolean operator for abort condition using the result (ignored)
    template <typename input, typename result_t>
    bool
    operator()(const result_t&, input& cache, double& stepMax) const
    {
      return operator()(cache, stepMax);
    }

    /// boolean operator for abort condition without using the result
    /// @param cache The propagation cache
    /// @param stepMax Maximum step for the propagation cache it might
    ///        be adapted to the remainim path length
    template <typename input>
    bool
    operator()(input& cache, double& stepMax) const
    {
      // calculate the distance to the surface
      const double distance
          = surface->intersectionEstimate(cache.position(), cache.direction())
                .pathLength;
      // Adjust the step size so that we cannot cross the target surface
      if (std::abs(stepMax) > std::abs(distance)) stepMax = distance;
      return (std::abs(distance) <= tolerance);
    }
  };

}  // namespace detail
}  // namespace Acts

#endif
