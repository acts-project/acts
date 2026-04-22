// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/helix.hpp"
#include "detray/utils/root_finding.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t>
struct helix_intersector_impl;

/// @brief Intersection implementation for helical trajectories with planar
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <detray::concepts::aos algebra_t>
struct helix_intersector_impl<cartesian2D<algebra_t>, algebra_t> {
 private:
  using scalar_t = dscalar<algebra_t>;

 public:
  using algebra_type = algebra_t;

  template <typename surface_descr_t>
  using intersection_type =
      intersection2D<surface_descr_t, algebra_t, intersection::contains_pos>;

  template <typename other_algebra_t>
  using trajectory_type = detail::helix<other_algebra_t>;

  // Maximum number of solutions this intersector can produce
  static constexpr std::uint8_t n_solutions{1u};

  using result_type = intersection_point_err<algebra_t>;

  /// Operator function to find intersections between a helix and a planar
  /// surface
  ///
  /// @param h is the input helix trajectory
  /// @param trf is the surface placement transform
  ///
  /// @return the intersection
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const trajectory_type<algebra_t> &h, const dtransform3D<algebra_t> &trf,
      const scalar_t = 0.f) const {
    using point3_t = dpoint3D<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    // Surface normal
    const vector3_t sn = trf.z();
    // Surface translation
    const point3_t st = trf.translation();

    // Starting point on the helix for the Newton iteration
    const vector3_t dist{st - h.pos(0.f)};
    scalar_t denom{vector::dot(sn, h.dir(0.5f * vector::norm(dist)))};
    scalar_t s_ini;
    if (denom == 0.f) {
      s_ini = vector::norm(dist);
    } else {
      s_ini = vector::dot(sn, dist) / denom;
    }

    /// Evaluate the function and its derivative at the point @param x
    auto plane_inters_func = [&h, &st, &sn](const scalar_t x) {
      // f(s) = sn * (h.pos(s) - st) == 0
      const scalar_t f_s{vector::dot(sn, (h.pos(x) - st))};
      // f'(s) = sn * h.dir(s)
      const scalar_t df_s{vector::dot(sn, h.dir(x))};

      return std::make_tuple(f_s, df_s);
    };

    // Run the root finding algorithm
    scalar_t s;
    scalar_t ds;
    if (run_rtsafe) {
      std::tie(s, ds) =
          newton_raphson_safe(plane_inters_func, s_ini, convergence_tolerance,
                              max_n_tries, max_path);
    } else {
      std::tie(s, ds) =
          newton_raphson(plane_inters_func, s_ini, convergence_tolerance,
                         max_n_tries, max_path);
    }

    return {s, h.pos(s), ds};
  }

  /// Tolerance for convergence
  scalar_t convergence_tolerance{1.f * unit<scalar_t>::um};
  // Guard against infinite loops
  std::size_t max_n_tries{1000u};
  // Early exit, if the intersection is too far away
  scalar_t max_path{5.f * unit<scalar_t>::m};
  // Complement the Newton algorithm with Bisection steps
  bool run_rtsafe{true};
};

template <detray::concepts::aos algebra_t>
struct helix_intersector_impl<polar2D<algebra_t>, algebra_t>
    : public helix_intersector_impl<cartesian2D<algebra_t>, algebra_t> {};

}  // namespace detray
