// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/line2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/helix.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/root_finding.hpp"

namespace detray {

template <typename frame_t, concepts::algebra algebra_t>
struct helix_intersector_impl;

/// @brief Intersection implementation for helical trajectories with line
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <detray::concepts::aos algebra_t>
struct helix_intersector_impl<line2D<algebra_t>, algebra_t> {
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

  /// Operator function to find intersections between a helix and a line
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

    // line axis direction
    const vector3_t l = getter::vector<3>(trf.matrix(), 0u, 2u);

    // line center
    const point3_t c = trf.translation();

    // initial track direction
    const vector3_t t0 = h.dir(0.f);

    // initial track position
    const point3_t r0 = h.pos(0.f);

    // Projection of line to track direction
    const scalar_t lt0{vector::dot(l, t0)};

    const scalar_t denom{1.f - (lt0 * lt0)};

    // Case for wire is parallel to track
    // @NOTE We might not have to call this which is meant to be for ray
    // intersection...
    if (denom < 1e-5f) {
      DETRAY_ERROR_HOST("Helix line intersector encountered invalid value!");
      return {};
    }

    // vector from track position to line center
    const vector3_t D = c - r0;

    // D projection on line direction
    const scalar_t P{vector::dot(D, l)};

    // D projection on track direction
    const scalar_t Q{vector::dot(D, t0)};

    // Path length to the point of closest approach on the track
    // @NOTE Ray intersection algorithm is used for the initial guess on
    // the path length
    scalar_t s_ini{1.f / denom * (Q - P * lt0)};

    /// Evaluate the function and its derivative at the point @param x
    auto line_inters_func = [&h, &c, &l](const scalar_t x) {
      // track direction
      const vector3_t t = h.dir(x);

      // track position
      const point3_t r = h.pos(x);

      // Projection of (track position - center) to the line
      const scalar_t A = vector::dot(r - c, l);

      // Vector orthogonal to the line and passing the track position
      // w = r - (c + ((r - c) * l)l)
      const vector3_t w = r - (c + A * l);

      // f(s) = t * w = 0
      const scalar_t f = vector::dot(t, w);

      // dtds = d^2r/ds^2 = qop * (t X b_field)
      const vector3_t dtds = h.qop() * vector::cross(t, h.b_field());
      // dwds = t - (t * l)l
      const vector3_t dwds = t - vector::dot(t, l) * l;

      // f'(s) = dtds * w + t * dwds
      const scalar_t dfds = vector::dot(dtds, w) + vector::dot(t, dwds);

      return std::make_tuple(f, dfds);
    };

    // Run the root finding algorithm
    scalar_t s;
    scalar_t ds;
    if (run_rtsafe) {
      std::tie(s, ds) =
          newton_raphson_safe(line_inters_func, s_ini, convergence_tolerance,
                              max_n_tries, max_path);
    } else {
      std::tie(s, ds) =
          newton_raphson(line_inters_func, s_ini, convergence_tolerance,
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

}  // namespace detray
