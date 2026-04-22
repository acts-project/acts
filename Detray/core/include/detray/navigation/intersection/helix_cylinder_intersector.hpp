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
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/geometry/coordinates/cylindrical2D.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_cylinder_intersector.hpp"
#include "detray/tracks/helix.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/root_finding.hpp"

namespace detray {

template <typename frame_t, concepts::algebra algebra_t>
struct helix_intersector_impl;

/// @brief Intersection implementation for cylinder surfaces using helical
/// trajectories.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
/// @note Don't use for low p_t tracks!
template <detray::concepts::aos algebra_t>
struct helix_intersector_impl<cylindrical2D<algebra_t>, algebra_t>
    : public ray_intersector_impl<cylindrical2D<algebra_t>, algebra_t,
                                  intersection::contains_pos> {
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
  static constexpr std::uint8_t n_solutions{2u};

  using result_type = darray<intersection_point_err<algebra_t>, n_solutions>;

  /// Operator function to find intersections between a helix and a cylinder
  /// surface
  ///
  /// @param h is the input helix trajectory
  /// @param trf is the surface placement transform
  /// @param mask is the input mask that defines the surface extent
  ///
  /// @return the intersection
  template <typename mask_t>
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const trajectory_type<algebra_t> &h, const dtransform3D<algebra_t> &trf,
      const mask_t &mask, const scalar_t = 0.f) const {
    using point3_t = dpoint3D<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    result_type ret{};

    // Cylinder z axis
    const vector3_t sz = trf.z();
    // Cylinder centre
    const point3_t sc = trf.translation();

    // Starting point on the helix for the Newton iteration
    // The mask is a cylinder -> it provides its radius as the first
    // value
    const scalar_t r{mask[cylinder2D::e_r]};

    // Try to guess the best starting positions for the iteration

    // Direction of the track at the helix origin
    const auto h_dir = h.dir(0.5f * r);
    // Default starting path length for the Newton iteration (assumes
    // concentric cylinder)
    const scalar_t default_s{r * vector::perp(h_dir)};

    // Initial helix path length parameter
    darray<scalar_t, 2> paths{default_s, default_s};

    // try to guess good starting path by calculating the intersection
    // path of the helix tangential with the cylinder. This only has a
    // chance of working for tracks with reasonably high p_T !
    detail::ray<algebra_t> t{h.pos(), h.time(), h_dir, h.qop()};
    const auto qe = this->solve_intersection(t, mask, trf);

    // Obtain both possible solutions by looping over the (different)
    // starting positions
    auto n_runs{static_cast<unsigned int>(qe.solutions())};

    // Note: the default path length might be smaller than either
    // solution
    switch (qe.solutions()) {
      case 2:
        paths[1] = qe.larger();
        // If there are two solutions, reuse the case for a single
        // solution to setup the intersection with the smaller path
        // in ret[0]
        [[fallthrough]];
      case 1: {
        paths[0] = qe.smaller();
        break;
      }
      default: {
        n_runs = 2u;
        paths[0] = r;
        paths[1] = -r;
      }
    }

    /// Evaluate the function and its derivative at the point @param x
    auto cyl_inters_func = [&h, &r, &sz, &sc](const scalar_t x) {
      const vector3_t crp = vector::cross(h.pos(x) - sc, sz);

      // f(s) = ((h.pos(s) - sc) x sz)^2 - r^2 == 0
      const scalar_t f_s{(vector::dot(crp, crp) - r * r)};
      // f'(s) = 2 * ( (h.pos(s) - sc) x sz) * (h.dir(s) x sz) )
      const scalar_t df_s{2.f * vector::dot(crp, vector::cross(h.dir(x), sz))};

      return std::make_tuple(f_s, df_s);
    };

    for (unsigned int i = 0u; i < n_runs; ++i) {
      const scalar_t &s_ini = paths[i];

      // Run the root finding algorithm
      scalar_t s;
      scalar_t ds;
      if (run_rtsafe) {
        std::tie(s, ds) =
            newton_raphson_safe(cyl_inters_func, s_ini, convergence_tolerance,
                                max_n_tries, max_path);
      } else {
        std::tie(s, ds) =
            newton_raphson(cyl_inters_func, s_ini, convergence_tolerance,
                           max_n_tries, max_path);
      }

      ret[i] = {s, h.pos(s), ds};
    }

    return ret;
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

template <concepts::algebra algebra_t>
struct helix_intersector_impl<concentric_cylindrical2D<algebra_t>, algebra_t>
    : public helix_intersector_impl<cylindrical2D<algebra_t>, algebra_t> {
  using algebra_type = algebra_t;

  template <typename surface_descr_t>
  using intersection_type =
      intersection2D<surface_descr_t, algebra_t, intersection::contains_pos>;

  template <typename other_algebra_t>
  using trajectory_type = detail::helix<other_algebra_t>;

  // Maximum number of solutions this intersector can produce
  static constexpr std::uint8_t n_solutions{1u};

  using result_type = intersection_point_err<algebra_t>;

  /// Operator function to find intersections between helix and a
  /// concentric cylinder surface
  ///
  /// @param h is the input helix trajectory
  /// @param trf is the surface placement transform
  /// @param mask is the input mask that defines the surface extent
  ///
  /// @return the intersection
  template <typename mask_t>
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const trajectory_type<algebra_t> &h, const dtransform3D<algebra_t> &trf,
      const mask_t &mask, const dscalar<algebra_t> = 0.f) const {
    using base_t = helix_intersector_impl<cylindrical2D<algebra_t>, algebra_t>;

    // Array of two solutions
    typename base_t::result_type results =
        base_t::point_of_intersection(h, trf, mask, 0.f);

    // For portals, only take the solution along the helix direction,
    // not the one in the opposite direction
    return math::signbit(results[0].path) ? results[1] : results[0];
  }
};

}  // namespace detray
