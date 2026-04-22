// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/detail/shape_utils.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <cstdint>
#include <limits>
#include <ostream>

namespace detray {

/// Base class for an intersector: Provides mask based call interfaces for
/// convenience
template <typename intersector_t>
struct intersector_base : public intersector_t {
  /// linear algebra types
  /// @{
  using algebra_type = typename intersector_t::algebra_type;
  using value_type = dvalue<algebra_type>;
  using scalar_type = dscalar<algebra_type>;
  using point2_type = dpoint2D<algebra_type>;
  using point3_type = dpoint3D<algebra_type>;
  using vector3_type = dvector3D<algebra_type>;
  using transform3_type = dtransform3D<algebra_type>;
  /// @}

  template <typename surface_descr_t>
  using intersection_type =
      typename intersector_t::template intersection_type<surface_descr_t>;

  template <typename other_algebra_t>
  using trajectory_type =
      typename intersector_t::template trajectory_type<other_algebra_t>;

  // Maximum number of solutions this intersector can produce
  static constexpr std::uint8_t n_solutions{intersector_t::n_solutions};

  /// Always includes the intersection position, in order to resolve the mask
  using result_type = typename intersector_t::result_type;

  /// Operator function to find intersections between trajectory and a surface
  ///
  /// @tparam surface_descr_t is the type of surface handle
  /// @tparam mask_t is the input mask type
  ///
  /// @param ray is the input ray trajectory
  /// @param sf the surface handle the mask is associated with
  /// @param mask is the input mask that defines the surface extent
  /// @param trf is the surface placement transform
  /// @param mask_tolerance is the tolerance for mask edges
  /// @param overstep_tol negative cutoff for the path
  ///
  /// @returns the intersection
  /// @{

  /// Shapes that produce one solution
  template <typename surface_descr_t, typename mask_t,
            concepts::algebra other_algebra_t, std::uint8_t N = n_solutions>
    requires(N == 1u)
  DETRAY_HOST_DEVICE constexpr intersection_type<surface_descr_t> operator()(
      const trajectory_type<other_algebra_t> &traj, const surface_descr_t &sf,
      const mask_t &mask, const transform3_type &trf,
      const intersection::config &cfg = {},
      const scalar_type external_mask_tol = 0.f) const {
    result_type result =
        call_intersector(traj, mask, trf, cfg.overstep_tolerance);

    intersection_type<surface_descr_t> is;
    resolve_mask(is, traj, result, sf, mask, trf, cfg, external_mask_tol);

    return is;
  }

  /// Shapes that produce multiple solutions
  template <typename surface_descr_t, typename mask_t,
            concepts::algebra other_algebra_t, std::uint8_t N = n_solutions>
    requires(N > 1u)
  DETRAY_HOST_DEVICE constexpr darray<intersection_type<surface_descr_t>,
                                      n_solutions>
  operator()(const trajectory_type<other_algebra_t> &traj,
             const surface_descr_t &sf, const mask_t &mask,
             const transform3_type &trf, const intersection::config &cfg = {},
             const scalar_type external_mask_tol = 0.f) const {
    // One or both of these solutions might be invalid
    result_type result =
        call_intersector(traj, mask, trf, cfg.overstep_tolerance);

    darray<intersection_type<surface_descr_t>, n_solutions> ret;

    for (std::size_t i = 0u; i < n_solutions; ++i) {
      if (detray::detail::any_of(result[i].is_valid())) {
        resolve_mask(ret[i], traj, result[i], sf, mask, trf, cfg,
                     external_mask_tol);
      }
    }

    // Even if there are two geometrically valid solutions, the smaller one
    // might not be passed on if it is below the overstepping tolerance
    return ret;
  }
  /// @}

  /// Interface that uses fixed mask tolerance
  template <typename surface_descr_t, typename mask_t,
            concepts::algebra other_algebra_t>
  DETRAY_HOST_DEVICE constexpr decltype(auto) operator()(
      const trajectory_type<other_algebra_t> &traj, const surface_descr_t &sf,
      const mask_t &mask, const transform3_type &trf,
      const scalar_type mask_tolerance,
      const scalar_type overstep_tol = 0.f) const {
    const intersection::config intr_cfg{
        .min_mask_tolerance = static_cast<float>(mask_tolerance),
        .max_mask_tolerance = static_cast<float>(mask_tolerance),
        .mask_tolerance_scalor = 0.f,
        .overstep_tolerance = static_cast<float>(overstep_tol)};

    return this->operator()(traj, sf, mask, trf, intr_cfg, 0.f);
  }

 private:
  /// Call the underlying intersector
  template <concepts::algebra other_algebra_t, typename mask_t>
  constexpr result_type call_intersector(
      const trajectory_type<other_algebra_t> &traj,
      [[maybe_unused]] const mask_t &mask, const transform3_type &trf,
      const scalar_type overstep_tol = 0.f) const {
    if constexpr (concepts::cylindrical<mask_t>) {
      return intersector_t{}.point_of_intersection(traj, trf, mask,
                                                   overstep_tol);
    } else {
      return intersector_t{}.point_of_intersection(traj, trf, overstep_tol);
    }
  }
};

/// @brief Fill an intersection with the result of the intersection alg.
///
/// @param [out] sfi the surface intersection
/// @param [in] traj the test trajectory that intersects the surface
/// @param [in] s path length to the intersection point
/// @param [in] mask the mask of the surface
/// @param [in] trf the transform of the surface
/// @param [in] mask_tol_scalor scale factor for adaptive mask tolerance
/// @param [in] mask_tolerance minimal and maximal mask tolerance
template <typename intersection_t, typename trajectory_t,
          concepts::algebra algebra_t, concepts::point point_t,
          typename surface_descr_t, typename mask_t,
          concepts::transform3D transform3_t, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr void resolve_mask(
    intersection_t &is, const trajectory_t &traj,
    const intersection_point<algebra_t, point_t, intersection::contains_pos>
        &ip,
    const surface_descr_t sf_desc, const mask_t &mask, const transform3_t &trf,
    const intersection::config &cfg = {},
    const scalar_t external_mask_tolerance = 0.f) {
  // Mask out solutions that don't meet the overstepping tolerance (SoA)
  if constexpr (concepts::soa<algebra_t>) {
    using status_t = typename intersection_t::status_type;

    is.status()(is.path() < cfg.overstep_tolerance) =
        static_cast<status_t>(intersection::status::e_outside);
  } else {
    is.set_status(intersection::status::e_outside);
  }

  // Build intersection struct from test trajectory, if the distance is valid
  if (detray::detail::none_of(ip.path >= cfg.overstep_tolerance)) {
    // Not a valid intersection
    return;
  }

  // Save local position for debug evaluation
  if constexpr (intersection_t::contains_pos()) {
    // Global position on the surface
    dpoint3D<algebra_t> glob_pos;
    if constexpr (concepts::soa<algebra_t>) {
      // The trajectory is given in AoS layout. Translate...
      const auto &origin = traj.pos();
      const auto &dir = traj.dir();

      // Broadcast
      const dvector3D<algebra_t> ro{origin[0], origin[1], origin[2]};
      const dvector3D<algebra_t> rd{dir[0], dir[1], dir[2]};

      glob_pos = ro + ip.path * rd;
    } else {
      // Works for any parameterized trajectory
      glob_pos = traj.pos(ip.path);
    }
    is.set_local(mask_t::to_local_frame3D(trf, glob_pos, traj.dir(ip.path)));
  }

  scalar_t base_tol = 0.f;
  scalar_t ext_tol = 0.f;

  // Tol.: scale with distance of surface to account for track bending
  if (!sf_desc.is_portal() ||
      cfg.min_mask_tolerance == std::numeric_limits<float>::max()) {
    ext_tol = external_mask_tolerance;
    base_tol =
        math::max(static_cast<scalar_t>(cfg.min_mask_tolerance),
                  math::min(static_cast<scalar_t>(cfg.max_mask_tolerance),
                            static_cast<scalar_t>(cfg.mask_tolerance_scalor) *
                                math::fabs(ip.path)));
  }

  // Mask check results with and without external tolerance
  typename mask_t::result_type mask_check{};

  // Intersector provides specialized local point
  if constexpr (std::same_as<point_t, dpoint2D<algebra_t>>) {
    mask_check = mask.resolve(ip.point, base_tol, ext_tol);
  } else {
    // Otherwise, let the shape transform the point to local
    mask_check = mask.resolve(trf, ip.point, base_tol, ext_tol);
  }

  // Set the less strict status first, then overwrite with more strict
  is.set_status_if(intersection::status::e_edge,
                   detray::get<check_type::e_with_edge>(mask_check));
  is.set_status_if(intersection::status::e_inside,
                   detray::get<check_type::e_precise>(mask_check));

  is.set_path(ip.path);
  is.set_surface(sf_desc);
  is.set_direction(!detail::signbit(ip.path));
  is.set_volume_link(mask.volume_link());
}

}  // namespace detray
