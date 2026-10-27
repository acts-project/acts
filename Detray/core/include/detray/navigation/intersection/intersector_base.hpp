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
  /// Scalar value type for the configured algebra layout.
  using value_type = dvalue<algebra_type>;
  /// Scalar calculation type for the configured algebra layout.
  using scalar_type = dscalar<algebra_type>;
  /// Two-dimensional point type for the configured algebra layout.
  using point2_type = dpoint2D<algebra_type>;
  /// Three-dimensional point type for the configured algebra layout.
  using point3_type = dpoint3D<algebra_type>;
  /// Three-dimensional vector type for the configured algebra layout.
  using vector3_type = dvector3D<algebra_type>;
  /// Transform type used to place the surface in 3D.
  using transform3_type = dtransform3D<algebra_type>;
  /// @}

  /// Intersection result type associated with a surface descriptor.
  template <typename surface_descr_t>
  using intersection_type =
      typename intersector_t::template intersection_type<surface_descr_t>;

  /// Trajectory representation accepted for the requested algebra layout.
  template <typename other_algebra_t>
  using trajectory_type =
      typename intersector_t::template trajectory_type<other_algebra_t>;

  /// Maximum number of intersection solutions returned by the intersector.
  static constexpr std::uint8_t n_solutions{intersector_t::n_solutions};

  /// Underlying geometric result retained before mask resolution.
  using result_type = typename intersector_t::result_type;

  /// Operator function to find intersections between trajectory and a surface
  ///
  /// @tparam surface_descr_t is the type of surface handle
  /// @tparam mask_t is the input mask type
  ///
  /// @{

  /// Shapes that produce one solution
  /// @param traj Trajectory to intersect with the surface.
  /// @param sf Surface descriptor associated with the mask.
  /// @param mask Surface mask used to classify the candidate point.
  /// @param trf Transform from the surface frame to the global frame.
  /// @param cfg Intersection tolerances and path cutoff.
  /// @param external_mask_tol Additional tolerance for the mask boundary.
  /// @returns The resolved intersection, including mask status.
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
  /// @param traj Trajectory to intersect with the surface.
  /// @param sf Surface descriptor associated with the mask.
  /// @param mask Surface mask used to classify each candidate point.
  /// @param trf Transform from the surface frame to the global frame.
  /// @param cfg Intersection tolerances and path cutoff.
  /// @param external_mask_tol Additional tolerance for the mask boundary.
  /// @returns The resolved intersection candidates.
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
  /// @param traj Trajectory to intersect with the surface.
  /// @param sf Surface descriptor associated with the mask.
  /// @param mask Mask defining the finite extent of the surface.
  /// @param trf Transform from the surface frame to the global frame.
  /// @param mask_tolerance Fixed tolerance used to classify the mask boundary.
  /// @param overstep_tol Minimum accepted signed path length.
  /// @returns The resolved intersection result or results.
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
/// @param is Intersection object to populate.
/// @param traj Trajectory used to evaluate the global intersection position.
/// @param ip Geometric intersection candidate returned by the intersector.
/// @param sf_desc Surface descriptor associated with the mask.
/// @param mask Mask used to classify the candidate's local position.
/// @param trf Transform from the surface frame to the global frame.
/// @param cfg Tolerances and path cutoff used for classification.
/// @param external_mask_tolerance Additional tolerance for the mask boundary.
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
  is.set_direction(!math::signbit(ip.path));
  is.set_volume_link(mask.volume_link());
}

}  // namespace detray
