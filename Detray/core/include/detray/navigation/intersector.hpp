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
#include "detray/geometry/concepts.hpp"
#include "detray/navigation/intersection/helix_intersector.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"

namespace detray {

/// @brief Intersection interface for detector surfaces.
///
/// Composes the different intersector options into a unifyed interface
template <typename shape_t, concepts::algebra algebra_t,
          bool resolve_pos = false>
struct intersector {
  /// Algebra policy used by the intersection calculations.
  using algebra_type = algebra_t;
  /// Scalar type associated with the algebra policy.
  using scalar_type = dscalar<algebra_t>;
  /// Transform type used to place surfaces in 3D.
  using transform3_type = dtransform3D<algebra_t>;

  /// How to intersect surfaces with rays
  using ray_intersector_type = ray_intersector<shape_t, algebra_t, resolve_pos>;

  /// How to intersect surfaces with helices
  using helix_intersector_type = helix_intersector<shape_t, algebra_t>;

  // Test with int as dummy surface descriptor type
  static_assert(
      std::same_as<
          typename ray_intersector_type::template intersection_type<int>,
          typename helix_intersector_type::template intersection_type<int>>);

  // Maximum number of solutions this intersector can produce
  /// Largest number of solutions produced by either trajectory intersector.
  static constexpr std::uint8_t n_solutions{math::max(
      ray_intersector_type::n_solutions, helix_intersector_type::n_solutions)};

  // Take the helix intersector result type, as it needs more data
  /// Common result representation, which accommodates helix error estimates.
  using result_type = typename helix_intersector_type::result_type;

  /// Operator function to find intersections between ray and planar mask
  ///
  /// @param ray Straight trajectory to test against the surface.
  /// @param trf Transform placing the surface in global coordinates.
  /// @param overstep_tol Minimum accepted signed path length.
  ///
  /// @return the intersection result
  /// @{
  template <typename S = shape_t>
    requires(!concepts::cylindrical_shape<S, algebra_type>)
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const detail::ray<algebra_t> &ray, const transform3_type &trf,
      const scalar_type overstep_tol = 0.f) const {
    return to_result_type(
        ray_intersector_type{}.point_of_intersection(ray, trf, overstep_tol));
  }

  /// Intersect a ray with a cylindrical surface using its radius mask.
  /// @param ray Straight trajectory to test.
  /// @param trf Transform placing the surface in global coordinates.
  /// @param mask Mask providing the cylinder radius.
  /// @param overstep_tol Minimum accepted signed path length.
  /// @returns The geometric intersection result.
  template <typename mask_t, typename S = shape_t>
    requires concepts::cylindrical_shape<S, algebra_type>
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const detail::ray<algebra_t> &ray, const transform3_type &trf,
      const mask_t &mask, const scalar_type overstep_tol = 0.f) const {
    return to_result_type(ray_intersector_type{}.point_of_intersection(
        ray, trf, mask, overstep_tol));
  }

  /// Intersect a helix with a planar or line surface.
  /// @param h Helical trajectory to test.
  /// @param trf Transform placing the surface in global coordinates.
  /// @returns The geometric intersection result.
  template <typename S = shape_t>
    requires(!concepts::cylindrical_shape<S, algebra_type>)
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const detail::helix<algebra_t> &h, const transform3_type &trf,
      const scalar_type /*overstep_tol*/ = 0.f) const {
    return helix_intersector_type{}.point_of_intersection(h, trf);
  }

  /// Intersect a helix with a cylindrical surface using its radius mask.
  /// @param h Helical trajectory to test.
  /// @param trf Transform placing the surface in global coordinates.
  /// @param mask Mask providing the cylinder radius.
  /// @returns The geometric intersection result.
  template <typename mask_t, typename S = shape_t>
    requires concepts::cylindrical_shape<S, algebra_type>
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const detail::helix<algebra_t> &h, const transform3_type &trf,
      const mask_t &mask, const scalar_type /*overstep_tol*/ = 0.f) const {
    return helix_intersector_type{}.point_of_intersection(h, trf, mask);
  }
  /// @}

  /// Resolve a ray intersection against the surface mask.
  /// @param ray Straight trajectory to test.
  /// @param sf Surface descriptor associated with the mask.
  /// @param mask Mask defining the surface extent.
  /// @param trf Transform placing the surface in global coordinates.
  /// @param cfg Tolerances and path cutoff for resolving the result.
  /// @param external_mask_tol Additional tolerance applied at the mask boundary.
  /// @returns The resolved intersection result or results.
  template <typename surface_descr_t, typename mask_t>
  DETRAY_HOST_DEVICE inline decltype(auto) operator()(
      const detail::ray<algebra_t> &ray, const surface_descr_t &sf,
      const mask_t &mask, const transform3_type &trf,
      const intersection::config &cfg = {},
      const scalar_type external_mask_tol = 0.f) const {
    return ray_intersector_type{}(ray, sf, mask, trf, cfg, external_mask_tol);
  }

  /// Resolve a helix intersection against the surface mask.
  /// @param h Helical trajectory to test.
  /// @param sf Surface descriptor associated with the mask.
  /// @param mask Mask defining the surface extent.
  /// @param trf Transform placing the surface in global coordinates.
  /// @param cfg Tolerances and path cutoff for resolving the result.
  /// @returns The resolved intersection result or results.
  template <typename surface_descr_t, typename mask_t>
  DETRAY_HOST_DEVICE inline decltype(auto) operator()(
      const detail::helix<algebra_t> &h, const surface_descr_t &sf,
      const mask_t &mask, const transform3_type &trf,
      const intersection::config &cfg = {},
      const scalar_type /*unused*/ = 0.f) const {
    return helix_intersector_type{}(h, sf, mask, trf, cfg);
  }

 private:
  /// Translate ray intersector results to the common intersector results
  DETRAY_HOST_DEVICE
  constexpr result_type to_result_type(
      const ray_intersector_type::result_type &ray_results) const {
    result_type result;
    if constexpr (n_solutions > 1) {
      for (std::size_t i = 0u; i < n_solutions; ++i) {
        result[i] = typename result_type::value_type{ray_results[i]};
      }
    } else {
      result = result_type{ray_results};
    }
    return result;
  }
};

}  // namespace detray
