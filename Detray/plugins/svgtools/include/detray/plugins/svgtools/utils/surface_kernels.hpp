// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/detail/radius_getter.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/utils/concepts.hpp"

// System include(s)
#include <cassert>
#include <optional>

namespace detray::svgtools::utils {

/// @brief Functor to calculate a suitable starting point for displaying the
/// link arrow.
struct link_start_getter {
 public:
  template <typename mask_group_t, concepts::index index_t,
            concepts::transform3D transform3_t>
  DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                     const index_t& index,
                                     const transform3_t& transform,
                                     const std::size_t /*unused*/ = 0) const {
    return link_start(mask_group.at(index), transform);
  }

  template <typename mask_group_t, concepts::interval idx_range_t,
            concepts::transform3D transform3_t>
  DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                     const idx_range_t& idx_range,
                                     const transform3_t& transform,
                                     const std::size_t mask_idx) const {
    assert(mask_idx < idx_range.size());
    return link_start(mask_group.at(idx_range.lower() + mask_idx), transform);
  }

 private:
  // Calculates the link starting location of the remaining shapes.
  template <typename mask_t, concepts::transform3D transform3_t>
  typename mask_t::point3_type inline link_start(
      const mask_t& mask, const transform3_t& transform) const {
    return transform.point_to_global(mask.centroid());
  }

  // Calculates the (optimal) link starting point for rings.
  template <concepts::transform3D transform3_t, concepts::algebra algebra_t>
  dpoint3D<algebra_t> inline link_start(
      const detray::mask<detray::ring2D, algebra_t>& mask,
      const transform3_t& transform) const {
    using shape_t = detray::ring2D;
    using mask_t = detray::mask<shape_t, algebra_t>;
    using point3_t = typename mask_t::point3_type;
    using scalar_t = typename mask_t::scalar_type;

    const scalar_t r{0.5f *
                     (mask[shape_t::e_inner_r] + mask[shape_t::e_outer_r])};
    const scalar_t phi{detray::constant<scalar_t>::pi_2};

    return mask_t::to_global_frame(transform, point3_t{r, phi, 0.f});
  }

  // Calculates the (optimal) link starting point for annuluses.
  template <concepts::transform3D transform3_t, concepts::algebra algebra_t>
  dpoint3D<algebra_t> inline link_start(
      const detray::mask<detray::annulus2D, algebra_t>& mask,
      const transform3_t& transform) const {
    using shape_t = detray::annulus2D;
    using mask_t = detray::mask<shape_t, algebra_t>;
    using point3_t = typename mask_t::point3_type;
    using scalar_t = typename mask_t::scalar_type;

    const scalar_t r{(mask[shape_t::e_min_r] + mask[shape_t::e_max_r]) / 2.f};
    const scalar_t phi{mask[shape_t::e_average_phi]};

    return mask_t::to_global_frame(transform, point3_t{r, phi, 0.f});
  }

  // Calculates the (optimal) link starting point for concentric cylinders
  template <concepts::transform3D transform3_t, concepts::algebra algebra_t>
  dpoint3D<algebra_t> inline link_start(
      const detray::mask<concentric_cylinder2D, algebra_t>& mask,
      const transform3_t& transform) const {
    using mask_t = detray::mask<concentric_cylinder2D, algebra_t>;
    using point3_t = typename mask_t::point3_type;
    using scalar_t = typename mask_t::scalar_type;

    const scalar_t r{mask[concentric_cylinder2D::e_r]};
    const scalar_t phi{detray::constant<scalar_t>::pi_2};
    // Shift the center to the actual cylinder bounds
    const scalar_t z{mask.centroid()[2]};

    return mask_t::to_global_frame(transform, point3_t{phi, z, r});
  }

  // Calculates the (optimal) link starting point for cylinders (2D).
  template <concepts::transform3D transform3_t, concepts::algebra algebra_t>
  dpoint3D<algebra_t> inline link_start(
      const detray::mask<cylinder2D, algebra_t>& mask,
      const transform3_t& transform) const {
    using mask_t = detray::mask<cylinder2D, algebra_t>;
    using point3_t = typename mask_t::point3_type;
    using scalar_t = typename mask_t::scalar_type;

    const scalar_t r{mask[cylinder2D::e_r]};
    const scalar_t phi{detray::constant<scalar_t>::pi_2};
    // Shift the center to the actual cylinder bounds
    const scalar_t z{mask.centroid()[2]};

    return mask_t::to_global_frame(transform, point3_t{r * phi, z, r});
  }

  // Calculates the (optimal) link starting point for cylinders (3D).
  template <concepts::transform3D transform3_t, concepts::algebra algebra_t>
  dpoint3D<algebra_t> inline link_start(
      const detray::mask<detray::cylinder3D, algebra_t>& mask,
      const transform3_t& transform) const {
    using shape_t = detray::cylinder3D;
    using mask_t = detray::mask<shape_t, algebra_t>;
    using point3_t = typename mask_t::point3_type;
    using scalar_t = typename mask_t::scalar_type;

    const scalar_t r{mask[shape_t::e_max_r]};
    const scalar_t phi{0.5f *
                       (mask[shape_t::e_max_phi] + mask[shape_t::e_max_phi])};
    const scalar_t z{mask.centroid()[2]};

    return mask_t::to_global_frame(transform, point3_t{r, phi, z});
  }
};

/// @brief Functor to calculate a suitable end point for displaying the link
/// arrow.
struct link_end_getter {
 public:
  template <typename mask_group_t, concepts::index index_t, typename detector_t,
            concepts::point3D point3_t, concepts::vector3D vector3_t,
            concepts::scalar scalar_t>
  DETRAY_HOST inline vector3_t operator()(
      const mask_group_t& mask_group, const index_t& index,
      const detector_t& detector,
      const detray::tracking_volume<detector_t>& volume,
      const point3_t& surface_point, const vector3_t& surface_normal,
      const scalar_t& link_length, const std::size_t /*unused*/ = 0) const {
    return link_dir(mask_group.at(index), detector, volume, surface_point,
                    surface_normal) *
               link_length +
           surface_point;
  }

  template <typename mask_group_t, concepts::interval idx_range_t,
            typename detector_t, concepts::point3D point3_t,
            concepts::vector3D vector3_t, concepts::scalar scalar_t>
  DETRAY_HOST inline vector3_t operator()(
      const mask_group_t& mask_group, const idx_range_t& idx_range,
      const detector_t& detector,
      const detray::tracking_volume<detector_t>& volume,
      const point3_t& surface_point, const vector3_t& surface_normal,
      const scalar_t& link_length, const std::size_t mask_idx) const {
    assert(mask_idx < idx_range.size());
    return link_dir(mask_group.at(idx_range.lower() + mask_idx), detector,
                    volume, surface_point, surface_normal) *
               link_length +
           surface_point;
  }

 private:
  /// @brief Calculates the direction of the link for remaining shapes.
  template <typename detector_t, typename mask_t, concepts::point3D point3_t,
            concepts::vector3D vector3_t>
  inline vector3_t link_dir(const mask_t& /*mask*/,
                            const detector_t& /*detector*/,
                            const detray::tracking_volume<detector_t>& volume,
                            const point3_t& surface_point,
                            const vector3_t& surface_normal) const {
    const auto dir = volume.center() - surface_point;
    const auto dot_prod = vector::dot(dir, surface_normal);

    // Should geometrically not happen with a local point 'surface_point'
    assert(dot_prod != 0.f);

    return math::copysign(1.f, dot_prod) * surface_normal;
  }

  /// @brief Calculates the direction of the link for cylinders (2D)
  template <typename detector_t, concepts::point3D point3_t,
            concepts::vector3D vector3_t, typename shape_t>
    requires std::is_same_v<shape_t, cylinder2D> ||
             std::is_same_v<shape_t, concentric_cylinder2D>
  inline vector3_t link_dir(
      const detray::mask<shape_t, typename detector_t::algebra_type>& mask,
      const detector_t& detector,
      const detray::tracking_volume<detector_t>& volume,
      const point3_t& /*surface_point*/,
      const vector3_t& surface_normal) const {
    for (const auto& desc : volume.portals()) {
      const detray::geometry::surface surface{detector, desc};

      if (auto r =
              surface
                  .template visit_mask<detray::detail::outer_radius_getter>()) {
        if (*r > mask[shape_t::e_r]) {
          return surface_normal;
        }
      }
    }
    return vector3_t{-1.f * surface_normal};
  }
};

}  // namespace detray::svgtools::utils
