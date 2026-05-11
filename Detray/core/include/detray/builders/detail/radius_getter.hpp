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
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/utils/concepts.hpp"

// System include(s)
#include <cassert>
#include <optional>

namespace detray::detail {

/// @brief Functor to calculate the outermost radius of a cylinder shape.
/// If the shape is not defined by a radius, then null option is returned.
struct outer_radius_getter {
 public:
  template <typename mask_group_t, concepts::index index_t>
  DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                     const index_t& index) const {
    return outer_radius(mask_group.at(index));
  }

  template <typename mask_group_t, concepts::interval idx_range_t>
  DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                     const idx_range_t& idx_range) const {
    // All masks on the same cylinder surface have the same radius
    return outer_radius(mask_group.at(idx_range.lower()));
  }

 private:
  // Struct to access the radius of a surface
  template <typename mask_t>
  DETRAY_HOST std::optional<typename mask_t::scalar_type> inline outer_radius(
      const mask_t& /*mask*/) const {
    return std::nullopt;
  }

  // Calculates the outer radius for cylinders (2D).
  template <concepts::algebra algebra_t>
  DETRAY_HOST inline auto outer_radius(
      const detray::mask<detray::cylinder2D, algebra_t>& mask) const {
    return std::optional(mask[cylinder2D::e_r]);
  }

  // Calculates the outer radius for concentric cylinders (2D).
  template <concepts::algebra algebra_t>
  DETRAY_HOST inline auto outer_radius(
      const detray::mask<detray::concentric_cylinder2D, algebra_t>& mask)
      const {
    return std::optional(mask[concentric_cylinder2D::e_r]);
  }

  // Calculates the outer radius for cylinders (3D).
  template <concepts::algebra algebra_t>
  DETRAY_HOST inline auto outer_radius(
      const detray::mask<detray::cylinder3D, algebra_t>& mask) const {
    return std::optional(mask[cylinder3D::e_max_r]);
  }
};

}  // namespace detray::detail
