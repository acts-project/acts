// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/utils/invalid_values.hpp"

// Plugin include(s)
#include "detray/plugins/svgtools/utils/surface_kernels.hpp"

// System include(s)
#include <cassert>
#include <tuple>

namespace detray::svgtools::utils {

/// @brief Checks if the detray surface has a volume link.
template <typename detector_t>
inline auto is_not_world_portal(
    const detray::geometry::surface<detector_t>& d_portal) {
  const auto d_link_idx = d_portal.volume_links();
  bool is_world_pt{false};
  for (const auto vol_link : d_link_idx) {
    is_world_pt = is_world_pt || detray::detail::is_invalid_value(vol_link);
  }
  return !is_world_pt;
}

/// @note expects that the detray surface has a volume link.
/// @returns the volume link of the detray surface.
template <typename detector_t>
inline auto get_linked_volume(
    const detector_t& detector,
    const detray::geometry::surface<detector_t>& d_portal,
    [[maybe_unused]] const std::size_t mask_idx) {
  assert(is_not_world_portal(d_portal));
  const auto d_link_idx = d_portal.volume_links();

  return tracking_volume{detector, d_link_idx.at(mask_idx)};
}

/// @brief Calculates the start and end point of the link.
/// @note The detray surface must have a volume link.
/// @returns (start, end).
template <typename detector_t>
inline auto link_points(const typename detector_t::geometry_context& context,
                        const detector_t& detector,
                        const detray::geometry::surface<detector_t>& d_portal,
                        typename detector_t::vector3_type dir,
                        const double link_length,
                        const std::size_t mask_idx = 0u) {
  assert(is_not_world_portal(d_portal));

  // Calculating the start position:
  const auto start = d_portal.template visit_mask<link_start_getter>(
      d_portal.transform(context), mask_idx);

  // Calculating the end position:
  const auto n =
      d_portal.normal(context, d_portal.global_to_local(context, start, dir));
  const auto volume = get_linked_volume(detector, d_portal, mask_idx);
  const auto end = d_portal.template visit_mask<link_end_getter>(
      detector, volume, start, n, link_length, mask_idx);

  return std::make_tuple(start, end);
}

}  // namespace detray::svgtools::utils
