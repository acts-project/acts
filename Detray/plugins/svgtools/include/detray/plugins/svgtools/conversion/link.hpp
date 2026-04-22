// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"

// Plugin include(s)
#include "detray/plugins/svgtools/utils/link_utils.hpp"
#include "detray/plugins/svgtools/utils/surface_kernels.hpp"

// Actsvg includes(s)
#include "actsvg/proto/portal.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools::conversion {

/// @returns The proto link calculated using the surface normal vector.
template <typename detector_t>
inline auto links(const typename detector_t::geometry_context& context,
                  const detector_t& detector,
                  const detray::geometry::surface<detector_t>& d_portal) {
  using point3_container_t = std::vector<typename detector_t::point3_type>;
  using p_link_t = typename actsvg::proto::portal<point3_container_t>::link;

  std::vector<p_link_t> pt_links{};
  typename detector_t::vector3_type dir{};

  // Length of link arrow is currently hardcoded.
  constexpr double link_length = 4.;

  for (std::size_t i = 0u; i < d_portal.n_masks(); ++i) {
    const auto [start, end] = svgtools::utils::link_points(
        context, detector, d_portal, dir, link_length, i);

    p_link_t p_link;
    p_link._start = start;
    p_link._end = end;
    pt_links.push_back(p_link);
  }

  return pt_links;
}

}  // namespace detray::svgtools::conversion
