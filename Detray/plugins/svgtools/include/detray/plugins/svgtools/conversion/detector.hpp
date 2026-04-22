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
#include "detray/geometry/tracking_volume.hpp"
#include "detray/plugins/svgtools/conversion/volume.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"

// Actsvg include(s)
#include "actsvg/proto/detector.hpp"

namespace detray::svgtools::conversion {

/// @brief Generates the proto detector object
///
/// @param context The geometry context.
/// @param detector The detector object
/// @param style the style settings.
/// @param hide_portals whether to display portals.
/// @param hide_passives whether to display passive surfaces.
/// @param hide_grids whether to display the volume surface grids.
///
/// @returns An actsvg proto detector representing
template <typename detector_t, typename view_t>
auto detector(const typename detector_t::geometry_context& context,
              const detector_t& detector, const view_t& view,
              const styling::detector_style& style =
                  styling::tableau_colorblind::detector_style,
              bool hide_portals = false, bool hide_passives = false,
              bool hide_grids = false) {
  using point3_container_t = std::vector<typename detector_t::point3_type>;
  actsvg::proto::detector<point3_container_t> p_detector;

  for (const auto& vol_desc : detector.volumes()) {
    auto [p_volume, gr_type] = svgtools::conversion::volume(
        context, detector, tracking_volume{detector, vol_desc}, view,
        style._volume_style, hide_portals, hide_passives, hide_grids);

    p_detector._volumes.push_back(p_volume);
  }

  return p_detector;
}

}  // namespace detray::svgtools::conversion
