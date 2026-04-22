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
#include "detray/geometry/tracking_volume.hpp"

// Plugin include(s)
#include "detray/plugins/svgtools/conversion/portal.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/conversion/surface_grid.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"

// Actsvg include(s)
#include "actsvg/display/geometry.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"

// System include(s)
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace detray::svgtools::conversion {

/// @brief Generates the proto volume of a detray volume.
///
/// @param context The geometry context.
/// @param detector The detractor the volume belongs to.
/// @param d_volume The detray volume.
/// @param style the style settings
/// @param hide_portals whether to display the volumes portals.
/// @param hide_passives whether to display the contained passive surfaces.
/// @param hide_grids whether to display the contained surface grid.
/// @param search_window neighborhood search window for the grid.
///
/// @returns An actsvg proto volume representing the volume.
template <typename detector_t, typename view_t>
auto volume(const typename detector_t::geometry_context& context,
            const detector_t& detector,
            const detray::tracking_volume<detector_t>& d_volume,
            const view_t& view,
            const styling::volume_style& style =
                styling::tableau_colorblind::volume_style,
            bool hide_portals = false, bool hide_passives = false,
            bool hide_grids = false, bool hide_material = true,
            const std::array<dindex, 2>& search_window = {2u, 2u}) {
  using point3_container_t = std::vector<typename detector_t::point3_type>;

  actsvg::proto::volume<point3_container_t> p_volume;
  p_volume._index = d_volume.index();

  // Prepare surfaces to be displayed in grids
  std::vector<actsvg::proto::surface<point3_container_t>> p_sensitves;

  // Convert grid, if present
  auto [p_grid, grid_type] = svgtools::conversion::surface_grid(
      detector, p_volume._index, view, style._grid_style);

  for (const auto& desc : d_volume.surfaces()) {
    const auto sf = detray::geometry::surface<detector_t>{detector, desc};

    if (sf.is_portal()) {
      if (!hide_portals) {
        auto p_portals = svgtools::conversion::portal(context, detector, sf,
                                                      view, style._portal_style,
                                                      false, hide_material);

        std::ranges::copy(p_portals, std::back_inserter(p_volume._portals));
      }
    } else if (!(sf.is_passive() && hide_passives)) {
      const auto& sf_style = sf.is_sensitive() ? style._sensitive_surface_style
                                               : style._passive_surface_style;

      auto p_surfaces = svgtools::conversion::surface(
          context, detector, sf, view, sf_style, hide_material);
      assert(p_surfaces.size() == sf.n_masks());

      for (std::size_t i = 0u; i < p_surfaces.size(); ++i) {
        auto& p_surface = p_volume._v_surfaces.emplace_back(p_surfaces.at(i));

        std::string sf_info{"* index " + std::to_string(sf.index()) +
                            "* mask " + std::to_string(i)};

        p_surface._aux_info["module_info"] = {sf_info};
        p_surface._aux_info["grid_info"] = {sf_info};

        // Put the sensitive surfaces in the grids
        if (sf.is_sensitive()) {
          p_sensitves.push_back(p_surface);
        }
      }
    }
  }

  // Add the proto grid to the proto volume and find bin associations
  if (!hide_grids && p_grid.has_value()) {
    p_volume._surface_grid = *p_grid;
    p_volume._grid_associations = {
        get_bin_association(detector, d_volume, search_window)};
  }

  p_volume._surfaces = {std::move(p_sensitves)};

  return std::tuple(p_volume, grid_type);
}

}  // namespace detray::svgtools::conversion
