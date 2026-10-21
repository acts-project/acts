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
#include "detray/utils/ranges.hpp"

// Plugin include(s)
#include "detray/plugins/svgtools/conversion/grid.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/meta.hpp"
#include "actsvg/proto/volume.hpp"

namespace detray::svgtools::meta::display {

/// @brief Converts a proto volume to its surface grid SVG.
///
/// @note the @param grid_svg has to be a non-const reference in order to select
/// the correct overload of the @c actsvg::connectors::connect_action function
template <concepts::point3D point3_t, typename view_t>
inline auto surface_grid(
    const std::string& id,
    const actsvg::proto::volume<std::vector<point3_t>>& p_volume,
    const detray::svgtools::conversion::detail::grid_type& gr_type,
    actsvg::svg::object& grid_svg, view_t view) {
  // Only display grid, if the correct view is requested
  if constexpr (std::is_same_v<view_t, actsvg::views::x_y>) {
    if (gr_type != detray::svgtools::conversion::detail::grid_type::e_endcap) {
      return actsvg::svg::object{};
    }
  } else if constexpr (std::is_same_v<view_t, actsvg::views::z_phi> ||
                       std::is_same_v<view_t, actsvg::views::z_rphi>) {
    if (gr_type != detray::svgtools::conversion::detail::grid_type::e_barrel) {
      return actsvg::svg::object{};
    }
  } else {
    return actsvg::svg::object{};
  }

  if (p_volume._surfaces.empty()) {
    DETRAY_DEBUG_HOST(
        "-> No sensitive surfaces present - skipping surface grid");
    return actsvg::svg::object{};
  }

  if (p_volume._grid_associations.empty()) {
    DETRAY_DEBUG_HOST(
        "-> No bin associations found between surfaces and grid - skipping "
        "surface grid");
    return actsvg::svg::object{};
  }

  // Draw the surface grid with bin highlighting
  const auto& p_surfaces = p_volume._surfaces.front();
  const auto& bin_associations = p_volume._grid_associations.front();

  // Draw the volume sensitive surfaces
  std::vector<actsvg::svg::object> sf_svgs;
  for (const auto& p_surface : p_surfaces) {
    std::string sf_id{p_surface._name + "_" + view._axis_names.at(0) +
                      view._axis_names.at(1)};
    if constexpr (std::is_same_v<view_t, actsvg::views::z_rphi>) {
      view._fixed_r = p_surface._radii.at(0u);
    }
    sf_svgs.push_back(
        actsvg::display::surface(std::move(sf_id), p_surface, view));
  }

  DETRAY_VERBOSE_HOST("-> Drawing surface grid SVG...");

  DETRAY_DEBUG_HOST("-> Found bin-surface associations for "
                    << bin_associations.size() << " bins:");
  assert(!bin_associations.empty());

  for ([[maybe_unused]] const auto& [gbin, assoc] :
       detray::views::enumerate(bin_associations)) {
    DETRAY_DEBUG_HOST("--> Glob bin: " << gbin);

    for ([[maybe_unused]] const std::size_t sf_idx : assoc) {
      DETRAY_DEBUG_HOST("   - surface " << sf_idx);
    }
  }

  actsvg::svg::object sf_grid_svg;
  sf_grid_svg._id = p_volume._name + "_surface_grid";
  sf_grid_svg._tag = "g";

  // Create the surface highlighting
  actsvg::connectors::connect_action(grid_svg._sub_objects, sf_svgs,
                                     bin_associations);

  DETRAY_VERBOSE_HOST("--> Connected surfaces to grid bins");

  sf_grid_svg.add_objects(sf_svgs);

  DETRAY_VERBOSE_HOST("--> Added " << sf_svgs.size()
                                   << " surfaces to grid SVG");

  // Create the association info boxes
  auto x_max = sf_grid_svg._x_range.at(1u);
  for (auto [ig, g_tile] : detray::views::enumerate(grid_svg._sub_objects)) {
    // Target surface text
    std::vector<std::string> bin_text;
    bin_text.push_back("Bin " + std::to_string(ig) + ":");

    for (const auto [is, sis] :
         detray::views::enumerate(bin_associations.at(ig))) {
      const auto& p_surface = p_surfaces.at(sis);
      std::string info{"- surface: " + std::to_string(sis)};
      if (p_surface._aux_info.contains("center")) {
        for (const auto& ci : p_surface._aux_info.at("center")) {
          info += ci;
        }
      }
      bin_text.push_back(info);
    }
    // Make the connected text
    std::string ctext_id{id + "_ct_" + std::to_string(ig)};
    auto ctext = actsvg::draw::connected_text(
        ctext_id, {static_cast<actsvg::scalar>(1.1 * x_max), 0}, bin_text,
        actsvg::style::font{}, actsvg::style::transform{}, g_tile);
    sf_grid_svg.add_object(ctext);
  }

  // Add the grid
  sf_grid_svg.add_object(grid_svg);

  DETRAY_DEBUG_HOST("-> Surface grid SVG ID: " << sf_grid_svg._id);

  return sf_grid_svg;
}

}  // namespace detray::svgtools::meta::display
