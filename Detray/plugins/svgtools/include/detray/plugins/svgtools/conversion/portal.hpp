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

// Plugin include(s)
#include "detray/plugins/svgtools/conversion/link.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"
#include "detray/plugins/svgtools/utils/link_utils.hpp"

// Actsvg includes(s)
#include "actsvg/proto/portal.hpp"

namespace detray::svgtools::conversion {

/// @returns An actsvg proto portal representing the portal.
/// @note detray portal is_portal() should be true.
template <typename detector_t, typename view_t>
auto portal(const typename detector_t::geometry_context& context,
            const detector_t& detector,
            const detray::geometry::surface<detector_t>& d_portal,
            const view_t& view,
            const styling::portal_style& style =
                styling::tableau_colorblind::portal_style,
            bool hide_links = false, bool hide_material = false) {
  assert(d_portal.is_portal());

  using point3_container_t = std::vector<typename detector_t::point3_type>;
  using p_portal_t = actsvg::proto::portal<point3_container_t>;
  using p_surface_t = actsvg::proto::surface<point3_container_t>;

  const auto p_surfaces = svgtools::conversion::surface(
      context, detector, d_portal, view, style._surface_style, hide_material);

  std::vector<typename p_portal_t::link> pt_links{};
  const bool display_links{!hide_links &&
                           svgtools::utils::is_not_world_portal(d_portal)};
  if (display_links) {
    pt_links = svgtools::conversion::links(context, detector, d_portal);
  }

  std::vector<p_portal_t> p_portals{};
  for (std::size_t i = 0u; i < p_surfaces.size(); ++i) {
    const auto& p_sub_portal = p_surfaces.at(i);

    auto& p_portal = p_portals.emplace_back();
    p_portal._name =
        "portal_" + std::to_string(d_portal.index()) + "_" + std::to_string(i);
    p_portal._surface = p_sub_portal;
    p_portal._surface._sf_type = p_surface_t::sf_type::e_portal;

    if (display_links) {
      p_portal._volume_links = {pt_links.at(i)};
    }

    svgtools::styling::apply_style(p_portal, style);
  }

  return p_portals;
}

}  // namespace detray::svgtools::conversion
