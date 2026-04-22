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
#include "detray/plugins/svgtools/meta/proto/information_section.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <string>

namespace detray::svgtools::meta::display {

template <concepts::point3D point3_t, concepts::point2D point2_t,
          typename view_t>
inline auto information_section(
    const std::string& id,
    const svgtools::meta::proto::information_section<point3_t>& is,
    const view_t& view, const point2_t& screen_offset,
    const actsvg::svg::object& connected_object,
    const bool use_relative_offset = false) {
  // Title style
  actsvg::style::fill title_fill;
  title_fill._fc._rgb = is._color;
  title_fill._fc._opacity = 0.8f;
  actsvg::style::font title_font;

  title_font._size = 16u;
  title_font._fc = actsvg::style::color{{255, 255, 255}};

  // Info text style
  actsvg::style::fill info_fill;
  info_fill._fc._rgb = is._color;
  info_fill._fc._opacity = 0.4f;

  actsvg::style::font info_font;
  info_font._size = 12;

  // Box stroke
  actsvg::style::stroke stroke;

  const auto position =
      use_relative_offset
          ? (screen_offset + view(std::vector{is._position}).at(0))
          : screen_offset;

  return actsvg::draw::connected_info_box(
      id, position, is._title, title_fill, title_font, is._info, info_fill,
      info_font, stroke, connected_object, {"mousedown", "mouseout"});
}
}  // namespace detray::svgtools::meta::display
