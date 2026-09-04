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
#include "detray/plugins/svgtools/conversion/landmark.hpp"
#include "detray/plugins/svgtools/meta/proto/intersection.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools::conversion {

/// @returns The proto intersection of a detray intersection.
template <typename detector_t, typename intersection_t>
inline auto intersection(const detector_t& detector,
                         const std::vector<intersection_t>& intersections,
                         const typename detector_t::vector3_type& dir = {},
                         const typename detector_t::geometry_context& gctx = {},
                         const dindex_range highlight_idx =
                             {detray::detail::invalid_value<dindex>(),
                              detray::detail::invalid_value<dindex>()},
                         const styling::landmark_style& style =
                             styling::svg_default::intersection_style) {
  using point3_t = typename detector_t::point3_type;
  using point2_t = typename detector_t::point2_type;
  using p_intersection_t = svgtools::meta::proto::intersection<point3_t>;
  p_intersection_t p_ir;

  for (const auto& intr : intersections) {
    const detray::geometry::surface sf{detector, intr.surface()};
    if (sf.identifier().is_invalid()) {
      continue;
    }

    const point2_t bound{intr.local()[0], intr.local()[1]};
    const auto position = sf.local_to_global(gctx, bound, dir);
    const auto p_lm = svgtools::conversion::landmark(position, style);

    p_ir._landmarks.push_back(p_lm);
  }

  svgtools::styling::apply_style(p_ir, style);

  // Place highlights
  assert(highlight_idx[0] <= highlight_idx[1]);
  if (highlight_idx[1] - highlight_idx[0] != 0u) {
    styling::landmark_style highlight_style{style};

    constexpr actsvg::scalar opacity{1.f};
    highlight_style._fill_color = {styling::colors::red, opacity};

    for (std::size_t i = highlight_idx[0]; i <= highlight_idx[1]; ++i) {
      styling::apply_style(p_ir._landmarks.at(i), highlight_style);
    }
  }

  return p_ir;
}

}  // namespace detray::svgtools::conversion
