// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/concepts.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <string>
#include <vector>

namespace detray::svgtools::meta::proto {

/// @brief A proto landmark class as a simple translation layer from a
/// description of a point.
template <concepts::point3D point3_t>
struct trajectory {
  /// Ordered global positions forming the trajectory polyline.
  std::vector<point3_t> _points;
  /// Identifier used for the trajectory display object.
  std::string _name;
  /// Stroke style used to draw the trajectory.
  actsvg::style::stroke _stroke = actsvg::style::stroke();
};
}  // namespace detray::svgtools::meta::proto
