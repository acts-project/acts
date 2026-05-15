// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <string>

namespace detray::svgtools::meta::proto {

/// @brief A proto landmark class as a simple translation layer from a
/// description of a point.
template <concepts::point3D point3_t>
struct landmark {
  point3_t _position{0.f, 0.f, 0.f};
  std::string _name{"unknown landmark"};
  actsvg::style::marker _marker{"x", 1.f};
};

}  // namespace detray::svgtools::meta::proto
