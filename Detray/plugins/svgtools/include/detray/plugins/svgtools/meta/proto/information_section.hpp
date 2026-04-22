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
#include <vector>

namespace detray::svgtools::meta::proto {

/// @brief A proto intersection record class as a simple translation layer from
/// a intersection record description.
template <concepts::point3D point3_t>
struct information_section {
  std::vector<std::string> _info{};
  std::string _title{};
  std::array<int, 3> _color{};
  point3_t _position{};
};

}  // namespace detray::svgtools::meta::proto
