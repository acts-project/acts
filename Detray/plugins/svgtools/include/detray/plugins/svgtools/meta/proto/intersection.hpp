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
#include "detray/plugins/svgtools/meta/proto/landmark.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <string>
#include <vector>

namespace detray::svgtools::meta::proto {

/// @brief A proto intersection class as a simple translation layer from
/// an intersection description.
template <concepts::point3D point3_t>
struct intersection {
  using landmark_type = svgtools::meta::proto::landmark<point3_t>;
  std::vector<landmark_type> _landmarks;
  std::string _name;
};

}  // namespace detray::svgtools::meta::proto
