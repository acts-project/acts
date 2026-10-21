// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <span>
#include <tuple>

#include <nlohmann/json.hpp>

/// @namespace Acts::DetrayJsonHelper
/// @ingroup json_plugin
/// @todo Move this to the @ref ActsPlugins namespace

namespace Acts::DetrayJsonHelper {

/// @addtogroup json_plugin
/// @{

/// Determine the acceleration link from a grid
///
///
///   brute_force = 0u,      // try all
///   cartesian2_grid = 1u,  // rectangle, trapezoid, (triangle) grids
///   cuboid3_grid = 2u,     // cuboid grid
///   polar2_grid = 3u,      // ring/disc, annulus grids
///   cylinder2_grid = 4u,   // 2D cylinder grid
///   cylinder3_grid = 5u,   // 3D cylinder grid
///
/// @param casts are the grid axes cast types
///
/// @return the acceleration link idnetifier
std::size_t accelerationLink(std::span<const AxisDirection> casts);

/// @}
}  // namespace Acts::DetrayJsonHelper
