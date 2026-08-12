// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <limits>

namespace Acts {

/// Tolerance for being numerical equal for geometry building
static constexpr double s_epsilon = 3 * std::numeric_limits<double>::epsilon();

/// Tolerance for being on Surface
///
/// @note This is intentionally given w/o an explicit unit to avoid having
///       to include the units header unnecessarily. With the native length
///       unit of mm this corresponds to 0.1um.
static constexpr double s_onSurfaceTolerance = 1e-4;

/// Tolerance in radians for a phi sector to count as the full azimuth
///
/// @note Geometry descriptions state phi edges with about ten significant
///       digits, so a full sector can miss its nominal value by ~1e-9. Even
///       at a metre radius that is a nanometre of arc, while the narrowest
///       real sectors are milliradians wide.
static constexpr double s_fullAzimuthTolerance = 1e-9;

/// Tolerance for not being within curvilinear projection
/// this allows using the same curvilinear frame to eta = 6,
/// validity tested with IntegrationTests/PropagationTest
static constexpr double s_curvilinearProjTolerance = 0.999995;

}  // namespace Acts
