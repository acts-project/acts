// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {

class DetectorVolume;
class Portal;

namespace DetrayJsonHelper {

/// Helper method to make suitable portal surfaces for detray
///
/// @param gctx the geometry context
/// @param volume the volume in question
/// @param logLevel is a screen logging level
///
/// @return new portal surfaces
std::vector<std::shared_ptr<Portal>> clipAndPickPortals(
    const Acts::GeometryContext& gctx,
    const Acts::Experimental::DetectorVolume& volume,
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// Clipping method for scalars and volumes
///
/// @param vBoundaries the unclipped volume boundaries
/// @param volumes the unclipped vector of volumes
/// @param clipRange is the new range to which these are to be clipeed to
/// @param logLevel is a screen logging level
///
/// @return the cluipped value vector vector and clipped volume vector
std::tuple<std::vector<ActsScalar>, std::vector<const DetectorVolume*>> clip(
    const std::vector<ActsScalar>& vBoundaries,
    const std::vector<const DetectorVolume*>& volumes,
    const std::array<ActsScalar, 2u>& clipRange,
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

}  // namespace DetrayJsonHelper
}  // namespace Experimental
}  // namespace Acts