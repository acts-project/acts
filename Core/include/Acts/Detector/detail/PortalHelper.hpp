// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {

class DetectorVolume;

/// @brief Definition of a portal replacement when building proto containers
///
/// It consists of the new portal, the index, the direction, the parameters
/// gathered from the sub volumes, the binning description
using PortalReplacement =
    std::tuple<std::shared_ptr<Experimental::Portal>, unsigned int, Direction,
               std::vector<ActsScalar>, BinningValue>;

namespace detail {
namespace PortalHelper {

/// @brief Create and attach the multi link updator, the portal will get
/// a volume updator attached, that points to the different sub volumes
/// depending on the global position and binning
///
/// @param gctx the geometry context
/// @param volumes are the volumes that are pointed to
/// @param pReplacements are the portal replacements that are newly connected
///
void attachDetectorVolumeUpdators(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    std::vector<PortalReplacement>& pReplacements);

/// @brief Method that strips out attached volumes from portals and
/// provides them back to the caller.
///
/// @note it throws an exception if both sides are already taken
///
/// @param portal the portal to be resolved
///
/// @return a vector of attached volumes
std::vector<std::shared_ptr<DetectorVolume>> attachedDetectorVolumes(
    Portal& portal) noexcept(false);

}  // namespace PortalHelper
}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
