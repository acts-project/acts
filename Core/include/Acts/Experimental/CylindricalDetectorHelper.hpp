// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <limits>
#include <memory>
#include <tuple>
#include <vector>

namespace Acts {

namespace Experimental {

class DetectorVolume;
class Portal;

using ProtoContainer = std::map<unsigned int, std::shared_ptr<Portal>>;

/// @brief Connect in R when having fully prepared input
///
/// @param gctx The geometry context
/// @param volumes the volumes
/// @param selectedOnly switch only selected boundaries
///
/// @note no checking for consistency is done at this stage
///
/// @returns the proto container surfaces of a Proto container
ProtoContainer connectVolumesInR(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly = {});

/// @brief Connect in Z when having fully prepared input
///
/// @param gctx The geometry context
/// @param volumes the volumes
/// @param selectedOnly switch only selected boundaries
///
/// @note no checking for consistency is done at this stage
ProtoContainer connectVolumesInZ(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly = {});

/// @brief Connect in phi when having fully prepared input
///
/// @param gctx The geometry context
/// @param volumes the volumes
/// @param selectedOnly switch only selected boundaries
///
/// @note no checking for consistency is done at this stage
ProtoContainer connectVolumesInPhi(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly = {});

/// @brief Connect containers in R when having fully prepared input
///
/// @param gctx The geometry context
/// @param containers the containers
/// @param selectedOnly switch only selected boundaries
///
/// @note no checking for consistency is done at this stage
///
/// @returns the proto container surfaces of a Proto container
ProtoContainer connectContainersInR(
    const GeometryContext& gctx,
    const std::vector<ProtoContainer>& containers,
    const std::vector<unsigned int>& selectedOnly = {});


/// @brief Connect containers in Z when having fully prepared input
///
/// @param gctx The geometry context
/// @param containers the containers
/// @param selectedOnly switch only selected boundaries
///
/// @note no checking for consistency is done at this stage
///
/// @returns the proto container surfaces of a Proto container
ProtoContainer connectContainersInZ(
    const GeometryContext& gctx,
    const std::vector<ProtoContainer>& containers,
    const std::vector<unsigned int>& selectedOnly = {});

}  // namespace Experimental
}  // namespace Acts