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
/// @param logLevel is the screen logging level
///
/// @note no checking for consistency is done at this stage
///
/// @returns the proto container surfaces of a Proto container
ProtoContainer connectVolumesInR(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Connect in Z when having fully prepared input
///
/// @param gctx The geometry context
/// @param volumes the volumes
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note no checking for consistency is done at this stage
ProtoContainer connectVolumesInZ(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Connect in phi when having fully prepared input
///
/// @param gctx The geometry context
/// @param volumes the volumes
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note no checking for consistency is done at this stage
ProtoContainer connectVolumesInPhi(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Connect containers in R when having fully prepared input
///
/// @param gctx The geometry context
/// @param containers the containers
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note no checking for consistency is done at this stage
///
/// @returns the proto container surfaces of a Proto container
ProtoContainer connectContainersInR(
    const GeometryContext& gctx,
    const std::vector<ProtoContainer>& containers,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);


/// @brief Connect containers in Z when having fully prepared input
///
/// @param gctx The geometry context
/// @param containers the containers
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note no checking for consistency is done at this stage
///
/// @returns the proto container surfaces of a Proto container
ProtoContainer connectContainersInZ(
    const GeometryContext& gctx,
    const std::vector<ProtoContainer>& containers,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Helper method to extract r,z,phi boundaries for 
/// eventual grid volume search
///
/// @param gctx the geometry context of the call 
/// @param volumes the volumes at input
/// @param logLevel is the screen logging level
///
/// @return extracted boundary values
std::array<std::vector<ActsScalar>, 3u> rzphiBoundaries(
    const GeometryContext& gctx,
    const std::vector<const DetectorVolume*>& volumes,
    Acts::Logging::Level logLevel = Acts::Logging::INFO);


/// @brief  Create a grid finder for this detector
///
/// @param gctx the geometry context of the call 
/// @param volumes the volumes at input
/// @param logLevel is the screen logging level
/// 
void attachGridVolumeFinder(const GeometryContext& gctx, std::shared_ptr<Detector>& detector,
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

}  // namespace Experimental
}  // namespace Acts