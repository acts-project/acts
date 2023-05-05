// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <limits>
#include <map>
#include <memory>
#include <vector>

namespace Acts {

namespace Experimental {

class DetectorVolume;
class Portal;

namespace detail {
namespace CylindricalDetectorHelper {

/// @brief Connect detector volumes in R
///
/// @param gctx The geometry context
/// @param volumes the volumes
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note a fair amount of consistency checking is done,
/// and exceptions are thrown if any of the tests fail
///
/// @return a proto container with the outside portals
DetectorComponent::PortalContainer connectInR(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Connect detector volumes in Z
///
/// @param gctx The geometry context
/// @param volumes the volumes
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note a fair amount of consistency checking is done,
/// and exceptions are thrown if any of the tests fail
///
/// @return a proto container with the outside portals
DetectorComponent::PortalContainer connectInZ(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Connect detector volumes in phi
///
/// @param gctx The geometry context
/// @param volumes the volumes
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note a fair amount of consistency checking is done,
/// and exceptions are thrown if any of the tests fail
///
/// @return a proto container with the outside portals
DetectorComponent::PortalContainer connectInPhi(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Wrap detector volumes in R,Z
///
/// @param gctx The geometry context
/// @param volumes the volumes
/// @param logLevel is the screen logging level
///
/// @note a fair amount of consistency checking is done,
/// and exceptions are thrown if any of the tests fail
///
/// @return a proto container with the outside portals
DetectorComponent::PortalContainer wrapInZR(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Connect containers in R
///
/// @param gctx The geometry context
/// @param containers the containers
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note not much checking is done anymore, as the DetectorComponent::PortalContainer
/// are assumed to come properly formed from the prior methods
///
/// @return a proto container with the outside portals
DetectorComponent::PortalContainer connectInR(
    const GeometryContext& gctx,
    const std::vector<DetectorComponent::PortalContainer>& containers,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Connect containers in Z
///
/// @param gctx The geometry context
/// @param containers the containers
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note not much checking is done anymore, as the DetectorComponent::PortalContainer
/// are assumed to come properly formed from the prior methods
///
/// @return a proto container with the outside portals
DetectorComponent::PortalContainer connectInZ(
    const GeometryContext& gctx,
    const std::vector<DetectorComponent::PortalContainer>& containers,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Connect containers in Phi
///
/// @param gctx The geometry context
/// @param containers the containers
/// @param selectedOnly switch only selected boundaries
/// @param logLevel is the screen logging level
///
/// @note not much checking is done anymore, as the DetectorComponent::PortalContainer
/// are assumed to come properly formed from the prior methods
///
/// @return a proto container with the outside portals
DetectorComponent::PortalContainer connectInPhi(
    const GeometryContext& gctx,
    const std::vector<DetectorComponent::PortalContainer>& containers,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// @brief Wrap container in R,Z - this uses the cutout cylinder bounds
///
/// @param gctx The geometry context
/// @param containers the containers, i.e. the inner volume and the wrapping container
/// @param logLevel is the screen logging level
///
/// @note not much checking is done anymore, as the DetectorComponent::PortalContainer
/// are assumed to come properly formed from the prior methods
///
/// @return a proto container with the outside portals
DetectorComponent::PortalContainer wrapInZR(
    const GeometryContext& gctx,
    const std::vector<DetectorComponent::PortalContainer>& containers,
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

}  // namespace CylindricalDetectorHelper
}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
