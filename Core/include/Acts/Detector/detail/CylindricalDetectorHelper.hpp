// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
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
/// @tparam volume_container_t the type of the container
///
/// @param gctx the geometry context of the call
/// @param volumes the volumes at input
/// @param logLevel is the screen logging level
///
/// @return extracted boundary values
template <typename volume_container_t>
std::array<std::vector<ActsScalar>, 3u> rzphiBoundaries(
    const GeometryContext& gctx, const volume_container_t& volumes,
    Acts::Logging::Level logLevel = Acts::Logging::INFO) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Estimate R/z/phi boundaries of  " << volumes.size()
                                                << " volumes.");

  // The return boundaries
  std::array<std::vector<ActsScalar>, 3u> boundaries;

  // The map for collecting
  std::array<std::map<ActsScalar, size_t>, 3u> valueMaps;
  auto& rMap = valueMaps[0u];
  auto& zMap = valueMaps[1u];
  auto& phiMap = valueMaps[2u];

  auto fillMap = [&](std::map<ActsScalar, size_t>& map,
                     const std::array<ActsScalar, 2u>& values) {
    for (auto v : values) {
      if (map.find(v) != map.end()) {
        ++map[v];
      } else {
        map[v] = 1u;
      }
    }
  };

  // Loop over the volumes and collect boundaries
  for (const auto& v : volumes) {
    if (v->volumeBounds().type() == VolumeBounds::BoundsType::eCylinder) {
      auto bValues = v->volumeBounds().values();
      // The min/max values
      ActsScalar rMin = bValues[CylinderVolumeBounds::BoundValues::eMinR];
      ActsScalar rMax = bValues[CylinderVolumeBounds::BoundValues::eMaxR];
      ActsScalar zCenter = v->transform(gctx).translation().z();
      ActsScalar zHalfLength =
          bValues[CylinderVolumeBounds::BoundValues::eHalfLengthZ];
      ActsScalar zMin = zCenter - zHalfLength;
      ActsScalar zMax = zCenter + zHalfLength;
      ActsScalar phiCenter =
          bValues[CylinderVolumeBounds::BoundValues::eAveragePhi];
      ActsScalar phiSector =
          bValues[CylinderVolumeBounds::BoundValues::eHalfPhiSector];
      ActsScalar phiMin = phiCenter - phiSector;
      ActsScalar phiMax = phiCenter + phiSector;
      // Fill the maps
      fillMap(rMap, {rMin, rMax});
      fillMap(zMap, {zMin, zMax});
      fillMap(phiMap, {phiMin, phiMax});
    }
  }

  for (auto [im, map] : enumerate(valueMaps)) {
    for (auto [key, value] : map) {
      boundaries[im].push_back(key);
    }
    std::sort(boundaries[im].begin(), boundaries[im].end());
  }

  ACTS_VERBOSE("- did yield " << boundaries[0u].size() << " boundaries in R.");
  ACTS_VERBOSE("- did yield " << boundaries[1u].size() << " boundaries in z.");
  ACTS_VERBOSE("- did yield " << boundaries[2u].size()
                              << " boundaries in phi.");

  return boundaries;
}

}  // namespace CylindricalDetectorHelper
}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
