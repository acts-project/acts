// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {

/// Definition of a portal replacement when building proto containers
/// It consists of the new portal, the index, the direction, the parameters
/// gathered from the sub volumes, the binning description
using PortalReplacement =
    std::tuple<std::shared_ptr<Experimental::Portal>, unsigned int,
               NavigationDirection, std::vector<ActsScalar>, BinningValue>;

/// @brief Create and attach the multi link updator
///
/// @param gctx the geometry context
/// @param volumes are the volumes that are pointed to
/// @param pReplacements are the portal replacements that are newly connected
///
///
void attachDetectorVolumeUpdators(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    std::vector<PortalReplacement>& pReplacements) {
  // Unpack to navigation bare points
  auto cVolumes = unpack_shared_const_vector(volumes);
  // Set to the contructed portals (p), at index (i), in direction (d)
  // using boundaries and binning
  for (auto& [p, i, dir, boundaries, binning] : pReplacements) {
    // Check if the boundaries need a transform
    const auto pTransform = p->surface().transform(gctx);
    // Creating a link to the mother
    auto volumes1D = std::make_unique<const BoundVolumesGrid1Impl>(
        boundaries, binning, cVolumes, pTransform.inverse());
    DetectorVolumeUpdator dVolumeUpdator;
    dVolumeUpdator.connect<&BoundVolumesGrid1Impl::update>(
        std::move(volumes1D));
    p->assignDetectorVolumeUpdator(dir, std::move(dVolumeUpdator), volumes);
  }
}

/// @brief  Method that strips out attached volumes from portals
///
/// @note it throws an exception if both sides are already taken
///
/// @param portal the portal to be resolved
///
/// @return a vector of attached volumes
std::vector<std::shared_ptr<DetectorVolume>> attachedDetectorVolumes(
    Portal& portal) noexcept(false) {
  auto& attachedVolumes = portal.attachedDetectorVolumes();
  if (not attachedVolumes[0u].empty() and not attachedVolumes[1u].empty()) {
    throw std::invalid_argument(
        "PortalHelper: trying to get attachedVolumes from already populated "
        "portal.");
  }
  unsigned int iu = attachedVolumes[0u].empty() ? 1u : 0u;
  return attachedVolumes[iu];
}

}  // namespace Experimental
}  // namespace Acts
