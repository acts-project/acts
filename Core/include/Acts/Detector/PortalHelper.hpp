// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <exception>
#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {

/// @brief Calls the portal generation and adds registration to sub portals
///
/// This code is split off the PortalGenerator code in order to allow
/// unit testing of the portal generation wihtout detector volume construction
///
/// @param dTransform a contextually resolved transform
/// @param dBounds the detecor volume bounds
/// @param dVolume the reference to the detector volume which generates this volume
///
/// @return a vector of newly created portals with registered inside volume
inline static std::vector<std::shared_ptr<Portal>>
generatePortalsUpdateInternals(
    const Transform3& dTransform, const VolumeBounds& dBounds,
    const std::shared_ptr<DetectorVolume>& dVolume) noexcept(false) {
  if (dVolume == nullptr) {
    throw std::runtime_error(
        "generatePortalsUpdateInternals: no detector volume provided.");
  }

  // Setting link to the mother volume to all sub volumes of this volume
  for (auto& vPtr : dVolume->volumePtrs()) {
    for (auto& pPtr : vPtr->portalPtrs()) {
      // Creating a link to the mother
      auto motherLinkImpl =
          std::make_unique<const SingleDetectorVolumeImpl>(dVolume.get());
      DetectorVolumeUpdator motherLink;
      motherLink.connect<&SingleDetectorVolumeImpl::update>(
          std::move(motherLinkImpl));
      pPtr->assignDetectorVolumeUpdator(std::move(motherLink), {dVolume});
    }
  }
  // Return from the standard generator
  return generatePortals(dTransform, dBounds, dVolume);
}

/// Create a default portal generator that connects to the
/// static method.
///
/// @note parameters are ignored in this case
inline static Delegate<std::vector<std::shared_ptr<Portal>>(
    const Transform3&, const VolumeBounds&,
    const std::shared_ptr<DetectorVolume>&)>
defaultPortalAndSubPortalGenerator() {
  Delegate<std::vector<std::shared_ptr<Portal>>(
      const Transform3&, const VolumeBounds&,
      const std::shared_ptr<DetectorVolume>&)>
      pGenerator;
  pGenerator.connect<&generatePortalsUpdateInternals>();
  return pGenerator;
}

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
