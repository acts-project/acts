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
#include "Acts/Experimental/detail/DetectorVolumeLinks.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {
namespace detail {

using PortalReplacement =
    std::tuple<std::shared_ptr<Experimental::Portal>, unsigned int,
               NavigationDirection, std::vector<ActsScalar>, BinningValue>;

/// @brief Create and attach the multi links
///
/// @param gctx the geometry context
/// @param volumes are the volumes that are pointed to
/// @param pReplacements are the portal replacements that are newly connected
///
/// @note this helper tries to creates a multi link with least transformation
/// involved
///
void attachMultiLinks(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    std::vector<PortalReplacement>& pReplacements) {
  // unpack to navigation bare points
  auto constVolumes = unpack_shared_const_vector(volumes);
  // Set to the contructed portals (p), at index (i), in direction (d)
  // using boundaries and binning
  for (auto& [p, i, dir, boundaries, binning] : pReplacements) {
    DetectorVolumeLink volumeLink;
    ManagedDetectorVolumeLink managedLink;
    // Check if the boundaries need a transform
    std::vector<ActsScalar> cBoundaries = boundaries;
    const auto pTransform = p->surface().transform(gctx);
    const auto pZ = pTransform.rotation().col(2);
    bool aligned = pTransform.isApprox(Transform3::Identity());
    bool alongZ =
        pZ.isApprox(Vector3(0., 0., 1.)) or pZ.isApprox(Vector3(0., 0., -1.));
    // Different checks for different binnings
    switch (binning) {
      case binR: {
        aligned =
            alongZ and Acts::VectorHelpers::perp(pTransform.translation()) <
                           std::numeric_limits<ActsScalar>::epsilon();
      } break;
      case binZ: {
        // z shift can be compensated
        if (not aligned and alongZ) {
          aligned = true;
          // Apply shift in z
          ActsScalar zOffset = pTransform.translation().z();
          std::for_each(cBoundaries.begin(), cBoundaries.end(),
                        [&](ActsScalar& z) { z += zOffset; });
        }
      } break;
      default:
        break;
    }

    // No transform necessary or already corrected
    if (aligned) {
      // Multi link implementation
      auto multiLinkImpl = std::make_shared<detail::MultiLink1DImpl>(
          constVolumes, cBoundaries, binning);
      // Create the delegate and its managed object
      volumeLink.connect<&detail::MultiLink1DImpl::targetVolume>(
          multiLinkImpl.get());
      managedLink = Acts::Experimental::ManagedDetectorVolumeLink{
          std::move(volumeLink), multiLinkImpl};
    } else {
      // Multilink with transform necessary
      detail::MultiLink1DImpl multiLinkImpl(constVolumes, cBoundaries, binning);
      auto transformedMuliLinkImpl =
          std::make_shared<detail::TransformedMultiLink1DImpl>(
              multiLinkImpl, pTransform.inverse());
      // Create the delegate and its managed object
      volumeLink.connect<&detail::TransformedMultiLink1DImpl::targetVolume>(
          transformedMuliLinkImpl.get());
      managedLink = Acts::Experimental::ManagedDetectorVolumeLink{
          std::move(volumeLink), transformedMuliLinkImpl};
    }
    // Now update the link
    p->updateVolumeLink(dir, std::move(managedLink), volumes);
  }
}

/// @brief  Method that strips out attached volumes from portals
///
/// @note it throws an exception if both sides are already taken
///
/// @param portal the portal to be resolved
///
/// @return a vector of attached volumes
std::vector<std::shared_ptr<DetectorVolume>> attachedVolumes(
    Portal& portal) noexcept(false) {
  auto& attachedVolumes = portal.attachedVolumes();
  if (not attachedVolumes[0u].empty() and not attachedVolumes[1u].empty()) {
    throw std::invalid_argument(
        "PortalHelper: trying to get attachedVolumes from already populated "
        "portal.");
  }
  unsigned int iu = attachedVolumes[0u].empty() ? 1u : 0u;
  return attachedVolumes[iu];
}

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts