// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Geometry/detail/DetectorVolumeUpdators.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <exception>
#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {
namespace detail {

/// @brief Generator function for creation of portal surfaces
/// for a cylindrical volume
///
/// @param dTransform a contextually resolved transform
/// @param dBounds the detecor volume bounds
/// @param dVolume the reference to the detector volume which generates this volume
///
/// @return a vector of newly created portals with registered inside volume
inline static std::vector<std::shared_ptr<Portal>> generatePortals(
    const Transform3& dTransform, const VolumeBounds& dBounds,
    const std::shared_ptr<DetectorVolume>& dVolume) noexcept(false) {
  if (dVolume == nullptr) {
    throw std::runtime_error("PortalsGenerator: no detector volume provided.");
  }
  // Get the oriented boundary surfaces
  auto orientedSurfaces = dBounds.orientedSurfaces(dTransform);

  // Create and fill the portal return vector
  std::vector<std::shared_ptr<Portal>> portals;
  for (auto [i, oSurface] : enumerate(orientedSurfaces)) {
    // Create a portal from the surface
    auto portal = Portal::makeShared(oSurface.first);
    // Create a shared link instance & delegate
    auto singleLinkImpl =
        std::make_unique<const SingleDetectorVolumeImpl>(dVolume.get());
    DetectorVolumeUpdator singleLink;
    singleLink.connect<&SingleDetectorVolumeImpl::update>(
        std::move(singleLinkImpl));
    // Update the volume link and the store
    NavigationDirection insideDir = oSurface.second;
    portal->assignDetectorVolumeUpdator(insideDir, std::move(singleLink),
                                        {dVolume});
    // Portal is prepared
    portals.push_back(std::move(portal));
  }

  // The portals are returned
  return portals;
}

/// Create a default portal generator that connects to the
/// static method.
///
inline static Delegate<std::vector<std::shared_ptr<Portal>>(
    const Transform3&, const VolumeBounds&,
    const std::shared_ptr<DetectorVolume>&)>
defaultPortalGenerator() {
  Delegate<std::vector<std::shared_ptr<Portal>>(
      const Transform3&, const VolumeBounds&,
      const std::shared_ptr<DetectorVolume>&)>
      pGenerator;
  pGenerator.connect<&generatePortals>();
  return pGenerator;
}

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts