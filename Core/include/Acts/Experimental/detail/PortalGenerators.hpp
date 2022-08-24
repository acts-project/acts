// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Experimental/detail/DetectorVolumeLinks.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {
namespace Experimental {
namespace detail {

/// @brief Generator function for creation of portal surfaces
/// for a cylindrical volume
///
/// @param dTransform the context-resolved transform of the detector volume
/// @param dBounds the detecor volume bounds
/// @param dVolume the reference to the detector volume
///
/// @return a vector of newly created portals with registered inside volume
inline static std::vector<std::shared_ptr<Portal>> portals(
    const Transform3& dTransform, const VolumeBounds& dBounds,
    const DetectorVolume& dVolume) {
  // Get the oriented boundary surfaces
  auto orientedSurfaces = dBounds.orientedSurfaces(dTransform);

  // Create and fill the portal return vector
  std::vector<std::shared_ptr<Portal>> portals;
  for (auto [i, oSurface] : enumerate(orientedSurfaces)) {
    // Create a portal from the surface
    auto portal = Portal::makeShared(oSurface.first);
    // Create a shared link instance & delegate
    auto singleLinkStored =
        std::make_shared<SingleLinkImpl>(SingleLinkImpl{&dVolume});
    DetectorVolumeLink singleLink;
    singleLink.connect<&SingleLinkImpl::targetVolume>(singleLinkStored.get());
    // Update the volume link and the store
    Acts::NavigationDirection insideDir = oSurface.second;
    portal->updateVolumeLink(insideDir, singleLink, singleLinkStored, true);
    // Create a null link
    DetectorVolumeLink nullLink;
    nullLink.connect<&nullVolumeLink>();
    Acts::NavigationDirection outsideDir =
        (insideDir == Acts::NavigationDirection::Forward)
            ? Acts::NavigationDirection::Backward
            : Acts::NavigationDirection::Forward;
    portal->updateVolumeLink(outsideDir, nullLink);
    // Portal is prepared
    portals.push_back(portal);
  }
  return portals;
};

/// Create a default portal generator that connects to the
/// static method.
///
inline static Delegate<std::vector<std::shared_ptr<Portal>>(
    const Transform3&, const VolumeBounds&, const DetectorVolume&)>
defaultPortalGenerator() {
  Delegate<std::vector<std::shared_ptr<Portal>>(
      const Transform3&, const VolumeBounds&, const DetectorVolume&)>
      pGenerator;
  pGenerator.connect<&portals>();
  return pGenerator;
}

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts