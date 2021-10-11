// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Experimental/DetectorEnvironment.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace Acts {

class DetectorVolume;

/// Helper function to get the portal candidates from a volume
///
/// @param gctx is the current geometry conbtext
/// @param volume is the volume for which the portals are intersected
/// @param position is the position at the query
/// @param direction is the direction at the query
///
/// @note onSurface solutions are ranked last
inline std::vector<PortalIntersection> portalCandidates(
    const GeometryContext& gctx, const DetectorVolume& volume,
    const Vector3& position, const Vector3& direction) {
  // The portal intersections
  std::vector<PortalIntersection> pIntersections;
  // Get all the portals
  const auto& portals = volume.portals();
  pIntersections.reserve(portals.size());
  // Loop over portals an intersect
  for (const auto& p : portals) {
    pIntersections.push_back(p->intersect(gctx, position, direction));
  }
  // Sort and push negative and on surface solutions to the end
  std::sort(pIntersections.begin(), pIntersections.end(),
            [](const auto& a, const auto& b) {
              if (a.intersection.pathLength < s_onSurfaceTolerance) {
                return false;
              }
              if (b.intersection.pathLength < s_onSurfaceTolerance) {
                return true;
              }
              return (a < b);
            });
  // Return the sorted solutions
  return pIntersections;
};

/// A portal link to a single volume
struct SinglePortalLink {
  /// The simple link to the detector volume
  const DetectorVolume* volume = nullptr;

  /// The link to the (optional) surfaces contained
  SurfaceLinks surfaces = VoidSurfaceLink{};

  /// Fullfills the call std::function call structure
  ///
  /// @param gctx the current geometry context
  /// @param portal the portal at this request (ignored)
  /// @param position the current position
  /// @param direction the current direction
  /// @param bCheck boundary check for surface search
  /// @param provideAll flag for trail&error navigation
  ///
  /// @return a new environment for the navigation
  DetectorEnvironment operator()(const GeometryContext& gctx,
                                 const Portal& /*portal*/,
                                 const Vector3& position,
                                 const Vector3& direction,
                                 const BoundaryCheck& bCheck,
                                 bool provideAll = false) const {
    // The portals one-time intersected
    std::vector<PortalIntersection> pCandidates = {};
    if (volume != nullptr) {
      // Get the portals coordinates (ordered)
      pCandidates = portalCandidates(gctx, *volume, position, direction);

      // Maximum path length
      ActsScalar maximumPath = pCandidates[0].intersection.pathLength;

      // The surface candidates one-time intersected & ordered
      std::vector<SurfaceIntersection> sCandidates = surfaces(
          gctx, *volume, position, direction, bCheck, maximumPath, provideAll);
      // Return the new environment
      return DetectorEnvironment{volume, sCandidates, pCandidates};
    }
    return DetectorEnvironment{};
  }
};

/// A portal link into multiple volumes
struct MultiplePortalLink {
  /// The volume links
  std::vector<PortalLink> portalLinks = {};
  /// The volume association
  VolumeLink index = VoidVolumeLink{};

  /// Fullfills the call std::function call structure
  ///
  /// @param gctx the current geometry context
  /// @param portal the portal at this request
  /// @param position the current position
  /// @param direction the current direction
  /// @param bCheck boundary check for surface search
  /// @param provideAll flag for trail&error navigation
  ///
  /// @return a new environment for the navigation
  DetectorEnvironment operator()(const GeometryContext& gctx,
                                 const Portal& portal, const Vector3& position,
                                 const Vector3& direction,
                                 const BoundaryCheck& bCheck,
                                 bool provideAll = false) const {
    // The portals one-time intersected
    const Surface& surface = portal.surfaceRepresentation();
    unsigned int vIndex = index(surface.transform(gctx), position);
    if (vIndex < portalLinks.size()) {
      return portalLinks[vIndex](gctx, portal, position, direction, bCheck,
                                 provideAll);
    }
    return DetectorEnvironment{};
  }
};

}  // namespace Acts
