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

#include <array>
#include <functional>
#include <memory>
#include <vector>

namespace Acts {

class DetectorVolume;

/// Helper function to get the portal candidates from a volume
///
/// @param gctx is the current geometry conbtext
/// @param portals is the volume for which the portals are intersected
/// @param position is the position at the query
/// @param direction is the direction at the query
/// @param pathRange is the allowed path range for this search
///
/// @note onSurface solutions are ranked last
inline std::vector<PortalIntersection> portalCandidates(
    const GeometryContext& gctx, const std::vector<const Portal*>& portals,
    const Vector3& position, const Vector3& direction,
    const std::array<ActsScalar, 2>& pathRange = {
        0., std::numeric_limits<ActsScalar>::infinity()}) {
  // The portal intersections
  std::vector<PortalIntersection> pIntersections;
  // Get all the portals
  pIntersections.reserve(portals.size());
  // Loop over portals an intersect
  for (const auto& p : portals) {
    // Get the intersection
    auto pIntersection = p->intersect(gctx, position, direction);
    // Re-order if necessary
    if (pIntersection.intersection.pathLength + s_onSurfaceTolerance <
            pathRange[0] and
        pIntersection.alternative.pathLength + s_onSurfaceTolerance >
            pathRange[0] and
        pIntersection.alternative.status >= Intersection3D::Status::reachable) {
      // Let's swap the solutions
      pIntersection.swapSolutions();
    }
    // Exclude on-portal solution
    if (std::abs(pIntersection.intersection.pathLength) < s_onSurfaceTolerance){
      continue;
    }
    pIntersections.push_back(pIntersection);
  }
  // Sort and non-allowed solutions to the end
  std::sort(
      pIntersections.begin(), pIntersections.end(),
      [&](const auto& a, const auto& b) {
        if (a.intersection.pathLength + s_onSurfaceTolerance < pathRange[0]) {
          return false;
        } else if (b.intersection.pathLength + s_onSurfaceTolerance <
                   pathRange[0]) {
          return true;
        }
        return a.intersection.pathLength < b.intersection.pathLength;
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
                                 const Portal& portal, const Vector3& position,
                                 const Vector3& direction,
                                 const BoundaryCheck& bCheck,
                                 bool provideAll = false) const {
    // The portals one-time intersected
    std::vector<PortalIntersection> pCandidates = {};
    if (volume != nullptr) {
      // Get the portals coordinates (ordered)
      pCandidates =
          portalCandidates(gctx, volume->portals(), position, direction);

      // Maximum path length - if the intersection is on portal it needs to go
      ActsScalar maximumPath =
          std::abs(pCandidates[0].intersection.pathLength) >
                  s_onSurfaceTolerance
              ? pCandidates[0].intersection.pathLength
              : pCandidates[1].intersection.pathLength;

      // The surface candidates one-time intersected & ordered
      std::vector<SurfaceIntersection> sCandidates =
          surfaces(gctx, *volume, position, direction, bCheck,
                   {0., maximumPath}, provideAll);
      // Return the new environment
      DetectorEnvironment environment{volume, sCandidates, pCandidates};
      environment.currentSurface = &(portal.surfaceRepresentation());
      environment.status = DetectorEnvironment::eOnPortal;
      return environment;
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
    unsigned int vIndex = index(position);
    if (vIndex < portalLinks.size()) {
      return portalLinks[vIndex](gctx, portal, position, direction, bCheck,
                                 provideAll);
    }
    return DetectorEnvironment{};
  }
};

}  // namespace Acts
