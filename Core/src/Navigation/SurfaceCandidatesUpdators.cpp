// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"

namespace {

/// Helper method to update the candidates (portals/surfaces),
/// this can be called for initial surface/portal estimation,
/// but also during the navigation to update the current list
/// of candidates.
///
/// @param gctx is the Geometry context of this call
/// @param nState [in,out] is the navigation state to be updated
///
/// @todo for surfaces skip the non-reached ones, while keep for portals
void updateCandidates(const Acts::GeometryContext& gctx,
                      Acts::Experimental::NavigationState& nState) {
  const auto& position = nState.position;
  const auto& direction = nState.direction;
  auto& nCandidates = nState.surfaceCandidates;

  for (auto& c : nCandidates) {
    // Get the surface representation: either native surfcae of portal
    const Acts::Surface& sRep =
        (c.surface != nullptr) ? (*c.surface) : (c.portal->surface());

    // Get the intersection @todo make a templated intersector
    // TODO surface tolerance
    auto sIntersection = sRep.intersect(
        gctx, position, direction, c.boundaryCheck, Acts::s_onSurfaceTolerance);
    // Re-order and swap if necessary
    if (sIntersection.intersection.pathLength + Acts::s_onSurfaceTolerance <
            nState.overstepTolerance and
        sIntersection.alternative.status >=
            Acts::Intersection3D::Status::reachable) {
      sIntersection.swapSolutions();
    }
    c.objectIntersection = sIntersection;
  }
  // Sort and stuff non-allowed solutions to the end
  std::sort(
      nCandidates.begin(), nCandidates.end(),
      [&](const auto& a, const auto& b) {
        // The two path lengths
        Acts::ActsScalar pathToA = a.objectIntersection.intersection.pathLength;
        Acts::ActsScalar pathToB = b.objectIntersection.intersection.pathLength;
        if (pathToA + Acts::s_onSurfaceTolerance < nState.overstepTolerance or
            std::abs(pathToA) < Acts::s_onSurfaceTolerance) {
          return false;
        } else if (pathToB + Acts::s_onSurfaceTolerance <
                       nState.overstepTolerance or
                   std::abs(pathToB) < Acts::s_onSurfaceTolerance) {
          return true;
        }
        return pathToA < pathToB;
      });
  // Set the surface candidate
  nState.surfaceCandidate = nCandidates.begin();
}

}  // namespace

void Acts::Experimental::AllPortalsImpl::update(const GeometryContext& gctx,
                                                NavigationState& nState) const {
  if (nState.currentVolume == nullptr) {
    throw std::runtime_error(
        "AllPortalsImpl: no detector volume set to navigation state.");
  }
  // Retrieval necessary
  if (nState.surfaceCandidates.empty()) {
    // Fill internal portals if existing
    for (const auto v : nState.currentVolume->volumes()) {
      const auto& iPortals = v->portals();
      PortalsFiller::fill(nState, iPortals);
    }
    // Filling the new portal candidates
    const auto& portals = nState.currentVolume->portals();
    PortalsFiller::fill(nState, portals);
  }
  // Sort and update
  updateCandidates(gctx, nState);
}

void Acts::Experimental::AllPortalsAndSurfacesImpl::update(
    const GeometryContext& gctx, NavigationState& nState) const {
  if (nState.currentDetector == nullptr) {
    throw std::runtime_error(
        "AllPortalsAndSurfacesImpl: no detector volume set to navigation "
        "state.");
  }
  // A volume switch has happened, update list of portals & surfaces
  if (nState.surfaceCandidates.empty()) {
    // Fill internal portals if existing
    for (const auto v : nState.currentVolume->volumes()) {
      const auto& iPortals = v->portals();
      PortalsFiller::fill(nState, iPortals);
    }
    // Assign the new volume
    const auto& portals = nState.currentVolume->portals();
    const auto& surfaces = nState.currentVolume->surfaces();
    PortalsFiller::fill(nState, portals);
    SurfacesFiller::fill(nState, surfaces);
  }
  // Update internal candidates
  updateCandidates(gctx, nState);
}
