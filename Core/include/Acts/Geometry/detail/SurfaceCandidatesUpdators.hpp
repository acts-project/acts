// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationState.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/detail/NavigationStateUpdators.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <tuple>

namespace Acts {
namespace Experimental {
namespace detail {

/// Helper method to update the candidates (portals/surfaces),
/// this can be called for initial surface/portal estimation,
/// but also during the navigation to update the current list
/// of candidates.
///
/// @param gctx is the Geometry context of this call
/// @param nState [in,out] is the navigation state to be updated
///
/// @todo for surfaces skip the non-reached ones, while keep for portals
inline static void updateCandidates(const GeometryContext& gctx,
                                    NavigationState& nState) {
  const auto& position = nState.position;
  const auto& direction = nState.direction;
  auto& nCandidates = nState.surfaceCandidates;

  for (auto& c : nCandidates) {
    // Get the surface reprensentation: either native surfcae of portal
    const Surface& sRep =
        (c.surface != nullptr) ? (*c.surface) : (c.portal->surface());

    // Get the intersection @todo make a templated intersector
    auto sIntersection = sRep.intersect(gctx, position, direction, c.bCheck);
    // Re-order and swap if necessary
    if (sIntersection.intersection.pathLength + s_onSurfaceTolerance <
            nState.overstepTolerance and
        sIntersection.alternative.status >= Intersection3D::Status::reachable) {
      sIntersection.swapSolutions();
    }
    c.objectIntersection = sIntersection;
  }
  // Sort and stuff non-allowed solutions to the end
  std::sort(
      nCandidates.begin(), nCandidates.end(),
      [&](const auto& a, const auto& b) {
        // The two path lengths
        ActsScalar pathToA = a.objectIntersection.intersection.pathLength;
        ActsScalar pathToB = b.objectIntersection.intersection.pathLength;
        if (pathToA + s_onSurfaceTolerance < nState.overstepTolerance or
            std::abs(pathToA) < s_onSurfaceTolerance) {
          return false;
        } else if (pathToB + s_onSurfaceTolerance < nState.overstepTolerance or
                   std::abs(pathToB) < s_onSurfaceTolerance) {
          return true;
        }
        return pathToA < pathToB;
      });
  // Set the surface candidate
  nState.surfaceCandidate = nCandidates.begin();
}

/// A ordered portal provider
///
/// @param gctx is the Geometry context of this call
/// @param nState is the navigation state to be updated
///
/// @note that the intersections are ordered, such that the
/// smallest intersection pathlength >= overstep tolerance is the lowest
///
/// @return an ordered list of portal candidates
inline static void portalCandidates(const GeometryContext& gctx,
                                    NavigationState& nState) noexcept(false) {
  if (nState.currentVolume == nullptr) {
    throw std::runtime_error(
        "PortalCandidates: no detector volume set to navigation state.");
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

/// A ordered portal and surface provider
///
/// @param gctx is the Geometry context of this call
/// @param nState is the navigation state to be updated
///
/// @note that the intersections are ordered, such that the
/// smallest intersection pathlength >= overstep tolerance is the lowest
///
/// @return an ordered list of portal candidates
inline static void portalAndSurfaceCandidates(const GeometryContext& gctx,
                                              NavigationState& nState) {
  if (nState.currentDetector == nullptr) {
    throw std::runtime_error(
        "PortalAndSurfaceCandidates: no detector volume set to navigation "
        "state.");
  }
  // A volu
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

/// Generate a provider for all portals
///
/// @return a connected navigationstate updator
inline static ManagedSurfaceCandidatesUpdator allPortals() {
  ManagedSurfaceCandidatesUpdator managedUpdator;
  SurfaceCandidatesUpdator nStateUpdator;
  nStateUpdator.connect<&portalCandidates>();
  managedUpdator.delegate = std::move(nStateUpdator);
  managedUpdator.implementation = nullptr;
  return managedUpdator;
}

/// Generate a provider for all portals and Surfacess
///
/// @note this is a try-and error navigation, not recommended for production
/// setup with many surfaces
///
/// @return a connected navigationstate updator
inline static ManagedSurfaceCandidatesUpdator allPortalsAndSurfaces() {
  ManagedSurfaceCandidatesUpdator managedUpdator;
  SurfaceCandidatesUpdator nStateUpdator;
  nStateUpdator.connect<&portalAndSurfaceCandidates>();
  managedUpdator.delegate = std::move(nStateUpdator);
  managedUpdator.implementation = nullptr;
  return managedUpdator;
}

/// @brief This holds and extracts a collection of surfaces without much
/// checking, this could be e.g. support surfaces for layer structures,
/// e.g.
///
struct AdditionalSurfacesImpl : public IDelegateImpl {
  /// The volumes held by this collection
  std::vector<const Surface*> surfaces = {};

  /// Extract the sub volumes from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call (ignored)
  /// @param nState is the navigation state
  ///
  inline void update([[maybe_unused]] const GeometryContext& gctx,
                     NavigationState& nState) const {
    SurfacesFiller::fill(nState, surfaces);
  }
};

/// @brief  An indexed surface implementation access
///
/// @tparam grid_type is the grid type used for this
template <typename grid_type>
using IndexedSurfacesImpl =
    IndexedUpdatorImpl<grid_type, IndexedSurfacesExtractor, SurfacesFiller>;

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
