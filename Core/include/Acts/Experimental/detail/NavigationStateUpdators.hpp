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
#include "Acts/Experimental/NavigationState.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <tuple>

namespace Acts {
namespace Experimental {
namespace detail {

/// A candidate updator
///
/// @param nState is the navigation state to be updated
/// @param gctx is the Geometry context of this call
/// @param position is the position of the update call
/// @param direction is the direction of the update call
/// @param absMomentum the absolute momentum
/// @param charge the charge parameters
inline static void updateCandidates(NavigationState& nState,
                                    const GeometryContext& gctx,
                                    const Vector3& position,
                                    const Vector3& direction,
                                    [[maybe_unused]] ActsScalar absMomentum,
                                    [[maybe_unused]] ActsScalar charge) {
  for (auto& c : nState.surfaceCandidates) {
    // Get the surface reprensentation
    const Surface* sRep =
        (c.surface != nullptr) ? c.surface : &(c.portal->surface());

    // Get the intersection @todo make a templated intersector
    auto sIntersection = sRep->intersect(gctx, position, direction, c.bCheck);
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
      nState.surfaceCandidates.begin(), nState.surfaceCandidates.end(),
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
  nState.surfaceCandidate = nState.surfaceCandidates.begin();
}

/// A ordered portal provider
///
/// @param nState is the navigation state to be updated
/// @param volume is the detector volume
/// @param gctx is the Geometry context of this call
/// @param position is the position of the update call
/// @param direction is the direction of the update call
/// @param absMomentum the absolute momentum
/// @param charge the charge parameters
/// @param clear is a boolean telling if you should clear the canidates
///
/// @note that the intersections are ordered, such that the
/// smallest intersection pathlength >= overstep tolerance is the lowest
///
/// @return an ordered list of portal candidates
inline static void portalCandidates(NavigationState& nState,
                                    const DetectorVolume& volume,
                                    const GeometryContext& gctx,
                                    const Vector3& position,
                                    const Vector3& direction,
                                    ActsScalar absMomentum, ActsScalar charge) {
  // A volume switch has happened, update list of portals
  if (nState.currentVolume != &volume) {
    // Assign the new volume
    nState.currentVolume = &volume;
    const auto& portals = volume.portals();
    nState.surfaceCandidates.clear();
    nState.surfaceCandidate = nState.surfaceCandidates.end();
    nState.surfaceCandidates.reserve(portals.size());
    for (const auto* p : portals) {
      nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
          ObjectIntersection<Surface>{}, nullptr, p, true});
    }
    return;
  }
  // Update internal candidates
  updateCandidates(nState, gctx, position, direction, absMomentum, charge);
}

/// Generate a default portal provider
///
/// @return a connected navigationstate updator
inline static NavigationStateUpdator defaultPortalProvider() {
  NavigationStateUpdator nStateUpdator;
  nStateUpdator.connect<&portalCandidates>();
  return nStateUpdator;
}

struct AllSurfacesAttacher {
  /// @brief An attacher of all surface candidates in a volume
  ///
  /// @param nState the navigation state to which the surfaces are attached
  /// @param volume is the detector volume
  /// @param gctx is the Geometry context of this call
  /// @param position is the position of the update call
  /// @param direction is the direction of the update call
  /// @param absMomentum the absolute momentum
  /// @param charge the charge parameters
  ///
  /// @note this is attaching objects without intersecting nor checking
  void operator()(NavigationState& nState, const DetectorVolume& volume,
                  [[maybe_unused]] const GeometryContext& gctx,
                  [[maybe_unused]] const Vector3& position,
                  [[maybe_unused]] const Vector3& direction,
                  [[maybe_unused]] ActsScalar absMomentum,
                  [[maybe_unused]] ActsScalar charge) const {
    for (const auto* s : volume.surfaces()) {
      nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
          ObjectIntersection<Surface>{}, s, nullptr,
          nState.surfaceBoundaryCheck});
    }
  }
};

struct AllPortalsAttacher {
  /// @brief An attacher of all portal surface candidates of internal volumes
  ///
  /// @param nState the navigation state to which the surfaces are attached
  /// @param volume is the detector volume
  /// @param gctx is the Geometry context of this call
  /// @param position is the position of the update call
  /// @param direction is the direction of the update call
  /// @param absMomentum the absolute momentum
  /// @param charge the charge parameters
  ///
  /// @note this is attaching objects without intersecting nor checking
  void operator()(NavigationState& nState, const DetectorVolume& volume,
                  [[maybe_unused]] const GeometryContext& gctx,
                  [[maybe_unused]] const Vector3& position,
                  [[maybe_unused]] const Vector3& direction,
                  [[maybe_unused]] ActsScalar absMomentum,
                  [[maybe_unused]] ActsScalar charge) const {
    for (const auto* v : volume.volumes()) {
      for (const auto* p : v->portals()) {
        nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
            ObjectIntersection<Surface>{}, nullptr, p, true});
      }
    }
  }
};

// This struct allows to combine the portals from the volume itself
// with  attachers
template <typename... candidate_attachers_t>
struct NavigationStateUpdator {
  std::tuple<candidate_attachers_t...> m_candidateAttachers;

  NavigationStateUpdator(const std::tuple<candidate_attachers_t...>& attachers)
      : m_candidateAttachers(attachers) {}

  /// A combined portal provider
  ///
  /// @param nState is the navigation state to be updated
  /// @param volume is the detector volume
  /// @param gctx is the Geometry context of this call
  /// @param position is the position of the update call
  /// @param direction is the direction of the update call
  /// @param absMomentum the absolute momentum
  /// @param charge the charge parameters
  ///
  /// @note that the intersections are ordered, such that the
  /// smallest intersection pathlength >= overstep tolerance is the lowe
  void update(NavigationState& nState, const DetectorVolume& volume,
              const GeometryContext& gctx, const Vector3& position,
              const Vector3& direction, ActsScalar absMomentum,
              ActsScalar charge) const {
    // The portal candidates
    if (nState.currentVolume != &volume) {
      // Let's add the portal candidates of this volume
      portalCandidates(nState, volume, gctx, position, direction, absMomentum,
                       charge);
      // Unfold the tuple and add the attachers
      std::apply(
          [&](auto&&... attacher) {
            ((attacher(nState, volume, gctx, position, direction, absMomentum,
                       charge)),
             ...);
          },
          m_candidateAttachers);
    }
    // Update internal candidates
    updateCandidates(nState, gctx, position, direction, absMomentum, charge);
  }
};

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
