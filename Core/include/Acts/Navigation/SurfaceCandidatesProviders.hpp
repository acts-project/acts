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
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationDelegateHelpers.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <tuple>

namespace Acts {
namespace Experimental {

/// Helper method to update the candidates (portals/surfaces),
/// this can be called for initial surface/portal estimation,
/// but also during the navigation to update the current list
/// of candidates.
///
/// @param gctx is the Geometry context of this call
/// @param nState [in] is the navigation state
/// @param candidates [out] are the candidates to be updated
///
/// @todo for surfaces skip the non-reached ones, while keep for portals
inline static void updateCandidates(
    const GeometryContext& gctx, const NavigationState& nState,
    NavigationState::SurfaceCandidates& candidates) {
  const auto& position = nState.position;
  const auto& direction = nState.direction;

  for (auto& c : candidates) {
    // Get the surface representation: either native surfcae of portal
    const Surface& sRep =
        (c.surface != nullptr) ? (*c.surface) : (c.portal->surface());

    // Get the intersection @todo make a templated intersector
    // TODO surface tolerance
    auto sIntersection = sRep.intersect(gctx, position, direction,
                                        c.boundaryCheck, s_onSurfaceTolerance);
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
      candidates.begin(), candidates.end(), [&](const auto& a, const auto& b) {
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
}

struct AllPortalsImpl : public ISurfaceCandidatesProvider {
  /// A ordered portal provider
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState is the navigation state to be updated
  ///
  /// @note that the intersections are ordered, such that the
  /// smallest intersection pathlength >= overstep tolerance is the lowest
  ///
  /// @return an ordered list of portal candidates
  inline void update(
      const GeometryContext& gctx, const NavigationState& nState,
      NavigationState::SurfaceCandidates& candidates) const override {
    auto currentVolume = nState.currentVolume;

    if (currentVolume == nullptr) {
      throw std::runtime_error("no detector volume set");
    }

    // Fill child volume portals
    for (const auto v : currentVolume->volumes()) {
      fillSurfaceCandidates(candidates, v->portals());
    }

    // Fill portals from the current volume
    fillSurfaceCandidates(candidates, currentVolume->portals());

    // Sort and update
    updateCandidates(gctx, nState, candidates);
  }
};

struct AllPortalsAndSurfacesImpl : public ISurfaceCandidatesProvider {
  /// An ordered list of portals and surfaces provider
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState is the navigation state to be updated
  ///
  /// @note that the intersections are ordered, such that the
  /// smallest intersection pathlength >= overstep tolerance is the lowest
  ///
  /// @return an ordered list of portal and surface candidates
  inline void update(
      const GeometryContext& gctx, const NavigationState& nState,
      NavigationState::SurfaceCandidates& candidates) const override {
    auto currentVolume = nState.currentVolume;

    if (currentVolume == nullptr) {
      throw std::runtime_error("no detector volume set");
    }

    // Fill child volume portals
    for (const auto v : currentVolume->volumes()) {
      fillSurfaceCandidates(candidates, v->portals());
    }

    // Fill surfaces and portals from the current volume
    fillSurfaceCandidates(candidates, currentVolume->portals());
    fillSurfaceCandidates(candidates, currentVolume->surfaces(),
                          nState.surfaceBoundaryCheck);

    // Update internal candidates
    updateCandidates(gctx, nState, candidates);
  }
};

/// Generate a provider for all portals
///
/// @return a connected navigationstate updator
inline static SurfaceCandidatesProvider tryAllPortals() {
  auto ap = std::make_unique<const AllPortalsImpl>();
  SurfaceCandidatesProvider nStateUpdator;
  nStateUpdator.connect<&AllPortalsImpl::update>(std::move(ap));
  return nStateUpdator;
}

/// Generate a provider for all portals and Surfacess
///
/// @note this is a try-and error navigation, not recommended for production
/// setup with many surfaces
///
/// @return a connected navigationstate updator
inline static SurfaceCandidatesProvider tryAllPortalsAndSurfaces() {
  auto aps = std::make_unique<const AllPortalsAndSurfacesImpl>();
  SurfaceCandidatesProvider nStateUpdator;
  nStateUpdator.connect<&AllPortalsAndSurfacesImpl::update>(std::move(aps));
  return nStateUpdator;
}

/// @brief This holds and extracts a collection of surfaces without much
/// checking, this could be e.g. support surfaces for layer structures,
/// e.g.
///
struct AdditionalSurfacesImpl : public ISurfaceCandidatesProvider {
  /// The volumes held by this collection
  std::vector<const Surface*> surfaces = {};

  /// Extract the sub volumes from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call (ignored)
  /// @param nState is the navigation state
  ///
  inline void update([[maybe_unused]] const GeometryContext& gctx,
                     NavigationState& nState) const override {
    SurfacesFiller::fill(nState, surfaces);
  }
};

/// @brief  An indexed surface implementation access
///
/// @tparam grid_type is the grid type used for this indexed lookup
template <typename grid_type>
using IndexedSurfacesImpl =
    IndexedUpdatorImpl<grid_type, IndexedSurfacesExtractor, SurfacesFiller>;

/// @brief An indexed surface implementation with portal access
///
///@tparam inexed_updator is the updator for the indexed surfaces
template <typename grid_type, template <typename> class indexed_updator>
using IndexedSurfacesAllPortalsImpl =
    ChainedUpdatorImpl<AllPortalsImpl, indexed_updator<grid_type>>;

}  // namespace Experimental
}  // namespace Acts
