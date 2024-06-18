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
#include "Acts/Navigation/MultiLayerNavigation.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Navigation/NavigationStateUpdaters.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <tuple>

namespace Acts::Experimental {

struct AllPortalsNavigation : public IInternalNavigation {
  /// A ordered portal provider
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState is the navigation state to be updated
  ///
  /// @note that the intersections are ordered, such that the
  /// smallest intersection pathlength >= overstep tolerance is the lowest
  ///
  /// @return an ordered list of portal candidates
  inline void update(const GeometryContext& gctx,
                     NavigationState& nState) const {
    if (nState.currentVolume == nullptr) {
      throw std::runtime_error(
          "AllPortalsNavigation: no detector volume set to navigation state.");
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
};

struct AllPortalsAndSurfacesNavigation : public IInternalNavigation {
  /// An ordered list of portals and surfaces provider
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState is the navigation state to be updated
  ///
  /// @note that the intersections are ordered, such that the
  /// smallest intersection pathlength >= overstep tolerance is the lowest
  ///
  /// @return an ordered list of portal and surface candidates
  inline void update(const GeometryContext& gctx,
                     NavigationState& nState) const {
    if (nState.currentDetector == nullptr) {
      throw std::runtime_error(
          "AllPortalsAndSurfacesNavigation: no detector volume set to "
          "navigation "
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
};

struct FrustumNavigation : public IInternalNavigation {
  /// The frustum navigation
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState is the navigation state to be updated
  ///
  /// @note that the intersections are ordered, such that the
  /// smallest intersection pathlength >= overstep tolerance is the lowest
  ///
  /// @return an ordered list of portal and surface candidates
  inline void update(const GeometryContext& gctx,
                     NavigationState& nState) const {
    if (nState.currentDetector == nullptr) {
      throw std::runtime_error(
          "FrustumNavigation: no detector volume set to navigation state.");
    }
    // Create the frustum if it's not set
    if (nState.frustum == nullptr) {
      nState.frustum = std::make_shared<Frustum3>(nState.position,
                                                  nState.direction, M_PI / 4);
    }
    // Check if we leave the frustum, reset candidates if so
    const auto& normals = nState.frustum->normals();
    const bool outside = std::any_of(
        normals.begin(), normals.end(), [&nState](const auto& normal) {
          return (nState.position - nState.frustum->origin()).dot(normal) >= 0;
        });
    if (outside) {
      nState.surfaceCandidates.clear();
      nState.frustum = std::make_shared<Frustum3>(nState.position,
                                                  nState.direction, M_PI / 4);
    }
    // A volume switch has happened, or we left the frustum, update list of
    // portals & surfaces
    if (nState.surfaceCandidates.empty()) {
      // Fill internal portals if existing
      auto topBoxCopy = nState.topBox;
      while (topBoxCopy) {
        if (topBoxCopy->intersect(*nState.frustum)) {
          if (topBoxCopy->hasEntity()) {
            const auto& portals = topBoxCopy->entity()->portals();
            const auto& surfaces = topBoxCopy->entity()->surfaces();
            PortalsFiller::fill(nState, portals);
            SurfacesFiller::fill(nState, surfaces);
            topBoxCopy = topBoxCopy->getSkip();
          } else {
            topBoxCopy = topBoxCopy->getLeftChild();
          }
        } else {
          topBoxCopy = topBoxCopy->getSkip();
        }
      }
    }
    // Assign the new volume
    const auto& portals = nState.currentVolume->portals();
    const auto& surfaces = nState.currentVolume->surfaces();
    PortalsFiller::fill(nState, portals);
    SurfacesFiller::fill(nState, surfaces);
    // Update internal candidates
    updateCandidates(gctx, nState);
  }
};

/// Generate a provider for all portals
///
/// @return a connected navigationstate updator
inline static InternalNavigationDelegate tryAllPortals() {
  auto ap = std::make_unique<const AllPortalsNavigation>();
  InternalNavigationDelegate nStateUpdater;
  nStateUpdater.connect<&AllPortalsNavigation::update>(std::move(ap));
  return nStateUpdater;
}

/// Generate a provider for all portals and Surfacess
///
/// @note this is a try-and error navigation, not recommended for production
/// setup with many surfaces
///
/// @return a connected navigationstate updator
inline static InternalNavigationDelegate tryAllPortalsAndSurfaces() {
  auto aps = std::make_unique<const AllPortalsAndSurfacesNavigation>();
  InternalNavigationDelegate nStateUpdater;
  nStateUpdater.connect<&AllPortalsAndSurfacesNavigation::update>(
      std::move(aps));
  return nStateUpdater;
}

/// Generate a provider for portals and surfaces inside a frustum
///
/// @return a connected navigationstate updator
inline static InternalNavigationDelegate tryFrustumPortalsAndSurfaces() {
  auto fr = std::make_unique<const FrustumNavigation>();
  InternalNavigationDelegate nStateUpdater;
  nStateUpdater.connect<&FrustumNavigation::update>(std::move(fr));
  return nStateUpdater;
}

/// @brief This holds and extracts a collection of surfaces without much
/// checking, this could be e.g. support surfaces for layer structures,
/// e.g.
///
struct AdditionalSurfacesNavigation : public IInternalNavigation {
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
/// @tparam grid_type is the grid type used for this indexed lookup
template <typename grid_type>
using IndexedSurfacesNavigation =
    IndexedGridNavigation<IInternalNavigation, grid_type,
                          IndexedSurfacesExtractor, SurfacesFiller>;

/// @brief  An indexed multi layer surface implementation access
///
/// @tparam grid_type is the grid type used for this indexed lookup
template <typename grid_type>
using MultiLayerSurfacesNavigation =
    MultiLayerNavigation<grid_type, PathGridSurfacesGenerator>;

/// @brief An indexed surface implementation with portal access
///
///@tparam inexed_updator is the updator for the indexed surfaces
template <typename grid_type, template <typename> class indexed_updator>
using IndexedSurfacesAllPortalsNavigation =
    ChainedNavigation<IInternalNavigation, AllPortalsNavigation,
                      indexed_updator<grid_type>>;

}  // namespace Acts::Experimental
