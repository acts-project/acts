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

/// A ordered portal and surface provider
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
inline static void portalAndSurfaceCandidates(
    NavigationState& nState, const DetectorVolume& volume,
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, ActsScalar absMomentum, ActsScalar charge) {
  // A volume switch has happened, update list of portals & surfaces
  if (nState.currentVolume != &volume) {
    // Assign the new volume
    nState.currentVolume = &volume;
    const auto& portals = volume.portals();
    const auto& surfaces = volume.surfaces();
    nState.surfaceCandidates.clear();
    nState.surfaceCandidate = nState.surfaceCandidates.end();
    // Reserve portals and surfaces size
    nState.surfaceCandidates.reserve(portals.size() + surfaces.size());
    // Add portal candidates
    for (const auto* p : portals) {
      nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
          ObjectIntersection<Surface>{}, nullptr, p, true});
    }
    // Add surface candidates
    for (const auto* s : surfaces) {
      nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
          ObjectIntersection<Surface>{}, s, nullptr,
          nState.surfaceBoundaryCheck});
    }
  }
  // Update internal candidates
  updateCandidates(nState, gctx, position, direction, absMomentum, charge);
}

/// Generate a provider for all portals
///
/// @return a connected navigationstate updator
inline static ManagedNavigationStateUpdator allPortals() {
  ManagedNavigationStateUpdator managedUpdator;
  NavigationStateUpdator nStateUpdator;
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
inline static ManagedNavigationStateUpdator allPortalsAndSurfaces() {
  ManagedNavigationStateUpdator managedUpdator;
  NavigationStateUpdator nStateUpdator;
  nStateUpdator.connect<&portalAndSurfaceCandidates>();
  managedUpdator.delegate = std::move(nStateUpdator);
  managedUpdator.implementation = nullptr;
  return managedUpdator;
}

class AllSurfacesAttacher : public IDelegateImpl {
 public:
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

class AllPortalsAttacher : public IDelegateImpl {
 public:
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

/// @brief  This is a templated grid surface attacher, it relies on
/// an indexed based grid lookup and a range interation
///
/// @tparam grid_type the grid type
///
template <typename grid_type>
class GridSurfaceAttacher : public IDelegateImpl {
 public:
  /// @brief The grid
  grid_type grid;

  std::array<BinningValue, grid_type::DIM> casts;

  /// @brief  Constructor for a grid based surface attacher
  /// @param igrid the grid that is moved into this attacher
  /// @param icasts is the cast values array
  GridSurfaceAttacher(grid_type&& igrid, const std::array<BinningValue, grid_type::DIM>& icasts) 
  : grid(std::move(igrid)), casts(icasts) {}


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
                  const Vector3& position,
                  [[maybe_unused]] const Vector3& direction,
                  [[maybe_unused]] ActsScalar absMomentum,
                  [[maybe_unused]] ActsScalar charge) const {
    // Get the surfaces containerd in this volume
    const auto& surfaces = volume.surfaces();
    // Retrieve the grid indices and fill the surface candidates
    const auto& indices = grid.atPosition(castPosition(position));
    for (const auto& idx : indices) {
      nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
          ObjectIntersection<Surface>{}, surfaces[idx], nullptr,
          nState.surfaceBoundaryCheck});
    }
  }

  /// Cast into a lookup position
  ///
  /// @param position is the position of the update call
  std::array<ActsScalar, grid_type::DIM> castPosition(
      const Vector3& position) const {
    std::array<ActsScalar, grid_type::DIM> casted;
    fillCasts(position, casted,
              std::make_integer_sequence<std::size_t, grid_type::DIM>{});
    return casted;
  }

 private:
  /// Unroll cast loop
  /// @param position is the position of the update call
  /// @param a is the array to be filled
  template <typename Array, std::size_t... idx>
  void fillCasts(const Vector3& position, Array& a,
                 std::index_sequence<idx...>) const {
    ((a[idx] = VectorHelpers::cast(position, casts[idx])), ...);
  }
};

// This struct allows to combine the portals from the volume itself
// with  attachers
template <typename... candidate_attachers_t>
class NavigationStateUpdatorImpl : public IDelegateImpl {
 public:
  std::tuple<candidate_attachers_t...> candidateAttachers;

  NavigationStateUpdatorImpl(const std::tuple<candidate_attachers_t...>& attachers)
      : candidateAttachers(attachers) {}

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
          candidateAttachers);
    }
    // Update internal candidates
    updateCandidates(nState, gctx, position, direction, absMomentum, charge);
  }
};

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
