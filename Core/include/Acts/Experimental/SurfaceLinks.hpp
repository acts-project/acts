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
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Enumerate.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <functional>
#include <memory>
#include <set>
#include <vector>

namespace Acts {

/// Helper function to get the surface candidates from a list, can and should
/// be reuesed by all surface finders as it guaraqntees consistent path,
/// overstep and alternative solutions handling
///
/// @tparam surface_container_t is an iterable surface container
///
/// @param gctx is the current geometry conbtext
/// @param surfaces is the iterable surface container to draw candidates from
/// @param position is the position at the query
/// @param direction is the direction at the query
/// @param bCheck is the BoundaryCheck
/// @param pathRange is the allowed path range for the intersection
/// @param guessedNumber is the guessed number of intersections for
/// reserving
///
/// @note onSurface solutions are ranked last
template <typename surface_container_t>
std::vector<SurfaceIntersection> surfaceCandidates(
    const GeometryContext& gctx, const surface_container_t& surfaces,
    const Vector3& position, const Vector3& direction,
    const BoundaryCheck& bCheck, const std::array<ActsScalar, 2>& pathRange,
    size_t guessedNumber) {
  // Return object and reserving
  std::vector<SurfaceIntersection> sIntersections;
  sIntersections.reserve(guessedNumber);
  for (const auto& s : surfaces) {
    auto sIntersection = s->intersect(gctx, position, direction, bCheck);
    // Swap in case of two solutions and conflict with overstep tolerance
    if (sIntersection.intersection.pathLength + s_onSurfaceTolerance <
            pathRange[0] and
        sIntersection.alternative.pathLength + s_onSurfaceTolerance >
            pathRange[0] and
        sIntersection.alternative.status >= Intersection3D::Status::reachable) {
      sIntersection.swapSolutions();
    }
    // The intersection needs to be valid and within range
    if (sIntersection.intersection.status >=
            Intersection3D::Status::reachable and
        sIntersection.intersection.pathLength + s_onSurfaceTolerance >=
            pathRange[0] and
        sIntersection.intersection.pathLength - s_onSurfaceTolerance <
            pathRange[1]) {
      sIntersections.push_back(sIntersection);
    }
  }
  // Sort and return
  std::sort(sIntersections.begin(), sIntersections.end());
  return sIntersections;
}

struct AllSurfaces {
  /// A guess for the number of candidates (for vector reserve)
  size_t guessedNumberOfCandidates = 25;

  /// Call operator for this struct to provide the list of surfaces
  ///
  /// @param gctx the current geometry context
  /// @param volume the currenct volume
  /// @param position the current position for intersection start
  /// @param direction the current momentum for intersecton start
  /// @param bCheck the boundary check for this search
  /// @param pathRange the allowed path range for this search
  /// @note provideAll is (obvisouly) ignored in this context
  ///
  /// @note the onSurface solution will be given first
  ///
  /// @return the surface intersections (ordered)
  std::vector<SurfaceIntersection> operator()(
      const GeometryContext& gctx, const DetectorVolume& volume,
      const Vector3& position, const Vector3& direction,
      const BoundaryCheck& bCheck, const std::array<ActsScalar, 2>& pathRange,
      bool) const {
    return surfaceCandidates(gctx, volume.surfaces(), position, direction,
                             bCheck, pathRange, guessedNumberOfCandidates);
  }
};

struct SingleSurface {};

template <typename grid_t>
struct SingleGridSurfaces {
  /// The grid with access indices
  grid_t accessGrid;

  /// The parameter casts from local into grid point definition
  std::vector<BinningValue> parameterCasts = {};

  /// The transform into grid local frame
  Transform3 toLocal = Transform3::Identity();

  /// Constructor
  ///
  /// @param accessGrid_ the access grid for indices
  /// @param parameterCasts_ the binning value list for casting parameters
  ///        into the grid point definition
  /// @param toLocal_ the transform to local for the grid access
  ///
  SingleGridSurfaces(grid_t&& accessGrid_,
                     const std::vector<BinningValue>& parameterCasts_,
                     const Transform3& toLocal_ = Transform3::Identity())
      : accessGrid(std::move(accessGrid_)),
        parameterCasts(parameterCasts_),
        toLocal(toLocal_) {}

  /// Call operator for this struct to provide the list of surfaces
  ///
  /// @param gctx the current geometry context
  /// @param volume the currenct volume
  /// @param position the current position for intersection start
  /// @param direction the current momentum for intersecton start
  /// @param bCheck the boundary check for this search
  /// @param pathRange the allowed path range for this search
  /// @note provideAll returns all surfaces from the volume for intersection
  ///
  /// @note the onSurface solution will be given first
  ///
  /// @return the surface intersections (ordered)
  std::vector<SurfaceIntersection> operator()(
      const GeometryContext& gctx, const DetectorVolume& volume,
      const Vector3& position, const Vector3& direction,
      const BoundaryCheck& bCheck, const std::array<ActsScalar, 2>& pathRange,
      bool provideAll) const {
    // Validation fallback to provide all candidates
    if (provideAll) {
      return surfaceCandidates(gctx, volume.surfaces(), position, direction,
                               bCheck, pathRange, 8u);
    }

    // Bring the position into local frame & cast into the grid point definition
    Vector3 posInFrame = toLocal * position;
    typename decltype(accessGrid)::point_t castedPosition;
    for (auto [i, castValue] : enumerate(parameterCasts)) {
      castedPosition[i] = VectorHelpers::cast(posInFrame, castValue);
    }

    // Set of unique candidates (tried set/unordered set but was slower)
    std::vector<const Surface*> localCandidates;
    localCandidates.reserve(10u);
    const auto& surfaces = volume.surfaces();

    // Get central Index and neighbors around and fill into the unique set
    auto centralIndex = accessGrid.localBinsFromPosition(castedPosition);
    auto neighborIndices = accessGrid.neighborHoodIndices(centralIndex, 1u);
    auto surfaceIndex = accessGrid.atLocalBins(centralIndex);
    if (surfaceIndex >= 0 and
        static_cast<size_t>(surfaceIndex) < surfaces.size()) {
      localCandidates.push_back(surfaces[surfaceIndex]);
    }
    for (const auto& neighborIndex : neighborIndices) {
      surfaceIndex = accessGrid.at(neighborIndex);
      if (surfaceIndex >= 0 and
          static_cast<size_t>(surfaceIndex) < surfaces.size()) {
        localCandidates.push_back(surfaces[surfaceIndex]);
      }
    }
    
    std::sort(localCandidates.begin(), localCandidates.end()); 
    auto last = std::unique(localCandidates.begin(), localCandidates.end());
    localCandidates.erase(last, localCandidates.end());

    // Pre-intersect the unique candidates and return to the portal
    return surfaceCandidates(gctx, localCandidates, position, direction, bCheck,
                             pathRange, 8u);
  }
};

}  // namespace Acts
